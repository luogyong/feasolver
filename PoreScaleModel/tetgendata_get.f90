module tetgendata
    use MeshDs,only:lowcase,strtoint,skipcomment
    use penf
    use vtk_fortran, only : vtk_file
    use VTK_CELLTYPE
    use quicksort
    implicit none
    
    public::tetgendata_tydef
    
    private
    
    
    integer::order=1
    character(512)::filepath
 
    logical::isfirstcall=.true.

    real(8)::eps=1.d-6
    
    !type ar2d_tydef3
    !    integer::nnode=0
    !    integer,allocatable::node(:)
    !endtype   
    
    type::tetgen_node_tydef
        integer::na=0,marker=0,surfin=0,a1=0,uid=0      
        !for vnode_tg, marker>0 inside the con and is the nodeid for tecplot out
        !if(node(i).uid<i) then the node is a duplicated node, which is the copy of node(uid) 
        !for tetgen node ,marker=-1,means the node is out the model box;and >0 ,nodal id
        real(8)::x(3)
        real(8),allocatable::at(:)
    endtype
    
    type adjlist_tydef
        integer::nnode=0,nedge=0,nface=0,nelt=0
        integer,allocatable::node(:),edge(:),face(:),elt(:)
        integer,allocatable::esubid(:),fsubid(:),eltsubid(:) 
    contains
        procedure::push=>tet_adjlist_push_element
    endtype
    
 
    
    type tetgen_element_tydef
        integer::nnode=0,marker=0,nedge=0,isclose=0,nvcnode=0
        !isclose,for vedge,vface and vcell =0 open/infinite. >0 close/finite,=-1,some vertex are outside the modelbox
        !for vedge,marker=0 outside container;marker/=0,inside the con;=2,clipped;=3,clipped and coplan of the con; =4,clipped and not coplan;
        !for vface,isclose>0, close and (all nodes) inside the con,and it will be output to the tec file; =-1,close and some node outside the cont; =-2 close and all node outside the con;
        !for vface,iclose=0, open
        !for vcell, isclose>0,close and (all nodes) inside the con,and it will output to the the tec file.
        !for vface, marker==0,not clipped face; =3,clipped on only the generated edge is on one boundary; =4,clopped and the generated egde is on different boundary 
        
        integer,allocatable::node(:)
        integer::cell(2)=-1 !for voronoi face, they are cell sharing the face, for vcell, cell(1) is the corresponding node id (node(cell(1))) for the cell
        real(8),allocatable::V(:)  
        !for edge ray, its direction, for face, its unit normal
        !for cell, its center.
        integer,allocatable::edge(:) !for voronoi face only
        !integer::celoc=0,surfin=-1 !for vface, celoc=clipped edge id; surfin is the boudary id where the vertex  of the clipped edge laying. 
        integer,allocatable::vcnode(:) !nodes of the vcell. the node(:) is actually the vface ids. 
        integer,allocatable::vtu_vcface(:) !data for vtu output,the format is according the vtu polyhedron: [nf,[[nnode,ni for j=1,nnode] for i=1,nf]]
        real(8)::property  !for faces,the area; for edges,the length.
        
        logical::isbigballcell=.false.
           
    end type  
    
 
    type tetgendata_tydef
        integer::nnode,nelt,nface,nedge,nneigh,isremoveduplicates=1
        integer::nvnode,nvcell,nvface,nvedge
        type(tetgen_node_tydef),allocatable::node(:),vnode(:)
        type(tetgen_element_tydef),allocatable::edge(:),face(:),elt(:), vface(:),&
                vedge(:),vcell(:),neigh(:)
        integer,allocatable::t2e(:,:),t2f(:,:),f2e(:,:),t2t(:,:) !t2t为单元的邻接单元,相当于neigh
        type(adjlist_tydef),allocatable::nadjlist(:),eadjlist(:),vnadjlist(:) !nodal,and edge list
        character(512)::file        
    contains
        procedure::readin=>read_tetgen_file
        procedure::setadjlist=>tetgen_setup_adjlist
        procedure::out_mesh_vtu=>output_tetgen_mesh_vtu
        procedure::out_voro_vtu=>output_tetgen_voro_mesh_vtu
        procedure,private::vfnode_setup,vcnode_setup,vtu_vcface_setup
        
    endtype
    

    
    INTERFACE ENLARGE_AR
        MODULE PROCEDURE I_ENLARGE_AR,NODE_ENLARGE_AR,ELEMENT_ENLARGE_AR                        
    END INTERFACE 
    
    contains
    
     
        
    subroutine output_tetgen_voro_mesh_vtu(self,vtkformat,box)
    !only output the closed cell info.
        implicit none
        class(tetgendata_tydef)::self
        character(len=*)::vtkformat ![ascii,raw,binary]
        real(8)::box(6)
        !logical,optional::isremoveduplicates
        integer::i,j,k,nelt1
        real(R8P),allocatable::x1(:)
        integer(I4P),ALLOCATABLE::IX1(:)
        type(vtk_file)                :: a_vtk_file                                                    !< A VTK file.
        integer(I1P), allocatable, dimension(:)    :: cell_type
        integer(I4P), allocatable, dimension(:)    :: offset
        integer(I4P), allocatable, dimension(:)    :: connect
        integer(I4P), allocatable, dimension(:)    :: face
        integer(I4P), allocatable, dimension(:)    :: faceoffset
        integer(I4P)                               :: error
        INTEGER :: CSTAT, ESTAT
        CHARACTER(100) :: CMSG
        integer::nptcle1,n1,n2,n3,n4
        character(8)::ch1     
        logical::tof1
        
        
        call self.vfnode_setup()
        call self.vcnode_setup()
        call self.vtu_vcface_setup()
        
        nelt1=0;n1=0;n2=0
        do i=1,self.nvcell  
            if(self.vcell(i).isclose<1) then
                if(self.vcell(i).isclose==0) self.node(i).marker=-2 !mark it
                cycle
            endif
            
            !如果voro的节点在box之外,也不输出,
            if(any(self.vnode(self.vcell(i).vcnode).marker<0)) then
                if(self.vcell(i).isclose==1) then
                    self.vcell(i).isclose=-1 
                    self.node(i).marker=-1 !close but outof the model box, mark it as =-1
                    !node与vcell一一对应
                else
                    if(self.vcell(i).isclose==0) self.node(i).marker=-2 !unclose and outside the box mark it
                endif
                               
                cycle
            endif
            
            nelt1=nelt1+1
            n1=n1+self.vcell(i).nvcnode
            n2=n2+size(self.vcell(i).vtu_vcface)
        enddo

        
        
        allocate(cell_type(nelt1),offset(nelt1),connect(n1),face(n2),faceoffset(nelt1))
        !cell_type=3;offset=2;
        n1=0;nelt1=0;n3=0
        do i=1,self.nvcell
            if(self.vcell(i).isclose<1) cycle
            nelt1=nelt1+1
            n2=self.vcell(i).nvcnode
            if(self.isremoveduplicates/=0) then
                connect(n1+1:n1+n2)=self.vnode(self.vnode(self.vcell(i).vcnode(1:n2)).uid).marker-1 !vtu is 0-based
            else                
                connect(n1+1:n1+n2)=self.vcell(i).vcnode(1:n2)-1 !vtu is 0-based，duplicates are allowed
            endif
            n1=n1+n2
            cell_type(nelt1)=VTK_POLYHEDRON
            !face
            n4=size(self.vcell(i).vtu_vcface)
            face(n3+1:n3+n4)=self.vcell(i).vtu_vcface(1:n4)
            n3=n3+n4
                        
            if(nelt1==1) then
                offset(nelt1)=n2
                faceoffset(nelt1)=n4
            else
                offset(nelt1)=offset(nelt1-1)+n2
                faceoffset(nelt1)=faceoffset(nelt1-1)+n4
            endif
            
        enddo

        ! ascii
        if(self.isremoveduplicates/=0) then
            ix1=pack([1:self.nvnode],self.vnode(1:self.nvnode).marker>0)
        else
            ix1=pack([1:self.nvnode],.true.)
        endif
        n2=size(ix1,dim=1)
        
        error = a_vtk_file%initialize(format=trim(vtkformat), filename=trim(self.file)//'_voro.vtu', mesh_topology='UnstructuredGrid')
        error = a_vtk_file%xml_writer%write_piece(np=n2, nc=nelt1)
        error = a_vtk_file%xml_writer%write_geo(np=n2, nc=nelt1, x=self.vnode(ix1).x(1), y=self.vnode(ix1).x(2), z=self.vnode(ix1).x(3))
        error = a_vtk_file%xml_writer%write_connectivity(nc=nelt1, connectivity=connect, offset=offset, cell_type=cell_type,face=face,faceoffset=faceoffset)
    
        error = a_vtk_file%xml_writer%write_dataarray(location='node', action='open') 
        do i=1,self.vnode(1).na
            if(i==1) then
                !注意,在此假定(默认)第一个属性值为半径的平方.
                allocate(x1(n2))
                do j=1,n2
                    x1(j)=sqrt(self.vnode(ix1(j)).at(i))
                enddo            
                error = a_vtk_file%xml_writer%write_dataarray(data_name='Rt(正交球半径)', x=x1)
                deallocate(x1)
            else
                write(ch1,'(A)') i
                allocate(x1(n2))
                do j=1,n2
                    x1(j)=self.vnode(ix1(j)).at(i)
                enddo 
                error = a_vtk_file%xml_writer%write_dataarray(data_name='at'//trim(adjustl(ch1)), x=x1)
                deallocate(x1)
            endif
        enddo
        error = a_vtk_file%xml_writer%write_dataarray(data_name='Marker', x=self.vnode(ix1).marker)
    
        error = a_vtk_file%xml_writer%write_dataarray(location='node', action='close')
             
        error = a_vtk_file%xml_writer%write_dataarray(location='cell', action='open')
        
        if(allocated(ix1)) deallocate(ix1)
        allocate(ix1(nelt1))    
        n1=0
        do i=1,self.nvcell
            if(self.vcell(i).isclose<1) cycle
            n1=n1+1
            ix1(n1)=self.vcell(i).marker
        enddo
        error = a_vtk_file%xml_writer%write_dataarray(data_name='Marker', x=ix1) 
        deallocate(ix1)
        !error = a_vtk_file%xml_writer%write_dataarray(data_name='Rh_t(喉水力半径)', x=self.elt(1:self.nelt).rht)
        !error = a_vtk_file%xml_writer%write_dataarray(data_name='V_t(喉体积)', x=self.elt(1:self.nelt).pv)
        !error = a_vtk_file%xml_writer%write_dataarray(data_name='Ap_t(单元颗粒表面积)', x=self.elt(1:self.nelt).pa)
        !error = a_vtk_file%xml_writer%write_dataarray(data_name='Req_t(喉面积等效半径)', x=self.elt(1:self.nelt).req)
        !error = a_vtk_file%xml_writer%write_dataarray(data_name='Reff_t2(喉有效半径)', x=self.elt(1:self.nelt).reff)
        !error = a_vtk_file%xml_writer%write_dataarray(data_name='Ve_t(喉总体积)', x=self.elt(1:self.nelt).v)
        !error = a_vtk_file%xml_writer%write_dataarray(data_name='e_t(喉单元孔隙比)', x=self.elt(1:self.nelt).e)
        !
        error = a_vtk_file%xml_writer%write_dataarray(location='cell', action='close')
        error = a_vtk_file%xml_writer%write_piece()
        error = a_vtk_file%finalize()
    
        !因为生产xml文件包含中文,convert to utf8,方便paraview读入,也可以外部手动转
        CALL EXECUTE_COMMAND_LINE ('powershell -command "Get-Content '//trim(self.file)//'_voro.vtu'//' | Set-Content -Encoding utf8 '//trim(self.file)//'_voro-utf8.vtu'//'"', EXITSTAT=ESTAT,CMDSTAT=CSTAT, CMDMSG=CMSG)
        
        IF (CSTAT > 0) THEN
            PRINT *, "Tetgen Command execution failed with error ", TRIM(CMSG)
            pause
        ELSE IF (CSTAT < 0) THEN
            PRINT *, "Tetgen Command execution not supported"
            pause
        !ELSE
        !    PRINT *, "Command completed with status ", ESTAT
        END IF
        !delete the origin file
        CALL EXECUTE_COMMAND_LINE ('del "'//trim(self.file)//'_voro.vtu"', EXITSTAT=ESTAT,CMDSTAT=CSTAT, CMDMSG=CMSG)
        
        IF (CSTAT > 0) THEN
            PRINT *, "Tetgen Command execution failed with error ", TRIM(CMSG)
            pause
        ELSE IF (CSTAT < 0) THEN
            PRINT *, "Tetgen Command execution not supported"
            pause
        !ELSE
        !    PRINT *, "Command completed with status ", ESTAT
        END IF
    
        !! raw
        !error = a_vtk_file%initialize(format='raw', filename='XML_UNST-raw.vtu', mesh_topology='UnstructuredGrid')
        !call write_data
        !error = a_vtk_file%finalize()
        !! binary
        !error = a_vtk_file%initialize(format='binary', filename='XML_UNST-binary.vtu', mesh_topology='UnstructuredGrid')
        !call write_data
        !error = a_vtk_file%finalize()

        !test_passed = .true. ! nothing to test yet

        !print "(A,L1)", new_line('a')//'Are all tests passed? ', all(test_passed)
    
        deallocate(cell_type,offset,connect,face,faceoffset)        
        
       
        
    endsubroutine
    
    subroutine output_tetgen_mesh_vtu(self,vtkformat,isoutbigball,nptcle)
        implicit none
        class(tetgendata_tydef)::self
        character(len=*)::vtkformat ![ascii,raw,binary]
        integer,optional::isoutbigball
        integer,optional::nptcle
        integer::i,j,k,nelt1
        real(R8P),allocatable::x1(:)
        integer(I4P),ALLOCATABLE::IX1(:)
        type(vtk_file)                :: a_vtk_file                                                    !< A VTK file.
        integer(I1P), allocatable, dimension(:)    :: cell_type
        integer(I4P), allocatable, dimension(:)    :: offset
        integer(I4P), allocatable, dimension(:)    :: connect
        integer(I4P)                               :: error
        INTEGER :: CSTAT, ESTAT
        CHARACTER(100) :: CMSG
        integer::isoutbb1=0
        integer::nptcle1,n1
        character(8)::ch1
        
        if(present(nptcle)) then
            nptcle1=nptcle
        else
            nptcle1=self.nnode    
        endif
        
        isoutbb1=0
        if(present(isoutbigball)) then
            isoutbb1=isoutbigball
        endif
        
        nelt1=self.nelt
        if(isoutbb1==0) then
            !含最后6个节点的单元不输出
            do i=1,self.nelt
                if(any(self.elt(i).node>nptcle1)) then
                    self.elt(i).isbigballcell=.true.
                    nelt1=nelt1-1
                else
                    self.elt(i).isbigballcell=.false.
                endif            
            enddo
        else
            nptcle1=self.nnode
        endif
        
    
    
    
    allocate(cell_type(nelt1),offset(nelt1),connect(4*nelt1))
    !cell_type=3;offset=2;
    n1=0
    do i=1,self.nelt
        if(self.elt(i).isbigballcell) cycle
        n1=n1+1
        connect(4*(n1-1)+1:4*n1)=self.elt(i).node(1:self.elt(i).nnode)-1 !vtu is 0-based
        if(n1==1) then
            offset(n1)=self.elt(i).nnode
        else
            offset(n1)=offset(n1-1)+self.elt(i).nnode
        endif
        cell_type(n1)=10
    enddo
    !offset(n1+1)=offset(n1-1)+self.elt(n1).nnode
    ! ascii
        
    error = a_vtk_file%initialize(format=trim(vtkformat), filename=trim(self.file)//'_tet.vtu', mesh_topology='UnstructuredGrid')
    error = a_vtk_file%xml_writer%write_piece(np=nptcle1, nc=nelt1)
    error = a_vtk_file%xml_writer%write_geo(np=nptcle1, nc=nelt1, x=self.node(1:nptcle1).x(1), y=self.node(1:nptcle1).x(2), z=self.node(1:nptcle1).x(3))
    error = a_vtk_file%xml_writer%write_connectivity(nc=nptcle1, connectivity=connect, offset=offset, cell_type=cell_type)
    
    error = a_vtk_file%xml_writer%write_dataarray(location='node', action='open') 
    do i=1,self.node(1).na
        if(i==1) then
            !注意,在此假定(默认)第一个属性值为半径的平方.
            allocate(x1(nptcle1))
            do j=1,nptcle1
                x1(j)=sqrt(self.node(j).at(i))
            enddo            
            error = a_vtk_file%xml_writer%write_dataarray(data_name='Rptcle(颗粒半径)', x=x1)
            deallocate(x1)
        else
            write(ch1,'(A)') i
            allocate(x1(nptcle1))
            do j=1,nptcle1
                x1(j)=self.node(j).at(i)
            enddo 
            error = a_vtk_file%xml_writer%write_dataarray(data_name='at'//trim(adjustl(ch1)), x=x1)
            deallocate(x1)
        endif
    enddo
    error = a_vtk_file%xml_writer%write_dataarray(data_name='Marker', x=self.node(1:nptcle1).marker)
    
    error = a_vtk_file%xml_writer%write_dataarray(location='node', action='close')
             
    error = a_vtk_file%xml_writer%write_dataarray(location='cell', action='open')
    
    allocate(ix1(nelt1))    
    n1=0
    do i=1,self.nelt
        if(self.elt(i).isbigballcell) cycle
        n1=n1+1
        ix1(n1)=self.elt(i).marker
    enddo
    error = a_vtk_file%xml_writer%write_dataarray(data_name='Marker', x=ix1) 
    deallocate(ix1)
    !error = a_vtk_file%xml_writer%write_dataarray(data_name='Rh_t(喉水力半径)', x=self.elt(1:self.nelt).rht)
    !error = a_vtk_file%xml_writer%write_dataarray(data_name='V_t(喉体积)', x=self.elt(1:self.nelt).pv)
    !error = a_vtk_file%xml_writer%write_dataarray(data_name='Ap_t(单元颗粒表面积)', x=self.elt(1:self.nelt).pa)
    !error = a_vtk_file%xml_writer%write_dataarray(data_name='Req_t(喉面积等效半径)', x=self.elt(1:self.nelt).req)
    !error = a_vtk_file%xml_writer%write_dataarray(data_name='Reff_t2(喉有效半径)', x=self.elt(1:self.nelt).reff)
    !error = a_vtk_file%xml_writer%write_dataarray(data_name='Ve_t(喉总体积)', x=self.elt(1:self.nelt).v)
    !error = a_vtk_file%xml_writer%write_dataarray(data_name='e_t(喉单元孔隙比)', x=self.elt(1:self.nelt).e)
    !
    error = a_vtk_file%xml_writer%write_dataarray(location='cell', action='close')
    error = a_vtk_file%xml_writer%write_piece()
    error = a_vtk_file%finalize()
    
    !因为生产xml文件包含中文,convert to utf8,方便paraview读入,也可以外部手动转
    CALL EXECUTE_COMMAND_LINE ('powershell -command "Get-Content '//trim(self.file)//'_tet.vtu'//' | Set-Content -Encoding utf8 '//trim(self.file)//'_tet-utf8.vtu'//'"', EXITSTAT=ESTAT,CMDSTAT=CSTAT, CMDMSG=CMSG)
        
    IF (CSTAT > 0) THEN
        PRINT *, "Tetgen Command execution failed with error ", TRIM(CMSG)
        pause
    ELSE IF (CSTAT < 0) THEN
        PRINT *, "Tetgen Command execution not supported"
        pause
    !ELSE
    !    PRINT *, "Command completed with status ", ESTAT
    END IF
    !delete the origin file
    CALL EXECUTE_COMMAND_LINE ('del "'//trim(self.file)//'_tet.vtu"', EXITSTAT=ESTAT,CMDSTAT=CSTAT, CMDMSG=CMSG)
        
    IF (CSTAT > 0) THEN
        PRINT *, "Tetgen Command execution failed with error ", TRIM(CMSG)
        pause
    ELSE IF (CSTAT < 0) THEN
        PRINT *, "Tetgen Command execution not supported"
        pause
    !ELSE
    !    PRINT *, "Command completed with status ", ESTAT
    END IF
    
    !! raw
    !error = a_vtk_file%initialize(format='raw', filename='XML_UNST-raw.vtu', mesh_topology='UnstructuredGrid')
    !call write_data
    !error = a_vtk_file%finalize()
    !! binary
    !error = a_vtk_file%initialize(format='binary', filename='XML_UNST-binary.vtu', mesh_topology='UnstructuredGrid')
    !call write_data
    !error = a_vtk_file%finalize()

    !test_passed = .true. ! nothing to test yet

    !print "(A,L1)", new_line('a')//'Are all tests passed? ', all(test_passed)
    
    deallocate(cell_type,offset,connect)        
        
        
        
    endsubroutine
    
    
    subroutine tetgen_setup_adjlist(self)
        implicit none
        class(tetgendata_tydef)::self
        integer::i,j,k,n1
        
        allocate(self.nadjlist(self.nnode),self.eadjlist(self.nedge),self.vnadjlist(self.nvnode))
        
        associate(alist=>self.nadjlist,elist=>self.eadjlist,edge=>self.edge,f2e=>self.f2e,&
            t2e=>self.t2e,face=>self.face,elt=>self.elt)
            do i=1,self.nedge            
                do j=1,2
                    call alist(edge(i).node(j)).push(0,edge(i).node(mod(j,2)+1),0)
                    call alist(edge(i).node(j)).push(1,i,j)
                enddo
            end do
            do i=1,self.nface
                do j=1,3
                     call alist(face(i).node(j)).push(2,i,j)
                     call elist(f2e(j,i)).push(2,i,j)
                enddo
            enddo
            do i=1,self.nelt
                do j=1,6
                     if(j<5) call alist(elt(i).node(j)).push(3,i,j)
                     call elist(t2e(j,i)).push(3,i,j)
                enddo
            enddo
        end associate
        associate(alist=>self.vnadjlist,edge=>self.vedge)
            do i=1,self.nvedge
                if(any(edge(i).node(1:2)<1)) cycle !infinite lines are skipped.
                do j=1,2                    
                    call alist(edge(i).node(j)).push(0,edge(i).node(mod(j,2)+1),0)
                    call alist(edge(i).node(j)).push(1,i,j)
                enddo
            end do
            !do i=1,self.nface
            !    do j=1,3
            !         call alist(face(i).node(j)).push(2,i,j)
            !         call elist(f2e(j,i)).push(2,i,j)
            !    enddo
            !enddo
            !do i=1,self.nelt
            !    do j=1,6
            !         if(j<5) call alist(elt(i).node(j)).push(3,i,j)
            !         call elist(t2e(j,i)).push(3,i,j)
            !    enddo
            !enddo
        end associate        
        
    endsubroutine
  
    subroutine tet_adjlist_push_element(self,itype,ielt,isub)
        implicit none
        class(adjlist_tydef)::self
        integer,intent(in)::itype,ielt,isub
        
        select case(itype)
        case(0) !node
            call push(self.node,self.nnode,ielt)            
        case(1) !edge
            call push(self.edge,self.nedge,ielt,self.esubid,isub)  
        case(2) !face
            call push(self.face,self.nface,ielt,self.fsubid,isub) 
        case(3) !tet
            call push(self.elt,self.nelt,ielt,self.eltsubid,isub) 
        case default
            error stop 'no such type. sub=tet_adjlist_push_element'
        endselect
    contains
        subroutine push(ia,nelt,ielt,ib,subid)
            implicit none
            integer,allocatable::ia(:)
            integer,allocatable,optional::ib(:)
            integer::nelt
            integer,intent(in)::ielt
            integer,optional::subid
            
            if(allocated(ia)) then
                if(any(ia-ielt==0)) return
            else
                allocate(ia(5))
                if(present(ib)) allocate(ib(5))
            endif
            nelt=nelt+1
            if(nelt>size(ia,dim=1)) then
                call enlarge_ar(ia,5)
                if(present(ib)) call enlarge_ar(ib,5)
            else
                ia(nelt)=ielt
                if(present(ib)) ib(nelt)=subid
            endif
        end
        
    endsubroutine

    
    subroutine read_tetgen_file(self,unit,fext,box)
        use dflib
        USE IFPORT
        implicit none
        class(tetgendata_tydef)::self
        integer,intent(in),optional::unit
        character(len=*),optional,intent(in)::fext(:)
        real(8),intent(in),optional::box(:)
        CHARACTER(3)        drive
	    CHARACTER(512)      dir
	    CHARACTER(512)      name,file1,FILEPATH
	    CHARACTER(16)      ext
        CHARACTER(len=:),allocatable::ext1(:)
	    integer(4)::length,msg
        logical::isexist
        integer::unit1,i,j,hasread
        
        if(present(unit)) then
            inquire(UNIT,name=file1)
        else
            file1=self.file
        endif
		length = SPLITPATHQQ(file1, drive, dir, name, ext)
		!tec_title=trim(name)
        FILEPATH=trim(drive)//trim(dir)
        msg = CHDIR(FILEPATH)
        FILEPATH=trim(drive)//trim(dir)//trim(name)
        
        if(trim(adjustl(ext))=='.cell') then
            FILEPATH=filepath(:len_trim(filepath)-2)
        endif

        if(present(unit)) close(unit)
        
        if(present(fext)) then
            ext1=fext
        else
            allocate(character(16)::ext1(12))
            ext1(1:12)=['node','ele','face','edge','neigh','t2e','t2f','f2e','v.node','v.edge','v.face','v.cell']
        endif
        
        
        do i=1,size(ext1)
            file1=trim(filepath)//'.'//trim(adjustl(ext1(i)))
            inquire(file=file1,exist=isexist)
            if(isexist) then
                unit1=10
                open(unit=unit1,file=file1,status='old')
                hasread=1
                select case(trim(adjustl(ext1(i))))
                case('node')
                    call read_tetgen_node(unit1,trim(adjustl(ext1(i))),self.node,self.nnode)
                    !natr=node(1).na
                case('ele')
                    call read_tetgen_element(unit1,trim(adjustl(ext1(i))),self.elt,self.nelt)
                    !ntetnode=element(1).nnode
                case('neigh')
                    call read_tetgen_element(unit1,trim(adjustl(ext1(i))),self.neigh,self.nneigh)                    
                case('face')
                    call read_tetgen_element(unit1,trim(adjustl(ext1(i))),self.face,self.nface)
                    !cal length
                    do j=1,self.nface
                        call triangle_area_3d (reshape([self.node(self.face(j).node(1)).x, &
                            self.node(self.face(j).node(2)).x,&
                            self.node(self.face(j).node(3)).x],[3,3]), self.face(j).property )                        
                    enddo
                case('edge')
                    call read_tetgen_element(unit1,trim(adjustl(ext1(i))),self.edge,self.nedge)
                    !cal length
                    do j=1,self.nedge
                        self.edge(j).property=norm2(self.node(self.edge(j).node(1)).x-self.node(self.edge(j).node(2)).x)
                    enddo
                case('v.node')
                    call read_tetgen_node(unit1,trim(adjustl(ext1(i))),self.vnode,self.nvnode)
                    call find_duplicated_node(self.vnode,box)
                case('v.face')
                    call read_tetgen_element(unit1,trim(adjustl(ext1(i))),self.vface,self.nvface)                                      
                case('v.edge')
                    call read_tetgen_element(unit1,trim(adjustl(ext1(i))),self.vedge,self.nvedge)
                    do j=1,self.nvedge
                        if(self.vedge(j).node(2)>0) then
                        self.vedge(j).property=norm2(self.vnode(self.vedge(j).node(1)).x-self.vnode(self.vedge(j).node(2)).x)
                        else
                            self.vedge(j).property=1.e10
                        endif
                    enddo
                case('v.cell')
                    call read_tetgen_element(unit1,trim(adjustl(ext1(i))),self.vcell,self.nvcell)  
                case('t2e')
                    call read_adj_table(unit1,trim(adjustl(ext1(i))),self.t2e,self.nelt) 
                case('t2f')
                    call read_adj_table(unit1,trim(adjustl(ext1(i))),self.t2f,self.nelt)  
                case('f2e')
                    call read_adj_table(unit1,trim(adjustl(ext1(i))),self.f2e,self.nface)                       
                case default
                    hasread=0
                    print *, 'No such file type=',trim(ext1(i))
                end select
                if(hasread>0) print *, 'Done in reading file=',trim(file1)
            else
                print *, 'file is not exist and skipped. file=',trim(file1)
            endif
        enddo    
        print *,'Done in readin tetgen data.begin to setup adjacent table.'
        call self.setadjlist()
        allocate(self.t2t(4,self.nelt)) 
        do i=1,self.nelt
            self.t2t(:,i)=self.neigh(i).node(1:4)
        enddo
        print *, 'Done in adjacent table setup.'
    
    end subroutine
        
    subroutine read_tetgen_element(unit,ftype,element,nelt)
        integer,intent(in)::unit
        character(len=*)::ftype
        integer,intent(out)::nelt
        type(tetgen_element_tydef),allocatable::element(:)
        
        integer::na1=0,ismarker1=0,n1=0,n2,nelt1,nmax,maxset,i
        integer::nread,nset,nneed,nnode1,n3
        integer::iar1(10)
        
        parameter(nmax=200)
	    parameter(maxset=200)
	
	    real(8)::linedata(nmax),ar1(nmax)
	    character(32)::set(maxset)
        !call skipcomment(unit)
        !read(unit,*) nnode_tg,ndim_tg,na1,ismarker1
        nneed=nmax
        
        call strtoint(unit,linedata,nmax,nread,nneed)
        
        nelt=int(linedata(1))
        
        call lowcase(ftype)
        select case(trim(adjustl(ftype)))
        case('ele')
            n1=sum(linedata(2:nread))
            nnode1=int(linedata(2))
            if(nnode1>4) order=2
            ismarker1=linedata(nread)
        case('face')
            n1=-1 !order=1 by default.
            if(order==2) then
                nnode1=6
            else
                nnode1=3
            endif
            
            ismarker1=linedata(nread)
        case('edge')
            n1=-1 !order=1 by default.
            if(order==2) then
                nnode1=3
            else
                nnode1=2
            endif
            ismarker1=linedata(nread)
        case('neigh')
            n1=4
            nnode1=4
            ismarker1=0
        case('v.cell')
            !if(nread>1) container.type=int(linedata(2))
            !if(nread>2) then
            !    if(container.type==0) then
            !        container.box=reshape(linedata(3:8),([3,2]))
            !    else
            !        container.cylinder.p1=linedata(3:5)
            !        container.cylinder.p2=linedata(6:8)
            !        container.cylinder.r=linedata(9)
            !    endif
            !    container.isread=.true.
            !endif
            
            n1=-1
            nnode1=-1
            ismarker1=0
        case default
            n1=-1
            nnode1=-1
            ismarker1=0
        end select
              
        if(allocated(element))deallocate(element)
        allocate(element(nelt))
        
        select case(trim(adjustl(ftype)))
        case('ele','neigh')
            do i=1,nelt
                read(unit,*) n2,ar1(1:n1)
                element(n2).nnode=nnode1
                element(n2).node=ar1(1:nnode1)
                if(ismarker1>0) element(n2).marker=int(ar1(n1))
            end do
        case('face','edge')
            call strtoint(unit,linedata,nmax,nread,nneed)
       
            n1=nread-1
            n2=int(linedata(1))
            element(n2).node=int(linedata(2:1+nnode1))
            if(ismarker1>0) element(n2).marker=int(linedata(2+nnode1))
            element(n2).nnode=nnode1
            n3=1+nnode1+ismarker1
            if(nread>n3) element(n2).cell(1:(nread-n3))=int(linedata(n3+1:nread))
            
            do i=2,nelt
                read(unit,*) n2,ar1(1:n1)
                element(n2).nnode=nnode1
                element(n2).node=ar1(1:nnode1)
                if(ismarker1>0) element(n2).marker=int(ar1(nnode1+1))
                element(n2).cell(1:(n1-(nnode1+ismarker1)))=int(ar1(nnode1+ismarker1+1:n1))
            end do
!5.2.11  .v.node, .v.edge, .v.face, .v.cell
!
!A .v.node file contains a list of vertices of the Voronoi diagram or the power diagram ( -w switch). Each Voronoi vertex is the circumcenter (or the orthocenter) of a Delaunay (or weighted Delaunay) tetrahedron. The format of .v.node is the same as that of a .node file.
!A .v.edge file contains a list of edges of the Voronoi diagram or the power diagram ( -w switch). Each edge corresponds to a face of the Delaunay (or weighted Delaunay) tetrahedralization. Each Voronoi edge is either a line segment connecting two Voronoi vertices or a ray starting from a Voronoi vertex. The file format of a .v.edge file is
!
!  First line:  <# of edges>
!  Following lines list # of edges:
!    <edge #> <vertex 1> <vertex 2> <V_x> <V_y> <V_z>
!    ...
!<vertex 1> and <vertex 2> are two indices pointing to the list of Voronoi vertices. <vertex 1> must be non-negative, while <vertex 2> may be ?1 which means it is a ray. In this case, the unit vector of this ray is given by <V_x>, <V_y>, and <V_z>.
!A .v.face file contains a list of faces of the Voronoi diagram or the power diagram ( -w switch). Each face corresponds to an edge of the Delaunay (or weighted Delaunay) tetrahedralization. It is formed by a list of Voronoi edges, which may not be closed. The file format of a .v.face file is
!
!  First line:  <# of faces>
!  Following lines list # of faces:
!    <face #> <cell 1> <cell 2> <# of edges> <edge 1> <edge 2> ...
!    ...
!<cell 1> and <cell 2> are two indices pointing into the list of Voronoi cells, i.e., the two cells share this face. <edge 1>, <edge 2> ..., are indices pointing into the edge list, i.e., they are the edges of this face, there are total <# of edges>. If the face is not closed, the index of the last edge of this face is ?1.
!A .v.cell file contains a list of cells of the Voronoi diagram or the power diagram ( -w switch). Each cell corresponds to a vertex of the Delaunay (or weighted Delaunay) tetrahedralization. A cell is formed by a list of Voronoi faces, which may not be closed. The file format of a .v.cell file is
!
!  First line:  <# of cells>
!  Following lines list # of cells:
!    <cell #> <# of faces> <face 1> <face 2> ...
!    ...
!<face 1>, <face 2>, ... are indices pointing into the Voronoi face list. There are total <# of faces> faces. If the cell is not closed, the index of the last face in this cell is ?1.
       case('v.edge')
            do i=1,nelt
                call strtoint(unit,linedata,nmax,nread,nneed)
                n2=int(linedata(1))
                element(n2).node=int(linedata(2:3))
                element(n2).nnode=2
                if(nread>3) then
                    if(nread/=6) error stop 'error in readin numbers for vedge. subroutine=read_tetgen_element'
                    element(n2).v=linedata(4:6)
                endif
                if(element(n2).node(2)>0) element(n2).isclose=1
            enddo
       case('v.face')           
            do i=1,nelt
                call strtoint(unit,linedata,nmax,nread,nneed)
                n2=int(linedata(1))
                element(n2).nedge=int(linedata(4))
                if(nread/=4+element(n2).nedge) error stop 'error in readin numbers for vface. subroutine=read_tetgen_element'
                element(n2).cell=int(linedata(2:3))
                element(n2).edge=int(linedata(5:nread))
                if(element(n2).edge(element(n2).nedge)>0) element(n2).isclose=1                
            enddo  
       case('v.cell')
           !我改过tetgenv.cell的输出格式(在编号后增加此cell对应的颗粒号.),此处只能读入改后版本输出的v.cell文件
            do i=1,nelt
                call strtoint(unit,linedata,nmax,nread,nneed)
                n2=int(linedata(1))
                !the corresponding node for the cell is  node_tg(element(n2).cell(1))                     
                element(n2).cell(1)=int(linedata(2))
                element(n2).nnode=int(linedata(3))
                if(nread/=3+element(n2).nnode) then
                    error stop 'error in readin numbers for vcell. subroutine=read_tetgen_element'
                endif
                element(n2).node=int(linedata(4:nread))
                if(element(n2).node(element(n2).nnode)>0) element(n2).isclose=1    
            enddo             
        end select
        
        close(unit)
    end subroutine
    
    subroutine vtu_vcface_setup(self)
    !!if needed ,call after vfnode_setup
        implicit none
        class(tetgendata_tydef)::self
        integer::i,j,k,n1,n2,if1
        integer::ia1(500) !assume data size <=500
        do i=1,self.nvcell
            if(self.vcell(i).isclose==0) cycle
            n1=1
            ia1(n1)=self.vcell(i).nnode !the data in vcell.node is vface id 
            do j=1,self.vcell(i).nnode
                if1=self.vcell(i).node(j)
                n2=self.vface(if1).nnode
                if(self.isremoveduplicates/=0) then
                    ia1(n1+1:n1+n2+1)=[n2,self.vnode(self.vnode(self.vface(if1).node(1:n2)).uid).marker-1]  !vtu,0-based
                else                    
                    ia1(n1+1:n1+n2+1)=[n2,self.vface(if1).node(1:n2)-1] 
                endif
                n1=n1+n2+1
            enddo
            self.vcell(i).vtu_vcface=ia1(1:n1)
        enddo        
        
    endsubroutine
    
    subroutine vfnode_setup(self)
        implicit none
        class(tetgendata_tydef)::self   
        real(8)::av1(3,3)    
        integer::i,j,e1,v1(2),n1,n2,n3,ic=0,if1
        integer::node1(100),ischeck1(100),edge1(100)
          
           
        !set face node
        do if1=1,self.nvface
            !skip the unclosed faces.
            if(self.vface(if1).isclose==0) cycle
            ischeck1=0
            ic=0
            do while(any(ischeck1(:self.vface(if1).nedge)==0))
                ic=ic+1
                if(ic>self.vface(if1).nedge**2) then
                    error stop 'failed to order the node.sub=update_face'                
                endif
                i=mod(ic-1,self.vface(if1).nedge)+1
                if(ischeck1(i)==1) cycle
            
                if(self.vface(if1).edge(i)<1) then
                    ischeck1(i)=1
                    cycle
                endif
                v1=self.vedge(self.vface(if1).edge(i)).node
                if(i==1) then
                    node1(1:2)=v1
                    n1=2
                    ischeck1(i)=1
                    edge1(1)=self.vface(if1).edge(i)
                else
                    do j=1,2
                        if(node1(n1)==v1(j)) then
                            n2=mod(j,2)+1
                            n1=n1+1
                            node1(n1)=v1(n2)
                            if(n2==2) then
                                edge1(n1-1)=self.vface(if1).edge(i)
                            else
                                edge1(n1-1)=-self.vface(if1).edge(i)
                            endif
                            ischeck1(i)=1
                            exit
                        elseif(node1(1)==v1(j)) then
                            n2=mod(j,2)+1
                            n1=n1+1
                            node1(n1:2:-1)=node1(n1-1:1:-1)
                            node1(1)=v1(n2)
                            edge1(n1-1:2:-1)=edge1(n1-2:1:-1)
                            if(n2==2) then
                                edge1(1)=-self.vface(if1).edge(i)
                            else
                                edge1(1)=self.vface(if1).edge(i)
                            endif
                            ischeck1(i)=1
                            exit
                        endif
                    enddo
                                           
                endif
            enddo
            if(node1(1)==node1(n1)) n1=n1-1 
            self.vface(if1).nnode=n1 
            self.vface(if1).node=node1(1:n1)
        enddo
            
    endsubroutine    

    subroutine vcnode_setup(self)
        implicit none
        class(tetgendata_tydef)::self   
        real(8)::av1(3,3)    
        integer::i,j,e1,v1(2),n1,n2,n3,ic1=0,if1
        !integer::node1(200),n2unode1(200)
        integer,allocatable::unode1(:)
          
           
        !set vcell node
        allocate(unode1(self.nvnode))
        unode1=0
        do ic1=1,self.nvcell
            !skip the unclosed faces.
            if(self.vcell(ic1).isclose==0) cycle
            n1=0
            do i=1,self.vcell(ic1).nnode
                if1=self.vcell(ic1).node(i)
                n2=self.vface(if1).nnode
                unode1(self.vface(if1).node(1:n2))=ic1
                !node1(n1+1:n1+n2)=self.vface(if1).node(1:n2)
                !n1=n1+n2
                !if(n1>200) error stop 'vcnode_setup. vcell node number is > 200.'
            enddo
            !call quick_sort(node1(1:n1),unode1,n2unode1(1:n1))            
            self.vcell(ic1).vcnode=pack([1:self.nvnode],unode1==ic1)
            self.vcell(ic1).nvcnode=size(self.vcell(ic1).vcnode,dim=1)
        enddo
        
        deallocate(unode1)
            
    endsubroutine        
    
    subroutine read_tetgen_node(unit,ftype,element,nelt)
        integer,intent(in)::unit
        character(len=*)::ftype
        integer,intent(out)::nelt
        type(tetgen_node_tydef),allocatable::element(:) 
        
        integer::na1=0,ismarker1=0,n1=0,n2,i,ndim_tg
        real(8),allocatable::ar1(:)
        
        call skipcomment(unit)
        read(unit,*) nelt,ndim_tg,na1,ismarker1
        n1=3+na1+ismarker1
        if(allocated(element)) deallocate(element)
        allocate(element(nelt),ar1(n1))
        !if(na1>0) allocate(node_tg.at(na1))
        do i=1,nelt
            read(unit,*) n2,ar1(1:n1)            
            element(n2).x=ar1(1:3)
            if(na1>0) element(n2).at=ar1(4:3+na1)
            element(n2).na=na1
            if(ismarker1>0) element(n2).marker=ar1(3+na1+1)
        end do
        close(unit)
        deallocate(ar1)
    end subroutine
    
    subroutine read_adj_table(unit,ftype,element,nelt)
        integer,intent(in)::unit,nelt
        character(len=*)::ftype
        integer,allocatable::element(:,:)
        
        integer::i,nnode1,n1
        
        select case(trim(adjustl(ftype)))
        case('f2e')            
            nnode1=3
        case('t2e')
            nnode1=6
        case('t2f')
            nnode1=4        
        end select
        
        if(allocated(element))deallocate(element)
        allocate(element(nnode1,nelt))
        
        read(unit,*) ((n1,element(:,n1)),i=1,nelt)
        
        close(unit)
    end subroutine
    
    
    subroutine find_duplicated_node(node1,box)
    !remove duplicated node and stored in uvertex
        type(tetgen_node_tydef)::node1(:)
        integer::i,j,k,n1,n2
        real(8),intent(in),optional::box(:)
        
        !allocate(ver2node(nver))
        n1=size(node1)
        n2=0
        do i=1,n1
            !check whether the node is inside the box
            if(present(box)) then
                if(node1(i).x(1)<box(1).or.node1(i).x(1)>box(2) &
                    .or.node1(i).x(2)<box(3).or.node1(i).x(2)>box(4) &
                    .or.node1(i).x(3)<box(5).or.node1(i).x(3)>box(6)) then
                    node1(i).marker=-1 !outside the box                    
                endif
            endif
            
            if(node1(i).uid>0) cycle
            do j=i+1,n1
                if(node1(j).uid>0) cycle
                if(node1(j).x(1)<node1(i).x(1)-eps) cycle
                if(node1(j).x(1)>node1(i).x(1)+eps) cycle
                if(node1(j).x(2)<node1(i).x(2)-eps) cycle
                if(node1(j).x(2)>node1(i).x(2)+eps) cycle
                if(node1(j).x(3)<node1(i).x(3)-eps) cycle
                if(node1(j).x(3)>node1(i).x(3)+eps) cycle                
                node1(j).uid=i
            enddo
            node1(i).uid=i
            if(node1(i).marker/=-1)  then
                n2=n2+1
                node1(i).marker=n2 !real nodal id
            endif

        enddo

    end subroutine
       

    SUBROUTINE NODE_ENLARGE_AR(AVAL,DSTEP)
        TYPE(tetgen_node_tydef),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(tetgen_node_tydef),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
    
        LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1,SOURCE=AVAL)
        DEALLOCATE(AVAL)
        ALLOCATE(AVAL(LB1:UB1+DSTEP))
        AVAL(LB1:UB1)=VAL1
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
    END SUBROUTINE

    SUBROUTINE ELEMENT_ENLARGE_AR(AVAL,DSTEP)
        TYPE(tetgen_element_tydef),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(tetgen_element_tydef),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
    
        LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1,SOURCE=AVAL)
        DEALLOCATE(AVAL)
        ALLOCATE(AVAL(LB1:UB1+DSTEP))
        AVAL(LB1:UB1)=VAL1
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
    END SUBROUTINE

    SUBROUTINE I_ENLARGE_AR(AVAL,DSTEP)
        INTEGER,ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        INTEGER,ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0,istat
    
        LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1,SOURCE=AVAL)
        DEALLOCATE(AVAL,STAT=ISTAT)
        ALLOCATE(AVAL(LB1:UB1+DSTEP),STAT=ISTAT)
        AVAL(LB1:UB1)=VAL1
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1,STAT=ISTAT)
    END SUBROUTINE

    subroutine triangle_area_3d ( t, area )

    !*****************************************************************************80
    !
    !! TRIANGLE_AREA_3D computes the area of a triangle in 3D.
    !
    !  Discussion:
    !
    !    This routine uses the fact that the norm of the cross product 
    !    of two vectors is the area of the parallelogram they form.  
    !
    !    Therefore, the area of the triangle is half of that value.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    27 December 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Adrian Bowyer, John Woodwark,
    !    A Programmer's Geometry,
    !    Butterworths, 1983,
    !    ISBN: 0408012420.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) T(3,3), the triangle vertices.
    !
    !    Output, real ( kind = 8 ) AREA, the area of the triangle.
    !
      implicit none

      integer ( kind = 4 ), parameter :: dim_num = 3

      real ( kind = 8 ) area
      real ( kind = 8 ) cross(dim_num)
      real ( kind = 8 ) t(dim_num,3)
    !
    !  Compute the cross product vector.
    !
      cross(1) = ( t(2,2) - t(2,1) ) * ( t(3,3) - t(3,1) ) &
               - ( t(3,2) - t(3,1) ) * ( t(2,3) - t(2,1) )

      cross(2) = ( t(3,2) - t(3,1) ) * ( t(1,3) - t(1,1) ) &
               - ( t(1,2) - t(1,1) ) * ( t(3,3) - t(3,1) )

      cross(3) = ( t(1,2) - t(1,1) ) * ( t(2,3) - t(2,1) ) &
               - ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) )

      area = 0.5D+00 * sqrt ( sum ( cross(1:3)**2 ) )

      return
    end

    
end module