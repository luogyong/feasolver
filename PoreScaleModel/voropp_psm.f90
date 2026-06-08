module voropp
    use meshDS,only:strtoint,ENLARGE_AR
    use penf
    use vtk_fortran, only : vtk_file
    use VTK_CELLTYPE 
    implicit none   
    public::voropp_data_typdef
    
    private
    
    type voropp_face_tydef
        integer::ne,adjc(2)=0,id=0 !this face=globle vccpface(id) 
        !for vppcell.face, when read adjc(1)=adjancent particle id,adjc(2)= the face subid in the faces of the vcell(this.p2c(adjc(1)) 
        !for vppface,adjc(1) right cell id; adjc(2) left cell id
        integer,allocatable::e(:),v(:)    
        !for vppcell.face.v, the id is local.
        !for vppface.v,the id is global
        real(8)::area=-999.,perimeters=-999.,normal(3)=0
    end type    
    type voropp_cell_tydef
        integer::id=0,nv=0,ne=0,nf=0        
        !logical::isP=.true.,isC=.true.
        real(8)::vol,area,center(3),x(3),r=0,perimeters=0.0
        real(8),allocatable::v(:,:)
        type(voropp_face_tydef),allocatable::face(:)
        integer,allocatable::vid(:),vtu_fnode(:) !note that ,the vtu_fnode is the global id
        !character(256)::informat
    contains
        procedure::set_vtu_face_node
    end type
    
    type voropp_data_typdef
        integer::nc=0 !mcell
        integer::nv=0 !nvertex
        integer::nf=0 !nface       
        character(512)::voroformat
        real(8),allocatable::vertex(:,:)
        type(voropp_cell_tydef),allocatable::vppCell(:)                 
        type(voropp_face_tydef),allocatable::vppface(:)   
        integer,allocatable::p2c(:) !particle id to corresponding cell id
        !logical::ismerged=.true. !merge identical node
    contains
        procedure::gen_vppdata
        procedure::read=>voropp_read
        procedure::out=>voropp_out2vtu
        procedure::vppface_set
        procedure,private::vertex_unique_insert        
    endtype
    
    !integer,allocatable::ver2node(:)
    !character(512)::tec_title,FILEPATH,VarString,VarLocation
    !character(len=64)::MeshPassVar,VoroPassVar
    !integer::nvartec=0,nvoropassvar=0
    !real(8),parameter::Pi=3.141592653589793
    real(8)::eps=1.0d-9
    
    contains
    
    subroutine set_vtu_face_node(this)
        implicit none
        class(voropp_cell_tydef)::this
        integer::nsize1,i,n1,n2
        
        nsize1=1+this.nf+sum(this.face.ne)
        allocate(this.vtu_fnode(nsize1))
        n1=1
        this.vtu_fnode(n1)=this.nf
        do i=1,this.nf
            n2=this.face(i).ne
            this.vtu_fnode(n1+1:n1+n2+1)=[n2,this.vid(this.face(i).v)-1]
            n1=n1+n2+1
        enddo
        
    endsubroutine
    
    subroutine gen_vppdata(this,cmd)
        implicit none
        class(voropp_data_typdef)::this
        character(*)::cmd
        INTEGER :: CSTAT, ESTAT
        CHARACTER(100) :: CMSG
            
        print *,'run voro++ to generate data...'
        CALL EXECUTE_COMMAND_LINE (trim(cmd), EXITSTAT=ESTAT,CMDSTAT=CSTAT, CMDMSG=CMSG)
        
        IF (CSTAT > 0) THEN
            PRINT *, "voro++ Command execution failed with error ", TRIM(CMSG)
            pause
        ELSE IF (CSTAT < 0) THEN
            PRINT *, "voro++ Command execution not supported"
            pause
        ELSE
            PRINT *, "Command completed with status ", ESTAT
        END IF  
        
    endsubroutine
    
    subroutine voropp_out2vtu(this,file,vtkformat)
        implicit none
        class(voropp_data_typdef)::this
        character(*),intent(in)::file,vtkformat
        
        integer::i,j,k,nelt1,nnode1
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
        
        
        
        nelt1=0;n1=0;n2=0
        do i=1,this.nc
            nelt1=nelt1+1
            n1=n1+this.vppcell(i).nv
            if(.not.allocated(this.vppcell(i).vtu_fnode)) call this.vppcell(i).set_vtu_face_node()
            n2=n2+size(this.vppcell(i).vtu_fnode)
        enddo

        
        
        allocate(cell_type(nelt1),offset(nelt1),connect(n1),face(n2),faceoffset(nelt1))
        !cell_type=3;offset=2;
        n1=0;nelt1=0;n3=0;
        do i=1,this.nc           
            nelt1=nelt1+1
            n2=this.vppcell(i).nv
            
            connect(n1+1:n1+n2)=this.vppcell(i).vid-1 !vtu is 0-based
            
            n1=n1+n2
            cell_type(nelt1)=VTK_POLYHEDRON
            !face
            n4=size(this.vppcell(i).vtu_fnode)
            face(n3+1:n3+n4)=this.vppcell(i).vtu_fnode(1:n4)
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
        
        ix1=[1:this.nv]
   
        n2=size(ix1,dim=1)
        
        error = a_vtk_file%initialize(format=trim(vtkformat), filename=trim(file)//'_voropp.vtu', mesh_topology='UnstructuredGrid')
        error = a_vtk_file%xml_writer%write_piece(np=n2, nc=nelt1)
        error = a_vtk_file%xml_writer%write_geo(np=n2, nc=nelt1, x=this.vertex(1,ix1), y=this.vertex(2,ix1), z=this.vertex(3,ix1))
        error = a_vtk_file%xml_writer%write_connectivity(nc=nelt1, connectivity=connect, offset=offset, cell_type=cell_type,face=face,faceoffset=faceoffset)
    
        !error = a_vtk_file%xml_writer%write_dataarray(location='node', action='open') 
        !do i=1,this.vertex(1).na
        !    if(i==1) then
        !        !注意,在此假定(默认)第一个属性值为半径的平方.
        !        allocate(x1(n2))
        !        do j=1,n2
        !            x1(j)=sqrt(this.vertex(ix1(j)).at(i))
        !        enddo            
        !        error = a_vtk_file%xml_writer%write_dataarray(data_name='Rt(正交球半径)', x=x1)
        !        deallocate(x1)
        !    else
        !        write(ch1,'(A)') i
        !        allocate(x1(n2))
        !        do j=1,n2
        !            x1(j)=this.vertex(ix1(j)).at(i)
        !        enddo 
        !        error = a_vtk_file%xml_writer%write_dataarray(data_name='at'//trim(adjustl(ch1)), x=x1)
        !        deallocate(x1)
        !    endif
        !enddo
        !error = a_vtk_file%xml_writer%write_dataarray(data_name='Marker', x=this.vertex(ix1).marker)
        !
        !error = a_vtk_file%xml_writer%write_dataarray(location='node', action='close')
             
        error = a_vtk_file%xml_writer%write_dataarray(location='cell', action='open')
        
        !if(allocated(ix1)) deallocate(ix1)
        !allocate(ix1(nelt1))    
        !n1=0
        !do i=1,this.nc
        !
        !    n1=n1+1
        !    ix1(n1)=this.vppcell(i).marker
        !enddo
        error = a_vtk_file%xml_writer%write_dataarray(data_name='Vol', x=this.vppcell(:this.nc).vol) 
  
        !n1=0
        !do i=1,this.nc
        !
        !    n1=n1+1
        !    ix1(n1)=i
        !enddo
        error = a_vtk_file%xml_writer%write_dataarray(data_name='iParticle', x=this.vppcell(:this.nc).id) 
        !deallocate(ix1)        
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
        CALL EXECUTE_COMMAND_LINE ('powershell -command "Get-Content '//trim(file)//'_voropp.vtu'//' | Set-Content -Encoding utf8 '//trim(file)//'_voropp-utf8.vtu'//'"', EXITSTAT=ESTAT,CMDSTAT=CSTAT, CMDMSG=CMSG)
        
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
        CALL EXECUTE_COMMAND_LINE ('del "'//trim(file)//'_voropp.vtu"', EXITSTAT=ESTAT,CMDSTAT=CSTAT, CMDMSG=CMSG)
        
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
     
    subroutine vppface_set(this)
        implicit none
        class(voropp_data_typdef)::this
        integer::i,j,k,n1,n2
        real(8)::vec1(3),t1
        
        !A list of the neighboring particle or wall IDs corresponding to each face. The list can contain negative numbers. For the non-periodic case these 
        !correspond to when the particles have faces created by the edges of the computational region. The numbers -1 to -6 correspond to the minimum x, 
        !maximum x, minimum y, maximum y, minimum z, and maximum z walls respectively. For periodic boundary conditions, negative numbers correspond to the 
        !cases when a face of the Voronoi cell is created by the periodic image of the current particle.
        !
        !In general, the neighbor information will be symmetric, so that if particle A reports particle B as a neighbor, then particle B will report particle 
        !A as a neighbor. However, due to the fact that Voro++ computes each Voronoi cell individually, it does not provide an explicit guarantee that the 
        !neighbor information will always be symmetric. Suppose there is a very small Voronoi face connecting A to B – it may be the case that due to roundoff 
        !error, the Voronoi cell computed for particle A has a face connecting it to B, but the cell computed for particle B does not have a face connecting 
        !it to A. If the user requires perfectly symmetric neighbor information, this can be achieved by scanning the output for any one-sided connections, 
        !and either deleting them or adding in the reverse connections. The face areas output from “%f” can also be used to remove connections between 
        !particles that only have a very small face between them.          
        
        if(.not.allocated(this.vppface)) allocate(this.vppface(sum(this.vppCell.nf)))
        this.nf=0
        do i=1,this.nc
            do j=1,this.vppCell(i).nf
                if(this.vppCell(i).face(j).id==0) then
                    this.nf=this.nf+1
                    this.vppCell(i).face(j).id=this.nf
                    n1=this.vppCell(i).face(j).adjc(1)                    
                    !check pairs
                    if(n1>0) then
                        n1=this.p2c(n1)
                        !this.vppCell(i).face(j).adjc(1)=n1 !change it to celll id
                        
                        do k=1,this.vppCell(n1).nf
                            if(this.vppCell(n1).face(k).adjc(1)==this.vppCell(i).id) then
                                this.vppCell(n1).face(k).id=-this.nf 
                                !- ,the face normal is negetive to the face in vppface.
                                this.vppCell(n1).face(k).adjc(2)=j
                                this.vppCell(i).face(j).adjc(2)=k
                                !this.vppCell(n1).face(k).adjc(1)=i !change it to cell id
                                exit
                            endif    
                        enddo
                        if(k>this.vppCell(n1).nf) then
                            write(*,10) i,n1,this.vppCell(i).face(j).area
                            error stop "error stop at sub=vppface_set"
                        endif
                    endif
                    this.vppface(this.nf)=this.vppCell(i).face(j)
                    this.vppface(this.nf).v=this.vppCell(i).vid(this.vppface(this.nf).v)
                    
                    !for vppface,adjc(1) right cell id; adjc(2) left cell
                    !vec1=this.vppCell(i).center-this.vppCell(i).v(:,vppface(this.nf).v(1))
                    !t1=dot_product(vec1,this.vppCell(i).face(j).normal)                    
                    !if(t1<0.d0) then
                        this.vppface(this.nf).adjc(1)=max(n1,0) 
                        this.vppface(this.nf).adjc(2)=i
                    !else
                    !    this.vppface(this.nf).adjc(1)=i
                    !    this.vppface(this.nf).adjc(2)=n1                       
                    !                       
                    !endif
                endif
            enddo
        enddo
10  format('No face pair found between cell ',i7, ' and cell ',i7,'the face area is ',g)
    endsubroutine
    
    subroutine vertex_unique_insert(this,point,aindex)
    !insert point(:,:)  to this.vertex if they are not in this.vertex and 
    !return their indices by aindex
        use omp_lib
        implicit none
        class(voropp_data_typdef)::this
        real(8),intent(in)::point(3,*)
        integer,intent(out)::aindex(:)
        
        integer::i,j,n1,n2
        real(8)::t1
        
        n1=size(aindex,dim=1)
        if(.not.allocated(this.vertex)) allocate(this.vertex(3,10000))
        do i=1,n1
            n2=0
            do j=1,this.nv
                if(point(1,i)<this.vertex(1,j)-eps) cycle
                if(point(1,i)>this.vertex(1,j)+eps) cycle
                if(point(2,i)<this.vertex(2,j)-eps) cycle
                if(point(2,i)>this.vertex(2,j)+eps) cycle
                if(point(3,i)<this.vertex(3,j)-eps) cycle
                if(point(3,i)>this.vertex(3,j)+eps) cycle                
                !t1=norm2(point(:,i)-this.vertex(:,j))
                !if(abs(t1)<eps) then
                    aindex(i)=j
                    n2=1
                    exit
                !endif
            enddo
            if(n2==0) then
                this.nv=this.nv+1
                if(this.nv>size(this.vertex,dim=2)) call enlarge_ar(this.vertex,100)
                this.vertex(:,this.nv)=point(:,i)
                aindex(i)=this.nv
            endif
        enddo
    end subroutine
    
    subroutine voropp_read(this,iunit,file)
        !class(voropp_cell_tydef)::this
        !use dflib
        USE IFPORT
        implicit none
        class(voropp_data_typdef)::this
        integer,intent(in),optional::iunit
        CHARACTER(len=*),intent(in),optional::file
        CHARACTER(3)        drive
	    CHARACTER(512)      dir
	    CHARACTER(512)      name,file1
	    CHARACTER(16)      ext,ext1(12)
	    integer(4)::length,msg
        logical::isexist,isin1
        integer::unit,i,hasread

               
        integer::na1=0,ismarker1=0,n1=0,n2,nelt1,j,k
        integer::nread,nset,nneed,nnode1,n3
        integer::iar1(10),ef
        integer,parameter::nmax=1000
	    integer,parameter::maxset=100
        
	    real(8)::linedata(nmax),ar1(nmax)
	    character(256)::set(maxset)
        
        if(present(iunit)) then
            inquire(iUNIT,name=file1,EXIST=isexist)
            unit=iunit
            if(.not.isexist) file1=file
        else
            unit=10
            file1=file
        endif
        
		!length = SPLITPATHQQ(file1, drive, dir, name, ext)
		!tec_title=trim(name)
  !      FILEPATH=trim(drive)//trim(dir)
  !      msg = CHDIR(FILEPATH)
  !      FILEPATH=trim(drive)//trim(dir)//trim(name)
        
        !if(trim(adjustl(ext))=='.cell') then
        !    FILEPATH=filepath(:len_trim(filepath)-2)
        !endif
        close(unit)        
        
        !call skipcomment(unit)
        !read(unit,*) nnode_tg,ndim_tg,na1,ismarker1
        
        open(unit,file=trim(file1),status='old')
        nneed=nmax
        ef=0
        call strtoint(unit,linedata,nmax,nread,nneed,set,maxset,nset)
        this.nc=int(linedata(1))
        allocate(this.vppCell(this.nc),this.p2c(this.nc))
        this.p2c=0;
        
        this.voroformat=''
        do i=1,nset
            n1=index(set(i),'%')
            if(n1>0) then
                this.voroformat=trim(this.voroformat)//set(i)(n1+1:n1+1)
            endif            
        enddo
        !input order check
        n1=index(trim(this.voroformat),'i')
        if(n1/=1) then
            error stop "outputformat: %i should be in the first place."
        endif
        n1=maxval([index(trim(this.voroformat),'q'),index(trim(this.voroformat),'x') ,index(trim(this.voroformat),'y') ,index(trim(this.voroformat),'z')])        
        n2=index(trim(this.voroformat),'w')        
        n3=index(trim(this.voroformat),'p')
        if(n3>0.and.(n1>n3.or.n1==0)) then
            write(*,40) 'p'
            error stop
        endif        
        if(n3>0.and.(n2>n3.or.n2==0)) then
            write(*,30) 'w','p'
            error stop
        endif
        n3=index(trim(this.voroformat),'P')
        if(n3>0.and.(n2>n3.or.n2==0)) then
            write(*,30) 'w','P'
            error stop
        endif
        n3=index(trim(this.voroformat),'c')
        if(n3>0.and.(n1>n3.or.n1==0)) then
            write(*,40) 'c'
            error stop
        endif          
        n2=index(trim(this.voroformat),'s')        
        n3=minval([index(trim(this.voroformat),'e'),index(trim(this.voroformat),'f'),index(trim(this.voroformat),'a'),index(trim(this.voroformat),'t'), &
           & index(trim(this.voroformat),'l'),index(trim(this.voroformat),'n')])
        if(n3>0.and.(n2>n3.or.n2==0)) then
            write(*,50) 's'
            error stop
        endif  
        n2=index(trim(this.voroformat),'a')        
        n3=index(trim(this.voroformat),'t')
        if(n3>0.and.(n2>n3.or.n2==0)) then
            write(*,30) 'a','t'
            error stop
        endif        
        
        n1=0;
        do i=1,this.nc     
            call strtoint(unit,linedata,nmax,nread,nneed,set,maxset,nset,ef)
            if(ef<0) exit
            n1=n1+1
            n2=1
            n3=int(linedata(n2))
            if(n3>size(this.p2c)) call enlarge_ar(this.p2c,100)
            this.p2c(n3)=n1
            isin1=.false.
            do j=1,len_trim(this.voroformat)
               
                select case(this.voroformat(j:j))
                !Particle-related:
                !  %i The particle ID number
                !  %x The x coordinate of the particle
                !  %y The y coordinate of the particle
                !  %z The z coordinate of the particle
                !  %q The position vector of the particle, short for "%x %y %z"
                !  %r The radius of the particle (only printed if -r enabled)                    
                case('i')
                    this.vppCell(n1).id=int(linedata(n2))
                    n2=n2+1
                case('x')
                    this.vppCell(n1).x(1)=linedata(n2)
                    n2=n2+1
                    isin1=.true.
                case('y')
                    this.vppCell(n1).x(2)=linedata(n2)
                    n2=n2+1 
                    isin1=.true.
                case('z')
                    this.vppCell(n1).x(3)=linedata(n2)
                    n2=n2+1 
                    isin1=.true.
                case('q')
                    this.vppCell(n1).x=linedata(n2:n2+2)
                    n2=n2+3
                    isin1=.true.
                case('r')
                    this.vppCell(n1).r=linedata(n2)
                    n2=n2+1 
                !this.vertex-related:
                !%w The number of vertices in the Voronoi cell
                !%p A list of the vertices of the Voronoi cell in the format (x,y,z),
                !   relative to the particle center
                !%P A list of the vertices of the Voronoi cell in the format (x,y,z),
                !   relative to the global coordinate system
                !%o A list of the orders of each this.vertex
                !%m The maximum radius squared of a this.vertex position, relative to the
                !   particle center
                case('w')
                    this.vppCell(n1).nv=int(linedata(n2))
                    allocate(this.vppCell(n1).v(3,this.vppCell(n1).nv))
                    n2=n2+1
                case('p')
                    !this.vppCell(n1).isP=.false.
                 
                    this.vppCell(n1).v=reshape(linedata(n2:n2+3*this.vppCell(n1).nv-1),([3,this.vppCell(n1).nv]))
                    n2=n2+3*this.vppCell(n1).nv
                    do k=1,this.vppCell(n1).nv
                        this.vppCell(n1).v(:,k)=this.vppCell(n1).v(:,k)+this.vppCell(n1).x
                    enddo                    
                    allocate(this.vppCell(n1).vid(this.vppCell(n1).nv))
                    call this.vertex_unique_insert(this.vppCell(n1).v,this.vppCell(n1).vid)
                case('P')
                      
                    this.vppCell(n1).v=reshape(linedata(n2:n2+3*this.vppCell(n1).nv-1),([3,this.vppCell(n1).nv]))
                    n2=n2+3*this.vppCell(n1).nv
                    allocate(this.vppCell(n1).vid(this.vppCell(n1).nv))
                    call this.vertex_unique_insert(this.vppCell(n1).v,this.vppCell(n1).vid)                    
                case('o')
                    !not used yet.skipped 
                    n2=n2+this.vppCell(n1).nv
                case('m')
                    !not used yet. skipped
                    n2=n2+1
                !Edge-related:
                !  %g The number of edges of the Voronoi cell
                !  %E The total edge distance
                !  %e A list of perimeters of each face
                case('g')
                    this.vppCell(n1).ne=int(linedata(n2))
                    n2=n2+1
                case('E')
                    this.vppCell(n1).perimeters=linedata(n2)
                    n2=n2+1
                case('e')
                    this.vppCell(n1).face.perimeters=linedata(n2:n2+this.vppCell(n1).nf-1)
                    n2=n2+this.vppCell(n1).nf
                !Face-related:
                !  %s The number of faces of the Voronoi cell
                !  %F The total surface area of the Voronoi cell
                !  %A A frequency table of the number of edges for each face
                !  %a A list of the number of edges for each face
                !  %f A list of areas of each face
                !  %t A list of bracketed sequences of vertices that make up each face
                !  %l A list of normal vectors for each face
                !  %n A list of neighboring particle or wall IDs corresponding to each face   
                !A list of the neighboring particle or wall IDs corresponding to each face. The list can contain negative numbers. For the non-periodic case these 
                !correspond to when the particles have faces created by the edges of the computational region. The numbers -1 to -6 correspond to the minimum x, 
                !maximum x, minimum y, maximum y, minimum z, and maximum z walls respectively. For periodic boundary conditions, negative numbers correspond to the 
                !cases when a face of the Voronoi cell is created by the periodic image of the current particle.
                !
                !In general, the neighbor information will be symmetric, so that if particle A reports particle B as a neighbor, then particle B will report particle 
                !A as a neighbor. However, due to the fact that Voro++ computes each Voronoi cell individually, it does not provide an explicit guarantee that the 
                !neighbor information will always be symmetric. Suppose there is a very small Voronoi face connecting A to B – it may be the case that due to roundoff 
                !error, the Voronoi cell computed for particle A has a face connecting it to B, but the cell computed for particle B does not have a face connecting 
                !it to A. If the user requires perfectly symmetric neighbor information, this can be achieved by scanning the output for any one-sided connections, 
                !and either deleting them or adding in the reverse connections. The face areas output from “%f” can also be used to remove connections between 
                !particles that only have a very small face between them.                    
                case('s')
                    this.vppCell(n1).nf=int(linedata(n2))
                    allocate(this.vppCell(n1).face(this.vppCell(n1).nf))
                    n2=n2+1
                case('F')
                    this.vppCell(n1).area=linedata(n2)
                    n2=n2+1
                case('f')
                    this.vppCell(n1).face.area=linedata(n2:n2+this.vppCell(n1).nf-1)
                    n2=n2+this.vppCell(n1).nf
                case('a')
                    this.vppCell(n1).face.ne=int(linedata(n2:n2+this.vppCell(n1).nf-1))
                    n2=n2+this.vppCell(n1).nf
                    
                case('t')
                    do k=1,this.vppCell(n1).nf
                     
                        this.vppCell(n1).face(k).v=int(linedata(n2:n2+this.vppCell(n1).face(k).ne-1))+1 !base 0 
                        n2=n2+this.vppCell(n1).face(k).ne
                    enddo
                case('l')
                    do k=1,this.vppCell(n1).nf
                        !allocate(this.vppCell(n1).face(k).v(this.vppCell(n1).face(k).ne))
                        this.vppCell(n1).face(k).normal=linedata(n2:n2+2)
                        n2=n2+3
                    enddo                    
                case('n')
                    do k=1,this.vppCell(n1).nf
                        !allocate(this.vppCell(n1).face(k).v(this.vppCell(n1).face(k).ne))
                        this.vppCell(n1).face(k).adjc(1)=int(linedata(n2))
                        n2=n2+1
                    enddo
                !
                !Volume-related:
                !  %v The volume of the Voronoi cell
                !  %c The centroid of the Voronoi cell, relative to the particle center
                !  %C The centroid of the Voronoi cell, in the global coordinate system   
                case('v')
                    this.vppCell(n1).vol=linedata(n2)
                    n2=n2+1
                case('c')
                    this.vppCell(n1).center=linedata(n2:n2+2)+this.vppCell(n1).x
                    n2=n2+3
                case('C')
                    !this.vppCell(n1).isC=.true.
                    this.vppCell(n1).center=linedata(n2:n2+2)
                    n2=n2+3                    
                case default                
                    write(*, 10), this.voroformat(j:j)
                    stop
                end select
            enddo
            if(nread/=n2-1) then
                write(*,20) n2-1,nread
                stop
            endif
                
        enddo
        this.nc=n1
        
        close(unit)    
10      format('Custom output format: %',A1,' is not supported yet.')
20      format("The number of readin data is expected to be ",i4,".but it is ",i4,'.sub=voropp_read')
30      format('Custom output format: %',A1,' should be before: %',A1,'.')
40      format('Custom output format: %q/(%x,%y,%z) should be before: %',A1,'.')
50      format('Custom output format: %e/%f/%a/%t/%l/%n should be after: %',A1,'.')        
    end subroutine
 
  
    
end module 