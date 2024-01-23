module triangle_io
    use ds_t,only:arr_t,inpn
    use meshds,only:SEG,NSEG,CSL,CLN,Find_Point_Inside_Polygon_2D,ZONE,ZNUM, &
        & write_help,readline,property,pro_num,Err_msg,path_name,title, &
        & point_tydef,element_tydef,Edge_tydef,adjlist_tydef,addadjlist,v2edge,&
        & constrained_edge_tydef,soillayer,strtoint,SETUP_SEARCH_ZONE_2D,POINTlOC_2D,&
        & xmin,xmax,ymin,ymax,xyscale,SEARCHZONE_TYDEF,isnorefined
    use sample_pdf,only:example_pdf
    use tetgen
    use, intrinsic :: iso_c_binding
    implicit none
    private
    public::libtriangle,bgmesh,triangle_tydef

    ! type mesh_lgy_tydef
    !     integer::nnode=0,nelt=0,nedge=0,ncedge=0,NSZONE=0
    !     type(point_tydef),allocatable::node(:)
    !     type(element_tydef),allocatable::element(:)
    !     type(edge_tydef),allocatable::edge(:)
    !     type(adjlist_tydef),allocatable::adjlist(:)
    !     type(constrained_edge_tydef),allocatable::cedge(:)
    ! contains
    !     procedure::tri2mesh_lgy
    ! end type


! struct triangulateio {
!   REAL *pointlist;                                               /* In / out */
!   REAL *pointattributelist;                                      /* In / out */
!   int *pointmarkerlist;                                          /* In / out */
!   int numberofpoints;                                            /* In / out */
!   int numberofpointattributes;                                   /* In / out */

!   int *trianglelist;                                             /* In / out */
!   REAL *triangleattributelist;                                   /* In / out */
!   REAL *trianglearealist;                                         /* In only */
!   int *neighborlist;                                             /* Out only */
!   int numberoftriangles;                                         /* In / out */
!   int numberofcorners;                                           /* In / out */
!   int numberoftriangleattributes;                                /* In / out */

!   int *segmentlist;                                              /* In / out */
!   int *segmentmarkerlist;                                        /* In / out */
!   int numberofsegments;                                          /* In / out */

!   REAL *holelist;                        /* In / pointer to array copied out */
!   int numberofholes;                                      /* In / copied out */

!   REAL *regionlist;                      /* In / pointer to array copied out */
!   int numberofregions;                                    /* In / copied out */

!   int *edgelist;                                                 /* Out only */
!   int *edgemarkerlist;            /* Not used with Voronoi diagram; out only */
!   REAL *normlist;                /* Used only with Voronoi diagram; out only */
!   int numberofedges;                                             /* Out only */
! };
    type,bind(c)::triangulateio_tydef

        type(c_ptr)::pointlist=c_null_ptr;                                             ! /* In / out */
        type(c_ptr)::pointattributelist;                                    !  /* In / out */
        type(c_ptr)::pointmarkerlist;                                        !  /* In / out */
        integer(c_int):: numberofpoints=0;                                          !  /* In / out */
        integer(c_int):: numberofpointattributes=0;                                 !  /* In / out */

        type(c_ptr)::trianglelist;                                           !  /* In / out */
        type(c_ptr)::triangleattributelist;                                 !  /* In / out */
        type(c_ptr)::trianglearealist;                                      !   /* In only */
        type(c_ptr)::neighborlist;                                           !  /* Out only */
        integer(c_int):: numberoftriangles;                                       !  /* In / out */
        integer(c_int):: numberofcorners;                                         !  /* In / out */
        integer(c_int):: numberoftriangleattributes;                              !  /* In / out */

        type(c_ptr)::segmentlist;                                            ! /* In / out */
        type(c_ptr)::segmentmarkerlist;                                      !  /* In / out */
        integer(c_int):: numberofsegments;                                        !  /* In / out */

        type(c_ptr)::holelist;                        !/* In / pointer to array copied out */
        integer(c_int):: numberofholes=0;                     !                 /* In / copied out */

        type(c_ptr)::regionlist;                      !/* In / pointer to array copied out */
        integer(c_int):: numberofregions=0;                   !                 /* In / copied out */

        type(c_ptr)::edgelist;                         !                        /* Out only */
        type(c_ptr)::edgemarkerlist;           !/* Not used with Voronoi diagram; out only */
        type(c_ptr)::normlist;                !/* Used only with Voronoi diagram; out only */
        integer(c_int):: numberofedges;             !                                /* Out only */

    end type
    
    ! type triangulateio_fortran_tydef
    !     !type(triangulateio_tydef)::tridata

    ! contains
    !     procedure,private,nopass::triangulateio_real_array_setvalue,triangulateio_int_array_setvalue,&
    !                        triangulateio_real_array_getvalue,triangulateio_int_array_getvalue
    !     generic::setvalue=>triangulateio_real_array_setvalue,triangulateio_int_array_setvalue
    !     generic::getvalue=>triangulateio_real_array_getvalue,triangulateio_int_array_getvalue
    !     procedure,nopass::free=>triangulateio_free_array        
    !     procedure::triangulate=>triangulateio_triangulate
    ! endtype


    logical::isbgready=.false.
    
    type triangle_tydef
        integer::method=-1,nnode=0,nelt=0,nedge=0,ncedge=0,NSZONE=0,islib=0,dim=2        
        character(1024)::path='',cmd='-penI',polyfile=''
        type(point_tydef),allocatable::node(:)
        type(element_tydef),allocatable::element(:)
        type(edge_tydef),allocatable::edge(:)
        type(adjlist_tydef),allocatable::adjlist(:)
        type(constrained_edge_tydef),allocatable::cedge(:)
        TYPE(SEARCHZONE_TYDEF),ALLOCATABLE::SEARCHZONE(:)
        type(triangulateio_tydef)::in,out,vorout
        type(tetgenio)::tetin,tetout
        character(8192):: helpstring= &
            &"相对于程序自带算法,triangle能够实现用更少的节点实现区域的划分，比如生成geo三模型时，可选用此选项. \n &
            & triangle的输入格式为: \n &
            & triangle,method=I;path=string,cmd=string,poly=string,islib=N0Y1,dim=2|3 \n &
            & NOTE: \n &
            & 0)method=0,根据模型在当前目录下生成poly文件,然后根据cmd命令进行网格划分,并读入网络数据文件(node,ele,edge和neigh文件); \n &
            &   method=1,利用当前目录下已有poly文件,根据cmd命令进行网格划分,并读入网络数据; \n &
            &   method=2,根据模型在当前目录下生成poly文件,然后退出; \n &
            &   method=3,仅读入triangle生成的网格数据文件; \n &
            &   method=4,根据模型在当前目录下生成poly文件，然后根据cmd命令进行网格划分后退出; \n &
            &   method=5,执行method=0的命令后,更新各节点node.atrrib(1)为颗粒大小随机场,然后按命令-rpqu重新划分后，再次读入网格数据文件; \n &
            &   method=-1(默认),令本命令无效(即不执行上述任一任务). \n &
            & 1)path=triangle.exe所在的目录.比如'c:\triangle.exe'  不输入时默认为当前目录或已设环境变量.\n & 
            & 2)cmd=triangle的命令参数.比如'cmd='-penI'' 不输入时默认为'-penI'. \n &
            & 3)polyfile=triangle的输入文件.比如'polyfile='.\a.poly'. 不输入时默认为当前目录下与mes文件同名的poly的文件。\n &
            & 4)当采用第三方程序Triangle进行网格划分且要进行地层插值时，所有节点都必须输入确定的地层信息。\N &
            & 5)islib=1时,采用triangle.lib提供的io的方式进行处理,而不是通过文件进行数据交换,默认为0。\N &
            & 6)dim=3(默认为2)时,读入poly采用tetgen进行三维网格划分,此时poly必须输入。\N &            
            & 7)triangle命令参数可参考:https://www.cs.cmu.edu/~quake/triangle.html. 常用的命令参数为\n &
            & -p Triangulates a Planar Straight Line Graph (.poly file). \n &
            & -r Refines a previously generated mesh. \n &
            & -q Quality mesh generation with no angles smaller than 20 degrees. An alternate minimum angle may be specified after the `q'. \n &
            & -a Imposes a maximum triangle area constraint. A fixed area constraint (that applies to every triangle) may be specified after the `a', or varying area constraints may be read from a .poly file or .area file. \n &
            & -u Imposes a user-defined constraint on triangle size.\n &
            & -A Assigns a regional attribute to each triangle that identifies what segment-bounded region it belongs to. \n &
            & -c Encloses the convex hull with segments.\n &
            & -D Conforming Delaunay: use this switch if you want all triangles in the mesh to be Delaunay, and not just constrained Delaunay; or if you want to ensure that all Voronoi vertices lie within the triangulation.\n &
            & -j Jettisons vertices that are not part of the final triangulation from the output .node file (including duplicate input vertices and vertices ``eaten'' by holes).\n &
            & -e Outputs (to an .edge file) a list of edges of the triangulation.\n &
            & -v Outputs the Voronoi diagram associated with the triangulation. Does not attempt to detect degeneracies, so some Voronoi vertices may be duplicated.\n &
            & -n Outputs (to a .neigh file) a list of triangles neighboring each triangle.\n &
            & -g Outputs the mesh to an Object File Format (.off) file, suitable for viewing with the Geometry Center's Geomview package.\n &
            & -I Suppresses mesh iteration numbers.\n &
            & -O Suppresses holes: ignores the holes in the .poly file.\n &
            & -Y Prohibits the insertion of Steiner points on the mesh boundary. If specified twice (-YY), it prohibits the insertion of Steiner points on any segment, including internal segments.\n &
            & -S Specifies the maximum number of added Steiner points. \n &
            & -i Uses the incremental algorithm for Delaunay triangulation, rather than the divide-and-conquer algorithm.\n &
            & -F Uses Steven Fortune's sweepline algorithm for Delaunay triangulation, rather than the divide-and-conquer algorithm.\n &
            & -Q Quiet: Suppresses all explanation of what Triangle is doing, unless an error occurs.\n &
            & -V Verbose: Gives detailed information about what Triangle is doing. Add more `V's for increasing amount of detail. `-V' gives information on algorithmic progress and detailed statistics.\n &
            & -h Help: Displays complete instructions. \n &
            & 8) tetgen 的命令参数可参考：https://wias-berlin.de/software/tetgen/1.5/doc/manual/manual005.html 。常用的命令参数为：\n &
            & -p Tetrahedralizes a piecewise linear complex (PLC).\n &
            & -Y Preserves the input surface mesh (does not modify it).\n &
            & -r Reconstructs a previously generated mesh.\n &
            & -q Refines mesh (to improve mesh quality).\n &
            & -R Mesh coarsening (to reduce the mesh elements).\n &
            & -A Assigns attributes to tetrahedra in different regions.\n &
            & -a Applies a maximum tetrahedron volume constraint.\n &
            & -m Applies a mesh sizing function.\n &
            & -i Inserts a list of additional points.\n &
            & -O Specifies the level of mesh optimization.\n &
            & -S Specifies maximum number of added points.\n &
            & -T Sets a tolerance for coplanar test (default 10-8).\n &
            & -X Suppresses use of exact arithmetic.\n &
            & -M No merge of coplanar facets or very close vertices.\n &
            & -w Generates weighted Delaunay (regular) triangulation.\n &
            & -c Retains the convex hull of the PLC.\n &
            & -d Detects self-intersections of facets of the PLC.\n &
            & -z Numbers all output items starting from zero.\n &
            & -f Outputs all faces to .face file.\n &
            & -e Outputs all edges to .edge file.\n &
            & -n Outputs tetrahedra neighbors to .neigh file.\n &
            & -v Outputs Voronoi diagram to files.\n &
            & -g Outputs mesh to .mesh file for viewing by Medit.\n &
            & -k Outputs mesh to .vtk file for viewing by Paraview.\n &
            & -J No jettison of unused vertices from output .node file.\n &
            & -B Suppresses output of boundary information.\n &
            & -N Suppresses output of .node file.\n &
            & -E Suppresses output of .ele file.\n &
            & -F Suppresses output of .face and .edge file.\n &
            & -I Suppresses mesh iteration numbers.\n &
            & -C Checks the consistency of the final mesh.\n &
            & -Q Quiet: No terminal output except errors.\n &
            & -V Verbose: Detailed information, more terminal output.\n &
            & -h Help: A brief instruction for using TetGen.\n &
            & "C        
    contains
        procedure,nopass::help=>write_help
        procedure::readin=>triangle_option_read
        procedure::outpoly=>write_poly_file
        procedure::readmesh=>read_triangle_file
        procedure::meshing=>meshbytriangle
        procedure::exe=>execute_triangle_cmd
        procedure::getdata=>get_triangle_meshdata
        procedure,nopass::set_searchzone=>SETUP_SEARCH_ZONE_2D
        procedure::getattrib=>interpolate_attrib_bgmesh
        procedure::out_node=>triangle_write_node_file
        !use triangle.lib to mesh
        procedure,private,nopass::triangulateio_real_array_setvalue,triangulateio_int_array_setvalue,&
                           triangulateio_real_array_getvalue,triangulateio_int_array_getvalue
        generic::setvalue=>triangulateio_real_array_setvalue,triangulateio_int_array_setvalue
        generic::getvalue=>triangulateio_real_array_getvalue,triangulateio_int_array_getvalue
        procedure::triangulate=>triangulateio_triangulate
        procedure,nopass::free=>triangulateio_free_array  
        procedure::set_indata=>triangulateio_set_indata
        procedure::tri2mesh_lgy=>triangulateio_tri2mesh_lgy
    endtype
    type(triangle_tydef)::libtriangle,bgmesh
    
    CHARACTER(1024)::POLY2D_FILE_FORMAT= &
        &"#2D POLY FILE FORMAT \n &        
        &#First line: <# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)> \n &
        &#Following lines: <vertex #> <x> <y> [attributes] [boundary marker] \n &
        &#One line: <# of segments> <# of boundary markers (0 or 1)> \n &
        &#Following lines: <segment #> <endpoint> <endpoint> [boundary marker] \n &
        &#One line: <# of holes> \n &
        &#Following lines: <hole #> <x> <y> \n &
        &#Optional line: <# of regional attributes and/or area constraints> \n &
        &#Optional following lines: <region #> <x> <y> <attribute> <maximum area> \n &
        &"C        
    
    interface
        function real_allocate(n) bind(c)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_size_t),value::n
            type(c_ptr)::real_allocate
        end function

        function int_allocate(n) bind(c)
            use, intrinsic :: iso_c_binding
            implicit none
            integer(c_size_t),value::n
            type(c_ptr)::int_allocate
        end function

        SUBROUTINE real_deallocate(p) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR
            IMPLICIT NONE
            TYPE(C_PTR), INTENT(IN), VALUE :: p
        END SUBROUTINE real_deallocate

        SUBROUTINE int_deallocate(p) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR
            IMPLICIT NONE
            TYPE(C_PTR), INTENT(IN), VALUE :: p
        END SUBROUTINE int_deallocate
        SUBROUTINE triangulate(cmd,in,out,vorout) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR,C_CHAR
            IMPORT::triangulateio_tydef
            IMPLICIT NONE            
            CHARACTER(KIND=C_CHAR), INTENT(IN)::cmd
            TYPE(triangulateio_tydef)::in,out,vorout
        END SUBROUTINE triangulate


    end interface
        
        
        
    contains
    
    subroutine triangulateio_tri2mesh_lgy(this,tridata)
        implicit none
        class(triangle_tydef)::this
        type(triangulateio_tydef)::tridata
        real(c_double),pointer::ar1(:)
        integer(c_int),pointer::iar1(:)
        integer::err,i,j,iar2(3),n1,n2,n3

        this.nnode=tridata.numberofpoints
        
        if(allocated(this.node)) deallocate(this.node)
        allocate(this.node(this.nnode),stat=err)
       
        call c_f_pointer(tridata.pointlist,ar1,[2*this.nnode])
        this.node.x=ar1(1:2*this.nnode:2)  
        this.node.y=ar1(2:2*this.nnode:2)
        this.node.number=[1:this.nnode]
        if(tridata.numberofpointattributes>0) then
            call c_f_pointer(tridata.pointattributelist,ar1,[tridata.numberofpointattributes*this.nnode])
            do i=1,this.nnode
                this.node(i).nat=tridata.numberofpointattributes
                this.node(i).at=ar1((i-1)*tridata.numberofpointattributes+1:i*tridata.numberofpointattributes)
            enddo
        endif
        if(C_associated(tridata.pointmarkerlist)) then
            call c_f_pointer(tridata.pointmarkerlist,iar1,[this.nnode])
            this.node.marker=iar1(1:this.nnode)
        endif

        !elt
        this.nelt=tridata.numberoftriangles
        if(allocated(this.element)) deallocate(this.element)
        allocate(this.element(this.nelt),stat=err)        
        call c_f_pointer(tridata.trianglelist,iar1,[tridata.numberofcorners*this.nelt])
        do i=1,this.nelt
            this.element(i).NNUM=tridata.numberofcorners
            this.element(i).node(1:tridata.numberofcorners)=iar1((i-1)*tridata.numberofcorners+1:i*tridata.numberofcorners)
        enddo
        if(tridata.numberoftriangleattributes>0) then
            call c_f_pointer(tridata.triangleattributelist,ar1,[tridata.numberoftriangleattributes*this.nelt])
            do i=1,this.nelt
                this.element(i).nat=tridata.numberoftriangleattributes
                this.element(i).at=ar1((i-1)*tridata.numberoftriangleattributes+1:i*tridata.numberoftriangleattributes)
            enddo
        endif

        !neigbor
        call c_f_pointer(tridata.neighborlist,iar1,[3*this.nelt])
        do i=1,this.nelt
            iar2=iar1((i-1)*3+1:i*3)
            this.element(i).adj(1:3)=iar2([3,1,2])
        enddo
        this.nedge=tridata.numberofedges        
        if(allocated(this.edge)) deallocate(this.edge)
        allocate(this.edge(this.nedge),stat=err)
        if(allocated(this.adjlist)) deallocate(this.adjlist)
        allocate(this.adjlist(this.nnode))
        call c_f_pointer(tridata.edgelist,iar1,[2*this.nedge])
        do i=1,this.nedge
            this.edge(i).v=iar1((i-1)*2+1:i*2)
            call addadjlist(this.adjlist,this.edge(i).v(1),this.edge(i).v(2),i)
        enddo
        if(C_associated(tridata.edgemarkerlist)) then
            call c_f_pointer(tridata.edgemarkerlist,iar1,[this.nedge])
            this.edge.marker=iar1
        endif
        this.ncedge=count(this.edge.marker/=0)
        if(allocated(this.cedge)) deallocate(this.cedge)
        allocate(this.cedge(this.ncedge))
        n1=0
        do j=1,this.nedge
            if(this.edge(j).marker==0) cycle
            n1=n1+1
            this.cedge(n1).v=this.edge(j).v
            this.cedge(n1).cl=this.edge(j).marker
            this.cedge(n1).edge=j
            this.edge(j).iscedge=n1
            !marker==-1 model outer bonndaries.
            if(this.edge(j).marker==-1) this.cedge(n1).cl=0                        
        enddo
               
        
        !initialize element.edge
        
        do i=1,this.nelt
            do j=1,3
                n1=this.element(i).node(j)
                n2=this.element(i).node(mod(j,3)+1)
                n3=v2edge(this.adjlist,n1,n2)
                this.element(i).edge(j)=n3
                if(this.edge(n3).e(1)==-1) then
                    this.edge(n3).e(1)=i
                else
                    this.edge(n3).e(2)=i
                endif
            enddo
        enddo 

    endsubroutine

    
    subroutine tetgenio_read_indata(tetdata,basefilename,filetype)
        implicit none
        !!class(triangle_tydef)::self
        type(tetgenio)::tetdata
        character(len=*),intent(in)::basefilename
        integer(tetgenbehavior_objecttype),intent(in)::filetype
        logical::isok
        character(8)::ext1

        !tetdata=tetgenio()

        select case(filetype)
        case(tetgenbehavior_POLY)
            isok=tetdata.load_poly(trim(adjustL(basefilename)))
            ext1='.poly'
        case(tetgenbehavior_NODES)
            isok=tetdata.load_node(trim(adjustL(basefilename)))
            ext1='.node'
        case(tetgenbehavior_OFF)
            isok=tetdata.load_off(trim(adjustL(basefilename)))
            ext1='.off'
        case(tetgenbehavior_PLY)
           isok=tetdata.load_ply(trim(adjustL(basefilename)))
           ext1='.ply'
        case(tetgenbehavior_STL)
            isok=tetdata.load_stl(trim(adjustL(basefilename)))
            ext1='.stl'
        case(tetgenbehavior_MEDIT)
            isok=tetdata.load_medit(trim(adjustL(basefilename)),1)
            ext1='.medit'
        case(tetgenbehavior_VTK)
            isok=tetdata.load_vtk(trim(adjustL(basefilename)))
            ext1='.vtk'
        case(tetgenbehavior_MESH)  
            isok=tetdata.load_tetmesh(trim(adjustL(basefilename)),tetgenbehavior_MESH)
            ext1='.tetmesh'
        case(tetgenbehavior_NEU_MESH)

        case default
            print *, 'error.sub=tetgenio_read_indata.no such filetype:',filetype
            error stop
        end select

        if(.not.isok) then
            print *, 'fail to load file:'//trim(basefilename)//trim(ext1)
        endif

    endsubroutine

    subroutine triangulateio_set_indata(self,tridata,polyfile)
        implicit none
        class(triangle_tydef)::self
        type(triangulateio_tydef)::tridata
        character(len=*),intent(in),optional::polyfile
        real(c_double),allocatable::ar1(:,:)
        integer(c_int),allocatable::iar1(:,:)
        integer::unit,i,j,n1,n2,n3,ismarker1,err
        integer::nread,nneed,nmax,unit1        
        parameter(nmax=100)
	    real(8)::linedata(nmax)


        if(present(polyfile)) then
            if(len_trim(adjustL(polyfile))>0) then
                unit1=10
                open(unit,file=polyfile,status='old')
                nneed=nmax

                !point data
                call strtoint(unit1,linedata,nmax,nread,nneed)
                tridata.numberofpoints=int(linedata(1))
                if(tridata.numberofpoints==0) then
                    print *, 'no point data in the polyfile.'
                    print *, 'please add point data in it.'
                    error stop  "sub=triangulateio_indata"
                endif
                tridata.numberofpointattributes=int(linedata(3))
                ismarker1=int(linedata(4))
                deallocate(ar1,stat=err)
                allocate(ar1(2+tridata.numberofpointattributes+ismarker1,tridata.numberofpoints),stat=err)
                do i=1,tridata.numberofpoints
                    call strtoint(unit1,linedata,nmax,nread,nneed)
                    n1=int(linedata(1))
                    ar1(:,n1)=linedata(2:nread)                    
                enddo
                call self.setvalue(tridata.pointlist,reshape(AR1(1:2,:),([2*tridata.numberofpoints])))
                if(tridata.numberofpointattributes>0) &
                call self.setvalue(tridata.pointattributelist,reshape(AR1(3:2+tridata.numberofpointattributes,:),([tridata.numberofpointattributes*tridata.numberofpoints])))
                if(ismarker1==1)  call self.setvalue(tridata.pointmarkerlist,reshape(int(AR1(nread,:)),([tridata.numberofpoints])))
                
                !seg data
                call strtoint(unit1,linedata,nmax,nread,nneed)
                tridata.numberofsegments=int(linedata(1))
                ismarker1=int(linedata(2))
                
                if(tridata.numberofsegments>0) then
                    deallocate(iar1,stat=err)
                    allocate(iar1(2+ismarker1,tridata.numberofsegments),stat=err)
                    do i=1,tridata.numberofsegments
                        call strtoint(unit1,linedata,nmax,nread,nneed)
                        n1=int(linedata(1))
                        iar1(:,n1)=linedata(2:nread) 
                    enddo
                    call self.setvalue(tridata.segmentlist,reshape(IAR1(1:2,:),([2*tridata.numberofsegments])))
                    if(ismarker1==1)  call self.setvalue(tridata.segmentmarkerlist,reshape(int(IAR1(nread,:)),([tridata.numberofsegments])))
                endif

                !hole
                call strtoint(unit1,linedata,nmax,nread,nneed)
                tridata.numberofholes=int(linedata(1))
                if(tridata.numberofholes>0) then
                    deallocate(ar1,stat=err)
                    allocate(ar1(2,tridata.numberofholes),stat=err)
                    do i=1,tridata.numberofholes
                        call strtoint(unit1,linedata,nmax,nread,nneed)
                        n1=int(linedata(1))
                        ar1(:,n1)=linedata(2:nread) 
                    enddo
                    call self.setvalue(tridata.holelist,reshape(AR1(1:2,:),([2*tridata.numberofholes])))
                endif

                !region
                call strtoint(unit1,linedata,nmax,nread,nneed)
                tridata.numberofregions=int(linedata(1))
                if(tridata.numberofregions>0) then
                    deallocate(ar1,stat=err)
                    allocate(ar1(4,tridata.numberofregions),stat=err)
                    do i=1,tridata.numberofregions
                        call strtoint(unit1,linedata,nmax,nread,nneed)
                        n1=int(linedata(1))
                        ar1(1:nread-1,n1)=linedata(2:nread)
                        if(nread<5) ar1(4,n1)=-1.d0 
                    enddo
                    call self.setvalue(tridata.regionlist,reshape(AR1(1:4,:),([4*tridata.numberofregions])))
                endif                

                close(unit1)

                return
            endif
        endif

       N1=SOILLAYER
        IF(N1>0) THEN
            N1=N1+1
        ENDIF
        tridata.numberofpoints=inpn
        tridata.numberofpointattributes=1+n1
        n2=3+n1
        deallocate(ar1,stat=err)
        allocate(ar1(n2,inpn),stat=err)
        
        DO I=1,INPN
            if(n1>0) then
                AR1(:,i)=[ ARR_T(I).X,ARR_T(I).Y,ARR_T(I).S,ARR_T(I).SOILDATA(0:SOILLAYER)]
            ELSE
                 AR1(:,i)=[ ARR_T(I).X,ARR_T(I).Y,ARR_T(I).S]
            ENDIF
        enddo
        call self.setvalue(tridata.pointlist,reshape(AR1(1:2,:),([2*inpn])))
        call self.setvalue(tridata.pointattributelist,reshape(AR1(3:n2,:),([(n2-2)*inpn])))
        !iar2=arr_t(1:inpn).marker
        call self.setvalue(tridata.pointmarkerlist,arr_t(1:inpn).marker)
        
        tridata.numberofsegments=NSEG
        deallocate(iar1,stat=err)
        ALLOCATE(IAR1(3,NSEG),stat=err)
        DO I=1,NSEG
            call seg(i).getparas(sv=n1,ev=n2,icl=n3)
            if(n3==0) n3=-1 !marker==-1 model outer bonndaries.
            IAR1(:,I)=[N1,N2,N3]
        ENDDO
        if(nseg>0) then
            call self.setvalue(tridata.segmentlist,reshape(iAR1(1:2,:),([2*nseg])))
            call self.setvalue(tridata.segmentmarkerlist,reshape(iAR1(3,:),([nseg])))
        endif

        !HOLE
        N1=COUNT(CSL(0:CLN).HOLE>0)
        tridata.numberofholes=n1
        if(n1>0) then
            deallocate(ar1,stat=err)
            allocate(ar1(2,n1),stat=err)
        
            N1=0
            DO I=1,CLN
                IF(CSL(I).HOLE==0) CYCLE  
                N1=N1+1          
                ar1(:,n1)=Find_Point_Inside_Polygon_2D(CSL(I).conpoint(1:2,:))           
            ENDDO
        
            call self.setvalue(tridata.holelist,reshape(AR1(1:2,:),([2*n1])))
        endif

        N1=COUNT(ZONE(0:ZNUM).FORMAT>0)
        IF(N1>0) THEN
            N1=0           
            DO I=0,ZNUM
                IF(ZONE(I).FORMAT==0) CYCLE
                N1=N1+ZONE(I).NUM
            ENDDO
            tridata.numberofregions=N1
            deallocate(ar1,stat=err)
            allocate(ar1(4,n1),stat=err)            
            N1=0
            DO I=0,ZNUM
                IF(ZONE(I).FORMAT==0) CYCLE
                DO J=1,ZONE(I).NUM
                    N1=N1+1
                    ar1(:,j)=[ZONE(I).POINT(1:2,J),real(I,8),ZONE(I).EAREA]
                ENDDO
            ENDDO
            call self.setvalue(tridata.regionlist,reshape(AR1,([4*n1])))          
        ENDIF




    endsubroutine




    subroutine triangulateio_triangulate(self,cmd,in,out,vorout)
        implicit none
        class(triangle_tydef)::self
        character(len=*,kind=c_char),intent(in)::cmd
        type(triangulateio_tydef),intent(in)::in
        type(triangulateio_tydef),intent(out)::out
        type(triangulateio_tydef),intent(out),optional::vorout

        
        call triangulate(cmd,in,out,vorout)
        


    endsubroutine

    subroutine triangulateio_free_array(tridata)
        implicit none
        !class(triangle_tydef)::self
        type(triangulateio_tydef)::tridata

        call real_deallocate(tridata.pointlist)
        call real_deallocate(tridata.pointattributelist)
        call real_deallocate(tridata.triangleattributelist)
        call real_deallocate(tridata.trianglearealist)
        call real_deallocate(tridata.holelist)
        call real_deallocate(tridata.regionlist)                
        call real_deallocate(tridata.normlist)

        call int_deallocate(tridata.pointmarkerlist)
        call int_deallocate(tridata.trianglelist)
        call int_deallocate(tridata.neighborlist)
        call int_deallocate(tridata.segmentlist)
        call int_deallocate(tridata.segmentmarkerlist)
        call int_deallocate(tridata.edgelist)
        call int_deallocate(tridata.edgemarkerlist)
    end

    subroutine triangulateio_real_array_getvalue(cptr,array,n)
        implicit none
        !class(triangle_tydef)::self
        type(C_PTR)::cptr
        integer(c_size_t),intent(in)::n
        real(c_double),intent(out),pointer::array(:)
   
        call c_f_pointer(cptr,array,[n])        

 
    end

    subroutine triangulateio_int_array_getvalue(cptr,array,n)
        implicit none
        !class(triangle_tydef)::self
        !type(triangulateio_tydef)::tridata
        type(C_PTR)::cptr
        integer(c_size_t)::n
        integer(c_int),intent(out),pointer::array(:)

        call c_f_pointer(cptr,array,[n])        

    end
    
    
    
    subroutine triangulateio_real_array_setvalue(cptr,array)
        implicit none
        type(C_PTR)::cptr
        real(c_double),intent(in),target,value::array(:)
        real(c_double),pointer::art1(:)
        integer(c_size_t)::n

        n=size(array)
        cptr=real_allocate(n)
        call c_f_pointer(cptr,art1,[n]) 
        art1=array
        !isloc1=C_associated(self.tridata.pointlist,c_loc(art1))
        !isloc1=C_associated(self.tridata.pointlist,c_loc(array))
        
    end

    subroutine triangulateio_int_array_setvalue(cptr,array)
        implicit none
        !class(triangle_tydef)::self
        !type(triangulateio_tydef)::tridata
        !character(len=*),intent(in)::name
        type(C_PTR)::cptr
        integer(c_int),intent(in),target,value::array(:)
        integer(c_int),pointer::art1(:)
        integer(c_size_t)::n

        n=size(array)
        cptr=int_allocate(n)
        call c_f_pointer(cptr,art1,[n]) 
        art1=array

    end

    subroutine get_triangle_meshdata(this,element,node,edge,adjlist,cedge)
        implicit none
        class(triangle_tydef)::this
        type(point_tydef),allocatable,optional::node(:)
        type(element_tydef),allocatable,optional::element(:)
        type(edge_tydef),allocatable,optional::edge(:)
        type(adjlist_tydef),allocatable,optional::adjlist(:)
        type(constrained_edge_tydef),allocatable,optional::cedge(:)
        if(present(element)) element=this.element
        if(present(node)) node=this.node
        if(present(edge)) edge=this.edge
        if(present(adjlist)) adjlist=this.adjlist
        if(present(cedge)) cedge=this.cedge
    endsubroutine
    
    subroutine read_triangle_file(this,filepath,fext)
        use dflib
        USE IFPORT
        implicit none
        class(triangle_tydef)::this
        character(len=*),intent(in),optional::filepath
        character(len=*),optional,intent(in)::fext(:)
        CHARACTER(3)        drive
	    CHARACTER(512)      dir
	    CHARACTER(512)      name,file1
	    CHARACTER(16)      ext
        CHARACTER(len=:),allocatable::ext1(:)
	    integer(4)::length,msg
        logical::isexist
        integer::unit1,i,j,hasread,offset1=0

        integer::na1=0,ismarker1=0,n1=0,n2,nelt1,nmax
        integer::nread,nneed,nnode1,n3
        integer::iar1(10)
        
        parameter(nmax=100)
	    !parameter(maxset=100)
	
	    real(8)::linedata(nmax),ar1(nmax)
	    !character(32)::set(maxset)



        if(this.islib==1) then
            call this.tri2mesh_lgy(this.out)
            return
        endif


        nneed=nmax

        if(present(fext)) then
            ext1=fext
        else
            allocate(character(16)::ext1(4))
            ext1(1:4)=['node','ele','edge','neigh'] !,'face','t2e','t2f','f2e','v.node','v.edge','v.face','v.cell']
        endif
        
        
        do i=1,size(ext1)
            file1=trim(filepath)//'.'//trim(adjustl(ext1(i)))
            inquire(file=file1,exist=isexist)
            if(isexist) then
                unit1=10
                open(unit=unit1,file=file1,status='old')
                hasread=1
                
                call strtoint(unit1,linedata,nmax,nread,nneed)
                select case(trim(adjustl(ext1(i))))
                case('node')
                    
                    !read(unit,*) nelt,ndim_tg,na1,ismarker1
                    this.nnode=int(linedata(1))
                    n1=sum(linedata(2:4))
                    na1=int(linedata(3))
                    ismarker1=int(linedata(4))
                    if(allocated(this.node)) deallocate(this.node)
                    allocate(this.node(this.nnode))
                    !if(na1>0) allocate(node_tg.at(na1))
                    offset1=0
                    if(index(trim(this.cmd),'u')/=0) offset1=1
                    do j=1,this.nnode
                        read(unit1,*) n2,ar1(1:n1)            
                        this.node(n2).x=ar1(1);this.node(n2).y=ar1(2);
                        if(na1>0) this.node(n2).at=ar1(3:2+na1)
                        
                        if(soillayer>0) then
                            allocate(this.node(n2).elevation(0:soillayer))
                            this.node(n2).elevation=ar1(offset1+3:3+soillayer+offset1)
                            this.node(n2).havesoildata=int(ar1(4+soillayer+offset1))
                        endif
                        this.node(n2).nat=na1
                        this.node(n2).number=n2
                        if(ismarker1>0) then
                            this.node(n2).marker=int(ar1(n1))
                        endif
                    end do
                    
                    !call read_tetgen_node(unit1,trim(adjustl(ext1(i))),this.node,this.nnode)
                    !natr=node_tg(1).na
                case('ele')
                    !call read_tetgen_element(unit1,trim(adjustl(ext1(i))),this.element,this.nelt)
                    !ntetnode=element_tg(1).nnode
                    this.nelt=int(linedata(1))
                    if(allocated(this.element)) deallocate(this.element)
                    allocate(this.element(this.nelt))
                    nnode1=int(linedata(2))
                    ismarker1=int(linedata(3))
                    n1=nnode1+ismarker1
                    do j=1,this.nelt
                        read(unit1,*) n2,ar1(1:n1)
                        this.element(n2).nnum=nnode1
                        this.element(n2).node(1:nnode1)=int(ar1(1:nnode1))
                        if(ismarker1>0) this.element(n2).at=ar1(nnode1+1:n1)
                        this.element(n2).nat=ismarker1
                    end do                    
                case('neigh')
                    !call read_tetgen_element(unit1,trim(adjustl(ext1(i))),neigh_tg,nneigh_tg)
                    n1=3
                    do j=1,int(linedata(1))
                        read(unit1,*) n2,ar1(1:3)
                        this.element(n2).adj(1:3)=int(ar1([3,1,2]))
                        !
                    end do  
                !case('face')
                !    call read_tetgen_element(unit1,trim(adjustl(ext1(i))),face_tg,nface_tg)
                case('edge')
                    this.nedge=int(linedata(1))
                    if(allocated(this.edge)) deallocate(this.edge)
                    allocate(this.edge(this.nedge))
                    ismarker1=int(linedata(2))
                    n1=2+ismarker1
                    if(allocated(this.adjlist)) deallocate(this.adjlist)
                    allocate(this.adjlist(this.nnode))
                    do j=1,this.nedge
                        read(unit1,*) n2,ar1(1:n1)
                        this.edge(n2).v=int(ar1(1:2))
                        if(ismarker1>0) this.edge(n2).marker=int(ar1(n1))
                        call addadjlist(this.adjlist,int(ar1(1)),int(ar1(2)),n2)
                    end do 
                    this.ncedge=count(this.edge.marker/=0)
                    if(allocated(this.cedge)) deallocate(this.cedge)
                    allocate(this.cedge(this.ncedge))
                    n1=0
                    do j=1,this.nedge
                        if(this.edge(j).marker==0) cycle
                        n1=n1+1
                        this.cedge(n1).v=this.edge(j).v
                        this.cedge(n1).cl=this.edge(j).marker
                        this.cedge(n1).edge=j
                        this.edge(j).iscedge=n1
                        !marker==-1 model outer bonndaries.
                        if(this.edge(j).marker==-1) this.cedge(n1).cl=0                        
                    enddo
                    
                        
                    !call read_tetgen_element(unit1,trim(adjustl(ext1(i))),edge_tg,nedge_tg)
                !case('v.node')
                !    call read_tetgen_node(unit1,trim(adjustl(ext1(i))),vnode_tg,nvnode_tg)
                !case('v.face')
                !    call read_tetgen_element(unit1,trim(adjustl(ext1(i))),vface_tg,nvf_tg)                                      
                !case('v.edge')
                !    call read_tetgen_element(unit1,trim(adjustl(ext1(i))),vedge_tg,nve_tg)
                !case('v.cell')
                !    call read_tetgen_element(unit1,trim(adjustl(ext1(i))),vcell_tg,nvc_tg)  
                !case('t2e')
                !    call read_adj_table(unit1,trim(adjustl(ext1(i))),t2e,nelt_tg) 
                !case('t2f')
                !    call read_adj_table(unit1,trim(adjustl(ext1(i))),t2f,nelt_tg)  
                !case('f2e')
                !    call read_adj_table(unit1,trim(adjustl(ext1(i))),f2e,nface_tg)                       
                case default
                    hasread=0
                    print *, 'No such file type=',trim(ext1(i))
                end select
                if(hasread>0) print *, 'Done in reading file=',trim(file1)
                close(unit1)
            else
                print *, 'file is not exist and skipped. file=',trim(file1)
            endif
        enddo    
        
        !initialize element.edge
        
        do i=1,this.nelt
            do j=1,3
                n1=this.element(i).node(j)
                n2=this.element(i).node(mod(j,3)+1)
                n3=v2edge(this.adjlist,n1,n2)
                this.element(i).edge(j)=n3
                if(this.edge(n3).e(1)==-1) then
                    this.edge(n3).e(1)=i
                else
                    this.edge(n3).e(2)=i
                endif
            enddo
        enddo
        


    end subroutine
    
    
    
    subroutine triangle_option_read(this,unit)
        use dflib
        implicit none
        class(triangle_tydef)::this
        integer,intent(in)::unit
        integer::oldcolor,i,n1
        logical::isexist
        
        !type(triangle_tydef)::ftri
        !
        !call ftri.triangulate()
        !stop
        
        
        print *,'Reading mesh_by_triangle data...'
        oldcolor = SETTEXTCOLOR(INT2(10))
        call this.help(this.helpstring)        
        oldcolor = SETTEXTCOLOR(INT2(15))
        do i=1, pro_num
            select case(property(i).name)
            case('method')
                this.method=int(property(i).value)
            case('path')
	             this.path=trim(adjustl(property(i).cvalue))
            case('cmd')
	             this.cmd=trim(adjustl(property(i).cvalue)) 
            case('polyfile','poly')
	             this.polyfile=trim(adjustl(property(i).cvalue))    
            case('islib')
	             this.islib=int(property(i).value) 
            case('dim')
	             this.dim=int(property(i).value)                                                     
            case default
	            call Err_msg(property(i).name)
            end select
        end do 
        
        if(len_trim(this.polyfile)==0) this.polyfile=trim(adjustl(path_name))//'.poly'
        if(len_trim(this.path)==0) then
            if(this.dim==2) then
                this.path='.\triangle.exe'
            else
                this.path='.\tetgen.exe'
            endif
            inquire(file=this.path,exist=isexist)
            if(.not.isexist) then
                if(this.dim==2) then
                    this.path='triangle'
                else
                    this.path='tetgen'
                endif
            endif
        endif
        if(this.dim==2) then
        !为保持与本程序的数据结构兼容,参数'en'为必输参数         
            if(index(trim(this.cmd),'e')==0) this.cmd=trim(this.cmd)//'e'
            if(index(trim(this.cmd),'n')==0) this.cmd=trim(this.cmd)//'n'
            if(index(trim(this.cmd),'g')==0) this.cmd=trim(this.cmd)//'g'
            if(index(trim(this.cmd),'Q')==0.and.this.islib==1) this.cmd=trim(this.cmd)//'Q'     
            n1=index(trim(this.cmd),'u')
            if(n1==0.and.isnorefined/=1) this.cmd=trim(this.cmd)//'u'
            if(n1/=0.and.isnorefined==1) then
                this.cmd=this.cmd(1:n1-1)//this.cmd(n1+1:len_trim(this.cmd))
            endif
        else
            if(index(trim(this.cmd),'e')==0) this.cmd=trim(this.cmd)//'e'
            if(index(trim(this.cmd),'nn')==0.or.index(trim(this.cmd),'n')==0) this.cmd=trim(this.cmd)//'nn'
            if(index(trim(this.cmd),'f')==0) this.cmd=trim(this.cmd)//'f'
            if(index(trim(this.cmd),'g')==0) this.cmd=trim(this.cmd)//'g'
            if(index(trim(this.cmd),'k')==0) this.cmd=trim(this.cmd)//'k'
            !if(index(trim(this.cmd),'m')==0) this.cmd=trim(this.cmd)//'m'
            if(index(trim(this.cmd),'Q')==0.and.this.islib==1) this.cmd=trim(this.cmd)//'Q'
        endif
    endsubroutine
    
    subroutine execute_triangle_cmd(this,filepath)
    
        use ifport
        class(triangle_tydef)::this
        character(len=*),intent(in)::filepath
        integer::unit1,i,err,n1,n2
        INTEGER :: CSTAT, ESTAT,result
        CHARACTER(100) :: CMSG
        character(1024)::file1,cmd1
        character(16)::ftype1(4)
        character(:),allocatable::filepath1
        real(c_double),pointer::ar1(:)
        type(triangulateio_tydef)::out1
        character(256)::basefilename
		integer(4)::length
		character(1024)::nme
		CHARACTER(3)        drive
		CHARACTER(1024)      dir
		CHARACTER(1024)      name
		CHARACTER(256)      ext
        type(SWIGTYPE_p_double)::dptr
        real(c_double),pointer::ar2(:)
        logical::tof1
        character(16)::ch1
        if(this.method<0) return

        if(this.dim==3) then
            length = SPLITPATHQQ(filepath, drive, dir, basefilename, ext)
            this.tetin=tetgenio()
            call this.tetin.set_basefilename(trim(adjustl(basefilename)))
            basefilename=trim(adjustl(basefilename))//C_NULL_CHAR            
            call tetgenio_read_indata(this.tetin,basefilename,tetgenbehavior_POLY)
            if(index(trim(this.cmd),'m')/=0) then
                !让pointmtrlist=pointattributes(1)
                dptr=this.tetin.get_pointattributelist()
                n1=this.tetin.get_numberofpointattributes()
                call c_f_pointer(dptr.swigdata.cptr,ar1,[n1*this.tetin.get_numberofpoints()])
                call this.tetin.set_numberofpointmtrs(1)
                ar2=ar1(1::n1)
                dptr=this.tetin.get_pointmtrlist()
                dptr.swigdata.cptr=c_loc(ar2)
                call this.tetin.set_pointmtrlist(dptr)
            endif
            if(index(trim(this.cmd),'r')/=0) then
                call tetgenio_read_indata(this.tetin,basefilename,tetgenbehavior_MESH)
            endif
            call this.meshing()
            if(this.islib==1) then
                
                !n1=len_trim(basefilename)-1 !remove '/0'
                !basefilename=basefilename(1:n1)
                !n2=1
                !if(n1>2) then
                !    tof1=ichar(basefilename(n1:n1))>=ichar('0').and.ichar(basefilename(n1:n1))<=ichar('9')
                !    tof1=tof1.and.basefilename(n1-1:n1-1)=='.'
                !    if(tof1) then
                !        read(basefilename(n1:n1),*) n2
                !        n2=n2+1
                !        basefilename=basefilename(1:n1-2)
                !    endif                
                !endif
                !write(ch1,'(i2)') n2
                !basefilename=trim(adjustl(basefilename))//'.'//trim(adjustl(ch1))//C_NULL_CHAR
                !call this.tetout.save_poly(basefilename)
                !call this.tetout.save_nodes(basefilename)
                !call this.tetout.save_elements(basefilename)
                !call this.tetout.save_neighbors(basefilename)
                !call this.tetout.save_faces(basefilename)
                !call this.tetout.save_edges(basefilename)
                call this.tetout.clean_memory()
                !call this.tetin.clean_memory() ！会出错,因为有些单元数组的空间不是用c++的new分配的.比如上述pointmtrlist.
            endif 
            
            stop
            !return

        endif

        !write poly file
        if(this.method==0.or.this.method==2.or.this.method==4.or.this.method==5) then
            !outpoly
            call this.outpoly()
            if(this.method==2) stop            
        endif
        !call triangle
        if(this.method==0.or.this.method==1.or.this.method==4.or.this.method==5) then
            if(this.islib/=0) THEN
                if(this.method/=1) then
                    call this.set_indata(this.in)
                else
                    call this.set_indata(this.in,this.polyfile)
                endif                
            endif
            call this.meshing() 
            if(this.method==4) stop
        endif
        
        if(index(trim(this.cmd),'I')==0) then
            filepath1=trim(filepath)//'.1'
        else
            filepath1=trim(filepath)
        endif


        if(this.method==5) then
            !refine mesh by generated s random size field based on pdf_sample
            if(this.islib==0) then
                call this.readmesh(filepath1,['node'])
                example_pdf.nval=this.nnode
                ar2=example_pdf.sample()
                do i=1,this.nnode
                    this.node(i).at(1)=ar2(i)
                enddo
                call this.out_node(filepath1)
                if(index(trim(this.cmd),'I')==0) then
                    cmd1='triangle '//'-rqpueng '//filepath1
                    filepath1=trim(filepath)//'.2'
                else
                    cmd1='triangle '//'-rpqueng '//filepath1
                    filepath1=trim(filepath)//'.1'
                endif
                call this.meshing(trim(cmd1))

            else
                this.in=this.out
                example_pdf.nval=this.in.numberofpoints
                if(this.in.numberofpointattributes>0) then 
                    call c_f_pointer(this.in.pointattributelist,ar1,[this.in.numberofpointattributes*this.in.numberofpoints])
                    ar1(1:size(ar1):this.in.numberofpointattributes)=example_pdf.sample()
                else
                    this.in.numberofpointattributes=1
                    ar1=example_pdf.sample()
                    call this.setvalue(this.in.pointattributelist,ar1)
                endif
                !call this.free(this.out)
                this.out=out1
                this.cmd='-rpquengQ'
                call this.meshing()
            endif
            
        endif
           
        call this.readmesh(filepath1)

        deallocate(filepath1,stat=err)
        if(this.islib) then
            call this.free(this.in)
            call this.free(this.out)
            call this.free(this.vorout)
        endif
        
!10  format(i7,3(e24.16,','),e24.16)        
    end subroutine  
    subroutine triangle_write_node_file(this,filepath,fext)
        implicit none
        class(triangle_tydef)::this
        character(len=*),intent(in)::filepath
        character(len=*),intent(in),optional::fext
        character(len=:),allocatable::file1
        character(16)::fext1
        integer::unit=10,nat1,ismarker1,i

        if(present(fext)) then
            fext1=trim(adjustl(fext))
        else
            fext1='node'
        endif
        file1=trim(adjustl(filepath))//'.'//trim(adjustl(fext1))

        open(unit,file=file1,status='replace')
        nat1=this.node(1).nat;ismarker1=1
        write(unit,10) this.nnode,2,nat1,ismarker1 
        do i=1,this.nnode
            write(unit,20),i,this.node(i).x,this.node(i).y,this.node(i).at,this.node(i).marker
        enddo
        close(unit)

10  format(4(i,x))
20  format(i,x,2(f,x),<nat1>(f,x),<ismarker1>(i,x))
    endsubroutine
    subroutine meshbytriangle(this,cmd)
        class(triangle_tydef)::this
        character(*),intent(in),optional::cmd
        INTEGER :: CSTAT, ESTAT
        CHARACTER(100) :: CMSG
        character(len=:),allocatable:: cmd1
        
        cmd1=' '
        if(present(cmd)) then
            cmd1=trim(adjustl(cmd))
        endif

        if(len_trim(adjustL(cmd1))==0) then
            cmd1=trim(this.path)//' '//trim(this.cmd)//' '//trim(this.polyfile)
            cmd1=trim(adjustl(cmd1))
        endif
        !result = SYSTEMQQ ("tetgen -wefnNI pack_six_cube_poly.node")
        if(this.method==4.or.this.islib==0) then
            CALL EXECUTE_COMMAND_LINE (cmd1, EXITSTAT=ESTAT,CMDSTAT=CSTAT, CMDMSG=CMSG)
            IF (CSTAT > 0) THEN
                PRINT *, "Triangle Command execution failed with error ", TRIM(CMSG)
                pause
            ELSE IF (CSTAT < 0) THEN
                PRINT *, "Triangle Command execution not supported"
                pause
            ELSE
                PRINT *, "Command completed with status ", ESTAT
            END IF
        else
            if(this.dim==2) then
                call triangulate(trim(this.cmd)//C_NULL_CHAR,this.in,this.out,this.vorout)
            else
                this.tetout=tetgenio()
                call tetrahedralize(trim(this.cmd)//C_NULL_CHAR,this.tetin,this.tetout)                
            endif
        endif
        
    endsubroutine
    
    
    subroutine write_poly_file(this,isbgmesh)
        implicit none
!!2D POLY FILE FORMAT        
!!First line: <# of vertices> <dimension (must be 2)> <# of attributes> <# of boundary markers (0 or 1)>
!!Following lines: <vertex #> <x> <y> [attributes] [boundary marker]
!!One line: <# of segments> <# of boundary markers (0 or 1)>
!!Following lines: <segment #> <endpoint> <endpoint> [boundary marker]
!!One line: <# of holes>
!!Following lines: <hole #> <x> <y>
!!Optional line: <# of regional attributes and/or area constraints>
!!Optional following lines: <region #> <x> <y> <attribute> <maximum area>
        class(triangle_tydef)::this
        integer,optional,intent(in)::isbgmesh
        !character(len=*),intent(in)::path_name
        integer::len1,UNIT,I,N1,N2,N3,J,isbgmesh1=0,nnode1,nseg1,n4,offset1=0
        REAL(8)::PT(2)
        character(1024)::outfile
        
        !outfile=trim(path_name)//'.poly'
        !if(len_trim(this.polyfile)==0.or.this.method==0) then
        !    this.polyfile=trim(outfile)
        !endif
        isbgmesh1=0
        if(present(isbgmesh)) isbgmesh1=isbgmesh
        outfile=trim(this.polyfile)
        nnode1=inpn;nseg1=nseg
        if(isbgmesh1/=0) then
            !outfile=trim(path_name)//'_bg.poly'
            nseg1=0
            DO I=1,NSEG
                call seg(i).getparas(sv=n1,ev=n2,icl=n3)
                if(arr_t(n1).havesoildata==2.and.arr_t(n2).havesoildata==2) then
                    nseg1=nseg1+1
                    call seg(i).setparas(isbg=nseg1)
                endif          
            ENDDO             
            nnode1=0
            do i=1,inpn
                if(ARR_T(i).havesoildata==2) then
                    nnode1=nnode1+1
                    ARR_T(i).bgnum=nnode1
                endif
            enddo
        endif
        
        UNIT=50
        OPEN(UNIT=UNIT,FILE=outfile,STATUS='REPLACE')
        len1=len_trim(POLY2D_FILE_FORMAT)
	    WRITE(UNIT,10) trim(POLY2D_FILE_FORMAT)
        
        N1=SOILLAYER
        IF(N1>0) THEN
            N1=N1+1
        ENDIF
        offset1=0
        if(index(trim(this.cmd),'u')/=0.and.isbgmesh1==0) offset1=1 !size info in the first attributes
        
        WRITE(UNIT,20) nnode1,N1+1+offset1
        if(offset1==0) then
            DO I=1,INPN            
                if(isbgmesh1/=0) then
                    if(ARR_T(I).havesoildata/=2) cycle
                    WRITE(UNIT,30) ARR_T(I).bgnum,ARR_T(I).X,ARR_T(I).Y,ARR_T(I).SOILDATA,ARR_T(I).havesoildata,0
                else
                    if(ARR_T(I).havesoildata>0) THEN
                        WRITE(UNIT,30) I,ARR_T(I).X,ARR_T(I).Y,ARR_T(I).SOILDATA,ARR_T(I).havesoildata,0
                    ELSE
                        WRITE(UNIT,31) I,ARR_T(I).X,ARR_T(I).Y,ARR_T(I).havesoildata,0
                    ENDIF
                endif
            ENDDO
        else
            DO I=1,INPN            
                if(isbgmesh1/=0) then
                    if(ARR_T(I).havesoildata/=2) cycle
                    WRITE(UNIT,30) ARR_T(I).bgnum,ARR_T(I).X,ARR_T(I).Y,ARR_T(I).SOILDATA,ARR_T(I).havesoildata,0 !bgmesh no need to refine.
                else
                    if(ARR_T(I).havesoildata>0) THEN
                        WRITE(UNIT,32) I,ARR_T(I).X,ARR_T(I).Y,ARR_T(I).S,ARR_T(I).SOILDATA,ARR_T(I).havesoildata,0
                    ELSE
                        WRITE(UNIT,33) I,ARR_T(I).X,ARR_T(I).Y,ARR_T(I).S,ARR_T(I).havesoildata,0
                    ENDIF
                endif
            ENDDO            
        endif
        
        
        WRITE(UNIT,40) NSEG1
        
        if(isbgmesh1/=0) then
            DO I=1,NSEG
                call seg(i).getparas(sv=n1,ev=n2,icl=n3,isbg=n4)
                if(n4>0) then
                    if(n3==0) n3=-1 !marker==-1 model outer bonndaries.
                    WRITE(UNIT,41) n4,arr_t(N1).bgnum,arr_t(N2).bgnum,N3
                endif
            ENDDO
            !HOLE
            WRITE(UNIT,'(I)') 0
        else
            DO I=1,NSEG
                call seg(i).getparas(sv=n1,ev=n2,icl=n3)
                if(n3==0) n3=-1 !marker==-1 model outer bonndaries.
                WRITE(UNIT,41) I,N1,N2,N3
            ENDDO
        endif
        
        if(isbgmesh1==0) then
        
            !HOLE
            N1=COUNT(CSL(0:CLN).HOLE>0)
            WRITE(UNIT,'(I)') N1
            N1=0
            DO I=1,CLN
                IF(CSL(I).HOLE==0) CYCLE            
                PT=Find_Point_Inside_Polygon_2D(CSL(I).conpoint(1:2,:))
                N1=N1+1
                WRITE(UNIT,50) N1,PT
            ENDDO
        
            N1=COUNT(ZONE(0:ZNUM).FORMAT>0)
            IF(N1>0) THEN
                N1=0           
                DO I=0,ZNUM
                    IF(ZONE(I).FORMAT==0) CYCLE
                    N1=N1+ZONE(I).NUM
                ENDDO
                WRITE(UNIT,'(I)') N1
                N1=0
                DO I=0,ZNUM
                    IF(ZONE(I).FORMAT==0) CYCLE
                    DO J=1,ZONE(I).NUM
                        N1=N1+1
                        WRITE(UNIT,60) N1,ZONE(I).POINT(1:2,J),I
                    ENDDO
                ENDDO          
            ENDIF
        endif
        
        CLOSE(UNIT)
        
10      FORMAT(A<len1>)  
20      FORMAT(I7,1X,'2',1X,I3,1X,'1')
30      FORMAT(I7,1X,<2+N1>(F15.7,1X),2(I3,1X))
31      FORMAT(I7,1X,<2>(F15.7,1X),<N1>('-999.0',1X),2(I3,1X)) 
32      FORMAT(I7,1X,<2+N1+1>(F15.7,1X),2(I3,1X))
33      FORMAT(I7,1X,<3>(F15.7,1X),2(I3,1X))         
40      FORMAT(I7,1X,'1')  
41      FORMAT(4(I7,1X)) 
50      FORMAT(I7,2(F15.7,1X))
60      FORMAT(I7,2(F15.7,1X),I)       
    end subroutine
    
    subroutine gen_backgroundmesh()
        implicit none
        logical::isexist
        
        if(len_trim(bgmesh.path)==0) then
            bgmesh.path='.\triangle.exe'
            inquire(file=bgmesh.path,exist=isexist)
            if(.not.isexist) then
                bgmesh.path='triangle'
            endif
        endif
        bgmesh.cmd='-pencgI'
        bgmesh.polyfile=trim(path_name)//'_bgm.poly'
        
        call bgmesh.outpoly(isbgmesh=1)
        call bgmesh.meshing()
        call bgmesh.readmesh(trim(path_name)//'_bgm')
        call bgmesh.set_searchzone(bgmesh.searchzone,bgmesh.nszone,(Xmax-Xmin)/xyscale,0.d0,(Ymax-Ymin)/xyscale,0.d0, &
            BGMESH.NODE,bgmesh.element)
 
        isbgready=.true.
    end subroutine
    
    function interpolate_attrib_bgmesh(this,Pt) result(attrib)
        implicit none
        class(triangle_tydef)::this
        real(8),intent(in)::pt(2)
        real(8),allocatable::attrib(:)
        real(8)::tri1(2,3),xsi1(3),ra1(3)
        integer::i,ielt1,nat1,j
        

        
        if(.not.isbgready) call gen_backgroundmesh()
        
        nat1=this.node(1).nat !all nodes nat is tha same
        if(nat1<1) return        
        ielt1=POINTlOC_2D(PT,-1,this.searchzone,this.node,this.element)
        if(ielt1<1) error stop 'Cannot locate the Point in the mesh.fun=interpolate_attrib_bgmesh'
        do i=1,3
            tri1(1,i)=this.node(this.element(ielt1).node(i)).x
            tri1(2,i)=this.node(this.element(ielt1).node(i)).y
        enddo
        xsi1= trishafun( tri1, pt(1), pt(2) )
        allocate(attrib(nat1))
        do i=1,nat1
            do j=1,3
                ra1(j)=this.node(this.element(ielt1).node(j)).at(i)
            enddo
            attrib(i)=dot_product(ra1,xsi1)
        enddo
    end function
    
    function trishafun(tri,x,y) result(shafun)
    implicit none
    real(8),intent(in)::tri(2,3),x,y
    real(8)::shafun(3)
    real(8)::a(3),b(3),c(3),t1
    integer::i
    
    
    a(1)=tri(1,2)*tri(2,3)-tri(1,3)*tri(2,2)
    a(2)=tri(1,3)*tri(2,1)-tri(1,1)*tri(2,3)
    a(3)=tri(1,1)*tri(2,2)-tri(1,2)*tri(2,1)
    b(1)=tri(2,2)-tri(2,3)
    b(2)=tri(2,3)-tri(2,1)
    b(3)=tri(2,1)-tri(2,2)
    c(1)=tri(1,3)-tri(1,2)
    c(2)=tri(1,1)-tri(1,3)
    c(3)=tri(1,2)-tri(1,1)
    t1=sum(a)
    if(abs(t1)<1d-14) then
        print *, 'the element area is zero. function=trishafun'
        stop
    endif
    shafun=(a+b*x+c*y)/t1
    return
endfunction
    
    
subroutine triangle_barycentric_2d ( t, p, xsi )

!*****************************************************************************80
!
!! TRIANGLE_BARYCENTRIC_2D finds the barycentric coordinates of a point in 2D.
!
!  Discussion:
!
!    The barycentric coordinate of point P related to vertex A can be
!    interpreted as the ratio of the area of the triangle with 
!    vertex A replaced by vertex P to the area of the original 
!    triangle.
!
!    This routine assumes that the triangle vertices are given in
!    counter clockwise order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!    The vertices should be given in counter clockwise order.
!
!    Input, real ( kind = 8 ) P(2), the point to be checked.
!
!    Output, real ( kind = 8 ) XSI(3), the barycentric coordinates of P
!    with respect to the triangle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  integer ( kind = 4 ), parameter :: rhs_num = 1

  real ( kind = 8 ) a(dim_num,dim_num+rhs_num)
  integer ( kind = 4 ) info
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) t(dim_num,3)
  real ( kind = 8 ) xsi(dim_num+1)
!
!  Set up the linear system
!
!    ( X2-X1  X3-X1 ) XSI(1)  = X-X1
!    ( Y2-Y1  Y3-Y1 ) XSI(2)    Y-Y1
!
!  which is satisfied by the barycentric coordinates of P.
!
  a(1,1) = t(1,2) - t(1,1)
  a(1,2) = t(1,3) - t(1,1)
  a(1,3) = p(1)   - t(1,1)

  a(2,1) = t(2,2) - t(2,1)
  a(2,2) = t(2,3) - t(2,1)
  a(2,3) = p(2)   - t(2,1)
!
!  Solve the linear system.
!
  call r8mat_solve ( dim_num, rhs_num, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TRIANGLE_BARYCENTRIC_2D - Fatal error!'
    write ( *, '(a)' ) '  The linear system is singular.'
    write ( *, '(a)' ) '  The input data does not form a proper triangle.'
    stop 1
  end if

  xsi(1) = a(1,3)
  xsi(2) = a(2,3)
  xsi(3) = 1.0D+00 - xsi(1) - xsi(2)

  return
end

subroutine r8mat_solve ( n, rhs_num, a, info )

!*****************************************************************************80
!
!! R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, integer ( kind = 4 ) RHS_NUM, the number of right hand sides.  
!    RHS_NUM must be at least 0.
!
!    Input/output, real ( kind = 8 ) A(N,N+rhs_num), contains in rows and
!    columns 1 to N the coefficient matrix, and in columns N+1 through
!    N+rhs_num, the right hand sides.  On output, the coefficient matrix
!    area has been destroyed, while the right hand sides have
!    been overwritten with the corresponding solutions.
!
!    Output, integer ( kind = 4 ) INFO, singularity flag.
!    0, the matrix was not singular, the solutions were computed;
!    J, factorization failed on step J, and the solutions could not
!    be computed.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) rhs_num

  real ( kind = 8 ) a(n,n+rhs_num)
  real ( kind = 8 ) apivot
  real ( kind = 8 ) factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot
  integer ( kind = 4 ) j

  info = 0

  do j = 1, n
!
!  Choose a pivot row.
!
    ipivot = j
    apivot = a(j,j)

    do i = j+1, n
      if ( abs ( apivot ) < abs ( a(i,j) ) ) then
        apivot = a(i,j)
        ipivot = i
      end if
    end do

    if ( apivot == 0.0D+00 ) then
      info = j
      return
    end if
!
!  Interchange.
!
    do i = 1, n + rhs_num
      call r8_swap ( a(ipivot,i), a(j,i) )
    end do
!
!  A(J,J) becomes 1.
!
    a(j,j) = 1.0D+00
    a(j,j+1:n+rhs_num) = a(j,j+1:n+rhs_num) / apivot
!
!  A(I,J) becomes 0.
!
    do i = 1, n

      if ( i /= j ) then

        factor = a(i,j)
        a(i,j) = 0.0D+00
        a(i,j+1:n+rhs_num) = a(i,j+1:n+rhs_num) - factor * a(j,j+1:n+rhs_num)

      end if

    end do

  end do

  return
end
  
subroutine r8_swap ( x, y )

!*****************************************************************************80
!
!! R8_SWAP switches two R8's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  implicit none

  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  z = x
  x = y
  y = z

  return
end    
    
end module 