module triangle_io
    use ds_t,only:arr_t,inpn
    use meshds,only:SEG,NSEG,CSL,CLN,Find_Point_Inside_Polygon_2D,ZONE,ZNUM, &
        & write_help,readline,property,pro_num,Err_msg,path_name,title, &
        & point_tydef,element_tydef,Edge_tydef,adjlist_tydef,addadjlist,v2edge,&
        & constrained_edge_tydef,soillayer,strtoint,SETUP_SEARCH_ZONE_2D,POINTlOC_2D,&
        & xmin,xmax,ymin,ymax,xyscale,SEARCHZONE_TYDEF,isnorefined
    implicit none
    private
    public::libtriangle,bgmesh

    logical::isbgready=.false.
    
    type triangle_tydef
        integer::method=-1,nnode=0,nelt=0,nedge=0,ncedge=0,NSZONE=0        
        character(1024)::path='',cmd='-penI',polyfile=''
        type(point_tydef),allocatable::node(:)
        type(element_tydef),allocatable::element(:)
        type(edge_tydef),allocatable::edge(:)
        type(adjlist_tydef),allocatable::adjlist(:)
        type(constrained_edge_tydef),allocatable::cedge(:)
        TYPE(SEARCHZONE_TYDEF),ALLOCATABLE::SEARCHZONE(:)
        character(4096):: helpstring= &
            &"相对于程序自带算法,triangle能够实现用更少的节点实现区域的划分，比如生成geo三模型时，可选用此选项. \n &
            & triangle的输入格式为: \n &
            & triangle,method=I;path=string,cmd=string,poly=string \n &
            & NOTE: \n &
            & 0)method=0,根据模型在当前目录下生成poly文件,然后根据cmd命令进行网格划分,并读入网络数据文件(node,ele,edge和neigh文件); \n &
            &   method=1,利用当前目录下已有poly文件,根据cmd命令进行网格划分,并读入网络数据; \n &
            &   method=2,根据模型在当前目录下生成poly文件,然后退出; \n &
            &   method=3,仅读入triangle生成的网格数据文件; \n &
            &   method=4,根据模型在当前目录下生成poly文件，然后根据cmd命令进行网格划分后退出; \n &
            &   method=-1(默认),令本命令无效(即不执行上述任一任务). \n &
            & 1)path=triangle.exe所在的目录.比如'c:\triangle.exe'  不输入时默认为当前目录或已设环境变量.\n & 
            & 2)cmd=triangle的命令参数.比如'cmd='-penI'' 不输入时默认为'-penI'. \n &
            & 3)polyfile=triangle的输入文件.比如'polyfile='.\a.poly'. 不输入时默认为当前目录下与mes文件同名的poly的文件。\n &
            & 4)当采用第三方程序Triangle进行网格划分且要进行地层插值时，所有节点都必须输入确定的地层信息。\N &
            & 5)triangle命令参数可参考:https://www.cs.cmu.edu/~quake/triangle.html. 常用的命令参数为\n &
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
    
    contains
    
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
        character(len=*),intent(in)::filepath
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
                    allocate(this.edge(this.nedge))
                    ismarker1=int(linedata(2))
                    n1=2+ismarker1
                    if(.not.allocated(this.adjlist)) allocate(this.adjlist(this.nnode))
                    do j=1,this.nedge
                        read(unit1,*) n2,ar1(1:n1)
                        this.edge(n2).v=int(ar1(1:2))
                        if(ismarker1>0) this.edge(n2).marker=int(ar1(n1))
                        call addadjlist(this.adjlist,int(ar1(1)),int(ar1(2)),n2)
                    end do 
                    this.ncedge=count(this.edge.marker/=0)
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
                        !if(this.edge(j).marker==-1) this.cedge(n1).cl=0                        
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
            case('polyfile')
	             this.polyfile=trim(adjustl(property(i).cvalue))                  
            case default
	            call Err_msg(property(i).name)
            end select
        end do 
        
        if(len_trim(this.polyfile)==0) this.polyfile=trim(adjustl(path_name))//'.poly'
        if(len_trim(this.path)==0) then
            this.path='.\triangle.exe'
            inquire(file=this.path,exist=isexist)
            if(.not.isexist) then
                this.path='triangle'
            endif
        endif
        !为保持与本程序的数据结构兼容,参数'en'为必输参数         
        if(index(trim(this.cmd),'e')==0) this.cmd=trim(this.cmd)//'e'
        if(index(trim(this.cmd),'n')==0) this.cmd=trim(this.cmd)//'n'
        if(index(trim(this.cmd),'I')==0) this.cmd=trim(this.cmd)//'I'     
        n1=index(trim(this.cmd),'u')
        if(n1==0.and.isnorefined/=1) this.cmd=trim(this.cmd)//'u'
        if(n1/=0.and.isnorefined==1) then
            this.cmd=this.cmd(1:n1-1)//this.cmd(n1+1:len_trim(this.cmd))
        endif
    endsubroutine
    
    subroutine execute_triangle_cmd(this,filepath)
    
        use ifport
        class(triangle_tydef)::this
        character(len=*),intent(in)::filepath
        integer::unit1,i
        INTEGER :: CSTAT, ESTAT,result
        CHARACTER(100) :: CMSG
        character(1024)::file1,cmd1
        character(16)::ftype1(4)
        
        if(this.method<0) return
        !write poly file
        if(this.method==0.or.this.method==2.or.this.method==4) then
            !outpoly
            call this.outpoly()
            if(this.method==2) stop
        endif
        !call triangle
        if(this.method==0.or.this.method==1.or.this.method==4) then
            call this.meshing() 
            if(this.method==4) stop
        endif
        
           
        call this.readmesh(filepath)

        
        
!10  format(i7,3(e24.16,','),e24.16)        
    end subroutine  
    
    subroutine meshbytriangle(this)
        class(triangle_tydef)::this
        INTEGER :: CSTAT, ESTAT
        CHARACTER(100) :: CMSG
        character(1024):: cmd1
        
        cmd1=trim(this.path)//' '//trim(this.cmd)//' '//trim(this.polyfile)
        !result = SYSTEMQQ ("tetgen -wefnNI pack_six_cube_poly.node")

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