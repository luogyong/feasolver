module voropp
    use solverds,only:strtoint,incount,node,nnum,nvo,vo,element,outvar,ndimension,tdisp,&
        head,PoreSize,Phead,discharge,throatsize,throatRe,throatFriction,throatQ
    use tetgen_io,only:read_tetgen_file,element_tg,nelt_tg
    implicit none   
    public::voropp_handle
    
    private
    
    INTERFACE ENLARGE_AR
        MODULE PROCEDURE I_ENLARGE_AR,R_ENLARGE_AR,R_ENLARGE_AR2D
    END INTERFACE     
    
    !type::voropp_node_tydef
    !    integer::id
    !    real(8)::x(3)
    !endtype
    type voropp_face_tydef
        integer::ne,adjc(2)=0,id=0 !this face=globle vccpface(id) 
        !for voroppcell.face, when read adjc(1)=adjancent particle id,adjc(2)= the face indice of the voroppcell.face corresponding to the the surface
        !for vppface,adjc(1) right cell id; adjc(2) left cell id
        integer,allocatable::e(:),v(:)        
        real(8)::area=-999.,perimeters=-999.,normal(3)=0
    end type    
    type voropp_cell_tydef
        integer::id=0,nv=0,ne=0,nf=0        
        !logical::isP=.true.,isC=.true.
        real(8)::vol,area,center(3),x(3),r=0,perimeters=0.0
        real(8),allocatable::v(:,:)
        type(voropp_face_tydef),allocatable::face(:)
        integer,allocatable::vid(:)
        !character(256)::informat
    contains
        !procedure::readin=>voropp_read
    end type
    
    type(voropp_cell_tydef),allocatable::VoroppCell(:)
    integer,allocatable::index_p2c(:) !particle id to corresponding cell id
    integer::nvppc=0
    
    real(8),allocatable::vertex(:,:)    
    integer::nver=0
    integer,allocatable::ver2node(:)
    type(voropp_face_tydef),allocatable::vppface(:)
    integer::nvppf=0
    
    character(512)::tec_title,FILEPATH,VarString,voroformat,VarLocation
    character(len=64)::MeshPassVar,VoroPassVar
    integer::nvartec=0,nvoropassvar=0
    real(8),parameter::Pi=3.141592653589793
    real(8)::eps=1.0d-6
    
    contains
    
    
    subroutine voropp_handle(unit,volfile)
        integer,intent(in),optional::unit
        character(len=*),intent(in),optional::volfile
        logical::isfirstcall,isexist
        
        

        
        
        isfirstcall=.true.
        if(present(unit)) then
            inquire(unit,EXIST=isexist)
            if(isexist) call voropp_read(iunit=unit)
        else
            inquire(file=volfile,EXIST=isexist)
            if(isexist) call voropp_read(volfile=volfile)
        endif
        if(.not.isexist) then
            print *, 'no voro++ vol file found.'
            return
        endif
        
        if(nnum/=nver) then
            print *,'the numbers of the voro vertexes and the nodes are not identical. And '
            print *, 'no voro cells are out put.'
            return
        endif 
        
        call vppface_set()
        call vertex2node()
        !call gen_weightedTet()
        call tec_variable_string()
        call voropp_to_tecplot(isfirstcall)
        
    end subroutine
    

    subroutine vertex2node()
    !find the corresponding node for each voro vertex
        
        integer::i,j,k,n1,n2
        logical::isf1(nnum)
        real(8)::t1
        
        allocate(ver2node(nver))
        ver2node=0
        isf1=.false.
        do i=1,nver            
            do j=1,nnum
                if(isf1(j)) cycle
                if(node(j).coord(1)<vertex(1,i)-eps) cycle
                if(node(j).coord(1)>vertex(1,i)+eps) cycle
                if(node(j).coord(2)<vertex(2,i)-eps) cycle
                if(node(j).coord(2)>vertex(2,i)+eps) cycle
                if(node(j).coord(3)<vertex(3,i)-eps) cycle
                if(node(j).coord(3)>vertex(3,i)+eps) cycle                
                ver2node(i)=j
                isf1(j)=.true.
                !t1=norm2(node(j).coord-vertex(:,i))
                !if(abs(t1)<eps) then
                !    ver2node(i)=j                    
                !endif
            enddo
        enddo
        
        if(any(ver2node==0)) then            
            error stop 'failed to match between voro vertex and node.'
        endif
    end subroutine

    
    
    subroutine gen_weightedTet()
        use ifport
        
        integer::unit1,i
        INTEGER :: CSTAT, ESTAT,result
        CHARACTER(100) :: CMSG
        character(512)::file1,cmd1
        character(16)::ftype1(2)
        
        !out .node file
        unit1=20
        file1=trim(filepath)//'.node'
        open(unit1,file=file1,status='replace')
        
        write(unit1,*) nvppc,3,1,0
        do i=1,nvppc
            write(unit1,10) i,voroppcell(i).x,voroppcell(i).r**2
        end do
        close(unit1)
        
        
        cmd1="tetgen -wefnnNI "//trim(file1)
        !result = SYSTEMQQ ("tetgen -wefnNI pack_six_cube_poly.node")

        CALL EXECUTE_COMMAND_LINE (cmd1, EXITSTAT=ESTAT,CMDSTAT=CSTAT, CMDMSG=CMSG)
        IF (CSTAT > 0) THEN
            PRINT *, "Command execution failed with error ", TRIM(CMSG)
        ELSE IF (CSTAT < 0) THEN
            PRINT *, "Command execution not supported"
        ELSE
            PRINT *, "Command completed with status ", ESTAT
            open(unit1,file=file1,status='old')
            ftype1(1)='node'
            ftype1(2)='ele'
            call read_tetgen_file(unit1,ftype1)
        END IF
        
        
10  format(i7,3(e24.16,','),e24.16)        
    end subroutine
    
    subroutine voropp_to_tecplot(isfirstcall)
        logical,intent(inout)::isfirstcall
        integer::unit1,nc1,i,j,natr1,nc2,nv1
        !integer,external::incount 
        
        
        unit1=10

        if(isfirstcall) then
            open(unit1,file=trim(filepath)//'_voro++.tec',status='replace')
            write(unit1,10) trim(tec_title)//'_voro++'
            
            write(unit1,*) trim(VarString)
            isfirstcall=.false.
        else
            open(unit1,file=trim(filepath)//'_voro++.tec',status='old',access='append')
        endif
        
        
        
        if(nver>0) then
            write(unit1,30) nver,nvppc,nvppf,trim(adjustl(varlocation)),sum(vppface(:nvppf).ne)
            
            call write_data(unit1,0)
            nc1=20
            !node count per face
            write(unit1,70)
            write(unit1,50) vppface(:nvppf).ne
            !nodes per face
            write(unit1,71)
            do i=1,nvppf
                nc2=vppface(i).ne
                write(unit1,51) vppface(i).v
            enddo
            !left element per face
            write(unit1,72)
            write(unit1,50) vppface(:nvppf).adjc(2)
            !right element per face
            write(unit1,73)
            write(unit1,50) vppface(:nvppf).adjc(1)
            
            !no interpolate between cell. No duplicated points and faces are merged.
            call write_no_merged_zone() 
            
        endif
        
        
        !call write_mesh_zone() 
   
        
        close(unit1)
        
10  format('Title="',a<len_trim(tec_title)+7>,'"')
20  format('Variables="X","Y","Z",',<natr1>('"ATR',i<incount(i)>,'","Rn"')) 
21  format('Variables="X","Y","Z","Rn"')   
30  format('Zone,T=Voro++,zonetype=FEPOLYHEDRON,N=',I7,',E=',I7,',datapacking=BLOCK,Faces=',i7,',' &
    A,',TOTALNUMFACENODES=',i7, ',NUMCONNECTEDBOUNDARYFACES=0,TOTALNUMBOUNDARYCONNECTIONS=0') 
40  format(<nc1>(E24.16,1X),i7)
41  format(<nc1>(E24.16,','))
50  format(<nc1>(I7,1X))
51  format(<nc2>(I7,1X))
52  format('# faceid:',<nc2>(I7,1X))
60  format('#node count per face for cell=',i7)
61  format('#face nodes for cell=',i7)
62  format('#face left elements for cell=',i7)
63  format('#face right elements for cell=',i7)
70  format('#node count per face')
71  format('#nodes per face')
72  format('#left elements per face')
73  format('#right elements per face')
80  format(i<incount(nc1)>,'*0.')
    
    contains
    
    
    subroutine write_no_merged_zone()
        integer::nver1,nvppf1,nfv1
        integer::i,j,k,nc1,n1
        !integer,external::incount
        
        nver1=sum(voroppcell(1:nvppc).nv)
        nvppf1=sum(voroppcell(1:nvppc).nf)
        nfv1=0
        do i=1,nvppc
            nfv1=nfv1+sum(voroppcell(i).face.ne)
        enddo
        
        if(nver1>0) then
            write(unit1,30) nver1,nvppc,nvppf1,trim(adjustl(varlocation)),nfv1
            
            call write_data(unit1,2)
            nc1=20           
            
       
            !node count per face
            write(unit1,70)
            write(unit1,50) (voroppcell(k).face.ne,k=1,nvppc)
            !nodes per face
            write(unit1,71)
            n1=0
            do i=1,nvppc
                do j=1,voroppcell(i).nf
                    nc2=voroppcell(i).face(j).ne
                    write(unit1,51) voroppcell(i).face(j).v+n1
                enddo
                n1=voroppcell(i).nv+n1
            enddo
            !left element per face
            write(unit1,72)
            write(unit1,50) ((i,j=1,voroppcell(i).nf),i=1,nvppc)
            !write(unit1,82) ((voroppcell(i).nf,i),i=1,nvppc)
            !right element per face
            write(unit1,73)
            n1=nvppf1
            write(unit1,81) n1

            
        endif
30  format('Zone,T=Voro++_NoMerge,zonetype=FEPOLYHEDRON,N=',I7,',E=',I7,',datapacking=BLOCK,Faces=',i7,',' &
    A,',TOTALNUMFACENODES=',i7, ',NUMCONNECTEDBOUNDARYFACES=0,TOTALNUMBOUNDARYCONNECTIONS=0') 
40  format(<nc1>(E24.16,1X),i7)
41  format(<nc1>(E24.16,','))
50  format(<nc1>(I7,1X))
51  format(<nc2>(I7,1X))
52  format('# faceid:',<nc2>(I7,1X))
60  format('#node count per face for cell=',i7)
61  format('#face nodes for cell=',i7)
62  format('#face left elements for cell=',i7)
63  format('#face right elements for cell=',i7)
70  format('#node count per face')
71  format('#nodes per face')
72  format('#left elements per face')
73  format('#right elements per face')      
81  format(i<incount(n1)>,'*0')
82  format(20(i<incount(voroppcell(i).nf)>,'*',i<incount(i)>))    
    end subroutine
    
    subroutine write_mesh_zone()
        
        integer::i,j,k,nc1
        
        if(nvppc<1) return
        
        MeshPassVar=adjustl(MeshPassVar)
        nc1=len_trim(MeshPassVar)
        if(nc1>0) then
            MeshPassVar='PassiveVarList=['//trim(MeshPassVar(1:nc1-1))//']'
        endif
        
        write(unit1,30) nvppc,nelt_tg,trim(MeshPassVar)
        
        call write_data(unit1,1)   
        
        
        nc1=4
        do i=1,nelt_tg
            write(unit1,50) element_tg(i).node
        enddo

30  format('Zone,T=Particle,zonetype=FETETRAHEDRON,N=',I7,',E=',I7,',datapacking=BLOCK',2X,A<len_trim(MeshPassVar)>) 
40  format(<nc1>(E24.16,1X),i7)
41  format(<nc1>(E24.16,','))    
50  format(<nc1>(I7,1X))    
51  format(<nc1>(I7,','))         
60  FORMAT('PassiveVarList=[',I<INCOUNT(nc1)>,'-',I<INCOUNT(nc2)>,']') 
61  FORMAT('PassiveVarList=[',I<INCOUNT(nc1)>,']')         
!81  format(i<incount(n1)>,'*0')        
    end subroutine
   
    
    subroutine write_data(unit1,ztype)
        integer,intent(in)::unit1
        integer,intent(in)::ztype !=0,voro,=1,zont_particle,=2,no_merged
        integer::i,j,k,n1,nc1,idof
        
        
        nc1=20
        do i=1,len_trim(adjustl(voroformat))
            select case(voroformat(i:i))
            case('q')
                if(ztype==2) then 
                    do j=1,3
                        write(unit1,41) (voroppcell(k).v(j,:),k=1,nvppc)
                    enddo
                elseif(ztype==1) then
                    do j=1,3
                        write(unit1,41) (voroppcell(k).x(j),k=1,nvppc)
                    enddo
                elseif(ztype==0) then
                    do j=1,3
                        write(unit1,41) vertex(j,1:nver)
                    enddo
                endif
            case('x')
                if(ztype==2) then 
                    write(unit1,41) (voroppcell(k).v(1,:),k=1,nvppc)
                elseif(ztype==1) then
                    write(unit1,41) voroppcell.x(1)
                elseif(ztype==0) then
                    write(unit1,41) vertex(1,1:nver)
                endif
            case('y')
                if(ztype==2) then     
                    write(unit1,41) (voroppcell(k).v(2,:),k=1,nvppc)
                elseif(ztype==1) then
                    write(unit1,41) voroppcell.x(2)
                elseif(ztype==0) then
                    write(unit1,41) vertex(2,1:nver)
                endif    
            case('z')
                if(ztype==2) then 
                    write(unit1,41) (voroppcell(k).v(3,:),k=1,nvppc)  
                elseif(ztype==1) then
                    write(unit1,41) voroppcell.x(3)
                elseif(ztype==0) then
                    write(unit1,41) vertex(3,1:nver)
                endif    
            case('r')
                 
                write(unit1,41) voroppcell(1:nvppc).r
                write(unit1,41) 2.*voroppcell(1:nvppc).r
          
            case('v')
                if(ztype==0.or.ztype==2) then 
                    write(unit1,41) voroppcell(1:nvppc).vol
                elseif(ztype==1) then
                    write(unit1,41) 4.0/3.0*PI*voroppcell(1:nvppc).r**3
                endif    
            case('F')
                if(ztype==0.or.ztype==2) then 
                    write(unit1,41) voroppcell(1:nvppc).area
                elseif(ztype==1) then
                    write(unit1,41) 4.0*PI*voroppcell(1:nvppc).r**2
                endif
            case('E')
                if(ztype==0.or.ztype==2) then 
                    write(unit1,41) voroppcell(1:nvppc).perimeters
                elseif(ztype==1) then
                    n1=nvppc
                    write(unit1,82) n1
                endif
            end select
        enddo  
        
        if(ztype==2) then
            i=1
	        do while(i<=nvo)
		        select case(vo(i))
     
			        case(head)
				        idof=4                    
				        write(unit1,41)  (tdisp(node(ver2node(voroppcell(j).vid)).dof(idof)),j=1,nvppc)
			        case(PoreSize)
				        write(unit1,41)  (node(ver2node(voroppcell(j).vid)).Poresize,j=1,nvppc)                
			        case(Phead)
				        idof=4
				        write(unit1,41) ((tdisp(node(ver2node(voroppcell(j).vid)).dof(idof))-node(ver2node(voroppcell(j).vid)).coord(ndimension)),j=1,nvppc)
				    case(discharge)
				        write(unit1,41) (node(ver2node(voroppcell(j).vid)).q,j=1,nvppc)
                    !case(throatsize)
                    !    write(unit1,41)  element.property(2)
                    !case(throatRe)
                    !    write(unit1,41)  element.property(3)
                    !case(throatQ      )
                    !    write(unit1,41)  abs(element.flux(1))
                    !case(throatFriction)
                    !    write(unit1,41)  1.0/element.km(1,1)
            
                
		        end select		
		        i=i+outvar(vo(i)).nval
            end do 
        elseif(ztype==0) then
            i=1
	        do while(i<=nvo)
		        select case(vo(i))
			        case(head)
				        idof=4                    
				        write(unit1,41)  tdisp(node(ver2node).dof(idof))
			        case(PoreSize)
				 
				        write(unit1,41)  node(ver2node).Poresize                
			        case(Phead)
				        idof=4
				        write(unit1,41) tdisp(node(ver2node).dof(idof))-node(ver2node).coord(ndimension)

			        case(discharge)
				        write(unit1,41) node(ver2node).q
                    !case(throatsize)
                    !    write(unit1,41)  element.property(2)
                    !case(throatRe)
                    !    write(unit1,41)  element.property(3)
                    !case(throatQ)
                    !    write(unit1,41)  abs(element.flux(1))
                    !case(throatFriction)
                    !    write(unit1,41)  1.0/element.km(1,1)
            
                
		        end select		
		        i=i+outvar(vo(i)).nval
            end do             
        endif
        
        
41  format(<nc1>(E24.16,',')) 
82  format(i<incount(n1)>,'*0.')         
    end subroutine
 
    
end subroutine
    

    
    subroutine tec_variable_string()
        !character(len=*)::VarString,MeshPassVar,VoroPassVar
        !integer::nvoropassvar
        integer::i,nc1,nc2,nmeshvar
        character(len=64)::VarTec(100),str1
                
        VarString='Variables='
        VarLocation='VarLocation=(['
        nvartec=0
        MeshPassVar=''
        
        do i=1,len_trim(adjustl(voroformat))
            select case(voroformat(i:i))
            case('q','x','y','z')
                Varstring=trim(Varstring)//'"X","Y","Z",'
                nvartec=nvartec+3
            case('r')
                Varstring=trim(Varstring)//'"Radius","Diameter"'
                nvartec=nvartec+1
                write(str1,'(I<INCOUNT(nvartec)>)') nvartec                
                VarLocation=trim(VarLocation)//TRIM(ADJUSTL(STR1))//','
                nvartec=nvartec+1
                write(str1,'(I<INCOUNT(nvartec)>)') nvartec                
                VarLocation=trim(VarLocation)//TRIM(ADJUSTL(STR1))//','                
            case('v')
                Varstring=trim(Varstring)//'"Vol",'
                nvartec=nvartec+1
                write(str1,'(I<INCOUNT(nvartec)>)') nvartec                
                VarLocation=trim(VarLocation)//TRIM(ADJUSTL(STR1))//','                
            case('F')
                Varstring=trim(Varstring)//'"Area",'
                nvartec=nvartec+1
                write(str1,'(I<INCOUNT(nvartec)>)') nvartec                
                VarLocation=trim(VarLocation)//TRIM(ADJUSTL(STR1))//','                
            case('E')
                Varstring=trim(Varstring)//'"Perimeter",'
                nvartec=nvartec+1
                write(str1,'(I<INCOUNT(nvartec)>)') nvartec                
                VarLocation=trim(VarLocation)//TRIM(ADJUSTL(STR1))//','
                MeshPassVar=trim(MeshPassVar)//TRIM(ADJUSTL(STR1))//','
            end select
        enddo
        i=1
	    do while(i<=nvo)
		    select case(vo(i))

			    case(head,PoreSize,Phead,discharge)
                    Varstring=trim(Varstring)//'"'//trim(outvar(vo(i)).name)//'",'
                    nvartec=nvartec+1
			    !case(throatsize,throatRe,throatFriction,throatQ)
				   ! Varstring=trim(Varstring)//'"'//trim(outvar(vo(i)).name)//'",'
       !             nvartec=nvartec+1
       !             write(str1,'(I<INCOUNT(nvartec)>)') nvartec                
       !             VarLocation=trim(VarLocation)//TRIM(ADJUSTL(STR1))//','
		    end select		
		    i=i+outvar(vo(i)).nval
        end do        

        
        
        
        VarLocation=trim(VarLocation(:len_trim(VarLocation)-1))//']=CELLCENTERED)'
        
        !nvartec=4
        !VarTec(1:4)=['"X"','"Y"','"Z"','"MARKER"']
        !do i=1,node_tg(1).na
        !    WRITE(VarTec(nvartec+i),10) I
        !enddo
        !
        !nvartec=nvartec+node_tg(1).na
        !nmeshvar=nvartec
        !
        !do i=1,vnode_tg(1).na
        !    WRITE(VarTec(nvartec+i),11) I
        !enddo
        !
        !nvartec=nvartec+vnode_tg(1).na
        !!rn
        !write(VarTec(nvartec+1),*) '"VORO_Rn"'        
        !nvartec=nvartec+1
        !
        !VarString='Variables='
        !do i=1,nVarTec
        !    VarString=trim(adjustl(VarString))//trim(VarTec(i))//','
        !enddo
        !
        !nc1=nmeshvar+1;nc2=nvartec
        !if(nc2>nc1) then
        !    write(MeshPassVar,20) nc1,nc2
        !else
        !     write(MeshPassVar,21) nc1
        !endif
        !VoroPassVar=''
        !if(nmeshvar>3) then
        !    nc1=4;nc2=nmeshvar
        !    nvoropassvar=nc2-nc1+1
        !    if(nc2>nc1) then
        !        write(VoroPassVar,20) nc1,nc2
        !    else
        !        write(VoroPassVar,21) nc1
        !    endif
        !endif
        
10      FORMAT('"MESH_ATR',I<INCOUNT(I)>,'"')
11      FORMAT('"VORO_ATR',I<INCOUNT(I)>,'"')  
20      FORMAT('PassiveVarList=[',I<INCOUNT(nc1)>,'-',I<INCOUNT(nc2)>,']') 
21      FORMAT('PassiveVarList=[',I<INCOUNT(nc1)>,']') 

    end subroutine    
    
    
    subroutine vppface_set()
    
        integer::i,j,k,n1,n2
        real(8)::vec1(3),t1
        
        if(.not.allocated(vppface)) allocate(vppface(sum(voroppcell.nf)))
        nvppf=0
        do i=1,nvppc
            do j=1,voroppcell(i).nf
                if(voroppcell(i).face(j).id==0) then
                    nvppf=nvppf+1
                    voroppcell(i).face(j).id=nvppf
                    n1=voroppcell(i).face(j).adjc(1)                    
                    !check pairs
                    if(n1>0) then
                        n1=index_p2c(n1)
                        !voroppcell(i).face(j).adjc(1)=n1 !change it to celll id
                        
                        do k=1,voroppcell(n1).nf
                            if(voroppcell(n1).face(k).adjc(1)==voroppcell(i).id) then
                                voroppcell(n1).face(k).id=-nvppf 
                                !- ,the face normal is negetive to the face in vppface.
                                voroppcell(n1).face(k).adjc(2)=j
                                voroppcell(i).face(j).adjc(2)=k
                                !voroppcell(n1).face(k).adjc(1)=i !change it to cell id
                                exit
                            endif    
                        enddo
                        if(k>voroppcell(n1).nf) then
                            write(*,10) i,n1,voroppcell(i).face(j).area
                            error stop "error stop at sub=vppface_set"
                        endif
                    endif
                    vppface(nvppf)=voroppcell(i).face(j)
                    vppface(nvppf).v=voroppcell(i).vid(vppface(nvppf).v)
                    
                    !for vppface,adjc(1) right cell id; adjc(2) left cell
                    !vec1=voroppcell(i).center-voroppcell(i).v(:,vppface(nvppf).v(1))
                    !t1=dot_product(vec1,voroppcell(i).face(j).normal)                    
                    !if(t1<0.d0) then
                        vppface(nvppf).adjc(1)=max(n1,0) 
                        vppface(nvppf).adjc(2)=i
                    !else
                    !    vppface(nvppf).adjc(1)=i
                    !    vppface(nvppf).adjc(2)=n1                       
                    !                       
                    !endif
                endif
            enddo
        enddo
10  format('No face pair found between cell ',i7, ' and cell ',i7,'the face area is ',g)
    endsubroutine
    
    subroutine vertex_unique_insert(point,aindex)
    !insert point(:,:)  to vertex if they are not in vertex and 
    !return their indices by aindex
        real(8),intent(in)::point(3,*)
        integer,intent(out)::aindex(:)
        
        integer::i,j,n1,n2
        real(8)::t1
        
        n1=size(aindex,dim=1)
        if(.not.allocated(vertex)) allocate(vertex(3,10000))
        do i=1,n1
            n2=0
            do j=1,nver
                if(point(1,i)<vertex(1,j)-eps) cycle
                if(point(1,i)>vertex(1,j)+eps) cycle
                if(point(2,i)<vertex(2,j)-eps) cycle
                if(point(2,i)>vertex(2,j)+eps) cycle
                if(point(3,i)<vertex(3,j)-eps) cycle
                if(point(3,i)>vertex(3,j)+eps) cycle                
                t1=norm2(point(:,i)-vertex(:,j))
                if(abs(t1)<eps) then
                    aindex(i)=j
                    n2=1
                endif
            enddo
            if(n2==0) then
                nver=nver+1
                if(nver>size(vertex,dim=2)) call enlarge_ar(vertex,100)
                vertex(:,nver)=point(:,i)
                aindex(i)=nver
            endif
        enddo
    end subroutine
    
    subroutine voropp_read(iunit,volfile)
        !class(voropp_cell_tydef)::this
        use dflib
        USE IFPORT
        implicit none
        integer,intent(in),optional::iunit
        CHARACTER(len=*),intent(in),optional::volfile
        CHARACTER(3)        drive
	    CHARACTER(512)      dir
	    CHARACTER(512)      name,file1
	    CHARACTER(16)      ext,ext1(12)
	    integer(4)::length,msg
        logical::isexist,isin1
        integer::unit,i,hasread

               
        integer::na1=0,ismarker1=0,n1=0,n2,nelt1,nmax,maxset,j,k
        integer::nread,nset,nneed,nnode1,n3
        integer::iar1(10),ef
        
        parameter(nmax=1000)
	    parameter(maxset=100)
	
	    real(8)::linedata(nmax),ar1(nmax)
	    character(256)::set(maxset)
        
        if(present(iunit)) then
            inquire(iUNIT,name=file1,EXIST=isexist)
            unit=iunit
            if(.not.isexist) file1=volfile
        else
            unit=10
            file1=volfile
        endif
        
		length = SPLITPATHQQ(file1, drive, dir, name, ext)
		tec_title=trim(name)
        FILEPATH=trim(drive)//trim(dir)
        msg = CHDIR(FILEPATH)
        FILEPATH=trim(drive)//trim(dir)//trim(name)
        
        !if(trim(adjustl(ext))=='.cell') then
        !    FILEPATH=filepath(:len_trim(filepath)-2)
        !endif
        close(unit)        
        
        !call skipcomment(unit)
        !read(unit,*) nnode_tg,ndim_tg,na1,ismarker1
        
        open(unit,file=file1,status='old')
        nneed=nmax
        ef=0
        call strtoint(unit,linedata,nmax,nread,nneed,set,maxset,nset)
        nvppc=int(linedata(1))
        allocate(voroppcell(nvppc),Index_p2c(nvppc))
        Index_p2c=0;
        
        voroformat=''
        do i=1,nset
            n1=index(set(i),'%')
            if(n1>0) then
                voroformat=trim(voroformat)//set(i)(n1+1:n1+1)
            endif            
        enddo
        !input order check
        n1=index(trim(voroformat),'i')
        if(n1/=1) then
            error stop "outputformat: %i should be in the first place."
        endif
        n1=maxval([index(trim(voroformat),'q'),index(trim(voroformat),'x') ,index(trim(voroformat),'y') ,index(trim(voroformat),'z')])        
        n2=index(trim(voroformat),'w')        
        n3=index(trim(voroformat),'p')
        if(n3>0.and.(n1>n3.or.n1==0)) then
            write(*,40) 'p'
            error stop
        endif        
        if(n3>0.and.(n2>n3.or.n2==0)) then
            write(*,30) 'w','p'
            error stop
        endif
        n3=index(trim(voroformat),'P')
        if(n3>0.and.(n2>n3.or.n2==0)) then
            write(*,30) 'w','P'
            error stop
        endif
        n3=index(trim(voroformat),'c')
        if(n3>0.and.(n1>n3.or.n1==0)) then
            write(*,40) 'c'
            error stop
        endif          
        n2=index(trim(voroformat),'s')        
        n3=minval([index(trim(voroformat),'e'),index(trim(voroformat),'f'),index(trim(voroformat),'a'),index(trim(voroformat),'t'), &
           & index(trim(voroformat),'l'),index(trim(voroformat),'n')])
        if(n3>0.and.(n2>n3.or.n2==0)) then
            write(*,50) 's'
            error stop
        endif  
        n2=index(trim(voroformat),'a')        
        n3=index(trim(voroformat),'t')
        if(n3>0.and.(n2>n3.or.n2==0)) then
            write(*,30) 'a','t'
            error stop
        endif        
        
        n1=0;
        do i=1,nvppc     
            call strtoint(unit,linedata,nmax,nread,nneed,set,maxset,nset,ef)
            if(ef<0) exit
            n1=n1+1
            n2=1
            n3=int(linedata(n2))
            if(n3>size(index_p2c)) call enlarge_ar(index_p2c,100)
            index_p2c(n3)=n1
            isin1=.false.
            do j=1,len_trim(voroformat)
               
                select case(voroformat(j:j))
                !Particle-related:
                !  %i The particle ID number
                !  %x The x coordinate of the particle
                !  %y The y coordinate of the particle
                !  %z The z coordinate of the particle
                !  %q The position vector of the particle, short for "%x %y %z"
                !  %r The radius of the particle (only printed if -r enabled)                    
                case('i')
                    voroppcell(n1).id=int(linedata(n2))
                    n2=n2+1
                case('x')
                    voroppcell(n1).x(1)=linedata(n2)
                    n2=n2+1
                    isin1=.true.
                case('y')
                    voroppcell(n1).x(2)=linedata(n2)
                    n2=n2+1 
                    isin1=.true.
                case('z')
                    voroppcell(n1).x(3)=linedata(n2)
                    n2=n2+1 
                    isin1=.true.
                case('q')
                    voroppcell(n1).x=linedata(n2:n2+2)
                    n2=n2+3
                    isin1=.true.
                case('r')
                    voroppcell(n1).r=linedata(n2)
                    n2=n2+1 
                !Vertex-related:
                !%w The number of vertices in the Voronoi cell
                !%p A list of the vertices of the Voronoi cell in the format (x,y,z),
                !   relative to the particle center
                !%P A list of the vertices of the Voronoi cell in the format (x,y,z),
                !   relative to the global coordinate system
                !%o A list of the orders of each vertex
                !%m The maximum radius squared of a vertex position, relative to the
                !   particle center
                case('w')
                    voroppcell(n1).nv=int(linedata(n2))
                    allocate(voroppcell(n1).v(3,voroppcell(n1).nv))
                    n2=n2+1
                case('p')
                    !voroppcell(n1).isP=.false.
                 
                    voroppcell(n1).v=reshape(linedata(n2:n2+3*voroppcell(n1).nv-1),([3,voroppcell(n1).nv]))
                    n2=n2+3*voroppcell(n1).nv
                    do k=1,voroppcell(n1).nv
                        voroppcell(n1).v(:,k)=voroppcell(n1).v(:,k)+voroppcell(n1).x
                    enddo                    
                    allocate(voroppcell(n1).vid(voroppcell(n1).nv))
                    call vertex_unique_insert(voroppcell(n1).v,voroppcell(n1).vid)
                case('P')
                      
                    voroppcell(n1).v=reshape(linedata(n2:n2+3*voroppcell(n1).nv-1),([3,voroppcell(n1).nv]))
                    n2=n2+3*voroppcell(n1).nv
                    allocate(voroppcell(n1).vid(voroppcell(n1).nv))
                    call vertex_unique_insert(voroppcell(n1).v,voroppcell(n1).vid)                    
                case('o')
                    !not used yet.skipped 
                    n2=n2+voroppcell(n1).nv
                case('m')
                    !not used yet. skipped
                    n2=n2+1
                !Edge-related:
                !  %g The number of edges of the Voronoi cell
                !  %E The total edge distance
                !  %e A list of perimeters of each face
                case('g')
                    voroppcell(n1).ne=int(linedata(n2))
                    n2=n2+1
                case('E')
                    voroppcell(n1).perimeters=linedata(n2)
                    n2=n2+1
                case('e')
                    voroppcell(n1).face.perimeters=linedata(n2:n2+voroppcell(n1).nf-1)
                    n2=n2+voroppcell(n1).nf
                !Face-related:
                !  %s The number of faces of the Voronoi cell
                !  %F The total surface area of the Voronoi cell
                !  %A A frequency table of the number of edges for each face
                !  %a A list of the number of edges for each face
                !  %f A list of areas of each face
                !  %t A list of bracketed sequences of vertices that make up each face
                !  %l A list of normal vectors for each face
                !  %n A list of neighboring particle or wall IDs corresponding to each face                    
                case('s')
                    voroppcell(n1).nf=int(linedata(n2))
                    allocate(voroppcell(n1).face(voroppcell(n1).nf))
                    n2=n2+1
                case('F')
                    voroppcell(n1).area=linedata(n2)
                    n2=n2+1
                case('f')
                    voroppcell(n1).face.area=linedata(n2:n2+voroppcell(n1).nf-1)
                    n2=n2+voroppcell(n1).nf
                case('a')
                    voroppcell(n1).face.ne=int(linedata(n2:n2+voroppcell(n1).nf-1))
                    n2=n2+voroppcell(n1).nf
                    
                case('t')
                    do k=1,voroppcell(n1).nf
                     
                        voroppcell(n1).face(k).v=int(linedata(n2:n2+voroppcell(n1).face(k).ne-1))+1 !base 0 
                        n2=n2+voroppcell(n1).face(k).ne
                    enddo
                case('l')
                    do k=1,voroppcell(n1).nf
                        !allocate(voroppcell(n1).face(k).v(voroppcell(n1).face(k).ne))
                        voroppcell(n1).face(k).normal=linedata(n2:n2+2)
                        n2=n2+3
                    enddo                    
                case('n')
                    do k=1,voroppcell(n1).nf
                        !allocate(voroppcell(n1).face(k).v(voroppcell(n1).face(k).ne))
                        voroppcell(n1).face(k).adjc(1)=int(linedata(n2))
                        n2=n2+1
                    enddo
                !
                !Volume-related:
                !  %v The volume of the Voronoi cell
                !  %c The centroid of the Voronoi cell, relative to the particle center
                !  %C The centroid of the Voronoi cell, in the global coordinate system   
                case('v')
                    voroppcell(n1).vol=linedata(n2)
                    n2=n2+1
                case('c')
                    voroppcell(n1).center=linedata(n2:n2+2)+voroppcell(n1).x
                    n2=n2+3
                case('C')
                    !voroppcell(n1).isC=.true.
                    voroppcell(n1).center=linedata(n2:n2+2)
                    n2=n2+3                    
                case default                
                    write(*, 10), voroformat(j:j)
                    stop
                end select
            enddo
            if(nread/=n2-1) then
                write(*,20) n2-1,nread
                stop
            endif
                
        enddo
        nvppc=n1
        
        close(unit)    
10      format('Custom output format: %',A1,' is not supported yet.')
20      format("The number of readin data is expected to be ",i4,".but it is ",i4,'.sub=voropp_read')
30      format('Custom output format: %',A1,' should be before: %',A1,'.')
40      format('Custom output format: %q/(%x,%y,%z) should be before: %',A1,'.')
50      format('Custom output format: %e/%f/%a/%t/%l/%n should be after: %',A1,'.')        
    end subroutine
    
     
    
    SUBROUTINE I_ENLARGE_AR(AVAL,DSTEP)
        INTEGER,ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        INTEGER,ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
    
        LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1,SOURCE=AVAL)
        DEALLOCATE(AVAL)
        ALLOCATE(AVAL(LB1:UB1+DSTEP))
        AVAL(LB1:UB1)=VAL1
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
    END SUBROUTINE

    SUBROUTINE R_ENLARGE_AR(AVAL,DSTEP)
        REAL(8),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        REAL(8),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
    
        LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1,SOURCE=AVAL)
        DEALLOCATE(AVAL)
        ALLOCATE(AVAL(LB1:UB1+DSTEP))
        AVAL(LB1:UB1)=VAL1
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
    END SUBROUTINE
    
    SUBROUTINE R_ENLARGE_AR2D(AVAL,DSTEP)
        REAL(8),ALLOCATABLE,INTENT(INOUT)::AVAL(:,:)
        INTEGER,INTENT(IN)::DSTEP
        REAL(8),ALLOCATABLE::VAL1(:,:)
        INTEGER::LB1=0,UB1=0,I,LB2,UB2
    
        LB1=LBOUND(AVAL,DIM=2);UB1=UBOUND(AVAL,DIM=2)
        ALLOCATE(VAL1,SOURCE=AVAL)
        DEALLOCATE(AVAL)
        LB2=LBOUND(AVAL,DIM=1);UB2=UBOUND(AVAL,DIM=1)
        ALLOCATE(AVAL(LB2:UB2,LB1:UB1+DSTEP))
        DO I=LB2,UB2
            AVAL(I,LB1:UB1)=VAL1(I,LB1:UB1)
        ENDDO
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
    END SUBROUTINE    
    
end module 