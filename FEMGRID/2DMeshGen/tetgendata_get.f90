module tetgendata
    use MeshDs,only:lowcase,strtoint,skipcomment
    implicit none
    
    public::tetgendata_tydef
    
    private
    
    
    integer::order=1
    character(512)::filepath
 
    logical::isfirstcall=.true.

    real(8)::eps=1.d-6
    
    
    type::tetgen_node_tydef
        integer::na=0,marker=0,surfin=0,a1=0,uid=0      
        !for vnode_tg, marker>0 inside the con and is the nodeid for tecplot out
        !if(node(i).uid<i) then the node is a duplicated node, which is the copy of node(uid) 
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
        integer::nnode=0,marker=0,nedge=0,isclose=0
        !isclose,for vedge,vface and vcell =0 open/infinite. >0 close/finite
        !for vedge,marker=0 outside container;marker/=0,inside the con;=2,clipped;=3,clipped and coplan of the con; =4,clipped and not coplan;
        !for vface,isclose>0, close and (all nodes) inside the con,and it will be output to the tec file; =-1,close and some node outside the cont; =-2 close and all node outside the con;
        !for vface,iclose=0, open
        !for vcell, isclose>0,close and (all nodes) inside the con,and it will output to the the tec file.
        !for vface, marker==0,not clipped face; =3,clipped on only the generated edge is on one boundary; =4,clopped and the generated egde is on different boundary 
        !
        integer,allocatable::node(:)
        integer::cell(2)=-1 !for voronoi face, they are cell sharing the face, for vcell, cell(1) is the corresponding node id (node(cell(1))) for the cell
        real(8),allocatable::V(:)  
        !for edge ray, its direction, for face, its unit normal
        !for cell, its center.
        integer,allocatable::edge(:) !for voronoi face only
        !integer::celoc=0,surfin=-1 !for vface, celoc=clipped edge id; surfin is the boudary id where the vertex  of the clipped edge laying. 
        real(8)::property  !for faces,the area; for edges,the length.
    end type  
    
 
    type tetgendata_tydef
        integer::nnode,nelt,nface,nedge,nneigh
        integer::nvnode,nvcell,nvface,nvedge
        type(tetgen_node_tydef),allocatable::node(:),vnode(:)
        type(tetgen_element_tydef),allocatable::edge(:),face(:),elt(:), vface(:),&
                vedge(:),vcell(:),neigh(:)
        integer,allocatable::t2e(:,:),t2f(:,:),f2e(:,:)
        type(adjlist_tydef),allocatable::nadjlist(:),eadjlist(:),vnadjlist(:) !nodal,and edge list
        character(512)::file
    contains
        procedure::readin=>read_tetgen_file
        procedure::setadjlist=>tetgen_setup_adjlist
    endtype
    

    
    INTERFACE ENLARGE_AR
        MODULE PROCEDURE I_ENLARGE_AR,NODE_ENLARGE_AR,ELEMENT_ENLARGE_AR                        
    END INTERFACE 
    
    contains
    
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

    
    subroutine read_tetgen_file(self,unit,fext)
        use dflib
        USE IFPORT
        implicit none
        class(tetgendata_tydef)::self
        integer,intent(in),optional::unit
        character(len=*),optional,intent(in)::fext(:)
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
                    call find_duplicated_node(self.vnode)
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
        
        parameter(nmax=100)
	    parameter(maxset=100)
	
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
    
    
    subroutine find_duplicated_node(node1)
    !remove duplicated node and stored in uvertex
        type(tetgen_node_tydef)::node1(:)
        integer::i,j,k,n1,n2
        
        !allocate(ver2node(nver))
        n1=size(node1)
        do i=1,n1
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