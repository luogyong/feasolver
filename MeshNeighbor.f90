MODULE MESHNEIGHBOR

IMPLICIT NONE

PRIVATE

PUBLIC::find_neighbors_mesh,find_mesh_edge,tet_edge_tydef

type tet_edge_tydef
    integer::nnode=0,nsize=0  
    integer,allocatable::node(:) !1,ele1,2,edge of ele1;3,ele2,4, edge of ele2;...
contains
    procedure::push=>ar2d_node_fill
    procedure,private,nopass::enlarge_node=>I_ENLARGE_AR
endtype




CONTAINS


    subroutine find_neighbors_mesh(nnode,nelt,node,adj,geo,ndim)
        implicit none
        integer,intent(in)::nnode,nelt,ndim,node(nnode,nelt)
        integer,intent(out)::adj(nnode,nelt)
        integer,allocatable,intent(out)::geo(:,:) !faces or edges in the model.
        !geo(1:4,i)=[elt1,index1,elt2,index2] !elements and its face/edge id share the UNIFIED face/edge
        !node(1:10,nelt):node(9,:)=element node number;node(10,:)=eshape,!eshape:102=line,203=triangle,204=quadrilateral,304=tet,308=hex,306=prism

        IF(ndim==2) THEN
            CALL triangulation_neighbor_triangles (nnode,nelt,node,adj,geo)
        ELSE
            CALL tet_mesh_neighbor_tets (nnode,nelt,node,adj,geo)
        ENDIF

        WHERE(ADJ<0) ADJ=0

    endsubroutine

    subroutine find_mesh_edge( triangle_order, triangle_num, &
    triangle_node,AEDGE )
    
    !base the subroutine triangulation_neighbor_triangles ( triangle_order, triangle_num, &
    !triangle_node, triangle_neighbor,AEDGE )
    !triangle_node(1:8,:)=nodes;triangle_node(9,:)=nedge;triangle_node(10,:)=eshape;
    implicit none

    integer ( kind = 4 ) ,intent(in)::triangle_num
    integer ( kind = 4 ) ,intent(in)::triangle_order    
    integer ( kind = 4 ) ,intent(in)::triangle_node(triangle_order,triangle_num)
    !integer ( kind = 4 ) ,intent(out)::triangle_neighbor(3,triangle_num)
    type(tet_edge_tydef),allocatable,intent(out)::AEDGE(:)
    type(tet_edge_tydef),allocatable::AEDGE1(:)
    integer ( kind = 4 ),allocatable:: col(:,:)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) icol
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    integer ( kind = 4 ) side1
    integer ( kind = 4 ) side2    
    integer ( kind = 4 ) tri    
    integer ( kind = 4 ) tri1
    integer ( kind = 4 ) tri2
    INTEGER::NE1,err,maxedge1,v1(2,12),iedge1
    !
    !  Step 1.
    !  From the list of nodes for triangle T, of the form: (I,J,K)
    !  construct the three neighbor relations:
    !
    !    (I,J,3,T) or (J,I,3,T),
    !    (J,K,1,T) or (K,J,1,T),
    !    (K,I,2,T) or (I,K,2,T)
    !
    !  where we choose (I,J,1,T) if I < J, or else (J,I,1,T)
    !
    maxedge1=SUM(triangle_node(9,:))
    ALLOCATE(col(4,maxedge1),STAT=ERR)
    maxedge1=0
    do tri = 1, triangle_num
        
        SELECT CASE(triangle_node(10,TRI))

        CASE(102)
            v1(:,1)=triangle_node([1,2],tri)           
        CASE(203)            
            v1(:,1:3)=reshape(triangle_node([1,2,2,3,3,1],tri),([2,3]))
        case(204)
            v1(:,1:4)=reshape(triangle_node([1,2,2,3,3,4,4,1],tri),([2,4]))
        case(304)
            v1(:,1:6)=reshape(triangle_node([1,2,2,3,3,1,1,4,2,4,3,4],tri),([2,6])) 
        case(306)
            v1(:,1:9)=reshape(triangle_node([1,2,2,3,3,1,4,5,5,6,6,4,1,4,2,5,3,6],tri),([2,9])) 
        case(308)
            v1(:,1:12)=reshape(triangle_node([1,2,2,3,3,4,4,1,&
									        5,6,6,7,7,8,8,5,&
									        1,5,2,6,3,7,4,8],tri),([2,12]))                                    
        CASE DEFAULT
            error stop "no such eshape. error in=triangulation_neighbor_triangles."
        END SELECT

        do iedge1=1,triangle_node(9,TRI)
            maxedge1=maxedge1+1
            if ( v1(1,iedge1) < v1(2,iedge1) ) then
                i=1;j=2
            else
                i=2;j=1
            end if
            col(1:4,maxedge1) = (/ v1(i,iedge1), v1(j,iedge1), iedge1, tri /)
        enddo
    end do
    !
    !  Step 2. Perform an ascending dictionary sort on the neighbor relations.
    !  We only intend to sort on rows 1 and 2; the routine we call here
    !  sorts on rows 1 through 4 but that won't hurt us.
    !
    !  What we need is to find cases where two triangles share an edge.
    !  Say they share an edge defined by the nodes I and J.  Then there are
    !  two columns of COL that start out ( I, J, ?, ? ).  By sorting COL,
    !  we make sure that these two columns occur consecutively.  That will
    !  make it easy to notice that the triangles are neighbors.
    !
    call i4col_sort_a ( 4, maxedge1, col )
    !
    !  Step 3. Neighboring triangles show up as consecutive columns with
    !  identical first two entries.  Whenever you spot this happening,
    !  make the appropriate entries in TRIANGLE_NEIGHBOR.
    !

    !triangle_neighbor(1:3,1:triangle_num) = -1
    allocate(aedge1(maxedge1),stat=err)

    icol = 1
    NE1=0
    do
        
        if ( maxedge1 <= icol ) then
            IF(maxedge1 == icol)  THEN
                NE1=NE1+1
                call AEDGE1(NE1).push(col(4:3:-1,icol))
            ENDIF
            exit
        end if

        NE1=NE1+1

        if ( col(1,icol) /= col(1,icol+1) .or. col(2,icol) /= col(2,icol+1) ) then
            call AEDGE1(NE1).push(col(4:3:-1,icol))
            icol = icol + 1            
            cycle
        end if

        ! side1 = col(3,icol)
        ! tri1 = col(4,icol)
        
        ! side2 = col(3,icol+1)
        ! tri2 = col(4,icol+1)

        ! ! triangle_neighbor(side1,tri1) = tri2
        ! ! triangle_neighbor(side2,tri2) = tri1
        ! !AEDGE1(:,NE1)=[TRI1,SIDE1,TRI2,SIDE2]
        call AEDGE1(NE1).push(col(4:3:-1,icol))
        ! call AEDGE1(NE1).push(col(4:3:-1,icol+1))
        !icol = icol + 2
       
        !when this method is used to generate unified edge for 3d faces of a tet mesh , one edge may be shared by more than two faces.
        do j=icol+1,maxedge1
             if( all(col(1:2,icol)==col(1:2,j))) then
                call AEDGE1(NE1).push(col(4:3:-1,j))
            else
                exit
            endif 
        end do

        icol=j

    end do
    aedge=aedge1(1:ne1)

    deallocate(aedge1,stat=err)

    return
    end

    subroutine triangulation_neighbor_triangles ( triangle_order, triangle_num, &
    triangle_node, triangle_neighbor,AEDGE )
    !2023-5-11,extend to support element(line,qua) other than tri

    !*****************************************************************************80
    !
    !! TRIANGULATION_NEIGHBOR_TRIANGLES determines triangle neighbors.
    !
    !  Discussion:
    !
    !    A triangulation of a set of nodes can be completely described by
    !    the coordinates of the nodes, and the list of nodes that make up
    !    each triangle.  However, in some cases, it is necessary to know
    !    triangle adjacency information, that is, which triangle, if any,
    !    is adjacent to a given triangle on a particular side.
    !
    !    This routine creates a data structure recording this information.
    !
    !    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
    !    data items.
    !
    !    Note that ROW is a work array allocated dynamically inside this
    !    routine.  It is possible, for very large values of TRIANGLE_NUM,
    !    that the necessary amount of memory will not be accessible, and the
    !    routine will fail.  This is a limitation of the implementation of
    !    dynamic arrays in FORTRAN90.  One way to get around this would be
    !    to require the user to declare ROW in the calling routine
    !    as an allocatable array, get the necessary memory explicitly with
    !    an ALLOCATE statement, and then pass ROW into this routine.
    !
    !    Of course, the point of dynamic arrays was to make it easy to
    !    hide these sorts of temporary work arrays from the poor user!
    !
    !    This routine was revised to store the edge data in a column
    !    array rather than a row array.
    !
    !  Example:
    !
    !    The input information from TRIANGLE_NODE:
    !
    !    Triangle   Nodes
    !    --------   ---------------
    !     1         3      4      1
    !     2         3      1      2
    !     3         3      2      8
    !     4         2      1      5
    !     5         8      2     13
    !     6         8     13      9
    !     7         3      8      9
    !     8        13      2      5
    !     9         9     13      7
    !    10         7     13      5
    !    11         6      7      5
    !    12         9      7      6
    !    13        10      9      6
    !    14         6      5     12
    !    15        11      6     12
    !    16        10      6     11
    !
    !    The output information in TRIANGLE_NEIGHBOR:
    !
    !    Triangle  Neighboring Triangles
    !    --------  ---------------------
    !
    !     1        -1     -1      2
    !     2         1      4      3
    !     3         2      5      7
    !     4         2     -1      8
    !     5         3      8      6
    !     6         5      9      7
    !     7         3      6     -1
    !     8         5      4     10
    !     9         6     10     12
    !    10         9      8     11
    !    11        12     10     14
    !    12         9     11     13
    !    13        -1     12     16
    !    14        11     -1     15
    !    15        16     14     -1
    !    16        13     15     -1
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    14 February 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) TRIANGLE_ORDER, the order of the triangles.
    !
    !    Input, integer ( kind = 4 ) TRIANGLE_NUM, the number of triangles.
    !
    !    Input, integer ( kind = 4 ) TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM), 
    !    the nodes that make up each triangle.
    !
    !    Output, integer ( kind = 4 ) TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the three
    !    triangles that are direct neighbors of a given triangle.  
    !    TRIANGLE_NEIGHBOR(1,I) is the index of the triangle which touches side 1, 
    !    defined by nodes 2 and 3, and so on.  TRIANGLE_NEIGHBOR(1,I) is negative 
    !    if there is no neighbor on that side.  In this case, that side of the 
    !    triangle lies on the boundary of the triangulation.
    !
    implicit none

    integer ( kind = 4 ) ,intent(in)::triangle_num
    integer ( kind = 4 ) ,intent(in)::triangle_order    
    integer ( kind = 4 ) ,intent(in)::triangle_node(triangle_order,triangle_num)
    integer ( kind = 4 ) ,intent(out)::triangle_neighbor(triangle_order,triangle_num)
    integer,allocatable,intent(out)::AEDGE(:,:)

    integer ( kind = 4 ),allocatable:: col(:,:),AEDGE1(:,:)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) icol
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    integer ( kind = 4 ) side1
    integer ( kind = 4 ) side2    
    integer ( kind = 4 ) tri    
    integer ( kind = 4 ) tri1
    integer ( kind = 4 ) tri2
    INTEGER::NE1,maxedge1,ERR,m,v1(2,4),iedge1
    !
    !  Step 1.
    !  From the list of nodes for triangle T, of the form: (I,J,K)
    !  construct the three neighbor relations:
    !
    !    (I,J,3,T) or (J,I,3,T),
    !    (J,K,1,T) or (K,J,1,T),
    !    (K,I,2,T) or (I,K,2,T)
    !
    !  where we choose (I,J,1,T) if I < J, or else (J,I,1,T)
    !

    !count the face
    maxedge1=SUM(triangle_node(9,:))
    ALLOCATE(col(4,maxedge1),AEDGE1(4,maxedge1),STAT=ERR)
    maxedge1=0
    do tri = 1, triangle_num
        
        SELECT CASE(triangle_node(10,TRI))

        CASE(102)
            v1(:,1)=triangle_node([1,2],tri)           
        CASE(203)            
            v1(:,1:3)=reshape(triangle_node([1,2,2,3,3,1],tri),([2,3]))
        case(204)
            v1=reshape(triangle_node([1,2,2,3,3,4,4,1],tri),([2,4]))
        CASE DEFAULT
            error stop "no such eshape. error in=triangulation_neighbor_triangles."
        END SELECT

        do iedge1=1,triangle_node(9,TRI)
            maxedge1=maxedge1+1
            if ( v1(1,iedge1) < v1(2,iedge1) ) then
                i=1;j=2
            else
                i=2;j=1
            end if
            col(1:4,maxedge1) = (/ v1(i,iedge1), v1(j,iedge1), iedge1, tri /)
        enddo
    end do
    !
    !  Step 2. Perform an ascending dictionary sort on the neighbor relations.
    !  We only intend to sort on rows 1 and 2; the routine we call here
    !  sorts on rows 1 through 4 but that won't hurt us.
    !
    !  What we need is to find cases where two triangles share an edge.
    !  Say they share an edge defined by the nodes I and J.  Then there are
    !  two columns of COL that start out ( I, J, ?, ? ).  By sorting COL,
    !  we make sure that these two columns occur consecutively.  That will
    !  make it easy to notice that the triangles are neighbors.
    !
    call i4col_sort_a ( 4, maxedge1, col )
    !
    !  Step 3. Neighboring triangles show up as consecutive columns with
    !  identical first two entries.  Whenever you spot this happening,
    !  make the appropriate entries in TRIANGLE_NEIGHBOR.
    !
    triangle_neighbor(1:4,1:triangle_num) = -1

    icol = 1
    NE1=0
    do
        
        if ( maxedge1 <= icol ) then
            IF(maxedge1 == icol)  THEN
                NE1=NE1+1
                AEDGE1(:,NE1)=[col(4,icol),col(3,icol),0,0]
            ENDIF
            exit
        end if

        NE1=NE1+1

        if ( col(1,icol) /= col(1,icol+1) .or. col(2,icol) /= col(2,icol+1) ) then
            AEDGE1(:,NE1)=[col(4,icol),col(3,icol),0,0]
            icol = icol + 1            
            cycle
        end if

        side1 = col(3,icol)
        tri1 = col(4,icol)
        
        side2 = col(3,icol+1)
        tri2 = col(4,icol+1)

        triangle_neighbor(side1,tri1) = tri2
        triangle_neighbor(side2,tri2) = tri1
        AEDGE1(:,NE1)=[TRI1,SIDE1,TRI2,SIDE2]
        
        icol = icol + 2
    end do

    AEDGE=AEDGE1(:,1:NE1)

    return
    end


    subroutine tet_mesh_neighbor_tets ( tetra_order, tetra_num, tetra_node, &
    tetra_neighbor,Aface)
    !2023-05-12,modified by lgy£¬extended to support elements(tri,qua,prism,hex) other than tet.
    !*****************************************************************************80
    !
    !! TET_MESH_NEIGHBOR_TETS determines tetrahedron neighbors.
    !
    !  Discussion:
    !
    !    A tet mesh of a set of nodes can be completely described by
    !    the coordinates of the nodes, and the list of nodes that make up
    !    each tetrahedron.  In the most common case, four nodes are used.
    !    There is also a 10 node case, where nodes are also placed on
    !    the midsides of the tetrahedral edges.
    !
    !    This routine can handle 4 or 10-node tetrahedral meshes.  The
    !    10-node case is handled simply by ignoring the six midside nodes,
    !    which are presumed to be listed after the vertices.
    !
    !    The tetrahedron adjacency information records which tetrahedron
    !    is adjacent to a given tetrahedron on a particular face.
    !
    !    This routine creates a data structure recording this information.
    !
    !    The primary amount of work occurs in sorting a list of 4 * TETRA_NUM
    !    data items.
    !
    !    The neighbor tetrahedrons are indexed by the face they share with
    !    the tetrahedron.
    !
    !    Each face of the tetrahedron is indexed by the node which is NOT
    !    part of the face.  That is:
    
    !    * Neighbor 1 shares face 1 defined by nodes 2, 3, 4.  for compatable with elements other than triangle, the rule is modified by lgy,and no longer applied.
    !    * Neighbor 2 shares face 2 defined by nodes 1, 3, 4;
    !    * Neighbor 3 shares face 3 defined by nodes 1, 2, 4;
    !    * Neighbor 4 shares face 4 defined by nodes 1, 2, 3.
    
    !    For instance, if the (transposed) TETRA_NODE array was:
    !
    !    Row       1      2      3      4
    !    Col
    !
    !      1       4      3      5      1
    !      2       4      2      5      1
    !      3       4      7      3      5
    !      4       4      7      8      5
    !      5       4      6      2      5
    !      6       4      6      8      5
    !
    !    then the (transposed) TETRA_NEIGHBOR array should be:
    !
    !    Row       1      2      3      4
    !    Col
    !
    !      1      -1      2     -1      3
    !      2      -1      1     -1      5
    !      3      -1      1      4     -1
    !      4      -1      6      3     -1
    !      5      -1      2      6     -1
    !      6      -1      4      5     -1
    !
    !    09 February 2006: Jeff Borggaard reported that the code
    !    was failing when TETRA_NUM hit 10,000, but not at 9,999.
    !    This sounded like the effect of some odd internal compiler limit.
    !    He changed the internal FACES array from being implicitly
    !    allocated to being explicitly allocated, and the problem 
    !    went away.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    09 February 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) TETRA_ORDER, the order of the tetrahedrons.
    !
    !    Input, integer ( kind = 4 ) TETRA_NUM, the number of tetrahedrons.
    !
    !    Input, integer ( kind = 4 ) TETRA_NODE(TETRA_ORDER,TETRA_NUM), the nodes 
    !    that make up each tetrahedron.
    !
    !    Output, integer ( kind = 4 ) TETRA_NEIGHBOR(4,TETRA_NUM), the four 
    !    tetrahedrons that are direct neighbors of a given tetrahedron.  If there 
    !    is no neighbor sharing a given face, the index is set to -1.
    !
    implicit none

    integer ( kind = 4 ),intent(in)::tetra_num
    integer ( kind = 4 ),intent(in)::tetra_order
    integer ( kind = 4 ),intent(in):: tetra_node(tetra_order,tetra_num)
    integer ( kind = 4 ),intent(out):: tetra_neighbor(tetra_order,tetra_num)
    integer,allocatable,intent(out)::Aface(:,:)
    integer ( kind = 4 ) a
    integer ( kind = 4 ) b
    integer ( kind = 4 ) c,d
    integer ( kind = 4 ) face
    integer ( kind = 4 ) face1
    integer ( kind = 4 ) face2
    integer ( kind = 4 ), allocatable, dimension ( :, : ) :: faces,Aface1
    integer ( kind = 4 ) i
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    integer ( kind = 4 ) l
    integer ( kind = 4 ) tetra    
    integer ( kind = 4 ) tetra1
    integer ( kind = 4 ) tetra2
    integer::err,NF1,maxnf1,v1(4,6),if1
    integer::ihuge,ifun1

    maxnf1=sum(tetra_node(9,:))
    allocate ( faces(1:6,maxnf1), AFACE1(1:4,maxnf1),stat=err)
    !
    !  Step 1.
    !  From the list of nodes for tetrahedron T, of the form: (I,J,K,L)
    !  construct the four face relations:
    !
    !    (J,K,L,1,T)
    !    (I,K,L,2,T)
    !    (I,J,L,3,T)
    !    (I,J,K,4,T)
    !
    !  In order to make matching easier, we reorder each triple of nodes
    !  into ascending order.
    !
    maxnf1=0
    ihuge=10**9
    do tetra = 1, tetra_num

        !i = tetra_node(1,tetra)
        !j = tetra_node(2,tetra)
        !k = tetra_node(3,tetra)
        !l = tetra_node(4,tetra)
        !if(tetra_node(10,tetra)>304) then
        !    m=tetra_node(5,tetra)
        !    n=tetra_node(6,tetra)
        !endif
        !if(tetra_node(10,tetra)>306) then
        !    p=tetra_node(7,tetra)
        !    q=tetra_node(8,tetra)
        !endif        
      
        select case(tetra_node(10,tetra))

        case(304)
            v1(1:3,1:4)=RESHAPE(tetra_node([2,1,3,1,2,4,2,3,4,3,1,4],tetra),([3,4]))
            ifun1=3
        case(203)
            v1(1:3,1)=tetra_node([1,2,3],tetra)
            ifun1=3
        case(204)
            v1(1:4,1)=tetra_node([1,2,3,4],tetra)
            ifun1=4
        case(306)
            v1(1:4,1:5)=RESHAPE(tetra_node([2,1,3,3,4,5,6,6,1,2,5,4,2,3,6,5,3,1,4,6],tetra),([4,5]))
            v1(4,1:2)=ihuge
            ifun1=4
        case(308)
            v1=RESHAPE(tetra_node([1,2,3,4,5,6,7,8,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8],tetra),([4,6]))
            ifun1=4
        case(102)
            ifun1=2
          
        end select

        
        do if1=1,tetra_node(9,tetra)
            maxnf1=maxnf1+1
            if(ifun1==4) then
                call IFour_sort_a ( v1(1,if1), v1(2,if1), v1(3,if1),v1(4,if1), a, b, c,d )
                faces(1:6,maxnf1) = (/ a, b, c,d,if1, tetra /)
            elseif(ifun1==3) then
                call i4i4i4_sort_a ( v1(1,if1), v1(2,if1), v1(3,if1), a, b, c )
                faces(1:6,maxnf1) = (/ a, b, c,ihuge,if1, tetra /)
            else
                if(tetra_node(1,tetra)<tetra_node(2,tetra)) then
                    a=1;b=2
                else
                    a=2;b=1
                endif
                faces(1:6,maxnf1) = (/ tetra_node(a,tetra), tetra_node(b,tetra), ihuge,ihuge,1, tetra /)  
            endif
        enddo 


    end do
    !
    !  Step 2. Perform an ascending dictionary sort on the neighbor relations.
    !  We only intend to sort on rows 1:3; the routine we call here
    !  sorts on rows 1 through 5 but that won't hurt us.
    !
    !  What we need is to find cases where two tetrahedrons share a face.
    !  By sorting the columns of the FACES array, we will put shared faces
    !  next to each other.
    !
    call i4col_sort_a ( 6, maxnf1, faces )
    !
    !  Step 3. Neighboring tetrahedrons show up as consecutive columns with
    !  identical first three entries.  Whenever you spot this happening,
    !  make the appropriate entries in TETRA_NEIGHBOR.
    !
    tetra_neighbor(1:6,1:tetra_num) = -1

    face = 1
    NF1=0
    do

        if ( maxnf1 <= face ) then
            IF(maxnf1 == face) THEN
                NF1=NF1+1
                AFACE1(:,NF1)=[faces(6,face),faces(5,face),0,0]
            ENDIF
            exit
        end if
        NF1=NF1+1
        if ( all ( faces(1:4,face) == faces(1:4,face+1) ) ) then
            face1 = faces(5,face)
            tetra1 = faces(6,face)
            face2 = faces(5,face+1)
            tetra2 = faces(6,face+1)
            tetra_neighbor(face1,tetra1) = tetra2
            tetra_neighbor(face2,tetra2) = tetra1
            AFACE1(:,NF1)=[TETRA1,FACE1,TETRA2,FACE2]
            face = face + 2          
        else
            AFACE1(:,NF1)=[faces(6,face),faces(5,face),0,0]
            face = face + 1
        end if

    end do
    AFACE=AFACE1(:,1:NF1)

    deallocate ( faces,AFACE1,STAT=ERR )

    return
    end

    subroutine i4i4i4_sort_a ( i1, i2, i3, j1, j2, j3 )

    !*****************************************************************************80
    !
    !! I4I4I4_SORT_A ascending sorts a triple of I4's.
    !
    !  Discussion:
    !
    !    The program allows the reasonable call:
    !
    !      call i4i4i4_sort_a ( i1, i2, i3, i1, i2, i3 )
    !
    !    and this will return the reasonable result.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    13 October 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) I1, I2, I3, the values to sort.
    !
    !    Output, integer ( kind = 4 ) J1, J2, J3, the sorted values.
    !
    implicit none

    integer ( kind = 4 ) i1
    integer ( kind = 4 ) i2
    integer ( kind = 4 ) i3
    integer ( kind = 4 ) j1
    integer ( kind = 4 ) j2
    integer ( kind = 4 ) j3
    integer ( kind = 4 ) k1
    integer ( kind = 4 ) k2
    integer ( kind = 4 ) k3
    !
    !  Copy arguments, so that the user can make "reasonable" calls like:
    !
    !    call i4i4i4_sort_a ( i1, i2, i3, i1, i2, i3 )
    !
    k1 = i1
    k2 = i2
    k3 = i3

    j1 = min ( min ( k1, k2 ), min ( k2, k3 ) )
    j2 = min ( max ( k1, k2 ), &
        min ( max ( k2, k3 ), max ( k3, k1 ) ) )
    j3 = max ( max ( k1, k2 ), max ( k2, k3 ) )

    return
    end


   subroutine IFour_sort_a ( i1, i2, i3,i4, j1, j2, j3,j4 )
    !2023-5-12: sort four integer
    !*****************************************************************************80
    !
    !! I4I4I4_SORT_A ascending sorts a triple of I4's.
    !
    !  Discussion:
    !
    !    The program allows the reasonable call:
    !
    !      call i4i4i4_sort_a ( i1, i2, i3, i1, i2, i3 )
    !
    !    and this will return the reasonable result.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    13 October 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) I1, I2, I3, the values to sort.
    !
    !    Output, integer ( kind = 4 ) J1, J2, J3, the sorted values.
    !
    implicit none

    integer ( kind = 4 ) i1
    integer ( kind = 4 ) i2
    integer ( kind = 4 ) i3
    integer ( kind = 4 ) i4
    integer ( kind = 4 ) j1
    integer ( kind = 4 ) j2
    integer ( kind = 4 ) j3
    integer ( kind = 4 ) j4
    integer ( kind = 4 ) ka1(4)

    !
    !  Copy arguments, so that the user can make "reasonable" calls like:
    !
    !    call i4i4i4_sort_a ( i1, i2, i3, i1, i2, i3 )
    !
    ka1=[i1,i2,i3,i4]
    j1 = minval(ka1)
    j4 = maxval(ka1)
    j3=  maxval(ka1,ka1<j4)
    j2=  minval(ka1,ka1>j1)


    return
    end
    
  subroutine i4col_compare ( m, n, a, i, j, isgn )

    !*****************************************************************************80
    !
    !! I4COL_COMPARE compares columns I and J of an I4COL.
    !
    !  Example:
    !
    !    Input:
    !
    !      M = 3, N = 4, I = 2, J = 4
    !
    !      A = (
    !        1  2  3  4
    !        5  6  7  8
    !        9 10 11 12 )
    !
    !    Output:
    !
    !      ISGN = -1
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    12 June 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
    !
    !    Input, integer ( kind = 4 ) A(M,N), an array of N columns of vectors 
    !    of length M.
    !
    !    Input, integer ( kind = 4 ) I, J, the columns to be compared.
    !    I and J must be between 1 and N.
    !
    !    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
    !    -1, column I < column J,
    !     0, column I = column J,
    !    +1, column J < column I.
    !
    implicit none

    integer ( kind = 4 ) m
    integer ( kind = 4 ) n

    integer ( kind = 4 ) a(m,n)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) isgn
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    !
    !  Check.
    !
    if ( i < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
        write ( *, '(a,i8,a)' ) '  Column index I = ', i, ' is less than 1.'
        stop
    end if

    if ( n < i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
        write ( *, '(a,i8,a)' ) '  N = ', n, ' is less than column index I = ', i
        stop
    end if

    if ( j < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
        write ( *, '(a,i8,a)' ) '  Column index J = ', j, ' is less than 1.'
        stop
    end if

    if ( n < j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4COL_COMPARE - Fatal error!'
        write ( *, '(a,i8,a)' ) '  N = ', n, ' is less than column index J = ', j
        stop
    end if

    isgn = 0

    if ( i == j ) then
        return
    end if

    k = 1

    do while ( k <= m )

        if ( a(k,i) < a(k,j) ) then
        isgn = -1
        return
        else if ( a(k,j) < a(k,i) ) then
        isgn = +1
        return
        end if

        k = k + 1

    end do

    return
    end
    subroutine i4col_sort_a ( m, n, a )

    !*****************************************************************************80
    !
    !! I4COL_SORT_A ascending sorts an I4COL of columns.
    !
    !  Discussion:
    !
    !    In lexicographic order, the statement "X < Y", applied to two real
    !    vectors X and Y of length M, means that there is some index I, with
    !    1 <= I <= M, with the property that
    !
    !      X(J) = Y(J) for J < I,
    !    and
    !      X(I) < Y(I).
    !
    !    In other words, the first time they differ, X is smaller.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    25 September 2001
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, the number of rows of A, and the length of
    !    a vector of data.
    !
    !    Input, integer ( kind = 4 ) N, the number of columns of A.
    !
    !    Input/output, integer ( kind = 4 ) A(M,N).
    !    On input, the array of N columns of M-vectors.
    !    On output, the columns of A have been sorted in ascending
    !    lexicographic order.
    !
    implicit none

    integer ( kind = 4 ) m
    integer ( kind = 4 ) n

    integer ( kind = 4 ) a(m,n)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) indx
    integer ( kind = 4 ) isgn
    integer ( kind = 4 ) j

    if ( m <= 0 ) then
        return
    end if

    if ( n <= 1 ) then
        return
    end if
    !
    !  Initialize.
    !
    i = 0
    indx = 0
    isgn = 0
    j = 0
    !
    !  Call the external heap sorter.
    !
    do

        call sort_heap_external ( n, indx, i, j, isgn )
    !
    !  Interchange the I and J objects.
    !
        if ( 0 < indx ) then

            call i4col_swap ( m, n, a, i, j )
    !
    !  Compare the I and J objects.
    !
        else if ( indx < 0 ) then

            call i4col_compare ( m, n, a, i, j, isgn )

        else if ( indx == 0 ) then

            exit

        end if

    end do

    return
    end
    subroutine i4col_swap ( m, n, a, j1, j2 )

    !*****************************************************************************80
    !
    !! I4COL_SWAP swaps columns J1 and J2 of a integer array of column data.
    !
    !  Example:
    !
    !    Input:
    !
    !      M = 3, N = 4, J1 = 2, J2 = 4
    !
    !      A = (
    !        1  2  3  4
    !        5  6  7  8
    !        9 10 11 12 )
    !
    !    Output:
    !
    !      A = (
    !        1  4  3  2
    !        5  8  7  6
    !        9 12 11 10 )
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    04 April 2001
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
    !    in the array.
    !
    !    Input/output, integer ( kind = 4 ) A(M,N), an array of N columns 
    !    of length M.
    !
    !    Input, integer ( kind = 4 ) J1, J2, the columns to be swapped.
    !
    implicit none

    integer ( kind = 4 ) m
    integer ( kind = 4 ) n

    integer ( kind = 4 ) a(m,n)
    integer ( kind = 4 ) col(m)
    integer ( kind = 4 ) j1
    integer ( kind = 4 ) j2

    if ( j1 < 1 .or. n < j1 .or. j2 < 1 .or. n < j2 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4COL_SWAP - Fatal error!'
        write ( *, '(a)' ) '  J1 or J2 is out of bounds.'
        write ( *, '(a,i8)' ) '  J1 =    ', j1
        write ( *, '(a,i8)' ) '  J2 =    ', j2
        write ( *, '(a,i8)' ) '  N =     ', n
        stop

    end if

    if ( j1 == j2 ) then
        return
    end if

    col(1:m)  = a(1:m,j1)
    a(1:m,j1) = a(1:m,j2)
    a(1:m,j2) = col(1:m)

    return
    end

  subroutine sort_heap_external ( n, indx, i, j, isgn )

    !*****************************************************************************80
    !
    !! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
    !
    !  Discussion:
    !
    !    The actual list of data is not passed to the routine.  Hence this
    !    routine may be used to sort integers, reals, numbers, names,
    !    dates, shoe sizes, and so on.  After each call, the routine asks
    !    the user to compare or interchange two items, until a special
    !    return value signals that the sorting is completed.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    05 February 2004
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    A Nijenhuis and H Wilf,
    !    Combinatorial Algorithms,
    !    Academic Press, 1978, second edition,
    !    ISBN 0-12-519260-6.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of items to be sorted.
    !
    !    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
    !
    !    The user must set INDX to 0 before the first call.
    !    Thereafter, the user should not change the value of INDX until
    !    the sorting is done.
    !
    !    On return, if INDX is
    !
    !      greater than 0,
    !      * interchange items I and J;
    !      * call again.
    !
    !      less than 0,
    !      * compare items I and J;
    !      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
    !      * call again.
    !
    !      equal to 0, the sorting is done.
    !
    !    Output, integer ( kind = 4 ) I, J, the indices of two items.
    !    On return with INDX positive, elements I and J should be interchanged.
    !    On return with INDX negative, elements I and J should be compared, and
    !    the result reported in ISGN on the next call.
    !
    !    Input, integer ( kind = 4 ) ISGN, results of comparison of elements I 
    !    and J.
    !    (Used only when the previous call returned INDX less than 0).
    !    ISGN <= 0 means I is less than or equal to J;
    !    0 <= ISGN means I is greater than or equal to J.
    !
    implicit none

    integer ( kind = 4 ) i
    integer ( kind = 4 ), save :: i_save = 0
    integer ( kind = 4 ) indx
    integer ( kind = 4 ) isgn
    integer ( kind = 4 ) j
    integer ( kind = 4 ), save :: j_save = 0
    integer ( kind = 4 ), save :: k = 0
    integer ( kind = 4 ), save :: k1 = 0
    integer ( kind = 4 ) n
    integer ( kind = 4 ), save :: n1 = 0
    !
    !  INDX = 0: This is the first call.
    !
    if ( indx == 0 ) then

        i_save = 0
        j_save = 0
        k = n / 2
        k1 = k
        n1 = n
    !
    !  INDX < 0: The user is returning the results of a comparison.
    !
    else if ( indx < 0 ) then

        if ( indx == -2 ) then

        if ( isgn < 0 ) then
            i_save = i_save + 1
        end if

        j_save = k1
        k1 = i_save
        indx = -1
        i = i_save
        j = j_save
        return

        end if

        if ( 0 < isgn ) then
        indx = 2
        i = i_save
        j = j_save
        return
        end if

        if ( k <= 1 ) then

        if ( n1 == 1 ) then
            i_save = 0
            j_save = 0
            indx = 0
        else
            i_save = n1
            n1 = n1 - 1
            j_save = 1
            indx = 1
        end if

        i = i_save
        j = j_save
        return

        end if

        k = k - 1
        k1 = k
    !
    !  0 < INDX, the user was asked to make an interchange.
    !
    else if ( indx == 1 ) then

        k1 = k

    end if

    do

        i_save = 2 * k1

        if ( i_save == n1 ) then
        j_save = k1
        k1 = i_save
        indx = -1
        i = i_save
        j = j_save
        return
        else if ( i_save <= n1 ) then
        j_save = i_save + 1
        indx = -2
        i = i_save
        j = j_save
        return
        end if

        if ( k <= 1 ) then
        exit
        end if

        k = k - 1
        k1 = k

    end do

    if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
        i = i_save
        j = j_save
    else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
        i = i_save
        j = j_save
    end if

    return
    end

    SUBROUTINE I_ENLARGE_AR(AVAL,DSTEP)
        IMPLICIT NONE
        INTEGER,ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        INTEGER,ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0,ALLSTAT
        
        
        ALLOCATE(VAL1(dstep))
        AVAL=[AVAL,VAL1]
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
    END SUBROUTINE
    
    subroutine ar2d_node_fill(self,value)
        implicit none
        class(tet_edge_tydef)::self
        integer,intent(in)::value(2)
        integer::dsize1

        self.nnode=self.nnode+2       
        if(self.nnode>self.nsize) then
            if(self.nsize<1) then
                dsize1=10
            else
                dsize1=self.nsize
            endif
            CALL self.enlarge_node(self.node,dsize1)
            self.nsize=self.nsize+dsize1
        endif
        self.node(self.nnode-1:self.nnode)=value
    end subroutine
    

ENDMODULE