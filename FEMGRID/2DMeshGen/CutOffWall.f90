module CutoffWall
USE meshDS,ONLY:strtoint,seg,segindex,edge,nedge,node,elt,nnode,nelt,ENLARGE_AR, &
    adjlist,soillayer,zone,Removeadjlist,addadjlist,ar2d_tydef
use ds_t,only:arr_t

implicit none
private
real ( kind = 8 ), parameter :: Pi = 3.141592653589793D+00
public::cowall,ncow,xzone,nxzone

type cutoffwall_type
    integer::NCP=0,mat=1     
    real(8)::thick=0 !墙厚
    integer,allocatable::cp(:) !控制点在输入数组中编号
    real(8),allocatable::be(:),te(:) !底高程,顶高程
    integer,allocatable::node(:),edge(:),element(:) !node=分离前防渗墙线上的母节点（注意不是防渗墙单元的节点）,element=生成的防渗墙防渗墙单元
    !edge=防渗墙的母线段    
    !integer,allocatable::wbedge(:)
    type(ar2d_tydef)::face(4)   !face(1),all faces stretched by wall lines;face(2)=faces on the real walls;FACE(3)=BOUNDARY EDGE OF FACE(1);FACE(4)=BOUNDARY EDGE OF FACE(2)
    character(512):: helpstring= &
    'CUTOFFWALL的输入格式为: \n &
        & 1)NCOW(防渗墙个数); \n & 
        & 2.1) NCP(边界控制点数),MAT(材料号),THICKNESS(墙厚) \n &
        & 2.2) CP(1:NCP)(各控制点号) \n &
        & 2.3) BE(1:NCP)(各控制点的墙底高程)，输入be(i)=-999表示i处墙底高程为模型地底高程. \n &
        & 2.4) TE(1:NCP)(各控制点的墙顶高程)，输入te(i)=-999表示i处墙顶高程为模型地表高程. \n &
        & 注意: \n &
        &     a) 如为闭合的防渗墙，则令首尾节点相同. \n &
        &     b) 各控制点定义的边必须已在控制线(CL)定义. 'C
        
        
    contains
        procedure,nopass::help=>write_help
        PROCEDURE::READIN=>COW_read
        PROCEDURE::GEN_ELEMENT=>GEN_COW_ELEMENT
        PROCEDURE::wall_elevation=>update_wall_elevation
        procedure::setboundary=>set_cow_boundary_node_edge
        !PROCEDURE::OUTPUT=>COW_write            
endtype
type(cutoffwall_type),allocatable::cowall(:)
integer::ncow=0

type xzone_tydef
    integer::izone,isx=1,ncutf=0
    real(8)::te
    integer,allocatable::cutface(:) !if isx==1,facets on the cut face,isx==2,cutedge around the cut face.
    character(1024)::helpstring= &
        &"xzone的作用是实现切割和开挖。如isx==1(默认),则将处于izone内部的高于te的土体挖掉，==0仅仅进行切割；==2仅对zone的边界线进行切割\n &
        & \n xZone的输入格式为:\n &
        & 1)nxzone //结构体数  \n &	
        & 2) izone,elevation,isx //分别为区域号及其底高程。\n &
        & "C        
contains
    procedure,nopass::help=>write_help
    procedure::readin=>xzone_readin
    procedure::set=>xzone_set_elevation
    !procedure::cut=>xzone_cut_gedge
endtype
type(xzone_tydef),allocatable::xzone(:)
integer::nxzone=0

    contains   
    
    subroutine xzone_readin(this,unit)
        class(xzone_tydef)::this
        integer,intent(in)::unit
        INTEGER::I,J,N1,DN=0,NINC1
        INTEGER,PARAMETER::DNMAX=1000
        REAL(8)::AR(DNMAX) 
        
 	    call strtoint(unit,ar,dnmax,dn,dnmax)
	    !n1=I
        THIS.izone=int(ar(1)) 
        !THIS.ndim=int(ar(2))
        THIS.te=ar(2)
        if(dn>2) THIS.isx=int(ar(3))
     
    endsubroutine    

    subroutine xzone_set_elevation(this)
        class(xzone_tydef)::this
        INTEGER::I,J,k,N1,ielt1
        integer::node1(nnode)
        
        node1=0
        do i=1,zone(this.izone).ntrie3n
            ielt1=zone(this.izone).trie3n(i)
            do j=1,3
                n1=elt(ielt1).node(j)
                if(node1(n1)==0) then
                    do k=soillayer,0,-1
                        if(node(n1).elevation(k)>this.te) then
                            node(n1).elevation(k)=this.te
                        else
                            exit !elevation已经排序
                        endif
                    enddo
                    node1(n1)=1
                endif
            enddo
        enddo
     
    endsubroutine
    
    !subroutine xzone_cut_gedge()
    !    class(xzone_tydef)::this
    !endsubroutine
    
subroutine write_help(helpstring,unit)
    USE IFQWIN
    implicit none
    character(*),intent(in)::helpstring
    !class(ZoneBC_type)::this
    integer,intent(in),optional::unit 
    integer(4)::oldcolor
    
    oldcolor = SETTEXTCOLOR(INT2(10))
	write(*,'(A)') trim(helpstring)
	oldcolor = SETTEXTCOLOR(INT2(15)) 


endsubroutine  

subroutine  COW_read(this,unit)
    implicit none
    class(cutoffwall_type)::this
    integer,intent(in)::unit
    INTEGER::I,J,N1,DN=0,NINC1
    INTEGER,PARAMETER::DNMAX=300
    REAL(8)::AR(DNMAX)

	call strtoint(unit,ar,dnmax,dn,dnmax)
	!n1=I
    THIS.NCP=int(ar(1)) !izone
    THIS.MAT=int(ar(2))
    THIS.THICK=ar(3)
    
    !ALLOCATE(THIS.CP(THIS.NCP),THIS.BE(THIS.NCP))
    call strtoint(unit,ar,dnmax,dn,dnmax)
    THIS.CP=INT(AR(1:THIS.NCP))
    call strtoint(unit,ar,dnmax,dn,dnmax)
    THIS.BE=AR(1:THIS.NCP)
    call strtoint(unit,ar,dnmax,dn,dnmax)
    THIS.TE=AR(1:THIS.NCP)    
	!end do

endsubroutine

subroutine set_cow_boundary_node_edge(this)
    implicit none
    class(cutoffwall_type)::this
    integer::i,iseg1,nnum1,nnum2,n1
    integer,allocatable::node1(:),node2(:),edge1(:),edge2(:)
    logical::isclose
    
    if(allocated(edge2)) deallocate(edge2)
    
    do i=1,this.ncp-1
        iseg1=segindex(this.cp(i),this.cp(i+1))
        !n1=n1+seg(iseg1).nnum-1
        !call seg(iseg1).getparas(nnum=nnum1,nedge=nedge1)

        node1=seg(iseg1).get_node([arr_t(this.cp(i)).x,arr_t(this.cp(i)).y])
        nnum1=size(node1)
        edge1=seg(iseg1).get_edge([arr_t(this.cp(i)).x,arr_t(this.cp(i)).y])
        nnum2=size(edge1)
        n1=nnum1-1
        if(.not.isclose.and.i==this.ncp-1) then !not close. add the last one
            n1=nnum1
        endif        
        node2=[node2,node1(1:n1)]
        
        edge2=[edge2,edge1]
        !if(size(node2)/=size(edge2)) then
        !    print *, size(node2),size(edge2)
        !endif
      
    enddo 
    
    this.node=node2
    this.edge=edge2
end subroutine

subroutine GEN_COW_ELEMENT(this)
    implicit none
    class(cutoffwall_type)::this
    integer::i,j,k,K1,n1,n2,n3,N4,iseg1,nnum1,nedge1,ielt1,ielt2,iflag,ED1,ED2,v1(2),ALLOC_ERR,nnum2
    integer,allocatable::node2(:),edge1(:),elt1(:),edge2(:),elt2(:),edge3(:),edge4(:)
    logical::isclose
    real(8)::val1(2)
    
    !gen new nodes
    n1=0;n2=0
    node.subbw=0 !借用
    node.layer=0
    isclose=this.cp(1)==this.cp(this.ncp)
    !if(allocated(edge2)) deallocate(edge2)
    if(.not.allocated(this.node)) then
    !do i=1,this.ncp-1
    !    iseg1=segindex(this.cp(i),this.cp(i+1))
    !    !n1=n1+seg(iseg1).nnum-1
    !    !call seg(iseg1).getparas(nnum=nnum1,nedge=nedge1)
    !
    !    node1=seg(iseg1).get_node([arr_t(this.cp(i)).x,arr_t(this.cp(i)).y])
    !    nnum1=size(node1)
    !    edge1=seg(iseg1).get_edge([arr_t(this.cp(i)).x,arr_t(this.cp(i)).y])
    !    nnum2=size(edge1)
    !    n1=nnum1-1
    !    if(.not.isclose.and.i==this.ncp-1) then !not close. add the last one
    !        n1=nnum1
    !    endif        
    !    node2=[node2,node1(1:n1)]
    !    
    !    edge2=[edge2,edge1]
    !    !if(size(node2)/=size(edge2)) then
    !    !    print *, size(node2),size(edge2)
    !    !endif
    !  
    !enddo
        call this.setboundary()
    endif
    
    node2=this.node
    edge2=this.edge      
    n1=size(node2)
    n2=size(edge2)

    
    if(size(node)<nnode+n1) call enlarge_ar(node,n1)
    do j=1,n1
        nnode=nnode+1        
        node(nnode)=node(node2(j))
        node(node2(j)).subbw=nnode
        node(nnode).subbw=node2(j)
    enddo
    
           
    !split elements
    !n1=size(edge2);n2=n1
    !if(.not.isclose) n2=n2-1
    !EDGE1=SPREAD(0,1,n1) !借用
    if(allocated(edge1)) deallocate(edge1)
    allocate(edge1(0:n1))
    edge1=0
    !n1=0
    allocate(edge3(0:n1))
    
    if(isclose) then
        edge3=abs([edge2(n1),edge2])
    else
        edge3=[0,abs(edge2),0]
    endif
    
    do j=1,n1

        elt1=FIND_ELT_ON_HEADING_LEFT(node2(j),edge3(j-1),edge3(j))
 
        elt2=[elt2,elt1] !

        EDGE1(J)=EDGE1(J-1)+SIZE(ELT1)
        
    enddo
    
    !UPDATE EDGE
    edge4=merge(edge3,0,.false.)    
    DO J=1,N1
        
        IF(J==1) THEN
            N4=0
        ELSE
            N4=EDGE1(J-1)
        ENDIF
        DO K=N4+1,EDGE1(J)
            IF(ELT2(K)<1) CYCLE
            DO K1=1,ELT(ELT2(K)).NNUM
                IF(ELT(ELT2(K)).NODE(K1)==NODE2(J)) THEN
                    ELT(ELT2(K)).NODE(K1)=NODE(NODE2(J)).SUBBW                    
                ENDIF
                ED1=ELT(ELT2(K)).EDGE(K1)
                IF(ANY(EDGE3==ED1)) THEN
                !GEN A NEW EDGE
                    nedge=nedge+1
                    if(nedge>size(edge,dim=1)) call enlarge_ar(edge,100)
                    edge(nedge)=edge(ED1)
                    EDGE(NEDGE).E(1)=ELT2(K)
                    EDGE(NEDGE).E(2)=-1
                    !EDGE(ED1).E=MERGE(EDGE(ED1).E,-1,EDGE(ED1).E/=ELT2(K))                     
                    ELT(ELT2(K)).EDGE(K1)=NEDGE
                    WHERE(EDGE3==ED1) EDGE4=NEDGE
                    ED1=NEDGE
                    
                ENDIF
                IF(ED1/=NEDGE.AND.ANY(EDGE(ED1).V==NODE2(J))) THEN
                    CALL Removeadjlist(adjList,EDGE(ED1).V(1),EDGE(ED1).V(2))
                ENDIF
                WHERE(EDGE(ED1).V==NODE2(J))
                    EDGE(ED1).V=NODE(NODE2(J)).SUBBW                    
                ENDWHERE
                IF(ANY(EDGE(ED1).V==NODE(NODE2(J)).SUBBW )) THEN
                    CALL addadjlist(adjList,EDGE(ED1).V(1),EDGE(ED1).V(2),ED1)
                ENDIF                
            ENDDO
        ENDDO
    
    ENDDO
    
        
    !gen zerothickness element
    n2=size(edge2)
    if(size(elt)<nelt+n2) call enlarge_ar(elt,n2)
    this.element=nelt+[1:n2]
    do i=1,n2
        nelt=nelt+1
        if(i==1) n4=nelt
        ed1=abs(edge2(i))
        if(edge2(i)>0) then
            v1=edge(ed1).v
        else
            v1=edge(ed1).v(2:1:-1)
        endif
        elt(nelt).et=-1 
        elt(nelt).nnum=4
        elt(nelt).node(1:4)=[v1,node(v1(2:1:-1)).subbw]
        elt(nelt).mat=this.mat
                
                
        ielt1=elt2(edge1(i))
        
        elt(nelt).adj(1)=ielt1
        elt(nelt).zn=0
        
        if(edge(ed1).e(1)/=ielt1) then
            ielt2=edge(ed1).e(1)
            edge(ed1).e(1)=nelt
        else
            ielt2=edge(ed1).e(2)
            edge(ed1).e(2)=nelt
        endif
        elt(nelt).adj(3)=ielt2
        
        
        
        !把单元L中原来与T相连的边更新为与EPT1相连
        !iflag==1 ept1=>null()
        IF(IELT2>0) call EDG(ielt1,ielt2,nelt,IFlag) 
        IF(IELT1>0) call edg(ielt2,ielt1,nelt,iflag)

        
        elt(nelt).adj(2)=nelt+1
        elt(nelt).adj(4)=nelt-1
       
        if(isclose) then
            if(i==1) elt(nelt).adj(4)=n4+n2-1
            if(i==n2) elt(nelt).adj(2)=n4
        else
            if(i==1) elt(nelt).adj(4)=-1
            if(i==n2) elt(nelt).adj(2)=-1                
        endif
       
        
        !update edges
        
        nedge=nedge+1
        if(nedge>size(edge,dim=1)) call enlarge_ar(edge,100)
        !edge(nedge)=edge(ed1)
        edge(nedge).v=elt(nelt).node([4,1])
        edge(nedge).isxyoverlap=1
        edge(nedge).e=[nelt,elt(nelt).adj(4)]
        if(i==1) n3=nedge
        call addadjlist(adjList,edge(nedge).v(1),edge(nedge).v(2),nedge)
        if((.not.isclose).and.(i==n2)) then
            nedge=nedge+1
            if(nedge>size(edge,dim=1)) call enlarge_ar(edge,100)
            !edge(nedge)=edge(ed1)
            edge(nedge).v=elt(nelt).node(2:3)
            edge(nedge).isxyoverlap=1
            edge(nedge).e=[nelt,elt(nelt).adj(2)]
            call addadjlist(adjList,edge(nedge).v(1),edge(nedge).v(2),nedge)            
        end if

        if(isclose.and.(i==n2)) then
            elt(nelt).edge=[ed1,n3,edge4(i),nedge]
        else
            elt(nelt).edge=[ed1,nedge+1,edge4(i),nedge]
        endif
    enddo
    
    
    deallocate(node2,edge1,elt1,edge2,elt2,edge3,edge4,STAT = ALLOC_ERR)
 
    node.subbw=0
    node.layer=0
endsubroutine

subroutine update_wall_elevation(this)
    implicit none
    class(cutoffwall_type)::this
    integer::i,j,k,K1,n1,iseg1,nnum1
    integer,allocatable::node1(:)
    logical::isclose
    
    
    isclose=this.cp(1)==this.cp(this.ncp)
    do i=1,this.ncp-1
        iseg1=segindex(this.cp(i),this.cp(i+1))
        !n1=n1+seg(iseg1).nnum-1
        !call seg(iseg1).getparas(nnum=nnum1,nedge=nedge1)

        node1=seg(iseg1).get_node([arr_t(this.cp(i)).x,arr_t(this.cp(i)).y])
        nnum1=size(node1)
        
        n1=nnum1-1
        if(.not.isclose.and.i==this.ncp-1) then !not close. add the last one
            n1=nnum1
        endif  
        !interpolate wall elevation
        if(abs(this.te(i)+999.)<1e-6) this.te(i)=node(node1(1)).elevation(soillayer)
        if(abs(this.te(i+1)+999.)<1e-6) this.te(i+1)=node(node1(nnum1)).elevation(soillayer)
        if(abs(this.be(i)+999.)<1e-6) this.be(i)=node(node1(1)).elevation(0)
        if(abs(this.be(i+1)+999.)<1e-6) this.be(i+1)=node(node1(nnum1)).elevation(0)
        
        do j=1,n1
            if(.not.allocated(node(node1(j)).we)) allocate(node(node1(j)).we(2))
            node(node1(j)).we(1)=linearint([arr_t(this.cp(i)).x,arr_t(this.cp(i)).y],this.te(i),[arr_t(this.cp(i+1)).x,arr_t(this.cp(i+1)).y],this.te(i+1),[node(node1(j)).x,node(node1(j)).y])
            node(node1(j)).we(2)=linearint([arr_t(this.cp(i)).x,arr_t(this.cp(i)).y],this.be(i),[arr_t(this.cp(i+1)).x,arr_t(this.cp(i+1)).y],this.be(i+1),[node(node1(j)).x,node(node1(j)).y])
        enddo        
    enddo
end subroutine


FUNCTION find_ELT_ON_HEADING_LEFT(V,BE1,BE2) RESULT(ELT1)
!返回以V为顶点，BE1-BE2前进方向左边的所有单元,BE1,BE2是以V公共顶点的边
!如果没有，返回ELT1(1)=-1
!如果BE1<1,这令BE1为与BE2夹角最大的边(最平行的边)
!如果BE2<1,这令BE2为与BE1夹角最最小的边(最平行的边)
IMPLICIT NONE
INTEGER,INTENT(IN)::V
INTEGER,INTENT(INOUT)::BE1,BE2
INTEGER,ALLOCATABLE::ELT1(:)

INTEGER::I,J,N1,n2,IELT1,ELT2(100),NC1,IELT2=0,IEDGE0,V0,e1(2)
REAL(8)::T1,x1(2,3),t2,t3=0

E1=[be1,be2]

if(all(e1<1)) then
    print *,'no edge supplied.'
    elt1=[-1]    
    return
endif
    

if(any(E1<1)) then
    if(E1(1)<1) then
        n2=2
    else
        n2=1
    endif
    n1=minloc(abs(adjlist(v).edge(1:adjlist(v).count)-E1(n2)),dim=1)
    x1(:,3)=[node(N1).x,node(N1).y]
    t1=1e6;
    do i=1,adjlist(v).count
        if(adjlist(v).edge(i)==E1(n2)) cycle
        n1=adjlist(v).node(i)
        x1(:,1)=[node(n1).x,node(n1).y]
        x1(:,2)=[node(v).x,node(v).y]
        t2=abs(angle_rad_2d (X1(:,1),X1(:,2),X1(:,3))-Pi)
        if(t2<t1) then
            t1=t2
            E1(mod(n2,2)+1)=adjlist(v).edge(i)
        endif        
    enddo
endif    
BE1=E1(1);BE2=E1(2)

!FOUND THE FIRST ELEMENT
IF(EDGE(E1(1)).V(1)==V) THEN
    V0=EDGE(E1(1)).V(2)
ELSEIF(EDGE(E1(1)).V(2)==V) THEN
    V0=EDGE(E1(1)).V(1)
ELSE
    PRINT *,'THE EDGE E1(1) SHOULD INCLUDE THE VETEX V. E1(1),V=',E1(1),V
    STOP
ENDIF

IELT2=0
DO I=1,2
    IELT1=EDGE(E1(1)).E(I)
    IF(IELT1>0.AND.IELT2==0) THEN
        DO J=1,ELT(IELT1).NNUM
            !T1=((NODE(ELT(IELT1).NODE(J)).X-NODE(V0).X)**2+(NODE(ELT(IELT1).NODE(J)).Y-NODE(V0).Y)**2)**0.5 !V0很可能已经分离出两个重合的节点。
        
            IF(ELT(IELT1).NODE(J)==V0.AND.ELT(IELT1).NODE(MOD(J,3)+1)==V) THEN
                IELT2=IELT1
                IEDGE0=J
                EXIT               
            ENDIF
        ENDDO
    ENDIF
ENDDO  

NC1=0
IF(IELT2>0) THEN
  
    DO WHILE(.TRUE.)
        
        NC1=NC1+1        
        IF(NC1>100) THEN
            ERROR STOP  'ERROR IN FUNCTION find_ELT_ON_HEADING_LEFT'
        ENDIF        
        ELT2(NC1)=IELT2 
        
        N1=MOD(IEDGE0,ELT(IELT2).NNUM)+1
        IEDGE0=ELT(IELT2).EDGE(N1)
        IF(IEDGE0==E1(2)) EXIT
        IELT2=ELT(IELT2).ADJ(N1)
        IF(IELT2<1) EXIT        
        DO I=1,ELT(IELT2).NNUM
            IF(ELT(IELT2).EDGE(I)==IEDGE0) THEN
                IEDGE0=I
                EXIT
            ENDIF
        ENDDO
        
    ENDDO
    ELT1=ELT2(1:NC1)
ELSE
    ELT1=[-1]
ENDIF



endfunction


function angle_rad_2d ( p1, p2, p3 )

!*****************************************************************************80
!
!! ANGLE_RAD_2D returns the angle in radians swept out between two rays in 2D.
!
!  Discussion:
!    Clockwise is positive
!    Except for the zero angle case, it should be true that
!
!      ANGLE_RAD_2D ( P1, P2, P3 ) + ANGLE_RAD_2D ( P3, P2, P1 ) = 2 * PI
!
!        P1
!        /
!       /    
!      /     
!     /  
!    P2--------->P3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(2), P2(2), P3(2), define the rays
!    P1 - P2 and P3 - P2 which define the angle.
!
!    Output, real ( kind = 8 ) ANGLE_RAD_2D, the angle swept out by the rays,
!    in radians.  0 <= ANGLE_RAD_2D < 2 * PI.  If either ray has zero
!    length, then ANGLE_RAD_2D is set to 0.
!
  implicit none

  real ( kind = 8 ) angle_rad_2d
  real ( kind = 8 ) p(2)
  real ( kind = 8 ) p1(2)
  real ( kind = 8 ) p2(2)
  real ( kind = 8 ) p3(2)
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00

  p(1) = ( p3(1) - p2(1) ) * ( p1(1) - p2(1) ) &
       + ( p3(2) - p2(2) ) * ( p1(2) - p2(2) )

  p(2) = ( p3(1) - p2(1) ) * ( p1(2) - p2(2) ) &
       - ( p3(2) - p2(2) ) * ( p1(1) - p2(1) )

  if ( all ( p(1:2) == 0.0D+00)  ) then
    angle_rad_2d = 0.0D+00
    return
  end if

  angle_rad_2d = atan2 ( p(2), p(1) )

  if ( angle_rad_2d < 0.0D+00 ) then
    angle_rad_2d = angle_rad_2d + 2.0D+00 * r8_pi
  end if

  return
end

real(8) function linearint(x1,v1,x2,v2,xi)
    real(8),intent(in)::x1(:),v1,x2(:),v2,xi(:)
    real(8)::dx(size(x1))
    real(8)::k
    integer::i,n1
    
    if(abs(v1-v2)<1e-7) then
        linearint=v1
        return
    endif
    
    n1=size(x1,dim=1)
    dx=x2-x1
    
    do i=1,n1
        if(abs(dx(i))>1e-6) then
            linearint=(v2-v1)/dx(i)*(xi(i)-x1(i))+v1            
            exit
        endif
    enddo
    
    if(i>n1) linearint=1/0.d0
    
    
end function


end module
    
    
    