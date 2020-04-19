module CutoffWall
USE meshDS,ONLY:strtoint,seg,segindex,edge,nedge,node,elt,nnode,nelt,ENLARGE_AR
use ds_t,only:arr_t

implicit none
private
    
type cutoffwall_type
    integer::NCP=0,mat=1 !icl为防渗墙所在的控制线编号
    
    real(8)::thick=0 !墙厚
    integer,allocatable::cp(:) !控制点在输入数组中编号
    real(8),allocatable::be(:) !底高程
        
        
    character(512):: helpstring= &
    'CUTOFFWALL的输入格式为: \n &
        & 1)NCOW(防渗墙个数); \n & 
        & 2.1) NCP(边界控制点数),MAT(材料号),THICKNESS(墙厚) \n &
        & 2.2) CP(1:NCP)(各控制点号) \n &
        & 2.3) BE(1:NCP)(各控制点的墙底高程) \n'C    
    contains
        procedure,nopass::help=>write_help
        PROCEDURE::READIN=>COW_read
        PROCEDURE::GEN_ELEMENT=>GEN_COW_ELEMENT
        !PROCEDURE::OUTPUT=>COW_write            
endtype
    
contains    
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
    
    ALLOCATE(THIS.CP(THIS.NCP),THIS.BE(THIS.NCP))
    call strtoint(unit,ar,dnmax,dn,dnmax)
    THIS.CP=INT(AR(1:SIZE(THIS.CP,DIM=1)))
    call strtoint(unit,ar,dnmax,dn,dnmax)
    THIS.BE=AR(1:SIZE(THIS.BE,DIM=1))
    
	!end do

endsubroutine


subroutine GEN_COW_ELEMENT(this)
    implicit none
    class(cutoffwall_type)::this
    integer::i,j,k,K1,n1,n2,n3,N4,iseg1,nnum1,nedge1,ielt1,ielt2,iflag,ED1,ED2
    integer,allocatable::node1(:),node2(:),edge1(:),elt1(:),edge2(:),elt2(:)
    logical::isclose
    
    !gen new nodes
    n1=0
    node.subbw=0 !借用
    node.layer=0
    isclose=this.cp(1)==this.cp(this.ncp)
    do i=1,this.ncp-1
        iseg1=segindex(this.cp(i),this.cp(i+1))
        !n1=n1+seg(iseg1).nnum-1
        !call seg(iseg1).getparas(nnum=nnum1,nedge=nedge1)

        node1=seg(iseg1).get_node([arr_t(this.cp(i)).x,arr_t(this.cp(i)).y])
        if(i<this.ncp-1) then
            n1=nnum1-1
        elseif(.not.isclose) then !not close. add the last one
            n1=nnum1
        endif        
        node2=[node2,node1(1:n1)]
        edge1=seg(iseg1).get_edge([arr_t(this.cp(i)).x,arr_t(this.cp(i)).y])
        edge2=[edge2,edge1]
    enddo
      
    n1=size(node2)
    if(size(node)<nnode+n1) call enlarge_ar(node,n1)
    do j=1,n1
        nnode=nnode+1
        node(nnode)=node(node2(j))
        node(node1(j)).subbw=nnode
        node(nnode).subbw=node2(j)
    enddo        
        
    !split elements
    n1=size(edge2);n2=n1
    if(.not.isclose) n2=n2-1
    EDGE1=SPREAD(0,1,N2) !借用   
    !n1=0
    do j=1,n2
        n3=mod(j,n2)+1
        elt1=FIND_ELT_ON_HEADING_LEFT(node2(n3),edge2(j),edge2(n3))
        elt2=[elt2,elt1] !
        IF(J==1) THEN
            edge1(J)=size(elt1)
        ELSE
            EDGE1(J)=EDGE1(J-1)+SIZE(ELT1)
        ENDIF
    enddo
    
    !UPDATE EDGE
    DO J=1,N2
        n3=mod(j,n2)+1
        IF(J==1) THEN
            N4=0
        ELSE
            N4=EDGE1(J-1)
        ENDIF
        DO K=N4+1,EDGE1(J)
            IF(ELT2(K)<1) CYCLE
            DO K1=1,ELT(ELT2(K)).NNUM
                IF(ELT(ELT2(K)).NODE(K1)==NODE2(N3)) THEN
                    ELT(ELT2(K)).NODE(K1)=NODE(NODE2(N3)).SUBBW
                ENDIF
                ED1=ELT(ELT2(K)).EDGE(K1)
                IF(ANY(EDGE2==ED1)) THEN
                !GEN A NEW EDGE
                    nedge=nedge+1
                    if(nedge>size(edge,dim=1)) call enlarge_ar(edge,100)
                    edge(nedge)=edge(ED1)
                    EDGE(NEDGE).E(1)=ELT2(K)
                    EDGE(NEDGE).E(2)=-1
                    EDGE(ED1).E=MERGE(EDGE(ED1).E,-1,EDGE(ED1).E/=ELT2(K))                    
                    ELT(ELT2(K)).EDGE(K1)=NEDGE
                    ED1=NEDGE
                ENDIF
                
                WHERE(EDGE(ED1).V==NODE2(N3))
                    EDGE(ED1).V=NODE(NODE2(N3)).SUBBW
                ENDWHERE
                
                
            ENDDO
        ENDDO
    
    ENDDO
    
        
    !gen zerothickness element
    if(size(elt)<nelt+n2) call enlarge_ar(elt,n2)
    do i=1,n2
        nelt=nelt+1

        elt(nelt).et=-1 
        elt(nelt).nnum=4
        elt(nelt).node(1:4)=[edge(edge2(i)).v,node(edge(edge2(i)).v(2:1:-1)).subbw]
        ielt1=elt2(i)
        elt(nelt).adj(3)=ielt1
            
        if(edge(edge2(i)).e(1)/=ielt1) then
            ielt2=edge(edge2(i)).e(1)
        else
            ielt2=edge(edge2(i)).e(2)
        endif
        elt(nelt).adj(1)=ielt2
        
        
        
        !把单元L中原来与T相连的边更新为与EPT1相连
        !iflag==1 ept1=>null()
        IF(IELT2>0) call EDG(ielt1,ielt2,nelt,IFlag) 
        IF(IELT1>0) call edg(ielt2,ielt1,nelt,iflag)

        if(i>1.and.i<n2) then
            elt(nelt).adj(2)=nelt+1
            elt(nelt).adj(4)=nelt-1
        else
            if(isclose) then
                if(i==1) elt(nelt).adj(4)=n2
                if(i==n1) elt(nelt).adj(2)=1
            else
                if(i==1) elt(nelt).adj(4)=-1
                if(i==n1) elt(nelt).adj(2)=-1                
            endif
        endif
        
        !update edges
        
        nedge=nedge+1
        if(nedge>size(edge,dim=1)) call enlarge_ar(edge,100)
        edge(nedge)=edge(edge2(i))
        edge(nedge).v=node(edge(edge2(i)).v).subbw
        edge(nedge).e=[ielt1,ielt2]
        
            
    enddo
        
    
 
    node.subbw=0
    node.layer=0
endsubroutine


FUNCTION find_ELT_ON_HEADING_LEFT(V,BE1,BE2) RESULT(ELT1)
!返回以V为顶点，BE1-BE2前进方向左边的所有单元,BE1,BE2是以V公共顶点的边
!如果没有，返回ELT1(1)=-1
IMPLICIT NONE
INTEGER,INTENT(IN)::V,BE1,BE2
INTEGER,ALLOCATABLE::ELT1(:)

INTEGER::I,J,N1,IELT1,ELT2(100),NC1,IELT2=0,IEDGE0,V0
REAL(8)::T1
!FOUND THE FIRST ELEMENT
IF(EDGE(BE1).V(1)==V) THEN
    V0=EDGE(BE1).V(2)
ELSEIF(EDGE(BE1).V(2)==V) THEN
    V0=EDGE(BE1).V(1)
ELSE
    PRINT *,'THE EDGE BE1 SHOULD INCLUDE THE VETEX V. BE1,V=',BE1,V
    STOP
ENDIF

IELT2=0
DO I=1,2
    IELT1=EDGE(BE1).E(I)
    IF(IELT1>0.AND.IELT2==0) THEN
        DO J=1,3
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
    NC1=NC1+1
    ELT2(NC1)=IELT2   
    DO WHILE(.TRUE.)
        N1=MOD(IEDGE0,3)+1
        IEDGE0=ELT(IELT2).EDGE(N1)
        IF(IEDGE0==BE2) EXIT
        IELT2=ELT(IELT2).ADJ(N1)
        IF(IELT2<1) EXIT
    ENDDO
    ELT1=ELT2(1:NC1)
ELSE
    ELT1=[-1]
ENDIF



endfunction

end module
    
    
    