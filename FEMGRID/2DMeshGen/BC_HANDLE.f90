MODULE BC_HANDLE
    USE meshDS, ONLY: NNODE,TNODE,ZONE,ELT,NODE,modeldimension,strtoint,seg,segindex,edge,enlarge_ar
    USE ds_t, ONLY:arr_t
    implicit none

    private
    PUBLIC::ZoneBC_type,LineBC_type,zonebc,nzbc,linebc,nlbc
    
    !INTERFACE ENLARGE_AR
    !    MODULE PROCEDURE I_ENLARGE_AR,R_ENLARGE_AR
    !END INTERFACE 
    
	Type linebc_type
		integer::nnum,LOC,ILAYER,BCTYPE,NDIM,DOF,SF,ISFIELD=0 !control point number
		integer::NELT=0
		integer,allocatable::NBC(:) 
        real(8),allocatable::vbc(:) !vbc(nvb/4)
		integer,allocatable::ELT(:,:),ESEG(:) !line elt(2,nelt),ESEG(nelt),单元所属的线段。
        integer::ISWELLCONDITION=0 !是否为井的边界，或出溢面
		real(8),allocatable::lincof(:)
        INTEGER::NNODE=0
        integer,allocatable::NODE(:)
		REAL(8),ALLOCATABLE::NVAL(:)
        character(512):: helpstring= &
        'LINEBC的输入格式为: \n &
            & 1)线边界数(nLBC); \n & 
            & 2.1) 控制参数[NCP(边界控制点数),ILAYER,LOC,BCTYPE,NDIM,DOF,STEPFUNC(0),ISFIELD(0) \n &
            & 2.2) 控制点号(CP)*ncp \n &
            & 2.3) 边界值(VBC)*(MAX(ncp,4)) \\各点之间默认线插,ISFIELD==1,则输入A,B,C,D四个数即可。 \n &
            & NOTE: \n &
            &   a) NDIM=[0,1(11,12,13),2(21,22,23)] 分别表示点、线、面(由此线拉伸形成的面)的约束. 因为边界条件（如位移，水头等）不具叠加性，所以NDIM取值对边界条件无影响. \n &
            &      (21,22,23)分别表示面单元沿x,y,z平面的投影面积,2表全面积。(11,12,13)分别表示线单元沿x,y,z轴的投影长度，1表全长度。 \n &
            &   b) DOF=4(水头)，且VALUE=-999,表示节点的水头边界值为节点高程值。(模拟暴雨工况。) \n &
            &   c) DOF(I)=1,2,...,7,分别表示约束X,Y,Z,H,MX,MY,MZ. \n &
            &   d) [A,B,C,D] 平面场方程计算参数，边界值VALUE=A*X+B*Y+C*Z+D. \n &
            &   e) BCTYPE=[0,1,2] 分别表示位移类，力类，出溢面类边界。 \n &
            &   f) LOC=[1,0] 分别表示ILAYER指向的是第ilayer高程面的节点，否则指向的是第ilayer单元层节点(或面)。'C 
    contains
        procedure,nopass::help=>write_help
        procedure::initialize=>LineBC_initialize
        PROCEDURE::GETVALUE=>LineBC_interpolate
        PROCEDURE::READIN=>LineBC_read
        PROCEDURE::Set_BC=>LineBC_load_translate
        PROCEDURE::OUTPUT=>LineBC_write
	end type 
    type(linebc_type),allocatable::LineBC(:)
    integer::nLBC=0

    
    type zonebc_type
		integer::izone,ISFIELD=0
		INTEGER::ndim=0 !=0,point load; =1,line load; =2, planar load; =3,volume load;
		integer::dof
		integer::sf=0
		real(8)::value=0  !当输入seepageface时，value=1,2,3 分别表示节点的水头值等于坐标x,y,z.
		integer::n1=0,n2=0 !for spgface output
        real(8)::LFC(4)=0.d0 !FIELD=AX+BY+CY+D LFC()=[A,B,C,D]
        integer::ISWELLCONDITION=0 !是否为井的边界，或出溢面
        integer::ilayer=0,BCtype=-1!BCTYPE=0,bc,=1;load;=2,spgface
        INTEGER::LOC=1 !IF LOCATION=1, ILAYER指的是第ilayer高程面，否则指的是第ilayer单元层
        !具体的边界节点及其对应的值
        integer::nnode=0
        integer,allocatable::node(:)
        real(8),allocatable::nval(:)
        character(512):: helpstring= &
        'ZoneBC的输入格式为: \n &
            & 1)组数(nbc); \n & 
            & 2)[IZONE,ILAYER,LOC,BCTYPE,NDIM,DOF,STEPFUNC(0),VALUE,[A,B,C,D]]*nbc \n &
            & NOTE: \n &
            &   a) NDIM=[0,1(11,12,13),2(21,22,23),3] 分别表示点、线、面、体的约束. 因为边界条件（如位移，水头等）不具叠加性，所以NDIM取值对边界条件无影响. \n &
            &      (21,22,23)分别表示面单元沿x,y,z平面的投影面积,2表全面积。(11,12,13)分别表示线单元沿x,y,z轴的投影长度，1表全长度。 \n &
            &   b) DOF=4(水头)，且VALUE=-999,表示节点的水头边界值为节点高程值。(模拟暴雨工况。) \n &
            &   c) DOF(I)=1,2,...,7,分别表示约束X,Y,Z,H,MX,MY,MZ. \n &
            &   d) [A,B,C,D] 平面场方程计算参数，边界值VALUE=A*X+B*Y+C*Z+D. \n &
            &   e) BCTYPE=[0,1,2] 分别表示位移类，力类，出溢面类边界。 \n &
            &   f) IZONE, 为zone的下标(ZONE(IGROUP) \n &
            &   g) LOC=[1,0] 分别表示ILAYER指向的是第ilayer高程面，否则指向的是第ilayer单元层。'C
    CONTAINS
        procedure,nopass::help=>write_help
        PROCEDURE::GETVALUE=>LINEARFIELDCAL
        PROCEDURE::READIN=>ZoneBC_read
        PROCEDURE::set_bc=>ZoneBC_load_translate
        PROCEDURE::OUTPUT=>ZoneBC_write
	end type
	type(ZoneBC_type),allocatable::ZoneBC(:)
	integer::nZBC=0
    
	type et_type
		integer::nnode=0,nedge=0,nface=0		
		character(512)::description
		integer,allocatable::edge(:,:),face(:,:),FaceEdge(:,:)
		!edge(2,nedge),
		!face: use node index to represent face. face(0:4,nface),face(0,:)==3,triangular face,==4, quadrilateral face
		!FaceEdge: use edge index to represent face. FaceEdge(0:4,nface),FaceEdge(0,:)==3,triangular face, ==4, quadrilateral face
		real(8),allocatable::weight(:,:) !分布单元荷载各节点分布系数, weight(:,1) 平面均布荷载；weight(:,2) 平面三角形荷载;weight(:,3) 轴对称均布荷载；weight(:,4) 轴对称三角形荷载，!目前只能处理平面应变的两种情况
	contains
        procedure::GETGMSHET=>GETGMSHET
        procedure::GETWeight=>not_nodal_force_weight
	end type
	!type(et_type)::gmsh_elttype(100)    
    
    

CONTAINS

subroutine ZoneBC_write(this,unit) 
    implicit none
    class(ZoneBC_type)::this
    integer,intent(in)::unit   
    integer::i
    
    select case(this.bctype)
    case(0) !displacement type
        if(this.ISWELLCONDITION<1) then
            write(unit,130) this.nnode
        else
            write(unit,133) this.nnode
        endif
        write(unit, 132) 
		do i=1,this.nnode
			write(unit,131) node(this.node(i)).number,this.dof,this.nval(i),this.sf
		end do 
    case(1) !force
        write(unit,140) this.nnode
		write(unit, 142) 
		do i=1,this.nnode
			write(unit,141) node(this.node(i)).number,this.dof,this.nval(i),this.sf
		end do    
    case(2) !seeapge face
        if(this.ISWELLCONDITION<1) then
            write(unit,150) this.nnode,this.sf
        else
            write(unit,153) this.nnode,this.sf,this.ISWELLCONDITION
        endif
       
		write(unit, 152) 
		write(unit,151) (node(this.node(i)).number,i=1,this.nnode)        
    endselect

    
130 FORMAT(/'BC,NUM=',I7,',ISINC=0') 
131 FORMAT(I7,1X,I2,1X,E15.7,1X,I4)
132 FORMAT("// ","NODE DOF VALUE [STEPFUNC.]")

133 FORMAT(/'BC,NUM=',I7,',ISINC=0,ISWELLHEAD=1') 


140 FORMAT(/'LOAD,NUM=',I7,',ISINC=0')
141 FORMAT(I7,1X,I2,1X,E15.7,1X,I4)
142 FORMAT("// ","NODE DOF VALUE [STEPFUNC.] ")

150 FORMAT(/'SEEPAGE FACE,NUM=',I7,', sf=',I7) 
151 FORMAT(10(I7,1X))
152 FORMAT("// NODE")
153 FORMAT(/'SEEPAGE FACE,NUM=',I7,', sf=',I7,",ISWELLBORE=",I7)     
endsubroutine

subroutine lineBC_write(this,unit) 
    implicit none
    class(linebc_type)::this
    integer,intent(in)::unit   
    integer::i
    
    select case(this.bctype)
    case(0) !displacement type
        if(this.ISWELLCONDITION<1) then
            write(unit,130) this.nnode
        else
            write(unit,133) this.nnode
        endif
        write(unit, 132) 
		do i=1,this.nnode
			write(unit,131) node(this.node(i)).number,this.dof,this.nval(i),this.sf
		end do 
    case(1) !force
        write(unit,140) this.nnode
		write(unit, 142) 
		do i=1,this.nnode
			write(unit,141) node(this.node(i)).number,this.dof,this.nval(i),this.sf
		end do    
    case(2) !seeapge face
        if(this.ISWELLCONDITION<1) then
            write(unit,150) this.nnode,this.sf
        else
            write(unit,153) this.nnode,this.sf,this.ISWELLCONDITION
        endif
       
		write(unit, 152) 
		
		write(unit,151) (node(this.node(i)).number,i=1,this.nnode)
		        
    endselect

    
130 FORMAT(/'BC,NUM=',I7,',ISINC=0') 
131 FORMAT(I7,1X,I2,1X,E15.7,1X,I4)
132 FORMAT("// ","NODE DOF VALUE [STEPFUNC.]")

133 FORMAT(/'BC,NUM=',I7,',ISINC=0,ISWELLHEAD=1') 


140 FORMAT(/'LOAD,NUM=',I7,',ISINC=0')
141 FORMAT(I7,1X,I2,1X,E15.7,1X,I4)
142 FORMAT("// ","NODE DOF VALUE [STEPFUNC.] ")

150 FORMAT(/'SEEPAGE FACE,NUM=',I7,', sf=',I7) 
151 FORMAT(10(I7,1X))
152 FORMAT("// NODE")
153 FORMAT(/'SEEPAGE FACE,NUM=',I7,', sf=',I7,",ISWELLBORE=",I7)     
endsubroutine

subroutine LineBC_initialize(this)
implicit none
class(linebc_type)::this
integer::i,j,k,n1,n2,nedge1
integer,allocatable::edge1(:)
real(8)::t2

    !gen line elements
	n1=0
	do j=1,this.nnum-1
        call seg(segindex(this.NBC(j),this.NBC(j+1))).getparas(nedge=nedge1)
		this.nelt=this.nelt+nedge1 !element count
	end do
	allocate(this.elt(2,this.nelt),this.eseg(this.nelt))
	n2=0
	do j=1,this.nnum-1
		n1=segindex(this.NBC(j),this.NBC(j+1))
        edge1=seg(n1).get_edge()
        call seg(n1).getparas(nedge=nedge1)
		do i=1,nedge1
            n2=n2+1
            this.elt(:,n2)=edge(edge1(i)).v
            this.eseg(n2)=j
        enddo
	end do 

    !cal linof
    
    
	t2=0
	do j=1,this.nnum-1
		k=j+1
		n1=this.NBC(j)
		n2=this.NBC(k)
		t2=((arr_t(n1).x-arr_t(n2).x)**2+ &
					(arr_t(n1).y-arr_t(n2).y)**2)**0.5
		if(abs(t2)<1e-7) then
			print *, 'Error in LineBC_initialize'
			pause
		end if
		
		if(j==1) allocate(this.lincof(this.nnum-1))
		this.lincof(j)=(this.vbc(k)-this.vbc(j))/t2
	end do	

endsubroutine



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

!10 FORMAT('BC的输入格式为:' \ &
!            '1)组数(nbc);' \ & 
!            '2)[BASE,IGROUP,ILAYER,LOC,BCTYPE,NDIM,DOF,STEPFUNC(0),VALUE,[A,B,C,D]]*nbc' \ &
!            'NOTE:' \ &
!            '   a) NDIM=[0,1,2,3] 分别表示点、线、面、体的约束. 因为边界条件（如位移，水头等）不具叠加性，所以NDIM取值对边界条件无影响.' \ &
!            '   b) DOF=4(水头)，且VALUE=-999,表示节点的水头边界值为节点高程值。(模拟暴雨工况。)' \ &
!            '   c) DOF(I)=1,2,...,7,分别表示约束X,Y,Z,H,MX,MY,MZ.' \ &
!            '   d) [A,B,C,D] 平面场方程计算参数，边界值VALUE=A*X+B*Y+C*Z+D.' \ &
!            '   e) BCTYPE=[0,1,2] 分别表示位移类，力类，出溢面类边界。' \ &
!            '   f) BASE=[1,0] 分别表示边界是定义在ZONE上，还是在CL (control line)上。' \ &
!            '   g) LOC=[1,0] 分别表示ILAYER指向的是第ilayer高程面，否则指向的是第ilayer单元层。' \ &
!            '   h) IGROUP, IGROUP为zone的下标(ZONE(IGROUP))')
endsubroutine

subroutine  LineBC_read(this,unit)
    implicit none
    class(LineBC_type)::this
    integer,intent(in)::unit
    INTEGER::I,J,N1,DN=0,NINC1
    INTEGER,PARAMETER::DNMAX=300
    REAL(8)::AR(DNMAX)
	!print *,'Reading bc data...'
 !   call this.help(this.helpstring)       
	!call skipcomment(unit)
	!read(unit,*) nZBC
	!allocate(ZoneBC(nZBC))
	!do i=1,nZBC
	call strtoint(unit,ar,dnmax,dn,dnmax)
	!n1=I
    THIS.nnum=int(ar(1)) !izone
    THIS.ilayer=int(ar(2))
    THIS.LOC=int(ar(3))
    THIS.BCtype=int(ar(4))
	THIS.ndim=int(ar(5))
	THIS.dof=int(ar(6))
	THIS.sf=int(ar(7))
	THIS.isfield=int(ar(8))
    ALLOCATE(THIS.NBC(THIS.NNUM))
    IF(THIS.ISFIELD==0) THEN
        ALLOCATE(THIS.VBC(THIS.NNUM))
        
    ELSE
        ALLOCATE(THIS.VBC(MAX(THIS.NNUM,4)))
    ENDIF
    call strtoint(unit,ar,dnmax,dn,dnmax)
    THIS.NBC=AR(1:SIZE(THIS.NBC,DIM=1))
    
    call strtoint(unit,ar,dnmax,dn,dnmax)
    THIS.VBC=AR(1:SIZE(THIS.VBC,DIM=1))
    
	!end do

endsubroutine


subroutine  ZoneBC_read(this,unit)
    implicit none
    class(ZoneBC_type)::this
    integer,intent(in)::unit
    INTEGER::I,J,N1,DN=0,NINC1
    INTEGER,PARAMETER::DNMAX=300
    REAL(8)::AR(DNMAX)
	!print *,'Reading bc data...'
 !   call this.help(this.helpstring)       
	!call skipcomment(unit)
	!read(unit,*) nZBC
	!allocate(ZoneBC(nZBC))
	!do i=1,nZBC
	call strtoint(unit,ar,dnmax,dn,dnmax)
	!n1=I
    THIS.IZONE=int(ar(1)) !izone
    THIS.ilayer=int(ar(2))
    THIS.LOC=int(ar(3))
    THIS.BCtype=int(ar(4))
	THIS.ndim=int(ar(5))
	THIS.dof=int(ar(6))
	THIS.sf=int(ar(7))
	THIS.value=ar(8)
    if(dn>8) then
        IF(dn==12) THEN
            THIS.LFC(1:4)=AR(9:12)
            THIS.ISFIELD=1
        ELSE
            PRINT *, "INPUT NUMBERS IS NOT AS EXPECTED(8 OR 12). IZONE=",THIS.IZONE
            STOP
        ENDIF
    endif
	!end do

endsubroutine

subroutine ZoneBC_load_translate(this)
	
	implicit none
    class(ZoneBC_type)::this
	integer::i,j,k,n1,n2,nc,iar1(5),n3,NELT1,NINC1
	real(8)::LAV1,t1
	integer,allocatable::nodalload1(:),node1(:)
	real(8),allocatable::load1(:)
	INTEGER,POINTER::ELT1(:)
    TYPE(ET_TYPE)::ETWEIGHT1
	! 荷载具可叠加，边界条件不具有可加性。分布力转化为节点力时要积分，但分布位移转化为节点位移时不需积分，节点位移等于分布位移。
    
   
	allocate(nodalload1(tnode))
	if(.NOT.allocateD(load1)) allocate(load1(tnode))
    
	IF(THIS.LOC==1) THEN !ELEVATION SURFACE
        NELT1=ZONE(THIS.IZONE).NTRIE3N
        ELT1=>ZONE(THIS.IZONE).TRIE3N
        NINC1=NNODE*THIS.ILAYER
    ELSE
        NINC1=0
        IF(INDEX(ZONE(THIS.IZONE).SOLVER_ET(THIS.ILAYER),'prm')>0) THEN
            NELT1=ZONE(THIS.IZONE).NPRM(THIS.ILAYER)
            ELT1=>ZONE(THIS.IZONE).PRM(:,THIS.ILAYER)            
        ELSEIF(INDEX(ZONE(THIS.IZONE).SOLVER_ET(THIS.ILAYER),'tet')>0) THEN
            NELT1=ZONE(THIS.IZONE).ntet(THIS.ILAYER)
            ELT1=>ZONE(THIS.IZONE).tet(:,THIS.ILAYER)             
        ENDIF
    ENDIF
    
    CALL ETWEIGHT1.GETWEIGHT(ETWEIGHT1.GETGMSHET(ELT(elt1(1)).ET))
    
    SELECT CASE(THIS.BCTYPE) 
    CASE(1)        
		select case(THIS.ndim)
			case(0) !ndim=0,是点荷载，各个单元中的同一节点，不具可加性，均指同一值。
				nodalload1=0
				do j=1,nelt1
					n1=elt1(j)
					do k=1,elt(n1).nnum
						n2=elt(n1).node(k)+NINC1
						if(nodalload1(n2)==0) then							
                            if(.not.allocated(this.node)) then
                                allocate(this.node(100),THIS.nval(100))                                
                                this.nnode=0
                            endif
                            this.nnode=this.nnode+1
							if(this.nnode>size(this.node,dim=1)) THEN
                                call enlarge_AR(this.node,100)
                                call enlarge_AR(this.NVAL,100)
                            ENDIF
							this.node(this.nnode)=n2
							this.NVAL(this.nnode)=THIS.GETVALUE(node(n2).X,node(n2).Y,node(n2).Z)
							nodalload1(n2)=1
						end if
					end do
				end do
								
			case default
				load1=0
				do j=1,NELT1
					n1=elt1(j)
					call Element_LAV_Cal(n1,LAV1,this.ndim)
					t1=0.d0
					!equivalent uniform load
					do k=1,elt(n1).nnum
						n2=elt(n1).node(k)+ninc1
						t1=t1+THIS.getvalue(node(n2).X,node(n2).Y,node(n2).Z)								 				
					end do
					t1=t1/elt(n1).nnum
					LAV1=LAV1*t1
					do k=1,elt(n1).nnum
						n2=elt(n1).node(k)+NINC1
						load1(n2)=load1(n2)+ETWEIGHT1.WEIGHT(K,1)*LAV1 !均布荷载                        								 				
					end do
					
				end do				
				do j=1,tnode
					if(abs(load1(j))>1e-10) then
                        if(.not.allocated(this.node)) then
                            allocate(this.node(100),THIS.nval(100))                                
                            this.nnode=0
                        endif
                        this.nnode=this.nnode+1
						if(this.nnode>size(this.node,dim=1)) THEN
                            call enlarge_AR(this.node,100)
                            call enlarge_AR(this.NVAL,100)
                        ENDIF
						this.node(this.nnode)=j
						this.NVAL(this.nnode)=load1(j)  
                        
					end if
				end do

		end select
	CASE(0,2)
        nodalload1=0
		do j=1,nelt1
			n1=elt1(j)
			do k=1,elt(n1).nnum
				n2=elt(n1).node(k)+NINC1
				if(nodalload1(n2)==0) then							
                    if(.not.allocated(this.node)) then
                        allocate(this.node(100),THIS.nval(100))                                
                        this.nnode=0
                    endif
                    this.nnode=this.nnode+1
					if(this.nnode>size(this.node,dim=1)) THEN
                        call enlarge_AR(this.node,100)
                        call enlarge_AR(this.NVAL,100)
                    ENDIF
					this.node(this.nnode)=n2
					this.NVAL(this.nnode)=THIS.GETVALUE(node(n2).X,node(n2).Y,node(n2).Z)
					nodalload1(n2)=1
				end if
			end do
		end do
    CASE DEFAULT
        PRINT *, 'NONE SUCH BC TYPE IS EXPECTED. subroutine ZoneBC_load_translate()'
        STOP
    END SELECT
        

	if(allocated(nodalload1)) deallocate(nodalload1)
	if(allocated(load1)) deallocate(load1)
	
end subroutine

subroutine LineBC_load_translate(this)
	
	implicit none
    class(LineBC_type)::this
	integer::i,j,k,n1,n2,nc,iar1(5),n3,NELT1,NINC1,ET1,nnum1,node1(4)
	real(8)::LAV1,t1
	integer,allocatable::nodalload1(:)
	real(8),allocatable::load1(:)
	INTEGER,POINTER::ELT1(:)
    TYPE(ET_TYPE)::ETWEIGHT1
	! 荷载具可叠加，边界条件不具有可加性。分布力转化为节点力时要积分，但分布位移转化为节点位移时不需积分，节点位移等于分布位移。
    
   
	allocate(nodalload1(tnode))
	if(.NOT.allocateD(load1)) allocate(load1(tnode))
    
    NELT1=this.nelt
    NINC1=NNODE*THIS.ILAYER
	IF(THIS.LOC==1) THEN !ELEVATION SURFACE
        ET1=21 !一维线单元
        nnum1=2
        !ELT1=>ZONE(THIS.IZONE).TRIE3N
    ELSE
        ET1=42 !二维四边形单元
        Nnum1=4
    ENDIF
    
    CALL ETWEIGHT1.GETWEIGHT(ETWEIGHT1.GETGMSHET(et1))
    
    SELECT CASE(THIS.BCTYPE) 
    CASE(1)        
		select case(THIS.ndim)
			case(0) !ndim=0,是点荷载，各个单元中的同一节点，不具可加性，均指同一值。
				nodalload1=0
				do j=1,nelt1
                
					node1(1:2)=this.elt(:,j)+ninc1
                    if(et1==42) node1(3:4)=[node1(2),NODE1(1)]-NNODE
                    
					do k=1,nnum1
                    
						n2=node1(k)
						if(nodalload1(n2)==0) then							
                            if(.not.allocated(this.node)) then
                                allocate(this.node(100),THIS.nval(100))                                
                                this.nnode=0
                            endif
                            this.nnode=this.nnode+1
							if(this.nnode>size(this.node,dim=1)) THEN
                                call enlarge_AR(this.node,100)
                                call enlarge_AR(this.NVAL,100)
                            ENDIF
							this.node(this.nnode)=n2
							this.NVAL(this.nnode)=THIS.GETVALUE(this.eseg(j),node(n2).X,node(n2).Y,node(n2).Z)
							nodalload1(n2)=1
						end if
					end do
				end do
								
			case default
				load1=0
				do j=1,NELT1
					node1(1:2)=this.elt(:,j)+ninc1
                    if(et1==42) node1(3:4)=[node1(2),NODE1(1)]-NNODE
					call Element_LAV_Cal(n1,LAV1,this.ndim,node1,et1,nnum1)
					t1=0.d0
					!equivalent uniform load
					do k=1,nnum1
						n2=node1(k)
						t1=t1+THIS.getvalue(this.eseg(j),node(n2).X,node(n2).Y,node(n2).Z)								 				
					end do
					t1=t1/nnum1
					LAV1=LAV1*t1
					do k=1,nnum1
						n2=node1(k)
						load1(n2)=load1(n2)+ETWEIGHT1.WEIGHT(K,1)*LAV1 !均布荷载                        								 				
					end do
					
				end do				
				do j=1,tnode
					if(abs(load1(j))>1e-10) then
                        if(.not.allocated(this.node)) then
                            allocate(this.node(100),THIS.nval(100))                                
                            this.nnode=0
                        endif
                        this.nnode=this.nnode+1
						if(this.nnode>size(this.node,dim=1)) THEN
                            call enlarge_AR(this.node,100)
                            call enlarge_AR(this.NVAL,100)
                        ENDIF
						this.node(this.nnode)=j
						this.NVAL(this.nnode)=load1(j)  
                        
					end if
				end do

		end select
	CASE(0,2)
        nodalload1=0
		do j=1,nelt1
			node1(1:2)=this.elt(:,j)+ninc1
            if(et1==42) node1(3:4)=[node1(2),NODE1(1)]-NNODE
			do k=1,nnum1
				n2=node1(k)
				if(nodalload1(n2)==0) then							
                    if(.not.allocated(this.node)) then
                        allocate(this.node(100),THIS.nval(100))                                
                        this.nnode=0
                    endif
                    this.nnode=this.nnode+1
					if(this.nnode>size(this.node,dim=1)) THEN
                        call enlarge_AR(this.node,100)
                        call enlarge_AR(this.NVAL,100)
                    ENDIF
					this.node(this.nnode)=n2
					this.NVAL(this.nnode)=THIS.GETVALUE(this.eseg(j),node(n2).X,node(n2).Y,node(n2).Z)
					nodalload1(n2)=1
				end if
			end do
		end do
    CASE DEFAULT
        PRINT *, 'NONE SUCH BC TYPE IS EXPECTED. subroutine LineBC_load_translate()'
        STOP
    END SELECT
        

	if(allocated(nodalload1)) deallocate(nodalload1)
	if(allocated(load1)) deallocate(load1)
	
end subroutine



!subroutine SortByPath(igroup,ispgroup,Local_Node,NLNDE)
!	use DS_Gmsh2Solver
!	implicit none
!	integer,intent(in)::igroup,ispgroup,NLNDE
!	integer,intent(out)::Local_node(NLNDE)
!	integer::i,j,k,iar1(5),n1,n2
!	integer,allocatable::element1(:)
!	
!	
!	allocate(element1(ZONE(igroup).nel))
!	element1=0
!	Local_Node(1)=elt(ZONE(ispgroup).elt(1)).node(1)	
!	j=1	
!	do while(j<=NLNDE-1)
!		do k=1,ZONE(igroup).nel
!			if(element1(k)/=0) cycle
!			
!			n1=elt(ZONE(igroup).elt(k)).nnode
!			iar1(1:n1)=elt(ZONE(igroup).elt(k)).node(1:n1)
!			n1=iar1(1)-Local_Node(j)
!			n2=iar1(2)-Local_Node(j)
!			if(n1*n2/=0) cycle
!			
!			element1(k)=1
!			
!			if(n1==0) then
!				select case(ZONE(igroup).ET_GMSH)
!				case(1)
!					Local_Node(j+1)=iar1(2)
!					j=j+1
!				case(8)
!					Local_Node(j+1)=iar1(3)
!					Local_Node(j+2)=iar1(2)
!					j=j+2
!				case(27)
!					Local_Node(j+1)=iar1(3)
!					Local_Node(j+2)=iar1(4)
!					Local_Node(j+3)=iar1(5)
!					Local_Node(j+4)=iar1(2)
!					j=j+4
!				case default
!					print *, "The element is expected to be a line element. datapoint(i),i=",i
!					stop
!				end select
!			else
!				select case(ZONE(igroup).ET_GMSH)
!				case(1)
!					Local_Node(j+1)=iar1(1)
!					j=j+1
!				case(8)
!					Local_Node(j+1)=iar1(3)
!					Local_Node(j+2)=iar1(1)
!					j=j+2
!				case(27)
!					Local_Node(j+1)=iar1(5)
!					Local_Node(j+2)=iar1(4)
!					Local_Node(j+3)=iar1(3)
!					Local_Node(j+4)=iar1(1)
!					j=j+4
!				case default
!					print *, "The element is expected to be a line element. datapoint(i),i=",i
!					stop
!				end select						
!			end if
!			exit
!		end do
!	enddo	
!
!	deallocate(element1)
!
!end subroutine

subroutine Element_LAV_Cal(ienum,LAV,ndim,node1,et1,nnum1) !计算单元的长度、面积或体积
	!use DS_Gmsh2Solver
	implicit none	
	integer,intent(in)::ienum
    integer,intent(in),optional::ndim,node1(:),et1,nnum1
	real(8),intent(out)::LAV
	!REAL(8),EXTERNAL::TRIAREA,TETVOL
	integer::i,j,ia1(4,3),et2,node2(30),nnum2
	real(8)::a1(6,3)=0
    
    if(present(et1)) then
        et2=et1
        nnum2=nnum1
    else
        et2=elt(ienum).et
        nnum2=elt(ienum).nnum
    endif
    if(present(node1)) then
        node2(1:nnum2)=node1(1:nnum2)
    else
        node2(1:nnum2)=elt(ienum).node(1:nnum2)
    endif
    
	
	select case(et2)
		case(21) !Length for a Line element 
            do i=1,2
		        a1(i,:)=[node(node2(i)).x,node(node2(i)).y,node(node2(i)).z]
	        end do
            if(ndim==11) then
                a1(:,2)=0;a1(:,3)=0
            elseif(ndim==12) then
                a1(:,1)=0;a1(:,3)=0
            elseif(ndim==13) then
                a1(:,1)=0;a1(:,2)=0    
            endif
			LAV=norm2(a1(1,:)-a1(2,:))
		case(0,6,15) !三角形面积
						
			do i=1,3
				a1(i,:)=[node(node2(i)).x,node(node2(i)).y,node(node2(i)).z]
			end do
            if(ndim==21) then
                a1(:,1)=0
            elseif(ndim==22) then
                a1(:,2)=0
            elseif(ndim==23) then
                a1(:,3)=0   
            endif            
			LAV=TRIAREA(A1(1:3,1:3))
		case(42) !四边形面积
			do i=1,4
				a1(i,:)=[node(node2(i)).x,node(node2(i)).y,node(node2(i)).z]
			end do
            if(ndim==21) then
                a1(:,1)=0
            elseif(ndim==22) then
                a1(:,2)=0
            elseif(ndim==23) then
                a1(:,3)=0   
            endif  
            
			LAV=TRIAREA(A1(1:3,1:3))
			
			A1(2,:)=A1(1,:)
			LAV=LAV+TRIAREA(A1(2:4,1:3))
			
		case(43,103) !四面体体积

                do i=1,4
				    a1(i,:)=[node(node2(i)).x,node(node2(i)).y,node(node2(i)).z]
			    enddo
			    LAV=TETVOL(A1(1:4,1:3))
		    
		case(63,153)
            ia1=reshape([1,2,3,4,2,5,3,4,3,5,6,4],(/4,3/))
            LAV=0
			do j=1,3
                do i=1,4
				    a1(i,:)=[node(node2(ia1(i,j))).x,node(node2(ia1(i,j))).y,node(node2(ia1(i,j))).z]
			    enddo
			    LAV=LAV+TETVOL(A1(1:4,1:3))
		    end do            

		case default
			print *, "NONE SUCH AN ELEMENT TYPE WAS EXPECTED. TO BE IMPROVED. SUB Element_LAV_Cal"
			STOP
		
	end select

end subroutine

!xyz(1:3,:) 三角形三个节点的空间坐标,第一维为节点。
real(8) function TriArea(xyz)  
	implicit none
	real(8),intent(in)::xyz(3,3)
	integer::i,j
	real(8)::a1(3,3)
	
	a1=xyz
	do i=2,3
		do j=1,3
			a1(i,j)=xyz(i,j)-xyz(1,j)
		end do
	end do
	TriArea=(a1(2,2)*a1(3,3)-a1(3,2)*a1(2,3))**2+ &
		(a1(2,1)*a1(3,3)-a1(3,1)*a1(2,3))**2+ &
		(a1(2,1)*a1(3,2)-a1(3,1)*a1(2,2))**2
	TriArea=0.5*abs(TriArea)**0.5

end function

!xyz(1:4,:) 四面体四个节点的空间坐标,第一维为节点。
real(8) function TetVOL(xyz)  
	implicit none
	real(8),intent(in)::xyz(4,3)
	integer::i,j
	real(8)::a1(3,3)
		
	do i=2,4
		do j=1,3
			a1(i-1,j)=xyz(i,j)-xyz(1,j)
		end do
	end do
	TetVOL=a1(1,1)*(a1(2,2)*a1(3,3)-a1(3,2)*a1(2,3))- &
			 a1(1,2)*(a1(2,1)*a1(3,3)-a1(3,1)*a1(2,3))+ &
			 a1(1,3)*(a1(2,1)*a1(3,2)-a1(3,1)*a1(2,2))
	TetVOL=1./6.*ABS(TetVOL)
end function


REAL(8) FUNCTION LINEARFIELDCAL(BCG,X,Y,Z)
        IMPLICIT NONE
        CLASS(ZoneBC_type),INTENT(IN)::BCG
        REAL(8),INTENT(IN),OPTIONAL::X,Y,Z
        REAL(8)::X1=0.D0,Y1=0.D0,Z1=0.D0
        
        IF(BCG.ISFIELD==0) THEN
            LINEARFIELDCAL=BCG.VALUE
            IF(BCG.DOF==4) THEN
                IF(ABS(BCG.VALUE+999.0D0)<1.D-7) THEN
                    IF(modeldimension==3) THEN
                        LINEARFIELDCAL=Z
                    ELSE
                        LINEARFIELDCAL=Y
                    ENDIF
                ENDIF
            ENDIF
        ELSE
            IF(PRESENT(X)) THEN
                X1=X 
            ELSE
                X1=0.D0
            ENDIF
            IF(PRESENT(Y)) THEN
                Y1=Y 
            ELSE
                Y1=0.D0
            ENDIF
            IF(PRESENT(Z)) THEN
                Z1=Z 
            ELSE
                Z1=0.D0
            ENDIF       
            
            LINEARFIELDCAL=BCG.LFC(1)*X1+BCG.LFC(2)*Y1+BCG.LFC(3)*Z1+BCG.LFC(4)
        ENDIF
ENDFUNCTION

FUNCTION LineBC_interpolate(THIS,ISEG,X,Y,Z) RESULT(IVAL)
        IMPLICIT NONE
        CLASS(LINEBC_type),INTENT(IN)::THIS
        REAL(8),INTENT(IN)::X,Y,Z
        INTEGER,INTENT(IN)::ISEG
        REAL(8)::ival,T2
        INTEGER::N1
        IF(this.ISFIELD==0) THEN
            n1=this.NBC(ISEG)
		    t2=((arr_t(n1).x-X)**2+ &
					(arr_t(n1).y-Y)**2)**0.5
            IVAL=this.LINCOF(ISEG)*T2+THIS.VBC(ISEG)
            IF(this.DOF==4) THEN
                IF(ABS(IVAL+999.0D0)<1.D-7) THEN
                    IF(modeldimension==3) THEN
                        IVAL=Z
                    ELSE
                        IVAL=Y
                    ENDIF
                ENDIF
            ENDIF
        ELSE

            IVAL=this.VBC(1)*X+this.VBC(2)*Y+this.VBC(3)*Z+this.VBC(4)
        ENDIF
ENDFUNCTION



subroutine not_nodal_force_weight(this,et)
	
	implicit none
    class(et_type)::this
	integer,intent(in)::et
	
	select case(et)
		case(1)
			if(.not.allocated(this.weight)) allocate(this.weight(this.nnode,2))
			this.weight(:,1)=0.5
			this.weight(1,2)=1./6.
			this.weight(2,2)=1./3.
		case(2,3,4,5,6)
			if(.not.allocated(this.weight)) allocate(this.weight(this.nnode,1))
			this.weight=1./this.nnode
		case(8)
			if(.not.allocated(this.weight)) allocate(this.weight(this.nnode,2))
			this.weight(1:2,1)=1./6.
			this.weight(3,1)=2./3.
			this.weight(1,2)=0.
			this.weight(2,2)=1./6.
			this.weight(3,2)=1./3.
		case(9)
			if(.not.allocated(this.weight)) allocate(this.weight(this.nnode,1))
			this.weight=0.
			this.weight(4:6,1)=1./3.
		case(11)	
			if(.not.allocated(this.weight)) allocate(this.weight(this.nnode,1))
			this.weight=0.
			this.weight(1:4,1)=-1./20.
			this.weight(5:10,1)=1./5.
		case(16)
			if(.not.allocated(this.weight)) allocate(this.weight(this.nnode,1))
			this.weight=0.
			this.weight(1:4,1)=-1./12.
			this.weight(5:8,1)=1./3.
		case(23) !CPE15
			allocate(this.weight(this.nnode,1))
			this.weight=0.D0
			this.weight(4:6,1)=-1./45.
			this.weight(7:12,1)=4./45.
			this.weight(13:15,1)=8./45.
		case(27)
			if(.not.allocated(this.weight)) allocate(this.weight(this.nnode,2))
			this.weight(1,1)=7./90.
			this.weight(2,1)=7./90.
			this.weight(3,1)=16./45.
			this.weight(5,1)=16./45.
			this.weight(4,1)=2/15.
			this.weight(1,2)=0
			this.weight(2,2)=7./90.
			this.weight(3,2)=4/45.
			this.weight(5,2)=4./15.
			this.weight(4,2)=1/15.
		case default
			print *, 'To be improved. sub=not_nodal_force_weight.et=',et
			stop
		end select
    
end subroutine    

INTEGER FUNCTION GETGMSHET(THIS,ET)
    
    IMPLICIT NONE
    CLASS(ET_TYPE)::THIS
    INTEGER,INTENT(IN)::ET
    
    SELECT CASE(ET)
    !LINE,2-NODE LINE
    CASE(21)
        THIS.NNODE=2
        GETGMSHET=1
    CASE(0,32)
        THIS.NNODE=3
        GETGMSHET=2
    CASE(42)
        THIS.NNODE=4
        GETGMSHET=3 
    CASE(43)
        THIS.NNODE=4
        GETGMSHET=4
    CASE(83)
        THIS.NNODE=8
        GETGMSHET=5
    CASE(63)
        THIS.NNODE=6
        GETGMSHET=6
    CASE(6,62)
        THIS.NNODE=6
        GETGMSHET=9       
    CASE(103)
        THIS.NNODE=10
        GETGMSHET=11
    CASE(10)
        THIS.NNODE=1
        GETGMSHET=15
    CASE(82)
        THIS.NNODE=8
        GETGMSHET=16
    CASE(15,152)
        THIS.NNODE=15
        GETGMSHET=23
    CASE(153)
        THIS.NNODE=15
        GETGMSHET=18        
    CASE DEFAULT
        PRINT *, 'NO SUCH ELEMENT TYPE. FUNC=GETGMSHET.'
        PAUSE
    ENDSELECT
        
            
    
ENDFUNCTION
    
END MODULE