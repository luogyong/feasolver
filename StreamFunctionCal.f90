MODULE STREAMFUNCTION
!假定：
!1)渗透系数的主轴与坐标轴平行
!2)承压
!3)单联通域，无奇异点
    
    USE solverds, ONLY:ELEMENT,ENUM,NODE,NNUM,ECP,SPG2D,CAX_SPG,TDISP,SOLVER_CONTROL,DOFADJL_TYDEF,&
        UM,isIniSEdge,PI,flownet_file,title,esetid,eset,neset,timestep,cpe8_spg,cps8_spg,CAX8_spg,cpe8r_spg,cps8r_spg,CAX8R_spg,&
        cpe6_spg,cps6_spg,CAX6_spg,cpe15_spg,cps15_spg,CAX15_spg,bar,bar2d,beam,beam2d,ssp2d,sf,bc_disp,bc_load,BD_NUM,BL_NUM
    USE MESHGEO,ONLY:sedge,nsedge,snadjl,EDGE_TYDEF,NODE_ADJ_TYDEF,I_ENLARGE_AR,GETGMSHET,ELTTYPE
    USE SolverMath,ONLY:ANGLE2D
    USE MKL_DSS	
    IMPLICIT NONE
    
    PRIVATE
    INTEGER::MAXBLOOPS=100
    INTEGER::nzone_tec_sf=0,varsharezone_tec_sf=0
    logical::isfirstcall=.true.
    PUBLIC::STREAMFUNCTIONCAL   
    
    
    integer, parameter :: RPRE = 8

    !INTEGER,ALLOCATABLE::ISBCELT1(:),ISBCDOF1(:)
    
    !TYPE SF_KQ
    !    LOGICAL::ISINI=.FALSE.
    !    REAL(KIND=RPRE),ALLOCATABLE::K(:,:),QC(:,:)
    !CONTAINS
    !    PROCEDURE,PASS::INITIALIZE=>KQ_INI
    !END TYPE
    !TYPE(SF_KQ),ALLOCATABLE::SFKQ(:)
    !REAL(KIND=RPRE),ALLOCATABLE::Q1(:)
    
    type node_sf_tydef
        integer::node=0    
    end type
    type element_sf_tydef        
        integer::NEDGE=0
        LOGICAL::ISINI=.FALSE.
		integer,allocatable::node(:)
        integer,allocatable::EDGE(:),ADJELT(:) 
        REAL(KIND=RPRE),ALLOCATABLE::K(:,:),QC(:,:) 
    CONTAINS
        PROCEDURE,PASS::KQ_CAL=>SF_KQ_INI
        PROCEDURE::FORM_DOFADJL=>SF_FORM_DOFADJL
    end type    
    
    TYPE BLOOPS_TYDEF
        INTEGER::NNODE=0,NEDGE=0,ISOUTBC=0
        INTEGER::PCUT(2) !主边界与此边界相连路径的的端点
        INTEGER,ALLOCATABLE::NODE(:),EDGE(:),EDGECUT(:)  !注意，node仅存储边界的主节点，不包含单元边的中间节点。  
        REAL(8)::BBOX(2,3)    
    ENDTYPE

    
    TYPE M2S_DOMAIN_TYDEF
        INTEGER::nnum,enum,NEDGE,NBLOOPS=0
        type(node_sf_tydef),allocatable::node(:)
        type(element_sf_tydef),allocatable::element(:)
        TYPE(EDGE_TYDEF),ALLOCATABLE::EDGE(:)
        TYPE(NODE_ADJ_TYDEF),ALLOCATABLE::NADJL(:)
        TYPE(BLOOPS_TYDEF),ALLOCATABLE::BLOOPS(:)
        !REAL(KIND=RPRE),ALLOCATABLE::Q(:),H(:)
    contains
        procedure::initialize=>m2s_initialize
        procedure,private::split_edge=>m2s_split_edge
        PROCEDURE,private::get_bc_loop=>m2s_get_bc_loop  
        procedure,private::m2sdomain=>m2sdomain_transform
        PROCEDURE,private::GET_PATH=>M2S_GETPATH
        
    END TYPE
    
    type(M2S_DOMAIN_TYDEF)::model_sf
        
    TYPE MKL_SOLVER_TYDEF
        INTEGER::NNZ=0,NDOF=0
        INTEGER,ALLOCATABLE::ROWINDEX(:),JCOL(:),diaglkmloc(:),ADOF(:)
        TYPE(DOFADJL_TYDEF),ALLOCATABLE::dofadjl(:)
        REAL(KIND=RPRE),ALLOCATABLE::KM(:),Q(:),H(:) !Q一开始是方程右端项,后面为求解结果。
        TYPE(MKL_DSS_HANDLE) :: handle
    CONTAINS
        PROCEDURE::SOLVER_INI=>MKL_SOLVER_INI
        PROCEDURE::assemble_km=>assemble_km_SF
        PROCEDURE::bc=>bc_SF
        PROCEDURE::solve=>solve_SF
        PROCEDURE::Q_CAL=>SF_Q_CAL
    ENDTYPE
    TYPE(MKL_SOLVER_TYDEF)::SOLVER_SF
    
    CONTAINS
    
    FUNCTION STREAMFUNCTIONCAL(iincs,isubts) RESULT(VSF)
        INTEGER,INTENT(IN)::IINCS,ISUBTS
        INTEGER::I,J,K,OUTBC1,N1
        REAL(KIND=RPRE),ALLOCATABLE::VSF(:)

        
        if(.not.isIniSEdge) CALL Model_MESHTOPO_INI()
        
        CALL MODEL_SF.INITIALIZE()
        
        ALLOCATE(VSF(MODEL_SF.NNUM))
        
        CALL SOLVER_SF.SOLVER_INI(MODEL_SF)
        CALL SOLVER_SF.assemble_km(MODEL_SF)
        CALL SOLVER_SF.q_cal(MODEL_SF)
        CALL SOLVER_SF.bc(MODEL_SF)
        CALL SOLVER_SF.solve()
        
        CALL outdata_flownet(IINCS,ISUBTS) 
        
        IF(MODEL_SF.NBLOOPS>1) then
            VSF=0.D0
        ELSE
            WHERE(NODE.DOF(4)>0)  VSF(NODE.DOF(4))=SOLVER_SF.Q(MODEL_SF.NODE.NODE)
        endif
    
        
        
        
    END FUNCTION
    

    

SUBROUTINE SF_Q_CAL(SELF,MODEL)
    class(MKL_SOLVER_TYDEF)::self
    CLASS(M2S_DOMAIN_TYDEF)::MODEL
    INTEGER::I
        
    IF(.NOT.ALLOCATED(self.Q)) ALLOCATE(self.Q(SELF.NDOF))
    IF(.NOT.ALLOCATED(self.H)) ALLOCATE(self.H(SELF.NDOF))
    self.Q=0.D0
    self.h=TDISP(MODEL.node.node)
    DO I=1,MODEL.ENUM
        if(element(I).isactive==0) CYCLE              
        SELF.Q(MODEL.ELEMENT(I).NODE)=SELF.Q(MODEL.ELEMENT(I).NODE)+ &
            MATMUL(MODEL.ELEMENT(I).QC,SELF.H(MODEL.ELEMENT(I).NODE))            
    ENDDO
    
    
END SUBROUTINE

    
 !!convert multiple connected domain to simple connected domain.
SUBROUTINE m2sdomain_transform(self)
    class(M2S_DOMAIN_TYDEF)::self
    INTEGER::I,J,OUTBC1,N1,N2,N3,N4,N5,N6,N7,N8
    REAL(KIND=RPRE)::XC1(3)
    INTEGER,ALLOCATABLE::BEDGE1(:),IA1(:)
        
    OUTBC1=MAXLOC(SELF.BLOOPS(1:SELF.NBLOOPS).ISOUTBC,DIM=1)
    !find internal singularities
    ALLOCATE(IA1(SELF.NNUM))
    IA1=0
    DO I=1,SELF.NBLOOPS
        IA1(SELF.NODE(SELF.BLOOPS(I).NODE).NODE)=-1
        DO J=1,SELF.BLOOPS(I).NEDGE
            N1=SELF.BLOOPS(I).EDGE(J)
            IF(SELF.EDGE(N1).NMIDPNT>0) THEN
                IA1(SELF.NODE(SELF.EDGE(N1).MIDPNT).NODE)=-1
            ENDIF
        ENDDO
    ENDDO
    IF(BD_NUM>0) then
        WHERE(IA1(BC_DISP.NODE)==0.and.abs(node(BC_DISP.NODE).q)>1.e-6) IA1(BC_DISP.NODE)=BC_DISP.NODE
    ENDIF
    IF(BL_NUM>0) then
        WHERE(IA1(BC_LOAD.NODE)==0.and.abs(node(BC_LOAD.NODE).q)>1.e-6) IA1(BC_LOAD.NODE)=BC_LOAD.NODE
    ENDIF
    
    IA1=PACK(IA1,IA1>0)    
    DO J=1,SIZE(IA1)        
        SELF.NBLOOPS=SELF.NBLOOPS+1
        SELF.BLOOPS(SELF.NBLOOPS).NNODE=1
        SELF.BLOOPS(SELF.NBLOOPS).NODE=[IA1(J)]
    ENDDO
    
    IF(SELF.NBLOOPS>1) THEN
        DO I=1,SELF.NBLOOPS
            IF(I==OUTBC1) CYCLE
            IF(SELF.BLOOPS(I).NNODE>1) THEN
                DO J=1,2
                    XC1(J)=SUM(NODE(SELF.BLOOPS(I).NODE(1:SELF.BLOOPS(I).NNODE-1)).COORD(J))/(SELF.BLOOPS(I).NNODE-1)
                ENDDO
                N8=0
            ELSE
                XC1(1:2)=NODE(SELF.BLOOPS(I).NODE(1)).COORD(1:2)
                N8=1
            ENDIF
            
            
            N1=MINLOC(NORM2(RESHAPE([NODE(SELF.NODE(SELF.BLOOPS(0).NODE).NODE).COORD(1)-XC1(1),&
                NODE(SELF.NODE(SELF.BLOOPS(0).NODE).NODE).COORD(2)-XC1(2)],([SELF.BLOOPS(0).NNODE,2])),DIM=2),DIM=1)
            N2=MINLOC(NORM2(RESHAPE([NODE(SELF.NODE(SELF.BLOOPS(I).NODE).NODE).COORD(1)-NODE(SELF.NODE(SELF.BLOOPS(0).NODE(N1)).NODE).COORD(1),&
                NODE(SELF.NODE(SELF.BLOOPS(I).NODE).NODE).COORD(2)-NODE(SELF.NODE(SELF.BLOOPS(0).NODE(N1)).NODE).COORD(2)],([SELF.BLOOPS(I).NNODE,2])),DIM=2),DIM=1)
            SELF.BLOOPS(I).PCUT=[SELF.BLOOPS(0).NODE(N1),SELF.BLOOPS(I).NODE(N2)]
            SELF.BLOOPS(I).EDGECUT=SELF.GET_PATH(SELF.BLOOPS(I).PCUT(1),SELF.BLOOPS(I).PCUT(2))
            
            DO j=0,SIZE(SELF.BLOOPS(I).EDGECUT)-N8
                IF(j<1) THEN
                    N6=N1-1
                    IF(N6<1) N3=SELF.BLOOPS(0).NEDGE
                ELSE
                    N6=N6+1
                ENDIF
                N3=SELF.BLOOPS(0).EDGE(N6)
                 
                IF(j<SIZE(SELF.BLOOPS(I).EDGECUT)) THEN
                    N4=SELF.BLOOPS(I).EDGECUT(j+1)
                ELSE
                    N4=SELF.BLOOPS(I).EDGE(N2)
                ENDIF
                
                
                BEDGE1=SELF.split_edge(N3,N4)
                IF(NORM2(REAL(BEDGE1))>0) THEN
                    N5=SIZE(BEDGE1)
                    IF(N5==1) THEN
                        SELF.BLOOPS(0).EDGE=[SELF.BLOOPS(0).EDGE(:N6),BEDGE1,ABS(N4),SELF.BLOOPS(0).EDGE(N6+1:)]
                        
                        SELF.BLOOPS(0).NEDGE=SELF.BLOOPS(0).NEDGE+2
                    ELSEIF(N5==2) THEN
                        SELF.BLOOPS(0).EDGE=[SELF.BLOOPS(0).EDGE(:N6-1),BEDGE1,ABS(N4),ABS(N3),SELF.BLOOPS(0).EDGE(N6+1:)]
                        SELF.BLOOPS(0).NEDGE=SELF.BLOOPS(0).NEDGE+3
                    ELSE
                        STOP "UNEXPECTED RESULT. SUB=m2sdomain_transform"
                    ENDIF
                    N7=BEDGE1(N5)
                ENDIF
                
                
                
            ENDDO
            
            !ADD BLOOPS(I).EDGE TO BLOOPS(0).EDGE
            IF(SELF.BLOOPS(I).NEDGE>0) THEN
                N3=MINLOC(ABS(SELF.BLOOPS(0).EDGE-N7),DIM=1)
                SELF.BLOOPS(0).EDGE=[SELF.BLOOPS(0).EDGE(:N3),SELF.BLOOPS(I).EDGE,SELF.BLOOPS(0).EDGE(N3+1:)]
                SELF.BLOOPS(0).NEDGE=SELF.BLOOPS(0).NEDGE+SELF.BLOOPS(I).NEDGE
            ENDIF
            !UPDATE BLOOPS(0).NODE
            SELF.BLOOPS(0).NNODE=SELF.BLOOPS(0).NEDGE+1
            SELF.BLOOPS(0).NODE=[SELF.EDGE(SELF.BLOOPS(0).EDGE).V(1),SELF.EDGE(SELF.BLOOPS(0).EDGE(1)).V(1)]
                
            
            
                
        ENDDO
        
    ENDIF 
    
    DEALLOCATE(IA1)
    
    RETURN
    
END SUBROUTINE
!
!subroutine bcedge_merge(ibloop)
!    integer,intent(in)::ibloop
!    
!    integer::i,j,k,n1
!    
!    n1=minloc(abs(bloops(0).node-bloops(ibloop).pcut(1)),dim=1)
!    
!    
!endsubroutine
    

FUNCTION M2S_GETPATH(SELF,S1,E1) RESULT(PATH)
!返回连接节点S1和E1之间的单元边。
    class(M2S_DOMAIN_TYDEF)::self
    INTEGER,INTENT(IN)::S1,E1
    INTEGER,ALLOCATABLE::PATH(:)
    REAL(8)::ALPHA1
    INTEGER::N1,IEDGE1,N2,MAXLENGTH1=1000,I
    REAL(8),ALLOCATABLE::ANGLE1(:)
    INTEGER,ALLOCATABLE::PATH1(:)
    N1=S1
    
    DO WHILE(.TRUE.)
        ALPHA1=angle2d(NODE(self.node(n1).node).COORD(1:2),NODE(self.node(E1).node).COORD(1:2))
        ANGLE1=self.edge(ABS(Self.NADJL(N1).EDGE)).ANGLE
        DO I=1,SELF.NADJL(N1).NNUM
            IF(Self.NADJL(N1).EDGE(I)<0) THEN 
                IF(ANGLE1(I)>=0.D0) THEN
                    ANGLE1(I)=ANGLE1(I)-PI()
                ELSE
                    ANGLE1(I)=ANGLE1(I)+PI()
                ENDIF
            ENDIF
         END DO   
        ANGLE1=ALPHA1-ANGLE1
        IEDGE1=MOD(minloc(ABS([ANGLE1,ANGLE1+2*PI(),ANGLE1-2*PI()]),DIM=1)-1,SELF.NADJL(N1).NNUM)+1
        N2=Self.NADJL(N1).EDGE(IEDGE1)
        PATH1=[PATH1,N2]
        IF(N2>0) THEN
            N1=self.edge(N2).V(2)
        ELSE
            N1=Self.EDGE(-N2).V(1)
        ENDIF
        IF(N1==E1) EXIT      
        IF(SIZE(PATH1)>MAXLENGTH1) THEN
            PRINT *, 'PATH SEEMS TOO LONG.PLEASE CONFIRM.'
            MAXLENGTH1=MAXLENGTH1+100
            PAUSE
        ENDIF
    ENDDO
    
    PATH=PATH1
    
END FUNCTION   
    
    
    subroutine m2s_initialize(self)
        class(M2S_DOMAIN_TYDEF)::self
        integer::i
        self.nnum=nnum
        self.enum=enum
        self.nedge=nsedge
        allocate(self.node(nnum),self.element(enum))
        do i=1,enum
            self.element(i).node=element(i).node
            self.element(i).edge=element(i).edge
            self.element(i).nedge=size(self.element(i).edge)
            !self.element(i).g=element(i).g
            self.element(i).ADJELT=element(i).ADJELT
            call self.element(i).KQ_CAL(i)
        enddo
        do i=1,nnum
            self.node(i).node=i
            !self.node(i).dof=node(i).dof(4)
        enddo
        
        self.edge=sedge(:nsedge)
        self.nadjl=snadjl(:nnum)
        
        CALL SELF.get_bc_loop()
        CALL SELF.m2sdomain()
        

        
    end subroutine  
    
    FUNCTION m2s_split_edge(self,eg1,eg2) RESULT(BEDGE)
        !更新EG1-EG2前进方向左边单元节点V为新节点Vc,v为EG1-EG2的共有节点，并更新网格的关系
        !返回新生成的由eg2分裂而成的新的edge。
        class(M2S_DOMAIN_TYDEF)::self
        integer,intent(in)::eg1,eg2
        INTEGER,ALLOCATABLE::BEDGE(:)
        integer::v,elt1,edge1,n1,i,j,n2,sign1,n3,v0,be1,be2,keynode1,ET1
        integer,allocatable::elt2(:),edge2(:),ar1(:),ar2(:)
        
        be1=abs(eg1);be2=abs(eg2)
        if(eg1<0) then
            v=self.edge(be1).v(1)
            v0=self.edge(be1).v(2)
            sign1=-1
        else
            v=self.edge(be1).v(2)
            v0=self.edge(be1).v(1)
            sign1=1
        endif
        !new a node
        self.nnum=self.nnum+1
        keynode1=self.nnum
        self.node=[self.node,self.node(v)]
        self.nadjl=[self.nadjl,self.nadjl(v)]
        self.nadjl(keynode1).NNUM=0;self.nadjl(keynode1).ENUM=0
	    if(allocated(self.nadjl(keynode1).NODE)) deallocate(self.nadjl(keynode1).NODE)
        if(allocated(self.nadjl(keynode1).element)) deallocate(self.nadjl(keynode1).element)
        if(allocated(self.nadjl(keynode1).subid)) deallocate(self.nadjl(keynode1).subid)
        if(allocated(self.nadjl(keynode1).edge)) deallocate(self.nadjl(keynode1).edge)

        
        !找到以be1为边的左边单元
        if(self.edge(be1).enum>0.and.self.edge(be1).element(1)*sign1>0) then
            elt1=self.edge(be1).element(1)
            n1=self.edge(be1).subid(1)
            n2=self.edge(be1).subid(2)
        elseif(self.edge(be1).enum>1.and.self.edge(be1).element(2)*sign1>0) then
            elt1=self.edge(be1).element(2)
            n1=self.edge(be1).subid(2)
            n2=self.edge(be1).subid(1)
        else
            stop "No element-on-the-heading-left attached to edge BE1"
        endif
        
        edge1=be1
        
        !如果be1不是边界，则分离之，令其成为边界
        if(self.edge(be1).enum==2) then
            n3=self.element(elt1).adjelt(n1)
            self.element(n3).adjelt(n2)=0
            self.element(elt1).adjelt(n1)=0
            self.edge(be1).enum=1
            self.edge(be1).element(1)=n3
            self.edge(be1).subid(1)=n2
            self.edge(be1).v=[v,v0]
            
            !new an edge
            self.edge=[self.edge(:self.nedge),self.edge(be1)]
            self.nedge=self.nedge+1
            self.edge(self.nedge).v=[v0,v]
            self.element(elt1).edge(n1)=self.nedge
            self.edge(self.nedge).element(1)=elt1
            self.edge(self.nedge).subid(1)=n1
            self.edge(self.nedge).enum=1 
            
            self.nadjl(v0).edge=[self.nadjl(v0).edge(:self.nadjl(v0).nnum),self.nedge]
            self.nadjl(v).edge=[self.nadjl(v).edge(:self.nadjl(v).nnum),-self.nedge]
            self.nadjl(v0).node=[self.nadjl(v0).node(:self.nadjl(v0).nnum),v]
            self.nadjl(v).node=[self.nadjl(v).node(:self.nadjl(v).nnum),v0]
            self.nadjl(v).nnum=self.nadjl(v).nnum+1
            
            edge1=self.nedge
            BEDGE=ADD_ELT(BEDGE,EDGE1)
            
            call update_edge_midpnt()
            
        elseif(self.edge(be1).enum>2) then
            stop "For a 2D case , one side has at most 2 elements."  
        endif
        
        
        do while(.true.)

            !转化时，单元的数量不变，除少数单元的节点和自由度编号有变外，单元的其他属性保持不变。
            call update_edge1()
            
            where(self.element(elt1).node==v) self.element(elt1).node=keynode1
            self.nadjl(keynode1).element=ADD_ELT(self.nadjl(keynode1).element,elt1)
            self.nadjl(keynode1).subid=ADD_ELT(self.nadjl(keynode1).subid,minloc(abs(self.element(elt1).node-keynode1),dim=1))

            
            
            if(self.edge(edge1).element(1)==elt1) then
                n2=1
            else
                n2=2
            endif
            n1=mod(self.edge(edge1).subid(n2),self.element(elt1).nedge)+1 !
            edge1=abs(self.element(elt1).edge(n1))
            
            if(edge1==be2) then
               !update
                if(self.edge(be2).enum>1) then
                   !new an edge
                    self.edge=[self.edge(:self.nedge),self.edge(be2)]
                    self.nedge=self.nedge+1               
                    if(self.edge(self.nedge).v(1)/=v) then
                        self.edge(self.nedge).v(2)=self.edge(self.nedge).v(1)
                    endif
                    self.edge(self.nedge).v(1)=keynode1
                
                    self.element(elt1).edge(n1)=self.nedge
                    self.edge(self.nedge).element(1)=elt1
                    self.edge(self.nedge).subid(1)=n1
                    self.edge(self.nedge).enum=1
                    BEDGE=ADD_ELT(BEDGE,SELF.NEDGE)
                    self.nadjl(keynode1).edge=ADD_ELT(self.nadjl(keynode1).edge,self.nedge)                
                    self.nadjl(keynode1).node=ADD_ELT(self.nadjl(keynode1).node,self.edge(self.nedge).v(2))
                
                    n2=self.edge(self.nedge).v(2)
                    self.nadjl(n2).node=[self.nadjl(n2).node(:self.nadjl(n2).nnum),keynode1]
                    self.nadjl(n2).edge=[self.nadjl(n2).edge(:self.nadjl(n2).nnum),-self.nedge]
                    self.nadjl(n2).nnum=self.nadjl(n2).nnum+1
                    self.nadjl(keynode1).nnum=size(self.nadjl(keynode1).node)
                    self.nadjl(keynode1).enum=size(self.nadjl(keynode1).element)
                    
                    call update_edge_midpnt()
                else
                    call update_edge1()
                endif
               
                !update nadjl(v)
                if(allocated(ar1)) deallocate(ar1)
                n2=max(self.nedge,self.enum,keynode1)
                allocate(ar1(n2))
                ar1=0
                ar1(abs(self.nadjl(keynode1).edge(1:self.nadjl(keynode1).nnum)))=1
                where(ar1(abs(self.nadjl(v).edge(1:self.nadjl(v).nnum)))/=1) ar1(abs(self.nadjl(v).edge(1:self.nadjl(v).nnum)))=2
                self.nadjl(v).edge=pack(self.nadjl(v).edge(:self.nadjl(v).nnum),ar1(abs(self.nadjl(v).edge(:self.nadjl(v).nnum)))==2)
                
                ar1(self.nadjl(keynode1).element(1:self.nadjl(keynode1).enum))=3
                where(ar1(self.nadjl(v).element(1:self.nadjl(v).enum))/=3)  ar1(self.nadjl(v).element(1:self.nadjl(v).enum))=4
                ar2=pack([1:self.nadjl(v).enum],ar1(self.nadjl(v).element(1:self.nadjl(v).enum))==4)
                
                self.nadjl(v).element=self.nadjl(v).element(ar2)
                self.nadjl(v).subid=self.nadjl(v).subid(ar2)
                self.nadjl(v).nnum=size(self.nadjl(v).edge)
                self.nadjl(v).enum=size(self.nadjl(v).element)
                if(allocated(self.nadjl(v).node)) deallocate(self.nadjl(v).node)
                allocate(self.nadjl(v).node(self.nadjl(v).nnum))
                do i=1,self.nadjl(v).nnum
                    n2=self.nadjl(v).edge(i)
                    if(n2<0) then
                        self.nadjl(v).node(i)=self.edge(-n2).v(1)
                    else
                        self.nadjl(v).node(i)=self.edge(n2).v(2)
                    endif
                enddo
                
                
                
                !update the old edge Be2
                if(self.element(elt1).adjelt(n1)/=0)then
                    if(self.edge(be2).element(1)==elt1) then
                        self.edge(be2).element(1)=self.edge(be2).element(2)
                        self.edge(be2).subid(1)=self.edge(be2).subid(2)
                        !确保边界边(二维情况下只含一个单元的边)的方向与单元边的方向一致
                        !if(self.edge(be2).element(1)<0) then
                        !    n3=self.edge(be2).v(1);self.edge(be2).v(1)=self.edge(be2).v(2);self.edge(be2).v(2)=n3
                        !    self.edge(be2).element(1)=-self.edge(be2).element(1)
                        !endif                        
                    endif
                    self.edge(be2).element(2)=0
                    self.edge(be2).enum=1
                    self.element(self.edge(be2).element(1)).adjelt(self.edge(be2).subid(1))=0
                    self.element(elt1).adjelt(n1)=0                   
                endif
               
                exit     
            endif
            
            elt1=self.element(elt1).adjelt(n1)
        enddo
        
        IF(.NOT.ALLOCATED(BEDGE)) BEDGE=[0]
    contains
    
    subroutine update_edge1()
        integer::n1,n2,sign1        
        if(self.edge(edge1).v(1)==v) then
            n1=1;n2=2;sign1=1
        else
            n1=2;n2=1;sign1=-1           
        endif
        
        self.edge(edge1).v(n1)=keynode1
        if(allocated(self.nadjl(keynode1).edge)) then
            self.nadjl(keynode1).edge=[self.nadjl(keynode1).edge(1:self.nadjl(keynode1).nnum),sign1*edge1]
            self.nadjl(keynode1).node=[self.nadjl(keynode1).node(1:self.nadjl(keynode1).nnum),self.edge(edge1).v(n2)]
        else
            self.nadjl(keynode1).edge=[sign1*edge1]
            self.nadjl(keynode1).node=[self.edge(edge1).v(n2)]
        endif
        self.nadjl(keynode1).nnum=1+self.nadjl(keynode1).nnum       
        where(self.nadjl(self.edge(edge1).v(n2)).node==v) self.nadjl(self.edge(edge1).v(n2)).node=keynode1  
    
    endsubroutine
    
    subroutine update_edge_midpnt()
    
        if(self.edge(self.nedge).nmidpnt>0) then
            self.node=[self.node,self.node(self.edge(self.nedge).midpnt)]
            self.nadjl=[self.nadjl,self.nadjl(self.edge(self.nedge).midpnt)] !只更新空间，中间节点信息暂不更新。
            self.edge(self.nedge).midpnt=[self.nnum+1:self.nnum+self.edge(self.nedge).nmidpnt]
            ET1=GETGMSHET(element(elt1).ET)
		    SELF.ELEMENT(ELT1).NODE(ELTTYPE(ET1).MIDPNT(:,N1))=self.edge(self.nedge).midpnt
            self.nnum=self.nnum+self.edge(self.nedge).nmidpnt            
        endif
        
    end subroutine
      
    end FUNCTION
    
SUBROUTINE m2s_get_bc_loop(self)
    class(M2S_DOMAIN_TYDEF)::self
    INTEGER::I,J,N1,N2,N3,INODE1,OUTBC1=0,NSE_BC1
    INTEGER,ALLOCATABLE::IA1(:,:),ISACTIVE1(:),LOOP1(:),IA2(:,:),ELOOP1(:),SEDGE_BC1(:)
    REAL(8)::X1(2,4),XC1(3),T1
    
    NSE_BC1=COUNT(SELF.EDGE.ENUM==1)
    ALLOCATE(SEDGE_BC1(NSE_BC1),ISACTIVE1(NSE_BC1),LOOP1(NSE_BC1+1),ELOOP1(NSE_BC1),IA1(2,SELF.NNUM),IA2(2,SELF.NNUM))
    NSE_BC1=0;IA1=0;IA2=0
    ALLOCATE(SELF.BLOOPS(0:MAXBLOOPS))
    
    DO I=1,SELF.NEDGE
        IF(SELF.EDGE(I).ENUM==1) THEN
            NSE_BC1=NSE_BC1+1
            SEDGE_BC1(NSE_BC1)=I
            DO J=1,2
                N1=MINLOC(IA1(:,SELF.EDGE(I).V(J)),MASK=IA1(:,SELF.EDGE(I).V(J))==0,DIM=1)
                IF(N1==0) THEN
                    STOP "SIMPLE LOOP IS ASSUMED"
                ENDIF
                IA1(N1,SELF.EDGE(I).V(J))=SELF.EDGE(I).V(MOD(J,2)+1)
                !IA1(I,J)=与边界j节点相邻的两个边界主节点(不考虑边中间的节点)，
                !IA2(I,J)=与边界j节点相邻的两个边界边(不考虑边中间的节点)，
                IF(J==1) THEN
                    IA2(N1,SELF.EDGE(I).V(J))=NSE_BC1
                ELSE
                    IA2(N1,SELF.EDGE(I).V(J))=-NSE_BC1
                ENDIF
            ENDDO
        ENDIF
    ENDDO

    ISACTIVE1=1
    T1=1.D20
    DO I=1,NSE_BC1
        IF(ISACTIVE1(I)/=1) CYCLE
        INODE1=2
        IF(INODE1<3) THEN
            LOOP1(1:2)=SELF.EDGE(SEDGE_BC1(I)).V
            ELOOP1(1)=SEDGE_BC1(I)
            ISACTIVE1(I)=0
        ENDIF
        DO WHILE(.TRUE.) 
            N2=LOOP1(INODE1);N1=LOOP1(INODE1-1)        
            N3=MINLOC(IA1(:,N2)-N1,MASK=ABS(IA1(:,N2)-N1)>0,DIM=1)
            ELOOP1(INODE1)=SIGN(SEDGE_BC1(ABS(IA2(N3,N2))),IA2(N3,N2))
            LOOP1(INODE1+1)=IA1(N3,N2)
            INODE1=INODE1+1
            ISACTIVE1(ABS(IA2(N3,N2)))=0
            IF(LOOP1(1)==LOOP1(INODE1)) THEN   !假定边界都是闭合的环             
                SELF.NBLOOPS=SELF.NBLOOPS+1
                IF(SELF.NBLOOPS>MAXBLOOPS) THEN
                    PRINT *, "MAX SELF.BLOOPS IS ",MAXBLOOPS
                    STOP
                ENDIF
                SELF.BLOOPS(SELF.NBLOOPS).NODE=LOOP1(1:INODE1)
                SELF.BLOOPS(SELF.NBLOOPS).EDGE=ELOOP1(1:INODE1-1)
                !DO J=1,INODE1-1
                !    SELF.BLOOPS(SELF.NBLOOPS).NODE=[SELF.BLOOPS(SELF.NBLOOPS).NODE,LOOP1(J)]
                !ENDDO
                !SELF.BLOOPS(SELF.NBLOOPS).NODE=ADD_ELT(SELF.BLOOPS(SELF.NBLOOPS).NODE,LOOP1(1)) !首尾节点相同。
                !IF(ELOOP1(1)>0) THEN
                !    SELF.BLOOPS(SELF.NBLOOPS).ISOUTBC=1
                !    OUTBC1=SELF.NBLOOPS                    
                !ENDIF
                SELF.BLOOPS(SELF.NBLOOPS).NNODE=INODE1
                SELF.BLOOPS(SELF.NBLOOPS).NEDGE=INODE1-1
                DO J=1,3
                    SELF.BLOOPS(SELF.NBLOOPS).BBOX(1,J)=MINVAL(NODE(LOOP1(1:INODE1)).COORD(J))
                    SELF.BLOOPS(SELF.NBLOOPS).BBOX(2,J)=MAXVAL(NODE(LOOP1(1:INODE1)).COORD(J))
                    !最外围的边界
                    if(j==1.AND.SELF.BLOOPS(SELF.NBLOOPS).BBOX(1,J)<T1) THEN
                        T1=SELF.BLOOPS(SELF.NBLOOPS).BBOX(1,J)
                        OUTBC1=SELF.NBLOOPS 
                    ENDIF
                ENDDO
                EXIT
            ENDIF
        ENDDO
    ENDDO
    
    SELF.BLOOPS(OUTBC1).ISOUTBC=1
    SELF.BLOOPS(0)=SELF.BLOOPS(OUTBC1) !存多连通区域合并为单连通域的边界
    
    

    !CHECK    
    !IF(NODE_LOOP_BC(1)/=NODE_LOOP_BC(NNLBC)) THEN    
    !    STOP "ERROR.THE FIRST AND THE LAST NODE SHOULD BE IDENTICAL.SUB=SETUP_EDGE_BC"
    !ENDIF
    
    
    DEALLOCATE(IA1,IA2,LOOP1,ISACTIVE1)
    
END SUBROUTINE

    

    
    
    
    
    SUBROUTINE SF_KQ_INI(SELF,IELT)
        CLASS(element_sf_tydef)::SELF
        INTEGER,INTENT(IN)::IELT
        REAL(KIND=RPRE)::KM1(2,2),r1,DET1,t1
        REAL(KIND=RPRE),ALLOCATABLE::B1(:,:)
        INTEGER::I,J,K,NDOF1
    

        if(SELF.ISINI.OR.element(IELT).isactive==0) RETURN 
        if(element(iELT).ec==spg2d.or.element(iELT).ec==cax_spg)  then
            NDOF1=SIZE(SELF.NODE)
            allocate(SELF.K(NDOF1,NDOF1),SELF.QC(NDOF1,NDOF1))
            allocate(B1(2,NDOF1))

		    r1=1.0
		    if(element(ielt).ec==CAX_SPG) then
			    r1=abs(element(ielt).xygp(1,i))
			    if(abs(r1)<1e-7) r1=1.0e-7
		    end if            
            
            
          
            self.k=0.0D0
            self.qc=0.0d0
            
	        do i=1,element(IELT).ngp
		        self.k=self.k+ecp(element(IELT).et).weight(I)*element(IELT).DETJAC(I)*R1* MATMUL( &
				        TRANSPOSE(element(IELT).b(:,:,I)),element(IELT).b(:,:,I))
                
                
                KM1=ELEMENT(ielt).D*element(ielt).kr(i)
                !DET1=KM1(1,1)*KM1(2,2)-KM1(1,2)*KM1(2,1)
                !KM1=KM1/DET1
                T1=km1(1,1);km1(1,1)=km1(2,2);km1(2,2)=t1;
                B1(1,:)=element(IELT).b(2,:,I)
                B1(2,:)=-element(IELT).b(1,:,I)                
		        self.QC=self.QC+ecp(element(IELT).et).weight(I)*element(IELT).DETJAC(I)*R1* MATMUL( &
				        MATMUL(TRANSPOSE(element(IELT).b(:,:,I)),KM1),B1)	                
            ENDDO
            
        endif
        
        SELF.ISINI=.TRUE.
    
    ENDSUBROUTINE
    
    subroutine assemble_km_SF(SELF,MODEL)
    !************************************************************
    !MODEL.ELEMENT仅包含几何关系，单元的其他属性依赖solver中的element
    !model.element与solver.element一一对应
    !************************************************************
	    CLASS(MKL_SOLVER_TYDEF)::SELF
	    TYPE(M2S_DOMAIN_TYDEF)::MODEL    
	    integer::i,j,k,nj,n1=0
	    integer::loc,rdof,cdof,loc1
        !integer,allocatable::SELF.ADOF(:)
	    real(8)::var1,t1,t2
	    
        
	    IF(.NOT.ALLOCATED(SELF.KM)) ALLOCATE(SELF.KM(SELF.NNZ))
        IF(.NOT.ALLOCATED(SELF.ADOF)) ALLOCATE(SELF.ADOF(MODEL.NNUM))
	    n1=0
        SELF.ADOF=1
        self.km=0.d0
	    do i=1,MODEL.enum
            if(element(i).isactive==0) then
                SELF.ADOF(MODEL.ELEMENT(I).NODE)=0
                cycle 
            endif
		    !renew elements stiffness
		    t1=1.d0
			
		    !assemble the total matrix
		    do j=1,element(i).NNUM !row
			
			    rdof=MODEL.element(i).NODE(j)

			    do k=1,element(i).NNUM  !colume
				    cdof=MODEL.element(i).NODE(k)
				    if(solver_control.issym) then
					    if(solver_control.ismkl) then
						    !mkl solver, stored upper part
						    if(rdof>cdof) cycle						
					    else
						    !default solver, stored lower part
						    if(rdof<cdof) cycle  !if colume>row, skip for it is in upper part
					    end if
                    end if
                    
				    if(solver_control.ismkl) then
				
                        do loc=SELF.rowindex(rdof),SELF.rowindex(rdof+1)-1
                            if(self.jcol(loc)==cdof) exit
                        enddo
                    else
                        PRINT *, 'KML IS REQUIRED.'
                        
                    endif

				
				    SELF.km(loc)=SELF.km(loc)+t1*MODEL.ELEMENT(I).K(j,k)
                    if(isnan(SELF.km(loc))) then
                        print *, 'NAN VALUE IN ELEMENT(I).KM(J,K),I,J,K',I,J,K
                        STOP
                    endif
				
			    end do
		    end do
	    end do
	
	    !为避免对角线元素为零，保持正定，令没激活的自由度的对角线元素等于1，不会影响计算结果。
	    do i=1,MODEL.NNUM
		    if(SELF.ADOF(i)==0) then
			    SELF.km(SELF.diaglkmloc(i))=1.d0
		    endif	
	    enddo
	
        !DEALLOCATE(ADOF1)
        
    end subroutine    
    
    SUBROUTINE BC_SF(SELF,MODEL)
    	CLASS(MKL_SOLVER_TYDEF)::SELF
	    TYPE(M2S_DOMAIN_TYDEF)::MODEL 
        integer:: j,dof1,N1,n2,N3,outbc1
        real,allocatable::ra1(:)
        
       
        outbc1=0
        do j=1,model.bloops(outbc1).nnode-1
		    !if(bc_disp(j).isdead==1) cycle
            !if(na1(model_sf.bloops(1).node(j))/=0) cycle
		    dof1=model.bloops(outbc1).node(j)
            if(SELF.adof(dof1)==0) cycle
            
		    self.km(self.diaglkmloc(dof1))=UM
            !令原第一个水头边界节点的流函数的值为0
            SELF.Q(dof1)=0.d0
      !      if(solver_control.solver==N_R.OR.solver_control.solver==INISTIFF) then
			   ! load1(dof1)=(0.0d0-stepdis(dof1))*UM	!		
		    !else
			   ! load1(dof1)=0.0d0*UM			
		    !end if
            
            exit !一个即可,发现只施加一个节点的边界要比施加所有定流函数边界所得流网的质量要好。
	    end do
    ENDSUBROUTINE
    
    
    SUBROUTINE SOLVE_SF(self)
        
        CLASS(MKL_SOLVER_TYDEF)::SELF
	    TYPE(M2S_DOMAIN_TYDEF)::MODEL
		INTEGER::ERROR
        REAL(KIND=RPRE),ALLOCATABLE::SOL1(:)
        
        
        IF(.NOT.ALLOCATED(SOL1)) ALLOCATE(SOL1(SELF.NDOF))
        ! Factor the matrix.！MKL_DSS_POSITIVE_DEFINITE
        error = DSS_FACTOR_REAL( self.handle, MKL_DSS_POSITIVE_DEFINITE,self.km)
        IF (error /= MKL_DSS_SUCCESS) THEN
            !PRINT *, 'MKL_DSS_POSITIVE_DEFINITE FAILED.TRY MKL_DSS_INDEFINITE OPTION.'
			error = DSS_FACTOR_REAL( self.handle, MKL_DSS_INDEFINITE,self.km)
        ENDIF
                        
		IF (error /= MKL_DSS_SUCCESS) call mkl_solver_error(error)	
    
        
		error = DSS_SOLVE_REAL(self.handle, MKL_DSS_DEFAULTS, SELF.Q,1,SOL1 )
		IF (error /= MKL_DSS_SUCCESS) call mkl_solver_error(error)
        
		SELF.Q=SOL1 
        
        DEALLOCATE(SOL1)
    
    ENDSUBROUTINE
    
    FUNCTION ADD_ELT(AR1,ELT) RESULT(AR2)
    INTEGER,ALLOCATABLE,INTENT(IN)::AR1(:)
    INTEGER,INTENT(IN)::ELT
    INTEGER,ALLOCATABLE::AR2(:)
    
    IF(ALLOCATED(AR1)) THEN
        AR2=[AR1,ELT]
    ELSE
        AR2=[ELT]
    ENDIF
    END FUNCTION
    
    subroutine MKL_SOLVER_INI(self,MODEL)
        use quicksort        
        class(MKL_SOLVER_TYDEF)::self
        TYPE(M2S_DOMAIN_TYDEF)::MODEL
        integer::i,err,ndof1,n2,N1,J,error,perm(1)
        !integer,allocatable::ar1(:)
         

        
        self.nnz=0;SELF.NDOF=MODEL.NNUM
        IF(.NOT.ALLOCATED(SELF.DOFADJL)) ALLOCATE(SELF.DOFADJL(SELF.NDOF))
        DO I=1,MODEL.ENUM
            CALL MODEL.ELEMENT(I).FORM_DOFADJL(SELF.DOFADJL)
        ENDDO
        self.nnz=sum(self.dofadjl.ndof)
        !do i=1,SELF.NDOF
        !    call quick_sort(MODEL.NADJL(i).node(:MODEL.NADJL(i).nnum))
        !    SELF.DOFADJL(I).NDOF=count(MODEL.NADJL(i).node>=i)+1
        !    SELF.DOFADJL(I).DOF=[I,PACK(MODEL.NADJL(i).node,MODEL.NADJL(i).node>I)]
        !    self.nnz=self.nnz+SELF.DOFADJL(I).NDOF
        !enddo
	    IF(self.NNZ<0) THEN
            PRINT *, "THE NUMBER NONZERO ENTRIES IN KM IS TOO LARGE TO BE HANDEL BY INTEGER*4."
            STOP
        ENDIF

        !diaglkmloc(irow(i))=n1 !FOR SETTING THE BOUNDARY CONDITIONS 
	    if(.not.allocated(self.diaglkmloc))	allocate(self.diaglkmloc(MODEL.NNUM))	!假定一个节点只有一个自由度

        !rowindex for mkl format

        ndof1=MODEL.NNUM
	    allocate(self.jcol(self.nnz),self.rowindex(ndof1+1),STAT=ERR)
        n1=0;    
        DO I=1,ndof1
            DO J=1,SELF.DOFADJL(I).NDOF
                N1=N1+1
                IF(J==1) self.rowindex(I)=N1
                self.JCOL(N1)=SELF.DOFADJL(I).DOF(J)            
                IF(self.JCOL(N1)==I) self.diaglkmloc(I)=N1
            ENDDO                   
        ENDDO
        
	    self.rowindex(ndof1+1)=self.rowindex(ndof1)+SELF.DOFADJL(NDOF1).NDOF 
        
        

		! Initialize the solver.
		error = DSS_CREATE( SELF.handle, MKL_DSS_MSG_LVL_INFO + MKL_DSS_TERM_LVL_ERROR )
		IF (error /= MKL_DSS_SUCCESS) call mkl_solver_error(error)
		! Define the non-zero structure of the matrix.        
		if(solver_control.issym) then
			error = DSS_DEFINE_STRUCTURE( SELF.handle, MKL_DSS_SYMMETRIC, self.rowIndex, self.NDOF, self.NDOF, self.jcol, self.nnz )
		else
			error = DSS_DEFINE_STRUCTURE( SELF.handle, MKL_DSS_NON_SYMMETRIC, self.rowIndex, self.NDOF, self.NDOF, self.jcol, self.nnz )
		end if
		IF (error /= MKL_DSS_SUCCESS) call mkl_solver_error(error)
		! Reorder the matrix.
		error = DSS_REORDER( SELF.handle, MKL_DSS_DEFAULTS, perm )
		IF (error /= MKL_DSS_SUCCESS) call mkl_solver_error(error)	

        
        
        
	    !deallocate(ar1)
        
        
        
    end subroutine   
    
!according element connection, calculate the bandwidth.
subroutine SF_form_dofadjl(self,dofadjl)
	USE quicksort
    class(element_sf_tydef)::self
    type(dofadjl_tydef)::dofadjl(:)   

	integer::n1,n2,n3,n4,j,DOF1(100),I,ndof1
    

	!n1=minval(element(ienum).g)
	!n3=maxval(element(ienum).g)
    ndof1=size(self.node)
    DOF1(1:NDOF1)=self.node
    
    CALL quick_sort(DOF1(1:ndof1))
    n1=dof1(1);n3=dof1(ndof1)
	do j=1,ndof1
		if(solver_control.issym) then			
			if(.not.solver_control.ismkl) then
				!for default solver,存下三角,bw(i)为第i行下角形带宽。				
                DO I=1,J
                    CALL dofadjl(DOF1(J)).ADD_ITEM(DOF1(I))
                ENDDO
                
			else
				!for mkl solver, 存上三角，bw(i)为第i行上角形带宽
		
				 DO I=J,ndof1
                    CALL dofadjl(DOF1(J)).ADD_ITEM(DOF1(I))
                 ENDDO
			end if

		else
            DO I=1,ndof1
                CALL dofadjl(DOF1(J)).ADD_ITEM(DOF1(I))                
            ENDDO
		end if

    end do

    

end subroutine

subroutine outdata_flownet(iincs,isubts)	

	integer::iincs,isubts,file_unit,ieset
	logical::isbarfamily,isset1
	integer::i,j,nc,n1,k,k1,iset1
	character(1024)::cstring=''
	real(8)::t1


	
	file_unit=15
	if(isfirstcall) then 
		open(unit=file_unit,file=flownet_file,status='replace')
		!if(anybarfamily) open(unit=file_diagram,file=resultfile21,status='replace')		
		cstring='TITLE = '//'"'//trim(title)//'_flownet"'
		write(file_unit,'(a1024)') cstring
		!call tecplot_variables(cstring)
		write(file_unit,'(a)') 'Variables="X","Y","H","VSF"'
		isfirstcall=.false.
	else
		open(unit=file_unit,file=flownet_file,status='old',access='append')
	endif	

    
	call tecplot_zonetitle_flownet()
	isset1=.false.
	
	do i=1,neset
		iset1=esetid(i)
		if(sf(eset(iset1).sf).factor(iincs)==0) cycle
        IF(ESET(ISET1).COUPLESET>0 &            
            .AND.ESET(ISET1).COUPLESET<ISET1) CYCLE !附属单元不输出 !!!!!!!
		
		
		isbarfamily=eset(iset1).et==bar.or.eset(iset1).et==bar2d.or.eset(iset1).et==beam.or.eset(iset1).et==beam2d.or.eset(iset1).et==ssp2d
		
		if((.not.isbarfamily).and.(.not.isset1)) then
            write(file_unit,'(a1024)') eset(iset1).zonetitle_sf
            nc=4
			do j=1,model_sf.nnum
                n1=model_sf.node(j).node	
		        write(file_unit,999) node(n1).coord(1:2),solver_sf.h(j),solver_sf.Q(j)
            enddo

			isset1=.true.
			
		end if
		
		if(.not.eset(iset1).out_mesh_sf) cycle
		
		do j=eset(iset1).enums,eset(iset1).enume
			
			select case(eset(iset1).et)
				case(cpe8_spg,cps8_spg,CAX8_spg,cpe8r_spg,cps8r_spg,CAX8R_spg) 
					nc=4
					write(file_unit,9999) model_sf.element(j).node(8),model_sf.element(j).node(1),&
										model_sf.element(j).node(5),model_sf.element(j).node(5)
					write(file_unit,9999) model_sf.element(j).node(5),model_sf.element(j).node(2),&
										model_sf.element(j).node(6),model_sf.element(j).node(6)
					write(file_unit,9999) model_sf.element(j).node(6),model_sf.element(j).node(3),&
										model_sf.element(j).node(7),model_sf.element(j).node(7)
					write(file_unit,9999) model_sf.element(j).node(7),model_sf.element(j).node(4),&
										model_sf.element(j).node(8),model_sf.element(j).node(8)
					write(file_unit,9999) model_sf.element(j).node(5),model_sf.element(j).node(6),&
										model_sf.element(j).node(7),model_sf.element(j).node(8)								
					case(cpe6_spg,cps6_spg,CAX6_spg)
						nc=3
						write(file_unit,9999) model_sf.element(j).node(1),model_sf.element(j).node(4),&
											model_sf.element(j).node(6)
						write(file_unit,9999) model_sf.element(j).node(2),model_sf.element(j).node(5),&
											model_sf.element(j).node(4)
						write(file_unit,9999) model_sf.element(j).node(3),model_sf.element(j).node(6),&
											model_sf.element(j).node(5)
						write(file_unit,9999) model_sf.element(j).node(4),model_sf.element(j).node(5),&
											model_sf.element(j).node(6)
					case( cpe15_spg,cps15_spg,CAX15_spg)
						nc=3
						write(file_unit,9999) model_sf.element(j).node(3),model_sf.element(j).node(11),model_sf.element(j).node(10)
						write(file_unit,9999) model_sf.element(j).node(11),model_sf.element(j).node(6),model_sf.element(j).node(14)
						write(file_unit,9999) model_sf.element(j).node(11),model_sf.element(j).node(14),model_sf.element(j).node(10)
						write(file_unit,9999) model_sf.element(j).node(10),model_sf.element(j).node(14),model_sf.element(j).node(5)
						write(file_unit,9999) model_sf.element(j).node(6),model_sf.element(j).node(12),model_sf.element(j).node(15)
						write(file_unit,9999) model_sf.element(j).node(6),model_sf.element(j).node(15),model_sf.element(j).node(14)
						write(file_unit,9999) model_sf.element(j).node(14),model_sf.element(j).node(15),model_sf.element(j).node(13)
						write(file_unit,9999) model_sf.element(j).node(14),model_sf.element(j).node(13),model_sf.element(j).node(5)
						write(file_unit,9999) model_sf.element(j).node(5),model_sf.element(j).node(13),model_sf.element(j).node(9)
						write(file_unit,9999) model_sf.element(j).node(12),model_sf.element(j).node(1),model_sf.element(j).node(7)
						write(file_unit,9999) model_sf.element(j).node(12),model_sf.element(j).node(7),model_sf.element(j).node(15)
						write(file_unit,9999) model_sf.element(j).node(15),model_sf.element(j).node(7),model_sf.element(j).node(4)
						write(file_unit,9999) model_sf.element(j).node(15),model_sf.element(j).node(4),model_sf.element(j).node(13)
						write(file_unit,9999) model_sf.element(j).node(13),model_sf.element(j).node(4),model_sf.element(j).node(8)
						write(file_unit,9999) model_sf.element(j).node(13),model_sf.element(j).node(8),model_sf.element(j).node(9)
						write(file_unit,9999) model_sf.element(j).node(9),model_sf.element(j).node(8),model_sf.element(j).node(2)
                    case default
                        nc=size(model_sf.element(j).node)
						write(file_unit,9999) model_sf.element(j).node(1:nc)
			end select
			
		end do			
		
    end do
    

		
	close(file_unit)


	999     format(<nc>(E24.15,1X))
    9999    format(<nc>I7)
contains

subroutine tecplot_zonetitle_flownet()
	!use solverds
	integer::i,j,k,nc,n1,iset1
    real(kind=rpre)::t1=0.d0
	logical::isset1=.false.
	character(1024)::cstring='',cstring2='',cstring3=''
	character(48)::cword1='',cword2='',cword3='',cword4='',cword5='',cword6=''
	character(48)::cword7='',cword8='',cword9=''
	
	isset1=.false.
	do i=1,neset
        
		
        iset1=esetid(i)
        if(sf(eset(iset1).sf).factor(iincs)==0) cycle
		nzone_tec_sf=nzone_tec_sf+1
		write(cword1,*) i
		write(cword2,*) model_sf.nnum
		
		
		
		
		select case(eset(iset1).et)
			case(cpe6_spg,cps6_spg,CAX6_spg)
				n1=(eset(iset1).enume-eset(iset1).enums+1)*4
			case(cpe15_spg,cps15_spg,CAX15_spg)
				n1=(eset(iset1).enume-eset(iset1).enums+1)*16
			case(cpe8_spg,cps8_spg,CAX8_spg,cpe8r_spg,cps8r_spg,CAX8R_spg) 
				n1=(eset(iset1).enume-eset(iset1).enums+1)*5

			case default
				n1=eset(iset1).enume-eset(iset1).enums+1
		end select
		write(cword3,*) n1
		
        write(cword4,*) 'point'
		
		if(len_trim(eset(iset1).grouptitle)==0) then
			eset(iset1).zonetitle_sf ='ZONE,T=ESET'//trim(adjustL(cword1))//',N='//trim(adjustL(cword2))//',E=' &
					//trim(adjustL(cword3))//',ZONETYPE=' &
					//trim(adjustL(eset(iset1).stype))//',DATAPACKING='//trim(cword4) 
		else
			eset(iset1).zonetitle_sf ='ZONE,T='//trim(adjustL(eset(iset1).grouptitle))//',N='//trim(adjustL(cword2))//',E=' &
					//trim(adjustL(cword3))//',ZONETYPE=' &
					//trim(adjustL(eset(iset1).stype))//',DATAPACKING='//trim(cword4) 		
		endif


        t1=0.d0
        do j=1,iincs
            if(j<iincs) then
                n1=timestep(j).nsubts
            else
                n1=isubts
            end if
            do k=1,n1
                t1=t1+timestep(j).subts(k)
            end do			
        end do
		
		write(cword7,'(E15.7)') t1
		write(cword8,'(i7)') ISET1
		eset(iset1).zonetitle_sf=trim(eset(iset1).zonetitle_sf)//',StrandID='//trim(adjustL(cword8))//',Solutiontime='//trim(adjustL(cword7))
		if(eset(iset1).mesh_share_id_sf<1) then
			eset(iset1).mesh_share_id_sf=nzone_tec_sf
			eset(iset1).out_mesh_sf=.true.
		else
			write(cword1,*) eset(iset1).mesh_share_id_sf
			eset(iset1).zonetitle_sf=trim(eset(iset1).zonetitle_sf)//',connectivitysharezone='//trim(adjustL(cword1))
			eset(iset1).out_mesh_sf=.false.
		end if
		
		if(.not.isset1) then
			varsharezone_tec_sf=nzone_tec_sf
			isset1=.true.
		else
			write(cword1,*) 4
			write(cword2,*) varsharezone_tec_sf
			eset(iset1).zonetitle_sf=trim(eset(iset1).zonetitle_sf)//',VARSHARELIST=([1-'//trim(adjustL(cword1))//']=' &
												//trim(adjustL(cword2))//')'
		
		end if
		
	end do
end subroutine

end subroutine 
  
    
    
END MODULE
    