module DS_Gmsh2Solver
    

    INTERFACE ENLARGE_AR
        MODULE PROCEDURE I_ENLARGE_AR,R_ENLARGE_AR,NODE_ENLARGE_AR,ELEMENT_ENLARGE_AR,BCGROUP_ENLARGE_AR,&
                        FACE_TYDEF_ENLARGE_AR,EDGE_TYDEF_ENLARGE_AR,BCTYPE_ENLARGE_AR
    END INTERFACE    
    
    INTEGER::NLAYER=1
    INTEGER,ALLOCATABLE::SUBLAYER(:)
    LOGICAL::ISGENLAYER=.FALSE.
    INTEGER::NOPOPUP=1
    TYPE EDGE_TYDEF
        LOGICAL::ISINI=.FALSE.
        INTEGER::V(2)=0
        INTEGER::HKEY=-1
        CHARACTER(64)::CKEY=""
        INTEGER::NEL=0
        INTEGER,ALLOCATABLE::ELEMENT(:),SUBID(:) !SUBID IS INDEX WHICH EDGE OF THE ELEMENT IS THE EDGE
    ENDTYPE
    TYPE(EDGE_TYDEF),ALLOCATABLE::EDGE(:)
    TYPE FACE_TYDEF
        LOGICAL::ISINI=.FALSE.
        INTEGER::SHAPE=3
        INTEGER::V(4)=0,EDGE(4)=0
        INTEGER::HKEY=-1
        CHARACTER(64)::CKEY="" 		
        INTEGER::ISTRISURFACE=0 !<-1 与face反向
        REAL(8)::UNORMAL(3)=0.D0		
        INTEGER::NEL=0
        INTEGER,ALLOCATABLE::ELEMENT(:),SUBID(:) !SUBID IS INDEX WHICH FACE OF THE ELEMENT IS THE FACE       
    ENDTYPE
    TYPE(FACE_TYDEF),ALLOCATABLE::FACE(:)    
    INTEGER::NEDGE=0,NFACE=0,MAXNEDGE=10000,MAXNFACE=10000
    
	type node_type
		integer::inode=0,NEL=0
		real(8)::xy(3)=0
		INTEGER::ISDEAD=0
        INTEGER::N1=0 !FOR TEMPORARY USE ONLY.
        INTEGER,ALLOCATABLE::ELEMENT(:)
	end type
	type(node_type),allocatable::node(:)
	integer::nnode=0,TNNODE=0,gnnode=0
	integer,allocatable::Noutputorder(:),g2n(:)

	type element_type
        integer::et=-1
		integer::ntag=0
		integer,allocatable::tag(:)
		!number-of-tags
		!gives the number of integer tags that follow for the n-th element. 
		!By default, the first tag is the number of the physical entity to which the element belongs;
		!the second is the number of the elementary geometrical entity to which the element belongs; 
		!the third is the number of mesh partitions to which the element belongs, 
		!followed by the partition ids (negative partition ids indicate ghost cells). A zero tag is equivalent to no tag. 
        integer::ilayer=1
		integer::nnode=0,nedge=0,nface=0
		integer,allocatable::node(:)
        INTEGER,ALLOCATABLE::EDGE(:),FACE(:)
    !CONTAINS
    !    PROCEDURE::GET_EDGE=>SET_ELEMENT_EDGE
    !    PROCEDURE::GET_FACE=>SET_ELEMENT_FACE
	end type
	type(element_type),allocatable::element(:)
	integer::nel=0
	
	
	type physicalGroup_type
		integer::isini=0,SF=0
		integer::ndim,COUPLESET=0 !COUPLESET>0 AND <>itself，表此单元组的单元与physicalgroup(COUPLESET)的单元相同。 by default it was set to be itself
		logical::ismodel=.FALSE.,ISMASTER=.TRUE. !phgpnum(:) OUPLESET>0 AND <>itself, ismaster=.false. then let elemnt=coupleset.element
		integer::ET_GMSH=0 
		character(32)::name=''
		integer::nel=0
		integer,allocatable::element(:)
		integer::mat(50)=0
		character(32)::ET="ELT_BC_OR_LOAD"
		INTEGER::LAYERGROUP(50)=0	
        INTEGER,ALLOCATABLE::TRISURFACE(:),QUASURFACE(:),TRINODEG2L(:),TRINODEL2G(:) !global to local
        INTEGER::NTRISURFACE=0,NQUASURFACE=0,NTRINODEL2G=0
        REAL(8)::PROPERTY(3)=0.D0
	end type
	integer,parameter::maxphgp=10100
	type(physicalGroup_type)::physicalgroup(maxphgp)
	integer,allocatable::phgpnum(:)
	integer::nphgp=0
	
	type et_type
		integer::nnode=0,nedge=0,nface=0		
		character(512)::description
		integer,allocatable::edge(:,:),face(:,:),FaceEdge(:,:)
		!edge(2,nedge),
		!face: use node index to represent face. face(0:4,nface),face(0,:)==3,triangular face,==4, quadrilateral face
		!FaceEdge: use edge index to represent face. FaceEdge(0:4,nface),FaceEdge(0,:)==3,triangular face, ==4, quadrilateral face
		real(8),allocatable::weight(:,:) !分布单元荷载各节点分布系数, weight(:,1) 平面均布荷载；weight(:,2) 平面三角形荷载;weight(:,3) 轴对称均布荷载；weight(:,4) 轴对称三角形荷载，
													!目前只能处理平面应变的两种情况。
	end type
	type(et_type)::elttype(100)
	
	type bcgroup_type
		integer::group,ISFIELD=0,ISINPUT=0
		integer::ndim=0 !=0,point load; =1,line load; =2, planar load; =3,volume load;
		integer::dof
		integer::sf=0
		real(8)::value=0  !当输入seepageface时，value=1,2,3 分别表示节点的水头值等于坐标x,y,z.
		integer::n1=0,n2=0 !for spgface output
        real(8)::LFC(4)=0.d0 !FIELD=AX+BY+CY+D LFC()=[A,B,C,D]
        integer::ISWELLCONDITION=0 !是否为井的边界，或出溢面
        integer::spg_isdual=0
    CONTAINS
        PROCEDURE::GETVALUE=>LINEARFILEDCAL
	end type
	type(bcgroup_type),allocatable::elt_bc(:),elt_load(:),elt_spgface(:)
	integer::nelt_bc,nelt_load,nelt_spgface
	
	type bc_tydef
		integer::node,dof	!when body force is input, node equels to element	
		integer::sf=0 !step function for this load
		integer::isdead=0 !for seepage face iterative, =1,the condition is no longer on work. 
		real(8)::value=0
        integer::spg_isdual=0
	end type
	type(bc_tydef),allocatable::nodalLoad(:),nodalBC(:),spgface(:),Wellhead(:),WellSpgface(:)
	integer::nnodalLoad=0,maxnnodalLoad=1000
	integer::nnodalBC=0,maxnnodalBC=1000
	integer::nspgface=0,maxnspgface=1000
    integer::nwellhead=0,maxnwellhead=0
    integer::nwellspgface=0,maxnwellspgface=0
    
    type wsp_typdef
        integer::Group,chgroup,spgroup !chgroup, characteristic point.
        integer::xdirection=1
        integer::nnode=0,nchnode=0
        integer,allocatable::node(:)
		integer,allocatable::chnode(:)
    end type
    type(wsp_typdef),allocatable::wsp(:)
    integer::nwsp=0
    
    type out_data_typdef
        integer::Group
        integer::SPGroup !Startpoint Group,这个Group只有一个点。
		integer::order=0
        integer::issumq=0 !/=0, 表示仅输出各节点的流量和。
        integer::nnode=0
        integer,allocatable::node(:)        
    end type
    type(out_data_typdef),allocatable::DataPoint(:)
    integer::NDataPoint=0
	
	TYPE OFF_MESH_TYPDEF
		CHARACTER(512)::FILE=''
		INTEGER::PGPID=0
		INTEGER::NNODE=0,NEL=0,NEDGE=0
		integer,ALLOCATABLE::ELEMENT(:)
	ENDTYPE
	integer,parameter::MAXOFF=100
	integer,allocatable::OFFINDEX(:)
	integer::NOFF=0,NOFFNODE=0,NOFFFACE=0
	TYPE(OFF_MESH_TYPDEF)::OFFMESH(MAXOFF)
	REAL(8),ALLOCATABLE::OFFNODE(:,:)
	INTEGER,ALLOCATABLE::OFFFACE(:,:)
	REAL(8)::XYZTOL=1.D-4
	
	integer,allocatable::modelgroup(:)
	integer::nmodelgroup=0
		
	character(512)::resultfile,title,INCLUDEFILE(0:100),COPYFILE(100),FILEPATH,meshstructurefile
	INTEGER::NINCLUDEFILE=0,NCOPYFILE=0
		
	integer,allocatable::adjL(:,:)
	integer::maxadj=100
	
	integer::modeldimension=2
    integer::IsoutMS=0
    REAL(8),ALLOCATABLE::ELEVATION(:,:)
    
    !INTEGER::EDGE
    
    TYPE TRISURFACE_TYDEF
        INTEGER::NV,NF
        INTEGER,ALLOCATABLE::V(:),TRI(:)
    ENDTYPE
    TYPE(TRISURFACE_TYDEF),ALLOCATABLE::TRISURFACE(:)
    INTEGER::NTRISURFACE
    INTEGER::ISGENTRISURFACE=0,ISGEO=0,ISSTL=0,ISOFF=0
    
    TYPE WELLBORE_TYDEF
        INTEGER::IGP=0,WELLNODE=0,NWELLBORESEG=-1,BCTYPE=0,NWNODE=0,PIPEFLOW=-1,SPHERICALFLOW=0,NSEMI_SF=0
        INTEGER::NSPG_FACE=0,MAT=0
        REAL(8)::R,VALUE=0.D0
        INTEGER,ALLOCATABLE::SEMI_SF_IPG(:),SPG_FACE(:),SINK_NODE_SPG_FACE(:)
        REAL(8),ALLOCATABLE::DIR_VECTOR(:,:)        
    CONTAINS
        PROCEDURE::INITIALIZE=>WELLBORE_INITIALIZE
        
    ENDTYPE
    TYPE(WELLBORE_TYDEF),ALLOCATABLE::WELLBORE(:)
    INTEGER::NWELLBORE=0
    

    CONTAINS
    
    REAL(8) FUNCTION LINEARFILEDCAL(BCG,X,Y,Z)
        IMPLICIT NONE
        CLASS(bcgroup_type),INTENT(IN)::BCG
        REAL(8),INTENT(IN),OPTIONAL::X,Y,Z
        REAL(8)::X1=0.D0,Y1=0.D0,Z1=0.D0
        
        IF(BCG.ISFIELD==0) THEN
            LINEARFILEDCAL=BCG.VALUE
            IF(BCG.DOF==4) THEN
                IF(ABS(BCG.VALUE+999.0D0)<1.D-7) THEN
                    IF(modeldimension==3) THEN
                        LINEARFILEDCAL=Z
                    ELSE
                        LINEARFILEDCAL=Y
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
            
            LINEARFILEDCAL=BCG.LFC(1)*X1+BCG.LFC(2)*Y1+BCG.LFC(3)*Z1+BCG.LFC(4)
        ENDIF
    ENDFUNCTION
    
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
    SUBROUTINE FACE_TYDEF_ENLARGE_AR(AVAL,DSTEP)
        TYPE(FACE_TYDEF),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(FACE_TYDEF),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
    
        LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1,SOURCE=AVAL)
        DEALLOCATE(AVAL)
        ALLOCATE(AVAL(LB1:UB1+DSTEP))
        AVAL(LB1:UB1)=VAL1
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
    END SUBROUTINE 
    SUBROUTINE EDGE_TYDEF_ENLARGE_AR(AVAL,DSTEP)
        TYPE(EDGE_TYDEF),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(EDGE_TYDEF),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
    
        LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1,SOURCE=AVAL)
        DEALLOCATE(AVAL)
        ALLOCATE(AVAL(LB1:UB1+DSTEP))
        AVAL(LB1:UB1)=VAL1
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
    END SUBROUTINE
    SUBROUTINE NODE_ENLARGE_AR(AVAL,DSTEP)
        TYPE(node_type),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(node_type),ALLOCATABLE::VAL1(:)
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
        TYPE(element_type),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(element_type),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
    
        LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1,SOURCE=AVAL)
        DEALLOCATE(AVAL)
        ALLOCATE(AVAL(LB1:UB1+DSTEP))
        AVAL(LB1:UB1)=VAL1
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
    END SUBROUTINE
    SUBROUTINE BCGROUP_ENLARGE_AR(AVAL,DSTEP)
        TYPE(bcgroup_type),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(bcgroup_type),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0,ERR
    
        LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1,SOURCE=AVAL,STAT=ERR)
        DEALLOCATE(AVAL,STAT=ERR)
        ALLOCATE(AVAL(LB1:UB1+DSTEP),STAT=ERR)
        AVAL(LB1:UB1)=VAL1
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1,STAT=ERR)
    END SUBROUTINE
    SUBROUTINE BCTYPE_ENLARGE_AR(AVAL,DSTEP,UBAVAL)
        TYPE(bc_tydef),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        INTEGER,INTENT(IN OUT),OPTIONAL::UBAVAL
        TYPE(bc_tydef),ALLOCATABLE::VAL1(:)
        
        INTEGER::LB1=0,UB1=0,ERR
    
        LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
    
        ALLOCATE(VAL1,SOURCE=AVAL,STAT=ERR)
        IF(ALLOCATED(AVAL)) DEALLOCATE(AVAL)
        ALLOCATE(AVAL(LB1:UB1+DSTEP))
        AVAL(LB1:UB1)=VAL1
        IF(PRESENT(UBAVAL)) UBAVAL=UB1+DSTEP
        !AVAL(UB1+1:UB1+10)=0
        IF(ALLOCATED(VAL1)) DEALLOCATE(VAL1)
    END SUBROUTINE   
    
    SUBROUTINE WELLBORE_INITIALIZE(SELF)
        CLASS(WELLBORE_TYDEF)::SELF
        INTEGER::I,J,K,N1,N2,N3,N4,INODE1,IEL1,II1
        INTEGER::NODE1(10)
        INTEGER,ALLOCATABLE::ISWBE1(:)
        INTEGER,ALLOCATABLE::IA1(:)
        CHARACTER(16),PARAMETER::CAR1(3)=["PIPE2","WELLBORE","WELLBORE_SPGFACE"]
        
        !generate four-noded well element
        ! 1 +-----------+ 2
        !   |           |
        ! 4 +           + 3
        ! Node 1 and 2: Line element simulating well flow along the well.
        ! Element 1-4 and 2-3 : virtual branch element 
        
        
        
        NODE.N1=0
        N3=2+SELF.NSPG_FACE
        ALLOCATE(IA1(N3))
        IF(SELF.NSPG_FACE>0) THEN
            IA1=[SELF.IGP,SELF.PIPEFLOW,SELF.SPG_FACE]
        ELSE
            IA1=[SELF.IGP,SELF.PIPEFLOW]
        ENDIF

        DO IJ1=1,N3
            
            N1=IA1(IJ1)
            
            IF(N1<1) CYCLE
            
            IF(IJ1==1) THEN
                PHYSICALGROUP(N1).ET="WELLBORE"
            ELSEIF(IJ1==2) THEN
                PHYSICALGROUP(N1).ET="PIPE2"
            ELSE
                PHYSICALGROUP(N1).ET="WELLBORE_SPGFACE"
            ENDIF
            
            PHYSICALGROUP(N1).ISMODEL=.TRUE.
            PHYSICALGROUP(N1).MAT(1)=SELF.MAT
            
            IF(PHYSICALGROUP(N1).MAT(1)==0.AND.IA1(1)>0) THEN
                PHYSICALGROUP(N1).MAT(1)=PHYSICALGROUP(IA1(1)).MAT(1)  !!!!          
            ENDIF
            IF(PHYSICALGROUP(N1).MAT(1)==0.AND.IA1(2)>0) THEN
                PHYSICALGROUP(N1).MAT(1)=PHYSICALGROUP(IA1(2)).MAT(1)            
            ENDIF            
            IF(PHYSICALGROUP(N1).MAT(1)==0) THEN
                PRINT *,'PLEASE SIGN A MATID TO PHYSICALGROUP(N1) THROUGH "GROUPPARAMETER". N1=',N1
                STOP
            ENDIF
                        
            DO J=1,PHYSICALGROUP(N1).NEL
            
                IEL1=PHYSICALGROUP(N1).ELEMENT(J)
            
                DO K=1,ELEMENT(IEL1).NNODE
                    INNODE1=ELEMENT(IEL1).NODE(K)
                    IF(NODE(INNODE1).N1==0) THEN
                        NNODE=NNODE+1
                        NODE(INNODE1).N1=NNODE
                        IF(NNODE>UBOUND(NODE,DIM=1)) THEN
                            CALL ENLARGE_AR(NODE,1000)       
                        ENDIF
                        NODE(NNODE)=NODE(INNODE1)
                        NODE(NNODE).N1=INNODE1
                    ENDIF                
                ENDDO
            
                ELEMENT(IEL1).NODE=NODE(ELEMENT(IEL1).NODE).N1
                
                IF(IJ1==2) THEN
                    PHYSICALGROUP(N1).ET_GMSH=1
                ELSE
                    !滤管井单元,井流出溢面单元              
                    ELEMENT(IEL1).NNODE=2*ELEMENT(IEL1).NNODE
                    NODE1(1:ELEMENT(IEL1).NNODE)=[ELEMENT(IEL1).NODE,NODE(ELEMENT(IEL1).NODE).N1]
                    DEALLOCATE(ELEMENT(IEL1).NODE)
                    N2=ELEMENT(IEL1).NNODE/2
                    ALLOCATE(ELEMENT(IEL1).NODE,SOURCE=NODE1([1:N2,ELEMENT(IEL1).NNODE:N2+1:-1]))
                    PHYSICALGROUP(N1).ET_GMSH=3                
                ENDIF
                
                !SELECT CASE(IJ1)
                !CASE(1)                   
                ! !滤管井单元           
                !    ELEMENT(IEL1).NNODE=2*ELEMENT(IEL1).NNODE
                !    NODE1(1:ELEMENT(IEL1).NNODE)=[ELEMENT(IEL1).NODE,NODE(ELEMENT(IEL1).NODE).N1]
                !    DEALLOCATE(ELEMENT(IEL1).NODE)
                !    N2=ELEMENT(IEL1).NNODE/2
                !    ALLOCATE(ELEMENT(IEL1).NODE,SOURCE=NODE1([1:N2,ELEMENT(IEL1).NNODE:N2+1:-1]))
                !    PHYSICALGROUP(N1).ET_GMSH=3   
                !CASE(2) !管流单元
                !    !NOTHING TO DO
                !    PHYSICALGROUP(N1).ET_GMSH=1
                !CASE DEFAULT !井流出溢面单元                    
                !    ! :>>>>>>>>>>>>>>>>>>>>>>:
                !    ! 1...........2>>>>>>>>>>5(SinkNode(滤管最上面的边界节点))   (井壁侧)  !实际为1-4和2-3两个线单元，第5个节点只是为统计井流量用。
                !    ! |           |
                !    ! 4           3                                                
                !    ELEMENT(IEL1).NNODE=2*ELEMENT(IEL1).NNODE+1
                !    N4=ELEMENT(PHYSICALGROUP(SELF.SINK_NODE_SPG_FACE(IJ1-2)).ELEMENT(1)).NODE(1) !此集只含一个单节点的单元
                !    N4=NODE(N4).N1
                !    NODE1(1:ELEMENT(IEL1).NNODE)=[ELEMENT(IEL1).NODE,NODE(ELEMENT(IEL1).NODE).N1,N4]
                !    DEALLOCATE(ELEMENT(IEL1).NODE)
                !    N2=(ELEMENT(IEL1).NNODE-1)/2
                !    ALLOCATE(ELEMENT(IEL1).NODE,SOURCE=NODE1([1:N2,ELEMENT(IEL1).NNODE-1:N2+1:-1,ELEMENT(IEL1).NNODE]))
                !    PHYSICALGROUP(N1).ET_GMSH=3
                !ENDSELECT
            ENDDO
        
        ENDDO

     
        
        IF(SELF.SPHERICALFLOW>0) THEN
            N2=SELF.SPHERICALFLOW
            PHYSICALGROUP(N2).ISMODEL=.TRUE.
            IF(PHYSICALGROUP(N2).MAT(1)==0.AND.IA1(1)>0) THEN
                PHYSICALGROUP(N2).MAT=PHYSICALGROUP(IA1(1)).MAT
            ENDIF
            IF(PHYSICALGROUP(N2).MAT(1)==0.AND.IA1(2)>0) THEN
                PHYSICALGROUP(N2).MAT=PHYSICALGROUP(IA1(2)).MAT
            ENDIF  
            IF(PHYSICALGROUP(N2).MAT(1)==0) THEN
                PRINT *,'PLEASE SIGN A MATID TO PHYSICALGROUP(N2) THROUGH "GROUPPARAMETER". N2=',N2
                STOP
            ENDIF            
            PHYSICALGROUP(N2).ET="SPHFLOW"
            PHYSICALGROUP(N2).ET_GMSH=1
            DO I=1,PHYSICALGROUP(N2).NEL
                IEL1=PHYSICALGROUP(N2).ELEMENT(I)
                ELEMENT(IEL1).NNODE=2
                IF(NODE(ELEMENT(IEL1).NODE(1)).N1==0) THEN
                    STOP "GHOST NODE N1 SHOULD BE >0. Error in sub WELLBORE_INITIALIZE."
                ENDIF
                NODE1(1:ELEMENT(IEL1).NNODE)=[NODE(ELEMENT(IEL1).NODE).N1,ELEMENT(IEL1).NODE]
                DEALLOCATE(ELEMENT(IEL1).NODE)
                ALLOCATE(ELEMENT(IEL1).NODE,SOURCE=NODE1(1:2))                           
            ENDDO
        ENDIF
        
        DO I=1,SELF.NSEMI_SF            
            N2=SELF.SEMI_SF_IPG(I)
            PHYSICALGROUP(N2).ISMODEL=.TRUE.
            IF(PHYSICALGROUP(N2).MAT(1)==0.AND.IA1(1)>0) THEN
                PHYSICALGROUP(N2).MAT=PHYSICALGROUP(IA1(1)).MAT
            ENDIF
            IF(PHYSICALGROUP(N2).MAT(1)==0.AND.IA1(2)>0) THEN
                PHYSICALGROUP(N2).MAT=PHYSICALGROUP(IA1(2)).MAT
            ENDIF  
            IF(PHYSICALGROUP(N2).MAT(1)==0) THEN
                PRINT *,'PLEASE SIGN A MATID TO PHYSICALGROUP(N2) THROUGH "GROUPPARAMETER". N2=',N2
                STOP
            ENDIF 
            PHYSICALGROUP(N2).ET="SEMI_SPHFLOW"
            PHYSICALGROUP(N2).ET_GMSH=1
            PHYSICALGROUP(N2).PROPERTY(1:3)=SELF.DIR_VECTOR(:,I)
            DO J=1,PHYSICALGROUP(N2).NEL
                IEL1=PHYSICALGROUP(N2).ELEMENT(J)
                ELEMENT(IEL1).NNODE=2
                IF(NODE(ELEMENT(IEL1).NODE(1)).N1==0) THEN
                    STOP "GHOST NODE N1 SHOULD BE >0. Error in sub WELLBORE_INITIALIZE."
                ENDIF
                NODE1(1:ELEMENT(IEL1).NNODE)=[NODE(ELEMENT(IEL1).NODE).N1,ELEMENT(IEL1).NODE]
                DEALLOCATE(ELEMENT(IEL1).NODE)
                ALLOCATE(ELEMENT(IEL1).NODE,SOURCE=NODE1(1:2))                           
            ENDDO             
        ENDDO
        
        !BC
        N1=SELF.WELLNODE
        
        DO I=1,PHYSICALGROUP(N1).NEL 
            IEL1=PHYSICALGROUP(N1).ELEMENT(I)
            ELEMENT(IEL1).NODE=NODE(ELEMENT(IEL1).NODE).N1
        ENDDO
        IF(SELF.BCTYPE/=1) THEN
            
            NELT_BC=NELT_BC+1
            IF(NELT_BC>SIZE(ELT_BC,DIM=1)) CALL ENLARGE_AR(ELT_BC,10)
            
            ELT_BC(NELT_BC).GROUP=SELF.WELLNODE
            ELT_BC(NELT_BC).NDIM=0
            ELT_BC(NELT_BC).DOF=4
            ELT_BC(NELT_BC).VALUE=SELF.VALUE
            IF(SELF.BCTYPE==0) ELT_BC(NELT_BC).ISWELLCONDITION=1 !自流减压井水头边界，其出水量只出不进
        ELSE            
            NELT_LOAD=NELT_LOAD+1
            IF(NELT_LOAD>SIZE(ELT_LOAD,DIM=1)) CALL ENLARGE_AR(ELT_LOAD,10)
            ELT_LOAD(NELT_LOAD).GROUP=SELF.WELLNODE
            ELT_LOAD(NELT_LOAD).NDIM=0
            ELT_LOAD(NELT_LOAD).DOF=4
            ELT_LOAD(NELT_LOAD).VALUE=SELF.VALUE            
        ENDIF        
    
        DO I=1,SELF.NSPG_FACE
            IF(nelt_spgface+SELF.NSPG_FACE>SIZE(ELT_SPGFACE,DIM=1)) CALL ENLARGE_AR(ELT_SPGFACE,SELF.NSPG_FACE)
            nelt_spgface=nelt_spgface+1
            
            ELT_SPGFACE(NELT_SPGFACE).GROUP=SELF.SPG_FACE(I)
            ELT_SPGFACE(NELT_SPGFACE).NDIM=0
            ELT_SPGFACE(NELT_SPGFACE).DOF=4
            ELT_SPGFACE(NELT_SPGFACE).VALUE=SELF.VALUE            
            ELT_SPGFACE(NELT_SPGFACE).ISWELLCONDITION=SELF.SINK_NODE_SPG_FACE(I)
        ENDDO
        
        IF(ALLOCATED(ISWBE1)) DEALLOCATE(ISWBE1)
        IF(ALLOCATED(IA1)) DEALLOCATE(IA1)
    ENDSUBROUTINE
    
  !  SUBROUTINE SET_ELEMENT_EDGE(THIS)
  !      CLASS(element_type)::THIS
		!INTEGER::I
		!
		!DO I=1,elttype(THIS.ET).NEDGE
		!	
		!ENDDO
  !      
  !  END SUBROUTINE
  !  SUBROUTINE SET_ELEMENT_FACE(THIS)
  !      CLASS(element_type)::THIS
  !      
  !  END SUBROUTINE 
!translate all the characters in term into lowcase character string
    subroutine lowcase(term,iterm)
	    use dflib
	    implicit none
	    integer i,in,nA,nZ,nc,nd
        integer,optional::iterm
	    character(1)::ch
	    character::term*(*)
	
	    term=adjustl(trim(term))
	    nA=ichar('A')
	    nZ=ichar('Z')
	    nd=ichar('A')-ichar('a')
	    in=len_trim(term)
	    do i=1,in
		    ch=term(i:i)
		    nc=ichar(ch)
		    if(nc>=nA.and.nc<=nZ) then
			    term(i:i)=char(nc-nd)
		    end if
	    end do
    end subroutine    
    
    
end module
