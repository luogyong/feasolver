module DS_Gmsh2Solver
    
    INTEGER::NLAYER=1
    INTEGER,ALLOCATABLE::SUBLAYER(:)
    LOGICAL::ISGENLAYER=.FALSE.
    INTEGER::NOPOPUP=0
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
		integer::inode=0
		real(8)::xy(3)=0
		!INTEGER::CHILD=0
	end type
	type(node_type),allocatable::node(:)
	integer::nnode=0,TNNODE=0
	integer,allocatable::Noutputorder(:)

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
		integer::isini=0
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
		integer::group
		integer::ndim=0 !=0,point load; =1,line load; =2, planar load; =3,volume load;
		integer::dof
		integer::sf=0
		real(8)::value=0  !当输入seepageface时，value=1,2,3 分别表示节点的水头值等于坐标x,y,z.
		integer::n1=0,n2=0 !for spgface output
	end type
	type(bcgroup_type),allocatable::elt_bc(:),elt_load(:),elt_spgface(:)
	integer::nelt_bc,nelt_load,nelt_spgface
	
	type bc_tydef
		integer::node,dof	!when body force is input, node equels to element	
		integer::sf=0 !step function for this load
		integer::isdead=0 !for seepage face iterative, =1,the condition is no longer on work. 
		real(8)::value=0
	end type
	type(bc_tydef),allocatable::nodalLoad(:),nodalBC(:),spgface(:)
	integer::nnodalLoad=0,maxnnodalLoad=1000
	integer::nnodalBC=0,maxnnodalBC=1000
	integer::nspgface=0,maxnspgface=1000
    
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
	
    CONTAINS
    
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
end module
