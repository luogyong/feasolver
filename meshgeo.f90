MODULE MESHGEO

!USE solverds
USE POS_IO
USE hashtbl

!PRIVATE
!PUBLIC::EDGE,NEDGE,FACE,NFACE,SETUP_EDGE_TBL_TET,SETUP_FACE_TBL_TET,TET,NTET,&
!        SETUP_SUB_TET4_ELEMENT,&
!        MEDGE,NMEDGE,MFACE,NMFACE,SETUP_EDGE_TBL,SETUP_FACE_TBL,GETVAL,ProbeatPoint

LOGICAL::ISINI_GMSHET=.FALSE.

TYPE POINT_TYDEF
    REAL(8)::X(3)=0.D0
ENDTYPE

TYPE SEGMENT
    TYPE(POINT_TYDEF)::V1,V2
ENDTYPE


TYPE NODE_ADJ_TYDEF
    LOGICAL::ISINI=.FALSE.
    INTEGER::NNUM=0,ENUM=0
	INTEGER::ISDEAD=0
    INTEGER,ALLOCATABLE::NODE(:),ELEMENT(:),SUBID(:) !SUBID IS INDEX WHICH EDGE OF THE ELEMENT IS THE NODE
ENDTYPE
TYPE(NODE_ADJ_TYDEF),ALLOCATABLE::SNADJL(:)

TYPE EDGE_TYDEF
    LOGICAL::ISINI=.FALSE.
    INTEGER::V(2)=0
    INTEGER::HKEY=-1
    CHARACTER(64)::CKEY=""
    REAL(8)::DIS=0.D0
    INTEGER::ENUM=0
	INTEGER::ISDEAD=0
    INTEGER,ALLOCATABLE::ELEMENT(:),SUBID(:) !SUBID IS INDEX WHICH EDGE OF THE ELEMENT IS THE EDGE
ENDTYPE
TYPE(EDGE_TYDEF),ALLOCATABLE::EDGE(:),MEDGE(:),SEDGE(:)
TYPE FACE_TYDEF
    LOGICAL::ISINI=.FALSE.
    INTEGER::SHAPE=3
    INTEGER::V(4)=0,EDGE(4)=0 !IF EDGE(I)<0 MEANING ITS ORDER IS REVERSE .
    INTEGER::HKEY=-1
    CHARACTER(64)::CKEY="" 		
    INTEGER::ISTRISURFACE=0 !<-1 ��face����
    REAL(8)::BBOX(2,3)=0.D0 !UNORMAL(3)=0.D0,
    INTEGER::ENUM=0
	INTEGER::ISDEAD=0
    INTEGER,ALLOCATABLE::ELEMENT(:),SUBID(:) !SUBID IS INDEX WHICH FACE OF THE ELEMENT IS THE FACE
ENDTYPE
TYPE(FACE_TYDEF),ALLOCATABLE::FACE(:),MFACE(:)    
INTEGER::NEDGE=0,NFACE=0,MAXNEDGE=10000,MAXNFACE=10000
INTEGER::NMEDGE=0,NMFACE=0,MAXNMEDGE=10000,MAXNMFACE=10000,NSEDGE=0,MAXNSEDGE=10000 !M FOR MODEL,S FOR SOLVER MODEL 

!split all element into simple 3-noded-triangle , 4-noded-tetrahedron,point and 2-noded-line element
TYPE TET_TYDEF
	INTEGER::MOTHER=0,GMET=0,DIM=-1,NV=0,NE=0,NF=0,SF=0,ISDEAD=0,ISET=0
	INTEGER::V(4)=0,E(6)=0,F(4)=0	!F������û��ά��������ͨ��elttype(4).face��ȷ����
	REAL(8)::BBOX(2,3)=0.D0 !MIN,MAX OF X,Y AND Z
    REAL(8)::STEPSIZE=0.1D0
    INTEGER::ADJELT(4)=0  
ENDTYPE
TYPE(TET_TYDEF),ALLOCATABLE::TET(:)
INTEGER::NTET=0

INTEGER,ALLOCATABLE::EDGE_BC(:),FACE_BC(:),NODE_LOOP_BC(:)
INTEGER::NE_BC=0,NF_BC=0,NNLBC=0






type et_type
	integer::nnode=0,nedge=0,nface=0,ntet=0		!NOTE THAT THE NNODE IS ONLY THE END NODE NUMBER 
	character(512)::description
	integer,allocatable::edge(:,:),face(:,:),FaceEdge(:,:)
	integer,allocatable::tet(:,:) !,tetEDGE(:,:),TETFACE(:,:)
	INTEGER::DIM=-1
	!edge(2,nedge),
	!face: use node index to represent face. face(0:4,nface),face(0,:)==3,triangular face,==4, quadrilateral face
	!FaceEdge: use edge index to represent face. FaceEdge(0:4,nface),FaceEdge(0,:)==3,triangular face, ==4, quadrilateral face
	real(8),allocatable::weight(:,:) !�ֲ���Ԫ���ظ��ڵ�ֲ�ϵ��, weight(:,1) ƽ��������أ�weight(:,2) ƽ�������κ���;weight(:,3) ��Գƾ������أ�weight(:,4) ��Գ������κ��أ�
												!Ŀǰֻ�ܴ���ƽ��Ӧ������������
end type
type(et_type)::elttype(100)

REAL(8)::VTOL=1.D-7

INTEGER::LASTLOC=1
INTEGER,ALLOCATABLE::ISSEARCH(:)

TYPE SEARCHZONE_TYDEF
    INTEGER::NEL=0,NDX(3)=1 !NDX=DIVISION IN X Y AND Y AXIS
    REAL(8)::BBOX(2,3) !MIN,MAX
    INTEGER,ALLOCATABLE::ELEMENT(:)
ENDTYPE
TYPE(SEARCHZONE_TYDEF),ALLOCATABLE::SEARCHZONE(:)
INTEGER::NSZONE=1

CONTAINS 


SUBROUTINE SETUP_SUBZONE_TET()
    IMPLICIT NONE
    REAL(8)::DX,DY,DZ
    INTEGER::NDX,NDY,NDZ,I,J,K,N1
    
    NDX=5;NDY=5;NDZ=5
    DX=(POSDATA.MAXX-POSDATA.MINX)/NDX
    IF(ABS(DX)<1.D-7) NDX=1
    DY=(POSDATA.MAXY-POSDATA.MINY)/NDY
    IF(ABS(DY)<1.D-7) NDY=1
    DZ=(POSDATA.MAXZ-POSDATA.MINZ)/NDZ
    IF(ABS(DZ)<1.D-7) NDZ=1
    
    NSZONE=NDX*NDY*NDZ
    
    ALLOCATE(SEARCHZONE(NSZONE))
    
    

    N1=0
    DO I=1,NDZ
        DO J=1,NDY
            DO K=1,NDX
                N1=N1+1
                SEARCHZONE(N1).BBOX(1,1)=POSDATA.MINX+DX*(K-1)
                SEARCHZONE(N1).BBOX(2,1)=POSDATA.MINX+DX*K
                SEARCHZONE(N1).BBOX(1,2)=POSDATA.MINY+DY*(J-1)
                SEARCHZONE(N1).BBOX(2,2)=POSDATA.MINY+DY*J 
                SEARCHZONE(N1).BBOX(1,3)=POSDATA.MINZ+DZ*(I-1)
                SEARCHZONE(N1).BBOX(2,3)=POSDATA.MINZ+DZ*I
                IF(.NOT.ALLOCATED(SEARCHZONE(N1).ELEMENT)) ALLOCATE(SEARCHZONE(N1).ELEMENT(NTET))
                SEARCHZONE(N1).ELEMENT=0
                SEARCHZONE(N1).NDX=[NDX,NDY,NDZ]
            ENDDO
        ENDDO
    ENDDO
    
    DO I=1,NTET
        DO J=1,NSZONE
            
            IF(NDX>1) THEN
                IF(SEARCHZONE(J).BBOX(1,1)>TET(I).BBOX(2,1)) CYCLE !MIN>MAX
                IF(SEARCHZONE(J).BBOX(2,1)<TET(I).BBOX(1,1)) CYCLE !MAX<MIN
            ENDIF
            IF(NDY>1) THEN
                IF(SEARCHZONE(J).BBOX(1,2)>TET(I).BBOX(2,2)) CYCLE !MIN>MAX
                IF(SEARCHZONE(J).BBOX(2,2)<TET(I).BBOX(1,2)) CYCLE !MAX<MIN 
            ENDIF
            IF(NDZ>1) THEN
                IF(SEARCHZONE(J).BBOX(1,3)>TET(I).BBOX(2,3)) CYCLE !MIN>MAX
                IF(SEARCHZONE(J).BBOX(2,3)<TET(I).BBOX(1,3)) CYCLE !MAX<MIN 
            ENDIF
            SEARCHZONE(J).NEL=SEARCHZONE(J).NEL+1
            SEARCHZONE(J).ELEMENT(SEARCHZONE(J).NEL)=I
        ENDDO
    ENDDO
    
    
    

ENDSUBROUTINE

SUBROUTINE SETUP_ADJACENT_ELEMENT_TET()

    IMPLICIT NONE
    INTEGER::I,J
    
    !CAN'T NOT HANDLE 2D AND 3D COMBINATION MODEL, TO BE IMPROVED.
    
    IF(POSDATA.NDIM==2) THEN
        DO I=1,NEDGE
            IF(EDGE(I).ENUM==2) THEN
                TET(EDGE(I).ELEMENT(1)).ADJELT(EDGE(I).SUBID(1))=EDGE(I).ELEMENT(2)
                TET(EDGE(I).ELEMENT(2)).ADJELT(EDGE(I).SUBID(2))=EDGE(I).ELEMENT(1)
            ELSEIF(EDGE(I).ENUM>2) THEN
                PRINT *, 'MORE THEN 2 ELEMENTS SHARING ONE COMMON EDGE. EDGE=',I
            ENDIF    
        ENDDO 
        
    ELSEIF(POSDATA.NDIM>2) THEN
        DO I=1,NFACE
            IF(FACE(I).ENUM==2) THEN            
                TET(ABS(FACE(I).ELEMENT(1))).ADJELT(FACE(I).SUBID(1))=ABS(FACE(I).ELEMENT(2))
                TET(ABS(FACE(I).ELEMENT(2))).ADJELT(FACE(I).SUBID(2))=ABS(FACE(I).ELEMENT(1))
            ELSEIF(FACE(I).ENUM>2) THEN
                PRINT *, 'MORE THEN 2 ELEMENTS SHARING ONE COMMON FACE. FACE=',I
            ENDIF    
        ENDDO
    ENDIF
 
    
ENDSUBROUTINE



SUBROUTINE SETUP_SUB_TET4_ELEMENT(ELEMENT,NEL,ESET,NESET,NODE,NNODE)
	IMPLICIT NONE
    
	INTEGER,INTENT(IN)::NEL,NESET,NNODE
	TYPE(ELEMENT_TYDEF)::ELEMENT(NEL)    
    TYPE(ESET_TYDEF)::ESET(NESET)
    TYPE(NODE_TYDEF)::NODE(NNODE)
    
    INTEGER::i,J,K,GMET1
    !IF(.NOT.ISINI_GMSHET) THEN
    !    IF(.NOT.ALLOCATED(EDGE)) ALLOCATE(EDGE(MAXNEDGE),FACE(MAXNFACE))
    !    CALL Initialize_et2numNodes()
    !    CALL ET_GMSH_EDGE_FACE()
    !    ISINI_GMSHET=.TRUE.
    !ENDIF
	
	NTET=0
	DO I=1,NEL
        IF(ESET(ELEMENT(I).ISET).COUPLESET>0.AND.ESET(ELEMENT(I).ISET).COUPLESET<ELEMENT(I).ISET) CYCLE !Ӱ�ӵ�Ԫ������
        
		GMET1=GETGMSHET(ESET(ELEMENT(I).ISET).ET)
		SELECT CASE(ELTTYPE(GMET1).DIM)
		CASE(0,1)
			NTET=NTET+1	
			ELEMENT(I).NTET=1
		CASE(2,3)
			NTET=NTET+Elttype(GMET1).NTET
			ELEMENT(I).NTET=Elttype(GMET1).NTET
		CASE DEFAULT
			STOP "ERROR IN ELEMENT DIMENSION. SUB=SETUP_SUB_TET4_ELEMENT"
		END SELECT
		ALLOCATE(ELEMENT(I).TET(ELEMENT(I).NTET))
	ENDDO
	
	IF(ALLOCATED(TET)) DEALLOCATE(TET)
	ALLOCATE(TET(NTET))
	NTET=0
	DO I=1,NEL
        IF(ESET(ELEMENT(I).ISET).COUPLESET>0.and.(ESET(ELEMENT(I).ISET).COUPLESET<ELEMENT(I).ISET)) CYCLE !Ӱ�ӵ�Ԫ������
        
		GMET1=GETGMSHET(ESET(ELEMENT(I).ISET).ET)
		ELEMENT(I).TET=[NTET+1:NTET+ELEMENT(I).NTET]
		TET(ELEMENT(I).TET).ISET=ELEMENT(I).ISET
		SELECT CASE(ELTTYPE(GMET1).DIM)		
		CASE(0,1)
			NTET=NTET+1
			TET(NTET).MOTHER=I
			TET(NTET).NV=ELEMENT(I).NNUM
			TET(NTET).V(1:TET(NTET).NV)=ELEMENT(I).NODE
			IF(ELTTYPE(GMET1).DIM==0) THEN
				TET(NTET).GMET=15
				TET(NTET).DIM=0
			ELSE
				TET(NTET).GMET=1
				TET(NTET).DIM=1
			ENDIF
			DO K=1,3
				TET(NTET).BBOX(1,K)=MINVAL(NODE(TET(NTET).V(1:TET(NTET).NV)).COORD(K))
				TET(NTET).BBOX(2,K)=MAXVAL(NODE(TET(NTET).V(1:TET(NTET).NV)).COORD(K))
			ENDDO
		CASE(2)
			DO J=1,Elttype(GMET1).NTET				
				NTET=NTET+1
				TET(NTET).MOTHER=I
				TET(NTET).NV=3;TET(NTET).GMET=2;TET(NTET).DIM=2
				TET(NTET).V(1:3)=ELEMENT(I).NODE(Elttype(GMET1).TET(1:3,J))
				DO K=1,3
					TET(NTET).BBOX(1,K)=MINVAL(NODE(TET(NTET).V(1:TET(NTET).NV)).COORD(K))
					TET(NTET).BBOX(2,K)=MAXVAL(NODE(TET(NTET).V(1:TET(NTET).NV)).COORD(K))
				ENDDO
			ENDDO
		CASE(3)
			DO J=1,Elttype(GMET1).NTET				
				NTET=NTET+1
				TET(NTET).MOTHER=I
				TET(NTET).NV=4;TET(NTET).GMET=4;TET(NTET).DIM=3
				TET(NTET).V(1:4)=ELEMENT(I).NODE(Elttype(GMET1).TET(:,J))
				DO K=1,3
					TET(NTET).BBOX(1,K)=MINVAL(NODE(TET(NTET).V(1:TET(NTET).NV)).COORD(K))
					TET(NTET).BBOX(2,K)=MAXVAL(NODE(TET(NTET).V(1:TET(NTET).NV)).COORD(K))
				ENDDO
			ENDDO
		END SELECT
		
	ENDDO	
	
END SUBROUTINE

    INTEGER FUNCTION GETGMSHET(ET)
        USE solverds
        IMPLICIT NONE
        INTEGER,INTENT(IN)::ET
    
        SELECT CASE(ET)
        !LINE,2-NODE LINE
        CASE(BAR,BEAM,BAR2D,BEAM2D,FELINESEG,PIPE2,SEMI_SPHFLOW,SPHFLOW)    
            GETGMSHET=1
        CASE(FETRIANGLE,CPE3,CPE3_SPG,CPE3_CPL,CPS3,SHELL3,DKT3,&
             CAX3,CAX3_SPG,CAX3_CPL)
             GETGMSHET=2
        CASE(FEQUADRILATERAL,CPE4,CPE4_SPG,CPE4_CPL,CPS4,&
             CAX4,CAX4_SPG,CAX4_CPL,&
             CPE4R,CPE4R_SPG,CPE4R_CPL,CPS4R,&
             CAX4R,CAX4R_SPG,CAX4R_CPL,ZT4_SPG,WELLBORE)
             GETGMSHET=3 
        CASE(FETETRAHEDRON,TET4,TET4_SPG,TET4_CPL)         
            GETGMSHET=4
        CASE(FEBRICK)
            GETGMSHET=5
        CASE(PRM6,PRM6_SPG,PRM6_CPL,ZT6_SPG)
            GETGMSHET=6
        CASE(CPE6,CPE6_SPG,CPE6_CPL,CPS6,&
             CAX6,CAX6_SPG,CAX6_CPL)
            GETGMSHET=9       
        CASE(TET10,TET10_SPG,TET10_CPL)         
            GETGMSHET=11
        CASE(SPRINGX,SPRINGY,SPRINGZ,&
             SPRINGMX,SPRINGMY,SPRINGMZ,&
             SOILSPRINGX,SOILSPRINGY,SOILSPRINGZ)
             GETGMSHET=15
        CASE(CPE8,CPE8_SPG,CPE8_CPL,CPS8,&
             CAX8,CAX8_SPG,CAX8_CPL,&
             CPE8R,CPE8R_SPG,CPE8R_CPL,CPS8R,&
             CAX8R,CAX8R_SPG,CAX8R_CPL)
             GETGMSHET=16
        CASE(CPE15,CPE15_SPG,CPE15_CPL,CPS15)
             GETGMSHET=23
        CASE(PRM15,PRM15_SPG,PRM15_CPL)
             GETGMSHET=18        
        CASE DEFAULT
            PRINT *, 'NO SUCH ELEMENT TYPE. FUNC=GETGMSHET.'
            PAUSE
        ENDSELECT
        
            
    
    ENDFUNCTION  

    
!SUBROUTINE SETUP_EDGE_TBL_TET()
!    
!    IMPLICIT NONE
!    INTEGER::I,J,N1(2),ET1,NEDGE1,TBL_LEN
!	CHARACTER(LEN=:),ALLOCATABLE::KEY1
!    TYPE(DICT_DATA)::VAL1
!    CHARACTER(64)::CKEY1
!    INTEGER::HKEY1
!    
!    !IF(.NOT.ISINI_GMSHET) THEN
!    !    
!    !    CALL Initialize_et2numNodes()
!    !    CALL ET_GMSH_EDGE_FACE()
!    !    ISINI_GMSHET=.TRUE.
!    !ENDIF
!	
!    IF(.NOT.ALLOCATED(EDGE)) ALLOCATE(EDGE(MAXNEDGE))
!	
!	TBL_LEN=2*POSDATA.NNODE	
!    CALL EDGE_TBL.INIT(TBL_LEN)
!	DO I=1,NTET
!		ET1=TET(I).GMET
!		NEDGE1=ELTTYPE(ET1).NEDGE
!		TET(I).NE=NEDGE1
!        !ALLOCATE(TET(I).E(NEDGE1))
!		DO J=1,NEDGE1
!            N1=TET(I).V(ELTTYPE(ET1).EDGE(:,J))
!            CALL I2C_KEY_hash_tbl_sll(KEY1,N1(:),2)
!            CKEY1=TRIM(KEY1)
!            HKEY1 = MOD(ABS(HASH_DJB(KEY1)),EDGE_TBL%vec_len)
!            VAL1.IEL=I;VAL1.ISE=J;VAL1.IITEM=EDGE_TBL.TBL_ID
!			CALL EDGE_TBL.PUT(TRIM(KEY1),VAL1,HKEY1)
!			TET(I).E(J)=VAL1.IITEM
!            IF(VAL1.IITEM>MAXNEDGE) THEN
!                CALL EDGE_TYDEF_ENLARGE_AR(EDGE,1000)
!                MAXNEDGE=MAXNEDGE+1000
!            ENDIF
!            IF(.NOT.EDGE(VAL1.IITEM).ISINI) THEN
!                EDGE(VAL1.IITEM).CKEY=CKEY1
!                EDGE(VAL1.IITEM).HKEY=HKEY1
!                EDGE(VAL1.IITEM).V=TET(I).V(ELTTYPE(ET1).EDGE(:,J))
!                EDGE(VAL1.IITEM).DIS=NORM2(POSDATA.NODE(EDGE(VAL1.IITEM).V(1)).COORD-POSDATA.NODE(EDGE(VAL1.IITEM).V(2)).COORD)
!                EDGE(VAL1.IITEM).ISINI=.TRUE.
!                NEDGE=VAL1.IITEM
!            ENDIF
!            IF(EDGE(VAL1.IITEM).ENUM==0) ALLOCATE(EDGE(VAL1.IITEM).ELEMENT(5),EDGE(VAL1.IITEM).SUBID(5))
!            EDGE(VAL1.IITEM).ENUM=EDGE(VAL1.IITEM).ENUM+1
!            IF(EDGE(VAL1.IITEM).ENUM>SIZE(EDGE(VAL1.IITEM).ELEMENT,DIM=1)) THEN
!				CALL I_ENLARGE_AR(EDGE(VAL1.IITEM).ELEMENT,5)
!				CALL I_ENLARGE_AR(EDGE(VAL1.IITEM).SUBID,5)
!			ENDIF
!            EDGE(VAL1.IITEM).ELEMENT(EDGE(VAL1.IITEM).ENUM)=I			
!            EDGE(VAL1.IITEM).SUBID(EDGE(VAL1.IITEM).ENUM)=J
!		ENDDO
!    END DO
!    RETURN
!ENDSUBROUTINE
!
!SUBROUTINE SETUP_FACE_TBL_TET()
!    IMPLICIT NONE
!   
!    INTEGER::I,J,N1(4),ET1,NFACE1,TBL_LEN,N2,K
!    REAL(8)::V1(3),V2(3),NORMAL1(3),T1
!	CHARACTER(LEN=:),ALLOCATABLE::KEY1
!    TYPE(DICT_DATA)::VAL1
!    CHARACTER(64)::CKEY1
!    INTEGER::HKEY1
!    
!    !IF(.NOT.ISINI_GMSHET) THEN
!    !    
!    !    CALL Initialize_et2numNodes()
!    !    CALL ET_GMSH_EDGE_FACE()
!    !    ISINI_GMSHET=.TRUE.
!    !ENDIF 	
!
!	IF(.NOT.ALLOCATED(FACE)) ALLOCATE(FACE(MAXNFACE))
!	
!	TBL_LEN=2*POSDATA.NNODE	
!    CALL FACE_TBL.INIT(TBL_LEN)
!	DO I=1,NTET
!		ET1=TET(I).GMET
!		NFACE1=ELTTYPE(ET1).NFACE
!		TET(I).NF=NFACE1
!        !ALLOCATE(TET(I).F(NFACE1))
!		DO J=1,NFACE1
!            N2=ELTTYPE(ET1).FACE(0,J)
!            N1(1:N2)=TET(I).V(ELTTYPE(ET1).FACE(1:N2,J))
!            CALL FACE_TBL.KEY(KEY1,N1(1:N2),ELTTYPE(ET1).FACE(0,J))
!            CKEY1=TRIM(KEY1)
!            HKEY1 = MOD(ABS(HASH_DJB(KEY1)),FACE_TBL%vec_len)            
!            VAL1.IEL=I;VAL1.ISE=J;VAL1.IITEM=FACE_TBL.TBL_ID
!			CALL FACE_TBL.PUT(TRIM(KEY1),VAL1,HKEY1)
!            TET(I).F(J)=VAL1.IITEM
!            
!            IF(VAL1.IITEM>MAXNFACE) THEN
!                CALL FACE_TYDEF_ENLARGE_AR(FACE,1000)
!                MAXNFACE=MAXNFACE+1000
!            ENDIF
!            IF(.NOT.FACE(VAL1.IITEM).ISINI) THEN
!                FACE(VAL1.IITEM).CKEY=CKEY1 
!                FACE(VAL1.IITEM).HKEY=HKEY1
!                FACE(VAL1.IITEM).SHAPE=ELTTYPE(ET1).FACE(0,J)
!                FACE(VAL1.IITEM).V(1:FACE(VAL1.IITEM).SHAPE)=TET(I).V(ELTTYPE(ET1).FACE(1:FACE(VAL1.IITEM).SHAPE,J))
!				FACE(VAL1.IITEM).EDGE(1:FACE(VAL1.IITEM).SHAPE)=TET(I).E(ABS(ELTTYPE(ET1).FACEEDGE(1:FACE(VAL1.IITEM).SHAPE,J)))
!				!CHECK LOOP ORDER
!				DO K=1,FACE(VAL1.IITEM).SHAPE
!					IF(EDGE(FACE(VAL1.IITEM).EDGE(K)).V(1)/=FACE(VAL1.IITEM).V(K)) FACE(VAL1.IITEM).EDGE(K)=-FACE(VAL1.IITEM).EDGE(K)
!				ENDDO
!                V1=POSDATA.NODE(FACE(VAL1.IITEM).V(2)).COORD-POSDATA.NODE(FACE(VAL1.IITEM).V(1)).COORD
!                V2=POSDATA.NODE(FACE(VAL1.IITEM).V(3)).COORD-POSDATA.NODE(FACE(VAL1.IITEM).V(1)).COORD
!                call r8vec_cross_3d ( v1, v2, NORMAL1(1:3) )
!                T1 = sqrt ( sum ( ( NORMAL1 )**2 ) )
!                if ( T1 /= 0.0D+00 ) then
!                  NORMAL1 = NORMAL1/T1
!                end if
!                FACE(VAL1.IITEM).UNORMAL=NORMAL1
!                FACE(VAL1.IITEM).ISINI=.TRUE.
!                NFACE=VAL1.IITEM
!				DO K=1,3
!					FACE(NFACE).BBOX(1,K)=MINVAL(POSDATA.NODE(FACE(NFACE).V(1:FACE(NFACE).SHAPE)).COORD(K))-VTOL
!					FACE(NFACE).BBOX(2,K)=MAXVAL(POSDATA.NODE(FACE(NFACE).V(1:FACE(NFACE).SHAPE)).COORD(K))+VTOL
!				ENDDO
!            ENDIF
!            IF(FACE(VAL1.IITEM).ENUM==0) ALLOCATE(FACE(VAL1.IITEM).ELEMENT(2),FACE(VAL1.IITEM).SUBID(2))
!            FACE(VAL1.IITEM).ENUM=FACE(VAL1.IITEM).ENUM+1
!            IF(FACE(VAL1.IITEM).ENUM>SIZE(FACE(VAL1.IITEM).ELEMENT,DIM=1)) THEN
!                CALL I_ENLARGE_AR(FACE(VAL1.IITEM).ELEMENT,5)
!                CALL I_ENLARGE_AR(FACE(VAL1.IITEM).SUBID,5)
!            ENDIF
!            
!            IF(FACE(VAL1.IITEM).ENUM>1) THEN
!                DO K=1,FACE(VAL1.IITEM).SHAPE
!                    IF(FACE(VAL1.IITEM).V(K)==TET(I).V(ELTTYPE(TET(I).GMET).FACE(1,J))) THEN
!                        IF(FACE(VAL1.IITEM).V(MOD(K,FACE(VAL1.IITEM).SHAPE)+1)/= &
!                            TET(I).V(ELTTYPE(TET(I).GMET).FACE(2,J))) THEN
!                            FACE(VAL1.IITEM).ELEMENT(FACE(VAL1.IITEM).ENUM)=-I !I��Ԫ����ķ�����face�ķ�����,Ϊ-1
!                        ELSE
!                            FACE(VAL1.IITEM).ELEMENT(FACE(VAL1.IITEM).ENUM)=I
!                        ENDIF                         
!                        EXIT
!                    ENDIF
!                ENDDO
!            ELSE
!                FACE(VAL1.IITEM).ELEMENT(FACE(VAL1.IITEM).ENUM)=I !ͬ�� =1
!            ENDIF
!            FACE(VAL1.IITEM).SUBID(FACE(VAL1.IITEM).ENUM)=J
!                
!            
!		ENDDO
!    END DO    
!    
!    
!ENDSUBROUTINE
!
!SUBROUTINE SETUP_EDGE_TBL(ELEMENT,NEL,ESET,NESET,NODE,NNODE)
!    
!    IMPLICIT NONE
!	INTEGER,INTENT(IN)::NEL,NESET,NNODE
!	TYPE(ELEMENT_TYDEF)::ELEMENT(NEL)    
!    TYPE(ESET_TYDEF)::ESET(NESET)
!    TYPE(NODE_TYDEF)::NODE(NNODE)
!    
!    INTEGER::I,J,N1(2),ET1,NEDGE1,TBL_LEN
!	CHARACTER(LEN=:),ALLOCATABLE::KEY1
!    TYPE(DICT_DATA)::VAL1
!    CHARACTER(64)::CKEY1
!    INTEGER::HKEY1
!    
!    !IF(.NOT.ISINI_GMSHET) THEN        
!    !    CALL Initialize_et2numNodes()
!    !    CALL ET_GMSH_EDGE_FACE()
!    !    ISINI_GMSHET=.TRUE.
!    !ENDIF
!	
!    IF(.NOT.ALLOCATED(MEDGE)) ALLOCATE(MEDGE(MAXNMEDGE))
!	
!	TBL_LEN=NNODE	
!    CALL MEDGE_TBL.INIT(TBL_LEN)
!	DO I=1,NEL
!		ET1=GETGMSHET(ESET(ELEMENT(I).ISET).ET)
!		NEDGE1=ELTTYPE(ET1).NEDGE
!		ELEMENT(I).NEDGE=NEDGE1
!        ALLOCATE(ELEMENT(I).EDGE(NEDGE1))
!		DO J=1,NEDGE1
!            N1=ELEMENT(I).NODE(ELTTYPE(ET1).EDGE(:,J))
!            CALL I2C_KEY_hash_tbl_sll(KEY1,N1(:),2)
!            CKEY1=TRIM(KEY1)
!            HKEY1 = MOD(ABS(HASH_DJB(KEY1)),MEDGE_TBL%vec_len)
!            VAL1.IEL=I;VAL1.ISE=J;VAL1.IITEM=MEDGE_TBL.TBL_ID
!			CALL MEDGE_TBL.PUT(TRIM(KEY1),VAL1,HKEY1)
!			ELEMENT(I).EDGE(J)=VAL1.IITEM
!            IF(VAL1.IITEM>MAXNMEDGE) THEN
!                CALL EDGE_TYDEF_ENLARGE_AR(MEDGE,1000)
!                MAXNMEDGE=MAXNMEDGE+1000
!            ENDIF
!            IF(.NOT.MEDGE(VAL1.IITEM).ISINI) THEN
!                MEDGE(VAL1.IITEM).CKEY=CKEY1
!                MEDGE(VAL1.IITEM).HKEY=HKEY1
!                MEDGE(VAL1.IITEM).V=ELEMENT(I).NODE(ELTTYPE(ET1).EDGE(:,J))
!                MEDGE(VAL1.IITEM).DIS=NORM2(NODE(MEDGE(VAL1.IITEM).V(1)).COORD-NODE(MEDGE(VAL1.IITEM).V(2)).COORD)
!                MEDGE(VAL1.IITEM).ISINI=.TRUE.
!                NMEDGE=VAL1.IITEM
!            ENDIF
!            IF(MEDGE(VAL1.IITEM).ENUM==0) ALLOCATE(MEDGE(VAL1.IITEM).ELEMENT(5),MEDGE(VAL1.IITEM).SUBID(5))
!            MEDGE(VAL1.IITEM).ENUM=MEDGE(VAL1.IITEM).ENUM+1
!            IF(MEDGE(VAL1.IITEM).ENUM>SIZE(MEDGE(VAL1.IITEM).ELEMENT,DIM=1)) THEN
!				CALL I_ENLARGE_AR(MEDGE(VAL1.IITEM).ELEMENT,5)
!				CALL I_ENLARGE_AR(MEDGE(VAL1.IITEM).SUBID,5)
!			ENDIF
!            MEDGE(VAL1.IITEM).ELEMENT(MEDGE(VAL1.IITEM).ENUM)=I			
!            MEDGE(VAL1.IITEM).SUBID(MEDGE(VAL1.IITEM).ENUM)=J
!		ENDDO
!    END DO
!    RETURN
!ENDSUBROUTINE
!
!SUBROUTINE SETUP_FACE_TBL(ELEMENT,NEL,ESET,NESET,NODE,NNODE)
!    IMPLICIT NONE
!    INTEGER,INTENT(IN)::NEL,NESET,NNODE
!	TYPE(ELEMENT_TYDEF)::ELEMENT(NEL)    
!    TYPE(ESET_TYDEF)::ESET(NESET)
!    TYPE(NODE_TYDEF)::NODE(NNODE)
!    INTEGER::I,J,N1(4),ET1,NFACE1,TBL_LEN,N2,K
!    REAL(8)::V1(3),V2(3),NORMAL1(3),T1
!	CHARACTER(LEN=:),ALLOCATABLE::KEY1
!    TYPE(DICT_DATA)::VAL1
!    CHARACTER(64)::CKEY1
!    INTEGER::HKEY1
!    
!    !IF(.NOT.ISINI_GMSHET) THEN
!    !    CALL Initialize_et2numNodes()
!    !    CALL ET_GMSH_EDGE_FACE()
!    !    ISINI_GMSHET=.TRUE.
!    !ENDIF    
!	
!    IF(.NOT.ALLOCATED(MFACE)) ALLOCATE(MFACE(MAXNMFACE))
!	
!	TBL_LEN=NNODE	
!    CALL MFACE_TBL.INIT(TBL_LEN)
!	DO I=1,NEL
!		ET1=GETGMSHET(ESET(ELEMENT(I).ISET).ET)
!		NFACE1=ELTTYPE(ET1).NFACE
!		ELEMENT(I).NFACE=NFACE1
!        ALLOCATE(ELEMENT(I).FACE(NFACE1))
!		DO J=1,NFACE1
!            N2=ELTTYPE(ET1).FACE(0,J)
!            N1(1:N2)=ELEMENT(I).NODE(ELTTYPE(ET1).FACE(1:N2,J))
!            CALL MFACE_TBL.KEY(KEY1,N1(1:N2),ELTTYPE(ET1).FACE(0,J))
!            CKEY1=TRIM(KEY1)
!            HKEY1 = MOD(ABS(HASH_DJB(KEY1)),MFACE_TBL%vec_len)            
!            VAL1.IEL=I;VAL1.ISE=J;VAL1.IITEM=MFACE_TBL.TBL_ID
!			CALL MFACE_TBL.PUT(TRIM(KEY1),VAL1,HKEY1)
!            ELEMENT(I).FACE(J)=VAL1.IITEM
!            
!            IF(VAL1.IITEM>MAXNMFACE) THEN
!                CALL FACE_TYDEF_ENLARGE_AR(MFACE,1000)
!                MAXNMFACE=MAXNMFACE+1000
!            ENDIF
!            IF(.NOT.MFACE(VAL1.IITEM).ISINI) THEN
!                MFACE(VAL1.IITEM).CKEY=CKEY1 
!                MFACE(VAL1.IITEM).HKEY=HKEY1
!                MFACE(VAL1.IITEM).SHAPE=ELTTYPE(ET1).FACE(0,J)
!                MFACE(VAL1.IITEM).V(1:MFACE(VAL1.IITEM).SHAPE)=ELEMENT(I).NODE(ELTTYPE(ET1).FACE(1:MFACE(VAL1.IITEM).SHAPE,J))
!				MFACE(VAL1.IITEM).EDGE(1:MFACE(VAL1.IITEM).SHAPE)=ELEMENT(I).EDGE(ABS(ELTTYPE(ET1).FACEEDGE(1:MFACE(VAL1.IITEM).SHAPE,J)))
!				!CHECK LOOP ORDER
!				DO K=1,MFACE(VAL1.IITEM).SHAPE
!					IF(MEDGE(MFACE(VAL1.IITEM).EDGE(K)).V(1)/=MFACE(VAL1.IITEM).V(K)) MFACE(VAL1.IITEM).EDGE(K)=-MFACE(VAL1.IITEM).EDGE(K)
!				ENDDO
!                V1=NODE(MFACE(VAL1.IITEM).V(2)).COORD-NODE(MFACE(VAL1.IITEM).V(1)).COORD
!                V2=NODE(MFACE(VAL1.IITEM).V(3)).COORD-NODE(MFACE(VAL1.IITEM).V(1)).COORD
!                call r8vec_cross_3d ( v1, v2, NORMAL1(1:3) )
!                T1 = sqrt ( sum ( ( NORMAL1 )**2 ) )
!                if ( T1 /= 0.0D+00 ) then
!                  NORMAL1 = NORMAL1/T1
!                end if
!                MFACE(VAL1.IITEM).UNORMAL=NORMAL1
!                MFACE(VAL1.IITEM).ISINI=.TRUE.
!                NMFACE=VAL1.IITEM
!            ENDIF
!            IF(MFACE(VAL1.IITEM).ENUM==0) ALLOCATE(MFACE(VAL1.IITEM).ELEMENT(2),MFACE(VAL1.IITEM).SUBID(2))
!            MFACE(VAL1.IITEM).ENUM=MFACE(VAL1.IITEM).ENUM+1
!            IF(MFACE(VAL1.IITEM).ENUM>SIZE(MFACE(VAL1.IITEM).ELEMENT,DIM=1)) THEN
!                CALL I_ENLARGE_AR(MFACE(VAL1.IITEM).ELEMENT,5)
!                CALL I_ENLARGE_AR(MFACE(VAL1.IITEM).SUBID,5)
!            ENDIF
!            IF(MFACE(VAL1.IITEM).ENUM>1) THEN
!                DO K=1,MFACE(VAL1.IITEM).SHAPE
!                    IF(MFACE(VAL1.IITEM).V(K)==ELEMENT(I).NODE(ELTTYPE(GETGMSHET(ESET(ELEMENT(I).ISET).ET)).FACE(1,J))) THEN
!                        IF(MFACE(VAL1.IITEM).V(MOD(K,MFACE(VAL1.IITEM).SHAPE)+1)/= &
!                            ELEMENT(I).NODE(ELTTYPE(GETGMSHET(ESET(ELEMENT(I).ISET).ET)).FACE(2,J))) THEN
!                            MFACE(VAL1.IITEM).ELEMENT(MFACE(VAL1.IITEM).ENUM)=-I !I?????????????face???????,?-1
!                        ELSE
!                            MFACE(VAL1.IITEM).ELEMENT(MFACE(VAL1.IITEM).ENUM)=I
!                        ENDIF                  
!                        EXIT
!                    ENDIF
!                ENDDO
!            ELSE
!                MFACE(VAL1.IITEM).ELEMENT(MFACE(VAL1.IITEM).ENUM)=I !??? =1
!            ENDIF
!            MFACE(VAL1.IITEM).SUBID(MFACE(VAL1.IITEM).ENUM)=J
!                
!            
!		ENDDO
!    END DO    
!    
!    
!ENDSUBROUTINE


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

SUBROUTINE ET_GMSH_EDGE_FACE()

	IMPLICIT NONE
	INTEGER::I,ET,AET1(33)
	
    AET1=[1:31,92,93]
	DO I=1,33
        ET=AET1(I)
		SELECT CASE(ET)
			CASE(1) !LINE
                ELTTYPE(ET).NNODE=2
				Elttype(ET).NEDGE=1;Elttype(ET).NFACE=0;Elttype(ET).NTET=0;ELTTYPE(ET).DIM=1
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,1)=[1,2]
			CASE(2) !TRIANGLE
                ELTTYPE(ET).NNODE=3
				Elttype(ET).NEDGE=3;Elttype(ET).NFACE=1;Elttype(ET).NTET=1;ELTTYPE(ET).DIM=2
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE),Elttype(ET).TET(4,Elttype(ET).NTET))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,1],(/2,3/))
				Elttype(ET).FACE(:,1)=[3,1,2,3,0]
				Elttype(ET).FACEEDGE=Elttype(ET).FACE
				Elttype(ET).TET(:,1)=[1,2,3,0]
			CASE(9) !6-NODED-TRIANGLE
                ELTTYPE(ET).NNODE=3
				Elttype(ET).NEDGE=3;Elttype(ET).NFACE=1;Elttype(ET).NTET=4;ELTTYPE(ET).DIM=2
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE),Elttype(ET).TET(4,Elttype(ET).NTET))
				
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,1],(/2,3/))
				Elttype(ET).FACE(:,1)=[3,1,2,3,0]
				Elttype(ET).FACEEDGE=Elttype(ET).FACE
				Elttype(ET).TET(:,:)=RESHAPE([1,4,6,0,&
                                               4,2,5,0,&
                                               5,3,6,0,&
                                               4,5,6,0],&
                                               (/4,4/))				
			CASE(23) !15-NODED-TRIANGLE
                ELTTYPE(ET).NNODE=3
				Elttype(ET).NEDGE=3;Elttype(ET).NFACE=1;Elttype(ET).NTET=16;ELTTYPE(ET).DIM=2
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE),Elttype(ET).TET(4,Elttype(ET).NTET))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,1],(/2,3/))
				Elttype(ET).FACE(:,1)=[3,1,2,3,0]
				Elttype(ET).FACEEDGE=Elttype(ET).FACE
				Elttype(ET).TET(:,:)=RESHAPE([ 1,7,12,0,&
											   7,4,15,0,&
											   15,6,12,0,&
											   7,15,12,0,&
											   
											   4,8,13,0,&
											   8,2,9,0,&
											   9,5,13,0,&
											   8,9,13,0,&
											   
											   5,10,14,0,&
											   10,3,11,0,&
											   11,6,14,0,&
											   10,11,14,0,&
											   
											   4,13,15,0,&
											   13,5,14,0,&
											   14,6,15,0,&
											   13,14,15,0],&
											   (/4,16/))
				
				
			CASE(3) !QUADRANGLE
                ELTTYPE(ET).NNODE=4
				Elttype(ET).NEDGE=4;Elttype(ET).NFACE=1;Elttype(ET).NTET=2;ELTTYPE(ET).DIM=2
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE),Elttype(ET).TET(4,Elttype(ET).NTET))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,4,4,1],(/2,4/))                
				Elttype(ET).FACE(:,:)=RESHAPE([4,1,2,3,4],&
											  (/5,1/))
				Elttype(ET).FACEEDGE=Elttype(ET).FACE
				Elttype(ET).TET(:,:)=RESHAPE([1,2,3,0,&
                                              3,4,1,0],(/4,2/))				
                
			CASE(16) !8-noded-QUADRANGLE
                ELTTYPE(ET).NNODE=4
				Elttype(ET).NEDGE=4;Elttype(ET).NFACE=1;Elttype(ET).NTET=6;ELTTYPE(ET).DIM=2				
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE),Elttype(ET).TET(4,Elttype(ET).NTET))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,4,4,1],(/2,4/))                
				Elttype(ET).FACE(:,:)=RESHAPE([4,1,2,3,4],&
											  (/5,1/))
				Elttype(ET).FACEEDGE=Elttype(ET).FACE
				Elttype(ET).TET(:,:)=RESHAPE([1,5,8,0,&
											   5,2,6,0,&
											   6,3,7,0,&
											   7,4,8,0,&
											   5,6,8,0,&
											   6,7,8,0],&
											   (/4,6/))				
				
			CASE(4) !TETRAHEDRON
                ELTTYPE(ET).NNODE=4
				Elttype(ET).NEDGE=6;Elttype(ET).NFACE=4;Elttype(ET).NTET=1;ELTTYPE(ET).DIM=3
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE),Elttype(ET).TET(4,Elttype(ET).NTET))	
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,1,1,4,2,4,3,4],(/2,6/))
				Elttype(ET).FACE(:,:)=RESHAPE([3,2,1,3,0,&
											   3,1,2,4,0,&
											   3,2,3,4,0,&
											   3,3,1,4,0],(/5,4/))
				Elttype(ET).FACEEDGE(:,:)=RESHAPE([3,-1,-3,-2,0,&
											   3,1,5,-4,0,&
											   3,2,6,-5,0,&
											   3,3,4,-6,0],(/5,4/))
                Elttype(ET).TET(:,:)=RESHAPE([1,2,3,4], (/4,1/))

			CASE(11) !10-noded-TETRAHEDRON
                ELTTYPE(ET).NNODE=4
				Elttype(ET).NEDGE=6;Elttype(ET).NFACE=4;Elttype(ET).NTET=8;ELTTYPE(ET).DIM=3
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE),Elttype(ET).TET(4,Elttype(ET).NTET))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,1,1,4,2,4,3,4],(/2,6/))
				Elttype(ET).FACE(:,:)=RESHAPE([3,2,1,3,0,&
											   3,1,2,4,0,&
											   3,2,3,4,0,&
											   3,3,1,4,0],(/5,4/))
				Elttype(ET).FACEEDGE(:,:)=RESHAPE([3,-1,-3,-2,0,&
											   3,1,5,-4,0,&
											   3,2,6,-5,0,&
											   3,3,4,-6,0],(/5,4/))
				Elttype(ET).TET(:,:)=RESHAPE([7,1,5,8,&
												7,3,10,6,&
												10,4,8,9,&
												5,2,6,9,&
												10,9,5,6,&
												10,5,7,6,&
												5,9,10,8,&
												7,5,10,8],&
                                                (/4,8/))

			CASE(5) !HEXAHEDRON
                ELTTYPE(ET).NNODE=8
				Elttype(ET).NEDGE=12;Elttype(ET).NFACE=6;Elttype(ET).NTET=6;ELTTYPE(ET).DIM=3
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE),Elttype(ET).TET(4,Elttype(ET).NTET))				
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,4,4,1,&
											   5,6,6,7,7,8,8,5,&
											   1,5,2,6,3,7,4,8],(/2,12/))
				Elttype(ET).FACE(:,:)=RESHAPE([4,4,3,2,1,&
											   4,5,6,7,8,&
											   4,1,2,6,5,&
											   4,2,3,7,6,&
											   4,3,4,8,7,&
											   4,4,1,5,8],(/5,6/))
											   	
				Elttype(ET).FACEEDGE(:,:)=RESHAPE([4,-3,-2,-1,-4,&
												   4,5,6,7,8,&
												   4,1,10,-5,-9,&
												   4,2,11,6,-10,&
												   4,3,4,12,-7,&
												   4,4,9,-8,-12],(/5,18/))
				Elttype(ET).TET(:,:)=RESHAPE([1,2,3,6,&
												1,7,5,6,&
												1,3,7,6,&
												1,3,4,8,&
												1,7,3,8,&
												1,5,7,8],&
                                                (/4,6/))
			
			CASE(6) !6-NODE PRISM
                ELTTYPE(ET).NNODE=6
				Elttype(ET).NEDGE=9;Elttype(ET).NFACE=5;Elttype(ET).NTET=3;ELTTYPE(ET).DIM=3
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE),Elttype(ET).TET(4,Elttype(ET).NTET))				
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,1,&
											   4,5,5,6,6,4,&
											   1,4,2,5,3,6],(/2,9/))
				Elttype(ET).FACE(:,:)=RESHAPE([3,2,1,3,0,&
											   3,4,5,6,0,&
											   4,1,2,5,4,&
											   4,2,3,6,5,&
											   4,3,1,4,6],(/5,5/))
				Elttype(ET).FACEEDGE(:,:)=RESHAPE([3,-1,-3,-2,0,&
											   3,4,5,6,0,&
											   4,1,8,-4,-7,&
											   4,2,9,-5,-8,&
											   4,3,7,-6,-9],(/5,5/))
				Elttype(ET).TET(:,:)=RESHAPE([1,2,3,5,&
												1,3,6,5,&
												1,6,4,5],(/4,3/))

			CASE(18) !15-NODE PRISM
                ELTTYPE(ET).NNODE=6
				Elttype(ET).NEDGE=9;Elttype(ET).NFACE=5;Elttype(ET).NTET=14;ELTTYPE(ET).DIM=3
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE),Elttype(ET).TET(4,Elttype(ET).NTET))			
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,1,&
											   4,5,5,6,6,4,&
											   1,4,2,5,3,6],(/2,9/))
				Elttype(ET).FACE(:,:)=RESHAPE([3,2,1,3,0,&
											   3,4,5,6,0,&
											   4,1,2,5,4,&
											   4,2,3,6,5,&
											   4,3,1,4,6],(/5,5/))
				Elttype(ET).FACEEDGE(:,:)=RESHAPE([3,-1,-3,-2,0,&
											   3,4,5,6,0,&
											   4,1,8,-4,-7,&
											   4,2,9,-5,-8,&
											   4,3,7,-6,-9],(/5,5/))
				Elttype(ET).TET(:,:)=RESHAPE([15,13,8,14,&
												13,15,11,14,&
												4,13,12,10,&
												9,7,13,1,&
												13,15,12,11,&
												9,15,13,8,&
												8,13,7,14,&
												10,13,11,14,&
												13,12,10,11,&
												9,7,8,13,&
												15,9,3,8,&
												15,6,12,11,&
												7,8,14,2,&
												5,10,11,14],&
                                                (/4,14/))
			CASE(7) !5-node pyramid
                ELTTYPE(ET).NNODE=5
				Elttype(ET).NEDGE=8;Elttype(ET).NFACE=5;Elttype(ET).NTET=2;ELTTYPE(ET).DIM=3
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE),Elttype(ET).TET(4,Elttype(ET).NTET))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,4,4,1,&
											   1,5,2,5,3,5,4,5],(/2,8/))
				Elttype(ET).FACE(:,:)=RESHAPE([4,4,3,2,1,&
											   3,1,2,5,0,&
											   3,2,3,5,0,&
											   3,3,4,5,0,&
											   3,4,1,5,0],(/5,5/))
				Elttype(ET).FACEEDGE(:,:)=RESHAPE([4,-3,-2,-1,-4,&
											   3,1,6,-5,0,&
											   3,2,7,-6,0,&
											   3,3,8,-7,0,&
											   3,4,5,-8,0],(/5,5/))	
				Elttype(ET).TET(:,:)=RESHAPE([1,2,3,5,&
												1,3,4,5],&
                                                (/4,2/))
				
				
			CASE(15) !POINTS
                ELTTYPE(ET).NNODE=1
				Elttype(ET).NEDGE=0;Elttype(ET).NFACE=0
				ELTTYPE(ET).DIM=0
            !CASE DEFAULT
            !    PRINT *, 'NO SUCH GMSH ET OR STILL UNDER WORK.'
            !    STOP
		ENDSELECT
	ENDDO
	
ENDSUBROUTINE

subroutine Initialize_et2numNodes()
	implicit none
	integer::et1
	
	Elttype(1).nnode=2
	Elttype(1).description="2-node line."
	et1=1
	call not_nodal_force_weight(et1)

	
	Elttype(2).nnode=3
	Elttype(2).description="3-node triangle."
	et1=2
	call not_nodal_force_weight(et1)
	
	Elttype(3).nnode=4
	Elttype(3).description="4-node quadrangle."
	et1=3
	call not_nodal_force_weight(et1)
	
	Elttype(4).nnode=4
	Elttype(4).description="4-node tetrahedron."
	et1=4
	call not_nodal_force_weight(et1)
	
	Elttype(5).nnode=8
	Elttype(5).description="8-node hexahedron."
	et1=5
	call not_nodal_force_weight(et1)
	
	Elttype(6).nnode=6
	Elttype(6).description="6-node prism."
	et1=6
	call not_nodal_force_weight(et1)	
	
	Elttype(7).nnode=5
	Elttype(7).description="5-node pyramid."
	
	Elttype(8).nnode=3
	Elttype(8).description="3-node second order line (2 nodes associated with the vertices and 1 with the edge)."
	et1=8
	call not_nodal_force_weight(et1)	
	
	Elttype(9).nnode=6
	Elttype(9).description="6-node second order triangle (3 nodes associated with the vertices and 3 with the edges)."
	et1=9
	call not_nodal_force_weight(et1)	
	
	Elttype(10).nnode=9
	Elttype(10).description="9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face)."
		
	Elttype(11).nnode=10
	Elttype(11).description="10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges)."
	et1=11
	call not_nodal_force_weight(et1)	
	
	Elttype(12).nnode=27
	Elttype(12).description="27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume)."
	Elttype(13).nnode=18
	Elttype(13).description="18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces)."
	Elttype(14).nnode=14
	Elttype(14).description="14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face)."
	Elttype(15).nnode=1
	Elttype(15).description="1-node point."
	
	Elttype(16).nnode=8
	Elttype(16).description="8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges)."
	et1=16
	call not_nodal_force_weight(et1)	
	
	
	Elttype(17).nnode=20
	Elttype(17).description="20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges)."
	
	Elttype(18).nnode=15
	Elttype(18).description="15-node second order prism (6 nodes associated with the vertices and 9 with the edges)."
	
	Elttype(19).nnode=13
	Elttype(19).description="13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges)."
	Elttype(20).nnode=9
	Elttype(20).description="9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)"
	Elttype(21).nnode=10
	Elttype(21).description="10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)"
	Elttype(22).nnode=12
	Elttype(22).description="12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)"
	
	Elttype(23).nnode=15
	Elttype(23).description="15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)"
	et1=23
	call not_nodal_force_weight(et1)	
	

	
	Elttype(24).nnode=15	
	Elttype(24).description="15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)"
	Elttype(25).nnode=21
	Elttype(25).description="21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)"
	Elttype(26).nnode=4
	Elttype(26).description="4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge)"
	
	Elttype(27).nnode=5
	Elttype(27).description="5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)"
	et1=27
	call not_nodal_force_weight(et1)	

	
	Elttype(28).nnode=6
	Elttype(28).description="6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge)"
	Elttype(29).nnode=20
	Elttype(29).description="20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)"
	Elttype(30).nnode=35
	Elttype(30).description="35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)"
	Elttype(31).nnode=56
	Elttype(31).description="56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)"
	Elttype(92).nnode=64
	Elttype(92).description="64-node third order hexahedron (8 nodes associated with the vertices, 24 with the edges, 24 with the faces, 8 in the volume)"
	Elttype(93).nnode=125
	Elttype(93).description="125-node fourth order hexahedron (8 nodes associated with the vertices, 36 with the edges, 54 with the faces, 27 in the volume)"

	
end subroutine


subroutine not_nodal_force_weight(et)
	implicit none
	integer,intent(in)::et
	
	select case(et)
		case(1)
			allocate(Elttype(et).weight(Elttype(et).nnode,2))
			Elttype(et).weight(:,1)=0.5
			Elttype(et).weight(1,2)=1./6.
			Elttype(et).weight(2,2)=1./3.
		case(2,3,4,5,6)
			allocate(Elttype(et).weight(Elttype(et).nnode,1))
			Elttype(et).weight=1./Elttype(et).nnode
		case(8)
			allocate(Elttype(et).weight(Elttype(et).nnode,2))
			Elttype(et).weight(1:2,1)=1./6.
			Elttype(et).weight(3,1)=2./3.
			Elttype(et).weight(1,2)=0.
			Elttype(et).weight(2,2)=1./6.
			Elttype(et).weight(3,2)=1./3.
		case(9)
			allocate(Elttype(et).weight(Elttype(et).nnode,1))
			Elttype(et).weight=0.
			Elttype(et).weight(4:6,1)=1./3.
		case(11)	
			allocate(Elttype(et).weight(Elttype(et).nnode,1))
			Elttype(et).weight=0.
			Elttype(et).weight(1:4,1)=-1./20.
			Elttype(et).weight(5:10,1)=1./5.
		case(16)
			allocate(Elttype(et).weight(Elttype(et).nnode,1))
			Elttype(et).weight=0.
			Elttype(et).weight(1:4,1)=-1./12.
			Elttype(et).weight(5:8,1)=1./3.
		case(23) !CPE15
			allocate(Elttype(et).weight(Elttype(et).nnode,1))
			Elttype(et).weight=0.D0
			Elttype(et).weight(4:6,1)=-1./45.
			Elttype(et).weight(7:12,1)=4./45.
			Elttype(et).weight(13:15,1)=8./45.
		case(27)
			allocate(Elttype(et).weight(Elttype(et).nnode,2))
			Elttype(et).weight(1,1)=7./90.
			Elttype(et).weight(2,1)=7./90.
			Elttype(et).weight(3,1)=16./45.
			Elttype(et).weight(5,1)=16./45.
			Elttype(et).weight(4,1)=2/15.
			Elttype(et).weight(1,2)=0
			Elttype(et).weight(2,2)=7./90.
			Elttype(et).weight(3,2)=4/45.
			Elttype(et).weight(5,2)=4./15.
			Elttype(et).weight(4,2)=1/15.
		end select
end subroutine

SUBROUTINE I_ENLARGE_AR(AVAL,DSTEP)
    IMPLICIT NONE
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
    IMPLICIT NONE
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




 

ENDMODULE

 



subroutine r8vec_cross_3d ( v1, v2, v3 )

!*****************************************************************************80
!
!! R8VEC_CROSS_3D computes the cross product of two vectors in 3D.
!
!  Discussion:
!
!    The cross product in 3D can be regarded as the determinant of the
!    symbolic matrix:
!
!          |  i  j  k |
!      det | x1 y1 z1 |
!          | x2 y2 z2 |
!
!      = ( y1 * z2 - z1 * y2 ) * i
!      + ( z1 * x2 - x1 * z2 ) * j
!      + ( x1 * y2 - y1 * x2 ) * k
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), the two vectors.
!
!    Output, real ( kind = 8 ) V3(3), the cross product vector.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) v1(dim_num)
  real ( kind = 8 ) v2(dim_num)
  real ( kind = 8 ) v3(dim_num)

  v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
  v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
  v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

  return
end




SUBROUTINE R2_ENLARGE_AR(AVAL,DSTEP,DIM1)
    IMPLICIT NONE
    REAL(8),ALLOCATABLE,INTENT(INOUT)::AVAL(:,:)
    INTEGER,INTENT(IN)::DSTEP,DIM1
    REAL(8),ALLOCATABLE::VAL1(:,:)
    INTEGER::LB1=0,UB1=0
    
    LB1=LBOUND(AVAL,DIM=2);UB1=UBOUND(AVAL,DIM=2)
    ALLOCATE(VAL1,SOURCE=AVAL)
    DEALLOCATE(AVAL)
    ALLOCATE(AVAL(DIM1,LB1:UB1+DSTEP))
    AVAL(:,LB1:UB1)=VAL1
    !AVAL(UB1+1:UB1+10)=0
    DEALLOCATE(VAL1)
END SUBROUTINE

SUBROUTINE I2_ENLARGE_AR(AVAL,DSTEP,DIM1)
    IMPLICIT NONE
    INTEGER,ALLOCATABLE,INTENT(INOUT)::AVAL(:,:)
    INTEGER,INTENT(IN)::DSTEP,DIM1
    INTEGER,ALLOCATABLE::VAL1(:,:)
    INTEGER::LB1=0,UB1=0
    
    LB1=LBOUND(AVAL,DIM=2);UB1=UBOUND(AVAL,DIM=2)
    ALLOCATE(VAL1,SOURCE=AVAL)
    DEALLOCATE(AVAL)
    ALLOCATE(AVAL(DIM1,LB1:UB1+DSTEP))
    AVAL(:,LB1:UB1)=VAL1
    !AVAL(UB1+1:UB1+10)=0
    DEALLOCATE(VAL1)
END SUBROUTINE


!!SET EDGE FOR TRIANGLE AND TET. 
!SUBROUTINE SETUP_EDGE_TBLi(EDGE_TBL_L,TBL_SIZE_L,EDGE_L,NEDGE_L,ELEMENT_L,NEL_L,NODE_L,NNODE_L)
!    USE MESHGEO
!    USE hashtbl
!    IMPLICIT NONE
!	INTEGER,INTENT(IN)::TBL_SIZE_L,NEL_L,NNODE_L
!    INTEGER::NEDGE_L
!	TYPE(hash_tbl_sll)::EDGE_TBL_L
!	TYPE(EDGE_TYDEF),ALLOCATABLE::EDGE_L(:)
!	type(TET_TYDEF)::ELEMENT_L(NEL_L)
!	REAL(8),INTENT(IN)::NODE_L(3,NNODE_L)
!    INTEGER::I,J,N1(2),ET1,NEDGE1
!	CHARACTER(LEN=:),ALLOCATABLE::KEY1
!    TYPE(DICT_DATA)::VAL1
!    CHARACTER(64)::CKEY1
!    INTEGER::HKEY1
!    
!     
!    
!	IF(.NOT.ALLOCATED(EDGE_L)) ALLOCATE(EDGE_L(TBL_SIZE_L))
!	CALL EDGE_TBL_L.FREE
!    IF(.NOT.EDGE_TBL_L.IS_INIT) CALL EDGE_TBL_L.INIT(TBL_SIZE_L)
!	DO I=1,NEL_L
!		ET1=ELEMENT_L(I).GMET
!		NEDGE1=ELTTYPE(ET1).NEDGE
!        ELEMENT_L(I).NE=NEDGE1
!        !IF(.NOT.ALLOCATED(ELEMENT_L(I).EDGE)) ALLOCATE(ELEMENT_L(I).EDGE(NEDGE1))
!		DO J=1,NEDGE1
!            N1=ELEMENT_L(I).V(ELTTYPE(ET1).EDGE(:,J))
!            CALL I2C_KEY_hash_tbl_sll(KEY1,N1(:),2)
!            CKEY1=TRIM(KEY1)
!            HKEY1 = MOD(ABS(HASH_DJB(KEY1)),EDGE_TBL_L%vec_len)
!            VAL1.IEL=I;VAL1.ISE=J;VAL1.IITEM=EDGE_TBL_L.TBL_ID
!			CALL EDGE_TBL_L.PUT(TRIM(KEY1),VAL1,HKEY1)
!			ELEMENT_L(I).E(J)=VAL1.IITEM
!            IF(VAL1.IITEM>SIZE(EDGE_L,DIM=1)) THEN
!                CALL EDGE_TYDEF_ENLARGE_AR(EDGE_L,1000)
!                !MAXNEDGE=MAXNEDGE+1000
!            ENDIF
!            IF(.NOT.EDGE_L(VAL1.IITEM).ISINI) THEN
!                EDGE_L(VAL1.IITEM).CKEY=CKEY1
!                EDGE_L(VAL1.IITEM).HKEY=HKEY1
!                EDGE_L(VAL1.IITEM).V=ELEMENT_L(I).V(ELTTYPE(ET1).EDGE(:,J))
!                EDGE_L(VAL1.IITEM).DIS=NORM2(NODE_L(:,EDGE_L(VAL1.IITEM).V(1))-NODE_L(:,EDGE_L(VAL1.IITEM).V(2)))
!                EDGE_L(VAL1.IITEM).ISINI=.TRUE.
!                NEDGE_L=VAL1.IITEM
!            ENDIF
!            IF(EDGE_L(VAL1.IITEM).ENUM==0) ALLOCATE(EDGE_L(VAL1.IITEM).ELEMENT(5),EDGE_L(VAL1.IITEM).SUBID(5))
!            EDGE_L(VAL1.IITEM).ENUM=EDGE_L(VAL1.IITEM).ENUM+1
!            IF(EDGE_L(VAL1.IITEM).ENUM>SIZE(EDGE_L(VAL1.IITEM).ELEMENT,DIM=1)) THEN
!				CALL I_ENLARGE_AR(EDGE_L(VAL1.IITEM).ELEMENT,5)
!				CALL I_ENLARGE_AR(EDGE_L(VAL1.IITEM).SUBID,5)
!			ENDIF
!            EDGE_L(VAL1.IITEM).ELEMENT(EDGE_L(VAL1.IITEM).ENUM)=I			
!            EDGE_L(VAL1.IITEM).SUBID(EDGE_L(VAL1.IITEM).ENUM)=J
!		ENDDO
!    END DO
!    RETURN
!ENDSUBROUTINE
!
!
!SUBROUTINE SETUP_FACE_TBLi(FACE_TBL_L,TBL_SIZE_L,FACE_L,NFACE_L,EDGE_L,NEDGE_L,ELEMENT_L,NEL_L,NODE_L,NNODE_L)
!    use MESHGEO
!    use hashtbl
!    IMPLICIT NONE
!    INTEGER,INTENT(IN)::TBL_SIZE_L,NEL_L,NEDGE_L,NNODE_L
!    INTEGER::NFACE_L
!	TYPE(hash_tbl_sll)::FACE_TBL_L
!	TYPE(EDGE_TYDEF),INTENT(IN)::EDGE_L(NEDGE_L)
!    TYPE(FACE_TYDEF),ALLOCATABLE::FACE_L(:)
!	type(TET_TYDEF)::ELEMENT_L(NEL_L)
!    REAL(8),INTENT(IN)::NODE_L(3,NNODE_L)
!    INTEGER::I,J,N1(4),ET1,NFACE1,TBL_LEN,N2,K
!    REAL(8)::V1(3),V2(3),NORMAL1(3),T1
!	CHARACTER(LEN=:),ALLOCATABLE::KEY1
!    TYPE(DICT_DATA)::VAL1
!    CHARACTER(64)::CKEY1
!    INTEGER::HKEY1
!    
!    IF(.NOT.ISINI_GMSHET) THEN
!        CALL Initialize_et2numNodes()
!        CALL ET_GMSH_EDGE_FACE()
!        ISINI_GMSHET=.TRUE.
!    ENDIF 	
!
!	IF(ALLOCATED(FACE_L)) DEALLOCATE(FACE_L)
!	ALLOCATE(FACE_L(TBL_SIZE_L))
!    CALL FACE_TBL_L.FREE
!	IF(.NOT.FACE_TBL_L.IS_INIT) CALL FACE_TBL_L.INIT(TBL_SIZE_L)
!
!	DO I=1,NEL_L
!		ET1=ELEMENT_L(I).GMET
!		NFACE1=ELTTYPE(ET1).NFACE
!		ELEMENT_L(I).NF=NFACE1
!        !ALLOCATE(ELEMENT_L(I).F(NFACE1))
!		DO J=1,NFACE1
!            N2=ELTTYPE(ET1).FACE(0,J)
!            N1(1:N2)=ELEMENT_L(I).V(ELTTYPE(ET1).FACE(1:N2,J))
!            CALL FACE_TBL_L.KEY(KEY1,N1(1:N2),ELTTYPE(ET1).FACE(0,J))
!            CKEY1=TRIM(KEY1)
!            HKEY1 = MOD(ABS(HASH_DJB(KEY1)),FACE_TBL_L%vec_len)            
!            VAL1.IEL=I;VAL1.ISE=J;VAL1.IITEM=FACE_TBL_L.TBL_ID
!			CALL FACE_TBL_L.PUT(TRIM(KEY1),VAL1,HKEY1)
!            ELEMENT_L(I).F(J)=VAL1.IITEM
!            
!            IF(VAL1.IITEM>SIZE(FACE_L,DIM=1)) THEN
!                CALL FACE_TYDEF_ENLARGE_AR(FACE_L,1000)
!                !MAXNFACE=MAXNFACE+1000
!            ENDIF
!            IF(.NOT.FACE_L(VAL1.IITEM).ISINI) THEN
!                FACE_L(VAL1.IITEM).CKEY=CKEY1 
!                FACE_L(VAL1.IITEM).HKEY=HKEY1
!                FACE_L(VAL1.IITEM).SHAPE=ELTTYPE(ET1).FACE(0,J)
!                FACE_L(VAL1.IITEM).V(1:FACE_L(VAL1.IITEM).SHAPE)=ELEMENT_L(I).V(ELTTYPE(ET1).FACE(1:FACE_L(VAL1.IITEM).SHAPE,J))
!				FACE_L(VAL1.IITEM).EDGE(1:FACE_L(VAL1.IITEM).SHAPE)=ELEMENT_L(I).E(ABS(ELTTYPE(ET1).FACEEDGE(1:FACE_L(VAL1.IITEM).SHAPE,J)))
!				!CHECK LOOP ORDER
!				DO K=1,FACE_L(VAL1.IITEM).SHAPE
!					IF(EDGE_L(FACE_L(VAL1.IITEM).EDGE(K)).V(1)/=FACE_L(VAL1.IITEM).V(K)) FACE_L(VAL1.IITEM).EDGE(K)=-FACE_L(VAL1.IITEM).EDGE(K)
!				ENDDO
!                V1=NODE_L(:,FACE_L(VAL1.IITEM).V(2))-NODE_L(:,FACE_L(VAL1.IITEM).V(1))
!                V2=NODE_L(:,FACE_L(VAL1.IITEM).V(3))-NODE_L(:,FACE_L(VAL1.IITEM).V(1))
!                call r8vec_cross_3d ( v1, v2, NORMAL1(1:3) )
!                T1 = sqrt ( sum ( ( NORMAL1 )**2 ) )
!                if ( T1 /= 0.0D+00 ) then
!                  NORMAL1 = NORMAL1/T1
!                end if
!                FACE_L(VAL1.IITEM).UNORMAL=NORMAL1
!                FACE_L(VAL1.IITEM).ISINI=.TRUE.
!                NFACE_L=VAL1.IITEM
!				DO K=1,3
!					FACE_L(NFACE_L).BBOX(1,K)=MINVAL(NODE_L(K,FACE_L(NFACE_L).V(1:FACE_L(NFACE_L).SHAPE)))-VTOL
!					FACE_L(NFACE_L).BBOX(2,K)=MAXVAL(NODE_L(K,FACE_L(NFACE_L).V(1:FACE_L(NFACE_L).SHAPE)))+VTOL
!				ENDDO
!            ENDIF
!            IF(FACE_L(VAL1.IITEM).ENUM==0) ALLOCATE(FACE_L(VAL1.IITEM).ELEMENT(2))
!            FACE_L(VAL1.IITEM).ENUM=FACE_L(VAL1.IITEM).ENUM+1
!            IF(FACE_L(VAL1.IITEM).ENUM>SIZE(FACE_L(VAL1.IITEM).ELEMENT,DIM=1)) CALL I_ENLARGE_AR(FACE_L(VAL1.IITEM).ELEMENT,5)
!            
!            IF(FACE_L(VAL1.IITEM).ENUM>1) THEN
!                DO K=1,FACE_L(VAL1.IITEM).SHAPE
!                    IF(FACE_L(VAL1.IITEM).V(K)==ELEMENT_L(I).V(ELTTYPE(ELEMENT_L(I).GMET).FACE(1,J))) THEN
!                        IF(FACE_L(VAL1.IITEM).V(MOD(K,FACE_L(VAL1.IITEM).SHAPE)+1)/= &
!                            ELEMENT_L(I).V(ELTTYPE(ELEMENT_L(I).GMET).FACE(2,J))) THEN
!                            FACE_L(VAL1.IITEM).ELEMENT(FACE_L(VAL1.IITEM).ENUM)=-I !I��Ԫ����ķ�����face�ķ�����,Ϊ-1
!                        ELSE
!                            FACE_L(VAL1.IITEM).ELEMENT(FACE_L(VAL1.IITEM).ENUM)=I
!                        ENDIF                  
!                        EXIT
!                    ENDIF
!                ENDDO
!            ELSE
!                FACE_L(VAL1.IITEM).ELEMENT(FACE_L(VAL1.IITEM).ENUM)=I !ͬ�� =1
!            ENDIF
!            
!                
!            
!		ENDDO
!    END DO    
!    
!    
!ENDSUBROUTINE


