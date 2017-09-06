MODULE MESHGEO

USE solverds
USE hashtbl

PRIVATE
PUBLIC::EDGE,NEDGE,FACE,NFACE,SETUP_EDGE_TBL_TET,SETUP_FACE_TBL_TET,TET,NTET,&
        SETUP_SUB_TET4_ELEMENT

LOGICAL::ISINI_GMSHET=.FALSE.

TYPE EDGE_TYDEF
    LOGICAL::ISINI=.FALSE.
    INTEGER::V(2)=0
    INTEGER::HKEY=-1
    CHARACTER(64)::CKEY=""
    INTEGER::ENUM=0
    INTEGER,ALLOCATABLE::ELEMENT(:),SUBID(:) !SUBID IS INDEX WHICH EDGE OF THE ELEMENT IS THE EDGE
ENDTYPE
TYPE(EDGE_TYDEF),ALLOCATABLE::EDGE(:)
TYPE FACE_TYDEF
    LOGICAL::ISINI=.FALSE.
    INTEGER::SHAPE=3
    INTEGER::V(4)=0,EDGE(4)=0 !IF EDGE(I)<0 MEANING ITS ORDER IS REVERSE .
    INTEGER::HKEY=-1
    CHARACTER(64)::CKEY="" 		
    INTEGER::ISTRISURFACE=0 !<-1 与face反向
    REAL(8)::UNORMAL(3)=0.D0,BBOX(2,3)=0.D0		
    INTEGER::ENUM=0
    INTEGER,ALLOCATABLE::ELEMENT(:),SUBID(:) !SUBID IS INDEX WHICH FACE OF THE ELEMENT IS THE FACE
ENDTYPE
TYPE(FACE_TYDEF),ALLOCATABLE::FACE(:)    
INTEGER::NEDGE=0,NFACE=0,MAXNEDGE=10000,MAXNFACE=10000

!split all element into simple 3-noded-triangle , 4-noded-tetrahedron,point and 2-noded-line element
TYPE TET_TYDEF
	INTEGER::MOTHER=0,GMET=0,DIM=-1,NV=0,NE=0,NF=0 
	INTEGER::V(4)=0,E(6)=0,F(4)=0	!F的正反没有维护，可以通过elttype(4).face来确定。
	REAL(8)::BBOX(2,3)=0.D0 !MIN,MAX OF X,Y AND Z 
ENDTYPE
TYPE(TET_TYDEF),ALLOCATABLE::TET(:)
INTEGER::NTET=0

type et_type
	integer::nnode=0,nedge=0,nface=0,ntet=0		
	character(512)::description
	integer,allocatable::edge(:,:),face(:,:),FaceEdge(:,:)
	integer,allocatable::tet(:,:),tetEDGE(:,:),TETFACE(:,:)
	INTEGER::DIM=-1
	!edge(2,nedge),
	!face: use node index to represent face. face(0:4,nface),face(0,:)==3,triangular face,==4, quadrilateral face
	!FaceEdge: use edge index to represent face. FaceEdge(0:4,nface),FaceEdge(0,:)==3,triangular face, ==4, quadrilateral face
	real(8),allocatable::weight(:,:) !分布单元荷载各节点分布系数, weight(:,1) 平面均布荷载；weight(:,2) 平面三角形荷载;weight(:,3) 轴对称均布荷载；weight(:,4) 轴对称三角形荷载，
												!目前只能处理平面应变的两种情况。
end type
type(et_type)::elttype(100)

REAL(8)::VTOL=1.D-3

    
CONTAINS 

SUBROUTINE SETUP_SUB_TET4_ELEMENT()
	IMPLICIT NONE
	INTEGER::i,J,K,GMET1
	
    IF(.NOT.ISINI_GMSHET) THEN
        IF(.NOT.ALLOCATED(EDGE)) ALLOCATE(EDGE(MAXNEDGE),FACE(MAXNFACE))
        CALL Initialize_et2numNodes()
        CALL ET_GMSH_EDGE_FACE()
        ISINI_GMSHET=.TRUE.
    ENDIF
	
	NTET=0
	DO I=1,ENUM
        IF(ESET(ELEMENT(I).SET).COUPLESET<ELEMENT(I).SET) CYCLE !影子单元不加入
        
		GMET1=GETGMSHET(ELEMENT(I).ET)
		SELECT CASE(ELTTYPE(GMET1).DIM)
		CASE(0,1)
			NTET=NTET+1	
			ELEMENT(I).NTET=1
		CASE(2)
			NTET=NTET+Elttype(GMET1).NFACE
			ELEMENT(I).NTET=Elttype(GMET1).NFACE
		CASE(3)
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
	DO I=1,ENUM
        IF(ESET(ELEMENT(I).SET).COUPLESET<ELEMENT(I).SET) CYCLE !影子单元不加入
        
		GMET1=GETGMSHET(ELEMENT(I).ET)
		ELEMENT(I).TET=[NTET+1:NTET+ELEMENT(I).NTET]
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
				TET(NTET).BBOX(1,K)=MINVAL(NODE(TET(NTET).V(1:TET(NTET).NV)).COORD(K))-VTOL
				TET(NTET).BBOX(2,K)=MAXVAL(NODE(TET(NTET).V(1:TET(NTET).NV)).COORD(K))+VTOL
			ENDDO
		CASE(2)
			DO J=1,Elttype(GMET1).NFACE				
				NTET=NTET+1
				TET(NTET).MOTHER=I
				TET(NTET).NV=3;TET(NTET).GMET=2;TET(NTET).DIM=2
				TET(NTET).V(1:3)=ELEMENT(I).NODE(Elttype(GMET1).FACE(1:3,J))
				DO K=1,3
					TET(NTET).BBOX(1,K)=MINVAL(NODE(TET(NTET).V(1:TET(NTET).NV)).COORD(K))-VTOL
					TET(NTET).BBOX(2,K)=MAXVAL(NODE(TET(NTET).V(1:TET(NTET).NV)).COORD(K))+VTOL
				ENDDO
			ENDDO
		CASE(3)
			DO J=1,Elttype(GMET1).NTET				
				NTET=NTET+1
				TET(NTET).MOTHER=I
				TET(NTET).NV=4;TET(NTET).GMET=4;TET(NTET).DIM=3
				TET(NTET).V(1:4)=ELEMENT(I).NODE(Elttype(GMET1).TET(:,J))
				DO K=1,3
					TET(NTET).BBOX(1,K)=MINVAL(NODE(TET(NTET).V(1:TET(NTET).NV)).COORD(K))-VTOL
					TET(NTET).BBOX(2,K)=MAXVAL(NODE(TET(NTET).V(1:TET(NTET).NV)).COORD(K))+VTOL
				ENDDO
			ENDDO
		END SELECT
		
	ENDDO	
	
END SUBROUTINE

INTEGER FUNCTION GETGMSHET(ET)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ET
    
    SELECT CASE(ET)
    !LINE,2-NODE LINE
    CASE(BAR,BEAM,BAR2D,BEAM2D)    
        GETGMSHET=1
    CASE(CPE3,CPE3_SPG,CPE3_CPL,CPS3,SHELL3,DKT3,&
         CAX3,CAX3_SPG,CAX3_CPL)
         GETGMSHET=2
    CASE(CPE4,CPE4_SPG,CPE4_CPL,CPS4,&
         CAX4,CAX4_SPG,CAX4_CPL,&
         CPE4R,CPE4R_SPG,CPE4R_CPL,CPS4R,&
         CAX4R,CAX4R_SPG,CAX4R_CPL)
         GETGMSHET=3 
    CASE(TET4,TET4_SPG,TET4_CPL)         
        GETGMSHET=4         
    CASE(PRM6,PRM6_SPG,PRM6_CPL)
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

    
SUBROUTINE SETUP_EDGE_TBL_TET()
    
    IMPLICIT NONE
    INTEGER::I,J,N1(2),ET1,NEDGE1,TBL_LEN
	CHARACTER(LEN=:),ALLOCATABLE::KEY1
    TYPE(DICT_DATA)::VAL1
    CHARACTER(64)::CKEY1
    INTEGER::HKEY1
    
    IF(.NOT.ISINI_GMSHET) THEN
        IF(.NOT.ALLOCATED(EDGE)) ALLOCATE(EDGE(MAXNEDGE),FACE(MAXNFACE))
        CALL Initialize_et2numNodes()
        CALL ET_GMSH_EDGE_FACE()
        ISINI_GMSHET=.TRUE.
    ENDIF
    
	TBL_LEN=2*NNUM	
    CALL EDGE_TBL.INIT(TBL_LEN)
	DO I=1,NTET
		ET1=TET(I).GMET
		NEDGE1=ELTTYPE(ET1).NEDGE
		TET(I).NE=NEDGE1
        !ALLOCATE(TET(I).E(NEDGE1))
		DO J=1,NEDGE1
            N1=TET(I).V(ELTTYPE(ET1).EDGE(:,J))
            CALL I2C_KEY_hash_tbl_sll(KEY1,N1(:),2)
            CKEY1=TRIM(KEY1)
            HKEY1 = MOD(ABS(HASH_DJB(KEY1)),EDGE_TBL%vec_len)
            VAL1.IEL=I;VAL1.ISE=J;VAL1.IITEM=EDGE_TBL.TBL_ID
			CALL EDGE_TBL.PUT(TRIM(KEY1),VAL1,HKEY1)
			TET(I).E(J)=VAL1.IITEM
            IF(VAL1.IITEM>MAXNEDGE) THEN
                CALL EDGE_TYDEF_ENLARGE_AR(EDGE,1000)
                MAXNEDGE=MAXNEDGE+1000
            ENDIF
            IF(.NOT.EDGE(VAL1.IITEM).ISINI) THEN
                EDGE(VAL1.IITEM).CKEY=CKEY1
                EDGE(VAL1.IITEM).HKEY=HKEY1
                EDGE(VAL1.IITEM).V=TET(I).V(ELTTYPE(ET1).EDGE(:,J))
                EDGE(VAL1.IITEM).ISINI=.TRUE.
                NEDGE=VAL1.IITEM
            ENDIF
            IF(EDGE(VAL1.IITEM).ENUM==0) ALLOCATE(EDGE(VAL1.IITEM).ELEMENT(5),EDGE(VAL1.IITEM).SUBID(5))
            EDGE(VAL1.IITEM).ENUM=EDGE(VAL1.IITEM).ENUM+1
            IF(EDGE(VAL1.IITEM).ENUM>SIZE(EDGE(VAL1.IITEM).ELEMENT,DIM=1)) THEN
				CALL I_ENLARGE_AR(EDGE(VAL1.IITEM).ELEMENT,5)
				CALL I_ENLARGE_AR(EDGE(VAL1.IITEM).SUBID,5)
			ENDIF
            EDGE(VAL1.IITEM).ELEMENT(EDGE(VAL1.IITEM).ENUM)=I			
            EDGE(VAL1.IITEM).SUBID(EDGE(VAL1.IITEM).ENUM)=J
		ENDDO
    END DO
    RETURN
ENDSUBROUTINE

SUBROUTINE SETUP_FACE_TBL_TET()
    IMPLICIT NONE
    INTEGER::I,J,N1(4),ET1,NFACE1,TBL_LEN,N2,K
    REAL(8)::V1(3),V2(3),NORMAL1(3),T1
	CHARACTER(LEN=:),ALLOCATABLE::KEY1
    TYPE(DICT_DATA)::VAL1
    CHARACTER(64)::CKEY1
    INTEGER::HKEY1
    
    IF(.NOT.ISINI_GMSHET) THEN
        IF(.NOT.ALLOCATED(EDGE)) ALLOCATE(EDGE(MAXNEDGE),FACE(MAXNFACE))
        CALL Initialize_et2numNodes()
        CALL ET_GMSH_EDGE_FACE()
        ISINI_GMSHET=.TRUE.
    ENDIF    
    
	TBL_LEN=2*NNUM	
    CALL FACE_TBL.INIT(TBL_LEN)
	DO I=1,NTET
		ET1=TET(I).GMET
		NFACE1=ELTTYPE(ET1).NFACE
		TET(I).NF=NFACE1
        !ALLOCATE(TET(I).F(NFACE1))
		DO J=1,NFACE1
            N2=ELTTYPE(ET1).FACE(0,J)
            N1(1:N2)=TET(I).V(ELTTYPE(ET1).FACE(1:N2,J))
            CALL FACE_TBL.KEY(KEY1,N1(1:N2),ELTTYPE(ET1).FACE(0,J))
            CKEY1=TRIM(KEY1)
            HKEY1 = MOD(ABS(HASH_DJB(KEY1)),FACE_TBL%vec_len)            
            VAL1.IEL=I;VAL1.ISE=J;VAL1.IITEM=FACE_TBL.TBL_ID
			CALL FACE_TBL.PUT(TRIM(KEY1),VAL1,HKEY1)
            TET(I).F(J)=VAL1.IITEM
            
            IF(VAL1.IITEM>MAXNFACE) THEN
                CALL FACE_TYDEF_ENLARGE_AR(FACE,1000)
                MAXNFACE=MAXNFACE+1000
            ENDIF
            IF(.NOT.FACE(VAL1.IITEM).ISINI) THEN
                FACE(VAL1.IITEM).CKEY=CKEY1 
                FACE(VAL1.IITEM).HKEY=HKEY1
                FACE(VAL1.IITEM).SHAPE=ELTTYPE(ET1).FACE(0,J)
                FACE(VAL1.IITEM).V(1:FACE(VAL1.IITEM).SHAPE)=TET(I).V(ELTTYPE(ET1).FACE(1:FACE(VAL1.IITEM).SHAPE,J))
				FACE(VAL1.IITEM).EDGE(1:FACE(VAL1.IITEM).SHAPE)=TET(I).E(ABS(ELTTYPE(ET1).FACEEDGE(1:FACE(VAL1.IITEM).SHAPE,J)))
				!CHECK LOOP ORDER
				DO K=1,FACE(VAL1.IITEM).SHAPE
					IF(EDGE(FACE(VAL1.IITEM).EDGE(K)).V(1)/=FACE(VAL1.IITEM).V(K)) FACE(VAL1.IITEM).EDGE(K)=-FACE(VAL1.IITEM).EDGE(K)
				ENDDO
                V1=NODE(FACE(VAL1.IITEM).V(2)).COORD-NODE(FACE(VAL1.IITEM).V(1)).COORD
                V2=NODE(FACE(VAL1.IITEM).V(3)).COORD-NODE(FACE(VAL1.IITEM).V(1)).COORD
                call r8vec_cross_3d ( v1, v2, NORMAL1(1:3) )
                T1 = sqrt ( sum ( ( NORMAL1 )**2 ) )
                if ( T1 /= 0.0D+00 ) then
                  NORMAL1 = NORMAL1/T1
                end if
                FACE(VAL1.IITEM).UNORMAL=NORMAL1
                FACE(VAL1.IITEM).ISINI=.TRUE.
                NFACE=VAL1.IITEM
				DO K=1,3
					FACE(NFACE).BBOX(1,K)=MINVAL(NODE(FACE(NFACE).V(1:FACE(NFACE).SHAPE)).COORD(K))-VTOL
					FACE(NFACE).BBOX(2,K)=MAXVAL(NODE(FACE(NFACE).V(1:FACE(NFACE).SHAPE)).COORD(K))+VTOL
				ENDDO
            ENDIF
            IF(FACE(VAL1.IITEM).ENUM==0) ALLOCATE(FACE(VAL1.IITEM).ELEMENT(2))
            FACE(VAL1.IITEM).ENUM=FACE(VAL1.IITEM).ENUM+1
            IF(FACE(VAL1.IITEM).ENUM>SIZE(FACE(VAL1.IITEM).ELEMENT,DIM=1)) CALL I_ENLARGE_AR(FACE(VAL1.IITEM).ELEMENT,5)
            
            IF(FACE(VAL1.IITEM).ENUM>1) THEN
                DO K=1,FACE(VAL1.IITEM).SHAPE
                    IF(FACE(VAL1.IITEM).V(K)==TET(I).V(ELTTYPE(TET(I).GMET).FACE(1,J))) THEN
                        IF(FACE(VAL1.IITEM).V(MOD(K,FACE(VAL1.IITEM).SHAPE)+1)/= &
                            TET(I).V(ELTTYPE(TET(I).GMET).FACE(2,J))) THEN
                            FACE(VAL1.IITEM).ELEMENT(FACE(VAL1.IITEM).ENUM)=-I !I单元此面的方向与face的方向反向,为-1
                        ELSE
                            FACE(VAL1.IITEM).ELEMENT(FACE(VAL1.IITEM).ENUM)=I
                        ENDIF                  
                        EXIT
                    ENDIF
                ENDDO
            ELSE
                FACE(VAL1.IITEM).ELEMENT(FACE(VAL1.IITEM).ENUM)=I !同向 =1
            ENDIF
            
                
            
		ENDDO
    END DO    
    
    
ENDSUBROUTINE

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
				Elttype(ET).NEDGE=1;Elttype(ET).NFACE=0;Elttype(ET).NTET=0;ELTTYPE(ET).DIM=1
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,1)=[1,2]
			CASE(2) !TRIANGLE
				Elttype(ET).NEDGE=3;Elttype(ET).NFACE=1;Elttype(ET).NTET=0;ELTTYPE(ET).DIM=2
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,1],(/2,3/))
				Elttype(ET).FACE(:,1)=[3,1,2,3,0]
				Elttype(ET).FACEEDGE=Elttype(ET).FACE
			CASE(9) !6-NODED-TRIANGLE
				Elttype(ET).NEDGE=9;Elttype(ET).NFACE=4;Elttype(ET).NTET=0;ELTTYPE(ET).DIM=2
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,4,4,2,2,5,5,3,3,6,6,1,4,5,5,6,6,4],(/2,9/))
				Elttype(ET).FACE(:,:)=RESHAPE([3,1,4,6,0,&
                                               3,4,2,5,0,&
                                               3,5,3,6,0,&
                                               3,4,5,6,0],&
                                               (/5,4/))
				Elttype(ET).FACEEDGE(:,:)=RESHAPE([3,1,-9,6,0,&
                                              3,2,3,-7,0,&
                                              3,4,5,-8,0,&
                                              3,7,8,9,0],&
                                              (/5,4/))
			CASE(23) !15-NODED-TRIANGLE
				Elttype(ET).NEDGE=30;Elttype(ET).NFACE=16;Elttype(ET).NTET=0;ELTTYPE(ET).DIM=2
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,7,7,4,4,15,15,6,6,12,12,1,7,15,15,12,12,7,&
                                               4,8,8,2,2,9,9,5,5,13,13,4,8,9,9,13,13,8,&
                                               5,10,10,3,3,11,11,6,6,14,14,5,10,11,11,14,14,10,&
                                               13,14,14,15,15,13],&
                                               (/2,30/))
				Elttype(ET).FACE(:,:)=RESHAPE([3,1,7,12,0,&
                                               3,7,4,15,0,&
                                               3,15,6,12,0,&
                                               3,7,15,12,0,&
                                               
                                               3,4,8,13,0,&
                                               3,8,2,9,0,&
                                               3,9,5,13,0,&
                                               3,8,9,13,0,&
                                               
                                               3,5,10,14,0,&
                                               3,10,3,11,0,&
                                               3,11,6,14,0,&
                                               3,10,11,14,0,&
                                               
                                               3,4,13,15,0,&
                                               3,13,5,14,0,&
                                               3,14,6,15,0,&
                                               3,13,14,15,0],&
                                               (/5,16/))
                                               
				Elttype(ET).FACEEDGE(:,:)=RESHAPE([3,1,-9,6,0,&
                                                   3,2,3,-7,0,&
                                                   3,4,5,-8,0,&
                                                   3,7,8,9,0,&
                                                   
                                                   3,10,-18,15,0,&
                                                   3,11,12,-16,0,&
                                                   3,13,14,-17,0,&
                                                   3,16,17,18,0,&
                                                   
                                                   3,19,-27,24,0,&
                                                   3,20,21,-25,0,&
                                                   3,22,23,-26,0,&
                                                   3,25,26,27,0,&
                                                   
                                                   3,-15,-30,-3,0,&
                                                   3,-14,-24,-28,0,&
                                                   3,-23,-4,-29,0,&
                                                   3,28,29,30,0],&
                                                    (/5,16/))                                              
			CASE(3) !QUADRANGLE
				Elttype(ET).NEDGE=5;Elttype(ET).NFACE=2;Elttype(ET).NTET=0;ELTTYPE(ET).DIM=2
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,4,4,1,3,1],(/2,5/))                
				Elttype(ET).FACE(:,:)=RESHAPE([3,1,2,3,0,&
                                               3,3,4,1,0],(/5,2/))
				Elttype(ET).FACEEDGE=RESHAPE([3,1,2,5,0,&
                                              3,3,4,-5,0],(/5,2/))
                
			CASE(16) !8-noded-QUADRANGLE
				Elttype(ET).NEDGE=13;Elttype(ET).NFACE=6;Elttype(ET).NTET=0;ELTTYPE(ET).DIM=2
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,5,5,2,2,6,6,3,3,7,7,4,4,8,8,1,&
                                               5,6,6,7,7,8,8,5,8,6], &
                                               (/2,13/))                
				Elttype(ET).FACE(:,:)=RESHAPE([3,1,5,8,0,&
                                               3,5,2,6,0,&
                                               3,6,3,7,0,&
                                               3,7,4,8,0,&
                                               3,5,6,8,0,&
                                               3,6,7,8,0],&
                                               (/5,6/))
				Elttype(ET).FACEEDGE=RESHAPE([3,1,-12,8,0,&
                                              3,2,3,-9,0,&
                                              3,4,5,-10,0,&
                                              3,6,7,-11,0,&
                                              3,9,-13,12,0,&
                                              3,10,11,13,0],&
                                              (/5,6/))                                             
			CASE(4) !TETRAHEDRON
				Elttype(ET).NEDGE=6;Elttype(ET).NFACE=4;Elttype(ET).NTET=1;ELTTYPE(ET).DIM=3
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE),Elttype(ET).TET(4,Elttype(ET).NTET),&
						 Elttype(ET).TETEDGE(6,Elttype(ET).NTET),Elttype(ET).TETFACE(4,Elttype(ET).NTET))	
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
				Elttype(ET).TETEDGE(:,:)=RESHAPE([1,2,3,4,5,6], (/6,1/))
				Elttype(ET).TETFACE(:,:)=Elttype(ET).TET(:,:)
			CASE(11) !10-noded-TETRAHEDRON
				Elttype(ET).NEDGE=25;Elttype(ET).NFACE=24;Elttype(ET).NTET=8;ELTTYPE(ET).DIM=3
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE),Elttype(ET).TET(4,Elttype(ET).NTET),&
						 Elttype(ET).TETEDGE(6,Elttype(ET).NTET),Elttype(ET).TETFACE(4,Elttype(ET).NTET))
				Elttype(ET).EDGE(:,:)=RESHAPE([7,1,&
												1,5,&
												5,7,&
												7,8,&
												1,8,&
												5,8,&
												7,3,&
												3,10,&
												10,7,&
												7,6,&
												3,6,&
												10,6,&
												10,4,&
												4,8,&
												8,10,&
												10,9,&
												4,9,&
												8,9,&
												5,2,&
												2,6,&
												6,5,&
												5,9,&
												2,9,&
												6,9,&
												5,10],&
                                               (/2,25/))
				Elttype(ET).FACE(:,:)=RESHAPE([3,1,7,5,0,&
												3,7,1,8,0,&
												3,1,5,8,0,&
												3,5,7,8,0,&
												3,3,7,10,0,&
												3,7,3,6,0,&
												3,3,10,6,0,&
												3,10,7,6,0,&
												3,4,10,8,0,&
												3,10,4,9,0,&
												3,4,8,9,0,&
												3,8,10,9,0,&
												3,2,5,6,0,&
												3,5,2,9,0,&
												3,2,6,9,0,&
												3,6,5,9,0,&
												3,9,10,5,0,&
												3,10,9,6,0,&
												3,5,10,6,0,&
												3,5,10,7,0,&
												3,5,7,6,0,&
												3,5,9,8,0,&
												3,10,5,8,0,&
												3,10,7,8,0],&
											    (/5,24/))
				Elttype(ET).FACEEDGE(:,:)=RESHAPE([3,-1,-3,-2,0,&
													3,1,5,-4,0,&
													3,2,6,-5,0,&
													3,3,4,-6,0,&
													3,-7,-9,-8,0,&
													3,7,11,-10,0,&
													3,8,12,-11,0,&
													3,9,10,-12,0,&
													3,-13,-15,-14,0,&
													3,13,17,-16,0,&
													3,14,18,-17,0,&
													3,15,16,-18,0,&
													3,-19,-21,-20,0,&
													3,19,23,-22,0,&
													3,20,24,-23,0,&
													3,21,22,-24,0,&
													3,-16,-25,22,0,&
													3,16,-24,-12,0,&
													3,25,12,21,0,&
													3,25,9,-3,0,&
													3,3,10,21,0,&
													3,22,-18,-6,0,&
													3,-25,6,15,0,&
													3,9,4,15,0],&
                                                    (/5,24/))
				Elttype(ET).TET(:,:)=RESHAPE([7,1,5,8,&
												7,3,10,6,&
												10,4,8,9,&
												5,2,6,9,&
												10,9,5,6,&
												10,5,7,6,&
												5,9,10,8,&
												7,5,10,8],&
                                                (/4,8/))
				Elttype(ET).TETEDGE(:,:)=RESHAPE([1,2,3,4,5,6,&
												7,8,9,10,11,12,&
												13,14,15,16,17,18,&
												19,20,21,22,23,24,&
												16,22,25,12,24,21,&
												25,3,9,12,21,10,&
												22,16,25,6,18,15,&
												3,25,9,4,6,15],&
                                                (/6,8/))				
				Elttype(ET).TETFACE(:,:)=RESHAPE([1,2,3,4,&
												5,6,7,8,&
												9,10,11,12,&
												13,14,15,16,&
												17,18,16,19,&
												20,19,21,8,&
												17,22,12,23,&
												20,4,23,24],&
                                                (/4,8/))
			CASE(5) !HEXAHEDRON
				Elttype(ET).NEDGE=19;Elttype(ET).NFACE=18;Elttype(ET).NTET=6;ELTTYPE(ET).DIM=3
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE),Elttype(ET).TET(4,Elttype(ET).NTET),&
						 Elttype(ET).TETEDGE(6,Elttype(ET).NTET),Elttype(ET).TETFACE(4,Elttype(ET).NTET))				
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,&
												2,3,&
												3,1,&
												1,6,&
												2,6,&
												3,6,&
												1,7,&
												7,5,&
												5,1,&
												7,6,&
												5,6,&
												3,7,&
												3,4,&
												4,1,&
												1,8,&
												3,8,&
												4,8,&
												7,8,&
												5,8],(/2,19/))
				Elttype(ET).FACE(:,:)=RESHAPE([3,2,1,3,0,&
												3,1,2,6,0,&
												3,2,3,6,0,&
												3,3,1,6,0,&
												3,7,1,5,0,&
												3,1,7,6,0,&
												3,7,5,6,0,&
												3,5,1,6,0,&
												3,3,1,7,0,&
												3,3,7,6,0,&
												3,3,1,4,0,&
												3,1,3,8,0,&
												3,3,4,8,0,&
												3,4,1,8,0,&
												3,1,7,8,0,&
												3,7,3,8,0,&
												3,1,5,8,0,&
												3,5,7,8,0],&
												(/5,18/))
				Elttype(ET).FACEEDGE(:,:)=RESHAPE([3,-1,-3,-2,0,&
													3,1,5,-4,0,&
													3,2,6,-5,0,&
													3,3,4,-6,0,&
													3,-7,-9,-8,0,&
													3,7,10,-4,0,&
													3,8,11,-10,0,&
													3,9,4,-11,0,&
													3,3,7,-12,0,&
													3,12,10,-6,0,&
													3,3,-14,-13,0,&
													3,-3,16,-15,0,&
													3,13,17,-16,0,&
													3,14,15,-17,0,&
													3,7,18,-15,0,&
													3,-12,16,-18,0,&
													3,-9,19,-15,0,&
													3,-8,18,-19,0],&
													(/5,18/))
				Elttype(ET).TET(:,:)=RESHAPE([1,2,3,6,&
												1,7,5,6,&
												1,3,7,6,&
												1,3,4,8,&
												1,7,3,8,&
												1,5,7,8],&
                                                (/4,6/))
				Elttype(ET).TETEDGE(:,:)=RESHAPE([1,2,3,4,5,6,&
													7,8,9,4,10,11,&
													3,12,7,4,6,10,&
													3,13,14,15,16,17,&
													7,12,3,15,18,16,&
													9,8,7,15,19,18],&
                                                (/6,6/))
				Elttype(ET).TETFACE(:,:)=RESHAPE([1,2,3,4,&
												5,6,7,8,&
												9,4,10,6,&
												11,12,13,14,&
												9,15,16,12,&
												5,17,18,15],&
                                                (/4,6/))				
			CASE(6) !6-NODE PRISM
				Elttype(ET).NEDGE=12;Elttype(ET).NFACE=10;Elttype(ET).NTET=3;ELTTYPE(ET).DIM=3
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE),Elttype(ET).TET(4,Elttype(ET).NTET),&
						 Elttype(ET).TETEDGE(6,Elttype(ET).NTET),Elttype(ET).TETFACE(4,Elttype(ET).NTET))				
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,&
												2,3,&
												3,1,&
												1,5,&
												2,5,&
												3,5,&
												3,6,&
												6,1,&
												6,5,&
												6,4,&
												4,1,&
												4,5],(/2,12/))
				Elttype(ET).FACE(:,:)=RESHAPE([3,2,1,3,0,&
												3,1,2,5,0,&
												3,2,3,5,0,&
												3,3,1,5,0,&
												3,3,1,6,0,&
												3,3,6,5,0,&
												3,6,1,5,0,&
												3,6,1,4,0,&
												3,6,4,5,0,&
												3,4,1,5,0],(/5,10/))
				Elttype(ET).FACEEDGE(:,:)=RESHAPE([3,-1,-3,-2,0,&
												3,1,5,-4,0,&
												3,2,6,-5,0,&
												3,3,4,-6,0,&
												3,3,-8,-7,0,&
												3,7,9,-6,0,&
												3,8,4,-9,0,&
												3,8,-11,-10,0,&
												3,10,12,-9,0,&
												3,11,4,-12,0],(/5,10/))
				Elttype(ET).TET(:,:)=RESHAPE([1,2,3,5,&
												1,3,6,5,&
												1,6,4,5],(/4,3/))
				Elttype(ET).TETEDGE(:,:)=RESHAPE([1,2,3,4,5,6,&
												3,7,8,4,6,9,&
												8,10,11,4,9,12],&
                                                (/6,3/))
				Elttype(ET).TETFACE(:,:)=RESHAPE([1,2,3,4,&
												5,4,6,7,&
												8,7,9,10],&
												(/4,3/))
			CASE(18) !15-NODE PRISM
				Elttype(ET).NEDGE=41;Elttype(ET).NFACE=41;Elttype(ET).NTET=14;ELTTYPE(ET).DIM=3
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE),Elttype(ET).TET(4,Elttype(ET).NTET),&
						 Elttype(ET).TETEDGE(6,Elttype(ET).NTET),Elttype(ET).TETFACE(4,Elttype(ET).NTET))			
				Elttype(ET).EDGE(:,:)=RESHAPE([15,13,&
												13,8,&
												8,15,&
												15,14,&
												13,14,&
												8,14,&
												15,11,&
												11,13,&
												11,14,&
												4,13,&
												13,12,&
												12,4,&
												4,10,&
												13,10,&
												12,10,&
												9,7,&
												7,13,&
												13,9,&
												9,1,&
												7,1,&
												13,1,&
												15,12,&
												12,11,&
												9,15,&
												9,8,&
												7,8,&
												7,14,&
												11,10,&
												10,14,&
												9,3,&
												3,15,&
												3,8,&
												15,6,&
												6,12,&
												6,11,&
												7,2,&
												8,2,&
												14,2,&
												5,10,&
												11,5,&
												5,14],&
											   (/2,41/))
				Elttype(ET).FACE(:,:)=RESHAPE([3,13,15,8,0,&
												3,15,13,14,0,&
												3,13,8,14,0,&
												3,8,15,14,0,&
												3,15,13,11,0,&
												3,15,11,14,0,&
												3,11,13,14,0,&
												3,13,4,12,0,&
												3,4,13,10,0,&
												3,13,12,10,0,&
												3,12,4,10,0,&
												3,7,9,13,0,&
												3,9,7,1,0,&
												3,7,13,1,0,&
												3,13,9,1,0,&
												3,15,13,12,0,&
												3,15,12,11,0,&
												3,12,13,11,0,&
												3,15,9,13,0,&
												3,9,15,8,0,&
												3,13,9,8,0,&
												3,13,8,7,0,&
												3,13,7,14,0,&
												3,7,8,14,0,&
												3,13,10,11,0,&
												3,10,13,14,0,&
												3,11,10,14,0,&
												3,12,10,11,0,&
												3,7,9,8,0,&
												3,9,15,3,0,&
												3,9,3,8,0,&
												3,3,15,8,0,&
												3,6,15,12,0,&
												3,15,6,11,0,&
												3,6,12,11,0,&
												3,7,8,2,0,&
												3,8,14,2,0,&
												3,14,7,2,0,&
												3,10,5,11,0,&
												3,5,10,14,0,&
												3,11,5,14,0],&
											   (/5,41/))
                                               
				Elttype(ET).FACEEDGE(:,:)=RESHAPE([3,-1,-3,-2,0,&
													3,1,5,-4,0,&
													3,2,6,-5,0,&
													3,3,4,-6,0,&
													3,1,-8,-7,0,&
													3,7,9,-4,0,&
													3,8,5,-9,0,&
													3,-10,-12,-11,0,&
													3,10,14,-13,0,&
													3,11,15,-14,0,&
													3,12,13,-15,0,&
													3,-16,-18,-17,0,&
													3,16,20,-19,0,&
													3,17,21,-20,0,&
													3,18,19,-21,0,&
													3,1,11,-22,0,&
													3,22,23,-7,0,&
													3,-11,-8,-23,0,&
													3,-24,-18,-1,0,&
													3,24,-3,-25,0,&
													3,18,25,-2,0,&
													3,2,-26,17,0,&
													3,-17,27,-5,0,&
													3,26,6,-27,0,&
													3,14,-28,8,0,&
													3,-14,5,-29,0,&
													3,28,29,-9,0,&
													3,15,-28,-23,0,&
													3,-16,25,-26,0,&
													3,24,-31,-30,0,&
													3,30,32,-25,0,&
													3,31,-3,-32,0,&
													3,-33,22,-34,0,&
													3,33,35,-7,0,&
													3,34,23,-35,0,&
													3,26,37,-36,0,&
													3,6,38,-37,0,&
													3,-27,36,-38,0,&
													3,-39,-40,28,0,&
													3,39,29,-41,0,&
													3,40,41,-9,0],&
											       (/5,41/))
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
				Elttype(ET).TETEDGE(:,:)=RESHAPE([1,2,3,4,5,6,&
													1,7,8,5,4,9,&
													10,11,12,13,14,15,&
													16,17,18,19,20,21,&
													1,22,11,8,7,23,&
													24,1,18,25,3,2,&
													2,17,26,6,5,27,&
													14,8,28,29,5,9,&
													11,15,14,8,23,28,&
													16,26,25,18,17,2,&
													24,30,31,3,25,32,&
													33,34,22,7,35,23,&
													26,6,27,36,37,38,&
													39,28,40,41,29,9],&
													(/6,14/))
					Elttype(ET).TETFACE(:,:)=RESHAPE([1,2,3,4,&
													5,2,6,7,&
													8,9,10,11,&
													12,13,14,15,&
													16,5,17,18,&
													19,20,1,21,&
													22,3,23,24,&
													25,26,7,27,&
													10,18,28,25,&
													29,12,22,21,&
													30,20,31,32,&
													33,34,35,17,&
													24,36,37,38,&
													39,40,27,41],&
													(/4,14/))
			CASE(7) !5-node pyramid
				Elttype(ET).NEDGE=9;Elttype(ET).NFACE=7;Elttype(ET).NTET=2;ELTTYPE(ET).DIM=3
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE),Elttype(ET).TET(4,Elttype(ET).NTET),&
						 Elttype(ET).TETEDGE(6,Elttype(ET).NTET),Elttype(ET).TETFACE(4,Elttype(ET).NTET))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,4,4,1,&
											   1,5,2,5,3,5,4,5,&
											   1,3],(/2,9/))
				Elttype(ET).FACE(:,:)=RESHAPE([3,2,1,3,0,&
											   3,1,4,3,0,&
											   3,1,2,5,0,&
											   3,2,3,5,0,&
											   3,3,4,5,0,&
											   3,4,1,5,0,&
											   3,1,3,5,0],(/5,7/))
				Elttype(ET).FACEEDGE(:,:)=RESHAPE([3,-1,9,-2,0,&
												   3,-4,-3,-9,0,&
												   3,1,6,-5,0,&
												   3,2,7,-6,0,&
												   3,3,8,-7,0,&
												   3,4,5,-8,0,&
												   3,9,7,-5],(/5,7/))
				Elttype(ET).TET(:,:)=RESHAPE([1,2,3,5,&
												1,3,4,5],&
                                                (/4,2/))
				Elttype(ET).TETEDGE(:,:)=RESHAPE([1,2,9,5,6,7,&
												9,3,4,5,7,8],&
                                                (/6,2/))
				Elttype(ET).TETFACE(:,:)=RESHAPE([1,3,4,7,&
												  2,5,6,7],&
												(/4,2/))				
				
			CASE(15) !POINTS
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

INTEGER FUNCTION POINTlOC(PT)
    USE solverds
	USE MESHGEO
	IMPLICIT NONE
	REAL(8),INTENT(IN)::PT(3)
	INTEGER::I
	LOGICAL,EXTERNAL::PtInTri,PtInTET
	
	!I=INT(1+(RANDOM(1))(NFACE-1))

	do i=1,NTET
		
		IF(TET(I).BBOX(2,1)>=PT(1).and.TET(I).BBOX(1,1)<=PT(1)) then
			IF(TET(I).BBOX(2,2)>=PT(2).and.TET(I).BBOX(1,2)<=PT(2)) then
				IF(NDIMENSION>2) THEN
					IF(TET(I).BBOX(2,3)<PT(3).and.TET(I).BBOX(1,3)>PT(3)) CYCLE
				ENDIF
				IF(TET(I).DIM==2) THEN
					IF(PtInTri(PT, NODE(TET(I).V(1)).COORD, NODE(TET(I).V(2)).COORD, NODE(TET(I).V(3)).COORD)) EXIT
				ELSEIF(TET(I).DIM==3) THEN
					IF(PtInTET(PT, NODE(TET(I).V(1)).COORD, NODE(TET(I).V(2)).COORD, NODE(TET(I).V(3)).COORD,NODE(TET(I).V(4)).COORD)) EXIT
				ENDIF
			ENDIF
		endif
	ENDDO
	
	IF(I>NTET) I=0
	POINTlOC=I
	RETURN
ENDFUNCTION

INTEGER FUNCTION PTINTRIlOC(PT)
    USE solverds
	USE MESHGEO
	IMPLICIT NONE
	REAL(8),INTENT(IN)::PT(3)
	INTEGER::I
	LOGICAL,EXTERNAL::PtInTri,PtInTET
    REAL(8)::TOF1
	
	!I=INT(1+(RANDOM(1))(NFACE-1))
    
	do i=1,NFACE
		
		IF(FACE(I).BBOX(2,1)>=PT(1).and.FACE(I).BBOX(1,1)<=PT(1)) then
			IF(FACE(I).BBOX(2,2)>=PT(2).and.FACE(I).BBOX(1,2)<=PT(2)) then
				IF(NDIMENSION>2) THEN
					IF(FACE(I).BBOX(2,3)<PT(3).and.FACE(I).BBOX(1,3)>PT(3)) CYCLE
				ENDIF
				IF(PtInTri(PT, NODE(FACE(I).V(1)).COORD, NODE(FACE(I).V(2)).COORD, NODE(FACE(I).V(3)).COORD)) EXIT
				
			ENDIF
		endif
	ENDDO
	
	IF(I>NFACE) I=0
	PTINTRIlOC=I
	RETURN
ENDFUNCTION

logical function PtInTri (pt, v1, v2, v3)
	implicit none
	integer b1, b2
	real(8),intent(in)::pt(3),v1(3),v2(3),v3(3)
    integer,EXTERNAL::isacw
    
	PtInTri=.false.
    b1 = isacw(pt(1),pt(2),pt(3),v1(1),v1(2),v1(3),v2(1),v2(2),v2(3)) ;
    if(b1==2) then
        PtInTri=.true.
        return
    endif
    b2 = isacw(pt(1),pt(2),pt(3),v2(1),v2(2),v2(3),v3(1),v3(2),v3(3)) ;
    if(b2==2) then
        PtInTri=.true.
        return
    endif    
	if(b1/=b2) return
    b1 = isacw(pt(1),pt(2),pt(3),v3(1),v3(2),v3(3),v1(1),v1(2),v1(3)) ;
    if(b1==2) then
        PtInTri=.true.
        return
    endif    
    if(b1/=b2) return
    
	PtInTri=.true.	
end function

logical function PtInTet(pt, v1, v2, v3,v4)
    !假定四面体四个面的正向一致，都为正向(节点逆时针)或负向
	implicit none
	integer:: b1, b2
	real(8),intent(in)::pt(3),v1(3),v2(3),v3(3),v4(3)
	integer,EXTERNAL::ISFRONT
    
	PtInTet=.false.
    b1 = Isfront([v2,v1,v3,pt]) 
    if(b1==2) then
        PtInTet=.true.
        return
    endif
    b2 = Isfront([v1,v2,v4,pt]) 
    if(b2==2) then
        PtInTet=.true.
        return
    endif    
	if(b1/=b2) return
    b1 = Isfront([v2,v3,v4,pt])
    if(b1==2) then
        PtInTet=.true.
        return
    endif
	if(b1/=b2) return	
	b1 = Isfront([v3,v1,v4,pt])
    if(b1==2) then
        PtInTet=.true.
        return
    endif    
	if(b1/=b2) return
    
    PtInTet=.true.

end function

integer function Isfront(V)
!v(:,1-3) are face, and v(:,4) is point to be tested
    
	implicit none
	real(8),intent(in)::V(3,4)
	real(8)::V1(3,3)
	integer::I
    real(8)::T1
    logical,external::PtInTri
    real(8),EXTERNAL::determinant
    
    Isfront=0
	do i=1,3
		v1(:,i)=v(:,i+1)-v(:,1)
	enddo

	if(determinant(v1)>0.d0) Isfront=1
    
    if(abs(t1)<1e-10) then
        if(PtInTri (v(:,4), v(:,1), v(:,2), v(:,3))) Isfront=2 !on the surface        
    endif
    

end function

integer function isacw(x1,y1,z1,x2,y2,z2,x3,y3,z3)
	implicit none
	real(8),intent(in)::x1,y1,z1,x2,y2,z2,x3,y3,z3
    real(8)::t1,yn2,xn2,zn2,yn3,xn3,zn3
	
	isacw=0
    
	yn2=y2-y1
	xn2=x2-x1
    zn2=z2-z1
	yn3=y3-y1
	xn3=x3-x1
    zn3=z3-z1
	t1=(xn2*yn3-yn2*xn3)+(yn2*zn3-zn2*yn3)-(xn2*zn3-zn2*xn3)
    if(t1>0.d0) isacw=1
    
    if(abs(t1)<1e-10) then
        if(t1<min(x1,x2,x3)) return
        if(t1>max(x1,x2,x3)) return
        if(t1<min(y1,y2,y3)) return
        if(t1>max(y1,y2,y3)) return 
        if(t1<min(z1,z2,z3)) return
        if(t1>max(z1,z2,z3)) return 
        isacw=2 !on the edge
    endif
	

end function



