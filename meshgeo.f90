MODULE MESHGEO

USE solverds
USE hashtbl

PRIVATE
PUBLIC::EDGE,NEDGE,FACE,NFACE,SETUP_EDGE_TBL,SETUP_FACE_TBL

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
    INTEGER::V(4)=0,EDGE(4)=0
    INTEGER::HKEY=-1
    CHARACTER(64)::CKEY="" 		
    INTEGER::ISTRISURFACE=0 !<-1 与face反向
    REAL(8)::UNORMAL(3)=0.D0		
    INTEGER::ENUM=0
    INTEGER,ALLOCATABLE::ELEMENT(:),SUBID(:) !SUBID IS INDEX WHICH FACE OF THE ELEMENT IS THE FACE       
ENDTYPE
TYPE(FACE_TYDEF),ALLOCATABLE::FACE(:)    
INTEGER::NEDGE=0,NFACE=0,MAXNEDGE=10000,MAXNFACE=10000

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
    
CONTAINS 



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

    
SUBROUTINE SETUP_EDGE_TBL()
    
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
    
	TBL_LEN=NNUM	
    CALL EDGE_TBL.INIT(TBL_LEN)
	DO I=1,ENUM
		ET1=GETGMSHET(ELEMENT(I).ET)
		NEDGE1=ELTTYPE(ET1).NEDGE
        ALLOCATE(ELEMENT(I).EDGE(NEDGE1))
		DO J=1,NEDGE1
            N1=ELEMENT(I).NODE(ELTTYPE(ET1).EDGE(:,J))
            CALL I2C_KEY_hash_tbl_sll(KEY1,N1(:),2)
            CKEY1=TRIM(KEY1)
            HKEY1 = MOD(ABS(HASH_DJB(KEY1)),EDGE_TBL%vec_len)
            VAL1.IEL=I;VAL1.ISE=J;VAL1.IITEM=EDGE_TBL.TBL_ID
			CALL EDGE_TBL.PUT(TRIM(KEY1),VAL1,HKEY1)
			ELEMENT(I).EDGE(J)=VAL1.IITEM
            IF(VAL1.IITEM>MAXNEDGE) THEN
                CALL EDGE_TYDEF_ENLARGE_AR(EDGE,1000)
                MAXNEDGE=MAXNEDGE+1000
            ENDIF
            IF(.NOT.EDGE(VAL1.IITEM).ISINI) THEN
                EDGE(VAL1.IITEM).CKEY=CKEY1
                EDGE(VAL1.IITEM).HKEY=HKEY1
                EDGE(VAL1.IITEM).V=ELEMENT(I).NODE(ELTTYPE(ET1).EDGE(:,J))
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

SUBROUTINE SETUP_FACE_TBL()
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
    
	TBL_LEN=NNUM	
    CALL FACE_TBL.INIT(TBL_LEN)
	DO I=1,ENUM
		ET1=GETGMSHET(ELEMENT(I).ET)
		NFACE1=ELTTYPE(ET1).NFACE
        ALLOCATE(ELEMENT(I).FACE(NFACE1))
		DO J=1,NFACE1
            N2=ELTTYPE(ET1).FACE(0,J)
            N1(1:N2)=ELEMENT(I).NODE(ELTTYPE(ET1).FACE(1:N2,J))
            CALL FACE_TBL.KEY(KEY1,N1(1:N2),ELTTYPE(ET1).FACE(0,J))
            CKEY1=TRIM(KEY1)
            HKEY1 = MOD(ABS(HASH_DJB(KEY1)),FACE_TBL%vec_len)            
            VAL1.IEL=I;VAL1.ISE=J;VAL1.IITEM=FACE_TBL.TBL_ID
			CALL FACE_TBL.PUT(TRIM(KEY1),VAL1,HKEY1)
            ELEMENT(I).FACE(J)=VAL1.IITEM
            
            IF(VAL1.IITEM>MAXNFACE) THEN
                CALL FACE_TYDEF_ENLARGE_AR(FACE,1000)
                MAXNFACE=MAXNFACE+1000
            ENDIF
            IF(.NOT.FACE(VAL1.IITEM).ISINI) THEN
                FACE(VAL1.IITEM).CKEY=CKEY1 
                FACE(VAL1.IITEM).HKEY=HKEY1
                FACE(VAL1.IITEM).SHAPE=ELTTYPE(ET1).FACE(0,J)
                FACE(VAL1.IITEM).V(1:FACE(VAL1.IITEM).SHAPE)=ELEMENT(I).NODE(ELTTYPE(ET1).FACE(1:FACE(VAL1.IITEM).SHAPE,J))
				FACE(VAL1.IITEM).EDGE(1:FACE(VAL1.IITEM).SHAPE)=ELEMENT(I).EDGE(ABS(ELTTYPE(ET1).FACEEDGE(1:FACE(VAL1.IITEM).SHAPE,J)))
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
            ENDIF
            IF(FACE(VAL1.IITEM).ENUM==0) ALLOCATE(FACE(VAL1.IITEM).ELEMENT(2))
            FACE(VAL1.IITEM).ENUM=FACE(VAL1.IITEM).ENUM+1
            IF(FACE(VAL1.IITEM).ENUM>SIZE(FACE(VAL1.IITEM).ELEMENT,DIM=1)) CALL I_ENLARGE_AR(FACE(VAL1.IITEM).ELEMENT,5)
            
            IF(FACE(VAL1.IITEM).ENUM>1) THEN
                DO K=1,FACE(VAL1.IITEM).SHAPE
                    IF(FACE(VAL1.IITEM).V(K)==ELEMENT(I).NODE(ELTTYPE(GETGMSHET(ELEMENT(I).ET)).FACE(1,J))) THEN
                        IF(FACE(VAL1.IITEM).V(MOD(K,FACE(VAL1.IITEM).SHAPE)+1)/= &
                            ELEMENT(I).NODE(ELTTYPE(GETGMSHET(ELEMENT(I).ET)).FACE(2,J))) THEN
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

SUBROUTINE ET_GMSH_EDGE_FACE()

	IMPLICIT NONE
	INTEGER::I,ET,AET1(33)
	
    AET1=[1:31,92,93]
	DO I=1,33
        ET=AET1(I)
		SELECT CASE(ET)
			CASE(1) !LINE
				Elttype(ET).NEDGE=1;Elttype(ET).NFACE=0
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,1)=[1,2]
			CASE(2) !TRIANGLE
				Elttype(ET).NEDGE=3;Elttype(ET).NFACE=1
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,1],(/2,3/))
				Elttype(ET).FACE(:,1)=[3,1,2,3,0]
				Elttype(ET).FACEEDGE=Elttype(ET).FACE
			CASE(9) !6-NODED-TRIANGLE
				Elttype(ET).NEDGE=9;Elttype(ET).NFACE=4
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
				Elttype(ET).NEDGE=30;Elttype(ET).NFACE=16
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
				Elttype(ET).NEDGE=4;Elttype(ET).NFACE=1
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,4,4,1],(/2,4/))                
				Elttype(ET).FACE(:,:)=RESHAPE([4,1,2,3,4],(/5,1/))
				Elttype(ET).FACEEDGE=RESHAPE([4,1,2,3,4],(/5,1/))
                
			CASE(16) !8-noded-QUADRANGLE
				Elttype(ET).NEDGE=12;Elttype(ET).NFACE=5
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,5,5,2,2,6,6,3,3,7,7,4,4,8,8,1,&
                                               5,6,6,7,7,8,8,5], &
                                               (/2,12/))                
				Elttype(ET).FACE(:,:)=RESHAPE([3,1,5,8,0,&
                                               3,5,2,6,0,&
                                               3,6,3,7,0,&
                                               3,7,4,8,0,&
                                               4,5,6,7,8],&
                                               (/5,5/))
				Elttype(ET).FACEEDGE=RESHAPE([3,1,-12,8,0,&
                                              3,2,3,-9,0,&
                                              3,4,5,-10,0,&
                                              3,6,7,-11,0,&
                                              4,9,10,11,12],&
                                              (/5,5/))                                             
			CASE(4) !TETRAHEDRON
				Elttype(ET).NEDGE=6;Elttype(ET).NFACE=4
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,1,1,4,2,4,3,4],(/2,6/))
				Elttype(ET).FACE(:,:)=RESHAPE([3,2,1,3,0,&
											   3,1,2,4,0,&
											   3,2,3,4,0,&
											   3,3,1,4,0],(/5,4/))
				Elttype(ET).FACEEDGE(:,:)=RESHAPE([3,-1,-3,-2,0,&
											   3,1,5,-4,0,&
											   3,2,6,-5,0,&
											   3,3,4,-6,0],(/5,4/))
			CASE(11) !10-noded-TETRAHEDRON
				Elttype(ET).NEDGE=24;Elttype(ET).NFACE=16
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,7,7,3,3,6,6,2,2,5,5,1,&
                                               1,8,8,4,&
                                               2,9,9,4,&
                                               3,10,10,4,&
                                               5,7,7,6,6,5,&
                                               5,9,9,8,8,5,&
                                               6,10,10,9,9,6,&
                                               7,8,8,10,10,7],&
                                               (/2,24/))
				Elttype(ET).FACE(:,:)=RESHAPE([3,1,7,5,0,&
											   3,2,5,6,0,&
											   3,3,6,7,0,&
											   3,5,7,6,0,&
                                               
                                               3,1,5,8,0,&
											   3,2,9,5,0,&
											   3,4,8,9,0,&
											   3,5,9,8,0,&
                                               
                                               3,2,6,9,0,&
											   3,3,10,6,0,&
											   3,4,9,10,0,&
											   3,6,10,9,0,&
                                               
                                               3,1,8,7,0,&
											   3,3,7,10,0,&
											   3,4,10,8,0,&
											   3,7,8,10,0],&
                                               (/5,16/))
				Elttype(ET).FACEEDGE(:,:)=RESHAPE([3,1,-13,6,0,&
											       3,5,-15,4,0,&
											       3,3,-14,2,0,&
											       3,13,14,15,0,&
                                               
                                                   3,-6,-18,-7,0,&
											       3,9,-16,-5,0,&
											       3,-8,-17,10,0,&
											       3,16,17,18,0,&
                                               
                                                   3,-4,-21,-9,0,&
											       3,11,-19,-13,0,&
											       3,-10,-20,12,0,&
											       3,19,20,21,0,&
                                               
                                                   3,7,-22,-1,0,&
											       3,-2,-24,-11,0,&
											       3,-12,-23,8,0,&
											       3,22,23,24,0],&
                                                   (/5,16/))                                             
			CASE(5) !HEXAHEDRON
				Elttype(ET).NEDGE=12;Elttype(ET).NFACE=6
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE))				
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
											   4,2,11,-6,-10,&
											   4,3,12,-7,-11,&
											   4,4,9,-8,-12],(/5,6/))			
			CASE(6) !6-NODE PRISM
				Elttype(ET).NEDGE=9;Elttype(ET).NFACE=5
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE))				
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
			CASE(18) !15-NODE PRISM
				Elttype(ET).NEDGE=39;Elttype(ET).NFACE=26
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE))				
				Elttype(ET).EDGE(:,:)=RESHAPE([1,9,9,3,3,8,8,2,2,7,7,1, &
											   4,10,10,5,5,11,11,6,6,12,12,4,&
                                               1,13,13,4,&
                                               2,14,14,5,&
                                               3,15,15,6,&
                                               7,9,9,8,8,7,&
                                               10,11,11,12,12,10,&
                                               7,14,14,10,10,13,13,7,13,14,&
                                               8,15,15,11,11,14,14,8,14,15,&
                                               9,13,13,12,12,15,15,9,15,13],&
											   (/2,39/))
				Elttype(ET).FACE(:,:)=RESHAPE([3,1,9,7,0,&
											   3,2,7,8,0,&
											   3,3,8,9,0,&
											   3,7,9,8,0,&
                                               
                                               3,4,10,12,0,&
											   3,5,11,10,0,&
											   3,6,12,11,0,&
											   3,10,11,12,0,&
                                               
                                               3,1,7,13,0,&
											   3,2,14,7,0,&
											   3,5,10,14,0,&
											   3,4,13,10,0,&
                                               3,7,14,13,0,&
											   3,10,13,14,0,&
                                               
                                               3,2,8,14,0,&
											   3,3,15,8,0,&
											   3,6,11,15,0,&
											   3,5,14,11,0,&
                                               3,8,15,14,0,&
											   3,11,14,15,0,&                                               
                                               
                                               3,3,9,15,0,&
											   3,1,13,9,0,&
											   3,4,12,13,0,&
											   3,6,15,12,0,&
                                               3,9,13,15,0,&
											   3,12,15,13,0],&                                               
											   (/5,26/))
                                               
				Elttype(ET).FACEEDGE(:,:)=RESHAPE([3,1,-19,6,0,&
											       3,5,-21,4,0,&
											       3,3,-20,2,0,&
											       3,19,20,21,0,&
                                               
                                                   3,7,-24,12,0,&
											       3,9,-22,8,0,&
											       3,11,-23,10,0,&
											       3,22,23,24,0,&
                                               
                                                   3,-6,-28,-13,0,&
											       3,15,-25,-5,0,&
											       3,-8,-26,16,0,&
											       3,-14,-27,-7,0,&
                                                   3,25,-29,28,0,&
											       3,27,29,26,0,&
                                               
                                                   3,-4,-33,-15,0,&
											       3,17,-30,-3,0,&
											       3,-10,-21,18,0,&
											       3,-16,-32,-9,0,&
                                                   3,30,-34,33,0,&
											       3,32,34,31,0,&                                               
                                               
                                                   3,-2,-38,-17,0,&
											       3,13,-35,-1,0,&
											       3,-12,-36,14,0,&
											       3,-18,-37,-11,0,&
                                                   3,35,-39,38,0,&
											       3,37,39,36,0],&                                               
											       (/5,26/))		
			CASE(7) !5-node pyramid
				Elttype(ET).NEDGE=8;Elttype(ET).NFACE=5
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE))
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
			CASE(15) !POINTS
				Elttype(ET).NEDGE=0;Elttype(ET).NFACE=0
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