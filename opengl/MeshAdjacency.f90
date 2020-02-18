MODULE MESHADJ
USE solverds,ONLY:MAX_NODE_ADJ,MAX_FACE_ADJ   
USE MESHGEO
USE quicksort

INTERFACE SETUP_EDGE_ADJL
    MODULE PROCEDURE SETUP_EDGE_ADJL_MODEL,SETUP_EDGE_ADJL_TET,SETUP_EDGE_ADJL_SOLVER
END INTERFACE
INTERFACE SETUP_FACE_ADJL
    MODULE PROCEDURE SETUP_FACE_ADJL_MODEL,SETUP_FACE_ADJL_TET,SETUP_FACE_ADJL_SOLVER
END INTERFACE

PRIVATE::MAX_NODE_ADJ,MAX_FACE_ADJ
    
CONTAINS

SUBROUTINE SETUP_EDGE_BC()
    IMPLICIT NONE
    INTEGER::I,J,N1,N2,N3
    INTEGER,ALLOCATABLE::IA1(:,:)
    
    NE_BC=COUNT(EDGE.ENUM==1)
    ALLOCATE(EDGE_BC(NE_BC),NODE_LOOP_BC(NE_BC+1),IA1(2,POSDATA.NNODE))
    NE_BC=0;IA1=0
    DO I=1,NEDGE
        IF(EDGE(I).ENUM==1) THEN
            NE_BC=NE_BC+1
            EDGE_BC(NE_BC)=I
            DO J=1,2
                N1=MINLOC(IA1(:,EDGE(I).V(J)),MASK=IA1(:,EDGE(I).V(J))==0,DIM=1)
                IF(N1==0) THEN
                    STOP "SIMPLE LOOP IS ASSUMED"
                ENDIF
                IA1(N1,EDGE(I).V(J))=EDGE(I).V(MOD(J,2)+1)
            ENDDO
        ENDIF
    ENDDO
    N1=EDGE(EDGE_BC(1)).V(1)
    N2=EDGE(EDGE_BC(1)).V(2)
    NODE_LOOP_BC(1:2)=[N1,N2]
    DO I=2,NE_BC
        N2=NODE_LOOP_BC(I);N1=NODE_LOOP_BC(I-1)        
        N3=MINLOC(IA1(:,N2)-N1,MASK=ABS(IA1(:,N2)-N1)>0,DIM=1)
        NODE_LOOP_BC(I+1)=IA1(N3,N2)        
    ENDDO
    NNLBC=NE_BC
    !CHECK    
    !IF(NODE_LOOP_BC(1)/=NODE_LOOP_BC(NNLBC)) THEN    
    !    STOP "ERROR.THE FIRST AND THE LAST NODE SHOULD BE IDENTICAL.SUB=SETUP_EDGE_BC"
    !ENDIF
    
    
    DEALLOCATE(IA1)
    
  
    
END SUBROUTINE


!SET EDGE FOR SOLVER MODEL ELEMENT. 
SUBROUTINE SETUP_EDGE_ADJL_SOLVER(EDGE_L,NEDGE_L,NADJL)
    USE solverds,ONLY:ELEMENT_L=>ELEMENT,NODE_L=>NODE,NNODE_L=>NNUM,NEL_L=>ENUM
    IMPLICIT NONE
    INTEGER::NEDGE_L
	TYPE(EDGE_TYDEF),ALLOCATABLE::EDGE_L(:)
    TYPE(NODE_ADJ_TYDEF),ALLOCATABLE::NADJL(:)
    INTEGER::MAXADJL1
    INTEGER::ADJL1(MAX_NODE_ADJ,NNODE_L),E_ADJL1(MAX_NODE_ADJ,NNODE_L),NADJL1(NNODE_L) 
    INTEGER::VMAX1,VMIN1,IEDGE1,ET1,NEDGE1,N1,N2,V1(2)
    INTEGER::I,J    
  
    MAXADJL1=MAX_NODE_ADJ
    
    IF(.NOT.ISINI_GMSHET) THEN
        
        CALL Initialize_et2numNodes()
        CALL ET_GMSH_EDGE_FACE()
        ISINI_GMSHET=.TRUE.
    
    ENDIF
    
    
	IF(ALLOCATED(EDGE_L)) DEALLOCATE(EDGE_L)
	ALLOCATE(EDGE_L(NNODE_L))
    
    IF(ALLOCATED(NADJL)) DEALLOCATE(NADJL)
    ALLOCATE(NADJL(NNODE_L))
    
    ADJL1=0;E_ADJL1=0;NADJL1=0
	DO I=1,NEL_L
		ET1=GETGMSHET(ELEMENT_L(I).ET)
		NEDGE1=ELTTYPE(ET1).NEDGE
        ELEMENT_L(I).NEDGE=NEDGE1
        ALLOCATE(ELEMENT_L(I).EDGE(NEDGE1))
		DO J=1,NEDGE1
            V1=ELEMENT_L(I).NODE(ELTTYPE(ET1).EDGE(:,J))
            IF(V1(1)>V1(2)) THEN
                VMAX1=V1(1);VMIN1=V1(2)
            ELSE
                VMAX1=V1(2);VMIN1=V1(1)
            ENDIF
 
            N2=NADJL1(VMAX1)
            
            N1=MINLOC(ABS(ADJL1(1:N2,VMAX1)-VMIN1),MASK=ADJL1(1:N2,VMAX1)-VMIN1==0,DIM=1)
            if(N1>0) THEN
                IEDGE1=E_ADJL1(N1,VMAX1)            
            ELSE
                NADJL1(VMAX1)=N2+1
                IF(N2+1>MAXADJL1) THEN
                     print *, 'MAX_NODE_ADJ NEEDS TO BE  ENLARGED BY APPENDING "MAX_NODE_ADJ=XXX" TO THE KEYWORD "SolverControl".MAX_NODE_ADJ=',MAX_NODE_ADJ
                     stop
                ENDIF             
                ADJL1(N2+1,VMAX1)=VMIN1
                NEDGE_L=NEDGE_L+1
                IEDGE1=NEDGE_L
                E_ADJL1(N2+1,VMAX1)=NEDGE_L                
            ENDIF
     

            ELEMENT_L(I).EDGE(J)=IEDGE1
            IF(IEDGE1>SIZE(EDGE_L,DIM=1)) THEN
                CALL EDGE_TYDEF_ENLARGE_AR(EDGE_L,1000)
                !MAXNEDGE=MAXNEDGE+1000
            ENDIF
            IF(.NOT.EDGE_L(IEDGE1).ISINI) THEN
                EDGE_L(IEDGE1).V=ELEMENT_L(I).NODE(ELTTYPE(ET1).EDGE(:,J))
                EDGE_L(IEDGE1).DIS=NORM2(NODE_L(EDGE_L(IEDGE1).V(1)).COORD-NODE_L(EDGE_L(IEDGE1).V(2)).COORD)
                EDGE_L(IEDGE1).ISINI=.TRUE.
                NEDGE_L=IEDGE1
                
                NADJL(EDGE_L(IEDGE1).V).NNUM=NADJL(EDGE_L(IEDGE1).V).NNUM+1
            ENDIF
            IF(EDGE_L(IEDGE1).ENUM==0) ALLOCATE(EDGE_L(IEDGE1).ELEMENT(5),EDGE_L(IEDGE1).SUBID(5))
            EDGE_L(IEDGE1).ENUM=EDGE_L(IEDGE1).ENUM+1
            IF(EDGE_L(IEDGE1).ENUM>SIZE(EDGE_L(IEDGE1).ELEMENT,DIM=1)) THEN
				CALL I_ENLARGE_AR(EDGE_L(IEDGE1).ELEMENT,5)
				CALL I_ENLARGE_AR(EDGE_L(IEDGE1).SUBID,5)
			ENDIF
            EDGE_L(IEDGE1).ELEMENT(EDGE_L(IEDGE1).ENUM)=I			
            EDGE_L(IEDGE1).SUBID(EDGE_L(IEDGE1).ENUM)=J
		ENDDO
        
        NADJL(ELEMENT_L(I).NODE).ENUM=NADJL(ELEMENT_L(I).NODE).ENUM+1
    END DO
    
    DO I=1, NNODE_L
        ALLOCATE(NADJL(I).NODE(NADJL(I).NNUM),NADJL(I).ELEMENT(NADJL(I).ENUM),NADJL(I).SUBID(NADJL(I).ENUM))
        NADJL(I).ENUM=0;NADJL(I).NNUM=0
    ENDDO
    
    DO I=1,NEL_L
        NADJL(ELEMENT_L(I).NODE).ENUM=NADJL(ELEMENT_L(I).NODE).ENUM+1
        DO J=1,ELEMENT_L(I).NNUM
            N1=ELEMENT_L(I).NODE(J)
            NADJL(N1).ELEMENT(NADJL(N1).ENUM)=I
            NADJL(N1).SUBID(NADJL(N1).ENUM)=J
        ENDDO
    ENDDO
    
     DO I=1,NEDGE_L
        NADJL(EDGE_L(I).V).NNUM=NADJL(EDGE_L(I).V).NNUM+1
        DO J=1,2
            N1=EDGE_L(I).V(J)
            NADJL(N1).NODE(NADJL(N1).NNUM)=EDGE_L(I).V(MOD(J,2)+1)            
        ENDDO
    ENDDO   
    
    RETURN
ENDSUBROUTINE

!SET EDGE FOR MODEL ELEMENT. 
SUBROUTINE SETUP_EDGE_ADJL_MODEL(EDGE_L,NEDGE_L,ELEMENT_L,NEL_L,NODE_L,NNODE_L,ESET,NESET)
    
    IMPLICIT NONE
	INTEGER,INTENT(IN)::NEL_L,NNODE_L,NESET
    INTEGER::NEDGE_L
	TYPE(EDGE_TYDEF),ALLOCATABLE::EDGE_L(:)
	type(ELEMENT_TYDEF)::ELEMENT_L(NEL_L)
	REAL(8),INTENT(IN)::NODE_L(3,NNODE_L)
    TYPE(ESET_TYDEF),INTENT(IN)::ESET(NESET)
    INTEGER::MAXADJL1
    INTEGER::ADJL1(MAX_NODE_ADJ,NNODE_L),E_ADJL1(MAX_NODE_ADJ,NNODE_L),NADJL1(NNODE_L) 
    INTEGER::VMAX1,VMIN1,IEDGE1,ET1,NEDGE1,N1,N2,V1(2)
    INTEGER::I,J
    
    
    MAXADJL1=MAX_NODE_ADJ
    
	IF(ALLOCATED(EDGE_L)) DEALLOCATE(EDGE_L)
	ALLOCATE(EDGE_L(NNODE_L))
    
    ADJL1=0;E_ADJL1=0;NADJL1=0
	DO I=1,NEL_L
		ET1=GETGMSHET(ESET(ELEMENT_L(I).ISET).ET)
		NEDGE1=ELTTYPE(ET1).NEDGE
        ELEMENT_L(I).NEDGE=NEDGE1
        ALLOCATE(ELEMENT_L(I).EDGE(NEDGE1))
		DO J=1,NEDGE1
            V1=ELEMENT_L(I).NODE(ELTTYPE(ET1).EDGE(:,J))
            IF(V1(1)>V1(2)) THEN
                VMAX1=V1(1);VMIN1=V1(2)
            ELSE
                VMAX1=V1(2);VMIN1=V1(1)
            ENDIF
            !VMAX1=MAXVAL(ELEMENT_L(I).V(ELTTYPE(ET1).EDGE(:,J)))    
            !VMIN1=MINVAL(ELEMENT_L(I).V(ELTTYPE(ET1).EDGE(:,J)))   
            N2=NADJL1(VMAX1)
            !N2=MINLOC(ADJL1(:,VMAX1),MASK=ADJL1(:,VMAX1)==0,DIM=1)-1
            !IF(N2<0) THEN
            !     STOP 'ERROR IN SUB SETUP_EDGE_TBLI2. THE SIZE OF ADJL1 MAY BE TOO SMALL.'
            !ENDIF
            
            N1=MINLOC(ABS(ADJL1(1:N2,VMAX1)-VMIN1),MASK=ADJL1(1:N2,VMAX1)-VMIN1==0,DIM=1)
            if(N1>0) THEN
                IEDGE1=E_ADJL1(N1,VMAX1)            
            ELSE
                NADJL1(VMAX1)=N2+1
                IF(N2+1>MAXADJL1) THEN
                     print *, 'MAX_NODE_ADJ NEEDS TO BE  ENLARGED BY APPENDING "MAX_NODE_ADJ=XXX" TO THE KEYWORD "SolverControl".MAX_NODE_ADJ=',MAX_NODE_ADJ
                     stop
                ENDIF             
                ADJL1(N2+1,VMAX1)=VMIN1
                NEDGE_L=NEDGE_L+1
                IEDGE1=NEDGE_L
                E_ADJL1(N2+1,VMAX1)=NEDGE_L                
            ENDIF
     

            ELEMENT_L(I).EDGE(J)=IEDGE1
            IF(IEDGE1>SIZE(EDGE_L,DIM=1)) THEN
                CALL EDGE_TYDEF_ENLARGE_AR(EDGE_L,1000)
                !MAXNEDGE=MAXNEDGE+1000
            ENDIF
            IF(.NOT.EDGE_L(IEDGE1).ISINI) THEN
                EDGE_L(IEDGE1).V=ELEMENT_L(I).NODE(ELTTYPE(ET1).EDGE(:,J))
                EDGE_L(IEDGE1).DIS=NORM2(NODE_L(:,EDGE_L(IEDGE1).V(1))-NODE_L(:,EDGE_L(IEDGE1).V(2)))
                EDGE_L(IEDGE1).ISINI=.TRUE.
                NEDGE_L=IEDGE1
            ENDIF
            IF(EDGE_L(IEDGE1).ENUM==0) ALLOCATE(EDGE_L(IEDGE1).ELEMENT(5),EDGE_L(IEDGE1).SUBID(5))
            EDGE_L(IEDGE1).ENUM=EDGE_L(IEDGE1).ENUM+1
            IF(EDGE_L(IEDGE1).ENUM>SIZE(EDGE_L(IEDGE1).ELEMENT,DIM=1)) THEN
				CALL I_ENLARGE_AR(EDGE_L(IEDGE1).ELEMENT,5)
				CALL I_ENLARGE_AR(EDGE_L(IEDGE1).SUBID,5)
			ENDIF
            EDGE_L(IEDGE1).ELEMENT(EDGE_L(IEDGE1).ENUM)=I			
            EDGE_L(IEDGE1).SUBID(EDGE_L(IEDGE1).ENUM)=J
		ENDDO
    END DO
    
    RETURN
ENDSUBROUTINE

!SET EDGE FOR MODEL ELEMENT. 
SUBROUTINE SETUP_EDGE_ADJL_TET(EDGE_L,NEDGE_L,ELEMENT_L,NEL_L,NODE_L,NNODE_L,ESET,NESET)
    
    IMPLICIT NONE
	INTEGER,INTENT(IN)::NEL_L,NNODE_L,NESET
    INTEGER::NEDGE_L
	TYPE(EDGE_TYDEF),ALLOCATABLE::EDGE_L(:)
	type(TET_TYDEF)::ELEMENT_L(NEL_L)
	REAL(8),INTENT(IN)::NODE_L(3,NNODE_L)
    TYPE(ESET_TYDEF),INTENT(IN)::ESET(NESET)
    INTEGER::MAXADJL1
    INTEGER::ADJL1(MAX_NODE_ADJ,NNODE_L),E_ADJL1(MAX_NODE_ADJ,NNODE_L),NADJL1(NNODE_L) 
    INTEGER::VMAX1,VMIN1,IEDGE1,ET1,NEDGE1,N1,N2,V1(2)
    INTEGER::I,J
    
    
    MAXADJL1=MAX_NODE_ADJ
    
    IF(ALLOCATED(EDGE_L)) DEALLOCATE(EDGE_L)
	ALLOCATE(EDGE_L(NNODE_L))
    ADJL1=0;E_ADJL1=0;NADJL1=1
	DO I=1,NEL_L
		ET1=ELEMENT_L(I).GMET
		NEDGE1=ELTTYPE(ET1).NEDGE
        ELEMENT_L(I).NE=NEDGE1
        
		DO J=1,NEDGE1
            V1=ELEMENT_L(I).V(ELTTYPE(ET1).EDGE(:,J))
            IF(V1(1)>V1(2)) THEN
                VMAX1=V1(1);VMIN1=V1(2)
            ELSE
                VMAX1=V1(2);VMIN1=V1(1)
            ENDIF
            !VMAX1=MAXVAL(ELEMENT_L(I).V(ELTTYPE(ET1).EDGE(:,J)))    
            !VMIN1=MINVAL(ELEMENT_L(I).V(ELTTYPE(ET1).EDGE(:,J))) 
            N2=NADJL1(VMAX1)
            !N2=MINLOC(ADJL1(:,VMAX1),MASK=ADJL1(:,VMAX1)==0,DIM=1)-1
            !IF(N2<0) THEN
            !     STOP 'ERROR IN SUB SETUP_EDGE_TBLI2. THE SIZE OF ADJL1 MAY BE TOO SMALL.'
            !ENDIF
            
            N1=MINLOC(ABS(ADJL1(1:N2,VMAX1)-VMIN1),MASK=ADJL1(1:N2,VMAX1)-VMIN1==0,DIM=1)
            if(N1>0) THEN
                IEDGE1=E_ADJL1(N1,VMAX1)            
            ELSE
                NADJL1(VMAX1)=N2+1
                IF(N2+1>MAXADJL1) THEN
                     print *, 'MAX_NODE_ADJ NEEDS TO BE  ENLARGED BY APPENDING "MAX_NODE_ADJ=XXX" TO THE KEYWORD "SolverControl".MAX_NODE_ADJ=',MAX_NODE_ADJ
                     stop
                ENDIF             
                ADJL1(N2+1,VMAX1)=VMIN1
                NEDGE_L=NEDGE_L+1
                IEDGE1=NEDGE_L
                E_ADJL1(N2+1,VMAX1)=NEDGE_L                
            ENDIF
     

            ELEMENT_L(I).E(J)=IEDGE1
            IF(IEDGE1>SIZE(EDGE_L,DIM=1)) THEN
                CALL EDGE_TYDEF_ENLARGE_AR(EDGE_L,1000)
                !MAXNEDGE=MAXNEDGE+1000
            ENDIF
            IF(.NOT.EDGE_L(IEDGE1).ISINI) THEN
                EDGE_L(IEDGE1).V=ELEMENT_L(I).V(ELTTYPE(ET1).EDGE(:,J))
                EDGE_L(IEDGE1).DIS=NORM2(NODE_L(:,EDGE_L(IEDGE1).V(1))-NODE_L(:,EDGE_L(IEDGE1).V(2)))
                EDGE_L(IEDGE1).ISINI=.TRUE.
                NEDGE_L=IEDGE1
            ENDIF
            IF(EDGE_L(IEDGE1).ENUM==0) ALLOCATE(EDGE_L(IEDGE1).ELEMENT(5),EDGE_L(IEDGE1).SUBID(5))
            EDGE_L(IEDGE1).ENUM=EDGE_L(IEDGE1).ENUM+1
            IF(EDGE_L(IEDGE1).ENUM>SIZE(EDGE_L(IEDGE1).ELEMENT,DIM=1)) THEN
				CALL I_ENLARGE_AR(EDGE_L(IEDGE1).ELEMENT,5)
				CALL I_ENLARGE_AR(EDGE_L(IEDGE1).SUBID,5)
			ENDIF
            EDGE_L(IEDGE1).ELEMENT(EDGE_L(IEDGE1).ENUM)=I			
            EDGE_L(IEDGE1).SUBID(EDGE_L(IEDGE1).ENUM)=J
		ENDDO
    END DO
    
    
    
    RETURN
ENDSUBROUTINE


SUBROUTINE SETUP_FACE_ADJL_TET(FACE_L,NFACE_L,EDGE_L,NEDGE_L,ELEMENT_L,NEL_L,NODE_L,NNODE_L,ESET,NESET)
    
    IMPLICIT NONE
    INTEGER,INTENT(IN)::NEL_L,NEDGE_L,NNODE_L,NESET
    INTEGER::NFACE_L
	TYPE(EDGE_TYDEF),INTENT(IN)::EDGE_L(NEDGE_L)
    TYPE(FACE_TYDEF),ALLOCATABLE::FACE_L(:)
	type(TET_TYDEF)::ELEMENT_L(NEL_L)
    REAL(8),INTENT(IN)::NODE_L(3,NNODE_L)
    TYPE(ESET_TYDEF),INTENT(IN)::ESET(NESET)
    INTEGER::I,J,N1,ET1,NFACE1,N2,K,K1,NV1,IFACE1,N2ORDER1(10),NODE1(10),O2NODE1(10),ORDER1(10),NODE_SORT1(10),N2ORDER2(10)
    REAL(8)::V1(3),V2(3),NORMAL1(3),T1
    INTEGER::MAXADJL1
    INTEGER::ADJL1(3,MAX_FACE_ADJ,NNODE_L),F_ADJL1(MAX_FACE_ADJ,NNODE_L),NADJL1(NNODE_L) 

    MAXADJL1=MAX_FACE_ADJ
    
    IF(.NOT.ISINI_GMSHET) THEN
        CALL Initialize_et2numNodes()
        CALL ET_GMSH_EDGE_FACE()
        ISINI_GMSHET=.TRUE.
    ENDIF 	

	IF(ALLOCATED(FACE_L)) DEALLOCATE(FACE_L)
	ALLOCATE(FACE_L(NNODE_L))
    ADJL1=0;F_ADJL1=0;NADJL1=0
	DO I=1,NEL_L
		ET1=ELEMENT_L(I).GMET
		NFACE1=ELTTYPE(ET1).NFACE
		ELEMENT_L(I).NF=NFACE1
        NV1=ELTTYPE(ET1).NNODE
        NODE1(1:NV1)=ELEMENT_L(I).V(1:NV1)
        !NODE_SORT1(1:NV1)=NODE1(1:NV1)
        O2NODE1(1:NV1)=[1:NV1]
        DO J=1,NV1-1
            !N1=MINLOC(NODE1(1:NV1),MASK=NODE1(1:NV1)>0,DIM=1)
            DO K=J+1,NV1
                IF(NODE1(K)<NODE1(J)) THEN
                    N1=NODE1(J);NODE1(J)=NODE1(K);NODE1(K)=N1;
                    N1=O2NODE1(J);O2NODE1(J)=O2NODE1(K);O2NODE1(K)=N1
                    !N2ORDER1(J)=O2NODE1(J);
                    !N2ORDER1(N1)=J;NODE1(N1)=-1
                    !O2NODE1(J)=N1
                ENDIF
            ENDDO
        ENDDO
        N2ORDER1(O2NODE1(1:NV1))=[1:NV1]
        !N2ORDER2(1:NV1)=[1:NV1]
        !CALL QUICK_SORT(NODE_SORT1(1:NV1),N2ORDER2(1:NV1))
        
 		DO J=1,NFACE1
            N2=ELTTYPE(ET1).FACE(0,J)            
            ORDER1(1:N2)=N2ORDER1(ELTTYPE(ET1).FACE(1:N2,J))
            NODE1(1:NV1)=0
            NODE1(ORDER1(1:N2))=ELEMENT_L(I).V(O2NODE1(ORDER1(1:N2)))
            N1=0
            DO K=1,NV1
                IF(NODE1(K)>0) THEN
                    N1=N1+1
                    IF(N1<K) NODE1(N1)=NODE1(K)
                ENDIF
            ENDDO
            !N1=MINLOC(F_ADJL1(:,NODE1(N2)),MASK=F_ADJL1(:,NODE1(N2))==0,DIM=1)-1
            !IF(N1<0) THEN                
            !    STOP "PLEASE ENLARGE ARRAY SIZE. SUB=SETUP_FACE_ADJL_TET"
            !ENDIF
            N1=NADJL1(NODE1(N2))
            !N2<=4
            IFACE1=0
            !if(i==23.and.j==3) then
            !    print *,i
            !endif
            DO K=1,N1
                DO K1=1,N2-1
                    if(ADJL1(K1,K,NODE1(N2))-NODE1(K1)/=0) EXIT
                ENDDO
                IF(K1==N2) THEN
                    IFACE1=F_ADJL1(K,NODE1(N2))
                    EXIT
                ENDIF 
            ENDDO
            IF(IFACE1==0) THEN
                NADJL1(NODE1(N2))=N1+1
                IF(N1+1>MAXADJL1) THEN                
                    print *, 'MAX_FACE_ADJ NEEDS TO BE  ENLARGED BY APPENDING "MAX_FACE_ADJ=XXX" TO THE KEYWORD "SolverControl".MAX_FACEE_ADJ=',MAX_FACE_ADJ
                    stop
                ENDIF
                ADJL1(1:N2-1,N1+1,NODE1(N2))=NODE1(1:N2-1)
                NFACE_L=NFACE_L+1
                F_ADJL1(N1+1,NODE1(N2))=NFACE_L
                IFACE1=NFACE_L
            ENDIF

            ELEMENT_L(I).F(J)=IFACE1
            
            IF(IFACE1>SIZE(FACE_L,DIM=1)) THEN
                CALL FACE_TYDEF_ENLARGE_AR(FACE_L,1000)
                !MAXNFACE=MAXNFACE+1000
            ENDIF
            IF(.NOT.FACE_L(IFACE1).ISINI) THEN
                FACE_L(IFACE1).SHAPE=ELTTYPE(ET1).FACE(0,J)
                FACE_L(IFACE1).V(1:FACE_L(IFACE1).SHAPE)=ELEMENT_L(I).V(ELTTYPE(ET1).FACE(1:FACE_L(IFACE1).SHAPE,J))
				FACE_L(IFACE1).EDGE(1:FACE_L(IFACE1).SHAPE)=ELEMENT_L(I).E(ABS(ELTTYPE(ET1).FACEEDGE(1:FACE_L(IFACE1).SHAPE,J)))
				!CHECK LOOP ORDER
				DO K=1,FACE_L(IFACE1).SHAPE
					IF(EDGE_L(FACE_L(IFACE1).EDGE(K)).V(1)/=FACE_L(IFACE1).V(K)) FACE_L(IFACE1).EDGE(K)=-FACE_L(IFACE1).EDGE(K)
				ENDDO
                !V1=NODE_L(:,FACE_L(IFACE1).V(2))-NODE_L(:,FACE_L(IFACE1).V(1))
                !V2=NODE_L(:,FACE_L(IFACE1).V(3))-NODE_L(:,FACE_L(IFACE1).V(1))
                !call r8vec_cross_3d ( v1, v2, NORMAL1(1:3) )
                !T1 = sqrt ( sum ( ( NORMAL1 )**2 ) )
                !if ( T1 /= 0.0D+00 ) then
                !  NORMAL1 = NORMAL1/T1
                !end if
                !FACE_L(IFACE1).UNORMAL=NORMAL1
                FACE_L(IFACE1).ISINI=.TRUE.
                NFACE_L=IFACE1
				DO K=1,3
					FACE_L(NFACE_L).BBOX(1,K)=MINVAL(NODE_L(K,FACE_L(NFACE_L).V(1:FACE_L(NFACE_L).SHAPE)))-VTOL
					FACE_L(NFACE_L).BBOX(2,K)=MAXVAL(NODE_L(K,FACE_L(NFACE_L).V(1:FACE_L(NFACE_L).SHAPE)))+VTOL
				ENDDO
            ENDIF
            IF(FACE_L(IFACE1).ENUM==0) THEN
                ALLOCATE(FACE_L(IFACE1).ELEMENT(2))
                ALLOCATE(FACE_L(IFACE1).SUBID(2))
            ENDIF
            FACE_L(IFACE1).ENUM=FACE_L(IFACE1).ENUM+1
            IF(FACE_L(IFACE1).ENUM>SIZE(FACE_L(IFACE1).ELEMENT,DIM=1)) THEN
                CALL I_ENLARGE_AR(FACE_L(IFACE1).ELEMENT,5)
                CALL I_ENLARGE_AR(FACE_L(IFACE1).SUBID,5)
            ENDIF
            
            IF(FACE_L(IFACE1).ENUM>1) THEN
                DO K=1,FACE_L(IFACE1).SHAPE
                    IF(FACE_L(IFACE1).V(K)==ELEMENT_L(I).V(ELTTYPE(ELEMENT_L(I).GMET).FACE(1,J))) THEN
                        IF(FACE_L(IFACE1).V(MOD(K,FACE_L(IFACE1).SHAPE)+1)/= &
                            ELEMENT_L(I).V(ELTTYPE(ELEMENT_L(I).GMET).FACE(2,J))) THEN
                            FACE_L(IFACE1).ELEMENT(FACE_L(IFACE1).ENUM)=-I !I��Ԫ����ķ�����face�ķ�����,Ϊ-1
                        ELSE
                            FACE_L(IFACE1).ELEMENT(FACE_L(IFACE1).ENUM)=I
                        ENDIF                  
                        EXIT
                    ENDIF
                ENDDO
            ELSE
                FACE_L(IFACE1).ELEMENT(FACE_L(IFACE1).ENUM)=I 
            ENDIF
            FACE_L(IFACE1).SUBID(FACE_L(IFACE1).ENUM)=J
                
            
		ENDDO
    END DO    
    
    
ENDSUBROUTINE

SUBROUTINE SETUP_FACE_ADJL_MODEL(FACE_L,NFACE_L,EDGE_L,NEDGE_L,ELEMENT_L,NEL_L,NODE_L,NNODE_L,ESET,NESET)
    
    IMPLICIT NONE
    INTEGER,INTENT(IN)::NEL_L,NEDGE_L,NNODE_L,NESET
    INTEGER::NFACE_L
	TYPE(EDGE_TYDEF),INTENT(IN)::EDGE_L(NEDGE_L)
    TYPE(FACE_TYDEF),ALLOCATABLE::FACE_L(:)
	type(ELEMENT_TYDEF)::ELEMENT_L(NEL_L)    
    REAL(8),INTENT(IN)::NODE_L(3,NNODE_L)
    TYPE(ESET_TYDEF),INTENT(IN)::ESET(NESET)
    INTEGER::I,J,N1,ET1,NFACE1,N2,K,K1,NV1,IFACE1,N2ORDER1(10),NODE1(10),O2NODE1(10),ORDER1(10)
    REAL(8)::V1(3),V2(3),NORMAL1(3),T1
    INTEGER::MAXADJL1
    INTEGER::ADJL1(3,MAX_FACE_ADJ,NNODE_L),F_ADJL1(MAX_FACE_ADJ,NNODE_L),NADJL1(NNODE_L) 

    MAXADJL1=MAX_FACE_ADJ

    
    IF(.NOT.ISINI_GMSHET) THEN
        CALL Initialize_et2numNodes()
        CALL ET_GMSH_EDGE_FACE()
        ISINI_GMSHET=.TRUE.
    ENDIF 	

	IF(ALLOCATED(FACE_L)) DEALLOCATE(FACE_L)
	ALLOCATE(FACE_L(NNODE_L))
    ADJL1=0;F_ADJL1=0;NADJL1=0
	DO I=1,NEL_L
		ET1=GETGMSHET(ESET(ELEMENT_L(I).ISET).ET)
		NFACE1=ELTTYPE(ET1).NFACE
		ELEMENT_L(I).NFACE=NFACE1
        NV1=ELTTYPE(ET1).NNODE
        NODE1(1:NV1)=ELEMENT_L(I).NODE(1:NV1)
        ALLOCATE(ELEMENT_L(I).FACE(NFACE1))
        !DO J=1,NV1
        !    N1=MINLOC(NODE1(1:NV1),MASK=NODE1(1:NV1)>0,DIM=1)
        !    N2ORDER1(N1)=J;NODE1(N1)=-1
        !    O2NODE1(J)=N1
        !ENDDO
        O2NODE1(1:NV1)=[1:NV1]
        DO J=1,NV1-1
            !N1=MINLOC(NODE1(1:NV1),MASK=NODE1(1:NV1)>0,DIM=1)
            DO K=J+1,NV1
                IF(NODE1(K)<NODE1(J)) THEN
                    N1=NODE1(J);NODE1(J)=NODE1(K);NODE1(K)=N1;
                    N1=O2NODE1(J);O2NODE1(J)=O2NODE1(K);O2NODE1(K)=N1
                    !N2ORDER1(J)=O2NODE1(J);
                    !N2ORDER1(N1)=J;NODE1(N1)=-1
                    !O2NODE1(J)=N1
                ENDIF
            ENDDO
        ENDDO
        N2ORDER1(O2NODE1(1:NV1))=[1:NV1]
        
 		DO J=1,NFACE1
            N2=ELTTYPE(ET1).FACE(0,J)            
            ORDER1(1:N2)=N2ORDER1(ELTTYPE(ET1).FACE(1:N2,J))
            NODE1(1:NV1)=0
            NODE1(ORDER1(1:N2))=ELEMENT_L(I).NODE(O2NODE1(ORDER1(1:N2)))
            N1=0
            DO K=1,NV1
                IF(NODE1(K)>0) THEN
                    N1=N1+1
                    IF(N1<K) NODE1(N1)=NODE1(K)
                ENDIF
            ENDDO

            
            !N1=MINLOC(F_ADJL1(:,NODE1(N2)),MASK=F_ADJL1(:,NODE1(N2))==0,DIM=1)-1
            !IF(N1<0) THEN                
            !    STOP "PLEASE ENLARGE ARRAY SIZE. SUB=SETUP_FACE_ADJL_MODEL"
            !ENDIF
            !!N2<=4
            !IFACE1=0
            !DO K=1,N1
            !    if(SUM(ABS(ADJL1(1:N2-1,K,NODE1(N2))-NODE1(1:N2-1)))==0) THEN
            !        IFACE1=F_ADJL1(K,NODE1(N2))
            !        EXIT
            !    ENDIF                
            !ENDDO
            !IF(IFACE1==0) THEN

            N1=NADJL1(NODE1(N2))
            !N2<=4
            IFACE1=0
            !if(i==23.and.j==3) then
            !    print *,i
            !endif
            DO K=1,N1
                DO K1=1,N2-1
                    if(ADJL1(K1,K,NODE1(N2))-NODE1(K1)/=0) EXIT
                ENDDO
                IF(K1==N2) THEN
                    IFACE1=F_ADJL1(K,NODE1(N2))
                    EXIT
                ENDIF 
            ENDDO
            IF(IFACE1==0) THEN
                NADJL1(NODE1(N2))=N1+1
                IF(N1+1>MAXADJL1) THEN                
                    print *, 'MAX_FACE_ADJ NEEDS TO BE  ENLARGED BY APPENDING "MAX_FACE_ADJ=XXX" TO THE KEYWORD "SolverControl".MAX_FACEE_ADJ=',MAX_FACE_ADJ
                    stop
                ENDIF
                ADJL1(1:N2-1,N1+1,NODE1(N2))=NODE1(1:N2-1)
                NFACE_L=NFACE_L+1
                F_ADJL1(N1+1,NODE1(N2))=NFACE_L
                IFACE1=NFACE_L
                !IF(NODE1(N2)==50) THEN
                !    PRINT *, I,J,NODE1(1:N2)
                !ENDIF            
            ENDIF

            ELEMENT_L(I).FACE(J)=IFACE1
            
            IF(IFACE1>SIZE(FACE_L,DIM=1)) THEN
                CALL FACE_TYDEF_ENLARGE_AR(FACE_L,1000)
                !MAXNFACE=MAXNFACE+1000
            ENDIF
            IF(.NOT.FACE_L(IFACE1).ISINI) THEN
                FACE_L(IFACE1).SHAPE=ELTTYPE(ET1).FACE(0,J)
                FACE_L(IFACE1).V(1:FACE_L(IFACE1).SHAPE)=ELEMENT_L(I).NODE(ELTTYPE(ET1).FACE(1:FACE_L(IFACE1).SHAPE,J))
				FACE_L(IFACE1).EDGE(1:FACE_L(IFACE1).SHAPE)=ELEMENT_L(I).EDGE(ABS(ELTTYPE(ET1).FACEEDGE(1:FACE_L(IFACE1).SHAPE,J)))
				!CHECK LOOP ORDER
				DO K=1,FACE_L(IFACE1).SHAPE
					IF(EDGE_L(FACE_L(IFACE1).EDGE(K)).V(1)/=FACE_L(IFACE1).V(K)) FACE_L(IFACE1).EDGE(K)=-FACE_L(IFACE1).EDGE(K)
				ENDDO
                !V1=NODE_L(:,FACE_L(IFACE1).V(2))-NODE_L(:,FACE_L(IFACE1).V(1))
                !V2=NODE_L(:,FACE_L(IFACE1).V(3))-NODE_L(:,FACE_L(IFACE1).V(1))
                !call r8vec_cross_3d ( v1, v2, NORMAL1(1:3) )
                !T1 = sqrt ( sum ( ( NORMAL1 )**2 ) )
                !if ( T1 /= 0.0D+00 ) then
                !  NORMAL1 = NORMAL1/T1
                !end if
                !FACE_L(IFACE1).UNORMAL=NORMAL1
                FACE_L(IFACE1).ISINI=.TRUE.
                NFACE_L=IFACE1
				DO K=1,3
					FACE_L(NFACE_L).BBOX(1,K)=MINVAL(NODE_L(K,FACE_L(NFACE_L).V(1:FACE_L(NFACE_L).SHAPE)))-VTOL
					FACE_L(NFACE_L).BBOX(2,K)=MAXVAL(NODE_L(K,FACE_L(NFACE_L).V(1:FACE_L(NFACE_L).SHAPE)))+VTOL
				ENDDO
            ENDIF
            IF(FACE_L(IFACE1).ENUM==0) THEN
                ALLOCATE(FACE_L(IFACE1).ELEMENT(2))
                ALLOCATE(FACE_L(IFACE1).SUBID(2))
            ENDIF
            FACE_L(IFACE1).ENUM=FACE_L(IFACE1).ENUM+1
            IF(FACE_L(IFACE1).ENUM>SIZE(FACE_L(IFACE1).ELEMENT,DIM=1)) THEN
                CALL I_ENLARGE_AR(FACE_L(IFACE1).ELEMENT,5)
                CALL I_ENLARGE_AR(FACE_L(IFACE1).SUBID,5)
            ENDIF
            
            IF(FACE_L(IFACE1).ENUM>1) THEN
                DO K=1,FACE_L(IFACE1).SHAPE
                    IF(FACE_L(IFACE1).V(K)==ELEMENT_L(I).NODE(ELTTYPE(ET1).FACE(1,J))) THEN
                        IF(FACE_L(IFACE1).V(MOD(K,FACE_L(IFACE1).SHAPE)+1)/= &
                            ELEMENT_L(I).NODE(ELTTYPE(ET1).FACE(2,J))) THEN
                            FACE_L(IFACE1).ELEMENT(FACE_L(IFACE1).ENUM)=-I !I��Ԫ����ķ�����face�ķ�����,Ϊ-1
                        ELSE
                            FACE_L(IFACE1).ELEMENT(FACE_L(IFACE1).ENUM)=I
                        ENDIF                  
                        EXIT
                    ENDIF
                ENDDO
            ELSE
                FACE_L(IFACE1).ELEMENT(FACE_L(IFACE1).ENUM)=I 
            ENDIF
            FACE_L(IFACE1).SUBID(FACE_L(IFACE1).ENUM)=J
                
            
		ENDDO
    END DO    
    
    
ENDSUBROUTINE

SUBROUTINE SETUP_FACE_ADJL_Solver(FACE_L,NFACE_L,EDGE_L,NEDGE_L)
    USE solverds,ONLY:ELEMENT_L=>ELEMENT,NODE_L=>NODE,NNODE_L=>NNUM,NEL_L=>ENUM,ESET=>eset,NESET=>NESET,ESETID
    IMPLICIT NONE
    INTEGER,INTENT(IN)::NEDGE_L
    INTEGER::NFACE_L
	TYPE(EDGE_TYDEF),INTENT(IN)::EDGE_L(NEDGE_L)
    TYPE(FACE_TYDEF),ALLOCATABLE::FACE_L(:)
	!type(ELEMENT_TYDEF)::ELEMENT_L(NEL_L)    
    !REAL(8),INTENT(IN)::NODE_L(3,NNODE_L)
    !TYPE(ESET_TYDEF),INTENT(IN)::ESET(NESET)
    INTEGER::I,J,N1,ET1,NFACE1,N2,K,K1,NV1,IFACE1,N2ORDER1(10),NODE1(10),O2NODE1(10),ORDER1(10)
    REAL(8)::V1(3),V2(3),NORMAL1(3),T1
    INTEGER::MAXADJL1
    INTEGER::ADJL1(3,MAX_FACE_ADJ,NNODE_L),F_ADJL1(MAX_FACE_ADJ,NNODE_L),NADJL1(NNODE_L) 

    MAXADJL1=MAX_FACE_ADJ

    
    IF(.NOT.ISINI_GMSHET) THEN
        CALL Initialize_et2numNodes() !根据gmsh单元类型和格式，确定每个单元的节点数
        CALL ET_GMSH_EDGE_FACE() !根据gmsh单元类型和格式，确定每个边、面数
        ISINI_GMSHET=.TRUE.
    ENDIF 	

	IF(ALLOCATED(FACE_L)) DEALLOCATE(FACE_L)
	ALLOCATE(FACE_L(NNODE_L))
    ADJL1=0;F_ADJL1=0;NADJL1=0
	DO I=1,NEL_L
		ET1=GETGMSHET(ESET(ELEMENT_L(I).SET).ET)
		NFACE1=ELTTYPE(ET1).NFACE
		ELEMENT_L(I).NFACE=NFACE1
        NV1=ELTTYPE(ET1).NNODE
        NODE1(1:NV1)=ELEMENT_L(I).NODE(1:NV1)
        ALLOCATE(ELEMENT_L(I).FACE(NFACE1),ELEMENT_L(I).ADJELT(NFACE1))
        IF(ALLOCATED(ELEMENT_L(I).ADJELT)) ELEMENT_L(I).ADJELT=0
        !DO J=1,NV1
        !    N1=MINLOC(NODE1(1:NV1),MASK=NODE1(1:NV1)>0,DIM=1)
        !    N2ORDER1(N1)=J;NODE1(N1)=-1
        !    O2NODE1(J)=N1
        !ENDDO
        O2NODE1(1:NV1)=[1:NV1] !order2node
        DO J=1,NV1-1
            !N1=MINLOC(NODE1(1:NV1),MASK=NODE1(1:NV1)>0,DIM=1)
            DO K=J+1,NV1
                IF(NODE1(K)<NODE1(J)) THEN
                    N1=NODE1(J);NODE1(J)=NODE1(K);NODE1(K)=N1;
                    N1=O2NODE1(J);O2NODE1(J)=O2NODE1(K);O2NODE1(K)=N1
                    !N2ORDER1(J)=O2NODE1(J);
                    !N2ORDER1(N1)=J;NODE1(N1)=-1
                    !O2NODE1(J)=N1
                ENDIF
            ENDDO
        ENDDO
        N2ORDER1(O2NODE1(1:NV1))=[1:NV1] !node2order
        
 		DO J=1,NFACE1
            N2=ELTTYPE(ET1).FACE(0,J)            
            ORDER1(1:N2)=N2ORDER1(ELTTYPE(ET1).FACE(1:N2,J))
            NODE1(1:NV1)=0
            NODE1(ORDER1(1:N2))=ELEMENT_L(I).NODE(O2NODE1(ORDER1(1:N2)))
            N1=0
            DO K=1,NV1
                IF(NODE1(K)>0) THEN
                    N1=N1+1
                    IF(N1<K) NODE1(N1)=NODE1(K)
                ENDIF
            ENDDO

            
            !N1=MINLOC(F_ADJL1(:,NODE1(N2)),MASK=F_ADJL1(:,NODE1(N2))==0,DIM=1)-1
            !IF(N1<0) THEN                
            !    STOP "PLEASE ENLARGE ARRAY SIZE. SUB=SETUP_FACE_ADJL_MODEL"
            !ENDIF
            !!N2<=4
            !IFACE1=0
            !DO K=1,N1
            !    if(SUM(ABS(ADJL1(1:N2-1,K,NODE1(N2))-NODE1(1:N2-1)))==0) THEN
            !        IFACE1=F_ADJL1(K,NODE1(N2))
            !        EXIT
            !    ENDIF                
            !ENDDO
            !IF(IFACE1==0) THEN

            N1=NADJL1(NODE1(N2))
            !N2<=4
            IFACE1=0
            !if(i==23.and.j==3) then
            !    print *,i
            !endif
            DO K=1,N1
                DO K1=1,N2-1
                    if(ADJL1(K1,K,NODE1(N2))-NODE1(K1)/=0) EXIT
                ENDDO
                IF(K1==N2) THEN
                    IFACE1=F_ADJL1(K,NODE1(N2))
                    EXIT
                ENDIF 
            ENDDO
            IF(IFACE1==0) THEN
                NADJL1(NODE1(N2))=N1+1
                IF(N1+1>MAXADJL1) THEN                
                    print *, 'MAX_FACE_ADJ NEEDS TO BE  ENLARGED BY APPENDING "MAX_FACE_ADJ=XXX" TO THE KEYWORD "SolverControl".MAX_FACEE_ADJ=',MAX_FACE_ADJ
                    stop
                ENDIF
                ADJL1(1:N2-1,N1+1,NODE1(N2))=NODE1(1:N2-1)
                NFACE_L=NFACE_L+1
                F_ADJL1(N1+1,NODE1(N2))=NFACE_L
                IFACE1=NFACE_L
                !IF(NODE1(N2)==50) THEN
                !    PRINT *, I,J,NODE1(1:N2)
                !ENDIF            
            ENDIF

            ELEMENT_L(I).FACE(J)=IFACE1
            
            IF(IFACE1>SIZE(FACE_L,DIM=1)) THEN
                CALL FACE_TYDEF_ENLARGE_AR(FACE_L,1000)
                !MAXNFACE=MAXNFACE+1000
            ENDIF
            IF(.NOT.FACE_L(IFACE1).ISINI) THEN
                FACE_L(IFACE1).SHAPE=ELTTYPE(ET1).FACE(0,J)
                FACE_L(IFACE1).V(1:FACE_L(IFACE1).SHAPE)=ELEMENT_L(I).NODE(ELTTYPE(ET1).FACE(1:FACE_L(IFACE1).SHAPE,J))
				FACE_L(IFACE1).EDGE(1:FACE_L(IFACE1).SHAPE)=ELEMENT_L(I).EDGE(ABS(ELTTYPE(ET1).FACEEDGE(1:FACE_L(IFACE1).SHAPE,J)))
				!CHECK LOOP ORDER
				DO K=1,FACE_L(IFACE1).SHAPE
					IF(EDGE_L(FACE_L(IFACE1).EDGE(K)).V(1)/=FACE_L(IFACE1).V(K)) FACE_L(IFACE1).EDGE(K)=-FACE_L(IFACE1).EDGE(K)
				ENDDO
                !V1=NODE_L(:,FACE_L(IFACE1).V(2))-NODE_L(:,FACE_L(IFACE1).V(1))
                !V2=NODE_L(:,FACE_L(IFACE1).V(3))-NODE_L(:,FACE_L(IFACE1).V(1))
                !call r8vec_cross_3d ( v1, v2, NORMAL1(1:3) )
                !T1 = sqrt ( sum ( ( NORMAL1 )**2 ) )
                !if ( T1 /= 0.0D+00 ) then
                !  NORMAL1 = NORMAL1/T1
                !end if
                !FACE_L(IFACE1).UNORMAL=NORMAL1
                FACE_L(IFACE1).ISINI=.TRUE.
                NFACE_L=IFACE1
				DO K=1,3
					FACE_L(NFACE_L).BBOX(1,K)=MINVAL(NODE_L(FACE_L(NFACE_L).V(1:FACE_L(NFACE_L).SHAPE)).COORD(K))-VTOL
					FACE_L(NFACE_L).BBOX(2,K)=MAXVAL(NODE_L(FACE_L(NFACE_L).V(1:FACE_L(NFACE_L).SHAPE)).COORD(K))+VTOL
				ENDDO
            ENDIF
            IF(FACE_L(IFACE1).ENUM==0) THEN
                ALLOCATE(FACE_L(IFACE1).ELEMENT(2))
                ALLOCATE(FACE_L(IFACE1).SUBID(2))
            ENDIF
            FACE_L(IFACE1).ENUM=FACE_L(IFACE1).ENUM+1
            IF(FACE_L(IFACE1).ENUM>SIZE(FACE_L(IFACE1).ELEMENT,DIM=1)) THEN
                CALL I_ENLARGE_AR(FACE_L(IFACE1).ELEMENT,5)
                CALL I_ENLARGE_AR(FACE_L(IFACE1).SUBID,5)
            ENDIF
            
            IF(FACE_L(IFACE1).ENUM>1) THEN
                DO K=1,FACE_L(IFACE1).SHAPE
                    IF(FACE_L(IFACE1).V(K)==ELEMENT_L(I).NODE(ELTTYPE(ET1).FACE(1,J))) THEN
                        IF(FACE_L(IFACE1).V(MOD(K,FACE_L(IFACE1).SHAPE)+1)/= &
                            ELEMENT_L(I).NODE(ELTTYPE(ET1).FACE(2,J))) THEN
                            FACE_L(IFACE1).ELEMENT(FACE_L(IFACE1).ENUM)=-I !I��Ԫ����ķ�����face�ķ�����,Ϊ-1
                        ELSE
                            FACE_L(IFACE1).ELEMENT(FACE_L(IFACE1).ENUM)=I
                        ENDIF                  
                        EXIT
                    ENDIF
                ENDDO
            ELSE
                FACE_L(IFACE1).ELEMENT(FACE_L(IFACE1).ENUM)=I 
            ENDIF
            FACE_L(IFACE1).SUBID(FACE_L(IFACE1).ENUM)=J
                
            
		ENDDO
    END DO    
    
    
ENDSUBROUTINE


END MODULE
    