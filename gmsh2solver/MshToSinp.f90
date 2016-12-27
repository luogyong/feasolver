program GmshToSolver
	use ifcore
	use ifqwin
	use DS_Gmsh2Solver
	implicit none
	integer::i
	character(1)::key
	type(qwinfo) winfo
	LOGICAL(4)::tof,pressed
	integer,allocatable::IPERM(:)
	
	winfo%TYPE = QWIN$MAX
	tof=SETWSIZEQQ(QWIN$FRAMEWINDOW, winfo)
	tof=SETWSIZEQQ(0, winfo)  
	
	print *, 'LGY Works. Msh2Sinp'
	
    IF(NOPOPUP==0) THEN
	    write(*, 10) 
	    key=getcharqq()
	    if(ichar(key)==ichar('h').or.ichar(key)==ichar('H')) then		
		    call write_readme_gmsh2sinp()			
		    stop
        end if
    ENDIF

	
	print *, 'Begin to read in data...'
	call readin()
	
	!统计每个物理组单元的个数及单元数组下标。
	do i=1,maxphgp
		if(physicalgroup(i).nel>0) allocate(physicalgroup(i).element(physicalgroup(i).nel))
	end do
	
	physicalgroup.nel=0
	do i=1,nel
		physicalgroup(element(i).tag(1)).nel=physicalgroup(element(i).tag(1)).nel+1
		physicalgroup(element(i).tag(1)).element(physicalgroup(element(i).tag(1)).nel)=i		
	end do	
	
	CALL GEN_ZEROTHICKNESS_ELEMENT()
	IF(ISGENLAYER) CALL GEN_LAYERED_ELEMENT()
	IF(ISGENTRISURFACE>0) THEN
        ALLOCATE(EDGE(MAXNEDGE),FACE(MAXNFACE))
		CALL SETUP_EDGE_TBL()
		CALL SETUP_FACE_TBL()
        CALL GEN_TRISURFACE()
        STOP "DONE IN GENERATING TRISURFACE."
	ENDIF
	! reordering nodal number
	allocate(Noutputorder(nnode))
	ALLOCATE(IPERM(nnode))
	call setup_adjList()
	call reorder_nodal_number(IPERM,nnode,adjL,maxadj)
	do i=1,nnode
		noutputorder(IPERM(i))=i	
		node(i).inode=IPERM(i)		
	end do
	DEALLOCATE(IPERM)
	

	
	CALL elt_bc_load_translate()
	
	call Tosolver()
	
	Print *, 'GMSHTOSOLVER IS DONE. PLEASE ADD OTHER PARAMETERS BY HAND.'
    IF(NOPOPUP/=0) I =setexitqq(QWIN$EXITNOPERSIST)
10 format("Press 'H' to write a keyword help file named 'D:\README_GMSH2SINP.TXT'. Any other key to read an Msh file.")		
end program

SUBROUTINE SETUP_EDGE_TBL()
    USE DS_Gmsh2Solver
    USE hashtbl
    IMPLICIT NONE
    INTEGER::I,J,N1(2),ET1,NEDGE1,TBL_LEN
	CHARACTER(LEN=:),ALLOCATABLE::KEY1
    TYPE(DICT_DATA)::VAL1
    CHARACTER(64)::CKEY1
    INTEGER::HKEY1
    
	TBL_LEN=NNODE	
    CALL EDGE_TBL.INIT(TBL_LEN)
	DO I=1,NEL
		ET1=ELEMENT(I).ET
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
            IF(EDGE(VAL1.IITEM).NEL==0) ALLOCATE(EDGE(VAL1.IITEM).ELEMENT(5),EDGE(VAL1.IITEM).SUBID(5))
            EDGE(VAL1.IITEM).NEL=EDGE(VAL1.IITEM).NEL+1
            IF(EDGE(VAL1.IITEM).NEL>SIZE(EDGE(VAL1.IITEM).ELEMENT,DIM=1)) THEN
				CALL I_ENLARGE_AR(EDGE(VAL1.IITEM).ELEMENT,5)
				CALL I_ENLARGE_AR(EDGE(VAL1.IITEM).SUBID,5)
			ENDIF
            EDGE(VAL1.IITEM).ELEMENT(EDGE(VAL1.IITEM).NEL)=I			
            EDGE(VAL1.IITEM).SUBID(EDGE(VAL1.IITEM).NEL)=J
		ENDDO
    END DO
    RETURN
ENDSUBROUTINE

SUBROUTINE SETUP_EDGE_TBLi(EDGE_TBL_L,TBL_SIZE_L,EDGE_L,NEDGE_L,ELEMENT_L,NEL_L)
    USE DS_Gmsh2Solver
    USE hashtbl
    IMPLICIT NONE
	INTEGER::TBL_SIZE_L,NEDGE_L,NEL_L
	TYPE(hash_tbl_sll)::EDGE_TBL_L
	TYPE(EDGE_TYDEF),ALLOCATABLE::EDGE_L(:)
	type(element_type)::ELEMENT_L(NEL_L)
	
    INTEGER::I,J,N1(2),ET1,NEDGE1
	CHARACTER(LEN=:),ALLOCATABLE::KEY1
    TYPE(DICT_DATA)::VAL1
    CHARACTER(64)::CKEY1
    INTEGER::HKEY1
    
	IF(.NOT.ALLOCATED(EDGE_L)) ALLOCATE(EDGE_L(TBL_SIZE_L))
	
    IF(.NOT.EDGE_TBL_L.IS_INIT) CALL EDGE_TBL_L.INIT(TBL_SIZE_L)
	DO I=1,NEL_L
		ET1=ELEMENT_L(I).ET
		NEDGE1=ELTTYPE(ET1).NEDGE
        IF(.NOT.ALLOCATED(ELEMENT_L(I).EDGE)) ALLOCATE(ELEMENT_L(I).EDGE(NEDGE1))
		DO J=1,NEDGE1
            N1=ELEMENT_L(I).NODE(ELTTYPE(ET1).EDGE(:,J))
            CALL I2C_KEY_hash_tbl_sll(KEY1,N1(:),2)
            CKEY1=TRIM(KEY1)
            HKEY1 = MOD(ABS(HASH_DJB(KEY1)),EDGE_TBL_L%vec_len)
            VAL1.IEL=I;VAL1.ISE=J;VAL1.IITEM=EDGE_TBL_L.TBL_ID
			CALL EDGE_TBL_L.PUT(TRIM(KEY1),VAL1,HKEY1)
			ELEMENT_L(I).EDGE(J)=VAL1.IITEM
            IF(VAL1.IITEM>SIZE(EDGE_L,DIM=1)) THEN
                CALL EDGE_TYDEF_ENLARGE_AR(EDGE_L,1000)
                !MAXNEDGE=MAXNEDGE+1000
            ENDIF
            IF(.NOT.EDGE_L(VAL1.IITEM).ISINI) THEN
                EDGE_L(VAL1.IITEM).CKEY=CKEY1
                EDGE_L(VAL1.IITEM).HKEY=HKEY1
                EDGE_L(VAL1.IITEM).V=ELEMENT_L(I).NODE(ELTTYPE(ET1).EDGE(:,J))
                EDGE_L(VAL1.IITEM).ISINI=.TRUE.
                NEDGE_L=VAL1.IITEM
            ENDIF
            IF(EDGE_L(VAL1.IITEM).NEL==0) ALLOCATE(EDGE_L(VAL1.IITEM).ELEMENT(5),EDGE_L(VAL1.IITEM).SUBID(5))
            EDGE_L(VAL1.IITEM).NEL=EDGE_L(VAL1.IITEM).NEL+1
            IF(EDGE_L(VAL1.IITEM).NEL>SIZE(EDGE_L(VAL1.IITEM).ELEMENT,DIM=1)) THEN
				CALL I_ENLARGE_AR(EDGE_L(VAL1.IITEM).ELEMENT,5)
				CALL I_ENLARGE_AR(EDGE_L(VAL1.IITEM).SUBID,5)
			ENDIF
            EDGE_L(VAL1.IITEM).ELEMENT(EDGE_L(VAL1.IITEM).NEL)=I			
            EDGE_L(VAL1.IITEM).SUBID(EDGE_L(VAL1.IITEM).NEL)=J
		ENDDO
    END DO
    RETURN
ENDSUBROUTINE


SUBROUTINE SETUP_FACE_TBL()
    USE DS_Gmsh2Solver
    USE hashtbl
    IMPLICIT NONE
    INTEGER::I,J,N1(4),ET1,NFACE1,TBL_LEN,N2,K
    REAL(8)::V1(3),V2(3),NORMAL1(3),T1
	CHARACTER(LEN=:),ALLOCATABLE::KEY1
    TYPE(DICT_DATA)::VAL1
    CHARACTER(64)::CKEY1
    INTEGER::HKEY1
    
	TBL_LEN=NNODE	
    CALL FACE_TBL.INIT(TBL_LEN)
	DO I=1,NEL
		ET1=ELEMENT(I).ET
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
                V1=NODE(FACE(VAL1.IITEM).V(2)).XY-NODE(FACE(VAL1.IITEM).V(1)).XY
                V2=NODE(FACE(VAL1.IITEM).V(3)).XY-NODE(FACE(VAL1.IITEM).V(1)).XY
                call r8vec_cross_3d ( v1, v2, NORMAL1(1:3) )
                T1 = sqrt ( sum ( ( NORMAL1 )**2 ) )
                if ( T1 /= 0.0D+00 ) then
                  NORMAL1 = NORMAL1/T1
                end if
                FACE(VAL1.IITEM).UNORMAL=NORMAL1
                FACE(VAL1.IITEM).ISINI=.TRUE.
                NFACE=VAL1.IITEM
            ENDIF
            IF(FACE(VAL1.IITEM).NEL==0) ALLOCATE(FACE(VAL1.IITEM).ELEMENT(2))
            FACE(VAL1.IITEM).NEL=FACE(VAL1.IITEM).NEL+1
            IF(FACE(VAL1.IITEM).NEL>SIZE(FACE(VAL1.IITEM).ELEMENT,DIM=1)) CALL I_ENLARGE_AR(FACE(VAL1.IITEM).ELEMENT,5)
            
            IF(FACE(VAL1.IITEM).NEL>1) THEN
                DO K=1,FACE(VAL1.IITEM).SHAPE
                    IF(FACE(VAL1.IITEM).V(K)==ELEMENT(I).NODE(ELTTYPE(ELEMENT(I).ET).FACE(1,J))) THEN
                        IF(FACE(VAL1.IITEM).V(MOD(K,FACE(VAL1.IITEM).SHAPE)+1)/= &
                            ELEMENT(I).NODE(ELTTYPE(ELEMENT(I).ET).FACE(2,J))) THEN
                            FACE(VAL1.IITEM).ELEMENT(FACE(VAL1.IITEM).NEL)=-I !I单元此面的方向与face的方向反向,为-1
                        ELSE
                            FACE(VAL1.IITEM).ELEMENT(FACE(VAL1.IITEM).NEL)=I
                        ENDIF                  
                        EXIT
                    ENDIF
                ENDDO
            ELSE
                FACE(VAL1.IITEM).ELEMENT(FACE(VAL1.IITEM).NEL)=I !同向 =1
            ENDIF
            
                
            
		ENDDO
    END DO 
    
ENDSUBROUTINE     
    
    
SUBROUTINE GEN_TRISURFACE()
    USE DS_Gmsh2Solver
    USE hashtbl
    IMPLICIT NONE
    INTEGER::I,J,K,IGP1,NFACE1,IEL1,NTRI1=0,N1=0,NQUA1=0,K1=0,N2,AR1(3)
    INTEGER,ALLOCATABLE::NODE1(:),TRI1(:,:),FACE_ORDER1(:),NODEID1(:)
	REAL(8),ALLOCATABLE::XYZ1(:,:),FACENOR1(:,:)
    TYPE(DICT_DATA), ALLOCATABLE:: val1(:)
    LOGICAL::ISTRI1
	CHARACTER(512)::STLFILE1,OFFFILE1
	
    ALLOCATE(NODE1(NNODE))
    
	DO CONCURRENT (I=1:NPHGP)
		IGP1=phgpnum(i)
        NFACE1=ELTTYPE(PHYSICALGROUP(IGP1).ET_GMSH).NFACE
        NTRI1=0;NQUA1=0;N1=0     

		DO CONCURRENT (J=1:PHYSICALGROUP(IGP1).NEL)
            IEL1=PHYSICALGROUP(IGP1).ELEMENT(J)
			DO CONCURRENT (K=1:NFACE1)
                !2D
                ISTRI1=(PHYSICALGROUP(IGP1).NDIM==2.AND.PHYSICALGROUP(IGP1).ISMODEL)
                !3D
                IF(PHYSICALGROUP(IGP1).NDIM==3) THEN
                    ISTRI1=ISTRI1.OR.(FACE(ELEMENT(IEL1).FACE(K)).NEL==1)
                    IF(FACE(ELEMENT(IEL1).FACE(K)).NEL==2) ISTRI1=ISTRI1.OR. &
                        (ELEMENT(ABS(FACE(ELEMENT(IEL1).FACE(K)).ELEMENT(1))).TAG(1)/= &
                        ELEMENT(ABS(FACE(ELEMENT(IEL1).FACE(K)).ELEMENT(2))).TAG(1))
                ENDIF
                
                IF(ISTRI1) THEN
                    
                    IF(ELEMENT(ABS(FACE(ELEMENT(IEL1).FACE(K)).ELEMENT(1))).TAG(1)==IGP1)  THEN
                        N2=SIGN(1,FACE(ELEMENT(IEL1).FACE(K)).ELEMENT(1))
                    ELSE
                        N2=SIGN(1,FACE(ELEMENT(IEL1).FACE(K)).ELEMENT(2))
                    ENDIF
                    FACE(ELEMENT(IEL1).FACE(K)).ISTRISURFACE=IGP1*N2
                    
                    IF(.NOT.ALLOCATED(PHYSICALGROUP(IGP1).TRINODEG2L)) THEN
                        ALLOCATE(PHYSICALGROUP(IGP1).TRINODEG2L(NNODE))
                        PHYSICALGROUP(IGP1).TRINODEG2L=0
                    ENDIF
                    DO CONCURRENT (K1=1:FACE(ELEMENT(IEL1).FACE(K)).SHAPE)
                        IF(PHYSICALGROUP(IGP1).TRINODEG2L(FACE(ELEMENT(IEL1).FACE(K)).V(K1))==0)THEN
                            PHYSICALGROUP(IGP1).NTRINODEL2G=PHYSICALGROUP(IGP1).NTRINODEL2G+1
                            PHYSICALGROUP(IGP1).TRINODEG2L(FACE(ELEMENT(IEL1).FACE(K)).V(K1))=PHYSICALGROUP(IGP1).NTRINODEL2G
                            NODE1(PHYSICALGROUP(IGP1).NTRINODEL2G)=FACE(ELEMENT(IEL1).FACE(K)).V(K1)
                        ENDIF
                    ENDDO
                    IF(FACE(ELEMENT(IEL1).FACE(K)).SHAPE==3) THEN 
                        NTRI1=NTRI1+1
                    ELSE
                        NQUA1=NQUA1+1
                    ENDIF
                ENDIF                  
                    
            ENDDO
                
        ENDDO

        IF(NTRI1>0.OR.NQUA1>0) THEN
            ALLOCATE(PHYSICALGROUP(IGP1).TRINODEL2G(PHYSICALGROUP(IGP1).NTRINODEL2G))
            PHYSICALGROUP(IGP1).TRINODEL2G=NODE1(1:PHYSICALGROUP(IGP1).NTRINODEL2G)
            ALLOCATE(PHYSICALGROUP(IGP1).TRISURFACE(NTRI1),PHYSICALGROUP(IGP1).QUASURFACE(NQUA1))
            DO CONCURRENT (J=1:NFACE)
                IF(ABS(FACE(J).ISTRISURFACE)==IGP1) THEN   
					IF (FACE(J).ISTRISURFACE<0) THEN
						N2=-J
					ELSE
						N2=J
					ENDIF
                    IF(FACE(J).SHAPE==3) THEN
						PHYSICALGROUP(IGP1).NTRISURFACE=PHYSICALGROUP(IGP1).NTRISURFACE+1
						PHYSICALGROUP(IGP1).TRISURFACE(PHYSICALGROUP(IGP1).NTRISURFACE)=N2
					ELSE
						PHYSICALGROUP(IGP1).NQUASURFACE=PHYSICALGROUP(IGP1).NQUASURFACE+1
						PHYSICALGROUP(IGP1).QUASURFACE(PHYSICALGROUP(IGP1).NQUASURFACE)=N2
                    ENDIF
                ENDIF                    
			ENDDO
			
			IF (isstl/=0.or.isoff/=0) THEN

				IF(ALLOCATED(XYZ1)) DEALLOCATE(XYZ1)
				IF(ALLOCATED(TRI1)) DEALLOCATE(TRI1)
				IF(ALLOCATED(FACENOR1)) DEALLOCATE(FACENOR1)
				N1=PHYSICALGROUP(IGP1).NTRISURFACE+PHYSICALGROUP(IGP1).NQUASURFACE*2
				ALLOCATE(XYZ1(3,PHYSICALGROUP(IGP1).NTRINODEL2G),TRI1(3,N1),FACENOR1(3,N1))
				DO CONCURRENT (J=1:PHYSICALGROUP(IGP1).NTRINODEL2G)
					XYZ1(:,J)=NODE(PHYSICALGROUP(IGP1).TRINODEL2G(J)).XY
				ENDDO
				
				N1=0
				DO CONCURRENT (J=1:PHYSICALGROUP(IGP1).NTRISURFACE)
					N1=N1+1
					N2=ABS(PHYSICALGROUP(IGP1).TRISURFACE(J))
					IF(FACE(N2).ISTRISURFACE>0) THEN
						TRI1(:,N1)=PHYSICALGROUP(IGP1).TRINODEG2L(FACE(N2).V(1:3))
					ELSE
						TRI1(:,N1)=PHYSICALGROUP(IGP1).TRINODEG2L(FACE(N2).V(3:1:-1))
					ENDIF
					FACENOR1(:,N1)=FACE(N2).UNORMAL*SIGN(1,FACE(N2).ISTRISURFACE)
				ENDDO
				DO CONCURRENT (J=1:PHYSICALGROUP(IGP1).NQUASURFACE)
					N2=ABS(PHYSICALGROUP(IGP1).QUASURFACE(J))
					N1=N1+1
					TRI1(:,N1)=PHYSICALGROUP(IGP1).TRINODEG2L(FACE(N2).V(1:3))
									
					N1=N1+1
					TRI1(1:2,N1)=PHYSICALGROUP(IGP1).TRINODEG2L(FACE(N2).V(3:4))
					TRI1(3,N1)=PHYSICALGROUP(IGP1).TRINODEG2L(FACE(N2).V(1))
					DO K=-1,0
						IF(FACE(N2).ISTRISURFACE<0) THEN
							AR1=TRI1(3:1:-1,N1+K)
							TRI1(:,N1+K)=AR1
						ENDIF                
						FACENOR1(:,N1+K)=FACE(N2).UNORMAL*SIGN(1,FACE(N2).ISTRISURFACE)
					ENDDO
				ENDDO
				!COMPUTE FACE NORMAL
				

				if(isstl/=0) then
					!OUTPUT TO A FILE IN STL FORMAT
					STLFILE1=TRIM(FILEPATH)//'_'//TRIM(PHYSICALGROUP(IGP1).NAME)//'.STL'
					CALL stla_write (STLFILE1, PHYSICALGROUP(IGP1).NTRINODEL2G,N1, XYZ1(:,:), TRI1(:,:), FACENOR1(:,:) )
				endif
				if(isoff/=0) then
					!OUTPUT TO A FILE IN off FORMAT
                    IF(ALLOCATED(FACE_ORDER1)) DEALLOCATE(FACE_ORDER1)
                    ALLOCATE(FACE_ORDER1(N1))
                    FACE_ORDER1=3
					OFFFILE1=TRIM(FILEPATH)//'_'//TRIM(PHYSICALGROUP(IGP1).NAME)//'.OFF'
					call off_write ( XYZ1(:,:), PHYSICALGROUP(IGP1).NTRINODEL2G, PHYSICALGROUP(IGP1).NTRINODEL2G, TRI1(:,:), n1, n1, &
									FACE_ORDER1, OFFFILE1, 3 )
					
				endif				
			ENDIF
        ENDIF       

		
		
	ENDDO
	!WRITE THE WHOLE MODLE TO OFF FILE
    IF(ISOFF) THEN
	    IF (ALLOCATED(NODEID1)) DEALLOCATE(NODEID1)
        IF(ALLOCATED(TRI1)) DEALLOCATE(TRI1)
        IF(ALLOCATED(FACE_ORDER1)) DEALLOCATE(FACE_ORDER1)
        
        ALLOCATE(NODEID1(NNODE))
	    NODEID1=0
	    !PREPARE THE OUTPUT NODES AND EDGES
	    DO CONCURRENT (I=1:NFACE)
		    IF(FACE(I).ISTRISURFACE==0) CYCLE
		    NODEID1(FACE(I).V(1:FACE(I).SHAPE))=1		    		
	    ENDDO
        N1=0
        DO I=1,NNODE
            IF (NODEID1(I)/=0)THEN
                N1=N1+1
                NODEID1(I)=N1
            ENDIF
        ENDDO
	    IF(ALLOCATED(XYZ1)) DEALLOCATE(XYZ1)
        ALLOCATE(XYZ1(3,N1))
        DO CONCURRENT (I=1:NNODE)
            IF (NODEID1(I)>0)  XYZ1(:,NODEID1(I))=NODE(I).XY
        ENDDO
        !SURFACE        
	    N2=0
	    DO I=1,NFACE
		    IF (FACE(I).ISTRISURFACE/=0)  N2=N2+1
	    ENDDO
        ALLOCATE(TRI1(4,N2),FACE_ORDER1(N2))
        N2=0
        DO CONCURRENT (I=1:NFACE)
            IF (FACE(I).ISTRISURFACE/=0) THEN
                N2=N2+1
                TRI1(:,N2)=FACE(I).V
                FACE_ORDER1(N2)=FACE(I).SHAPE
            ENDIF                
        ENDDO

	    OFFFILE1=TRIM(FILEPATH)//'_WHOLEMODEL.OFF'
        call off_write ( XYZ1(:,:), N1, N1, TRI1(:,:), n2, n2, &
				FACE_ORDER1, OFFFILE1, 4 )
    ENDIF
    
	IF (ISGEO) THEN
		CALL WRITE2GMSH_GEOFILE()
	ENDIF
	
	IF(ALLOCATED(NODE1)) DEALLOCATE(NODE1)
	IF(ALLOCATED(XYZ1)) DEALLOCATE(XYZ1)
	IF(ALLOCATED(TRI1)) DEALLOCATE(TRI1)
	IF(ALLOCATED(FACENOR1)) DEALLOCATE(FACENOR1)
	IF(ALLOCATED(FACE_ORDER1)) DEALLOCATE(FACE_ORDER1)
	IF (ALLOCATED(NODEID1)) DEALLOCATE(NODEID1)
ENDSUBROUTINE




SUBROUTINE GEN_LAYERED_ELEMENT()
    USE DS_Gmsh2Solver
    IMPLICIT NONE
    INTEGER::I,J,K,N1,N2,N3,N4,IGP1,NA1(100),ALLOC_ERR,NSUBLAYER1
    REAL(8)::T1
    INTEGER,ALLOCATABLE::SUMLAYER1(:)
    REAL(8),ALLOCATABLE::AR1(:)
	character(32)::CH1
    type(element_type),allocatable::element1(:)
    type(node_type),allocatable::node1(:)
    
    ALLOCATE(SUMLAYER1(0:NLAYER))
    SUMLAYER1(0)=0
    DO I=1,NLAYER
        SUMLAYER1(I)=SUBLAYER(I)+SUMLAYER1(I-1)
    ENDDO
    NSUBLAYER1=SUMLAYER1(NLAYER)
    ALLOCATE(ELEMENT1(NSUBLAYER1*NEL),STAT=ALLOC_ERR)
 	
    DO CONCURRENT (I=1:NEL)
        N1=(I-1)*NSUBLAYER1
        N3=ELEMENT(I).NNODE
        N2=2*N3
        NA1(1:N3)=(ELEMENT(I).NODE-1)*(NSUBLAYER1+1)
        DO CONCURRENT (J=1:NLAYER)
            DO CONCURRENT (K=1:SUBLAYER(J))
            N4=N1+SUMLAYER1(J-1)+K
            ELEMENT1(N4).NNODE=N2
            ELEMENT1(N4).ILAYER=J
            ELEMENT1(N4).NTAG=ELEMENT(I).NTAG
			SELECT CASE(ELEMENT(I).ET)
				CASE(1)	
					ELEMENT1(N4).ET=3					
				CASE(2)
					ELEMENT1(N4).ET=6
				CASE(3)
					ELEMENT1(N4).ET=5
				CASE DEFAULT
					STOP "ETTYPE UNEXPECTED. SUB=GEN_LAYERED_ELEMENT()."
			ENDSELECT
            ALLOCATE(ELEMENT1(N4).NODE(N2),STAT=ALLOC_ERR)
            ALLOCATE(ELEMENT1(N4).TAG,SOURCE=ELEMENT(I).TAG,STAT=ALLOC_ERR)
			IF(J>1) THEN
				ELEMENT1(N4).TAG(1)=ELEMENT(I).TAG(1)*100+J				
				ELEMENT1(N4).TAG(2)=ELEMENT(I).TAG(2)*100+J								
			ENDIF
			PHYSICALGROUP(ELEMENT(I).TAG(1)).LAYERGROUP(J)=ELEMENT1(N4).TAG(1)
			
            ELEMENT1(N4).NODE(1:N3)=NA1(1:N3)+SUMLAYER1(J-1)+K
            IF(ELEMENT1(N4).ET/=3) THEN
                ELEMENT1(N4).NODE(N3+1:N2)=ELEMENT1(N4).NODE(1:N3)+1
            ELSE
                ELEMENT1(N4).NODE(N2:N3+1:-1)=ELEMENT1(N4).NODE(1:N3)+1 !QUA ELEMENT.
            ENDIF
             
            ENDDO
        ENDDO
    ENDDO
    
    DEALLOCATE(ELEMENT)
    ALLOCATE(ELEMENT,SOURCE=ELEMENT1,STAT=ALLOC_ERR)
    NEL=NEL*NSUBLAYER1
    DEALLOCATE(ELEMENT1)
    
    !UPDATE NODE
    ALLOCATE(NODE1(NNODE*(NSUBLAYER1+1)),STAT=ALLOC_ERR)
    ALLOCATE(AR1(NSUBLAYER1),STAT=ALLOC_ERR)
    
    DO CONCURRENT (I=1:NNODE)
        N1=(I-1)*(NSUBLAYER1+1)        
        DO CONCURRENT (J=1:NLAYER)
            T1=(ELEVATION(J+1,I)-ELEVATION(J,I))/SUBLAYER(J)
            N3=SUBLAYER(J)
            IF(J==NLAYER) N3=N3+1
            DO CONCURRENT (K=1:N3)
                N2=N1+SUMLAYER1(J-1)+K
                NODE1(N2)=NODE(I)            
                NODE1(N2).XY(3)=ELEVATION(J,I)+T1*(K-1)
            ENDDO
        ENDDO 
    ENDDO
    
    DEALLOCATE(NODE)
    ALLOCATE(NODE,SOURCE=NODE1,STAT=ALLOC_ERR)
    NNODE=NNODE*(NSUBLAYER1+1)
    DEALLOCATE(NODE1)
    DEALLOCATE(ELEVATION)
    
	!update physical group
	!统计每个物理组单元的个数及单元数组下标。
	N3=0
	do i=1,nphgp
		IGP1=phgpnum(i)
		physicalgroup(IGP1).NDIM=physicalgroup(IGP1).NDIM+1
        
		SELECT CASE(physicalgroup(IGP1).ET_GMSH)
			CASE(1)	
				physicalgroup(IGP1).ET_GMSH=3					
			CASE(2)
				physicalgroup(IGP1).ET_GMSH=6
			CASE(3)
				physicalgroup(IGP1).ET_GMSH=5
			CASE DEFAULT
				STOP "ET_GMSH UNEXPECTED. SUB=GEN_LAYERED_ELEMENT()."
		ENDSELECT		
		DO J=2,NLAYER
			N2=IGP1*100+J
			physicalgroup(N2)=physicalgroup(IGP1)
            IF(SUBLAYER(J)>1) THEN
                IF(ALLOCATED(physicalgroup(N2).ELEMENT)) DEALLOCATE(physicalgroup(N2).ELEMENT)
                physicalgroup(N2).NEL=physicalgroup(IGP1).NEL*SUBLAYER(J)
                ALLOCATE(physicalgroup(N2).ELEMENT(physicalgroup(N2).NEL),STAT=ALLOC_ERR)
            ENDIF
            
			physicalgroup(N2).MAT(1)=physicalgroup(IGP1).MAT(J)
			WRITE(CH1,'(I3)') J
			physicalgroup(N2).name=TRIM(physicalgroup(IGP1).NAME)//'_L'//TRIM(ADJUSTL(CH1))
			N3=N3+1
			PHGPNUM(NPHGP+N3)=N2            
		ENDDO
		physicalgroup(IGP1).name=TRIM(physicalgroup(IGP1).NAME)//'_L1'
        IF(SUBLAYER(1)>1) THEN
            DEALLOCATE(physicalgroup(IGP1).ELEMENT)
            physicalgroup(IGP1).NEL=physicalgroup(IGP1).NEL*SUBLAYER(1)
            ALLOCATE(physicalgroup(IGP1).ELEMENT(physicalgroup(IGP1).NEL),STAT=ALLOC_ERR)
        ENDIF
	end do
	NPHGP=NPHGP+N3
	
	physicalgroup(PHGPNUM(1:NPHGP)).nel=0
	do i=1,nel
		physicalgroup(element(i).tag(1)).nel=physicalgroup(element(i).tag(1)).nel+1
		physicalgroup(element(i).tag(1)).element(physicalgroup(element(i).tag(1)).nel)=i		
	end do
	
END SUBROUTINE

    
SUBROUTINE GEN_ZEROTHICKNESS_ELEMENT()
	USE DS_Gmsh2Solver
	USE HASHTBL
    USE IFPORT
	IMPLICIT NONE
	INTEGER::I,J,K,IPGROUP1,NNODE1,NODE1(50),NVAL1,N1,N2,IEL1,IEL2
	CHARACTER(LEN=:),ALLOCATABLE::KEY1
	TYPE(DICT_DATA),ALLOCATABLE::VAL1(:)
	REAL(8)::XY1(3)
	
	DO I=1,nphgp
		IPGROUP1=PHGPNUM(I)
		CALL LOWCASE(physicalgroup(IPGROUP1).ET,LEN_TRIM(physicalgroup(IPGROUP1).ET))
		IF(physicalgroup(IPGROUP1).ET=='zt_cpe4_spg'.OR.physicalgroup(IPGROUP1).ET=='zt_prm6_spg') THEN
			
			!IF(.NOT.EDGE_TBL.IS_INIT) CALL SETUP_EDGE_TBL()
			!IF(physicalgroup(IPGROUP1).NDIM/=1) STOP "DIMENSION ERROR.ERRORS IN SUB GEN_ZEROTHICKNESS_ELEMENT()."
			!IF(physicalgroup(IPGROUP1).ET_GMSH/=1) STOP "ET ERROR.ERRORS IN SUB GEN_ZEROTHICKNESS_ELEMENT()."
			IF(MOD(physicalgroup(IPGROUP1).NEL,2)/=0) STOP "ELEMENT NUMBER ERROR.ERRORS IN SUB GEN_ZEROTHICKNESS_ELEMENT()."
				
			N1=physicalgroup(IPGROUP1).NEL/2
				
			DO J=1,N1
				IEL1=physicalgroup(IPGROUP1).ELEMENT(J)
				IEL2=physicalgroup(IPGROUP1).ELEMENT(N1+J)
                NODE1(1:ELEMENT(IEL1).NNODE)=ELEMENT(IEL1).NODE
				NODE1(ELEMENT(IEL1).NNODE+1:ELEMENT(IEL1).NNODE+ELEMENT(IEL2).NNODE)=ELEMENT(IEL2).NODE
				ELEMENT(IEL1).NNODE=ELEMENT(IEL1).NNODE*2
				DEALLOCATE(ELEMENT(IEL1).NODE)
				ALLOCATE(ELEMENT(IEL1).NODE,SOURCE=NODE1(1:ELEMENT(IEL1).NNODE))
			ENDDO
				
			physicalgroup(IPGROUP1).NEL=N1
				
			!RAMDOM CHECK MATCH
			DO J=1,10
				N2=MOD(IRAND(1),physicalgroup(IPGROUP1).NEL)+1
				N2=physicalgroup(IPGROUP1).ELEMENT(N2)
				DO K=1,ELEMENT(N2).NNODE/2
					XY1=NODE(ELEMENT(N2).NODE(K)).XY-NODE(ELEMENT(N2).NODE(K+ELEMENT(N2).NNODE/2)).XY
					IF(NORM2(XY1)>1E-6) STOP "MATCH ERROR IN SUB GEN_ZEROTHICKNESS_ELEMENT()."
				ENDDO
			ENDDO
			!REORDER NODAL NUMBER zt_cpe4_spg
			IF(physicalgroup(IPGROUP1).ET=='zt_cpe4_spg') THEN
				DO J=1,physicalgroup(IPGROUP1).NEL
					N2=ELEMENT(physicalgroup(IPGROUP1).ELEMENT(J)).NODE(3)
					ELEMENT(physicalgroup(IPGROUP1).ELEMENT(J)).NODE(3)=ELEMENT(physicalgroup(IPGROUP1).ELEMENT(J)).NODE(4)
					ELEMENT(physicalgroup(IPGROUP1).ELEMENT(J)).NODE(4)=N2
				ENDDO
			ENDIF
		ENDIF
	ENDDO
	
ENDSUBROUTINE

SUBROUTINE ENLARGE_NODE()
	USE DS_Gmsh2Solver
	IMPLICIT NONE
	TYPE(node_type),ALLOCATABLE::NODE1(:)
	INTEGER::I
	TNNODE=NNODE+1000
	ALLOCATE(NODE1(TNNODE))
	NODE1(1:NNODE)=NODE
	DEALLOCATE(NODE)
	ALLOCATE(NODE,SOURCE=NODE1)
	DEALLOCATE(NODE1)
ENDSUBROUTINE




subroutine elt_bc_load_translate()
	use DS_Gmsh2Solver
	implicit none
	integer::i,j,k,n1,n2,nc,iar1(5)
	real(8)::LAV1
	integer,allocatable::nodalload1(:),node1(:)
	real(8),allocatable::load1(:)
	
	! 荷载具可叠加，边界条件不具有可加性。分布力转化为节点力时要积分，但分布位移转化为节点位移时不需积分，节点位移等于分布位移。
	if(nelt_load>0) allocate(nodalLoad(maxnnodalLoad))
	if(nelt_bc>0) allocate(nodalBC(maxnnodalBC))
	if(nelt_spgface>0) allocate(spgface(maxnspgface))
	allocate(nodalload1(nnode),node1(nnode))
	if(nelt_load>0) allocate(load1(nnode))
	do i=1, nelt_load
		select case(elt_load(i).ndim)
			case(0) !ndim=0,是点荷载，各个单元中的同一节点，不具可加性，均指同一值。
				nodalload1=0
				do j=1,physicalgroup(elt_load(i).group).nel
					n1=physicalgroup(elt_load(i).group).element(j)
					do k=1,element(n1).nnode
						n2=element(n1).node(k)
						if(nodalload1(n2)==0) then
							nnodalLoad=nnodalLoad+1 
							if(nnodalLoad>maxnnodalLoad) call enlargenodalLoad()
							nodalLoad(nnodalLoad).node=n2
							nodalLoad(nnodalLoad).dof=elt_load(i).dof
							nodalLoad(nnodalLoad).sf=elt_load(i).sf
							nodalLoad(nnodalLoad).value=elt_load(i).value
							nodalload1(n2)=1
						end if
					end do
				end do
								
			case default
				load1=0
				do j=1,physicalgroup(elt_load(i).group).nel
					n1=physicalgroup(elt_load(i).group).element(j)
					call Element_LAV_Cal(n1,LAV1) 
					LAV1=LAV1*elt_load(i).value
					do k=1,element(n1).nnode
						n2=element(n1).node(k)
						load1(n2)=load1(n2)+Elttype(element(n1).et).weight(k,1)*LAV1 !均布荷载		 				
					end do					
				end do				
				do j=1,nnode
					if(abs(load1(j))>1e-10) then
						nnodalLoad=nnodalLoad+1 
						if(nnodalLoad>maxnnodalLoad) call enlargenodalLoad()
						nodalLoad(nnodalLoad).node=j
						nodalLoad(nnodalLoad).dof=elt_load(i).dof
						nodalLoad(nnodalLoad).sf=elt_load(i).sf
						nodalLoad(nnodalLoad).value=load1(j)
					end if
				end do

		end select
		
	end do
	
	do i=1, nelt_bc
		!**********目前认为边界条件都没有可加性,NDIM=1,2,3的情况与NDIM=0的情况一样**********。
		nodalload1=0
		do j=1,physicalgroup(elt_bc(i).group).nel
			n1=physicalgroup(elt_bc(i).group).element(j)
			do k=1,element(n1).nnode
				n2=element(n1).node(k)
				if(nodalload1(n2)==0) then
					nnodalBC=nnodalBC+1 
					if(nnodalBC>maxnnodalBC) call enlargenodalBC()
					nodalBC(nnodalBC).node=n2
					nodalBC(nnodalBC).dof=elt_bc(i).dof
					nodalBC(nnodalBC).sf=elt_bc(i).sf
					nodalBC(nnodalBC).value=elt_bc(i).value
					nodalload1(n2)=1 !多次出现时以第一次出现的值为准。
				end if
			end do			
		end do		
	end do 
	
	do i=1, nelt_spgface
		!**********目前认为边界条件都没有可加性**********。
		nodalload1=0
		elt_spgface(i).n1=nspgface+1
		do j=1,physicalgroup(elt_spgface(i).group).nel
			n1=physicalgroup(elt_spgface(i).group).element(j)
			do k=1,element(n1).nnode
				n2=element(n1).node(k)
				if(nodalload1(n2)==0) then
					nspgface=nspgface+1 
					if(nspgface>maxnspgface) call enlargespgface()
					spgface(nspgface).node=n2
					spgface(nspgface).dof=elt_spgface(i).dof
					spgface(nspgface).sf=elt_spgface(i).sf
					spgface(nspgface).value=elt_spgface(i).value
					nodalload1(n2)=1 !多次出现时以第一次出现的值为准。
				end if
			end do			
		end do
		elt_spgface(i).n2=nspgface
	end do 
    
    do i=1,nwsp
        nodalload1=0
        !node1=0
        do j=1,physicalgroup(wsp(i).group).nel
            n1=physicalgroup(wsp(i).group).element(j)
            do k=1,element(n1).nnode
                n2=element(n1).node(k)
                if(nodalload1(n2)==0)then
                    wsp(i).nnode= wsp(i).nnode+1
                    node1(wsp(i).nnode)=n2
                    nodalload1(n2)=1
                end if                
            end do
        end do
		
        wsp(i).nchnode=physicalgroup(wsp(i).chgroup).nel
        allocate(wsp(i).node(wsp(i).nnode),wsp(i).chnode(wsp(i).nchnode))
		
		call SortByPath(wsp(i).group,wsp(i).spgroup,wsp(i).node,wsp(i).nnode)
        !do j=1,wsp(i).nnode
        !    
        !    do k=j+1,wsp(i).nnode
        !        !从小到大
        !        if(node(node1(k)).xy(1)<node(node1(j)).xy(1)) then
        !            n1=node1(k)
        !            node1(k)=node1(j)
        !            node1(j)=n1
        !        end if
        !    end do
        !end do
		
		
		
        
        !if(wsp(i).xdirection==1) then
        !    wsp(i).node=node1(1:wsp(i).nnode)
        !else
        !    wsp(i).node=node1(wsp(i).nnode:1:-1)
        !endif
        
		do j=1,wsp(i).nchnode
			do k=1,wsp(i).nnode
				if(wsp(i).node(k)==element(physicalgroup(wsp(i).chgroup).element(j)).node(1)) then
					wsp(i).chnode(j)=k
					exit
				end if
			end do
        end do
		
        do j=1,wsp(i).nchnode
            do k=j+1,wsp(i).nchnode
               if(wsp(i).chnode(j)>wsp(i).chnode(k)) then
                   n1=wsp(i).chnode(j)
                   wsp(i).chnode(j)=wsp(i).chnode(k)
                   wsp(i).chnode(k)=n1
               end if
            end do
        end do
		
    end do
	
    do i=1,nDataPoint
        nodalload1=0
        !node1=0
        do j=1,physicalgroup(DataPoint(i).group).nel
            n1=physicalgroup(DataPoint(i).group).element(j)
            do k=1,element(n1).nnode
                n2=element(n1).node(k)
                if(nodalload1(n2)==0)then
                    DataPoint(i).nnode= DataPoint(i).nnode+1
                    node1(DataPoint(i).nnode)=n2
                    nodalload1(n2)=1
                end if                
            end do
        end do
		
        allocate(DataPoint(i).node(DataPoint(i).nnode))
		
		IF(DATAPOINT(I).ORDER==0 .or. physicalgroup(DataPoint(i).group).ndim/=1) THEN
			DataPoint(i).node(1:DataPoint(i).nnode)=NODE1(1:DataPoint(i).nnode)
		ELSE
			call SortByPath(DataPoint(i).group,DataPoint(i).spgroup,DataPoint(i).node,DataPoint(i).nnode)
		END IF		

		
    end do	
	
	if(allocated(nodalload1)) deallocate(nodalload1)
	if(allocated(load1)) deallocate(load1)
	if(allocated(node1)) deallocate(node1)
end subroutine


subroutine SortByPath(igroup,ispgroup,Local_Node,NLNDE)
	use DS_Gmsh2Solver
	implicit none
	integer,intent(in)::igroup,ispgroup,NLNDE
	integer,intent(out)::Local_node(NLNDE)
	integer::i,j,k,iar1(5),n1,n2
	integer,allocatable::element1(:)
	
	
	allocate(element1(physicalgroup(igroup).nel))
	element1=0
	Local_Node(1)=element(physicalgroup(ispgroup).element(1)).node(1)	
	j=1	
	do while(j<=NLNDE-1)
		do k=1,physicalgroup(igroup).nel
			if(element1(k)/=0) cycle
			
			n1=element(physicalgroup(igroup).element(k)).nnode
			iar1(1:n1)=element(physicalgroup(igroup).element(k)).node(1:n1)
			n1=iar1(1)-Local_Node(j)
			n2=iar1(2)-Local_Node(j)
			if(n1*n2/=0) cycle
			
			element1(k)=1
			
			if(n1==0) then
				select case(physicalgroup(igroup).ET_GMSH)
				case(1)
					Local_Node(j+1)=iar1(2)
					j=j+1
				case(8)
					Local_Node(j+1)=iar1(3)
					Local_Node(j+2)=iar1(2)
					j=j+2
				case(27)
					Local_Node(j+1)=iar1(3)
					Local_Node(j+2)=iar1(4)
					Local_Node(j+3)=iar1(5)
					Local_Node(j+4)=iar1(2)
					j=j+4
				case default
					print *, "The element is expected to be a line element. datapoint(i),i=",i
					stop
				end select
			else
				select case(physicalgroup(igroup).ET_GMSH)
				case(1)
					Local_Node(j+1)=iar1(1)
					j=j+1
				case(8)
					Local_Node(j+1)=iar1(3)
					Local_Node(j+2)=iar1(1)
					j=j+2
				case(27)
					Local_Node(j+1)=iar1(5)
					Local_Node(j+2)=iar1(4)
					Local_Node(j+3)=iar1(3)
					Local_Node(j+4)=iar1(1)
					j=j+4
				case default
					print *, "The element is expected to be a line element. datapoint(i),i=",i
					stop
				end select						
			end if
			exit
		end do
	enddo	

	deallocate(element1)

end subroutine

subroutine Element_LAV_Cal(ienum,LAV) !计算单元的长度、面积或体积
	use DS_Gmsh2Solver
	implicit none	
	integer,intent(in)::ienum
	real(8),intent(out)::LAV
	REAL(8),EXTERNAL::TRIAREA,TETVOL
	integer::i,j
	real(8)::a1(4,3)=0
	
	select case(element(ienum).et)
		case(1,8,27) !Length for a Line element 
			LAV=((node(element(ienum).node(1)).xy(1)-node(element(ienum).node(2)).xy(1))**2+ &
			(node(element(ienum).node(1)).xy(2)-node(element(ienum).node(2)).xy(2))**2+ &
			(node(element(ienum).node(1)).xy(3)-node(element(ienum).node(2)).xy(3))**2)**0.5
		case(2,9,23) !三角形面积
						
			do i=1,3
				a1(i,:)=node(element(ienum).node(i)).xy
			end do
			LAV=TRIAREA(A1(1:3,1:3))
		case(3,16) !四边形面积
			do i=1,4
				a1(i,:)=node(element(ienum).node(i)).xy
			end do
			LAV=TRIAREA(A1(1:3,1:3))
			
			A1(2,:)=A1(1,:)
			LAV=LAV+TRIAREA(A1(2:4,1:3))
			
		case(4,11) !四面体体积	
			do i=1,4
				a1(i,:)=node(element(ienum).node(i)).xy
			end do
			LAV=TETVOL(A1(1:4,1:3))
		
		case(5,6)
			print *, "Still under work."
			pause
		case default
			print *, "Still under work."
			pause
		
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
	TriArea=0.5*TriArea**0.5

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


!Enlarge Array ELT() by an increment of 1000
subroutine enlargenodalLoad()
	use DS_Gmsh2Solver
	implicit none
	
	type(bc_tydef)::nodalLoad1(maxnnodalLoad)
	
	nodalLoad1=nodalLoad(1:maxnnodalLoad)
	deallocate(nodalLoad)	
	allocate(nodalLoad(maxnnodalLoad+1000))
	nodalLoad(1:maxnnodalLoad)=nodalLoad1
	maxnnodalLoad=maxnnodalLoad+1000

end subroutine


subroutine enlargenodalBC()
	use DS_Gmsh2Solver
	implicit none

	type(bc_tydef)::nodalBC1(maxnnodalBC)
	
	nodalBC1=nodalBC(1:maxnnodalBC)
	deallocate(nodalBC)	
	allocate(nodalBC(maxnnodalBC+1000))
	nodalBC(1:maxnnodalBC)=nodalBC1
	maxnnodalBC=maxnnodalBC+1000

end subroutine

subroutine enlargespgface()
	use DS_Gmsh2Solver
	implicit none

	type(bc_tydef)::spgface1(maxnspgface)
	
	spgface1=spgface(1:maxnspgface)
	deallocate(spgface)	
	allocate(spgface(maxnspgface+1000))
	spgface(1:maxnspgface)=spgface1
	maxnspgface=maxnspgface+1000

end subroutine

SUBROUTINE find_duplicate_nodes(NODE,DIM,NNODE,NEWNODE,NODE2NN,NNEW,TOL)
	IMPLICIT NONE
	INTEGER,INTENT(IN)::DIM,NNODE
	REAL(8),INTENT(IN)::NODE(DIM,NNODE),TOL
	INTEGER::NODE2NN(NNODE),NEWNODE(NNODE),NNEW 
	!NODE2NN,INDEX THE NODE(J)=NEWNODE(NODE2NN(J))
	!NEWNODE:NODES WITHOUT DUPLICATES,INDEX,NEWNODE(I)=NODE(NEWNODE(I))
	INTEGER::I,J
	
	
	NODE2NN=0;NEWNODE=0;NNEW=0
	DO I=1,NNODE
		DO J=I+1,NNODE
			IF(NODE2NN(J)>0) CYCLE
			IF(ABS(NODE(1,I)-NODE(1,J))>TOL) CYCLE
			IF(DIM>1) THEN
                IF(ABS(NODE(2,I)-NODE(2,J))>TOL) CYCLE
            ENDIF
			IF(DIM>2) THEN
                IF(ABS(NODE(3,I)-NODE(3,J))>TOL) CYCLE
            ENDIF                
			NODE2NN(J)=I !MEAN NODE(J)=NODE(I)
		ENDDO
		IF(NODE2NN(I)==0) THEN
		!NO DUPLICATE NODE.
			NNEW=NNEW+1
			NEWNODE(NNEW)=I
			NODE2NN(I)=NNEW
		ELSE
			NODE2NN(I)=NODE2NN(NODE2NN(I))
		ENDIF
		
	ENDDO	
	
ENDSUBROUTINE

SUBROUTINE FIND_DUPLICATE_FACES(FACE,NORDER,NFACE,NF2FACE,FACE2NF,NNF)
!only can find the elements have identical noded id.
	IMPLICIT NONE
	INTEGER,INTENT(IN)::NFACE,NORDER
	INTEGER,INTENT(IN)::FACE(NORDER,NFACE)
	INTEGER::NF2FACE(NFACE),FACE2NF(NFACE),NNF
	INTEGER,ALLOCATABLE::FACE1(:,:)
	INTEGER::I,J,K,N1,N2
	
	ALLOCATE(FACE1,SOURCE=FACE)
	
	NNF=0;NF2FACE=0;FACE2NF=0
	!SORT FACE NODE
	DO CONCURRENT (I=1:NFACE)
		N1=3
		IF (FACE1(4,I)>0) N1=4		
		DO J=1,N1-1
			DO K=J+1,N1
				IF(FACE1(K,I)<FACE1(J,I)) THEN
					N2=FACE1(J,I);FACE1(J,I)=FACE1(K,I);FACE1(K,I)=N2
				ENDIF
			ENDDO
		ENDDO		
	END DO    
	
	DO I=1,NFACE
		DO J=I+1,NFACE			
			IF(FACE2NF(J)>0) CYCLE
			IF(FACE1(1,I)-FACE1(1,J)/=0) CYCLE
			IF(FACE1(2,I)-FACE1(2,J)/=0) CYCLE
			IF(FACE1(3,I)-FACE1(3,J)/=0) CYCLE
			IF(FACE(4,I)>0.AND.FACE1(4,J)>0) THEN
				IF(FACE1(4,I)-FACE1(4,J)/=0) CYCLE
			ENDIF
			
			FACE2NF(J)=I !MEAN FACE1(J)=FACE1(I)
		ENDDO
		IF(FACE2NF(I)==0) THEN
		!NO DUPLICATE FACE1.
			NNF=NNF+1
			NF2FACE(NNF)=I
			FACE2NF(I)=NNF
		ELSE
			FACE2NF(I)=FACE2NF(FACE2NF(I))
		ENDIF		
	ENDDO
	
	IF(ALLOCATED(FACE1)) DEALLOCATE(FACE1)
	
ENDSUBROUTINE


SUBROUTINE ET_EDGE_FACE()
	USE DS_Gmsh2Solver
	IMPLICIT NONE
	INTEGER::I,ET,AET1(33)
	
    AET1=[1:31,92,93]
	DO I=1,33
        ET=AET1(I)
		SELECT CASE(ET)
			CASE(1,8,26,27,28) !LINE
				Elttype(ET).NEDGE=1;Elttype(ET).NFACE=0
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,1)=[1,2]
			CASE(2,9,20,21,22,23,24,25) !TRIANGLE
				Elttype(ET).NEDGE=3;Elttype(ET).NFACE=1
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,1],(/2,3/))
				Elttype(ET).FACE(:,1)=[3,1,2,3,0]
				Elttype(ET).FACEEDGE=Elttype(ET).FACE
			CASE(3,10,16) !QUADRANGLE
				Elttype(ET).NEDGE=4;Elttype(ET).NFACE=1
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE),&
						 Elttype(ET).FACEEDGE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,4,4,1],(/2,4/))
				Elttype(ET).FACE(:,1)=[4,1,2,3,4]
				Elttype(ET).FACEEDGE=Elttype(ET).FACE
			CASE(4,11,29,30,31) !TETRAHEDRON
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
			CASE(5,12,17,92,93) !HEXAHEDRON
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
			CASE(6,13,18) !6-NODE PRISM
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
			
			CASE(7,14,19) !5-node pyramid
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
		ENDSELECT
	ENDDO
	
ENDSUBROUTINE




   !n1为node()数组中全域节点总数
   subroutine reorder_nodal_number(IPERM,nnum,adjL,maxadj)

	  implicit none
	  integer,intent(in)::nnum,maxadj
	  integer,intent(in)::adjL(maxadj,nnum)
	  integer,intent(inout)::IPERM(nnum)
	  integer::i,j,k,nnz,ITYPE,IFLAG,LP,MP,n1,n2,n3
	  REAL(8)::IPROF(2)
	  integer,allocatable::IRN(:),JCN(:),IW(:),ICPTR(:),art(:)
	  integer::ar(15),allo_err
	  real(8)::t1
	  COMMON /MC40I/ LP,MP

 

!	  allocate(IPERM(nnum))
	  allocate(ICPTR(nnum+1))
	  allocate(IW(3*nnum+2))

 	  
	  !统计adjL()中非零元素的个数。
	  NNZ=count(adjL>0)
	  allocate(IRN(2*NNZ))
	  allocate(JCN(NNZ))
	
	  n2=0
	  do i=2,nnum
		 do j=1,maxadj
			if(adjL(j,i)/=0) then
			   n2=n2+1
			   IRN(n2)=i
			   jCN(n2)=adjL(j,i)
			else
			   exit
			end if
		 end do	
	  end do
	
	  ITYPE=1
	  ICPTR=0
	  IW=0
	  IPERM=0

	  CALL MC40AD(ITYPE,nnum,NNZ,IRN,JCN,ICPTR,IPERM,IW,IPROF,IFLAG)
	  ! Check for an error return
	  if (IFLAG.LT.0) then
		 print *, 'S.W. Sloan method failed. Try another normal method to reordering nodal number. Please Wait...'
		 !call nodecoding(nnum)
		 return
	  end if
!	  C Write out the profile
	  WRITE (6,210) IPROF(1)
	  IF (IPROF(1).EQ.IPROF(2))THEN
		 WRITE (6,200)
	  ELSE
		 WRITE (6,220) IPROF(2)
	  END IF
	  !pause

!
	  200 FORMAT(/5X,'The algorithm did not reduce the profile')
	  210 FORMAT(/5X,'The profile initially is',F15.0)
	  220 FORMAT(/5X,'The reduced profile is',F15.0)


	  !call csb()

	  !deallocate(adjl)
	  deallocate(IRN)
	  deallocate(JCN)	
	  deallocate(IW,stat=allo_err)
	  !deallocate(IPERM)
	  deallocate(ICPTR,stat=allo_err)

	  return	
   end subroutine
   
   
   	!对于每个节点，存储和其相邻(定义为总刚中，这两个节点对应的元素不为零。)且编号小于该节点自身编号的节点编号
	!这里假定一个节点最多有maxadj个与之相邻且编号小于该节点自身编号的节点。
	
   subroutine setup_adjList()
		use DS_Gmsh2Solver
		implicit none
		integer::i,j,iel,n1,ar(200),n2
		
		allocate(adjL(maxadj,nnode))
		adjL=0
		do iel=1,nel
	  	
	  		n1=element(iel).nnode
	  		ar(1:n1)=element(iel).node
	  		!从大到小
  			do i=1,n1
			  do j=i+1,n1
				 if(ar(i)<ar(j)) then
					n2=ar(i)
					ar(i)=ar(j)
					ar(j)=n2
				 end if
			  end do
			end do
			
			do i=1,n1-1
			  do j=i+1,n1
				 call addtoadjL(ar(j),adjL(:,ar(i)),maxadj)
			  end do
			end do					

		end do
   end subroutine

   !如果iar(:)中没有n1,则把n1加入到iar(:)中。
   subroutine addtoadjL(n1,iar,maxadj)
	  implicit none
	  integer,intent(in)::n1,maxadj
	  integer::iar(maxadj),i
	
		  
	  if(any(iar==n1)) return
	
	  if(any(iar==0)) then
		 do i=1,maxadj
		    if(iar(i)==0) then
			   iar(i)=n1
			   exit
			end if
		 end do
	  else
		 print *, '与节点相邻且编号小于该节点自身编号的节点数在>maxadj,adjL的原定空间不足.'
		 stop
	  end if

	
    end subroutine
    
logical function isacw(x1,y1,z1,x2,y2,z2,x3,y3,z3)
	implicit none
	real(8)::x1,y1,z1,x2,y2,z2,x3,y3,z3
    real(8)::t1
	
	isacw=.false.
	y2=y2-y1
	x2=x2-x1
    z2=z2-z1
	y3=y3-y1
	x3=x3-x1
    z3=z3-z1
	t1=(x2*y3-y2*x3)+(y2*z3-z2*y3)-(x2*z3-z2*x3)
	if(t1>0) isacw=.true.

end function
