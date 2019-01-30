module geomodel
	USE IFQWIN 
    use meshDS, only : node,nnode,elt,nelt,edge,nedge,adjlist,zone,znum, &
						bnode,nbnode,cedge,ncedge,soillayer,&
						path_name,ENUMBER,xyscale,xmin,ymin
    implicit none
    INTEGER,ALLOCATABLE::BEDGE(:) !区边界inode对应的每个节点生成的edge的整体编号,
contains

    subroutine out_volume_to_gmsh()
		character*256 term
		integer(4)::msg
		character(512)::outfile
		
        call st_membrance()
        call node_on_boundary()
        call find_zone_boudary()
        CALL CHECK_ORIENT()
		outfile=trim(path_name)//'_geomodel.geo'
        call WRITE2GMSH_GEOFILEi(outfile)
		
		term="THE VOLUME GEO FILE HAS BEEN OUTPUT! Click Yes to Exit,No to Continue"
	 	term=trim(term)
     	msg = MESSAGEBOXQQ(trim(term),'COMPLETED'C,MB$ICONINFORMATION.OR.MB$YESNO.OR.MB$DEFBUTTON1)
     	if(msg==MB$IDYES) msg=clickmenuqq(loc(WINEXIT))		
		
    endsubroutine
    
	subroutine node_on_boundary()
		integer::i
		do i=1,ncedge
			node(cedge(I).v).onbdy=1		
		enddo
		nbnode=count(node.onbdy==1)
		allocate(bnode(nbnode),BEDGE(NNODE))
		bnode=pack(node(1:nnode).number,node(1:nnode).onbdy==1)
	end subroutine
	
    subroutine find_zone_boudary()
		integer::izone
        integer::i,j,k,ielt,iedge,n1
		integer::bedge1(1000),nbe=0
		
        do izone=1,znum
		    nbe=0
            n1=zone(izone).ntrie3n
		    do i=1,zone(izone).ntrie3n
			    ielt=zone(izone).trie3n(i)			
			    do j=1,3
				    iedge=elt(ielt).edge(j)
				    if(edge(iedge).e(1)+1==0.or.edge(iedge).e(2)+1==0) then
					    nbe=nbe+1
					    bedge1(nbe)=iedge
				    else
					    if(elt(edge(iedge).e(1)).zn/=elt(edge(iedge).e(2)).zn) then
						    nbe=nbe+1
						    bedge1(nbe)=iedge
					    endif
				    endif
			    enddo
		    enddo
		    if(nbe>0) then
			    allocate(zone(izone).bedge,source=bedge1(1:nbe))
                ZONE(IZONE).NBE=NBE
		    endif
	    enddo
    endsubroutine
	
	SUBROUTINE WRITE2GMSH_GEOFILEi(FILE_L)
		
		IMPLICIT NONE

		CHARACTER(LEN=*),INTENT(IN)::FILE_L
		
		INTEGER::I,J,K,K1,N1,N2,N3,N4,NV1,INC,LEN1,INC2,AN1(4)
		INTEGER::NEDGE1=0
		REAL(8)::T1,AT1(3)
		REAL(8),ALLOCATABLE::MINDIS1(:)
		!INTEGER,ALLOCATABLE::NODEID1(:),EDGEID1(:),FACEID1(:) !LOCAL IDS FOR NODES AND EDGES 
		CHARACTER(512)::GEOFILE1
		CHARACTER(32)::CH1,CH2,CH3		
		CHARACTER(4096)::STR1
		

		
		OPEN(UNIT=50,FILE=FILE_L,STATUS='REPLACE')
		WRITE(50,10)

		!MINIMAL DISTANCE BETWEEN NODES
		ALLOCATE(MINDIS1(NNODE))
		MINDIS1=1.D20
		DO I=1,NNODE			
			DO J=I+1,NNODE				
				AT1=[(NODE(I).X-NODE(J).X)*XYSCALE,(NODE(I).Y-NODE(J).Y)*XYSCALE,NODE(I).Z-NODE(J).Z]
				T1=MAX(NORM2(AT1),0.1D0) !SET MINDIS>0.1
				IF (MINDIS1(I)>T1) MINDIS1(I)=T1
				IF (MINDIS1(J)>T1) MINDIS1(J)=T1
			ENDDO
		ENDDO
		
		DO J=0,SOILLAYER
			N1=NNODE*J
			DO I=1,NNODE
				N2=N1+I
				INC=INCOUNT(N2)
				WRITE(50,20) N2,NODE(N2).X*XYSCALE+XMIN,NODE(N2).Y*XYSCALE+YMIN,NODE(N2).Z,MINDIS1(I)
			ENDDO
		ENDDO
		
		N1=0;NEDGE1=0
		DO J=0,SOILLAYER
			N2=NNODE*J
			DO I=1,NEDGE
				IF(edge(i).v(1)*edge(i).v(2)==0) cycle
				IF(J==0) THEN
					NEDGE1=NEDGE1+1					
				ENDIF
				N1=N1+1
				INC=INCOUNT(N1)
				WRITE(50,30), N1,EDGE(I).V+N2			
			ENDDO			
		ENDDO
		!VERTICAL BOUNDARY EDGE 
		DO J=1,SOILLAYER
			N3=NNODE*(J-1);N4=NNODE*J
			DO I=1,NBNODE
				N1=N1+1
				IF(J==1) BEDGE(BNODE(I))=N1
				INC=INCOUNT(N1)
				
				WRITE(50,30), N1,BNODE(I)+N3,BNODE(I)+N4			
			ENDDO
		ENDDO
		
		!honrizontal triangle
		N1=0
		DO J=0,SOILLAYER
			N2=NEDGE1*J;
			DO I=1,NELT
				IF (ELT(I).ISDEL) CYCLE
				N1=N1+1			
				INC=INCOUNT(N1)
				NV1=2
				WRITE(50,40) N1,(EDGE(ELT(I).EDGE).NUM+N2)*ELT(I).ORIENT
				WRITE(50,50) N1,N1
			ENDDO		
		ENDDO
		!vertial boundary rectangle
		DO J=1,SOILLAYER
			N2=NEDGE1*(J-1);N3=NBNODE*(J-1)
			DO I=1,NCEDGE
				N1=N1+1			
				INC=INCOUNT(N1)
				NV1=3
                
				WRITE(50,40) N1,EDGE(CEDGE(I).EDGE).NUM+N2,BEDGE(EDGE(CEDGE(I).EDGE).V(2))+N3,&
								-(EDGE(CEDGE(I).EDGE).NUM+N2+NEDGE1),-(BEDGE(EDGE(CEDGE(I).EDGE).V(1))+N3)
				WRITE(50,50) N1,N1
				IF(J==1) EDGE(CEDGE(I).EDGE).BRECT=N1
			ENDDO		
		ENDDO		
		
		!VOLUME
		N1=0
		DO I=1,ZNUM
			DO J=0,SOILLAYER-1
				N1=N1+1
				INC=INCOUNT(N1)
				WRITE(CH1,'(I4)') I
				WRITE(CH2,'(I4)') J+1
				CH3='Z'//TRIM(ADJUSTL(CH1))//'_L'//TRIM(ADJUSTL(CH2))
				LEN1=LEN(TRIM(ADJUSTL(CH3)))	
			    !n4=ZONE(I).NTRIE3N
                STR1=""
				DO K=1,ZONE(I).NTRIE3N
					N3=ELT(ZONE(I).TRIE3N(K)).NUMBER
					DO K1=1,2
						IF(K1==1) THEN
							N2=ENUMBER*J+N3;
						ELSE
							N2=ENUMBER*(J+1)+N3;
						ENDIF						
						INC2=INCOUNT(N2)
						CH1=""
						WRITE(CH1,90) N2
						STR1=TRIM(ADJUSTL(STR1))//TRIM(ADJUSTL(CH1))	
					ENDDO
				ENDDO
				DO K=1,ZONE(I).NBE
					N2=NCEDGE*(J)+EDGE(ZONE(I).BEDGE(K)).BRECT
					INC2=INCOUNT(N2)
					CH1=""
					WRITE(CH1,90) N2
					STR1=TRIM(STR1)//TRIM(ADJUSTL(CH1))
				ENDDO			
				!WRITE(50,60) I,SIGN(FACEID1(ABS(PHYSICALGROUP(IPG1).TRISURFACE)),PHYSICALGROUP(IPG1).TRISURFACE), &
				!               SIGN(FACEID1(ABS(PHYSICALGROUP(IPG1).QUASURFACE)),PHYSICALGROUP(IPG1).QUASURFACE)
				INC2=LEN_TRIM(STR1)-1 !-1,get rid off ','
				IF(ZONE(I).NBE>0) THEN
					WRITE(50,61) N1,STR1(1:INC2)
					WRITE(50,70) N1,N1
					WRITE(50,80) TRIM(CH3),N1,N1
				ENDIF
			
			ENDDO
		ENDDO
		
		CLOSE(50)
		
		

	10  FORMAT("meshscale=1;")
	20	FORMAT("Point(",I<INC>,")={",3(E15.7,","),E15.7,"*meshscale};")
	30	FORMAT("Line(",I<INC>,")={",I6,",",I6,"};")
	40  FORMAT("Line Loop(",I<INC>,")={",<NV1>(I6,","),I6,"};")
	50	FORMAT("Plane Surface(",I<INC>,")={",I<INC>,"};")
	60	FORMAT("Surface Loop(",I<INC>,")={",<J>(I6,","),I6,"};")
	61	FORMAT("Surface Loop(",I<INC>,")={",A<INC2>"};")
	70	FORMAT("Volume(",I<INC>,")={",I<INC>,"};")
	80	FORMAT('Physical Volume("',A<LEN1>,'",',I<INC>,")={",I<INC>,"};")
	90 	FORMAT(I<INC2>,",")
	100 FORMAT('Physical Surface("',A<LEN1>,'",',I<INC>,")={",A<INC2>"};")
    
    CONTAINS
        !FUNCTION ORIENT_EDGE()
        FUNCTION EDGE_NODE(IEDGE,ILAYER) RESULT(V)
            INTEGER,INTENT(IN)::IEDGE,ILAYER !O<=ILAYER<=SOILLAYER
            INTEGER::V(2)           
            
            V=EDGE(IEDGE).V+NNODE*(ILAYER)
        
    
        ENDFUNCTION
        FUNCTION VEDGE_NODE(INODE,ISOIL) RESULT(V)
            INTEGER,INTENT(IN)::INODE,ISOIL !ISOIL>=1
            INTEGER::V(2)           
            
            V=[INODE+NNODE*(ISOIL-1),INODE+NNODE*(ISOIL)]
        
    
        ENDFUNCTION
	ENDSUBROUTINE	 
    
    INTEGER FUNCTION INCOUNT(N)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    REAL(8)::T1
    INTEGER::I

    T1=ABS(N)
    INCOUNT=INT(LOG10(T1))+1
	IF(N<0) INCOUNT=INCOUNT+1
    

	END FUNCTION
    
    
    SUBROUTINE CHECK_ORIENT()
        INTEGER::I,J,V1,V2
        
        DO I=1,NELT
            IF(ELT(I).ISDEL) CYCLE
            V2=EDGE(ELT(I).EDGE(1)).V(2)
            DO J=2,3
                V1=EDGE(ELT(I).EDGE(J)).V(1)
                IF(V2/=V1)THEN
                    ELT(I).ORIENT(J)=-1
                    V2=EDGE(ELT(I).EDGE(J)).V(1)
                ELSE
                    ELT(I).ORIENT(J)=1;V2=EDGE(ELT(I).EDGE(J)).V(2)
                ENDIF
            ENDDO
        ENDDO
        
    ENDSUBROUTINE
    
    
end module geomodel
