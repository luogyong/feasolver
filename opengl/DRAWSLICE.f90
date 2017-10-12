
SUBROUTINE GEN_SLICE_SURFACE_DATA()
    USE solverds
    USE MESHGEO
    USE view_modifier
    USE function_plotter
    IMPLICIT NONE

    real(8)::SLICEVAL1(NEDGE,NVO),SLICEPOINT1(3,NEDGE)
    INTEGER::NTRI1
    INTEGER,ALLOCATABLE::TRI1(:,:)
     INTEGER::I,AT1(NEDGE),J,K,N1,N2,ISLICE
    
    interface
    subroutine GENSLICE(AXIS,X,VAL,SLICEPOINT,TRI,NTRI)
        use MESHGEO
        use function_plotter
        implicit none
        integer,intent(in)::AXIS
        real(8),intent(in)::X
        integer,intent(out)::NTRI
        INTEGER,allocatable,intent(out)::TRI(:,:)
        real(8),intent(out)::SLICEPOINT(3,NEDGE),VAL(NEDGE,NVO)    
    end subroutine  
    
    endinterface    
    
    
     
    DO ISLICE=1,NSLICE
        CALL GENSLICE(SLICE(ISLICE).PLANE,SLICE(ISLICE).X,SLICEVAL1,SLICEPOINT1,TRI1,NTRI1)
        !PACKING
        AT1=0;N1=0
        DO I=1,NTRI1       
            DO J=1,3
                N2=TRI1(J,I)
                IF(AT1(N2)==0) THEN
                    N1=N1+1
                    AT1(N2)=N1
                ENDIF
                TRI1(J,I)=AT1(N2)
            ENDDO
        ENDDO
        SLICE(ISLICE).NV=N1;SLICE(ISLICE).NTRI=NTRI1
        IF(ALLOCATED(SLICE(ISLICE).V)) DEALLOCATE(SLICE(ISLICE).V,SLICE(ISLICE).VAL,SLICE(ISLICE).TRI)
        ALLOCATE(SLICE(ISLICE).V(3,SLICE(ISLICE).NV),SLICE(ISLICE).VAL(SLICE(ISLICE).NV,NVO))
        ALLOCATE(SLICE(ISLICE).TRI,SOURCE=TRI1)
        DO I=1,NEDGE
            IF(AT1(I)/=0) THEN
                SLICE(ISLICE).V(:,AT1(I))=SLICEPOINT1(:,I)
                SLICE(ISLICE).VAL(AT1(I),:)=SLICEVAL1(I,:)
            ENDIF
        ENDDO 
    

    ENDDO
    
    DEALLOCATE(TRI1)    
   
END SUBROUTINE    


SUBROUTINE GEN_SLICE_ISOLINE_DATA()
    USE solverds
    USE MESHGEO
    USE view_modifier
    USE function_plotter
    IMPLICIT NONE

    INTEGER::NISOLINE1,NISOLP1,NISOTRI1,NEDGE1,NFACE1
    INTEGER,ALLOCATABLE::ISOLINE1(:,:),ISOTRI1(:,:) 
    TYPE(TET_TYDEF),ALLOCATABLE::ELT1(:)
    TYPE(FACE_TYDEF),ALLOCATABLE::FACE1(:)
    TYPE(EDGE_TYDEF),allocatable::EDGE1(:)
    real(8),allocatable::PTISOLINE1(:,:)
    TYPE(hash_tbl_sll)::EDGE_TBL1,FACE_TBL1
    INTEGER::I,J,K,N1,N2,ISLICE
    INTEGER,ALLOCATABLE::AT1(:)
    interface
   
   SUBROUTINE SETUP_EDGE_TBLi(EDGE_TBL_L,TBL_SIZE_L,EDGE_L,NEDGE_L,ELEMENT_L,NEL_L,NODE_L,NNODE_L)
        USE MESHGEO
        USE hashtbl
        IMPLICIT NONE
	    INTEGER,INTENT(IN)::TBL_SIZE_L,NEL_L,NNODE_L
        INTEGER::NEDGE_L
	    TYPE(hash_tbl_sll)::EDGE_TBL_L
	    TYPE(EDGE_TYDEF),ALLOCATABLE::EDGE_L(:)
	    type(TET_tydef)::ELEMENT_L(NEL_L)
        REAL(8),INTENT(IN)::NODE_L(3,NNODE_L)
    ENDSUBROUTINE 
    
    SUBROUTINE SETUP_FACE_TBLi(FACE_TBL_L,TBL_SIZE_L,FACE_L,NFACE_L,EDGE_L,NEDGE_L,ELEMENT_L,NEL_L,NODE_L,NNODE_L)
        use MESHGEO
        use hashtbl
        IMPLICIT NONE
        INTEGER,INTENT(IN)::TBL_SIZE_L,NEL_L,NEDGE_L,NNODE_L
        INTEGER::NFACE_L
        REAL(8),INTENT(IN)::NODE_L(3,NNODE_L)
	    TYPE(hash_tbl_sll)::FACE_TBL_L
	    TYPE(EDGE_TYDEF),INTENT(IN)::EDGE_L(NEDGE_L)
        TYPE(FACE_TYDEF),ALLOCATABLE::FACE_L(:)
	    type(TET_TYDEF)::ELEMENT_L(NEL_L)
    ENDSUBROUTINE
    
    subroutine ContourLine(EDGE_L,NEDGE_L,FACE_L,NFACE_L,TET_L,NTET_L,XYZ,&
                            VA,nVA,VC,LINE,NLINE,CONTOURPOINT,TRI,NTRI)
		use MESHGEO
		implicit none
        integer,intent(in)::nVA,NEDGE_L,NFACE_L,NTET_L
        real(8),intent(in)::VA(nVA),VC,XYZ(3,NVA)
        TYPE(EDGE_TYDEF),INTENT(IN)::EDGE_L(NEDGE_L)
        TYPE(FACE_TYDEF),INTENT(IN)::FACE_L(NFACE_L)
        TYPE(TET_TYDEF),INTENT(IN)::TET_L(NTET_L)
		integer,intent(out)::NLINE,NTRI
		INTEGER,allocatable,intent(out)::LINE(:,:),TRI(:,:)
		real(8),intent(out)::CONTOURPOINT(3,NEDGE_L)
	end subroutine
    
    endinterface    
    
    
     
    DO ISLICE=1,NSLICE

    
        IF(ALLOCATED(ELT1)) DEALLOCATE(ELT1)
        ALLOCATE(ELT1(SLICE(ISLICE).NTRI))
        ELT1.NV=3;ELT1.GMET=2
        DO I=1,SLICE(ISLICE).NTRI
            ELT1(I).V(1:3)=SLICE(ISLICE).TRI(:,I)
        ENDDO
    
        CALL SETUP_EDGE_TBLi(EDGE_TBL1,2*SLICE(ISLICE).NTRI,EDGE1,NEDGE1,ELT1,SLICE(ISLICE).NTRI,SLICE(ISLICE).V,SLICE(ISLICE).NV) 
        CALL SETUP_FACE_TBLi(FACE_TBL1,SLICE(ISLICE).NTRI,FACE1,NFACE1,EDGE1,NEDGE1,ELT1,SLICE(ISLICE).NTRI,SLICE(ISLICE).V,SLICE(ISLICE).NV)
        CALL FACE_TBL1.FREE
        CALL EDGE_TBL1.FREE
        N1=0
        DO I=1,NEDGE1
            IF(EDGE1(I).ENUM==1) N1=N1+1        
        ENDDO
        IF(ALLOCATED(SLICE(ISLICE).BCEDGE)) DEALLOCATE(SLICE(ISLICE).BCEDGE)
        ALLOCATE(SLICE(ISLICE).BCEDGE(2,N1))
        SLICE(ISLICE).NBCE=N1
        N1=0
        DO I=1,NEDGE1
            IF(EDGE1(I).ENUM==1) THEN
                N1=N1+1
                SLICE(ISLICE).BCEDGE(:,N1)=EDGE1(I).V
            ENDIF
        ENDDO
        
		IF(CONTOURBAR.IVARPLOT/=OUTVAR(SLICE_PLOT_VARIABLE).IVO) CALL initialize_contourplot(OUTVAR(SLICE_PLOT_VARIABLE).IVO)        
        SLICE(ISLICE).NISOLINE=CONTOURBAR.NVAL
        IF(ALLOCATED(SLICE(ISLICE).ISOLINE)) DEALLOCATE(SLICE(ISLICE).ISOLINE)
        ALLOCATE(SLICE(ISLICE).ISOLINE(SLICE(ISLICE).NISOLINE))
		
        IF(ALLOCATED(PTISOLINE1)) DEALLOCATE(PTISOLINE1)
        ALLOCATE(PTISOLINE1(3,NEDGE1))
        DO I=1,CONTOURBAR.NVAL
            SLICE(ISLICE).ISOLINE(I).VAL=CONTOURBAR.VAL(I)
            CALL ContourLine(EDGE1,NEDGE1,FACE1,NFACE1,ELT1,SLICE(ISLICE).NTRI,SLICE(ISLICE).V,&
                          SLICE(ISLICE).VAL(:,CONTOURBAR.IVARPLOT),SLICE(ISLICE).NV,CONTOURBAR.VAL(I),ISOLINE1,NISOLINE1,PTISOLINE1,ISOTRI1,NISOTRI1)
            !!!ASSUME NEDGE1<=NEDGE...
            IF(ALLOCATED(AT1)) DEALLOCATE(AT1)
            ALLOCATE(AT1(NEDGE1))
            AT1=0;N1=0
            DO K=1,NISOLINE1       
                DO J=1,2
                    N2=ISOLINE1(J,K)
                    IF(AT1(N2)==0) THEN
                        N1=N1+1
                        AT1(N2)=N1
                    ENDIF
                    ISOLINE1(J,K)=AT1(N2)
                ENDDO
            ENDDO
            SLICE(ISLICE).ISOLINE(I).NV=N1;SLICE(ISLICE).ISOLINE(I).NE=NISOLINE1;
            IF(ALLOCATED(SLICE(ISLICE).ISOLINE(I).V)) DEALLOCATE(SLICE(ISLICE).ISOLINE(I).V,SLICE(ISLICE).ISOLINE(I).EDGE)        
            ALLOCATE(SLICE(ISLICE).ISOLINE(I).V(3,N1))
            ALLOCATE(SLICE(ISLICE).ISOLINE(I).EDGE,SOURCE=ISOLINE1)
            DO J=1,NEDGE1
                IF(AT1(J)/=0) THEN
                    SLICE(ISLICE).ISOLINE(I).V(:,AT1(J))=PTISOLINE1(:,J) 
                ENDIF
            ENDDO 
        ENDDO 
    ENDDO
    
    DEALLOCATE(ISOLINE1,ELT1,FACE1,EDGE1,PTISOLINE1) 
    IF(ALLOCATED(ISOTRI1)) DEALLOCATE(ISOTRI1)
    
     
END SUBROUTINE

SUBROUTINE INPUTSLICELOCATION()
    USE function_plotter 
    IMPLICIT NONE
    
    info.str='Set/Reset Slice Location(Format:x=1,2;y=0.3;z=1). Press Enter to complet and q to exit'C
    info.color=green;info.qkey=.TRUE.
    INFO.ISNEEDINPUT=.TRUE.
    INFO.INTERPRETER=STR2REALARRAY
    info.func_id=FUNC_ID_GETSLICELOCATION
    info.inputstr=''
    ISPLOTSLICE=.TRUE.
	ISPLOTSLICESURFACE=.TRUE.

ENDSUBROUTINE

SUBROUTINE GETSLICELOCATION()
    use function_plotter
    implicit none
    INTEGER::NX1=0,NY1=0,NZ1=0,I,J,N1
    REAL(8)::X1(10),Y1(10),Z1(10)
    

    
    NSLICE=0
    DO I=1,INFO.NINPUTVAR
        SELECT CASE(TRIM(ADJUSTL(INFO.INPUTVAR(I).NAME)))
        CASE('x','X')
            N1=1
        CASE('y','Y')
            N1=2
        CASE('z','Z')
            N1=3
        END SELECT
        
        DO J=1,INFO.INPUTVAR(I).NVAL
            NSLICE=NSLICE+1
            SLICE(NSLICE).PLANE=N1
            SLICE(NSLICE).X=INFO.INPUTVAR(I).VAL(J)
        ENDDO
        
    ENDDO
    !do while(INFO.ISNEEDINPUT)
    !    
    !    
    !enddo
    
    !WRITE(*,10)
    !WRITE(*,20), 'X'
    !READ(*,*) NX1,X1(1:NX1)
    !WRITE(*,20), 'Y'
    !READ(*,*) NY1,Y1(1:NY1)
    !WRITE(*,20), 'Z'
    !READ(*,*) NZ1,Z1(1:NZ1)
    
    

    !IF(ALLOCATED(SLICE)) DEALLOCATE(SLICE)
    !NSLICE=NX1+NY1+NZ1
    IF(NSLICE>30) STOP 'NSLICE MUST BE <=30.'
    
    IF(NSLICE<1) RETURN
    
    !DO I=1,NX1
    !    SLICE(I).PLANE=1;SLICE(I).X=X1(I)
    !ENDDO
    !DO I=1,NY1
    !    SLICE(NX1+I).PLANE=2;SLICE(NX1+I).X=Y1(I)
    !ENDDO    
    !DO I=1,NZ1
    !    SLICE(NX1+NY1+I).PLANE=3;SLICE(NX1+NY1+I).X=Z1(I)
    !ENDDO
	
	CALL GEN_SLICE_SURFACE_DATA()
    CALL GEN_SLICE_ISOLINE_DATA()
    CALL SLICEPLOT()
	
10 FORMAT('INPUT FORMAT:NSLICES,LOC1,LOC2,....(NSLICES<=10)')
20 FORMAT('SET SLICES PARALLEL TO ', A, ' PLANE:')
!30 FORMAT(I2,<N1>F14.6)

ENDSUBROUTINE



SUBROUTINE SLICEPLOT()
    use solverds
    use opengl_gl
    use function_plotter
    use MESHGEO
    use view_modifier
    implicit none

    real(GLFLOAT),allocatable::vcolor(:,:)
    integer :: i,j,k,n1,MAT1(10),ISLICE  



    n1=outvar(SLICE_PLOT_VARIABLE).ivo
    IF(CONTOURBAR.IVARPLOT/=n1) call initialize_contourplot(n1)

    call glDeleteLists(SLICELIST, 1_glsizei)
    call reset_view
    
    call glNewList(SLICELIST, gl_compile_and_execute)

    


    call glPolygonMode(gl_front_and_back, gl_fill)
    call glenable(gl_polygon_offset_fill)
    call glPolygonoffset(1._glfloat,1._glfloat)
	
	
    
    call glBegin(gl_triangles)
	
	DO ISLICE=1,NSLICE
        if (ISPLOTSLICESURFACE) then
		    if(allocated(vcolor)) deallocate(vcolor)
		    allocate(vcolor(4,SLICE(ISLICE).NV))
 
		    do i=1,SLICE(ISLICE).NV
			    call get_rainbow(SLICE(ISLICE).VAL(I,N1),contourbar.val(1),contourbar.val(contourbar.nval),vcolor(:,i))
			    if(isTransparency) vcolor(4,i)=0.6
		    enddo
        ENDIF
        
		do i=1,SLICE(ISLICE).NTRI

			do j=1,3            
				IF(ISPLOTSLICESURFACE)THEN
					call glcolor4fv(vcolor(:,SLICE(ISLICE).TRI(j,I)))
				ELSE
					call glcolor4fv(mycolor(:,gray))
				ENDIF
				call glvertex3dv(SLICE(ISLICE).V(:,SLICE(ISLICE).TRI(J,I)))            
			enddo    
		enddo
   ENDDO
   
   call glEnd
   
   CALL glBegin(gl_LINES)
   
   call glcolor4fv(mycolor(:,BLACK))
   DO I=1,SLICE(ISLICE).NBCE
        DO J=1,2
		    call glvertex3dv(SLICE(ISLICE).V(:,SLICE(ISLICE).BCEDGE(J,I))) 
		ENDDO
   ENDDO
   
   CALL GLEND
   
   


    IF(ISPLOTSLICEISOLINE) THEN

	
	
	    call glBegin(gl_LINES)
	
	    CALL glcolor4fv(mycolor(:,BLACK))
	    DO ISLICE=1,NSLICE
	
		    DO I=1,SLICE(ISLICE).NISOLINE
			    DO J=1,SLICE(ISLICE).ISOLINE(I).NE
				    DO K=1,2
					    call glvertex3dv(SLICE(ISLICE).ISOLINE(I).V(:,SLICE(ISLICE).ISOLINE(I).EDGE(K,J))) 
				    ENDDO
			    ENDDO
		    ENDDO
		
	    ENDDO
	
	    CALL GLEND()

    ENDIF
 


    call glEndList

    call gldisable(gl_polygon_offset_fill)


    call glutPostRedisplay

if(allocated(vcolor)) deallocate(vcolor)
!if(allocated(contour_value)) deallocate(contour_value)
	

ENDSUBROUTINE
    
subroutine GENSLICE(AXIS,X,VAL,CONTOURPOINT,TRI,NTRI)
    use solverds
    use MESHGEO
    use function_plotter
    implicit none
    integer,intent(in)::AXIS
    real(8),intent(in)::X
    integer,intent(out)::NTRI
    INTEGER,allocatable,intent(out)::TRI(:,:)
    real(8),intent(out)::CONTOURPOINT(3,NEDGE),VAL(NEDGE,NVO)
 
    real(8)::V1,V2,T1,T2,T3,T4,NTE1,PT1(3,6)
    INTEGER::ISPVC1(NEDGE),E1(4)=0,E2(6)=0,E3(6)=0,IPT1(6),IPT2(6),IPT3(6)
    integer::i,j,MAXNTE1=1000,N1,MAXNTRI1=1000,N2,N3,K,VAR1,NFACE1,NCONTRI1=0
    REAL(GLDOUBLE)::VEC1(3,NNUM),SCALE1,X1(3),X2(3) 
	REAL(8),EXTERNAL::VECTOR_SCALE
    
    
    interface
    
    SUBROUTINE I2_ENLARGE_AR(AVAL,DSTEP,DIM1)
        INTEGER,ALLOCATABLE,INTENT(INOUT)::AVAL(:,:)
        INTEGER,INTENT(IN)::DSTEP,DIM1
    end subroutine
    
   
    endinterface
	
	
    
    IF(ALLOCATED(TRI)) DEALLOCATE(TRI)
	ALLOCATE(TRI(3,MAXNTRI1))	
	
	VEC1=0;SCALE1=1.D0
	IF(IsContour_In_DeformedMesh.AND.OUTVAR(DISX).VALUE>0.AND.OUTVAR(DISY).VALUE>0) THEN
		VEC1(1,:)=NODALQ(:,OUTVAR(DISX).IVO,STEPPLOT.ISTEP)
		VEC1(2,:)=NODALQ(:,OUTVAR(DISY).IVO,STEPPLOT.ISTEP)
		IF(NDIMENSION>2) VEC1(3,:)=NODALQ(:,OUTVAR(DISZ).IVO,STEPPLOT.ISTEP)
		scale1=Scale_Deformed_Grid*STEPPLOT.VSCALE(1)
	ENDIF	
	
	!VAR1=outvar(CONTOUR_PLOT_VARIABLE).ivo
    ISPVC1=0
    do i=1,nedge
		IF(EDGE(I).ISDEAD==1) CYCLE
        X1=NODE(EDGE(I).V(1)).COORD+VEC1(:,EDGE(I).V(1))*SCALE1
		X2=NODE(EDGE(I).V(2)).COORD+VEC1(:,EDGE(I).V(2))*SCALE1
        
        T1=X-X1(AXIS);T2=X-X2(AXIS);T3=X2(AXIS)-X1(AXIS)
        IF(ABS(T3)<1.D-7)THEN
            !IF(ABS(T1)<1.D-14) THEN                
            !    ISPVC1(I)=2 !两个节点都是
            !    PVC1(:,I)=NODE(V2).COORD
            !ENDIF        
        ELSE
            IF(T1*T2<=0) THEN
                ISPVC1(I)=1
                T4=MIN(MAX(0.D0,T1/T3),1.D0)
                CONTOURPOINT(:,I)=X1+t4*(X2-X1)
                DO J=1,NVO
					V1=NODALQ(EDGE(I).V(1),J,STEPPLOT.ISTEP);V2=NODALQ(EDGE(I).V(2),J,STEPPLOT.ISTEP);
					VAL(I,J)=V1+t4*(V2-V1)
                ENDDO
            ENDIF   
        ENDIF 
    enddo

    
	NTRI=0
	DO I=1,NTET
		if(tet(i).isdead==1) cycle
        IF(TET(I).DIM/=3) CYCLE
        
		E2=TET(I).E
		IF(SUM(ISPVC1(E2))<3) CYCLE
        N1=0
        
        IPT1=0
        DO J=1,TET(J).NE
            IF(ISPVC1(E2(J))==1) THEN
                N1=N1+1
                PT1(:,N1)=CONTOURPOINT(:,E2(J))
                IPT1(N1)=J
            ENDIF
        ENDDO
        IPT2=0
        CALL removeduplicatedpoint(Pt1(:,1:N1),IPT2(1:N1),n1)
        IF(SUM(IPT2)<3) CYCLE
        
        IPT3=0
        DO J=1,N1
            IF(IPT2(J)==1) IPT3(IPT1(J))=1
        ENDDO
        
        N1=0
        NTRI=NTRI+1
		IF(NTRI+1>MAXNTRI1) THEN !至少有两个空间
            CALL I2_ENLARGE_AR(TRI,500,3)
            MAXNTRI1=MAXNTRI1+500
        ENDIF
        E3=0
        DO J=1,6      
            IF(IPT3(J)==1) THEN              
                N1=N1+1
                IF(N1<4) THEN
					TRI(N1,NTRI)=E2(J)
					E3(N1)=J                    
				ELSEIF(N1==4) THEN
					NTRI=NTRI+1
					TRI(1,NTRI)=E2(J)
					IF(J==5) N2=3
					IF(J==4) N2=2
					IF(J==6) N2=1
					N3=1
					DO K=1,3
						IF(E3(K)/=N2) THEN
							N3=N3+1
							TRI(N3,NTRI)=E2(E3(K))							
						ENDIF
					ENDDO
					
				ENDIF
            ENDIF
       ENDDO 		
	ENDDO
    

    
	
end subroutine


subroutine removeduplicatedpoint(Pt,unify,ns)
    implicit none
    integer,intent(in)::ns
    integer,intent(out)::unify(ns)
    real(8),intent(in)::Pt(3,ns)
    integer::i,j
    
    
    unify=1
    do i=1,ns-1
        if(unify(i)==0) cycle
        do j=i+1,ns
            if(unify(j)==0) cycle
            if(norm2(pt(:,i)-pt(:,j))<1.e-10) unify(j)=0
        enddo
    enddo
    
endsubroutine

    
