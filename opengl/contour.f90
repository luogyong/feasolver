subroutine DrawSurfaceContour()
!use solverds
use opengl_gl
use function_plotter
!use MESHGEO
!use view_modifier
implicit none
REAL(GLDOUBLE)::SCALE1
REAL(8),ALLOCATABLE::VEC1(:,:)
REAL(8),EXTERNAL::VECTOR_SCALE
real(GLFLOAT),allocatable::vcolor(:,:)


integer :: i,j,k,n1,MAT1(10),ie1  !,nfrac

! prepare to make a new display list

ALLOCATE(VEC1(3,POSDATA.NNODE))
VEC1=0.D0
if (draw_surface_solid) then

    SCALE1=1.D0
    IF(IsContour_In_DeformedMesh.AND.POSDATA.IDISX>0.AND.POSDATA.IDISY>0) THEN

		VEC1(1,:)=POSDATA.NODALQ(:,POSDATA.IDISX,STEPPLOT.ISTEP)
		VEC1(2,:)=POSDATA.NODALQ(:,POSDATA.IDISY,STEPPLOT.ISTEP)
		IF(POSDATA.NDIM>2) VEC1(3,:)=POSDATA.NODALQ(:,POSDATA.IDISZ,STEPPLOT.ISTEP)

        scale1=STEPPLOT.VSCALE(1)*Scale_Deformed_Grid
    ENDIF

    n1=CONTOUR_PLOT_VARIABLE

    call glDeleteLists(contourlist, 1_glsizei)


    call reset_view
    call glNewList(contourlist, gl_compile_and_execute)

! draw the solid surface

    call glPolygonMode(gl_front_and_back, gl_fill)
    call glenable(gl_polygon_offset_fill)
    call glPolygonoffset(0.5_glfloat,0.5_glfloat)
    call glBegin(gl_triangles)

   
    if (surface_color == rainbow_surface) then
        if(allocated(vcolor)) deallocate(vcolor)
        allocate(vcolor(4,POSDATA.NNODE))

        
        do i=1,POSDATA.NNODE
            call get_rainbow(POSDATA.NODALQ(i,n1,STEPPLOT.ISTEP),contourbar.val(1),contourbar.val(contourbar.nval),vcolor(:,i))
            if(isTransparency) vcolor(4,i)=0.6
        enddo
 
    
    endif

    do i=1,nface
        !if(face(i).shape/=3) cycle
		!只对外边界及材料边界进行渲染。
		IF(FACE(I).ISDEAD==1) CYCLE
		MAT1(1:FACE(I).ENUM)=POSDATA.ELEMENT(TET(ABS(FACE(I).ELEMENT)).MOTHER).ISET !
        IF(FACE(I).ENUM>1) THEN			
			IF(all(MAT1(1:FACE(I).ENUM)-MAT1(1)==0)) CYCLE		
		ENDIF
        do j=1,3
            !call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,vcolor(:,face(i).v(j)))
            IF(surface_color == white_surface)THEN
                call glcolor4fv(mycolor(:,light_gray))
            ELSEIF(surface_color==color_by_set) THEN
                call glcolor4fv(DISTINCT_COLOR(:,MAX(mod(mat1(1),45),1)))
            ELSE
                call glcolor4fv(vcolor(:,face(i).v(j)))
            ENDIF
            call glvertex3dv(POSDATA.node(face(i).v(j)).coord+VEC1(:,face(i).v(j))*SCALE1)            
        enddo    
    enddo


   call glEnd

   !call glBegin(GL_QUADS)
   !
   !
   ! !do i=1,enum
   ! !    IF(ESET(ELEMENT(I).SET).COUPLESET<ELEMENT(I).SET) CYCLE
   ! !    do j=1,3            
   ! !        call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,vcolor(:,element(i).node(j)))
   ! !        call glvertex3dv(node(element(i).node(j)).coord)
   ! !    enddo
   ! !enddo
   !
   ! do i=1,nface
   !     if(face(i).shape/=4) cycle
   !     do j=1,4
   !         !call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,vcolor(:,face(i).v(j)))
   !         IF(surface_color == white_surface)THEN
   !             call glcolor4fv(GRAY)
   !         ELSE
   !             call glcolor4fv(vcolor(:,face(i).v(j)))
   !         ENDIF
   !         call glvertex3dv(node(face(i).v(j)).coord+VEC1(:,face(i).v(j))*SCALE1)            
   !     enddo    
   ! enddo
   !
   !
   !call glEnd   
   
   
   !if (surface_color == rainbow_surface) call Color_Bar(contour_value,num_contour,NFRAC,color_bar_title)


    call glEndList

    call gldisable(gl_polygon_offset_fill)
    
    
    call glutPostRedisplay

    
    
endif ! draw_surface_solid

if(allocated(vcolor)) deallocate(vcolor)

DEALLOCATE(VEC1)
!if(allocated(contour_value)) deallocate(contour_value)

endsubroutine  


subroutine DrawLineContour()
!use solverds
use opengl_gl
use function_plotter
!use MESHGEO
!use view_modifier
implicit none
real(GLDOUBLE),ALLOCATABLE::CPNT1(:,:)    
INTEGER,allocatable::CLINE1(:,:),TRI1(:,:)
INTEGER::NCL1=0,NTRI1=0,AT1(NEDGE)
REAL(GLDOUBLE),ALLOCATABLE::VEC1(:,:),XYZ(:,:) 
REAL(GLDOUBLE)::SCALE1


integer :: i,j,k,n1,N2  !,nfrac

! prepare to make a new display list

interface

	subroutine ContourLine(EDGE_L,NEDGE_L,FACE_L,NFACE_L,TET_L,NTET_L,XYZ,&
                            VA,nVA,VC,LINE,NLINE,CONTOURPOINT,TRI,NTRI)
		use MESHGEO
		implicit none
        integer,intent(in)::nVA,NEDGE_L,NFACE_L,NTET_L
        real(8),intent(in)::VA(nVA),VC,XYZ(3,nVA)
        TYPE(EDGE_TYDEF),INTENT(IN)::EDGE_L(NEDGE_L)
        TYPE(FACE_TYDEF),INTENT(IN)::FACE_L(NFACE_L)
        TYPE(TET_TYDEF),INTENT(IN)::TET_L(NTET_L)
		integer,intent(out)::NLINE,NTRI
		INTEGER,allocatable,intent(out)::LINE(:,:),TRI(:,:)
		real(8),intent(out)::CONTOURPOINT(3,NEDGE_L)
	end subroutine

endinterface


ALLOCATE(CPNT1(3,NEDGE),VEC1(3,POSDATA.NNODE),XYZ(3,POSDATA.NNODE))    

if (draw_contour) then

    SCALE1=1.D0;VEC1=0.D0
    IF(IsContour_In_DeformedMesh.AND.POSDATA.IDISX>0.AND.POSDATA.IDISY>0) THEN
	    VEC1(1,:)=POSDATA.NODALQ(:,POSDATA.IDISX,STEPPLOT.ISTEP)
	    VEC1(2,:)=POSDATA.NODALQ(:,POSDATA.IDISY,STEPPLOT.ISTEP)
	    IF(POSDATA.NDIM>2) VEC1(3,:)=POSDATA.NODALQ(:,POSDATA.IDISZ,STEPPLOT.ISTEP)
	    scale1=Scale_Deformed_Grid*STEPPLOT.VSCALE(1)
    ENDIF
    DO I=1,POSDATA.NNODE
        XYZ(:,I)=POSDATA.NODE(I).COORD+VEC1(:,I)*SCALE1    
    ENDDO
    IF(ALLOCATED(ISOLINE)) DEALLOCATE(ISOLINE)
    NISOLINE=CONTOURBAR.NVAL
    ALLOCATE(ISOLINE(NISOLINE))

    do i=1,contourbar.nval
        CPNT1=0.0d0
        
	    call ContourLine(EDGE,NEDGE,FACE,NFACE,TET,NTET,XYZ,POSDATA.NODALQ(:,CONTOUR_PLOT_VARIABLE,STEPPLOT.ISTEP),POSDATA.NNODE, &
                            contourbar.VAL(I),CLINE1,NCL1,CPNT1,TRI1,NTRI1)
        AT1=0;N1=0
        DO J=1,NCL1
            DO K=1,2
                N2=CLINE1(K,J)
                IF(AT1(N2)==0) THEN
                    N1=N1+1
                    AT1(N2)=N1
                ENDIF
                CLINE1(K,J)=AT1(N2)
            ENDDO
        ENDDO
        DO J=1,NTRI1
            DO K=1,3
                N2=TRI1(K,J)
                IF(AT1(N2)==0) THEN
                    N1=N1+1
                    AT1(N2)=N1
                ENDIF
                TRI1(K,J)=AT1(N2)
            ENDDO
        ENDDO
        ISOLINE(I).NV=N1;ISOLINE(I).NE=NCL1;ISOLINE(I).NTRI=NTRI1
        ISOLINE(I).VAL=contourbar.val(I)
        ALLOCATE(ISOLINE(I).V(3,N1))
        DO J=1,NEDGE
            IF(AT1(J)/=0) ISOLINE(I).V(:,AT1(J))=CPNT1(:,J)
        ENDDO
        IF(NCL1>0) ALLOCATE(ISOLINE(I).EDGE,SOURCE=CLINE1)
        IF(NTRI1>0) ALLOCATE(ISOLINE(I).TRI,SOURCE=TRI1)    
    ENDDO


 


    call reset_view 
    call glDeleteLists(contourLinelist, 1_glsizei)
    call glNewList(contourLinelist, gl_compile_and_execute)
    
    call gldisable(GL_CULL_FACE);    
    call glPolygonMode(gl_front_and_back, gl_fill)  
  
    
	do i=1,NISOLINE
  !      CPNT1=0.0d0
		!call ContourLine(EDGE,NEDGE,FACE,NFACE,TET,NTET,XYZ,NODALQ(:,N1),POSDATA.NNODE,contourbar.VAL(I),CLINE1,NCL1,CPNT1,TRI1,NTRI1)
		call glcolor4fv(mycolor(:,black))		
		CALL GLBEGIN(GL_LINES)
			DO J=1,ISOLINE(I).NE
                DO K=1,2
				call glvertex3dv(ISOLINE(I).V(:,ISOLINE(I).EDGE(K,J))) 
				ENDDO  
			ENDDO
		CALL GLEND()
        
        call glcolor4fv(contourbar.Color(:,i))       
        
        CALL GLBEGIN(GL_TRIANGLES)
			DO J=1,ISOLINE(I).NTRI
                DO K=1,3
				    call glvertex3dv(ISOLINE(I).V(:,ISOLINE(I).TRI(K,J)))
                ENDDO
			ENDDO		
        
        CALL GLEND()
	enddo
    
    call glEndList
    
    call glutPostRedisplay
endif

DEALLOCATE(CPNT1,VEC1,XYZ) 

endsubroutine     

subroutine ContourLine(EDGE_L,NEDGE_L,FACE_L,NFACE_L,TET_L,NTET_L,XYZ,&
                    VA,nVA,VC,LINE,NLINE,CONTOURPOINT,TRI,NTRI)

    !use MESHGEO
    use function_plotter
    implicit none
    integer,intent(in)::nVA,NEDGE_L,NFACE_L,NTET_L
    real(8),intent(in)::VA(nVA),VC,XYZ(3,nVA)
    TYPE(EDGE_TYDEF),INTENT(IN)::EDGE_L(NEDGE_L)
    TYPE(FACE_TYDEF),INTENT(IN)::FACE_L(NFACE_L)
    TYPE(TET_TYDEF),INTENT(IN)::TET_L(NTET_L)
    integer,intent(out)::NLINE,NTRI
    INTEGER,allocatable,intent(out)::LINE(:,:),TRI(:,:)
    real(8),intent(out)::CONTOURPOINT(3,NEDGE_L)
    real(8)::V1,V2,T1,T2,T3,T4,X1(3),X2(3),PT1(3,6) 
    INTEGER::ISPVC1(NEDGE_L),E1(4)=0,E2(6)=0,E3(6)=0,IPT1(6),IPT2(6),IPT3(6)
    integer::i,j,MAXNLINE1=1000,N1,MAXNTRI1=1000,N2,N3,K
    
	
    
    interface
    SUBROUTINE I2_ENLARGE_AR(AVAL,DSTEP,DIM1)
        INTEGER,ALLOCATABLE,INTENT(INOUT)::AVAL(:,:)
        INTEGER,INTENT(IN)::DSTEP,DIM1
    end subroutine
    endinterface
	
	
    
    IF(ALLOCATED(LINE)) DEALLOCATE(LINE)
	ALLOCATE(LINE(2,MAXNLINE1))
    IF(ALLOCATED(TRI)) DEALLOCATE(TRI)
	ALLOCATE(TRI(3,MAXNTRI1))	
	
	
	
	
    ISPVC1=0
    do i=1,NEDGE_L
		IF(EDGE_L(I).ISDEAD==1) CYCLE
        V1=VA(EDGE_L(I).V(1));V2=VA(EDGE_L(I).V(2))
        T1=VC-V1;T2=VC-V2;T3=V2-V1
        IF(ABS(T3)<1.D-7)THEN
            !IF(ABS(T1)<1.D-14) THEN                
            !    ISPVC1(I)=2 !两个节点都是
            !    PVC1(:,I)=NODE(V2).COORD
            !ENDIF        
        ELSE
            IF(T1*T2<=0) THEN
                ISPVC1(I)=1
                T4=MIN(MAX(0.D0,T1/T3),1.D0)
				X1=XYZ(:,EDGE_L(I).V(1))
				X2=XYZ(:,EDGE_L(I).V(2))
                CONTOURPOINT(:,I)=X1+t4*(X2-X1)
            ENDIF   
        ENDIF 
    enddo
    NLINE=0	
    DO I=1,NFACE_L
		!IF(ANY(TET(ABS(FACE_L(I).ELEMENT(1:FACE_L(I).ENUM))).MOTHER>0)) CYCLE !NOT PLANE ELEMENT
		IF(FACE_L(I).ISDEAD==1) CYCLE
		IF(FACE_L(I).ENUM>1) CYCLE
		E1=ABS(FACE_L(I).EDGE)		
        !!!!!!!!!!!!!!!!!!!!
		IF(SUM(ISPVC1(E1(1:FACE_L(I).SHAPE)))/=2) CYCLE
        N1=0
        NLINE=NLINE+1
		IF(NLINE>MAXNLINE1) THEN
            CALL I2_ENLARGE_AR(LINE,500,2)
            MAXNLINE1=MAXNLINE1+500
        ENDIF        
        DO J=1,FACE_L(I).SHAPE        
            IF(ISPVC1(E1(J))==1) THEN              
                N1=N1+1
                LINE(N1,NLINE)=E1(J)       
            ENDIF
       ENDDO    
    ENDDO
    
	NTRI=0
	DO I=1,NTET_L
		IF(TET_L(I).ISDEAD==1) CYCLE
        IF(TET_L(I).DIM/=3) CYCLE
        
		E2=TET_L(I).E
		IF(SUM(ISPVC1(E2))<3) CYCLE
        N1=0        
        IPT1=0
        DO J=1,TET_L(J).NE
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
        
        NTRI=NTRI+1
		IF(NTRI+1>MAXNTRI1) THEN !至少有两个空间
            CALL I2_ENLARGE_AR(TRI,500,3)
            MAXNTRI1=MAXNTRI1+500
        ENDIF
        
        
        E3=0;N1=0
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

subroutine initialize_contourplot(IVARPLOT)
    !use solverds
    use function_plotter
    implicit none
	INTEGER,INTENT(IN)::IVARPLOT
    integer::i
    real(GLDOUBLE) :: graphmin,graphmax,graphstep
    !real(GLDOUBLE),allocatable :: contour_value(:)
    character(128)::str1,str2  
    
    CONTOURBAR.IVARPLOT=IVARPLOT
    minv=minval(POSDATA.NODALQ(:,IVARPLOT,:))
    maxv=maxval(POSDATA.NODALQ(:,IVARPLOT,:))
    WRITE(STR1,'(G13.3)') MAXV;WRITE(STR2,'(G13.3)') MINV;
        
    contourbar.title=trim(adjustl(POSDATA.OUTVAR(IVARPLOT).name))//',MAX='//TRIM(ADJUSTL(STR1))//',MIN='//TRIM(ADJUSTL(STR2))
    IF(ABS(MINV)<1E-7) MINV=SIGN(1E-7,MINV)
    call loose_label(minv, maxv,init_num_contour,graphmin,graphmax,graphstep,contourbar.nfrac)
    graphstep=graphstep/scale_contour_num
    contourbar.nval=MAX(floor((graphmax-graphmin)/graphstep)+1,2)
    if(allocated(contourbar.VAL)) deallocate(contourbar.VAL)
    if(allocated(contourbar.color)) deallocate(contourbar.color)
    allocate(contourbar.VAL(contourbar.nval),contourbar.color(4,contourbar.nval))
    do i=1,contourbar.nval
        contourbar.VAL(i)=graphmin+(i-1)*graphstep
        call get_rainbow(contourbar.VAL(i),graphmin,graphmax,contourbar.color(:,i))            
    enddo
end subroutine

    
!---+----3----+----2----+----1----+---<>---+----1----+----2----+----3----+----4
!--------------------------------  Color_Bar  ---------------------------------
 
subroutine Color_Bar()
    use opengl_gl
    use opengl_glut
	use function_plotter
    implicit none
    integer i;
    !real(GLFLOAT):: color(4,num_contour)
    real(GLDOUBLE)::left1(2,contourbar.nval),right1(2,contourbar.nval),SPACING,orig(3),PPM1,scale1
    CHARACTER(16)::STR1    
    integer(glCint),dimension(4)::viewport1

           
    !! Set up coordinate system to position color bar near bottom of window.
    call glgetintegerv(gl_viewport,viewport1)
    CALL glMatrixMode(GL_PROJECTION);
    CALL glPushMatrix()
    CALL glLoadIdentity();
    CALL glOrtho(0., real(viewport1(3),gldouble), 0., real(viewport1(4),gldouble), -1._GLDOUBLE, 1._GLDOUBLE);
    CALL glMatrixMode(GL_MODELVIEW);
    CALL glPushMatrix()
    CALL glLoadIdentity();
    
    call glPushAttrib(GL_ALL_ATTRIB_BITS)
	!call glPushAttrib(GL_LIGHTING_BIT .or. GL_CURRENT_BIT); ! lighting and color mask
	call glDisable(GL_LIGHTING);     ! need to disable lighting for proper text color
	call glDisable(GL_TEXTURE_2D);
    call glDisable(GL_CULL_FACE);
	call glDepthFunc(GL_ALWAYS);
    
    
    
    orig(1)=7./8.*viewport1(3);orig(2)=1./5.*viewport1(4);orig(3)=0.
    SPACING=3./5.*viewport1(4)/(contourbar.nval-1)
    
    left1(1,:)=orig(1);right1(1,:)=orig(1)+1./60.*viewport1(3);
    DO I=1,contourbar.nval
        left1(2,I)=orig(2)+(I-1)*SPACING
        right1(2,I)=left1(2,I)
        !CALL get_rainbow(contour_value(I),contour_value(1),contour_value(contourbar.nval),color(:,I))
    ENDDO    
    
    ! Use Quad strips to make color bar.

    call glPolygonMode(gl_front_and_back, gl_fill)
    CALL glBegin(GL_QUAD_STRIP)
        do i = 1,contourbar.nval
            CALL glColor3fv(contourbar.color(:,i));
            CALL glVertex2Dv(left1(:,i));
            CALL glVertex2Dv(right1(:,i));
        enddo
    CALL glEnd()
 
    ! Label ends of color bar.
 
    CALL glColor3f(0._GLFLOAT, 0._GLFLOAT, 0._GLFLOAT)
    
    call glPolygonMode(gl_front_and_back, gl_line)
    
    CALL glBegin(GL_QUAD_STRIP)
    do i = 1,contourbar.nval
        !CALL glColor3fv(color(:,i));
        CALL glVertex2Dv(left1(:,i));
        CALL glVertex2Dv(right1(:,i));
    enddo
    CALL glEnd()
    
    
 
    !CALL bitmap_output(-5, 7, 0, "Min_H", GLUT_BITMAP_9_BY_15);
    !CALL bitmap_output(95, 7, 0, "Max_H", GLUT_BITMAP_9_BY_15);
    DO I=1,contourbar.nval
        WRITE(STR1,10) contourbar.val(I)
        call output(REAL(right1(1,1)+viewport1(3)/200.,GLFLOAT),REAL(right1(2,I),GLFLOAT),TRIM(ADJUSTL(STR1)))
    ENDDO
    
    !CALL GLROTATED(90.,0.,0.,1.)
    !104.76/viewport1(3)
    !call drawStrokeText(ORIG(1), ORIG(2),ORIG(3),TRIM(ADJUSTL(BAR_TITLE)))
    !call output(0._GLFLOAT, 0._glfloat,TRIM(ADJUSTL(BAR_TITLE)))
    !assume window size equals to viewport, the unit is pixel
    PPM1=glutget(GLUT_SCREEN_WIDTH)/REAL(glutget(GLUT_SCREEN_WIDTH_MM)) !PIXELS PER MM
    !10 pound 
    scale1=PPM1*3.527777778/119.05*0.75
    call glLineWidth(1.0_glfloat)
    call drawStrokeText(90.0,ORIG(1)-viewport1(3)/250., ORIG(2),ORIG(3),scale1,TRIM(ADJUSTL(contourbar.title)))
    
    CALL glPopMatrix()
    CALL glMatrixMode(GL_PROJECTION)
    CALL glPopMatrix()
    
 !  	CALL glEnable(GL_TEXTURE_2D);
	!CALL glEnable(GL_LIGHTING);
 !   call glEnable(GL_CULL_FACE);
	!CALL glDepthFunc(GL_LEQUAL);
    call glPolygonMode(gl_front_and_back, gl_FILL)
	CALL glPopAttrib();
    
    call glutPostRedisplay
    
10  FORMAT(F16.<contourbar.nfrac>)    
end subroutine

subroutine get_rainbow(val,minval,maxval,c)
use opengl_gl
implicit none
real(GLDOUBLE), intent(in) :: val,maxval,minval
real(GLFLOAT), intent(out) :: c(4)

real(GLFLOAT) :: f

if (maxval > minval) then
   f = (val-minval)/(maxval-minval)
else ! probably maxval==minval
   f = 0.5_glfloat
endif

if (f < .25) then
   c(1) = 0.0_glfloat
   c(2) = 4.0_glfloat * f
   c(3) = 1.0_glfloat
   c(4) = 1.0_glfloat
elseif (f < .5) then
   c(1) = 0.0_glfloat
   c(2) = 1.0_glfloat
   c(3) = 2.0_glfloat - 4.0_glfloat*f
   c(4) = 1.0_glfloat
elseif (f < .75) then
   c(1) = 4.0_glfloat * f - 2.0_glfloat
   c(2) = 1.0_glfloat
   c(3) = 0.0_glfloat
   c(4) = 1.0_glfloat
else
   c(1) = 1.0_glfloat
   c(2) = 4.0_glfloat - 4.0_glfloat*f
   c(3) = 0.0_glfloat
   c(4) = 1.0_glfloat
endif

end subroutine get_rainbow


subroutine loose_label(minv, maxv,ntick,graphmin,graphmax,graphstep,nfrac)
    implicit none
    integer,intent(in)::ntick
    real(8),intent(in)::minv,maxv
    integer,intent(out)::nfrac
    real(8),intent(out)::graphmin,graphmax,graphstep !/* graph range min and max and spacing */

    character::str(6), temp(20);
    !real(8) graphstep;				!/* tick mark spacing */
    !real(8) graphmin, graphmax;		
    real(8) range, x;
    
    real(8):: nicenum
    external nicenum
    
    !/* we expect min!=max */
    range = nicenum(maxv-minv, 0);
    graphstep = nicenum(range/(NTICK-1), 1);
    graphmin = floor(minv/graphstep)*graphstep;
    graphmax = ceilING(maxv/graphstep)*graphstep;
    !NFRAC=MAX(1,2)
    nfrac = MAX(-floor(log10(graphstep)),0)	!/* # of fractional digits to show */
    
    !sprintf(str, "%%.%df", nfrac);	!/* simplest axis labels */
    !WRITE(STR,"F3.1") NFRAC
    !printf("graphmin=%g graphmax=%g increment=%g\n", graphmin, graphmax, graphstep);

end subroutine

!/*
! * nicenum: find a "nice" number approximately equal to x.
! * Round the number if round=1, take ceiling if round=0
! */

function nicenum(x, round)
implicit none
real(8) nicenum,x;
integer round;

    integer expv;				!/* exponent of x */
    real(8) f;				!/* fractional part of x */
    real(8) nf;				!/* nice, rounded fraction */

    expv = floor(log10(x));
    f = x/10.**expv;		!/* between 1 and 10 */
    if (round) THEN
	    if (f<1.5) THEN
            nf = 1.;
	    else if (f<3.) THEN
            nf = 2.;
	    else if (f<7.) THEN
            nf = 5.;
	    else 
            nf = 10.;
        ENDIF
    else
	    if (f<=1.) THEN
            nf = 1.;
	    else if (f<=2.) THEN
            nf = 2.;
	    else if (f<=5.) THEN
            nf = 5.;
	    else 
            nf = 10.;
        ENDIF
    ENDIF
    
    NICENUM= nf*10.**expv;
end function



