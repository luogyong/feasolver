subroutine DrawSurfaceContour()
use solverds
use opengl_gl
use function_plotter
use MESHGEO
use view_modifier
implicit none
REAL(GLDOUBLE)::VEC1(3,NNUM),SCALE1 
REAL(8),EXTERNAL::VECTOR_SCALE
real(GLFLOAT),allocatable::vcolor(:,:)


integer :: i,j,k,n1,MAT1(10)  !,nfrac

! prepare to make a new display list




VEC1=0;SCALE1=1.D0
IF(IsContour_In_DeformedMesh.AND.OUTVAR(DISX).VALUE>0.AND.OUTVAR(DISY).VALUE>0) THEN
    VEC1(1,:)=NODALQ(:,OUTVAR(DISX).IVO)
    VEC1(2,:)=NODALQ(:,OUTVAR(DISY).IVO)
    IF(NDIMENSION>2) VEC1(3,:)=NODALQ(:,OUTVAR(DISZ).IVO)
    scale1=VECTOR_SCALE(VEC1,NNUM,Scale_Deformed_Grid)
ENDIF

n1=outvar(CONTOUR_PLOT_VARIABLE).ivo

call glDeleteLists(contourlist, 1_glsizei)

if (draw_surface_solid) then
    call reset_view
    call glNewList(contourlist, gl_compile_and_execute)

! draw the solid surface

    call glPolygonMode(gl_front_and_back, gl_fill)
    call glBegin(gl_triangles)

   
    if (surface_color == rainbow_surface) then
        if(allocated(vcolor)) deallocate(vcolor)
        allocate(vcolor(4,NNUM))

        
        do i=1,nnum
            call get_rainbow(nodalq(i,n1),contourbar.val(1),contourbar.val(contourbar.nval),vcolor(:,i))
            if(isTransparency) vcolor(4,i)=0.6
        enddo
 
    
    endif

    do i=1,nface
        !if(face(i).shape/=3) cycle
		!只对外边界及材料边界进行渲染。
        IF(FACE(I).ENUM>1) THEN
			MAT1(1:FACE(I).ENUM)=ELEMENT(TET(ABS(FACE(I).ELEMENT)).MOTHER).MAT !
			IF(all(MAT1(1:FACE(I).ENUM)-MAT1(1)==0)) CYCLE
		ENDIF
        do j=1,3
            !call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,vcolor(:,face(i).v(j)))
            IF(surface_color == white_surface)THEN
                call glcolor4fv(GRAY)
            ELSE
                call glcolor4fv(vcolor(:,face(i).v(j)))
            ENDIF
            call glvertex3dv(node(face(i).v(j)).coord+VEC1(:,face(i).v(j))*SCALE1)            
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
endif ! draw_surface_solid

call glEndList

call glutPostRedisplay

if(allocated(vcolor)) deallocate(vcolor)
!if(allocated(contour_value)) deallocate(contour_value)

endsubroutine  


subroutine DrawLineContour()
use solverds
use opengl_gl
use function_plotter
use MESHGEO
use view_modifier
implicit none
real(GLDOUBLE)::CPNT1(3,NEDGE)    
INTEGER,allocatable::CLINE1(:,:),TRI1(:,:)
INTEGER::NCL1=0,NTRI1=0

integer :: i,j,k,n1  !,nfrac

! prepare to make a new display list

interface

	subroutine ContourLine(VA,nVA,VC,LINE,NLINE,CONTOURPOINT,TRI,NTRI)
		use MESHGEO
		implicit none
		integer,intent(in)::nVA
		real(8),intent(in)::VA(nVA),VC
		integer,intent(out)::NLINE,NTRI
		INTEGER,allocatable,intent(out)::LINE(:,:),TRI(:,:)
		real(8),intent(out)::CONTOURPOINT(3,NEDGE)
	end subroutine

endinterface



n1=outvar(CONTOUR_PLOT_VARIABLE).ivo

call glDeleteLists(contourLinelist, 1_glsizei)

if (draw_contour) then
    call reset_view
    call glNewList(contourLinelist, gl_compile_and_execute)
    call gldisable(GL_CULL_FACE);
	do i=1,contourbar.nval
        CPNT1=0.0d0
		call ContourLine(NODALQ(:,N1),NNUM,contourbar.VAL(I),CLINE1,NCL1,CPNT1,TRI1,NTRI1)
		call glcolor4fv(BLACK)		
		CALL GLBEGIN(GL_LINES)
			DO J=1,NCL1
				call glvertex3dv(CPNT1(:,CLINE1(1,J))) 
				call glvertex3dv(CPNT1(:,CLINE1(2,J))) 
			ENDDO
		CALL GLEND()
        call glcolor4fv(contourbar.Color(:,i))
        CALL GLBEGIN(GL_TRIANGLES)
			DO J=1,NTRI1
				call glvertex3dv(CPNT1(:,tri1(1,J))) 
				call glvertex3dv(CPNT1(:,TRI1(2,J)))
                call glvertex3dv(CPNT1(:,TRI1(3,J))) 
			ENDDO		
        
        CALL GLEND()
	enddo

endif


call glEndList

call glutPostRedisplay

endsubroutine 
    


subroutine initialize_contourplot()
    use solverds
    use function_plotter
    implicit none
    integer::i,n1
    real(GLDOUBLE) :: graphmin,graphmax,graphstep
    !real(GLDOUBLE),allocatable :: contour_value(:)
    character(128)::str1,str2  
    
    n1=outvar(CONTOUR_PLOT_VARIABLE).ivo
    minv=minval(nodalq(:,n1))
    maxv=maxval(nodalq(:,n1))
    WRITE(STR1,'(G13.3)') MAXV;WRITE(STR2,'(G13.3)') MINV;
        
    contourbar.title=trim(adjustl(outvar(CONTOUR_PLOT_VARIABLE).name))//',MAX='//TRIM(ADJUSTL(STR1))//',MIN='//TRIM(ADJUSTL(STR2))
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
    scale1=PPM1*3.527777778/119.05*0.8
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


subroutine ContourLine(VA,nVA,VC,LINE,NLINE,CONTOURPOINT,TRI,NTRI)
    use solverds
    use MESHGEO
    use function_plotter
    implicit none
    integer,intent(in)::nVA
    real(8),intent(in)::VA(nVA),VC
    integer,intent(out)::NLINE,NTRI
    INTEGER,allocatable,intent(out)::LINE(:,:),TRI(:,:)
    real(8),intent(out)::CONTOURPOINT(3,NEDGE)
    real(8)::V1,V2,T1,T2,T3,T4
    INTEGER::ISPVC1(NEDGE),E1(4)=0,E2(6)=0,E3(6)=0
    integer::i,j,MAXNLINE1=1000,N1,MAXNTRI1=1000,N2,N3,K
    REAL(GLDOUBLE)::VEC1(3,NNUM),SCALE1,X1(3),X2(3) 
	REAL(8),EXTERNAL::VECTOR_SCALE
    
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
	
	VEC1=0;SCALE1=1.D0
	IF(IsContour_In_DeformedMesh.AND.OUTVAR(DISX).VALUE>0.AND.OUTVAR(DISY).VALUE>0) THEN
		VEC1(1,:)=NODALQ(:,OUTVAR(DISX).IVO)
		VEC1(2,:)=NODALQ(:,OUTVAR(DISY).IVO)
		IF(NDIMENSION>2) VEC1(3,:)=NODALQ(:,OUTVAR(DISZ).IVO)
		scale1=VECTOR_SCALE(VEC1,NNUM,Scale_Deformed_Grid)
	ENDIF	
	
	
    ISPVC1=0
    do i=1,nedge
        V1=VA(EDGE(I).V(1));V2=VA(EDGE(I).V(2))
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
				X1=NODE(EDGE(I).V(1)).COORD+VEC1(:,EDGE(I).V(1))*SCALE1
				X2=NODE(EDGE(I).V(2)).COORD+VEC1(:,EDGE(I).V(2))*SCALE1
                CONTOURPOINT(:,I)=X1+t4*(X2-X1)
            ENDIF   
        ENDIF 
    enddo
    NLINE=0
	
    DO I=1,NFACE
		!IF(ANY(TET(ABS(FACE(I).ELEMENT(1:FACE(I).ENUM))).MOTHER>0)) CYCLE !NOT PLANE ELEMENT
		IF(FACE(I).ENUM>1) CYCLE
		E1=ABS(FACE(I).EDGE)		
        !!!!!!!!!!!!!!!!!!!!
		IF(SUM(ISPVC1(E1(1:FACE(I).SHAPE)))/=2) CYCLE
        N1=0
        NLINE=NLINE+1
		IF(NLINE>MAXNLINE1) THEN
            CALL I2_ENLARGE_AR(LINE,500,2)
            MAXNLINE1=MAXNLINE1+500
        ENDIF        
        DO J=1,FACE(I).SHAPE        
            IF(ISPVC1(E1(J))==1) THEN              
                N1=N1+1
                LINE(N1,NLINE)=E1(J)       
            ENDIF
       ENDDO    
    ENDDO
    
	NTRI=0
	DO I=1,NTET
        IF(TET(I).DIM/=3) CYCLE
        
		E2=TET(I).E
		IF(SUM(ISPVC1(E2))<3) CYCLE
        N1=0
        NTRI=NTRI+1
		IF(NTRI+1>MAXNTRI1) THEN !至少有两个空间
            CALL I2_ENLARGE_AR(TRI,500,3)
            MAXNTRI1=MAXNTRI1+500
        ENDIF
        E3=0
        DO J=1,6      
            IF(ISPVC1(E2(J))==1) THEN              
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
