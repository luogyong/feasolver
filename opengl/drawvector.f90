subroutine DrawVector()
use POS_IO
use opengl_gl
use function_plotter
implicit none 
INTEGER::I,J
REAL(8)::GRIDSCALE1,SCALE1,MAXV1,minv1,V1(3)
REAL(8),EXTERNAL::Vector_SCALE
REAL(8),ALLOCATABLE::VEC1(:,:)

!MAXV1=MAXVAL(NORM2(VEC,DIM=1))
!scale1=modelr/40./maxv1*Scale_Vector_len


call glDeleteLists(VectorList, 1_glsizei)

if(IsDrawVector) then
	ALLOCATE(VEC1(3,POSDATA.NNODE))
    GRIDSCALE1=1.D0;VEC1=0.D0
    IF(IsContour_In_DeformedMesh.AND.POSDATA.IDISX>0.AND.POSDATA.IDISY>0) THEN

		VEC1(1,:)=POSDATA.NODALQ(:,POSDATA.IDISX,STEPPLOT.ISTEP)
		VEC1(2,:)=POSDATA.NODALQ(:,POSDATA.IDISY,STEPPLOT.ISTEP)
		IF(POSDATA.NDIM>2) VEC1(3,:)=POSDATA.NODALQ(:,POSDATA.IDISZ,STEPPLOT.ISTEP)

        GRIDSCALE1=STEPPLOT.VSCALE(1)*Scale_Deformed_Grid
    ENDIF


	SCALE1=Scale_Vector_len*STEPPLOT.VSCALE(VECTOR_PLOT_GROUP)
    call reset_view
    
    call glNewList(VectorList, gl_compile_and_execute)

    DO I=1,POSDATA.NNODE
        IF(POSDATA.NODE(I).ISDEAD==1) CYCLE
		V1=POSDATA.node(i).coord+VEC1(:,I)*GRIDSCALE1
        call drawArrow(V1,V1+POSDATA.VEC(:,i)*scale1)    
    ENDDO
    
    !drawVectorLegend.
    VabsMax=STEPPLOT.VMAX(VECTOR_PLOT_GROUP);VabsMin=STEPPLOT.VMIN(VECTOR_PLOT_GROUP)
    Vscale=scale1
    

endif

call glEndList

IF(ALLOCATED(VEC1)) DEALLOCATE(VEC1)

call glutPostRedisplay

endsubroutine

SUBROUTINE VEC_PLOT_DATA()
	USE POS_IO
	USE function_plotter
	IMPLICIT NONE
	
	select case(VECTOR_PLOT_GROUP)

	case(VECTOR_GROUP_DIS)
		POSDATA.VEC(1,:)=POSDATA.NODALQ(:,POSDATA.IDISX,STEPPLOT.ISTEP)
		POSDATA.VEC(2,:)=POSDATA.NODALQ(:,POSDATA.IDISY,STEPPLOT.ISTEP)
		IF(POSDATA.NDIM>2) THEN
			POSDATA.VEC(3,:)=POSDATA.NODALQ(:,POSDATA.IDISZ,STEPPLOT.ISTEP)
		ELSE
			POSDATA.VEC(3,:)=0.D0
		ENDIF	
		VectorPairName='DIS.'
	case(VECTOR_GROUP_SEEPAGE_VEC)       
		POSDATA.VEC(1,:)=POSDATA.NODALQ(:,POSDATA.IVX,STEPPLOT.ISTEP)
		POSDATA.VEC(2,:)=POSDATA.NODALQ(:,POSDATA.IVY,STEPPLOT.ISTEP)
		IF(POSDATA.NDIM>2) THEN
			POSDATA.VEC(3,:)=POSDATA.NODALQ(:,POSDATA.IVZ,STEPPLOT.ISTEP)
		ELSE
			POSDATA.VEC(3,:)=0.D0
		ENDIF
		VectorPairName='SEEP.V'
	case(VECTOR_GROUP_SEEPAGE_GRAD)

		POSDATA.VEC(1,:)=POSDATA.NODALQ(:,POSDATA.IGRADX,STEPPLOT.ISTEP)
		POSDATA.VEC(2,:)=POSDATA.NODALQ(:,POSDATA.IGRADY,STEPPLOT.ISTEP)
		IF(POSDATA.NDIM>2) THEN
			POSDATA.VEC(3,:)=POSDATA.NODALQ(:,POSDATA.IGRADZ,STEPPLOT.ISTEP)
		ELSE
			POSDATA.VEC(3,:)=0.D0
		ENDIF    
		VectorPairName='SEEP.I'
	case(VECTOR_GROUP_SFR)

		POSDATA.VEC(1,:)=POSDATA.NODALQ(:,POSDATA.ISFR_SFRX,STEPPLOT.ISTEP)
		POSDATA.VEC(2,:)=POSDATA.NODALQ(:,POSDATA.ISFR_SFRY,STEPPLOT.ISTEP)
		POSDATA.VEC(3,:)=0.D0       
		VectorPairName='SFR'
	end select

ENDSUBROUTINE

subroutine drawVectorLegend(Vtmax,Vtmin,Scale,VTITLE)
use opengl_gl
use opengl_glut
use view_modifier
implicit none
real(8),intent(in)::Vtmax,Vtmin,Scale
CHARACTER(128),INTENT(IN)::VTITLE
integer(glCint),dimension(4)::viewport1
real(gldouble),dimension(16)::model
real(glfloat)::len1
real(gldouble)::t1
character(128)::str1,str2,STR3
    
    
    call glgetintegerv(gl_viewport,viewport1)
	call glViewport(10,100,80,80);

	call glMatrixMode(GL_PROJECTION);
    call glpushmatrix()
	call glLoadIdentity();
    t1=Vtmax*Scale
    len1=real(t1)
	call glOrtho(-t1,t1,-t1,t1,-1.0,1.0);

	!// Strip translation
    call glMatrixMode(GL_MODELVIEW);
    call glpushmatrix()
    call glLoadIdentity();
    !call glGetDoublev(GL_MODELVIEW_MATRIX, model)
    !!Matrix4 mv = view;
    !model(13:15)=0.0d0
    !
    !if(.not.isperspect) then
    !    model(1:3)=model(1:3)/xscale_factor
    !    model(5:7)=model(5:7)/yscale_factor
    !    model(9:11)=model(9:11)/zscale_factor
    !endif
    !!call glMatrixMode(GL_MODELVIEW);
    !call glLoadMatrixd(reshape(model,(/4,4/)));

	! // Axes
    call glPushAttrib(GL_ALL_ATTRIB_BITS)
	!call glPushAttrib(GL_CURRENT_BIT);
	!call glPushAttrib(GL_LINE_BIT);
	!call glPushAttrib(GL_VIEWPORT_BIT)
    call glDisable(GL_LIGHTING);
    call glDisable(GL_CULL_FACE);
	call glLineWidth(3.0_glfloat);
    
    call drawArrow([-0.5*len1,0.,0.],[0.5*len1,0.,0.])
    


	!// Axes labels
    call glColor4d(0.,0.,0.0,1.0);
	call glRasterPos3f(-len1,-len1/2._glfloat,0._glfloat);
    WRITE(STR2,'(G13.5)') VTMAX
    WRITE(STR3,'(G13.5)') VTMIN
    STR1=TRIM(ADJUSTL(VTITLE))
    call output(-len1,-len1/4._glfloat,str1)
    STR1='Mag:'
    call output(-len1,-len1/4._glfloat*2._glfloat,str1) 
    STR1=' Max:'//TRIM(ADJUSTL(str2))
    call output(-len1,-len1/4._glfloat*3._glfloat,str1) 
    STR1=' Min:'//TRIM(ADJUSTL(str3))
    call output(-len1,-len1/4._glfloat*4._glfloat,str1)     
!	call glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,ichar('x'));


	call glPopAttrib();


	!// Restore viewport, projection and modelview matrices
	!call glViewport(viewport1(1),viewport1(2),viewport1(3),viewport1(4));
    call glViewport(0,0,GLUTGET(GLUT_WINDOW_WIDTH),GLUTGET(GLUT_WINDOW_HEIGHT));
	call glMatrixMode(GL_PROJECTION);
	call glpopmatrix()
	call glMatrixMode(GL_MODELVIEW);
	call glpopmatrix()    

    call glutPostRedisplay
    
endsubroutine


subroutine drawVectorLegend2(Vtmax,Vtmin,Scale,VTITLE)
    use opengl_gl
    use opengl_glut
    use view_modifier
    implicit none
    real(8),intent(in)::Vtmax,Vtmin,Scale
    CHARACTER(128),INTENT(IN)::VTITLE
    integer(glCint),dimension(4)::viewport1
    real(gldouble),dimension(16)::model
    real(gldouble)::len1,ppm1
    real(gldouble)::t1,orig(3),dest(3),winp1(3)
    character(128)::str1,str2,STR3


    
    !! Set up coordinate system to position color bar near bottom of window
    
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
    call glPolygonMode(gl_front_and_back, gl_line)
    
    orig(1)=7./8.*viewport1(3);orig(2)=9./10.*viewport1(4);orig(3)=0.
        !find scale!************i don't know how to scale 
    !t1=Vtmax*Scale 
    PPM1=glutget(GLUT_SCREEN_WIDTH)/REAL(glutget(GLUT_SCREEN_WIDTH_MM))
    t1=ppm1*10
    !t1=viewport1(3)/16; 
    
    dest=orig;dest(1)=dest(1)+t1
    call glLineWidth(3.0_glfloat);
    call glColor4d(1.,0.,0.0,1.0);
    !call glScaled(xscale_factor,yscale_factor,zscale_factor)
    call drawArrow(orig,dest)
    


	!// Axes labels
    call glColor4d(0.,0.,0.0,1.0);
	!call glRasterPos3d(orig(1),orig(2)-viewport1(4)/40.0,orig(3));
    WRITE(STR2,'(G13.5)') VTMAX
    WRITE(STR3,'(G13.5)') VTMIN
    STR1=TRIM(ADJUSTL(VTITLE))
    call output(real(orig(1),glfloat),real(orig(2)-viewport1(4)/60.0,glfloat),str1)
    STR1='Mag:'
    call output(real(orig(1),glfloat),real(orig(2)-2.*viewport1(4)/60.0,glfloat),str1) 
    STR1=' Max:'//TRIM(ADJUSTL(str2))
    call output(real(orig(1),glfloat),real(orig(2)-3.*viewport1(4)/60.0,glfloat),str1) 
    STR1=' Min:'//TRIM(ADJUSTL(str3))
    call output(real(orig(1),glfloat),real(orig(2)-4.*viewport1(4)/60.0,glfloat),str1)     
!	call glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,ichar('x'));
    
    CALL glPopMatrix()
    CALL glMatrixMode(GL_PROJECTION)
    CALL glPopMatrix()
	CALL glPopAttrib();
    
    call glutPostRedisplay
 
end subroutine

REAL(8) FUNCTION VECTOR_SCALE(VEC1,NV1,NSTEP1,SCALE1)
    USE POS_IO
    IMPLICIT NONE
    INTEGER,INTENT(IN)::NV1,NSTEP1
    REAL(8),INTENT(IN)::VEC1(3,NV1,NSTEP1),SCALE1
    REAL(8)::MAXV1
    MAXV1=MAXVAL(NORM2(VEC1,DIM=1))
    VECTOR_SCALE=POSDATA.modelr/40./maxv1*SCALE1
ENDFUNCTION



    
SUBROUTINE drawAxes()
use opengl_gl
use opengl_glut
use view_modifier
implicit none
integer(glCint),dimension(4)::viewport1
real(gldouble),dimension(16)::model
real(glfloat)::len1 = 0.5;
real(glfloat),dimension(3)::axesOrigin=[0.0,0.0,0.0];
    
    call glgetintegerv(gl_viewport,viewport1)
	call glViewport(10,10,80,80);

	call glMatrixMode(GL_PROJECTION);
    call glpushmatrix()
	call glLoadIdentity();
	call glOrtho(-0.5,0.5,-0.5,0.5,-1.0,1.0);

	!// Strip translation
    call glMatrixMode(GL_MODELVIEW);
    call glpushmatrix()
    call glGetDoublev(GL_MODELVIEW_MATRIX, model)
    !Matrix4 mv = view;
    model(13:15)=0.0d0
    
    if(.not.isperspect) then
        model(1:3)=model(1:3)/xscale_factor
        model(5:7)=model(5:7)/yscale_factor
        model(9:11)=model(9:11)/zscale_factor
    endif
    !call glMatrixMode(GL_MODELVIEW);
    call glLoadMatrixd(reshape(model,(/4,4/)));

	! // Axes
    call glPushAttrib(GL_ALL_ATTRIB_BITS)
	!call glPushAttrib(GL_CURRENT_BIT);
	!call glPushAttrib(GL_LINE_BIT);
	!call glPushAttrib(GL_VIEWPORT_BIT)
    call glDisable(GL_LIGHTING);
    call glDisable(GL_CULL_FACE);
	call glLineWidth(3.0_glfloat);
    

    
	call glBegin(GL_LINES);
    call glColor3d(1.,0.,0.);
	call glVertex3fv(axesOrigin);
	call glVertex3f(axesOrigin(1) + len1,axesOrigin(2),axesOrigin(3));
    call glColor3d(0.,1.0,0.);
	call glVertex3fv(axesOrigin);    
	call glVertex3f(axesOrigin(1),axesOrigin(2) + len1,axesOrigin(3));
    call glColor3d(0.,0.,1.);
	call glVertex3fv(axesOrigin);
	call glVertex3f(axesOrigin(1),axesOrigin(2),axesOrigin(3) + len1);
	call glEnd();

	!// Axes labels
	call glRasterPos3f(axesOrigin(1) + len1,axesOrigin(2),axesOrigin(3));
	call glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,ichar('x'));
	call glRasterPos3f(axesOrigin(1),axesOrigin(2) + len1,axesOrigin(3));
	call glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,ichar('y'));
	call glRasterPos3f(axesOrigin(1),axesOrigin(2),axesOrigin(3) + len1);
	call glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,ichar('z'));

	call glPopAttrib();
	!call glPopAttrib();
    !call glEnable(GL_LIGHTING);
    !call glEnable(GL_CULL_FACE);

	!// Restore viewport, projection and modelview matrices
	call glViewport(viewport1(1),viewport1(2),viewport1(3),viewport1(4));
	call glMatrixMode(GL_PROJECTION);
	call glpopmatrix()
	call glMatrixMode(GL_MODELVIEW);
	call glpopmatrix()
    
    call glutPostRedisplay
    
end subroutine

subroutine drawArrow(orig,dest)
use opengl_gl
use opengl_glut
implicit none
real(gldouble),intent(in)::orig(3),dest(3)
real(gldouble)::radius,height,length,PI,angle,rot(3),v1(3)
integer(glint)::slices=8, stacks=1

v1=dest-orig
length=norm2(v1)
IF(LENGTH<1E-10) RETURN
PI=ATAN(1.0D0)*4.0D0
call glPushAttrib(GL_ALL_ATTRIB_BITS)
call glColor3d(1.,0.,0.);
CALL GLLINEWIDTH(1._GLFLOAT)
call glbegin(gl_lines)
call glvertex3dv(orig)
call glvertex3dv(dest)
call glend()

!cone
height=length/4.0
radius=height*tan(30/180.*PI)
call r8vec_cross_3d ( [0.,0.,1.0],v1,rot )
angle=asin(norm2(rot)/length)/PI*180.
CALL glPolygonMode(gl_front_and_back, gl_fill)
call glColor4d(1.,0.,0.,0.8);
CALL glMatrixMode(GL_MODELVIEW);
CALL glPushMatrix()
CALL GLTRANSLATED(DEST(1),DEST(2),DEST(3))
CALL GLROTATED(ANGLE,ROT(1),ROT(2),ROT(3))
CALL glutSolidCone(radius,height,slices,stacks)
CALL glPopMatrix() 
call glPopAttrib()
    
endsubroutine

