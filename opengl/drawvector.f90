subroutine DrawVector()
use solverds
use opengl_gl
use function_plotter
implicit none 
INTEGER::I,J
REAL(8)::SCALE1,MAXV1

MAXV1=MAXVAL(NORM2(VEC,DIM=1))
scale1=modelr/40./maxv1*Scale_Vector_len

call glDeleteLists(VectorList, 1_glsizei)
if(IsDrawVector) then
    call reset_view
    call glNewList(VectorList, gl_compile_and_execute)

    DO I=1,NNUM
        call drawArrow(node(i).coord,node(i).coord+vec(:,i)*scale1)    
    ENDDO

endif

call glEndList

call glutPostRedisplay

endsubroutine
    
    
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

