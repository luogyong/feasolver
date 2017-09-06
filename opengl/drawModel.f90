subroutine drawGrid()
use solverds
use opengl_gl
use function_plotter
use MESHGEO
implicit none    
integer::i,j
real(gldouble)::VEC1(3,NNUM),maxv1,scale1
logical::isdg1=.false.
REAL(8),EXTERNAL::VECTOR_SCALE

call glDeleteLists(gridlist, 1_glsizei)

if (draw_surface_grid) then

    call reset_view 
    call glNewList(gridlist, gl_compile_and_execute)

    call glPolygonMode(gl_front_and_back, gl_line)
    VEC1=0;SCALE1=1.D0
    
    IF(ISDEFORMEDMESH.AND.OUTVAR(DISX).VALUE>0.AND.OUTVAR(DISY).VALUE>0) THEN
        VEC1(1,:)=NODALQ(:,OUTVAR(DISX).IVO)
        VEC1(2,:)=NODALQ(:,OUTVAR(DISY).IVO)
        IF(NDIMENSION>2) VEC1(3,:)=NODALQ(:,OUTVAR(DISZ).IVO)
        scale1=VECTOR_SCALE(VEC1,NNUM,Scale_Deformed_Grid)
        !MAXV1=MAXVAL(NORM2(VEC1,DIM=1))
        !scale1=modelr/40./maxv1*Scale_Deformed_Grid
    ENDIF
	
	

    call glBegin(gl_lines)

        call glcolor4fv(black)

		do i=1,nface

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
		
        !do i=1,nedge
        !    do j=1,2
        !        !call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,vcolor(:,face(i).v(j)))
        !        call glvertex3dv(node(edge(i).v(j)).coord+VEC1(:,edge(i).v(j))*scale1)            
        !    enddo    
        !enddo

    call glEnd

    call glEndList

endif

call glutPostRedisplay

endsubroutine
