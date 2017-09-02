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

        do i=1,nedge
            do j=1,2
                !call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,vcolor(:,face(i).v(j)))
                call glvertex3dv(node(edge(i).v(j)).coord+VEC1(:,edge(i).v(j))*scale1)            
            enddo    
        enddo

    call glEnd

    call glEndList

endif

call glutPostRedisplay

endsubroutine