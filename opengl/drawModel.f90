subroutine drawGrid()
use solverds
use opengl_gl
use function_plotter
use MESHGEO
implicit none    
integer::i,j

call glDeleteLists(gridlist, 1_glsizei)

if (.not.draw_surface_grid) return

call reset_view 
call glNewList(gridlist, gl_compile_and_execute)

call glPolygonMode(gl_front_and_back, gl_line)
call glBegin(gl_lines)

    call glcolor4fv(black)

    do i=1,nedge
        do j=1,2
            !call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,vcolor(:,face(i).v(j)))
            call glvertex3dv(node(edge(i).v(j)).coord)            
        enddo    
    enddo

call glEnd

call glEndList

call glutPostRedisplay

endsubroutine