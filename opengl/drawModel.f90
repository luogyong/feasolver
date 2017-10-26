subroutine drawGrid()
!use solverds
use opengl_gl
use function_plotter
use MESHGEO
implicit none    
integer::i,j,k,MAT1(10)
REAL(8),ALLOCATABLE::VEC1(:,:)
real(gldouble)::maxv1,scale1
logical::isdg1=.false.
REAL(8),EXTERNAL::VECTOR_SCALE

call glDeleteLists(gridlist, 1_glsizei)
call reset_view 

call glPushAttrib(GL_ALL_ATTRIB_BITS)

ALLOCATE(VEC1(3,POSDATA.NNODE))
VEC1=0;SCALE1=1.D0
IF(ISDEFORMEDMESH.AND.POSDATA.IDISX>0.AND.POSDATA.IDISY>0) THEN
	VEC1(1,:)=POSDATA.NODALQ(:,POSDATA.IDISX,STEPPLOT.ISTEP)
	VEC1(2,:)=POSDATA.NODALQ(:,POSDATA.IDISY,STEPPLOT.ISTEP)
	IF(POSDATA.NDIM>2) VEC1(3,:)=POSDATA.NODALQ(:,POSDATA.IDISZ,STEPPLOT.ISTEP)
	scale1=Scale_Deformed_Grid*STEPPLOT.VSCALE(1)
	!MAXV1=MAXVAL(NORM2(VEC1,DIM=1))
	!scale1=modelr/40./maxv1*Scale_Deformed_Grid
ENDIF

call glNewList(gridlist, gl_compile_and_execute)
	if (draw_surface_grid) then

		do k=1,2
		    if(k==1) cycle
			if(k==1) then
                call glPolygonMode(gl_front_and_back, gl_fill)

            endif
			if(k==2) then
                call glPolygonMode(gl_front_and_back, gl_line)
                !call glenable(gl_polygon_offset_line)
                !call glPolygonoffset(-1.0_glfloat,-1.0_glfloat)
            endif
		
			call glBegin(gl_triangles)
				
	 
				do i=1,nmface
					if(mface(i).isdead==1) cycle
					if(mface(i).shape/=3) cycle
					!只对外边界及材料边界进行渲染。
					MAT1(1:MFACE(I).ENUM)=POSDATA.ELEMENT(ABS(MFACE(I).ELEMENT)).ISET !
					IF(MFACE(I).ENUM>1) THEN				    
						IF(all(MAT1(1:MFACE(I).ENUM)-MAT1(1)==0)) CYCLE
					ENDIF
					
					if(k==1) then
						call glcolor4fv(mycolor(:,max(mod(mat1(1)+2,139),3)))
					else
						call glcolor4fv(mycolor(:,forest_green))
					endif			
					do j=1,MFACE(I).SHAPE
						!call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,vcolor(:,MFACE(I).v(j)))
					
						call glvertex3dv(POSDATA.node(MFACE(I).v(j)).coord+VEC1(:,MFACE(I).v(j))*SCALE1)            
					enddo
					
	            enddo
			call glEnd
		
			call glBegin(GL_QUADS)

				do i=1,nmface
					if(mface(i).isdead==1) cycle
					if(mface(i).shape/=4) cycle
					!只对外边界及材料边界进行渲染。
					MAT1(1:MFACE(I).ENUM)=POSDATA.ELEMENT(ABS(MFACE(I).ELEMENT)).ISET !
					IF(MFACE(I).ENUM>1) THEN
						IF(all(MAT1(1:MFACE(I).ENUM)-MAT1(1)==0)) CYCLE
					ENDIF
					if(k==1) then
						call glcolor4fv(mycolor(:,max(mod(mat1(1)+2,139),3)))
					else
						call glcolor4fv(mycolor(:,black))
					endif
					do j=1,MFACE(I).SHAPE
						!call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,vcolor(:,MFACE(I).v(j)))
					
						call glvertex3dv(POSDATA.node(MFACE(I).v(j)).coord+VEC1(:,MFACE(I).v(j))*SCALE1)            
					enddo    
				enddo
			
			call glEnd
        enddo
    end if
	
    
    
	if(show_edge) then
		call glPolygonMode(gl_front_and_back, gl_line)
		call glcolor4fv(mycolor(:,Navy))

		call glbegin(gl_lines)
			do i=1,nmedge
				if(medge(I).isdead==1) cycle
				do j=1,2
					!call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,vcolor(:,face(i).v(j)))
					call glvertex3dv(POSDATA.node(medge(i).v(j)).coord+VEC1(:,medge(i).v(j))*scale1)            
				enddo    
			enddo
		
		call glend

	endif
	
	if(show_node) then
	
		call glcolor4fv(mycolor(:,orange))
		call glPointSize(4.0_glfloat)
		call glbegin(gl_points)
		do i=1,POSDATA.NNODE
			if(POSDATA.NODE(I).isdead==1) cycle
			call glvertex3dv(POSDATA.node(i).coord+VEC1(:,i)*scale1)            
		enddo
			
		call glend

	endif

call glEndList
    
call glPopAttrib();

call glutPostRedisplay

DEALLOCATE(VEC1)

endsubroutine
