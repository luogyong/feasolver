subroutine DRAW_NADAL_VALUE()
    USE solverds
    use opengl_gl
    use opengl_glut
	use function_plotter
    implicit none
    integer i
    CHARACTER(16)::STR1  
	
	IF(ISSHOWNODALVALUE) THEN
		call glPushAttrib(GL_ALL_ATTRIB_BITS)
		!call glPushAttrib(GL_LIGHTING_BIT .or. GL_CURRENT_BIT); ! lighting and color mask
		call glDisable(GL_LIGHTING);     ! need to disable lighting for proper text color
		call glDisable(GL_TEXTURE_2D);
		call glDisable(GL_CULL_FACE);
        call glEnable(GL_DEPTH_TEST)
		call glDepthFunc(GL_ALWAYS);
		
		CALL GLCOLOR4FV(MYCOLOR(:,BLACK))
        
        
		
		DO I=1,NNUM
			IF(VISDEAD(I)==1) CYCLE
			WRITE(STR1,10) NODALQ(I,OUTVAR(SHOW_NODAL_VALUE).IVO,STEPPLOT.ISTEP)
			CALL output3D(NODE(I).COORD(1),NODE(I).COORD(2),NODE(I).COORD(3),TRIM(ADJUSTL(STR1))) 
		ENDDO
		
	   
		
		CALL glPopAttrib();
    ENDIF
    call glutPostRedisplay
    
10  FORMAT(G16.4)    
end subroutine
