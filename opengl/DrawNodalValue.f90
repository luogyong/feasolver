subroutine DRAW_NADAL_VALUE()
    USE POS_IO
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
        
        
		
		DO I=1,POSDATA.NNODE
			IF(POSDATA.NODE(I).ISDEAD==1) CYCLE
			WRITE(STR1,10) POSDATA.NODALQ(I,SHOW_NODAL_VALUE,STEPPLOT.ISTEP)
			CALL output3D(POSDATA.NODE(I).COORD(1),POSDATA.NODE(I).COORD(2),POSDATA.NODE(I).COORD(3),TRIM(ADJUSTL(STR1))) 
		ENDDO
		
	   
		
		CALL glPopAttrib();
    ENDIF
    call glutPostRedisplay
    
10  FORMAT(G16.4)    
end subroutine
