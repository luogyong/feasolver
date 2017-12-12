subroutine DRAW_NADAL_VALUE()
    USE POS_IO
    use opengl_gl
    use opengl_glut
	use function_plotter
    implicit none
    integer i
    CHARACTER(16)::STR1 
    REAL(8)::X1,Y1,Z1
	
	IF(ISSHOWNODALVALUE) THEN
		call glPushAttrib(GL_ALL_ATTRIB_BITS)
		!call glPushAttrib(GL_LIGHTING_BIT .or. GL_CURRENT_BIT); ! lighting and color mask
		call glDisable(GL_LIGHTING);     ! need to disable lighting for proper text color
		call glDisable(GL_TEXTURE_2D);
		call glDisable(GL_CULL_FACE);
        call glEnable(GL_DEPTH_TEST)
		call glDepthFunc(GL_ALWAYS);
		
		CALL GLCOLOR4FV(MYCOLOR(:,BLACK))
        
        
		IF(SHOW_NODAL_VALUE/=-2) THEN
		    DO I=1,POSDATA.NNODE
			    IF(POSDATA.NODE(I).ISDEAD==1) CYCLE
                IF(SHOW_NODAL_VALUE>0) THEN
			        WRITE(STR1,10) POSDATA.NODALQ(I,SHOW_NODAL_VALUE,STEPPLOT.ISTEP)
                ELSEIF(SHOW_NODAL_VALUE==-1) THEN
                    WRITE(STR1,20) I !NODE ID
                ENDIF
			    CALL output3D(POSDATA.NODE(I).COORD(1),POSDATA.NODE(I).COORD(2),POSDATA.NODE(I).COORD(3),TRIM(ADJUSTL(STR1))) 
		    ENDDO
        ELSE
            !ELEMENT ID
		    DO I=1,POSDATA.NEL
                IF(POSDATA.ESET(POSDATA.ELEMENT(I).ISET).STEPSTATUS(STEPPLOT.ISTEP)==0) CYCLE 
			    
                WRITE(STR1,20) I !NODE ID
                X1=SUM(POSDATA.NODE(POSDATA.ELEMENT(I).NODE).COORD(1))/POSDATA.ELEMENT(I).NNUM
                Y1=SUM(POSDATA.NODE(POSDATA.ELEMENT(I).NODE).COORD(2))/POSDATA.ELEMENT(I).NNUM
                Z1=0.D0
                IF(POSDATA.NDIM>2) Z1=SUM(POSDATA.NODE(POSDATA.ELEMENT(I).NODE).COORD(3))/POSDATA.ELEMENT(I).NNUM                
			    CALL output3D(X1,Y1,Z1,TRIM(ADJUSTL(STR1))) 
		    ENDDO            
            
        ENDIF
		
	   
		
		CALL glPopAttrib();
    ENDIF
    call glutPostRedisplay
    
10  FORMAT(G16.4)
20  FORMAT(I8)
end subroutine
