 subroutine ProbeatPoint(x,y)

     use function_plotter
    
	implicit none    
    integer(kind=glcint),intent(in) ::  x, y
    real(8)::pt1(3)
    integer::iel

!    INTEGER,EXTERNAL::PTINTRIlOC,POINTlOC
    
    
    isProbeState=.TRUE.	
    call glutSetCursor(GLUT_CURSOR_CROSSHAIR) 
    

    call GetOGLPos(x, y,Pt1)

    !PT1(1)=19.2507603482789
    !PT1(2)=20.2225219827355
    !PT1(3)=0.D0
    
    iel=POINTlOC(pt1,0)
	
    IF(iel>0) then
        call ProbeShow(Pt1,iel)
        PT_PROBE=PT1;IEL_PROBE=IEL
    else
        info.str='the picked location is out of zone.\nPrss q to exit'C
        info.color=red;info.qkey=.true.
    endif


endsubroutine

subroutine ProbeatPhyscialspace(pt,val,iel)

    use function_plotter
    
	implicit none    
    real(8),intent(in)::pt(3)
    integer,intent(in out)::iel
    real(8),intent(in out)::val(POSDATA.NVAR)    
!    INTEGER,EXTERNAL::PTINTRIlOC,POINTlOC,POINTlOC_BC
   
    !iel=0
    iel=POINTlOC_BC(pt,iel)
	val=0.d0
    IF(iel>0) then
        call getval(Pt,iel,val)
    else        
        info.str='the Point is out of zone.\nPrss q to exit'C
        info.color=red;info.qkey=.true.
    endif


endsubroutine


SUBROUTINE ProbeShow(Pt1,iel)
    
    use function_plotter
    implicit none
    real(8),intent(in)::pt1(3)
    integer,intent(in)::iel
    integer::I,NTI1
    real(8)::val1(POSDATA.NVAR),POS1(2)
    CHARACTER(128)::TITLE1(POSDATA.NVAR+20),STR1
	CHARACTER(32)::CWORD1,CWORD2
    
    NTI1=0
    call getval(Pt1,iel,val1)
    VAL_PROBE(1:POSDATA.NVAR)=VAL1
	NTI1=NTI1+1
	WRITE(CWORD1,*) STEPPLOT.ISTEP
	WRITE(CWORD2,*) STEPPLOT.NSTEP
	TITLE1(NTI1)="ISTEP/NSTEP:"//TRIM(ADJUSTL(CWORD1))//'/'//TRIM(ADJUSTL(CWORD2))																
	WRITE(CWORD1,'(F10.5)') STEPPLOT.TIME(STEPPLOT.ISTEP)
	WRITE(CWORD2,'(F10.5)') STEPPLOT.TIME(STEPPLOT.NSTEP)
	NTI1=NTI1+1
	TITLE1(NTI1)="ITIME/NTIME:"//TRIM(ADJUSTL(CWORD1))//'/'//TRIM(ADJUSTL(CWORD2))
    NTI1=NTI1+1
    TITLE1(NTI1)="PROBE AT POINT(X,Y):"
    IF(POSDATA.NDIM>2) TITLE1(NTI1)="PROBE AT POINT(X,Y,Z):"
    NTI1=NTI1+1
    DO I=1,POSDATA.NDIM
        WRITE(STR1,10) PT1(I)
        IF(I==1) THEN
            TITLE1(NTI1)='('//TRIM(ADJUSTL(STR1))
        ELSE
            TITLE1(NTI1)=TRIM(ADJUSTL(TITLE1(NTI1)))//','//TRIM(ADJUSTL(STR1))
        ENDIF
        IF(I==POSDATA.NDIM) TITLE1(NTI1)=TRIM(ADJUSTL(TITLE1(NTI1)))//')'
    ENDDO
    DO I=1,POSDATA.NVAR
		WRITE(STR1,10) VAL1(I)
		NTI1=NTI1+1
        TITLE1(NTI1)=TRIM(ADJUSTL(POSDATA.OUTVAR(I).NAME))//':'//TRIM(ADJUSTL(STR1))
    ENDDO
    NTI1=NTI1+1
    WRITE(STR1,20) TET(IEL).MOTHER
    TITLE1(NTI1)="ELEMENT ID:"//':'//TRIM(ADJUSTL(STR1))
    NTI1=NTI1+1
    WRITE(STR1,20) POSDATA.ELEMENT(TET(IEL).MOTHER).ISET
    TITLE1(NTI1)="ELEMENT SET:"//':'//TRIM(ADJUSTL(STR1))
        
	POS1(1)=0.01;POS1(2)=0.95
    
    CALL SHOW_MTEXT(TITLE1,NTI1,POS1,BLACK,ProbeValuelist)
    
10 FORMAT(G16.3)
20 FORMAT(I7)
ENDSUBROUTINE




 subroutine getval(Pt,itet,val)
    !use solverds
    !use MESHGEO
    USE function_plotter
    implicit none
    integer,intent(in)::itet
    real(8),intent(in)::Pt(3)
    real(8)::val(POSDATA.NVAR)
    real(8)::shpfun(4),val1(4)
    integer::i,n1
    
    call tetshapefun(Pt,itet,shpfun)
    n1=tet(itet).nv    
    do i=1,POSDATA.NVAR
        val1(1:n1)=POSDATA.nodalq(tet(itet).v(1:n1),i,stepplot.istep)
        val(i)=dot_product(val1(1:n1),shpfun(1:n1))
    enddo
    
 endsubroutine
    
subroutine GetOGLPos(x, y,object) 
    use opengl_gl
    use opengl_glu
    implicit none
    integer(glint),intent(in)::x,y
    real(gldouble)::object(3)
    integer(GLint),dimension(4):: viewport;
    real(GLdouble),dimension(16):: modelview;
    real(GLdouble),dimension(16):: projection;
    real(gldouble)::winX, winY, winZ=0.d0,clip_Z
    real(GLdouble)::posX, posY, posZ,nf(2),far_z,near_z
    integer(GLint)::i
    REAL(GLFLOAT)::WINZ1
 
    call glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
    call glGetDoublev( GL_PROJECTION_MATRIX, projection );
    call glGetIntegerv( GL_VIEWPORT, viewport );
    !print *, viewport
    
    winX = real(x);
    winY = real(viewport(4)) - real(y);
    
    call f9y1glReadPixels( x, int(winY), 1, 1 , GL_DEPTH_COMPONENT, GL_FLOAT, WINZ1 );

    WINZ=WINZ1 !!!NOTE THAT f9y1glReadPixels CANN'T ACCEPT GL_DOUBLE. it took me two days to find it. 
    
    i= gluUnProject( winX, winY, winZ, modelview, projection, viewport, posX, posY, posZ);
    object=[posX,posY,posZ]
    return
    !return CVector3(posX, posY, posZ);
end subroutine


subroutine GetWinXYZ(objx,objy,objz,WinP) 
    use opengl_gl
    use opengl_glu
    implicit none
    real(gldouble),intent(in)::objx,objy,objz
    real(gldouble)::WinP(3),winX,winY,winZ
    integer(GLint),dimension(4):: viewport;
    real(GLdouble),dimension(16):: modelview;
    real(GLdouble),dimension(16):: projection;
    integer(GLint)::i
 
    call glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
    call glGetDoublev( GL_PROJECTION_MATRIX, projection );
    call glGetIntegerv( GL_VIEWPORT, viewport );

    i= gluProject( objX, objY, objZ, modelview, projection, viewport, winX, winY, winZ);
    WinP=[winX, winY, winZ]
    return
    !return CVector3(posX, posY, posZ);
end subroutine
 
 
 


    


