 subroutine ProbeatPoint(x,y)

     use function_plotter
    
	implicit none    
    integer(kind=glcint),intent(in) ::  x, y
    real(8)::pt1(3)
    integer::iel

    INTEGER,EXTERNAL::PTINTRIlOC,POINTlOC
    
    
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
    INTEGER,EXTERNAL::PTINTRIlOC,POINTlOC,POINTlOC_BC
   
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


subroutine tetshapefun(Pt,itet,shpfun) 
    !use solverds
    use MESHGEO
    use SolverMath
    implicit none
    integer,intent(in)::itet
    real(8),intent(in)::Pt(3)
    real(8),dimension(4)::shpfun
    real(8)::Va1(3,4),v1(3),v2(3),v3(3),vol1
    integer::i,vi1(3,4),n1
    !real(8),EXTERNAL::determinant
     
    do i=1,tet(itet).nv
        Va1(:,i)=POSDATA.NODE(tet(itet).v(i)).coord
    enddo
    if(tet(itet).dim==2) then
        v1=va1(:,2)-va1(:,1)
        v2=va1(:,3)-va1(:,1) 
        call r8vec_cross_3d ( v1, v2, v3 )
        n1=maxloc(abs(v3),dim=1)
        vol1=v3(n1)
        if(abs(vol1)<1e-7) then
            print *, "Error. the tet(itet) Area is 0. sub=tetshapefun,itet=",itet
            stop
        endif
        shpfun(3)=1.0d0
        do i=1,2
            v1=va1(:,mod(i,3)+1)-Pt
            v2=va1(:,mod(i+1,3)+1)-Pt
            call r8vec_cross_3d ( v1, v2, v3 )
            !shpfun(i)=max(min(norm2(v3)/vol1,1.d0),0.d0) !0-1
            shpfun(i)=v3(n1)/vol1
            shpfun(3)=shpfun(3)-shpfun(i)
        enddo
        
        
    elseif(tet(itet).dim==3) then
		vi1=reshape([2,4,3,1,3,4,1,4,2,1,2,3],([3,4]))
	
        v1=va1(:,2)-va1(:,1)
        v2=va1(:,3)-va1(:,1) 
        v3=va1(:,4)-va1(:,1) 
        vol1=determinant(reshape([v1,v2,v3],([3,3])))
        shpfun(4)=1.d0
        do i=1,3
            v1=va1(:,vi1(1,i))-Pt
            v2=va1(:,vi1(2,i))-Pt
            v3=va1(:,vi1(3,i))-Pt
            !shpfun(i)=min(max(abs(determinant(reshape([v1,v2,v3],([3,3])))/vol1),0.d0),1.d0)
            shpfun(i)=determinant(reshape([v2,v1,v3],([3,3])))/vol1
            shpfun(4)=shpfun(4)-shpfun(i)
        enddo
    endif
    
    
    
endsubroutine

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
    do i=1,size(POSDATA.nodalq,dim=2)
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
 
 
 

!search method:element-by-element，it is robust but may be comparily slow.
!if tryiel>0, it the element and its nerghborhood will be searched firstly.
INTEGER FUNCTION POINTlOC(PT,TRYIEL)
	USE MESHGEO
	IMPLICIT NONE
	REAL(8),INTENT(IN)::PT(3)
    INTEGER,INTENT(IN)::TRYIEL
	INTEGER::I,J,K,ELT1(-10:0)
	LOGICAL,EXTERNAL::PtInTri,PtInTET
    LOGICAL::ISFOUND=.FALSE.
	!INTEGER::SEARCHLIST(NTET)
	!I=INT(1+(RANDOM(1))(NFACE-1))
    
    ISFOUND=.FALSE.    
    
    DO J=1,NSZONE
        IF(SEARCHZONE(J).NEL<1) CYCLE
        IF(SEARCHZONE(J).NDX(1)>1) THEN
            IF(PT(1)<SEARCHZONE(J).BBOX(1,1)) CYCLE
            IF(PT(1)>SEARCHZONE(J).BBOX(2,1)) CYCLE
        ENDIF
        IF(SEARCHZONE(J).NDX(2)>1) THEN
            IF(PT(2)<SEARCHZONE(J).BBOX(1,2)) CYCLE
            IF(PT(2)>SEARCHZONE(J).BBOX(2,2)) CYCLE
        ENDIF        
        IF(SEARCHZONE(J).NDX(3)>1) THEN
            IF(PT(3)<SEARCHZONE(J).BBOX(1,3)) CYCLE
            IF(PT(3)>SEARCHZONE(J).BBOX(2,3)) CYCLE
        ENDIF     
        
        IF(TRYIEL>0) THEN
            ELT1(-TET(TRYIEL).NV)=TRYIEL
            DO K=1,TET(TRYIEL).NV
                ELT1(-TET(TRYIEL).NV+K)=TET(TRYIEL).ADJELT(K)
            ENDDO
            K=-TET(TRYIEL).NV
        ELSE
            K=1
        ENDIF
        
	    do K=K,SEARCHZONE(J).NEL
            IF(K>0) THEN
                I=SEARCHZONE(J).ELEMENT(K) 
            ELSE
                I=ELT1(K)
            ENDIF
            IF(I<1) CYCLE
		    IF(TET(I).ISDEAD==1) CYCLE
		    IF(TET(I).BBOX(2,1)>=PT(1).and.TET(I).BBOX(1,1)<=PT(1)) then
			    IF(TET(I).BBOX(2,2)>=PT(2).and.TET(I).BBOX(1,2)<=PT(2)) then
				    IF(POSDATA.NDIM>2) THEN
					    IF(TET(I).BBOX(2,3)<PT(3).and.TET(I).BBOX(1,3)>PT(3)) CYCLE
				    ENDIF
				    IF(TET(I).DIM==2) THEN
					    IF(PtInTri(PT, POSDATA.NODE(TET(I).V(1)).COORD, POSDATA.NODE(TET(I).V(2)).COORD, POSDATA.NODE(TET(I).V(3)).COORD)) THEN
                            ISFOUND=.TRUE.
                            EXIT
                        ENDIF
				    ELSEIF(TET(I).DIM==3) THEN
					    IF(PtInTET(PT, POSDATA.NODE(TET(I).V(1)).COORD, POSDATA.NODE(TET(I).V(2)).COORD, POSDATA.NODE(TET(I).V(3)).COORD,POSDATA.NODE(TET(I).V(4)).COORD)) THEN
                            ISFOUND=.TRUE.
                            EXIT
                        ENDIF
				    ENDIF
			    ENDIF
		    endif
	    ENDDO
        

    ENDDO
	
	IF(.NOT.ISFOUND) I=0
    LASTLOC=I
    POINTlOC=I
	RETURN
ENDFUNCTION




INTEGER FUNCTION PTINTRIlOC(PT)
    !USE solverds
	USE MESHGEO
	IMPLICIT NONE
	REAL(8),INTENT(IN)::PT(3)
	INTEGER::I
	LOGICAL,EXTERNAL::PtInTri,PtInTET
    REAL(8)::TOF1
	
	!I=INT(1+(RANDOM(1))(NFACE-1))
    
	do i=1,NFACE
		IF(FACE(I).ISDEAD==1) CYCLE
		IF(FACE(I).BBOX(2,1)>=PT(1).and.FACE(I).BBOX(1,1)<=PT(1)) then
			IF(FACE(I).BBOX(2,2)>=PT(2).and.FACE(I).BBOX(1,2)<=PT(2)) then
				IF(POSDATA.NDIM>2) THEN
					IF(FACE(I).BBOX(2,3)<PT(3).and.FACE(I).BBOX(1,3)>PT(3)) CYCLE
				ENDIF
				IF(PtInTri(PT, POSDATA.NODE(FACE(I).V(1)).COORD, POSDATA.NODE(FACE(I).V(2)).COORD, POSDATA.NODE(FACE(I).V(3)).COORD)) EXIT
				
			ENDIF
		endif
	ENDDO
	
	IF(I>NFACE) I=0
	PTINTRIlOC=I
	RETURN
ENDFUNCTION

! 
! Barycentric coordinates search method.
!if fail, recheck by searching one-by-one method (POINTlOC)
integer function POINTlOC_BC(Pt,TRYiel)
    use MESHGEO
    USE IFPORT
    implicit none
    real(8),intent(in)::pt(3)
    integer,intent(in)::TRYiel
    real(8)::shpfun(4)
    integer::IEL,n1,N2,I
    INTEGER,EXTERNAL::POINTlOC
    integer::N2E1(1:4,3:4)=reshape([2,3,1,0,3,4,2,1],[4,2])
  
    N2=0
    POINTlOC_BC=0
    IEL=0
    if(TRYiel<1) THEN
        DO I=1,NSZONE
            IF(SEARCHZONE(I).NEL<1) CYCLE
            IF(SEARCHZONE(I).NDX(1)>1) THEN
                IF(PT(1)<SEARCHZONE(I).BBOX(1,1)) CYCLE
                IF(PT(1)>SEARCHZONE(I).BBOX(2,1)) CYCLE
            ENDIF
            IF(SEARCHZONE(I).NDX(2)>1) THEN
                IF(PT(2)<SEARCHZONE(I).BBOX(1,2)) CYCLE
                IF(PT(2)>SEARCHZONE(I).BBOX(2,2)) CYCLE
            ENDIF        
            IF(SEARCHZONE(I).NDX(3)>1) THEN
                IF(PT(3)<SEARCHZONE(I).BBOX(1,3)) CYCLE
                IF(PT(3)>SEARCHZONE(I).BBOX(2,3)) CYCLE
            ENDIF                
        
            iel=SEARCHZONE(I).ELEMENT(MAX(INT(rand(1)*SEARCHZONE(I).NEL),1))

            EXIT
        ENDDO
        IF(IEL==0) THEN
            RETURN
        ENDIF
        
        !iel=(MAX(INT(rand(1)*NTET),1))
    ELSE
        IEL=TRYiel
    ENDIF
    
    do while(.true.)
        !N2=N2+1
        call tetshapefun(Pt,iel,shpfun)
        n2=tet(iel).nv
        
        n1=minloc(shpfun(1:n2),dim=1)
        if(shpfun(n1)<0.d0) then            
            !iel=tet(iel).adjelt(mod(n1,3)+1)
            iel=tet(iel).adjelt(N2E1(n1,n2))
            !PRINT *, IEL,N2
            if(iel==0) then
                !RECHECK BY THE FINAL METHOD.
                POINTlOC_BC=POINTlOC(PT,TRYiel)                
                exit
            endif
            
        else
            POINTlOC_BC=iel
            exit
        endif
        !ISSEARCH(IEL)=IEL
    enddo
    
    !if(.not.isfound) iel=0
    
end function



logical function PtInTri (pt, v1, v2, v3)
	implicit none
	integer b1, b2
	real(8),intent(in)::pt(3),v1(3),v2(3),v3(3)
    integer,EXTERNAL::isacw
    
	PtInTri=.false.
    b1 = isacw(pt(1),pt(2),pt(3),v1(1),v1(2),v1(3),v2(1),v2(2),v2(3)) ;
    if(b1==2) then
        PtInTri=.true.
        return
    endif
    b2 = isacw(pt(1),pt(2),pt(3),v2(1),v2(2),v2(3),v3(1),v3(2),v3(3)) ;
    if(b2==2) then
        PtInTri=.true.
        return
    endif    
	if(b1/=b2) return
    b1 = isacw(pt(1),pt(2),pt(3),v3(1),v3(2),v3(3),v1(1),v1(2),v1(3)) ;
    if(b1==2) then
        PtInTri=.true.
        return
    endif    
    if(b1/=b2) return
    
	PtInTri=.true.	
end function



logical function PtInTet(pt, v1, v2, v3,v4)
    !假定四面体四个面的正向一致，都为正向(节点逆时针)或负向
	implicit none
	integer:: b1, b2
	real(8),intent(in)::pt(3),v1(3),v2(3),v3(3),v4(3)
	integer,EXTERNAL::ISFRONT
    
	PtInTet=.false.
    b1 = Isfront([v2,v1,v3,pt]) 
    if(b1==2) then
        PtInTet=.true.
        return
    endif
    b2 = Isfront([v1,v2,v4,pt]) 
    if(b2==2) then
        PtInTet=.true.
        return
    endif    
	if(b1/=b2) return
    b1 = Isfront([v2,v3,v4,pt])
    if(b1==2) then
        PtInTet=.true.
        return
    endif
	if(b1/=b2) return	
	b1 = Isfront([v3,v1,v4,pt])
    if(b1==2) then
        PtInTet=.true.
        return
    endif    
	if(b1/=b2) return
    
    PtInTet=.true.

end function

integer function Isfront(V)
!v(:,1-3) are face, and v(:,4) is point to be tested
    use SolverMath
	implicit none
	real(8),intent(in)::V(3,4)
	real(8)::V1(3,3)
	integer::I
    real(8)::T1
    logical,external::PtInTri
    !real(8),EXTERNAL::determinant
    
    Isfront=0
	do i=1,3
		v1(:,i)=v(:,i+1)-v(:,1)
	enddo
    
    t1=determinant(v1)
    
	if(t1>0.d0) Isfront=1
    
    if(abs(t1)<1e-10) then
        ISFRONT=3 !4点共面
        if(PtInTri (v(:,4), v(:,1), v(:,2), v(:,3))) Isfront=2 !on the TRIsurface        
    endif
    

end function

integer function isacw(x1,y1,z1,x2,y2,z2,x3,y3,z3)
	implicit none
	real(8),intent(in)::x1,y1,z1,x2,y2,z2,x3,y3,z3
    real(8)::t1,yn2,xn2,zn2,yn3,xn3,zn3,norm1(3),t2
	
	isacw=0
    !isacw=1,=2,onedge,=3,coline;
    
	yn2=y2-y1
	xn2=x2-x1
    zn2=z2-z1
	yn3=y3-y1
	xn3=x3-x1    
    zn3=z3-z1
    norm1=[yn2*zn3-zn2*yn3,-(xn2*zn3-zn2*xn3),xn2*yn3-yn2*xn3]
    t2=norm2(norm1)
    if(t2<1e-10) then !共线
        ISACW=3
        if(x1<min(x3,x2)) return
        if(X1>max(x3,x2)) return
        if(Y1<min(y3,y2)) return
        if(Y1>max(y3,y2)) return 
        if(Z1<min(z3,z2)) return
        if(Z1>max(z3,z2)) return 
        isacw=2 !on the edge
    else
        !t1=(xn2*yn3-yn2*xn3)+(yn2*zn3-zn2*yn3)-(xn2*zn3-zn2*xn3)
        !if(abs(t1)<1.d-10) then 
            !与(1,1,1)垂直            
        if(abs(norm1(3))>1e-10) then
            isacw=sign(1.,norm1(3)) !从+z看
        elseif(abs(norm1(1))>1e-10) then
            isacw=sign(1.,norm1(1)) !从+x看
        elseif(abs(norm1(2))>1e-10) then
            isacw=sign(1.,norm1(2)) !从+y看
        endif
        !else
            !从(1,1,1)方向看
        !    isacw=sign(1.,t1)
        !endif
    endif
	

end function
    


