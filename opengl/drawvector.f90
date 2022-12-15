subroutine DrawVector()
!use POS_IO
use opengl_gl
use function_plotter
implicit none 

INTEGER::NNODE
TYPE(node_tydef),ALLOCATABLE::NODE(:)
REAL(8),ALLOCATABLE::NODEDATA(:,:),VECTOR(:,:,:)


INTEGER::I,J
REAL(8)::GRIDSCALE1,SCALE1,MAXV1,minv1,V1(3),V2(3),T1
REAL(8),EXTERNAL::Vector_SCALE
REAL(8),ALLOCATABLE::VEC1(:,:)



!MAXV1=MAXVAL(NORM2(VEC,DIM=1))
!scale1=modelr/40./maxv1*Scale_Vector_len


call glDeleteLists(VectorList, 1_glsizei)

if(IsDrawVector) then
	IF(VECTORLOC==VECTORLOC_NODE) THEN
		NNODE=POSDATA.NNODE
		ALLOCATE(NODE,SOURCE=POSDATA.NODE)
		ALLOCATE(NODEDATA(POSDATA.NVAR,NNODE))
		DO I=1,NNODE
			NODEDATA(:,I)=POSDATA.NODALQ(I,:,STEPPLOT.ISTEP)
		ENDDO
		ALLOCATE(VECTOR,SOURCE=POSDATA.VEC)
	ELSE
		NNODE=POSDATA.NGRIDNODE
		ALLOCATE(NODE,SOURCE=POSDATA.GRIDNODE)
		ALLOCATE(NODEDATA,SOURCE=POSDATA.GRIDDATA)
		ALLOCATE(VECTOR,SOURCE=POSDATA.GRIDVEC)		
	ENDIF
	

	ALLOCATE(VEC1(3,NNODE))
    GRIDSCALE1=1.D0;VEC1=0.D0
    IF(ISDEFORMEDMESH.AND.POSDATA.IDISX>0.AND.POSDATA.IDISY>0) THEN

		VEC1(1,:)=NODEDATA(POSDATA.IDISX,:)
		VEC1(2,:)=NODEDATA(POSDATA.IDISY,:)
		IF(POSDATA.NDIM>2) VEC1(3,:)=NODEDATA(POSDATA.IDISZ,:)		
        GRIDSCALE1=STEPPLOT.VSCALE(1)*Scale_Deformed_Grid
    ENDIF


	
    call reset_view
    
    call glNewList(VectorList, gl_compile_and_execute)
	DO J=1,NVECTORPAIR
		IF(ACTIVE_VECTOR_GROUP(J)) THEN
			SCALE1=Scale_Vector_len*STEPPLOT.VSCALE(J)
            IF(ISVECTORUNIFY) 	SCALE1=SCALE1*STEPPLOT.VMAX(J)		
			DO I=1,NNODE,VectorFrequency
				IF(NODE(I).ISDEAD==1) CYCLE
                V2=VECTOR(:,I,J)
                IF(ISVECTORUNIFY) THEN
                    T1=NORM2(V2)                    
                    IF(ABS(T1)>1.D-10) THEN
                        V2=V2/T1
                    ELSE
                        CYCLE
                    ENDIF                    
                ENDIF
				V1=node(i).coord+VEC1(:,I)*GRIDSCALE1
				call drawArrow(V1,V1+V2*scale1,MYCOLOR(:,COLOR_TOP_TEN(J)),isarrow,vectorbase)    
			ENDDO
		ENDIF
    ENDDO
    !drawVectorLegend.
    !VabsMax=STEPPLOT.VMAX(VECTOR_PLOT_GROUP);VabsMin=STEPPLOT.VMIN(VECTOR_PLOT_GROUP)
    !Vscale=scale1
    

endif

call glEndList

call glutPostRedisplay

IF(ALLOCATED(VEC1)) DEALLOCATE(VEC1)
IF(ALLOCATED(NODE))	DEALLOCATE(NODE)
IF(ALLOCATED(NODEDATA))	DEALLOCATE(NODEDATA)
IF(ALLOCATED(VECTOR))	DEALLOCATE(VECTOR)	



endsubroutine

SUBROUTINE VEC_PLOT_DATA()
	!USE POS_IO
	USE function_plotter
	IMPLICIT NONE
	INTEGER::I,J,N1,N2,IEL1
    REAL(8)::PI1=3.141592653589793,T1,VAL1(POSDATA.NVAR),PT(3)
	!INTEGER,EXTERNAL::POINTlOC
	
	IF(VECTORLOC==2) THEN
		IEL1=0
		DO I=1,POSDATA.NGRIDNODE
			PT=POSDATA.GRIDNODE(I).COORD
			IF(POSDATA.GRIDNODE(I).IEL<1) THEN
				iel1=POINTlOC(PT,iel1)
			ELSE
                IEL1=POSDATA.GRIDNODE(I).IEL
				IF(TET(IEL1).ISDEAD==1) IEL1=0				
			ENDIF
			POSDATA.GRIDDATA(:,I)=0.d0
			IF(iel1>0) then
				call getval(Pt,iel1,POSDATA.GRIDDATA(:,I))
				POSDATA.GRIDNODE(I).ISDEAD=0				
			else        
				POSDATA.GRIDNODE(I).ISDEAD=1
			endif
		ENDDO
	
	ENDIF
	
	DO I=1,NVECTORPAIR
		IF(VECTORLOC==1) THEN !AT NODES
			IF (ACTIVE_VECTOR_GROUP(I)) THEN
				select case(I)
					
				case(VECTOR_GROUP_DIS)
					
					POSDATA.VEC(1,:,I)=POSDATA.NODALQ(:,POSDATA.IDISX,STEPPLOT.ISTEP)
					POSDATA.VEC(2,:,I)=POSDATA.NODALQ(:,POSDATA.IDISY,STEPPLOT.ISTEP)
					IF(POSDATA.NDIM>2) THEN
						POSDATA.VEC(3,:,I)=POSDATA.NODALQ(:,POSDATA.IDISZ,STEPPLOT.ISTEP)
					ELSE
						POSDATA.VEC(3,:,I)=0.D0
					ENDIF	
					!VectorPairName='DIS.'
				case(VECTOR_GROUP_SEEPAGE_VEC)       
					POSDATA.VEC(1,:,I)=POSDATA.NODALQ(:,POSDATA.IVX,STEPPLOT.ISTEP)
					POSDATA.VEC(2,:,I)=POSDATA.NODALQ(:,POSDATA.IVY,STEPPLOT.ISTEP)
					IF(POSDATA.NDIM>2) THEN
						POSDATA.VEC(3,:,I)=POSDATA.NODALQ(:,POSDATA.IVZ,STEPPLOT.ISTEP)
					ELSE
						POSDATA.VEC(3,:,I)=0.D0
					ENDIF
					!VectorPairName='SEEP.V'
				case(VECTOR_GROUP_SEEPAGE_GRAD)

					POSDATA.VEC(1,:,I)=POSDATA.NODALQ(:,POSDATA.IGRADX,STEPPLOT.ISTEP)
					POSDATA.VEC(2,:,I)=POSDATA.NODALQ(:,POSDATA.IGRADY,STEPPLOT.ISTEP)
					IF(POSDATA.NDIM>2) THEN
						POSDATA.VEC(3,:,I)=POSDATA.NODALQ(:,POSDATA.IGRADZ,STEPPLOT.ISTEP)
					ELSE
						POSDATA.VEC(3,:,I)=0.D0
					ENDIF    
					!VectorPairName='SEEP.I'
				case(VECTOR_GROUP_SFR)

					POSDATA.VEC(1,:,I)=POSDATA.NODALQ(:,POSDATA.ISFR_SFRX,STEPPLOT.ISTEP)
					POSDATA.VEC(2,:,I)=POSDATA.NODALQ(:,POSDATA.ISFR_SFRY,STEPPLOT.ISTEP)
                    IF(POSDATA.NDIM==2) THEN
                        POSDATA.VEC(3,:,I)=0.D0
                    ELSE                        
					    POSDATA.VEC(3,:,I)=POSDATA.NODALQ(:,POSDATA.ISFR_SFRZ,STEPPLOT.ISTEP) 
                    ENDIF
					!VectorPairName='SFR'
				case(VECTOR_GROUP_PSIGMA1)
					POSDATA.VEC(1,:,I)=POSDATA.NODALQ(:,POSDATA.IXPS1,STEPPLOT.ISTEP)
					POSDATA.VEC(2,:,I)=POSDATA.NODALQ(:,POSDATA.IYPS1,STEPPLOT.ISTEP)
					POSDATA.VEC(3,:,I)=POSDATA.NODALQ(:,POSDATA.IZPS1,STEPPLOT.ISTEP)
				case(VECTOR_GROUP_PSIGMA2)
					POSDATA.VEC(1,:,I)=POSDATA.NODALQ(:,POSDATA.IXPS2,STEPPLOT.ISTEP)
					POSDATA.VEC(2,:,I)=POSDATA.NODALQ(:,POSDATA.IYPS2,STEPPLOT.ISTEP)
					POSDATA.VEC(3,:,I)=POSDATA.NODALQ(:,POSDATA.IZPS2,STEPPLOT.ISTEP)
				case(VECTOR_GROUP_PSIGMA3)
					POSDATA.VEC(1,:,I)=POSDATA.NODALQ(:,POSDATA.IXPS3,STEPPLOT.ISTEP)
					POSDATA.VEC(2,:,I)=POSDATA.NODALQ(:,POSDATA.IYPS3,STEPPLOT.ISTEP)
					POSDATA.VEC(3,:,I)=POSDATA.NODALQ(:,POSDATA.IZPS3,STEPPLOT.ISTEP)                    
				end select
				
			ELSE
				POSDATA.VEC(:,:,I)=0.D0
			ENDIF
		ELSE
			!AT STRUCTURAL GRIDS
			IF (ACTIVE_VECTOR_GROUP(I)) THEN
				
				select case(I)
					
				case(VECTOR_GROUP_DIS)
					
					POSDATA.GRIDVEC(1,:,I)=POSDATA.GRIDDATA(POSDATA.IDISX,:)
					POSDATA.GRIDVEC(2,:,I)=POSDATA.GRIDDATA(POSDATA.IDISY,:)
					IF(POSDATA.NDIM>2) THEN
						POSDATA.GRIDVEC(3,:,I)=POSDATA.GRIDDATA(POSDATA.IDISZ,:)
					ELSE
						POSDATA.GRIDVEC(3,:,I)=0.D0
					ENDIF	
					!VectorPairName='DIS.'
				case(VECTOR_GROUP_SEEPAGE_VEC)       
					POSDATA.GRIDVEC(1,:,I)=POSDATA.GRIDDATA(POSDATA.IVX,:)
					POSDATA.GRIDVEC(2,:,I)=POSDATA.GRIDDATA(POSDATA.IVY,:)
					IF(POSDATA.NDIM>2) THEN
						POSDATA.GRIDVEC(3,:,I)=POSDATA.GRIDDATA(POSDATA.IVZ,:)
					ELSE
						POSDATA.GRIDVEC(3,:,I)=0.D0
					ENDIF
					!VectorPairName='SEEP.V'
				case(VECTOR_GROUP_SEEPAGE_GRAD)

					POSDATA.GRIDVEC(1,:,I)=POSDATA.GRIDDATA(POSDATA.IGRADX,:)
					POSDATA.GRIDVEC(2,:,I)=POSDATA.GRIDDATA(POSDATA.IGRADY,:)
					IF(POSDATA.NDIM>2) THEN
						POSDATA.GRIDVEC(3,:,I)=POSDATA.GRIDDATA(POSDATA.IGRADZ,:)
					ELSE
						POSDATA.GRIDVEC(3,:,I)=0.D0
					ENDIF    
					!VectorPairName='SEEP.I'
				case(VECTOR_GROUP_SFR)

					POSDATA.GRIDVEC(1,:,I)=POSDATA.GRIDDATA(POSDATA.ISFR_SFRX,:)
					POSDATA.GRIDVEC(2,:,I)=POSDATA.GRIDDATA(POSDATA.ISFR_SFRY,:)
                    IF(POSDATA.NDIM>2) THEN
						POSDATA.GRIDVEC(3,:,I)=POSDATA.GRIDDATA(POSDATA.ISFR_SFRZ,:)
					ELSE
						POSDATA.GRIDVEC(3,:,I)=0.D0
					ENDIF       
					!VectorPairName='SFR'
				case(VECTOR_GROUP_PSIGMA1)
					POSDATA.GRIDVEC(1,:,I)=POSDATA.GRIDDATA(POSDATA.IXPS1,:)
					POSDATA.GRIDVEC(2,:,I)=POSDATA.GRIDDATA(POSDATA.IYPS1,:)
                    POSDATA.GRIDVEC(3,:,I)=POSDATA.GRIDDATA(POSDATA.IZPS1,:)
				case(VECTOR_GROUP_PSIGMA2)
					POSDATA.GRIDVEC(1,:,I)=POSDATA.GRIDDATA(POSDATA.IXPS2,:)
					POSDATA.GRIDVEC(2,:,I)=POSDATA.GRIDDATA(POSDATA.IYPS2,:)
                    POSDATA.GRIDVEC(3,:,I)=POSDATA.GRIDDATA(POSDATA.IZPS2,:)
                case(VECTOR_GROUP_PSIGMA3)
					POSDATA.GRIDVEC(1,:,I)=POSDATA.GRIDDATA(POSDATA.IXPS3,:)
					POSDATA.GRIDVEC(2,:,I)=POSDATA.GRIDDATA(POSDATA.IYPS3,:)
                    POSDATA.GRIDVEC(3,:,I)=POSDATA.GRIDDATA(POSDATA.IZPS3,:)    
				end select
				

			ELSE
				POSDATA.GRIDVEC(:,:,I)=0.D0
			ENDIF
		
		ENDIF
	ENDDO
ENDSUBROUTINE



!subroutine drawVectorLegend(Vtmax,Vtmin,Scale,VTITLE)
!use opengl_gl
!use opengl_glut
!use view_modifier
!implicit none
!real(8),intent(in)::Vtmax,Vtmin,Scale
!CHARACTER(128),INTENT(IN)::VTITLE
!integer(glCint),dimension(4)::viewport1
!real(gldouble),dimension(16)::model
!real(glfloat)::len1
!real(gldouble)::t1
!character(128)::str1,str2,STR3
!    
!    
!    call glgetintegerv(gl_viewport,viewport1)
!	call glViewport(10,100,80,80);

!	call glMatrixMode(GL_PROJECTION);
!    call glpushmatrix()
!	call glLoadIdentity();
!    t1=Vtmax*Scale
!    len1=real(t1)
!	call glOrtho(-t1,t1,-t1,t1,-1.0,1.0);

!	!// Strip translation
!    call glMatrixMode(GL_MODELVIEW);
!    call glpushmatrix()
!    call glLoadIdentity();
!    !call glGetDoublev(GL_MODELVIEW_MATRIX, model)
!    !!Matrix4 mv = view;
!    !model(13:15)=0.0d0
!    !
!    !if(.not.isperspect) then
!    !    model(1:3)=model(1:3)/xscale_factor
!    !    model(5:7)=model(5:7)/yscale_factor
!    !    model(9:11)=model(9:11)/zscale_factor
!    !endif
!    !!call glMatrixMode(GL_MODELVIEW);
!    !call glLoadMatrixd(reshape(model,(/4,4/)));

!	! // Axes
!    call glPushAttrib(GL_ALL_ATTRIB_BITS)
!	!call glPushAttrib(GL_CURRENT_BIT);
!	!call glPushAttrib(GL_LINE_BIT);
!	!call glPushAttrib(GL_VIEWPORT_BIT)
!    call glDisable(GL_LIGHTING);
!    call glDisable(GL_CULL_FACE);
!	call glLineWidth(3.0_glfloat);
!    
!    call drawArrow([-0.5*len1,0.,0.],[0.5*len1,0.,0.])
!    


!	!// Axes labels
!    call glColor4d(0.,0.,0.0,1.0);
!	call glRasterPos3f(-len1,-len1/2._glfloat,0._glfloat);
!    WRITE(STR2,'(G13.5)') VTMAX
!    WRITE(STR3,'(G13.5)') VTMIN
!    STR1=TRIM(ADJUSTL(VTITLE))
!    call output(-len1,-len1/4._glfloat,str1)
!    STR1='Mag:'
!    call output(-len1,-len1/4._glfloat*2._glfloat,str1) 
!    STR1=' Max:'//TRIM(ADJUSTL(str2))
!    call output(-len1,-len1/4._glfloat*3._glfloat,str1) 
!    STR1=' Min:'//TRIM(ADJUSTL(str3))
!    call output(-len1,-len1/4._glfloat*4._glfloat,str1)     
!!	call glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,ichar('x'));


!	call glPopAttrib();


!	!// Restore viewport, projection and modelview matrices
!	!call glViewport(viewport1(1),viewport1(2),viewport1(3),viewport1(4));
!    call glViewport(0,0,GLUTGET(GLUT_WINDOW_WIDTH),GLUTGET(GLUT_WINDOW_HEIGHT));
!	call glMatrixMode(GL_PROJECTION);
!	call glpopmatrix()
!	call glMatrixMode(GL_MODELVIEW);
!	call glpopmatrix()    

!    call glutPostRedisplay
!    
!endsubroutine


subroutine drawVectorLegend2(Vtmax,Vtmin,Scale,ACTIVE_VECTOR_PAIR,VTITLE,NVG)
    use opengl_gl
    use opengl_glut
    use view_modifier
	USE INDEXCOLOR
    implicit none
	INTEGER,INTENT(IN)::NVG
	LOGICAL,INTENT(IN)::ACTIVE_VECTOR_PAIR(NVG)
    real(8),intent(in)::Vtmax(NVG),Vtmin(NVG),Scale(NVG)
    CHARACTER(128),INTENT(IN)::VTITLE(NVG)
    integer(glCint),dimension(4)::viewport1
    real(gldouble),dimension(16)::model
    real(gldouble)::len1,ppm1
    real(gldouble)::t1,orig(3),dest(3),winp1(3)
    character(128)::str1,str2,STR3
	INTEGER::I,N1

    
    !! Set up coordinate system to position color bar near bottom of window
    
    call glgetintegerv(gl_viewport,viewport1)
    CALL glMatrixMode(GL_PROJECTION);
    CALL glPushMatrix()
    CALL glLoadIdentity();
    CALL glOrtho(0._GLDOUBLE, REAL(viewport1(3),gldouble), 0._GLDOUBLE, REAL(viewport1(4),gldouble), -1._GLDOUBLE, 1._GLDOUBLE);
    CALL glMatrixMode(GL_MODELVIEW);
    CALL glPushMatrix()
   
    
    
    CALL glLoadIdentity();
    
    call glPushAttrib(GL_ALL_ATTRIB_BITS)
	!call glPushAttrib(GL_LIGHTING_BIT .or. GL_CURRENT_BIT); ! lighting and color mask
	call glDisable(GL_LIGHTING);     ! need to disable lighting for proper text color
	call glDisable(GL_TEXTURE_2D);
    call glDisable(GL_CULL_FACE);
	call glDepthFunc(GL_ALWAYS);    
    call glPolygonMode(gl_front_and_back, gl_line)
        !find scale!************i don't know how to scale 
    !t1=Vtmax*Scale 
	!t1=viewport1(3)/16;
    PPM1=glutget(GLUT_SCREEN_WIDTH)/REAL(glutget(GLUT_SCREEN_WIDTH_MM))
    t1=ppm1*10    
	N1=-1
	DO I=1,NVG
		IF(.NOT.ACTIVE_VECTOR_PAIR(I)) CYCLE
		N1=N1+1
 		orig(1)=7./8.*viewport1(3);orig(2)=(36.-N1)/40.*viewport1(4);orig(3)=0.
	
		dest=orig;dest(1)=dest(1)+t1
		call glLineWidth(3.0_glfloat);
		call glColor4FV(MYCOLOR(:,COLOR_TOP_TEN(I)));
		!call glScaled(xscale_factor,yscale_factor,zscale_factor)
		call drawArrow(orig,dest,MYCOLOR(:,COLOR_TOP_TEN(I)),.true.,1)
		


		!// Axes labels
		!call glColor4d(DISTINCT_COLOR(I));
		!call glRasterPos3d(orig(1),orig(2)-viewport1(4)/40.0,orig(3));
		WRITE(STR2,'(G13.5)') VTMAX(I)
		WRITE(STR3,'(G13.5)') VTMIN(I)
		STR1=TRIM(ADJUSTL(VTITLE(I)))//': MAX='//TRIM(STR2)//',MIN='//TRIM(STR3)
		call output(real(orig(1),glfloat),real(orig(2)-viewport1(4)/60.0,glfloat),str1)
		!STR1=STR1//
		!call output(real(orig(1),glfloat),real(orig(2)-2.*viewport1(4)/60.0,glfloat),str1) 
		!STR1=' Max:'//TRIM(ADJUSTL(str2))
		!call output(real(orig(1),glfloat),real(orig(2)-3.*viewport1(4)/60.0,glfloat),str1) 
		!STR1=' Min:'//TRIM(ADJUSTL(str3))
		!call output(real(orig(1),glfloat),real(orig(2)-4.*viewport1(4)/60.0,glfloat),str1)     
	!	call glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,ichar('x'));
    ENDDO
    CALL glPopMatrix()
    CALL glMatrixMode(GL_PROJECTION)
    CALL glPopMatrix()
	CALL glPopAttrib();
    
    call glutPostRedisplay
 
end subroutine

REAL(8) FUNCTION VECTOR_SCALE(VEC1,NV1,NSTEP1,SCALE1)
    USE POS_IO
    IMPLICIT NONE
    INTEGER,INTENT(IN)::NV1,NSTEP1
    REAL(8),INTENT(IN)::VEC1(3,NV1,NSTEP1),SCALE1
    REAL(8)::MAXV1
    MAXV1=MAXVAL(NORM2(VEC1,DIM=1))
    VECTOR_SCALE=POSDATA.modelr/40./maxv1*SCALE1
ENDFUNCTION



    
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
	call glOrtho(-0.5D0,0.5D0,-0.5D0,0.5D0,-1.0D0,1.0D0);

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
    call glColor3d(1._GLDOUBLE,0._GLDOUBLE,0._GLDOUBLE);
	call glVertex3fv(axesOrigin);
	call glVertex3f(axesOrigin(1) + len1,axesOrigin(2),axesOrigin(3));
    call glColor3d(0._GLDOUBLE,1.0_GLDOUBLE,0._GLDOUBLE);
	call glVertex3fv(axesOrigin);    
	call glVertex3f(axesOrigin(1),axesOrigin(2) + len1,axesOrigin(3));
    call glColor3d(0._GLDOUBLE,0._GLDOUBLE,1._GLDOUBLE);
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
    
    call glutPostRedisplay
    
end subroutine

subroutine drawArrow(orig,dest,color,isarrow,base)
use opengl_gl
use opengl_glut
implicit none
logical,intent(in)::isarrow
integer,intent(in)::base !=1 orig,=2,dest,=3,center
real(gldouble),intent(in)::orig(3),dest(3)
real(gldouble)::radius,height,length,PI,angle,rot(3),v1(3),orig1(3),dest1(3)
REAL(GLFLOAT)::color(4)
integer(glint)::slices=8, stacks=1


v1=dest-orig
length=norm2(v1)
IF(LENGTH<1E-10) RETURN

select case(base)
case(2)
    orig1=orig-v1
    dest1=dest-v1
case(3)
    orig1=orig-v1/2.
    dest1=dest-v1/2.
case default
    orig1=orig;dest1=dest
endselect




PI=ATAN(1.0D0)*4.0D0
call glPushAttrib(GL_ALL_ATTRIB_BITS)
call glColor4FV(color);
!call glColor3d(1.,0.,0.);
CALL GLLINEWIDTH(1._GLFLOAT)
call glbegin(gl_lines)
call glvertex3dv(orig1)
call glvertex3dv(dest1)
call glend()

if(isarrow) then
    !cone
    height=length/4.0
    radius=height*tan(30/180.*PI)

    call r8vec_cross_3d ( [0.D0,0.D0,1.0D0],v1,rot )
    angle=asin(norm2(rot)/length)/PI*180.
    CALL glPolygonMode(gl_front_and_back, gl_fill)
    !call glColor4d(1.,0.,0.,0.8);
    CALL glMatrixMode(GL_MODELVIEW);
    CALL glPushMatrix()
    CALL GLTRANSLATED(dest1(1),dest1(2),dest1(3))
    CALL GLROTATED(ANGLE,ROT(1),ROT(2),ROT(3))
    CALL glutSolidCone(radius,height,slices,stacks)
    CALL glPopMatrix() 
endif

call glPopAttrib()
    
endsubroutine

