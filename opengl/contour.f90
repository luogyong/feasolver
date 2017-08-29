subroutine DrawContour()
use solverds
use opengl_gl
use function_plotter
use MESHGEO
use view_modifier
implicit none
    
    
real(GLFLOAT),allocatable::vcolor(:,:)

integer :: i,j,k,n1,nfrac
real(GLDOUBLE) :: frac,graphmin,graphmax,graphstep
real(GLDOUBLE),allocatable :: contour_value(:)
character(128)::color_bar_title,str1,str2
! prepare to make a new display list



call glDeleteLists(contourlist, 1_glsizei)

if (draw_surface_solid) then
    call reset_view
    call glNewList(contourlist, gl_compile_and_execute)

! draw the solid surface



    call glPolygonMode(gl_front_and_back, gl_fill)
    call glBegin(gl_triangles)

   
    if (surface_color == rainbow_surface) then
        if(allocated(vcolor)) deallocate(vcolor)
        allocate(vcolor(4,NNUM))
        n1=outvar(CONTOUR_PLOT_VARIABLE).ivo
        minv=minval(nodalq(:,n1))
        maxv=maxval(nodalq(:,n1))
        WRITE(STR1,'(G13.3)') MAXV;WRITE(STR2,'(G13.3)') MINV;
        
        color_bar_title=trim(adjustl(outvar(CONTOUR_PLOT_VARIABLE).name))//',MAX='//TRIM(ADJUSTL(STR1))//',MIN='//TRIM(ADJUSTL(STR2))
        IF(ABS(MINV)<1E-7) MINV=SIGN(1E-7,MINV)
        call loose_label(minv, maxv,init_num_contour,graphmin,graphmax,graphstep,nfrac)
        graphstep=graphstep/scale_contour_num
        num_contour=MAX(floor((graphmax-graphmin)/graphstep)+1,2)
        if(allocated(contour_value)) deallocate(contour_value)
        allocate(contour_value(num_contour))
        do i=1,num_contour
            contour_value(i)=graphmin+(i-1)*graphstep
        enddo
        
        do i=1,nnum
            call get_rainbow(nodalq(i,n1),graphmin,graphmax,vcolor(:,i))
            if(isTransparency) vcolor(4,i)=0.6
        enddo
 
    
    endif

    do i=1,nface
        if(face(i).shape/=3) cycle
        do j=1,3
            !call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,vcolor(:,face(i).v(j)))
            IF(surface_color == white_surface)THEN
                call glcolor4fv(GRAY)
            ELSE
                call glcolor4fv(vcolor(:,face(i).v(j)))
            ENDIF
            call glvertex3dv(node(face(i).v(j)).coord)            
        enddo    
    enddo


   call glEnd

   call glBegin(GL_QUADS)

   
    !do i=1,enum
    !    IF(ESET(ELEMENT(I).SET).COUPLESET<ELEMENT(I).SET) CYCLE
    !    do j=1,3            
    !        call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,vcolor(:,element(i).node(j)))
    !        call glvertex3dv(node(element(i).node(j)).coord)
    !    enddo
    !enddo

    do i=1,nface
        if(face(i).shape/=4) cycle
        do j=1,4
            !call glMaterialfv(gl_front_and_back,gl_ambient_and_diffuse,vcolor(:,face(i).v(j)))
            IF(surface_color == white_surface)THEN
                call glcolor4fv(GRAY)
            ELSE
                call glcolor4fv(vcolor(:,face(i).v(j)))
            ENDIF
            call glvertex3dv(node(face(i).v(j)).coord)            
        enddo    
    enddo


   call glEnd   
   
   
   if (surface_color == rainbow_surface) call Color_Bar(contour_value,num_contour,NFRAC,color_bar_title)
endif ! draw_surface_solid

call glEndList

call glutPostRedisplay

if(allocated(vcolor)) deallocate(vcolor)
if(allocated(contour_value)) deallocate(contour_value)

endsubroutine    
    
    
    
!---+----3----+----2----+----1----+---<>---+----1----+----2----+----3----+----4
!--------------------------------  Color_Bar  ---------------------------------
 
subroutine Color_Bar(contour_value,num_contour,nfrac,BAR_TITLE)
    use opengl_gl    
    implicit none
    integer,intent(in)::num_contour,NFRAC
    real(GLDOUBLE),intent(in)::contour_value(num_contour)
    character(128),intent(in)::BAR_TITLE
    integer i;
    real(GLFLOAT):: color(4,num_contour)
    real(GLDOUBLE)::bot(2,num_contour),top(2,num_contour),SPACING
    CHARACTER(16)::STR1    
    
    
    SPACING=100/(NUM_CONTOUR-1)
    
    BOT(2,:)=0;TOP(2,:)=5;
    DO I=1,NUM_CONTOUR
        BOT(1,I)=0+(I-1)*SPACING
        TOP(1,I)=BOT(1,I)
        CALL get_rainbow(contour_value(I),contour_value(1),contour_value(NUM_CONTOUR),color(:,I))
    ENDDO
    
    
           
    !! Set up coordinate system to position color bar near bottom of window.
 
    CALL glMatrixMode(GL_PROJECTION);
    CALL glPushMatrix()
    CALL glLoadIdentity();
    CALL glOrtho(-40._GLDOUBLE, 140._GLDOUBLE, -20._GLDOUBLE, 200._GLDOUBLE, -1._GLDOUBLE, 1._GLDOUBLE);
    CALL glMatrixMode(GL_MODELVIEW);
    CALL glPushMatrix()
    CALL glLoadIdentity();
    
    call glPushAttrib(GL_ALL_ATTRIB_BITS)
	!call glPushAttrib(GL_LIGHTING_BIT .or. GL_CURRENT_BIT); ! lighting and color mask
	call glDisable(GL_LIGHTING);     ! need to disable lighting for proper text color
	call glDisable(GL_TEXTURE_2D);
    call glDisable(GL_CULL_FACE);
	call glDepthFunc(GL_ALWAYS);    
    
    ! Use Quad strips to make color bar.

    
    CALL glBegin(GL_QUAD_STRIP)
        do i = 1,NUM_CONTOUR
            CALL glColor3fv(color(:,i));
            CALL glVertex2Dv(bot(:,i));
            CALL glVertex2Dv(top(:,i));
        enddo
    CALL glEnd()
 
    ! Label ends of color bar.
 
    CALL glColor3f(0._GLFLOAT, 0._GLFLOAT, 0._GLFLOAT)
    
    call glPolygonMode(gl_front_and_back, gl_line)
    
    CALL glBegin(GL_QUAD_STRIP)
    do i = 1,NUM_CONTOUR
        !CALL glColor3fv(color(:,i));
        CALL glVertex2Dv(bot(:,i));
        CALL glVertex2Dv(top(:,i));
    enddo
    CALL glEnd()
    
    call glPolygonMode(gl_front_and_back, gl_FILL)
 
    !CALL bitmap_output(-5, 7, 0, "Min_H", GLUT_BITMAP_9_BY_15);
    !CALL bitmap_output(95, 7, 0, "Max_H", GLUT_BITMAP_9_BY_15);
    DO I=1,NUM_CONTOUR
        WRITE(STR1,10) CONTOUR_VALUE(I)
        call output(REAL(TOP(1,I),GLFLOAT), -3._glfloat,TRIM(ADJUSTL(STR1)))
    ENDDO
    
    call output(0._GLFLOAT, 7._glfloat,TRIM(ADJUSTL(BAR_TITLE)))
    
    CALL glPopMatrix()
    CALL glMatrixMode(GL_PROJECTION)
    CALL glPopMatrix()
    
 !  	CALL glEnable(GL_TEXTURE_2D);
	!CALL glEnable(GL_LIGHTING);
 !   call glEnable(GL_CULL_FACE);
	!CALL glDepthFunc(GL_LEQUAL);
	CALL glPopAttrib();
    
10  FORMAT(F16.<nfrac>)    
end subroutine

subroutine get_rainbow(val,minval,maxval,c)
use opengl_gl
implicit none
real(GLDOUBLE), intent(in) :: val,maxval,minval
real(GLFLOAT), intent(out) :: c(4)

real(GLFLOAT) :: f

if (maxval > minval) then
   f = (val-minval)/(maxval-minval)
else ! probably maxval==minval
   f = 0.5_glfloat
endif

if (f < .25) then
   c(1) = 0.0_glfloat
   c(2) = 4.0_glfloat * f
   c(3) = 1.0_glfloat
   c(4) = 1.0_glfloat
elseif (f < .5) then
   c(1) = 0.0_glfloat
   c(2) = 1.0_glfloat
   c(3) = 2.0_glfloat - 4.0_glfloat*f
   c(4) = 1.0_glfloat
elseif (f < .75) then
   c(1) = 4.0_glfloat * f - 2.0_glfloat
   c(2) = 1.0_glfloat
   c(3) = 0.0_glfloat
   c(4) = 1.0_glfloat
else
   c(1) = 1.0_glfloat
   c(2) = 4.0_glfloat - 4.0_glfloat*f
   c(3) = 0.0_glfloat
   c(4) = 1.0_glfloat
endif

end subroutine get_rainbow


subroutine loose_label(minv, maxv,ntick,graphmin,graphmax,graphstep,nfrac)
    implicit none
    integer,intent(in)::ntick
    real(8),intent(in)::minv,maxv
    integer,intent(out)::nfrac
    real(8),intent(out)::graphmin,graphmax,graphstep !/* graph range min and max and spacing */

    character::str(6), temp(20);
    !real(8) graphstep;				!/* tick mark spacing */
    !real(8) graphmin, graphmax;		
    real(8) range, x;
    
    real(8):: nicenum
    external nicenum
    
    !/* we expect min!=max */
    range = nicenum(maxv-minv, 0);
    graphstep = nicenum(range/(NTICK-1), 1);
    graphmin = floor(minv/graphstep)*graphstep;
    graphmax = ceilING(maxv/graphstep)*graphstep;
    !NFRAC=MAX(1,2)
    nfrac = MAX(-floor(log10(graphstep)),0)	!/* # of fractional digits to show */
    
    !sprintf(str, "%%.%df", nfrac);	!/* simplest axis labels */
    !WRITE(STR,"F3.1") NFRAC
    !printf("graphmin=%g graphmax=%g increment=%g\n", graphmin, graphmax, graphstep);

end subroutine

!/*
! * nicenum: find a "nice" number approximately equal to x.
! * Round the number if round=1, take ceiling if round=0
! */

function nicenum(x, round)
implicit none
real(8) nicenum,x;
integer round;

    integer expv;				!/* exponent of x */
    real(8) f;				!/* fractional part of x */
    real(8) nf;				!/* nice, rounded fraction */

    expv = floor(log10(x));
    f = x/10.**expv;		!/* between 1 and 10 */
    if (round) THEN
	    if (f<1.5) THEN
            nf = 1.;
	    else if (f<3.) THEN
            nf = 2.;
	    else if (f<7.) THEN
            nf = 5.;
	    else 
            nf = 10.;
        ENDIF
    else
	    if (f<=1.) THEN
            nf = 1.;
	    else if (f<=2.) THEN
            nf = 2.;
	    else if (f<=5.) THEN
            nf = 5.;
	    else 
            nf = 10.;
        ENDIF
    ENDIF
    
    NICENUM= nf*10.**expv;
end function