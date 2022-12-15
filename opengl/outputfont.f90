subroutine output(x,y,s)
    use opengl_gl
    use opengl_glut
    implicit none
    
  real(glfloat) :: x,y
  character :: s*(*)
  character :: c
  integer :: i,lenc
  
  call glrasterpos2f(x,y)
  lenc = len(trim(adjustl(s)))
  do i=1,lenc
     c = s(i:i)
     call glutbitmapcharacter(GLUT_BITMAP_HELVETICA_10, &
          ichar(c))
  end do
  
  call glutPostRedisplay
end subroutine output

subroutine output3D(x,y,z,s)
    use opengl_gl
    use opengl_glut
    implicit none
    
  real(gldouble),intent(in) :: x,y,z
  character ,intent(in):: s*(*)
  character :: c
  integer :: i,lenc
  
  call glrasterpos3d(x,y,z)
  lenc = len(s)
  do i=1,lenc
     c = s(i:i)
     call glutbitmapcharacter(GLUT_BITMAP_HELVETICA_10, &
          ichar(c))
  end do
  
  call glutPostRedisplay
end subroutine output3D


subroutine drawStrokeText(angle,x,y,z,scale,s)
    use opengl_gl
    use opengl_glut
    use pos_io,only:posdata
    implicit none
    
  real(GLDOUBLE),INTENT(in) :: ANGLE,x,y,z,scale
  character,intent(in) :: s*(*)
  character :: c
  integer :: i,lenc
  
  call glPushMatrix();
  
  call glTranslated(x, y,z);
  if(posdata.ndim==2) then
      call glrotated(ANGLE,0.,0.,1.)
  else
      call glrotated(ANGLE+180.,0.,1.,0.)
  endif
  
  call glScaled(scale, scale, scale);
  lenc = len(s)
  do i=1,lenc
     c = s(i:i)
     call glutStrokeCharacter(GLUT_STROKE_ROMAN, &
          ichar(c))
  end do
  call glPopMatrix();
  
  call glutPostRedisplay
end subroutine drawStrokeText

subroutine drawStrokeText2d(angle,x,y,z,scale,s)
    use opengl_gl
    use opengl_glut
    implicit none
    
  real(GLDOUBLE),INTENT(in) :: ANGLE,x,y,z,scale
  character,intent(in) :: s*(*)
  character :: c
  integer :: i,lenc
  
  call glPushMatrix();
  
  call glTranslated(x, y,z);

  call glrotated(ANGLE,0.,0.,1.)

  
  call glScaled(scale, scale, scale);
  lenc = len(s)
  do i=1,lenc
     c = s(i:i)
     call glutStrokeCharacter(GLUT_STROKE_ROMAN, &
          ichar(c))
  end do
  call glPopMatrix();
  
  call glutPostRedisplay
end subroutine drawStrokeText2d
    
    
SUBROUTINE showinfo(scolor)
    use opengl_gl
    use opengl_glut
	use function_plotter
    implicit none
    integer,intent(in)::scolor
    !character(*),intent(in) :: str(:)
    
    integer i;
    !real(GLFLOAT):: color(4,num_contour)
    real(GLDOUBLE)::orig(3),rowheight1
  
    integer(glCint),dimension(4)::viewport1

    
    !! Set up coordinate system to position color bar near bottom of window.
    call glgetintegerv(gl_viewport,viewport1)
    CALL glMatrixMode(GL_PROJECTION);
    CALL glPushMatrix()
    CALL glLoadIdentity();
    CALL glOrtho(0., real(viewport1(3),gldouble), 0., real(viewport1(4),gldouble), -1._GLDOUBLE, 1._GLDOUBLE);
    CALL glMatrixMode(GL_MODELVIEW);
    CALL glPushMatrix()
    CALL glLoadIdentity();
    
    call glPushAttrib(GL_ALL_ATTRIB_BITS)
	!call glPushAttrib(GL_LIGHTING_BIT .or. GL_CURRENT_BIT); ! lighting and color mask
	call glDisable(GL_LIGHTING);     ! need to disable lighting for proper text color
	call glDisable(GL_TEXTURE_2D);
    call glDisable(GL_CULL_FACE);
	call glDepthFunc(GL_ALWAYS);
    
    
    ORIG(1)=viewport1(3)/4.; ORIG(2)=0.95*viewport1(4); ORIG(3)=0.0;
    rowheight1=viewport1(4)/40.
    !ORIG=0.0

    CALL GLCOLOR4FV(MYCOLOR(:,SCOLOR))

    call output(REAL(ORIG(1),GLFLOAT),REAL(ORIG(2),GLFLOAT),TRIM(ADJUSTL(INFO.STR)))
    call output(REAL(ORIG(1),GLFLOAT),REAL(ORIG(2)-rowheight1,GLFLOAT),TRIM(ADJUSTL(INFO.INPUTSTR)))
    


    
    CALL glPopMatrix()
    CALL glMatrixMode(GL_PROJECTION)
    CALL glPopMatrix()
    
 !  	CALL glEnable(GL_TEXTURE_2D);
	!CALL glEnable(GL_LIGHTING);
 !   call glEnable(GL_CULL_FACE);
	!CALL glDepthFunc(GL_LEQUAL);
    !call glPolygonMode(gl_front_and_back, gl_FILL)
    CALL glPopAttrib();
    
    call glutPostRedisplay
    
10  FORMAT(G16.4)   
    

END SUBROUTINE

SUBROUTINE SHOW_MTEXT(STR,NSTR,POS,STRCOLOR,SHOWLIST)
!POS(widthfraction,heightfaction)
    use opengl_gl
    use opengl_glut
	use function_plotter
    implicit none
    INTEGER,INTENT(IN)::NSTR,STRCOLOR,SHOWLIST
	REAL(8),INTENT(IN)::POS(2)
    CHARACTER(*),INTENT(IN)::STR(NSTR)
    integer i;
    !real(GLFLOAT):: color(4,num_contour)
    real(GLDOUBLE)::SPACING,orig(3)
    CHARACTER(32)::STR1,STR2    
    integer(glCint),dimension(4)::viewport1

    call glDeleteLists(SHOWLIST, 1_glsizei)
    
    call glNewList(SHOWLIST, gl_compile_and_execute)   
    
    !! Set up coordinate system to position color bar near bottom of window.
    call glgetintegerv(gl_viewport,viewport1)
    CALL glMatrixMode(GL_PROJECTION);
    CALL glPushMatrix()
    CALL glLoadIdentity();
    CALL glOrtho(0., real(viewport1(3),gldouble), 0., real(viewport1(4),gldouble), -1._GLDOUBLE, 1._GLDOUBLE);
    CALL glMatrixMode(GL_MODELVIEW);
    CALL glPushMatrix()
    CALL glLoadIdentity();
    
    call glPushAttrib(GL_ALL_ATTRIB_BITS)
	!call glPushAttrib(GL_LIGHTING_BIT .or. GL_CURRENT_BIT); ! lighting and color mask
	call glDisable(GL_LIGHTING);     ! need to disable lighting for proper text color
	call glDisable(GL_TEXTURE_2D);
    call glDisable(GL_CULL_FACE);
	call glDepthFunc(GL_ALWAYS);
    
    
    ORIG(1)=viewport1(3)*pos(1); ORIG(2)=POS(2)*viewport1(4); ORIG(3)=0.0;
    !ORIG=0.0

    CALL GLCOLOR4FV(MYCOLOR(:,STRCOLOR))
    DO I=1,NSTR
        call output(REAL(ORIG(1),GLFLOAT),REAL(ORIG(2)-viewport1(4)/50.*I,GLFLOAT),TRIM(ADJUSTL(STR(I))))
    ENDDO
    


    
    CALL glPopMatrix()
    CALL glMatrixMode(GL_PROJECTION)
    CALL glPopMatrix()
    
 !  	CALL glEnable(GL_TEXTURE_2D);
	!CALL glEnable(GL_LIGHTING);
 !   call glEnable(GL_CULL_FACE);
	!CALL glDepthFunc(GL_LEQUAL);
    call glPolygonMode(gl_front_and_back, gl_FILL)
    CALL glPopAttrib();
    
call glendlist
    
	
    
call glutPostRedisplay
    
10  FORMAT(G16.4)   
    

END SUBROUTINE

subroutine string_interpreter(str,interpreter)
    use view_modifier
    use strings
    implicit none
    integer,intent(in)::interpreter
    character(*)::str
    character(len(str))::substr1(100),substr2(100)
    integer::i,j,nsubstr1,nsubstr2


    select case(interpreter)
    
    case default
        call parse(str,';',substr1,nsubstr1)
        if(allocated(str2vals)) deallocate(str2vals)
        allocate(str2vals(nsubstr1))
        nstr2vals=nsubstr1
        do i=1,nsubstr1            
            if(len_trim(substr1(i))>0) then
                call parse(substr1(i),'= ,',substr2,nsubstr2)
                if(nsubstr2>1) then
                    str2vals(i).name=trim(adjustl(substr2(1)))
                    str2vals(i).nval=nsubstr2-1
                    allocate(str2vals(i).val(str2vals(i).nval))
                    do j=1,str2vals(i).nval
                        read(substr2(j+1),*) str2vals(i).val(j)
                    enddo
                endif
            endif
        enddo
    end select        
endsubroutine
