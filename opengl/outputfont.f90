subroutine output(x,y,s)
    use opengl_gl
    use opengl_glut
    implicit none
    
  real(glfloat) :: x,y
  character :: s*(*)
  character :: c
  integer :: i,lenc
  
  call glrasterpos2f(x,y)
  lenc = len(s)
  do i=1,lenc
     c = s(i:i)
     call glutbitmapcharacter(GLUT_BITMAP_HELVETICA_10, &
          ichar(c))
  end do
end subroutine output


subroutine drawStrokeText(angle,x,y,z,scale,s)
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
end subroutine drawStrokeText

