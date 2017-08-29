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