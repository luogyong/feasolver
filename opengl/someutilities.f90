function GetOGLPos2(x, y) result(object)
    use opengl_gl
    use opengl_glu

    implicit none
    integer(glint),intent(in)::x,y
    real(gldouble)::object(3)
    integer(GLint),dimension(4):: viewport;
    real(GLdouble),dimension(16):: modelview;
    real(GLdouble),dimension(16):: projection;
    real(gldouble) winX, winY, winZ;
    real(GLdouble) posX, posY, posZ;
    integer(GLint)::i
 
    call glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
    call glGetDoublev( GL_PROJECTION_MATRIX, projection );
    call glGetIntegerv( GL_VIEWPORT, viewport );
 
    winX = real(x);
    !winY = real(Y);
    winY = real(viewport(4)) - real(y);
    call f9y1glreadpixels( x, int(winY), 1, 1 , GL_DEPTH_COMPONENT, GL_double, winZ );
 
    i= gluUnProject( winX, winY, winZ, modelview, projection, viewport, posX, posY, posZ);
    object=[posX,posY,posZ]
    return
    !return CVector3(posX, posY, posZ);
end function