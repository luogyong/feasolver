
!/*
! * Project an x,y pair onto a sphere of radius r OR a hyperbolic sheet
! * if we are away from the center of the sphere.
! */
real(8) function tb_project_to_sphere( r,  x,  y)
implicit none
    real(8),intent(in)::r,x,y
    real(8) d, t, z


    d = sqrt(x*x + y*y);
    if (d < r * 0.70710678118654752440) then     
        !/* Inside sphere */
        z = sqrt(r*r - d*d);
    else 
        !/* On hyperbola */
        t = r / 1.41421356237309504880;
        z = t*t / d;
    endif
    tb_project_to_sphere=z;
endfunction

!/*
! * Ok, simulate a track-ball.  Project the points onto the virtual
! * trackball, then figure out the axis of rotation, which is the cross
! * product of P1 P2 and O P1 (O is the center of the ball, 0,0,0)
! * Note:  This is a deformed trackball-- is a trackball in the center,
! * but is deformed into a hyperbolic sheet of rotation away from the
! * center.  This particular function was chosen after trying out
! * several variations.
! *
! * It is assumed that the arguments to this routine are in the range
! * (-1.0 ... 1.0)
! */
subroutine trackball(p1x, p1y, p2x, p2y,a,phi)
use opengl_glut
implicit none
    real(8),intent(in)::p1x, p1y, p2x, p2y
    real(8),intent(out):: a(3); !/* Axis of rotation */
    real(8),intent(out):: phi;  !/* how much to rotate about axis */
    
    real(8) p1(3), p2(3), d(3),x1,y1,x2,y2,ww,wh
    real(8) t,TRACKBALLSIZE,tb_project_to_sphere
    external tb_project_to_sphere
    
    TRACKBALLSIZE=0.8d0
    
    !a=0.d0
    !a(3)=1.0
    if ((p1x == p2x) .and. (p1y == p2y)) then
        !/* Zero rotation */
        phi=0.d0
        return;
    endif

    !/*
    ! * First, figure out z-coordinates for projection of P1 and P2 to
    ! * deformed sphere
    ! */
    ww=real(glutget(GLUT_WINDOW_WIDTH),8);wh=real(glutget(GLUT_WINDOW_HEIGHT),8)
    x1=(2.0*p1x-ww)/ww;y1=(wh-2.0*p1y)/wh;
    x2=(2.0*p2x-ww)/ww;y2=(wh-2.0*p2y)/wh;

    p1=[x1,y1,tb_project_to_sphere(TRACKBALLSIZE,x1,y1)]
    p2=[x2,y2,tb_project_to_sphere(TRACKBALLSIZE,x2,y2)]
    !vset(p1,p1x,p1y,tb_project_to_sphere(TRACKBALLSIZE,p1x,p1y));
    !vset(p2,p2x,p2y,tb_project_to_sphere(TRACKBALLSIZE,p2x,p2y));

    !/*
    ! *  Now, we want the cross product of P1 and P2
    ! */
    !vcross(p2,p1,a);
    
    call r8vec_cross_3d ( p1, p2, a )
    a=a/norm2(a)

    !/*
    ! *  Figure out how much to rotate around that axis.
    ! */
    !vsub(p1,p2,d);
    t = norm2(p1-p2) / (2.0*TRACKBALLSIZE);

    !/*
    ! * Avoid problems with out-of-control values...
    ! */
    if (t > 1.0) t = 1.0;
    if (t < -1.0) t = -1.0;
    phi = 2.0 * asin(t);
    !print *,a, phi/3.1415926*180

    !axis_to_quat(a,phi,q);
endsubroutine
