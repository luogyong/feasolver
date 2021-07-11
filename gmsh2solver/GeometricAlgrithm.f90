module GeoMetricAlgorithm
    implicit none
    public::rayintbox,rayintcyl,NORMAL_TRIFACE,cylinder_point_inside_3d,box_contains_point_nd,coplane_points_hull_3d,Pi
    
    private
    
    real(8),parameter::Pi=3.141592653589793
    
    
    contains
    

subroutine coplane_points_hull_3d(node_num, node_xy, hull_num, hull)
    implicit none

    integer ( kind = 4 ),intent(in)::node_num

    real ( kind = 8 ),intent(in):: node_xy(3,node_num)
    integer ( kind = 4 ),intent(out)::hull(node_num)
    integer ( kind = 4 ),intent(out):: hull_num
    
    real(8)::diff1(3)
    integer::i,n1,ix1(3),ix2(2)    
    do i=1,3
        diff1(i)=maxval(node_xy(1,:))-minval(node_xy(1,:))
    enddo
    ix1=[1,2,3]
    n1=minloc(diff1,dim=1)
    ix2=pack(ix1,ix1/=n1)
    
    call points_hull_2d ( node_num, node_xy(ix2,:), hull_num, hull )
    
endsubroutine
    
    
subroutine points_hull_2d ( node_num, node_xy, hull_num, hull )

!*****************************************************************************80
!
!! POINTS_HULL_2D computes the convex hull of 2D points.
!
!  Discussion:
!
!    The work involved is N*log(H), where N is the number of points, and H is
!    the number of points that are on the hull.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Output, integer ( kind = 4 ) HULL_NUM, the number of nodes that lie on 
!    the convex hull.
!
!    Output, integer ( kind = 4 ) HULL(NODE_NUM).  Entries 1 through HULL_NUM 
!    contain the indices of the nodes that form the convex hull, in order.
!
  implicit none

  integer ( kind = 4 ) node_num

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle_max
  !real ( kind = 8 ) angle_rad_2d
  
  real ( kind = 8 ) di
  real ( kind = 8 ) dr
  integer ( kind = 4 ) first
  integer ( kind = 4 ) hull(node_num)
  integer ( kind = 4 ) hull_num
  integer ( kind = 4 ) i
  real ( kind = 8 ) node_xy(2,node_num)
  real ( kind = 8 ) p_xy(2)
  integer ( kind = 4 ) q
  real ( kind = 8 ) q_xy(2)
  integer ( kind = 4 ) r
  real ( kind = 8 ) r_xy(2)

  if ( node_num < 1 ) then
    hull_num = 0
    return
  end if
!
!  If NODE_NUM = 1, the hull is the point.
!
  if ( node_num == 1 ) then
    hull_num = 1
    hull(1) = 1
    return
  end if
!
!  If NODE_NUM = 2, then the convex hull is either the two distinct points,
!  or possibly a single (repeated) point.
!
  if ( node_num == 2 ) then

    if ( node_xy(1,1) /= node_xy(1,2) .or. node_xy(2,1) /= node_xy(2,2) ) then
      hull_num = 2
      hull(1) = 1
      hull(2) = 2
    else
      hull_num = 1
      hull(1) = 1
    end if

    return

  end if
!
!  Find the leftmost point and call it "Q".
!  In case of ties, take the bottom-most.
!
  q = 1
  do i = 2, node_num
    if ( node_xy(1,i) < node_xy(1,q) .or. &
       ( node_xy(1,i) == node_xy(1,q) .and. node_xy(2,i) < node_xy(2,q) ) ) then
      q = i
    end if
  end do

  q_xy(1:2) = node_xy(1:2,q)
!
!  Remember the starting point, so we know when to stop!
!
  first = q
  hull_num = 1
  hull(1) = q
!
!  For the first point, make a dummy previous point, 1 unit south,
!  and call it "P".
!
  p_xy(1) = q_xy(1)
  p_xy(2) = q_xy(2) - 1.0D+00
!
!  Now, having old point P, and current point Q, find the new point R
!  so the angle PQR is maximal.
!
!  Watch out for the possibility that the two nodes are identical.
!
  do

    r = 0
    angle_max = 0.0D+00

    do i = 1, node_num

      if ( i /= q .and. &
           ( node_xy(1,i) /= q_xy(1) .or. node_xy(2,i) /= q_xy(2) ) ) then

        angle = angle_rad_2d ( p_xy, q_xy, node_xy(1:2,i) )

        if ( r == 0 .or. angle_max < angle ) then

          r = i
          r_xy(1:2) = node_xy(1:2,r)
          angle_max = angle
!
!  In case of ties, choose the nearer point.
!
        else if ( r /= 0 .and. angle == angle_max ) then

          di = ( node_xy(1,i) - q_xy(1) )**2 + ( node_xy(2,i) - q_xy(2) )**2
          dr = ( r_xy(1)      - q_xy(1) )**2 + ( r_xy(2)      - q_xy(2) )**2

          if ( di < dr ) then
            r = i
            r_xy(1:2) = node_xy(1:2,r)
            angle_max = angle
          end if

        end if

      end if

    end do
!
!  We are done when we have returned to the first point on the convex hull.
!
    if ( r == first ) then
      exit
    end if

    hull_num = hull_num + 1

    if ( node_num < hull_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POINTS_HULL_2D - Fatal error!'
      write ( *, '(a)' ) '  The algorithm has failed.'
      stop 1
    end if
!
!  Add point R to convex hull.
!
    hull(hull_num) = r
!
!  Set P := Q, Q := R, and prepare to search for next point R.
!
    q = r

    p_xy(1:2) = q_xy(1:2)
    q_xy(1:2) = r_xy(1:2)

  end do

  return
  
contains

function angle_rad_2d ( p1, p2, p3 )

!*****************************************************************************80
!
!! ANGLE_RAD_2D returns the angle in radians swept out between two rays in 2D.
!
!  Discussion:
!
!    Except for the zero angle case, it should be true that
!
!      ANGLE_RAD_2D ( P1, P2, P3 ) + ANGLE_RAD_2D ( P3, P2, P1 ) = 2 * PI
!
!        P1
!        /
!       /    
!      /     
!     /  
!    P2--------->P3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(2), P2(2), P3(2), define the rays
!    P1 - P2 and P3 - P2 which define the angle.
!
!    Output, real ( kind = 8 ) ANGLE_RAD_2D, the angle swept out by the rays,
!    in radians.  0 <= ANGLE_RAD_2D < 2 * PI.  If either ray has zero
!    length, then ANGLE_RAD_2D is set to 0.
!
  implicit none

  real ( kind = 8 ) angle_rad_2d
  real ( kind = 8 ) p(2)
  real ( kind = 8 ) p1(2)
  real ( kind = 8 ) p2(2)
  real ( kind = 8 ) p3(2)
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00

  p(1) = ( p3(1) - p2(1) ) * ( p1(1) - p2(1) ) &
       + ( p3(2) - p2(2) ) * ( p1(2) - p2(2) )

  p(2) = ( p3(1) - p2(1) ) * ( p1(2) - p2(2) ) &
       - ( p3(2) - p2(2) ) * ( p1(1) - p2(1) )

  if ( all ( p(1:2) == 0.0D+00)  ) then
    angle_rad_2d = 0.0D+00
    return
  end if

  angle_rad_2d = atan2 ( p(2), p(1) )

  if ( angle_rad_2d < 0.0D+00 ) then
    angle_rad_2d = angle_rad_2d + 2.0D+00 * r8_pi
  end if

  return
end    

end subroutine    
    
    
    
    function box_contains_point_nd ( dim_num, p1, p2, p )

        !*****************************************************************************80
        !
        !! BOX_CONTAINS_POINT_ND determines if a point is inside a box in ND.
        !
        !  Discussion:
        !
        !    A box is a rectangle with sides aligned on coordinate
        !    axes.  It can be described by its low and high corners, P1 and P2
        !    as the set of points P satisfying:
        !
        !      P1(1:DIM_NUM) <= P(1:DIM_NUM) <= P2(1:DIM_NUM).
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license. 
        !
        !  Modified:
        !
        !    28 February 2005
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
        !
        !    Input, real ( kind = 8 ) P1(DIM_NUM), P2(DIM_NUM), the low and high 
        !    corners of the box.
        !
        !    Input, real ( kind = 8 ) P(DIM_NUM), the point to be checked.
        !
        !    Output, logical ( kind = 4 ) BOX_CONTAINS_POINT_ND, is TRUE if the point 
        !    is inside the box.
        !
            implicit none

            integer ( kind = 4 ) dim_num

            logical ( kind = 4 ) box_contains_point_nd
            integer ( kind = 4 ) i
            real ( kind = 8 ),intent(in):: p(dim_num)
            real ( kind = 8 ) p1(dim_num)
            real ( kind = 8 ) p2(dim_num)

            box_contains_point_nd = .false.

            do i = 1, dim_num
            if ( p(i) < p1(i) .or. p2(i) < p(i) ) then
                return
            end if
            end do

            box_contains_point_nd = .true.

            return
        end function
        
 
    
    
function cylinder_point_inside_3d ( p1, p2, r, p ) result(inside)

!*****************************************************************************80
!
!! CYLINDER_POINT_INSIDE_3D determines if a cylinder contains a point in 3D.
!
!  Discussion:
!
!    The surface and interior of a (right) (finite) cylinder in 3D is defined 
!    by an axis, which is the line segment from point P1 to P2, and a 
!    radius R.  The points contained in the volume include:
!    * points at a distance less than or equal to R from the line through P1
!      and P2, whose nearest point on the line through P1 and P2 is, in fact,
!      P1, P2, or any point between them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(3), P2(3), the first and last points
!    on the axis line of the cylinder.
!
!    Input, real ( kind = 8 ) R, the radius of the cylinder.
!
!    Input, real ( kind = 8 ) P(3), the point.
!
!    Output, logical ( kind = 4 ) INSIDE, is TRUE if the point is 
!    inside the cylinder.
!
    implicit none

    integer ( kind = 4 ), parameter :: dim_num = 3
    integer::n1,iax1(2)
    real ( kind = 8 ) axis(dim_num)
    real ( kind = 8 ) axis_length
    logical ( kind = 4 ) inside
    real ( kind = 8 ) off_axis_component
    real ( kind = 8 ) p(dim_num)
    real ( kind = 8 ) p_dot_axis
    real ( kind = 8 ) p_length
    real ( kind = 8 ) p1(dim_num)
    real ( kind = 8 ) p2(dim_num)
    real ( kind = 8 ) r
    !real ( kind = 8 ) r8vec_norm

    axis(1:dim_num) = p2(1:dim_num) - p1(1:dim_num)
    axis_length = norm2 ( axis )
    if ( axis_length == 0.0D+00 ) then
    inside = .false.
    return
    end if    
    
    !added by lgy,圆柱平行于坐标轴的简单情况
    if(count(abs(axis)<1.e-8)==2) then
        if(abs(axis(3))>0) then
            n1=3 !平行轴
            iax1=[1,2]
        elseif(abs(axis(2))>0) then
            n1=2
            iax1=[1,3]
        else
            n1=1
            iax1=[2,3]
        endif
        
        if(p(n1)<min(p1(n1),p2(n1)).or.p(n1)>max(p1(n1),p2(n1))) then
            inside=.false.
        else
            if(norm2(p1(iax1)-p(iax1))>r) then
                inside=.false.
            else
                inside=.true.
            endif
        endif
    endif

    axis(1:dim_num) = axis(1:dim_num) / axis_length

    p_dot_axis = dot_product ( p(1:dim_num) - p1(1:dim_num), axis )
    !
    !  If the point lies below or above the "caps" of the cylinder, we're done.
    !
    if ( p_dot_axis < 0.0D+00 .or. axis_length < p_dot_axis ) then

    inside = .false.
    !
    !  Otherwise, determine the distance from P to the axis.
    !
    else

    p_length = norm2 ( p(1:dim_num) - p1(1:dim_num) )

    off_axis_component = sqrt ( p_length**2 - p_dot_axis**2 )

    if ( off_axis_component <= r ) then
        inside = .true.
    else
        inside = .false.
    end if

    end if

    return
end function    
    
    
subroutine rayintbox(numdim,minB,maxB,origin,dir,isint,coord,Surfin)
!/* 
!Fast Ray-Box Intersection
!by Andrew Woo
!from "Graphics Gems", Academic Press, 1990
!modify by lgy
!*/
    implicit none
    integer,intent(in)::NUMDIM
    real(8),intent(in)::origin(NUMDIM),dir(NUMDIM),minB(NUMDIM),maxB(NUMDIM)
    logical,intent(out)::isint
    real(8),intent(out)::coord(NUMDIM)
    integer,intent(out)::Surfin !intpoint on which plane. 1=xmin,2=xmax,3=ymin...
    integer,parameter:: RIGHT=	2
    integer,parameter:: LEFT=	1
    integer,parameter:: MIDDLE=	0
    
    logical::inside

	integer quadrant(NUMDIM);
	integer i;
	integer whichPlane;
	real(8) maxT(NUMDIM),origin1(NUMDIM),dir1(numdim)
	real(8) candidatePlane(NUMDIM);

	!/* Find candidate planes; this loop can be avoided if
 !  	rays cast all from the eye(assume perpsective view) */
    inside = .TRUE.;
	do i=1, NUMDIM
		if(origin(i) < minB(i)) then
			quadrant(i) = LEFT;
			candidatePlane(i) = minB(i);
			inside = .FALSE.;
		else if (origin(i) > maxB(i)) then
			quadrant(i) = RIGHT;
			candidatePlane(i) = maxB(i);
			inside = .FALSE.;
		else	
			quadrant(i) = MIDDLE;
		endif
    end do
	!/* Ray origin inside bounding box,then move outside  */
	if(inside)	then
		origin1 = origin+dir*2*norm2(minB-maxB);
        dir1=-dir
                    
	    do i=1, NUMDIM
		    if(origin1(i) < minB(i)) then
			    quadrant(i) = LEFT;
			    candidatePlane(i) = minB(i);
			    inside = .FALSE.;
		    else if (origin1(i) > maxB(i)) then
			    quadrant(i) = RIGHT;
			    candidatePlane(i) = maxB(i);
			    inside = .FALSE.;
		    else	
			    quadrant(i) = MIDDLE;
		    endif
        end do

    else
        origin1=origin
        dir1=dir
	endif


	!/* Calculate T distances to candidate planes */
	do i = 1, NUMDIM
		if (quadrant(i) /= MIDDLE.and. dir1(i) /=0.) then
			maxT(i) = (candidatePlane(i)-origin1(i)) / dir1(i);
		else
			maxT(i) = -1.;
        endif
    enddo
    
	!/* Get largest of the maxT's for final choice of intersection */
	whichPlane = 1;
	do i = 2,  NUMDIM
		if (maxT(whichPlane) < maxT(i)) whichPlane = i
    enddo
    
	!/* Check final candidate actually inside box */
	if (maxT(whichPlane) < 0.)  then
        isint=.FALSE.;
        return
    endif
    
	do i = 1, NUMDIM
		if (whichPlane /= i) then
			coord(i) = origin1(i) + maxT(whichPlane) *dir1(i);
			if (coord(i) < minB(i) .or. coord(i) > maxB(i)) then                
                isint= (.FALSE.);
                return
            endif
		else 
			coord(i) = candidatePlane(i);
		endif
    enddo
    
    !surfi=1-6, xmin=1,xmax=2....
    Surfin=2*(whichPlane-1)+quadrant(whichPlane)
	isint= (.TRUE.);				!/* ray hits box */


end subroutine    
    
subroutine rayintcyl(raybase,raycos,base,axis,radius,isint,in,out,botplane,topplane,surfin,surfout)
    implicit none
!/*
! * ANSI C code from the article
! * "Intersecting a Ray with a Cylinder"
! * by Joseph M. Cychosz and Warren N. Waggenspack, Jr.,
! * (3ksnn64@ecn.purdue.edu, mewagg@mewnw.dnet.lsu.edu)
! * in "Graphics Gems IV", Academic Press, 1994
! */
!
!#include	"GraphicsGems.h"
!#include	<math.h>
!
!/* ---- intcyl - Intersect a ray with a cylinder. --------------------- */
!/*									*/
!/*									*/
!/*	Description:							*/
!/*	    Intcyl determines the intersection of a ray with a		*/
!/*	    cylinder.							*/
!/*									*/
!/*	On entry:							*/
!/*	    raybase = The base point of the intersecting ray.		*/
!/*	    raycos  = The direction cosines of the above ray. (unit)	*/
!/*	    base    = The base location of the cylinder.		*/
!/*	    axis    = The axis of symmetry for the cylinder.  (unit)	*/
!/*	    radius  = The radius of the cylinder.			*/
!/*									*/
!/*	On return:							*/
!/*	    in	    = The entering distance of the intersection.	*/
!/*	    out	    = The leaving  distance of the intersection.	*/
!/*									*/
!/*	Returns:  True if the ray intersects the cylinder.		*/
!/*									*/
!/*	Note:	  In and/or out may be negative indicating the		*/
!/*		  cylinder is located behind the origin of the ray.	*/
!/*									*/
!/* -------------------------------------------------------------------- */
	real(8),intent(in)::		raybase(3)	!/* Base of the intersection ray */
	real(8),intent(in)::		raycos(3);	!/* Direction cosines of the ray */
	real(8),intent(in)::		base(3);	!	/* Base of the cylinder		*/
	real(8),intent(in)::		axis(3);	!	/* Axis of the cylinder		*/
	real(8),intent(in)::		radius;		!/* Radius of the cylinder	*/
	real(8),optional,intent(in)::		botplane(4);		!/* Bottom end-cap plane		*//* Plane: ax + by + cz + d = 0	*/
	real(8),optional,intent(in)::		topplane(4);		!/* Top end-cap plane		*//* Plane: ax + by + cz + d = 0	*/
    integer,optional,intent(out)::		surfin;	 !/* Entering surface identifier	*/
	integer,optional,intent(out)::		surfout; !	/* Exiting  surface identifier	*/
    logical,intent(out)::       isint
	real(8),intent(out)::		in;		!/* Entering distance		*/
	real(8),intent(out)::		out;	!	/* Leaving distance		*/


	logical::		hit;		!/* True if ray intersects cyl	*/
	real(8)::		RC(3);		!/* Ray base to cylinder base	*/
	real(8)::		dis;		!/* Shortest distance between	*/
					!/*   the ray and the cylinder	*/
	real(8)	::	t, s;		!/* Distances along the ray	*/
	real(8)::		n(3), D(3), O(3);
	real(8)::		ln;
	real(8),parameter::	pinf = 1.d30;	!/* Positive infinity		*/
    integer,parameter::	SIDE=	0	!	/* Object surface		*/
    integer,parameter::BOT=	1		!/* Bottom end-cap surface	*/
    integer,parameter::TOP=	2		!/* Top	  end-cap surface	*/
    
        
	RC = raybase - base;

	n=NORMAL_TRIFACE(reshape([0.d0,0.d0,0.d0,raycos,axis],([3,3])),.FALSE.);
    ln=norm2(n)
	if  ( ln  == 0. ) then !	/* ray parallel to cyl	*/
	    dis	 = dot_product (RC,axis);
	    D	 = RC - dis*axis;

	    dis	 = norm2(D);
	    in	 = -pinf;
	    out =  pinf;
	    isint= (dis <= radius);		!/* true if ray is in cyl*/
	    
	else

	    n=n/ln !V3Normalize (&n);
	    dis    = abs (dot_product (RC,n));	!	/* shortest distance	*/
	    hit  = (dis <= radius);

	    if  (hit) then !				/* if ray hits cylinder */
	        O=NORMAL_TRIFACE (reshape([0.d0,0.d0,0.d0,RC,axis],([3,3])),.FALSE.);
	        t = - dot_product (O,n) / ln;
	        O=NORMAL_TRIFACE(reshape([0.d0,0.d0,0.d0,n,axis],([3,3]))); !V3Cross (&n,axis,&O);
            !O=O/norm2(O) !  V3Normalize (&O);
	        s = abs (sqrt(radius*radius - dis*dis) / dot_product (raycos,O));
	        in	 = t - s;			!/* entering distance	*/
	        out = t + s;			!/* exiting  distance	*/
	    endif

	    isint= (hit);
    
    endif
    
    if(isint.and.(present(botplane).or.present(topplane))) then
        call clipobj()    
    endif    
    
    contains
    
    subroutine clipobj()

    !/* ---- clipobj - Clip object with plane pair. ------------------------ */
    !/*									*/
    !/*									*/
    !/*	Description:							*/
    !/*	    Clipobj clips the supplied infinite object with two		*/
    !/*	    (a top and a bottom) bounding planes.			*/
    !/*									*/
    !/*	On entry:							*/
    !/*	    raybase = The base point of the intersecting ray.		*/
    !/*	    raycos  = The direction cosines of the above ray. (unit)	*/
    !/*	    bot	    = The normal and perpendicular distance of the	*/
    !/*		      bottom plane.					*/
    !/*	    top	    = The normal and perpendicular distance of the	*/
    !/*		      top plane.					*/
    !/*	    objin   = The entering distance of the intersection with	*/
    !/*		      the object.					*/
    !/*	    objout  = The exiting  distance of the intersection with	*/
    !/*		      the object.					*/
    !/*									*/
    !/*	On return:							*/
    !/*	    objin   = The entering distance of the intersection.	*/
    !/*	    objout  = The exiting  distance of the intersection.	*/
    !/*	    surfin  = The identifier for the entering surface.		*/
    !/*	    surfout = The identifier for the leaving surface.		*/
    !/*									*/
    !/*	Returns:  True if the ray intersects the bounded object.	*/
    !/*									*/
    !/* -------------------------------------------------------------------- */
	    !real(8),intent(in)::		raybase(3);	!/* Base of the intersection ray */
	    !real(8),intent(in)::		raycos(3);	 !/* Direction cosines of the ray */
	    !real(8),intent(in)::		botplane(4);		!/* Bottom end-cap plane		*//* Plane: ax + by + cz + d = 0	*/
	    !real(8),intent(in)::		topplane(4);		!/* Top end-cap plane		*//* Plane: ax + by + cz + d = 0	*/
	    !real(8),intent(inout)::		objin;		!/* Entering distance		*/
	    !real(8),intent(inout)::		objout;	!/* Exiting  distance		*/
	    !integer,intent(out)::		surfin;	 !/* Entering surface identifier	*/
	    !integer,intent(out)::		surfout; !	/* Exiting  surface identifier	*/
        !logical,intent(out)::       isint

	    real(8)	dc, dw, t;
	    !real(8)	in, out;		!/* Object  intersection dists.	*/
    


	    if(present(surfin)) surfin =  SIDE;
        if(present(surfout)) surfout = SIDE
	    !in  = objin;
	    !out = objout;

    !/*	Intersect the ray with the bottom end-cap plane.		*/
        
        if(present(botplane)) then
            dc = botplane(1)*raycos(1)  + botplane(2)*raycos(2)  + botplane(3)*raycos(3);
	        dw = botplane(1)*raybase(1) + botplane(2)*raybase(2) + botplane(3)*raybase(3) + botplane(4);

	        if  ( dc == 0.0 ) then !{		/* If parallel to bottom plane	*/
	            if	( dw > 0.d0 ) then
                     isint=.FALSE.;
                     return
                endif
	         else
	            t  = - dw / dc;
	            if	( dc >= 0.0 ) then !{			    /* If far plane	*/
		            if  ( t > in .and. t < out ) then
                         out = t; 
                         if(present(surfout)) surfout = BOT; 
                    endif
		            if  ( t < in  ) then
                        isint=.FALSE.;
                        return
                    endif
	      
                else 				    !/* If near plane	*/
		            if  ( t > in .and. t < out ) then 
                         in	= t; 
                         if(present(surfin)) surfin  = BOT; 
                    endif
		            if  ( t > out ) then
                        isint=.FALSE.;
                        return
                    endif
	            endif
             endif
         endif

    !/*	Intersect the ray with the top end-cap plane.			*/
        if(present(topplane)) then
	    
            dc = topplane(1)*raycos(1)  + topplane(2)*raycos(2)  + topplane(3)*raycos(3);
	        dw = topplane(1)*raybase(1) + topplane(2)*raybase(2) + topplane(3)*raybase(3) + topplane(4);

	        if  ( dc == 0.0 ) then !{		/* If parallel to top plane	*/
	            if	( dw > 0.d0 ) then
                     isint=.FALSE.;
                     return
                endif
	        else 
	            t  = - dw / dc;
	            if	( dc >= 0.0 ) then !{			    /* If far plane	*/
		            if  ( t > in .and. t < out ) then
                         out = t;
                         if(present(surfout)) surfout = TOP; 
                    endif
		            if  ( t < in  ) then
                        isint=.FALSE.;
                        return
                    endif
                else !{				    /* If near plane	*/
		            if  ( t > in .and. t < out ) then
                         in	= t; 
                         if(present(surfin)) surfin  = TOP;
                    endif
		            if  ( t > out ) then
                        isint=.FALSE.;
                        return
                    endif
	            endif!}
	        endif
        endif
	    !objin	= in;
	    !objout = out;
	    isint=(in < out);
        return
    end subroutine       
    end subroutine 
    
    function NORMAL_TRIFACE(V,isunify) result (Normal)

    !*****************************************************************************80
    !V, XY OF THE FACET.
    !CALCULATE THE NORMAL VECTOR OF A TRI-FACET. 
    !

    !
      implicit none

      REAL(8),INTENT(IN)::V(:,:) !3*3
      logical,optional::isunify
      real ( kind = 8 ) v1(SIZE(V,DIM=1))
      real ( kind = 8 ) v2(SIZE(V,DIM=1))
      real ( kind = 8 ) normal(SIZE(V,DIM=1))
      logical::isunify1=.true.
      
      if(.not.present(isunify)) then
          isunify1=.true.
      else
          isunify1=isunify
      endif
      
      V1=V(:,2)-V(:,1);V2=V(:,3)-V(:,1);

      normal(1) = v1(2) * v2(3) - v1(3) * v2(2)
      normal(2) = v1(3) * v2(1) - v1(1) * v2(3)
      normal(3) = v1(1) * v2(2) - v1(2) * v2(1)
      if(isunify1) normal=normal/norm2(normal)
      return
  
    end function     
    

    
end module