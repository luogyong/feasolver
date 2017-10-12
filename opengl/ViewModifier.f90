! This program plots a function of two variables.  The function, called
! func_to_plot, is an external procedure at the end of the file.  
! This begins with the same module used in the modview example, followed by
! another module for plotting the function, called function_plotter.
! You might want to change default initial settings in modules modview and
! function_plotter.

! William F. Mitchell
! william.mitchell@nist.gov
! Mathematical and Computational Sciences Division
! National Institute of Standards and Technology
! August, 1999

!---------------------------------------------------------------------------

module view_modifier
    
! This module provides facilities to modify the view in an OpenGL window.
! The mouse buttons and keyboard arrow keys can be used to zoom, pan,
! rotate and change the scale.  A menu or submenu can be used to select which
! buttons perform which function and to reset the view to the initial settings.
! This is limited to one window.

! William F. Mitchell
! william.mitchell@nist.gov
! Mathematical and Computational Sciences Division
! National Institute of Standards and Technology
! April, 1998

! To use this module:
!
! 1) put a USE view_modifier statement in any program unit that calls a
!    procedure in this module
!
! 2) set the initial operation assignments, view and scale below the
!    "Initial configuration" comment below
!
! 3) call view_modifier_init after glutCreateWindow
!    This is a function that returns integer(kind=glcint) menuid.  The menuid
!    is the ID returned by glutCreateMenu.  You can either use the view_modifier
!    menu as your menu by calling glutAttachMenu immediately after
!    view_modifier_init, as in
!       menuid = view_modifier_init()
!       call glutAttachMenu(GLUT_RIGHT_BUTTON)
!    or by using the menuid to attach a submenu to your own menu, as in
!       call glutAddSubMenu("View Modifier",menuid)
!
! 4) in any callback functions that update the display, put
!       call reset_view
!    as the first executable statement
!
! Note that view_modifier_init sets the callback functions for glutMouseFunc,
! glutMotionFunc and glutSpecialFunc, so don't call these yourself
!
! The menu allows you to select what operation is attached to the left and
! middle mouse buttons and arrow keys, reset to the initial view, and quit.
! The right mouse button should be used for the menu.

use opengl_gl
use opengl_glu
use opengl_glut
use IndexColor
!use MESHGEO
implicit none
!private
!public :: view_modifier_init, reset_view,init_lookat, init_lookfrom,view_from,VIEW_ZP,&
!          minx,miny,minz,maxx,maxy,maxz,model_radius,keyboardCB,info

integer(kind=glcint), parameter :: ZOOM = 1, PAN = 2, ROTATE = 3, SCALEX = 4, &
                      SCALEY = 5, SCALEZ = 6,POINTPICKUP=7,LB_DRAWLINE=8,STEP_KEY=9
integer(kind=glcint), parameter :: RESET = 10,QUIT = 11,PRO_SYS=12
integer(kind=glcint), parameter:: VIEW_ZP =13,VIEW_ZN=14,VIEW_XP=15,VIEW_XN=16,VIEW_YP=17,    VIEW_YN=18,View_Save=19,View_CAST=20
real(kind=gldouble), parameter :: PI = 3.141592653589793_gldouble
real(kind=gldouble)::minx,miny,minz,maxx,maxy,maxz,axis_rotate(3),phi_rotate=0.,model_radius
real(gldouble),dimension(16)::ModelViewSaved
LOGICAL,public::ISPERSPECT=.TRUE.
INTEGER,PUBLIC::INPUTKEY=-1
LOGICAL,PUBLIC::isPickforslice=.FALSE.,isPickforstreamline=.false.

INTEGER,PUBLIC,PARAMETER::str2realArray=1
type string2val_tydef
    character(64)::Name
    integer::nval
    real(8),allocatable::val(:)
endtype
type(string2val_tydef),public,allocatable::str2vals(:)
integer,public::nstr2vals=0

type infoshow_tydef
    character(256)::str='',inputstr=''
    LOGICAL::ISNEEDINPUT=.FALSE.
    integer::color=1,interpreter=str2realArray,func_id=0,ninputvar=0
    logical::qkey=.TRUE.
    type(string2val_tydef),allocatable::inputvar(:)
    
endtype
type(infoshow_tydef)::info
integer,public,parameter::FUNC_ID_GETSLICELOCATION=1



type :: cart2D ! 2D cartesian coordinates
   real(kind=gldouble) :: x, y
end type cart2D

type :: cart3D ! 3D cartesian coordinates
   real(kind=gldouble) :: x, y, z
end type cart3D

type :: sphere3D ! 3D spherical coordinates
   real(kind=gldouble) :: theta, phi, rho
end type sphere3D

type(cart2D), save :: angle
type(cart3D), save :: shift
real(kind=gldouble),public, save :: xscale_factor, yscale_factor, zscale_factor
logical, save :: moving_left, moving_middle
type(cart2D), save :: begin_left, begin_middle

interface operator(+)
   module procedure cart3D_plus_cart3D
end interface
interface operator(-)
   module procedure cart3D_minus_cart3D
end interface

! ------- Initial configuration -------

! Set the initial operation performed by each button and the arrow keys.
! The operations are ZOOM, PAN, ROTATE, SCALEX, SCALEY, and SCALEZ

integer, save ::   left_button_func = ROTATE, &
                 middle_button_func = ZOOM, &
                     arrow_key_func = STEP_KEY

! Set the initial view as the point you are looking at, the point you are
! looking from, and the scale factors

type(cart3D) :: init_lookat,init_lookfrom
   !init_lookat = cart3D(0.5_gldouble, 0.5_gldouble, 0.0_gldouble), &
   !init_lookfrom = cart3D(5.0_gldouble, 10.0_gldouble, 2.5_gldouble)
   !init_lookat = cart3D(10.0_gldouble, 10.0_gldouble, 0.0_gldouble), &
   !init_lookfrom = cart3D(0.0_gldouble, 0.0_gldouble, 50_gldouble)

real(kind=gldouble), parameter :: &
   init_xscale_factor = 1.0_gldouble, &
   init_yscale_factor = 1.0_gldouble, &
   init_zscale_factor = 1.0_gldouble

! -------- end of Initial configuration ------
TYPE LINE_TEMP_TYDEF   
    REAL(8)::V1(3)=0.D0,V2(3)=0.D0
    LOGICAL::SHOW=.FALSE.
    INTEGER::ICOLOR=RED
    REAL(GLFLOAT)::LW=1.5
    CONTAINS
    PROCEDURE::DRAW=>DRAW_LINE
ENDTYPE
TYPE(LINE_TEMP_TYDEF),PUBLIC::LINE_TEMP



contains



SUBROUTINE DRAW_LINE(LINE)
    IMPLICIT NONE
    CLASS(LINE_TEMP_TYDEF),INTENT(in):: LINE
    
    call glPushAttrib(GL_ALL_ATTRIB_BITS)
    
    CALL GLLINEWIDTH(LINE.LW)
    call glcolor4fv(mycolor(:,LINE.ICOLOR))		
	CALL GLBEGIN(GL_LINES)
        call glvertex3dv(LINE.V1)
        call glvertex3dv(LINE.V2)
    CALL GLEND()
    CALL GLPOPATTRIB()
    
END SUBROUTINE

!          -------------
subroutine reset_to_init
!          -------------

! This resets the view to the initial configuration

type(sphere3D) :: slookfrom

slookfrom = cart2sphere(init_lookfrom-init_lookat)
angle%x = -180.0_gldouble*slookfrom%theta/PI 
angle%y = -180.0_gldouble*slookfrom%phi/PI
!axis_rotate=[0.,0.,1.];phi_rotate=0.d0
shift%x = 0.0_gldouble
shift%y = 0.0_gldouble
shift%z = -slookfrom%rho
xscale_factor = init_xscale_factor
yscale_factor = init_yscale_factor
zscale_factor = init_zscale_factor

call glutPostRedisplay

return
end subroutine reset_to_init

!          ---------------
subroutine view_from(VALUE)
!          ---------------

! This sets the view to be from straight above
INTEGER,INTENT(IN)::VALUE
real(gldouble)::t1
type(sphere3D) :: slookfrom

t1=2.*model_radius
SELECT CASE(VALUE)
CASE(VIEW_ZP)
    slookfrom = cart2sphere(cart3D(0.0,0.0,t1))
CASE(VIEW_ZN)
    slookfrom = cart2sphere(cart3D(0.0,0.0,-t1))
CASE(VIEW_XP)
    slookfrom = cart2sphere(cart3D(0.0,t1,0.0))
CASE(VIEW_XN)
    slookfrom = cart2sphere(cart3D(0.0,-t1,0.0)) 
CASE(VIEW_YP)
    slookfrom = cart2sphere(cart3D(t1,0.0,0.0))
CASE(VIEW_YN)
    slookfrom = cart2sphere(cart3D(-t1,0.0,0.0))      
ENDSELECT

angle%x = -180.0_gldouble*slookfrom%theta/PI
angle%y = -180.0_gldouble*slookfrom%phi/PI

call glutPostRedisplay

return
end subroutine view_from

!          ----------
subroutine reset_view
!          ----------
real(gldouble),dimension(16)::objectXform
! This routine resets the view to the current orientation and scale

call glMatrixMode(GL_MODELVIEW)
!call glPopMatrix
!call glPushMatrix
!call glGetDoublev( GL_MODELVIEW_MATRIX, objectXform );
call glLoadIdentity()

!call glRotated(phi_rotate/pi*180.d0, axis_rotate(1), axis_rotate(2), axis_rotate(3))
!call glMultMatrixd( reshape(objectXform,(/4,4/) ))
call glTranslated(shift%x, shift%y, shift%z)
call glRotated(angle%x, 0.0_gldouble, 0.0_gldouble, 1.0_gldouble)
call glRotated(angle%y, cos(PI*angle%x/180.0_gldouble), &
               -sin(PI*angle%x/180.0_gldouble), 0.0_gldouble)
call glTranslated(-init_lookat%x, -init_lookat%y, -init_lookat%z)
call glScaled(xscale_factor,yscale_factor,zscale_factor)
!call glPopMatrix
return
end subroutine reset_view

subroutine menu_handler1(value)
!          ------------
integer(kind=glcint), intent(in out) :: value

! This routine handles the first level entries in the menu

select case(value)

case(RESET)
   call reset_to_init
   
case(VIEW_ZP,VIEW_ZN,VIEW_XP,VIEW_XN,VIEW_YP,VIEW_YN)
   call view_from(VALUE)
   
!CASE(VIEW_SAVE)
!    call glGetDoublev(GL_MODELVIEW_MATRIX, ModelViewSaved)
!CASE(VIEW_CAST)
!    CALL GLMATRIXMODE(GL_MODELVIEW_MATRIX) 
!    call glLoadMatrixd(reshape(ModelViewSaved,(/4,4/)));
!    call glutPostRedisplay
case(QUIT)
   stop
case(PRO_SYS)
    ISPERSPECT=.NOT.ISPERSPECT
    CALL myreshape(GLUTGET(GLUT_WINDOW_WIDTH),GLUTGET(GLUT_WINDOW_HEIGHT))
    call glutPostRedisplay
end select


return
end subroutine menu_handler1

!          ---------------
subroutine set_left_button(value)
!          ---------------
integer(kind=glcint), intent(in out) :: value

! This routine sets the function of the left button as given by menu selection

left_button_func = value

return
end subroutine set_left_button

!          -----------------
subroutine set_middle_button(value)
!          -----------------
integer(kind=glcint), intent(in out) :: value

! This routine sets the function of the middle button as given by menu selection

middle_button_func = value

return
end subroutine set_middle_button

!          --------------
subroutine set_arrow_keys(value)
!          --------------
integer(kind=glcint), intent(in out) :: value

! This routine sets the function of the arrow keys as given by menu selection

arrow_key_func = value

return
end subroutine set_arrow_keys

!        ------------------
function view_modifier_init() result(menuid)
!        ------------------
integer(kind=glcint) :: menuid

! This initializes the view modifier variables and sets initial view.
! It should be called immediately after glutCreateWindow

integer(kind=glcint) :: button_left, button_middle, arrow_keys

! set the callback functions


! create the menu


!button_left = glutCreateMenu(set_left_button)
!call glutAddMenuEntry("rotate",ROTATE)
!call glutAddMenuEntry("zoom",ZOOM)
!call glutAddMenuEntry("pan",PAN)
!call glutAddMenuEntry("scale x",SCALEX)
!call glutAddMenuEntry("scale y",SCALEY)
!call glutAddMenuEntry("scale z", SCALEZ)
!button_middle = glutCreateMenu(set_middle_button)
!call glutAddMenuEntry("rotate",ROTATE)
!call glutAddMenuEntry("zoom",ZOOM)
!call glutAddMenuEntry("pan",PAN)
!call glutAddMenuEntry("scale x",SCALEX)
!call glutAddMenuEntry("scale y",SCALEY)
!call glutAddMenuEntry("scale z", SCALEZ)
!arrow_keys = glutCreateMenu(set_arrow_keys)
!call glutAddMenuEntry("rotate",ROTATE)
!call glutAddMenuEntry("zoom",ZOOM)
!call glutAddMenuEntry("pan",PAN)
!call glutAddMenuEntry("scale x",SCALEX)
!call glutAddMenuEntry("scale y",SCALEY)
!call glutAddMenuEntry("scale z", SCALEZ)
menuid = glutCreateMenu(menu_handler1)
!call glutAddSubMenu("left mouse button",button_left)
!call glutAddSubMenu("middle mouse button",button_middle)
!call glutAddSubMenu("arrow keys",arrow_keys)
call glutAddMenuEntry("Perspect/Ortho",PRO_SYS)
call glutAddMenuEntry("reset to initial view",RESET)
call glutAddMenuEntry("view from Z+",VIEW_ZP)
call glutAddMenuEntry("view from Z-",VIEW_ZN)
call glutAddMenuEntry("view from X+",VIEW_XP)
call glutAddMenuEntry("view from X-",VIEW_XN)
call glutAddMenuEntry("view from Y+",VIEW_YP)
call glutAddMenuEntry("view from Y-",VIEW_YN)
!call glutAddMenuEntry("SaveCurrentView",View_Save)
!call glutAddMenuEntry("UseSavedView",View_Cast)
!call glutAddMenuEntry("quit",QUIT)

! set the perspective

!call glMatrixMode(GL_PROJECTION)
!call gluPerspective(10.0_gldouble, 1.0_gldouble, 0.1_gldouble, 200.0_gldouble)

! set the initial view

!call glMatrixMode(GL_MODELVIEW)
!call glPushMatrix
call reset_to_init

return
end function view_modifier_init

!        -----------
function sphere2cart(spoint) result(cpoint)
!        -----------
type(sphere3D), intent(in) :: spoint
type(cart3D) :: cpoint

! This converts a 3D point from spherical to cartesean coordinates

real(kind=gldouble) :: t,p,r

t=spoint%theta
p=spoint%phi
r=spoint%rho

cpoint%x = r*cos(t)*sin(p)
cpoint%y = r*sin(t)*sin(p)
cpoint%z = r*cos(p)

return
end function sphere2cart

!        -----------
function cart2sphere(cpoint) result(spoint)
!        -----------
type(cart3D), intent(in) :: cpoint
type(sphere3D) :: spoint

! This converts a 3D point from cartesean to spherical coordinates

real(kind=gldouble) :: x,y,z

x=cpoint%x
y=cpoint%y
z=cpoint%z

spoint%rho = sqrt(x*x+y*y+z*z)
if (x==0.0_gldouble .and. y==0.0_gldouble) then
   spoint%theta = 0.0_gldouble
else
   spoint%theta = atan2(y,x)
end if
if (spoint%rho == 0.0_gldouble) then
   spoint%phi = 0.0_gldouble
else
   spoint%phi = acos(z/spoint%rho)
endif

return
end function cart2sphere

!        ------------------
function cart3D_plus_cart3D(cart1,cart2) result(cart3)
!        ------------------
type(cart3D), intent(in) :: cart1, cart2
type(cart3D) :: cart3

! Compute the sum of two 3D cartesean points

cart3%x = cart1%x + cart2%x
cart3%y = cart1%y + cart2%y
cart3%z = cart1%z + cart2%z

return
end function cart3D_plus_cart3D

!        -------------------
function cart3D_minus_cart3D(cart1,cart2) result(cart3)
!        -------------------
type(cart3D), intent(in) :: cart1, cart2
type(cart3D) :: cart3

! Compute the difference of two 3D cartesean points

cart3%x = cart1%x - cart2%x
cart3%y = cart1%y - cart2%y
cart3%z = cart1%z - cart2%z

return
end function cart3D_minus_cart3D



end module view_modifier


