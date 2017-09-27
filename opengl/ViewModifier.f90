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
private
public :: view_modifier_init, reset_view,init_lookat, init_lookfrom,view_from,VIEW_ZP,&
          minx,miny,minz,maxx,maxy,maxz,model_radius,keyboardCB,info

integer(kind=glcint), parameter :: ZOOM = 1, PAN = 2, ROTATE = 3, SCALEX = 4, &
                      SCALEY = 5, SCALEZ = 6,POINTPICKUP=7
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



type, private :: cart2D ! 2D cartesian coordinates
   real(kind=gldouble) :: x, y
end type cart2D

type, private :: cart3D ! 3D cartesian coordinates
   real(kind=gldouble) :: x, y, z
end type cart3D

type, private :: sphere3D ! 3D spherical coordinates
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
                     arrow_key_func = PAN

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

contains

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

!          -----
subroutine mouse(button, state, x, y)
!          -----
integer(kind=glcint), intent(in out) :: button, state, x, y
!integer,external::POINTlOC,PTINTRIlOC
integer::iel
real(8)::Pt1(3)

! This gets called when a mouse button changes
 
  if (button == GLUT_LEFT_BUTTON) then
    moving_left = .false.
    if(state == GLUT_DOWN.and.glutGetModifiers() == GLUT_ACTIVE_CTRL) then
		call ProbeatPoint(x,y)
    else
        
        if(state == GLUT_DOWN.and.glutGetModifiers() == GLUT_ACTIVE_SHIFT) then
            call glutSetCursor(GLUT_CURSOR_LEFT_ARROW)
            moving_left = .true.
            begin_left = cart2D(x,y)
            left_button_func=PAN
        ELSEIF(state == GLUT_DOWN) THEN
            if(isPickforstreamline) then
                
                call PickPoint(x,y,Pt1,IEL)
                IF(iel==0) then
                    info.str='the picked location is out of zone.Please pick again.'C
                    info.color=red;info.qkey=.true.
                ELSE
                          
                    info.str='Pick one more or Press q to exit.'C
                    if(info.qkey==.false.) then
                        isPickforstreamline=.false. 
                        info.str=''
                        call glutSetCursor(GLUT_CURSOR_LEFT_ARROW)
                        return
                    endif
                    call streamline_integration(Pt1)
                    info.color=green;info.qkey=.true.                                        
                                      
                endif
                
                
            else
                call glutSetCursor(GLUT_CURSOR_LEFT_ARROW)
			    moving_left = .true.
			    begin_left = cart2D(x,y)
			    left_button_func=ROTATE
            endif
        ENDIF
        
    endif
  endif
  
  if (button == GLUT_MIDDLE_BUTTON ) then
    if(state == GLUT_DOWN) then
        call glutSetCursor(GLUT_CURSOR_LEFT_ARROW)
        moving_middle = .true.
        begin_middle = cart2D(x,y)
    else
        moving_middle = .false.
    endif
    
  endif
  

end subroutine mouse

!          ------
subroutine motion(x, y)
!          ------
integer(kind=glcint), intent(in out) :: x, y

! This gets called when the mouse moves

integer :: button_function
type(cart2D) :: begin
real(kind=gldouble) :: factor

! Determine and apply the button function

if (moving_left) then
   button_function = left_button_func
   begin = begin_left
else if(moving_middle) then
   button_function = middle_button_func
   begin = begin_middle
end if

select case(button_function)
case (ZOOM)
   if (y < begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(begin%y-y))
   else if (y > begin%y) then
      factor = 1.0_gldouble + .002_gldouble*(y-begin%y)
   else
      factor = 1.0_gldouble
   end if
   IF(ISPERSPECT) THEN
        shift%z = shift%z/factor
   ELSE
        xscale_factor = xscale_factor * factor
        yscale_factor = yscale_factor * factor
        zscale_factor = zscale_factor * factor
   ENDIF
   
   
case (PAN)
   shift%x = shift%x + .01*(x - begin%x)
   shift%y = shift%y - .01*(y - begin%y)
case (ROTATE)
   angle%x = angle%x + (x - begin%x)
   angle%y = angle%y + (y - begin%y)

   !call trackball(begin.x, begin.y, real(x,8), real(y,8),axis_rotate,phi_rotate)
   !print *, begin.x, begin.y,x,y
case (SCALEX)
   if (y < begin%y) then
      factor = 1.0_gldouble + .002_gldouble*(begin%y-y)
   else if (y > begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(y-begin%y))
   else
      factor = 1.0_gldouble
   end if
   xscale_factor = xscale_factor * factor
case (SCALEY)
   if (y < begin%y) then
      factor = 1.0_gldouble + .002_gldouble*(begin%y-y)
   else if (y > begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(y-begin%y))
   else
      factor = 1.0_gldouble
   end if
   yscale_factor = yscale_factor * factor
case (SCALEZ)
   if (y < begin%y) then
      factor = 1.0_gldouble + .002_gldouble*(begin%y-y)
   else if (y > begin%y) then
      factor = 1.0_gldouble/(1.0_gldouble + .002_gldouble*(y-begin%y))
   else
      factor = 1.0_gldouble
   end if
   zscale_factor = zscale_factor * factor
end select

! update private variables and redisplay

if (moving_left) then
   begin_left = cart2D(x,y)
else if(moving_middle) then
   begin_middle = cart2D(x,y)
endif

if (moving_left .or. moving_middle) then
   call glutPostRedisplay
endif

return
end subroutine motion

!          ------
subroutine arrows(key, x, y)
!          ------
integer(glcint), intent(in out) :: key, x, y

! This routine handles the arrow key operations

real(kind=gldouble) :: factor

!select case(arrow_key_func)
!case(ZOOM)
!   select case(key)
!   case(GLUT_KEY_DOWN)
!      factor = 1.0_gldouble + .02_gldouble
!   case(GLUT_KEY_UP)
!      factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
!   case default
!      factor = 1.0_gldouble
!   end select
!   shift%z = factor*shift%z
!case(PAN)
!   select case(key)
!   case(GLUT_KEY_LEFT)
!      shift%x = shift%x - .02
!   case(GLUT_KEY_RIGHT)
!      shift%x = shift%x + .02
!   case(GLUT_KEY_DOWN)
!      shift%y = shift%y - .02
!   case(GLUT_KEY_UP)
!      shift%y = shift%y + .02
!   end select
!case(ROTATE)
!   select case(key)
!   case(GLUT_KEY_LEFT)
!      angle%x = angle%x - 1.0_gldouble
!   case(GLUT_KEY_RIGHT)
!      angle%x = angle%x + 1.0_gldouble
!   case(GLUT_KEY_DOWN)
!      angle%y = angle%y + 1.0_gldouble
!   case(GLUT_KEY_UP)
!      angle%y = angle%y - 1.0_gldouble
!   end select
!case(SCALEX)
!   select case(key)
!   case(GLUT_KEY_DOWN)
!      factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
!   case(GLUT_KEY_UP)
!      factor = 1.0_gldouble + .02_gldouble
!   case default
!      factor = 1.0_gldouble
!   end select
!   xscale_factor = xscale_factor * factor
!case(SCALEY)
!   select case(key)
!   case(GLUT_KEY_DOWN)
!      factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
!   case(GLUT_KEY_UP)
!      factor = 1.0_gldouble + .02_gldouble
!   case default
!      factor = 1.0_gldouble
!   end select
!   yscale_factor = yscale_factor * factor
!case(SCALEZ)
!   select case(key)
!   case(GLUT_KEY_DOWN)
!      factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
!   case(GLUT_KEY_UP)
!      factor = 1.0_gldouble + .02_gldouble
!   case default
!      factor = 1.0_gldouble
!   end select
!   zscale_factor = zscale_factor * factor
!
!end select
   
call glutPostRedisplay

return
end subroutine arrows




!          ------------
subroutine menu_handler(value)
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
end subroutine menu_handler

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
call glutreshapefunc(myreshape)
call glutMouseFunc(mouse)
call glutMotionFunc(motion)
call glutSpecialFunc(arrows)

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
menuid = glutCreateMenu(menu_handler)
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

subroutine myreshape(w,h)

    implicit none
    integer(glcint)::w,h
    


	real(gldouble):: wL,wH,R,wR,clipAreaXLeft,clipAreaXRight,clipAreaYBottom,clipAreaYTop, &
                    R1,dis1,near,far,xc,yc

    !call reset_view
    
	if (h == 0) h = 1;
	call glViewport(0, 0, w, h);

	R = real(w) / real(h);
    r1=((minx-maxx)**2+(miny-maxy)**2+(minz-maxz)**2)**0.5/2.0
    dis1=((init_lookat.x-init_lookfrom.x)**2+(init_lookat.y-init_lookfrom.y)**2+(init_lookat.z-init_lookfrom.z)**2)**0.5 
    
    near=1.0;far=near+2*r1+dis1*10
    !write(*,'(2G13.6)') 'NEAR=',NEAR,'FAR=',FAR 
    
	call glMatrixMode(GL_PROJECTION);
	call glLoadIdentity();


    !glortho and gluperspective 的参数都是相对于eye坐标的。
    
	if(IsPerspect) then
        call gluPerspective(30.0_gldouble, R, near, far)
    else
        
	    !wL = maxx -minx;	
	    !wH = maxy -miny;
	    !if (wH <= 0) wH = 1.0;
	    !wR = wL / wH;
        !
	    !if (wR > R)	then
		    ! ! Projection clipping area
		    ! clipAreaXLeft = minx;
		    ! clipAreaXRight = maxx;
		    ! clipAreaYBottom = (miny + maxy) / 2 - wL / R / 2.0;
		    ! clipAreaYTop = (miny + maxy) / 2 + wL / R / 2.0;
	    !
	    !else 
		    ! clipAreaXLeft =( minx + maxx) / 2.0 - wH* R / 2.0;
		    ! clipAreaXRight = (minx + maxx) / 2.0 + wH* R / 2.0;
		    ! clipAreaYBottom = miny;
		    ! clipAreaYTop = maxy;
	    !endif
	    !call glOrtho(clipAreaXLeft-init_lookfrom.x, clipAreaXRight-init_lookfrom.x, &
        !             clipAreaYBottom-init_lookfrom.y, clipAreaYTop-init_lookfrom.y,near,far)
        if(R>1) then
            call glOrtho(-r1*R, r1*R, -r1, r1,near,far) 
        else
            call glOrtho(-r1, r1, -r1/R, r1/R,near,far) 
        endif
    endif
    call glMatrixMode(GL_MODELVIEW);
	call glLoadIdentity();


end subroutine

subroutine keyboardCB(key,  x,  y)
    use strings
    integer,intent(in)::key,x,y
    integer::nlen=0,nsubstr=0
    character(len(info.inputstr))::substr(50)
    
    
    INPUTKEY=KEY
   
    
    select case(key)
    case(ichar('q'),ichar('Q'))
        if(INFO.QKEY) then
            info.str=''
            info.inputstr=''           
            INFO.ISNEEDINPUT=.FALSE.
            info.qkey=.false.
        endif
    case DEFAULT
 
    end select
    
    if(info.isneedinput) then
        nlen=len_trim(adjustl(info.inputstr))
        if(nlen>=len(info.inputstr)) info.inputstr=''
        if(key==8) then
            info.inputstr=info.inputstr(1:nlen-1)            
        elseif(key==13) then !enter
            call string_interpreter(info.inputstr,str2realArray)
            if(allocated(info.inputvar)) deallocate(info.inputvar)
            allocate(info.inputvar,source=str2vals)
            info.ninputvar=nstr2vals
            select case(info.func_id)
            case(FUNC_ID_GETSLICELOCATION)
                call getslicelocation()
            end select
            
            info.str=''
            info.inputstr=''
            INFO.NINPUTVAR=0
            if(allocated(info.inputvar)) deallocate(info.inputvar)
            INFO.ISNEEDINPUT=.false.
        else
            info.inputstr=trim(adjustl(info.inputstr))//char(key)
        endif
    ELSE
        INFO.INPUTSTR=''
    endif 
    
    call glutPostRedisplay
    !switch(key)
    !{
    !case 27: // ESCAPE
    !    exit(0);
    !    break;
    !
    !case 'd': // switch rendering modes (fill -> wire -> point)
    !case 'D':
    !    drawMode = ++drawMode % 3;
    !    if(drawMode == 0)        // fill mode
    !    {
    !        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    !        glEnable(GL_DEPTH_TEST);
    !        glEnable(GL_CULL_FACE);
    !    }
    !    else if(drawMode == 1)  // wireframe mode
    !    {
    !        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    !        glDisable(GL_DEPTH_TEST);
    !        glDisable(GL_CULL_FACE);
    !    }
    !    else                    // point mode
    !    {
    !        glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
    !        glDisable(GL_DEPTH_TEST);
    !        glDisable(GL_CULL_FACE);
    !    }
    !    break;
    !
    !case 'r':
    !case 'R':
    !    // reset rotation
    !    quat.set(1, 0, 0, 0);
    !    break;
    !
    !case ' ':
    !    if(trackball.getMode() == Trackball::ARC)
    !        trackball.setMode(Trackball::PROJECT);
    !    else
    !        trackball.setMode(Trackball::ARC);
    !    break;
    !
    !default:
    !    ;
    !}
endsubroutine


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
