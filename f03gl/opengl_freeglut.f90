MODULE opengl_glut
USE opengl_kinds
IMPLICIT NONE
PRIVATE


!  freeglut.h
!  
!  The freeglut library include file
!  
!  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
!  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
!  PAWEL W. OLSZTA BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
!  IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
!  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


!  freeglut_std.h
!  
!  The GLUT-compatible part of the freeglut library include file
!  
!  Copyright (c) 1999-2000 Pawel W. Olszta. All Rights Reserved.
!  Written by Pawel W. Olszta, <olszta@sourceforge.net>
!  Creation date: Thu Dec 2 1999
!  
!  Permission is hereby granted, free of charge, to any person obtaining a
!  copy of this software and associated documentation files (the "Software"),
!  to deal in the Software without restriction, including without limitation
!  the rights to use, copy, modify, merge, publish, distribute, sublicense,
!  and/or sell copies of the Software, and to permit persons to whom the
!  Software is furnished to do so, subject to the following conditions:
!  
!  The above copyright notice and this permission notice shall be included
!  in all copies or substantial portions of the Software.
!  
!  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
!  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
!  PAWEL W. OLSZTA BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
!  IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
!  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


!  Non-Windows definition of FGAPI and FGAPIENTRY  

!  The freeglut and GLUT API versions
INTEGER(GLenum), PARAMETER, PUBLIC :: FREEGLUT = 1
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_API_VERSION = 4
INTEGER(GLenum), PARAMETER, PUBLIC :: FREEGLUT_VERSION_2_0 = 1
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_XLIB_IMPLEMENTATION = 13

!  GLUT API macro definitions -- the special key codes:
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F1 = INT(z'0001')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F2 = INT(z'0002')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F3 = INT(z'0003')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F4 = INT(z'0004')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F5 = INT(z'0005')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F6 = INT(z'0006')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F7 = INT(z'0007')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F8 = INT(z'0008')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F9 = INT(z'0009')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F10 = INT(z'000A')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F11 = INT(z'000B')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_F12 = INT(z'000C')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_LEFT = INT(z'0064')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_UP = INT(z'0065')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_RIGHT = INT(z'0066')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_DOWN = INT(z'0067')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_PAGE_UP = INT(z'0068')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_PAGE_DOWN = INT(z'0069')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_HOME = INT(z'006A')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_END = INT(z'006B')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_INSERT = INT(z'006C')

!  GLUT API macro definitions -- mouse state definitions
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_LEFT_BUTTON = INT(z'0000')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_MIDDLE_BUTTON = INT(z'0001')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_RIGHT_BUTTON = INT(z'0002')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DOWN = INT(z'0000')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_UP = INT(z'0001')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_LEFT = INT(z'0000')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ENTERED = INT(z'0001')

!  GLUT API macro definitions -- the display mode definitions
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_RGB = INT(z'0000')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_RGBA = INT(z'0000')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INDEX = INT(z'0001')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_SINGLE = INT(z'0000')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DOUBLE = INT(z'0002')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ACCUM = INT(z'0004')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ALPHA = INT(z'0008')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DEPTH = INT(z'0010')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_STENCIL = INT(z'0020')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_MULTISAMPLE = INT(z'0080')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_STEREO = INT(z'0100')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_LUMINANCE = INT(z'0200')

!  GLUT API macro definitions -- windows and menu related definitions
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_MENU_NOT_IN_USE = INT(z'0000')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_MENU_IN_USE = INT(z'0001')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NOT_VISIBLE = INT(z'0000')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VISIBLE = INT(z'0001')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HIDDEN = INT(z'0000')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_FULLY_RETAINED = INT(z'0001')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_PARTIALLY_RETAINED = INT(z'0002')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_FULLY_COVERED = INT(z'0003')

!  GLUT API macro definitions -- fonts definitions
!  
!  Steve Baker suggested to make it binary compatible with GLUT:
!  /*
!  I don't really know if it's a good idea... But here it goes:

!  Those pointers will be used by following definitions:

!  GLUT API macro definitions -- the glutGet parameters
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_X = INT(z'0064')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_Y = INT(z'0065')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_WIDTH = INT(z'0066')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_HEIGHT = INT(z'0067')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_BUFFER_SIZE = INT(z'0068')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_STENCIL_SIZE = INT(z'0069')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_DEPTH_SIZE = INT(z'006A')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_RED_SIZE = INT(z'006B')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_GREEN_SIZE = INT(z'006C')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_BLUE_SIZE = INT(z'006D')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_ALPHA_SIZE = INT(z'006E')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_ACCUM_RED_SIZE = INT(z'006F')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_ACCUM_GREEN_SIZE = INT(z'0070')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_ACCUM_BLUE_SIZE = INT(z'0071')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_ACCUM_ALPHA_SIZE = INT(z'0072')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_DOUBLEBUFFER = INT(z'0073')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_RGBA = INT(z'0074')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_PARENT = INT(z'0075')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_NUM_CHILDREN = INT(z'0076')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_COLORMAP_SIZE = INT(z'0077')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_NUM_SAMPLES = INT(z'0078')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_STEREO = INT(z'0079')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_CURSOR = INT(z'007A')

INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_SCREEN_WIDTH = INT(z'00C8')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_SCREEN_HEIGHT = INT(z'00C9')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_SCREEN_WIDTH_MM = INT(z'00CA')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_SCREEN_HEIGHT_MM = INT(z'00CB')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_MENU_NUM_ITEMS = INT(z'012C')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DISPLAY_MODE_POSSIBLE = INT(z'0190')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_WINDOW_X = INT(z'01F4')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_WINDOW_Y = INT(z'01F5')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_WINDOW_WIDTH = INT(z'01F6')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_WINDOW_HEIGHT = INT(z'01F7')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_DISPLAY_MODE = INT(z'01F8')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ELAPSED_TIME = INT(z'02BC')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_FORMAT_ID = INT(z'007B')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_INIT_STATE = INT(z'007C')

!  GLUT API macro definitions -- the glutDeviceGet parameters
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_KEYBOARD = INT(z'0258')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_MOUSE = INT(z'0259')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_SPACEBALL = INT(z'025A')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_DIAL_AND_BUTTON_BOX = INT(z'025B')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_TABLET = INT(z'025C')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NUM_MOUSE_BUTTONS = INT(z'025D')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NUM_SPACEBALL_BUTTONS = INT(z'025E')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NUM_BUTTON_BOX_BUTTONS = INT(z'025F')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NUM_DIALS = INT(z'0260')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NUM_TABLET_BUTTONS = INT(z'0261')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DEVICE_IGNORE_KEY_REPEAT = INT(z'0262')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DEVICE_KEY_REPEAT = INT(z'0263')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_JOYSTICK = INT(z'0264')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_OWNS_JOYSTICK = INT(z'0265')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_BUTTONS = INT(z'0266')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_AXES = INT(z'0267')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_POLL_RATE = INT(z'0268')

!  GLUT API macro definitions -- the glutLayerGet parameters
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_OVERLAY_POSSIBLE = INT(z'0320')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_LAYER_IN_USE = INT(z'0321')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_HAS_OVERLAY = INT(z'0322')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_TRANSPARENT_INDEX = INT(z'0323')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NORMAL_DAMAGED = INT(z'0324')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_OVERLAY_DAMAGED = INT(z'0325')

!  GLUT API macro definitions -- the glutVideoResizeGet parameters
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_POSSIBLE = INT(z'0384')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_IN_USE = INT(z'0385')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_X_DELTA = INT(z'0386')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_Y_DELTA = INT(z'0387')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_WIDTH_DELTA = INT(z'0388')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_HEIGHT_DELTA = INT(z'0389')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_X = INT(z'038A')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_Y = INT(z'038B')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_WIDTH = INT(z'038C')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VIDEO_RESIZE_HEIGHT = INT(z'038D')

!  GLUT API macro definitions -- the glutUseLayer parameters
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_NORMAL = INT(z'0000')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_OVERLAY = INT(z'0001')

!  GLUT API macro definitions -- the glutGetModifiers parameters
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ACTIVE_SHIFT = INT(z'0001')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ACTIVE_CTRL = INT(z'0002')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ACTIVE_ALT = INT(z'0004')

!  GLUT API macro definitions -- the glutSetCursor parameters
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_RIGHT_ARROW = INT(z'0000')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_LEFT_ARROW = INT(z'0001')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_INFO = INT(z'0002')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_DESTROY = INT(z'0003')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_HELP = INT(z'0004')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_CYCLE = INT(z'0005')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_SPRAY = INT(z'0006')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_WAIT = INT(z'0007')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_TEXT = INT(z'0008')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_CROSSHAIR = INT(z'0009')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_UP_DOWN = INT(z'000A')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_LEFT_RIGHT = INT(z'000B')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_TOP_SIDE = INT(z'000C')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_BOTTOM_SIDE = INT(z'000D')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_LEFT_SIDE = INT(z'000E')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_RIGHT_SIDE = INT(z'000F')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_TOP_LEFT_CORNER = INT(z'0010')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_TOP_RIGHT_CORNER = INT(z'0011')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_BOTTOM_RIGHT_CORNER = INT(z'0012')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_BOTTOM_LEFT_CORNER = INT(z'0013')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_INHERIT = INT(z'0064')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_NONE = INT(z'0065')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CURSOR_FULL_CROSSHAIR = INT(z'0066')

!  GLUT API macro definitions -- RGB color component specification definitions
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_RED = INT(z'0000')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GREEN = INT(z'0001')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_BLUE = INT(z'0002')

!  GLUT API macro definitions -- additional keyboard and joystick definitions
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_REPEAT_OFF = INT(z'0000')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_REPEAT_ON = INT(z'0001')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_KEY_REPEAT_DEFAULT = INT(z'0002')

INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_BUTTON_A = INT(z'0001')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_BUTTON_B = INT(z'0002')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_BUTTON_C = INT(z'0004')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_JOYSTICK_BUTTON_D = INT(z'0008')

!  GLUT API macro definitions -- game mode definitions
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_ACTIVE = INT(z'0000')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_POSSIBLE = INT(z'0001')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_WIDTH = INT(z'0002')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_HEIGHT = INT(z'0003')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_PIXEL_DEPTH = INT(z'0004')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_REFRESH_RATE = INT(z'0005')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_GAME_MODE_DISPLAY_CHANGED = INT(z'0006')

!  Initialization functions, see fglut_init.c
!  void    glutInit( int* pargc, char** argv )
PUBLIC glutInit
INTERFACE glutInit
MODULE PROCEDURE glutInit_f03
END INTERFACE glutInit
INTERFACE
SUBROUTINE glutInit_gl(pargc, argv) BIND(C,NAME="glutInit")
IMPORT
! INTEGER(GLint), DIMENSION(*) :: pargc
INTEGER(GLint) :: pargc
TYPE(C_PTR), INTENT(IN) :: argv
END SUBROUTINE glutInit_gl
END INTERFACE

!  void    glutInitWindowPosition( int x, int y )
PUBLIC glutInitWindowPosition
INTERFACE
SUBROUTINE glutInitWindowPosition(x, y) BIND(C,NAME="glutInitWindowPosition")
IMPORT
INTEGER(GLint), VALUE :: x, y
END SUBROUTINE glutInitWindowPosition
END INTERFACE

!  void    glutInitWindowSize( int width, int height )
PUBLIC glutInitWindowSize
INTERFACE
SUBROUTINE glutInitWindowSize(width, height) BIND(C,NAME="glutInitWindowSize")
IMPORT
INTEGER(GLint), VALUE :: width, height
END SUBROUTINE glutInitWindowSize
END INTERFACE

!  void    glutInitDisplayMode( unsigned int displayMode )
PUBLIC glutInitDisplayMode
INTERFACE
SUBROUTINE glutInitDisplayMode(displayMode) BIND(C,NAME="glutInitDisplayMode")
IMPORT
INTEGER(GLuint), VALUE :: displayMode
END SUBROUTINE glutInitDisplayMode
END INTERFACE

!  void    glutInitDisplayString( const char* displayMode )
PUBLIC glutInitDisplayString
INTERFACE
SUBROUTINE glutInitDisplayString(displayMode) BIND(C,NAME="glutInitDisplayString")
IMPORT
! CHARACTER, INTENT(IN) :: displayMode
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: displayMode
END SUBROUTINE glutInitDisplayString
END INTERFACE


!  Process loop function, see freeglut_main.c
!  void    glutMainLoop( void )
PUBLIC glutMainLoop
INTERFACE
SUBROUTINE glutMainLoop() BIND(C,NAME="glutMainLoop")
IMPORT
END SUBROUTINE glutMainLoop
END INTERFACE


!  Window management functions, see freeglut_window.c
!  int     glutCreateWindow( const char* title )
PUBLIC glutCreateWindow
INTERFACE
FUNCTION glutCreateWindow(title) BIND(C,NAME="glutCreateWindow")
IMPORT
INTEGER(GLint) :: glutCreateWindow
! CHARACTER, INTENT(IN) :: title
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: title
END FUNCTION glutCreateWindow
END INTERFACE

!  int     glutCreateSubWindow( int window, int x, int y, int width, int height )
PUBLIC glutCreateSubWindow
INTERFACE
FUNCTION glutCreateSubWindow(window, x, y, width, height) BIND(C,NAME="glutCreateSubWindow")
IMPORT
INTEGER(GLint) :: glutCreateSubWindow
INTEGER(GLint), VALUE :: window, x, y, width, height
END FUNCTION glutCreateSubWindow
END INTERFACE

!  void    glutDestroyWindow( int window )
PUBLIC glutDestroyWindow
INTERFACE
SUBROUTINE glutDestroyWindow(window) BIND(C,NAME="glutDestroyWindow")
IMPORT
INTEGER(GLint), VALUE :: window
END SUBROUTINE glutDestroyWindow
END INTERFACE

!  void    glutSetWindow( int window )
PUBLIC glutSetWindow
INTERFACE
SUBROUTINE glutSetWindow(window) BIND(C,NAME="glutSetWindow")
IMPORT
INTEGER(GLint), VALUE :: window
END SUBROUTINE glutSetWindow
END INTERFACE

!  int     glutGetWindow( void )
PUBLIC glutGetWindow
INTERFACE
FUNCTION glutGetWindow() BIND(C,NAME="glutGetWindow")
IMPORT
INTEGER(GLint) :: glutGetWindow
END FUNCTION glutGetWindow
END INTERFACE

!  void    glutSetWindowTitle( const char* title )
PUBLIC glutSetWindowTitle
INTERFACE
SUBROUTINE glutSetWindowTitle(title) BIND(C,NAME="glutSetWindowTitle")
IMPORT
! CHARACTER, INTENT(IN) :: title
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: title
END SUBROUTINE glutSetWindowTitle
END INTERFACE

!  void    glutSetIconTitle( const char* title )
PUBLIC glutSetIconTitle
INTERFACE
SUBROUTINE glutSetIconTitle(title) BIND(C,NAME="glutSetIconTitle")
IMPORT
! CHARACTER, INTENT(IN) :: title
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: title
END SUBROUTINE glutSetIconTitle
END INTERFACE

!  void    glutReshapeWindow( int width, int height )
PUBLIC glutReshapeWindow
INTERFACE
SUBROUTINE glutReshapeWindow(width, height) BIND(C,NAME="glutReshapeWindow")
IMPORT
INTEGER(GLint), VALUE :: width, height
END SUBROUTINE glutReshapeWindow
END INTERFACE

!  void    glutPositionWindow( int x, int y )
PUBLIC glutPositionWindow
INTERFACE
SUBROUTINE glutPositionWindow(x, y) BIND(C,NAME="glutPositionWindow")
IMPORT
INTEGER(GLint), VALUE :: x, y
END SUBROUTINE glutPositionWindow
END INTERFACE

!  void    glutShowWindow( void )
PUBLIC glutShowWindow
INTERFACE
SUBROUTINE glutShowWindow() BIND(C,NAME="glutShowWindow")
IMPORT
END SUBROUTINE glutShowWindow
END INTERFACE

!  void    glutHideWindow( void )
PUBLIC glutHideWindow
INTERFACE
SUBROUTINE glutHideWindow() BIND(C,NAME="glutHideWindow")
IMPORT
END SUBROUTINE glutHideWindow
END INTERFACE

!  void    glutIconifyWindow( void )
PUBLIC glutIconifyWindow
INTERFACE
SUBROUTINE glutIconifyWindow() BIND(C,NAME="glutIconifyWindow")
IMPORT
END SUBROUTINE glutIconifyWindow
END INTERFACE

!  void    glutPushWindow( void )
PUBLIC glutPushWindow
INTERFACE
SUBROUTINE glutPushWindow() BIND(C,NAME="glutPushWindow")
IMPORT
END SUBROUTINE glutPushWindow
END INTERFACE

!  void    glutPopWindow( void )
PUBLIC glutPopWindow
INTERFACE
SUBROUTINE glutPopWindow() BIND(C,NAME="glutPopWindow")
IMPORT
END SUBROUTINE glutPopWindow
END INTERFACE

!  void    glutFullScreen( void )
PUBLIC glutFullScreen
INTERFACE
SUBROUTINE glutFullScreen() BIND(C,NAME="glutFullScreen")
IMPORT
END SUBROUTINE glutFullScreen
END INTERFACE


!  Display-connected functions, see freeglut_display.c
!  void    glutPostWindowRedisplay( int window )
PUBLIC glutPostWindowRedisplay
INTERFACE
SUBROUTINE glutPostWindowRedisplay(window) BIND(C,NAME="glutPostWindowRedisplay")
IMPORT
INTEGER(GLint), VALUE :: window
END SUBROUTINE glutPostWindowRedisplay
END INTERFACE

!  void    glutPostRedisplay( void )
PUBLIC glutPostRedisplay
INTERFACE
SUBROUTINE glutPostRedisplay() BIND(C,NAME="glutPostRedisplay")
IMPORT
END SUBROUTINE glutPostRedisplay
END INTERFACE

!  void    glutSwapBuffers( void )
PUBLIC glutSwapBuffers
INTERFACE
SUBROUTINE glutSwapBuffers() BIND(C,NAME="glutSwapBuffers")
IMPORT
END SUBROUTINE glutSwapBuffers
END INTERFACE


!  Mouse cursor functions, see freeglut_cursor.c
!  void    glutWarpPointer( int x, int y )
PUBLIC glutWarpPointer
INTERFACE
SUBROUTINE glutWarpPointer(x, y) BIND(C,NAME="glutWarpPointer")
IMPORT
INTEGER(GLint), VALUE :: x, y
END SUBROUTINE glutWarpPointer
END INTERFACE

!  void    glutSetCursor( int cursor )
PUBLIC glutSetCursor
INTERFACE
SUBROUTINE glutSetCursor(cursor) BIND(C,NAME="glutSetCursor")
IMPORT
INTEGER(GLint), VALUE :: cursor
END SUBROUTINE glutSetCursor
END INTERFACE


!  Overlay stuff, see freeglut_overlay.c
!  void    glutEstablishOverlay( void )
PUBLIC glutEstablishOverlay
INTERFACE
SUBROUTINE glutEstablishOverlay() BIND(C,NAME="glutEstablishOverlay")
IMPORT
END SUBROUTINE glutEstablishOverlay
END INTERFACE

!  void    glutRemoveOverlay( void )
PUBLIC glutRemoveOverlay
INTERFACE
SUBROUTINE glutRemoveOverlay() BIND(C,NAME="glutRemoveOverlay")
IMPORT
END SUBROUTINE glutRemoveOverlay
END INTERFACE

!  void    glutUseLayer( GLenum layer )
PUBLIC glutUseLayer
INTERFACE
SUBROUTINE glutUseLayer(layer) BIND(C,NAME="glutUseLayer")
IMPORT
INTEGER(GLenum), VALUE :: layer
END SUBROUTINE glutUseLayer
END INTERFACE

!  void    glutPostOverlayRedisplay( void )
PUBLIC glutPostOverlayRedisplay
INTERFACE
SUBROUTINE glutPostOverlayRedisplay() BIND(C,NAME="glutPostOverlayRedisplay")
IMPORT
END SUBROUTINE glutPostOverlayRedisplay
END INTERFACE

!  void    glutPostWindowOverlayRedisplay( int window )
PUBLIC glutPostWindowOverlayRedisplay
INTERFACE
SUBROUTINE glutPostWindowOverlayRedisplay(window) BIND(C,NAME="glutPostWindowOverlayRedisplay")
IMPORT
INTEGER(GLint), VALUE :: window
END SUBROUTINE glutPostWindowOverlayRedisplay
END INTERFACE

!  void    glutShowOverlay( void )
PUBLIC glutShowOverlay
INTERFACE
SUBROUTINE glutShowOverlay() BIND(C,NAME="glutShowOverlay")
IMPORT
END SUBROUTINE glutShowOverlay
END INTERFACE

!  void    glutHideOverlay( void )
PUBLIC glutHideOverlay
INTERFACE
SUBROUTINE glutHideOverlay() BIND(C,NAME="glutHideOverlay")
IMPORT
END SUBROUTINE glutHideOverlay
END INTERFACE


!  Menu stuff, see freeglut_menu.c
!  int     glutCreateMenu( void (* callback)( int menu ) )
PUBLIC glutCreateMenu
INTERFACE
FUNCTION glutCreateMenu(proc) BIND(C,NAME="glutCreateMenu")
IMPORT
INTEGER(GLint) :: glutCreateMenu
!  void proc( int menu )
INTERFACE
SUBROUTINE proc(menu) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: menu
END SUBROUTINE proc
END INTERFACE
END FUNCTION glutCreateMenu
END INTERFACE

!  void    glutDestroyMenu( int menu )
PUBLIC glutDestroyMenu
INTERFACE
SUBROUTINE glutDestroyMenu(menu) BIND(C,NAME="glutDestroyMenu")
IMPORT
INTEGER(GLint), VALUE :: menu
END SUBROUTINE glutDestroyMenu
END INTERFACE

!  int     glutGetMenu( void )
PUBLIC glutGetMenu
INTERFACE
FUNCTION glutGetMenu() BIND(C,NAME="glutGetMenu")
IMPORT
INTEGER(GLint) :: glutGetMenu
END FUNCTION glutGetMenu
END INTERFACE

!  void    glutSetMenu( int menu )
PUBLIC glutSetMenu
INTERFACE
SUBROUTINE glutSetMenu(menu) BIND(C,NAME="glutSetMenu")
IMPORT
INTEGER(GLint), VALUE :: menu
END SUBROUTINE glutSetMenu
END INTERFACE

!  void    glutAddMenuEntry( const char* label, int value )
PUBLIC glutAddMenuEntry
INTERFACE
SUBROUTINE glutAddMenuEntry(label, value) BIND(C,NAME="glutAddMenuEntry")
IMPORT
! CHARACTER, INTENT(IN) :: label
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: label
INTEGER(GLint), VALUE :: value
END SUBROUTINE glutAddMenuEntry
END INTERFACE

!  void    glutAddSubMenu( const char* label, int subMenu )
PUBLIC glutAddSubMenu
INTERFACE
SUBROUTINE glutAddSubMenu(label, subMenu) BIND(C,NAME="glutAddSubMenu")
IMPORT
! CHARACTER, INTENT(IN) :: label
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: label
INTEGER(GLint), VALUE :: subMenu
END SUBROUTINE glutAddSubMenu
END INTERFACE

!  void    glutChangeToMenuEntry( int item, const char* label, int value )
PUBLIC glutChangeToMenuEntry
INTERFACE
SUBROUTINE glutChangeToMenuEntry(item, label, value) BIND(C,NAME="glutChangeToMenuEntry")
IMPORT
! CHARACTER, INTENT(IN) :: label
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: label
INTEGER(GLint), VALUE :: item, value
END SUBROUTINE glutChangeToMenuEntry
END INTERFACE

!  void    glutChangeToSubMenu( int item, const char* label, int value )
PUBLIC glutChangeToSubMenu
INTERFACE
SUBROUTINE glutChangeToSubMenu(item, label, value) BIND(C,NAME="glutChangeToSubMenu")
IMPORT
! CHARACTER, INTENT(IN) :: label
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: label
INTEGER(GLint), VALUE :: item, value
END SUBROUTINE glutChangeToSubMenu
END INTERFACE

!  void    glutRemoveMenuItem( int item )
PUBLIC glutRemoveMenuItem
INTERFACE
SUBROUTINE glutRemoveMenuItem(item) BIND(C,NAME="glutRemoveMenuItem")
IMPORT
INTEGER(GLint), VALUE :: item
END SUBROUTINE glutRemoveMenuItem
END INTERFACE

!  void    glutAttachMenu( int button )
PUBLIC glutAttachMenu
INTERFACE
SUBROUTINE glutAttachMenu(button) BIND(C,NAME="glutAttachMenu")
IMPORT
INTEGER(GLint), VALUE :: button
END SUBROUTINE glutAttachMenu
END INTERFACE

!  void    glutDetachMenu( int button )
PUBLIC glutDetachMenu
INTERFACE
SUBROUTINE glutDetachMenu(button) BIND(C,NAME="glutDetachMenu")
IMPORT
INTEGER(GLint), VALUE :: button
END SUBROUTINE glutDetachMenu
END INTERFACE


!  Global callback functions, see freeglut_callbacks.c
!  void    glutTimerFunc( unsigned int time, void (* callback)( int ), int value )
PUBLIC glutTimerFunc
INTERFACE glutTimerFunc
MODULE PROCEDURE glutTimerFunc_f03
END INTERFACE glutTimerFunc
INTERFACE
SUBROUTINE glutTimerFunc_gl(time, proc, value) BIND(C,NAME="glutTimerFunc")
IMPORT
INTEGER(GLint), VALUE :: value
INTEGER(GLuint), VALUE :: time
!  void proc( int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutTimerFunc_gl
END INTERFACE

!  void    glutIdleFunc( void (* callback)( void ) )
PUBLIC glutIdleFunc
INTERFACE glutIdleFunc
MODULE PROCEDURE glutIdleFunc_f03
END INTERFACE glutIdleFunc
INTERFACE
SUBROUTINE glutIdleFunc_gl(proc) BIND(C,NAME="glutIdleFunc")
IMPORT
!  void proc( void )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutIdleFunc_gl
END INTERFACE


!  Window-specific callback functions, see freeglut_callbacks.c
!  void    glutKeyboardFunc( void (* callback)( unsigned char, int, int ) )
PUBLIC glutKeyboardFunc
INTERFACE glutKeyboardFunc
MODULE PROCEDURE glutKeyboardFunc_f03
END INTERFACE glutKeyboardFunc
INTERFACE
SUBROUTINE glutKeyboardFunc_gl(proc) BIND(C,NAME="glutKeyboardFunc")
IMPORT
!  void proc( unsigned char, int, int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutKeyboardFunc_gl
END INTERFACE

!  void    glutSpecialFunc( void (* callback)( int, int, int ) )
PUBLIC glutSpecialFunc
INTERFACE glutSpecialFunc
MODULE PROCEDURE glutSpecialFunc_f03
END INTERFACE glutSpecialFunc
INTERFACE
SUBROUTINE glutSpecialFunc_gl(proc) BIND(C,NAME="glutSpecialFunc")
IMPORT
!  void proc( int, int, int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutSpecialFunc_gl
END INTERFACE

!  void    glutReshapeFunc( void (* callback)( int, int ) )
PUBLIC glutReshapeFunc
INTERFACE glutReshapeFunc
MODULE PROCEDURE glutReshapeFunc_f03
END INTERFACE glutReshapeFunc
INTERFACE
SUBROUTINE glutReshapeFunc_gl(proc) BIND(C,NAME="glutReshapeFunc")
IMPORT
!  void proc( int, int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutReshapeFunc_gl
END INTERFACE

!  void    glutVisibilityFunc( void (* callback)( int ) )
PUBLIC glutVisibilityFunc
INTERFACE glutVisibilityFunc
MODULE PROCEDURE glutVisibilityFunc_f03
END INTERFACE glutVisibilityFunc
INTERFACE
SUBROUTINE glutVisibilityFunc_gl(proc) BIND(C,NAME="glutVisibilityFunc")
IMPORT
!  void proc( int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutVisibilityFunc_gl
END INTERFACE

!  void    glutDisplayFunc( void (* callback)( void ) )
PUBLIC glutDisplayFunc
INTERFACE glutDisplayFunc
MODULE PROCEDURE glutDisplayFunc_f03
END INTERFACE glutDisplayFunc
INTERFACE
SUBROUTINE glutDisplayFunc_gl(proc) BIND(C,NAME="glutDisplayFunc")
IMPORT
!  void proc( void )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutDisplayFunc_gl
END INTERFACE

!  void    glutMouseFunc( void (* callback)( int, int, int, int ) )
PUBLIC glutMouseFunc
INTERFACE glutMouseFunc
MODULE PROCEDURE glutMouseFunc_f03
END INTERFACE glutMouseFunc
INTERFACE
SUBROUTINE glutMouseFunc_gl(proc) BIND(C,NAME="glutMouseFunc")
IMPORT
!  void proc( int, int, int, int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutMouseFunc_gl
END INTERFACE

!  void    glutMotionFunc( void (* callback)( int, int ) )
PUBLIC glutMotionFunc
INTERFACE glutMotionFunc
MODULE PROCEDURE glutMotionFunc_f03
END INTERFACE glutMotionFunc
INTERFACE
SUBROUTINE glutMotionFunc_gl(proc) BIND(C,NAME="glutMotionFunc")
IMPORT
!  void proc( int, int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutMotionFunc_gl
END INTERFACE

!  void    glutPassiveMotionFunc( void (* callback)( int, int ) )
PUBLIC glutPassiveMotionFunc
INTERFACE glutPassiveMotionFunc
MODULE PROCEDURE glutPassiveMotionFunc_f03
END INTERFACE glutPassiveMotionFunc
INTERFACE
SUBROUTINE glutPassiveMotionFunc_gl(proc) BIND(C,NAME="glutPassiveMotionFunc")
IMPORT
!  void proc( int, int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutPassiveMotionFunc_gl
END INTERFACE

!  void    glutEntryFunc( void (* callback)( int ) )
PUBLIC glutEntryFunc
INTERFACE glutEntryFunc
MODULE PROCEDURE glutEntryFunc_f03
END INTERFACE glutEntryFunc
INTERFACE
SUBROUTINE glutEntryFunc_gl(proc) BIND(C,NAME="glutEntryFunc")
IMPORT
!  void proc( int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutEntryFunc_gl
END INTERFACE


!  void    glutKeyboardUpFunc( void (* callback)( unsigned char, int, int ) )
PUBLIC glutKeyboardUpFunc
INTERFACE glutKeyboardUpFunc
MODULE PROCEDURE glutKeyboardUpFunc_f03
END INTERFACE glutKeyboardUpFunc
INTERFACE
SUBROUTINE glutKeyboardUpFunc_gl(proc) BIND(C,NAME="glutKeyboardUpFunc")
IMPORT
!  void proc( unsigned char, int, int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutKeyboardUpFunc_gl
END INTERFACE

!  void    glutSpecialUpFunc( void (* callback)( int, int, int ) )
PUBLIC glutSpecialUpFunc
INTERFACE glutSpecialUpFunc
MODULE PROCEDURE glutSpecialUpFunc_f03
END INTERFACE glutSpecialUpFunc
INTERFACE
SUBROUTINE glutSpecialUpFunc_gl(proc) BIND(C,NAME="glutSpecialUpFunc")
IMPORT
!  void proc( int, int, int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutSpecialUpFunc_gl
END INTERFACE

!  void    glutJoystickFunc( void (* callback)( unsigned int, int, int, int ), int pollInterval )
PUBLIC glutJoystickFunc
INTERFACE glutJoystickFunc
MODULE PROCEDURE glutJoystickFunc_f03
END INTERFACE glutJoystickFunc
INTERFACE
SUBROUTINE glutJoystickFunc_gl(proc, pollInterval) BIND(C,NAME="glutJoystickFunc")
IMPORT
INTEGER(GLint), VALUE :: pollInterval
!  void proc( unsigned int, int, int, int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutJoystickFunc_gl
END INTERFACE

!  void    glutMenuStateFunc( void (* callback)( int ) )
PUBLIC glutMenuStateFunc
INTERFACE glutMenuStateFunc
MODULE PROCEDURE glutMenuStateFunc_f03
END INTERFACE glutMenuStateFunc
INTERFACE
SUBROUTINE glutMenuStateFunc_gl(proc) BIND(C,NAME="glutMenuStateFunc")
IMPORT
!  void proc( int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutMenuStateFunc_gl
END INTERFACE

!  void    glutMenuStatusFunc( void (* callback)( int, int, int ) )
PUBLIC glutMenuStatusFunc
INTERFACE glutMenuStatusFunc
MODULE PROCEDURE glutMenuStatusFunc_f03
END INTERFACE glutMenuStatusFunc
INTERFACE
SUBROUTINE glutMenuStatusFunc_gl(proc) BIND(C,NAME="glutMenuStatusFunc")
IMPORT
!  void proc( int, int, int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutMenuStatusFunc_gl
END INTERFACE

!  void    glutOverlayDisplayFunc( void (* callback)( void ) )
PUBLIC glutOverlayDisplayFunc
INTERFACE glutOverlayDisplayFunc
MODULE PROCEDURE glutOverlayDisplayFunc_f03
END INTERFACE glutOverlayDisplayFunc
INTERFACE
SUBROUTINE glutOverlayDisplayFunc_gl(proc) BIND(C,NAME="glutOverlayDisplayFunc")
IMPORT
!  void proc( void )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutOverlayDisplayFunc_gl
END INTERFACE

!  void    glutWindowStatusFunc( void (* callback)( int ) )
PUBLIC glutWindowStatusFunc
INTERFACE glutWindowStatusFunc
MODULE PROCEDURE glutWindowStatusFunc_f03
END INTERFACE glutWindowStatusFunc
INTERFACE
SUBROUTINE glutWindowStatusFunc_gl(proc) BIND(C,NAME="glutWindowStatusFunc")
IMPORT
!  void proc( int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutWindowStatusFunc_gl
END INTERFACE


!  void    glutSpaceballMotionFunc( void (* callback)( int, int, int ) )
PUBLIC glutSpaceballMotionFunc
INTERFACE glutSpaceballMotionFunc
MODULE PROCEDURE glutSpaceballMotionFunc_f03
END INTERFACE glutSpaceballMotionFunc
INTERFACE
SUBROUTINE glutSpaceballMotionFunc_gl(proc) BIND(C,NAME="glutSpaceballMotionFunc")
IMPORT
!  void proc( int, int, int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutSpaceballMotionFunc_gl
END INTERFACE

!  void    glutSpaceballRotateFunc( void (* callback)( int, int, int ) )
PUBLIC glutSpaceballRotateFunc
INTERFACE glutSpaceballRotateFunc
MODULE PROCEDURE glutSpaceballRotateFunc_f03
END INTERFACE glutSpaceballRotateFunc
INTERFACE
SUBROUTINE glutSpaceballRotateFunc_gl(proc) BIND(C,NAME="glutSpaceballRotateFunc")
IMPORT
!  void proc( int, int, int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutSpaceballRotateFunc_gl
END INTERFACE

!  void    glutSpaceballButtonFunc( void (* callback)( int, int ) )
PUBLIC glutSpaceballButtonFunc
INTERFACE glutSpaceballButtonFunc
MODULE PROCEDURE glutSpaceballButtonFunc_f03
END INTERFACE glutSpaceballButtonFunc
INTERFACE
SUBROUTINE glutSpaceballButtonFunc_gl(proc) BIND(C,NAME="glutSpaceballButtonFunc")
IMPORT
!  void proc( int, int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutSpaceballButtonFunc_gl
END INTERFACE

!  void    glutButtonBoxFunc( void (* callback)( int, int ) )
PUBLIC glutButtonBoxFunc
INTERFACE glutButtonBoxFunc
MODULE PROCEDURE glutButtonBoxFunc_f03
END INTERFACE glutButtonBoxFunc
INTERFACE
SUBROUTINE glutButtonBoxFunc_gl(proc) BIND(C,NAME="glutButtonBoxFunc")
IMPORT
!  void proc( int, int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutButtonBoxFunc_gl
END INTERFACE

!  void    glutDialsFunc( void (* callback)( int, int ) )
PUBLIC glutDialsFunc
INTERFACE glutDialsFunc
MODULE PROCEDURE glutDialsFunc_f03
END INTERFACE glutDialsFunc
INTERFACE
SUBROUTINE glutDialsFunc_gl(proc) BIND(C,NAME="glutDialsFunc")
IMPORT
!  void proc( int, int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutDialsFunc_gl
END INTERFACE

!  void    glutTabletMotionFunc( void (* callback)( int, int ) )
PUBLIC glutTabletMotionFunc
INTERFACE glutTabletMotionFunc
MODULE PROCEDURE glutTabletMotionFunc_f03
END INTERFACE glutTabletMotionFunc
INTERFACE
SUBROUTINE glutTabletMotionFunc_gl(proc) BIND(C,NAME="glutTabletMotionFunc")
IMPORT
!  void proc( int, int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutTabletMotionFunc_gl
END INTERFACE

!  void    glutTabletButtonFunc( void (* callback)( int, int, int, int ) )
PUBLIC glutTabletButtonFunc
INTERFACE glutTabletButtonFunc
MODULE PROCEDURE glutTabletButtonFunc_f03
END INTERFACE glutTabletButtonFunc
INTERFACE
SUBROUTINE glutTabletButtonFunc_gl(proc) BIND(C,NAME="glutTabletButtonFunc")
IMPORT
!  void proc( int, int, int, int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutTabletButtonFunc_gl
END INTERFACE


!  State setting and retrieval functions, see freeglut_state.c
!  int     glutGet( GLenum query )
PUBLIC glutGet
INTERFACE
FUNCTION glutGet(query) BIND(C,NAME="glutGet")
IMPORT
INTEGER(GLint) :: glutGet
INTEGER(GLenum), VALUE :: query
END FUNCTION glutGet
END INTERFACE

!  int     glutDeviceGet( GLenum query )
PUBLIC glutDeviceGet
INTERFACE
FUNCTION glutDeviceGet(query) BIND(C,NAME="glutDeviceGet")
IMPORT
INTEGER(GLint) :: glutDeviceGet
INTEGER(GLenum), VALUE :: query
END FUNCTION glutDeviceGet
END INTERFACE

!  int     glutGetModifiers( void )
PUBLIC glutGetModifiers
INTERFACE
FUNCTION glutGetModifiers() BIND(C,NAME="glutGetModifiers")
IMPORT
INTEGER(GLint) :: glutGetModifiers
END FUNCTION glutGetModifiers
END INTERFACE

!  int     glutLayerGet( GLenum query )
PUBLIC glutLayerGet
INTERFACE
FUNCTION glutLayerGet(query) BIND(C,NAME="glutLayerGet")
IMPORT
INTEGER(GLint) :: glutLayerGet
INTEGER(GLenum), VALUE :: query
END FUNCTION glutLayerGet
END INTERFACE


!  Font stuff, see freeglut_font.c
!  void    glutBitmapCharacter( void* font, int character )
PUBLIC glutBitmapCharacter
INTERFACE
SUBROUTINE glutBitmapCharacter(font, character) BIND(C,NAME="glutBitmapCharacter")
IMPORT
INTEGER(GLint), VALUE :: character
TYPE(C_PTR), VALUE :: font
END SUBROUTINE glutBitmapCharacter
END INTERFACE

!  int     glutBitmapWidth( void* font, int character )
PUBLIC glutBitmapWidth
INTERFACE
FUNCTION glutBitmapWidth(font, character) BIND(C,NAME="glutBitmapWidth")
IMPORT
INTEGER(GLint) :: glutBitmapWidth
INTEGER(GLint), VALUE :: character
TYPE(C_PTR), VALUE :: font
END FUNCTION glutBitmapWidth
END INTERFACE

!  void    glutStrokeCharacter( void* font, int character )
PUBLIC glutStrokeCharacter
INTERFACE
SUBROUTINE glutStrokeCharacter(font, character) BIND(C,NAME="glutStrokeCharacter")
IMPORT
INTEGER(GLint), VALUE :: character
TYPE(C_PTR), VALUE :: font
END SUBROUTINE glutStrokeCharacter
END INTERFACE

!  int     glutStrokeWidth( void* font, int character )
PUBLIC glutStrokeWidth
INTERFACE
FUNCTION glutStrokeWidth(font, character) BIND(C,NAME="glutStrokeWidth")
IMPORT
INTEGER(GLint) :: glutStrokeWidth
INTEGER(GLint), VALUE :: character
TYPE(C_PTR), VALUE :: font
END FUNCTION glutStrokeWidth
END INTERFACE

!  int     glutBitmapLength( void* font, const unsigned char* string )
PUBLIC glutBitmapLength
INTERFACE
FUNCTION glutBitmapLength(font, string) BIND(C,NAME="glutBitmapLength")
IMPORT
INTEGER(GLint) :: glutBitmapLength
! CHARACTER, INTENT(IN) :: string
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: string
TYPE(C_PTR), VALUE :: font
END FUNCTION glutBitmapLength
END INTERFACE

!  int     glutStrokeLength( void* font, const unsigned char* string )
PUBLIC glutStrokeLength
INTERFACE
FUNCTION glutStrokeLength(font, string) BIND(C,NAME="glutStrokeLength")
IMPORT
INTEGER(GLint) :: glutStrokeLength
! CHARACTER, INTENT(IN) :: string
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: string
TYPE(C_PTR), VALUE :: font
END FUNCTION glutStrokeLength
END INTERFACE


!  Geometry functions, see freeglut_geometry.c
!  void    glutWireCube( GLdouble size )
PUBLIC glutWireCube
INTERFACE
SUBROUTINE glutWireCube(size) BIND(C,NAME="glutWireCube")
IMPORT
REAL(GLdouble), VALUE :: size
END SUBROUTINE glutWireCube
END INTERFACE

!  void    glutSolidCube( GLdouble size )
PUBLIC glutSolidCube
INTERFACE
SUBROUTINE glutSolidCube(size) BIND(C,NAME="glutSolidCube")
IMPORT
REAL(GLdouble), VALUE :: size
END SUBROUTINE glutSolidCube
END INTERFACE

!  void    glutWireSphere( GLdouble radius, GLint slices, GLint stacks )
PUBLIC glutWireSphere
INTERFACE
SUBROUTINE glutWireSphere(radius, slices, stacks) BIND(C,NAME="glutWireSphere")
IMPORT
INTEGER(GLint), VALUE :: slices, stacks
REAL(GLdouble), VALUE :: radius
END SUBROUTINE glutWireSphere
END INTERFACE

!  void    glutSolidSphere( GLdouble radius, GLint slices, GLint stacks )
PUBLIC glutSolidSphere
INTERFACE
SUBROUTINE glutSolidSphere(radius, slices, stacks) BIND(C,NAME="glutSolidSphere")
IMPORT
INTEGER(GLint), VALUE :: slices, stacks
REAL(GLdouble), VALUE :: radius
END SUBROUTINE glutSolidSphere
END INTERFACE

!  void    glutWireCone( GLdouble base, GLdouble height, GLint slices, GLint stacks )
PUBLIC glutWireCone
INTERFACE
SUBROUTINE glutWireCone(base, height, slices, stacks) BIND(C,NAME="glutWireCone")
IMPORT
INTEGER(GLint), VALUE :: slices, stacks
REAL(GLdouble), VALUE :: base, height
END SUBROUTINE glutWireCone
END INTERFACE

!  void    glutSolidCone( GLdouble base, GLdouble height, GLint slices, GLint stacks )
PUBLIC glutSolidCone
INTERFACE
SUBROUTINE glutSolidCone(base, height, slices, stacks) BIND(C,NAME="glutSolidCone")
IMPORT
INTEGER(GLint), VALUE :: slices, stacks
REAL(GLdouble), VALUE :: base, height
END SUBROUTINE glutSolidCone
END INTERFACE


!  void    glutWireTorus( GLdouble innerRadius, GLdouble outerRadius, GLint sides, GLint rings )
PUBLIC glutWireTorus
INTERFACE
SUBROUTINE glutWireTorus(innerRadius, outerRadius, sides, rings) BIND(C,NAME="glutWireTorus")
IMPORT
INTEGER(GLint), VALUE :: sides, rings
REAL(GLdouble), VALUE :: innerRadius, outerRadius
END SUBROUTINE glutWireTorus
END INTERFACE

!  void    glutSolidTorus( GLdouble innerRadius, GLdouble outerRadius, GLint sides, GLint rings )
PUBLIC glutSolidTorus
INTERFACE
SUBROUTINE glutSolidTorus(innerRadius, outerRadius, sides, rings) BIND(C,NAME="glutSolidTorus")
IMPORT
INTEGER(GLint), VALUE :: sides, rings
REAL(GLdouble), VALUE :: innerRadius, outerRadius
END SUBROUTINE glutSolidTorus
END INTERFACE

!  void    glutWireDodecahedron( void )
PUBLIC glutWireDodecahedron
INTERFACE
SUBROUTINE glutWireDodecahedron() BIND(C,NAME="glutWireDodecahedron")
IMPORT
END SUBROUTINE glutWireDodecahedron
END INTERFACE

!  void    glutSolidDodecahedron( void )
PUBLIC glutSolidDodecahedron
INTERFACE
SUBROUTINE glutSolidDodecahedron() BIND(C,NAME="glutSolidDodecahedron")
IMPORT
END SUBROUTINE glutSolidDodecahedron
END INTERFACE

!  void    glutWireOctahedron( void )
PUBLIC glutWireOctahedron
INTERFACE
SUBROUTINE glutWireOctahedron() BIND(C,NAME="glutWireOctahedron")
IMPORT
END SUBROUTINE glutWireOctahedron
END INTERFACE

!  void    glutSolidOctahedron( void )
PUBLIC glutSolidOctahedron
INTERFACE
SUBROUTINE glutSolidOctahedron() BIND(C,NAME="glutSolidOctahedron")
IMPORT
END SUBROUTINE glutSolidOctahedron
END INTERFACE

!  void    glutWireTetrahedron( void )
PUBLIC glutWireTetrahedron
INTERFACE
SUBROUTINE glutWireTetrahedron() BIND(C,NAME="glutWireTetrahedron")
IMPORT
END SUBROUTINE glutWireTetrahedron
END INTERFACE

!  void    glutSolidTetrahedron( void )
PUBLIC glutSolidTetrahedron
INTERFACE
SUBROUTINE glutSolidTetrahedron() BIND(C,NAME="glutSolidTetrahedron")
IMPORT
END SUBROUTINE glutSolidTetrahedron
END INTERFACE

!  void    glutWireIcosahedron( void )
PUBLIC glutWireIcosahedron
INTERFACE
SUBROUTINE glutWireIcosahedron() BIND(C,NAME="glutWireIcosahedron")
IMPORT
END SUBROUTINE glutWireIcosahedron
END INTERFACE

!  void    glutSolidIcosahedron( void )
PUBLIC glutSolidIcosahedron
INTERFACE
SUBROUTINE glutSolidIcosahedron() BIND(C,NAME="glutSolidIcosahedron")
IMPORT
END SUBROUTINE glutSolidIcosahedron
END INTERFACE


!  Teapot rendering functions, found in freeglut_teapot.c
!  void    glutWireTeapot( GLdouble size )
PUBLIC glutWireTeapot
INTERFACE
SUBROUTINE glutWireTeapot(size) BIND(C,NAME="glutWireTeapot")
IMPORT
REAL(GLdouble), VALUE :: size
END SUBROUTINE glutWireTeapot
END INTERFACE

!  void    glutSolidTeapot( GLdouble size )
PUBLIC glutSolidTeapot
INTERFACE
SUBROUTINE glutSolidTeapot(size) BIND(C,NAME="glutSolidTeapot")
IMPORT
REAL(GLdouble), VALUE :: size
END SUBROUTINE glutSolidTeapot
END INTERFACE


!  Game mode functions, see freeglut_gamemode.c
!  void    glutGameModeString( const char* string )
PUBLIC glutGameModeString
INTERFACE
SUBROUTINE glutGameModeString(string) BIND(C,NAME="glutGameModeString")
IMPORT
! CHARACTER, INTENT(IN) :: string
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: string
END SUBROUTINE glutGameModeString
END INTERFACE

!  int     glutEnterGameMode( void )
PUBLIC glutEnterGameMode
INTERFACE
FUNCTION glutEnterGameMode() BIND(C,NAME="glutEnterGameMode")
IMPORT
INTEGER(GLint) :: glutEnterGameMode
END FUNCTION glutEnterGameMode
END INTERFACE

!  void    glutLeaveGameMode( void )
PUBLIC glutLeaveGameMode
INTERFACE
SUBROUTINE glutLeaveGameMode() BIND(C,NAME="glutLeaveGameMode")
IMPORT
END SUBROUTINE glutLeaveGameMode
END INTERFACE

!  int     glutGameModeGet( GLenum query )
PUBLIC glutGameModeGet
INTERFACE
FUNCTION glutGameModeGet(query) BIND(C,NAME="glutGameModeGet")
IMPORT
INTEGER(GLint) :: glutGameModeGet
INTEGER(GLenum), VALUE :: query
END FUNCTION glutGameModeGet
END INTERFACE


!  Video resize functions, see freeglut_videoresize.c
!  int     glutVideoResizeGet( GLenum query )
PUBLIC glutVideoResizeGet
INTERFACE
FUNCTION glutVideoResizeGet(query) BIND(C,NAME="glutVideoResizeGet")
IMPORT
INTEGER(GLint) :: glutVideoResizeGet
INTEGER(GLenum), VALUE :: query
END FUNCTION glutVideoResizeGet
END INTERFACE

!  void    glutSetupVideoResizing( void )
PUBLIC glutSetupVideoResizing
INTERFACE
SUBROUTINE glutSetupVideoResizing() BIND(C,NAME="glutSetupVideoResizing")
IMPORT
END SUBROUTINE glutSetupVideoResizing
END INTERFACE

!  void    glutStopVideoResizing( void )
PUBLIC glutStopVideoResizing
INTERFACE
SUBROUTINE glutStopVideoResizing() BIND(C,NAME="glutStopVideoResizing")
IMPORT
END SUBROUTINE glutStopVideoResizing
END INTERFACE

!  void    glutVideoResize( int x, int y, int width, int height )
PUBLIC glutVideoResize
INTERFACE
SUBROUTINE glutVideoResize(x, y, width, height) BIND(C,NAME="glutVideoResize")
IMPORT
INTEGER(GLint), VALUE :: x, y, width, height
END SUBROUTINE glutVideoResize
END INTERFACE

!  void    glutVideoPan( int x, int y, int width, int height )
PUBLIC glutVideoPan
INTERFACE
SUBROUTINE glutVideoPan(x, y, width, height) BIND(C,NAME="glutVideoPan")
IMPORT
INTEGER(GLint), VALUE :: x, y, width, height
END SUBROUTINE glutVideoPan
END INTERFACE


!  Colormap functions, see freeglut_misc.c
!  void    glutSetColor( int color, GLfloat red, GLfloat green, GLfloat blue )
PUBLIC glutSetColor
INTERFACE
SUBROUTINE glutSetColor(color, red, green, blue) BIND(C,NAME="glutSetColor")
IMPORT
INTEGER(GLint), VALUE :: color
REAL(GLfloat), VALUE :: red, green, blue
END SUBROUTINE glutSetColor
END INTERFACE

!  GLfloat glutGetColor( int color, int component )
PUBLIC glutGetColor
INTERFACE
FUNCTION glutGetColor(color, component) BIND(C,NAME="glutGetColor")
IMPORT
REAL(GLfloat) :: glutGetColor
INTEGER(GLint), VALUE :: color, component
END FUNCTION glutGetColor
END INTERFACE

!  void    glutCopyColormap( int window )
PUBLIC glutCopyColormap
INTERFACE
SUBROUTINE glutCopyColormap(window) BIND(C,NAME="glutCopyColormap")
IMPORT
INTEGER(GLint), VALUE :: window
END SUBROUTINE glutCopyColormap
END INTERFACE


!  Misc keyboard and joystick functions, see freeglut_misc.c
!  void    glutIgnoreKeyRepeat( int ignore )
PUBLIC glutIgnoreKeyRepeat
INTERFACE
SUBROUTINE glutIgnoreKeyRepeat(ignore) BIND(C,NAME="glutIgnoreKeyRepeat")
IMPORT
INTEGER(GLint), VALUE :: ignore
END SUBROUTINE glutIgnoreKeyRepeat
END INTERFACE

!  void    glutSetKeyRepeat( int repeatMode )
PUBLIC glutSetKeyRepeat
INTERFACE
SUBROUTINE glutSetKeyRepeat(repeatMode) BIND(C,NAME="glutSetKeyRepeat")
IMPORT
INTEGER(GLint), VALUE :: repeatMode
END SUBROUTINE glutSetKeyRepeat
END INTERFACE

!  void    glutForceJoystickFunc( void )
PUBLIC glutForceJoystickFunc
INTERFACE
SUBROUTINE glutForceJoystickFunc() BIND(C,NAME="glutForceJoystickFunc")
IMPORT
END SUBROUTINE glutForceJoystickFunc
END INTERFACE


!  Misc functions, see freeglut_misc.c
!  int     glutExtensionSupported( const char* extension )
PUBLIC glutExtensionSupported
INTERFACE
FUNCTION glutExtensionSupported(extension) BIND(C,NAME="glutExtensionSupported")
IMPORT
INTEGER(GLint) :: glutExtensionSupported
! CHARACTER, INTENT(IN) :: extension
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: extension
END FUNCTION glutExtensionSupported
END INTERFACE

!  void    glutReportErrors( void )
PUBLIC glutReportErrors
INTERFACE
SUBROUTINE glutReportErrors() BIND(C,NAME="glutReportErrors")
IMPORT
END SUBROUTINE glutReportErrors
END INTERFACE


!  freeglut_ext.h
!  
!  The non-GLUT-compatible extensions to the freeglut library include file
!  
!  Copyright (c) 1999-2000 Pawel W. Olszta. All Rights Reserved.
!  Written by Pawel W. Olszta, <olszta@sourceforge.net>
!  Creation date: Thu Dec 2 1999
!  
!  Permission is hereby granted, free of charge, to any person obtaining a
!  copy of this software and associated documentation files (the "Software"),
!  to deal in the Software without restriction, including without limitation
!  the rights to use, copy, modify, merge, publish, distribute, sublicense,
!  and/or sell copies of the Software, and to permit persons to whom the
!  Software is furnished to do so, subject to the following conditions:
!  
!  The above copyright notice and this permission notice shall be included
!  in all copies or substantial portions of the Software.
!  
!  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
!  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
!  PAWEL W. OLSZTA BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
!  IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
!  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


!  GLUT API Extension macro definitions -- behaviour when the user clicks on an "x" to close a window
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ACTION_EXIT = 0
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ACTION_GLUTMAINLOOP_RETURNS = 1
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ACTION_CONTINUE_EXECUTION = 2

!  Create a new rendering context when the user opens a new window?
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_CREATE_NEW_CONTEXT = 0
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_USE_CURRENT_CONTEXT = 1

!  Direct/Indirect rendering context options (has meaning only in Unix/X11)
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_FORCE_INDIRECT_CONTEXT = 0
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ALLOW_DIRECT_CONTEXT = 1
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_TRY_DIRECT_CONTEXT = 2
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_FORCE_DIRECT_CONTEXT = 3

!  GLUT API Extension macro definitions -- the glutGet parameters
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_ACTION_ON_WINDOW_CLOSE = INT(z'01F9')

INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_BORDER_WIDTH = INT(z'01FA')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_WINDOW_HEADER_HEIGHT = INT(z'01FB')

INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_VERSION = INT(z'01FC')

INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_RENDERING_CONTEXT = INT(z'01FD')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_DIRECT_RENDERING = INT(z'01FE')

!  New tokens for glutInitDisplayMode.
!  Only one GLUT_AUXn bit may be used at a time.
!  Value 0x0400 is defined in OpenGLUT.
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_AUX1 = INT(z'1000')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_AUX2 = INT(z'2000')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_AUX3 = INT(z'4000')
INTEGER(GLenum), PARAMETER, PUBLIC :: GLUT_AUX4 = INT(z'8000')

!  Process loop function, see freeglut_main.c
!  void    glutMainLoopEvent( void )
PUBLIC glutMainLoopEvent
INTERFACE
SUBROUTINE glutMainLoopEvent() BIND(C,NAME="glutMainLoopEvent")
IMPORT
END SUBROUTINE glutMainLoopEvent
END INTERFACE

!  void    glutLeaveMainLoop( void )
PUBLIC glutLeaveMainLoop
INTERFACE
SUBROUTINE glutLeaveMainLoop() BIND(C,NAME="glutLeaveMainLoop")
IMPORT
END SUBROUTINE glutLeaveMainLoop
END INTERFACE


!  Window-specific callback functions, see freeglut_callbacks.c
!  void    glutMouseWheelFunc( void (* callback)( int, int, int, int ) )
PUBLIC glutMouseWheelFunc
INTERFACE glutMouseWheelFunc
MODULE PROCEDURE glutMouseWheelFunc_f03
END INTERFACE glutMouseWheelFunc
INTERFACE
SUBROUTINE glutMouseWheelFunc_gl(proc) BIND(C,NAME="glutMouseWheelFunc")
IMPORT
!  void proc( int, int, int, int )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutMouseWheelFunc_gl
END INTERFACE

!  void    glutCloseFunc( void (* callback)( void ) )
PUBLIC glutCloseFunc
INTERFACE glutCloseFunc
MODULE PROCEDURE glutCloseFunc_f03
END INTERFACE glutCloseFunc
INTERFACE
SUBROUTINE glutCloseFunc_gl(proc) BIND(C,NAME="glutCloseFunc")
IMPORT
!  void proc( void )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutCloseFunc_gl
END INTERFACE

!  void    glutWMCloseFunc( void (* callback)( void ) )
PUBLIC glutWMCloseFunc
INTERFACE glutWMCloseFunc
MODULE PROCEDURE glutWMCloseFunc_f03
END INTERFACE glutWMCloseFunc
INTERFACE
SUBROUTINE glutWMCloseFunc_gl(proc) BIND(C,NAME="glutWMCloseFunc")
IMPORT
!  void proc( void )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutWMCloseFunc_gl
END INTERFACE

!  A. Donev: Also a destruction callback for menus 
!  void    glutMenuDestroyFunc( void (* callback)( void ) )
PUBLIC glutMenuDestroyFunc
INTERFACE glutMenuDestroyFunc
MODULE PROCEDURE glutMenuDestroyFunc_f03
END INTERFACE glutMenuDestroyFunc
INTERFACE
SUBROUTINE glutMenuDestroyFunc_gl(proc) BIND(C,NAME="glutMenuDestroyFunc")
IMPORT
!  void proc( void )
TYPE(C_FUNPTR), VALUE :: proc
END SUBROUTINE glutMenuDestroyFunc_gl
END INTERFACE


!  State setting and retrieval functions, see freeglut_state.c
!  void    glutSetOption ( GLenum option_flag, int value )
PUBLIC glutSetOption
INTERFACE
SUBROUTINE glutSetOption(option_flag, value) BIND(C,NAME="glutSetOption")
IMPORT
INTEGER(GLenum), VALUE :: option_flag
INTEGER(GLint), VALUE :: value
END SUBROUTINE glutSetOption
END INTERFACE

!  A.Donev: User-data manipulation 
!  void*   glutGetWindowData( void )
PUBLIC glutGetWindowData
INTERFACE
FUNCTION glutGetWindowData() BIND(C,NAME="glutGetWindowData")
IMPORT
TYPE(C_PTR) :: glutGetWindowData
END FUNCTION glutGetWindowData
END INTERFACE

!  void    glutSetWindowData(void* data)
PUBLIC glutSetWindowData
INTERFACE
SUBROUTINE glutSetWindowData(data) BIND(C,NAME="glutSetWindowData")
IMPORT
TYPE(C_PTR), VALUE :: data
END SUBROUTINE glutSetWindowData
END INTERFACE

!  void*   glutGetMenuData( void )
PUBLIC glutGetMenuData
INTERFACE
FUNCTION glutGetMenuData() BIND(C,NAME="glutGetMenuData")
IMPORT
TYPE(C_PTR) :: glutGetMenuData
END FUNCTION glutGetMenuData
END INTERFACE

!  void    glutSetMenuData(void* data)
PUBLIC glutSetMenuData
INTERFACE
SUBROUTINE glutSetMenuData(data) BIND(C,NAME="glutSetMenuData")
IMPORT
TYPE(C_PTR), VALUE :: data
END SUBROUTINE glutSetMenuData
END INTERFACE


!  Font stuff, see freeglut_font.c
!  int     glutBitmapHeight( void* font )
PUBLIC glutBitmapHeight
INTERFACE
FUNCTION glutBitmapHeight(font) BIND(C,NAME="glutBitmapHeight")
IMPORT
INTEGER(GLint) :: glutBitmapHeight
TYPE(C_PTR), VALUE :: font
END FUNCTION glutBitmapHeight
END INTERFACE

!  GLfloat glutStrokeHeight( void* font )
PUBLIC glutStrokeHeight
INTERFACE
FUNCTION glutStrokeHeight(font) BIND(C,NAME="glutStrokeHeight")
IMPORT
REAL(GLfloat) :: glutStrokeHeight
TYPE(C_PTR), VALUE :: font
END FUNCTION glutStrokeHeight
END INTERFACE

!  void    glutBitmapString( void* font, const unsigned char *string )
PUBLIC glutBitmapString
INTERFACE
SUBROUTINE glutBitmapString(font, string) BIND(C,NAME="glutBitmapString")
IMPORT
! CHARACTER, INTENT(IN) :: string
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: string
TYPE(C_PTR), VALUE :: font
END SUBROUTINE glutBitmapString
END INTERFACE

!  void    glutStrokeString( void* font, const unsigned char *string )
PUBLIC glutStrokeString
INTERFACE
SUBROUTINE glutStrokeString(font, string) BIND(C,NAME="glutStrokeString")
IMPORT
! CHARACTER, INTENT(IN) :: string
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: string
TYPE(C_PTR), VALUE :: font
END SUBROUTINE glutStrokeString
END INTERFACE


!  Geometry functions, see freeglut_geometry.c
!  void    glutWireRhombicDodecahedron( void )
PUBLIC glutWireRhombicDodecahedron
INTERFACE
SUBROUTINE glutWireRhombicDodecahedron() BIND(C,NAME="glutWireRhombicDodecahedron")
IMPORT
END SUBROUTINE glutWireRhombicDodecahedron
END INTERFACE

!  void    glutSolidRhombicDodecahedron( void )
PUBLIC glutSolidRhombicDodecahedron
INTERFACE
SUBROUTINE glutSolidRhombicDodecahedron() BIND(C,NAME="glutSolidRhombicDodecahedron")
IMPORT
END SUBROUTINE glutSolidRhombicDodecahedron
END INTERFACE

!  void    glutWireSierpinskiSponge ( int num_levels, GLdouble offset[3], GLdouble scale )
PUBLIC glutWireSierpinskiSponge
INTERFACE
SUBROUTINE glutWireSierpinskiSponge(num_levels, offset, scale) BIND(C,NAME="glutWireSierpinskiSponge")
IMPORT
INTEGER(GLint), VALUE :: num_levels
REAL(GLdouble), DIMENSION(3) :: offset
REAL(GLdouble), VALUE :: scale
END SUBROUTINE glutWireSierpinskiSponge
END INTERFACE

!  void    glutSolidSierpinskiSponge ( int num_levels, GLdouble offset[3], GLdouble scale )
PUBLIC glutSolidSierpinskiSponge
INTERFACE
SUBROUTINE glutSolidSierpinskiSponge(num_levels, offset, scale) BIND(C,NAME="glutSolidSierpinskiSponge")
IMPORT
INTEGER(GLint), VALUE :: num_levels
REAL(GLdouble), DIMENSION(3) :: offset
REAL(GLdouble), VALUE :: scale
END SUBROUTINE glutSolidSierpinskiSponge
END INTERFACE

!  void    glutWireCylinder( GLdouble radius, GLdouble height, GLint slices, GLint stacks)
PUBLIC glutWireCylinder
INTERFACE
SUBROUTINE glutWireCylinder(radius, height, slices, stacks) BIND(C,NAME="glutWireCylinder")
IMPORT
INTEGER(GLint), VALUE :: slices, stacks
REAL(GLdouble), VALUE :: radius, height
END SUBROUTINE glutWireCylinder
END INTERFACE

!  void    glutSolidCylinder( GLdouble radius, GLdouble height, GLint slices, GLint stacks)
PUBLIC glutSolidCylinder
INTERFACE
SUBROUTINE glutSolidCylinder(radius, height, slices, stacks) BIND(C,NAME="glutSolidCylinder")
IMPORT
INTEGER(GLint), VALUE :: slices, stacks
REAL(GLdouble), VALUE :: radius, height
END SUBROUTINE glutSolidCylinder
END INTERFACE


!  Extension functions, see freeglut_ext.c
!  GLUTproc glutGetProcAddress( const char *procName )
PUBLIC glutGetProcAddress
INTERFACE
FUNCTION glutGetProcAddress(procName) BIND(C,NAME="glutGetProcAddress")
IMPORT
TYPE(C_FUNPTR) :: glutGetProcAddress
! CHARACTER, INTENT(IN) :: procName
CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: procName
END FUNCTION glutGetProcAddress
END INTERFACE



! Font variables in GLUT_fonts.c
TYPE(C_PTR), BIND(C), PUBLIC, PROTECTED :: GLUT_STROKE_ROMAN,         &
    GLUT_STROKE_MONO_ROMAN, GLUT_BITMAP_9_BY_15, GLUT_BITMAP_8_BY_13, &
    GLUT_BITMAP_TIMES_ROMAN_10, GLUT_BITMAP_TIMES_ROMAN_24,           &
    GLUT_BITMAP_HELVETICA_10, GLUT_BITMAP_HELVETICA_12,               &
    GLUT_BITMAP_HELVETICA_18

! A special callback function for compatibility with f90gl
TYPE(C_FUNPTR), PUBLIC, SAVE, BIND(C) :: GLUT_NULL_FUNC=C_NULL_FUNPTR

CONTAINS

SUBROUTINE glutInit_f03()
  INTEGER(C_INT) :: argcp=1
  TYPE(C_PTR), DIMENSION(1), TARGET :: argv=C_NULL_PTR
  CHARACTER(C_CHAR), DIMENSION(1), TARGET :: empty_string=C_NULL_CHAR

  ! A hack
  INTERFACE
   SUBROUTINE SetNullFunc() BIND(C,NAME='SetNullFunc')
   END SUBROUTINE
  END INTERFACE  

  argv(1)=C_LOC(empty_string)
  CALL glutInit_gl(argcp, C_LOC(argv))
  
END SUBROUTINE

SUBROUTINE glutTimerFunc_f03(time, proc, value)
INTEGER(GLint), VALUE :: value
INTEGER(GLuint), VALUE :: time
INTERFACE
SUBROUTINE proc(arg1) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutTimerFunc_gl(time, C_FUNLOC(proc), value)
ELSE
   CALL glutTimerFunc_gl(time, C_NULL_FUNPTR, value)
END IF
END SUBROUTINE glutTimerFunc_f03
SUBROUTINE glutIdleFunc_f03(proc)
INTERFACE
SUBROUTINE proc() BIND(C)
IMPORT
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutIdleFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutIdleFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutIdleFunc_f03
SUBROUTINE glutKeyboardFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3) BIND(C)
IMPORT
INTEGER(GLbyte), VALUE :: arg1
INTEGER(GLint), VALUE :: arg2, arg3
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutKeyboardFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutKeyboardFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutKeyboardFunc_f03
SUBROUTINE glutSpecialFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2, arg3
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutSpecialFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutSpecialFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutSpecialFunc_f03
SUBROUTINE glutReshapeFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1, arg2) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutReshapeFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutReshapeFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutReshapeFunc_f03
SUBROUTINE glutVisibilityFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutVisibilityFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutVisibilityFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutVisibilityFunc_f03
SUBROUTINE glutDisplayFunc_f03(proc)
INTERFACE
SUBROUTINE proc() BIND(C)
IMPORT
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutDisplayFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutDisplayFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutDisplayFunc_f03
SUBROUTINE glutMouseFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3, arg4) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2, arg3, arg4
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutMouseFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutMouseFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutMouseFunc_f03
SUBROUTINE glutMotionFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1, arg2) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutMotionFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutMotionFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutMotionFunc_f03
SUBROUTINE glutPassiveMotionFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1, arg2) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutPassiveMotionFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutPassiveMotionFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutPassiveMotionFunc_f03
SUBROUTINE glutEntryFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutEntryFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutEntryFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutEntryFunc_f03
SUBROUTINE glutKeyboardUpFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3) BIND(C)
IMPORT
INTEGER(GLbyte), VALUE :: arg1
INTEGER(GLint), VALUE :: arg2, arg3
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutKeyboardUpFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutKeyboardUpFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutKeyboardUpFunc_f03
SUBROUTINE glutSpecialUpFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2, arg3
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutSpecialUpFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutSpecialUpFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutSpecialUpFunc_f03
SUBROUTINE glutJoystickFunc_f03(proc, pollInterval)
INTEGER(GLint), VALUE :: pollInterval
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3, arg4) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg2, arg3, arg4
INTEGER(GLuint), VALUE :: arg1
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutJoystickFunc_gl(C_FUNLOC(proc), pollInterval)
ELSE
   CALL glutJoystickFunc_gl(C_NULL_FUNPTR, pollInterval)
END IF
END SUBROUTINE glutJoystickFunc_f03
SUBROUTINE glutMenuStateFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutMenuStateFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutMenuStateFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutMenuStateFunc_f03
SUBROUTINE glutMenuStatusFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2, arg3
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutMenuStatusFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutMenuStatusFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutMenuStatusFunc_f03
SUBROUTINE glutOverlayDisplayFunc_f03(proc)
INTERFACE
SUBROUTINE proc() BIND(C)
IMPORT
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutOverlayDisplayFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutOverlayDisplayFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutOverlayDisplayFunc_f03
SUBROUTINE glutWindowStatusFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutWindowStatusFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutWindowStatusFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutWindowStatusFunc_f03
SUBROUTINE glutSpaceballMotionFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2, arg3
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutSpaceballMotionFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutSpaceballMotionFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutSpaceballMotionFunc_f03
SUBROUTINE glutSpaceballRotateFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2, arg3
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutSpaceballRotateFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutSpaceballRotateFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutSpaceballRotateFunc_f03
SUBROUTINE glutSpaceballButtonFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1, arg2) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutSpaceballButtonFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutSpaceballButtonFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutSpaceballButtonFunc_f03
SUBROUTINE glutButtonBoxFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1, arg2) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutButtonBoxFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutButtonBoxFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutButtonBoxFunc_f03
SUBROUTINE glutDialsFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1, arg2) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutDialsFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutDialsFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutDialsFunc_f03
SUBROUTINE glutTabletMotionFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1, arg2) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutTabletMotionFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutTabletMotionFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutTabletMotionFunc_f03
SUBROUTINE glutTabletButtonFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3, arg4) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2, arg3, arg4
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutTabletButtonFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutTabletButtonFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutTabletButtonFunc_f03
SUBROUTINE glutMouseWheelFunc_f03(proc)
INTERFACE
SUBROUTINE proc(arg1, arg2, arg3, arg4) BIND(C)
IMPORT
INTEGER(GLint), VALUE :: arg1, arg2, arg3, arg4
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutMouseWheelFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutMouseWheelFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutMouseWheelFunc_f03
SUBROUTINE glutCloseFunc_f03(proc)
INTERFACE
SUBROUTINE proc() BIND(C)
IMPORT
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutCloseFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutCloseFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutCloseFunc_f03
SUBROUTINE glutWMCloseFunc_f03(proc)
INTERFACE
SUBROUTINE proc() BIND(C)
IMPORT
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutWMCloseFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutWMCloseFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutWMCloseFunc_f03
SUBROUTINE glutMenuDestroyFunc_f03(proc)
INTERFACE
SUBROUTINE proc() BIND(C)
IMPORT
END SUBROUTINE proc
END INTERFACE
OPTIONAL :: proc
IF(PRESENT(proc)) THEN
   CALL glutMenuDestroyFunc_gl(C_FUNLOC(proc))
ELSE
   CALL glutMenuDestroyFunc_gl(C_NULL_FUNPTR)
END IF
END SUBROUTINE glutMenuDestroyFunc_f03


END MODULE opengl_glut

