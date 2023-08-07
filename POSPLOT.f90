!---------------------------------------------------------------------------
subroutine plot_func(TECFILE)
use function_plotter
USE MESHADJ
USE IFPORT
USE solverds, ONLY: SOLVER_CONTROL,INPUTFILE
USE SLOPE_PSO, ONLY:SLOPE_OPTIM
!USE POS_IO
implicit none
CHARACTER(*),INTENT(IN)::TECFILE
character*8 char_time
INTEGER::I
REAL(8),ALLOCATABLE::NODE1(:,:)
CHARACTER(1024)::FILE1
CHARACTER(2)::CH1
integer :: winid, menuid, submenuid
logical :: exist

interface
    subroutine myreshape(w,h)
        integer,intent(in out)::w,h
    end subroutine
    subroutine keyboardCB(key,  x,  y)
        integer,intent(in out)::key,x,y
    end subroutine
    subroutine arrows(key, x, y)
       integer,intent(in out) :: key, x, y
    endsubroutine
    subroutine motion(x, y)
        integer, intent(in out) :: x, y
    endsubroutine
    subroutine mouse(button, state, x, y)
        integer, intent(in out) :: button, state, x, y
    endsubroutine
endinterface

IF(SOLVER_CONTROL.ISSLOPEPA==0) THEN
    CALL TIME(char_time)
    PRINT *, 'Begin to Render...Stage=glutInit',char_time
    !call glutInit()
    PRINT *, 'Begin to Render...Stage=glutInitDisplayMode',char_time
    call glutInitDisplayMode(ior(GLUT_DOUBLE,ior(GLUT_RGB,GLUT_DEPTH)))
    call glutInitWindowPosition(10_glcint,10_glcint)
    call glutInitWindowSize(800_glcint,600_glcint)

    PRINT *, 'Begin to Render...Stage=glutCreateWindow',char_time
    !winid = glutCreateWindow(trim(adjustl(POSDATA.title)))
    winid = glutCreateWindow('IFSOLVER')
    PRINT *, 'Create Window Successfully.'
    call glutInit()
ENDIF

CALL TIME(char_time)
PRINT *, 'importing data...',char_time
IF(LEN_TRIM(ADJUSTL(TECFILE))>0) THEN
    CALL TECDATA.READIN(TECFILE)
    CALL POSDATA.TEC2POSDATA(TECDATA)
    CALL TECDATA.FREE()
ELSE
    CALL POSDATA.SOLVER2POSDATA()
ENDIF


!initialize
CALL STEPPLOT.INITIALIZE(POSDATA.NSTEP,POSDATA.NSTEP,POSDATA.STEPTIME)

!mesh topology
IF(.NOT.ISINI_GMSHET) THEN
        
    CALL Initialize_et2numNodes()
    CALL ET_GMSH_EDGE_FACE()
    ISINI_GMSHET=.TRUE.
    
ENDIF

CALL TIME(char_time)
PRINT *, 'Setup mesh topology for the model mesh...',char_time
IF(ALLOCATED(NODE1)) DEALLOCATE(NODE1)
ALLOCATE(NODE1(3,POSDATA.NNODE))
DO I=1,POSDATA.NNODE
    NODE1(:,I)=POSDATA.NODE(I).COORD
ENDDO
CALL Setup_MODEL_MESHTOPO(MEDGE,NMEDGE,MFACE,NMFACE,POSDATA.ELEMENT,POSDATA.NEL,NODE1,POSDATA.NNODE,POSDATA.ESET,POSDATA.NESET)
!CALL SETUP_EDGE_ADJL(MEDGE,NMEDGE,POSDATA.ELEMENT,POSDATA.NEL,NODE1,POSDATA.NNODE,POSDATA.ESET,POSDATA.NESET)
!CALL SETUP_FACE_ADJL(MFACE,NMFACE,MEDGE,NMEDGE,POSDATA.ELEMENT,POSDATA.NEL,NODE1,POSDATA.NNODE,POSDATA.ESET,POSDATA.NESET)


CALL TIME(char_time)
PRINT *, 'Convert to Tri/Tet mesh...',char_time
call SETUP_SUB_TET4_ELEMENT(POSDATA.ELEMENT,POSDATA.NEL,POSDATA.ESET,POSDATA.NESET,POSDATA.NODE,POSDATA.NNODE)
CALL TIME(char_time)
PRINT *, 'Setup mesh topology for Tri/Tet mesh...',char_time
!CALL SETUP_EDGE_ADJL(EDGE,NEDGE,TET,NTET,NODE1,POSDATA.NNODE,POSDATA.ESET,POSDATA.NESET)
!CALL SETUP_FACE_ADJL(FACE,NFACE,EDGE,NEDGE,TET,NTET,NODE1,POSDATA.NNODE,POSDATA.ESET,POSDATA.NESET)
!CALL SETUP_ADJACENT_ELEMENT_TET(EDGE,FACE,TET,POSDATA.NDIM)
CALL Setup_TET_MESHTOPO(EDGE,NEDGE,FACE,NFACE,TET,NTET,NODE1,POSDATA.NNODE,POSDATA.ESET,POSDATA.NESET)
DEALLOCATE(NODE1)
IF(POSDATA.NDIM==2) THEN
    CALL SETUP_EDGE_BC()
ELSE
    CALL SETUP_FACE_BC()
ENDIF
CALL TIME(char_time)
PRINT *, 'Group element to subzones...',char_time
CALL SETUP_SUBZONE_TET(POSDATA.MAXX,POSDATA.MINX,POSDATA.MAXY,POSDATA.MINY,POSDATA.MAXZ,POSDATA.MINZ,TET)

IF(SOLVER_CONTROL.ISSLOPEPA/=0) THEN
    NOPLOT=.TRUE.
    PRINT *, 'DOING SLOPE PARAMETER ANALYSIS. NOPLOT WILL BE SHOW.'
    
    !EVOLUTION AOGORITHM PSO AND EA ETC...
    IF(SOLVER_CONTROL.ISSLOPEPA>10) THEN    
        CALL SLOPE_OPTIM()
        RETURN
    ENDIF
    
    CALL RESET_STREAMLINE()
    IF(STEPPLOT.NSTEP>1) THEN
        STEPPLOT.ISTEP=STEPPLOT.NSTEP
        !CALL STEPPLOT.UPDATE()
    ENDIF
    CALL SEARCH_MINIMAL_SF_SLOPE(SOLVER_CONTROL.slope_only_searchtop)
    WRITE(CH1,'(I2)') SOLVER_CONTROL.ISSLOPEPA
    FILE1=TRIM(INPUTFILE)//'_SlopePA='//TRIM(ADJUSTL(CH1))//'.dat'    
    CALL OUT_SHOWN_STREAMLINE(trim(file1),SF_SLOPE(1))
      
    FILE1=TRIM(INPUTFILE)//'_SlopePA='//TRIM(ADJUSTL(CH1))//'_FOS.dat'
    inquire(file=TRIM(FILE1), exist=exist)
    OPEN(UNIT=20,FILE=FILE1,STATUS='UNKNOWN',POSITION='APPEND')
    IF(.NOT.EXIST) WRITE(20,10)
    CALL TIME(char_time)
    WRITE(20,20) SOLVER_CONTROL.ISSLOPEPA,SOLVER_CONTROL.SLOPE_KBASE,SOLVER_CONTROL.SLOPE_KSCALE, &
                 SOLVER_CONTROL.SLOPE_KRATIO,STREAMLINE(SF_SLOPE(1)).SF_SLOPE,&
        &        STREAMLINE(SF_SLOPE(1)).V(1:2,1),STREAMLINE(SF_SLOPE(1)).V(1:2,STREAMLINE(SF_SLOPE(1)).NV), &
        &        MINVAL(STREAMLINE(SF_SLOPE(1)).V(1,:)),MAXVAL(STREAMLINE(SF_SLOPE(1)).V(1,:)), &
        &        MINVAL(STREAMLINE(SF_SLOPE(1)).V(2,:)),MAXVAL(STREAMLINE(SF_SLOPE(1)).V(2,:)), &
        &        CHAR_TIME,DATE()
    CLOSE(20)
    RETURN
ENDIF



model_radius=POSDATA.MODELR
init_lookat.x=(POSDATA.minx+POSDATA.maxx)/2.0
init_lookat.y=(POSDATA.miny+POSDATA.maxy)/2.0
init_lookat.z=(POSDATA.minz+POSDATA.maxz)/2.0
if(POSDATA.NDIM<3) then
    init_lookfrom.x=init_lookat.x
    init_lookfrom.y=init_lookat.y
    init_lookfrom.z=3.0*POSDATA.MODELR+init_lookat.z
else
    init_lookfrom.x=init_lookat.x+POSDATA.MODELR*3.0
    init_lookfrom.y=init_lookat.y+POSDATA.MODELR*3.0
    init_lookfrom.z=init_lookat.z+POSDATA.MODELR*3.0
endif

IF(POSDATA.IHEAD>0) THEN
    CONTOUR_PLOT_VARIABLE=POSDATA.IHEAD
ELSEIF(POSDATA.IDISZ>0) THEN
    CONTOUR_PLOT_VARIABLE=POSDATA.IDISZ
ELSEIF(POSDATA.IDISY>0) THEN
    CONTOUR_PLOT_VARIABLE=POSDATA.IDISY
ELSEIF(POSDATA.IDISX>0) THEN
    CONTOUR_PLOT_VARIABLE=POSDATA.IDISX
ELSEIF(POSDATA.IX>0) THEN
    CONTOUR_PLOT_VARIABLE=POSDATA.IX
ENDIF

!CONTOUR_PLOT_VARIABLE=MAX(HEAD,DISZ,DISY,DISX,LOCX)

SLICE_PLOT_VARIABLE=CONTOUR_PLOT_VARIABLE !INITIALIZATION

call initialize_contourplot(CONTOUR_PLOT_VARIABLE)



! initialize view_modifier, receiving the id for it's submenu
PRINT *, 'Make menu...',char_time
submenuid = view_modifier_init()

! create the menu

call make_menu(submenuid)

! Set the display callback
PRINT *, 'Shape...',char_time
call glutreshapefunc(myreshape)
PRINT *, 'Mouse...',char_time
call glutMouseFunc(mouse)
PRINT *, 'Motion...',char_time
call glutMotionFunc(motion)
PRINT *, 'SpecialFunc...',char_time
call glutSpecialFunc(arrows)
PRINT *, 'display...',char_time
call glutDisplayFunc(display)
!glutTimerFunc(33, timerCB, 33);             // redraw only every given millisec
!//glutIdleFunc(idleCB);                       // redraw whenever system is idle
!glutReshapeFunc(reshapeCB);
PRINT *, 'KeyboardFunc...',char_time
call glutKeyboardFunc(keyboardCB);
!glutPassiveMotionFunc(mousePassiveMotionCB);

call initGL(POSDATA.MODELR)
! Create the image
PRINT *, 'make data...',char_time
CALL STEPPLOT.UPDATE()

!call DrawSurfaceContour()
!call DrawLineContour() 
!call drawvector()
!call drawgrid()

! Let glut take over
PRINT *, 'MainLoop...',char_time
call glutMainLoop

call POSDATA.FREE()


10 FORMAT(1X,'PATYPE',3X,'KBASE',2X,'KSCALE',2X,'KRATIO',5X,'FOS',5X, &
        & 'XENTRY',2X,'YENTRY',2X,'XEXIT',3X,'YEXIT',3X,'XMIN',4X,'XMAX',4X,'YMIN',4X,'YMAX',4X,'TIME',6X,'DATE')
20 FORMAT(I7,12(1X,F7.3),2X,A8,2X,A8)

end subroutine plot_func

!---------------------------------------------------------------------------




!///////////////////////////////////////////////////////////////////////////////
!// initialize OpenGL
!// disable unused features
!///////////////////////////////////////////////////////////////////////////////
subroutine initGL(r)
use opengl_gl
implicit none
    real(gldouble),intent(in)::r
    real(glfloat)::white(4) = [1,1,1,1];

    call glShadeModel(GL_SMOOTH);                    !// shading mathod: GL_SMOOTH or GL_FLAT
    call glPixelStorei(GL_UNPACK_ALIGNMENT, 4);      !// 4-byte pixel alignment
    !
    !!// enable /disable features
    call glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    !!//glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    !!//glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    call glEnable(GL_DEPTH_TEST);
    !call glEnable(GL_LIGHTING);
    !
    !call glEnable(GL_TEXTURE_2D);
    call glEnable(GL_CULL_FACE);
    call glEnable(GL_BLEND);
    CALL glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    !// track material ambient and diffuse from surface color, call it before glEnable(GL_COLOR_MATERIAL)
    call glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    !//glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    call glEnable(GL_COLOR_MATERIAL);
	
    call glClearColor(0.9_glclampf, 0.9_glclampf, 0.9_glclampf, 1.0_glclampf)
    !call glClearColor(0._glclampf, 0._glclampf, 0._glclampf, 0._glclampf);                   !// background color
    call glClearStencil(0);                          !// clear stencil buffer
    call glClearDepth(1._GLclampd);                         !// 0 is near, 1 is far
    call glDepthFunc(GL_LEQUAL);
    

    call initLights(r);

    
    !call glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128._glfloat);
    !call glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, white);
end subroutine

!///////////////////////////////////////////////////////////////////////////////
!// initialize lights
!///////////////////////////////////////////////////////////////////////////////
subroutine initLights(r)
use opengl_gl
implicit none
    real(gldouble),intent(in)::r
    !// set up light colors (ambient, diffuse, specular)
    real(GLfloat):: lightKa(4) = [0.0_glfloat, 0.0_glfloat, 0.0_glfloat, 1.0_glfloat];  !// ambient light
    real(GLfloat):: lightKd(4) = [1.0_glfloat, 1.0_glfloat, 1.0_glfloat, 1.0_glfloat];  !// diffuse light
    real(GLfloat):: lightKs(4) = [1, 1, 1, 1];           !// specular light
        !// position the light
    real(GLfloat):: lightPos(4);
    
    lightPos(:)= [0.d0, r*3.d0, r*2.d0, 1.d0] 
    call glLightfv(GL_LIGHT0, GL_AMBIENT, lightKa);
    call glLightfv(GL_LIGHT0, GL_DIFFUSE, lightKd);
    call glLightfv(GL_LIGHT0, GL_SPECULAR, lightKs);


    call glLightfv(GL_LIGHT0, GL_POSITION, lightPos);

    call glEnable(GL_LIGHT0);                        !// MUST enable each light source after configuration
end subroutine


subroutine PickPoint(x,y,Pt1,IEL)

    use opengl_gl
    use opengl_glut
    !use solverds
    !use view_modifier
    use function_plotter
    !USE MESHGEO
	implicit none    
    integer(kind=glcint),intent(in) ::  x, y
    INTEGER,INTENT(OUT)::IEL
    real(8),intent(out)::pt1(3)
!    INTEGER,EXTERNAL::POINTlOC
   	
    

    call GetOGLPos(x, y,Pt1)

    iel=POINTlOC(pt1,0)


    
    
endsubroutine



subroutine mouse(button, state, x, y)
use function_plotter
implicit none
!          -----
integer(kind=glcint), intent(in out) :: button, state, x, y
!integer,external::POINTlOC,PTINTRIlOC
integer::iel,I
real(8)::Pt1(3)

! This gets called when a mouse button changes
  moving_left = .FALSE.
  moving_MIDDLE= .FALSE.
  left_button_func=-1
  if (button == GLUT_LEFT_BUTTON) then
    moving_left = .TRUE.
    select case(state)
    
    case(GLUT_DOWN)
        SELECT CASE(glutGetModifiers())
        CASE(GLUT_ACTIVE_CTRL)
            call ProbeatPoint(x,y)
            IF(isProbeState_SFS) CALL SLOPE_SFR_STATE_PARAMETER_SHOW( )
        CASE(GLUT_ACTIVE_SHIFT)
            call glutSetCursor(GLUT_CURSOR_LEFT_ARROW)
            begin_left = cart2D(x,y)
            left_button_func=PAN
        CASE DEFAULT
            
        
            if(isPickforstreamline.or.IsFilterLocalMinimalMode) then
                
                call PickPoint(x,y,Pt1,IEL)
                IF(iel==0) then
                    info.str='the picked location is out of zone.Please pick again.'C
                    info.color=red;info.qkey=.true.
                ELSE
                    left_button_func=LB_DRAWLINE
                    LINE_TEMP.V1=PT1
                    LINE_TEMP.V2=PT1
                    LINE_TEMP.SHOW=.TRUE.
                                     
                ENDIF
                
                
            else
                call glutSetCursor(GLUT_CURSOR_LEFT_ARROW)
			    begin_left = cart2D(x,y)
			    left_button_func=ROTATE
            endif
        END SELECT
        
    
    
    case(GLUT_UP)
        if(isPickforstreamline) then
            info.str='Click to Pick more or Press q to exit.'C
            if(info.qkey==.false.) then
                isPickforstreamline=.false.
                LINE_TEMP.SHOW=.FALSE.
                LINE_TEMP.V1=LINE_TEMP.V2
                left_button_func=ROTATE
                info.str=''
                call glutSetCursor(GLUT_CURSOR_LEFT_ARROW)
                return
            endif
			
            IF(NORM2(LINE_TEMP.V1-LINE_TEMP.V2)<1E-3) THEN
                call gen_new_streamline(LINE_TEMP.V1)
            ELSE
                DO I=1,11
                    PT1=LINE_TEMP.V1+(I-1)/10.0*(LINE_TEMP.V2-LINE_TEMP.V1)
                    call gen_new_streamline(PT1)
                ENDDO
            ENDIF
            LINE_TEMP.SHOW=.FALSE.
            LINE_TEMP.V1=LINE_TEMP.V2
            info.color=green;info.qkey=.true. 
        
        ENDIF

        if(IsFilterLocalMinimalMode) then
            info.str='Click to Pick more or Press q to exit.'C
            if(info.qkey==.false.) then
                IsFilterLocalMinimalMode=.false.
                LINE_TEMP.SHOW=.FALSE.
                LINE_TEMP.V1=LINE_TEMP.V2
                left_button_func=ROTATE
                info.str=''
                call glutSetCursor(GLUT_CURSOR_LEFT_ARROW)
                return
            endif
			
            IF(NORM2(LINE_TEMP.V1-LINE_TEMP.V2)>1E-3) THEN
                call Filterlocalminimalslope(LINE_TEMP.V1,LINE_TEMP.V2)
            ENDIF
            LINE_TEMP.SHOW=.FALSE.
            LINE_TEMP.V1=LINE_TEMP.V2
            info.color=green;info.qkey=.true. 
        
        ENDIF        
    
    end select
    
    
  endif
  
  if (button == GLUT_MIDDLE_BUTTON ) then
    if(state == GLUT_DOWN) then
        call glutSetCursor(GLUT_CURSOR_LEFT_ARROW)
        moving_middle = .true.
        begin_middle = cart2D(x,y)
        middle_button_func=ZOOM
    else
        moving_middle = .false.
    endif
    
  endif
  

end subroutine mouse

!          ------
subroutine motion(x, y)
use function_plotter
implicit none
!          ------
integer(kind=glcint), intent(in out) :: x, y
integer::iel
real(8)::Pt1(3)
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
CASE(LB_DRAWLINE)
    call PickPoint(x,y,Pt1,IEL)
    LINE_TEMP.V2=PT1
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
   shift%x = shift%x + .1*(x - begin%x)
   shift%y = shift%y - .1*(y - begin%y)
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

use function_plotter
implicit none
integer(glcint), intent(in out) :: key, x, y

! This routine handles the arrow key operations

real(kind=gldouble) :: factor

select case(glutGetModifiers())
case(GLUT_ACTIVE_CTRL) !PAN
   select case(key)
   case(GLUT_KEY_LEFT)
      shift%x = shift%x - .2
   case(GLUT_KEY_RIGHT)
      shift%x = shift%x + .2
   case(GLUT_KEY_DOWN)
      shift%y = shift%y - .2
   case(GLUT_KEY_UP)
      shift%y = shift%y + .2
   end select
case(GLUT_ACTIVE_SHIFT) !ROTATE
    select case(key)
    case(GLUT_KEY_LEFT)
        angle%x = angle%x - 1.0_gldouble
    case(GLUT_KEY_RIGHT)
        angle%x = angle%x + 1.0_gldouble
    case(GLUT_KEY_DOWN)
        angle%y = angle%y + 1.0_gldouble
    case(GLUT_KEY_UP)
        angle%y = angle%y - 1.0_gldouble
    end select    
CASE(GLUT_ACTIVE_ALT) !ZOOM
    select case(key)
    case(GLUT_KEY_DOWN)
        factor = 1.0_gldouble + .02_gldouble
    case(GLUT_KEY_UP)
        factor = 1.0_gldouble/(1.0_gldouble + .02_gldouble)
    case default
        factor = 1.0_gldouble
    end select
    IF(ISPERSPECT) THEN
        shift%z = shift%z/factor
    ELSE
        xscale_factor = xscale_factor * factor
        yscale_factor = yscale_factor * factor
        zscale_factor = zscale_factor * factor
    ENDIF
CASE DEFAULT
    select case(key)
    CASE(GLUT_KEY_DOWN,GLUT_KEY_RIGHT)
        STEPPLOT.ISTEP=MOD(STEPPLOT.ISTEP,STEPPLOT.NSTEP)+1        
    CASE(GLUT_KEY_UP,GLUT_KEY_LEFT)
        STEPPLOT.ISTEP=STEPPLOT.ISTEP-1
        IF(STEPPLOT.ISTEP<1) STEPPLOT.ISTEP=STEPPLOT.NSTEP
    ENDSELECT    
    CALL STEPPLOT.UPDATE()
end select

!select case(arrow_key_func)
!
!CASE(STEP_KEY)
!    select case(key)
!    CASE(GLUT_KEY_DOWN,GLUT_KEY_RIGHT)
!        STEPPLOT.ISTEP=MOD(STEPPLOT.ISTEP,STEPPLOT.NSTEP)+1        
!    CASE(GLUT_KEY_UP,GLUT_KEY_LEFT)
!        STEPPLOT.ISTEP=STEPPLOT.ISTEP-1
!        IF(STEPPLOT.ISTEP<1) STEPPLOT.ISTEP=STEPPLOT.NSTEP
!    ENDSELECT
!    
!    CALL STEPPLOT.UPDATE()
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

subroutine myreshape(w,h)
    use function_plotter    
    implicit none
    integer(glcint)::w,h
    


	real(gldouble):: wL,wH,R,wR,clipAreaXLeft,clipAreaXRight,clipAreaYBottom,clipAreaYTop, &
                    R1,dis1,near,far,xc,yc

    !call reset_view
    
	if (h == 0) h = 1;
	call glViewport(0, 0, w, h);

	R = real(w) / real(h);
    r1=POSDATA.MODELR/2.0
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
    use function_plotter    
    use strings
    implicit none
    integer,intent(in)::key,x,y
    real(kind=gldouble)::factor
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
            left_button_func=ROTATE
            isPickforstreamline=.FALSE.
            IsFilterLocalMinimalMode=.FALSE.
            call glutSetCursor(GLUT_CURSOR_LEFT_ARROW)
        endif
    case(ichar('+'),ichar('='))
           
      factor = 1._gldouble+0.02_gldouble
       IF(ISPERSPECT) THEN
            shift%z = shift%z/factor
       ELSE
            xscale_factor = xscale_factor * factor
            yscale_factor = yscale_factor * factor
            zscale_factor = zscale_factor * factor
       ENDIF
    case(ichar('_'),ichar('-'))
       factor = 1._gldouble/(1._gldouble+0.02_gldouble)
       IF(ISPERSPECT) THEN
            shift%z = shift%z/factor
       ELSE
            xscale_factor = xscale_factor * factor
            yscale_factor = yscale_factor * factor
            zscale_factor = zscale_factor * factor
       ENDIF       
        
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