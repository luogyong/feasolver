!---------------------------------------------------------------------------

module function_plotter
use opengl_gl
use opengl_glut
use view_modifier
use solverds
use MESHGEO
implicit none
!private
public :: display,menu_handler,make_menu,CONTOUR_PLOT_VARIABLE,VECTOR_PLOT_GROUP
private::SET_VARIABLE_SHOW,SET_VECTOR_PLOT_GROUP
! symbolic constants

!contour
integer,parameter:: Contour_surfsolid_toggle = 1, &
                    Contour_Line_toggle = 2, &
                    Contour_Densify=3,&
                    Contour_Sparsify=4,&
                    Contour_In_DeformedMesh=5
integer, parameter :: black_contour = 1, &
                      rainbow_contour = 2
integer, parameter :: white_surface = 1, &
                      red_surface = 2, &
                      rainbow_surface = 3 ,&
                      transparency=4
                      
real(8),public::Scale_Contour_Num=1.0
logical,public::istransparency=.false.,IsContour_In_DeformedMesh=.false. 

type contour_bar_tydef
    integer :: nfrac,nval
    real(GLDOUBLE),allocatable :: val(:)
    real(GLFLOAT),allocatable::COLOR(:,:)
    character(128)::TITLE
endtype
type(contour_bar_tydef)::ContourBar

!vector
INTEGER,PARAMETER:: VECTOR_GROUP_DIS=1,&
                    VECTOR_GROUP_SEEPAGE_VEC=2,&
                    VECTOR_GROUP_SEEPAGE_GRAD=3,&
                    VECTOR_GROUP_SFR=4
integer,parameter::Vector_toggle= 1, &
                   Vector_lengthen=2,&
                   Vector_shorten=3                   
                   
real(8),public::Scale_Vector_Len=1.0,VabsMax,VabsMin,Vscale                    
character(128)::VectorPairName
!model
integer, parameter :: surfgrid_toggle = 1, &
                      quit_selected = 4,&
                      DeformedMesh=5,&
                      Enlarge_Scale_DeformedMesh=6,&
                      Minify_Scale_DeformedMesh=7
                      
logical::IsDeformedMesh=.false.                     
real(8)::Scale_Deformed_Grid=1.d0

integer, parameter :: set_nx = 1, &
                      set_ny = 2, &
                      set_ncontour = 3, &
                      set_contour_val = 4, &
                      set_xrange = 5, &
                      set_yrange = 6, &
                      reset_params = 7



! Default initial settings

integer,parameter :: init_ngridx = 40, &
                      init_ngridy = 40, &                      
                      init_contour_color = black_contour, &
                      init_surface_color = rainbow_surface

real(GLDOUBLE), parameter :: init_minx = 0.0_GLDOUBLE, &
                             init_maxx = 1.0_GLDOUBLE, &
                             init_miny = 0.0_GLDOUBLE, &
                             init_maxy = 1.0_GLDOUBLE

logical, parameter :: init_draw_surface_grid = .false., &
                      init_draw_surface_solid = .true., &
                      init_draw_contour = .true.


! Current settings

integer :: ngridx = init_ngridx, &
           ngridy = init_ngridy, &
           num_contour = 20, &
           contour_color = init_contour_color, &
           surface_color = init_surface_color,&
           init_num_contour = 20

!real(GLDOUBLE) :: minx = init_minx, &
!                  maxx = init_maxx, &
!                  miny = init_miny, &
!                  maxy = init_maxy, &
!                  minz = 0.0_GLDOUBLE, &
!                  maxz = 0.0_GLDOUBLE

logical :: draw_surface_grid = init_draw_surface_grid, &
           draw_surface_solid = init_draw_surface_solid, &
           draw_contour = init_draw_contour, &
           contour_values_given = .false.,&
           IsDrawVector=.false.

real(GLDOUBLE), allocatable :: actual_contours(:)
real(GLDOUBLE) :: minv,maxv

integer::CONTOUR_PLOT_VARIABLE=0,VECTOR_PLOT_GROUP=0

integer,parameter,public::ContourList=1,&
                          VectorList=2,&
                          GridList=3,&
						  ContourLineList=4

real(GLFLOAT) :: red(4) = (/1.0,0.0,0.0,1.0/), &
                 black(4) = (/0.0,0.0,0.0,1.0/), &
                 white(4) = (/1.0,1.0,1.0,1.0/),&
                 GRAY(4)=(/0.82745098,0.82745098,0.82745098,1.0/)                          

contains

subroutine display

! This gets called when the display needs to be redrawn

call reset_view

call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))

if(glIsList(ContourList)) call glCallList(ContourList)
if(glIsList(ContourLineList)) call glCallList(ContourLineList)
if(glIsList(VectorList)) call glCallList(VectorList)
if(glIsList(GridList)) call glCallList(GridList)

call drawAxes()
if(IsDrawVector) call drawVectorLegend2(VabsMax,VabsMin,Vscale,VectorPairName)
if (surface_color == rainbow_surface) then
    call Color_Bar()
endif
call glutSwapBuffers

return
end subroutine display



function normcrossprod(x,y,z)
real(glfloat), dimension(3) :: normcrossprod
real(gldouble), dimension(3), intent(in) :: x,y,z
real(glfloat) :: t1(3),t2(3),norm
t1(1) = x(2) - x(1)
t1(2) = y(2) - y(1)
t1(3) = z(2) - z(1)
t2(1) = x(3) - x(1)
t2(2) = y(3) - y(1)
t2(3) = z(3) - z(1)
normcrossprod(1) = t1(2)*t2(3) - t1(3)*t2(2)
normcrossprod(2) = t1(3)*t2(1) - t1(1)*t2(3)
normcrossprod(3) = t1(1)*t2(2) - t1(2)*t2(1)
norm = sqrt(dot_product(normcrossprod,normcrossprod))
if (norm /= 0._glfloat) normcrossprod = normcrossprod/norm
end function normcrossprod

subroutine menu_handler(selection)
integer(kind=glcint), intent(in out) :: selection

select case (selection)

case (quit_selected)
   stop

end select

return
end subroutine menu_handler

SUBROUTINE set_variable_show(VALUE)
integer(kind=glcint), intent(in out) :: value 

CONTOUR_PLOT_VARIABLE=VALUE
call initialize_contourplot()
CALL DrawSurfaceContour()
call DrawLineContour()

RETURN
END SUBROUTINE

SUBROUTINE SET_VECTOR_PLOT_GROUP(VALUE)
integer(kind=glcint), intent(in out) :: value 

VECTOR_PLOT_GROUP=VALUE

IsDrawVector=.true.

select case(VECTOR_PLOT_GROUP)

case(VECTOR_GROUP_DIS)
    VEC(1,:)=NODALQ(:,OUTVAR(DISX).IVO)
    VEC(2,:)=NODALQ(:,OUTVAR(DISY).IVO)
    IF(NDIMENSION>2) THEN
        VEC(3,:)=NODALQ(:,OUTVAR(DISZ).IVO)
    ELSE
        VEC(3,:)=0.D0
    ENDIF
    VectorPairName='DIS.'
case(VECTOR_GROUP_SEEPAGE_VEC)       
    VEC(1,:)=NODALQ(:,OUTVAR(VX).IVO)
    VEC(2,:)=NODALQ(:,OUTVAR(VY).IVO)
    IF(NDIMENSION>2) THEN
        VEC(3,:)=NODALQ(:,OUTVAR(VZ).IVO)
    ELSE
        VEC(3,:)=0.D0
    ENDIF
    VectorPairName='SEEP.V'
case(VECTOR_GROUP_SEEPAGE_GRAD)

    VEC(1,:)=NODALQ(:,OUTVAR(GRADX).IVO)
    VEC(2,:)=NODALQ(:,OUTVAR(GRADY).IVO)
    IF(NDIMENSION>2) THEN
        VEC(3,:)=NODALQ(:,OUTVAR(GRADZ).IVO)
    ELSE
        VEC(3,:)=0.D0
    ENDIF    
    VectorPairName='SEEP.I'
case(VECTOR_GROUP_SFR)

    VEC(1,:)=NODALQ(:,OUTVAR(SFR_SFRX).IVO)
    VEC(2,:)=NODALQ(:,OUTVAR(SFR_SFRY).IVO)
    VEC(3,:)=0.D0       
    VectorPairName='SFR'
end select

CALL drawvector()

RETURN
END SUBROUTINE

subroutine Contour_handler(selection)
integer(kind=glcint), intent(in) :: selection
select case (selection)

case (Contour_surfsolid_toggle)
   draw_surface_solid = .not. draw_surface_solid

case (Contour_Line_toggle)
   draw_contour = .not. draw_contour

case(contour_densify)
    Scale_Contour_Num=Scale_Contour_Num*2
	call initialize_contourplot()
case(contour_sparsify)
    Scale_Contour_Num=Scale_Contour_Num/2
	call initialize_contourplot()
CASE(Contour_In_DeformedMesh)
    IsContour_In_DeformedMesh=.NOT.IsContour_In_DeformedMesh
end select

if(selection/=Contour_Line_toggle) call DrawSurfaceContour()
if(selection/=Contour_surfsolid_toggle) call DrawLineContour()

endsubroutine

subroutine Vector_handler(selection)
integer(kind=glcint), intent(in) :: selection
select case (selection)


case (Vector_toggle)
   IsDrawVector = .not. IsDrawVector   
case(Vector_lengthen)
   Scale_Vector_len=Scale_Vector_len*2.0
case(Vector_shorten)
   Scale_Vector_len=Scale_Vector_len/2.0
end select

call drawvector()

endsubroutine

subroutine Model_handler(selection)
integer(kind=glcint), intent(in) :: selection
select case (selection)


case (surfgrid_toggle)
   draw_surface_grid = .not. draw_surface_grid
case(DeformedMesh)
    draw_surface_grid=.true.
   IsDeformedMesh=.not.IsDeformedMesh
case(Enlarge_Scale_DeformedMesh)
    Scale_Deformed_Grid=Scale_Deformed_Grid*2.0
    if(IsContour_In_DeformedMesh) then
		call  DrawSurfaceContour()
		call  DrawLineContour()  
	endif
case(Minify_Scale_DeformedMesh)
    Scale_Deformed_Grid=Scale_Deformed_Grid/2.0  
    if(IsContour_In_DeformedMesh) then
		call  DrawSurfaceContour()
		call  DrawLineContour()  
	endif
end select

call drawgrid()

endsubroutine

subroutine param_handler(selection)
integer(kind=glcint), intent(in out) :: selection

select case (selection)

case (set_ncontour)
   print *,"Enter number of contour lines:"
   read *, init_num_contour   
   contour_values_given = .false.
   call DrawSurfaceContour()
   call DrawLineContour() 
case (set_contour_val)
   print *,"enter number of contours:"
   read *, num_contour
   if (allocated(actual_contours)) deallocate(actual_contours)
   allocate(actual_contours(num_contour))
   print *,"enter ",num_contour," contour values:"
   read *,actual_contours
   contour_values_given = .true.
   call DrawSurfaceContour()
   call DrawLineContour() 
case (reset_params)

   !num_contour = init_num_contour
   !contour_color = init_contour_color
   !surface_color = init_surface_color
   !minx = init_minx
   !maxx = init_maxx
   !miny = init_miny
   !maxy = init_maxy
   !draw_surface_grid = init_draw_surface_grid
   !draw_surface_solid = init_draw_surface_solid
   !draw_contour = init_draw_contour
   !call fem_draw

end select

end subroutine param_handler

subroutine contour_color_handler(selection)
integer(kind=glcint), intent(in out) :: selection

contour_color = selection
call DrawSurfaceContour()

end subroutine contour_color_handler

subroutine surface_color_handler(selection)
integer(kind=glcint), intent(in out) :: selection

select case(selection)
case(transparency)
    Istransparency=.not.Istransparency
case default
    surface_color = selection
end select
call DrawSurfaceContour()

end subroutine surface_color_handler

subroutine make_menu(submenuid)
integer, intent(in) :: submenuid
integer :: menuid, param_id, contour_color_menu, surface_color_menu
INTEGER::VSHOW_SPG_ID,VSHOW_STRESS_ID,VSHOW_STRAIN_ID,VSHOW_PSTRAIN_ID,&
        VSHOW_SFR_ID,VSHOW_DIS_ID,VSHOW_FORCE_ID,CONTOUR_PLOT_ID,VECTOR_PLOT_ID,&
        VECTOR_PAIR_ID,Model_ID
        
contour_color_menu = glutCreateMenu(contour_color_handler)
call glutAddMenuEntry("black",black_contour)
call glutAddMenuEntry("contour value",rainbow_contour)

surface_color_menu = glutCreateMenu(surface_color_handler)
call glutAddMenuEntry("white",white_surface)
call glutAddMenuEntry("rainbow",rainbow_surface)
call glutAddMenuEntry("transparency",transparency)

!param_id = glutCreateMenu(param_handler)
!call glutAddMenuEntry("number of x grid intervals",set_nx)
!call glutAddMenuEntry("number of y grid intervals",set_ny)
!call glutAddMenuEntry("reset to initial parameters",reset_params)

VSHOW_SPG_ID=glutCreateMenu(SET_VARIABLE_SHOW)
IF(OUTVAR(HEAD).VALUE>0) CALL GLUTADDMENUENTRY("HEAD",HEAD)
IF(OUTVAR(PHEAD).VALUE>0) CALL GLUTADDMENUENTRY("PHEAD",PHEAD)
IF(OUTVAR(discharge).VALUE>0) CALL GLUTADDMENUENTRY("Q",discharge) 
IF(OUTVAR(KR_SPG).VALUE>0) CALL GLUTADDMENUENTRY("Kr",KR_SPG)
IF(OUTVAR(MW_SPG).VALUE>0) CALL GLUTADDMENUENTRY("Mw",MW_SPG) 
IF(OUTVAR(GRADX).VALUE>0) CALL GLUTADDMENUENTRY("IX",GRADX) 
IF(OUTVAR(GRADY).VALUE>0) CALL GLUTADDMENUENTRY("IY",GRADY)
IF(OUTVAR(GRADZ).VALUE>0) CALL GLUTADDMENUENTRY("IZ",GRADZ) 
IF(OUTVAR(VX).VALUE>0)CALL GLUTADDMENUENTRY("VX",VX)
IF(OUTVAR(VY).VALUE>0) CALL GLUTADDMENUENTRY("VY",VY) 
IF(OUTVAR(VZ).VALUE>0) CALL GLUTADDMENUENTRY("VZ",VZ)

VSHOW_STRESS_ID=glutCreateMenu(SET_VARIABLE_SHOW)
IF(OUTVAR(SXX).VALUE>0) CALL GLUTADDMENUENTRY("SXX",SXX)
IF(OUTVAR(SYY).VALUE>0) CALL GLUTADDMENUENTRY("SYY",SYY)
IF(OUTVAR(SZZ).VALUE>0) CALL GLUTADDMENUENTRY("SZZ",SZZ) 
IF(OUTVAR(SXY).VALUE>0) CALL GLUTADDMENUENTRY("SXY",SXY)
IF(OUTVAR(SYZ).VALUE>0)  CALL GLUTADDMENUENTRY ("SYZ",SYZ) 
IF(OUTVAR(SZX).VALUE>0)  CALL GLUTADDMENUENTRY ("SZX",SZX) 
IF(OUTVAR(sigma_mises).VALUE>0)  CALL GLUTADDMENUENTRY ("MISES",sigma_mises) 




VSHOW_STRAIN_ID=glutCreateMenu(SET_VARIABLE_SHOW)
IF(OUTVAR(EXX).VALUE>0)  CALL GLUTADDMENUENTRY ("EXX",EXX)
IF(OUTVAR(EYY).VALUE>0)  CALL GLUTADDMENUENTRY ("EYY",EYY)
IF(OUTVAR(EZZ).VALUE>0)  CALL GLUTADDMENUENTRY ("EZZ",EZZ) 
IF(OUTVAR(EXY).VALUE>0)  CALL GLUTADDMENUENTRY ("EXY",EXY)
IF(OUTVAR(EYZ).VALUE>0)  CALL GLUTADDMENUENTRY ("EYZ",EYZ) 
IF(OUTVAR(EZX).VALUE>0)  CALL GLUTADDMENUENTRY ("EZX",EZX) 
IF(OUTVAR(EEQ).VALUE>0)  CALL GLUTADDMENUENTRY ("EEQ",EEQ) 


VSHOW_PSTRAIN_ID=glutCreateMenu(SET_VARIABLE_SHOW)
IF(OUTVAR(PEXX).VALUE>0)  CALL GLUTADDMENUENTRY ("PEXX",PEXX)
IF(OUTVAR(PEYY).VALUE>0)  CALL GLUTADDMENUENTRY ("PEYY",PEYY)
IF(OUTVAR(PEZZ).VALUE>0)  CALL GLUTADDMENUENTRY ("PEZZ",PEZZ) 
IF(OUTVAR(PEXY).VALUE>0)  CALL GLUTADDMENUENTRY ("PEXY",PEXY)
IF(OUTVAR(PEYZ).VALUE>0)  CALL GLUTADDMENUENTRY ("PEYZ",PEYZ) 
IF(OUTVAR(PEZX).VALUE>0)  CALL GLUTADDMENUENTRY ("PEZX",PEZX)
IF(OUTVAR(PEEQ).VALUE>0)  CALL GLUTADDMENUENTRY ("PEEQ",PEEQ) 

VSHOW_SFR_ID=glutCreateMenu(SET_VARIABLE_SHOW)
IF(OUTVAR(SFR).VALUE>0)  CALL GLUTADDMENUENTRY ("SFR",SFR)
IF(OUTVAR(SFR_SITA).VALUE>0)  CALL GLUTADDMENUENTRY ("SFR_SITA",SFR_SITA)
IF(OUTVAR(SFR_SN).VALUE>0)  CALL GLUTADDMENUENTRY ("Sn(Tension+)",SFR_SN) 
IF(OUTVAR(SFR_TN).VALUE>0)  CALL GLUTADDMENUENTRY ("Tn(CCW+)",SFR_TN)
IF(OUTVAR(SFR_SFRX).VALUE>0)  CALL GLUTADDMENUENTRY ("SFRX",SFR_SFRX) 
IF(OUTVAR(SFR_SFRY).VALUE>0)  CALL GLUTADDMENUENTRY ("SFRY",SFR_SFRY)

VSHOW_DIS_ID=glutCreateMenu(SET_VARIABLE_SHOW)
IF(OUTVAR(DISX).VALUE>0)  CALL GLUTADDMENUENTRY ("X",DISX)
IF(OUTVAR(DISY).VALUE>0)  CALL GLUTADDMENUENTRY ("Y",DISY)
IF(OUTVAR(DISZ).VALUE>0)  CALL GLUTADDMENUENTRY ("Z",DISZ) 
IF(OUTVAR(RX).VALUE>0)  CALL GLUTADDMENUENTRY ("RX",RX)
IF(OUTVAR(RY).VALUE>0)  CALL GLUTADDMENUENTRY ("RY",RY)
IF(OUTVAR(RZ).VALUE>0)  CALL GLUTADDMENUENTRY ("RZ",RZ)
IF(OUTVAR(HEAD).VALUE>0)  CALL GLUTADDMENUENTRY ("HEAD",HEAD)

VSHOW_FORCE_ID=glutCreateMenu(SET_VARIABLE_SHOW)
IF(OUTVAR(NFX).VALUE>0)  CALL GLUTADDMENUENTRY ("FX",NFX)
IF(OUTVAR(NFY).VALUE>0)  CALL GLUTADDMENUENTRY ("FY",NFY)
IF(OUTVAR(NFZ).VALUE>0)  CALL GLUTADDMENUENTRY ("FZ",NFZ) 
IF(OUTVAR(MX).VALUE>0)  CALL GLUTADDMENUENTRY ("MX",MX)
IF(OUTVAR(MY).VALUE>0)  CALL GLUTADDMENUENTRY ("MY",MY)
IF(OUTVAR(MZ).VALUE>0)  CALL GLUTADDMENUENTRY ("MZ",MZ)
IF(OUTVAR(DISCHARGE).VALUE>0)  CALL GLUTADDMENUENTRY ("Q",DISCHARGE)

CONTOUR_PLOT_ID=glutCreateMenu(contour_handler)
call glutAddMenuEntry("toggle contour surface",Contour_surfsolid_toggle)
call glutAddSubMenu("contour surface color",surface_color_menu)
call glutAddMenuEntry("toggle contour line",Contour_Line_toggle)
call glutAddSubMenu("contour line color",contour_color_menu)
call glutAddMenuEntry("Densify contour",contour_densify)
call glutAddMenuEntry("Sparsify contour",contour_Sparsify)
call glutAddMenuEntry("Plot In DeformedGrid",Contour_In_DeformedMesh)
call glutAddMenuEntry("contour values",set_contour_val)
call glutAddSubMenu("DISPLACE",VSHOW_DIS_ID)
call glutAddSubMenu("FORCE",VSHOW_FORCE_ID)
call glutAddSubMenu("SEEPAGE",VSHOW_SPG_ID)
call glutAddSubMenu("STRESS",VSHOW_STRESS_ID)
call glutAddSubMenu("STRAIN",VSHOW_STRAIN_ID)
call glutAddSubMenu("PSTRAIN",VSHOW_PSTRAIN_ID)
call glutAddSubMenu("SFR",VSHOW_SFR_ID)


VECTOR_PAIR_ID=glutCreateMenu(SET_VECTOR_PLOT_GROUP)
IF(OUTVAR(DISX).VALUE>0.AND.OUTVAR(DISY).VALUE>0) THEN
    CALL GLUTADDMENUENTRY ("DISPLACE",VECTOR_GROUP_DIS)   
ENDIF
IF(OUTVAR(VX).VALUE>0.AND.OUTVAR(VY).VALUE>0) THEN
    CALL GLUTADDMENUENTRY ("SEEPAGE_VELOCITY",VECTOR_GROUP_SEEPAGE_VEC)
ENDIF
IF(OUTVAR(GRADX).VALUE>0.AND.OUTVAR(GRADY).VALUE>0) THEN
    CALL GLUTADDMENUENTRY ("SEEPAGE_GRADIENT",VECTOR_GROUP_SEEPAGE_GRAD)
ENDIF
IF(OUTVAR(SFR_SFRX).VALUE>0.AND.OUTVAR(SFR_SFRY).VALUE>0) THEN
    CALL GLUTADDMENUENTRY ("STRESS_FAILUE_FACE",VECTOR_GROUP_SFR)
ENDIF

VECTOR_PLOT_ID=glutCreateMenu(VECTOR_HANDLER)
call glutAddMenuEntry("toggle Vector Plot",Vector_toggle)
call glutAddMenuEntry("Lengthen Vector",Vector_Lengthen)
call glutAddMenuEntry("Shorten Vector",Vector_Shorten)
call glutAddSubMenu("VectorPair",VECTOR_PAIR_ID)

Model_ID=glutCreateMenu(Model_HANDLER)
call glutAddMenuEntry("Grid Toggle",surfgrid_toggle)
call glutAddMenuEntry("DeformedGrid Toggle",DeformedMesh)
call glutAddMenuEntry("++DeformedGridScale",Enlarge_Scale_DeformedMesh)
call glutAddMenuEntry("--DeformedGridScale",Minify_Scale_DeformedMesh)


menuid = glutCreateMenu(menu_handler)
call glutAddSubMenu("Contour",CONTOUR_PLOT_ID)
call glutAddSubMenu("Vector",VECTOR_PLOT_ID)
call glutAddSubMenu("Model",Model_ID)
call glutAddSubMenu("View",submenuid)
!call glutAddSubMenu("plotting parameters",param_id)
call glutAddMenuEntry("quit",quit_selected)

call glutAttachMenu(GLUT_RIGHT_BUTTON)
end subroutine make_menu

end module function_plotter



!---------------------------------------------------------------------------

subroutine plot_func

use opengl_gl
use opengl_glut
use view_modifier
use function_plotter
use solverds
implicit none

integer :: winid, menuid, submenuid
!real(gldouble)::r1
interface
    subroutine myreshape(w,h)
        integer::w,h
    end subroutine    
endinterface

! Initializations

call glutInit
call glutInitDisplayMode(ior(GLUT_DOUBLE,ior(GLUT_RGB,GLUT_DEPTH)))
call glutInitWindowSize(800_glcint,600_glcint)

! Create a window



winid = glutCreateWindow(title)

minx=minval(node.coord(1));maxx=maxval(node.coord(1))
miny=minval(node.coord(2));maxy=maxval(node.coord(2))
minz=minval(node.coord(3));maxz=maxval(node.coord(3))

modelr=((minx-maxx)**2+(miny-maxy)**2+(minz-maxz)**2)**0.5/2.0
model_radius=modelr
init_lookat.x=(minx+maxx)/2.0
init_lookat.y=(miny+maxy)/2.0
init_lookat.z=(minz+maxz)/2.0
if(ndimension<3) then
init_lookfrom.x=init_lookat.x
init_lookfrom.y=init_lookat.y
init_lookfrom.z=3.0*modelr+init_lookat.z
else
init_lookfrom.x=modelr*3.0
init_lookfrom.y=modelr*3.0
init_lookfrom.z=modelr*3.0
endif


CONTOUR_PLOT_VARIABLE=MAX(HEAD,DISZ,DISY,DISX,LOCZ)
call initialize_contourplot()


! initialize view_modifier, receiving the id for it's submenu

submenuid = view_modifier_init()

! create the menu

call make_menu(submenuid)

! Set the display callback

call glutDisplayFunc(display)


call initGL(modelr)
! Create the image

call DrawSurfaceContour()
call DrawLineContour() 
call drawvector()
call drawgrid()

! Let glut take over

call glutMainLoop

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
    
    lightPos(:)= [0., r*3., r*2., 1.] 
    call glLightfv(GL_LIGHT0, GL_AMBIENT, lightKa);
    call glLightfv(GL_LIGHT0, GL_DIFFUSE, lightKd);
    call glLightfv(GL_LIGHT0, GL_SPECULAR, lightKs);


    call glLightfv(GL_LIGHT0, GL_POSITION, lightPos);

    call glEnable(GL_LIGHT0);                        !// MUST enable each light source after configuration
end subroutine


