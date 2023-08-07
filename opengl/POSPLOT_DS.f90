!---------------------------------------------------------------------------
    
module function_plotter
use opengl_gl
use opengl_glut
use view_modifier
use MESHGEO
use COLORMAP
use common_para_opengl
implicit none
!private
!public :: display,menu_handler,make_menu,CONTOUR_PLOT_VARIABLE,VECTOR_PLOT_GROUP
!REAL(8),PARAMETER::PI=3.141592653589793


private::SET_VARIABLE_SHOW,SET_VECTOR_PLOT_GROUP
! symbolic constants

!contour
integer,parameter:: Contour_surfsolid_toggle = 1, &
                    Contour_Line_toggle = 2, &
                    Contour_Densify=3,&
                    Contour_Sparsify=4,&
                    Contour_In_DeformedMesh=5,&
                    Coutour_MinMax=6
integer, parameter :: black_contour = 1, &
                      rainbow_contour = 2
integer, parameter :: white_surface = 1, &
                      red_surface = 2, &
                      rainbow_surface = 3 ,&
                      transparency=4,&
                      COLOR_BY_SET=5
                      
real(8),public::Scale_Contour_Num=1.0
logical,public::istransparency=.false.,IsContour_In_DeformedMesh=.false. 

type contour_bar_tydef
	INTEGER::IVARPLOT=0
    integer :: nfrac,nval
    real(GLDOUBLE),allocatable :: val(:)
    real(GLFLOAT),allocatable::COLOR(:,:)
    character(128)::TITLE
endtype
type(contour_bar_tydef)::ContourBar

TYPE ISOLINE_TYDEF
    INTEGER::NV=0,NE=0,NTRI=0
    REAL(8)::VAL
    REAL(8),ALLOCATABLE::V(:,:)
    INTEGER,ALLOCATABLE::EDGE(:,:),TRI(:,:)    
END TYPE
TYPE(ISOLINE_TYDEF),ALLOCATABLE::ISOLINE(:)
INTEGER::NISOLINE

TYPE SLICE_TYDEF
	INTEGER::IVARPLOT=0
    INTEGER::PLANE=-1
    INTEGER::NV=0,NTRI=0,NISOLINE=0,NBCE=0 
    REAL(8)::X=0
    REAL(8),ALLOCATABLE::V(:,:),VAL(:,:)
    INTEGER,ALLOCATABLE::TRI(:,:),BCEDGE(:,:)
    TYPE(ISOLINE_TYDEF),ALLOCATABLE::ISOLINE(:)
END TYPE
TYPE(SLICE_TYDEF)::SLICE(30)
INTEGER::NSLICE=0

INTEGER::IVO(3)=0 !FOR VECTOR PAIR LOCATION IN POSDATA.NODALQ
INTEGER::IEL_STREAMLINE=0 !当前积分点所在的单元
INTEGER,PARAMETER::streamline_location_click=1,Plot_streamline_CLICK=2,&
                   Reset_streamline_CLICK=3,Enlarger_StrokeFontSize_CLICK=4,&
                   Smaller_StrokeFontSize_CLICK=5,SHOW_STREAMLINE_NODE_CLICK=6,OUTPUT_SHOWNSTREAMLINE_CLICK=7
LOGICAL::isstreamlineinitialized=.false.
TYPE STREAMLINE_TYDEF
    INTEGER::NV=0,ISINPUTSLIP=0
    LOGICAL::SHOW=.TRUE.,ISLOCALSMALL=.FALSE.,ISCOMPATABLE=.TRUE.
    REAL(8)::PTstart(3),SF_SLOPE=HUGE(1.d0),SVEC(2)=0.D0 !SVEC=3D SLOPE SLIDE SURFACE DIRECTION 
    REAL(8),ALLOCATABLE::V(:,:),VAL(:,:)
    REAL(8),ALLOCATABLE::PARA_SFCAL(:,:) !SIGN,SIGT,RAD1,SIGTA,SEG_LENGTH,C,PHI,SFR,DCOS(3) !后面三个为滑动面方向矢量	
    INTEGER,ALLOCATABLE::IEL(:) !流线上每个节点所在的单元号
    REAL(8),ALLOCATABLE::SLIDE_STRIP(:,:,:) !以流线和第二主应力方向线生成的生成的滑动面（三维滑坡）,slide_strip为滑动面两侧边线的坐标,把两边线对应的节点连起来,即为一个四边形滑动面.
CONTAINS
    PROCEDURE::GET_SLIDE_STRIP=>CAL_SLIDE_STRIP
ENDTYPE
INTEGER::MAXNSTREAMLINE=0
TYPE(STREAMLINE_TYDEF),ALLOCATABLE::STREAMLINE(:)
INTEGER,ALLOCATABLE::SF_SLOPE(:)
INTEGER::NSTREAMLINE=0,NINPUTSLIP=0,INPUTSLIP(1000)=0
LOGICAL::IS_JUST_SHOW_TOP_TEN_SLOPE=.FALSE.,IS_JUST_SHOW_THE_MINIMAL_ONE_SLOPE=.FALSE.,&
        IS_SHOW_ALL_SLOPE=.TRUE.,SLOPE_CHECK_ADMISSIBILITY=.TRUE.
REAL(8),ALLOCATABLE::SLOPESURFACE(:,:)

TYPE STREAMLINENODE_IN_TET_TYDEF
    INTEGER::NNODE=0
    INTEGER,ALLOCATABLE::NODE(:) !单数,istreaamline,双数,iNODE
CONTAINS
    PROCEDURE::ENLARGE_NODE=>TET2SN_ENLARGE_NODE
ENDTYPE
TYPE(STREAMLINENODE_IN_TET_TYDEF),ALLOCATABLE::TET2SN(:)


!vector
TYPE VECTOR_PLOT_TYDEF
	INTEGER::GROUP
	REAL(8)::SCALE=1.0D0
	REAL(8),POINTER::VEC(:,:)=>NULL()
ENDTYPE

INTEGER,PARAMETER:: NVECTORPAIR=10,&
					VECTOR_GROUP_DIS=1,&
                    VECTOR_GROUP_SEEPAGE_VEC=2,&
                    VECTOR_GROUP_SEEPAGE_GRAD=3,&
                    VECTOR_GROUP_SFR=4,&
					VECTOR_GROUP_PSIGMA1=5,&
					VECTOR_GROUP_PSIGMA2=6,&
                    VECTOR_GROUP_PSIGMA3=7
integer,parameter::Vector_toggle= 1, &
                   Vector_lengthen=2,&
                   Vector_shorten=3,&
                   Vector_Unify=4,&
                   Vector_Less=5,&
                   Vector_More=6,Arrow_toggle=7,Vector_Base_Orig=8,Vector_Base_End=9,Vector_Base_Cent=10
                   
real(8),public::Scale_Vector_Len=1.0,VabsMax,VabsMin,Vscale                    
character(128)::VectorPairName(NVECTORPAIR)=['DIS.','SEEP.V','SEEP.I','SFR','PSIGMA1','PSIGMA2','PSIGMA3','','','']
INTEGER,PARAMETER::VECTORLOC_NODE=1,VECTORLOC_GRID=2
integer::VectorFrequency=1,VectorBase=1,VECTORLOC=VECTORLOC_GRID !AT STRUTRAL GRID.
logical::IsVectorUnify=.False.,IsArrow=.True.
!model
integer, parameter :: surfgrid_toggle = 1, &
                      quit_selected = 4,&
                      DeformedMesh=5,&
                      Enlarge_Scale_DeformedMesh=6,&
                      Minify_Scale_DeformedMesh=7,&
					  Edge_toggle=8,&
					  Node_toggle=9,&
					  probe_selected=10,&
                      SHOW_SET_TOGGLE=11,&
                      Probe_Get_data_click=12
					  
					  
!NODALVALUE
INTEGER,PARAMETER::SHOWNODALVALUE_TOGGLE=1
INTEGER::SHOW_NODAL_VALUE=0
LOGICAL::ISSHOWNODALVALUE=.FALSE.,isProbestate=.false.
                      
logical::IsDeformedMesh=.false.,show_edge=.false.,show_node=.false.,SHOW_SET=.FALSE.                     
real(8)::Scale_Deformed_Grid=1.d0


!SLICE
INTEGER,PARAMETER::SLICE_LOCATION_CLICK=1,PLOTSLICESURFACE_CLICK=2,PLOTSLICEISOLINE_CLICK=3,&
				PLOTSLICE_CLICK=4
	

!TYPE SLICE_TYDEF
!    
!
!ENDTYPE

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
           IsDrawVector=.false.,&
		   IsPlotSliceSurface=.false.,&
		   isPLotSliceIsoLIne=.false.,&
           isPlotStreamLine=.false.,&
		   ISPLOTSLICE=.FALSE.,&
           SHOW_STREAMLINE_NODE=.false.

real(GLDOUBLE), allocatable :: actual_contours(:)
real(GLDOUBLE) :: minv,maxv

integer::CONTOUR_PLOT_VARIABLE=0,VECTOR_PLOT_GROUP=VECTOR_GROUP_DIS,SLICE_PLOT_VARIABLE=0,STREAMLINE_VECTOR=VECTOR_GROUP_SEEPAGE_VEC
LOGICAL::ACTIVE_VECTOR_GROUP(NVECTORPAIR)=.FALSE.

integer,parameter,public::ContourList=1,&
                          VectorList=2,&
                          GridList=3,&
						  ContourLineList=4,&
                          ProbeValueList=5,&
						  SLICELIST=6,&
                          StreamLineList=7,&
						  STEPSTATUSLIST=8,&
                          SFS_ProbeValuelist=9 !,&
                          !PSO_SLIP_PLOT_LIST=10
                          
                          
integer,parameter::STREAMLINE_SLOPE_CLICK=1,ShowTopTen_SLOPE_CLICK=2,ShowMinimal_SLOPE_CLICK=3,&
ShowAll_SLOPE_CLICK=4,ReadSlipSurface_slope_click=5,ShowLocalMinimal_Slope_Click=6,&
                    SEARCH_SLOPE_CLICK=7,FilterLocalMinimal_Slope_Click=8,&
                    SFS_Slope_Click=9,CHECKADMISSIBILITY=10,ShowReadinSlip_Slope_CLICK=11,SEARCH_BY_PSO=12,&
                    ShowStreamline_AdmissiblityImprove_CLICK=13,ShowSFContour_click=14,ShowSSTRIP_click=15,& 
                    ShowFSLabel_click=16,OutSlideInfo_click=17
LOGICAL::ISSTREAMLINESLOPE=.FALSE.,IsFilterLocalMinimalMode=.false.,isProbeState_SFS=.FALSE.,&
         IsShowReadinSlip=.false.,IsShowSlideStrip=.false.,IsShowFSLabel=.false.

!real(GLFLOAT) :: red(4) = (/1.0,0.0,0.0,1.0/), &
!                 black(4) = (/0.0,0.0,0.0,1.0/), &
!                 white(4) = (/1.0,1.0,1.0,1.0/),&
!                 GRAY(4)=(/0.82745098,0.82745098,0.82745098,1.0/) 


TYPE TIMESTEPINFO_TYDEF
    INTEGER::ISTEP=1,NSTEP=1
    REAL(8),ALLOCATABLE::TIME(:)
	!INTEGER,ALLOCATABLE::CALSTEP(:) !当前步对应的计算步,注意,计算步与绘图步往往不一致。
    LOGICAL::ISSHOWN=.TRUE.
	REAL(8)::VSCALE(NVECTORPAIR)=1.0D0,VMIN(NVECTORPAIR),VMAX(NVECTORPAIR)
    CHARACTER(256)::INFO=''
    CONTAINS
    PROCEDURE::INITIALIZE=>STEP_INITIALIZE
    PROCEDURE::UPDATE=>STEP_UPDATE
	
ENDTYPE
TYPE(TIMESTEPINFO_TYDEF),PUBLIC::STEPPLOT

REAL(8)::PT_PROBE(3)=0.0D0,VAL_PROBE(150),STRIP_WIDTH=0.5
INTEGER::IEL_PROBE=0


INTEGER::NODALID_SHOW=-1,ELEMENTID_SHOW=-2



!TYPE COLORMAP_TYPDEF
!    CHARACTER(128)::NAME=''
!    INTEGER::IMAP,NCP=0
!    REAL(8),ALLOCATABLE::COTROL_COLOR(:) !R,G,B,VAL
!CONTAINS
!    !PROCEDURE::
!END TYPE

    contains

    SUBROUTINE CAL_SLIDE_STRIP(this,WIDTH)
    !计算3D边坡滑动带
    !滑动面由流线方向和第二主应力方向确定(近似)
    !
        USE SolverMath,ONLY:CS_VECTOR
        IMPLICIT NONE
        CLASS(STREAMLINE_TYDEF)::THIS
        REAL(8),OPTIONAL::WIDTH
        REAL(8)::DIR1(3),DX1(3),W1,T1,DIR2(3),w2,DX2(3),XC1,YC1
        INTEGER::I,J,N1,N2,N3

        IF(POSDATA.NDIM/=3) RETURN
        W1=0.5
        IF(PRESENT(WIDTH)) W1=WIDTH/2.0

        IF(ALLOCATED(THIS.SLIDE_STRIP)) DEALLOCATE(THIS.SLIDE_STRIP) 
        ALLOCATE(THIS.SLIDE_STRIP(3,THIS.nv,2))
        ! DO I=1,THIS.NV

        !     IF(I<THIS.NV) THEN
        !         DX1=THIS.V(:,I+1)-THIS.V(:,I)
        !         DX1=DX1/NORM2(DX1)
        !         !DIR2:第二主应力方向
        !         DIR2=CS_VECTOR(THIS.PARA_SFCAL(9:11,I),DX1,.true.)                
        !     ENDIF
        !     !dir1相邻滑动面交线向量
        !     IF(I==1.OR.I==THIS.NV) THEN 
        !         !假定第一个交线为z平面与第一个滑动面的交线
        !         DIR1=CS_VECTOR([0.,0.,1.0],THIS.PARA_SFCAL(9:11,I),.true.) 
        !     ELSE
        !         T1=DOT_PRODUCT(THIS.PARA_SFCAL(9:11,I),THIS.PARA_SFCAL(9:11,I-1))
        !         IF(ABS(ABS(T1)-1.D0)>0.001) THEN
        !             DIR1=CS_VECTOR(THIS.PARA_SFCAL(9:11,I),THIS.PARA_SFCAL(9:11,I-1),.true.)  
        !         ELSE
        !             DIR1=DIR2
        !         ENDIF              
        !     ENDIF

        !     T1=DOT_PRODUCT(DIR2,DIR1) !令条带的宽度近似相等
        !     W2=W1/ABS(T1)
        !     DX2=DIR1*W2

        !     !下面的运算是为了保证四边形的对边不相交
        !     IF(I>1) THEN
        !         !相邻对边的方向的方向应该一致
        !         DIR1=THIS.SLIDE_STRIP(:,I-1,1)-THIS.V(:,I-1)
        !         T1=DOT_PRODUCT(DIR1,DX2)
        !         IF(T1>0.D0) THEN
        !             T1=1.0D0
        !         ELSE
        !             T1=-1.D0
        !         ENDIF
        !     ENDIF
        !     THIS.SLIDE_STRIP(:,I,1)=THIS.V(:,I)+DX2*T1
        !     THIS.SLIDE_STRIP(:,I,2)=THIS.V(:,I)-DX2*T1


        ! ENDDO

        DO I=1,THIS.NV

            DX2=THIS.VAL([POSDATA.IXPS2,POSDATA.IYPS2,POSDATA.IZPS2],I)
            DX2=DX2/NORM2(DX2)

            !流线滑槽与坡表面的交线
            IF(I==1.OR.I==THIS.NV) THEN 
                !假定第一个交线为z平面与第一个滑动面的交线
                N1=this.iel(i)
                N2=COUNT(tet(n1).adjelt(1:TET(N1).NF)==0)
                N3=0
                DO J=1,tet(N1).NF
                    IF(tet(n1).adjelt(j)==0) THEN
                        IF(N2>1)THEN
                            IF(PtInTri (THIS.V(:,I), POSDATA.NODE(FACE(TET(N1).F(J)).V(1)).COORD, &
                                POSDATA.NODE(FACE(TET(N1).F(J)).V(2)).COORD, &
                                POSDATA.NODE(FACE(TET(N1).F(J)).V(3)).COORD)) THEN
                                N3=J
                                EXIT
                            ENDIF
                        ELSE
                            N3=J
                            EXIT
                        ENDIF
                    ENDIF
                ENDDO
                IF(N3>0) THEN
                    DIR1=CS_VECTOR(POSDATA.NODE(FACE(TET(N1).F(N3)).V(2)).COORD-POSDATA.NODE(FACE(TET(N1).F(N3)).V(1)).COORD,&
                    POSDATA.NODE(FACE(TET(N1).F(N3)).V(3)).COORD-POSDATA.NODE(FACE(TET(N1).F(N3)).V(1)).COORD,.true.)
                    T1=DOT_PRODUCT(DIR1,THIS.PARA_SFCAL(9:11,I))
                    IF(ABS(ABS(T1)-1.D0)>0.001) THEN
                        DX2=CS_VECTOR(DIR1,THIS.PARA_SFCAL(9:11,I),.true.)  
                    ENDIF
                ENDIF

            ENDIF 

            !下面的运算是为了保证四边形的对边不相交
            T1=W1
            IF(I>1) THEN
                !相邻对边的方向应该一致
                DIR1=THIS.SLIDE_STRIP(:,I-1,1)-THIS.SLIDE_STRIP(:,I-1,2)
                T1=DOT_PRODUCT(DIR1,DX2)
                IF(T1>0.D0) THEN
                    T1=W1
                ELSE
                    T1=-W1
                ENDIF
            ENDIF

            THIS.SLIDE_STRIP(:,I,1)=THIS.V(:,I)+DX2*T1
            THIS.SLIDE_STRIP(:,I,2)=THIS.V(:,I)-DX2*T1


        ENDDO

        XC1=(POSDATA.MINX+POSDATA.MAXX)/2
        YC1=(POSDATA.MINY+POSDATA.MAXY)/2
    
    
        !IF((.NOT.STREAMLINE(I).ISCOMPATABLE).OR.STREAMLINE(I).NV<2) CYCLE
        
        N1=1
            
        !reference vector
        DX1(1)=THIS.V(1,N1)-XC1
        DX1(2)=THIS.V(2,N1)-YC1
        T1=ATAN2(DX1(2),DX1(1))+1.57079632679
        DX1(1)=COS(T1);DX1(2)=SIN(T1)
        DX2=THIS.SLIDE_STRIP(:,N1,2)-THIS.SLIDE_STRIP(:,N1,1) 
        T1=NORM2(DX2(1:2))
        DX2=DX2/T1
        T1=DOT_PRODUCT(DX1(1:2),DX2(1:2))
        IF(T1<0.D0)  DX2=-DX2
        THIS.SVEC=DX2(1:2)


    ENDSUBROUTINE

    subroutine drawStrokeText(angle,x,y,z,scale,s,rotate_axis)
        use opengl_gl
        use opengl_glut
        
        implicit none
        
        real(GLDOUBLE),INTENT(in) :: ANGLE,x,y,z,scale
        real(GLDOUBLE),INTENT(in),optional :: rotate_axis(3)
        character,intent(in) :: s*(*)
        character :: c
        integer :: i,lenc
        real(GLDOUBLE)::raxis1(3)
        
        call glPushMatrix();
        !call glLoadIdentity()
        call glTranslated(x, y,z);
        if(present(rotate_axis)) then
            raxis1=rotate_axis
        else
            raxis1=[0.,0.,1.]
        endif
        call glrotated(ANGLE,raxis1(1),raxis1(2),raxis1(3))
            
        call glScaled(scale, scale, scale);
        lenc = len(s)
        do i=1,lenc
            c = s(i:i)
            call glutStrokeCharacter(GLUT_STROKE_ROMAN, &
                ichar(c))
        end do
        call glPopMatrix();
        
        call glutPostRedisplay
    end subroutine drawStrokeText

    SUBROUTINE TET2SN_ENLARGE_NODE(THIS,DSTEP)
        CLASS(STREAMLINENODE_IN_TET_TYDEF)::THIS
        INTEGER,INTENT(IN)::DSTEP
        INTEGER,ALLOCATABLE::VAL1(:)

        !LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1(2*dstep))
        THIS.NODE=[THIS.NODE,VAL1]
        DEALLOCATE(VAL1)
    END SUBROUTINE


SUBROUTINE SEARCH_MINIMAL_SF_SLOPE(IS_ONLY_SEARCHTOP)
    !use function_plotter
    use FindLocalEx1D
    implicit none
    INTEGER,INTENT(IN),OPTIONAL::IS_ONLY_SEARCHTOP
    INTEGER::I,J,K,IA1(NNLBC),NH1,IA2(NNLBC),IMSF1,N1,NSTREAMLINE1,IMSF2,N2
    LOGICAL::L1,L2
    REAL(8)::MSF1,PT1(3),MSF2,ASF1(NNLBC),PT2(3),MAXY1
    REAL(8),ALLOCATABLE::MAXTEM1(:,:),MINTEM1(:,:)
    character(16)::str1
    
    NH1=0;IA1=0
    MSF1=1E10
    isPlotStreamLine=.TRUE.
    STREAMLINE_VECTOR=VECTOR_GROUP_SEEPAGE_VEC
    NSTREAMLINE1=NSTREAMLINE+1
    ISSTREAMLINESLOPE=.TRUE.
   
    info.color=red;info.qkey=.true.
    MAXY1=-1.D20
    IF(PRESENT(IS_ONLY_SEARCHTOP)) THEN
        IF(IS_ONLY_SEARCHTOP/=0) MAXY1=MAXVAL(POSDATA.NODE(NODE_LOOP_BC(:NNLBC)).COORD(POSDATA.NDIM))
    ENDIF

    POSDATA.NODALQ(:,POSDATA.IFOS,stepplot.istep)=-1.0
    DO I=1,NNLBC 
        IF(POSDATA.IH_BC>0) THEN
            L1=(POSDATA.NODALQ(NODE_LOOP_BC(I),POSDATA.IH_BC,stepplot.istep)+4.0D0)<1E-7
            IF(POSDATA.NDIM==2) THEN
                L2=POSDATA.NODALQ(NODE_LOOP_BC(I),POSDATA.IQ,stepplot.istep)>0.d0   
            ELSE
                L2=.TRUE.
            ENDIF         
            IF(L1.AND.L2) THEN                
             
                PT1(3)=0.0D0
                PT1(1:POSDATA.NDIM)=POSDATA.NODE(NODE_LOOP_BC(I)).COORD(1:POSDATA.NDIM)
                IF(PT1(POSDATA.NDIM)<MAXY1) CYCLE
                !write(*,100) pt1(1:2)
                CALL gen_new_streamline(PT1)
                IF(STREAMLINE(NSTREAMLINE).ISCOMPATABLE) POSDATA.NODALQ(NODE_LOOP_BC(I),POSDATA.IFOS,stepplot.istep)=STREAMLINE(NSTREAMLINE).SF_SLOPE
                NH1=NH1+1
                IA1(NH1)=I
                IA2(NH1)=NSTREAMLINE
                ASF1(NH1)=STREAMLINE(NSTREAMLINE).SF_SLOPE
                STREAMLINE(NSTREAMLINE).SHOW=.FALSE.
                IF(STREAMLINE(NSTREAMLINE).SF_SLOPE<MSF1) THEN
                    MSF1=STREAMLINE(NSTREAMLINE).SF_SLOPE
                    IMSF1=NSTREAMLINE
                    N2=I                    
                ENDIF
            ENDIF
            
        ELSE
            STOP 'THERE SEEMS TO BE NO SUCH A VARIABLE "H_BC".'
        ENDIF
    
    ENDDO
    
    !REFINE
    
    MSF2=MSF1;IMSF2=IMSF1
    IF(POSDATA.NDIM==2) THEN
        DO K=1,2
            IF(K==1) THEN
                N1=N2-1
                IF(N1<1) N1=NNLBC
            ELSE
                N1=MOD(N2,NNLBC)+1   
            ENDIF
            L1=(POSDATA.NODALQ(NODE_LOOP_BC(N1),POSDATA.IH_BC,stepplot.istep)+4.0D0)<1E-7
            L2=POSDATA.NODALQ(NODE_LOOP_BC(N1),POSDATA.IQ,stepplot.istep)>0.d0
                
            IF(L1.AND.L2) THEN        
                PT2=POSDATA.NODE(NODE_LOOP_BC(N1)).COORD
            ELSE
                PT2(1)=POSDATA.NODE(NODE_LOOP_BC(N1)).COORD(1)
                PT2(2:3)=STREAMLINE(IMSF1).V(2:3,MIN(10,STREAMLINE(IMSF1).NV))            
            ENDIF
            
            DO J=2,10
                PT1=POSDATA.NODE(NODE_LOOP_BC(N2)).COORD+(J-1)/10.0*(PT2-POSDATA.NODE(NODE_LOOP_BC(N2)).COORD)
                call gen_new_streamline(PT1)
                STREAMLINE(NSTREAMLINE).SHOW=.FALSE.
                IF(STREAMLINE(NSTREAMLINE).SF_SLOPE<MSF2) THEN
                    MSF2=STREAMLINE(NSTREAMLINE).SF_SLOPE
                    IMSF2=NSTREAMLINE                  
                ENDIF
            ENDDO 
        
        ENDDO
    ELSE
        !TO BE IMPROVED  
    ENDIF
    
    STREAMLINE(IMSF2).SHOW=.TRUE.
    write(str1,'(f8.3)') msf2
    info.str='Search done...The MinimalFS='//trim(adjustl(str1))
    info.qkey=.true.
    !call peakdet(MAXTEM1,MINTEM1,NH1,ASF1(1:NH1),0.01)
    !
    !DO I=1,SIZE(MINTEM1,DIM=1)
    !    STREAMLINE(IA2(INT(MINTEM1(1,I)))).ISLOCALSMALL=.TRUE.    
    !ENDDO
    
    RETURN
    
100 format('To gen_new_streamline at point=',2F8.3)
110 format('Done in gen_new_streamline at point=',2F8.3)
ENDSUBROUTINE    

SUBROUTINE SET_SAFETY_FIELD()
    IMPLICIT NONE
    INTEGER::I,J,N1

    IF(.NOT.ALLOCATED(TET2SN)) ALLOCATE(TET2SN(NTET))

    DO I=1,NTET
        IF(ALLOCATED(TET2SN(I).NODE)) DEALLOCATE(TET2SN(I).NODE)
        ALLOCATE(TET2SN(I).NODE(10))
    ENDDO

    DO I=1,NSTREAMLINE
        IF(.NOT.STREAMLINE(I).ISCOMPATABLE) CYCLE
        DO J=1,STREAMLINE(I).NV
            N1=TET2SN(STREAMLINE(I).IEL(J)).NNODE            
            IF(2*TET2SN(STREAMLINE(I).IEL(J)).NNODE>SIZE(TET2SN(STREAMLINE(I).IEL(J)).NODE)) THEN
                CALL TET2SN(STREAMLINE(I).IEL(J)).ENLARGE_NODE(5)
            ENDIF
            TET2SN(STREAMLINE(I).IEL(J)).NODE(N1+1)=I
            TET2SN(STREAMLINE(I).IEL(J)).NODE(N1+2)=J
            TET2SN(STREAMLINE(I).IEL(J)).NNODE=TET2SN(STREAMLINE(I).IEL(J)).NNODE+1
        ENDDO
    ENDDO

    


ENDSUBROUTINE

SUBROUTINE getHeatMapColor(value,VMIN,VMAX,VCOLOR,BCOLORI)
!!http://andrewnoske.com/wiki/Code_-_heatmaps_and_color_gradients
  !int aR = 0;   int aG = 0; int aB=255;  ! RGB for our 1st color (blue in this case).
  !int bR = 255; int bG = 0; int bB=0;    ! RGB for our 2nd color (red in this case).
  !
  !red   = (float)(bR - aR) * value + aR;      ! Evaluated as -255*value + 255.
  !green = (float)(bG - aG) * value + aG;      ! Evaluates as 0.
  !blue  = (float)(bB - aB) * value + aB;      ! Evaluates as 255*value + 0.   REAL(8),INTENT(IN)::VALUE,VMIN,VMAX
    REAL(8),INTENT(IN)::VALUE,VMIN,VMAX   
    REAL(GLFLOAT),INTENT(OUT)::VCOLOR(3)
    INTEGER,INTENT(IN),OPTIONAL::BCOLORI(:)
    INTEGER::NBCOLOR
    INTEGER,ALLOCATABLE::DEFAULTBCOLORI(:)
    INTEGER:: IDX1,IDX2,I
    REAL(8)::FRACTBETWEEN=0,V1,T1
   
   
    IF(PRESENT(BCOLORI)) THEN
        ALLOCATE(DEFAULTBCOLORI,SOURCE=BCOLORI)
        NBCOLOR=SIZE(BCOLORI)
    ELSE
        ALLOCATE(DEFAULTBCOLORI(4))
        NBCOLOR=5
        DEFAULTBCOLORI=[blue,  Cyan, green,  yellow,  red]
    ENDIF
   
    IF(ABS(VMAX-VMIN)>1E-8) THEN
        V1=(VALUE-VMIN)/(VMAX-VMIN)
    ELSE
        V1=0.5D0
    ENDIF
    
    IF(V1<=0.D0) THEN
        IDX1=1;IDX2=1
    ELSEIF(V1>=1.D0) THEN
        IDX1=NBCOLOR;IDX2=IDX1
    ELSE
 !    value = value * (NUM_COLORS-1);        // Will multiply value by 3.
!    idx1  = floor(value);                  // Our desired color will be after this index.
!    idx2  = idx1+1;                        // ... and before this index (inclusive).
!    fractBetween = value - float(idx1);    // Distance between the two indexes (0-1).
        V1=V1*(NBCOLOR-1)+1
        IDX1=FLOOR(V1)
        IDX2=IDX1+1
        FRACTBETWEEN=(V1-REAL(IDX1))
    ENDIF

!  *red   = (color[idx2][0] - color[idx1][0])*fractBetween + color[idx1][0];
!  *green = (color[idx2][1] - color[idx1][1])*fractBetween + color[idx1][1];
!  *blue  = (color[idx2][2] - color[idx1][2])*fractBetween + color[idx1][2];
    DO I=1,3
        VCOLOR(I)=(MYCOLOR(I,DEFAULTBCOLORI(IDX2))-MYCOLOR(I,DEFAULTBCOLORI(IDX1)))*FRACTBETWEEN+MYCOLOR(I,DEFAULTBCOLORI(IDX1))
    ENDDO
   
   

ENDSUBROUTINE

!bool getHeatMapColor(float value, float *red, float *green, float *blue)
!{
!  const int NUM_COLORS = 4;
!  static float color[NUM_COLORS][3] = { {0,0,1}, {0,1,0}, {1,1,0}, {1,0,0} };
!    // A static array of 4 colors:  (blue,   green,  yellow,  red) using {r,g,b} for each.
! 
!  int idx1;        // |-- Our desired color will be between these two indexes in "color".
!  int idx2;        // |
!  float fractBetween = 0;  // Fraction between "idx1" and "idx2" where our value is.
! 
!  if(value <= 0)      {  idx1 = idx2 = 0;            }    // accounts for an input <=0
!  else if(value >= 1)  {  idx1 = idx2 = NUM_COLORS-1; }    // accounts for an input >=0
!  else
!  {
!    value = value * (NUM_COLORS-1);        // Will multiply value by 3.
!    idx1  = floor(value);                  // Our desired color will be after this index.
!    idx2  = idx1+1;                        // ... and before this index (inclusive).
!    fractBetween = value - float(idx1);    // Distance between the two indexes (0-1).
!  }
! 
!  *red   = (color[idx2][0] - color[idx1][0])*fractBetween + color[idx1][0];
!  *green = (color[idx2][1] - color[idx1][1])*fractBetween + color[idx1][1];
!  *blue  = (color[idx2][2] - color[idx1][2])*fractBetween + color[idx1][2];
!}


!ENLARGE STREAMLINE SF_slope by 100
SUBROUTINE ENLARGE_STREAMLINE()
    IMPLICIT NONE
    
    TYPE(STREAMLINE_TYDEF),ALLOCATABLE::STREAMLINE1(:)
	INTEGER,ALLOCATABLE::SF_SLOPE1(:)
    
	MAXNSTREAMLINE=MAXNSTREAMLINE+100
	Allocate(STREAMLINE1(MAXNSTREAMLINE),SF_SLOPE1(MAXNSTREAMLINE))
    
    IF(NSTREAMLINE>0) THEN
        STREAMLINE1(1:NSTREAMLINE)=STREAMLINE(1:NSTREAMLINE)
        SF_SLOPE1(1:NSTREAMLINE)=SF_SLOPE(1:NSTREAMLINE)
    ENDIF
	if(allocated(STREAMLINE)) deallocate(STREAMLINE,SF_SLOPE)
		
	allocate(STREAMLINE,SOURCE=STREAMLINE1)
    allocate(SF_SLOPE,SOURCE=SF_SLOPE1)
	deallocate(STREAMLINE1,SF_SLOPE1)    

ENDSUBROUTINE

SUBROUTINE VECTORGRIDGEN()
	IMPLICIT NONE
	 
	INTEGER,PARAMETER::NDIV1=100
	REAL(8)::DX1=0.0D0,X1(NDIV1+1),Y1(NDIV1+1),Z1(NDIV1+1),T1=0.D0
	INTEGER::I,J,K,NX1=0,NY1=0,NZ1=0,N1=0,IEL1=0
	

	IF(ALLOCATED(POSDATA.GRIDDATA)) DEALLOCATE(POSDATA.GRIDDATA)
	IF(ALLOCATED(POSDATA.GRIDNODE)) DEALLOCATE(POSDATA.GRIDNODE)
	IF(ALLOCATED(POSDATA.GRIDVEC)) DEALLOCATE(POSDATA.GRIDVEC)
	DX1=MAXVAL([POSDATA.maxx-POSDATA.minx,POSDATA.maxy-POSDATA.miny,POSDATA.maxz-POSDATA.minz])/NDIV1
	DO I=1,NDIV1+1
		T1=POSDATA.minx+DX1*(I-1)
		IF(T1<=POSDATA.MAXX) THEN
			X1(I)=T1
			NX1=I
		ENDIF
		T1=POSDATA.minY+DX1*(I-1)
		IF(T1<=POSDATA.MAXY) THEN
			Y1(I)=T1
			NY1=I
		ENDIF		
		T1=POSDATA.minZ+DX1*(I-1)
		IF(T1<=POSDATA.MAXZ) THEN
			Z1(I)=T1
			NZ1=I
		ENDIF		
	ENDDO
	NX1=MAX(NX1,1);NY1=MAX(NY1,1);NZ1=MAX(NZ1,1);
	POSDATA.NGRIDNODE=NX1*NY1*NZ1
	ALLOCATE(POSDATA.GRIDVEC(3,POSDATA.NGRIDNODE,NVECTORPAIR), &
			 POSDATA.GRIDNODE(POSDATA.NGRIDNODE), &
			 POSDATA.GRIDDATA(POSDATA.NVAR,POSDATA.NGRIDNODE))
    N1=0
	DO I=1,NZ1
		DO J=1,NY1
			DO K=1,NX1				
				N1=N1+1
				POSDATA.GRIDNODE(N1).COORD=[X1(K),Y1(J),Z1(I)]				
			ENDDO
		ENDDO
	ENDDO

	POSDATA.GRIDDATA=0.0d0
	POSDATA.GRIDVEC=0.0d0
	!POSDATA.VEC=0.D0
	
ENDSUBROUTINE

SUBROUTINE STEP_INITIALIZE(STEPINFO,ISTEP,NSTEP,TIME)
    IMPLICIT NONE    
    CLASS(TIMESTEPINFO_TYDEF),INTENT(in out):: STEPINFO
    INTEGER,INTENT(IN)::ISTEP,NSTEP
    REAL(8),INTENT(IN)::TIME(NSTEP)
	REAL(8),ALLOCATABLE::VEC1(:,:,:)
	INTEGER::I
    
    
    STEPINFO.ISTEP=ISTEP
    STEPINFO.NSTEP=NSTEP
    ALLOCATE(STEPINFO.TIME,SOURCE=TIME)
	!ALLOCATE(STEPINFO.CALSTEP,SOURCE=CALSTEP)
	ALLOCATE(VEC1(3,POSDATA.NNODE,NSTEP))
	
    IF(ALLOCATED(POSDATA.NODE)) DEALLOCATE(POSDATA.NODE)    
    ALLOCATE(POSDATA.NODE(POSDATA.NNODE))
	
	IF(ALLOCATED(POSDATA.VEC)) DEALLOCATE(POSDATA.VEC)
	ALLOCATE(POSDATA.VEC(3,POSDATA.NNODE,NVECTORPAIR))
	POSDATA.VEC=0.D0
	
    POSDATA.NODE.COORD(3)=0.D0
    DO I=1,POSDATA.NNODE
        IF(POSDATA.IX>0) POSDATA.NODE(I).COORD(1)=POSDATA.NODALQ(I,POSDATA.IX,ISTEP)
        IF(POSDATA.IY>0) POSDATA.NODE(I).COORD(2)=POSDATA.NODALQ(I,POSDATA.IY,ISTEP)
        IF(POSDATA.IZ>0) POSDATA.NODE(I).COORD(3)=POSDATA.NODALQ(I,POSDATA.IZ,ISTEP)
    ENDDO
    
    POSDATA.minx=minval(POSDATA.NODE.COORD(1));POSDATA.maxx=maxval(POSDATA.NODE.COORD(1))
    POSDATA.miny=minval(POSDATA.NODE.COORD(2));POSDATA.maxy=maxval(POSDATA.NODE.COORD(2))
    POSDATA.minz=minval(POSDATA.NODE.COORD(3));POSDATA.maxz=maxval(POSDATA.NODE.COORD(3))
    POSDATA.MODELR=((POSDATA.minx-POSDATA.maxx)**2+(POSDATA.miny-POSDATA.maxy)**2+(POSDATA.minz-POSDATA.maxz)**2)**0.5/2.0
    !二维模型有时会输出z,这时所有的z相等,令POSDATA.NDIM=2
    
    if(abs(POSDATA.minz-POSDATA.maxz)<1.e-6) POSDATA.NDIM=2

    CALL VECTORGRIDGEN()
	
	!SET UP VECSCALE
	IF(POSDATA.IDISX>0.AND.POSDATA.IDISY>0) THEN
		DO I=1,STEPPLOT.NSTEP
			VEC1(1,:,I)=POSDATA.NODALQ(:,POSDATA.IDISX,I)
			VEC1(2,:,I)=POSDATA.NODALQ(:,POSDATA.IDISY,I)
			IF(POSDATA.NDIM>2.AND.POSDATA.IDISZ>0) THEN
				VEC1(3,:,I)=POSDATA.NODALQ(:,POSDATA.IDISZ,I)
			ELSE
				VEC1(3,:,I)=0
			ENDIF
		ENDDO
		STEPINFO.VMAX(1)=MAX(MAXVAL(NORM2(VEC1,DIM=1)),1.0E-8)
		STEPINFO.VMIN(1)=MINVAL(NORM2(VEC1,DIM=1))
		STEPINFO.VSCALE(1)=POSDATA.modelr/40./STEPINFO.VMAX(1)
        
    ENDIF
	
	IF(POSDATA.IVX>0.AND.POSDATA.IVY>0) THEN
		DO I=1,STEPPLOT.NSTEP
			VEC1(1,:,I)=POSDATA.NODALQ(:,POSDATA.IVX,I)
			VEC1(2,:,I)=POSDATA.NODALQ(:,POSDATA.IVY,I)
			IF(POSDATA.NDIM>2) THEN
				VEC1(3,:,I)=POSDATA.NODALQ(:,POSDATA.IVZ,I)
			ELSE
				VEC1(3,:,I)=0
			ENDIF
		ENDDO
		STEPINFO.VMAX(2)=MAX(MAXVAL(NORM2(VEC1,DIM=1)),1.0E-8)
		STEPINFO.VMIN(2)=MINVAL(NORM2(VEC1,DIM=1))
		STEPINFO.VSCALE(2)=POSDATA.modelr/40./STEPINFO.VMAX(2)
        
    ENDIF
    
	IF(POSDATA.IGRADX>0.AND.POSDATA.IGRADY>0) THEN
		DO I=1,STEPPLOT.NSTEP
			VEC1(1,:,I)=POSDATA.NODALQ(:,POSDATA.IGRADX,I)
			VEC1(2,:,I)=POSDATA.NODALQ(:,POSDATA.IGRADY,I)
			IF(POSDATA.NDIM>2) THEN
				VEC1(3,:,I)=POSDATA.NODALQ(:,POSDATA.IGRADZ,I)
			ELSE
				VEC1(3,:,I)=0
			ENDIF
		ENDDO
		
		STEPINFO.VMAX(3)=MAX(MAXVAL(NORM2(VEC1,DIM=1)),1.0E-8)
		STEPINFO.VMIN(3)=MINVAL(NORM2(VEC1,DIM=1))
		STEPINFO.VSCALE(3)=POSDATA.modelr/40./STEPINFO.VMAX(3)
        
    ENDIF
	
	IF(POSDATA.ISFR_SFRX>0.AND.POSDATA.ISFR_SFRY>0) THEN
		DO I=1,STEPPLOT.NSTEP
			VEC1(1,:,I)=POSDATA.NODALQ(:,POSDATA.ISFR_SFRX,I)
			VEC1(2,:,I)=POSDATA.NODALQ(:,POSDATA.ISFR_SFRY,I)
			VEC1(3,:,I)=0
		ENDDO

		STEPINFO.VMAX(4)=MAX(MAXVAL(NORM2(VEC1,DIM=1)),1.0E-8)
		STEPINFO.VMIN(4)=MINVAL(NORM2(VEC1,DIM=1))
		STEPINFO.VSCALE(4)=POSDATA.modelr/40./STEPINFO.VMAX(4)
        
    ENDIF
	
	IF(POSDATA.IPSIGMA1>0) THEN
		DO I=1,STEPPLOT.NSTEP
			VEC1(1,:,I)=POSDATA.NODALQ(:,POSDATA.IPSIGMA1,I)
			VEC1(2,:,I)=0
			VEC1(3,:,I)=0
		ENDDO

		STEPINFO.VMAX(5)=MAX(MAXVAL(NORM2(VEC1,DIM=1)),1.0E-8)
		STEPINFO.VMIN(5)=MINVAL(NORM2(VEC1,DIM=1))
		STEPINFO.VSCALE(5)=POSDATA.modelr/40./STEPINFO.VMAX(5)
		
		DO I=1,STEPPLOT.NSTEP
			VEC1(1,:,I)=POSDATA.NODALQ(:,POSDATA.IPSIGMA3,I)
			VEC1(2,:,I)=0
			VEC1(3,:,I)=0
		ENDDO

		STEPINFO.VMAX(6)=MAX(MAXVAL(NORM2(VEC1,DIM=1)),1.0E-8)
		STEPINFO.VMIN(6)=MINVAL(NORM2(VEC1,DIM=1))
		STEPINFO.VSCALE(6)=POSDATA.modelr/40./STEPINFO.VMAX(6)		
        STEPINFO.VSCALE(5:6)=maxval(STEPINFO.VSCALE(5:6))
    ENDIF
	
	
	DEALLOCATE(VEC1)
    
ENDSUBROUTINE

SUBROUTINE STEP_UPDATE(STEPINFO)
    IMPLICIT NONE
    CLASS(TIMESTEPINFO_TYDEF),INTENT(in out):: STEPINFO
    CHARACTER(16)::CWORD1,CWORD2,CWORD3,CWORD4
    REAL(8)::POS1(2)
    INTEGER::I
	
    WRITE(CWORD1,*) STEPINFO.ISTEP
    WRITE(CWORD2,*) STEPINFO.NSTEP
    WRITE(CWORD3,'(F10.5)') STEPINFO.TIME(STEPINFO.ISTEP)
    WRITE(CWORD4,'(F10.5)') STEPINFO.TIME(STEPINFO.NSTEP)
    
    STEPINFO.info='ISTEP/NSTEP='//TRIM(ADJUSTL(CWORD1))//'/'//TRIM(ADJUSTL(CWORD2))//',ITIME/NTIME=' &    
              //TRIM(ADJUSTL(CWORD3))//'/'//TRIM(ADJUSTL(CWORD4))
	POS1(1)=0.25;POS1(2)=0.05
	CALL SHOW_MTEXT(STEPINFO.info,1,POS1,BLACK,STEPSTATUSLIST)		  
    !INFO.QKEY=.TRUE.;INFO.ISNEEDINPUT=.FALSE.;INFO.COLOR=GREEN
    
	!ACTIVATE FACE AND EDGE
	!IF(.NOT.ALLOCATED(VISDEAD)) ALLOCATE(VISDEAD(NNODE))
    POSDATA.NODE.COORD(3)=0.D0
    DO I=1,POSDATA.NNODE
        IF(POSDATA.IX>0) POSDATA.NODE(I).COORD(1)=POSDATA.NODALQ(I,POSDATA.IX,STEPINFO.ISTEP)
        IF(POSDATA.IY>0) POSDATA.NODE(I).COORD(2)=POSDATA.NODALQ(I,POSDATA.IY,STEPINFO.ISTEP)
        IF(POSDATA.IZ>0) POSDATA.NODE(I).COORD(3)=POSDATA.NODALQ(I,POSDATA.IZ,STEPINFO.ISTEP)
    ENDDO
    
	FACE.ISDEAD=1;EDGE.ISDEAD=1;TET.ISDEAD=1;POSDATA.NODE.ISDEAD=1	
	DO I=1,NTET
		IF(POSDATA.ESET(TET(I).ISET).STEPSTATUS(STEPINFO.ISTEP)==1) THEN
			TET(I).ISDEAD=0
			FACE(TET(I).F(1:TET(I).NF)).ISDEAD=0
			EDGE(TET(I).E(1:TET(I).NE)).ISDEAD=0
			POSDATA.NODE(TET(I).V(1:TET(I).NV)).ISDEAD=0
		ENDIF	
	ENDDO
    MFACE.ISDEAD=1;MEDGE.ISDEAD=1;
	DO I=1,POSDATA.NEL
		IF(POSDATA.ESET(POSDATA.ELEMENT(I).ISET).STEPSTATUS(STEPINFO.ISTEP)==1) THEN			
			MFACE(POSDATA.ELEMENT(I).FACE(1:POSDATA.ELEMENT(I).NFACE)).ISDEAD=0
			MEDGE(ABS(POSDATA.ELEMENT(I).EDGE(1:POSDATA.ELEMENT(I).NEDGE))).ISDEAD=0			
		ENDIF	
	ENDDO
	
    IF(IsDrawVector) THEN
		CALL VEC_PLOT_DATA()
		CALL DrawVector()
	ENDIF
	
    IF(draw_surface_solid) CALL DrawSurfaceContour()
	IF(draw_contour) CALL DrawLineContour()
	
	IF(ISPLOTSLICE) THEN
		CALL GEN_SLICE_SURFACE_DATA()
		CALL GEN_SLICE_ISOLINE_DATA()
		CALL SLICEPLOT()
	ENDIF
    
	IF(isPlotStreamLine) THEN		
		CALL streamline_update()	
	ENDIF
        
    if(draw_surface_grid.or.show_edge.or.show_node.OR.SHOW_SET) then
        call drawGrid()
    endif
    
    IF(isProbeState) THEN
        CALL ProbeShow(PT_PROBE,IEL_PROBE)    
    ENDIF
ENDSUBROUTINE


subroutine display

! This gets called when the display needs to be redrawn

call reset_view

call glClear(ior(GL_COLOR_BUFFER_BIT,GL_DEPTH_BUFFER_BIT))
if(glIsList(StepStatusList)) call glCallList(StepStatusList)
if(draw_surface_solid.AND.glIsList(ContourList)) call glCallList(ContourList)
if(draw_contour.AND.glIsList(ContourLineList)) call glCallList(ContourLineList)
if(glIsList(VectorList)) call glCallList(VectorList)
if(glIsList(GridList)) call glCallList(GridList)
if(isProbeState.AND.glIsList(ProbeValuelist)) call glCallList(ProbeValuelist)
if(glIsList(slicelist).AND.ISPLOTSLICE) call glcalllist(slicelist)
if(ISPLOTSTREAMLINE.AND.glIsList(STREAMLINElist)) call glcalllist(STREAMLINElist)
if(glIsList(PSO_SLIP_PLOT_LIST)) call glcalllist(PSO_SLIP_PLOT_LIST)

call drawAxes()
if(IsDrawVector) call drawVectorLegend2(STEPPLOT.VMAX,STEPPLOT.VMin,STEPPLOT.VSCALE,ACTIVE_VECTOR_GROUP,VectorPairName,NVECTORPAIR)
if ((surface_color==rainbow_surface.and.(draw_surface_solid.or.draw_Contour)).OR.ISPLOTSLICESURFACE) then
    call Color_Bar()
endif
IF(ISSHOWNODALVALUE) CALL DRAW_NADAL_VALUE()
IF(LEN_TRIM(ADJUSTL(INFO.STR))>0) CALL SHOWINFO(INFO.COLOR)
IF(LINE_TEMP.SHOW) CALL LINE_TEMP.DRAW()
if(isProbeState_SFS.AND.glIsList(SFS_ProbeValuelist)) THEN
    CALL SLOPE_SFR_STATE_PARAMETER_SHOW()
    call glCallList(SFS_ProbeValuelist) 
ENDIF

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




subroutine SLICE_handler(selection)
    implicit none
    integer(kind=glcint), intent(in out) :: selection

    select case (selection)

    case(slice_location_click)
	    !isPickforslice=.true.
		!CALL GETSLICELOCATION()
        CALL INPUTSLICELOCATION()
    CASE(PLOTSLICESURFACE_CLICK)
		ISPLOTSLICESURFACE=.NOT.ISPLOTSLICESURFACE
    CASE(PLOTSLICEISOLINE_CLICK)
		ISPLOTSLICEISOLINE=.NOT.ISPLOTSLICEISOLINE
	CASE(PLOTSLICE_CLICK)
		ISPLOTSLICE=.NOT.ISPLOTSLICE
        
		
    end select
    
	call SLICEPLOT()
	
end subroutine

subroutine streamline_handler(selection)
    implicit none
    integer(kind=glcint), intent(in out) :: selection
    INTEGER::I
    
    select case (selection)
	
	
    case(streamline_location_click)
        isPickforstreamline=.true.
        IsFilterLocalMinimalMode=.false.
        isPlotStreamLine=.true.
        info.qkey=.true.
        call glutSetCursor(GLUT_CURSOR_CROSSHAIR)  
    CASE(Plot_streamline_CLICK)
		isPlotStreamLine=.not.isPlotStreamLine
    CASE(SHOW_STREAMLINE_NODE_CLICK)
        SHOW_STREAMLINE_NODE=.NOT.SHOW_STREAMLINE_NODE
    case(Reset_streamline_CLICK)
        CALL RESET_STREAMLINE()
     CASE(OUTPUT_SHOWNSTREAMLINE_CLICK)
        CALL OUT_SHOWN_STREAMLINE()
     CASE(Enlarger_StrokeFontSize_CLICK) 
        stroke_fontsize=2.*stroke_fontsize
       
     CASE(Smaller_StrokeFontSize_CLICK)
        stroke_fontsize=0.5*stroke_fontsize
        
  !  CASE(PLOTSLICEISOLINE_CLICK)
		!ISPLOTSLICEISOLINE=.NOT.ISPLOTSLICEISOLINE
        
		
    end select
    
    CALL STREAMLINE_PLOT()
    !CALL STREAMLINE_INI()
	!call SLICEPLOT()
	
end subroutine

SUBROUTINE RESET_STREAMLINE()
    INTEGER::I
    
    DO I=1,NSTREAMLINE
        IF(ALLOCATED(STREAMLINE(I).V)) THEN
            DEALLOCATE(STREAMLINE(I).V,STREAMLINE(I).VAL)
        ENDIF
        IF(ALLOCATED(STREAMLINE(I).PARA_SFCAL)) DEALLOCATE(STREAMLINE(I).PARA_SFCAL)
        STREAMLINE(I).ISINPUTSLIP=0
        STREAMLINE(I).SHOW=.TRUE.
    ENDDO
    nstreamline=0

ENDSUBROUTINE

subroutine slope_handler(selection)
    USE SolverMath
    use solverds,only:strtoint
    USE SLOPE_PSO, ONLY:SLOPE_OPTIM,POLYLINE_FOS_CAL
    implicit none
    integer(kind=glcint), intent(in out) :: selection
    INTEGER::NSLIP1,I,IEL1,J,K,n1=0,N2=0,OFFSET1,N3=0
    INTEGER,SAVE::TRYIEL1=0
    !INTEGER,EXTERNAL::POINTlOC_BC
    real(8)::AR2D1(2,2000),AR2D2(2,2000)
    INTEGER::unit,NAR1=100,NTOREAD1=100,NSET1=10,IERR
    INTEGER::NREADIN1,NSETREADIN1,NNODE1,EF,FLAG1
    REAL(8)::AR1(100),DT1,FA1,FM1,PT1(3)
    CHARACTER(32)::SET1(10)
    CHARACTER(512)::HELPSTRING
    REAL(IWP),ALLOCATABLE::NODE1(:,:)
    
    select case (selection)
    !CASE(SEARCH_BY_PSO)
    !    CALL SLOPE_OPTIM()
    CASE(ShowStreamline_AdmissiblityImprove_CLICK)
        CALL Improve_StreamlineShowed_Admissiblity()
    case(SFS_SLOPE_CLICK)
        isProbeState_SFS=.not.isProbeState_SFS
    case(ShowReadinSlip_Slope_CLICK)
        isShowReadinSlip=.not.isShowReadinSlip
        DO I=1,NSTREAMLINE
            IF(STREAMLINE(I).ISINPUTSLIP==1) THEN
                STREAMLINE(I).SHOW=isShowReadinSlip
            ENDIF
        ENDDO
        CALL STREAMLINE_PLOT()
    CASE(CHECKADMISSIBILITY)
        SLOPE_CHECK_ADMISSIBILITY=.NOT.SLOPE_CHECK_ADMISSIBILITY
    case(STREAMLINE_SLOPE_CLICK)
        ISSTREAMLINESLOPE=.NOT.ISSTREAMLINESLOPE
        IF(ISSTREAMLINESLOPE) THEN
            call slopestability_streamline(0) 
        ENDIF
        CALL STREAMLINE_PLOT()
    CASE(ShowTopTen_SLOPE_CLICK)
        IS_JUST_SHOW_TOP_TEN_SLOPE=.TRUE.
        !IF(IS_JUST_SHOW_TOP_TEN_SLOPE) THEN
        !    IS_JUST_SHOW_THE_MINIMAL_ONE_SLOPE=.FALSE.
        !    IS_SHOW_ALL_SLOPE=.FALSE.
            STREAMLINE(1:NSTREAMLINE).SHOW=.FALSE.           
            !IF(IS_JUST_SHOW_TOP_TEN_SLOPE) 
            STREAMLINE(SF_SLOPE(1:MIN(10,NSTREAMLINE))).SHOW=.TRUE.            
        !ELSE
        !    IS_SHOW_ALL_SLOPE=.TRUE.
        !    STREAMLINE(1:NSTREAMLINE).SHOW=.TRUE.
        !ENDIF
        CALL STREAMLINE_PLOT()
    CASE(SHOWMINIMAL_SLOPE_CLICK)
        IS_JUST_SHOW_THE_MINIMAL_ONE_SLOPE=.TRUE.
        !
        !IF(IS_JUST_SHOW_THE_MINIMAL_ONE_SLOPE) THEN
        !    IS_JUST_SHOW_TOP_TEN_SLOPE=.FALSE.
        !    IS_SHOW_ALL_SLOPE=.FALSE.
            STREAMLINE(1:NSTREAMLINE).SHOW=.FALSE.
            STREAMLINE(SF_SLOPE(1)).SHOW=.TRUE.
        !ELSE
        !    IS_SHOW_ALL_SLOPE=.TRUE.
        !    STREAMLINE(1:NSTREAMLINE).SHOW=.TRUE.
        !ENDIF        
        CALL STREAMLINE_PLOT()
    CASE(SHOWALL_SLOPE_CLICK)
        IS_SHOW_ALL_SLOPE=.TRUE.
        !IF(IS_SHOW_ALL_SLOPE) THEN
        !    IS_JUST_SHOW_TOP_TEN_SLOPE=.FALSE.
        !    IS_JUST_SHOW_THE_MINIMAL_ONE_SLOPE=.FALSE.            
            WHERE(STREAMLINE(1:NSTREAMLINE).ISCOMPATABLE) STREAMLINE(1:NSTREAMLINE).SHOW=.TRUE.
        !ELSE
        !    STREAMLINE(1:NSTREAMLINE).SHOW=.FALSE.
        !    IF(IS_JUST_SHOW_THE_MINIMAL_ONE_SLOPE) STREAMLINE(SF_SLOPE(MIN(1,NSTREAMLINE))).SHOW=.TRUE.
        !    IF(IS_JUST_SHOW_TOP_TEN_SLOPE) STREAMLINE(SF_SLOPE(1:MIN(10,NSTREAMLINE))).SHOW=.TRUE.
        !ENDIF         
        CALL STREAMLINE_PLOT()
    CASE(ShowLocalMinimal_Slope_Click)
        call ShowLocalMinimalSlope()
        
        CALL STREAMLINE_PLOT()
        
    case(FilterLocalMinimal_Slope_Click)
        IsFilterLocalMinimalMode=.not.IsFilterLocalMinimalMode
        if(IsFilterLocalMinimalMode) then
            isPickforstreamline=.false. !avoiding conflict
            isPlotStreamLine=.true.
            info.qkey=.true.
            call glutSetCursor(GLUT_CURSOR_CROSSHAIR) 
        
        endif
    CASE(SEARCH_SLOPE_CLICK)
        info.color=red;info.qkey=.true.
        info.str='Searching.Please Wait...'
        call glutPostRedisplay
        CALL RESET_STREAMLINE()
        CALL SEARCH_MINIMAL_SF_SLOPE()
        CALL STREAMLINE_PLOT()
        
    CASE(ReadSlipSurface_slope_click)
      helpstring='输入格式:\n &
            & 0) nslip, [flag1]  &
            & // nslip=读入的滑弧数。 flag1:=1,表示仅输入滑弧的一个节点,其它节点由积分确定; =0,读入滑弧上的所有节点,滑弧不再进行积分确定(默认) \n  &
            & // 如果(flag==1) 则按 1） \n &
            & 1) x, y, [z]      //每行为一个滑弧上的一个节点坐标,共nlip行。\n &
            & // 否则  依次输入按 1.0）和 1.1）,共nslip组。\n &
            & 1.0) nnode, [dx]   //nnode=每滑弧上节点数。dx: 如果大于0,则当输入节点相邻间距大于dx时,对其做进一步的细分,dx为细分后的间距。\n &
            & 1.1) x,y,[z]   //共nnode行 \n'C

        WRITE(*,'(A)') HELPSTRING 

        EF=0
        open(10,file='',status='old',iostat=EF)
        IF(EF/=0) RETURN
        CALL strtoint(10,AR1,NAR1,NREADIN1,NTOREAD1,SET1,NSET1,NSETREADIN1)
        NSLIP1=INT(AR1(1))
        IF(NREADIN1>1) THEN
            FLAG1=INT(AR1(2))
        ELSE
            FLAG1=0 !THE STREAMLINE NODE IS INPUT AND NO NEED TO INTEGRATE IN
        ENDIF
        !
        IF(FLAG1/=1) THEN
            DO K=1,NSLIP1
                IF(NSTREAMLINE+1>MAXNSTREAMLINE) CALL ENLARGE_STREAMLINE()
                NSTREAMLINE=NSTREAMLINE+1
                CALL strtoint(10,AR1,NAR1,NREADIN1,NTOREAD1,SET1,NSET1,NSETREADIN1)
                N1=INT(AR1(1));AR2D1=0.d0;AR2D2=0.d0
                DO I=1,N1            
                    READ(10,*,IOSTAT=IERR) AR2D1(:,I)
                    IF (IERR > 0) THEN
                        WRITE(*,*) 'Check input.  Something was wrong,N=',I
                        CLOSE(10)
                        RETURN
                    ELSE IF (IERR < 0) THEN
                        WRITE(*,*)  'END_OF_FILE.,N=',I
                        CLOSE(10)
                        RETURN
                    ENDIF
                ENDDO
                IF(NREADIN1>1.AND.AR1(2)>1.D-7) THEN
                    DT1=AR1(2) !进一步细分
                    N2=0
                    DO I=1,N1-1
                        CALL EQUDIVIDE(AR2D1(:,I),AR2D1(:,I+1),DT1,NODE1,NNODE1)
                        !GET RID OFF THE REPEATED NODES
                        IF(I>1) then
                            OFFSET1=1 
                        ELSE
                            OFFSET1=0
                        ENDIF
                        IF(N2+NNODE1-OFFSET1>2000) THEN
                            STOP 'OVERSIZE.SUB=slope_handler'                        
                        ENDIF
                        AR2D2(:,N2+1:N2+NNODE1-OFFSET1)=NODE1(:,1+OFFSET1:NNODE1)
                        N2=N2+NNODE1-OFFSET1                    
                    ENDDO  
                    STREAMLINE(NSTREAMLINE).NV=N2                
                ELSE
                    STREAMLINE(NSTREAMLINE).NV=N1
                    AR2D2(:,1:N1)=AR2D1(:,1:N1)
                ENDIF
                
                
                IF(ALLOCATED(STREAMLINE(NSTREAMLINE).V)) DEALLOCATE(STREAMLINE(NSTREAMLINE).V)                
                ALLOCATE(STREAMLINE(NSTREAMLINE).V(3,STREAMLINE(NSTREAMLINE).NV))
                IF(ALLOCATED(STREAMLINE(NSTREAMLINE).VAL)) DEALLOCATE(STREAMLINE(NSTREAMLINE).VAL)
                ALLOCATE(STREAMLINE(NSTREAMLINE).VAL(1:POSDATA.NVAR,STREAMLINE(NSTREAMLINE).NV))
                STREAMLINE(NSTREAMLINE).VAL=0.D0
                STREAMLINE(NSTREAMLINE).V(3,:)=0.d0
                N3=0
                DO J=1,STREAMLINE(NSTREAMLINE).NV
                    !READ(10,*) STREAMLINE(NSTREAMLINE).V(1:2,J)
                    STREAMLINE(NSTREAMLINE).V(1:2,J)=AR2D2(:,J)
                    IEL1=POINTlOC_BC(STREAMLINE(NSTREAMLINE).V(:,j),TRYIEL1)
                    TRYIEL1=IEL1
                    IF(IEL1>0) THEN !输入的滑线头尾可能在区域外,去掉。
                        N3=N3+1
                        call getval(STREAMLINE(NSTREAMLINE).V(:,J),iel1,STREAMLINE(NSTREAMLINE).VAL(:,N3))
                        IF(J/=N3) STREAMLINE(NSTREAMLINE).V(:,N3)=STREAMLINE(NSTREAMLINE).V(:,J)
                    ENDIF
                ENDDO
                STREAMLINE(NSTREAMLINE).NV=N3
                CALL slopestability_streamline(NSTREAMLINE)
                STREAMLINE(NSTREAMLINE).ISINPUTSLIP=1
                isShowReadinSlip=.TRUE.
                STREAMLINE(NSTREAMLINE).SHOW=.TRUE.
                NINPUTSLIP=NINPUTSLIP+1
                INPUTSLIP(NINPUTSLIP)=NSTREAMLINE
                
                !CALL POLYLINE_FOS_CAL(STREAMLINE(NSTREAMLINE).V(1:2,:),STREAMLINE(NSTREAMLINE).NV,0.5,FM1,FA1)
                !!
                !PRINT *, "FM1,FA1,FOS=",FM1,FA1,FA1/FM1
                
            ENDDO
                
        ELSE
            DO I=1,NSLIP1
                READ(10,*) PT1(1:POSDATA.NDIM)
                N1=POINTlOC_BC(PT1,0)
                IF(N1<1.OR.N1>NTET) THEN
                    !由于读入的数据可能存在误差,导致点处于模型外,此时令z=z-0.005,企图将点拉入模型内(假定输入的点处于模型的上表面)
                    PT1(POSDATA.NDIM)=PT1(POSDATA.NDIM)-0.005
                    N1=POINTlOC_BC(PT1,0)
                ENDIF
                
                IF(N1>0.AND.N1<=NTET) THEN
                    call gen_new_streamline(PT1)
                    NINPUTSLIP=NINPUTSLIP+1
                    INPUTSLIP(NINPUTSLIP)=NSTREAMLINE
                    STREAMLINE(nstreamline).isinputslip=1
                    isShowReadinSlip=.TRUE.
                    STREAMLINE(NSTREAMLINE).SHOW=.TRUE.                    
                ELSE
                    PRINT *, 'THE INPUT POINT(I) SEEMS TO OUTSIDE THE MODEL ZONE. I= ', I
                ENDIF
            ENDDO
            
        ENDIF        
        PRINT *,'INPUT DONE!'
        
        CALL STREAMLINE_PLOT()
        close(10)
    CASE(ShowSFContour_click)
        N1=-1
        CALL SET_VARIABLE_SHOW(N1)
    CASE(ShowSSTRIP_click)        
        IsShowSlideStrip=.not.IsShowSlideStrip   
        CALL STREAMLINE_PLOT() 
     CASE(ShowFSLabel_click)        
        IsShowFSLabel=.not.IsShowFSLabel   
        CALL STREAMLINE_PLOT() 
    CASE(OutSlideInfo_click)
        CALL OUT_3D_SLOPE_SF_DIRECTION()          
    CASE default
		!isPlotStreamLine=.not.isPlotStreamLine
  !  CASE(PLOTSLICEISOLINE_CLICK)
		!ISPLOTSLICEISOLINE=.NOT.ISPLOTSLICEISOLINE
        
		
    end select
    
    !CALL STREAMLINE_INI()
	!call SLICEPLOT()
	
end subroutine

subroutine probe_handler(selection)
    use solverds,only:strtoint
    USE SolverMath
    USE IFPORT
    implicit none
    integer(kind=glcint), intent(in out) :: selection
    INTEGER::NSLIP1,I,IEL1,J,K,n1=0,N2=0,OFFSET1
    INTEGER,SAVE::TRYIEL1=0
    !INTEGER,EXTERNAL::POINTlOC_BC
     real(8)::AR2D1(3,500),AR2D2(3,2000)
    INTEGER::unit,NAR1=100,NTOREAD1=100,NSET1=10
    INTEGER::NREADIN1,NSETREADIN1,NNODE1,EF
    REAL(8)::AR1(100),DT1,VAL1(100)
    CHARACTER(32)::SET1(10)
    REAL(IWP),ALLOCATABLE::NODE1(:,:)
    
    character(1024)::FILE1
	CHARACTER(3)        drive
	CHARACTER(1024)      dir
	CHARACTER(1024)      name
	CHARACTER(8)      ext
    integer(4)::msg1
   
    
    select case(selection)
    case(probe_selected)
	    isProbeState=.not.isProbestate
        if(.not.isProbeState) call glutSetCursor(GLUT_CURSOR_LEFT_ARROW)
    case(probe_get_data_click)
        EF=0
        open(10,file='',status='old',iostat=EF)
        IF(EF/=0) RETURN
        inquire(10,name=FILE1)
        MSG1 = SPLITPATHQQ(FILE1, drive, dir, name, ext)
        FILE1=trim(drive)//trim(dir)//trim(name)//'_PointValue.dat'
        open(20,file=FILE1,status='replace')
        
        WRITE(20,10) (TRIM(POSDATA.OUTVAR(I).NAME),I=1,POSDATA.NVAR)
        
        READ(10,*) NSLIP1
        DO K=1,NSLIP1
            CALL strtoint(10,AR1,NAR1,NREADIN1,NTOREAD1,SET1,NSET1,NSETREADIN1)
            N1=INT(AR1(1));AR2D1=0.d0;AR2D2=0.d0
            READ(10,*) AR2D1(1:POSDATA.NDIM,1:N1)
            IF(NREADIN1>1.AND.AR1(2)>1.D-7) THEN
                DT1=AR1(2)
                N2=0
                DO I=1,N1-1
                    CALL EQUDIVIDE(AR2D1(1:POSDATA.NDIM,I),AR2D1(1:POSDATA.NDIM,I+1),DT1,NODE1,NNODE1)
                    !GET RID OFF THE REPEATED NODES
                    IF(I>1) then
                        OFFSET1=1 
                    ELSE
                        OFFSET1=0
                    ENDIF
                    IF(N2+NNODE1-OFFSET1>2000) THEN
                        STOP 'OVERSIZE.SUB=probe_handler'                        
                    ENDIF
                    AR2D2(1:POSDATA.NDIM,N2+1:N2+NNODE1-OFFSET1)=NODE1(1:POSDATA.NDIM,1+OFFSET1:NNODE1)
                    N2=N2+NNODE1-OFFSET1                    
                ENDDO  
                N1=N2                
            ELSE
                AR2D2(:,1:N1)=AR2D1(:,1:N1)
            ENDIF
            
            DO J=1,N1   
                IEL1=POINTlOC_BC(AR2D2(:,J),TRYIEL1)
                TRYIEL1=IEL1
                IF(IEL1>0) THEN
                    call getval(AR2D2(:,J),iel1,VAL1)
                    WRITE(20,20) K,J,VAL1(1:POSDATA.NVAR)
                ELSE
                    WRITE(20,30) K,J,AR2D2(1:POSDATA.NDIM,J)
                ENDIF
            ENDDO        
    
        END DO
    
        CLOSE(10);CLOSE(20)
        
        PRINT *, 'OUTDATA DONE!'
    
    end select

10  FORMAT('ILINE',1X,'  INODE',1X,<POSDATA.NVAR>(A24,1X))
20  FORMAT(I5,1X,I7,1X,<POSDATA.NVAR>(EN24.15,1X))
30  FORMAT(I5,1X,I7,1X,<POSDATA.NDIM>(EN24.15,1X),'*THE POINT IS OUT OF MODEL ZONE.*')
    
endsubroutine

SUBROUTINE set_variable_show(VALUE)
integer(kind=glcint), intent(in out) :: value 

CONTOUR_PLOT_VARIABLE=VALUE
IF(CONTOURBAR.IVARPLOT/=CONTOUR_PLOT_VARIABLE) call initialize_contourplot(CONTOUR_PLOT_VARIABLE)
CALL DrawSurfaceContour()
call DrawLineContour()

RETURN
END SUBROUTINE

SUBROUTINE SET_NODAL_VARIABLE_SHOW(VALUE)
integer(kind=glcint), intent(in out) :: value 

SHOW_NODAL_VALUE=VALUE
ISSHOWNODALVALUE=.TRUE.

CALL DRAW_NADAL_VALUE()

RETURN
END SUBROUTINE

SUBROUTINE SET_SLICE_VARIABLE_SHOW(VALUE)
integer(kind=glcint), intent(in out) :: value 

SLICE_PLOT_VARIABLE=VALUE
IF(NSLICE<1) CALL INPUTSLICELOCATION()

CALL SLICEPLOT()
RETURN
END SUBROUTINE


SUBROUTINE SHOW_NodalValue_HANDLER(selection)
integer(kind=glcint), intent(in out) :: selection

select case (selection)

case (SHOWNODALVALUE_TOGGLE)
   ISSHOWNODALVALUE = .not. ISSHOWNODALVALUE
end select

CALL DRAW_NADAL_VALUE()
    
ENDSUBROUTINE

SUBROUTINE SET_VECTOR_PLOT_GROUP(VALUE)
integer(kind=glcint), intent(in out) :: value 

VECTOR_PLOT_GROUP=VALUE
ACTIVE_VECTOR_GROUP(VALUE)=.NOT.ACTIVE_VECTOR_GROUP(VALUE)

IsDrawVector=.true.

CALL VEC_PLOT_DATA()

CALL drawvector()

!IF(isPlotStreamLine) THEN
!	
!	CALL streamline_update()

!ENDIF

RETURN
END SUBROUTINE

SUBROUTINE SET_STREAMLINE_VECTOR_GROUP(VALUE)
integer(kind=glcint), intent(in out) :: value 



IF(isPlotStreamLine) THEN
	IF(VALUE/=STREAMLINE_VECTOR) THEN
		STREAMLINE_VECTOR=VALUE
		CALL streamline_update()
	ENDIF
ENDIF

RETURN
END SUBROUTINE



subroutine Contour_handler(selection)
USE IFQWIN
integer(kind=glcint), intent(in out) :: selection
integer::res1
type(qwinfo) winfo

real(8)::minv1,maxv1

select case (selection)

case (Contour_surfsolid_toggle)
   draw_surface_solid = .not. draw_surface_solid
case (Contour_Line_toggle)
   draw_contour = .not. draw_contour
case(contour_densify)
    Scale_Contour_Num=Scale_Contour_Num*2
	call initialize_contourplot(CONTOUR_PLOT_VARIABLE)
    call DrawSurfaceContour()
    call DrawLineContour()
case(contour_sparsify)
    Scale_Contour_Num=Scale_Contour_Num/2
	call initialize_contourplot(CONTOUR_PLOT_VARIABLE)
    call DrawSurfaceContour()
    call DrawLineContour()
CASE(Contour_In_DeformedMesh)
    IsContour_In_DeformedMesh=.NOT.IsContour_In_DeformedMesh
    call DrawSurfaceContour()
    call DrawLineContour()
CASE(Coutour_MinMax)
    
    info.str='Please input the Min and Max for the contour in the console window.'
    info.color=green;info.qkey=.true.
    res1=FOCUSQQ(QWIN$FRAMEWINDOW)

    PRINT *, 'INPUT THE MIN AND MAX VALUE'
    READ(*,*) minv1,maxv1
	call initialize_contourplot(CONTOUR_PLOT_VARIABLE,minv1,maxv1)
    call DrawSurfaceContour()
    call DrawLineContour()    
end select

!if(selection/=Contour_Line_toggle) call DrawSurfaceContour()
!if(selection/=Contour_surfsolid_toggle) call DrawLineContour()

endsubroutine

subroutine Vector_Location_handler(selection)
integer(kind=glcint), intent(in out) :: selection

vectorloc=selection
call VEC_PLOT_DATA()
call drawvector()
end subroutine

subroutine Vector_Base_handler(selection)
integer(kind=glcint), intent(in out) :: selection

vectorbase=selection-7

call drawvector()
end subroutine

subroutine Vector_Format_handler(selection)
integer(kind=glcint), intent(in out) :: selection
select case (selection)
case(Arrow_toggle)
    IsArrow=.not.IsArrow 
end select
call drawvector()
end subroutine

subroutine Vector_handler(selection)
integer(kind=glcint), intent(in out) :: selection
select case (selection)


case (Vector_toggle)
   IsDrawVector = .not. IsDrawVector   
case(Vector_lengthen)
   Scale_Vector_len=Scale_Vector_len*2.0
case(Vector_shorten)
   Scale_Vector_len=Scale_Vector_len/2.0
case(Vector_Unify)
    IsVectorUnify=.not.IsVectorUnify
case(Vector_Less)
    VectorFrequency=max(VectorFrequency-1,1)
case(Vector_More)
    VectorFrequency=max(VectorFrequency+1,1)
   
end select

call drawvector()

endsubroutine

subroutine Model_handler(selection)
integer(kind=glcint), intent(in out) :: selection
select case (selection)


case (surfgrid_toggle)
   draw_surface_grid = .not. draw_surface_grid
case (edge_toggle)
   show_edge = .not. show_edge 
case (node_toggle)
   show_node = .not. show_node
CASE(SHOW_SET_TOGGLE)
    SHOW_SET=.NOT.SHOW_SET
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

subroutine Color_Map_handler(selection)
integer(kind=glcint), intent(in out) :: selection

iCOLORMAP=selection
surface_color=rainbow_surface
call Color_Bar()
call DrawSurfaceContour()

end subroutine Color_Map_handler


subroutine make_menu(submenuid)
USE SLOPE_PSO,ONLY:PSO_SLIP_PLOT_MENU
integer, intent(in) :: submenuid
integer :: menuid, param_id, contour_color_menu, surface_color_menu,I
INTEGER::VSHOW_SPG_ID,VSHOW_STRESS_ID,VSHOW_STRAIN_ID,VSHOW_PSTRAIN_ID,&
        VSHOW_SFR_ID,VSHOW_DIS_ID,VSHOW_FORCE_ID,CONTOUR_PLOT_ID,VECTOR_PLOT_ID,&
        VECTOR_PAIR_ID,Model_ID,&
        NODALVALSHOW_SPG_ID,&
        NODALVALSHOW_STRESS_ID,&
        NODALVALSHOW_STRAIN_ID,&
        NODALVALSHOW_PSTRAIN_ID,&
        NODALVALSHOW_SFR_ID,&
        NODALVALSHOW_DIS_ID,&
        NODALVALSHOW_FORCE_ID,&
        Show_NodalValue_ID,&
        SLICE_PLOT_ID,&
        SLICESHOW_SPG_ID,&
        SLICESHOW_STRESS_ID,&
        SLICESHOW_STRAIN_ID,&
        SLICESHOW_PSTRAIN_ID,&
        SLICESHOW_SFR_ID,&
        SLICESHOW_DIS_ID,&
        SLICESHOW_FORCE_ID,&
        STREAMLINE_PLOT_ID,&
        SLOPE_ID,Probe_ID,&
		STREAMLINE_VECTOR_ID,Vector_Format_ID,Vector_Base_ID,Vector_LOC_ID,&
        ColorMap_ID,EA_SLOPE_ID

        
contour_color_menu = glutCreateMenu(contour_color_handler)
call glutAddMenuEntry("black",black_contour)
call glutAddMenuEntry("contour value",rainbow_contour)

ColorMap_ID=glutCreateMenu(Color_Map_handler)
Call glutAddMenuEntry("Rainbow",CM_Rainbow)
Call glutAddMenuEntry("RainBow_DarkEnds",CM_RainBow_DarkEnds)
Call glutAddMenuEntry("RainBowLight",CM_RainBowLight)
Call glutAddMenuEntry("RainBowPastel",CM_RainBowPastel)
Call glutAddMenuEntry("Terrain",CM_Terrain)
Call glutAddMenuEntry("Gray",CM_Gray)
Call glutAddMenuEntry("HotMetal",CM_HotMetal)
Call glutAddMenuEntry("Doppler",CM_Doppler)
Call glutAddMenuEntry("Elevation",CM_Elevation)
Call glutAddMenuEntry("Magma",CM_Magma)
Call glutAddMenuEntry("Accent",CM_Accent)
Call glutAddMenuEntry("Red",CM_Red)
Call glutAddMenuEntry("Blue",CM_Blue)
Call glutAddMenuEntry("Green",CM_Green)
Call glutAddMenuEntry("BlueRed",CM_BlueRed)
Call glutAddMenuEntry("BlueYellowRed",CM_BlueYellowRed)
Call glutAddMenuEntry("HighPoint",CM_HighPoint)
Call glutAddMenuEntry("HighPoint2",CM_HighPoint2)
Call glutAddMenuEntry("YellowHigh",CM_YellowHigh)
Call glutAddMenuEntry("Soil",CM_Soil)
Call glutAddMenuEntry("Geology",CM_Geology)
Call glutAddMenuEntry("Land",CM_Land)


surface_color_menu = glutCreateMenu(surface_color_handler)
call glutAddMenuEntry("white",white_surface)
call glutAddSubMenu("ColorMap",ColorMap_ID)
call glutAddMenuEntry("Set/Material",color_by_set)
call glutAddMenuEntry("transparency",transparency)

!param_id = glutCreateMenu(param_handler)
!call glutAddMenuEntry("number of x grid intervals",set_nx)
!call glutAddMenuEntry("number of y grid intervals",set_ny)
!call glutAddMenuEntry("reset to initial parameters",reset_params)

VSHOW_SPG_ID=glutCreateMenu(SET_VARIABLE_SHOW)
DO I=1,POSDATA.NVAR
    CALL GLUTADDMENUENTRY(TRIM(ADJUSTL(POSDATA.OUTVAR(I).NAME)),I)
ENDDO


CONTOUR_PLOT_ID=glutCreateMenu(contour_handler)
call glutAddMenuEntry("toggle contour surface",Contour_surfsolid_toggle)
call glutAddSubMenu("contour surface color",surface_color_menu)
call glutAddMenuEntry("toggle contour line",Contour_Line_toggle)
call glutAddSubMenu("contour line color",contour_color_menu)
call glutAddMenuEntry("Densify contour",contour_densify)
call glutAddMenuEntry("Sparsify contour",contour_Sparsify)
call glutAddMenuEntry("Plot In DeformedGrid",Contour_In_DeformedMesh)
call glutAddMenuEntry("ResetContourMinMax",Coutour_MinMax)
call glutAddSubMenu("SHOW_VARS",VSHOW_SPG_ID)



VECTOR_PAIR_ID=glutCreateMenu(SET_VECTOR_PLOT_GROUP)
IF(POSDATA.IDISX>0.AND.POSDATA.IDISY>0) THEN
    CALL GLUTADDMENUENTRY ("DISPLACE",VECTOR_GROUP_DIS)   
ENDIF
IF(POSDATA.IVX>0.AND.POSDATA.IVY>0) THEN
    CALL GLUTADDMENUENTRY ("SEEPAGE_VELOCITY",VECTOR_GROUP_SEEPAGE_VEC)
ENDIF
IF(POSDATA.IGRADX>0.AND.POSDATA.IGRADY>0) THEN
    CALL GLUTADDMENUENTRY ("SEEPAGE_GRADIENT",VECTOR_GROUP_SEEPAGE_GRAD)
ENDIF
IF(POSDATA.ISFR_SFRX>0.AND.POSDATA.ISFR_SFRY>0) THEN
    CALL GLUTADDMENUENTRY ("STRESS_FAILUE_FACE",VECTOR_GROUP_SFR)
ENDIF
IF(POSDATA.IPSIGMA1>0) THEN
    CALL GLUTADDMENUENTRY ("PSIGMA1",VECTOR_GROUP_PSIGMA1)
ENDIF
IF(POSDATA.IPSIGMA2>0) THEN
    CALL GLUTADDMENUENTRY ("PSIGMA2",VECTOR_GROUP_PSIGMA2)
ENDIF
IF(POSDATA.IPSIGMA3>0) THEN
    CALL GLUTADDMENUENTRY ("PSIGMA3",VECTOR_GROUP_PSIGMA3)
ENDIF

Vector_Base_ID=glutCreateMenu(VECTOR_Base_HANDLER)
CALL glutAddMenuEntry("OrigBase",Vector_Base_Orig)
CALL glutAddMenuEntry("EndBase",Vector_Base_End)
CALL glutAddMenuEntry("CentBase",Vector_Base_Cent)

Vector_Format_ID=glutCreateMenu(VECTOR_FORMAT_HANDLER)
CALL glutAddMenuEntry("ToggleArrow",Arrow_toggle)
call glutAddSubMenu("VectorBase",Vector_Base_ID)

Vector_LOC_ID=glutCreateMenu(VECTOR_LOCATION_HANDLER)
CALL glutAddMenuEntry("NODES",VECTORLOC_NODE)
CALL glutAddMenuEntry("UniformGrid",VECTORLOC_Grid)

VECTOR_PLOT_ID=glutCreateMenu(VECTOR_HANDLER)
call glutAddSubMenu("VectorLocation",Vector_LOC_ID)
call glutAddMenuEntry("toggle Vector Plot",Vector_toggle)
call glutAddMenuEntry("Lengthen Vector",Vector_Lengthen)
call glutAddMenuEntry("Shorten Vector",Vector_Shorten)
call glutAddMenuEntry("VectorFrequency++",Vector_More)
call glutAddMenuEntry("VectorFrequency--",Vector_Less)
call glutAddMenuEntry("UnifyingVector",Vector_unify)
call glutAddSubMenu("VectorFormat",Vector_Format_ID)
call glutAddSubMenu("Toggle VectorPair",VECTOR_PAIR_ID)

Model_ID=glutCreateMenu(Model_HANDLER)
call glutAddMenuEntry("ShowMesh",surfgrid_toggle)
call glutAddMenuEntry("ShowEdge",edge_toggle)
call glutAddMenuEntry("ShowNode",node_toggle)
call glutAddMenuEntry("ShowSet",show_set_toggle)
call glutAddMenuEntry("DeformedGrid Toggle",DeformedMesh)
call glutAddMenuEntry("++DeformedGridScale",Enlarge_Scale_DeformedMesh)
call glutAddMenuEntry("--DeformedGridScale",Minify_Scale_DeformedMesh)


NodalValShow_SPG_ID=glutCreateMenu(SET_NODAL_VARIABLE_SHOW)
DO I=1,POSDATA.NVAR
    CALL GLUTADDMENUENTRY(TRIM(ADJUSTL(POSDATA.OUTVAR(I).NAME)),I)    
ENDDO
CALL GLUTADDMENUENTRY('NodalID',NODALID_SHOW)
CALL GLUTADDMENUENTRY('ElementID',ELEMENTID_SHOW)

Show_NodalValue_ID=glutCreateMenu(SHOW_NodalValue_HANDLER)
call glutAddSubMenu("SHOW_VARS",NodalValShow_SPG_ID)
CALL glutAddMenuEntry("Show Toggle",SHOWNODALVALUE_TOGGLE)

SliceShow_SPG_ID=glutCreateMenu(SET_SLICE_VARIABLE_SHOW)
DO I=1,POSDATA.NVAR
    CALL GLUTADDMENUENTRY(TRIM(ADJUSTL(POSDATA.OUTVAR(I).NAME)),I)
ENDDO


SLICE_PLOT_ID=glutCreateMenu(SLICE_handler)
call glutAddMenuEntry("SET LOCTIONS",SLICE_LOCATION_CLICK)
call glutAddMenuEntry("PlotSlice toggle",PLOTSLICE_CLICK)
call glutAddMenuEntry("SliceSurfacePlot toggle",PLOTSLICESURFACE_CLICK)
call glutAddMenuEntry("SliceIsoLinePlot toggle",PLOTSLICEISOLINE_CLICK)
call glutAddSubMenu("SHOW_VARS",SLICESHOW_SPG_ID)


STREAMLINE_VECTOR_ID=glutCreateMenu(SET_STREAMLINE_VECTOR_GROUP)
IF(POSDATA.IDISX>0.AND.POSDATA.IDISY>0) THEN
    CALL GLUTADDMENUENTRY ("DISPLACE",VECTOR_GROUP_DIS)   
ENDIF
IF(POSDATA.IVX>0.AND.POSDATA.IVY>0) THEN
    CALL GLUTADDMENUENTRY ("SEEPAGE_VELOCITY",VECTOR_GROUP_SEEPAGE_VEC)
ENDIF
IF(POSDATA.IGRADX>0.AND.POSDATA.IGRADY>0) THEN
    CALL GLUTADDMENUENTRY ("SEEPAGE_GRADIENT",VECTOR_GROUP_SEEPAGE_GRAD)
ENDIF
IF(POSDATA.ISFR_SFRX>0.AND.POSDATA.ISFR_SFRY>0) THEN
    CALL GLUTADDMENUENTRY ("STRESS_FAILUE_FACE",VECTOR_GROUP_SFR)
ENDIF
IF(POSDATA.IPSIGMA1>0) THEN
    CALL GLUTADDMENUENTRY ("PSIGMA1",VECTOR_GROUP_PSIGMA1)
ENDIF
IF(POSDATA.IPSIGMA3>0) THEN
    CALL GLUTADDMENUENTRY ("PSIGMA3",VECTOR_GROUP_PSIGMA3)
ENDIF

STREAMLINE_PLOT_ID=glutCreateMenu(STREAMLINE_handler)
call glutAddSubMenu("SetStreamlineVector",STREAMLINE_VECTOR_ID)
call glutAddMenuEntry("PICK START POINT",STREAMLINE_LOCATION_CLICK)
call glutAddMenuEntry("PlotStreamLine toggle",Plot_streamline_CLICK)
call glutAddMenuEntry("OutPutShownStreamline",OUTPUT_SHOWNSTREAMLINE_CLICK)
call glutAddMenuEntry("ResetStreamLine",Reset_streamline_CLICK)
call glutAddMenuEntry("ShowStreamlineNode",SHOW_STREAMLINE_NODE_CLICK)
call glutAddMenuEntry("++StrokeFontSize",Enlarger_StrokeFontSize_CLICK)
call glutAddMenuEntry("--StrokeFontSize",Smaller_StrokeFontSize_CLICK)

EA_SLOPE_ID=PSO_SLIP_PLOT_MENU()

SLOPE_ID=glutCreateMenu(SLOPE_handler)
CALL glutAddMenuEntry("CheckAdmissiblity Toggle",CHECKADMISSIBILITY)
call glutAddMenuEntry("Search...",SEARCH_SLOPE_CLICK)
CALL glutAddSubMenu("SearchByEA",EA_SLOPE_ID)
call glutAddMenuEntry("JustShowTopTenSlips",ShowTopTen_SLOPE_CLICK)
call glutAddMenuEntry("JustShowMinimalSlip",ShowMinimal_SLOPE_CLICK)
call glutAddMenuEntry("ShowAllSlips",ShowAll_SLOPE_CLICK)
call glutAddMenuEntry("ShowReadinSlips",ShowReadinSlip_Slope_CLICK)
call glutAddMenuEntry("FilterLocalSmallSlip(DelectLocalLarger)",filterLocalMinimal_Slope_Click)
call glutAddMenuEntry("ShowLocalMinimalSlip",ShowLocalMinimal_Slope_Click)
call glutAddMenuEntry("StressFailureState",SFS_SLOPE_CLICK)
Call glutAddMenuEntry('ReadSlips',ReadSlipSurface_slope_click)
if(POSDATA.ndim==3) then
    call glutAddMenuEntry('ShowSFContour',ShowSFContour_click)
    call glutAddMenuEntry('ShowSlideStrip',ShowSSTRIP_click)
    call glutAddMenuEntry('ShowFSLabel',ShowFSLabel_click)
    CALL glutAddMenuEntry('OutSlideInfo',OutSlideInfo_click)
endif
call glutAddMenuEntry("AdmissiblityImprove...",ShowStreamline_AdmissiblityImprove_CLICK)
call glutAddMenuEntry("STREAMLINE_METHOD",STREAMLINE_SLOPE_CLICK)

PROBE_ID=glutCreateMenu(PROBE_handler)
call glutAddMenuEntry("ShowProbeValue toggle",probe_selected)
call glutAddMenuEntry("GetData...",Probe_Get_data_click)

menuid = glutCreateMenu(menu_handler)
call glutAddSubMenu("Contour",CONTOUR_PLOT_ID)
call glutAddSubMenu("Vector",VECTOR_PLOT_ID)
call glutAddSubMenu("Slice",SLICE_PLOT_ID)
call glutAddSubMenu("Streamline",STREAMLINE_PLOT_ID)
call glutAddSubMenu("NodalValue",Show_NodalValue_ID)
call glutAddSubMenu("ProbeValue",PROBE_ID)
call glutAddSubMenu("Model",Model_ID)
CALL glutAddSubMenu("SlopeStability",SLOPE_ID)
call glutAddSubMenu("View",submenuid)
!call glutAddSubMenu("plotting parameters",param_id)


call glutAddMenuEntry("quit",quit_selected)


call glutAttachMenu(GLUT_RIGHT_BUTTON)
end subroutine make_menu

SUBROUTINE OUT_SHOWN_STREAMLINE(FILE,OUTID,INFO)    
    USE IFPORT
    IMPLICIT NONE
    CHARACTER(*),OPTIONAL,INTENT(IN)::FILE,INFO
    INTEGER,OPTIONAL,INTENT(IN)::OUTID
    INTEGER::I,J,EF,OUTID1=0
    CHARACTER(1024)::FILE1
    LOGICAL::EXIST
    CHARACTER(8)::char_time
    
    PRINT *, 'SELECT TO FILE TO SAVE THE STREAMLINE DATA...'
    
    FILE1=""
    OUTID1=0
    IF(PRESENT(FILE)) FILE1=TRIM(ADJUSTL(FILE))
    IF(PRESENT(OUTID)) OUTID1=OUTID
  !logical :: exist
  !
  !inquire(file=TRIM(FILE1), exist=exist)
  !if (exist) then
  !  open(20, file=TRIM(FILE1), status="UNKNOWN", position="append", action="write")
  !else
  !  open(2O, file=TRIM(FILE1), status="new", action="write")
  !end if
    
    open(20,file=TRIM(FILE1),status="UNKNOWN", position="append",iostat=EF)
    IF(EF/=0) RETURN
    call TIME(char_time)
    WRITE(20,10) (TRIM(POSDATA.OUTVAR(I).NAME),I=1,POSDATA.NVAR),DATE(),char_time
        
    IF(OUTID1<1) THEN 
        DO I=1,NSTREAMLINE
           IF((.NOT.STREAMLINE(I).SHOW).OR.STREAMLINE(I).NV<2) CYCLE
            
            DO J=1,STREAMLINE(I).NV 
                IF(ALLOCATED(STREAMLINE(I).PARA_SFCAL)) THEN
                    WRITE(20,20) I,J,STREAMLINE(I).SF_SLOPE,STREAMLINE(I).VAL(1:POSDATA.NVAR,J),STREAMLINE(I).PARA_SFCAL(:,J)
                ELSE
                    WRITE(20,20) I,J,STREAMLINE(I).SF_SLOPE,STREAMLINE(I).VAL(1:POSDATA.NVAR,J),0,0,0,0,0,0,0,0
                ENDIF
                IF(J==STREAMLINE(I).NV) WRITE(20,'(I5)') I
            ENDDO        
    
        END DO
    
    ELSE
        I=OUTID1
        DO J=1,STREAMLINE(I).NV 
            IF(ALLOCATED(STREAMLINE(I).PARA_SFCAL)) THEN
                WRITE(20,20) I,J,STREAMLINE(I).SF_SLOPE,STREAMLINE(I).VAL(1:POSDATA.NVAR,J),STREAMLINE(I).PARA_SFCAL(:,J)
            ELSE
                WRITE(20,20) I,J,STREAMLINE(I).SF_SLOPE,STREAMLINE(I).VAL(1:POSDATA.NVAR,J),0,0,0,0,0,0,0,0
            ENDIF
            IF(J==STREAMLINE(I).NV) WRITE(20,'(I5)') I
        ENDDO     
    ENDIF
    
    CLOSE(20)
        
    PRINT *, 'OUTDATA DONE!'
    

    
10  FORMAT('ILINE',1X,'  INODE',1X,22X,'FOS',1X,<POSDATA.NVAR>(A24,1X),&
        15X,'SLIPSEG_SN',14X,'SLIPSEG_TAU',5X,'SLIPSEG_ANGLE_IN_RAD',13X,'SLIPSEG_TAUF',&
        14X,'SLIPSEG_LEN',16X,'SLIPSEG_C',14X,'SLIPSEG_PHI',14X,'SLIPSEG_SFR',12X,'SLIPFACE_COS1',12X,'SLIPFACE_COS2',12X,'SLIPFACE_COS3',5X'//DATE=',A8,',TIME=',A8)
20  FORMAT(I5,1X,I7,1X,<POSDATA.NVAR+12>(EN24.15,1X))
         
30  FORMAT(I5,1X,I7,1X,<POSDATA.NDIM>(EN24.15,1X),'*THE POINT IS OUT OF MODEL ZONE.*')

ENDSUBROUTINE

SUBROUTINE OUT_3D_SLOPE_SF_DIRECTION(FILE,OUTID,INFO)    
    USE IFPORT
    IMPLICIT NONE
    CHARACTER(*),OPTIONAL,INTENT(IN)::FILE,INFO
    INTEGER,OPTIONAL,INTENT(IN)::OUTID
    INTEGER::I,J,EF,OUTID1=0,N1
    CHARACTER(1024)::FILE1
    LOGICAL::EXIST
    CHARACTER(8)::char_time
    REAL(8)::XC1,YC1,XY1(2),T1,XY2(3)
    REAL(8),PARAMETER::PI=3.141592653589793
    PRINT *, 'SELECT TO FILE TO SAVE THE STREAMLINE DATA...'
    
    FILE1=""
    OUTID1=0
    IF(PRESENT(FILE)) FILE1=TRIM(ADJUSTL(FILE))
    IF(PRESENT(OUTID)) OUTID1=OUTID
  !logical :: exist
  !
  !inquire(file=TRIM(FILE1), exist=exist)
  !if (exist) then
  !  open(20, file=TRIM(FILE1), status="UNKNOWN", position="append", action="write")
  !else
  !  open(2O, file=TRIM(FILE1), status="new", action="write")
  !end if
    
    open(20,file=TRIM(FILE1),status="REPLACE",iostat=EF)
    IF(EF/=0) RETURN
         
    WRITE(20,10)
    XC1=(POSDATA.MINX+POSDATA.MAXX)/2
    YC1=(POSDATA.MINY+POSDATA.MAXY)/2
    
    DO I=1,NSTREAMLINE
        IF((.NOT.STREAMLINE(I).ISCOMPATABLE).OR.STREAMLINE(I).NV<2) CYCLE
        !DO J=1,2
        N1=1
            !IF(J==2) N1=STREAMLINE(I).NV
        !reference vector
        XY1(1)=STREAMLINE(I).V(1,N1)-XC1
        XY1(2)=STREAMLINE(I).V(2,N1)-YC1
        T1=ATAN2(XY1(2),XY1(1))+PI/2
        XY1(1)=COS(T1);XY1(2)=SIN(T1)
        XY2=STREAMLINE(I).SLIDE_STRIP(:,N1,2)-STREAMLINE(I).SLIDE_STRIP(:,N1,1) 
        T1=NORM2(XY2(1:2))
        XY2=XY2/T1
        T1=DOT_PRODUCT(XY1,XY2(1:2))
        IF(T1<0.D0)  XY2=-XY2
        STREAMLINE(I).SVEC=XY2(1:2)
        WRITE(20,20) I,N1,STREAMLINE(I).V(:,N1),STREAMLINE(I).SF_SLOPE,XY2(1:2)
        !ENDDO
    END DO
    
    
    
    CLOSE(20)
        
    PRINT *, 'OUTDATA DONE!'
    

    
10  FORMAT('ILINE',10X,'NODE',10X,'X',10X,'Y',10X,'Z',10X,'FS',10X,'DX',10X,'DY')
20  FORMAT(2(I7,1X),6(EN24.15,1X))
    

ENDSUBROUTINE

subroutine initialize_contourplot(IVARPLOT,MINV1,MAXV1)
    !use solverds
    implicit none
	INTEGER,INTENT(IN)::IVARPLOT
    real(GLDOUBLE),OPTIONAL::MINV1,MAXV1
    integer::i
    real(GLDOUBLE) :: graphmin,graphmax,graphstep
    !real(GLDOUBLE),allocatable :: contour_value(:)
    character(128)::str1,str2  
    
    CONTOURBAR.IVARPLOT=IVARPLOT
    IF(PRESENT(MINV1)) THEN
        minv=MINV1
    ELSE
        minv=minval(POSDATA.NODALQ(:,IVARPLOT,:))
    ENDIF
    IF(PRESENT(MAXV1)) THEN
        MAXV=MAXV1
    ELSE
        MAXV=maxval(POSDATA.NODALQ(:,IVARPLOT,:))
    ENDIF    


    WRITE(STR1,'(G13.3)') MAXV;WRITE(STR2,'(G13.3)') MINV;

    contourbar.title=trim(adjustl(POSDATA.OUTVAR(IVARPLOT).name))//',MAX='//TRIM(ADJUSTL(STR1))//',MIN='//TRIM(ADJUSTL(STR2))
    

    
    IF(ABS(MINV)<1E-7) MINV=SIGN(1E-7,MINV)
    call loose_label(minv, maxv,init_num_contour,graphmin,graphmax,graphstep,contourbar.nfrac)
    graphstep=graphstep/scale_contour_num
    contourbar.nval=MAX(floor((graphmax-graphmin)/graphstep)+1,2)
    if(allocated(contourbar.VAL)) deallocate(contourbar.VAL)
    if(allocated(contourbar.color)) deallocate(contourbar.color)
    allocate(contourbar.VAL(contourbar.nval),contourbar.color(4,contourbar.nval))
    do i=1,contourbar.nval
        contourbar.VAL(i)=graphmin+(i-1)*graphstep
        !call get_rainbow(contourbar.VAL(i),graphmin,graphmax,contourbar.color(:,i))            
    enddo
end subroutine
end module function_plotter




