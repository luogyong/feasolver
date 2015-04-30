	!constants
	!INTEGER,PARAMETER::DPN=KIND(1.0D0)
	integer,parameter::inactive=-9999999
	integer,parameter::MNDOF=9
	integer,parameter::maximat=50,maxiet=500,maxilf=10,maxset=100
	integer,parameter::MaxNumRead=1000 !when input data, data numbers in each line is restricted to the MaxNumRead
												!Care should be taken when the following input way is used in SUB:strtoint: 1-5 or 1*100 
	!DEFINE PARAMETER 
	!element type
	integer,parameter::CONDUCT1D=1 ! One demenasional  CONDUCTIVITY element, the active DOF is x.
	!solid element
	integer,parameter::CPE3=2  !Continuum, plane strain, triangle element with 3 nodes, nodal number order is anticlockwise (1,2,3) 
	integer,parameter::CPE6=3  !Continuum Plane strain, triangle element with 6 nodes, nodal number order is anticlockwise (1,4,2,5,3,6) 
	integer,parameter::CPE4=4  !Continuum Plane strain, 4-node bilinear element, nodal number order is anticlockwise (1,4,2,5,3,6) 
	integer,parameter::CPE8=5  !Continuum Plane strain, 8-node biquadratic element, nodal number order is anticlockwise (1,4,2,5,3,6) 
	integer,parameter::CPE4R=6  !Continuum Plane strain, 4-node bilinear element, Reduced integeratoin,nodal number order is anticlockwise (1,4,2,5,3,6) 
	integer,parameter::CPE8R=7  !Continuum Plane strain, 8-node biquadratic element, Reduced integeratoin,nodal number order is anticlockwise (1,4,2,5,3,6) 
	integer,parameter::CPS3=8 ! plane stress
	integer,parameter::CPS6=9 ! plane stress
	integer,parameter::CPS8=10 !PLANE STRESS
	integer,parameter::CPS8R=11 !PLANE STRESS
	integer,parameter::CPS4=12 !PLANE STRESS
	integer,parameter::CPS4R=13 !PLANE STRESS	
	integer,parameter::CPE15=15	!PLANE STRAIN
	integer,parameter::CPS15=16	!PLANE STRESS
	integer,parameter::CAX3=17 !axissymmetric 3-noded triangle element
	integer,parameter::CAX6=18 !axissymmetric 6-noded element
	integer,parameter::CAX15=19 !axissymmetric 15-noded triangle element
	integer,parameter::CAX4=20 !axissymmetric 4-noded quadrilateral element
	integer,parameter::CAX4r=21 !axissymmetric 4-noded reduced quadrilateral element
	integer,parameter::CAX8=22 !axissymmetric 8-noded quadrilateral element
	integer,parameter::CAX8r=23 !axissymmetric 8-noded reduced quadrilateral element
	integer,parameter::PRM6=24  !6-noded prism element
	integer,parameter::PRM15=25  !15-noded prism element
	integer,parameter::PRM9=26 !6-noded prism element
	integer,parameter::tet4=27 !4-noded tetrahedron element
	integer,parameter::tet10=28 !10-noded tetrahedron element
	!structure element
	integer,parameter::BAR=301 !3D 2-noded bar element
	integer,parameter::BAR2D=311 !2D 2-noded bar element
	integer,parameter::BEAM=302 !3D 2-noded beam element
	integer,parameter::BEAM2D=312 !2D 2-noded beam element
	integer,parameter::DKT3=303 !3-noded Discrete Kirchihoff trianglual element
	integer,parameter::SHELL3=304 !18-dof flat triangular shell element introduced by C.A. Felippa .
	integer,parameter::SHELL3_KJB=305 !18-dof flat triangular shell element introduced by K.J. Bathe.
	integer,parameter::SPRINGX=321 !X方向的弹簧
	integer,parameter::SPRINGY=322
	integer,parameter::SPRINGZ=323
	integer,parameter::SPRINGMX=324
	integer,parameter::SPRINGMY=325
	integer,parameter::SPRINGMZ=326
	integer,parameter::soilspringx=327
	integer,parameter::soilspringy=328
	integer,parameter::soilspringz=329
	
	!钢板桩单元，目前存在假定如下：
	!1) 接触的两根钢板桩的材料，型号一样
	!2) 不考虑oblique bending
	integer,parameter::SSP2D=306 !2D 2-NODED STEEL SHEET PILE ELEMENT.
	integer,parameter::SSP2D1=307
	integer,parameter::pe_ssp2d=308

	
	!seepage element
	integer,parameter::CPE3_SPG=102  !Continuum, plane strain, triangle element with 3 nodes, nodal number order is anticlockwise (1,2,3) 
	integer,parameter::CPE6_SPG=103  !Continuum Plane strain, triangle element with 6 nodes, nodal number order is anticlockwise (1,4,2,5,3,6) 
	integer,parameter::CPE4_SPG=104  !Continuum Plane strain, 4-node bilinear element, nodal number order is anticlockwise (1,4,2,5,3,6) 
	integer,parameter::CPE8_SPG=105  !Continuum Plane strain, 8-node biquadratic element, nodal number order is anticlockwise (1,4,2,5,3,6) 
	integer,parameter::CPE4R_SPG=106  !Continuum Plane strain, 4-node bilinear element, Reduced integeratoin,nodal number order is anticlockwise (1,4,2,5,3,6) 
	integer,parameter::CPE8R_SPG=107  !Continuum Plane strain, 8-node biquadratic element, Reduced integeratoin,nodal number order is anticlockwise (1,4,2,5,3,6) 
	integer,parameter::CPS3_SPG=108 ! plane stress
	integer,parameter::CPS6_SPG=109 ! plane stress
	integer,parameter::CPS8_SPG=110 !PLANE STRESS
	integer,parameter::CPS8R_SPG=111 !PLANE STRESS
	integer,parameter::CPS4_SPG=112 !PLANE STRESS
	integer,parameter::CPS4R_SPG=113 !PLANE STRESS	
	integer,parameter::CPE15_SPG=115	!PLANE STRAIN
	integer,parameter::CPS15_SPG=116	!PLANE STRESS
	integer,parameter::CAX3_SPG=117 !axissymmetric 3-noded triangle element
	integer,parameter::CAX6_SPG=118 !axissymmetric 6-noded element
	integer,parameter::CAX15_SPG=119 !axissymmetric 15-noded triangle element
	integer,parameter::CAX4_SPG=120 !axissymmetric 4-noded quadrilateral element
	integer,parameter::CAX4r_SPG=121 !axissymmetric 4-noded reduced quadrilateral element
	integer,parameter::CAX8_SPG=122 !axissymmetric 8-noded quadrilateral element
	integer,parameter::CAX8r_SPG=123 !axissymmetric 8-noded reduced quadrilateral element
	integer,parameter::PRM6_SPG=124  !6-noded prism element
	integer,parameter::PRM15_SPG=125  !15-noded prism element	
	integer,parameter::tet4_spg=127 !4-noded tetrahedron element
	integer,parameter::tet10_spg=128 !10-noded tetrahedron element
	
	!pipe flow
	integer,parameter::pipe2=401 !2-noded line element for pipe flow simulation
	integer,parameter::ppipe2=402 !2-noded line element for perforated wellbore inflow simulation
	
	!Coupled element
	integer,parameter::CPE3_CPL=202  !Continuum, plane strain, triangle element with 3 nodes, nodal number order is anticlockwise (1,2,3) 
	integer,parameter::CPE6_CPL=203  !Continuum Plane strain, triangle element with 6 nodes, nodal number order is anticlockwise (1,4,2,5,3,6) 
	integer,parameter::CPE4_CPL=204  !Continuum Plane strain, 4-node bilinear element, nodal number order is anticlockwise (1,4,2,5,3,6) 
	integer,parameter::CPE8_CPL=205  !Continuum Plane strain, 8-node biquadratic element, nodal number order is anticlockwise (1,4,2,5,3,6) 
	integer,parameter::CPE4R_CPL=206  !Continuum Plane strain, 4-node bilinear element, Reduced integeratoin,nodal number order is anticlockwise (1,4,2,5,3,6) 
	integer,parameter::CPE8R_CPL=207  !Continuum Plane strain, 8-node biquadratic element, Reduced integeratoin,nodal number order is anticlockwise (1,4,2,5,3,6) 
	integer,parameter::CPS3_CPL=208 ! plane stress
	integer,parameter::CPS6_CPL=209 ! plane stress
	integer,parameter::CPS8_CPL=210 !PLANE STRESS
	integer,parameter::CPS8R_CPL=211 !PLANE STRESS
	integer,parameter::CPS4_CPL=212 !PLANE STRESS
	integer,parameter::CPS4R_CPL=213 !PLANE STRESS	
	integer,parameter::CPE15_CPL=215	!PLANE STRAIN
	integer,parameter::CPS15_CPL=216	!PLANE STRESS
	integer,parameter::CAX3_CPL=217 !axissymmetric 3-noded triangle element
	integer,parameter::CAX6_CPL=218 !axissymmetric 6-noded element
	integer,parameter::CAX15_CPL=219 !axissymmetric 15-noded triangle element
	integer,parameter::CAX4_CPL=220 !axissymmetric 4-noded quadrilateral element
	integer,parameter::CAX4r_CPL=221 !axissymmetric 4-noded reduced quadrilateral element
	integer,parameter::CAX8_CPL=222 !axissymmetric 8-noded quadrilateral element
	integer,parameter::CAX8r_CPL=223 !axissymmetric 8-noded reduced quadrilateral element
	integer,parameter::PRM6_CPL=224  !6-noded prism element
	integer,parameter::PRM15_CPL=225  !15-noded prism element	
	integer,parameter::tet4_cpl=227 !4-noded tetrahedron element
	integer,parameter::tet10_cpl=228 !10-noded tetrahedron element	
	
	integer,parameter::LBZT4=47	!4-node discontinuite element for low bound analysis
	integer,parameter::LB3=48	!3-node triangular element for low bound analysis 	
	integer,parameter::UBZT4=49     ! zero thickness discontinuite elements for upper bound analysis, the active dofs are vx vy
	integer,parameter::UB3=50  ! linear 2D trianglar elements for upper bound analysis, the active dofs are vx,vy.

	!integer,parameter::RCD1=501 !RIGID CONNECTION ELEMENT
	!integer,parameter::RCD2=502 !RIGID CONNECTION ELEMENT
	!integer,parameter::RCD3=503 !RIGID CONNECTION ELEMENT
	!integer,parameter::RCD4=504 !RIGID CONNECTION ELEMENT
	!integer,parameter::RCD5=505 !RIGID CONNECTION ELEMENT
	!integer,parameter::RCD6=506 !RIGID CONNECTION ELEMENT
	!integer,parameter::RCD7=507 !RIGID CONNECTION ELEMENT

	!coordinates for each node
	integer,parameter::LocX=1
	integer,parameter::LocY=2
	integer,parameter::LocZ=3
	!displacement for each node
	integer,parameter::DisX=4
	integer,parameter::DisY=5
	integer,parameter::DisZ=6
	!stress for each node
	integer,parameter::Sxx=7
	integer,parameter::Syy=8
	integer,parameter::Szz=9
	integer,parameter::Sxy=10
	integer,parameter::Syz=11
	integer,parameter::Szx=12
	!strain for each node	
	integer,parameter::Exx=13
	integer,parameter::Eyy=14
	integer,parameter::Ezz=15
	integer,parameter::Exy=16
	integer,parameter::Eyz=17
	integer,parameter::Ezx=18
	!plastic strain for each node
	integer,parameter::PExx=19
	integer,parameter::PEyy=20
	integer,parameter::PEzz=21
	integer,parameter::PExy=22
	integer,parameter::PEyz=23
	integer,parameter::PEzx=24
	!stress for each gaussian point
	integer,parameter::Sxxg=25
	integer,parameter::Syyg=26
	integer,parameter::Szzg=27
	integer,parameter::Sxyg=28
	integer,parameter::Syzg=29
	integer,parameter::Szxg=30
	!strain for each gaussian point
	integer,parameter::Exxg=31
	integer,parameter::Eyyg=32
	integer,parameter::Ezzg=33
	integer,parameter::Exyg=34
	integer,parameter::Eyzg=35
	integer,parameter::Ezxg=36
	!plastic strain for each gaussian point
	integer,parameter::PExxg=37
	integer,parameter::PEyyg=38
	integer,parameter::PEzzg=39
	integer,parameter::PExyg=40
	integer,parameter::PEyzg=41
	integer,parameter::PEzxg=42
	
	!plastic work for each element
	integer,parameter::PW=43
	integer,parameter::SIGMA_MISES=44	
	integer,parameter::EEQ=45
	integer,parameter::PEEQ=46
	integer,parameter::xf_out=47
	integer,parameter::yf_out=48
	integer,parameter::zf_out=49
	
	!ratation for each node
	integer,parameter::Rx=50
	integer,parameter::Ry=51
	integer,parameter::Rz=52
	!Moment for each node
	integer,parameter::Mx=53
	integer,parameter::My=54
	integer,parameter::Mz=55
	!shear force for each node
	integer,parameter::Qx=56
	integer,parameter::Qy=57
	integer,parameter::Qz=58
	
	!gradient for each node
	integer,parameter::Gradx=59
	integer,parameter::Grady=60
	integer,parameter::Gradz=61	
	!seepage velocity for each node
	integer,parameter::vx=62
	integer,parameter::vy=63
	integer,parameter::vz=64
	integer,parameter::Head=65
	integer,parameter::discharge=66
	integer,parameter::PHead=67
	integer,parameter::kr_spg=68
	integer,parameter::mw_spg=69
	
	integer::minet=10000,maxet=-10000   !the ultimate element type number. 
	!solver method
	integer,parameter::DIRECTI=1 !direct iterative algorithm, recompute the stiffness at each iteration
	integer,parameter::N_R=2 !Newton-Raphson method or tangential stiffness method,recompute the stiffness at each iteration
	integer,parameter::INISTIFF=3 !Initial stiffness method,computer the stiffness only at the beginning of the compution precedure
	integer,parameter::CINIATAN=4 !Combined initial and tangential stiffness approach, recompute the stiffness at 
								!the first iteration of each load increment only.
	integer,parameter::CINIATAN2=5 !Combined initial and tangential stiffness approach, recompute the stiffness at 
								!the second iteration of each load increment only.
	integer,parameter::LELASTIC=6 !LINEAR ELASTIC SOLVER
	integer,parameter::LA04=-2  !HSL LA04 SOLVER
	integer,parameter::LPSOLVER=-3 !USING LP_SOLVE TO SOLVE LIMIT ANALYSIS PROBLEM.
	integer,parameter::MOSEK=-4  !MOSEK OPTIMIZTION SOLVER.
	real(kind=DPN),parameter::UM=1.0d20

	!material type constants
	integer,parameter::ELASTIC=1 !elastic material
	integer,parameter::LA_MC=2 !Mohr_Coloumb material for limit analysis
	integer,parameter::CONDUCT=3 !
	integer,parameter::MISES=4
	integer,parameter::MC=5 !mohr-coloumb
	integer,parameter::CamClay=6 !modified Camclay
	integer,parameter::EIP_BAR=7 !弹理想塑性杆或弹簧。
	integer,parameter::EIP_BEAM=8 !弹理想塑性梁。
	
	!Seepage Parameters to control nonlinear conductivity 
	integer,parameter::step_spg=7 !step funtion >=0,kr=1; <0,kr=1e-6.
	integer,parameter::linear_spg=8 !default case. linear function, <-epsilon1,kr=1e-6;>epsilon2,kr=1; others, intepolating.
	integer,parameter::VG_SPG=9	!simulating the unsaturated soil behavior using the van genuchten model.
	integer,parameter::LR_SPG=10	!simulating the unsaturated soil behavior using the Leong and Rahardjo model.
	!Refer to Leong, E.C. and Rahardjo, H. (1997). “A review of soil-water characteristic curve equations”, J of Geotechnical and Geo-environment Engineering, 123(12), p1106-1117.
	integer,parameter::EXP_SPG=11	!simulating the unsaturated soil behavior using the exponent model.
	
	
	!node coordinates input format constants
	integer,parameter::POINT=1
	integer,parameter::BLOCK=2

	integer,parameter::YES=1
	integer,parameter::NO=0
	
	!element class
	integer,parameter::CPE=0	!plane strain solid
	integer,parameter::CPS=1	!plane stress solid
	integer,parameter::CAX=2	!axisymmetric solid
	integer,parameter::C3D=3	!3-d solid
	integer,parameter::CND=4 	!conduction element,conductive problem
	integer,parameter::SPG=5 	!3D SEEPAGE element,seepage problem
	integer,parameter::LMT=6 	!LIMIT ANALYSIS ELEMENT,limit analysis problem
	integer,parameter::SLD=7  !solid problem
	integer,parameter::CPL=8	!3D COUPLED ELEMENT,coupled problem
	integer,parameter::SPG2D=9 !2D SEEPAGE ELEMENT
	integer,parameter::CPL2D=10 !2D COUPLED ELEMENT
	integer,parameter::STRU=11 !STRUCTURE ELEMENT
	integer,parameter::CAX_SPG=12 !AXIAL SYMMETRICAL SEEPAGE PROBLEM.
	integer,parameter::PIPE=13 !PIPE FLOW ELEMENT
	integer,parameter::PE=14 !Penelty element
	integer,parameter::spring=15 !spring like element
	integer,parameter::soilspring=16 !soil spring
	integer,parameter::CAX_CPL=17
	
	!method for mapping the values in integration points to nodes
	integer,parameter::SHT=1
	integer,parameter::SPR=2
	integer,parameter::AVG=3
	integer,parameter::EXPO=4
	
	!parameters for plastic body force generation method
	integer,parameter::VISCOP=1 !viscoplasticity method
	integer,parameter::INISTRESS=2 !INITIAL STRESS MEHTOD 
	integer,parameter::consistent=3 !fully implicit constitutive integration scheme and consistent tangent modulus
	integer,parameter::continuum=4 !semi-implicit constitutive integration scheme and continuum tangent modulus
	integer,parameter::constant=5
	integer,parameter::iniflux=6 !for spg analysis, bathe method.
	integer,parameter::lacy=7	!for spg analysis, lacy and Jean method.
 	!parameters for geostatic initial stress calculation
	integer,parameter::CAL_GEO=1
	integer,parameter::KO_GEO=2
	!load type
!	integer,parameter::NF=1
!	integer,parameter::BF=2
	
	!special value
	integer,parameter::SV=-9999999
	
	!system constant
	integer,parameter::sys_cylinder=-999
	integer,parameter::sys_sphere=-9999
	integer,parameter::sys_local=-99 !for bar family element.
	
	
	!bc or load partern
	integer,parameter::ramp=1 !在一步内随时间线性施加。
	integer,parameter::step=2 !在步初，瞬间施加全部的量。
	integer,parameter::ReRamp=-1 !在一步内随时间线性减小。
	
	
	
