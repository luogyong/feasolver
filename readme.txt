!coordinate system
!plane strain (x,y), y upward is positive and x rightward is positive.

!1. Degrees of freedom

!Except for axisymmetric elements, the degrees of freedom are always referred to as follows:
!1 x-displacement /limit analysis vx
!2 y-displacement /limit analysis vy
!3 z-displacement 
!4 Pore pressure, hydrostatic fluid pressure, scale field variable
!5 Rotation about the x-axis, in radians
!6 Rotation about the y-axis, in radians
!7 Rotation about the z-axis, in radians
!Here the x-, y-, and z-directions coincide with the global X-, Y-, and Z-directions, respectively;

!2 Material.property
!case (CONDUCT)
!		property(1)=k0
!		property(2)=k1  
!		the conductivity k=k0+k1*field variable
!case (LA_MC)
!		property(1)=p
!		property(2)=phi
!		property(3)=c
!		property(4)=Gama, specified gravity of soil 
!		property(5)=pore pressure 
!case(Elastic)
!	property(1)=E
!	property(2)=v
!case(mises)
!	property(1)=E
!	property(2)=v
!	property(3)=sigma_y, yield stress under uniaxial stretch. It should be noted that , for clay under plain strain
!										sigma_y=3*Cu (Cu, undrained shear strength) and sigma_y=2*Cu for triaxial (axisymmetrical) state.
!	propery(4)=unit weight
!case(MC)
!	property(1)=E
!	property(2)=v
!	property(3)=cohesive strength
!	property(4)=friction angle
!	porperty(5)=dilation angle
!	property(6)=unit weight	
!case(camclay)
!	property(1)=M
!	property(2)=v
!	property(3)=lamda
!	property(4)=kapa
!	porperty(5)=
!	property(6)=
!case(spg)
!	property(1)=k1
!			(2)=k2
!			(3)=k3
!			(4)=INT(TransM) ! !for SPG, transform the material kij to globle kxy, L2G.
!	property(7)=alpha !for van genuchten model/Leong and Rahardjo model
!	property(8)=n !for van genuchten model/Leong and Rahardjo
!   property(9)=m !for Leong and Rahardjo model
!   property(10)=Mv !压缩系数
!	property(11)=sita_s !饱和体积含水量
!   property(12)=sita_r !残余体积含水量
!   property(13)=rw !水的重度
!	property(14)=zerothickness element thickness,the default value is 1.
!	Note that, for zt4_spg and zt6_spg, 实际为一维流单元，只用k(ndimension)表示其垂直于单元边或面的渗透系数，property(14)为防渗墙厚度。	
!case(bar)
!	property(1)=E
!	property(2)=A
!	property(3)=hy !梁截面y'(局部坐标)的高度，为后处理转化为六面体单元时所用
!	property(4)=hz !梁截面z'(局部坐标)的高度，为后处理转化为六面体单元时所用
!	(5)=MinN,(6)=MaxN !最大的轴向压力， For eip_bar

!case(beam)
!	property(1)=E
!	property(2)=A		
!	property(3)=u !poisson ratio
!	property(4)=J
!	property(5)=Iz !referred to the local coordinate system
!	property(6)=Iy !referred to the local coordinate system
!	property(7)=hy !梁截面y'(局部坐标)的高度，为后处理转化为六面体单元时所用
!	property(8)=hz !梁截面z'(局部坐标)的高度，为后处理转化为六面体单元时所用
!   (9)=MinN,(10)=MaxN,(11)=minMx,(12)=maxMx,(13)=minMy,(14)=maxMy,(15)=minMz,(16)=maxMz，(17)=C,(18)=PHI,(19)=yc  !对于弹理想塑性梁，相对于局部坐标的限值,暂不考虑剪力的限值。
!	当为单元SSP2D指定材料，输入参数(A,Iz,minN,maxN,minMz,maxMz)均为单根钢板桩的参数，,C和PHI分别为库仑接触摩擦定律的两个参数。.
!case(pipe2,ppipe2)
!	property(1)=r !管半径
!	property(2)=lamda !管壁摩阻系数
!	property(3)=epslon !管壁的绝对粗糙度
!	property(4)=v  !运动粘滞系数
!3 element.type
!	1,=conduct1d
!	-23,=ubtri2dl
!	-24, =ubzt2d
!	2,cpetri3n
!	3,cpetri6n
!	
!4 element.property
!   1) ub3: property(1)=2*Area_the element
!						.property(2)=dissipation work
!						.property(3)=centroid coordinates yc
!   2) ubzt4: property(1)=Length of the element  
!					.property(2)=dissipation work
!					.property(3)=centroid coordinates yc
!	3) plain strain element (spe series)
!		.property(1)=E*(1-v)/((1+v)*(1-2*v))
!		
!	4) plain stress element (sps series)
!		.property(1)=E/(1-v**2)
!	5) bar element
!	.propterty(1)=1. !
!    .property(2)=L
!	6) beam element,ssp2d
!	.property(1)=1.
!    .property(2)=L
!	7) DKT3,SHELL3,SHELL3_KJB
!	.PROPERTY(3)=THICKNESS
!	.property(1)=1.0 
!	.PROPERTY(2)=ELEMENT AREA
!	8) cpe3_spg
!		.property(1)=1.0
!		.property(2)=element area
!	9)pe_ssp2d: (1)=1,fixed(default); =2,slip;=0:free
!	10)soilspring: (1 to 6)=ko,ka,kp,ks,单元所涉及的长度,Pw。
!	11)spring： (2:5): minVal,maxVal,ks,以及单元所涉及的长度。
!	12)zt4_spg,ZT6_SPG: (1):AREA
