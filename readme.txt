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
!   property(10)=Mv !ѹ��ϵ��
!	property(11)=sita_s !���������ˮ��
!   property(12)=sita_r !���������ˮ��
!   property(13)=rw !ˮ���ض�
!	property(14)=zerothickness element thickness,the default value is 1.
!	Note that, for zt4_spg and zt6_spg, ʵ��Ϊһά����Ԫ��ֻ��k(ndimension)��ʾ�䴹ֱ�ڵ�Ԫ�߻������͸ϵ����property(14)Ϊ����ǽ��ȡ�	
!case(bar)
!	property(1)=E
!	property(2)=A
!	property(3)=hy !������y'(�ֲ�����)�ĸ߶ȣ�Ϊ����ת��Ϊ�����嵥Ԫʱ����
!	property(4)=hz !������z'(�ֲ�����)�ĸ߶ȣ�Ϊ����ת��Ϊ�����嵥Ԫʱ����
!	(5)=MinN,(6)=MaxN !��������ѹ���� For eip_bar

!case(beam)
!	property(1)=E
!	property(2)=A		
!	property(3)=u !poisson ratio
!	property(4)=J
!	property(5)=Iz !referred to the local coordinate system
!	property(6)=Iy !referred to the local coordinate system
!	property(7)=hy !������y'(�ֲ�����)�ĸ߶ȣ�Ϊ����ת��Ϊ�����嵥Ԫʱ����
!	property(8)=hz !������z'(�ֲ�����)�ĸ߶ȣ�Ϊ����ת��Ϊ�����嵥Ԫʱ����
!   (9)=MinN,(10)=MaxN,(11)=minMx,(12)=maxMx,(13)=minMy,(14)=maxMy,(15)=minMz,(16)=maxMz��(17)=C,(18)=PHI,(19)=yc  !���ڵ�����������������ھֲ��������ֵ,�ݲ����Ǽ�������ֵ��
!	��Ϊ��ԪSSP2Dָ�����ϣ��������(A,Iz,minN,maxN,minMz,maxMz)��Ϊ�����ְ�׮�Ĳ�����,C��PHI�ֱ�Ϊ���ؽӴ�Ħ�����ɵ�����������.
!case(pipe2,ppipe2)
!	property(1)=r !�ܰ뾶
!	property(2)=lamda !�ܱ�Ħ��ϵ��
!	property(3)=epslon !�ܱڵľ��Դֲڶ�
!	property(4)=v  !�˶�ճ��ϵ��
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
!	10)soilspring: (1 to 6)=ko,ka,kp,ks,��Ԫ���漰�ĳ���,Pw��
!	11)spring�� (2:5): minVal,maxVal,ks,�Լ���Ԫ���漰�ĳ��ȡ�
!	12)zt4_spg,ZT6_SPG: (1):AREA
