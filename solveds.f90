  
    !��������������ݽṹ
	!�ڵ�
module solverds

    
	include 'double.h'
	include 'constant.h'   
    
    CHARACTER(1024)::INPUTFILE=''
    
	type node_tydef
		real(kind=DPN)::coord(3)=0.0D0 !coordinates (x,y,z)
		integer::ndof=0 !�ڵ�����ɶȸ���
		integer::dof(MNDOF)=inactive !-9999999,inactive dof; >0:active dof. if dof()>0, then it indexes the dof number.						 		
		real(kind=DPN),allocatable::stress(:),strain(:),pstrain(:),SFR(:),PSIGMA(:)  !SFR=Ӧ���ƻ���
		real(kind=DPN),allocatable::FQ(:),M(:) !nodal shear forces and nodal bending moments.
		real(kind=DPN),allocatable::igrad(:),velocity(:)	!for spg element, gradients,velocities.
		real(kind=DPN)::angle=0.0D0 !the angle around the vertex.
		real(kind=DPN)::MISES=0.0D0,EEQ=0.0D0,PEEQ=0.0D0 !for poreflow model. Mises=clogging_poresize,eeq=free_cc
		real(kind=DPN)::Q=0.0D0,Kr=0.0D0,Mw=0.0D0,PORESIZE=0.D0,cc=0.d0  !nodal flux,relative permeability,slope of volumetric water content.		
		integer::nelist=0,NELIST_SPG=0,ISACTIVE=0
		!integer,allocatable::elist(:),ELIST_SPG(:) !elements sharing the node.
!		integer::Property=0 !for SPG, Property=1 suggesting that the node is on the seepage surface.

	end type
	integer::nnum !�ڵ���
	type(node_tydef),allocatable::node(:)
	real(DPN),ALLOCATABLE::GNODE(:,:) !GHOST NODE
	INTEGER::NGNODE=0
    
	type element_tydef
		integer::nnum,NEDGE=0,NFACE=0,NTET=0,ISTOPO=0 !node numbers of this element
		integer,allocatable::node(:),EDGE(:),FACE(:),TET(:),ADJELT(:) !��Ԫ�Ľڵ�,TETΪ��Ԫ��ϸ�ֵ�Ԫ����±ꡣ
        integer,allocatable::node2(:),NSPLOC(:) !Ϊ����Bar��Beam��Ԫ�ĺ���Ϊ��Ԫ���ڵĽڵ��ţ�����ת����ʵ�������嵥Ԫ�����,����ԪΪzt4_spg,��zt6_spgʱ��node2ָ��gnode��
        !��et=wellboreʱ��node2�洢��wellbore��Ԫ��3��4�ڵ�Ϊ�ߵ�3ά��Ԫ,NSPLOCΪ���ܲ��������ڵĵ�Ԫλ�� 
        !��ISTOPO==1ʱ���ɸߴε�Ԫ���ɵ�һ�ε�wellbore�൥Ԫ����������Ԫ����ԪΪ�ߴε�Ԫʱ����wellbore�൥ԪĿǰֻ��2�ڵ��һ�ε�Ԫ�����ԶԸߴ��ߵ�Ԫ�����˷ֽ⣬
        !node(nnum+1:nnum+2)��¼��ǰ��Ԫ��Ӧ�ĸߴε�Ԫ�Ķ˽ڵ㣬�Ա���к���������ڽӷ�����        
		integer::et  !��Ԫ����
		integer::mat,mattype  !material id and material type.the paramters is got from material(mat)
		!for et=soilspring, mat=-1,�����൥Ԫ��mat=-2,�����൥Ԫ
		integer::set=0 !element set
		integer::ndof !total dofs of the element
		integer::ngp !the number of gauss points		
		integer::isactive=1 !1Y0N��-1,COUPLEELEMET,�����м��㣬�������һ���Ľ����Ϊ�����Ľ����
        INTEGER::SF=0 !STEP FUNCTION .FACTOR=0,DEATICE;FACTOR=1,ACTIVE
		integer::sign=1 !for soilspring element .sign=1, Pa,Po,Pp,Pw=+; sign=-1, pa,po,pp,pw=-.
					!for pe_ssp2d,ngp=iָ��smnp(i).		
		integer::nd ! the dimension of the strain-stress matrix
		integer::id  !element id number in the set
		integer::ec=0,eshape=-1 !element class,shape
		integer::nspr=0 ! the patch number sharing the element
		integer::layer=1 !element layer, for horizontal seepage analysis; 
						!for sectional seepage analysis, =1, all nodes are under water table,=-1,all nodes are above water talbe,=0 cross the water table
		integer::ifreedof=0 !>0,��½ӵ�,ָ��freedof(ifreedof)
		                    !for poreflow element,particle number clogging in the throat. 
		!integer::status=1 !=1,fixed(default); =2,slip;=0:free.
		integer::referencestep=1 !=i,�����˵�Ԫ���μ����ʼ�ο�״̬���Ϊ��i��������λ�Ƴ���        
		integer::system=0 !local coordinate system, =0,global coordinate
		! local coordinate system for bar element:
		! x-axis is along the bar and the positive direction is from node 1 to node 2
		! y and z is trivial for bar element.
		! local coordinate systme for beam element:
		! x-axis is along the bar and the positive direction is from node 1 to node 2
		! y and z is defined by user but must be consistent with right hand rule.
		!and be consistent with Iy and Iz.
		real(kind=DPN)::property(6)=0.0D0,cc=0.d0  !for spg problem and iniflux is used, property(3)=low lamda,property(2)=up lamda
        !for wellbore element, property(1), element frictional resistance,(2) and (3) are geometrical resistance; (4) acceralated resistantce; (5)=WELL SKIN RESISTANCE , property(6) surround angle.
        !fore sphflow and semi_sphflow property(1)= geometrical resistance.for semi_sphflow,property(4-6)=domain direction vector
        !for pipe2/poreflow  element,property(1), element frictional resistance,(2)D2, throat diameter(for poreflow),(3) D1,(4)=Length of the throat.(5)=length of clogging volume,(6)K of the clogging volume
		real(kind=DPN),allocatable::angle(:)!internal angle for every nodes
        !��et=wellboreʱ,angle�洢node2��Ԫ��Ӧ�Ķ���ǡ�
		integer,allocatable::g(:) !��Ԫ�Ķ�λ����
		real(kind=DPN),allocatable::km(:,:) !��Ԫ�ĵ���,km(ndof,ndof)
		real(kind=DPN),allocatable::CMM(:,:) !Consistent mass matrix
		real(kind=DPN),allocatable::B(:,:,:)  !B matrix, B(nd,ndof,ngp+nnum)
		real(kind=DPN),allocatable::D(:,:)  !elastic strain-stress Matrix, D(nd,nd), its shape depends on the
										!for element cpe3_spg D stores the initial km
        real(kind=DPN),allocatable::MHK(:,:) !couple matrix
		real(kind=DPN),allocatable::detjac(:) ! the determinant of jacobian matrix in each gauss point. 
		real(kind=DPN),allocatable::stress(:,:),Dstress(:,:) !stress at gauss points and nodes, stress(6,ngp+nnum)
		real(kind=DPN),allocatable::strain(:,:),Dstrain(:,:) !strain at gauss points and nodes, strain(6,ngp+nnum)
														!for cpe3_spg, strain(2,2),strain(:,1)=gradient, strain(:,2)=velocity
		real(kind=DPN),allocatable::sfr(:,:) !1.stress failure ratio,2.failure plane angle respective to x, anticlockwise is +ve.
											 !3. sigman(tension is +ve ) 4 tn (clockwise is +ve),5. tnx, 6. tny,7, sfr_kR,8,SFR_KO 
		real(kind=DPN),allocatable::pstrain(:,:) !total plastic strain in gp and centroid
		real(kind=DPN),allocatable::evp(:,:) !plastic strain incremental of each incremental in gp and centroid.
		real(kind=DPN),allocatable::xygp(:,:) !globe co-ordinates of gausian points.
		real(kind=DPN),allocatable::ev(:) !evolution variable in each gausian point.
		real(kind=DPN),allocatable::e(:) !void ratio of each gausian point
        real(kind=DPN),allocatable::uw(:) !pore pressure at each gausian point
		real(kind=DPN),allocatable::GForce(:),DGforce(:) !the cumulative general force of each dof in globle system,Gf(ndof),useful for structure elements
														!the nodal discharges, for spg element.
                                                        !Dgforce(:),ÿ������������
		real(kind=DPN),allocatable::GForceILS(:) !the nodal force in the local system.
        !real(kind=DPN)::reference
 		real(kind=DPN),allocatable::G2L(:,:) !transformation matrix from global system to local system
		!For SPG elements
		real(kind=DPN),allocatable::igrad(:,:) !hydraulic gradient igrad(nd,ngp+nnum)
		real(kind=DPN),allocatable::velocity(:,:) !seepage velocity, velocity(nd,ngp+nnum)
		real(kind=DPN),allocatable::flux(:) !flux
		real(kind=DPN),allocatable::Kr(:) !relative permeability
		real(kind=DPN),allocatable::Mw(:) !slope of volumetric water content
		real(kind=DPN),allocatable::sita_fin(:),sita_ini(:) !volume water content in the end and start of the time step.  
        !INTEGER,ALLOCATABLE::WELL_SP3(:),WELL_SP4(:) !USE forlab WELLBORE ELEMENT
        !REAL(KIND=DPN),ALLOCATABLE::WELL_SP3_R(:),WELL_SP4_R(:)
        !INTEGER::NWSP3=0,NWSP4=0		
		
		!real(kind=DPN),allocatable::lamda(:) !lamda(ngp)=1.0(by default),relative k. 
		!additional variables for upper bound analysis
		real(kind=DPN),allocatable::A12(:,:)  ! for UBZT4 is A23;WELLBORE,SAMPLE POINTS
		real(kind=DPN),allocatable::X2(:)
        REAL(8),allocatable::PFP(:) !(1-5):D1,D2,LR,Vol1,Vol2�ֱ��Ӧ�ڵ�12��Բֱ̨�������λ��,(6,7)Lc1,Lc2�ֱ�Ϊ�ٶ���ĳ���,(8,9)npc1,npc2�ٶµĿ�����
        REAL(8)::BBOX(2,3)=0.D0 !MIN,MAX OF X,Y AND Z
        REAL(8)::fqw=1.0d0,FD=0.D0 !���ܵ�Ԫ��͸ϵ����������.FD=FRICTIONAL TERM FOR PIPE2/poreflow
		REAL(8)::cdt=1.0d0 !poreflow convective term flow integration time
		
	contains
		procedure::RIJ=>PNW_FLOW_RIJ
	end type
	integer::enum=0 !�ڵ���
	type(element_tydef),allocatable::element(:)
	
	
	!step function
	type stepfun_tydef
		real(kind=DPN),allocatable::factor(:)
        character(64)::title=""  
        integer::base=1 !��base=1ʱ��factor(0)=0����base=0,factor(0)=input by user.
		!ע��,������Ǹ������ػ�λ�Ʊ߽������.
		!��FACTOR(ISTEP)=-999ʱ,�˱߽������ڴ˲���ʧЧ�������ã�.	
        !��Ԫ����ʱ��1Ϊ��,����δ����
	endtype
	type(stepfun_tydef),allocatable::sf(:)
	integer::nsf=0
	integer::nstep=1	
	!the value of the default step function sf(0)
	!sf(0).factor(0)=0 !the 0 step is for geotechnical initial stress calculation.at this step
	!all  boundaries and loads are 0 exept for the gravity.
	!sf(0).factor(1:nstep)=1.0
	
	type subtimestep_tydef
		integer::nsubts=0
		real(kind=DPN),allocatable::subts(:)
		!real(kind=DPN),allocatable::subfactor(:)
	end type
	type(subtimestep_tydef),allocatable::TimeStep(:)
!	real(kind=DPN),allocatable::TimeStep(:) !timeStep(nstep) !incremental time for each step.
	integer::nts=1
    
    type Stepinfo_tydef
        integer::matherstep=0
        logical::issteady=.true.
        integer::bctype=step,loadtype=step
        
        !LOADTYPE(BCTYPE):Ϊ����(λ�Ʊ߽�)ʩ�ӷ�ʽ��
		!=1(ramp),��ʾ���ں���(λ�Ʊ߽�)��ʱ������ʩ�ӣ�
		!=2(step(default))����ʾ������(λ�Ʊ߽�)�ڲ���˲��ʩ�ӡ�
		!=-1(ReRamp) ��ʾ������(λ�Ʊ߽�)��ʱ�����ԴӴ��С��
	end type
    type(Stepinfo_tydef),allocatable::stepinfo(:)
    integer::nstepinfo=1
	
!	type IniValue_tydef
!		real(kind=DPN),allocatable::v(:) !��ֵ. initialvalue(nnode)
!	end type
!	type(IniValue_tydef)::InitialValue(MNDOF)
	
	
	!boundary condition data structure
	type bc_tydef
		integer::node,dof	!when body force is input, node equels to element	
		integer::sf=0,rtf=0 !step function for this load. 
		integer::isdead=0 !for =1,the condition is deactive in the current step.
		real(kind=DPN)::value=0
				
		!����bc_disp,���isdual==i(>0),���ʾ�����ɶȿ��������߽�Nseep(i)�ظ�������߽�ˮͷС��λ��ˮͷ�����Ϊ����߽硣
		!����Nseep, ���isdual==i(>0)�����ʾ�����ɶȿ�����ˮͷ�߽�BC_disp(i)�ظ�������߽�ˮͷ����λ��ˮͷ������ˮͷ�߽�Ϊ׼����ˮͷ�߽�����ȼ����ڳ���߽硣
		integer::isdual=0,iswellhead=0 
        !����bc_disp,iswellhead>0,��ʾ�˽ڵ�Ϊ��ѹ�������߽�ڵ㣬�ڵ�����ֻ�������������������������˱߽�ʧЧ
        !����Nseep,iswellhead>0,��ʾ�˽ڵ�Ϊ���ڳ����棬�Ҵ˳�����������������Ϊiswellhead�ľ�������
		integer::ssp_onepile=0 !ֻ��SSP��Ԫ�ڵ������ã���ʾ��������Ƿ�������������һ���ְ�׮��,�������������á�
								!/=0,�����ڵ����ְ�׮�ϣ�=0�����á�
		integer::isincrement=0 !���Ϊ1�������value��������������ȫ����Ĭ�ϣ�
        
	end type
	type(bc_tydef),allocatable::bc_disp(:),bc_load(:),bf(:),NSeep(:),IniValue(:),CFN(:),ccbc(:)
	integer::bd_num=0,bl_num=0,bfnum=0,NumNSeep=0,Niniv=0,NCFN=0,nccbc=0	
	real(kind=DPN),allocatable::iniValueDof(:) 
    
    TYPE QWELLNODE_TYDEF
        INTEGER::NNODE,DOF=4,NNODE2=0
        INTEGER,ALLOCATABLE::NODE(:),NODE2(:),ELEMENT(:) 
        REAL(kind=DPN),ALLOCATABLE::QAN(:,:) !���ڵ�node2(i)�Ľ�������qan(1,i)����ֵ����qan(2,i)
        !NODEΪÿ�ھ��ľ��ڽڵ㼰����ڵ㣬��ͳ�ƾ�������
        !node2Ϊÿ�ھ��ľ�����Ԫ�����������浥Ԫ�ĵ�3��4�Žڵ�  
        !ELEMENTΪ���ܸ�����ʵ�嵥Ԫ����״Ϊ203��304�ĵ�Ԫ��
        real(kind=DPN)::Q=0,QA=0,QN=0 !
    ENDTYPE
    TYPE(QWELLNODE_TYDEF),ALLOCATABLE::QWELLNODE(:),QWELLNODE1(:)
    REAL(kind=DPN),ALLOCATABLE::QWAN(:,:)   !���ڵ�node2(i)�Ľ�������qan(1,i)����ֵ����qan(2,i)
    INTEGER::NQWNODE=0
	LOGICAL::IS_WELL_GEO_INI=.FALSE.
    
    type hinge_typef
        integer::element,node,dof=7 !��Ԫ�ţ��ڵ�ţ����ɶ�
		integer::newnode=0 !ָ����ڵ�NodeΪ�غϵĽڵ�node(newnode)��
    endtype
    type(hinge_typef),allocatable::FreeDOF(:)
    integer::Nfreedof=0
	
	integer,allocatable::IsBCChange(:)
		

	!material 
	type mat_tydef
		integer::type=0,ISINPUT=0
		real(kind=DPN)::property(32)=0.0D0
		integer::sf(32)=0 !not considered yet.
		real(kind=DPN)::weight=0.0
!		real(kind=DPN)::dt !stress integration time step for viscoplasticity algorithm. 
		integer::FF1D(32)=0 !DEPENDENT FIELD FUNCTION,the default determined variable  is the centroid coordinates (Y)
		logical::isff=.false. ! whether the parameters are dependent on  the field function
		character(64)::name=""
        REAL(kind=DPN),ALLOCATABLE::D(:,:),DINV(:,:) !ELASTIC K AND INVERSE(K)
		CONTAINS
			PROCEDURE::GET_ARRAY=>GET_MAT_PROPERTY_ARRAY
			PROCEDURE::GET_SCALAR=>GET_MAT_PROPERTY_SCALAR
			GENERIC::GET=>GET_ARRAY,GET_SCALAR
			!
	end type
	type(mat_tydef)::material(-2:maximat) 
	
	
	!solution control
	type solver_tydef
		integer::type=SLD !problem type,=SLD,solid;=SPG,seepage;=CPL,coupled.
		integer::solver=N_R !solution method, =-2, upper bound analysis,the default solver is N_R
		integer::bfgm=continuum,bfgm_spg=continuum  !body force generation method(stress update algorithm) if the INISTIFF method is applied.
		integer::niteration=100  !Maximum number of iterations allowed for each increment
        integer::NYITER=20 !Maximum number of iterations allowed for STRESS RETURN
		real(kind=DPN)::tolerance=0.001D0 !Convergence tolenance factor.
		real(kind=DPN)::Yftol=1.D-3 !tolerance permitted for yiled criterion to be vilated
		real(kind=DPN)::force_tol=1.D-3 !tolerance permitted for unbalanced residual force to be vilated
		real(kind=DPN)::disp_tol=1.D-3 !tolerance permitted for residual displacement to be vilated
        real(kind=dpn)::slowtol=0.001 !when convergence is slow when convratio<slowtol,exit
		integer::output=1 !1:output for mod(increments,output)=0 iteration
		integer::i2ncal=4,I2NWEIGHT=1  !method mapping variables in integration points to nodes
			!SHT=1,spr=2,AVG=3,EXP=4
		logical::issym=.true.
!		logical::issteady=.true. 
		logical::datapacking=.true. !.t.=point,.f.=block
		logical::ismg=.false. !just for limit analysis, if ismg=.true. minimize soil specific gravity, it demands that no nonzero boundary are applied.
		logical::islaverify=.false. !.true. to verify a soluction is whether a feasible one.
		logical::ispg=.false.  ! whether the job is just to read in result file and out put to tecplot.
		logical::isfc=.true. !whether the residual force criteria is implemented in judge of the convergence.
        LOGICAL::ISCONVERGING=.TRUE.
		!real(kind=DPN)::epsilon_spg=0.1 !for spg only. control the conductivity slope
!		integer::Para_spg=linear_spg !step_spg,linear_spg,unsat_spg.
									!if Para_spg=unsat_spg, van genuchten model is used. and 
									!alpha=material(element().mat).property(7)
									!n=material(element().mat).property(8)
		logical::ISFU=.false.  !For seepage problem, =yes, use factorization update to update the boundary condition change on the seepage face during iteration.
											!=NO, refactorization completely.
		logical::ismkl=.true.
        integer::mur=1
		real(kind=DPN)::BarfamilyScale=15.D0
        !parameter for line search
        INTEGER::ISLS=0,ISACC=DANG  !ISACC=0,NO ACCELERATION, =1,SLOAN, =2,DANG(DEFAULT)
        INTEGER::NLS=10
        real(kind=DPN)::STOL=0.8D0,S0=0.D0,ALPHA=1.0D0 !ALPHA =INISTIFFNESS ACCELERATION PARAMETER
        integer::RF_EPP=0 !�����������ɿ�����ֵ�Ƿ�Ҫ������ʼ����ѹ����
        integer::RF_APP=0 !������������ѹ�����أ������������Ƿ񰴵������ۼ���
		integer::INIEPP=2  !�����������ɿ�����ֵ�Ƿ�Ҫ������ʼ����ѹ��,2=������ѹ����1=��ֹ��ѹ��
        integer::nopopup=0
        integer::isParasys=0,CaseID=0 !isParasys,�Ƿ�Ϊ���������Է���(must start form 1)
!		integer::isPostCal=0 !���е�δ֪����Ϊ��֪���ɱ߽��������룩�������к�����㡣
        !REAL(KIND=DPN),ALLOCATABLE::ETA(:),RATIO(:)
        integer::slidedirection=right,slope_isTensionCrack=1,ISSLOPEPA=0,slope_mko=1,slope_ONLY_SEARCHTOP=0
        real(kind=DPN)::slope_kscale=2.D0,slope_kbase=-0.2,slope_kratio=10,slope_sfrpeak=1.1,slope_guide_direction=0 
        !���slope_kbase<0,��kbȡ���ֵ��kb=abs(slope_kbase)*maxsfr;����kbȡ����ֵ��kb=slope_kbase
        !sfrpeak=a,��ʾ���ڵ�ĳ���sfr/maxsfr>=aʱ������͸ϵ����exp(2*sfr)���,ּ��������������������ƻ���,a>1�����Ǵ˼�ǿ,�����ֻᵽ�����������㲻�ȶ���������Ϊ1.1�������ã���
        !slope_mko>0,��sfrҪ������ʼkoӦ������sfrko.�ٶ�ko=v/1-v,sxx=ko*syy,szz=sxx,txy=txz=tyz=0
        !IS_ONLY_SEARCHTOP/=0,����streamline�������в�������ʱ������ʹ������������µĶ����������������ϲ���С������
        !if(slope_kratio>0),ky=kx/slope_kratio,else,ky=input value.
		!slope_guide_direction==0,������ȷ�����Ż�������;==1,��λ�Ƴ�ȷ��;==2,�����������������ο�����ͶӰ�����ߵļн�ȷ��.
        integer::wellmethod=3 !�������������ʱ��=0����ȡ����Ϊ��Ԫ�ڵ㣻=1��Ϊ��Ԫ����;=2,�ڵ�Ԫ�ܱ߾���3�������(���˼��м�),ÿ�Ų�������Ϊnspwell,��״����Ϊ��Ԫ�ڵ㡣>2,��Ϊ������
        integer::nspwell=12 !wellmethod=2ʱ,ÿ��(2*Pi)����������������ܽǷ�2PI,��Χ֮��ĵ㽫������
        integer::wellaniso=0 !ˮƽ�������ԣ�=0��directional K method�� =1, Charles R. Fitts method.(ת��Ϊ����ͬ�Բ��Ͻ���)��=2�������ٷ����������乫ʽ���м��� 
        real(kind=DPN)::disf_scale=1.0d0
        integer::len_unit=0 !0,m;3,mm;2,cm;1,dm;4,km;
        integer::time_unit=0 !0,day;1,sec;2,min;3,hour;
        real(kind=DPN)::cdt=1.d0 !concentration integerate time
		integer::pnw_clogging=0 !if>0
		integer::well_bottom_type=0 !=0,ƽ�׾�(Ĭ��);<>0,����Ϊ�����
        !integer::well_bottom_method=0 !=0,������Ԫ��(�ǵ���,Ĭ��);=1,����
    contains
        procedure::unit_factor=>unit_scaling_factor
        procedure::get_g=>get_gravity
	end type
	type(solver_tydef)::solver_control
	INTEGER::MAX_NODE_ADJ=50,MAX_FACE_ADJ=100
    
	type property_tydef
		character(128)::name=''
        integer::nval=1
		real(kind=DPN)::value=0.0
		character(512)::cvalue=''	 !character value
    contains
        !procedure split
	end type
	type(property_tydef)::property(maxset) !one control line has MAXSET property value at most.
	integer::pro_num 

	type eset_tydef !for output
		!integer::num=0 !the name of the set
		character(16)::Stype  !for tecplot FETRIANGLE, OR FEQUADRILATERAL
		integer::enums=0,enume=0 !the first and the last element number in the element() of the set. 
		integer::ec
		integer::et
        integer::eshape=-1 !0,point��101=line,203=triangle,204=quadrilateral,304=tet,308=hex,306=prism
        integer::system=0
		character(1024)::zonetitle="",zonetitle_sf="",zonetitle_pf=''
		character(64)::grouptitle=""
        integer::coupleset=-1
        integer::isini=0 !>0
        integer::sf=0 !stepfun
		integer::mesh_share_id=0,out_mesh=.true.!for output tecplot.
		integer::mesh_share_id_sf=0,out_mesh_sf=.true.
		integer::mesh_share_id_pf=0,out_mesh_pf=.true.
		!for bar ane beam element only.
		real(kind=DPN),allocatable::xyz_section(:,:) !��Ԫ����ͳһ�Ľ�����ĸ��ǵ�����꣬�ּٶ�һ����Ԫ���ڵ����и˵�Ԫ������Ԫ�ľֲ�����һ����
		integer,allocatable::outorder(:) !�������ʱ��ͬһ�ֲ������½ڵ�����˳����element.node2�ڵ�Ŷ�Ӧ��
        integer::noutorder=0
		integer,allocatable::elist(:,:) !��Ԫ���ڵ�elist(2,noutorder)��Ŀǰ�ٶ�һ���ڵ������2����ϵ��Ԫ��

    end type
	type(eset_tydef)::eset(maxset) !the maximum set number allowed is MAXSET.
	integer::esetid(maxset)=-1,neset=0  !set number in the model.	
	
	type ueset_type !for boundary conditions apply.
		integer::enum=0 !element number in the set
		character(32)::name=''
		integer,allocatable::element(:) !element index of the set
		integer::group=-1				
	end type
	type(ueset_type)::ueset(0:maxset) !ueset(0)=all elements in this model
	integer::nueset=1
	
	type unset_type
		integer::nnum=0 !node number in the set
		character(32)::name=''
		integer,allocatable::node(:) !node index of the set		
	end type
	type(unset_type)::unset(0:maxset) !unset(0)=all nodes in this model
	integer::nunset=1	
	
	!element class property
	type ECP_tydef
		integer::nshape,ngp,ndim,nnum,shtype  !number of the shape function,number of the gaussian integer point, space dimensions,number of nodes
		logical::isini=.false.  ! whether the element type is initialized or not.
		real(kind=DPN),allocatable::gp(:,:),weight(:) !gaussina point local coordinates and their corresponding weights
							!gp(ndim,ngp), the last nnum entries store the local coordinates in element nodal points.
		real(kind=DPN),allocatable::Lderiv(:,:,:)  ! Derivative of the shape functions with the local coordinates at gaussian points.
						  !Lderiv(nshape,ndim,ngp)
		real(kind=DPN),allocatable::Lshape(:,:)	! shape function value at gaussian points
							!Lshape(nshape,ngp)	
		real(kind=DPN),allocatable::expolating_Lshape(:,:) !for postprocedure only. storing the shape funciont values at nodes.
								!expolating_Lshape(nshape,nnum)
								!for cpe15,cpe15_spg, the matrix keeps the extrapolation matrix
		real(kind=DPN),allocatable::termval(:,:) !for exptrapolation, cpe15
	end type
	type(ECP_tydef)::ecp(maxiet)  !ecp(1:maxet), the class property of element(i) is stored in ecp(element(i).et)
	!element list sharing the same vertex, the mid- and insided nodes are excluded
	
	type geostatic_tydef
		logical::isgeo=.false.
		integer::method=1 ! cal=1,ko=2, the ko method only for horizontal ground surface
											!horizontal soil layer
		integer::nsoil=1 !the number of soil layers
		real(kind=DPN),allocatable::KO(:),weight(:),height(:)
        integer,allocatable::eset(:) !�����Ӧ��ʱ������ĵ�Ԫ��,Ĭ��ȫ����Ԫ�������
        integer::neset=0
	end type
	type(geostatic_tydef)::geostatic

	!output variables
	type outvar_tydef
		character(128)::name=''
		integer::value=0,ivo=0	!>0, variable output is required.,ivo=����������NodalQ(:,ivo)
		logical::iscentre=.false. !location, nodes or centroid of element.
		integer::system=0	!reference systerm,the default system is the globel sysytem
								!if system=-999, then the local system is a cylindrical system whose origin is along 
								!the central line of the cylinder.
								!if system=-9999,then the local system is a spherical syystem whose origin is located
								!at the center of the sphere
		integer::nval=1 !��������������ٸ�����
	end type
	type(outvar_tydef)::outvar(200)	
	integer::vo(200)=0,nvo=0
	
	type coordinate_tydef !global to local system transformation matrix 
		real(kind=DPN)::c(3,3)=0 !c(i,j)=cos(ei,Ej),ei and Ej are base vectors for new and old coordinates.
	end type	
	type(coordinate_tydef),allocatable::coordinate(:) !coordinate(0) is the whole coordinate system.
	integer::ncoord=0
    
	
	!ע�⣬master���ϵ�mdof���ɶ���Ŀǰ����ʩ��λ�Ʊ߽���������ʱ����ʱ�任һ��slave��master
	type slave_master_node_pairtydef		
		integer::slave=0,sdof=0		
		integer::master=0,mdof=0
		integer::nmbl=0 !������master�ڵ��ϵĺ��ظ���,Ŀǰ�������10����
		integer::mbl(10)=0
		integer::pe=0 !��Ӧ�ķ���Ԫ
		real(kind=DPN)::interforce=0.0d0,load=0.0d0
		real(kind=DPN)::aff=0.0d0 !allowablefrictionforce
	end type
	type(slave_master_node_pairtydef),allocatable::smnp(:)
	integer::nsmnp=0
	
	type out_data_typdef
        integer::nnode=0,isSUMQ=0,isStat=0 !isSumQ/=0,ouput the sum of nodal discharges in node(). 
        !isStat>0 ouput statistics(sum,max,min,mean,median,mad,std,kurtosis,skewness) of each variable in outputlist.
        integer,allocatable::node(:)
        real(kind=DPN),allocatable::stat(:,:) !stat([sum,max,min,mean,median,mad,std,kurtosis,skewness],[nval])
    end type
    type(out_data_typdef),allocatable::DataPoint(:)
    integer::NDataPoint=0
	
    TYPE DOFADJL_TYDEF
        INTEGER::NDOF=0,DOF_SIZE=20
        INTEGER,ALLOCATABLE::DOF(:) !SORTED BY COLUMNS 
    CONTAINS
        PROCEDURE::ADD_ITEM=>DOFADJL_ADD_ITEM
    ENDTYPE
    TYPE(DOFADJL_TYDEF),ALLOCATABLE::DOFADJL(:)
	
    !type rcd_set_tydef
    !    integer::rcd
    !    integer::enums,enume
    !end type
    !type(rcd_set_tydef),allocatable::rcdset(:)
    !integer::nrcdset=0
    !integer,allocatable::rcd(:,:)
    !integer::nrcd=0

	
	character(1024)::title,resultfile,resultfile1,resultfile2,resultfile3,resultfile21,resultfile22,EXCAMSGFILE,EXCAB_BEAMRES_FILE,&
					EXCAB_STRURES_FILE,EXCAB_EXTREMEBEAMRES_FILE,SLOPE_FILE,helpfile,well_file,flownet_file,volfile,poreflowfile
	INTEGER::DATAPOINT_UNIT=29,NSFR=0
	integer::datapacking=1	!=1,point format:{x1,y1,z1},{x2,y2,z2},..., . (Default Format)
						 != 2, block format, {x1,x2,...},{y1,y2,...},{z1,z2,...}
	integer::ndimension=2	!nodal dimension,
						!=1, input x
						!=2, input x and y
						!=3 input x,y and z
							
	integer::ndof=0 !total dof number
	INTEGER*8,allocatable::bw(:) 	!bw(i): firstly, it is the column number of the most left entry in the i row. and later, it is the bandwidth of the total matrix
	integer,allocatable::adof(:)	!adof=1,active or deactive												!!for default solver, finally it stores the diagonal address in the total stiffness matrix
	
	integer::bwmax=0 !the maximum value of the band width. 
	real(kind=DPN),allocatable::load(:),Tload(:),km(:),tdisp(:),bfload(:),stepload(:)	!km(:):Lower  trianglar part of total stiffness matrix(including diagonal elements)
	real(kind=DPN),allocatable::tdispInLS(:) ! Total displacement in a local system.
	real(kind=DPN),allocatable::NI_NodalForce(:) !�ɵ�Ԫ���ֶ��õĽڵ�����
	real(kind=DPN),allocatable::Tstepdis(:,:) !TstepDisp(ndof,nstep):��������ʱ�ĸ����ɶȵ����� 
	!additioinal variables defined when the linear system is unsymmetric.
!	real(kind=DPN),allocatable::ukm(:)	!ukm(:):Upper trianglar part of total stiffness matrix(including diagonal elements,but set equel to 0)
!	integer,allocatable::ubw(:)
	
!	real(kind=DPN),allocatable::a(:) !nonzero elments in km(:) and ukm(:)	

	INTEGER::nnz=0
	integer,allocatable::jcol(:),Lmre(:),adrn(:) !Lmre(i): the column number of the most rigth entry in the i row.
	integer,allocatable::ROWINDEX(:) !rowindex:�ܸ���ÿһ�е�һ������Ԫ�����ܸ������е�λ��;
											!bw: һ��ʼΪ�ܸ���ÿһ�е�һ������Ԫ�����ܸ��е��кţ����Ϊ�ܸ���ÿһ�����һ������Ԫ�����ܸ������е�λ��;
											!Lmre:�ܸ���ÿһ�����һ������Ԫ�����ܸ��е��кš�
											!ֵ��ע����ǣ���default solver�У��ܸպ�����Ԫ��������ÿ�е�һ������Ԫ�����һ������Ԫ֮�䣩������mkl������洢��ʽ�У��ܸղ�����Ԫ��
											!adrn(i):��default solver�е��ܸ������з���Ԫ i ��mkl��ʽ�ܸ������λ�á�
	
	
	real(kind=DPN),allocatable::DIAGLKM(:) !STORED DIAGONAL ELEMENT OF THE KM BEFORE FACTORIZATION FOR FACTOR UPDATING.
	integer,allocatable::diaglkmloc(:) !�ܸ��жԽ���Ԫ�����ܸ������е�λ�á�

	!one dimensional linear field function
	real(kind=DPN)::LF1D(0:maxilf,2)  !LF1D(:,1)=k,LF1D(:,2)=c. then y=k*x+c
	real(kind=dpn),allocatable:: slope_guide_line(:,:)
	integer::nsgline=0
!	integer::NLF1D=0 

	!additional variables for upper bound analysis
	integer::ncons=0 !number of constraints imposed by elements
	integer::nuvar=0 !uncontraint variables in linear programming
	integer::nzone_tec=0 !the number of zone in output file in tecplot format.
	integer::varsharezone_tec=0 !the zone number whose variables are shared with others in output file in tecplot format.
	
	real(kind=DPN)::NormL=0.0
	
	real(kind=DPN)::minNPH=1e20,MAXSFR=-1.D20,MAXSFR_LAST=-1.D20
	integer::isref_spg=0 !�Ƿ����·ֽ����
	real(kind=DPN)::eps1=1e20,eps2=1e20
	real(kind=DPN)::Origin_Sys_cylinder(3)=0.0  !in the order of x,y and z. the default value is 0.
	real(kind=DPN)::Qinput=0.0d0,Qstored=0.0d0,Qstorted_ini=0.0d0 !�����غ� 
	
	real(kind=DPN)::BARFAMILY_DIAGRAM_SCALE(12)=0.0D0 !�ȴ�������ֵ����滭ͼ�Ŵ�ϵ����DISX,DISY,DISZ,RX,RY,RZ,FX,FY,FZ,MX,MY,MZ
	real(kind=DPN)::barfamily_minxyz(3)=1.0d20,barfamily_maxxyz(3)=-1.0d20,MYFVAL=-1.0D20,SICR=0 ! SICR=STRESS INTEGRATION CONVERGE RATIO
    INTEGER::ISEXCA2D=0,ISHBEAM=0,ISSLOPE=0,NYITER(2)=0
	
	INTEGER::NDOFHEAD=0,NDOFMEC=0 !!ÿ�����������ɶ�������Ϣ���ɶ�����
	INTEGER,ALLOCATABLE::DOFHEAD(:),DOFMEC(:),CalStep(:) !head dofs in the model.
	real(kind=DPN),allocatable::NodalQ(:,:,:),RTime(:) !NODALQ(INODE,IVO,NNODALQ) !NNODALQ=SUM(TEIMSTEP.nsubts)
	INTEGER::NNODALQ=0,isporeflow=0
	integer,allocatable::bfgm_step(:)
    integer::mpi_rank = 0, mpi_size = 1, mpi_ierr
    LOGICAL::ISINISEDGE=.FALSE.,ISOUT_WELL_FILE=.FALSE.
    INTEGER,ALLOCATABLE::ISPF(:,:) 
    
    INTERFACE
          subroutine INVARIANT(stress,inv)
	        implicit none
	        real(8),intent(in)::stress(6)
	        real(8),intent(out)::inv(3)
	
        end subroutine 
		
        SUBROUTINE KSITA_MC_C2(KSITA,LODE,SITA,PHI)
            !LODE IN RAD,SITA IN RAD, TRANSITIONAL ANGLE.
	        IMPLICIT NONE
	        REAL(8),INTENT(IN)::LODE,SITA,PHI
	        REAL(8),INTENT(OUT)::KSITA(6)
        END SUBROUTINE
		subroutine yieldfun(vyf,mat,inv,ayf,ev,ISTEP)
			!use solverds	
			implicit none
			INTEGER,INTENT(IN)::MAT,AYF,ISTEP
			REAL(8),INTENT(IN)::INV(3),EV
			REAL(8),INTENT(OUT)::VYF
		END SUBROUTINE
		
		SUBROUTINE MC_KSITA(LODE,SITA,PHI,A,B,KSITA,D_KSITA)
			IMPLICIT NONE
			REAL(8),INTENT(IN)::LODE,SITA,A(2),B(2),PHI
			REAL(8),INTENT(OUT)::KSITA,D_KSITA
		END SUBROUTINE
		
		 subroutine deriv_yf(dyf,dywi,m)
			!use solverds
			implicit none
			REAL(8),INTENT(IN)::dywi(3),M(6,3)
			REAL(8),INTENT(OUT)::DYF(6)
		END SUBROUTINE
		
		subroutine deriv_sinv(m,stress)
			implicit none
			REAL(8),INTENT(IN)::STRESS(6)
			real(8),INTENT(OUT)::m(6,3)
		END SUBROUTINE
		
		subroutine deriv_yf_with_inv(dywi,inv,mat,ayf,ev,ISTEP)
			
			implicit none
			INTEGER,INTENT(IN)::MAT,AYF,ISTEP
			real(8),INTENT(IN)::inv(3),EV
			real(8),INTENT(OUT)::dywi(3)
			
		END SUBROUTINE
		
		subroutine deriv_qf_with_inv(dywi,inv,mat,ayf,ev,ISTEP)
		
			implicit none
			INTEGER,INTENT(IN)::MAT,AYF,ISTEP
			real(8),INTENT(IN)::inv(3),EV
			real(8),INTENT(OUT)::dywi(3)	
			
		END SUBROUTINE
        
        FUNCTION PARA_MC_CLAUSEN(MATID,ISTEP) RESULT(PARA)
            IMPLICIT NONE
            INTEGER,INTENT(IN)::MATID,ISTEP
            REAL(8)::PARA(3)
        END FUNCTION
        

		
	END INTERFACE    
   
    
    contains
    real(8) function get_gravity(this)
    !calculte the g value based on the model unit
        class(solver_tydef)::this
        real(8)::t1,t2
 
       
        select case(this.len_unit)
        case(0) !m
            t1=1.0
        case(1) !dm
            t1=10
        case(2) !cm
            t1=100
        case(3) !mm
            t1=1000
        case(4) !km
            t1=1.0/1000
        case default
            error stop 'no such length unit.'            
        end select
        
        
        select case(this.time_unit)
        case(0) !day
            t2=(24*3600)
        case(1) !sec
            t2=1.0
        case(2) !min
            t2=60
        case(3) !hour
            t2=3600
        case default
            error stop 'no such time unit.'            
        end select        
        t2=t2**2
        
        get_gravity=9.8*t1*t2
        
    end function
    real(8) function unit_scaling_factor(this,unit)
    !calculate the scale factor scaling the current model unit to the input unit.
        class(solver_tydef)::this
        character(len=*)::unit
        real(8)::t1,t2
        !to m
        t1=1.0
        select case(this.len_unit)
        case(0) !m
            t1=1.0
        case(1) !dm
            t1=0.1
        case(2) !cm
            t1=0.01
        case(3) !mm
            t1=0.001
        case(4) !km
            t1=1000
        case default
            error stop 'no such length unit.'
            
        end select
        
        !to day
        select case(this.time_unit)
        case(0) !day
            t2=1.0
        case(1) !sec
            t2=1.0/(24*3600)
        case(2) !min
            t2=1./(24*60)
        case(3) !hour
            t2=1./(24)
        case default
            error stop 'no such time unit.'            
        end select        
        
        select case(trim(adjustl(unit)))
        case('m')
            unit_scaling_factor=t1
        case('mm')
            unit_scaling_factor=t1*1000
        case('cm')
            unit_scaling_factor=t1*100
        case('dm')
            unit_scaling_factor=t1*10
        case('km')
            unit_scaling_factor=t1/1000
        case('day')
            unit_scaling_factor=t2
        case('hour','h')
            unit_scaling_factor=t2*24
        case('sec','s')
            unit_scaling_factor=t2*24*3600
        case('min')
            unit_scaling_factor=t2*24*60
        case default
            print *, 'no such unit=',unit
            error stop
        end select
        
    endfunction

    SUBROUTINE DOFADJL_ADD_ITEM(SELF,ITEM) 
        
        class(DOFADJL_TYDEF)::self
        INTEGER,INTENT(IN)::ITEM
        INTEGER::I,N1
    

        IF(.NOT.ALLOCATED(SELF.DOF)) THEN
            ALLOCATE(SELF.DOF(SELF.DOF_SIZE))
            SELF.DOF=0
        ENDIF
        IF(.NOT.ANY(SELF.DOF(1:SELF.NDOF)-ITEM==0)) THEN
            
            IF(SELF.NDOF+1>SELF.DOF_SIZE) THEN
                CALL I_ENLARGE_AR(SELF.DOF,10)
                SELF.DOF_SIZE=SELF.DOF_SIZE+10
            ENDIF
            IF(SELF.NDOF==0) THEN
                N1=SELF.NDOF+1
            ELSE
                IF(ITEM>SELF.DOF(SELF.NDOF)) THEN
                    N1=SELF.NDOF+1
                ELSEIF(ITEM<SELF.DOF(1)) THEN
                    N1=1
                    SELF.DOF(SELF.NDOF+1:2:-1)=SELF.DOF(SELF.NDOF:1:-1)
                ELSE
                    N1=MINVAL(PACK([1:SELF.NDOF],SELF.DOF(:SELF.NDOF)-ITEM>0))
                    SELF.DOF(SELF.NDOF+1:N1+1:-1)=SELF.DOF(SELF.NDOF:N1:-1)
                ENDIF
            ENDIF
            SELF.DOF(N1)=ITEM
            SELF.NDOF=SELF.NDOF+1
        ENDIF
    CONTAINS
        SUBROUTINE I_ENLARGE_AR(AVAL,DSTEP)
            IMPLICIT NONE
            INTEGER,ALLOCATABLE,INTENT(INOUT)::AVAL(:)
            INTEGER,INTENT(IN)::DSTEP
            INTEGER,ALLOCATABLE::VAL1(:)
            INTEGER::LB1=0,UB1=0
    
            LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
            ALLOCATE(VAL1,SOURCE=AVAL)
            DEALLOCATE(AVAL)
            ALLOCATE(AVAL(LB1:UB1+DSTEP))
            AVAL(LB1:UB1)=VAL1
            !AVAL(UB1+1:UB1+10)=0
            DEALLOCATE(VAL1)
        END SUBROUTINE

    ENDSUBROUTINE     

    
    FUNCTION GET_MAT_PROPERTY_ARRAY(MAT,IPARA,ISTEP) RESULT(VAL)
		CLASS(mat_tydef),INTENT(in):: MAT
		INTEGER,INTENT(IN)::IPARA(:),ISTEP
        REAL(DPN)::VAL(SIZE(IPARA))
        INTEGER::I
		VAL=MAT.property(ipara)*sf(MAT.sf(ipara)).factor(istep) 
        IF(MAT.TYPE==MC) THEN
            DO I=1,SIZE(IPARA)
                IF(IPARA(I)==4.OR.IPARA(I)==5) THEN
                    VAL(I)=TAN(MAT.property(IPARA(I))/180.*PI())*sf(MAT.sf(IPARA(I))).factor(istep) 
                    VAL(I)=ATAN(VAL(I))/PI()*180
                ENDIF
            ENDDO
        ENDIF
	END FUNCTION GET_MAT_PROPERTY_ARRAY

    FUNCTION GET_MAT_PROPERTY_SCALAR(MAT,IPARA,ISTEP) RESULT(VAL)
		CLASS(mat_tydef),INTENT(in):: MAT
		INTEGER,INTENT(IN)::IPARA,ISTEP
        REAL(DPN)::VAL
		VAL=MAT.property(ipara)*sf(MAT.sf(ipara)).factor(istep) 
        IF(MAT.TYPE==MC.AND.(IPARA==4.OR.IPARA==5)) THEN
            !SCALE TAN(ANGLE) BUT NOT ANGLE
           VAL=TAN(MAT.property(IPARA)/180.*PI())*sf(MAT.sf(IPARA)).factor(istep) 
           VAL=ATAN(VAL)/PI()*180
        ENDIF        
	END FUNCTION GET_MAT_PROPERTY_SCALAR
	
    INTEGER FUNCTION ESET_GETFREEID(ISET) 
		!CLASS(ESET_tydef),INTENT(in):: ESET
		INTEGER,OPTIONAL,INTENT(IN)::ISET
        INTEGER::I
        
        ESET_GETFREEID=-1
		IF(PRESENT(ISET)) THEN
			IF(ESET(ISET).ISINI==0) THEN
				ESET_GETFREEID=ISET
				ESET(ISET).ISINI=ISET
			ELSE
				PRINT *,'THE ESET(ID) HAS BEEN USED. PLEASE GIVE ANOTHER SETID.ID=',ISET
				STOP
			ENDIF
		ELSE
			DO I=1,MAXSET
				IF(ESET(I).ISINI==0) THEN
					ESET_GETFREEID=I
					ESET(I).ISINI=I
					EXIT
				ENDIF
			ENDDO
		ENDIF
        IF(ESET_GETFREEID==-1) THEN
            STOP 'exceed maxset limit. please increase it.'
        ENDIF
   
	END FUNCTION ESET_GETFREEID    
    
	 function pi()
		real(kind=DPN)::pi
		pi=1.d0
		pi=datan(pi)*4.0d0	
	end function
	
	real(DPN) function matproperty(imat,ipara,istep)
	
		implicit none
		integer,intent(in)::imat,ipara,istep
		matproperty=material(imat).property(ipara)*sf(material(imat).sf(ipara)).factor(istep)
	endfunction
	!
	real(kind=DPN) function materialfun(Fid,x)
		integer::fid
		real(kind=DPN)::x	
			
		materialfun=LF1D(FID,1)*X+LF1D(FID,2)
	end function


FUNCTION csproduct(a1,a2) RESULT(maxoa)
!
! This function returns the dyadic product of two arrays
! 
	 IMPLICIT NONE
     INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
	 real(kind=iwp),INTENT(IN)::a1(:),a2(:)
	 real(kind=iwp)::maxoa(size(a1),size(a2))
	 INTEGER::nrow,i,j,ncol
	 
	 nrow=UBOUND(a1,1)
	 ncol=UBOUND(a2,1)
	 DO i=1,nrow
		DO j=1,ncol
			maxoa(i,j)=a1(i)*a2(j)
		ENDDO
	ENDDO
	RETURN
END FUNCTION csproduct    


! computate the elastic compliance matrix [C] so that strain=matmul([C],stress)
!ND>=4
 function DINV(E,V,ND) result(C)
	implicit none
	INTEGER,INTENT(IN)::ND
	real(8),intent(in)::E,V
	real(8)::C(ND,ND)
	
	C=0.D0
	C(1:3,1:3)=-V
	C(1,1)=1.D0
	C(2,2)=1.D0
	C(3,3)=1.D0
	C(4,4)=2.D0*(1.D0+V)
	IF(ND>4) THEN
		C(5,5)=C(4,4)
		C(6,6)=C(4,4)
	ENDIF
    C=C/E
	
end function

FUNCTION KM_WELLBORE(A1,A2,A3,A4,A5) RESULT(KM)
    IMPLICIT NONE
    REAL(8),INTENT(IN)::A1,A2,A3
    REAL(8),INTENT(IN),OPTIONAL::A4,A5
    REAL(8)::A15,A25
    INTEGER::ND1=4
    REAL(8),ALLOCATABLE::KM(:,:)
    
    IF(PRESENT(A4)) THEN
        ND1=5
        A15=A4;A25=A5
    ELSE
        ND1=4
        A15=0.D0
        A25=0.D0
    ENDIF
    ALLOCATE(KM(ND1,ND1))
    
    KM=0.D0
    KM(1,1)=A1+A3+A15
    KM(2,2)=A1+A2+A25
    KM(3,3)=A2
    KM(4,4)=A3
    KM(1,2)=-A1
    KM(2,1)=-A1
    KM(2,3)=-A2
    KM(3,2)=-A2
    KM(1,4)=-A3
    KM(4,1)=-A3
    
    IF(PRESENT(A4)) THEN
        KM(5,5)=A15+A25
        KM(1,5)=-A15;KM(5,1)=-A15
        KM(2,5)=-A25;KM(5,2)=-A25
    ENDIF

END FUNCTION


REAL(8) FUNCTION fD_PF(RE,KR,MODEL,REW,POROSITY) !darcy-friction for pipe flow
    IMPLICIT NONE
!function: calculate darcy-friction for pipe flow,if ReW>0,then it is a porous pipe flow
!Re, Reynolds number (unitless);
!Kr, relative roughness (unitless)
!ReW, wall Reynolds number (unitless);
!REF:[1] Fang X, Xu Y, Zhou Z. New Correlations of Single-Phase Friction Factor for Turbulent Pipe Flow and Evaluation of Existing Single-Phase Friction Factor Correlations[J]. Nuclear Engineering and Design, 2011, 241(3): 897-902. 
![2] Ouyang L-B, Arbabi S, Aziz K. General Wellbore Flow Model for Horizontal, Vertical, and Slanted Well Completions[J]. SPE Journal 1998, 3(2): 124~133.
    REAL(8),INTENT(IN)::RE,KR
    REAL(8),INTENT(IN),OPTIONAL::REW,POROSITY
    INTEGER,INTENT(IN),OPTIONAL::MODEL
    REAL(8)::LAMDA1,REW1,FC1,A,B,C,PO1,FO1,KR1
    INTEGER::MODEL1
    
    !IF(ABS(RE)<1.D-7) THEN
    !    LAMDA1=1E7
    !    fD_PF=LAMDA1
    !    RETURN
    !ENDIF
    
    MODEL1=0
    IF(PRESENT(MODEL)) MODEL1=MODEL
    
    
    IF(MODEL1==1) THEN
        PO1=0.D0
        IF(PRESENT(POROSITY)) PO1=POROSITY
        KR1=KR+0.282*PO1**2.4
        FO1=0.0106*PO1**0.413
    ENDIF
    
    IF(RE<2260) THEN
        LAMDA1=64./MAX(RE,1.0D-10)    
    ELSEIF(RE<3400) THEN
        A=-2.*log10(12/Re+KR1/3.7)
        B=-2*log10(2.51*A/Re+KR1/3.7)
        C=-2*log10(2.51/Re+KR1/3.7)
        LAMDA1=A-(B-A)**2/(C-2*B-A)    
        LAMDA1=(1./LAMDA1)**2          
    ELSEIF(KR1>0.D0) THEN        
        LAMDA1=1.613*(LOG(0.234*(KR1)**1.1007-60.525/RE**1.1105+56.291/RE**1.0712))**(-2.0D0)
    ELSE
        LAMDA1=0.25*(LOG10(150.39/RE**0.98865-152.66/RE))**(-2.D0)
    ENDIF
    
    
    
    SELECT CASE(MODEL1)
    
    CASE(0) !DARCY,DEFAULT
        fD_PF=LAMDA1    
    CASE(1) !SIWON
        fD_PF=LAMDA1+FO1
    CASE(2) !OUYANG
        REW1=0.D0
        IF(PRESENT(REW))    REW1=REW
        FC1=1.D0
        IF(ABS(REW1)>1.D-7) THEN
            !INFLOW
            IF(REW1>0.D0) THEN
                IF(RE<3000) THEN
                !LAMINAR
                    FC1=1.0+0.04304*REW1**0.6142
                ELSE
                    !TURBULENT
                    FC1=1.0-0.0153*REW1**0.3978
                ENDIF
        
            ELSE
            !OUTFLOW
                IF(RE<3000) THEN
                !LAMINAR
                    !FC1=1.0-0.0625*(-REW1)**1.3056/(REW1+4.626)**-0.2724
                
                        !!!!!!!!!!!!
                    !IF(ISNAN(FC1)) FC1=1.0-0.04304*(-REW1)**0.6142
                    FC1=1.D0
                ELSE
                !TURBULENT
                    FC1=1.0-17.5*REW1/RE**0.75
                ENDIF        
        
            ENDIF
        ENDIF
    
        IF(FC1<0.D0) FC1=1.D0
    
        fD_PF=LAMDA1*FC1    
    CASE(3) !USER INPUT
        
    ENDSELECT
    
    
    
    
    
    IF(ISNAN(fD_PF)) THEN
        PRINT *, 'fD_PF IS NAN'
    ENDIF

ENDFUNCTION

REAL(8) FUNCTION Vk(T)
    IMPLICIT NONE
    !Kinematic viscosity, unit=m2/s,fitting range 0<=T<=80
    REAL(8),INTENT(IN)::T
!VK= -7.2087579E-10*T1**5 + 1.9555068E-7*T1**4 - 2.2103680E-5*T1**3 + 1.4147594E-3*T1**2 - 6.0104832E-2*T1 + 1.7867110
VK=0.040598559873*10**(183.079140630271/(273.15+T-161.737209235748)) 
VK=VK*1.D-6 !m2/s

ENDFUNCTION

REAL(8) FUNCTION Re_W(V,D,T)
!calculate Reynold number of water for pipe flow
!V,VELOCITY,m/s
!D,DIAMETER,m
!T,temperature,celsius ���϶�, fitting range 0<=T<=80
IMPLICIT NONE
REAL(8),INTENT(IN)::V,D !unit,v,m/s, D,m
REAL(8),INTENT(IN),OPTIONAL::T !unit= ���϶�
REAL(8)::T1,cf1,cf2

IF(.NOT.PRESENT(T)) THEN
    T1=25
ELSE
    T1=T
    IF(T1<1.0D-6) T1=25
ENDIF
!to sec
cf1=solver_control.unit_factor('s')
!convert to m
cf2=solver_control.unit_factor('m')



Re_W=ABS(V*cf2/cf1*D*cf2/Vk(T1))

    
ENDFUNCTION

SUBROUTINE PNW_FLOW_RIJ(THIS)
    !�����϶���絥Ԫ����������,��Rij=1/Kij
    !Kij=gij/rw,rwΪˮ���ض�
    !gij=Pi*D**4/(128*u*L),uΪˮ�ľ���ճ��ϵ��
    !Kij=Pi*D**4/(128*vk/g*L)

    IMPLICIT NONE
    CLASS(element_tydef)::THIS
	!real(8),optional::Lc,Kc !�ٶ���ĳ��ȼ���͸ϵ�� 
    real(8)::t1,t2
    !Rt,��뾶,Lt����,Ri,Rj�ֱ�Ϊ���˿׵İ뾶
    REAL(8)::lc1,kc1,lt,alpha1,d1
    integer::i,j
	
	
	THIS.property(1)=0.d0
	do i=1,2
		t1=this.pfp(3)
		if(i==2) t1=1-t1
		lc1=this.property(4)*t1-this.pfp(5+i)
		if(lc1>0.d0) then
			alpha1=(this.pfp(i)-this.property(2))/this.property(2)
			d1=this.property(2)+this.pfp(5+i)*alpha1			
			this.property(1)=this.property(1)+hagen_poiseuille_friction(this.mat,lc1,d1,this.pfp(i),node(this.node(i)).cc)
		endif
	enddo

	this.property(1)=this.property(1)+this.Property(5)+this.Property(6)
	if(abs(this.property(1))<1.d-10) this.property(1)=1.0d-10
    THIS.FD=THIS.PROPERTY(1)

ENDSUBROUTINE

REAL(8) FUNCTION KC_K(e,dp,temp)
	!����KC��ʽ��������͸ϵ��,unit=L/T
	!input:e=��϶��,dp=����ֱ��
	implicit none 
	REAL(8),intent(in)::e,dp,temp
	optional::temp
	real(8)::t1,t2,temp1
	temp1=20
	if(present(temp)) temp1=temp
	!t1=u/rw
	t1=vk(temp1)/9.8 !unit:(m.sec)
	t2=solver_control.unit_factor('m')*solver_control.unit_factor('s')
	t1=t1/t2 !to model unit
	!t1=rw/u
	t1=1./t1
	KC_K=t1*(dp**2)*(e**3)/(180.*(1-e)**2)

ENDFUNCTION

    real(8) function hagen_poiseuille_friction(imat,l,d1,d2,cc)
        !�����϶���絥Ԫ����������,��Rij=1/Kij
        !Kij=gij/rw,rwΪˮ���ض�
        !gij=Pi*D**4/(128*u*L),uΪˮ�ľ���ճ��ϵ��
        !Kij=Pi*D**4/(128*vk/g*L)
        implicit none
        integer,intent(in)::imat
        real(8),intent(in)::l
        real(8),optional::d1,d2,cc
        real(8)::t1,t2,d12,d22,ccmax1=0.66,cc1

        ![1] Hirabayashi S, Sato T, Mitsuhori K, Yamamoto Y. Microscopic numerical simulations of suspension with particle accumulation in porous media[J]. Powder Technology, 2012, 225: 143-148.
        !DOI: https://doi.org/10.1016/j.powtec.2012.04.001
        !consider the effect of suspension particles on the viscosity.
        cc1=0.d0
        if(present(cc)) cc1=cc
        if(cc1>0.d0) then
            t2=(1-cc1/ccmax1)**(-2.5*ccmax1) 
        else
            t2=1.0d0
        endif
        
        d12=material(imat).property(1)*2        
        if(present(d1)) d12=d1
        d22=d12
        if(present(d2)) d22=d2

        t1=vk(material(imat).property(2))/9.8*t2 !unit=m.sec
        t2=solver_control.unit_factor('m')*solver_control.unit_factor('s')
        t1=t1/t2 !to model unit 
        t2=(d12+d22)/2.0   
        t1=PI()*t2**4/(128.*t1*L)		
        hagen_poiseuille_friction=1./t1
	          

    endfunction 

real(8) function vol_cone(h,r1,r2)
	implicit none 
	real(8),intent(in)::h,r1,r2 
	vol_cone=1/3.0*PI()*h*(r1**2+r2**2+r1*r2)
endfunction

   !���ַ������൱�������ַ�(����������)ת��Ϊ��Ӧ������
   !�� '123'תΪ123,'14-10'תΪ14,13,12,11,10
   !string��ת���������������ar(n1)���أ�����,n1Ϊ�ַ��������ֵĸ���:(ע��1-3ת����Ϊ3�����֣�1,2,3)
   !nmaxΪ����ar�Ĵ�С,stringĬ���ַ�����Ϊ1024��
   !num_readΪҪ�������ݵĸ�����
   !unitΪ�ļ���
   !ÿ��ֻ����һ����Ч�У�����'/'��ͷ���У�
   !ÿ�к�����'/'��ʼ�ĺ�����ַ�����Ч�ġ�
   subroutine  strtoint(unit,ar,nmax,n1,num_read,set,maxset,nset,ef1,isall)
		implicit none
		INTEGER,INTENT(IN)::unit,nmax,num_read,maxset
		INTEGER,INTENT(INOUT)::N1,NSET
		INTEGER,OPTIONAL::EF1
		REAL(8),INTENT(INOUT)::ar(nmax)
		character(*)::set(maxset)
		logical::tof1,tof2
		integer::i,j,k,strl,ns,ne,n2,n3,n4,step,& 
				ef,n5,nsubs
		real(8)::t1	  
		character(20000)::string
		character(512)::substring(1000)
		character(16)::legalC,SC
		logical,optional::isall

		LegalC='0123456789.-+eE*'
		sc=',;() '//char(9)
		n1=0
		nset=0
		ar=0
		!set(1:maxset)=''
	  do while(.true.)
		 read(unit,'(a)',iostat=ef) string
         
		 if(ef<0) then
            if(present(ef1))  then
                ef1=ef
            else
			    print *, 'file ended unexpected. sub strtoint()'
			    stop
            endif
		 end if

		 string=adjustL(string)
		 strL=len_trim(string)
		 
		do i=1,strL !remove 'Tab'
			if(string(i:i)/=char(9)) exit
		end do
		string=string(i:strL)
		string=adjustl(string)
		strL=len_trim(string)
		if(strL==0) cycle

		 if(string(1:2)/='//'.and.string(1:1)/='#') then
			
			!ÿ�к�����'/'��ʼ�ĺ�����ַ�����Ч�ġ�
			if(index(string,'//')/=0) then
				strL=index(string,'//')-1
				string=string(1:strL)
				strL=len_trim(string)
			end if

			nsubs=0
			n5=1
			do i=2,strL+1
				if(index(sc,string(i:i))/=0.and.index(sc,string(i-1:i-1))==0) then
					nsubs=nsubs+1					
					substring(nsubs)=string(n5:i-1)					
				end if
				if(index(sc,string(i:i))/=0) n5=i+1
			end do
			
			do i=1, nsubs
				substring(i)=adjustl(substring(i))				
				n2=len_trim(substring(i))
				!the first character should not be a number if the substring is a set.
				if(index('0123456789-+.', substring(i)(1:1))==0) then
					!set
					nset=nset+1
					set(nset)=substring(i)
					cycle
				end if
				n3=index(substring(i),'-')
				n4=index(substring(i),'*')
				tof1=.false.
				if(n3>1) then
				    tof1=(substring(i)(n3-1:n3-1)/='e'.and.substring(i)(n3-1:n3-1)/='E')
				end if
				if(tof1) then !����������'1-5'��������ʽ�Ķ�������
					read(substring(i)(1:n3-1),'(i8)') ns
					read(substring(i)(n3+1:n2),'(i8)') ne
					if(ns>ne) then
						step=-1
					else
						step=1
					end if
					do k=ns,ne,step
						n1=n1+1
						ar(n1)=k
					end do				     	
				else
				     tof2=.false.
				     if(n4>1) then
				             tof2=(substring(i)(n4-1:n4-1)/='e'.and.substring(i)(n4-1:n4-1)/='E')
				     end if
					if(tof2) then !����������'1*5'(��ʾ5��1)��������ʽ�Ķ�������
						read(substring(i)(1:n4-1),*) t1
						read(substring(i)(n4+1:n2),'(i8)') ne
						ar((n1+1):(n1+ne))=t1
						n1=n1+ne
					else
						n1=n1+1
						read(substring(i),*) ar(n1)
					end if	
				end if			
			end do
		 else
			cycle
		 end if
		
		 if(n1<=num_read) then
			if(present(isall)) then
				if(isall.and.n1==num_read) exit
			 else
				exit
			 end if
		 else
		    if(n1>num_read)  print *, 'error!nt2>num_read. i=',n1
		 end if
	
	  end do	

   end subroutine 

    INTEGER FUNCTION INCOUNT(N)
        IMPLICIT NONE
        INTEGER,INTENT(IN)::N
        REAL(8)::T1
        INTEGER::I

        T1=ABS(N)
        INCOUNT=INT(LOG10(T1))+1
	    IF(N<0) INCOUNT=INCOUNT+1
    

    END FUNCTION   

    !translate all the characters in term into lowcase character string
    subroutine lowcase(term)
	    use dflib
	    implicit none
	    integer i,in,nA,nZ,nc,nd
	    character(1)::ch
	    character(*)::term
	
	    term=adjustl(trim(term))
	    nA=ichar('A')
	    nZ=ichar('Z')
	    nd=ichar('A')-ichar('a')
	    in=len_trim(term)
	    do i=1,in
		    ch=term(i:i)
		    nc=ichar(ch)
		    if(nc>=nA.and.nc<=nZ) then
			    term(i:i)=char(nc-nd)
		    end if
	    end do
    end subroutine

    integer function get_free_file_unit_number

		integer, parameter :: first_unit=  10
		integer, parameter :: last_unit =9999
		logical test

		do get_free_file_unit_number=first_unit,last_unit
			inquire(get_free_file_unit_number,opened=test)
			if (.not.test) return
		enddo
		get_free_file_unit_number=-1

	end function

	function get_gp_dis(istep,ielt,igp,Ddis) result(dis)                                                                                                                                                                                    
																																																											
		implicit none                                                                                                                                                                                                                         
		integer,intent(in)::istep,ielt,igp                                                                                                                                                                                                    
		real(8),intent(in)::Ddis(:)                                                                                                                                                                                                           
		real(8)::dis(ndimension)                                                                                                                                                                                                              
		integer::I,nnode1,ndof1                                                                                                                                                                                                               
		real(8)::xdis(50)                                                                                                                                                                                                                     
																																																											
		nnode1=element(ielt).nnum                                                                                                                                                                                                             
		ndof1=element(ielt).ndof                                                                                                                                                                                                              
		do i=1,ndimension                                                                                                                                                                                                                     
			xdis(1:nnode1)=Tstepdis(node(element(ielt).node(1:nnode1)).dof(i),istep)+ddis([i:ndof1:ndimension])                                                                                                                                 
			dis(i)=dot_product(xdis(1:nnode1),ecp(element(ielt).et).lshape(:,igp))                                                                                                                                                              
		enddo                                                                                                                                                                                                                                 
																																																											
	end function 

	function get_gp_direction(istep,ielt,igp) result(dis)                                                                                                                                                                                    
																																																											
		implicit none                                                                                                                                                                                                                         
		integer,intent(in)::istep,ielt,igp
		integer::i,n1                                                                                                                                                                                                    
 		real(8)::dis(ndimension)                                                                                                                                                                                                              
		real(8)::x1(ndimension),t,f,g,h,x2(ndimension),t2=0,t1  
        

        if(nsgline==0) then
			error stop 'Please input the parameters by keyword= slope_guide_line'
		endif

		x1=element(ielt).xygp(1:ndimension,igp)
		t2=1.d20
		if(nsgline<2.or.ndimension==2) then
			dis=slope_guide_line(1:ndimension,1)-x1
		else
			do i=1,nsgline-1
				f=slope_guide_line(1,i+1)-slope_guide_line(1,i)
				g=slope_guide_line(2,i+1)-slope_guide_line(2,i)
				h=slope_guide_line(3,i+1)-slope_guide_line(3,i)
				call line_par_point_near_3d ( f, g, h, slope_guide_line(1,i), slope_guide_line(2,i), slope_guide_line(3,i), &
				x1, x2,t )				
				if(t>1.0d0) then
					x2=slope_guide_line(:,i+1)
				elseif(t<0.0d0) then
					x2=slope_guide_line(:,i)
				endif
				t1=norm2(x2-x1)
				if(t1<t2) then
					t2=t1
					dis=(x2-x1)/max(t1,1.0d-8)
				endif 
			enddo 
		endif                                                                                                                                                                                                                     
																																																											
	end function
	subroutine line_par_point_near_3d ( f, g, h, x0, y0, z0, p, pn,t )

		!*****************************************************************************80
		!
		!! LINE_PAR_POINT_NEAR_3D: nearest point on parametric line to given point, 3D.
		!
		!  Discussion:
		!
		!    The parametric form of a line in 3D is:
		!
		!      X = X0 + F * T
		!      Y = Y0 + G * T
		!      Z = Z0 + H * T
		!
		!    We may normalize by choosing F*F + G*G + H*H = 1, and F nonnegative.
		!
		!  Licensing:
		!
		!    This code is distributed under the GNU LGPL license. 
		!
		!  Modified:
		!
		!    12 April 2013
		!
		!  Author:
		!
		!    John Burkardt
		!
		!  Reference:
		!
		!    Adrian Bowyer, John Woodwark,
		!    A Programmer's Geometry,
		!    Butterworths, 1983,
		!    ISBN: 0408012420.
		!
		!  Parameters:
		!
		!    Input, real ( kind = 8 ) F, G, H, X0, Y0, Z0, the parametric 
		!    line parameters.
		!
		!    Input, real ( kind = 8 ) P(3), the point whose distance from the line is
		!    to be measured.
		!
		!    Output, real ( kind = 8 ) PN(3), the point on the parametric line which
		!    is nearest to P.
		!
		implicit none

		integer ( kind = 4 ), parameter :: dim_num = 3

		real ( kind = 8 ) f
		real ( kind = 8 ) g
		real ( kind = 8 ) h
		real ( kind = 8 ) p(dim_num)
		real ( kind = 8 ) pn(dim_num)
		real ( kind = 8 ) t
		real ( kind = 8 ) x0
		real ( kind = 8 ) y0
		real ( kind = 8 ) z0

		t = ( f * ( p(1) - x0 ) + g * ( p(2) - y0 ) + h * ( p(3) - z0 ) ) &
			/ ( f * f + g * g + h * h )

		pn(1) = x0 + t * f
		pn(2) = y0 + t * g
		pn(3) = z0 + t * h

		return
	end subroutine  
end module

INCLUDE 'mkl_dss.f90' 
MODULE MKLDS

	USE MKL_DSS
	TYPE(MKL_DSS_HANDLE) :: handle ! Allocate storage for the solver handle.
	
END MODULE
    

MODULE SOLVERLIB

    USE solverds
    
    INTERFACE
        subroutine enlarge_bc(bc,nbc,enbc,ibc)
            USE solverds
	        integer,intent(in)::enbc
	        type(bc_tydef),intent(in out),allocatable::bc(:)
	        integer,intent(out)::ibc
	        integer,intent(in out)::nbc
	        !type(bc_tydef),ALLOCATABLE::bc1(:)
        ENDSUBROUTINE
        
        subroutine enlarge_element(EL,NEL,enel,iel)
        !����EL����,ͬʱupdate�ܵĵ�Ԫ��NEL=NEL+enel
        !Enel:����ĵ�Ԫ����
        !iel,:���ݲ��ֵ���λ
            USE solverds
	        integer,intent(in)::enel
            INTEGER,INTENT(IN OUT)::NEL
            type(element_tydef),INTENT(IN OUT),ALLOCATABLE::EL(:)
	        integer,intent(out)::iel        
        ENDSUBROUTINE

        subroutine enlarge_node(EL,NEL,enel,iel)
        !����EL����,ͬʱupdate�ܵĵ�Ԫ��NEL=NEL+enel
        !Enel:����ĵ�Ԫ����
        !iel,:���ݲ��ֵ���λ
            USE solverds
	        integer,intent(in)::enel
            INTEGER,INTENT(IN OUT)::NEL
            type(node_tydef),INTENT(IN OUT),ALLOCATABLE::EL(:)
	        integer,intent(out)::iel        
        ENDSUBROUTINE		
		
        !SUBROUTINE invert(matrix)
        !    !
        !    ! This subroutine inverts a small square matrix onto itself.
        !    !
        !     IMPLICIT NONE
        !     INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
        !     REAL(iwp),INTENT(IN OUT)::matrix(:,:)
        !
        !ENDSUBROUTINE  
    ENDINTERFACE
    
    CONTAINS
    
    

    ENDMODULE

    
subroutine enlarge_bc(bc,nbc,enbc,ibc)
    !����bc����,ͬʱupdate�ܵĵ�Ԫ��NBC=NBC+enBC
    !EnBC:����ĵ�Ԫ����
    !iBC:���ݲ��ֵĳ�ʼ��λ
	    use solverds	
	    implicit none
	    integer,intent(in)::enbc
	    type(bc_tydef),intent(in out),allocatable::bc(:)
	    integer,intent(out)::ibc
	    integer,intent(in out)::nbc
	    type(bc_tydef),ALLOCATABLE::bc1(:)
	
	    allocate(bc1(nbc+enbc))
	    bc1(1:nbc)=bc
	    if(allocated(bc)) deallocate(bc)
	    ibc=nbc+1
	    nbc=nbc+enbc
	    !nbc=jbc
	    allocate(bc,source=bc1)
	    !element=element1
	    deallocate(bc1)
	
    endsubroutine
    
subroutine enlarge_element(EL,NEL,enel,iel)
!����EL����,ͬʱupdate�ܵĵ�Ԫ��NEL=NEL+enel
!Enel:����ĵ�Ԫ����
!iel,:���ݲ��ֵ���λ
	use solverds	
	implicit none
	integer,intent(in)::enel
    INTEGER,INTENT(IN OUT)::NEL
    type(element_tydef),INTENT(IN OUT),ALLOCATABLE::EL(:)
	integer,intent(out)::iel
	type(element_tydef),ALLOCATABLE::element1(:)
	
	IF(ENEL<1) RETURN
	
	allocate(element1(NEL+enel))
	element1(1:NEL)=EL
	if(allocated(EL)) deallocate(EL)
	iel=NEL+1
	NEL=NEL+enel	
	allocate(EL,SOURCE=ELEMENT1)
	deallocate(element1)
	
endsubroutine    

subroutine enlarge_node(EL,NEL,enel,iel)
!����EL����,ͬʱupdate�ܵĵ�Ԫ��NEL=NEL+enel
!Enel:����ĵ�Ԫ����
!iel,:���ݲ��ֵ���λ
	use solverds	
	implicit none
	integer,intent(in)::enel
    INTEGER,INTENT(IN OUT)::NEL
    type(node_tydef),INTENT(IN OUT),ALLOCATABLE::EL(:)
	integer,intent(out)::iel
	type(node_tydef),ALLOCATABLE::element1(:)
	
	IF(ENEL<1) RETURN
	
	allocate(element1(NEL+enel))
	element1(1:NEL)=EL
	if(allocated(EL)) deallocate(EL)
	iel=NEL+1
	NEL=NEL+enel	
	allocate(EL,SOURCE=ELEMENT1)
	deallocate(element1)
	
endsubroutine

subroutine enlarge_Gnode(ENEL)
!ENLARGE GNODE(3,N1) TO GNODE(3,N1+ENEL)

	use solverds	
	implicit none
	integer,intent(in)::ENEL    
	REAL(DPN),ALLOCATABLE::element1(:,:)
	INTEGER::N1
	
	IF(ENEL<1) RETURN
	N1=SIZE(GNODE,DIM=2)
	
	allocate(element1(3,N1+enel))
	element1(:,1:N1)=GNODE
	if(allocated(GNODE)) deallocate(GNODE)
	allocate(GNODE,SOURCE=ELEMENT1)
	deallocate(element1)
	
endsubroutine

    




