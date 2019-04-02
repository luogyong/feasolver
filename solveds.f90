  
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
		real(kind=DPN)::MISES=0.0D0,EEQ=0.0D0,PEEQ=0.0D0
		real(kind=DPN)::Q=0.0D0,Kr=0.0D0,Mw=0.0D0  !nodal flux,relative permeability,slope of volumetric water content.		
		integer::nelist=0,NELIST_SPG=0,ISACTIVE=0
		!integer,allocatable::elist(:),ELIST_SPG(:) !elements sharing the node.
!		integer::Property=0 !for SPG, Property=1 suggesting that the node is on the seepage surface.

	end type
	integer::nnum !�ڵ���
	type(node_tydef),allocatable::node(:)
	real(DPN),ALLOCATABLE::GNODE(:,:) !GHOST NODE
	INTEGER::NGNODE=0
    
	type element_tydef
		integer::nnum,NEDGE=0,NFACE=0,NTET=0 !node numbers of this element
		integer,allocatable::node(:),EDGE(:),FACE(:),TET(:) !��Ԫ�Ľڵ�,TETΪ��Ԫ��ϸ�ֵ�Ԫ����±ꡣ
        integer,allocatable::node2(:) !Ϊ����Bar��Beam��Ԫ�ĺ���Ϊ��Ԫ���ڵĽڵ��ţ�����ת����ʵ�������嵥Ԫ�����,����ԪΪzt4_spg,��zt6_spgʱ��node2ָ��gnode��
        !��et=wellboreʱ��node2�洢��wellbore��Ԫ��3��4�ڵ�Ϊ�ߵ�3ά��Ԫ 
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
		integer::ec=0 !element class
		integer::nspr=0 ! the patch number sharing the element
		integer::layer=1 !element layer, for horizontal seepage analysis; 
						!for sectional seepage analysis, =1, all nodes are under water table,=-1,all nodes are above water talbe,=0 cross the water table
		integer::ifreedof=-1 !>0,��½ӵ�,ָ��freedof(ifreedof)
		
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
		real(kind=DPN)::property(6)=0.0D0  !for spg problem and iniflux is used, property(3)=low lamda,property(2)=up lamda
        !for wellbore element, property(1-3), element hydraulic conductivity, property(4) surround angle. (5)=WELL SKIN RESISTANCE
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
											 !3. sigman(tension is +ve ) 4 tn (clockwise is +ve),5. tnx, 6. tny 
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
        INTEGER,ALLOCATABLE::WELL_SP3(:),WELL_SP4(:) !USE forlab WELLBORE ELEMENT
        REAL(KIND=DPN),ALLOCATABLE::WELL_SP3_R(:),WELL_SP4_R(:)
        INTEGER::NWSP3=0,NWSP4=0		
		
		!real(kind=DPN),allocatable::lamda(:) !lamda(ngp)=1.0(by default),relative k. 
		!additional variables for upper bound analysis
		real(kind=DPN),allocatable::A12(:,:)  ! for UBZT4 is A23
		real(kind=DPN),allocatable::X2(:)
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
		integer::isdual=0
		integer::ssp_onepile=0 !ֻ��SSP��Ԫ�ڵ������ã���ʾ��������Ƿ�������������һ���ְ�׮��,�������������á�
								!/=0,�����ڵ����ְ�׮�ϣ�=0�����á�
		integer::isincrement=0 !���Ϊ1�������value��������������ȫ����Ĭ�ϣ�
	end type
	type(bc_tydef),allocatable::bc_disp(:),bc_load(:),bf(:),NSeep(:),IniValue(:),CFN(:)
	integer::bd_num=0,bl_num=0,bfnum=0,NumNSeep=0,Niniv=0,NCFN=0	
	real(kind=DPN),allocatable::iniValueDof(:) 
    
    type hinge_typef
        integer::element,node,dof=7 !��Ԫ�ţ��ڵ�ţ����ɶ�
		integer::newnode=0 !ָ����ڵ�NodeΪ�غϵĽڵ�node(newnode)��
    endtype
    type(hinge_typef),allocatable::FreeDOF(:)
    integer::Nfreedof=0
	
	integer,allocatable::IsBCChange(:)
		

	!material 
	type mat_tydef
		integer::type=0
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
		real(kind=DPN)::force_tol=1.D-2 !tolerance permitted for unbalanced residual force to be vilated
		real(kind=DPN)::disp_tol=1.D-3 !tolerance permitted for residual displacement to be vilated
		integer::output=1 !1:output for mod(increments,output)=0 iteration
		integer::i2ncal=4,I2NWEIGHT=1  !method mapping variables in integration points to nodes
			!SHT=1,spr=2,AVG=3,EXP=4
		logical::issym=.true.
!		logical::issteady=.true. 
		logical::datapaking=.true. !.t.=point,.f.=block
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
        integer::slidedirection=right,slope_isTensionCrack=1,ISSLOPEPA=0
        real(kind=DPN)::slope_kscale=1.D0,slope_kbase=1.D0,slope_kratio=10
	end type
	type(solver_tydef)::solver_control
	
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
        integer::system=0
		character(1024)::zonetitle=""
		character(64)::grouptitle=""
        integer::coupleset=-1
        integer::isini=0 !>0
        integer::sf=0 !stepfun
		integer::mesh_share_id=0,out_mesh=.true. !for output tecplot.
		
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
	type(outvar_tydef)::outvar(150)	
	integer::vo(150)=0,nvo=0
	
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
        integer::nnode=0,isSUMQ=0 !isSumQ/=0,ouput the sum of nodal discharges in node(). 
        integer,allocatable::node(:)        
    end type
    type(out_data_typdef),allocatable::DataPoint(:)
    integer::NDataPoint=0
	

	
    !type rcd_set_tydef
    !    integer::rcd
    !    integer::enums,enume
    !end type
    !type(rcd_set_tydef),allocatable::rcdset(:)
    !integer::nrcdset=0
    !integer,allocatable::rcd(:,:)
    !integer::nrcd=0

	
	character(1024)::title,resultfile,resultfile1,resultfile2,resultfile3,resultfile21,resultfile22,EXCAMSGFILE,EXCAB_BEAMRES_FILE,&
					EXCAB_STRURES_FILE,EXCAB_EXTREMEBEAMRES_FILE,SLOPE_FILE,helpfile
	INTEGER::DATAPOINT_UNIT=29
	integer::datapacking=1	!=1,point format:{x1,y1,z1},{x2,y2,z2},..., . (Default Format)
						 != 2, block format, {x1,x2,...},{y1,y2,...},{z1,z2,...}
	integer::ndimension=2	!nodal dimension,
						!=1, input x
						!=2, input x and y
						!=3 input x,y and z
							
	integer::ndof=0 !total dof number
	integer,allocatable::bw(:) 	!bw(i): firstly, it is the column number of the most left entry in the i row. and later, it is the bandwidth of the total matrix
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

	integer::nnz=0
	integer,allocatable::irow(:),jcol(:),Lmre(:),adrn(:) !Lmre(i): the column number of the most rigth entry in the i row.
	integer,allocatable::ROWINDEX(:) !rowindex:�ܸ���ÿһ�е�һ������Ԫ�����ܸ������е�λ��;
											!bw: һ��ʼΪ�ܸ���ÿһ�е�һ������Ԫ�����ܸ��е��кţ����Ϊ�ܸ���ÿһ�����һ������Ԫ�����ܸ������е�λ��;
											!Lmre:�ܸ���ÿһ�����һ������Ԫ�����ܸ��е��кš�
											!ֵ��ע����ǣ���default solver�У��ܸպ�����Ԫ��������ÿ�е�һ������Ԫ�����һ������Ԫ֮�䣩������mkl������洢��ʽ�У��ܸղ�����Ԫ��
											!adrn(i):��default solver�е��ܸ������з���Ԫ i ��mkl��ʽ�ܸ������λ�á�
	
	
	real(kind=DPN),allocatable::DIAGLKM(:) !STORED DIAGONAL ELEMENT OF THE KM BEFORE FACTORIZATION FOR FACTOR UPDATING.
	integer,allocatable::diaglkmloc(:) !�ܸ��жԽ���Ԫ�����ܸ������е�λ�á�

	!one dimensional linear field function
	real(kind=DPN)::LF1D(0:maxilf,2)  !LF1D(:,1)=k,LF1D(:,2)=c. then y=k*x+c

!	integer::NLF1D=0 

	!additional variables for upper bound analysis
	integer::ncons=0 !number of constraints imposed by elements
	integer::nuvar=0 !uncontraint variables in linear programming
	integer::nzone_tec=0 !the number of zone in output file in tecplot format.
	integer::varsharezone_tec=0 !the zone number whose variables are shared with others in output file in tecplot format.
	
	real(kind=DPN)::NormL=0.0
	
	real(kind=DPN)::minNPH=1e20
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
	INTEGER::NNODALQ=0
	integer,allocatable::bfgm_step(:)
    integer::mpi_rank = 0, mpi_size = 1, mpi_ierr
    LOGICAL::ISINISEDGE=.FALSE.
    
    INTERFACE
         PURE subroutine INVARIANT(stress,inv)
	        implicit none
	        real(8),intent(in)::stress(6)
	        real(8),intent(out)::inv(3)
	
        end subroutine 
		
        PURE SUBROUTINE KSITA_MC_C2(KSITA,LODE,SITA,PHI)
            !LODE IN RAD,SITA IN RAD, TRANSITIONAL ANGLE.
	        IMPLICIT NONE
	        REAL(8),INTENT(IN)::LODE,SITA,PHI
	        REAL(8),INTENT(OUT)::KSITA(6)
        END SUBROUTINE
		PURE subroutine yieldfun(vyf,mat,inv,ayf,ev,ISTEP)
			!use solverds	
			implicit none
			INTEGER,INTENT(IN)::MAT,AYF,ISTEP
			REAL(8),INTENT(IN)::INV(3),EV
			REAL(8),INTENT(OUT)::VYF
		END SUBROUTINE
		
		PURE SUBROUTINE MC_KSITA(LODE,SITA,PHI,A,B,KSITA,D_KSITA)
			IMPLICIT NONE
			REAL(8),INTENT(IN)::LODE,SITA,A(2),B(2),PHI
			REAL(8),INTENT(OUT)::KSITA,D_KSITA
		END SUBROUTINE
		
		PURE subroutine deriv_yf(dyf,dywi,m)
			!use solverds
			implicit none
			REAL(8),INTENT(IN)::dywi(3),M(6,3)
			REAL(8),INTENT(OUT)::DYF(6)
		END SUBROUTINE
		
		PURE subroutine deriv_sinv(m,stress)
			implicit none
			REAL(8),INTENT(IN)::STRESS(6)
			real(8),INTENT(OUT)::m(6,3)
		END SUBROUTINE
		
		PURE subroutine deriv_yf_with_inv(dywi,inv,mat,ayf,ev,ISTEP)
			
			implicit none
			INTEGER,INTENT(IN)::MAT,AYF,ISTEP
			real(8),INTENT(IN)::inv(3),EV
			real(8),INTENT(OUT)::dywi(3)
			
		END SUBROUTINE
		
		PURE subroutine deriv_qf_with_inv(dywi,inv,mat,ayf,ev,ISTEP)
		
			implicit none
			INTEGER,INTENT(IN)::MAT,AYF,ISTEP
			real(8),INTENT(IN)::inv(3),EV
			real(8),INTENT(OUT)::dywi(3)	
			
		END SUBROUTINE
        
        PURE FUNCTION PARA_MC_CLAUSEN(MATID,ISTEP) RESULT(PARA)
            IMPLICIT NONE
            INTEGER,INTENT(IN)::MATID,ISTEP
            REAL(8)::PARA(3)
        END FUNCTION
		
	END INTERFACE    
   
    
    contains
    
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
    
	PURE function pi()
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
pure function DINV(E,V,ND) result(C)
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

FUNCTION KM_WELLBORE(A1,A2,A3) RESULT(KM)
    IMPLICIT NONE
    REAL(8),INTENT(IN)::A1,A2,A3
    REAL(8)::KM(4,4)
    
    
    KM=0.D0
    KM(1,1)=A1+A3
    KM(2,2)=A1+A2
    KM(3,3)=A2
    KM(4,4)=A3
    KM(1,2)=-A1
    KM(2,1)=-A1
    KM(2,3)=-A2
    KM(3,2)=-A2
    KM(1,4)=-A3
    KM(4,1)=-A3

END FUNCTION



REAL(8) FUNCTION fD_PF(RE,KR,REW) !darcy-friction for pipe flow
    IMPLICIT NONE
!function: calculate darcy-friction for pipe flow,if ReW>0,then it is a porous pipe flow
!Re, Reynolds number (unitless);
!Kr, relative roughness (unitless)
!ReW, wall Reynolds number (unitless);
!REF:[1] Fang X, Xu Y, Zhou Z. New Correlations of Single-Phase Friction Factor for Turbulent Pipe Flow and Evaluation of Existing Single-Phase Friction Factor Correlations[J]. Nuclear Engineering and Design, 2011, 241(3): 897-902. 
![2] Ouyang L-B, Arbabi S, Aziz K. General Wellbore Flow Model for Horizontal, Vertical, and Slanted Well Completions[J]. SPE Journal 1998, 3(2): 124~133.
    REAL(8),INTENT(IN)::RE,KR
    REAL(8),INTENT(IN),OPTIONAL::REW
    REAL(8)::LAMDA1,REW1,FC1
    
    IF(RE<3000) THEN
        LAMDA1=64./RE    
    ELSEIF(KR>0.D0) THEN        
        LAMDA1=1.613*(LOG(0.234*(KR)**1.1007-60.525/RE**1.1105+56.291/RE**1.0712))**(-2.0D0)
    ELSE
        LAMDA1=0.25*(LOG10(150.39/RE**0.98865-152.66/RE))**(-2.D0)
    ENDIF
    
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
REAL(8)::T1

IF(.NOT.PRESENT(T)) THEN
    T1=25
ELSE
    T1=T
ENDIF


Re_W=V*D/Vk(T1)

    
ENDFUNCTION


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

    




