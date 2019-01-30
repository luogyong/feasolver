module meshDS

	integer::maxstk, &
					maxnnode=100000, &
					newNodemax=10000, &
					maxnedge=300000, &
					maxnadjlist=10, &
					maxncedge=1000, &
					maxnelement=210000, &
					soillayer=0, &
					maxntetelt=10000
    real(8)::um
	parameter(maxstk=1000,um=1e15)

	type point_tydef
     integer::number,bw !��ı��,semibandwidth������,�������ӱ�š�
	   real(8) x,y,z !�������
	   real(8)::s  !�õ㸽���ĵ�Ԫ��С
	   INTEGER::SUBBW  !=-999 ,the node is dead.
	   integer::layer=0,onbdy=0 !nodal layer
	   
	end type
	type(point_tydef),allocatable,target::node(:)
	type(point_tydef),allocatable::cp(:)  !control point
	integer::nnode=0 !��¼�ڵ�������Ϊ�Ż�ǰ����

	                            										
	type BP_tydef !bounded points
	   type(point_tydef),pointer::npt=>null()
	   type(BP_tydef),pointer::next=>null()  !previousָ��ǰһ����Ԫ��nextָ���һ����Ԫ
	end type

	type constrainline_tydef
	   real(8),pointer::conpoint(:,:)=>null() !�������ϸ�������꣨x,y��
	   integer::flag,hole  !�Ƿ�պ�,�Ƿ������� 0Ϊ��1Ϊ��
	   integer::num,material=0  !�����ߵ����Ŀ,materialֻ�Է���ǽ����Ϊ��Ϊ����ǽ�Ĳ��Ϻ�,
       integer,allocatable::point(:)
	end type
	

	type element_tydef
	   !type(point_tydef),pointer::xy1,xy2,xy3,xy4   !��Ԫ���������
	   logical::isdel=.false.
	   integer::node(15) !6-noded/15-noded triangleelement node in node().
	   integer::et=0           !element type, 3-noded:0, 6-noded triangle:6; 15-noded:15 
	                                      !6-noded prism element: 63; 15-noded prism element:153;
	                                      !4-noded tetrahedral element:43, second order tetrahedral element: 103
	   integer::zn=1,number
	   integer::edge(3)=0,ORIENT(3)=1 !the edge number in edge set.
	   real(8)::property(5)=0    !��Ԫ����ֵ
	   integer::kcd=0,Maxedge=0  !��Ԫ�Ŀɷֶ�,
	   integer::adj(3)=-1 !if =-1,no adjacent element
	end type
!	type(element_tydef),pointer::Ehead,Ept,element,Etail!ǰ�ĸ�ָ�򻮷ֺõĵ�Ԫ��������ָ����������������ε�Ԫ
	type(element_tydef),allocatable::ELT(:)
	integer::nelt=0,Pehead=-1,PEpt=-1,IEpt=-1,ept=-1
	
	type PrismElement_tydef !������Ԫ
		logical::isdel=.false.
	    integer::nnum
	    integer,allocatable::node(:)
	    integer::mother=0 !generated from the elt(mother)
	    integer::et=63   !6-noded prism element: 63; 15-noded prism element:153;
	                                 !4-noded tetrahedral element:43, second order tetrahedral element: 103
	    integer::matlayer=1 !material 
	    integer::nlayer(2)=0 !nlayer(1��2) �ֱ�Ϊ�����͵ײ��ڵ���
    
	end type 
	type(PrismElement_tydef),allocatable::PRMELT(:),TetELT(:)
	integer::NPRMELT=0,NTetelt=0
	
	type physicalGroup_type
		integer::ndim=0
		integer::icl=-1 !��cl�����߽��ж���
		integer::nvseg=-1 !��segment���ж����physicalgroup  
		integer::izone=-1 !��zone���ж���
		integer::ET_GMSH=0 
		character(32)::name,ET
		integer::nelt=0
		integer::ilayer=0 !-1,output all layers. 
		integer,allocatable::elt(:)
		logical,allocatable::isdel(:) 
		integer,allocatable::vseg(:)
	end type
	integer::npg=0
	type(physicalGroup_type),allocatable::physicalgroup(:)
	

	

	type size_point_tydef  !��Ԫ�ߴ���Ƶ�
	    real(8)::x,y,a,d !���꣬�Ȳ�ĵ�һ��Ȳ�
	end type

	type zone_tydef
	   integer::num,k
	   real(8),allocatable::point(:,:)
	   real(8)::xmin=1e15,ymin=1e15,xmax=-1e15,ymax=-1e15
	   integer,allocatable::item(:),trie3n(:),Dise4n(:), trie6n(:),trie15n(:)
	   integer::nitem=0,ntrie3n=0,ndise4n=0,ntrie6n=0,ntrie15n=0  !restriction:nitem=ntrie3n+ndise4n
	   integer,allocatable::bedge(:) !boundary edges id(point to edge)
	   integer::nbe=0 !number of bedge
	end type


	type bc_tydef
	   integer::node=0
	   real(8)::vbc=0     !��Ϊˮͷ�߽���Ϊˮͷֵ����Ϊ�����߽���Ϊ����ֵ(��Ϊ��״ʱ���洢���˵��ֵ�����ڸ����ϵ�ֵ�߾����߳��ȵȱ������в�ֵ����״��ֻ��v(1)�����ã�v(2)=v(1))
	   integer::bt=0 !=0, displacement-like(unknown);
							!=1, force-like(load); 
	                  !=2, signorini boundary(seepage exiting face) 					

	   integer::dof=0 !1Ϊ����,��ʾ�ñ߽������Ը�ǿ͸ˮ�����ã�0����֮��0,Ϊ���ǲ㣬��������Ϊǿ͸ˮ�㣬���Ϊ5��ǿ͸ˮ��

	end type

	type material_tydef
	  real(8)::kx=0,ky=0,u=0,angle=0 !�������������͸ϵ������ˮϵ��,angleΪ��͸ϵ�����������뼸������ļнǡ��Զ�����(�����תΪ����)��Ĭ��ͬ����ʱ�뷽��Ϊ����
	end type

	type ccl_tydef !circular control line
	   integer::point
       integer::hole  !�Ƿ������� 0Ϊ��1Ϊ��
	   real(8)::r
	   !real(8)::s
	end type

	type aupoint_tydef !���񻮷ָ�����
	   integer::point
	   real(8)::s
	end type

	type dataline_tydef
		integer::ni,nj,nnum
		real(8)::v(2,2)  !, v(:,1)=x,v(:,2)=y
		integer,pointer::node(:)
	end type
	type(dataline_tydef),allocatable::ninseg(:)

	type property_tydef
		character(32)::name
		real(8)::value
	end type
	type(property_tydef)::property(10) !one control line has 10 property value at most.
	integer::pro_num
	
	type structmesh_tydef
		integer::izone !zone number
		integer::csl !csl(csl) is the zone boudary 
		integer::nnum ! number of the node
		integer::enum ! number of element
		real(8)::v(4,2) ! 4 vertexes for the zone, Note that, for structure mesh, every zone has 4 vertexes.
		real(8),pointer::node(:,:) !element nodes in this zone
		integer,pointer::element(:,:) !nodes for each element
		integer,pointer::adjelem(:,:) ! adajency table for each element		
	end type
	type(structmesh_tydef),allocatable::structmesh(:) 
	integer::nsm=0
	
	type InsertPoint_3DModel
		real(8)::x,y,z
		real(8)::ratio=0 !�󽻵�ı���ϵ��
		integer::isMesh=0 !Is the apoint a existing node in the mesh? 0=No;
		!>0,Yes,the existing point=node(ismesh)
		integer::InsertingLayer(2)=0 !the two crossing layers generating the point
		integer::isRepeated=0 !�Ƿ�Ϊ�ص㣬0��NO��>0 ,�ص㣬=������ȫ�ĵ�һ������±꣬NIP3DM,
	end type
	type(InsertPoint_3DModel),pointer::InsertPoint_3DM(:)=>null()
	integer::nIP3DM=0
	
	type Edge_tydef	
		integer::num=0
		integer::V(2)=0 !two vetex for each element edge 
		integer::E(2)=-1  !the two elements sharing the element	
		INTEGER::BRECT=0 !zone boudary rectangle face id of out to gmsh volume. 		
	end type
	type(edge_tydef),allocatable::edge(:)
	integer::nedge=0

	type adjlist_tydef
		integer::count=0 !the number of adjacent node 
		integer,pointer::node(:)=>null() !adjacent node
		integer,pointer::edge(:)=>null() !adjacent edge
	end type
	type(adjlist_tydef),allocatable::adjlist(:) !adjacent table for nodes

	type constrained_edge_tydef
		integer::v(2) !edge vertex in node
		integer::cl=-1 !belonged to which constroled line.
		integer::edge=-1 !if existed, it points to edge(edge)
	end type
	type(constrained_edge_tydef),allocatable::cedge(:)
	integer::ncedge=0
	
	type constrained_line_tydef
		integer::nedge=0
		integer,allocatable::edge(:)
	end type
	type(constrained_line_tydef),allocatable::edge_cl(:)
	
	
	Type meminp_tydef
		integer::icl,nnum !cpphead(icl),control line number��control points of icl.
		integer::nvb !ϸ�ֺ�������ϵ�ĸ�����
		integer,allocatable::cp(:) !for meminp2 to input the control points,cp(nnum)
		integer,allocatable::nbc(:) !nbc(nvb)
		real(8),allocatable::elevation(:,:) !elevation(nnum,0:soillayer) ��������icl�ϸ���ĸ̡߳�
		real(8),allocatable::lincof(:,:,:) !lincof(2,n1,0:soillyger) !if csl(icl).flag==1,n1=nnum,==0,n1=nnum-1
		real(8),allocatable::vbc(:,:) !vbc(nvb,0:soillyger)
	end type
	type(meminp_tydef),allocatable::meminp(:),meminp2(:)
	integer::nmeminp=0,nmeminp2=0

	integer,allocatable::segindex(:,:) 
	!�洢ģ���߶εģ����ǵ�Ԫ�ߣ���seg�е�λ�ã���������ccl��ʽ����Ŀ����ߣ���Ϊ����������ϵĵ㲻�������������棩��
	type seg_type
		integer::icl !���߶���cpphead(icl)��λ��
		integer::sv=0,ev=0 !�ڵ������������е�λ�á�
		integer::isA2Z=1 !�Ƿ���˳��1Ϊ˳�򣬼���cpphead(icl)������Ϊsv-ev. 0Ϊ����ev-sv. 
		integer::isT2H=0 !�Ƿ���β�������ĶΡ�=1,yes.
		integer::nnum=0 !ϸ�ֺ���߶εĽڵ������
		integer,allocatable::node(:) !ϸ�ֺ���߶εĽڵ���tnode���±�,�����˽ڵ㡣
		type(BP_tydef),pointer::svp,evp !�˽ڵ���cpphead(icl)�е�λ��
	end type
	type(seg_type),allocatable::seg(:)
	integer::nseg=0,maxnseg=5000	
	
	type geology_point_tydef 
	    integer::node=0
	    integer::isini=0
	    real(8),allocatable::elevation(:)  !elevation(0:nelevation) ,nelevation=soillayer
	                                                             !elevation(i)=-999, indicating the elevation in the node is interpolted from other nodes.
	end type
	type(geology_point_tydef),allocatable::geology(:)
	integer::ngeo=0
	
	
	
	type solver_typdef
		character(8)::ProblemType='SPG'
		character(8)::SolverType='N_R'		
	end type
	type(solver_typdef)::SolverInfo

	type(BP_tydef),pointer::BNhead,BNpt,BP !ȷ��ģ����߽�Ĺؼ���
	type(constrainline_tydef),allocatable::csl(:)
	type(BP_tydef),pointer::cpp,cpptail
	type(BP_tydef),allocatable,target::cpphead(:) !ָ���γɿ����ߵĵ�Ķ��У�Ϊ���������Ƿ��������εı�����׼��
!	type(element_tydef),target::element
 	type(size_point_tydef),allocatable::s_P(:)
	type(zone_tydef),allocatable::zone(:)
	type(material_tydef),allocatable::material(:)
!	type(hqboundary_tydef),allocatable::hqb(:)
	type(ccl_tydef),allocatable::ccl(:)
 	type(aupoint_tydef),allocatable::aupoint(:),aupoint_m(:)
    !type(element_tydef),allocatable,target::sehead(:) !singular element head.
	
	integer::stack(maxstk)

   ! type(element_tydef)::CE(20)
	integer::isnorefined=0
   	integer::ENUMBER=0  !��Ԫ�ĸ���
!	integer::nnode !��Ч�Ľڵ���
    integer::keypn  !key point numbers
    integer::cnn   !control point numbers
	integer::cln  !constrainline numbers
	real(8)::winxul,winxlr,winyul,winylr !for graphic
	real(8)::triangle(2,3) !���������������������
	integer(4)::i4
	integer::sizepoint  !!��Ԫ�ߴ���Ƶ�ĸ�����
	integer::znum=1   !������,�����򻮷����������򻮷�ʱ��zone����Ϊ�ա�
	integer::mnum   !������
	integer::hqnum  !����ˮͷ�ı߽߱�
	integer::ccln   !Բ�ο����ߵĸ���
	integer::aupn
	character(256)::resultfile,checkfile,title,path_name	!����ļ����ļ���
	integer,allocatable::cclincsl(:) !coffwnum(:)��¼����ǽ��csl�����е�λ��,buildingincsl(:)��¼�ṹ������csl�����е�λ��
	integer,allocatable::Noutputorder(:) !��¼node(i)��������˳�����noutputorder(j)=i����ڵ�node(i)���˳��ΪΪj.
	integer::topstk=0,dln=0
	integer::tt=1,slrn  !tt,Ϊ�ܲ���,slrnΪslrn�ĸ���
!	integer,allocatable::edge(:,:) !store two end points of edges of all triangle elements.  edge(nedge,2)
!	integer::nedge=0 !the edges of all triangle element 
	
	integer::maxadj=100  !the maximum adjacent vertice 

	real(8),allocatable::STR(:) ! STEP-TIME RELATIONSHIP,��¼ÿһ���Ĳ�����ֻ��һ���������±��ǲ���
	real(8),allocatable::SLR(:,:)  ! STEP-LOAD RELATIONSHIP, ��¼��һ���ĺ�������.��һ�������±��ǲ������ڶ����±��ǹ�ϵ�������ж��slr��ϵ��Ӧ���ڲ�ͬ�ı߽��ϡ�
	integer,allocatable::SOR(:) ! STEP-OUT RELATIONSHIP ��ʾ�ļ�����Ҫ�����,1,���,0,�������SRT,SLR,SOR�ĵ�һ���±��СӦһ�¡�	
	real(8)::ic=0. !��ʼ������-999��������д��ļ����룬-9999,��ʾ��ʼ����Ϊ�����ʼˮλ��
	integer,allocatable::dataline(:)!ָ��Ҫ������ݵ���������cpphead(i)�е�i��Ҳ��ָ��csl(i)�е�i,��߽��i=0
	integer::issupel=0,bnum=0,inum=0,saupn=0,iskpm=0 !saupn,����Ԫ�ڲ������ڵ�ĸ���;iskbm,�߽���Ƿ�Ҫ��һ��ϸ�֣�0:yes,others:No. default value:yes.
	real(8)::precision=1e-7

	integer,allocatable::isLM(:) ! the islm(i)=a, zone is to be meshed into a limit analaysis grid.
	integer::nislm=0,nninseg=0,isoffset=0
	real(8)::xyscale=1.0D0
	!real(8),allocatable::GeoElevation(:,:) !geoelevation(nnode,nlayer),�洢���ڵ�ز�߳���Ϣ�����϶��´洢��
	!integer::nlayer=1 !soil layer
	real(8),allocatable::elevation(:,:)
	integer,allocatable::nlayer(:,:) !elevation data of each node. the value of 
	! elevation(i,node) is the elevation of nlayer(i,node) 
	integer::nelevation !number of elevation data for each node	
	real(8)::xmin,ymin,xmax,ymax !the extreme value of mesh zone
	integer,allocatable::bnode(:) 
	integer::nbnode=0

end module
