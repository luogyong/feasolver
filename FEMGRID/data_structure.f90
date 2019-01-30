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
     integer::number,bw !点的编号,semibandwidth（带宽）,荷载因子编号。
	   real(8) x,y,z !点的坐标
	   real(8)::s  !该点附近的单元大小
	   INTEGER::SUBBW  !=-999 ,the node is dead.
	   integer::layer=0,onbdy=0 !nodal layer
	   
	end type
	type(point_tydef),allocatable,target::node(:)
	type(point_tydef),allocatable::cp(:)  !control point
	integer::nnode=0 !记录节点数，作为优化前结点号

	                            										
	type BP_tydef !bounded points
	   type(point_tydef),pointer::npt=>null()
	   type(BP_tydef),pointer::next=>null()  !previous指向前一个单元，next指向后一个单元
	end type

	type constrainline_tydef
	   real(8),pointer::conpoint(:,:)=>null() !控制线上各点的坐标（x,y）
	   integer::flag,hole  !是否闭合,是否是留洞 0为否，1为真
	   integer::num,material=0  !控制线点的数目,material只对防渗墙有用为，为防渗墙的材料号,
       integer,allocatable::point(:)
	end type
	

	type element_tydef
	   !type(point_tydef),pointer::xy1,xy2,xy3,xy4   !单元顶点的坐标
	   logical::isdel=.false.
	   integer::node(15) !6-noded/15-noded triangleelement node in node().
	   integer::et=0           !element type, 3-noded:0, 6-noded triangle:6; 15-noded:15 
	                                      !6-noded prism element: 63; 15-noded prism element:153;
	                                      !4-noded tetrahedral element:43, second order tetrahedral element: 103
	   integer::zn=1,number
	   integer::edge(3)=0,ORIENT(3)=1 !the edge number in edge set.
	   real(8)::property(5)=0    !单元属性值
	   integer::kcd=0,Maxedge=0  !单元的可分度,
	   integer::adj(3)=-1 !if =-1,no adjacent element
	end type
!	type(element_tydef),pointer::Ehead,Ept,element,Etail!前四个指向划分好的单元，后两个指向包含插入点的三角形单元
	type(element_tydef),allocatable::ELT(:)
	integer::nelt=0,Pehead=-1,PEpt=-1,IEpt=-1,ept=-1
	
	type PrismElement_tydef !派生单元
		logical::isdel=.false.
	    integer::nnum
	    integer,allocatable::node(:)
	    integer::mother=0 !generated from the elt(mother)
	    integer::et=63   !6-noded prism element: 63; 15-noded prism element:153;
	                                 !4-noded tetrahedral element:43, second order tetrahedral element: 103
	    integer::matlayer=1 !material 
	    integer::nlayer(2)=0 !nlayer(1，2) 分别为顶部和底部节点层号
    
	end type 
	type(PrismElement_tydef),allocatable::PRMELT(:),TetELT(:)
	integer::NPRMELT=0,NTetelt=0
	
	type physicalGroup_type
		integer::ndim=0
		integer::icl=-1 !用cl控制线进行定义
		integer::nvseg=-1 !用segment进行定义的physicalgroup  
		integer::izone=-1 !用zone进行定义
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
	

	

	type size_point_tydef  !单元尺寸控制点
	    real(8)::x,y,a,d !坐标，等差的第一项，等差
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
	   real(8)::vbc=0     !若为水头边界则为水头值，若为流量边界则为流量值(当为线状时，存储两端点的值，处于该线上的值线据沿线长度等比例进行插值，面状，只有v(1)起作用，v(2)=v(1))
	   integer::bt=0 !=0, displacement-like(unknown);
							!=1, force-like(load); 
	                  !=2, signorini boundary(seepage exiting face) 					

	   integer::dof=0 !1为激活,表示该边界条件对该强透水层作用，0，反之。0,为覆盖层，其它各层为强透水层，最多为5层强透水层

	end type

	type material_tydef
	  real(8)::kx=0,ky=0,u=0,angle=0 !两个主方向的渗透系数和贮水系数,angle为渗透系数的主方向与几何坐标的夹角。以度输入(输入后转为弧度)，默认同向，逆时针方向为正。
	end type

	type ccl_tydef !circular control line
	   integer::point
       integer::hole  !是否是留洞 0为否，1为真
	   real(8)::r
	   !real(8)::s
	end type

	type aupoint_tydef !网格划分辅助点
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
		real(8)::ratio=0 !求交点的比例系数
		integer::isMesh=0 !Is the apoint a existing node in the mesh? 0=No;
		!>0,Yes,the existing point=node(ismesh)
		integer::InsertingLayer(2)=0 !the two crossing layers generating the point
		integer::isRepeated=0 !是否为重点，0，NO。>0 ,重点，=与其重全的第一个点的下标，NIP3DM,
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
		integer::icl,nnum !cpphead(icl),control line number，control points of icl.
		integer::nvb !细分后控制线上点的个数。
		integer,allocatable::cp(:) !for meminp2 to input the control points,cp(nnum)
		integer,allocatable::nbc(:) !nbc(nvb)
		real(8),allocatable::elevation(:,:) !elevation(nnum,0:soillayer) ！控制线icl上各点的高程。
		real(8),allocatable::lincof(:,:,:) !lincof(2,n1,0:soillyger) !if csl(icl).flag==1,n1=nnum,==0,n1=nnum-1
		real(8),allocatable::vbc(:,:) !vbc(nvb,0:soillyger)
	end type
	type(meminp_tydef),allocatable::meminp(:),meminp2(:)
	integer::nmeminp=0,nmeminp2=0

	integer,allocatable::segindex(:,:) 
	!存储模型线段的（不是单元边）在seg中的位置（不包括以ccl形式输入的控制线，因为这类控制线上的点不在输入数组里面）。
	type seg_type
		integer::icl !该线段中cpphead(icl)中位置
		integer::sv=0,ev=0 !节点在输入数组中的位置。
		integer::isA2Z=1 !是否是顺序，1为顺序，即在cpphead(icl)链表中为sv-ev. 0为反序，ev-sv. 
		integer::isT2H=0 !是否是尾首相连的段。=1,yes.
		integer::nnum=0 !细分后此线段的节点个数。
		integer,allocatable::node(:) !细分后此线段的节点在tnode的下标,包括端节点。
		type(BP_tydef),pointer::svp,evp !此节点在cpphead(icl)中的位置
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

	type(BP_tydef),pointer::BNhead,BNpt,BP !确定模型外边界的关键点
	type(constrainline_tydef),allocatable::csl(:)
	type(BP_tydef),pointer::cpp,cpptail
	type(BP_tydef),allocatable,target::cpphead(:) !指向形成控制线的点的队列，为检查控制线是否在三角形的边上做准备
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
   	integer::ENUMBER=0  !单元的个数
!	integer::nnode !有效的节点数
    integer::keypn  !key point numbers
    integer::cnn   !control point numbers
	integer::cln  !constrainline numbers
	real(8)::winxul,winxlr,winyul,winylr !for graphic
	real(8)::triangle(2,3) !外包三角形三个顶点坐标
	integer(4)::i4
	integer::sizepoint  !!单元尺寸控制点的个数。
	integer::znum=1   !区域数,分区域划分数。单区域划分时，zone不能为空。
	integer::mnum   !材料数
	integer::hqnum  !给定水头的边边界
	integer::ccln   !圆形控制线的个数
	integer::aupn
	character(256)::resultfile,checkfile,title,path_name	!输出文件的文件名
	integer,allocatable::cclincsl(:) !coffwnum(:)记录防渗墙在csl数组中的位置,buildingincsl(:)记录结构物在在csl数组中的位置
	integer,allocatable::Noutputorder(:) !记录node(i)数组的输出顺序，如果noutputorder(j)=i，则节点node(i)输出顺序为为j.
	integer::topstk=0,dln=0
	integer::tt=1,slrn  !tt,为总步数,slrn为slrn的个数
!	integer,allocatable::edge(:,:) !store two end points of edges of all triangle elements.  edge(nedge,2)
!	integer::nedge=0 !the edges of all triangle element 
	
	integer::maxadj=100  !the maximum adjacent vertice 

	real(8),allocatable::STR(:) ! STEP-TIME RELATIONSHIP,记录每一步的步长。只有一个。数组下标是步数
	real(8),allocatable::SLR(:,:)  ! STEP-LOAD RELATIONSHIP, 记录第一步的荷载因子.第一个数组下标是步数。第二个下标是关系数，可有多个slr关系，应用在不同的边界上。
	integer,allocatable::SOR(:) ! STEP-OUT RELATIONSHIP 标示哪几步是要输出的,1,输出,0,不输出。SRT,SLR,SOR的第一个下标大小应一致。	
	real(8)::ic=0. !初始条件。-999在求解器中从文件读入，-9999,表示初始条件为输入初始水位。
	integer,allocatable::dataline(:)!指出要输出数据的数据线在cpphead(i)中的i，也是指在csl(i)中的i,外边界的i=0
	integer::issupel=0,bnum=0,inum=0,saupn=0,iskpm=0 !saupn,超单元内部保留节点的个数;iskbm,边界点是否要进一步细分，0:yes,others:No. default value:yes.
	real(8)::precision=1e-7

	integer,allocatable::isLM(:) ! the islm(i)=a, zone is to be meshed into a limit analaysis grid.
	integer::nislm=0,nninseg=0,isoffset=0
	real(8)::xyscale=1.0D0
	!real(8),allocatable::GeoElevation(:,:) !geoelevation(nnode,nlayer),存储各节点地层高程信息，从上而下存储。
	!integer::nlayer=1 !soil layer
	real(8),allocatable::elevation(:,:)
	integer,allocatable::nlayer(:,:) !elevation data of each node. the value of 
	! elevation(i,node) is the elevation of nlayer(i,node) 
	integer::nelevation !number of elevation data for each node	
	real(8)::xmin,ymin,xmax,ymax !the extreme value of mesh zone
	integer,allocatable::bnode(:) 
	integer::nbnode=0

end module
