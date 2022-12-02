module meshDS
    use quicksort
	integer::maxstk, &
					maxnnode=100000, &
					newNodemax=10000, &
					maxnedge=300000, &
					maxnadjlist=10, &
					maxncedge=1000, &
					maxnelement=210000, &
					soillayer=0, &
					maxntetelt=10000,&
                    modeldimension=3,&
                    INPMethod=1,& !=1,linear, =0, membrance
                    Zorder=0,& !=0土层高程从下往上输入,=1,反之。
                    poly3d=1,&
                    ismerged=1,&
                    iscompounded=0,&
                    ISMESHSIZE=1
    integer,parameter:: ET_PRM=63,&
                    ET_TET=43,Linear=1,Membrance=0
    logical::ismeminpdone=.false.
    real(8)::um
    
	parameter(maxstk=2000,um=1e15)
    
    INTERFACE ENLARGE_AR
        MODULE PROCEDURE I_ENLARGE_AR,R_ENLARGE_AR,NODE_ENLARGE_AR,ELEMENT_ENLARGE_AR,&
                         edge_enlarge_ar,ADJLIST_ENLARGE_AR,AR2D_ENLARGE_AR,AR2D_ENLARGE_AR2
    END INTERFACE 

	type point_tydef
     integer::number,bw !点的编号,semibandwidth（带宽）,荷载因子编号。
	   real(8)::x,y,z=0.0D0 !点的坐标
	   real(8)::s  !该点附近的单元大小
	   INTEGER::SUBBW  !=-999 ,the node is dead.
	   integer::layer=0,onbdy=0,marker,mother=0 !nodal layer
       integer::havesoildata=0 ![0,1,2] 0,no;1,partially;2,completely yes.(all soillayers have elevation(no -999)) 
       integer::iptr=0,isb=0,nat=0 !iptr=i,表明此节点在高程上i节点重合；!isb=1,zone boundary node;=2,model boundary node;=0,zone inside node
       real(8),allocatable::elevation(:),we(:),at(:) !we(1:2)=防渗墙顶和底高程
	end type
	type(point_tydef),allocatable,target::node(:)
	type(point_tydef),allocatable::cp(:)  !control point
	integer::nnode=0,tnode=0 !单层节点数，总节点数
    integer,allocatable::kp(:)
	                            										
	type BP_tydef !bounded points
	   type(point_tydef),pointer::npt=>null()
	   type(BP_tydef),pointer::next=>null()  !previous指向前一个单元，next指向后一个单元
	end type

	type constrainline_tydef
	   real(8),pointer::conpoint(:,:)=>null() !控制线上各点的坐标（x,y）
	   integer::flag=0,hole=0  !是否闭合,是否是留洞 0为否，1为真
	   integer::num,material=0  !控制线点的数目,material只对防渗墙有用为，为防渗墙的材料号,
       integer,allocatable::point(:)
	end type
	
   
	type element_tydef
	   !type(point_tydef),pointer::xy1,xy2,xy3,xy4   !单元顶点的坐标
	   logical::isdel=.false.
       INTEGER::NNUM=3,nat=0
	   integer::node(15) !6-noded/15-noded triangleelement node in node().
	   integer::et=0           !element type, 3-noded:0, 6-noded triangle:6; 15-noded:15 
	                                      !6-noded prism element: 63; 15-noded prism element:153;
	                                      !4-noded tetrahedral element:43, second order tetrahedral element: 103
                                !21,line,element
                                !42,quadrangle element
                                !-1,不连续单元或防渗墙单元
	   integer::zn=0,number,ZN2=-1,mat=-1
	   integer::edge(4)=0,ORIENT(4)=1 !the edge number in edge set.
	   real(8)::property(5)=0    !单元属性值, (1)=AREA
	   integer::kcd=0,Maxedge=0  !单元的可分度,
	   integer::adj(4)=-1 !if =-1,no adjacent element
       INTEGER::MOTHER=0,iLAYER=0,ISMODEL=1
       integer,allocatable::selt(:) !selt(i)=0,no element gen at the ith layer.
       real(8),allocatable::bbox(:,:),at(:)
	end type
!	type(element_tydef),pointer::Ehead,Ept,element,Etail!前四个指向划分好的单元，后两个指向包含插入点的三角形单元
	type(element_tydef),target,allocatable::ELT(:)
	integer::nelt=0,Pehead=-1,PEpt=-1,IEpt=-1,ept=-1
	
	!type PrismElement_tydef !派生单元
	!	logical::isdel=.false.
	!    integer::nnum
	!    integer,allocatable::node(:)
	!    integer::mother=0 !generated from the elt(mother)
	!    integer::et=63   !6-noded prism element: 63; 15-noded prism element:153;
	!                                 !4-noded tetrahedral element:43, second order tetrahedral element: 103
	!    integer::mat !material 
	!    integer::nlayer(2)=0 !nlayer(1，2) 分别为顶部和底部节点层号
 !               
	!end type 
	!type(PrismElement_tydef),allocatable::PRMELT(:),TetELT(:)
	!integer::NPRMELT=0,NTetelt=0
	
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
	

	type modelgroup_type
        integer::izone=1,ilayer=0,mat=1,coupleset=-1,sf=0 
        !对于2D的面单元,ilayer=i,表示第i个高程面。
        !对于2D的面单元,ilayer=i,表示第i层单元。注意，n层单元，有n+1高程面，第i层单元由i-1及i层高程面组成(棱柱单元)。
        character(32)::et=''
        character(32)::name=''
    endtype
    type(modelgroup_type),allocatable::model(:)
    integer::nmgroup=0

	type size_point_tydef  !单元尺寸控制点
	    real(8)::x,y,a,d !坐标，等差的第一项，等差
	end type
    type ar2d_tydef
        integer::nnum=0,IVOL=0  !IVOL=GMSH几何体编号
        REAL(8)::INSIDEPT(3)
        integer,allocatable::node(:)
    endtype
    type ar2d_tydef2
        integer::nnum=0
        integer,allocatable::node(:),edge(:)
    endtype    
    TYPE ZONE_LAYGER_VOLUME_TYDEF
        INTEGER::NVOL=0        
        type(ar2d_tydef),allocatable::bface(:)        
    ENDTYPE
	type zone_tydef
	   integer::num,k,OutGmshType=1,iElevation=1 !OutGmshType=1 Physical Volume only; =2 Physical Surface only; =3, both;if OutGmshType=2/3, iElevation指定输出哪个高程的面。       
	   real(8),allocatable::point(:,:)
	   real(8)::xmin=1e15,ymin=1e15,xmax=-1e15,ymax=-1e15
	   integer,allocatable::item(:),trie3n(:),Dise4n(:), trie6n(:),trie15n(:)
	   integer::nitem=0,ntrie3n=0,ndise4n=0,ntrie6n=0,ntrie15n=0  !restriction:nitem=ntrie3n+ndise4n
	   integer,allocatable::bedge(:),bnode(:),mat(:),cp(:) !boundary edges id(point to edge),mat=每层土的材料号
	   integer::nbe=0,nbn=0
       integer::format=0 !=0，区域由边界点定义；=1，区域由内点定义
       !当生成PRM和TET时，以下各量会用到
       
       integer,allocatable::NPrm(:),NTet(:),PRM(:,:),TET(:,:),ISMODEL(:)
       !区内每层土的棱柱单元和四面体单元的个数,PRM,每层棱柱单元在PRM_ELT中的下标,TET类似。
       !ISmodel=1,yes output in element form. 
       type(ar2d_tydef),allocatable::bface(:),hface(:),vface(:) !分别为每区每层的所有边界面，外露的水平面和竖直面。
       type(ZONE_LAYGER_VOLUME_TYDEF),allocatable::vbface(:) !faces for each volume in each layer of each zone
       character(32),allocatable::solver_et(:),NAME(:)
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
        character(512)::cvalue=''
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
        INTEGER::ISZONEBC=-1,ISCEDGE=0
        integer::isxyoverlap=0,marker=0 !在xy平面上，线段的两端点是否重合,主要是无厚度防渗墙单元的[2,3]和[4,1]节点的边
	end type
	type(edge_tydef),allocatable::edge(:)
	integer::nedge=0

	type adjlist_tydef
		integer::count=0,ipt1=0 !the number of adjacent node 
		integer,pointer::node(:)=>null() !adjacent node
		integer,pointer::edge(:)=>null() !adjacent edge
        !integer,pointer::elt(:)=>null() !adjacent element
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
		integer,allocatable::nbc(:),niseg(:) !nbc(nvb)
		real(8),allocatable::elevation(:,:) !elevation(nnum,0:soillayer) ！控制线icl上各点的高程。
		real(8),allocatable::lincof(:,:,:) !lincof(2,n1,0:soillyger) !if csl(icl).flag==1,n1=nnum,==0,n1=nnum-1
		real(8),allocatable::vbc(:,:) !vbc(nvb,0:soillyger)
	end type
	type(meminp_tydef),allocatable::meminp(:),meminp2(:)
	integer::nmeminp=0,nmeminp2=0

	integer,allocatable::segindex(:,:) 
	!存储模型线段的（不是单元边）在seg中的位置（不包括以ccl形式输入的控制线，因为这类控制线上的点不在输入数组里面）。
	type seg_type
        private
		integer::icl,isbg=0 !该线段中cpphead(icl)中位置
		integer::sv=0,ev=0 !节点在输入数组中的位置。
        logical::isini=.false.
		!integer::isA2Z=1 !是否是顺序，1为顺序，即在cpphead(icl)链表中为sv-ev. 0为反序，ev-sv. 
		!integer::isT2H=0 !是否是尾首相连的段。=1,yes.
		integer::nnum=0,nedge=0 !细分后此线段的节点个数。
		integer,allocatable::node(:),edge(:) !细分后此线段的节点在tnode的下标,包括端节点。
		!type(BP_tydef),pointer::svp,evp !此节点在cpphead(icl)中的位置
    contains
        procedure::setparas=>seg_set_parameters
        procedure::getparas=>seg_get_parameters
        procedure::initialize=>seg_initialize
        procedure::get_node=>seg_get_ordered_node
        procedure::get_edge=>seg_get_ordered_edge
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
	
    TYPE::SEARCHZONE_TYDEF(dim)
        integer,len::dim=2
        INTEGER::NEL=0
        integer::NDX(dim) !=1 !NDX=DIVISION IN X Y AND Y AXIS
        REAL(8)::BBOX(2,dim) !MIN,MAX
        INTEGER,ALLOCATABLE::ELEMENT(:)
    ENDTYPE
    TYPE(SEARCHZONE_TYDEF),ALLOCATABLE::SEARCHZONE(:)    
    INTEGER::NSZONE=1
    logical::isszinit=.false.
	
	
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
	type(zone_tydef),TARGET,allocatable::zone(:)
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
	integer::znum=0   !区域数,分区域划分数。单区域划分时，zone不能为空。
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
	real(8)::xmin,ymin,xmax,ymax,zmin=1.d20,zmax=-1.d20 !the extreme value of mesh zone
	integer,allocatable::bnode(:) 
	integer::nbnode=0
    
    character(512)::COPYFILE(100)
	INTEGER::NCOPYFILE=0
    logical::issoilinterpolated=.false.
    
    contains
    
    function readline(nunitr) result(line)
        implicit none
    
        ! Reads line from unit=nunitr, ignoring blank lines
        ! and deleting comments beginning with an exclamation point(//)
        integer,intent(in)::nunitr
        character (len=:),allocatable:: line
        character(len=2048)::line1
        integer::ios,ipos
        do  
          read(nunitr,'(a)', iostat=ios) line1      ! read input line
          if(ios /= 0) return
          line1=adjustl(line1)
          ipos=index(line1,'//')
          if(ipos == 1) cycle
          if(ipos /= 0) line1=line1(:ipos-1)
          if(len_trim(line1) /= 0) exit
        end do
        line=trim(line1)
        return

    end function readline 
    
    subroutine write_help(helpstring,unit)
    USE IFQWIN
    implicit none
    character(*),intent(in)::helpstring
    !class(ZoneBC_type)::this
    integer,intent(in),optional::unit 
    integer(4)::oldcolor
    
    oldcolor = SETTEXTCOLOR(INT2(10))
	write(*,'(A)') trim(helpstring)
	oldcolor = SETTEXTCOLOR(INT2(15)) 

!10 FORMAT('BC的输入格式为:' \ &
!            '1)组数(nbc);' \ & 
!            '2)[BASE,IGROUP,ILAYER,LOC,BCTYPE,NDIM,DOF,STEPFUNC(0),VALUE,[A,B,C,D]]*nbc' \ &
!            'NOTE:' \ &
!            '   a) NDIM=[0,1,2,3] 分别表示点、线、面、体的约束. 因为边界条件（如位移，水头等）不具叠加性，所以NDIM取值对边界条件无影响.' \ &
!            '   b) DOF=4(水头)，且VALUE=-999,表示节点的水头边界值为节点高程值。(模拟暴雨工况。)' \ &
!            '   c) DOF(I)=1,2,...,7,分别表示约束X,Y,Z,H,MX,MY,MZ.' \ &
!            '   d) [A,B,C,D] 平面场方程计算参数，边界值VALUE=A*X+B*Y+C*Z+D.' \ &
!            '   e) BCTYPE=[0,1,2] 分别表示位移类，力类，出溢面类边界。' \ &
!            '   f) BASE=[1,0] 分别表示边界是定义在ZONE上，还是在CL (control line)上。' \ &
!            '   g) LOC=[1,0] 分别表示ILAYER指向的是第ilayer高程面，否则指向的是第ilayer单元层。' \ &
!            '   h) IGROUP, IGROUP为zone的下标(ZONE(IGROUP))')
    endsubroutine
    
   !把字符串中相当的数字字符(包括浮点型)转化为对应的数字
   !如 '123'转为123,'14-10'转为14,13,12,11,10
   !string中转化后的数字以数组ar(n1)返回，其中,n1为字符串中数字的个数:(注　1-3转化后为3个数字：1,2,3)
   !nmax为数组ar的大小,string默认字符长度为256。
   !num_read为要读入数据的个数。
   !unit为文件号
   !每次只读入一个有效行（不以'/'开头的行）
   !每行后面以'/'开始的后面的字符是无效的。
   subroutine  strtoint(unit,ar,nmax,n1,num_read,isall)
	  implicit none
	  integer::i,j,k,strl,ns,ne,n1,n2,n3,step,nmax,num_read,unit,ef,n4
	  logical::tof1,tof2
	  real(8)::ar(nmax),t1
	  character*256::string
	  character*16::substring
	  character*16::legalC
	  logical::isall
	  optional::isall

	  LegalC='0123456789.-+eE*'


	  n1=0
	  ar=0

	  do while(.true.)
		 read(unit,'(a256)',iostat=ef) string
		 if(ef<0) then
			print *, 'file ended unexpended. sub strtoint()'
			stop
		 end if
		 string=adjustL(string)
		 strL=len_trim(string)
		 if(strL==0) cycle


		 if(string(1:2)/='//') then
			
			!每行后面以'/'开始的后面的字符是无效的。
			if(index(string,'//')/=0) then
				strL=index(string,'//')-1
				string=string(1:strL)
				!string=adjustL(adjustR(string))
				strL=len_trim(string)
			end if
					
			i=1			
			do while(i<=strL)
			   if(index(legalc,string(i:i))/=0) then
				  ns=i
			
				  if(i<strL) then
					 do j=i+1,strl
						if(index(legalc,string(j:j))==0) then
						   ne=j-1
						   exit
						end if
						if(j==strL) ne=strL
					 end do
				  else
					 ne=i
					 j=i
				  end if

				  substring=string(ns:ne)
				  n2=len_trim(substring)
				  n3=index(substring,'-')
				  n4=index(substring,'*')
				    tof1=.false.
				    if(n3>1) then
				         tof1=(substring(n3-1:n3-1)/='e'.and.substring(n3-1:n3-1)/='E')
				    end if				  
				  if(tof1) then !处理类似于'1-5'这样的形式的读入数据
					 read(substring(1:n3-1),'(i8)') ns
					 read(substring(n3+1:n2),'(i8)') ne
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
				             tof2=(substring(n4-1:n4-1)/='e'.and.substring(n4-1:n4-1)/='E')
				         end if					 
					 if(tof2) then !处理类似于'1*5'(表示5个1)这样的形式的读入数据
						 read(substring(1:n4-1),*) t1
						 read(substring(n4+1:n2),'(i8)') ne
						 ar((n1+1):(n1+ne))=t1
						 n1=n1+ne
					 else
						n1=n1+1
						read(string(ns:ne),*) ar(n1)
					 end if
			
				  end if
				
				  i=j+1
			   else
				  i=i+1
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

subroutine Err_msg(cstring)
	use dflib
	implicit none
	character(*)::cstring
	character(46)::term
	integer(4)::msg

	term="No such Constant:  "//trim(cstring)
	msg = MESSAGEBOXQQ(term,'Mistake'C,MB$ICONSTOP.OR.MB$OK.OR.MB$DEFBUTTON1)
	if(msg==MB$IDOK) then
		stop
	end if	
	
end subroutine   
   
SUBROUTINE DO_COPYFILE(SRCFILE,IDESFILE)
	IMPLICIT NONE
	CHARACTER(512),INTENT(IN)::SRCFILE
	INTEGER,INTENT(IN)::IDESFILE
	integer::ef,ISRC1,ITERM
	parameter(iterm=512)
	character(iterm)::term
	
	ISRC1=23
	OPEN(ISRC1,FILE=SRCFILE,STATUS='OLD')
	ef=0	
	do while(ef==0)
		read(ISRC1,'(A<ITERM>)',iostat=ef) term	
		WRITE(IDESFILE,'(A<ITERM>)') TERM
		TERM=''
		if(ef<0) exit
	ENDDO 
	CLOSE(23)

ENDSUBROUTINE   


    SUBROUTINE I_ENLARGE_AR(AVAL,DSTEP)
        INTEGER,ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        INTEGER,ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0,ALLSTAT
        
        
        ALLOCATE(VAL1(dstep))
        AVAL=[AVAL,VAL1]
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
    END SUBROUTINE

    SUBROUTINE R_ENLARGE_AR(AVAL,DSTEP)
        REAL(8),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        REAL(8),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
    
        ALLOCATE(VAL1(dstep))
        AVAL=[AVAL,VAL1]
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
    END SUBROUTINE
    SUBROUTINE NODE_ENLARGE_AR(AVAL,DSTEP)
        TYPE(point_tydef),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(point_tydef),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
    
        ALLOCATE(VAL1(dstep))
        AVAL=[AVAL,VAL1]
        DEALLOCATE(VAL1)
    END SUBROUTINE   	
    SUBROUTINE ELEMENT_ENLARGE_AR(AVAL,DSTEP)
        TYPE(element_tydef),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(element_tydef),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
    
        ALLOCATE(VAL1(dstep))
        AVAL=[AVAL,VAL1]
        DEALLOCATE(VAL1)
    END SUBROUTINE   
    SUBROUTINE EDGE_ENLARGE_AR(AVAL,DSTEP)
        TYPE(EDGE_tydef),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(EDGE_tydef),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
        
        
        !LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1(dstep))
        AVAL=[AVAL,VAL1]
        DEALLOCATE(VAL1)
    END SUBROUTINE 
    SUBROUTINE AR2D_ENLARGE_AR(AVAL,DSTEP)
        TYPE(ar2d_tydef),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(ar2d_tydef),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
        
        
        !LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1(dstep))
        AVAL=[AVAL,VAL1]
        DEALLOCATE(VAL1)
    END SUBROUTINE     
    SUBROUTINE AR2D_ENLARGE_AR2(AVAL,DSTEP)
        TYPE(ar2d_tydef2),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(ar2d_tydef2),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
        
        
        !LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1(dstep))
        AVAL=[AVAL,VAL1]
        DEALLOCATE(VAL1)
    END SUBROUTINE    
    
    SUBROUTINE ADJLIST_ENLARGE_AR(AVAL,DSTEP)
        TYPE(ADJLIST_tydef),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(ADJLIST_tydef),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
        
        
        !LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1(dstep))
        AVAL=[AVAL,VAL1]
        DEALLOCATE(VAL1)
    END SUBROUTINE    
    
	subroutine seg_initialize(this)
		use ds_t,only:arr_t
		!use meshds
		implicit none
        class(seg_type)::this
		integer::i,j,k,n1,n2,n3,sign1
        integer,allocatable::nnode1(:),nar1(:)
		real(8)::t1,t2,xi,yi,xj,yj,x1,y1        
		logical::tof1
		
        
        if(this.isini) return
        
        allocate(nnode1(1000))
        
		!do i=1,nseg
            
        xi=arr_t(this.sv).x;yi=arr_t(this.sv).y
        xj=arr_t(this.ev).x;yj=arr_t(this.ev).y
        n1=0;nnode1=0
        
        node.subbw=0!借用
        node.layer=0
        do j=1,ncedge
            if(cedge(j).cl/=this.icl) cycle
            x1=sum(node(cedge(j).v).x)/2.0
            y1=sum(node(cedge(j).v).y)/2.0
            call vins2d(xi,yi,xj,yj,x1,y1,tof1)
            if(tof1) then
                n1=n1+1
                if(n1>size(nnode1)) then
                    allocate(nar1(2*size(nnode1)))
                    nar1(1:size(nnode1))=nnode1
                    deallocate(nnode1)
                    allocate(nnode1,source=nar1)
                    deallocate(nar1)
                endif
                nnode1(n1)=cedge(j).edge
                    
                                
                    
                do k=1,2
                    n2=edge(cedge(j).edge).v(k)
                    if(k==1) then
                        sign1=1 !first vetex
                    else
                        sign1=-1 !second vetex
                    endif
                    if(node(n2).subbw==0) then
                        node(n2).subbw=j*sign1
                    else
                        node(n2).layer=j*sign1
                    endif
                enddo
            endif
		enddo
            
        if(n1>0) then
			this.nnum=count(node(1:).subbw/=0)
			allocate(this.node(this.nnum))
            this.node=pack([1:nnode],node(1:).subbw/=0)
            this.nedge=n1
			allocate(this.edge(n1))
            this.edge=nnode1(1:n1)
                
                
            !sorted
            nnode1=0
            !find the head
            do j=1,this.nnum
                n2=this.node(j)
                if(node(n2).subbw*node(n2).layer/=0) cycle
                !if(node(n2).subbw==
                if(abs(xi-node(n2).x)>1e-6) cycle
                if(abs(yi-node(n2).y)>1e-6) cycle
                nnode1(1)=n2
                exit
            enddo
                
            !find the successor
            do j=2,this.nnum
                n2=nnode1(j-1)
                if(j==2) then
                    if(node(n2).subbw/=0) then
                        n3=node(n2).subbw
                    else
                        n3=node(n2).layer
                    endif
                else
                    if(abs(node(n2).subbw)==n3) then
                        n3=node(n2).layer
                    else
                        n3=node(n2).subbw
                    endif
                endif
                    
                if(n3>0) then
                    nnode1(j)=edge(cedge(n3).edge).v(2)
                    this.edge(j-1)=cedge(n3).edge
                else
                    n3=abs(n3)
                    nnode1(j)=edge(cedge(n3).edge).v(1)
                    this.edge(j-1)=-cedge(n3).edge
                endif
                !this.edge(j-1)=cedge(n3).edge
            enddo
                
            this.node=nnode1(1:this.nnum)
               
        endif
            
		!end do
        this.isini=.true.
        
        if(allocated(nnode1)) deallocate(nnode1)
        
        node.subbw=0
        node.layer=0
	
	end subroutine    
    
    
    function seg_get_ordered_node(this,x) result(INode)
    !return ordered nodes starting headnode on this seg. 
        use ds_t,only:arr_t
        implicit none
        class(seg_type)::this
        real(8),intent(in),optional::x(2)
        integer,target,allocatable::inode(:)
        real(8)::t1
        
        if(.not.this.isini) call this.initialize
        
        allocate(inode(this.nnum))
        
        if(present(x)) then
            t1=((x(1)-arr_t(this.sv).x)**2+(x(2)-arr_t(this.sv).y)**2)**0.5
            if(abs(t1)<1e-7) then
                inode=this.node
            else
                inode=this.node(this.nnum:1:-1)
            endif
        else 
            inode=this.node
        endif 
    
    endfunction
    function seg_get_ordered_edge(this,x) result(INode)
    !return ordered edges starting headnode on this seg. 
        use ds_t,only:arr_t
        implicit none
        class(seg_type)::this
        real(8),intent(in),optional::x(2)
        integer,target,allocatable::inode(:)
        real(8)::t1
        
        
        if(.not.this.isini) call this.initialize
        
        allocate(inode(this.nedge))
        
        if(present(x)) then
            t1=((x(1)-arr_t(this.sv).x)**2+(x(2)-arr_t(this.sv).y)**2)**0.5
            if(abs(t1)<1e-7) then
                inode=this.edge
            else
                inode=-this.edge(this.nedge:1:-1)
            endif
        else 
            inode=this.edge
        endif 
    
    endfunction
    
    subroutine seg_get_parameters(this,icl,sv,ev,isini,nnum,nedge,isbg)
        implicit none
        class(seg_type)::this
        integer,intent(out),optional::icl,sv,ev,nnum,nedge,isbg
        logical,intent(out),optional::isini
        
        if(present(icl)) icl=this.icl
        if(present(sv)) sv=this.sv
        if(present(ev)) ev=this.ev
        if(present(isini)) isini=this.isini
        if(present(isbg)) isbg=this.isbg
        if(present(nnum)) then
            if(.not.this.isini) call this.initialize
            nnum=this.nnum
        endif
        if(present(nedge)) then
            if(.not.this.isini) call this.initialize
            nedge=this.nedge
        endif
    endsubroutine
    
    subroutine seg_set_parameters(this,icl,sv,ev,isbg)
        implicit none
        class(seg_type)::this
        integer,intent(in),optional::icl,sv,ev,isbg

        
        if(present(icl)) this.icl=icl
        if(present(sv)) this.sv=sv
        if(present(ev)) this.ev=ev 
        if(present(isbg)) this.isbg=isbg 
    endsubroutine   
    
    integer function v2edge(adjlist,v1,v2)
        implicit none
        type(adjlist_tydef)::adjlist(:)
        integer,intent(in)::v1,v2
        integer::n1,i
        
        n1=adjlist(v1).count
        v2edge=0
        do i=1,n1
            if(adjlist(v1).node(i)==v2) then
                v2edge=adjlist(v1).edge(i)
                exit
            endif
        enddo
    endfunction
    
    function iar2str(ia) result(str)
        implicit none
        integer,intent(in)::ia(:)
        character(len=:),allocatable::str        
        character(32)::ch1
        integer::I,n1
        
        n1=size(ia)        
        allocate(character(len=n1*17)::str)
        str=""
        if(n1==0) return
        do i=1,n1
            write(ch1,"(I,',')") ia(i)
            str=trim(adjustl(str))//trim(adjustl(ch1))
            if(mod(i,100)==0) str=trim(adjustl(str))//new_line('A')
        enddo
        
        N1=LEN_TRIM(STR)                    
        IF(STR(N1:N1)==NEW_LINE('A')) THEN
            STR(N1:N1)=""
            N1=N1-1
        ENDIF
		IF(STR(N1:N1)==',') THEN
            STR(N1:N1)=""
            N1=N1-1 !-1,get rid off ','
        ENDIF
            
        
        
    endfunction
    
     !update the adjlist(i).node and adjlist(i).edge,where i=p and v.
     subroutine addadjlist(adjlist,v1,v2,iedge,jedge)
	    !use meshds,only:adjlist_tydef
	    implicit none
        type(adjlist_tydef),allocatable::adjlist(:)
	    integer,intent(in)::v1,v2,iedge
        integer,optional,intent(out)::jedge !如果iedge已经存在，则以jedge返回其位置,否则jedge=iedge
        integer::i,j,v(2),n1,n2,n3,nsize1,jedge1
	    integer,pointer::node1(:)=>null(),edge1(:)=>null()
        
	    if(v1<0.or.v2<0) return !the edge in SuperTri is not taken into the list.
	    v(1)=v1
	    v(2)=v2
        if(any(v>size(adjlist))) call enlarge_ar(adjlist,100)
	    do i=1,2
		    if(.not.associated(adjlist(v(i)).node)) then
			    allocate(adjlist(v(i)).node(maxnadjlist),adjlist(v(i)).edge(maxnadjlist))
			    adjlist(v(i)).node=-1
			    adjlist(v(i)).edge=-1
		    end if
        end do
        !jedge1=0
	    if(any(adjlist(v1).node==v2)) then
            n2=minloc(abs(adjlist(v1).node-v2),dim=1)
            jedge1=adjlist(v1).edge(n2)
        else            
		    do i=1,2
			    n1=v(i)
			    n2=v(mod(i,2)+1)
			    adjlist(n1).count=adjlist(n1).count+1
			    n3=adjlist(n1).count
			    nsize1=size(adjlist(n1).node)
			    if(adjlist(n1).count>nsize1) then				
				    allocate(node1(2*nsize1), edge1(2*nsize1))
				    node1(1:nsize1)=adjlist(n1).node
				    edge1(1:nsize1)=adjlist(n1).edge
				    deallocate(adjlist(n1)%node,adjlist(n1)%edge)
				    !allocate(adjlist(n1).node(2*nsize1),adjlist(n1).edge(2*nsize1))
				    adjlist(n1).node=>node1
				    adjlist(n1).edge=>edge1
				    adjlist(n1).node(nsize1+1:2*nsize1)=-1
				    adjlist(n1).edge(nsize1+1:2*nsize1)=-1
				    nullify(node1,edge1)
			    end if
			    adjlist(n1).node(n3)=n2
			    adjlist(n1).edge(n3)=iedge
            end do
            jedge1=iedge
        end if
        if(present(jedge)) jedge=jedge1
     end subroutine

    subroutine Removeadjlist(adjlist,v1,v2)	
	    !use meshds,only:adjlist_tydef
        implicit none
        type(adjlist_tydef)::adjlist(:)
	    
	    integer::v1,v2,i,j,v(2),n1,n2,nc1
	
	    if(v1<0.or.v2<0) return !the edge in SuperTri is not taken into the list.
	    v(1)=v1
	    v(2)=v2
	    do i=1,2
		    n1=v(i)
		    n2=v(mod(i,2)+1)
		    nc1=adjlist(n1).count
		    do j=1,nc1
			    if(adjlist(n1).node(j)==n2) then
				    !move the last entry at the place j
				    adjlist(n1).node(j)=adjlist(n1).node(nc1)
				    adjlist(n1).node(nc1)=-1
				    adjlist(n1).edge(j)=adjlist(n1).edge(nc1)
				    adjlist(n1).edge(nc1)=-1
				    adjlist(n1).count=nc1-1
				    exit
			    end if
		    end do
	    end do
    end subroutine   
    
    function v2elts(inode) result(elts)
    !return all the 2D element sharing the node inode.
        implicit none
        integer,intent(in)::inode
        integer,allocatable,dimension(:)::elts
        integer::i,elt1(2)
        !借用
        elt(1:nelt).kcd=0
        
        do i=1,adjlist(inode).count
            elt1=edge(adjlist(inode).edge(i)).e
            where(elt1>0) elt(elt1).kcd=-1            
        enddo
        elts=pack([1:nelt],elt(1:nelt).kcd==-1)
        
    endfunction
    

    SUBROUTINE SETUP_SEARCH_ZONE_2D(SZ,NSZ,MAXX,MINX,MAXY,MINY,NODE1,TET1)
    !USE solverds,ONLY:ELEMENT_TYDEF
    IMPLICIT NONE
    TYPE(SEARCHZONE_TYDEF),ALLOCATABLE,INTENT(INOUT)::SZ(:)
    INTEGER,INTENT(INOUT)::NSZ
    REAL(8),INTENT(IN)::MAXX,MINX,MAXY,MINY
    TYPE(ELEMENT_TYDEF)::TET1(:)
    type(point_tydef),INTENT(IN)::NODE1(:)
    REAL(8)::DX,DY,DZ,BBOX1(2,2)
    INTEGER::NDX,NDY,NDZ,I,J,K,N1,NTET1
    
    NDX=5;NDY=5 !;NDZ=5
    DX=(MAXX-MINX)/NDX
    IF(ABS(DX)<1.D-7) NDX=1
    DY=(MAXY-MINY)/NDY
    IF(ABS(DY)<1.D-7) NDY=1
    !DZ=(MAXZ-MINZ)/NDZ
    !IF(ABS(DZ)<1.D-7) NDZ=1
    
    NSZ=NDX*NDY  !*NDZ
    
    IF(ALLOCATED(SZ)) DEALLOCATE(SZ)
    ALLOCATE(SZ(NSZ))
    
    NTET1=SIZE(TET1)

    N1=0
    !DO I=1,NDZ
        DO J=1,NDY
            DO K=1,NDX
                N1=N1+1
                SZ(N1).BBOX(1,1)=MINX+DX*(K-1)
                SZ(N1).BBOX(2,1)=MINX+DX*K
                SZ(N1).BBOX(1,2)=MINY+DY*(J-1)
                SZ(N1).BBOX(2,2)=MINY+DY*J 
                !SZ(N1).BBOX(1,3)=MINZ+DZ*(I-1)
                !SZ(N1).BBOX(2,3)=MINZ+DZ*I
                IF(.NOT.ALLOCATED(SZ(N1).ELEMENT)) ALLOCATE(SZ(N1).ELEMENT(MAX(NTET1/NSZ,100)))
                SZ(N1).ELEMENT=0
                !SZ(N1).NDX=[NDX,NDY,NDZ]
                SZ(N1).NDX(1:2)=[NDX,NDY]
            ENDDO
        ENDDO
    !ENDDO
    

    
    DO I=1,NTET1
        
        if(tet1(i).isdel.or.TET1(I).et/=0) CYCLE
        IF(ANY(TET1(i).NODE(1:3)<1)) CYCLE
		bbox1(1,1)=MINVAL(NODE1(TET1(i).NODE(1:3)).X)
        bbox1(2,1)=MAXVAL(NODE1(TET1(i).NODE(1:3)).X)
		bbox1(1,2)=MINVAL(NODE1(TET1(i).NODE(1:3)).Y)
        bbox1(2,2)=MAXVAL(NODE1(TET1(i).NODE(1:3)).Y)

        tet1(i).bbox=bbox1
        
        
        DO J=1,NSZ
            
            IF(NDX>1) THEN
                IF(SZ(J).BBOX(1,1)>BBOX1(2,1)) CYCLE !MIN>MAX
                IF(SZ(J).BBOX(2,1)<BBOX1(1,1)) CYCLE !MAX<MIN
            ENDIF
            IF(NDY>1) THEN
                IF(SZ(J).BBOX(1,2)>BBOX1(2,2)) CYCLE !MIN>MAX
                IF(SZ(J).BBOX(2,2)<BBOX1(1,2)) CYCLE !MAX<MIN 
            ENDIF
            !IF(SZ(j).DIM>2.AND.NDZ>1) THEN
            !    IF(SZ(J).BBOX(1,3)>BBOX1(2,3)) CYCLE !MIN>MAX
            !    IF(SZ(J).BBOX(2,3)<BBOX1(1,3)) CYCLE !MAX<MIN 
            !ENDIF
            SZ(J).NEL=SZ(J).NEL+1
            IF(SIZE(SZ(J).ELEMENT)<SZ(J).NEL) THEN
                CALL I_ENLARGE_AR(SZ(J).ELEMENT,100)
            ENDIF
            SZ(J).ELEMENT(SZ(J).NEL)=I
        ENDDO
    ENDDO
    
    !isszinit=.true.
    
    ENDSUBROUTINE
    
    
!search method:element-by-element，it is robust but may be comparily slow.
!if tryiel>0, it the element and its nerghborhood will be searched firstly.
INTEGER FUNCTION POINTlOC_2D(PT,TRYIEL,SZ,NSCH,ESCH)

	IMPLICIT NONE
	REAL(8),INTENT(IN)::PT(2)
    INTEGER,INTENT(IN)::TRYIEL
    TYPE(SEARCHZONE_TYDEF),INTENT(IN)::SZ(:)
    type(point_tydef),INTENT(IN)::NSCH(:)
    type(element_tydef),INTENT(IN)::ESCH(:)
	INTEGER::I,J,K,ELT1(-10:0),NSZ
	!LOGICAL,EXTERNAL::PtInTri,PtInTET
    LOGICAL::ISFOUND=.FALSE.
    REAL(8)::TRI1(2,3)
	!INTEGER::SEARCHLIST(NTET)
	!I=INT(1+(RANDOM(1))(NFACE-1))
    
    ISFOUND=.FALSE.    
    POINTlOC_2D=0
    NSZ=SIZE(SZ)
    DO J=1,NSZ
        IF(SZ(J).NEL<1) CYCLE
        IF(SZ(J).NDX(1)>1) THEN
            IF(PT(1)<SZ(J).BBOX(1,1)) CYCLE
            IF(PT(1)>SZ(J).BBOX(2,1)) CYCLE
        ENDIF
        IF(SZ(J).NDX(2)>1) THEN
            IF(PT(2)<SZ(J).BBOX(1,2)) CYCLE
            IF(PT(2)>SZ(J).BBOX(2,2)) CYCLE
        ENDIF        
  
        
        IF(TRYIEL>0) THEN
            ELT1(-ESCH(TRYIEL).NNUM)=TRYIEL
            DO K=1,ESCH(TRYIEL).NNUM
                ELT1(-ESCH(TRYIEL).NNUM+K)=ESCH(TRYIEL).ADJ(K)
            ENDDO
            K=-ESCH(TRYIEL).NNUM
        ELSE
            K=1
        ENDIF
        
	    do K=K,SZ(J).NEL
            IF(K>0) THEN
                I=SZ(J).ELEMENT(K) 
            ELSE
                I=ELT1(K)
            ENDIF
            IF(I<1) CYCLE
		    IF(ESCH(I).ISDEL) CYCLE
            IF(ANY(ESCH(I).node(1:3)<1)) CYCLE 
		    IF(ESCH(I).BBOX(2,1)>=PT(1).and.ESCH(I).BBOX(1,1)<=PT(1)) then
			    IF(ESCH(I).BBOX(2,2)>=PT(2).and.ESCH(I).BBOX(1,2)<=PT(2)) then
				    !IF(POSDATA.NDIM>2) THEN
					   ! IF(ESCH(I).BBOX(2,3)<PT(3).and.ESCH(I).BBOX(1,3)>PT(3)) CYCLE
				    !ENDIF
                    TRI1(1,:)=NSCH(ESCH(I).NODE(1:3)).X
                    TRI1(2,:)=NSCH(ESCH(I).NODE(1:3)).Y
                    CALL triangle_contains_point_2d_2 ( TRI1, PT, ISFOUND )
					IF(ISFOUND) THEN
                        POINTlOC_2D=I
                        EXIT
                    ENDIF
			    ENDIF
		    endif
        ENDDO
        
        IF(ISFOUND) EXIT

    ENDDO
	   

	RETURN
ENDFUNCTION    
    
subroutine triangle_contains_point_2d_2 ( t, p, inside )

!*****************************************************************************80
!
!! TRIANGLE_CONTAINS_POINT_2D_2 finds if a point is inside a triangle in 2D.
!
!  Discussion:
!
!    The routine assumes that the vertices are given in counter clockwise
!    order.  If the triangle vertices are actually given in clockwise 
!    order, this routine will behave as though the triangle contains
!    no points whatsoever!
!
!    The routine determines if a point P is "to the right of" each of the lines
!    that bound the triangle.  It does this by computing the cross product
!    of vectors from a vertex to its next vertex, and to P.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!    The vertices should be given in counter clockwise order.
!
!    Input, real ( kind = 8 ) P(2), the point to be checked.
!
!    Output, logical ( kind = 4 ) INSIDE, is TRUE if the point is 
!    inside the triangle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  logical ( kind = 4 ) inside
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) t(dim_num,3)

  do j = 1, 3

    k = mod ( j, 3 ) + 1

    if ( 1.0D-7 < ( p(1) - t(1,j) ) * ( t(2,k) - t(2,j) ) &
                 - ( p(2) - t(2,j) ) * ( t(1,k) - t(1,j) ) ) then
      inside = .false.
      return
    end if

  end do

  inside = .true.

  return
end

function Find_Point_Inside_Polygon_2D(Poly) result(Pt)
    
    implicit none
    real(8),intent(in)::poly(:,:)
    real(8)::Pt(2)
    real(8)::x1,xmin1,xmax1,dx1
    real(8),allocatable::y1(:)
    integer::i,n1,n2
    
    n1=size(poly,dim=2)
    xmin1=minval(poly(1,:))
    xmax1=maxval(poly(1,:))
    dx1=(xmax1-xmin1)/200.
    x1=(xmax1+xmin1)/2.0
    do while(any(abs(poly(1,:)-x1) <1e-7))
        x1=x1+dx1
        if(x1>=xmax1) then
            print *, 'failed in getting a trail x in function Find_Point_Inside_Polygon_2D'
            stop
        endif
    enddo
        
    do i=1,n1
        n2=mod(i,n1)+1
        if((x1-poly(1,i))*(x1-poly(1,n2))<0) then
            y1=[y1,(poly(2,i)-poly(2,n2))/(poly(1,i)-poly(1,n2))*(x1-poly(1,i))+poly(2,i)]
        endif
    enddo

    call quick_sort(y1)
    pt=[x1,(y1(1)+y1(2))/2]
  
end function

function Find_Point_Inside_Polygon_2D_edge(bedge) result(Pt)
    
    implicit none
    integer,intent(in)::bedge(:)
    real(8)::Pt(2)
    real(8)::x1,xmin1,xmax1,dx1
    real(8),allocatable::y1(:),xt1(:)
    integer::i,n1,n2,v1(2)
    
    n1=size(bedge,dim=1)
	xt1=node(edge(bedge).v(1)).x !close boundaries, edge.v(2)与edge.v(1)重合
    xmin1=minval(xt1) 
    xmax1=maxval(xt1)
    dx1=(xmax1-xmin1)/200.
    x1=(xmax1+xmin1)/2.0
    do while(any(abs(xt1-x1) <1e-7))
        x1=x1+dx1
        if(x1>=xmax1) then
            print *, 'failed in getting a trail x in function Find_Point_Inside_Polygon_2D_edge'
            stop
        endif
    enddo
        
    do i=1,n1
        v1=edge(bedge(i)).v
        if((x1-node(v1(1)).x)*(x1-node(v1(2)).x)<0.) then
            y1=[y1,(node(v1(2)).y-node(v1(1)).y)/(node(v1(2)).x-node(v1(1)).x) &
		*(x1-node(v1(1)).x)+node(v1(1)).y]
        endif
    enddo

    call quick_sort(y1)
    pt=[x1,(y1(1)+y1(2))/2]
  
end function

function elt_bounded_by_edges(PT,ibedge) result(ielt)
!给出区域内的一点pt和区域边界的单元边edge(ibedge)
!返回区域内的所有单元
!区域可以是多连通区域，只要iedge包含所有边界边
    implicit none
    real(8),intent(in)::PT(2)
    integer,intent(in)::ibedge(:)
    integer,allocatable::ielt(:)
    integer::J,IELT1,IELT2,STACK1(10000),N1,N2,EDGE1(NEDGE),EDGE2(NEDGE),ELT1(1:NELT)           

    IELT1=POINTlOC_2D(PT,-1,SEARCHZONE,NODE(1:NNODE),ELT)
    IF(IELT1<1) THEN
        PRINT *, 'NO TRIANGLE ELEMENT IN THE REGION.'
    ENDIF
    ELT1=0
    STACK1(1)=IELT1
    ELT1(IELT1)=1
    N2=1
                
    edge2=0
    edge2(Ibedge)=1
                
    DO WHILE(N2>0)
        IELT1=STACK1(N2) 
			
        N2=N2-1
        DO J=1,ELT(IELT1).NNUM
            N1=ELT(IELT1).EDGE(J)
            IF(EDGE2(N1)==0.AND.EDGE1(N1)/=1) THEN
                EDGE1(N1)=1                            
                IELT2=ELT(IELT1).ADJ(J)
                N2=N2+1
                STACK1(N2)=IELT2
                ELT1(IELT2)=1
            ENDIF
        ENDDO
    ENDDO
    IELT=PACK([1:NELT],ELT1(1:NELT)==1)
end function
subroutine triangle_contains_point_2d_3 ( t, p, inside )

!*****************************************************************************80
!
!! TRIANGLE_CONTAINS_POINT_2D_3 finds if a point is inside a triangle in 2D.
!
!  Discussion:
!
!    This routine is the same as TRIANGLE_CONTAINS_POINT_2D_2, except
!    that it does not assume an ordering of the points.  It should
!    work correctly whether the vertices of the triangle are listed
!    in clockwise or counter clockwise order.
!
!    The routine determines if a point P is "to the right of" each of the lines
!    that bound the triangle.  It does this by computing the cross product
!    of vectors from a vertex to its next vertex, and to P.
!
!    The point is inside the triangle if it is to the right of all
!    the lines, or to the left of all the lines.
!
!    This version was suggested by Paulo Ernesto of Maptek Brasil.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!
!    Input, real ( kind = 8 ) P(2), the point to be checked.
!
!    Output, logical ( kind = 4 ) INSIDE, is TRUE if the point is 
!    inside the triangle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) dir_new
  real ( kind = 8 ) dir_old
  logical ( kind = 4 ) inside
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) t(dim_num,3)

  dir_old = 0.0D+00

  do j = 1, 3

    k = mod ( j, 3 ) + 1

    dir_new = ( p(1) - t(1,j) ) * ( t(2,k) - t(2,j) ) &
            - ( p(2) - t(2,j) ) * ( t(1,k) - t(1,j) )
    
    if ( dir_new * dir_old < 0.0D+00 ) then
      inside = .false.
      return
    end if

    if ( dir_new /= 0.0D+00 ) then
      dir_old = dir_new
    end if

  end do

  inside = .true.

  return
end

subroutine plane_imp_line_par_int_3d ( a, b, c, d, x0, y0, z0, f, g, h, &
  intersect, p ,t)

!*****************************************************************************80
!
!! PLANE_IMP_LINE_PAR_INT_3D: intersection ( impl plane, param line ) in 3D.
!
!  Discussion:
!
!    The implicit form of a plane in 3D is:
!
!      A * X + B * Y + C * Z + D = 0
!
!    The parametric form of a line in 3D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!      Z = Z0 + H * T
!
!    We normalize by always choosing F*F + G*G + H*H = 1, 
!    and F nonnegative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 January 2005
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
!    ISBN: 0408012420,
!    page 111.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, C, D, the implicit plane parameters.
!
!    Input, real ( kind = 8 ) X0, Y0, Z0, F, G, H, parameters that define the
!    parametric line.
!
!    Output, integer ( kind = 4 ) INTERSECT, the kind of intersection;
!    0, the line and plane seem to be parallel and separate;
!    1, the line and plane intersect at a single point;
!    2, the line and plane seem to be parallel and joined.
!
!    Output, real ( kind = 8 ) P(3), is a point of intersection of the line
!    and the plane, if INTERSECT is TRUE.
!     Output ,real t,intersection location, 0<=t<=1.0 the  intersect poin inside the edge.(add and modified by lgy )
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) denom
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) h
  integer ( kind = 4 ) intersect
  real ( kind = 8 ) norm1
  real ( kind = 8 ) norm2
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) t
  real ( kind = 8 ), parameter :: tol = 0.00001D+00
  real ( kind = 8 ) x0
  real ( kind = 8 ) y0
  real ( kind = 8 ) z0
  
!
!  Check.
!
  norm1 = sqrt ( a * a + b * b + c * c )

  if ( norm1 == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP_LINE_PAR_INT_3D - Fatal error!'
    write ( *, '(a)' ) '  The plane normal vector is null.'
    stop 1
  end if

  norm2 = sqrt ( f * f + g * g + h * h )

  if ( norm2 == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PLANE_IMP_LINE_PAR_INT_3D - Fatal error!'
    write ( *, '(a)' ) '  The line direction vector is null.'
    stop 1
  end if

  denom = a * f + b * g + c * h
!
!  The line and the plane may be parallel.
!
!  if ( abs ( denom ) < tol * norm1 * norm2 ) then
  if ( abs ( denom ) < 1.d-8 ) then !modified by lgy 
    if ( a * x0 + b * y0 + c * z0 + d == 0.0D+00 ) then
      intersect = 2
      p(1) = x0
      p(2) = y0
      p(3) = z0
      t=0.d0
    else
      intersect = 0
      p(1:dim_num) = 0.0D+00
    end if
!
!  If they are not parallel, they must intersect.
!
  else

    intersect = 1
    t = - ( a * x0 + b * y0 + c * z0 + d ) / denom
    p(1) = x0 + t * f
    p(2) = y0 + t * g
    p(3) = z0 + t * h

  end if

  return
  end

subroutine plane_exp2imp_3d ( p1, p2, p3, a, b, c, d )

!*****************************************************************************80
!
!! PLANE_EXP2IMP_3D converts an explicit plane to implicit form in 3D.
!
!  Discussion:
!
!    The explicit form of a plane in 3D is
!
!      the plane through P1, P2 and P3.
!
!    The implicit form of a plane in 3D is
!
!      A * X + B * Y + C * Z + D = 0
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2005
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
!    Input, real ( kind = 8 ) P1(3), P2(3), P3(3), three points on the plane.
!
!    Output, real ( kind = 8 ) A, B, C, D, coefficients which describe 
!    the plane.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) p3(dim_num)

  a = ( p2(2) - p1(2) ) * ( p3(3) - p1(3) ) &
    - ( p2(3) - p1(3) ) * ( p3(2) - p1(2) )

  b = ( p2(3) - p1(3) ) * ( p3(1) - p1(1) ) &
    - ( p2(1) - p1(1) ) * ( p3(3) - p1(3) )

  c = ( p2(1) - p1(1) ) * ( p3(2) - p1(2) ) &
    - ( p2(2) - p1(2) ) * ( p3(1) - p1(1) )

  d = - p2(1) * a - p2(2) * b - p2(3) * c

  return
end  
  
end module
