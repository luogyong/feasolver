module meshDS
    
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
                    Zorder=0 !=0����̴߳�����������,=1,��֮��
    integer,parameter:: ET_PRM=63,&
                    ET_TET=43,Linear=1,Membrance=0
    logical::ismeminpdone=.false.
    real(8)::um
    
	parameter(maxstk=2000,um=1e15)
    
    INTERFACE ENLARGE_AR
        MODULE PROCEDURE I_ENLARGE_AR,R_ENLARGE_AR,NODE_ENLARGE_AR,ELEMENT_ENLARGE_AR,&
                         edge_enlarge_ar
    END INTERFACE 

	type point_tydef
     integer::number,bw !��ı��,semibandwidth������,�������ӱ�š�
	   real(8) x,y,z !�������
	   real(8)::s  !�õ㸽���ĵ�Ԫ��С
	   INTEGER::SUBBW  !=-999 ,the node is dead.
	   integer::layer=0,onbdy=0 !nodal layer
       integer::havesoildata=0 ![0,1,2] 0,no;1,partially;2,completely yes.(all soillayers have elevation(no -999)) 
       !
       real(8),allocatable::elevation(:)
	end type
	type(point_tydef),allocatable,target::node(:)
	type(point_tydef),allocatable::cp(:)  !control point
	integer::nnode=0,tnode=0 !����ڵ������ܽڵ���

	                            										
	type BP_tydef !bounded points
	   type(point_tydef),pointer::npt=>null()
	   type(BP_tydef),pointer::next=>null()  !previousָ��ǰһ����Ԫ��nextָ���һ����Ԫ
	end type

	type constrainline_tydef
	   real(8),pointer::conpoint(:,:)=>null() !�������ϸ�������꣨x,y��
	   integer::flag=0,hole=0  !�Ƿ�պ�,�Ƿ������� 0Ϊ��1Ϊ��
	   integer::num,material=0  !�����ߵ����Ŀ,materialֻ�Է���ǽ����Ϊ��Ϊ����ǽ�Ĳ��Ϻ�,
       integer,allocatable::point(:)
	end type
	

	type element_tydef
	   !type(point_tydef),pointer::xy1,xy2,xy3,xy4   !��Ԫ���������
	   logical::isdel=.false.
       INTEGER::NNUM=3
	   integer::node(15) !6-noded/15-noded triangleelement node in node().
	   integer::et=0           !element type, 3-noded:0, 6-noded triangle:6; 15-noded:15 
	                                      !6-noded prism element: 63; 15-noded prism element:153;
	                                      !4-noded tetrahedral element:43, second order tetrahedral element: 103
                                !21,line,element
                                !42,quadrangle element
	   integer::zn=1,number,ZN2=-1,mat=-1
	   integer::edge(3)=0,ORIENT(3)=1 !the edge number in edge set.
	   real(8)::property(5)=0    !��Ԫ����ֵ
	   integer::kcd=0,Maxedge=0  !��Ԫ�Ŀɷֶ�,
	   integer::adj(4)=-1 !if =-1,no adjacent element
       INTEGER::MOTHER=0,iLAYER=0,ISMODEL=1
	end type
!	type(element_tydef),pointer::Ehead,Ept,element,Etail!ǰ�ĸ�ָ�򻮷ֺõĵ�Ԫ��������ָ����������������ε�Ԫ
	type(element_tydef),target,allocatable::ELT(:)
	integer::nelt=0,Pehead=-1,PEpt=-1,IEpt=-1,ept=-1
	
	!type PrismElement_tydef !������Ԫ
	!	logical::isdel=.false.
	!    integer::nnum
	!    integer,allocatable::node(:)
	!    integer::mother=0 !generated from the elt(mother)
	!    integer::et=63   !6-noded prism element: 63; 15-noded prism element:153;
	!                                 !4-noded tetrahedral element:43, second order tetrahedral element: 103
	!    integer::mat !material 
	!    integer::nlayer(2)=0 !nlayer(1��2) �ֱ�Ϊ�����͵ײ��ڵ���
 !               
	!end type 
	!type(PrismElement_tydef),allocatable::PRMELT(:),TetELT(:)
	!integer::NPRMELT=0,NTetelt=0
	
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
	

	type modelgroup_type
        integer::izone=1,ilayer=0,mat=1,coupleset=-1,sf=0 
        !����2D���浥Ԫ,ilayer=i,��ʾ��i���߳��档
        !����2D���浥Ԫ,ilayer=i,��ʾ��i�㵥Ԫ��ע�⣬n�㵥Ԫ����n+1�߳��棬��i�㵥Ԫ��i-1��i��߳������(������Ԫ)��
        character(32)::et=''
        character(32)::name=''
    endtype
    type(modelgroup_type),allocatable::model(:)
    integer::nmgroup=0

	type size_point_tydef  !��Ԫ�ߴ���Ƶ�
	    real(8)::x,y,a,d !���꣬�Ȳ�ĵ�һ��Ȳ�
	end type

	type zone_tydef
	   integer::num,k,OutGmshType=1,iElevation=1 !OutGmshType=1 Physical Volume only; =2 Physical Surface only; =3, both;if OutGmshType=2/3, iElevationָ������ĸ��̵߳��档       
	   real(8),allocatable::point(:,:)
	   real(8)::xmin=1e15,ymin=1e15,xmax=-1e15,ymax=-1e15
	   integer,allocatable::item(:),trie3n(:),Dise4n(:), trie6n(:),trie15n(:)
	   integer::nitem=0,ntrie3n=0,ndise4n=0,ntrie6n=0,ntrie15n=0  !restriction:nitem=ntrie3n+ndise4n
	   integer,allocatable::bedge(:),mat(:) !boundary edges id(point to edge),mat=ÿ�����Ĳ��Ϻ�
	   integer::nbe=0 !number of bedge
       !������PRM��TETʱ�����¸������õ�
       
       integer,allocatable::NPrm(:),NTet(:),PRM(:,:),TET(:,:),ISMODEL(:) 
       !����ÿ������������Ԫ�������嵥Ԫ�ĸ���,PRM,ÿ��������Ԫ��PRM_ELT�е��±�,TET���ơ�
       !ISmodel=1,yes output in element form. 
       character(32),allocatable::solver_et(:),NAME(:)
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
        INTEGER::ISZONEBC=-1,ISCEDGE=0
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
		integer,allocatable::nbc(:),niseg(:) !nbc(nvb)
		real(8),allocatable::elevation(:,:) !elevation(nnum,0:soillayer) ��������icl�ϸ���ĸ̡߳�
		real(8),allocatable::lincof(:,:,:) !lincof(2,n1,0:soillyger) !if csl(icl).flag==1,n1=nnum,==0,n1=nnum-1
		real(8),allocatable::vbc(:,:) !vbc(nvb,0:soillyger)
	end type
	type(meminp_tydef),allocatable::meminp(:),meminp2(:)
	integer::nmeminp=0,nmeminp2=0

	integer,allocatable::segindex(:,:) 
	!�洢ģ���߶εģ����ǵ�Ԫ�ߣ���seg�е�λ�ã���������ccl��ʽ����Ŀ����ߣ���Ϊ����������ϵĵ㲻�������������棩��
	type seg_type
        private
		integer::icl !���߶���cpphead(icl)��λ��
		integer::sv=0,ev=0 !�ڵ������������е�λ�á�
        logical::isini=.false.
		!integer::isA2Z=1 !�Ƿ���˳��1Ϊ˳�򣬼���cpphead(icl)������Ϊsv-ev. 0Ϊ����ev-sv. 
		!integer::isT2H=0 !�Ƿ���β�������ĶΡ�=1,yes.
		integer::nnum=0,nedge=0 !ϸ�ֺ���߶εĽڵ������
		integer,allocatable::node(:),edge(:) !ϸ�ֺ���߶εĽڵ���tnode���±�,�����˽ڵ㡣
		!type(BP_tydef),pointer::svp,evp !�˽ڵ���cpphead(icl)�е�λ��
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
	type(zone_tydef),TARGET,allocatable::zone(:)
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
    
    character(512)::COPYFILE(100)
	INTEGER::NCOPYFILE=0
contains
   !���ַ������൱�������ַ�(����������)ת��Ϊ��Ӧ������
   !�� '123'תΪ123,'14-10'תΪ14,13,12,11,10
   !string��ת���������������ar(n1)���أ�����,n1Ϊ�ַ��������ֵĸ���:(ע��1-3ת����Ϊ3�����֣�1,2,3)
   !nmaxΪ����ar�Ĵ�С,stringĬ���ַ�����Ϊ256��
   !num_readΪҪ�������ݵĸ�����
   !unitΪ�ļ���
   !ÿ��ֻ����һ����Ч�У�����'/'��ͷ���У�
   !ÿ�к�����'/'��ʼ�ĺ�����ַ�����Ч�ġ�
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


		 if(string(1:1)/='/') then
			
			!ÿ�к�����'/'��ʼ�ĺ�����ַ�����Ч�ġ�
			if(index(string,'/')/=0) then
				strL=index(string,'/')-1
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
				  if(tof1) then !����������'1-5'��������ʽ�Ķ�������
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
					 if(tof2) then !����������'1*5'(��ʾ5��1)��������ʽ�Ķ�������
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
        INTEGER::LB1=0,UB1=0
        
        LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1,SOURCE=AVAL)
        DEALLOCATE(AVAL)
        ALLOCATE(AVAL(LB1:UB1+DSTEP))
        AVAL(LB1:UB1)=VAL1
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
    END SUBROUTINE

    SUBROUTINE R_ENLARGE_AR(AVAL,DSTEP)
        REAL(8),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        REAL(8),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
    
        LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1,SOURCE=AVAL)
        DEALLOCATE(AVAL)
        ALLOCATE(AVAL(LB1:UB1+DSTEP))
        AVAL(LB1:UB1)=VAL1
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
    END SUBROUTINE
    SUBROUTINE NODE_ENLARGE_AR(AVAL,DSTEP)
        TYPE(point_tydef),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(point_tydef),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
    
        LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1,SOURCE=AVAL)
        DEALLOCATE(AVAL)
        ALLOCATE(AVAL(LB1:UB1+DSTEP))
        AVAL(LB1:UB1)=VAL1
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
    END SUBROUTINE   	
    SUBROUTINE ELEMENT_ENLARGE_AR(AVAL,DSTEP)
        TYPE(element_tydef),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(element_tydef),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
    
        LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1,SOURCE=AVAL)
        DEALLOCATE(AVAL)
        ALLOCATE(AVAL(LB1:UB1+DSTEP))
        AVAL(LB1:UB1)=VAL1
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
    END SUBROUTINE   
    SUBROUTINE EDGE_ENLARGE_AR(AVAL,DSTEP)
        TYPE(EDGE_tydef),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(EDGE_tydef),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
    
        LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1,SOURCE=AVAL)
        DEALLOCATE(AVAL)
        ALLOCATE(AVAL(LB1:UB1+DSTEP))
        AVAL(LB1:UB1)=VAL1
        !AVAL(UB1+1:UB1+10)=0
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
        
        node.subbw=0!����
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
                    n2=cedge(j).v(k)
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
                if(abs(xi-node(n2).x)>1e-7) cycle
                if(abs(yi-node(n2).y)>1e-7) cycle
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
                    nnode1(j)=cedge(n3).v(2)
                else
                    n3=abs(n3)
                    nnode1(j)=cedge(n3).v(1)
                endif
                this.edge(j-1)=cedge(n3).edge
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
                inode=this.edge(this.nedge:1:-1)
            endif
        else 
            inode=this.edge
        endif 
    
    endfunction
    
    subroutine seg_get_parameters(this,icl,sv,ev,isini,nnum,nedge)
        implicit none
        class(seg_type)::this
        integer,intent(out),optional::icl,sv,ev,nnum,nedge
        logical,intent(out),optional::isini
        
        if(present(icl)) icl=this.icl
        if(present(sv)) sv=this.sv
        if(present(ev)) ev=this.ev
        if(present(isini)) isini=this.isini
        if(present(nnum)) then
            if(.not.this.isini) call this.initialize
            nnum=this.nnum
        endif
        if(present(nedge)) then
            if(.not.this.isini) call this.initialize
            nedge=this.nedge
        endif
    endsubroutine
    
    subroutine seg_set_parameters(this,icl,sv,ev)
        implicit none
        class(seg_type)::this
        integer,intent(in),optional::icl,sv,ev

        
        if(present(icl)) this.icl=icl
        if(present(sv)) this.sv=sv
        if(present(ev)) this.ev=ev 
        
    endsubroutine    
end module
