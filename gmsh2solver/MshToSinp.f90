module DS_Gmsh2Solver

	type node_type
		integer::inode=0
		real(8)::xy(3)=0
		INTEGER::CHILD=0
	end type
	type(node_type),allocatable::node(:)
	integer::nnode=0,TNNODE=0
	integer,allocatable::Noutputorder(:)

	type element_type
		integer::et=-1
		integer::ntag=0
		integer,allocatable::tag(:)
		!number-of-tags
		!gives the number of integer tags that follow for the n-th element. 
		!By default, the first tag is the number of the physical entity to which the element belongs;
		!the second is the number of the elementary geometrical entity to which the element belongs; 
		!the third is the number of mesh partitions to which the element belongs, 
		!followed by the partition ids (negative partition ids indicate ghost cells). A zero tag is equivalent to no tag. 

		integer::nnode=0,ndege=0,nface=0
		integer,allocatable::node(:),edge(:),face(:)
    CONTAINS
        PROCEDURE::GET_EDGE=>SET_ELEMENT_EDGE
        PROCEDURE::GET_FACE=>SET_ELEMENT_FACE
	end type
	type(element_type),allocatable::element(:)
	integer::nel=0
	
	
	type physicalGroup_type
		integer::ndim
		logical::ismodel=.FALSE. !is a part of the model
		integer::ET_GMSH=0 
		character(32)::name
		integer::nel=0
		integer,allocatable::element(:)
		integer::mat=0
		character(32)::ET="ELT_BC_OR_LOAD"
	end type
	integer,parameter::maxphgp=100
	type(physicalGroup_type)::physicalgroup(maxphgp)
	integer,allocatable::phgpnum(:)
	integer::nphgp=0
	
	type et_type
		integer::nnode=0,nedge=0,nface=0		
		character(512)::description
		integer,allocatable::edge(:,:),face(:,:)
		!edge(2,nedge),face(0:4,nface),face(0,:)==3,triangular face,==4, quadrilateral face
		real(8),allocatable::weight(:,:) !分布单元荷载各节点分布系数, weight(:,1) 平面均布荷载；weight(:,2) 平面三角形荷载;weight(:,3) 轴对称均布荷载；weight(:,4) 轴对称三角形荷载，
													!目前只能处理平面应变的两种情况。
	end type
	type(et_type)::elttype(100)
	
	type bcgroup_type
		integer::group
		integer::ndim=0 !=0,point load; =1,line load; =2, planar load; =3,volume load;
		integer::dof
		integer::sf=0
		real(8)::value=0  !当输入seepageface时，value=1,2,3 分别表示节点的水头值等于坐标x,y,z.
		integer::n1=0,n2=0 !for spgface output
	end type
	type(bcgroup_type),allocatable::elt_bc(:),elt_load(:),elt_spgface(:)
	integer::nelt_bc,nelt_load,nelt_spgface
	
	type bc_tydef
		integer::node,dof	!when body force is input, node equels to element	
		integer::sf=0 !step function for this load
		integer::isdead=0 !for seepage face iterative, =1,the condition is no longer on work. 
		real(8)::value=0
	end type
	type(bc_tydef),allocatable::nodalLoad(:),nodalBC(:),spgface(:)
	integer::nnodalLoad=0,maxnnodalLoad=1000
	integer::nnodalBC=0,maxnnodalBC=1000
	integer::nspgface=0,maxnspgface=1000
    
    type wsp_typdef
        integer::Group,chgroup,spgroup !chgroup, characteristic point.
        integer::xdirection=1
        integer::nnode=0,nchnode=0
        integer,allocatable::node(:)
		integer,allocatable::chnode(:)
    end type
    type(wsp_typdef),allocatable::wsp(:)
    integer::nwsp=0
    
    type out_data_typdef
        integer::Group
        integer::SPGroup !Startpoint Group,这个Group只有一个点。
		integer::order=0
        integer::nnode=0
        integer,allocatable::node(:)        
    end type
    type(out_data_typdef),allocatable::DataPoint(:)
    integer::NDataPoint=0
	
	
	integer,allocatable::modelgroup(:)
	integer::nmodelgroup=0
		
	character(512)::resultfile,title,INCLUDEFILE(0:100),COPYFILE(100)
	INTEGER::NINCLUDEFILE=0,NCOPYFILE=0
		
	integer,allocatable::adjL(:,:)
	integer::maxadj=50
	
	integer::modeldimension=2
	
    CONTAINS
    

	
    SUBROUTINE SET_ELEMENT_EDGE(THIS)
        CLASS(element_type)::THIS
		INTEGER::I
		
		DO I=1,elttype(THIS.ET).NEDGE
			
		ENDDO
        
    END SUBROUTINE
    SUBROUTINE SET_ELEMENT_FACE(THIS)
        CLASS(element_type)::THIS
        
    END SUBROUTINE    
end module


program GmshToSolver
	use ifcore
	use ifqwin
	use DS_Gmsh2Solver
	implicit none
	integer::i
	character(1)::key
	type(qwinfo) winfo
	LOGICAL(4)::tof,pressed
	integer,allocatable::IPERM(:)
	
	winfo%TYPE = QWIN$MAX
	tof=SETWSIZEQQ(QWIN$FRAMEWINDOW, winfo)
	tof=SETWSIZEQQ(0, winfo)  
	
	print *, 'LGY Works. Msh2Sinp'
	
	write(*, 10) 
	key=getcharqq()
	if(ichar(key)==ichar('h').or.ichar(key)==ichar('H')) then		
		call write_readme_gmsh2sinp()			
		stop
	end if

	
	print *, 'Begin to read in data...'
	call readin()
	
	!统计每个物理组单元的个数及单元数组下标。
	do i=1,maxphgp
		if(physicalgroup(i).nel>0) allocate(physicalgroup(i).element(physicalgroup(i).nel))
	end do
	
	physicalgroup.nel=0
	do i=1,nel
		physicalgroup(element(i).tag(1)).nel=physicalgroup(element(i).tag(1)).nel+1
		physicalgroup(element(i).tag(1)).element(physicalgroup(element(i).tag(1)).nel)=i		
	end do	
	
	!CALL GEN_ZEROTHICKNESS_ELEMENT()
	
	! reordering nodal number
	allocate(Noutputorder(nnode))
	ALLOCATE(IPERM(nnode))
	call setup_adjList()
	call reorder_nodal_number(IPERM,nnode,adjL,maxadj)
	do i=1,nnode
		noutputorder(IPERM(i))=i	
		node(i).inode=IPERM(i)		
	end do
	DEALLOCATE(IPERM)
	

	
	CALL elt_bc_load_translate()
	
	call Tosolver()
	
	Print *, 'GMSHTOSOLVER IS DONE. PLEASE ADD OTHER PARAMETERS BY HAND.'
10 format("Press 'H' to write a keyword help file named 'D:\README_GMSH2SINP.TXT'. Any other key to read an Msh file.")		
end program

SUBROUTINE SETUP_EDGE_TBL()
    USE DS_Gmsh2Solver
    USE hashtbl
    IMPLICIT NONE
    INTEGER::I,J,N1(2),ET1,NEDGE1,TBL_LEN
	CHARACTER(LEN=:),ALLOCATABLE::KEY1
    TYPE(DICT_DATA)::VAL1
    
    
	TBL_LEN=NNODE	
    CALL EDGE_TBL.INIT(TBL_LEN)
	DO I=1,NEL
		ET1=ELEMENT(I).ET
		NEDGE1=ELTTYPE(ET1).NEDGE
		DO J=1,NEDGE1
            N1=ELEMENT(I).NODE(ELTTYPE(ET1).EDGE(:,J))
            CALL EDGE_TBL.KEY(KEY1,N1(:),2)
            VAL1.IEL=I;VAL1.ISE=J
			CALL EDGE_TBL.PUT(TRIM(KEY1),VAL1)
			
		ENDDO
	END DO 
ENDSUBROUTINE

SUBROUTINE SETUP_FACE_TBL()
    USE DS_Gmsh2Solver
    USE hashtbl
    IMPLICIT NONE
    INTEGER::I,J,N1(4),ET1,NFACE1,TBL_LEN
	CHARACTER(LEN=:),ALLOCATABLE::KEY1
    TYPE(DICT_DATA)::VAL1
    
	TBL_LEN=NNODE	
    CALL FACE_TBL.INIT(TBL_LEN)
	DO I=1,NEL
		ET1=ELEMENT(I).ET
		NFACE1=ELTTYPE(ET1).NFACE
		DO J=1,NFACE1
            N1=ELEMENT(I).NODE(ELTTYPE(ET1).FACE(1:,J))
            CALL FACE_TBL.KEY(KEY1,N1(:),ELTTYPE(ET1).FACE(0,J))
            VAL1.IEL=I;VAL1.ISE=J
			CALL FACE_TBL.PUT(TRIM(KEY1),VAL1)
			
		ENDDO
	END DO 
ENDSUBROUTINE    
    
    
SUBROUTINE GEN_ZEROTHICKNESS_ELEMENT()
	USE DS_Gmsh2Solver
	USE HASHTBL
    USE IFPORT
	IMPLICIT NONE
	INTEGER::I,J,K,IPGROUP1,NNODE1,NODE1(50),NVAL1,N1,N2,IEL1,IEL2
	CHARACTER(LEN=:),ALLOCATABLE::KEY1
	TYPE(DICT_DATA),ALLOCATABLE::VAL1(:)
	REAL(8)::XY1(3)
	
	DO I=1,nphgp
		IPGROUP1=PHGPNUM(I)
		CALL LOWCASE(physicalgroup(IPGROUP1).ET,LEN_TRIM(physicalgroup(IPGROUP1).ET))
		SELECT CASE(physicalgroup(IPGROUP1).ET)
			CASE('zt_cpe4_spg','zt_prm6_spg')
				!IF(.NOT.EDGE_TBL.IS_INIT) CALL SETUP_EDGE_TBL()
				!IF(physicalgroup(IPGROUP1).NDIM/=1) STOP "DIMENSION ERROR.ERRORS IN SUB GEN_ZEROTHICKNESS_ELEMENT()."
				!IF(physicalgroup(IPGROUP1).ET_GMSH/=1) STOP "ET ERROR.ERRORS IN SUB GEN_ZEROTHICKNESS_ELEMENT()."
				IF(MOD(physicalgroup(IPGROUP1).NEL,2)/=0) STOP "ELEMENT NUMBER ERROR.ERRORS IN SUB GEN_ZEROTHICKNESS_ELEMENT()."
				
				N1=physicalgroup(IPGROUP1).NEL/2
				
				DO J=1,N1
					IEL1=physicalgroup(IPGROUP1).ELEMENT(J)
					IEL2=physicalgroup(IPGROUP1).ELEMENT(N1+J)
                    NODE1(1:ELEMENT(IEL1).NNODE)=ELEMENT(IEL1).NODE
					NODE1(ELEMENT(IEL1).NNODE+1:ELEMENT(IEL1).NNODE+ELEMENT(IEL2).NNODE)=ELEMENT(IEL2).NODE
					ELEMENT(IEL1).NNODE=ELEMENT(IEL1).NNODE*2
					DEALLOCATE(ELEMENT(IEL1).NODE)
					ALLOCATE(ELEMENT(IEL1).NODE,SOURCE=NODE1(1:ELEMENT(IEL1).NNODE))
				ENDDO
				
				physicalgroup(IPGROUP1).NEL=N1
				
				!RAMDOM CHECK MATCH
				DO J=1,10
					N2=MOD(IRAND(1),physicalgroup(IPGROUP1).NEL)
					N2=physicalgroup(IPGROUP1).ELEMENT(N2)
					DO K=1,ELEMENT(N2).NNODE/2
						XY1=NODE(ELEMENT(N2).NODE(K)).XY-NODE(ELEMENT(N2).NODE(K+ELEMENT(N2).NNODE/2)).XY
						IF(NORM2(XY1)>1E-6) STOP "MATCH ERROR IN SUB GEN_ZEROTHICKNESS_ELEMENT()."
					ENDDO
				ENDDO
				!REORDER NODAL NUMBER zt_cpe4_spg
				IF(physicalgroup(IPGROUP1).ET=='zt_cpe4_spg') THEN
					DO J=1,physicalgroup(IPGROUP1).NEL
						N2=ELEMENT(physicalgroup(IPGROUP1).ELEMENT(J)).NODE(3)
						ELEMENT(physicalgroup(IPGROUP1).ELEMENT(J)).NODE(3)=ELEMENT(physicalgroup(IPGROUP1).ELEMENT(J)).NODE(4)
						ELEMENT(physicalgroup(IPGROUP1).ELEMENT(J)).NODE(4)=N2
					ENDDO
				ENDIF
				
			CASE DEFAULT
				print *, "TRY TO GEN SOME KIND ZEROTHICKNESS-ELEMENT. BUT IT IS STILL UNDER DEVELOPMENT.NOT COMPLETED."
				PAUSE
		ENDSELECT
	ENDDO
	
ENDSUBROUTINE

SUBROUTINE ENLARGE_NODE()
	USE DS_Gmsh2Solver
	IMPLICIT NONE
	TYPE(node_type),ALLOCATABLE::NODE1(:)
	INTEGER::I
	TNNODE=NNODE+1000
	ALLOCATE(NODE1(TNNODE))
	NODE1(1:NNODE)=NODE
	DEALLOCATE(NODE)
	ALLOCATE(NODE,SOURCE=NODE1)
	DEALLOCATE(NODE1)
ENDSUBROUTINE

subroutine write_readme_gmsh2sinp()	
	use DS_Gmsh2Solver
	use ifport
	implicit none
	integer::i,j,item
	LOGICAL(4)::tof,pressed
	!integer,external::ipp
	character(512)::readme(512)
	
	open(2,file='D:\README_GMSH2SINP.TXT',STATUS='REPLACE')
	
	I=0
	README(IPP(I))= "//THE PROGRAM USES AN MSH FILE GENERATED BY THE GMSH AS AN INPUT FILE."
	README(IPP(I))= "//THE OUTPUT IS AN SINP FILE WHICH IS PREPARED FOR THE FEASOLVER AS AN INPUT FILE."
	README(IPP(I))= "//TO OUTPUT THE SINP FILE, ADDITIONAL INFORMATION IS NEED, WHICH IS PROVIDED BY SOME SELFDEFINED KEYWORDS."
	README(IPP(I)) ="//THE KEYWORD STRUCTURE IS EXPLAINED HEREIN"
	README(IPP(I)) ="//NOTE THAT NOT ALL OF THEM ARE REQUIRED TO SET UP THE MODEL."
	README(IPP(I)) ="//AND THAT THE DEPTH DIRECTION IS THE Z AXIS, WHOSE POSITIVE DIRECTION IS UPWARD.(KEEP CONSISTENCY WITH FEASOLVER)"
	README(IPP(I)) ="//所有的边界条件及荷载均以点作用的形式输出.ET=ELT_BC_OR_LOAD的单元组不输出."

	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//$ELT_BC"
	README(IPP(I))=  "//"//'"'//"THE KEYWORD ELT_BC IS USED TO DEFINE THE BOUNDARY CONDITION APPLIED ON THE ELEMENT GROUP."//'"'
	README(IPP(I))="//NOTE THAT THE FINAL BC VALUE OF A NODE IS WHAT INPUT BY THE FIRST ELEMENT SHARING THE NODE."
	README(IPP(I))= "//{NELT_BC(I)} "
	README(IPP(I)) = "//{NO(I),GROUP(I),NDIM(I),DOF(I),STEPFUNC(I),VALUE(R)}  // NDIM=0,1,2,3 分别表示点、线、面、体的约束. 因为边界条件（如位移，水头等）不具叠加性，所以NDIM取值对边界条件无影响."
	README(IPP(I))=  "//                                                      // DOF(I)=1,2,...,7,分别表示约束X,Y,Z,H,MX,MY,MZ "
    README(IPP(I))=  "//{......}   //共NELT_BC行. "
	README(IPP(I)) = "//$ENDELT_BC"
	
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//$ELT_LOAD"
	README(IPP(I))=  "//"//'"'//"THE KEYWORD ELT_LOAD IS USED TO DEFINE THE LOAD APPLIED ON THE ELEMENT GROUP."//'"'
	README(IPP(I))="//NOTE THAT, IF THE NDIM=1,2 OR 3, THE FINAL BC VALUE OF A NODE IS THE SUM OF ALL THE ELEMENT SHARING THE NODE."
	README(IPP(I)) = "//{NELT_LOAD(I)} "
	README(IPP(I)) = "//{NO(I),GROUP(I),NDIM(I),DOF(I),STEPFUNC(I),VALUE(R)}  // NDIM=0,1,2,3 分别表示点、线、面、体的荷载."
	README(IPP(I)) =  "//{......}   //共NELT_LOAD行. "
	README(IPP(I)) ="//目前可处理:NDIM=0(点荷载); "
	README(IPP(I)) ="//			 NDIM=1的情况目前可处理作用在2节点、3节点以及5节点的单元边，且仅限于均布线荷载;"
	README(IPP(I)) ="//			 NDIM=2的情况目前可处理作用在3节点（三角形）、4节点（四边形）、6节点（二次三角形）、8节点（二次四边形）以及15节点（四次三角形）的单元面，且仅限于均布面荷载;"
	README(IPP(I)) ="//			 NDIM=3的情况目前可处理作用在4节点（四面体）、10节点（二次四面体）、6节点(Prism)、8节点(长方体)，且仅限于均布体荷载;"
	README(IPP(I)) = "//$ENDELT_LOAD"		
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//$ELT_SPGFACE"
	README(IPP(I))=  "//"//'"'//"THE KEYWORD ELT_SPGFACE IS USED TO DEFINE THE SEEPAGEFACE CONDITION APPLIED ON THE ELEMENT GROUP."//'"'
	README(IPP(I)) = "//{NELT_SPGFACE(I)} "
	README(IPP(I)) = "//{NO(I),GROUP(I),NDIM(I),DOF(I),STEPFUNC(I),VALUE(R)}  // NDIM=0,1,2,3 均可，VALUE=任意值（程序内部设为Y或Z）"
	README(IPP(I)) =  "//{......}   //共NELT_SPGFACE行. "
	README(IPP(I)) = "//$ENDELT_SPGFACE"	
	
	!README(IPP(I)) ="\N//******************************************************************************************************"C
	!README(IPP(I)) = "//$MODELGROUP"
	!README(IPP(I))=  "//"//'"'//"THE KEYWORD MODELGROUP IS USED TO DEFINE THE GROUPS DEFINING THE MODEL, EXCLUDING THE GROUPS DEFINING THE BOUNDARIES."//'"'
	!README(IPP(I)) = "//{NMODELGROUP(I)} "
	!README(IPP(I)) = "//{GROUPID1(I),GROUPID2(I),...,GROUPIDN(I)}  // N=NMODELGROUP "
	!README(IPP(I)) = "//$ENDMODELGROUP"	
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//$MODELDIMENSION"
	README(IPP(I))=  "//"//'"'//"THE KEYWORD MODELDIMENSION IS USED TO DEFINE THE MODEL DIMENSION."//'"'
	README(IPP(I)) = "//{MODELDIMENSION(I)}   // MODELDIMENSION=2, 3"  
	README(IPP(I)) = "//$ENDMODELDIMENSION"	
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//$GROUPPARAMETER"
	README(IPP(I))=  "//"//'"'//"THE KEYWORD GROUPPARAMETER IS USED TO DEFINE THE MATERIAL AND ELEMENT TYPE IN THE MODELGROUP ELEMENTS."//'"'
	README(IPP(I)) = "//{NGROUPPARAMETER(I)}"  
	README(IPP(I)) = "//{GROUDID(I),MATID(I),ET(A32)}  //ET=CPE3,CPE6,TET4,TET10,... "
	README(IPP(I)) = "//{......}  // 共NGROUPPARAMETER 个"     
	README(IPP(I)) = "//$ENDGROUPPARAMETER"	
    
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//$WSP"
	README(IPP(I))=  "//"//'"'//"THE KEYWORD WSP IS USED TO DEFINE THE BOUNDRAY LINE WHERE A WATER SURFACE PROFILE IS CALCULATED."//'"'
	README(IPP(I)) = "//{NWSP(I)}"  
	README(IPP(I)) = "//{GROUPID(I),CHARACTERPOINTGROUP(I),STARTPOINTGROUP(I)}*{NWSP(I)}  // 每行3个，共NWSP行，STARTPOINTGROUP为输出起点(这个GROUP只含一个点)，根据起点按顺序输出线上各点。"
	README(IPP(I)) = "//注：CHARACTERPOINTGROUP，为从上游到下游流道的特征点。这些特征点将流道分成不同性质的子流道段"
	README(IPP(I)) = "//$ENDWSP"	
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//$DATAPOINT"
	README(IPP(I))=  "//"//'"'//"THE DataPoint IS USED TO OUTPUT THE DATA POINTS ALONG A POLYLINE YOU WANT TO LIST THE COMPUTATED RESULT."//'"'
	README(IPP(I)) = "//{NDataPoint(I)}"  
	README(IPP(I)) = "//{GROUPID(I)[,ORDER(I)=0,STARTPOINTGROUP(I)]}*{NDataPoint(I)}  // 每行3个，共NDataPoint行.STARTPOINTGROUP为输出起点(这个GROUP只含一个点)，根据起点按顺序输出线上各点。"
	README(IPP(I)) = "//ORDER=0,表示输出散点的形式（不排序）.=1，表示以STARTPOINTGROUP为起点按点的顺序输出。"
	README(IPP(I)) = "//$ENDDATAPOINT"		
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//$INCLUDEFILE"
	README(IPP(I))=  "//"//'"'//"THE INCLUDEFILE IS USED TO READ IN DATA IN ANOTHER FILE."//'"'
	README(IPP(I)) = "//{NINCLUDEFILE(I)}"  
	README(IPP(I)) = "//{ABSOLUTED PATH OF THE FILES}*{NINCLUDEFILE(I)}  // 共NINCLUDEFILE行."
	README(IPP(I)) = "//$ENDINCLUDEFILE"    
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//$COPYFILE"
	README(IPP(I))=  "//"//'"'//"THE COPYFILE IS USED TO COPY DATA IN THE FILES TO THE OUTPUT FILE."//'"'
	README(IPP(I)) = "//{NCOPYFILE(I)}"  
	README(IPP(I)) = "//{ABSOLUTED PATH OF THE FILES}*{NCOPYFILE(I)}  // 共NCOPYFILE行."
	README(IPP(I)) = "//$ENDCOPYFILE" 
    
	do 	j=1,i
		item=len_trim(readme(j))
		write(2,20) readme(j)
	end do
	
	tof=system("D:\README_GMSH2SINP.TXT")	
	
20	format(a<item>)
	
	contains
	
	integer function ipp(i) 
		implicit none
		integer,intent(inout)::i
		
		i=i+1
		ipp=i
		if(ipp>512) then
			print *, "The size of README is 512."
			stop
		end if
				
	end function
	
end subroutine

subroutine Tosolver()
	use DS_Gmsh2Solver
	implicit none
	integer::unit,item,n1,i,j,width,ITEM1
	
	
	unit=20
	open(unit,file=resultfile,status='replace')
	item=len('Title')
	write(unit,100) 'Title'
	item=len_trim(title)
	write(unit,100) title
	
	write(unit,110) nnode,1,modeldimension
	if(modeldimension==2) write(unit,112) "X","Y"
	if(modeldimension==3) write(unit,112) "X","Y","Z"
	do i=1,nnode
		n1=Noutputorder(i)
		write(unit,111) node(n1).xy(1:modeldimension)
	end do
	do i=1,nphgp
		item=len_trim(adjustL(physicalgroup(phgpnum(i)).et))
		ITEM1=len_trim(adjustL(physicalgroup(phgpnum(i)).NAME))
		
		call lowcase(physicalgroup(phgpnum(i)).et,len(physicalgroup(phgpnum(i)).et))
		if(.NOT.physicalgroup(phgpnum(i)).ISMODEL) cycle  !ET=ELT_BC_OR_LOAD的单元组不输出
		
		write(unit,120) physicalgroup(phgpnum(i)).nel,phgpnum(i),physicalgroup(phgpnum(i)).et, &
						physicalgroup(phgpnum(i)).mat,physicalgroup(phgpnum(i)).NAME
		item=len_trim(elttype(physicalgroup(phgpnum(i)).et_gmsh).description)
		write(unit,122) elttype(physicalgroup(phgpnum(i)).et_gmsh).description
		item=elttype(physicalgroup(phgpnum(i)).et_gmsh).nnode
		write(unit,123) ((J),J=1,item)		
		do j=1,physicalgroup(phgpnum(i)).nel
			n1=physicalgroup(phgpnum(i)).element(j)
			item=element(n1).nnode
			write(unit,121) node(element(n1).node).inode
		end do
	end do
	
	if(nnodalBC>0) then
		write(unit,130) nnodalbc
		write(unit, 132) 
		do i=1,nnodalBC
			write(unit,131) node(nodalBC(i).node).inode,nodalBC(i).dof,nodalBC(i).value,nodalBC(i).sf
		end do
	end if
	
	if(nnodalload>0) then
		write(unit,140) nnodalload
		write(unit,142)
		do i=1,nnodalload
			write(unit,141) node(nodalload(i).node).inode,nodalload(i).dof,nodalload(i).value,nodalload(i).sf
		end do
	end if
	
	if(nspgface>0) then
		do i=1, nelt_spgface
			write(unit,150) elt_spgface(i).n2-elt_spgface(i).n1+1,elt_spgface(i).sf
			write(unit,152)
			
			write(unit,151) (node(SPGFACE(j).node).inode,j=elt_spgface(i).n1,elt_spgface(i).n2)
		end do
	end if
    
    if(nwsp>0) then
        WRITE(UNIT,160) NWSP
        do i=1, nwsp
			write(unit,162) wsp(i).nchnode-1,WSP(I).NNODE
			DO J=1,physicalgroup(WSP(I).CHGROUP).nel-1
				write(unit,163) wsp(i).chnode(j),wsp(i).chnode(j+1)
			END DO
			write(unit,161) node(WSP(I).node).inode
		end do
    end if
	
    if(ndatapoint>0) then
        WRITE(UNIT,170) NDATAPOINT
        do i=1, NDATAPOINT
			write(unit,171) DATAPOINT(I).NNODE
			write(unit,172) node(DATAPOINT(I).NODE).inode
		end do
    end if
	

	DO I=1,NCOPYFILE
		CALL DO_COPYFILE(COPYFILE(I),UNIT)			
	END DO
	
	
100 FORMAT(A<ITEM>)	
101 FORMAT('"',A<ITEM>,'"')

110 FORMAT(/,'NODE,NUM=',I7,',DATAPACKING=',I2,',DIMENSION=',I2)
111 FORMAT(<MODELDIMENSION>(F15.7,1X))
112 FORMAT("//",<MODELDIMENSION>(A15,1X))

120 FORMAT(/'ELEMENT,NUM=',I7,',SET=',I3,',ET=',A<ITEM>,',MATID=',I3,',TITLE=',A<ITEM1>)
121 FORMAT(<ITEM>(I7,1X))
122 FORMAT("//",A<ITEM>)
123 FORMAT("// ", <ITEM>("N",I<width(j)>,5X))

130 FORMAT(/'BC,NUM=',I7) 
131 FORMAT(I7,1X,I2,1X,E15.7,1X,I4)
132 FORMAT("// ","NODE DOF VALUE [STEPFUNC.]")

140 FORMAT(/'LOAD,NUM=',I7)
141 FORMAT(I7,1X,I2,1X,E15.7,1X,I4)
142 FORMAT("// ","NODE DOF VALUE [STEPFUNC.] ")

150 FORMAT(/'SEEPAGE FACE,NUM=',I7,', sf=',I7) 
151 FORMAT(10(I7,1X))
152 FORMAT("// NODE")

160 FORMAT(/"WSP,NUM=",I7,",CP=0|1") 
161 FORMAT(10(I7,1X))
162 FORMAT(/"//{NSEGMENT(I)=",I7,",NNODE(I)=",I7,",Q(R),B(R),UYBC,DYBC,[Q_SF,UYBC_SF,DYBC_SF,CALTYPE,G,Kn]}")
163 FORMAT("//{UNODE(I)=",I7,",DNODE(I)=",I7,",n(R),So(R),[Profileshape(1),Profileshape(2)]}")

170 FORMAT(/"DATAPOINT,NUM=",I7)
171	FORMAT(/I7)
172 FORMAT(10(I7,1X))

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

integer function width(j)
	implicit none
	integer::j
	if(j<10) then
		width=1
	else
		width=2
	end if
	
end function

subroutine elt_bc_load_translate()
	use DS_Gmsh2Solver
	implicit none
	integer::i,j,k,n1,n2,nc,iar1(5)
	real(8)::LAV1
	integer,allocatable::nodalload1(:),node1(:)
	real(8),allocatable::load1(:)
	
	! 荷载具可叠加，边界条件不具有可加性。分布力转化为节点力时要积分，但分布位移转化为节点位移时不需积分，节点位移等于分布位移。
	if(nelt_load>0) allocate(nodalLoad(maxnnodalLoad))
	if(nelt_bc>0) allocate(nodalBC(maxnnodalBC))
	if(nelt_spgface>0) allocate(spgface(maxnspgface))
	allocate(nodalload1(nnode),node1(nnode))
	if(nelt_load>0) allocate(load1(nnode))
	do i=1, nelt_load
		select case(elt_load(i).ndim)
			case(0) !ndim=0,是点荷载，各个单元中的同一节点，不具可加性，均指同一值。
				nodalload1=0
				do j=1,physicalgroup(elt_load(i).group).nel
					n1=physicalgroup(elt_load(i).group).element(j)
					do k=1,element(n1).nnode
						n2=element(n1).node(k)
						if(nodalload1(n2)==0) then
							nnodalLoad=nnodalLoad+1 
							if(nnodalLoad>maxnnodalLoad) call enlargenodalLoad()
							nodalLoad(nnodalLoad).node=n2
							nodalLoad(nnodalLoad).dof=elt_load(i).dof
							nodalLoad(nnodalLoad).sf=elt_load(i).sf
							nodalLoad(nnodalLoad).value=elt_load(i).value
							nodalload1(n2)=1
						end if
					end do
				end do
								
			case default
				load1=0
				do j=1,physicalgroup(elt_load(i).group).nel
					n1=physicalgroup(elt_load(i).group).element(j)
					call Element_LAV_Cal(n1,LAV1) 
					LAV1=LAV1*elt_load(i).value
					do k=1,element(n1).nnode
						n2=element(n1).node(k)
						load1(n2)=load1(n2)+Elttype(element(n1).et).weight(k,1)*LAV1 !均布荷载		 				
					end do					
				end do				
				do j=1,nnode
					if(abs(load1(j))>1e-10) then
						nnodalLoad=nnodalLoad+1 
						if(nnodalLoad>maxnnodalLoad) call enlargenodalLoad()
						nodalLoad(nnodalLoad).node=j
						nodalLoad(nnodalLoad).dof=elt_load(i).dof
						nodalLoad(nnodalLoad).sf=elt_load(i).sf
						nodalLoad(nnodalLoad).value=load1(j)
					end if
				end do

		end select
		
	end do
	
	do i=1, nelt_bc
		!**********目前认为边界条件都没有可加性,NDIM=1,2,3的情况与NDIM=0的情况一样**********。
		nodalload1=0
		do j=1,physicalgroup(elt_bc(i).group).nel
			n1=physicalgroup(elt_bc(i).group).element(j)
			do k=1,element(n1).nnode
				n2=element(n1).node(k)
				if(nodalload1(n2)==0) then
					nnodalBC=nnodalBC+1 
					if(nnodalBC>maxnnodalBC) call enlargenodalBC()
					nodalBC(nnodalBC).node=n2
					nodalBC(nnodalBC).dof=elt_bc(i).dof
					nodalBC(nnodalBC).sf=elt_bc(i).sf
					nodalBC(nnodalBC).value=elt_bc(i).value
					nodalload1(n2)=1 !多次出现时以第一次出现的值为准。
				end if
			end do			
		end do		
	end do 
	
	do i=1, nelt_spgface
		!**********目前认为边界条件都没有可加性**********。
		nodalload1=0
		elt_spgface(i).n1=nspgface+1
		do j=1,physicalgroup(elt_spgface(i).group).nel
			n1=physicalgroup(elt_spgface(i).group).element(j)
			do k=1,element(n1).nnode
				n2=element(n1).node(k)
				if(nodalload1(n2)==0) then
					nspgface=nspgface+1 
					if(nspgface>maxnspgface) call enlargespgface()
					spgface(nspgface).node=n2
					spgface(nspgface).dof=elt_spgface(i).dof
					spgface(nspgface).sf=elt_spgface(i).sf
					spgface(nspgface).value=elt_spgface(i).value
					nodalload1(n2)=1 !多次出现时以第一次出现的值为准。
				end if
			end do			
		end do
		elt_spgface(i).n2=nspgface
	end do 
    
    do i=1,nwsp
        nodalload1=0
        !node1=0
        do j=1,physicalgroup(wsp(i).group).nel
            n1=physicalgroup(wsp(i).group).element(j)
            do k=1,element(n1).nnode
                n2=element(n1).node(k)
                if(nodalload1(n2)==0)then
                    wsp(i).nnode= wsp(i).nnode+1
                    node1(wsp(i).nnode)=n2
                    nodalload1(n2)=1
                end if                
            end do
        end do
		
        wsp(i).nchnode=physicalgroup(wsp(i).chgroup).nel
        allocate(wsp(i).node(wsp(i).nnode),wsp(i).chnode(wsp(i).nchnode))
		
		call SortByPath(wsp(i).group,wsp(i).spgroup,wsp(i).node,wsp(i).nnode)
        !do j=1,wsp(i).nnode
        !    
        !    do k=j+1,wsp(i).nnode
        !        !从小到大
        !        if(node(node1(k)).xy(1)<node(node1(j)).xy(1)) then
        !            n1=node1(k)
        !            node1(k)=node1(j)
        !            node1(j)=n1
        !        end if
        !    end do
        !end do
		
		
		
        
        !if(wsp(i).xdirection==1) then
        !    wsp(i).node=node1(1:wsp(i).nnode)
        !else
        !    wsp(i).node=node1(wsp(i).nnode:1:-1)
        !endif
        
		do j=1,wsp(i).nchnode
			do k=1,wsp(i).nnode
				if(wsp(i).node(k)==element(physicalgroup(wsp(i).chgroup).element(j)).node(1)) then
					wsp(i).chnode(j)=k
					exit
				end if
			end do
        end do
		
        do j=1,wsp(i).nchnode
            do k=j+1,wsp(i).nchnode
               if(wsp(i).chnode(j)>wsp(i).chnode(k)) then
                   n1=wsp(i).chnode(j)
                   wsp(i).chnode(j)=wsp(i).chnode(k)
                   wsp(i).chnode(k)=n1
               end if
            end do
        end do
		
    end do
	
    do i=1,nDataPoint
        nodalload1=0
        !node1=0
        do j=1,physicalgroup(DataPoint(i).group).nel
            n1=physicalgroup(DataPoint(i).group).element(j)
            do k=1,element(n1).nnode
                n2=element(n1).node(k)
                if(nodalload1(n2)==0)then
                    DataPoint(i).nnode= DataPoint(i).nnode+1
                    node1(DataPoint(i).nnode)=n2
                    nodalload1(n2)=1
                end if                
            end do
        end do
		
        allocate(DataPoint(i).node(DataPoint(i).nnode))
		
		IF(DATAPOINT(I).ORDER==0 .or. physicalgroup(DataPoint(i).group).ndim/=1) THEN
			DataPoint(i).node(1:DataPoint(i).nnode)=NODE1(1:DataPoint(i).nnode)
		ELSE
			call SortByPath(DataPoint(i).group,DataPoint(i).spgroup,DataPoint(i).node,DataPoint(i).nnode)
		END IF		

		
    end do	
	
	if(allocated(nodalload1)) deallocate(nodalload1)
	if(allocated(load1)) deallocate(load1)
	if(allocated(node1)) deallocate(node1)
end subroutine


subroutine SortByPath(igroup,ispgroup,Local_Node,NLNDE)
	use DS_Gmsh2Solver
	implicit none
	integer,intent(in)::igroup,ispgroup,NLNDE
	integer,intent(out)::Local_node(NLNDE)
	integer::i,j,k,iar1(5),n1,n2
	integer,allocatable::element1(:)
	
	
	allocate(element1(physicalgroup(igroup).nel))
	element1=0
	Local_Node(1)=element(physicalgroup(ispgroup).element(1)).node(1)	
	j=1	
	do while(j<=NLNDE-1)
		do k=1,physicalgroup(igroup).nel
			if(element1(k)/=0) cycle
			
			n1=element(physicalgroup(igroup).element(k)).nnode
			iar1(1:n1)=element(physicalgroup(igroup).element(k)).node(1:n1)
			n1=iar1(1)-Local_Node(j)
			n2=iar1(2)-Local_Node(j)
			if(n1*n2/=0) cycle
			
			element1(k)=1
			
			if(n1==0) then
				select case(physicalgroup(igroup).ET_GMSH)
				case(1)
					Local_Node(j+1)=iar1(2)
					j=j+1
				case(8)
					Local_Node(j+1)=iar1(3)
					Local_Node(j+2)=iar1(2)
					j=j+2
				case(27)
					Local_Node(j+1)=iar1(3)
					Local_Node(j+2)=iar1(4)
					Local_Node(j+3)=iar1(5)
					Local_Node(j+4)=iar1(2)
					j=j+4
				case default
					print *, "The element is expected to be a line element. datapoint(i),i=",i
					stop
				end select
			else
				select case(physicalgroup(igroup).ET_GMSH)
				case(1)
					Local_Node(j+1)=iar1(1)
					j=j+1
				case(8)
					Local_Node(j+1)=iar1(3)
					Local_Node(j+2)=iar1(1)
					j=j+2
				case(27)
					Local_Node(j+1)=iar1(5)
					Local_Node(j+2)=iar1(4)
					Local_Node(j+3)=iar1(3)
					Local_Node(j+4)=iar1(1)
					j=j+4
				case default
					print *, "The element is expected to be a line element. datapoint(i),i=",i
					stop
				end select						
			end if
			exit
		end do
	enddo	

	deallocate(element1)

end subroutine

subroutine Element_LAV_Cal(ienum,LAV) !计算单元的长度、面积或体积
	use DS_Gmsh2Solver
	implicit none	
	integer,intent(in)::ienum
	real(8),intent(out)::LAV
	REAL(8),EXTERNAL::TRIAREA,TETVOL
	integer::i,j
	real(8)::a1(4,3)=0
	
	select case(element(ienum).et)
		case(1,8,27) !Length for a Line element 
			LAV=((node(element(ienum).node(1)).xy(1)-node(element(ienum).node(2)).xy(1))**2+ &
			(node(element(ienum).node(1)).xy(2)-node(element(ienum).node(2)).xy(2))**2+ &
			(node(element(ienum).node(1)).xy(3)-node(element(ienum).node(2)).xy(3))**2)**0.5
		case(2,9,23) !三角形面积
						
			do i=1,3
				a1(i,:)=node(element(ienum).node(i)).xy
			end do
			LAV=TRIAREA(A1(1:3,1:3))
		case(3,16) !四边形面积
			do i=1,4
				a1(i,:)=node(element(ienum).node(i)).xy
			end do
			LAV=TRIAREA(A1(1:3,1:3))
			
			A1(2,:)=A1(1,:)
			LAV=LAV+TRIAREA(A1(2:4,1:3))
			
		case(4,11) !四面体体积	
			do i=1,4
				a1(i,:)=node(element(ienum).node(i)).xy
			end do
			LAV=TETVOL(A1(1:4,1:3))
		
		case(5,6)
			print *, "Still under work."
			pause
		case default
			print *, "Still under work."
			pause
		
	end select

end subroutine

!xyz(1:3,:) 三角形三个节点的空间坐标,第一维为节点。
real(8) function TriArea(xyz)  
	implicit none
	real(8),intent(in)::xyz(3,3)
	integer::i,j
	real(8)::a1(3,3)
	
	a1=xyz
	do i=2,3
		do j=1,3
			a1(i,j)=xyz(i,j)-xyz(1,j)
		end do
	end do
	TriArea=(a1(2,2)*a1(3,3)-a1(3,2)*a1(2,3))**2+ &
		(a1(2,1)*a1(3,3)-a1(3,1)*a1(2,3))**2+ &
		(a1(2,1)*a1(3,2)-a1(3,1)*a1(2,2))**2
	TriArea=0.5*TriArea**0.5

end function

!xyz(1:4,:) 四面体四个节点的空间坐标,第一维为节点。
real(8) function TetVOL(xyz)  
	implicit none
	real(8),intent(in)::xyz(4,3)
	integer::i,j
	real(8)::a1(3,3)
		
	do i=2,4
		do j=1,3
			a1(i-1,j)=xyz(i,j)-xyz(1,j)
		end do
	end do
	TetVOL=a1(1,1)*(a1(2,2)*a1(3,3)-a1(3,2)*a1(2,3))- &
			 a1(1,2)*(a1(2,1)*a1(3,3)-a1(3,1)*a1(2,3))+ &
			 a1(1,3)*(a1(2,1)*a1(3,2)-a1(3,1)*a1(2,2))
	TetVOL=1./6.*ABS(TetVOL)
end function


!Enlarge Array ELT() by an increment of 1000
subroutine enlargenodalLoad()
	use DS_Gmsh2Solver
	implicit none
	
	type(bc_tydef)::nodalLoad1(maxnnodalLoad)
	
	nodalLoad1=nodalLoad(1:maxnnodalLoad)
	deallocate(nodalLoad)	
	allocate(nodalLoad(maxnnodalLoad+1000))
	nodalLoad(1:maxnnodalLoad)=nodalLoad1
	maxnnodalLoad=maxnnodalLoad+1000

end subroutine


subroutine enlargenodalBC()
	use DS_Gmsh2Solver
	implicit none

	type(bc_tydef)::nodalBC1(maxnnodalBC)
	
	nodalBC1=nodalBC(1:maxnnodalBC)
	deallocate(nodalBC)	
	allocate(nodalBC(maxnnodalBC+1000))
	nodalBC(1:maxnnodalBC)=nodalBC1
	maxnnodalBC=maxnnodalBC+1000

end subroutine

subroutine enlargespgface()
	use DS_Gmsh2Solver
	implicit none

	type(bc_tydef)::spgface1(maxnspgface)
	
	spgface1=spgface(1:maxnspgface)
	deallocate(spgface)	
	allocate(spgface(maxnspgface+1000))
	spgface(1:maxnspgface)=spgface1
	maxnspgface=maxnspgface+1000

end subroutine

subroutine readin()
	use DS_Gmsh2Solver
	use dflib
	implicit none
	integer:: unit1,i
	
	character(512) term,keyword
	!character(512)::nme
	CHARACTER(3)        drive
	CHARACTER(512)      dir
	CHARACTER(512)      name
	CHARACTER(16)      ext
	integer(4)::length,msg

	term="Msh Files(*.msh), *.msh; &
				Data Files(*.dat),*.dat; &
			  Prof.Cao Program files(*.z7),*.z7; &
			  All Files(*.*),*.*;"
	term=trim(term)
	call setmessageqq(term,QWIN$MSG_FILEOPENDLG)
	term=' '


	
	call Initialize_et2numNodes()
	CALL ET_EDGE_FACE()
    I=0
	DO WHILE(I<=NINCLUDEFILE)
		unit1=1
		IF(I==0) THEN
			open(UNIT1,file=' ',status='old' )
			inquire(UNIT1,name=INCLUDEFILE(I))
			length = SPLITPATHQQ(INCLUDEFILE(I), drive, dir, name, ext)
			title=trim(name)
			resultfile=trim(drive)//trim(dir)//trim(name)//'.sinp'
		ELSE
			open(UNIT1,file=INCLUDEFILE(I),status='old' )
		ENDIF
		
		call read_execute(unit1)
		print * ,'DONE IN READING THE FILE:'//TRIM(INCLUDEFILE(I))
		close(UNIT1)
        I=I+1
	END DO
	
	 

	do i=1, nmodelgroup
		physicalgroup(modelgroup(i)).ismodel=.true.
    end do
	
    
    
	return

 end subroutine



subroutine read_execute(unit)

	use DS_Gmsh2Solver
	implicit none
	integer::unit,ef,iterm,i,strL
	parameter(iterm=512)
	character(iterm)::term
	character(1)::ch
	
	ef=0
	
	do while(ef==0)
		term=''
		read(unit,999,iostat=ef) term	
		
		if(ef<0) exit
		term=adjustl(term)
		strL=len_trim(term)
		if(strL==0) cycle
		do i=1,strL !remove 'Tab'
			if(term(i:i)/=char(9)) exit
		end do
		term=term(i:strL)
		term=adjustl(term)
		strL=len_trim(term)
		if(strL==0) cycle		
		!write(ch,'(a1)') term
		if(term(1:1)=='$') then
			term=term(2:strL)		
			call lowcase(term,iterm)
			!call translatetoproperty(term)			
			term=adjustl(trim(term))
			call kwcommand(term,unit)			 	
		end if
	end do

999	format(a<iterm>)

end subroutine

subroutine kwcommand(term,unit)		
	use DS_Gmsh2Solver
	implicit none
	integer,intent(in)::unit	
	character(512),intent(in)::term
	integer::n1,n2,n3,nacw=0
	integer::en1(15)=0
	logical,external::isacw
    logical::ischeckacw=.false.
	integer::strL1,inode,i
	real(8)::at1(100),xy1(3,3)=0
	
	integer::nmax,nread,maxset,nset
	
    
	parameter(nmax=100)
	parameter(maxset=100)
	
	real(8)::ar(nmax)
	character(32)::set(maxset)	,str1
	
	strL1=len_trim(term)
	
	select case (term)
		case('nodes')
			read(unit,*) nnode
			allocate(node(nnode))
			do i=1,nnode
				call skipcomment(unit)				
				read(unit,*) at1(1:4)
				inode=int(at1(1))
				node(inode).xy(1:3)=at1(2:4)
			end do
		case('endnodes')
			strL1=5
			write(*, 20) "Nodes"	
							
		case('elements')
			call skipcomment(unit)
			read(unit,*) nel
			allocate(element(nel))
            nacw=-1
			do i=1,nel				
				call strtoint(unit,ar,nmax,nread,nmax,set,maxset,nset)
				n1=int(ar(1))
				element(n1).et=int(ar(2))
				
				element(n1).nnode=elttype(element(n1).et).nnode
				
				element(n1).ntag=int(ar(3))
				if(element(n1).ntag>0) then
					allocate(element(n1).tag(element(n1).ntag))
					element(n1).tag(1:element(n1).ntag)=int(ar(3+1:3+element(n1).ntag))
					physicalgroup(element(n1).tag(1)).nel=physicalgroup(element(n1).tag(1)).nel+1
					physicalgroup(element(n1).tag(1)).ET_GMSH=ELEMENT(N1).ET !假定同一GROUP的单元类型相同
                    
                    if(nacw/=element(n1).tag(2)) then                        
                        ischeckacw=.true.
                        nacw=element(n1).tag(2)
                    else
                        ischeckacw=.false.
                    end if  
				end if
	  
				allocate(element(n1).node(element(n1).nnode))
				element(n1).node(1:element(n1).nnode)=int(ar(3+element(n1).ntag+1:3+element(n1).ntag+element(n1).nnode))
				select case(element(n1).et)
                    case(2,9,23)
                        if(ischeckacw) then
                            xy1(1,1:3)=node(element(n1).node(1:3)).xy(1)
                            xy1(2,1:3)=node(element(n1).node(1:3)).xy(2)
                            xy1(3,1:3)=node(element(n1).node(1:3)).xy(3)
                            if(.not.isacw(xy1(1,1),xy1(2,1),xy1(3,1),xy1(1,2),xy1(2,2),xy1(3,2),xy1(1,3),xy1(2,3),xy1(3,3))) then
                                print *, 'Please reverse the normal direction of geometrical entity=',element(n1).tag(2)
                                pause
                                
                            end if
                        end if
                        if(element(n1).et==23) then
                            en1(4)=element(n1).node(5)
						    en1(5)=element(n1).node(8)
						    en1(6)=element(n1).node(11)
						    en1(7)=element(n1).node(4)
						    en1(8)=element(n1).node(6)
						    en1(9)=element(n1).node(7)
						    en1(10)=element(n1).node(9)
						    en1(11)=element(n1).node(10)
						    en1(12)=element(n1).node(12)
						    en1(13)=element(n1).node(14)
						    en1(14)=element(n1).node(15)
						    en1(15)=element(n1).node(13)
						    element(n1).node(4:15)=en1(4:15)                            
                        end if    
					case(11) !tet10最后两个节点的的顺序互换，以便和FEASOLVER一致。
						n2=element(n1).node(element(n1).nnode-1)
						element(n1).node(element(n1).nnode-1)=element(n1).node(element(n1).nnode)
						element(n1).node(element(n1).nnode)=n2
						
				end select			
			end do
		case('endelements')
			strL1=8
			write(*, 20) "Elements"
			
		case('elt_bc')
			call skipcomment(unit)
			read(unit,*) nelt_bc
			allocate(elt_bc(nelt_bc))
			do i=1,nelt_bc
				call strtoint(unit,ar,nmax,nread,nmax,set,maxset,nset)
				n1=int(ar(1))
				elt_bc(n1).group=int(ar(2))
				elt_bc(n1).ndim=int(ar(3))
				elt_bc(n1).dof=int(ar(4))
				elt_bc(n1).sf=int(ar(5))
				elt_bc(n1).value=ar(6)
			end do
		case('endelt_bc')
			strL1=LEN('ELT_BC')
			write(*, 20) "ELT_BC"		
		case('elt_load')
			call skipcomment(unit)
			read(unit,*) nelt_load
			allocate(elt_load(nelt_load))
			do i=1,nelt_load
				call strtoint(unit,ar,nmax,nread,nmax,set,maxset,nset)
				n1=int(ar(1))
				elt_load(n1).group=int(ar(2))
				elt_load(n1).ndim=int(ar(3))
				elt_load(n1).dof=int(ar(4))
				elt_load(n1).sf=int(ar(5))
				elt_load(n1).value=ar(6)
			end do
		case('endelt_load')
			strL1=len('elt_load')
			write(*, 20) "ELT_LOAD"
			
		case('elt_spgface')
			call skipcomment(unit)
			read(unit,*) nelt_spgface
			allocate(elt_spgface(nelt_spgface))
			do i=1,nelt_spgface
				call strtoint(unit,ar,nmax,nread,nmax,set,maxset,nset)
				n1=int(ar(1))
				elt_spgface(n1).group=int(ar(2))
				elt_spgface(n1).ndim=int(ar(3))
				elt_spgface(n1).dof=int(ar(4))
				elt_spgface(n1).sf=int(ar(5))
				elt_spgface(n1).value=ar(6)
			end do
		case('endelt_spgface')
			strL1=LEN('elt_spgface')
			write(*, 20) "ELT_SPGFACE"
            
		case('wsp')
            call skipcomment(unit)
			read(unit,*) nwsp
			allocate(wsp(nwsp))            
            do i=1,nwsp
				call strtoint(unit,ar,nmax,nread,nmax,set,maxset,nset)
                wsp(i).group=int(ar(1))
				wsp(i).CHgroup=int(ar(2))
                wsp(i).spgroup=int(ar(3))
            end do
		case('endwsp')
        	strL1=LEN('wsp')
			write(*, 20) "WSP"
			
		case('datapoint')
            call skipcomment(unit)
			read(unit,*) nDataPoint
			allocate(DataPoint(nDataPoint))            
            do i=1,nDataPoint
				call strtoint(unit,ar,nmax,nread,nmax,set,maxset,nset)
                DataPoint(i).group=int(ar(1))
				if(nread>1) DataPoint(i).order=int(ar(2))  
				if(nread>2) DataPoint(i).spgroup=int(ar(3))                
            end do
		case('enddatapoint')
        	strL1=LEN('DataPoint')
			write(*, 20) "DATAPOINT"
            
		case('modelgroup')
			call skipcomment(unit)
			read(unit,*) nmodelgroup
			allocate(modelgroup(nmodelgroup))
			call skipcomment(unit)
			read(unit,*) modelgroup(1:nmodelgroup)
		case('endmodelgroup')	
			strL1=LEN('modelgroup')
			write(*, 20) "modelgroup"		
		case('modeldimension')	
			call skipcomment(unit)
			read(unit,*) modeldimension
		case('endmodeldimension')	
			strL1=LEN('modeldimension')
			write(*, 20) "modeldimension"				
									
		case('physicalnames')
			call skipcomment(unit)
			read(unit,*) nphgp
			allocate(phgpnum(nphgp))
			!allocate(physicalgroup(nphgp))
			do i=1,nphgp
				call skipcomment(unit)
				read(unit,*) n1,n2,str1
				if(n2>maxphgp) then
					print *, 'Group Number >maxphgp. Please Renumber.,maxphgp=',maxphgp
					stop
				end if
				phgpnum(i)=n2
				physicalgroup(n2).ndim=n1
				physicalgroup(n2).name=str1
			end do
		
			
		
		case('endphysicalnames')
			strL1=LEN('PhysicalNames')
			write(*, 20) "PhysicalNames"		
		
		case('groupparameter')
			call skipcomment(unit)
			read(unit,*) N3
			
			do i=1,n3
				call skipcomment(unit)
				read(unit, *) n1,n2,str1
				physicalgroup(n1).mat=n2
				physicalgroup(n1).et=str1
				physicalgroup(n1).ismodel=.true.
			end do
		case('endgroupparameter')	
			strL1=LEN('groupparameter')
			write(*, 20) "GROUPPARAMETER"		
		CASE('includefile')
			call skipcomment(unit)
			read(unit,*) N3
			
			do i=NINCLUDEFILE+1,NINCLUDEFILE+n3
				call skipcomment(unit)
				READ(UNIT,'(A<512>)') INCLUDEFILE(I)
			end do			
			NINCLUDEFILE=NINCLUDEFILE+N3
		case('endincludefile')	
			strL1=LEN('includefile')
			write(*, 20) "INCLUDEFILE"
		CASE('copyfile')
			call skipcomment(unit)
			read(unit,*) N3
			
			do i=ncopyfile+1,ncopyfile+n3
				call skipcomment(unit)
				READ(UNIT,'(A<512>)') copyfile(I)
			end do			
			ncopyfile=ncopyfile+N3
		case('endcopyfile')	
			strL1=LEN('copyfile')
			write(*, 20) "COPYFILE"
		case('cutoffwall')
			
		case default
			strL1=len_trim(term)
			write(*, 10) term
			
	end select

10 format("The Keyword"," '",a<strL1>,"'"," cannot be recognized and be ignored.")	
20 format("The Content Introduced by the Keyord ","'",a<strL1>,"'"," is read in.")
end subroutine

subroutine skipcomment(unit)
	implicit none
	integer,intent(in)::unit
	integer::i,ef,strL
	character(512) string
	
	do while(.true.)
		read(unit,'(a512)',iostat=ef) string
		if(ef<0) then
			print *, 'file ended unexpected. sub strtoint()'
			stop
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

		if(string(1:1)/='/') then
			backspace(unit)
			exit
		end if
	end do
	
end subroutine

SUBROUTINE ET_EDGE_FACE()
	USE DS_Gmsh2Solver
	IMPLICIT NONE
	INTEGER::I,ET,AET1(33)
	
    AET1=[1:31,92,93]
	DO I=1,33
        ET=AET1(I)
		SELECT CASE(ET)
			CASE(1,8,26,27,28) !LINE
				Elttype(ET).NEDGE=1;Elttype(ET).NFACE=0
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,1)=[1,2]
			CASE(2,9,20,21,22,23,24,25) !TRIANGLE
				Elttype(ET).NEDGE=3;Elttype(ET).NFACE=1
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,1],(/2,3/))
				Elttype(ET).FACE(:,1)=[3,1,2,3,0]
			CASE(3,10,16) !QUADRANGLE
				Elttype(ET).NEDGE=4;Elttype(ET).NFACE=1
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,4,4,1],(/2,4/))
				Elttype(ET).FACE(:,1)=[4,1,2,3,4]
			
			CASE(4,11,29,30,31) !TETRAHEDRON
				Elttype(ET).NEDGE=6;Elttype(ET).NFACE=4
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,1,1,4,2,4,3,4],(/2,6/))
				Elttype(ET).FACE(:,:)=RESHAPE([3,1,2,3,0,&
											   3,1,2,4,0,&
											   3,2,3,4,0,&
											   3,3,1,4,0],(/5,4/))
			
			CASE(5,12,17,92,93) !HEXAHEDRON
				Elttype(ET).NEDGE=12;Elttype(ET).NFACE=6
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,4,4,1,&
											   5,6,6,7,7,8,8,5,&
											   1,5,2,6,3,7,4,8],(/2,12/))
				Elttype(ET).FACE(:,:)=RESHAPE([4,1,2,3,4,&
											   4,5,6,7,8,&
											   4,1,2,6,5,&
											   4,2,3,7,6,&
											   4,3,4,8,7,&
											   4,4,1,5,8],(/5,6/))
			
			CASE(6,13,18) !6-NODE PRISM
				Elttype(ET).NEDGE=9;Elttype(ET).NFACE=5
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,1,&
											   4,5,5,6,6,1,&
											   1,4,2,5,3,6],(/2,9/))
				Elttype(ET).FACE(:,:)=RESHAPE([3,1,2,3,0,&
											   3,4,5,6,0,&
											   4,1,2,5,4,&
											   4,2,3,6,5,&
											   4,3,1,4,6],(/5,5/))
												
			
			CASE(7,14,19) !5-node pyramid
				Elttype(ET).NEDGE=8;Elttype(ET).NFACE=5
				ALLOCATE(Elttype(ET).EDGE(2,Elttype(ET).NEDGE),Elttype(ET).FACE(0:4,Elttype(ET).NFACE))
				Elttype(ET).EDGE(:,:)=RESHAPE([1,2,2,3,3,4,4,1,&
											   1,5,2,5,3,5,4,5],(/2,8/))
				Elttype(ET).FACE(:,:)=RESHAPE([4,1,2,3,4,&
											   3,1,2,5,0,&
											   3,2,3,5,0,&
											   3,3,4,5,0,&
											   3,4,1,5,0],(/5,5/))
			CASE(15) !POINTS
				Elttype(ET).NEDGE=0;Elttype(ET).NFACE=0
		ENDSELECT
	ENDDO
	
ENDSUBROUTINE

subroutine Initialize_et2numNodes()
	use DS_Gmsh2Solver
	implicit none
	integer::et1
	
	Elttype(1).nnode=2
	Elttype(1).description="2-node line."
	et1=1
	call not_nodal_force_weight(et1)

	
	Elttype(2).nnode=3
	Elttype(2).description="3-node triangle."
	et1=2
	call not_nodal_force_weight(et1)
	
	Elttype(3).nnode=4
	Elttype(3).description="4-node quadrangle."
	et1=3
	call not_nodal_force_weight(et1)
	
	Elttype(4).nnode=4
	Elttype(4).description="4-node tetrahedron."
	et1=4
	call not_nodal_force_weight(et1)
	
	Elttype(5).nnode=8
	Elttype(5).description="8-node hexahedron."
	et1=5
	call not_nodal_force_weight(et1)
	
	Elttype(6).nnode=6
	Elttype(6).description="6-node prism."
	et1=6
	call not_nodal_force_weight(et1)	
	
	Elttype(7).nnode=5
	Elttype(7).description="5-node pyramid."
	
	Elttype(8).nnode=3
	Elttype(8).description="3-node second order line (2 nodes associated with the vertices and 1 with the edge)."
	et1=8
	call not_nodal_force_weight(et1)	
	
	Elttype(9).nnode=6
	Elttype(9).description="6-node second order triangle (3 nodes associated with the vertices and 3 with the edges)."
	et1=9
	call not_nodal_force_weight(et1)	
	
	Elttype(10).nnode=9
	Elttype(10).description="9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face)."
		
	Elttype(11).nnode=10
	Elttype(11).description="10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges)."
	et1=11
	call not_nodal_force_weight(et1)	
	
	Elttype(12).nnode=27
	Elttype(12).description="27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume)."
	Elttype(13).nnode=18
	Elttype(13).description="18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces)."
	Elttype(14).nnode=14
	Elttype(14).description="14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face)."
	Elttype(15).nnode=1
	Elttype(15).description="1-node point."
	
	Elttype(16).nnode=8
	Elttype(16).description="8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges)."
	et1=16
	call not_nodal_force_weight(et1)	
	
	
	Elttype(17).nnode=20
	Elttype(17).description="20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges)."
	
	Elttype(18).nnode=15
	Elttype(18).description="15-node second order prism (6 nodes associated with the vertices and 9 with the edges)."
	
	Elttype(19).nnode=13
	Elttype(19).description="13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges)."
	Elttype(20).nnode=9
	Elttype(20).description="9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)"
	Elttype(21).nnode=10
	Elttype(21).description="10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)"
	Elttype(22).nnode=12
	Elttype(22).description="12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)"
	
	Elttype(23).nnode=15
	Elttype(23).description="15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)"
	et1=23
	call not_nodal_force_weight(et1)	
	

	
	Elttype(24).nnode=15	
	Elttype(24).description="15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)"
	Elttype(25).nnode=21
	Elttype(25).description="21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)"
	Elttype(26).nnode=4
	Elttype(26).description="4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge)"
	
	Elttype(27).nnode=5
	Elttype(27).description="5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)"
	et1=27
	call not_nodal_force_weight(et1)	

	
	Elttype(28).nnode=6
	Elttype(28).description="6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge)"
	Elttype(29).nnode=20
	Elttype(29).description="20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)"
	Elttype(30).nnode=35
	Elttype(30).description="35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)"
	Elttype(31).nnode=56
	Elttype(31).description="56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)"
	Elttype(92).nnode=64
	Elttype(92).description="64-node third order hexahedron (8 nodes associated with the vertices, 24 with the edges, 24 with the faces, 8 in the volume)"
	Elttype(93).nnode=125
	Elttype(93).description="125-node fourth order hexahedron (8 nodes associated with the vertices, 36 with the edges, 54 with the faces, 27 in the volume)"

	
end subroutine

subroutine not_nodal_force_weight(et)
	use DS_Gmsh2Solver
	implicit none
	integer,intent(in)::et
	
	select case(et)
		case(1)
			allocate(Elttype(et).weight(Elttype(et).nnode,2))
			Elttype(et).weight(:,1)=0.5
			Elttype(et).weight(1,2)=1./6.
			Elttype(et).weight(2,2)=1./3.
		case(2,3,4,5,6)
			allocate(Elttype(et).weight(Elttype(et).nnode,1))
			Elttype(et).weight=1./Elttype(et).nnode
		case(8)
			allocate(Elttype(et).weight(Elttype(et).nnode,2))
			Elttype(et).weight(1:2,1)=1./6.
			Elttype(et).weight(3,1)=2./3.
			Elttype(et).weight(1,2)=0.
			Elttype(et).weight(2,2)=1./6.
			Elttype(et).weight(3,2)=1./3.
		case(9)
			allocate(Elttype(et).weight(Elttype(et).nnode,1))
			Elttype(et).weight=0.
			Elttype(et).weight(4:6,1)=1./3.
		case(11)	
			allocate(Elttype(et).weight(Elttype(et).nnode,1))
			Elttype(et).weight=0.
			Elttype(et).weight(1:4,1)=-1./20.
			Elttype(et).weight(5:10,1)=1./5.
		case(16)
			allocate(Elttype(et).weight(Elttype(et).nnode,1))
			Elttype(et).weight=0.
			Elttype(et).weight(1:4,1)=-1./12.
			Elttype(et).weight(5:8,1)=1./3.
		case(23)
			allocate(Elttype(et).weight(Elttype(et).nnode,1))
			Elttype(et).weight=0.D0
			Elttype(et).weight(4:6,1)=-1./45.
			Elttype(et).weight(7:12,1)=4./45.
			Elttype(et).weight(13:15,1)=8./45.
		case(27)
			allocate(Elttype(et).weight(Elttype(et).nnode,2))
			Elttype(et).weight(1,1)=7./90.
			Elttype(et).weight(2,1)=7./90.
			Elttype(et).weight(3,1)=16./45.
			Elttype(et).weight(5,1)=16./45.
			Elttype(et).weight(4,1)=2/15.
			Elttype(et).weight(1,2)=0
			Elttype(et).weight(2,2)=7./90.
			Elttype(et).weight(3,2)=4/45.
			Elttype(et).weight(5,2)=4./15.
			Elttype(et).weight(4,2)=1/15.
		end select
end subroutine

   !把字符串中相当的数字字符(包括浮点型)转化为对应的数字
   !如 '123'转为123,'14-10'转为14,13,12,11,10
   !string中转化后的数字以数组ar(n1)返回，其中,n1为字符串中数字的个数:(注　1-3转化后为3个数字：1,2,3)
   !nmax为数组ar的大小,string默认字符长度为512。
   !num_read为要读入数据的个数。
   !unit为文件号
   !每次只读入一个有效行（不以'/'开头的行）
   !每行后面以'/'开始的后面的字符是无效的。
   subroutine  strtoint(unit,ar,nmax,n1,num_read,set,maxset,nset)
	  implicit none
	  logical::tof1,tof2
	  integer::i,j,k,strl,ns,ne,n1,n2,n3,n4,step,nmax,& 
			num_read,unit,ef,n5,nsubs,maxset,nset
	  real(8)::ar(nmax),t1
		character(32)::set(maxset)
	  character(512)::string
	  character(32)::substring(100)
	  character(16)::legalC,SC

		LegalC='0123456789.-+eE*'
		sc=',; '//char(9)

		n1=0
		nset=0
		ar=0
		!set(1:maxset)=''
	  do while(.true.)
		 read(unit,'(a512)',iostat=ef) string
		 if(ef<0) then
			print *, 'file ended unexpected. sub strtoint()'
			stop
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

		 if(string(1:1)/='/') then
			
			!每行后面以'/'开始的后面的字符是无效的。
			if(index(string,'/')/=0) then
				strL=index(string,'/')-1
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
				if(tof1) then !处理类似于'1-5'这样的形式的读入数据
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
					if(tof2) then !处理类似于'1*5'(表示5个1)这样的形式的读入数据
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
		    exit
		 else
		    if(n1>num_read)  print *, 'error!nt2>num_read. i=',n1
		 end if
	
	  end do	

   end subroutine


!translate all the characters in term into lowcase character string
subroutine lowcase(term,iterm)
	use dflib
	implicit none
	integer i,in,nA,nZ,nc,nd,iterm
	character(1)::ch
	character(iterm)::term
	
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


   !n1为node()数组中全域节点总数
   subroutine reorder_nodal_number(IPERM,nnum,adjL,maxadj)

	  implicit none
	  integer,intent(in)::nnum,maxadj
	  integer,intent(in)::adjL(maxadj,nnum)
	  integer,intent(inout)::IPERM(nnum)
	  integer::i,j,k,nnz,ITYPE,IFLAG,LP,MP,n1,n2,n3
	  REAL(8)::IPROF(2)
	  integer,allocatable::IRN(:),JCN(:),IW(:),ICPTR(:),art(:)
	  integer::ar(15),allo_err
	  real(8)::t1
	  COMMON /MC40I/ LP,MP

 

!	  allocate(IPERM(nnum))
	  allocate(ICPTR(nnum+1))
	  allocate(IW(3*nnum+2))

 	  
	  !统计adjL()中非零元素的个数。
	  NNZ=count(adjL>0)
	  allocate(IRN(2*NNZ))
	  allocate(JCN(NNZ))
	
	  n2=0
	  do i=2,nnum
		 do j=1,maxadj
			if(adjL(j,i)/=0) then
			   n2=n2+1
			   IRN(n2)=i
			   jCN(n2)=adjL(j,i)
			else
			   exit
			end if
		 end do	
	  end do
	
	  ITYPE=1
	  ICPTR=0
	  IW=0
	  IPERM=0

	  CALL MC40AD(ITYPE,nnum,NNZ,IRN,JCN,ICPTR,IPERM,IW,IPROF,IFLAG)
	  ! Check for an error return
	  if (IFLAG.LT.0) then
		 print *, 'S.W. Sloan method failed. Try another normal method to reordering nodal number. Please Wait...'
		 !call nodecoding(nnum)
		 return
	  end if
!	  C Write out the profile
	  WRITE (6,210) IPROF(1)
	  IF (IPROF(1).EQ.IPROF(2))THEN
		 WRITE (6,200)
	  ELSE
		 WRITE (6,220) IPROF(2)
	  END IF
	  !pause

!
	  200 FORMAT(/5X,'The algorithm did not reduce the profile')
	  210 FORMAT(/5X,'The profile initially is',F15.0)
	  220 FORMAT(/5X,'The reduced profile is',F15.0)


	  !call csb()

	  !deallocate(adjl)
	  deallocate(IRN)
	  deallocate(JCN)	
	  deallocate(IW,stat=allo_err)
	  !deallocate(IPERM)
	  deallocate(ICPTR,stat=allo_err)

	  return	
   end subroutine
   
   
   	!对于每个节点，存储和其相邻(定义为总刚中，这两个节点对应的元素不为零。)且编号小于该节点自身编号的节点编号
	!这里假定一个节点最多有maxadj个与之相邻且编号小于该节点自身编号的节点。
	
   subroutine setup_adjList()
		use DS_Gmsh2Solver
		implicit none
		integer::i,j,iel,n1,ar(200),n2
		
		allocate(adjL(maxadj,nnode))
		adjL=0
		do iel=1,nel
	  	
	  		n1=element(iel).nnode
	  		ar(1:n1)=element(iel).node
	  		!从大到小
  			do i=1,n1
			  do j=i+1,n1
				 if(ar(i)<ar(j)) then
					n2=ar(i)
					ar(i)=ar(j)
					ar(j)=n2
				 end if
			  end do
			end do
			
			do i=1,n1-1
			  do j=i+1,n1
				 call addtoadjL(ar(j),adjL(:,ar(i)),maxadj)
			  end do
			end do					

		end do
   end subroutine

   !如果iar(:)中没有n1,则把n1加入到iar(:)中。
   subroutine addtoadjL(n1,iar,maxadj)
	  implicit none
	  integer,intent(in)::n1,maxadj
	  integer::iar(maxadj),i
	
		  
	  if(any(iar==n1)) return
	
	  if(any(iar==0)) then
		 do i=1,maxadj
		    if(iar(i)==0) then
			   iar(i)=n1
			   exit
			end if
		 end do
	  else
		 print *, '与节点相邻且编号小于该节点自身编号的节点数在>maxadj,adjL的原定空间不足.'
		 stop
	  end if

	
    end subroutine
    
logical function isacw(x1,y1,z1,x2,y2,z2,x3,y3,z3)
	implicit none
	real(8)::x1,y1,z1,x2,y2,z2,x3,y3,z3
    real(8)::t1
	
	isacw=.false.
	y2=y2-y1
	x2=x2-x1
    z2=z2-z1
	y3=y3-y1
	x3=x3-x1
    z3=z3-z1
	t1=(x2*y3-y2*x3)+(y2*z3-z2*y3)-(x2*z3-z2*x3)
	if(t1>0) isacw=.true.

end function
