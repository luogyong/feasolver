subroutine write_readme_feasolver()	
	use solverds
	use ifport
	implicit none
	integer::i,j,item
	LOGICAL(4)::tof,pressed
	!integer,external::ipp
	integer,parameter::nreadme=2048
	character(1024)::readme(nreadme)
	
    print *, "The help file is in d:\README_FEASOLVER.TXT."
	open(2,file=HELPFILE,STATUS='REPLACE')
	
	I=0
	README(IPP(I)) ="//THE KEYWORD STRUCTURE USED IN THE INPUT FILE SINP IS EXPLAINED HEREIN"
	README(IPP(I)) ="//THE [] MEANS OPTIONAL."
	README(IPP(I)) ="//THE | MEANS OR(或)."
	README(IPP(I)) ="//THE CHARACTER INSIDE () MEANS THE DATATYPE OR AN ARRAY."
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//TITLE"
	README(IPP(I))=  "//"//'"'//"THE KEYWORD TITLE IS USED TO INPUT THE INFOMATION OF THE MODEL."//'"'
	README(IPP(I))= "//{TITLE(A)} "
    
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//HELPFILE,EXIT=YES|NO"
	README(IPP(I))=  "//"//'"'//"THE KEYWORD TITLE IS USED TO OUTPUT A IFSOLVER HELPFILE IN THE CURRENT DIRECTION."//'"'
	README(IPP(I))= "//NO PARAMETER IS NEEDED.IS EIXT=YES, WHEN THE OUTPUT IS COMPLETED, STOP RUNNING AND EXIT."    
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//NODE,NUM=...(I) [,DATAPACKING=1|2] [,DEMENSION=1|2|3]   // NUM=节点数"
	README(IPP(I))=  "//"//'"'//"THE KEYWORD NODE IS USED TO INPUT THE NODAL INFOMATION."//'"'
	README(IPP(I)) = "//{X(R),Y(R)[,Z(R)])} //DATAPAKING=1.|"
	README(IPP(I)) = "//{X1(R),...,XNUM(R),Y1(R),...,YNUM(R)[,Z1(R),...,ZNUM(R),])} //DATAPAKING=2 "
	README(IPP(I)) =  "//{......}   //当DATAPAKING=1,共NUM行; 当DATAPAKING=2,共DEMENSION*NUM个数 "

    
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//ELEMENT,NUM=...(I),ET=...(A),MATID=...(I)[,SET=...(I),SYSTEM=1]  // NUM=单元个数,ET=单元类型,MATID=材料号,SET=单元集号,单元局部坐标。此关键词可重复出现。"
	README(IPP(I)) = "                                                                  //当ET=BEAM时,SYSTEM=局部坐标，为梁单元的y轴定向; | "
    README(IPP(I))=  "//"//'"'//"THE KEYWORD ELEMENT IS USED TO INPUT THE ELEMENT INFOMATION."//'"'
	README(IPP(I)) = "//{N1(I),N2(I),...,NN(I)} //N为单元节点号; |"
	README(IPP(I)) = "//{N1(I),N2(I),...,NN(I),H} //当ET=DKT3,SHELL3,SHELL3_KJB时, H=单元的厚度"
	README(IPP(I)) =  "//{......}   //共NUM行. "

	
	Do j=1,2
		README(IPP(I)) ="\N//******************************************************************************************************"C
		if(j==1) README(IPP(I)) = "//BC,NUM=(I)[,SF=(I),SSP_ONEPILE=(I),SPG_ISDUAL=(I),ISINC=(I)]   //NUM=节点约束个数,SF=时间因子，SSP_ONEPILE=是否只作用在一根钢板桩上(0,NO).此关键词可重复出现。ISINC=是否为增量(默认为全量)."  
		if(j==2) README(IPP(I)) = "//LOAD,NUM=(I)[,SF=(I),SSP_ONEPILE=(I),SPG_ISDUAL=(I),ISINC=(I)]   //NUM=节点荷载个数,SF=时间因子，SSP_ONEPILE=是否只作用在一根钢板桩上(0,NO).此关键词可重复出现。ISINC=是否为增量(默认为全量)."  
		README(IPP(I))=  "//"//'"'//"THE KEYWORD BC IS USED TO INPUT NODAL CONSTRAINS."//'"'
		README(IPP(I)) = "//{NODE(I),DOF(I),VALUE(R)[,SF(I),SPG_ISDUAL(I),SSP_ONEPILE(I),ISINC(I)]}  //ISDUAL:如果isdual==i(>0),则表示此自由度可能与出溢边界Nseep(i)重复，如果边界水头小于位置水头，则变为出溢边界。"  
		README(IPP(I)) = "//{...}  //共NUM行"	
		README(IPP(I)) = "//注意，对于水头边界(dof=4)，如果边界值小于其相应的位置水头，则认为该边界无效，不起作用。"
		README(IPP(I)) = "//注意，多步计算时，输入时要求每一步的水头均是总量，而不是增量."	
	end do
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//SEEPAGE FACE,NUM=...(I)，STEP FUNCTION=...(I)   //NUM=出溢点的个数,STEP FUNCTION=出溢点的时间步函数，默认为0"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD SEEPAGE FACE IS USED TO INPUT THE NODES CONSTRAINED BY A SEEPAGE FACE CONDITION."//'"'
	README(IPP(I)) = "//{N1(I),N2(I),...,NN(I)]} //节点号，共NUM个。此行可重复出现，直至节点数=NUM." 
	README(IPP(I)) = "//此关键词可重复出现，以输入不同STEP FUNCTION的出溢点."
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//BODY FORCE,NUM=...(I)  //NUM=单元荷载个数"
	README(IPP(I))=  "//"//'"'//"THE KEYWORD BODY FORCE IS USED TO INPUT ELEMENTAL LOADS."//'"'
	README(IPP(I)) = "//{E1,E2,...,EN,DOF,VALUE(R),STEPFUNC.}  //Ei为单元编号或单元组名。此行可重复出现，直至单元个数=NUM."  
	README(IPP(I)) = "// 当Ei为单元组名时，该单元组前面必须用关键词UESET定过."
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//ELT_LOAD,NUM=...(I)  //NUM=单元荷载行数 //此关键词可取代BODY FORCE. 此关键词目前不能使用，需进一步完善."
	README(IPP(I))=  "//"//'"'//"THE KEYWORD ELT_LOAD IS USED TO INPUT LOADS APPLIED ON AN ELEMENT GROUP."//'"'
	README(IPP(I)) = "//{ELTGROUP(I),DOF(I),VALUE(I),STEPFUNC.(I)}  //ELTGROUP=I, 即关键词ELEMENT中定义的SET号.亦即为ESET(I)的单元."  
	README(IPP(I)) = "//{......} //共NUM行. "
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//DATAPOINT,NUM=...(I)  //NUM=输出数据集的个数."
	README(IPP(I))=  "//"//'"'//"THE KEYWORD DATAPOINT IS USED TO OUTPUT VALUES IN SPECIFIC NODES ."//'"'
	README(IPP(I)) = "//{A:NNODE,ISSUMQ}  //NNODE=数据集的节点数，ISSUMQ/=0,表示仅输出此点集的各节点的流量和"  
	README(IPP(I)) = "//{A1:iNodes} //共NNODE节点. "
	README(IPP(I)) = "//{...} //共{A,A1}*NUM组. "

	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//TIME STEP,NUM=...(I)  //NUM=个数."
	README(IPP(I))=  "//"//'"'//"THE KEYWORD TIME STEP IS USED TO INPUT THE SUBINCREMENTAL TIME FOR EACH STEP IN A TRANSIENT ANALYSIS."//'"'
	README(IPP(I)) = "//{ISTEP(I),NSUBTIMESTEP(I),TIMESTEP(1)(R),TIMESTEP(2)(R),...,TIMESTEP(STEP)(R)}  //共NSUBTIMESTEP个时间,总的分析时间步长为各子步长之和."  
	README(IPP(I)) = "//{......} //共NUM行. "
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//INITIAL VALUE,NUM=...(I)   //NUM=初值个数"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD INITIAL VALUE IS USED TO INITIAL VALUE DATA."//'"'
	README(IPP(I)) = "//{NODE(I),DOF(I),VALUE(R)[,STEPFUNC.(I)]}"  
	README(IPP(I)) = "//{...}  //共NUM行"	
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//STEPINFO,NUM=...(I)   //NUM=步信息个数"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD STEPINFO IS USED TO STEP INFOMATION DATA."//'"'
	README(IPP(I)) = "//{ISTEP(I),MATHERSTEP(I),ISSTEADY(I),ISSTEADY(I)[,BCTYPE,LOADTYPE]} //ISTEP:STEP NUMBER; MATHERSTEP:此步的依托步; ISSTEADY:是否稳态分析, 1,Yes,others NO."  
	README(IPP(I)) = "//LOADTYPE(BCTYPE):为荷载(位移边界)施加方式，=1,表示步内荷载(位移边界)随时间线性施加，=2(default)，表示步荷载(位移边界)在步初瞬间施加。"
    README(IPP(I)) = "//{...}  //共NUM行"	
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//STEP FUNCTION,NUM=...(I),STEP=...(I)[,BASE=1]   //NUM=步方程的个数,STEP=步数.,BASE=0，请同时输入地应力步(第0步)的FACTOR，这是要输入STEP+1个factors"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD STEP FUNCTION IS USED TO STEP FUNCTION DATA."//'"'
	README(IPP(I)) = "//{[FACTOR(0)(R)],FACTOR(1)(R),FACTOR(2)(R),...,FACTOR(STEP)(R),TITLE(A)& 
                         \N// FACTOR(ISTEP)=第ISTEP步边界或荷载的系数. &
                         \N//当FACTOR(ISTEP)=-999时,此边界或荷载在此步中失效（无作用）. &
                         \N//!表单元生死时，0为死1为生。"C  
	README(IPP(I)) = "//{...}  //共NUM行"
	!README(IPP(I)) = "//注意,对于力学模型（SOLVER_CONTROL.TYPE=SLD），输入的是各步荷载或位移边界的增量;而对渗流模型（SOLVER_CONTROL.TYPE=SPG）,输入的是各步水头或流量边界的总量."

	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//SLAVE_MASTER_NODE_PAIR,NUM=...(I)   //NUM=约束节点对的个数."  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD SMNP IS USED TO SLAVE-MASTER NODE PAIR DATA."//'"'
	README(IPP(I)) = "//{SLAVE,SDOF,MASTER,MDOF} //SLAVE节点号，SLAVE自由度，MASTER节点号，MASTER自由度， "  
	README(IPP(I)) = "//{...}  //共NUM行.!注意，master点上的mdof自由度上目前不能施加位移边界条件，这时输入时请互换slave和master."
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//WSP,NUM=...(I)[CP=...(I),Method=...(I)]   //CP=1,仅为计算水面算,水面线计算完成后就退出。CP=2,渗流计算是边界按不考虑水跃作用取值。Mothed=1,2,3,4(kinds,WangXG,Ohtsu,Chow)"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD IS USED TO INPUT WATERSURFACEPROFILE CALCULATION DATA."//'"'
	README(IPP(I)) = "//{NSEGMENT(I),NNODE(I),Q(R),B(R),UYBC,DYBC,[Q_SF,UYBC_SF,DYBC_SF,caltype,g,kn]}" 
	README(IPP(I)) = "//流道数;流道剖分总节点数;过流量;流道宽;上游边界水深;下游边界水深; [Q的步长函数,UYBC的步长函数,DYBC的步长函数,计算控制参数;重力加速度;曼宁量纲系数]"
	README(IPP(I)) = "//{UNODE(I),DNODE(I),n(R),So(R),[Profileshape(1),Profileshape(2)]}  //起点号（局部编号，最上游节点编号为1，最下游编号为NNODE),终点号,糙率,坡率,[急流形式,缓流形式]。(流道参数行，从上游到下游依次输入,共NSEGMENT行)"
    README(IPP(I)) = "//{N1,N2,...,N(NNODE)}  //流道内节点（全局）编号,共NNODE个，从上游往下游输入。"
	README(IPP(I)) = "//注：Profileshape=M1,M2,M3,S1,S2,S3,H2,H3,A2,A3,C1,C3 "
	README(IPP(I)) = "//注：UYBC（DYBC）=-1,表示上游（下游）边界水深为临界水深，=-2，表示上游（下游）边界水深为正常水深。"
    README(IPP(I)) = "//注：CALTYPE=1，只算急流；=2，只算缓流；=0（默认），两者都算."
    README(IPP(I)) = "//注：对于/=4(Chow)，目前只能处理NSEGMENT=2的情况，即斜坡段+水平段，且上游在左下游在右.默认水平段与斜坡段的交点为第二段的第一个节点。"
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//SOILPROFILE,NUM=(I),[spmethod=(I),kmethod=(I),rf_epp=(I),rf_app=(I),iniepp=(I),soilspringmodel=0]   \n//spmethod=土压力计算方法，0，郎肯； &
        \n// kmethod=基床系数的计算方法，0，m法；1，(E,V)法；2,zhu;3 biot;4 vesic;-1,按直接输入. &
        \n// rf_epp=被动侧土弹簧抗力限值是否要减掉初始的主动土压力.(0N1Y). &
        \n// RF_APP=0 !主动侧主动土压力荷载，开挖面以下是否按倒三角折减.(0N1Y). & 
		\n// iniepp=2 !被动侧土弹簧土压力初始值(位移=0时的土压力),2=主动土压力，1=静止土压力. &
		\n// soilspringmodel !土弹簧模型：=HYPERBOLIC,为双曲模型;=EIP_SPRING,弹理想塑性(默认);=ELASTIC,为弹性"C
		
	README(IPP(I))=  "//"//'"'//"THE KEYWORD SOILPROFILE IS USED TO INPUT SOILPROFILE DATA."//'"'C
	README(IPP(I))=  "//A0:{TITLE(C)}  //土层剖面的名字"  
	README(IPP(I)) = "//A:{NASOIL(I),NPSOIL(I),BEAMID(I)}  //主动侧土层数(负数表主动土压力为负，被动为正，反之亦然。)，被动侧土层数，地基梁号." 
	README(IPP(I)) = "//B:{(Z1,Z2,MAT,WPMETHOD,STEPFUN [,PV])*NASOIL}   //层顶高程点号，层底高程点号，材料号，水压力考虑方法(0=合算，1=常规分算，2=分算，考虑渗透力,3=手动输入),步函数 [,超载]。共NASOIL行"  
	README(IPP(I)) = "//C:{(Z1,Z2,MAT,WPMETHOD,STEPFUN [,PV])*NPSOIL}   //层顶高程点号，层底高程点号，材料号，水压力考虑方法(0=合算，1=常规分算，2=分算，考虑渗透力,3=手动输入),步函数 [,超载]。共NPSOIL行"  
	README(IPP(I)) = "//D:{AWATERLEVEL,ASTEPFUN,PWATERLEVEL,PSTEPFUN}   //主动侧水位,主动侧水位步函数，被动侧水位,被动侧水位步函数，"  
	README(IPP(I)) = "//E:{ALoad,ALoadSTEPFUN,PLoad,PLoadSTEPFUN}   //主动侧超载,主动侧超载步函数，被动侧超载,被动侧超载步函数，" 
!	README(IPP(I)) = "//F:{NO_ACTION(I)}*NACTION   //约束号，共NACTION个" 
!	README(IPP(I)) = "//G:{NO_STRUT(I)}*NSTRUT   //支撑号，共NSTRUT个" 
	README(IPP(I)) = "//{A0,A,B,C,D,E}*NUM。   //共NUM组"
	README(IPP(I)) = "//注意：\n//1)每一时间步，土层按顺序从上而下输入，不同时间步间的土层可以重叠,同一时间步有效各土层不要出现重叠和空隙。 &
						\n//2)各时间步地下水位处要分层;桩顶(底)处要分层(如果桩在土里面);桩的材料分界处要分层。 &
						\n//3)考虑渗透力时，假定awL>pwL. &
						\n//4)如果水面高于地表，将水等效为为土层（令c,phi,模量均设为0，渗透系数<=0)"C
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//PILE,NUM=...(I)"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD PILE IS USED TO INPUT BEAM(RetainingStructure) DATA."//'"'
	README(IPP(I)) = "//A:{NSEG(I),[SYSTEM=0]}  //材料分段数，坐标号" 
	README(IPP(I)) = "//B:{Z(1:NSEG+1)}   //材料分段点,应从上往下（或从左往右）输入（方便单元寻址）"  
	README(IPP(I)) = "//C:{MAT(1:NSEG)}   //各段材料号"  
	README(IPP(I)) = "//{A,B,C}*NUM。   //共NUM组"	

	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//STRUT,NUM=...(I),[ISBAR=YES|NO]   //如果模拟两边同时开挖，则支撑用杆单元进行模拟，此时令isbar=yes,反之，isbar=no（默认）"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD STRUT IS USED TO INPUT STRUT(RetainingStructure) DATA."//'"'
	README(IPP(I)) = "//A:{Z,MAT,STEPFUN}  //点号，材料号，步函数" 
	README(IPP(I)) = "//{A}*NUM。   //共NUM组"
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//KPOINT,NUM=...(I)"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD KPOINT IS USED TO INPUT KeyPoint DATA."//'"'
	README(IPP(I)) = "//A:{NO(I),XY(1:NDIMENSION),[ELEMENTSIZE(R)]}  //点号，XY(1:NDIMENSION),[ELEMENT SIZE]" 
 	README(IPP(I)) = "//{A}*NUM。   //共NUM组"	
    
 	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//GEOLINE,NUM=...(I)"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD GeoLINE IS USED TO INPUT GEOMETRY LINE DATA."//'"'
	README(IPP(I)) = "//A:{NO(I),MATID(I),NPOINT(I),IPOINT(I),[TITLE(A)]}  //线号，材料号(投影)，控制点个数，控制点号(共NPOINT个)，名字" 
 	README(IPP(I)) = "//{A}*NUM。   //共NUM组"   
	README(IPP(I)) = "//注意：1）按X从左至右的顺序输入,既x1<=x2<=..<=xn" 
	
 	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//WATERLEVEL,NUM=...(I),FMT=0|1,VAR=1|2|3  //VAR=INTERPOLATION INDEPENDENT VARIABLE. 1/2/3 MEANS X/Y/Z"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD IS USED TO INPUT WALTERLEVEL LINE DATA."//'"'
    README(IPP(I)) = "//IF fmt=0(DEFAULT) then input the id number in the array KPOINT of the points. " 
    README(IPP(I)) = "//IF fmt=1 then INPUT (VAR,H) of the points. " 
	README(IPP(I)) = "//A:{P1,P2,...,P(NUM)]}  //(IF FMT=0)控制点号(共NUM个)" 
    README(IPP(I)) = "//A:{XI,HI}*NUM  //(IF FMT=1)(共NUM个)" 
	README(IPP(I)) = "//注意：1）按X从左至右的顺序输入,既x1<=x2<=..<=xn" 	
	
 	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//RIGHT TURN POINT,NUM=...(I)"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD  IS USED TO INPUT RIGHT TURN POINT DATA."//'"'
	README(IPP(I)) = "//A:{IPOINT1,IPOINT2,...IPOINTN}  //点号(共NUM个)" 
 	README(IPP(I)) = "//{A}*1。   "   
	README(IPP(I)) = "//此命令用于输入竖直坡面的下端点的点号(既同一X处有两个地表高程Y)，以解决同一X点有两个地表高程的问题。" 	
    
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//ACTION,NUM=...(I)"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD ACTION IS USED TO INPUT Laction DATA..."//'"'
	README(IPP(I))=  "//A0:{TITLE(C)}  //作用的名字"  
	README(IPP(I)) = "//A:{NKP,TYPE,DOF,NDIM,[SF,ISVALUESTEPFUN,ISEXVALUE]}  //控制点数，作用类型（0=力，1=位移，2=刚度），作用的自由度，作用的维度，生死步函数,值步函数开关(0N1Y)，极值开关(0N1Y)" 
	README(IPP(I)) = "//B:{KPOINT(1:NKP)}  //控制点号,应从上往下（或从左往右）输入（方便单元寻址）" 
	README(IPP(I)) = "//C:{VALUE(1:NKP)}  //控制点上作用的大小。注意单位统一(线作用的宽度取决于支护桩的间距)，当ndim=1,type=2时，单位：F/L^3;(ndim=1,type=1,L); (ndim=1,type=0,F/L^2) "
	README(IPP(I)) = "//D:{SF(1:NKP)}  //控制点上作用的步函数，如果isstepfun=1."
	README(IPP(I)) = "//E:{EXVALUE(1:NKP,1:2)}  //控制点上作用的上下限值(先各点下限(exvalue(:,1))，后各点上限(exvalue(:,2))) 如果isexvalue=1."
 	README(IPP(I)) = "//{A0,A,B,C,D,E}*NUM。   //共NUM组"
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//HINGE,NUM=...(I) //(此功能目前仅对beam2d及beam单元有效)"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD HINGE IS USED TO INPUT HINGE/FREEDOF DATA..."//'"'
	README(IPP(I)) = "//A:{ELEMENT,NODE,DOF}  //(要释放的自由度所在的)单元，节点和自由度编号。" 
 	README(IPP(I)) = "//{A}*NUM。   //共NUM组"
    
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//SLOPEPARAMETER,SLOPEMETHOD=0|1|2|3|4|5,OPTIMIZATIONMETHOD=1,SLIPSHAPE=CIRCULAR|NONCIRCULAR,SLICEWIDTH=1.0,..." 
	README(IPP(I))=  "//"//'"'//"THE KEYWORD HINGE IS USED TO INPUT HINGE/FREEDOF DATA..."//'"'
	README(IPP(I)) = "//SLOPEMETHOD:边坡分析方法，1，ordinary,2,bishop,3,spencer,4,janbu,5,gle; 0 for all." 
 	README(IPP(I)) = "//OMTIMIZATION,搜寻方法，1，grid;2,MONTE CARLO"	
    README(IPP(I)) = "//SLIPSHAPE,滑弧形状，1,circular;0,noncircular;"
    README(IPP(I)) = "//SLICEWIDTH,土条宽度,(1.0,bydefault)"
    README(IPP(I)) = "//XMIN_MC,XMAX_MC,SEARCH RANGE FOR MONTE CARLO TECHNIQUE"
    
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//MATERIAL,MATID=...(I),[TYPE=...(I)],[ISFF=YES|NO],[NAME=...(c)],ISSF=...(I)//MATID=材料号，TYPE=材料类型，ISFF=是否为依赖某场变量,NAME=材料文字注释.此关键词可重复出现.ISSF=是否输入材料参数的步函数(0N1Y)" 
	README(IPP(I))=  "//"//'"'//"THE KEYWORD MATERIAL IS USED TO INPUT MATERIAL INFOMATION."//'"'
	README(IPP(I)) = "//{A: MATERIAL.PROPERTIES}  //PLEASE REFER TO THE EXPLANATION IN THE END OF THE FILE."
	README(IPP(I)) = "//[B:{FIELD FUNCTION PARAMETERS}]  // 如ISFF=YES "	
	README(IPP(I)) = "//[C:{STEPFUNCTION(1:NPARAMETERS)}]  // 如ISSF=YES "	
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) ="//MATERIAL.PROPERTY (FYI)" 
	README(IPP(I)) ="//1. CASE (CONDUCT): PROPERTY(1)=K0 .(2)=K1  "
	README(IPP(I)) ="//THE CONDUCTIVITY K=K0+K1*FIELD VARIABLE\n"c
	README(IPP(I)) ="//2. CASE (LA_MC)  : PROPERTY(1)=P  .(2)=PHI .(3)=C	.(4)=GANMA	.(5)=Pore Pressure"
	README(IPP(I)) ="//GANMA=SPECIFIED GRAVITY OF SOIL\N"c
	README(IPP(I)) ="//3. CASE(ELASTIC) : PROPERTY(1)=E  .(2)=V\n"C
	README(IPP(I)) ="//4. CASE(MISES)   : PROPERTY(1)=E  .(2)=V   .(3)=SIGMA_Y  .(4)=GAMA"
	README(IPP(I)) ="//SIGMA_Y=YIELD STRESS UNDER UNIAXIAL STRETCH. FOR CLAY UNDER PLAIN STRAIN SIGMA_Y=3*CU (CU, UNDRAINED SHEAR STRENGTH) AND SIGMA_Y=2*CU FOR TRIAXIAL (AXISYMMETRICAL) STATE.\N"C
	README(IPP(I)) ="//5. CASE(MC)      : PROPERTY(1)=E  .(2)=V   .(3)=C  .(4)=PHI(DEG)	.(5)=DILATION ANGLE(DEG) .(6)=weight, .(7)=a .(8)=SitaT\N   "C
	README(IPP(I)) ="//a>=0 and SitaT are roundoff parameters,the default values are a=0.05*c*cot(phi),sitaT=25o.IF a=0. and sitaT=30. no rounded and Griffiths algorithm is used. \N "C
	README(IPP(I)) ="//Referred to: Abbo A J, Sloan S W. A Smooth Hyperbolic Approximation to the Mohr-Coulomb Yield Criterion[J]. Computers and Structures, 1995, 54(3): 427-441. )\N   "C
	README(IPP(I)) ="//6. CASE(CAMCLAY) : PROPERTY(1)=M  .(2)=V   .(3)=LAMDA \N  "C
	README(IPP(I)) ="//7. CASE(SPG)     : PROPERTY(1)=K1 .(2)=K2  .(3)=K3       .(4)=INT(TRANSM)    .(7)=ALPHA .(8)=N  .(9)=M .(10)=Mv	.(11)=Sita_s	.(12)=Sita_r    .(13)=rw 	.(14)=THICKNESS"
	README(IPP(I)) ="//TRANSM=L2G,FOR SPG, TRANSFORM THE MATERIAL KIJ TO GLOBLE KXY."
    README(IPP(I)) ="//ALPHA= A CURVE FITTING PARAMETER FOR THE VAN GENUCHTEN MODEL,LEONG AND RAHARDJO MODEL AND EXPONENT MODEL,UNIT=1/L." 
    README(IPP(I)) ="//N,M=CURVE FITTING DIMENSIONLESS PARAMETERS FOR THE VAN GENUCHTEN MODEL/LEONG AND RAHARDJO MODEL." 
    README(IPP(I)) ="//FOR THE VAN GENUCHTEN MODEL, THE DEFAULT M=1-1/N, AND M IS OMITTED FROM INPUT."
    README(IPP(I)) ="//Mv=COEFFICIENT OF COMPRESSIBILITY.(UNIT:L**2/F)"
    README(IPP(I)) ="//Sita_s=Saturated volumetric water content."
	README(IPP(I)) ="//Sita_r=Residual volumetric water content."
    README(IPP(I)) ="//rw=bulk Gravity of water.(UNIT:F/L**3)"C
	README(IPP(I)) ="//THICKNESS=PHYSICAL THICKNESS FOR A ZEROTHICKNESS ELEMENT\N"C
    README(IPP(I)) ="//注意单位的统一。尤其注意ALPHA的单位，在VG模型中为：1/L，在LR模型中为:L.\N"C
	README(IPP(I)) ="//8. CASE(BAR,BAR2D)     : PROPERTY(1)=E  .(2)=A	[.(3)=hy	.(4)=hz	 .(5)=minN(最大轴向压力)	.(6)=MaxN(最大轴向拉力)]"
	README(IPP(I)) ="//9. CASE(BEAM,BEAM2D,SSP2D)   : PROPERTY(1)=E  .(2)=A   .(3)=v	.(4)=J	.(5)=Iz	.(6)=Iy	[.(7)=hy	.(8)=hz]"
    README(IPP(I)) ="//                          [PROPERTY(9)=MinN .(10)=MaxN .(11)=minMx .(12)=maxMx .(13)=minMy .(14)=maxMy .(15)=minMz .(16)=maxMz]"
    README(IPP(I)) ="//                          [.(17)=yc ]"
	README(IPP(I)) ="//J,IY,IZ ARE REFERRED TO THE LOCAL COORDINATE SYSTEM."
	README(IPP(I)) ="//当为单元SSP2D指定材料，输入参数(A,Iz,minN,maxN,minMz,maxMz)均为单根钢板桩的参数,YC钢板桩形心轴距离锁口内沿的距离。"
	README(IPP(I)) ="//hy,hz分别为梁截面在y'和z'轴方向(局部坐标)的高度，为后处理转化为六面体单元时所用. \N"C
	README(IPP(I)) ="//当为支护桩时，hy,hz分别为桩径和桩距，计算作用于支护桩上的土压力和土弹簧的刚度用 \N"C
	README(IPP(I)) ="//10. CASE(PIPE2,PPIPE2) : PROPERTY(1)=R(管半径,L)  .(2)=LAMDA(管壁摩阻系数)	.(3)=EPSLON(管壁的绝对粗糙度,L)	.(4)=V(运动粘滞系数L**2/T)"C
	README(IPP(I)) ="//11. CASE(ExcavationSoil/SLOPESOIL) : PROPERTY(1:10)=黏聚力，摩擦角，天然/饱和重度，变形模量，泊松比,渗透系数，水平基床系数(F/L**3),墙土间摩擦角（度）,[层顶水压力PW1,层底水压力PW2]"C
	README(IPP(I)) ="//PW1，PW2仅当土层的水压力为手动输入时有用（WPMETHOD=3）。"
	README(IPP(I)) ="//12. CASE(spirng) : PROPERTY(1:3)=k,minV(发生负位移),maxV(发生正位移)，预加力，预加位移"C
	
						
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) ="//Coordinate system(FYI)" 
	README(IPP(I)) ="//对于2D问题,重力方向为Y轴方向，且向上为正,X轴向右为正."
	README(IPP(I)) ="//对于3D问题,重力方向为Z轴方向，且向上为正,水平面为XY平面,XYZ满足右手螺旋法则"
	README(IPP(I)) ="//弹簧类单元的局部坐标系与整体坐标系相同。"C
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) ="//DOF (FYI)" 
	README(IPP(I)) ="//Except for axisymmetric elements, the degrees of freedom are always referred to as follows:" 
	README(IPP(I)) ="//1 x-displacement /limit analysis vx" 
	README(IPP(I)) ="//2 y-displacement /limit analysis vy" 
	README(IPP(I)) ="//3 z-displacement " 
	README(IPP(I)) ="//4 Pore pressure, hydrostatic fluid pressure, scale field variable" 
	README(IPP(I)) ="//5 Rotation about the x-axis, in radians" 
	README(IPP(I)) ="//6 Rotation about the y-axis, in radians" 
	README(IPP(I)) ="//7 Rotation about the z-axis, in radians" 
	README(IPP(I)) ="//Here the x-, y-, and z-directions coincide with the global X-, Y-, and Z-directions, respectively." 
	
	do j=1,i
		item=len_trim(readme(j))
		write(2,20) readme(j)
	end do
	
	tof=system("D:\README_FEASOLVER.TXT")	
	
20	format(a<item>)
	
	contains
	
	integer function ipp(i) 
		implicit none
		integer,intent(inout)::i
		
		i=i+1
		ipp=i
		if(ipp>nreadme) then
			print *, "The length of README is", nreadme
			stop
		end if
				
	end function
	
end subroutine
