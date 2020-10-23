subroutine write_readme_feasolver()	
	use solverds
	use ifport
	implicit none
	integer::i,j,item
	LOGICAL(4)::tof,pressed
	!integer,external::ipp
	integer,parameter::nreadme=2048
	character(1024)::readme(nreadme)
	
    !print *, "The help file is in d:\README_FEASOLVER.TXT."
	open(20,file=HELPFILE,STATUS='REPLACE')
	
	I=0
	README(IPP(I)) ="//THE KEYWORD STRUCTURE USED IN THE INPUT FILE SINP IS EXPLAINED HEREIN"
	README(IPP(I)) ="//THE [] MEANS OPTIONAL."
	README(IPP(I)) ="//THE | MEANS OR(��)."
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
	README(IPP(I)) = "//SETTLEMENT_HEAD,NNODE=...(I),NHEAD=...(I)   //NNODE=ģ�͵Ľڵ�����,NHEAD=���ڵ������ˮͷ����"
	README(IPP(I))=  "//"//'"'//"THE KEYWORD SETTLEMENT_HEAD IS USED TO INPUT HEAD1 AND HEAD2 FOR SETTLEMENT INDUCED BY DOWNWATER."//'"'
    README(IPP(I))="//{SF(1),...,SF(NHEAD)}    //��ˮͷ��Ӧ��������ע�⣬��I����ѹ��ˮͷΪ����ĸ�ѹ��ˮͷ֮�ͣ���SUM(HEAD(I)*SF(I)-�߳�)�������Ǹ�ѹ�����" 
	README(IPP(I))="//{INODE,HEAD(1),...,HEAD(NHEAD)}*NNODE    //�ڵ�ţ�ˮͷ����Ӧ������"    	
    
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//NODE,NUM=...(I) [,DATAPACKING=1|2] [,DEMENSION=1|2|3]   // NUM=�ڵ���"
	README(IPP(I))=  "//"//'"'//"THE KEYWORD NODE IS USED TO INPUT THE NODAL INFOMATION."//'"'
	README(IPP(I)) = "//{X(R),Y(R)[,Z(R)])} //DATAPAKING=1.|"
	README(IPP(I)) = "//{X1(R),...,XNUM(R),Y1(R),...,YNUM(R)[,Z1(R),...,ZNUM(R),])} //DATAPAKING=2 "
	README(IPP(I)) =  "//{......}   //��DATAPAKING=1,��NUM��; ��DATAPAKING=2,��DEMENSION*NUM���� "

    
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//ELEMENT,NUM=...(I),ET=...(A),MATID=...(I)[,SET=...(I),SYSTEM=1]  // NUM=��Ԫ����,ET=��Ԫ����,MATID=���Ϻ�,SET=��Ԫ����,��Ԫ�ֲ����ꡣ�˹ؼ��ʿ��ظ����֡�"
	README(IPP(I)) = "                                                                  //��ET=BEAMʱ,SYSTEM=�ֲ����꣬Ϊ����Ԫ��y�ᶨ��; | "
    README(IPP(I))=  "//"//'"'//"THE KEYWORD ELEMENT IS USED TO INPUT THE ELEMENT INFOMATION."//'"'
	README(IPP(I)) = "//{N1(I),N2(I),...,NN(I)} //NΪ��Ԫ�ڵ��; |"
	README(IPP(I)) = "//{N1(I),N2(I),...,NN(I),H} //��ET=DKT3,SHELL3,SHELL3_KJBʱ, H=��Ԫ�ĺ��"
	README(IPP(I)) =  "//{......}   //��NUM��. "

	
	Do j=1,2
		README(IPP(I)) ="\N//******************************************************************************************************"C
		if(j==1) README(IPP(I)) = "//BC,NUM=(I)[,SF=(I),SSP_ONEPILE=(I),SPG_ISDUAL=(I),ISINC=(I)]   //NUM=�ڵ�Լ������,SF=ʱ�����ӣ�SSP_ONEPILE=�Ƿ�ֻ������һ���ְ�׮��(0,NO).�˹ؼ��ʿ��ظ����֡�ISINC=�Ƿ�Ϊ����(Ĭ��Ϊȫ��)."  
		if(j==2) README(IPP(I)) = "//LOAD,NUM=(I)[,SF=(I),SSP_ONEPILE=(I),SPG_ISDUAL=(I),ISINC=(I)]   //NUM=�ڵ���ظ���,SF=ʱ�����ӣ�SSP_ONEPILE=�Ƿ�ֻ������һ���ְ�׮��(0,NO).�˹ؼ��ʿ��ظ����֡�ISINC=�Ƿ�Ϊ����(Ĭ��Ϊȫ��)."  
		README(IPP(I))=  "//"//'"'//"THE KEYWORD BC IS USED TO INPUT NODAL CONSTRAINS."//'"'
		README(IPP(I)) = "//{NODE(I),DOF(I),VALUE(R)[,SF(I),SPG_ISDUAL(I),SSP_ONEPILE(I),ISINC(I)]}  //ISDUAL:���isdual==i(>0),���ʾ�����ɶȿ��������߽�Nseep(i)�ظ�������߽�ˮͷС��λ��ˮͷ�����Ϊ����߽硣"  
		README(IPP(I)) = "//{...}  //��NUM��"	
		README(IPP(I)) = "//ע�⣬����ˮͷ�߽�(dof=4)������߽�ֵС������Ӧ��λ��ˮͷ������Ϊ�ñ߽���Ч���������á�"
		README(IPP(I)) = "//ע�⣬�ಽ����ʱ������ʱҪ��ÿһ����ˮͷ��������������������."	
	end do
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//SEEPAGE FACE,NUM=...(I)��STEP FUNCTION=...(I)   //NUM=�����ĸ���,STEP FUNCTION=������ʱ�䲽������Ĭ��Ϊ0"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD SEEPAGE FACE IS USED TO INPUT THE NODES CONSTRAINED BY A SEEPAGE FACE CONDITION."//'"'
	README(IPP(I)) = "//{N1(I),N2(I),...,NN(I)]} //�ڵ�ţ���NUM�������п��ظ����֣�ֱ���ڵ���=NUM." 
	README(IPP(I)) = "//�˹ؼ��ʿ��ظ����֣������벻ͬSTEP FUNCTION�ĳ����."
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//BODY FORCE,NUM=...(I)  //NUM=��Ԫ���ظ���"
	README(IPP(I))=  "//"//'"'//"THE KEYWORD BODY FORCE IS USED TO INPUT ELEMENTAL LOADS."//'"'
	README(IPP(I)) = "//{E1,E2,...,EN,DOF,VALUE(R),STEPFUNC.}  //EiΪ��Ԫ��Ż�Ԫ���������п��ظ����֣�ֱ����Ԫ����=NUM."  
	README(IPP(I)) = "// ��EiΪ��Ԫ����ʱ���õ�Ԫ��ǰ������ùؼ���UESET����."
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//ELT_LOAD,NUM=...(I)  //NUM=��Ԫ�������� //�˹ؼ��ʿ�ȡ��BODY FORCE. �˹ؼ���Ŀǰ����ʹ�ã����һ������."
	README(IPP(I))=  "//"//'"'//"THE KEYWORD ELT_LOAD IS USED TO INPUT LOADS APPLIED ON AN ELEMENT GROUP."//'"'
	README(IPP(I)) = "//{ELTGROUP(I),DOF(I),VALUE(I),STEPFUNC.(I)}  //ELTGROUP=I, ���ؼ���ELEMENT�ж����SET��.�༴ΪESET(I)�ĵ�Ԫ."  
	README(IPP(I)) = "//{......} //��NUM��. "
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//DATAPOINT,NUM=...(I)  //NUM=������ݼ��ĸ���."
	README(IPP(I))=  "//"//'"'//"THE KEYWORD DATAPOINT IS USED TO OUTPUT VALUES IN SPECIFIC NODES ."//'"'
	README(IPP(I)) = "//{A:NNODE,ISSUMQ}  //NNODE=���ݼ��Ľڵ�����ISSUMQ/=0,��ʾ������˵㼯�ĸ��ڵ��������"  
	README(IPP(I)) = "//{A1:iNodes} //��NNODE�ڵ�. "
	README(IPP(I)) = "//{...} //��{A,A1}*NUM��. "

	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//TIME STEP,NUM=...(I)  //NUM=����."
	README(IPP(I))=  "//"//'"'//"THE KEYWORD TIME STEP IS USED TO INPUT THE SUBINCREMENTAL TIME FOR EACH STEP IN A TRANSIENT ANALYSIS."//'"'
	README(IPP(I)) = "//{ISTEP(I),NSUBTIMESTEP(I),TIMESTEP(1)(R),TIMESTEP(2)(R),...,TIMESTEP(STEP)(R)}  //��NSUBTIMESTEP��ʱ��,�ܵķ���ʱ�䲽��Ϊ���Ӳ���֮��."  
	README(IPP(I)) = "//{......} //��NUM��. "
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//INITIAL VALUE,NUM=...(I)   //NUM=��ֵ����"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD INITIAL VALUE IS USED TO INITIAL VALUE DATA."//'"'
	README(IPP(I)) = "//{NODE(I),DOF(I),VALUE(R)[,STEPFUNC.(I)]}"  
	README(IPP(I)) = "//{...}  //��NUM��"	
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//STEPINFO,NUM=...(I)   //NUM=����Ϣ����"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD STEPINFO IS USED TO STEP INFOMATION DATA."//'"'
	README(IPP(I)) = "//{ISTEP(I),MATHERSTEP(I),ISSTEADY(I),ISSTEADY(I)[,BCTYPE,LOADTYPE]} //ISTEP:STEP NUMBER; MATHERSTEP:�˲������в�; ISSTEADY:�Ƿ���̬����, 1,Yes,others NO."  
	README(IPP(I)) = "//LOADTYPE(BCTYPE):Ϊ����(λ�Ʊ߽�)ʩ�ӷ�ʽ��=1,��ʾ���ں���(λ�Ʊ߽�)��ʱ������ʩ�ӣ�=2(default)����ʾ������(λ�Ʊ߽�)�ڲ���˲��ʩ�ӡ�"
    README(IPP(I)) = "//{...}  //��NUM��"	
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//STEP FUNCTION,NUM=...(I),STEP=...(I)[,BASE=1,SCALE=1.0]   //NUM=�����̵ĸ���,STEP=����.,BASE=0����ͬʱ�����Ӧ����(��0��)��FACTOR����ʱҪ����STEP+1��factors"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD STEP FUNCTION IS USED TO STEP FUNCTION DATA."//'"'
	README(IPP(I)) = "//{[FACTOR(0)(R)],FACTOR(1)(R),FACTOR(2)(R),...,FACTOR(STEP)(R),TITLE(A)& 
                         \N// FACTOR(ISTEP)=��ISTEP���߽����ص�ϵ��. &
                         \N//��FACTOR(ISTEP)=-999ʱ,�˱߽������ڴ˲���ʧЧ�������ã�. &
                         \N//!��Ԫ����ʱ��0Ϊ��1Ϊ��.&
                         \N//scale,���ӵ���������,��factor=factor_input*scale.Ĭ��Ϊ1.����������"C  
	README(IPP(I)) = "//{...}  //��NUM��"


	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//SLAVE_MASTER_NODE_PAIR,NUM=...(I)   //NUM=Լ���ڵ�Եĸ���."  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD SMNP IS USED TO SLAVE-MASTER NODE PAIR DATA."//'"'
	README(IPP(I)) = "//{SLAVE,SDOF,MASTER,MDOF} //SLAVE�ڵ�ţ�SLAVE���ɶȣ�MASTER�ڵ�ţ�MASTER���ɶȣ� "  
	README(IPP(I)) = "//{...}  //��NUM��.!ע�⣬master���ϵ�mdof���ɶ���Ŀǰ����ʩ��λ�Ʊ߽���������ʱ����ʱ�뻥��slave��master."
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//WSP,NUM=...(I)[CP=...(I),Method=...(I)]   //CP=1,��Ϊ����ˮ����,ˮ���߼�����ɺ���˳���CP=2,���������Ǳ߽簴������ˮԾ����ȡֵ��Mothed=1,2,3,4(kinds,WangXG,Ohtsu,Chow)"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD IS USED TO INPUT WATERSURFACEPROFILE CALCULATION DATA."//'"'
	README(IPP(I)) = "//{NSEGMENT(I),NNODE(I),Q(R),B(R),UYBC,DYBC,[Q_SF,UYBC_SF,DYBC_SF,caltype,g,kn]}" 
	README(IPP(I)) = "//������;�����ʷ��ܽڵ���;������;������;���α߽�ˮ��;���α߽�ˮ��; [Q�Ĳ�������,UYBC�Ĳ�������,DYBC�Ĳ�������,������Ʋ���;�������ٶ�;��������ϵ��]"
	README(IPP(I)) = "//{UNODE(I),DNODE(I),n(R),So(R),[Profileshape(1),Profileshape(2)]}  //���ţ��ֲ���ţ������νڵ���Ϊ1�������α��ΪNNODE),�յ��,����,����,[������ʽ,������ʽ]��(���������У������ε�������������,��NSEGMENT��)"
    README(IPP(I)) = "//{N1,N2,...,N(NNODE)}  //�����ڽڵ㣨ȫ�֣����,��NNODE�������������������롣"
	README(IPP(I)) = "//ע��Profileshape=M1,M2,M3,S1,S2,S3,H2,H3,A2,A3,C1,C3 "
	README(IPP(I)) = "//ע��UYBC��DYBC��=-1,��ʾ���Σ����Σ��߽�ˮ��Ϊ�ٽ�ˮ�=-2����ʾ���Σ����Σ��߽�ˮ��Ϊ����ˮ�"
    README(IPP(I)) = "//ע��CALTYPE=1��ֻ�㼱����=2��ֻ�㻺����=0��Ĭ�ϣ������߶���."
    README(IPP(I)) = "//ע������/=4(Chow)��Ŀǰֻ�ܴ���NSEGMENT=2���������б�¶�+ˮƽ�Σ�������������������.Ĭ��ˮƽ����б�¶εĽ���Ϊ�ڶ��εĵ�һ���ڵ㡣"
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//SOILPROFILE,NUM=(I),[spmethod=(I),kmethod=(I),rf_epp=(I),rf_app=(I),iniepp=(I),soilspringmodel=0]   \n//spmethod=��ѹ�����㷽����0���ɿϣ� &
        \n// kmethod=����ϵ���ļ��㷽����0��m����1��(E,V)����2,zhu;3 biot;4 vesic;-1,��ֱ������. &
        \n// rf_epp=�����������ɿ�����ֵ�Ƿ�Ҫ������ʼ��������ѹ��.(0N1Y). &
        \n// RF_APP=0 !������������ѹ�����أ������������Ƿ񰴵������ۼ�.(0N1Y). & 
		\n// iniepp=2 !��������������ѹ����ʼֵ(λ��=0ʱ����ѹ��),2=������ѹ����1=��ֹ��ѹ��. &
		\n// soilspringmodel !������ģ�ͣ�=HYPERBOLIC,Ϊ˫��ģ��;=EIP_SPRING,����������(Ĭ��);=ELASTIC,Ϊ����"C
		
	README(IPP(I))=  "//"//'"'//"THE KEYWORD SOILPROFILE IS USED TO INPUT SOILPROFILE DATA."//'"'C
	README(IPP(I))=  "//A0:{TITLE(C)}  //�������������"  
	README(IPP(I)) = "//A:{NASOIL(I),NPSOIL(I),BEAMID(I)}  //������������(������������ѹ��Ϊ��������Ϊ������֮��Ȼ��)�����������������ػ�����." 
	README(IPP(I)) = "//B:{(Z1,Z2,MAT,WPMETHOD,STEPFUN [,PV])*NASOIL}   //�㶥�̵߳�ţ���׸̵߳�ţ����Ϻţ�ˮѹ�����Ƿ���(0=���㣬1=������㣬2=���㣬������͸��,3=�ֶ�����),������ [,����]����NASOIL��"  
	README(IPP(I)) = "//C:{(Z1,Z2,MAT,WPMETHOD,STEPFUN [,PV])*NPSOIL}   //�㶥�̵߳�ţ���׸̵߳�ţ����Ϻţ�ˮѹ�����Ƿ���(0=���㣬1=������㣬2=���㣬������͸��,3=�ֶ�����),������ [,����]����NPSOIL��"  
	README(IPP(I)) = "//D:{AWATERLEVEL,ASTEPFUN,PWATERLEVEL,PSTEPFUN}   //������ˮλ,������ˮλ��������������ˮλ,������ˮλ��������"  
	README(IPP(I)) = "//E:{ALoad,ALoadSTEPFUN,PLoad,PLoadSTEPFUN}   //�����೬��,�����೬�ز������������೬��,�����೬�ز�������" 
!	README(IPP(I)) = "//F:{NO_ACTION(I)}*NACTION   //Լ���ţ���NACTION��" 
!	README(IPP(I)) = "//G:{NO_STRUT(I)}*NSTRUT   //֧�źţ���NSTRUT��" 
	README(IPP(I)) = "//{A0,A,B,C,D,E}*NUM��   //��NUM��"
	README(IPP(I)) = "//ע�⣺\n//1)ÿһʱ�䲽�����㰴˳����϶������룬��ͬʱ�䲽�����������ص�,ͬһʱ�䲽��Ч�����㲻Ҫ�����ص��Ϳ�϶�� &
						\n//2)��ʱ�䲽����ˮλ��Ҫ�ֲ�;׮��(��)��Ҫ�ֲ�(���׮��������);׮�Ĳ��Ϸֽ紦Ҫ�ֲ㡣 &
						\n//3)������͸��ʱ���ٶ�awL>pwL. &
						\n//4)���ˮ����ڵر���ˮ��ЧΪΪ���㣨��c,phi,ģ������Ϊ0����͸ϵ��<=0)"C
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//PILE,NUM=...(I)"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD PILE IS USED TO INPUT BEAM(RetainingStructure) DATA."//'"'
	README(IPP(I)) = "//A:{NSEG(I),[SYSTEM=0]}  //���Ϸֶ����������" 
	README(IPP(I)) = "//B:{Z(1:NSEG+1)}   //���Ϸֶε�,Ӧ�������£���������ң����루���㵥ԪѰַ��"  
	README(IPP(I)) = "//C:{MAT(1:NSEG)}   //���β��Ϻ�"  
	README(IPP(I)) = "//{A,B,C}*NUM��   //��NUM��"	

	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//STRUT,NUM=...(I),[ISBAR=YES|NO]   //���ģ������ͬʱ���ڣ���֧���ø˵�Ԫ����ģ�⣬��ʱ��isbar=yes,��֮��isbar=no��Ĭ�ϣ�"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD STRUT IS USED TO INPUT STRUT(RetainingStructure) DATA."//'"'
	README(IPP(I)) = "//A:{Z,MAT,STEPFUN}  //��ţ����Ϻţ�������" 
	README(IPP(I)) = "//{A}*NUM��   //��NUM��"
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//KPOINT,NUM=...(I)"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD KPOINT IS USED TO INPUT KeyPoint DATA."//'"'
	README(IPP(I)) = "//A:{NO(I),XY(1:NDIMENSION),[ELEMENTSIZE(R)]}  //��ţ�XY(1:NDIMENSION),[ELEMENT SIZE]" 
 	README(IPP(I)) = "//{A}*NUM��   //��NUM��"	
    
 	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//GEOLINE,NUM=...(I)"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD GeoLINE IS USED TO INPUT GEOMETRY LINE DATA."//'"'
	README(IPP(I)) = "//A:{NO(I),MATID(I),NPOINT(I),IPOINT(I),[TITLE(A)]}  //�ߺţ����Ϻ�(ͶӰ)�����Ƶ���������Ƶ��(��NPOINT��)������" 
 	README(IPP(I)) = "//{A}*NUM��   //��NUM��"   
	README(IPP(I)) = "//ע�⣺1����X�������ҵ�˳������,��x1<=x2<=..<=xn" 
	
 	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//WATERLEVEL,NUM=...(I),FMT=0|1,VAR=1|2|3,SF=0  //VAR=INTERPOLATION INDEPENDENT VARIABLE. 1/2/3 MEANS X/Y/Z,SF=STEPFUNCTION"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD IS USED TO INPUT WALTERLEVEL LINE DATA."//'"'
    README(IPP(I)) = "//IF fmt=0(DEFAULT) then input the id number in the array KPOINT of the points. " 
    README(IPP(I)) = "//IF fmt=1 then INPUT (VAR,H) of the points. " 
	README(IPP(I)) = "//A:{P1,P2,...,P(NUM)]}  //(IF FMT=0)���Ƶ��(��NUM��)" 
    README(IPP(I)) = "//A:{XI,HI}*NUM  //(IF FMT=1)(��NUM��)" 
	README(IPP(I)) = "//ע�⣺1����X�������ҵ�˳������,��x1<=x2<=..<=xn" 	
	
 	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//RIGHT TURN POINT,NUM=...(I)"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD  IS USED TO INPUT RIGHT TURN POINT DATA."//'"'
	README(IPP(I)) = "//A:{IPOINT1,IPOINT2,...IPOINTN}  //���(��NUM��)" 
 	README(IPP(I)) = "//{A}*1��   "   
	README(IPP(I)) = "//����������������ֱ������¶˵�ĵ��(��ͬһX���������ر�߳�Y)���Խ��ͬһX���������ر�̵߳����⡣" 	
    
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//ACTION,NUM=...(I)"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD ACTION IS USED TO INPUT Laction DATA..."//'"'
	README(IPP(I))=  "//A0:{TITLE(C)}  //���õ�����"  
	README(IPP(I)) = "//A:{NKP,TYPE,DOF,NDIM,[SF,ISVALUESTEPFUN,ISEXVALUE]}  //���Ƶ������������ͣ�0=����1=λ�ƣ�2=�նȣ������õ����ɶȣ����õ�ά�ȣ�����������,ֵ����������(0N1Y)����ֵ����(0N1Y)" 
	README(IPP(I)) = "//B:{KPOINT(1:NKP)}  //���Ƶ��,Ӧ�������£���������ң����루���㵥ԪѰַ��" 
	README(IPP(I)) = "//C:{VALUE(1:NKP)}  //���Ƶ������õĴ�С��ע�ⵥλͳһ(�����õĿ��ȡ����֧��׮�ļ��)����ndim=1,type=2ʱ����λ��F/L^3;(ndim=1,type=1,L); (ndim=1,type=0,F/L^2) "
	README(IPP(I)) = "//D:{SF(1:NKP)}  //���Ƶ������õĲ����������isstepfun=1."
	README(IPP(I)) = "//E:{EXVALUE(1:NKP,1:2)}  //���Ƶ������õ�������ֵ(�ȸ�������(exvalue(:,1))�����������(exvalue(:,2))) ���isexvalue=1."
 	README(IPP(I)) = "//{A0,A,B,C,D,E}*NUM��   //��NUM��"
	
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//HINGE,NUM=...(I) //(�˹���Ŀǰ����beam2d��beam��Ԫ��Ч)"  
	README(IPP(I))=  "//"//'"'//"THE KEYWORD HINGE IS USED TO INPUT HINGE/FREEDOF DATA..."//'"'
	README(IPP(I)) = "//A:{ELEMENT,NODE,DOF}  //(Ҫ�ͷŵ����ɶ����ڵ�)��Ԫ���ڵ�����ɶȱ�š�" 
 	README(IPP(I)) = "//{A}*NUM��   //��NUM��"
    
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//SLOPEPARAMETER,SLOPEMETHOD=0|1|2|3|4|5,OPTIMIZATIONMETHOD=1,SLIPSHAPE=CIRCULAR|NONCIRCULAR,SLICEWIDTH=1.0,..." 
	README(IPP(I))=  "//"//'"'//"THE KEYWORD HINGE IS USED TO INPUT HINGE/FREEDOF DATA..."//'"'
	README(IPP(I)) = "//SLOPEMETHOD:���·���������1��ordinary,2,bishop,3,spencer,4,janbu,5,gle; 0 for all." 
 	README(IPP(I)) = "//OMTIMIZATION,��Ѱ������1��grid;2,MONTE CARLO"	
    README(IPP(I)) = "//SLIPSHAPE,������״��1,circular;0,noncircular;"
    README(IPP(I)) = "//SLICEWIDTH,�������,(1.0,bydefault)"
    README(IPP(I)) = "//XMIN_MC,XMAX_MC,SEARCH RANGE FOR MONTE CARLO TECHNIQUE"
    
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) = "//MATERIAL,MATID=...(I),[TYPE=...(I)],[ISFF=YES|NO],[NAME=...(c)],ISSF=...(I)//MATID=���Ϻţ�TYPE=�������ͣ�ISFF=�Ƿ�Ϊ����ĳ������,NAME=��������ע��.�˹ؼ��ʿ��ظ�����.ISSF=�Ƿ�������ϲ����Ĳ�����(0N1Y)" 
	README(IPP(I))=  "//"//'"'//"THE KEYWORD MATERIAL IS USED TO INPUT MATERIAL INFOMATION."//'"'
	README(IPP(I)) = "//{A: MATERIAL.PROPERTIES}  //PLEASE REFER TO THE EXPLANATION IN THE END OF THE FILE."
	README(IPP(I)) = "//[B:{FIELD FUNCTION PARAMETERS}]  // ��ISFF=YES "	
	README(IPP(I)) = "//[C:{STEPFUNCTION(1:NPARAMETERS)}]  // ��ISSF=YES "	
	
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
    README(IPP(I)) ="//ע�ⵥλ��ͳһ������ע��ALPHA�ĵ�λ����VGģ����Ϊ��1/L����LRģ����Ϊ:L.\N"C
	README(IPP(I)) ="//8. CASE(BAR,BAR2D)     : PROPERTY(1)=E  .(2)=A	[.(3)=hy	.(4)=hz	 .(5)=minN(�������ѹ��)	.(6)=MaxN(�����������)]"
	README(IPP(I)) ="//9. CASE(BEAM,BEAM2D,SSP2D)   : PROPERTY(1)=E  .(2)=A   .(3)=v	.(4)=J	.(5)=Iz	.(6)=Iy	[.(7)=hy	.(8)=hz]"
    README(IPP(I)) ="//                          [PROPERTY(9)=MinN .(10)=MaxN .(11)=minMx .(12)=maxMx .(13)=minMy .(14)=maxMy .(15)=minMz .(16)=maxMz]"
    README(IPP(I)) ="//                          [.(17)=yc ]"
	README(IPP(I)) ="//J,IY,IZ ARE REFERRED TO THE LOCAL COORDINATE SYSTEM."
	README(IPP(I)) ="//��Ϊ��ԪSSP2Dָ�����ϣ��������(A,Iz,minN,maxN,minMz,maxMz)��Ϊ�����ְ�׮�Ĳ���,YC�ְ�׮����������������صľ��롣"
	README(IPP(I)) ="//hy,hz�ֱ�Ϊ��������y'��z'�᷽��(�ֲ�����)�ĸ߶ȣ�Ϊ����ת��Ϊ�����嵥Ԫʱ����. \N"C
	README(IPP(I)) ="//��Ϊ֧��׮ʱ��hy,hz�ֱ�Ϊ׮����׮�࣬����������֧��׮�ϵ���ѹ���������ɵĸն��� \N"C
	README(IPP(I)) ="//10. CASE(wellbore) : PROPERTY(1)=R(���뾶). (2)Temp.\n &
                    &  .(3)=Kr(relative roughness of the  inner surface of the pipe ), \N&
                    &  .(4)=time unit (day=0(default),second=1), \N&
                    &  .(5)=g (gravity acc. =0(SET by default, 73156608000.00 m/day2. ), \N&
                    &  .(6)PIPE-FLOW MODEL(=0,Darcy(no porous effect,no friction effect,DEFAULT);=1,Siwon; =2,OUYang EFFECT;=3,Input by user;=4,(no porous effect,DEFAULT)); \N &
                    &  .(7)��Ƥ�ĺ�Ⱥ���͸ϵ���ı�L/K.=0(Ĭ��,�����Ǿ���)); &
                    &  .(8)=f(Darcy Friction factor(����Ħ��ϵ��),=0,�ɼ��㶨(Ĭ��)��һ��Ϊ0.02-0.03����,��model=3ʱ���롣);.(9)=POROSITY OF THE WELLBORE(��model=1ʱ����)"C
	README(IPP(I)) ="//11. CASE(ExcavationSoil/SLOPESOIL) : PROPERTY(1:10)=������Ħ���ǣ���Ȼ/�����ضȣ�����ģ�������ɱ�,��͸ϵ����ˮƽ����ϵ��(F/L**3),ǽ����Ħ���ǣ��ȣ�,[�㶥ˮѹ��PW1,���ˮѹ��PW2]"C
	README(IPP(I)) ="//PW1��PW2���������ˮѹ��Ϊ�ֶ�����ʱ���ã�WPMETHOD=3����"
	README(IPP(I)) ="//12. CASE(spring) : PROPERTY(1:3)=k,minV(������λ��),maxV(������λ��)��Ԥ������Ԥ��λ��"C
	
						
	README(IPP(I)) ="\N//******************************************************************************************************"C
	README(IPP(I)) ="//Coordinate system(FYI)" 
	README(IPP(I)) ="//����2D����,��������ΪY�᷽��������Ϊ��,X������Ϊ��."
	README(IPP(I)) ="//����3D����,��������ΪZ�᷽��������Ϊ��,ˮƽ��ΪXYƽ��,XYZ����������������"
	README(IPP(I)) ="//�����൥Ԫ�ľֲ�����ϵ����������ϵ��ͬ��"C
	
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
		write(20,20) readme(j)
	end do
	
	!tof=system("D:\README_FEASOLVER.TXT")	
	close(20)
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
