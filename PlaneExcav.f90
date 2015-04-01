Module ExcaDS
    
    INTEGER,PARAMETER::Double=KIND(1.0D0)
    INTEGER::showvalue,CGSTEP=0
    REAL(DOUBLE)::DEFAULTSIZE=0.2D0
    INTEGER::BGC=1,linecolor=0,IsMarker=0,IsThinLine=0 !绘图背景色：默认黑色，=0，白色
!    type soil_tydef
 !       real(double)::c=0.d0,fi=0.d0,g=0.d0,Eo=0.d0,v=0.d0,k=-1.d0,ks=0.d0 !黏聚力，摩擦角，重度，变形模量，泊松比,渗透系数，水平基床系数
              
 !   end type
    !integer::nsoil=0
    !type(soil_tydef),allocatable::soil(:)
    real(double),allocatable::kpoint(:,:)
    integer::nkp=0
	integer,allocatable::kpnode(:),KPELEMENT(:,:)
    
	
    type soillayer_tydef
        integer::z(2)=0,mat=1,sf=0,wpflag=1,nnode=0,nel=0  !土层号，步函数号，是否水土分算（1Yes0NO）,此段范围内土层所依附的支护梁的节点个数,单元数
        integer::IUNSET=-1,IUESET=-1 !此段范围内土层所依附的支护梁的节点集和单元集编号
		real(double)::ko=0.d0,ka=0.d0,kp=0.d0 !静止、主动、被动土压力系数 
        real(double)::SigmaVT(2),SigmaV(2),Sigmako(2),Sigmaka(2),Sigmakp(2),PW(2)=0.0d0,KS(2)=0.D0		
        !（竖向应力，有效竖向应力，静止土压力，主动土压力，被动土压力，孔压）。1表层顶的值，2为层底的值。
		
    endtype
    
    type soilprofile_tydef
        integer::nasoil=0,npsoil=0,beam=1,spm=0,naction=0,nstrut=0,kmethod=0 !spm=土压力计算方法,0=郎肯。 
		!kmethod,水平基床的计算方法: =0,m法(规范方法)；=1(E，v方法，);
        integer::sf_awL=0,sf_pwL=0,sf_aLoad=0,sf_pLoad=0 !主被侧水位的步长函数，主被动侧超载的步函数
		real(DOUBLE)::awL,pwL,aLoad=0,pLoad=0 !主被侧水位,主被动侧超载
		integer::za=0,zp=0 !各时间步主动侧和被动侧地表高程。
		integer::aside=1 !aside=1,表主动土压力正，被动为负。
        CHARACTER(64)::TITLE=""
		type(soillayer_tydef),allocatable::asoil(:),psoil(:)
		integer,allocatable::iaction(:),istrut(:)
        real(double),allocatable::stepload(:,:)
		
    end type
	type(soilprofile_tydef),allocatable::soilprofile(:)
	integer::nsoilprofile
    
    type beam_tydef
		integer::nseg=1,nnode=0,nel=0,system=0
        integer,allocatable::mat(:),kpoint(:)        
		integer,allocatable::node(:),element(:)
		integer,allocatable::element_sp(:,:) !节点的土弹簧单元号
		integer,allocatable::node2BCload(:)
        real(double),allocatable::Beamresult(:,:,:) !(node,13,istep),/dis,M,Q,Load/
		INTEGER::NVA=14
		real(double),allocatable::Nlength(:) !节点所辖的长度
		!element_sp(2,nnode),element_sp(1,i)=节点i主动侧弹簧单元在element的位置
		!element_sp(2,i)=节点被动侧弹簧单元在element的位置
    endtype
    
	type(beam_tydef),allocatable::pile(:)
	integer::npile=0
	
    type strut_tydef
        integer::z(2)=0,sf=0,mat=0,ELEMENT=0
		integer::isbar=0,ISINITIAL=0
        real(double)::PreLoad=0,PreDis=0
    end type
    type(strut_tydef),allocatable::strut(:)
	integer::nstrut=0
	
	type Action_tydef
		integer::nkp=0,type=0,dof=0,ndim=0,sf=0,IUNSET=-1,IUESET=-1,nstiffelement=0,nbcnode=0,nloadnode=0
		!iunset,iueset为此作用所依附的支护梁单元的节点组合单元组
		CHARACTER(64)::TITLE=""
		!type=0,force;type=1,displacement;=2,spring.
		!dof,施加在哪个自由度上。
        !Action的维度，0=point,1=line.
		integer,allocatable::kpoint(:),vsf(:),NODE(:) !控制点号,各控制点的步函数,控制点在node()中的编号
		real(double),allocatable::value(:) !控制点上的值。
		real(double),allocatable::exvalue(:,:) !如果为弹簧，为弹簧力的上下限值。
		integer::istiffelement=0 !IF(TYPE==2) 生成弹簧单元进行处理，在第一个节点unset(iunset).node(1)生成第一个单元，istiffelement存储第一个单元ELEMENT中的位置。
		integer::iBCNode=0   !IF(TYPE==1) 生成位移边界进行处理，在第一个节点unset(iunset).node(1)生成第一个边界，ibcnode存储第一个边界bc_disp中的位置。
		integer::iLoadNode=0 !IF(TYPE==0) 生成力进行处理，在第一个节点unset(iunset).node(1)生成第一个力，ibcnode存储第一个力在bc_load中的位置。
		integer,allocatable::node2stiffelement(:)
		integer,allocatable::node2bcdisp(:)
		integer,allocatable::node2bcload(:)
	endtype 
	type(action_tydef),allocatable::action(:)
	integer::naction=0
	
    contains
    


    
    end module
    
module ExcaDSLIB
    use ExcaDS
    interface
	    subroutine enlarge_strut(PROP,NP,EXN,SN)
        !扩大PROP数组,同时update总的单元数NP=NP+EXN
        !EXN:扩大的单元个数
        !SN,:扩容部分的起位
            USE EXCADS
	        integer,intent(in)::EXN
            INTEGER,INTENT(IN OUT)::NP
            type(strut_tydef),INTENT(IN OUT),ALLOCATABLE::PROP(:)
	        integer,intent(out)::SN        
        ENDSUBROUTINE
    end interface
end module 


subroutine Excavation(istep)
    use ExcaDS
    use solverds
    implicit none
    integer,intent(in)::istep
    integer::i,j,k,flag1,nc1,nc2,nc3,NC4,NC5,NC6
	type(soillayer_tydef)::soil1
	
	call EarthPressure(istep)
    !print *, "pass1"
    call soilspring_EXCA(istep)
    !print *, "pass2"
	CALL SOILLOAD_EXCA(ISTEP)
    !print *, "pass3"
	call ACTION_EXCA(istep)
	
    
	if(istep==1) then
		open(12,file=EXCAMSGFILE,status='replace')
		write(12,30)
		write(12,10)
	else
		open(12,file=EXCAMSGFILE,status='old',access='append')
	endif
	nc1=0
	do i=1,nsoilprofile
		do j=1,soilprofile(i).nasoil+soilprofile(i).npsoil
			if(j<=soilprofile(i).nasoil) then
				nc2=j
				nc3=1
				SOIL1=soilprofile(i).asoil(j)
				flag1=soilprofile(i).aside
			else
				nc2=j-soilprofile(i).nasoil
				nc3=-1
				SOIL1=soilprofile(i).psoil(nc2)
				flag1=-soilprofile(i).aside
			endif			
			if(sf(SOIL1.sf).factor(istep)==0) cycle
			do k=1,2
				nc1=nc1+1
				write(12,20) nc1,nc2,istep,i,nc3,kpoint(ndimension,soil1.z(k)), &
				flag1*soil1.sigmaVT(k),flag1*soil1.sigmaV(k),soil1.sigmako(k)+flag1*soil1.Pw(k),soil1.sigmaka(k)+flag1*soil1.Pw(k),& 
                    soil1.sigmakp(k)+flag1*soil1.Pw(k), & 
				flag1*soil1.Pw(k),soil1.ks(k),MATPROPERTY(SOIL1.MAT,1,ISTEP),MATPROPERTY(SOIL1.MAT,2,ISTEP), &
                    MATPROPERTY(SOIL1.MAT,3,ISTEP), MATPROPERTY(SOIL1.MAT,6,ISTEP)
				
			end do
		enddo
		
		nc2=soilprofile(i).beam
		DO K=1,2
			do J=1,pile(soilprofile(i).beam).nnode
				NC4=PILE(NC2).node(J)
				nc3=pile(nc2).element_sp(K,NC4)	
				IF(K==1) THEN
					NC5=1
				ELSE
					NC5=-1
                ENDIF
                IF(.NOT.ALLOCATED(ELEMENT(NC3).GFORCE)) ALLOCATE(ELEMENT(NC3).GFORCE(1))
				WRITE(12,40) NC3,ELEMENT(NC3).NODE(1),ISTEP,I,NC5, &
					node(element(nc3).node(1)).coord(1:2),ELEMENT(NC3).PROPERTY(1:3)+ELEMENT(NC3).PROPERTY(6), &
					ELEMENT(NC3).PROPERTY(6),-ELEMENT(NC3).GFORCE(1),ELEMENT(NC3).PROPERTY(4:5)				
			enddo
		ENDDO
		
	enddo
	

	close(12)
	
10  format(7X,"SOILINFO",13X,"NO",10X,"Layer",10X,"ISTEP",3X,"ISOILPROFILE",11X,"SIDE", &
             11X, "Z(L)",5X,"SvT(F/L^2)",6X,"Sv(F/L^2)",2X,"Sko+Pw(F/L^2)", &
               2X,"Ska+Pw(F/L^2)",2X,"Skp+Pw(F/L^2)",6X,"PW(F/L^2)",6X,"Ks(F/L^3)", &
            7X,"C(F/L^2)",9X,"PHI(O)",3X,"GAMMA(F/L^3)",9X,"K(L/T)")	
20	format(7X,"SOILINFO",5(X,I14),12(X,E14.7))
30	format(5X,"SOILSPRING",10X,"ELENO",7X,"ELE_NODE",10X,"ISTEP",3X,"ISOILPROFILE",11X,"SIDE", &
		   11X,"X(L)",11X,"Y(L)",5X,"FKO+FPw(F)",5X,"FKA+FPw(F)",5X,"FKP+FPw(F)",9X,"FPW(F)",4X,"F-SPRING(F)",8X,"KS(F/L)",10X,"eL(L)")
40	format(5X,"SOILSPRING",5(X,I14),9(X,E14.7))
50	format(7X,"BEAMINFO","X(L)","Y(L)","DISX(L)","Q(F)","LOAD(F)","M(F.L)")

end subroutine

SUBROUTINE SOILLOAD_EXCA(ISTEP)
	USE EXCADS
	USE SOLVERLIB
	IMPLICIT NONE
    INTEGER,INTENT(IN)::ISTEP
	INTEGER::I,J,K,IPILE,NS1,NS2,IELA1(2),INNODE1(2),NC1,NC2
	
	DO I=1,NSOILPROFILE
		IPILE=SOILPROFILE(I).BEAM
        
		NS1=LBOUND(PILE(IPILE).ELEMENT_SP,DIM=2)
		NS2=UBOUND(PILE(IPILE).ELEMENT_SP,DIM=2)
        
		BC_LOAD(PILE(IPILE).NODE2BCLOAD).VALUE=0
        NC2=PILE(IPILE).NODE2BCLOAD(NS2)
        NC1=PILE(IPILE).NODE2BCLOAD(NS1)
		IF(.NOT.ALLOCATED(SOILPROFILE(I).STEPLOAD)) THEN
            ALLOCATE(SOILPROFILE(I).STEPLOAD(NS1:NS2,NSTEP))            
        ENDIF
		if(istep==1) SOILPROFILE(I).STEPLOAD=0.D0
		DO J=NS1,NS2
			IELA1=PILE(IPILE).ELEMENT_SP(:,J)
            
			FORALL (K=1:2) INNODE1(K)=ELEMENT(IELA1(K)).NODE(1)
            
			!BC_LOAD(INNODE1).VALUE=BC_LOAD(INNODE1).VALUE+ELEMENT(IELA1).PROPERTY(6)
			!BC_LOAD(INNODE1(1)).VALUE=BC_LOAD(INNODE1(1)).VALUE+ELEMENT(IELA1(1)).PROPERTY(2)
			!BC_LOAD(INNODE1(2)).VALUE=BC_LOAD(INNODE1(2)).VALUE+ELEMENT(IELA1(2)).PROPERTY(1)
            
            SOILPROFILE(I).STEPLOAD(INNODE1,ISTEP)=SOILPROFILE(I).STEPLOAD(INNODE1,ISTEP)+ELEMENT(IELA1).PROPERTY(6)
            SOILPROFILE(I).STEPLOAD(INNODE1(1),ISTEP)=SOILPROFILE(I).STEPLOAD(INNODE1(1),ISTEP)+ELEMENT(IELA1(1)).PROPERTY(2)
			SOILPROFILE(I).STEPLOAD(INNODE1(2),ISTEP)=SOILPROFILE(I).STEPLOAD(INNODE1(2),ISTEP)+ELEMENT(IELA1(2)).PROPERTY(2)
        ENDDO
        
        !SOILPROFILE(I).STEPLOAD(BC_LOAD(NC1:NC2).NODE,ISTEP)=BC_LOAD(NC1:NC2).VALUE
        !求增量 
        IF(ISTEP>1) THEN        
            BC_LOAD(NC1:NC2).VALUE=SOILPROFILE(I).STEPLOAD(BC_LOAD(NC1:NC2).NODE,ISTEP) &
                -SOILPROFILE(I).STEPLOAD(BC_LOAD(NC1:NC2).NODE,ISTEP-1)
        ELSE
            BC_LOAD(NC1:NC2).VALUE=SOILPROFILE(I).STEPLOAD(BC_LOAD(NC1:NC2).NODE,ISTEP)
        ENDIF
        
	ENDDO

ENDSUBROUTINE

SUBROUTINE ACTION_EXCA(ISTEP)
    USE ExcaDS
    USE SOLVERLIB
    IMPLICIT NONE
    integer,intent(in)::istep	
	integer::I,j,k,k1,IAC1,nnode1(2),nc1,nc2,iel1,iel2(2),n1,ndim1
	real(double)::t1,d1,b1,bp1,v1(2)
	real(double),allocatable::ar1(:)
    real(double),external::interpolation
	
	!do i=1,nsoilprofile
	!	do j=1,soilprofile(i).naction
	!		IAC1=soilprofile(i).iaction(j)
	!	enddo    
	!end do
    
    if(isexca2D==1) ndim1=ndimension
    if(isexca2D==2) ndim1=1

	
	do IAC1=1,NACTION
	
		if(allocated(ar1)) deallocate(ar1)
		allocate(ar1,MOLD=action(iac1).value)

		forall (k=1:action(iac1).nkp) ar1(k)=action(iac1).value(k)*sf(action(iac1).vsf(k)).factor(istep)
	
		if(action(IAC1).NDIM==1) then
			nc1=2
			nc2=ueset(action(iac1).iueset).enum
		else
			nc1=1
			nc2=action(iac1).nkp				
		endif
		
		!intialization
		n1=size(action(iac1).node2stiffelement)
		if(allocated(action(iac1).node2stiffelement))  &
			forall (k1=1:n1,k=1:6,action(iac1).node2stiffelement(k1)>0) element(action(iac1).node2stiffelement(k1)).property(k)=0.d0
		where(action(iac1).node2bcdisp>0) bc_disp(action(iac1).node2bcdisp).value=0
		where(action(iac1).node2bcload>0) bc_load(action(iac1).node2bcload).value=0
		
		do k=1,nc2
			if(action(IAC1).NDIM==1) then
				iel1=ueset(action(iac1).iueset).element(k)
				nnode1=element(iel1).node	
				t1=(node(nnode1(1)).coord(1)-node(nnode1(2)).coord(1))**2+ &
					(node(nnode1(1)).coord(2)-node(nnode1(2)).coord(2))**2 + &
					(node(nnode1(1)).coord(3)-node(nnode1(2)).coord(3))**2
				t1=t1**0.5d0/2.0d0
				
				d1=matproperty(element(ieL1).mat,7,istep) !假定同一土层内的支护桩的参数保持一致。
				b1=matproperty(element(ieL1).mat,8,istep) !假定同一土层内的支护桩的参数保持一致。
			
				
			else
				nnode1(1)=action(iac1).node(k)
				t1=1.d0
                b1=1.d0
			endif
			
			
			do k1=1,nc1

				v1(k1)=b1*interpolation(kpoint(ndim1,action(iac1).kpoint(1:action(iac1).nkp)),&
							 ar1,action(iac1).nkp,node(nnode1(k1)).coord(ndim1))
			enddo
			if(nc1==1) v1(2)=v1(1)

			
			select case(action(iac1).type)
				case(2)
					iel2(1:nc1)=action(iac1).node2stiffelement(nnode1(1:nc1))    
					element(iel2(1:nc1)).property(4)=element(iel2(1:nc1)).property(4)+(v1(1)+v1(2))*t1/2.
					do k1=1,nc1
						if(.not.allocated(element(iel2(k1)).km)) allocate(element(iel2(k1)).km(1,1))
						element(iel2(k1)).km(1,1)=element(iel2(k1)).property(4)
					enddo
					do k1=1,nc1
						v1(k1)=interpolation(kpoint(ndim1,action(iac1).kpoint(1:action(iac1).nkp)),&
							 action(iac1).exvalue(1:action(iac1).nkp,1),action(iac1).nkp,node(nnode1(k1)).coord(ndim1))
					enddo
					if(nc1==1) v1(2)=v1(1)
					element(iel2(1:nc1)).property(2)=element(iel2(1:nc1)).property(2)+(v1(1)+v1(2))*t1/2.
					do k1=1,nc1
						v1(k1)=interpolation(kpoint(ndim1,action(iac1).kpoint(1:action(iac1).nkp)),&
							 action(iac1).exvalue(1:action(iac1).nkp,2),action(iac1).nkp,node(nnode1(k1)).coord(ndim1))
					enddo
					if(nc1==1) v1(2)=v1(1)
					element(iel2(1:nc1)).property(3)=element(iel2(1:nc1)).property(3)+(v1(1)+v1(2))*t1/2.
					element(iel2(1:nc1)).property(5)=element(iel2(1:nc1)).property(5)+t1*2
				case(1)
					bc_disp(action(iac1).node2bcdisp(nnode1(1:nc1))).value=v1(1:nc1) !位移边界没有叠加性						
				case(0)
					bc_load(action(iac1).node2bcload(nnode1(1:nc1))).value=bc_load(action(iac1).node2bcload(nnode1(1:nc1))).value+(v1(1)+v1(2))*t1/2.
			endselect				
		enddo		
	
	end do
ENDSUBROUTINE



subroutine GenElement_EXCA2() !STRUCTURAL MESH
	!use solverds
	use excaDS
    use SOLVERLIB
	implicit none
	integer::i,j,k,K1,IAR(500),NAR=500,NNODE1=0,NEL1=0,N1=0,N2=0,IDP,JDP,IPILE,ISP1,IAC1
	real(double)::DPOINT1(3,1000)=0,t1
	INTEGER::et1,nnum1,ndof1,ngp1,nd1,ec1
	CHARACTER(16)::STYPE,CH1,CH2
    CHARACTER(32)::TITLE1=""
    !integer,allocatable::kpelement(:)
	type(node_tydef),ALLOCATABLE::node1(:),NODE2(:)
    type(element_tydef),ALLOCATABLE::element1(:),ELEMENT2(:)
	type(bc_tydef),allocatable::bf1(:),bf2(:)
	
	
    allocate(kpelement(2,nkp)) 
	!单元定位表,记录与每个控制点相连在单元在element()中的编号，(1,i)为i点的上（左）边的单元，(2,i)则为i点的下（右）边的单元
	!假定梁从左往右或从上往下输入。
	kpnode=0
    kpelement=0
    
	do ipile=1,npile
	
		ALLOCATE(element1(10000))
		ALLOCATE(node1(NNUM+1:NNUM+10000))	
		NNODE1=NNUM
		NEL1=0
		
		ET1=BEAM2D		
		call ettonnum(et1,nnum1,ndof1,ngp1,nd1,stype,ec1)
		
		do j=1,PILE(IPILE).nseg	

			NAR=SIZE(IAR)
			call aovis2d(PILE(IPILE).kpoint(j),PILE(IPILE).kpoint(j+1),IAR,NAR)
			do k=1,NAR-1
                IDP=3
                JDP=1000
				call divideLine2D(kpoint(1,IAR(k)),kpoint(ndimension,IAR(k)),kpoint(ndimension+1,IAR(k)), &
                                  kpoint(1,IAR(k+1)),kpoint(ndimension,IAR(k+1)),kpoint(ndimension+1,IAR(k+1)),&
                                  DPOINT1,IDP,JDP)
                
				kpnode(iar(k))=nnode1+1
                kpnode(iar(k+1))=nnode1+JDP !!!!
                if(kpelement(2,iar(k))==0) then
					kpelement(2,iar(k))=enum+nel1+1
				else
					stop "ERROR #1 OCCURED IN SUB GenElement_EXCA"
				endif
				 if(kpelement(1,iar(k+1))==0) then
					kpelement(1,iar(k+1))=enum+nel1+JDP-1
				else
					stop "ERROR #2 OCCURED IN SUB GenElement_EXCA"
				endif
                
                
                N1=JDP-1
				IF(J==PILE(IPILE).NSEG.and.(K==NAR-1)) N1=JDP
                DO K1=1,NDIMENSION
                    NODE1((NNODE1+1):(NNODE1+N1)).COORD(K1)=DPOINT1(K1,1:N1)
				ENDDO
                DO K1=1,JDP-1
                    NEL1=NEL1+1
                    ELEMENT1(NEL1).NNUM=NNUM1
					ALLOCATE(ELEMENT1(NEL1).NODE(NNUM1))
                    ELEMENT1(NEL1).NODE(1)=NNODE1+K1
                    ELEMENT1(NEL1).NODE(2)=ELEMENT1(NEL1).NODE(1)+1
                    ELEMENT1(NEL1).MAT=PILE(IPILE).MAT(J)
					ELEMENT1(NEL1).ET=BEAM2D
					element1(NEL1).ndof=ndof1
					element1(NEL1).ngp=ngp1
					element1(NEL1).nd=nd1
					element1(NEL1).ec=ec1
					ELEMENT1(NEL1).SET=NESET+1 !!!!
					!ELEMENT(NEL1).SYSTEM=0  !!!!!
                ENDDO				
                NNODE1=NNODE1+N1
                
            enddo
			
		enddo
		
		neset=neset+1
		eset(neset).num=NESET
		eset(neset).stype=stype
		eset(neset).grouptitle="BEAM"
		eset(neset).et=et1
		eset(neset).ec=ec1
		eset(neset).system=0  !!!!!
		eset(neset).enums=enum+1
		
		allocate(element2(enum+NEL1),PILE(IPILE).ELEMENT(NEL1))
		element2(1:enum)=element(1:enum)
		element2(enum+1:ENUM+NEL1)=element1(1:NEL1)
		if(allocated(element))	deallocate(element)
		deallocate(element1)
		enum=enum+NEL1
		eset(neset).enume=enum
		allocate(element(enum))
		element=element2
		deallocate(element2)
        DO J=1,NEL1
            PILE(IPILE).ELEMENT(J)=eset(neset).ENUMS-1+J
        ENDDO
		pile(ipile).nel=nel1
        
		
		allocate(NODE2(NNODE1),PILE(IPILE).NODE(NNODE1-NNUM))
		NODE2(1:NNUM)=NODE(1:NNUM)
		NODE2(NNUM+1:NNODE1)=NODE1
        
        DO J=NNUM+1,NNODE1
            PILE(IPILE).NODE(J-NNUM)=J
        ENDDO
		pile(ipile).nnode=nnode1-nnum
		
		ALLOCATE(pile(ipile).Nlength(NNUM+1:NNODE1), &
				 pile(ipile).BeamResult(NNUM+1:NNODE1,pile(ipile).NVA,nstep))
		pile(ipile).BeamResult=0.d0
		
		if(nsoilprofile>0) then
			ALLOCATE(PILE(IPILE).NODE2BCLOAD(NNUM+1:NNODE1))
			CALL enlarge_bc(BC_LOAD,BL_NUM,pile(ipile).nnode,N1)
			FORALL (J=1:PILE(IPILE).NNODE)
				BC_LOAD(N1+J-1).NODE=NNUM+J
				PILE(IPILE).NODE2BCLOAD(NNUM+J)=N1+J-1
				BC_LOAD(N1+J-1).DOF=1
			ENDFORALL
		endif
		
		
		if(allocated(NODE))	deallocate(NODE)
		deallocate(NODE1)
		NNUM=NNODE1
		allocate(NODE,SOURCE=NODE2)
		deallocate(NODE2)
		
		pile(ipile).nlength=0.0d0
		do j=1,pile(ipile).nnode-1
			t1=0.d0
			do k1=1,ndimension
				t1=t1+(node(pile(ipile).node(j)).coord(k1)-node(pile(ipile).node(j+1)).coord(k1))**2
			enddo
			t1=t1**0.5
			pile(ipile).Nlength(pile(ipile).node(j:j+1))=pile(ipile).Nlength(pile(ipile).node(j:j+1))+t1/2.0d0
		enddo
		
	enddo
	
	
	do i=1,NSOILPROFILE
	
	    IPILE=SOILPROFILE(I).BEAM

		!防止竖向失稳，在底部生成1坚向弹簧单元
		CALL enlarge_element(ELEMENT,ENUM,1,N1)
		element(n1).et=springy
		element(n1).nnum=1
		allocate(element(n1).node(1))
		element(n1).node(1)=pile(ipile).node(pile(ipile).nnode)
		et1=springy
		call ettonnum(et1,nnum1,ndof1,ngp1,nd1,stype,ec1)
		element(n1).ndof=ndof1
		element(n1).ngp=ngp1
		element(n1).nd=nd1
		element(n1).ec=ec1
		element(n1).property(2)=-1.d20
		element(n1).property(3)=1.d20
		element(n1).property(4)=1.0d7
		
		
		n1=pile(ipile).node(pile(ipile).nnode)
		allocate(pile(ipile).element_sp(2,pile(ipile).node(1):n1))
		pile(ipile).element_sp=0
		nel1=2*pile(ipile).nnode
		CALL enlarge_element(ELEMENT,ENUM,NEL1,N1)
		
		et1=SOILSPRINGX
		call ettonnum(et1,nnum1,ndof1,ngp1,nd1,stype,ec1)
		element(N1:).nnum=1
		element(N1:).et=SOILSPRINGX
		element(N1:).ndof=ndof1
		element(N1:).ngp=ngp1
		element(N1:).nd=nd1
		element(N1:).ec=ec1
		element(N1:).isactive=0
		do j=1,pile(ipile).nnode						
			allocate(element(N1-1+2*j-1).node(1),element(N1-1+2*j).node(1))
			element(N1-1+2*j-1).node(1)=pile(ipile).node(j)			
			element(N1-1+2*j).node(1)=pile(ipile).node(j)
			pile(ipile).element_sp(1,pile(ipile).node(j))=N1-1+2*j-1 !!!!
			pile(ipile).element_sp(2,pile(ipile).node(j))=N1-1+2*j   !!!!
		enddo
		
		neset=neset+1
		eset(neset).num=NESET
		eset(neset).stype=stype
		eset(neset).grouptitle="soilspring"
		eset(neset).et=et1
		eset(neset).ec=ec1
		eset(neset).system=0  !!!!!
		eset(neset).enums=N1
	    eset(neset).enume=enum

		
		
        
        !ueset,unset
        do j=1, soilprofile(i).nasoil
            
			if(kpnode(soilprofile(i).asoil(J).z(1))>0) then
                N1=SOILPROFILE(I).ASOIL(J).Z(1)
                N2=SOILPROFILE(I).ASOIL(J).Z(2)
                WRITE(CH1,'(I3)') I
                WRITE(CH2,'(I3)') J
                TITLE1="SP_"//TRIM(CH1)//"ASOIL_"//TRIM(CH2)
				CALL Gen_USER_NSET_EXCA(N1,N2,TITLE1,SOILPROFILE(I).ASOIL(J).IUNSET)
				CALL Gen_USER_ESET_EXCA(N1,N2,TITLE1,SOILPROFILE(I).ASOIL(J).IUESET)							 
			else
				print *, "The soilprofile(i).asoil(j) has no corresponding elements.i,j=",i,j
					
			endif

        end do
        
	    do j=1, soilprofile(i).npsoil
            
			if(kpnode(soilprofile(i).psoil(J).z(1))>0) then
                N1=SOILPROFILE(I).psoil(J).Z(1)
                N2=SOILPROFILE(I).psoil(J).Z(2)
                WRITE(CH1,'(I3)') I
                WRITE(CH2,'(I3)') J
                TITLE1="SP_"//TRIM(CH1)//"psoil_"//TRIM(CH2)
				CALL Gen_USER_NSET_EXCA(N1,N2,TITLE1,SOILPROFILE(I).psoil(J).IUNSET)
				CALL Gen_USER_ESET_EXCA(N1,N2,TITLE1,SOILPROFILE(I).psoil(J).IUESET)							 
			else
				print *, "The soilprofile(i).psoil(j) has no corresponding elements.i,j=",i,j
					
			endif

        end do
		
		
    enddo
        

	do IAC1=1,naction
		if(action(IAC1).NDIM==1) then
			
			N1=action(IAC1).KPOINT(1)
			N2=action(IAC1).KPOINT(ACTION(IAC1).NKP)
			TITLE1=ACTION(IAC1).TITLE
			CALL Gen_USER_NSET_EXCA(N1,N2,TITLE1,ACTION(IAC1).IUNSET)
			CALL Gen_USER_ESET_EXCA(N1,N2,TITLE1,ACTION(IAC1).IUESET)
			NNODE1=UNSET(ACTION(IAC1).IUNSET).NNUM
			
		elseif(action(IAC1).NDIM==0) then
			
			DO K=1,action(IAC1).NKP
				action(IAC1).node(k)=kpnode(action(IAC1).kpoint(k))
			ENDDO	
			NNODE1=action(IAC1).NKP
			
		endif 
		
		select case(action(IAC1).TYPE) 
			case(2) !生成弹簧单元
				NEL1=NNODE1
				CALL enlarge_element(ELEMENT,ENUM,NEL1,N1)
				action(iac1).istiffelement=n1
				action(iac1).nstiffelement=nel1
				allocate(action(iac1).node2stiffelement(nnum))
				action(iac1).node2stiffelement=0
				
				element(N1:enum).sf=action(iac1).sf
				
				DO K=0,NEL1-1
					ALLOCATE(ELEMENT(N1+K).NODE(1))
					IF(ACTION(IAC1).NDIM==1) THEN
						ELEMENT(N1+K).NODE(1)=UNSET(ACTION(IAC1).iunset).NODE(K+1)
						
					ELSEIF(ACTION(IAC1).NDIM==0) THEN
						ELEMENT(N1+K).NODE(1)=ACTION(IAC1).NODE(K+1)
					ENDIF
					action(iac1).node2stiffelement(ELEMENT(N1+K).NODE(1))=n1+k	
				ENDDO
				
					
				!ALLOCATE(ACTION(IAC1).STIFFELEMENT(NEL1))
				!DO K=1,NEL1
				!	ACTION(IAC1).STIFFELEMENT(K)=N1+K-1
				!ENDDO
				IF(ACTION(IAC1).DOF<=3) THEN					
					et1=SPRINGX-1+ACTION(IAC1).DOF
				ELSE
					et1=SPRINGX+ACTION(IAC1).DOF
				ENDIF
				call ettonnum(et1,nnum1,ndof1,ngp1,nd1,stype,ec1)
				element(N1:).nnum=1
				element(N1:).et=et1
				element(N1:).ndof=ndof1
				element(N1:).ngp=ngp1
				element(N1:).nd=nd1
				element(N1:).ec=ec1
				
			case(1) !生成位移边界

				CALL enlarge_bc(BC_DISP,BD_NUM,NNODE1,N1)
				action(iac1).nbcnode=nnode1
				action(iac1).ibcnode=n1
				allocate(action(iac1).node2bcdisp(nnum))
				!BC_DISP.SF=
				IF(ACTION(IAC1).NDIM==1) THEN
					bc_disp(n1:).node=UNSET(ACTION(IAC1).iunset).NODE						
				ELSEIF(ACTION(IAC1).NDIM==0) THEN
					bc_disp(n1:).node=ACTION(IAC1).NODE
				ENDIF                    
				bc_disp(n1:).dof=action(iac1).dof
				
				!bc_disp(n1:).sf=action(iac1).sf
				
				forall (k=n1:bd_num) action(iac1).node2bcdisp(bc_disp(k).node)=k	
			case(0) !生成力边界

				CALL enlarge_bc(BC_LOAD,BL_NUM,NNODE1,N1)
				action(iac1).iloadnode=n1
				action(iac1).nloadnode=nnode1
				allocate(action(iac1).node2bcload(nnum))
				action(iac1).node2bcload=0
				IF(ACTION(IAC1).NDIM==1) THEN
					bc_LOAD(n1:).node=UNSET(ACTION(IAC1).iunset).NODE
				ELSEIF(ACTION(IAC1).NDIM==0) THEN
					bc_LOAD(n1:).node=ACTION(IAC1).NODE						
				ENDIF 
				forall (k=n1:bL_num) action(iac1).node2bcload(bc_load(k).node)=k
				BC_LOAD(N1:).dof=action(iac1).dof
				!BC_LOAD(N1:).sf=action(iac1).sf

		end select	
	enddo
	
	
	
	CALL enlarge_element(ELEMENT,ENUM,NSTRUT,N1)	
            
	DO J=1,NSTRUT
		
				
		IF(STRUT(J).ISBAR==0) THEN
			ET1=SPRINGX
			ALLOCATE(ELEMENT(N1+J-1).NODE(1))
			ELEMENT(N1+J-1).NODE(1)=KPNODE(STRUT(J).Z(1))
		ELSE
			ET1=BAR2D
			ALLOCATE(ELEMENT(N1+J-1).NODE(2))
			ELEMENT(N1+J-1).NODE=KPNODE(STRUT(J).Z)
		ENDIF
		ELEMENT(N1+J-1).SF=STRUT(J).SF
		ELEMENT(N1+J-1).MAT=STRUT(J).MAT
		STRUT(J).ELEMENT=N1+J-1

		call ettonnum(et1,nnum1,ndof1,ngp1,nd1,stype,ec1)
		element(N1+J-1).nnum=NNUM1
		element(N1+J-1).et=et1
		element(N1+J-1).ndof=ndof1
		element(N1+J-1).ngp=ngp1
		element(N1+J-1).nd=nd1
		element(N1+J-1).ec=ec1				
				
	ENDDO				

	
	
    call checkdata()
   

endsubroutine



subroutine Gen_USER_NSET_EXCA(N1,N2,NAME,ISET)
	use solverDS
	use excaDS
	integer,intent(in)::N1,N2
	CHARACTER(32),INTENT(IN)::NAME
	integer,intent(out)::ISET
	INTEGER::I
	
	ISET=NUNSET+1
	UNSET(ISET).NAME=NAME
	UNSET(ISET).NNUM=KPNODE(N2)-KPNODE(N1)+1
	ALLOCATE(UNSET(ISET).NODE(UNSET(ISET).NNUM))
	DO I=1,UNSET(ISET).NNUM
		UNSET(ISET).NODE(I)=KPNODE(N1)+I-1
	ENDDO
	NUNSET=NUNSET+1
	
endsubroutine

subroutine Gen_USER_ESET_EXCA(N1,N2,NAME,ISET)
	use solverDS
	use excaDS
	integer,intent(in)::N1,N2
	CHARACTER(32),INTENT(IN)::NAME
	integer,intent(out)::ISET
	INTEGER::I
	
	ISET=NUESET+1
	UESET(ISET).NAME=NAME
	UESET(ISET).ENUM=KPELEMENT(1,N2)-KPELEMENT(2,N1)+1
    IF(KPELEMENT(1,N2)<KPELEMENT(2,N1)) THEN
        PAUSE
    ENDIF
	ALLOCATE(UESET(ISET).ELEMENT(UESET(ISET).ENUM))
	DO I=1,UESET(ISET).ENUM
		UESET(ISET).ELEMENT(I)=KPELEMENT(2,N1)+I-1
	ENDDO
	NUESET=NUESET+1
	
endsubroutine

subroutine divideLine2D(xi,yi,si,xj,yj,sj,NODE1,IN1,NNODE1)
	use excaDS
	implicit none
    INTEGER,INTENT(IN OUT)::IN1,NNODE1
	REAL(DOUBLE),INTENT(IN)::xi,yi,si,xj,yj,sj
	real(double),INTENT(OUT)::node1(IN1,NNODE1)

	REAL(DOUBLE)::S1,L1,cos1,sin1,w1,w2
	
	
	
	nnode1=0
	nnode1=nnode1+1
	node1(1,nnode1)=xi
	node1(2,nnode1)=yi
	node1(3,nnode1)=si
	L1=((XI-XJ)**2+(YI-YJ)**2)**0.5
    IF(ABS(L1)<1E-7) THEN
        STOP "I(X,Y)=J(X,Y).ERROR IN SUB divideLine2D"
    END IF
    COS1=(XJ-XI)/L1
    SIN1=(YJ-YI)/L1
	!s1=2*si*sj/(si+sj)
	s1=1.05*max(si,sj)
	do while(L1>S1)
		nnode1=nnode1+1
		!fac1=node1(3,node1-1)/(sj+node1(3,node1-1))
		node1(1,nnode1)=node1(1,nnode1-1)+node1(3,nnode1-1)*cos1
		node1(2,nnode1)=node1(2,nnode1-1)+node1(3,nnode1-1)*sin1
        
		L1=L1-node1(3,nnode1-1)
		w1=1./node1(3,nnode1-1)
		w2=1./L1		
		node1(3,nnode1)=node1(3,nnode1-1)*w1/(w1+w2)+sj*w2/(w1+w2)
		S1=1.05*max(node1(3,nnode1),sj)		
    enddo
    nnode1=nnode1+1
	node1(1,nnode1)=xj
	node1(2,nnode1)=yj
	node1(3,nnode1)=sj
	
	
end subroutine


subroutine soilstiffness_cal(isp,istep,iz,soil) !iz,为开挖面高程点。
	use solverDS
	use excaDS
	implicit none

	integer,intent(in)::isp,istep,iz
	type(soillayer_tydef),intent(in out)::soil
	real(double)::fi1,C1,m1,z1,dz1,Xdis1,b1,d1,Es1,Ep1,mu1,Iz1,Gs1
	integer::i
	
	Es1=matproperty(soil.mat,4,istep)
	mu1=matproperty(soil.mat,5,istep)
	Ep1=matproperty(element(ueset(soil.iueset).element(1)).mat,1,istep) !假定同一土层内的支护桩的参数保持一致。
	Iz1=matproperty(element(ueset(soil.iueset).element(1)).mat,5,istep) !假定同一土层内的支护桩的参数保持一致。
	d1=matproperty(element(ueset(soil.iueset).element(1)).mat,7,istep) !假定同一土层内的支护桩的参数保持一致。
	b1=matproperty(element(ueset(soil.iueset).element(1)).mat,8,istep) !假定同一土层内的支护桩的参数保持一致。
	
	!考虑到支护桩为平面应变状态，每根桩所辖的范围应为桩距，所以，计算桩径应为桩距。			
	if(d1<=1.0d0) then
		d1=0.9*(1.5*d1+0.5)
	else
		d1=0.9*(d1+1)
	endif
	b1=min(d1,b1)	
	
	
	select case(soilprofile(isp).kmethod)
		case(0) !m法，ks的单位为:kN/m3
			fi1=matproperty(soil.mat,2,istep) 
			c1=matproperty(soil.mat,1,istep) !kPa			
			m1=0.2*(fi1**2-fi1+c1)/10.d0 !MN/m4
			do i=1,2
				dz1=kpoint(ndimension,iz)-kpoint(ndimension,soil.z(i)) !m
				soil.ks(i)=m1*1000*dz1
			enddo
			
		case(1) !E.V法 ks的单位为:kN/m3	
			soil.ks=Es1/(1-mu1**2)/0.88/b1
			
		case(2) !zhubitang ks的单位为:kN/m3
			Gs1=Es1/(2*(1+mu1))*(1+0.75*mu1)
			soil.ks=5.46*Gs1*(Ep1/Gs1)**(-1.d0/7.0d0)/b1
			
		case(3) !Biot ks的单位为:kN/m3
			soil.ks=0.95*Es1/(1-mu1**2)*(b1**4*Es1/((1-mu1**2)*Ep1*Iz1))**0.108
			soil.ks=soil.ks/b1
		case(4) !Vesic ks的单位为:kN/m3
			soil.ks=0.65*Es1/(1-mu1**2)*(b1**4*Es1/((1-mu1**2)*Ep1*Iz1))**(1./12.0)	
			soil.ks=soil.ks/b1
			
		case(-1) !直接输入，ks的单位为:kN/m3
			soil.ks=matproperty(soil.mat,7,istep)
		
			
			
	endselect
	
endsubroutine


subroutine soilspring_EXCA(istep)
	use solverDS
	use excaDS
	implicit none
	integer,intent(in)::istep
	integer::i,j
	
	do i=1,nsoilprofile 
		element(pile(soilprofile(i).beam).element_sp(1,:)).isactive=0
        element(pile(soilprofile(i).beam).element_sp(2,:)).isactive=0
        do j=1,6
		    element(pile(soilprofile(i).beam).element_sp(1,:)).property(j)=0.d0
		    element(pile(soilprofile(i).beam).element_sp(2,:)).property(j)=0.d0
        end do
		
		do j=1,soilprofile(i).nasoil
			call soilstiffness_cal(i,istep,soilprofile(i).za,soilprofile(i).asoil(j))
			call Initialize_soilspringEelement_EXCA(soilprofile(i).beam,istep,soilprofile(i).aside,1,soilprofile(i).asoil(j))
        enddo
        
		do j=1,soilprofile(i).npsoil
			call soilstiffness_cal(i,istep,soilprofile(i).zp,soilprofile(i).psoil(j))
			call Initialize_soilspringEelement_EXCA(soilprofile(i).beam,istep,soilprofile(i).aside,2,soilprofile(i).psoil(j))
		enddo
	end do
endsubroutine

subroutine Initialize_soilspringEelement_EXCA(ipile,istep,aside,aop,soil)
	use excaDS
	use solverDS
	implicit none
	integer,intent(in)::ipile,istep,aside,aop
	type(soillayer_tydef),intent(in)::soil	
	integer::i,j,ieL1,inode1(2),sign1,iel2(2)
	real(double)::zi1(2),vi1(2),t1,b1,d1,bp1
	real(double),external::interpolation
	

	
	if(sf(soil.sf).factor(istep)/=0) then 
	
	
		sign1=aside
		if(aop==2) sign1=-aside
		
	
		do i=1,ueset(soil.iueset).enum
			ieL1=ueset(soil.iueset).element(i)
			d1=matproperty(element(ieL1).mat,7,istep) !假定同一土层内的支护桩的参数保持一致。
			b1=matproperty(element(ieL1).mat,8,istep) !假定同一土层内的支护桩的参数保持一致。
			
			!考虑到支护桩为平面应变状态，每根桩所辖的范围应为桩距，所以，计算桩径应为桩距。			
			if(d1<=1.0d0) then
				d1=0.9*(1.5*d1+0.5)
			else
				d1=0.9*(d1+1)
			endif
			bp1=min(d1,b1)	

			
			inode1=element(ieL1).node
			if(element(ieL1).ifreedof>0) then
				if(inode1(1)==freedof(element(ieL1).ifreedof).newnode) inode1(1)=freedof(element(ieL1).ifreedof).node
				if(inode1(2)==freedof(element(ieL1).ifreedof).newnode) inode1(2)=freedof(element(ieL1).ifreedof).node
			endif
			iel2=pile(ipile).element_sp(aop,inode1)
			element(iel2).isactive=1
			
			element(iel2).mat=-aop
			element(iel2).sign=sign1
			element(iel2).sf=soil.sf
			
			zi1=node(inode1).coord(ndimension)
			t1=0
			do j=1,ndimension
				t1=t1+(node(element(ieL1).node(1)).coord(j)-node(element(ieL1).node(2)).coord(j))**2
			enddo
			t1=t1**0.5/2.0
			!ko,ka,kp,ks,L,Pw
			element(iel2).property(5)=element(iel2).property(5)+t1
			vi1(1)=interpolation(kpoint(ndimension,soil.z(1:2)),soil.sigmaKo,2,zi1(1))
			vi1(2)=interpolation(kpoint(ndimension,soil.z(1:2)),soil.sigmaKo,2,zi1(2))
			element(iel2).property(1)=element(iel2).property(1)+b1*(vi1(1)+vi1(2))*t1/2
			vi1(1)=interpolation(kpoint(ndimension,soil.z(1:2)),soil.sigmaKa,2,zi1(1))
			vi1(2)=interpolation(kpoint(ndimension,soil.z(1:2)),soil.sigmaKa,2,zi1(2))
			element(iel2).property(2)=element(iel2).property(2)+sign1*b1*(max(sign1*vi1(1),0.0)+max(sign1*vi1(2),0.0))*t1/2			
			vi1(1)=interpolation(kpoint(ndimension,soil.z(1:2)),soil.sigmaKp,2,zi1(1))
			vi1(2)=interpolation(kpoint(ndimension,soil.z(1:2)),soil.sigmaKp,2,zi1(2))
			element(iel2).property(3)=element(iel2).property(3)+bp1*(vi1(1)+vi1(2))*t1/2
			vi1(1)=interpolation(kpoint(ndimension,soil.z(1:2)),soil.Ks,2,zi1(1))
			vi1(2)=interpolation(kpoint(ndimension,soil.z(1:2)),soil.Ks,2,zi1(2))
			element(iel2).property(4)=element(iel2).property(4)+bp1*(vi1(1)+vi1(2))*t1/2
			vi1(1)=interpolation(kpoint(ndimension,soil.z(1:2)),soil.Pw,2,zi1(1))
			vi1(2)=interpolation(kpoint(ndimension,soil.z(1:2)),soil.Pw,2,zi1(2))
			element(iel2).property(6)=element(iel2).property(6)+b1*sign1*(vi1(1)+vi1(2))*t1/2
			
			if(istep==1) then
			!	!单元初力=ko+pw
                do j=1,2
				!if(.not.allocated(element(iel2(j)).gforce)) allocate(element(iel2(j)).gforce(1))
				if(.not.allocated(element(iel2(j)).km)) allocate(element(iel2(j)).km(1,1))
                element(iel2(j)).km(1,1)=element(iel2(j)).property(4)
				!element(iel2(j)).gforce(1) =element(iel2(j)).property(1)+element(iel2(j)).property(6)
			    enddo
            endif
			
		end do

	endif
	
	

endsubroutine

function interpolation(x,y,nx,xi)
!x,y must be in order.
	implicit none
	INTEGER,PARAMETER::Double=KIND(1.0D0)
    integer,intent(in)::nx
	real(double),intent(in)::x(nx),y(nx),xi
	real(double)::interpolation,t1
	integer::i
    

    if(nx==1) then
       interpolation=y(1)
       return
    endif
    do i=1,nx-1
        if((xi<=x(i+1).and.xi>=x(i)).or.(xi<=x(i).and.xi>=x(i+1))) then
	        t1=x(i+1)-x(i)
	        if(abs(t1)<1e-7) then
		        print *, "Warning! 分母=0,function=Interpolation()"
		        interpolation=(y(i)+y(i+1))/2.0d0
	        else
		        interpolation=(y(i+1)-y(i))/(t1)*(xi-x(i))+y(i)
            endif
            return
        endif
    enddo
    if(i==nx) then
        stop "xi is out side the range.function=Interpolation()"
    endif
    
    
endfunction

subroutine EarthPressure(istep)
	use solverds
	use ExcaDS
	implicit none
    integer,intent(in)::istep
    
	integer::i,j,k,ilayer1=0
	real(double)::fi1,gamma1,t1,WL1
	logical::tof1
	
	do i=1,nsoilprofile
		
		!water pressure calculation
		call waterpressure_Cal(i,istep)
		ilayer1=0
		tof1=.true.
		do j=1,soilprofile(i).nasoil
			if(sf(soilprofile(i).asoil(j).sf).factor(istep)==0) cycle
			ilayer1=ilayer1+1
			!竖向有效压力（如果是合算，则为总应力）
			if(tof1) then
                WL1=soilprofile(i).awL*sf(soilprofile(i).sf_awL).factor(istep)
                T1=MAX((wL1-kpoint(ndimension,soilprofile(i).asoil(j).z(1)))*GA,0.0)
				soilprofile(i).asoil(j).SigmaVT(1)=T1+soilprofile(i).aLoad*sf(soilprofile(i).sf_aload).factor(istep)
                soilprofile(i).asoil(j).SigmaV(1)=soilprofile(i).asoil(j).SigmaVT(1)-soilprofile(i).asoil(j).pw(1)				
				soilprofile(i).za=soilprofile(i).asoil(j).z(1)
				tof1=.false.
			else
				soilprofile(i).asoil(j).SigmaV(1)=soilprofile(i).asoil(j-1).SigmaV(2)
                soilprofile(i).asoil(j).SigmaVT(1)=soilprofile(i).asoil(j-1).SigmaVT(2)
			endif
			gamma1=matproperty(soilprofile(i).asoil(j).mat,3,istep)
			t1=kpoint(ndimension,soilprofile(i).asoil(j).z(1))-kpoint(ndimension,soilprofile(i).asoil(j).z(2))
            soilprofile(i).asoil(j).SigmaVT(2)=soilprofile(i).asoil(j).SigmaVT(1)+gamma1*(t1)
            soilprofile(i).asoil(j).SigmaV(2)=soilprofile(i).asoil(j).SigmaVT(2)-soilprofile(i).asoil(j).pw(2)
			
			!fi1=material(soilprofile(i).asoil(j).mat).property(2)/180.d0*PI()*SF(material(soilprofile(i).asoil(j).mat).sf(2)).factor(istep)
			
			!土压力系数
			call ep_Cal(i,istep,soilprofile(i).asoil(j),soilprofile(i).aside)									
												
			
		end do
		tof1=.true.
		do j=1,soilprofile(i).npsoil
			if(sf(soilprofile(i).psoil(j).sf).factor(istep)==0) cycle
			
			!竖向有效压力（如果是合算，则为总应力）
			
			if(tof1) then
                WL1=soilprofile(i).PwL*sf(soilprofile(i).sf_PwL).factor(istep)
                T1=MAX((wL1-kpoint(ndimension,soilprofile(i).Psoil(j).z(1)))*GA,0.d0)
				soilprofile(i).psoil(j).SigmaVT(1)=T1+soilprofile(i).pLoad*sf(soilprofile(i).sf_pload).factor(istep)
				soilprofile(i).psoil(j).SigmaV(1)=soilprofile(i).psoil(j).SigmaVT(1)-soilprofile(i).psoil(j).pw(1)
                soilprofile(i).zp=soilprofile(i).psoil(j).z(1)
				tof1=.false.
			else
				soilprofile(i).psoil(j).SigmaV(1)=soilprofile(i).psoil(j-1).SigmaV(2)
                soilprofile(i).psoil(j).SigmaVT(1)=soilprofile(i).psoil(j-1).SigmaVT(2)
			endif
			gamma1=matproperty(soilprofile(i).psoil(j).mat,3,istep)
            t1=kpoint(ndimension,soilprofile(i).psoil(j).z(1))-kpoint(ndimension,soilprofile(i).psoil(j).z(2))
			soilprofile(i).psoil(j).SigmaVT(2)=soilprofile(i).psoil(j).SigmaVT(1)+gamma1*(t1)
            soilprofile(i).psoil(j).SigmaV(2)=soilprofile(i).psoil(j).SigmaVT(2)-soilprofile(i).psoil(j).pw(2)
			
			!fi1=material(soilprofile(i).psoil(j).mat).property(2)/180.d0*PI()*SF(material(soilprofile(i).psoil(j).mat).sf(2)).factor(istep)
			
			!土压力系数
			call ep_Cal(i,istep,soilprofile(i).psoil(j),-soilprofile(i).aside)						
												
			
		end do	
		
	end do
	
10	format("THE SOIL PRESSURE INFO IS LISTED AS FOLLOWS:\N"C, &
		   "NO",13X,"Layer",10X,"ISTEP",10X,"ISOILPROFILE",3X,"SIDE",11X,"Z(L)",X,"Sv(F/L^2)",6X,"Sko(F/L^2)",5X,"Ska(F/L^2)",5X,"Skp(F/L^2)",5X,"PW(F/L^2)",6X,"Ks(F/L^3)",6X)	
20	format(5(I14,X),7(E14.7,X))
end subroutine

subroutine ep_Cal(isp,istep,soil,sign)
	use solverDS
	use excaDS
	implicit none
	integer,intent(in)::isp,istep,sign
	type(soillayer_tydef),intent(in out)::soil
	real(double)::fi1,C1,delta1=0,CDelta1=0,t1
	
	fi1=matproperty(soil.mat,2,istep)/180*PI()
	C1=matproperty(soil.mat,1,istep)
	delta1=matproperty(soil.mat,8,istep)/180*PI()
	if(abs(fi1)<1e-7) then
		CDelta1=0.d0
	else
		CDelta1=asin(sin(delta1)/sin(fi1))
	endif
	
	select case(soilprofile(isp).spm)
	
		case(0) !郎肯
			soil.ko=1-sin(fi1)
            
			soil.ka=(1-sin(fi1)*cos(CDelta1-delta1))/(1+sin(fi1))*exp(-(CDelta1-delta1)*tan(fi1))
			soil.kp=(1+sin(fi1)*cos(CDelta1+delta1))/(1-sin(fi1))*exp((CDelta1+delta1)*tan(fi1))
			soil.sigmaKo=soil.Ko*soil.sigmaV*sign
			soil.sigmaKa=(soil.Ka*soil.sigmaV-2*C1*soil.Ka**0.5)*sign
			soil.sigmaKp=(soil.Kp*soil.sigmaV+2*C1*soil.Kp**0.5)*sign
			
		case default
		
	end select
endsubroutine

subroutine waterpressure_cal(isoilprofile,istep)
	use solverds
	use excaDS
	implicit none
	integer,intent(in)::isoilprofile,istep
	
	integer::i,j,k
	
	real(double)::R1,R2,awL1,pwL1,h1,t1
	
	awL1=soilprofile(isoilprofile).awl*sf(soilprofile(isoilprofile).sf_awl).factor(istep)
	pwL1=soilprofile(isoilprofile).pwl*sf(soilprofile(isoilprofile).sf_pwl).factor(istep)
	
	!求总水阻
	R2=0.0d0
	do i=1,soilprofile(isoilprofile).nasoil
		if(sf(soilprofile(isoilprofile).asoil(i).sf).factor(istep)==0) cycle
		!简化考虑渗透力，墙身不透水，墙底两侧水压力相等。
		if(kpoint(ndimension,soilprofile(isoilprofile).asoil(i).z(2))<awl1) then!水位处必须分层
			!等效水阻
			!如果水面高于地表，将水等效为为土层（令c,phi,k,E均为0）
			if(matproperty(soilprofile(isoilprofile).asoil(i).mat,6,istep)>1e-10) then
				R2=R2+(kpoint(ndimension,soilprofile(isoilprofile).asoil(i).z(1))-kpoint(ndimension,soilprofile(isoilprofile).asoil(i).z(2))) &
						/matproperty(soilprofile(isoilprofile).asoil(i).mat,6,istep)	
			endif
		endif		
	enddo	
	do i=1,soilprofile(isoilprofile).npsoil
		if(sf(soilprofile(isoilprofile).psoil(i).sf).factor(istep)==0) cycle
		!简化考虑渗透力，墙身不透水，墙底两侧水压力相等。
		if(kpoint(ndimension,soilprofile(isoilprofile).psoil(i).z(2))<pwl1) then!水位处必须分层
			!等效水阻
			if(matproperty(soilprofile(isoilprofile).psoil(i).mat,6,istep)>1e-10) then
				R2=R2+(kpoint(ndimension,soilprofile(isoilprofile).psoil(i).z(1))-kpoint(ndimension,soilprofile(isoilprofile).psoil(i).z(2))) &
						/matproperty(soilprofile(isoilprofile).psoil(i).mat,6,istep)
			endif
		endif		
	end do	
	
	R1=0.0d0
	do i=1,soilprofile(isoilprofile).nasoil
		if(sf(soilprofile(isoilprofile).asoil(i).sf).factor(istep)==0) cycle
		
		select case(soilprofile(isoilprofile).asoil(i).wpflag) 
			case(1) !常规分算，不考虑渗透力
				soilprofile(isoilprofile).asoil(i).PW(1)=max(awL1-kpoint(ndimension,soilprofile(isoilprofile).asoil(i).z(1)),0.d0)*GA
				soilprofile(isoilprofile).asoil(i).PW(2)=max(awL1-kpoint(ndimension,soilprofile(isoilprofile).asoil(i).z(2)),0.d0)*GA   
				
			case(2) !简化考虑渗透力，墙身不透水，墙底两侧水压力相等。
				if(kpoint(ndimension,soilprofile(isoilprofile).asoil(i).z(2))<awl1) then!水位处必须分层
					!等效水阻
					h1=awl1-(awl1-pwl1)*r1/r2
					soilprofile(isoilprofile).asoil(i).pw(1)=(h1-kpoint(ndimension,soilprofile(isoilprofile).asoil(i).z(1)))*GA
					if(matproperty(soilprofile(isoilprofile).asoil(i).mat,6,istep)>1e-10) then
					R1=R1+(kpoint(ndimension,soilprofile(isoilprofile).asoil(i).z(1))-kpoint(ndimension,soilprofile(isoilprofile).asoil(i).z(2))) &
						/matproperty(soilprofile(isoilprofile).asoil(i).mat,6,istep)
					endif
						
					h1=awl1-(awl1-pwl1)*r1/r2
					soilprofile(isoilprofile).asoil(i).pw(2)=(h1-kpoint(ndimension,soilprofile(isoilprofile).asoil(i).z(2)))*GA
					
				endif
			case default
				soilprofile(isoilprofile).asoil(i).pw=0.0d0
		endselect		
	enddo
	
	do i=soilprofile(isoilprofile).npsoil,1,-1
		if(sf(soilprofile(isoilprofile).psoil(i).sf).factor(istep)==0) cycle
		
		select case(soilprofile(isoilprofile).psoil(i).wpflag) 
			case(1) !常规分算，不考虑渗透力
				soilprofile(isoilprofile).psoil(i).PW(1)=max(pwL1-kpoint(ndimension,soilprofile(isoilprofile).psoil(i).z(1)),0.d0)
				soilprofile(isoilprofile).psoil(i).PW(2)=max(pwL1-kpoint(ndimension,soilprofile(isoilprofile).psoil(i).z(2)),0.d0)   
				
			case(2) !简化考虑渗透力，墙身不透水，墙底两侧水压力相等。
				if(kpoint(ndimension,soilprofile(isoilprofile).psoil(i).z(2))<pwl1) then!水位处必须分层
					
					h1=awl1-(awl1-pwl1)*r1/r2
					soilprofile(isoilprofile).psoil(i).pw(2)=(h1-kpoint(ndimension,soilprofile(isoilprofile).psoil(i).z(2)))*GA
					!等效水阻
					if(matproperty(soilprofile(isoilprofile).psoil(i).mat,6,istep)>1e-10) then
					R1=R1+(kpoint(ndimension,soilprofile(isoilprofile).psoil(i).z(1))-kpoint(ndimension,soilprofile(isoilprofile).psoil(i).z(2))) &
						/matproperty(soilprofile(isoilprofile).psoil(i).mat,6,istep)
					endif	
					h1=awl1-(awl1-pwl1)*r1/r2
					soilprofile(isoilprofile).psoil(i).pw(1)=(h1-kpoint(ndimension,soilprofile(isoilprofile).psoil(i).z(1)))*GA
					
				endif
			case default
				soilprofile(isoilprofile).psoil(i).pw=0.0d0
		endselect		
	enddo

endsubroutine

!判断在arr_t中是否还其它点在线段xi,xj上，并形成数组ar,使ar包含vertex中在线段xi,xj上的所有点
!以n4返回该数组的大小，xi,xj为该数组的两端,之间的点按顺序排列。
subroutine aovis2d(N1,N2,IAR,NAR) !any other vertex inside the segment?
  use solverDS
  use ExcaDS
  !use ds_t
  implicit none
  integer,intent(in out)::NAR
  integer,intent(in)::N1,N2
  integer,intent(out)::IAR(NAR)
  real(double)::t1,t2
  integer::i,j,k
  logical::tof1,tof2
  
  NAR=0
  IAR=0
  NAR=NAR+1
  IAR(NAR)=N1
  do i=1,nkp
	 t1=(kpoint(1,i)-KPOINT(1,N1))**2+(kpoint(ndimension,i)-KPOINT(NDIMENSION,N1))**2
	 if(t1<1e-7) cycle
	 t1=(kpoint(1,i)-KPOINT(1,N2))**2+(kpoint(ndimension,i)-KPOINT(NDIMENSION,N2))**2
	 if(t1<1e-7) cycle
	 call vins2d(KPOINT(1,N1),KPOINT(NDIMENSION,N1),KPOINT(1,N2),KPOINT(NDIMENSION,N2),kpoint(1,i),kpoint(ndimension,i),tof1)
	 if(tof1) then
		tof2=.true.
		do j=1,NAR
			if(abs(KPOINT(1,IAR(J))-kpoint(1,i))>1e-6) cycle
			if(abs(KPOINT(NDIMENSION,IAR(J))-kpoint(ndimension,i))>1e-6) cycle
			tof2=.false.
			EXIT
		end do
		if(tof2) then
			NAR=NAR+1
			if(NAR>size(IAR)) then
			   print *,'sub aovis(),NAR>',size(IAR)
			   stop
			end if
			IAR(NAR)=I 
			!ar(NAR)=kpoint(ndimension,i)
			!ar(3,NAR)=arr_t(i).s	
		end if	   
	 end if
  end do
  NAR=NAR+1
  IAR(NAR)=N2
  !ar(3,NAR)=sj

  !按离ar(:,1)的由近而远的顺序重排
  do i=2,NAR-1
	 t1=(KPOINT(1,IAR(I))-KPOINT(1,IAR(1)))**2+(KPOINT(NDIMENSION,IAR(I))-KPOINT(NDIMENSION,IAR(1)))**2
	 do j=i+1,NAR-1
		t2=(KPOINT(1,IAR(J))-KPOINT(1,IAR(1)))**2+(KPOINT(NDIMENSION,IAR(J))-KPOINT(NDIMENSION,IAR(1)))**2
		if(t2<t1) then
		   K=IAR(J)
           IAR(J)=IAR(I)
           IAR(I)=K
		   t1=t2
		end if
	 end do
  end do

end subroutine

  !判断点x(),是否在线段xixj上。如果是，tof1=.t. or .f.
subroutine vins2d(xi,yi,xj,yj,x,y,tof1) !vertex in the segment or noe
  !use meshds
  implicit none
  logical::tof1
  real(8)::xi,yi,xj,yj,x,y,t1,t2,t3
  
  tof1=.false.
  if(x<min(xi,xj))  return  
  if(x>max(xi,xj))  return
  if(y<min(yi,yj))  return  
  if(y>max(yi,yj))  return

  
  t3=(xi-x)*(yj-y)-(xj-x)*(yi-y)
  if(abs(t3)>1e-7) return
  
  tof1=.true.	 

end subroutine


subroutine Beam_Result_EXCA(istep)
	use solverlib
	use excaDS
	implicit none
    integer,intent(in)::istep
	integer::i,j,k,iel1=0,nnode1(2)=0,ipile=0,n1,n2,nc1=0,NC2=0,iesp1(2)=0
    logical,allocatable::isout1(:)
	real(double)::t1=0
	
	if(istep==1) then
		open(unit=23,file=EXCAB_BEAMRES_FILE,status='replace')
		write(23,10)
		open(unit=24,file=EXCAB_STRURES_FILE,status='replace')
		write(24,20)		
	else
		open(unit=23,file=EXCAB_BEAMRES_FILE,status='old',access='append')
		open(unit=24,file=EXCAB_STRURES_FILE,status='old',access='append')
	endif
	
	do iPILE=1,NPILE
		!ipile=soilprofile(i).beam
		
        n1=Lbound(pile(ipile).beamresult,dim=1)
        n2=Ubound(pile(ipile).beamresult,dim=1)
        if(allocated(isout1)) deallocate(isout1)
        allocate(isout1(n1:n2))
        isout1=.false.
		do j=1,pile(ipile).nel
			
			iel1=pile(IPILE).element(j)
			nnode1=element(iel1).node
			if(element(ieL1).ifreedof>0) then
				if(nnode1(1)==freedof(element(ieL1).ifreedof).newnode) nnode1(1)=freedof(element(ieL1).ifreedof).node
				if(nnode1(2)==freedof(element(ieL1).ifreedof).newnode) nnode1(2)=freedof(element(ieL1).ifreedof).node
			endif            
			do k=1,2
				!displacement
                IF(ISEXCA2D==1) THEN
                    NC1=1
                    NC2=2
                ELSE
                    NC1=2
                    NC2=2
                ENDIF
				pile(ipile).beamresult(nnode1(k),1,istep)=Tstepdis(element(iel1).g(NC1+3*(k-1)),istep)
				!M
				pile(ipile).beamresult(nnode1(k),2,istep)=pile(ipile).beamresult(nnode1(k),2,istep)+ &
															-element(iel1).gforceILS(3+3*(k-1))*0.5d0
				!Q
				pile(ipile).beamresult(nnode1(k),3,istep)=pile(ipile).beamresult(nnode1(k),3,istep)+ &
															-element(iel1).gforceILS(NC2+3*(k-1))*0.5d0
				
                
                if(isout1(nnode1(k))) cycle
                
                if(nsoilprofile>0) then
                    iesp1(:)=pile(ipile).element_sp(:,nnode1(k))
				
                    if(element(iesp1(1)).isactive==1) then
				        !主动侧土荷载
				        pile(ipile).beamresult(nnode1(k),5,istep)=pile(ipile).beamresult(nnode1(k),5,istep) &
					        -(element(iesp1(1)).property(2)+element(iesp1(1)).property(6))/ &
					        element(iesp1(1)).property(5)
				        !主动侧土弹簧力									
				        pile(ipile).beamresult(nnode1(k),6,istep)=pile(ipile).beamresult(nnode1(k),6,istep) &
					        +element(iesp1(1)).gforceILS(1)/ &
                            element(iesp1(1)).property(5)  
                        !主动侧土弹簧力限值		
				        pile(ipile).beamresult(nnode1(k),7,istep)=pile(ipile).beamresult(nnode1(k),7,istep) &
					        -(element(iesp1(1)).property(3)- &
					          element(iesp1(1)).property(2))/ &
                            element(iesp1(1)).property(5) 
                        !主动侧土水压力		
				        pile(ipile).beamresult(nnode1(k),8,istep)=pile(ipile).beamresult(nnode1(k),8,istep) &
					        -element(iesp1(1)).property(6)/element(iesp1(1)).property(5)                      
				        !主动动侧土弹簧刚度
                        IF(pile(ipile).beamresult(nnode1(k),6,istep)<0) THEN
                            T1=-1
                        ELSE
                            T1=1
                        ENDIF
				        pile(ipile).beamresult(nnode1(k),13,istep)=pile(ipile).beamresult(nnode1(k),13,istep)+ &
					        T1*element(iesp1(1)).property(4)/ &
                            element(iesp1(1)).property(5) 
                    endif
                    if(element(iesp1(2)).isactive==1) then
				        !被动侧土荷载	
				        pile(ipile).beamresult(nnode1(k),9,istep)=pile(ipile).beamresult(nnode1(k),9,istep) &
					        -(element(iesp1(2)).property(1)+ &
                            element(iesp1(2)).property(6))/ &
					        element(iesp1(2)).property(5)
				        !被动侧土弹簧力									
				        pile(ipile).beamresult(nnode1(k),10,istep)=pile(ipile).beamresult(nnode1(k),10,istep) &
					        +element(iesp1(2)).gforceILS(1)/ &
                            element(iesp1(2)).property(5)
				        !被动侧土弹簧力限值		
				        pile(ipile).beamresult(nnode1(k),11,istep)=pile(ipile).beamresult(nnode1(k),11,istep) &
					        -(element(iesp1(2)).property(3)- &
					          element(iesp1(2)).property(2))/ &
                            element(iesp1(2)).property(5)
                        !被动侧土水压力		
				        pile(ipile).beamresult(nnode1(k),12,istep)=pile(ipile).beamresult(nnode1(k),12,istep) &
					        -element(iesp1(2)).property(6)/element(iesp1(1)).property(5)                    
				    
                        !被动侧土弹簧力刚度
                        IF(pile(ipile).beamresult(nnode1(k),10,istep)<0) THEN
                            T1=-1
                        ELSE
                            T1=1
                        ENDIF                    
				        pile(ipile).beamresult(nnode1(k),14,istep)=pile(ipile).beamresult(nnode1(k),14,istep)+ &
					        +T1*(element(iesp1(2)).property(4))/ &
                            element(iesp1(2)).property(5)   
                    endif
                
                endif
                !NODAL LOAD
                
                pile(ipile).beamresult(nnode1(k),4,istep)=-TLOAD(NODE(NNODE1(K)).DOF(NC1))/pile(ipile).Nlength(nnode1(k))
                
				isout1(nnode1(k))=.true.
			enddo		
			
		enddo
		
		do j=n1,n2
			write(23,11) IPILE,istep,J,node(j).coord(1:2),pile(ipile).beamresult(j,:,istep)
		enddo
		WRITE(23,12) IPILE,ISTEP
		

        
	enddo
	
	do j=1,nstrut
		nc1=strut(j).element
		IF(ELEMENT(NC1).ISACTIVE==1) THEN
		    write(24,21) istep,j,node(element(nc1).node(1)).coord(ndimension),element(nc1).gforceILS(1)
		    if(strut(j).isbar==1) write(24,21) istep,j,node(element(nc1).node(2)).coord(ndimension),element(nc1).gforceILS(2)
        ENDIF
	enddo
	
	
	close(23)
	CLOSE(24)
	
10 FORMAT(10X,"IPILE",10X,"ISTEP",10X,"INODE",11X,"X(L)",11X,"Y(L)",9X,"DIS(L)",4X,"MOMENT(F.M)",11X,"Q(F)",6X,"TFORCE(F)",3X,"PAA+PWA(F/L)", &
	   5X,"PA_SP(F/L)",X,"MAX_PA_SP(F/L)",7X,"PWA(F/L)",3X,"PAP+PWP(F/L)",5X,"PP_SP(F/L)",X,"MAX_PP_SP(F/L)",7X,"PWP(F/L)",4X,"KSA_SP(F/L)",4X,"KSP_SP(F/L)")	
11 FORMAT(3(X,I14),16(X,E14.7))	
12 FORMAT(2(X,I14))
20 FORMAT(10X,"ISTEP",9X,"ISTRUT",11X,"Y(L)",11X,"P(F)")
21 FORMAT(2(X,I14),2(X,E14.7))	   
endsubroutine

    subroutine enlarge_strut(PROP,NP,EXN,SN)
!扩大PROP数组,同时update总的单元数NP=NP+EXN
!EXN:扩大的单元个数
!SN,:扩容部分的起位
    USE EXCADS
	integer,intent(in)::EXN
    INTEGER,INTENT(IN OUT)::NP
    type(strut_tydef),INTENT(IN OUT),ALLOCATABLE::PROP(:)
	integer,intent(out)::SN
    type(strut_tydef),ALLOCATABLE::PROP1(:)
        
	IF(EXN<1) RETURN
	
	allocate(PROP1(NP+EXN))
	PROP1(1:NP)=PROP
	if(allocated(PROP)) deallocate(PROP)
	SN=NP+1
	NP=NP+EXN	
	allocate(PROP,SOURCE=PROP1)
	deallocate(PROP1)        
        
ENDSUBROUTINE  



  subroutine checkdata()
      use excaDS
      use solverds
	  use ifqwin
      use IFCORE
	  implicit none
      !integer,intent(in)::istep
      LOGICAL  statusmode,tof2,tof1
	  character(1)::butter
      real(8)::xmin,ymin,xmax,ymax
	  common /xy/ xmin,ymin,xmax,ymax
	  real(8)::t1,t2,size
	  type(wxycoord) wxy
	  type(qwinfo) winfo
	  TYPE (windowconfig):: thescreen 
	  common /c2/ thescreen
 	  EXTERNAL  shown,showY,SHOWBW,SHOWBC,SHOWELNUM,SHOWELMAT,SHOWCANCEL,CHECKSOIL,SHOWET
	  EXTERNAL	SHOWNA
      EXTERNAL  checkgraph
	  integer(4)          oldcolor,event,ix,iy,res,result,ikeystate,iunit
	  INTEGER(2)          status
      integer(4)   keystate,i4,ef
	  real(8)::scale
	
	ef=0
!    close(90)
	open(unit=90,file='user',title='checkdata',iostat=ef)
 
    

	i4=-1
	thescreen.numypixels   = -1
	thescreen.numxpixels   = -1
    thescreen.numtextcols  = -1
    thescreen.numtextrows  = -1
    thescreen.numcolors    = -1
	thescreen.fontsize = -1

	
	xmin=minval(kpoint(1,:))
	ymin=minval(kpoint(ndimension,:))
	xmax=maxval(kpoint(1,:))
	ymax=maxval(kpoint(ndimension,:))
    
	if(abs(xmax)<1e-6) xmax=ymax
	if(abs(ymax)<1e-6) ymax=xmax


    statusmode = SETWINDOWCONFIG(thescreen)
    iF(.NOT. statusmode) statusmode = SETWINDOWCONFIG(thescreen)
    statusmode = GETWINDOWCONFIG( thescreen )
	
	winfo%TYPE = QWIN$MAX
    result =SETWSIZEQQ(90, winfo)
	
	result = INSERTMENUQQ (6, 0, $MENUENABLED, 'Check'c, nul)
	result = INSERTMENUQQ (6, 1, $MENUENABLED, 'NODE_Num'c, Shown)
	result = INSERTMENUQQ (6, 2, $MENUENABLED, 'NODE_Y'c, ShowY)
	result = INSERTMENUQQ (6, 3, $MENUENABLED, 'NODE_BW'c, ShowBW)
	result = INSERTMENUQQ (6, 4, $MENUENABLED, 'NODE_BC'c, showbc)
	result = INSERTMENUQQ (6, 5, $MENUENABLED, 'ELEMENT_NUM'c, SHOWELNUM)
	result = INSERTMENUQQ (6, 6, $MENUENABLED, 'ELEMENT_MAT'c, SHOWELMAT)
	result = INSERTMENUQQ (6, 7, $MENUENABLED, 'ELEMENT_TYPE'c, SHOWET)
	result = INSERTMENUQQ (6, 8, $MENUENABLED, 'CHECKSOIL'c, CHECKSOIL)
	result = INSERTMENUQQ (6, 9, $MENUENABLED, 'NACTION'c, SHOWNA)
	result = INSERTMENUQQ (6, 10, $MENUENABLED, 'CANCEL'c, SHOWCANCEL)
								
	oldcolor=setbkcolorrgb(#00)
	oldcolor=setcolorrgb(#0000ff)
    CALL CLEARSCREEN( $GCLEARSCREEN )

	!CALL SETVIEWPORT( INT2(10), INT2(10), int2(940), int2(590) )
    thescreen.numxpixels=thescreen.numxpixels*0.98
	thescreen.numypixels=thescreen.numypixels*0.88

	CALL SETVIEWPORT( INT2(10), INT2(10), int2(thescreen.numxpixels), int2(thescreen.numypixels) )
    CALL SETLINEWIDTHQQ (3)
    status=RECTANGLE( $GBORDER,INT2(0),INT2(0), int2(thescreen.numxpixels-10), int2(thescreen.numypixels-10))
	!status=RECTANGLE( $GBORDER,INT2(1), INT2(1),int2(thescreen.numxpixels-11), int2(thescreen.numypixels-11) )
	
	CALL SETVIEWPORT( INT2(15), INT2(15), int2(thescreen.numxpixels-5), int2(thescreen.numypixels-5) )
	CALL CLEARSCREEN($GVIEWPORT)

	if(ymax-ymin>xmax-xmin) then	
	    t1=xmax
		ymax=ymax+(ymax-ymin)/20
        ymin=ymin-(ymax-ymin)/20
		xmax=(ymax-ymin)*(thescreen.numxpixels-31)/(thescreen.numypixels-31)+xmin
		xmin=xmin-(xmax-t1)/2
		xmax=xmax-(xmax-t1)/2
    else
	    t1=ymax
		xmax=xmax+(xmax-xmin)/20
        xmin=xmin-(xmax-xmin)/20		
		ymax=(xmax-xmin)*(thescreen.numypixels-31)/(thescreen.numxpixels-31)+ymin
		ymin=ymin-(ymax-t1)/2
		ymax=ymax-(ymax-t1)/2
	end if
	
	status=SETWINDOW( .true., xmin,ymin,xmax,ymax)
	oldcolor=setcolorrgb(#ff00)
	status=RECTANGLE_w( $GBORDER,xmin,ymin,xmax,ymax)
	
	!status=SETWINDOW( .false., xmin,ymin,xmax,ymax)
    
    iunit=90
	call checkgraph(iunit, event, ikeystate, ix, iy)

	event = MOUSE$LBUTTONDOWN
	event=  ior(event,MOUSE$LBUTTONUP)
    event = IOR (event, MOUSE$RBUTTONDOWN)
    event = IOR (event, MOUSE$LBUTTONDBLCLK)
    event = IOR (event, MOUSE$RBUTTONDBLCLK)
    event = IOR (event, MOUSE$MOVE)
	i4 = registermouseevent(90, event, checkgraph)
	
	

	result =  SETWSIZEQQ(90, winfo)

	do while(.true.)

	 !res = waitonmouseevent(event, keystate, ix, iy)
	 !if((MOUSE$KS_SHIFT .AND. res) == MOUSE$KS_SHIFT)   return
	 !i4 = registermouseevent(90, event, checkgraph)
	 !res = waitonmouseevent(MOUSE$move, keystate, ix, iy)

	 butter = GETCHARQQ()
	 tof1=(ichar(butter)==ichar('q')).or.(ichar(butter)==ichar('Q'))
	 tof2=(ichar(butter)==ichar('c')).or.(ichar(butter)==ichar('C'))

	 if(tof1) stop

	 if(tof2) then
	   !winfo%TYPE = QWIN$MIN
       !result =  SETWSIZEQQ(90, winfo)
	   close(90)
	   exit
	 end  if
	

	end do
	
	result = DELETEMENUQQ(6,0)	

	

   end subroutine

   subroutine checkgraph(iunit, ievent, ikeystate, ixpos, iypos)
      use ExcaDS
      use solverds
      USE DFLIB
	  implicit none
      
      integer(4)    oldcolor,iunit, ievent, ikeystate, ixpos, iypos
	  INTEGER(2)    status,result
      INTEGER(4)     res
	  character(256)  msg
	  character(16) CH,CH1
	  integer::i,j,k,n1,n2,N3,i4
      logical::tof,tof1,tof2
      real(8)::x1,y1,x2,y2,size,scale
	  real(8)::xmin,ymin,xmax,ymax
	  common /xy/ xmin,ymin,xmax,ymax
	  TYPE (windowconfig):: thescreen 
	  common /c2/ thescreen
	  TYPE (wxycoord)       wxy

	  scale=1

	  !if(i4>=0) then
         call getwindowcoord(ixpos-15, iypos-15,wxy)
		
		 if(ievent==MOUSE$MOVE) then
		    write(msg,*) 'Put C to close window and continue meshgen. Put Q to terminate program.'
            CALL SETMESSAGEQQ (msg, QWIN$MSG_INPUTPEND)
		    return
	     end if

		 if((MOUSE$KS_lbutton .AND. ikeystate) == MOUSE$KS_lbutton) scale=scale/2
	          	
		 if((MOUSE$KS_Rbutton .AND. ikeystate) == MOUSE$KS_Rbutton) scale=scale*2


		 xmin=wxy.wx-(wxy.wx-xmin)*scale
		 ymin=wxy.wy-(wxy.wy-ymin)*scale
		 xmax=wxy.wx+(xmax-wxy.wx)*scale
		 ymax=wxy.wy+(ymax-wxy.wy)*scale
         CALL SETVIEWPORT( INT2(15), INT2(15), int2(thescreen.numxpixels-5), int2(thescreen.numypixels-5) )
	     CALL CLEARSCREEN($GVIEWPORT)
		 oldcolor=setcolorrgb(#0000ff)
	     status=SETWINDOW( .TRUE., xmin,ymin,xmax,ymax)
	     oldcolor=setcolorrgb(#ff00)
	     status=RECTANGLE_w($GBORDER,xmin,ymin,xmax,ymax)
	
	 ! end if
  	
	  result = INITIALIZEFONTS()
      result = SETFONT('t''Arial''h15w8p')
	  SIZE=((XMIN-XMAX)**2+(YMAX-YMIN)**2)**0.5/1200
	  do i=1,enum
		
		n1=element(i).node(1)
		IF(ELEMENT(I).NNUM>1) THEN
            n2=element(i).node(2)
        ELSE
            N2=N1
        ENDIF    
		x1=node(n1).coord(1)
		y1=node(n1).coord(ndimension)
		x2=node(n2).coord(1)
		y2=node(n2).coord(ndimension)
		if(min(x1,x2)>xmax) cycle
		if(max(x1,x2)<xmin) cycle
		if(min(y1,y2)>ymax) cycle
		if(max(y1,y2)<ymin) cycle
		select case(showvalue)
			case(8)
				call setc(element(i).et)
			case default
				call setc(element(i).mat)
		end select

		call moveto_w(x1,y1,wxy)		
		status=lineto_w(x2,y2)
		
	  end do
	  
	  call  PLOTDATA(xmin,ymin,xmax,ymax,size)
  




   end subroutine

	subroutine showY()
		use EXCADS
		implicit none
		showvalue=2
	end subroutine

	subroutine shown()
		use EXCADS
		implicit none
		showvalue=1
	end subroutine
	subroutine showBW()
		use EXCADS
		implicit none
		showvalue=3
	end subroutine
	subroutine showBC()
		use EXCADS
		implicit none
		showvalue=4
	end subroutine

	subroutine SHOWELNUM()
		use EXCADS
		implicit none
		showvalue=5
	end subroutine
	subroutine SHOWELMAT()
		use EXCADS
		implicit none
		showvalue=6
	end subroutine	

	subroutine SHOWCANCEL()
		use EXCADS
		implicit none
		showvalue=-1
	end subroutine

	subroutine CHECKSOIL()
		use EXCADS
        USE SOLVERDS
		implicit none
        
        CGSTEP=CGSTEP+1
        CGSTEP=MOD(CGSTEP,NSTEP)
        IF(CGSTEP==0) CGSTEP=nstep
        call Excavation(cgstep)
        
		showvalue=10
	end subroutine
	subroutine SHOWET()
		use EXCADS
		implicit none
		showvalue=8
	end subroutine
	subroutine SHOWNA()
		use EXCADS
		implicit none
		showvalue=9
    end subroutine				
    
	subroutine PLOTDATA(xmin,ymin,xmax,ymax,size)
		USE solverds
        use EXCADS
		use dflib
		implicit none
		integer::n1,n2,i,J,flag
		integer(4)    oldcolor,status
		REAL(8)::X1,Y1,X2,Y2,size,xmin,ymin,xmax,ymax,T1
		character(32)::ch,ch1
		TYPE (wxycoord)       wxy
        type(soillayer_tydef)::SOIL1

		ch=''
		select case(showvalue)
		   case(1,2,3)
				do i=1,nnum
					x1=node(i).COORD(1)
					y1=node(i).COORD(NDIMENSION)
					if(x1>xmax) cycle
					if(x1<xmin) cycle
					if(y1>ymax) cycle
					if(y1<ymin) cycle
					call moveto_w(x1,y1,wxy)
					oldcolor=setcolorrgb(#C0C0C0 )
					status=ellipse_W($GFILLINTERIOR,x1-size,Y1-size,X1+size,Y1+size)
					if(showvalue==1) then
						write(ch,'(I6)') i
					else
						if(showvalue==2) then
                            
							write(ch,'(F10.3)') node(i).COORD(1)
                            write(ch1,'(F10.3)') node(i).COORD(2)
                            ch=trim(adjustL(ch))//","//trim(adjustL(ch1))
						else
							if(node(i).dof(1)/=0) then
								write(ch,'(I6)') node(i).dof(1)
							else
								write(ch,'(I6)') node(i).dof(2)
							end if
						end if
						
					end if
					oldcolor=setcolorrgb(#00FF00)
					call outgtext(trim(adjustL(CH)))

				end do
				
		   case(4)
			   do i=1,bd_num
					n1=BC_DISP(I).NODE
					x1=node(n1).COORD(1)
					y1=node(n1).COORD(NDIMENSION)
					if(x1>xmax) cycle
					if(x1<xmin) cycle
					if(y1>ymax) cycle
					if(y1<ymin) cycle
					ch=''

					write(ch,'(i1)')  BC_DISP(I).DOF

					call moveto_w(x1,y1,wxy)
					oldcolor=setcolorrgb(#00FF00)
					call outgtext(trim(adjustL(CH)))
				end do
			CASE(5,6,8)
				do i=1,enum
					x1=node(element(i).NODE(1)).COORD(1)
					y1=node(element(i).NODE(1)).COORD(NDIMENSION)
                    IF(ELEMENT(I).NNUM>1) THEN
					    x2=node(element(i).NODE(2)).COORD(1)
					    y2=node(element(i).NODE(2)).COORD(NDIMENSION)
                    ELSE
                        X2=X1
                        Y2=Y1
                    ENDIF
					if(min(x1,x2)>xmax) cycle
					if(max(x1,x2)<xmin) cycle
					if(min(y1,y2)>ymax) cycle
					if(max(y1,y2)<ymin) cycle
					call moveto_w(x1,y1,wxy)
					oldcolor=setcolorrgb(#C0C0C0 )
					status=ellipse_W($GFILLINTERIOR,x1-size,Y1-size,X1+size,Y1+size)
					call moveto_w(x2,y2,wxy)
					oldcolor=setcolorrgb(#C0C0C0 )
					status=ellipse_W($GFILLINTERIOR,x2-size,Y2-size,X2+size,Y2+size)
					call moveto_w((x1+x2)/2,(y1+y2)/2,wxy)
					SELECT CASE(SHOWVALUE)
						CASE(5)
							write(ch,'(I6)') I
						CASE(6)
							write(ch,'(I6)') ELEMENT(I).MAT
						CASE(8)
							write(ch,'(I6)') ELEMENT(I).ET
					END SELECT
					oldcolor=setcolorrgb(#00FF00)
					call outgtext(trim(adjustL(CH)))
                end do
			!CASE(7)
			!	do i=1,struts
			!		n1=nloc(strutdata(i).num,strutdata(i).num1)
			!		x1=node(n1).y
			!		y1=node(n1).z
			!		if(x1>xmax) cycle
			!		if(x1<xmin) cycle
			!		if(y1>ymax) cycle
			!		if(y1<ymin) cycle
			!		write(ch,'(i6)') i
			!		CALL moveto_w(x1,y1,wxy)
			!		oldcolor=setcolorrgb(#FFFF00)
			!		status=RECTANGLE_w($GFILLINTERIOR,x1-size*2,Y1-size*2,X1+size*2,Y1+size*2)
			!		oldcolor=setcolorrgb(#00FF00)
			!		call outgtext(trim(adjustL(CH)))			 
			!	end do
			!CASE(9)
			!	do i=1,na
			!		n1=nloc(naction(i).num,naction(i).num1)
			!		x1=node(n1).y
			!		y1=node(n1).z
			!		if(x1>xmax) cycle
			!		if(x1<xmin) cycle
			!		if(y1>ymax) cycle
			!		if(y1<ymin) cycle
			!		write(ch,'(f5.1)') naction(i).act(1)
			!		CALL moveto_w(x1,y1,wxy)
			!		oldcolor=setcolorrgb($LOGREEN)
			!		status=ellipse_w($GFILLINTERIOR,x1-size*2,Y1-size*2,X1+size*2,Y1+size*2)
			!		oldcolor=setcolorrgb(#00FF00)
			!		call outgtext(trim(adjustL(CH)))
			!	end do
            case(10) !soilprofile
                T1=abs(MAXVAL(KPOINT(NDIMENSION,:))-MINVAL(KPOINT(NDIMENSION,:)))/3
                do i=1,nsoilprofile
                    do j=1,soilprofile(i).nasoil+soilprofile(i).npsoil
                        if(j<=soilprofile(i).nasoil) then
                            SOIL1=soilprofile(i).asoil(j)
                            flag=1
                        else
                            SOIL1=soilprofile(i).psoil(j-soilprofile(i).nasoil)
                            flag=-1
                        endif
                        
                        if(sf(SOIL1.sf).factor(CGSTEP)==0) cycle
                        
                        call setc(SOIL1.mat)
                        x1=kpoint(1,SOIL1.z(1))
                        y1= kpoint(ndimension,SOIL1.z(1))
                        y2= kpoint(ndimension,SOIL1.z(2))
                        x2=x1-T1*soilprofile(i).aside*flag
                        status=RECTANGLE_w($GFILLINTERIOR,x1,y1,x2,y2)                
                            
                    enddo
                    
                enddo    
            CASE DEFAULT
                PRINT *, "TO BE IMPROVED.SUB=PLOTDATA"
				
            end select
            
10 Format("(",F7.3,",",F7.3,")")            
    end subroutine
    
    

	subroutine setc(n1)
		use dflib
		implicit none
		integer::n1
		integer(4)::oldcolor
		select case(mod(n1,15))
		   case(1)
              oldcolor=setcolorrgb($HIRED)
		   case(2)
              oldcolor=setcolorrgb($HIGREEN)
		   case(3)
              oldcolor=setcolorrgb($HIBLUE)
		   case(4)
              oldcolor=setcolorrgb($HIYELLOW)
		   case(5)
              oldcolor=setcolorrgb($HICYAN)
		   case(6)
              oldcolor=setcolorrgb($HIMAGENTA)
		   case(7)
              oldcolor=setcolorrgb($HIWHITE)
		   case(8)
			  oldcolor=setcolorrgb($Logray)
		   case(9)
              oldcolor=setcolorrgb($LORED)
		   case(10)
              oldcolor=setcolorrgb($LOGREEN)
		   case(11)
              oldcolor=setcolorrgb($LOBLUE)
		   case(12)
              oldcolor=setcolorrgb($LOBROWN)
		   case(13)
              oldcolor=setcolorrgb($LOCYAN)
		   case(14)
              oldcolor=setcolorrgb($LOMAGENTA)
		   case(0)
              oldcolor=setcolorrgb($Lowhite)			
!		   case default
!              oldcolor=setcolorrgb(#fff00)
		 end select
	end subroutine


    
