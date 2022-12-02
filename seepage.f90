SUBROUTINE SPG_Q_UPDATE(STEPDIS,bload,HHEAD,INIHEAD,DT,nbload,ienum,iiter,istep,iscon)
    USE solverds
    IMPLICIT NONE
    logical,intent(in)::iscon
	integer,intent(in)::iiter,ienum,istep,nbload
	real(kind=DPN),intent(in)::STEPDIS(NDOF),HHEAD(nbload),INIHEAD(nbload),DT
	real(kind=DPN),intent(out)::bload(nbload)
    INTEGER::J,N1,ND1
    REAL(DPN)::KT1(NDIMENSION,NDIMENSION),HJ,hj_ini,slope1,sita1,R1,A1,A2,A3,t1
    
    
    SELECT CASE(ELEMENT(IENUM).ET)
        
    CASE(PIPE2,WELLBORE,WELLBORE_SPGFACE)
        if(solver_control.wellmethod==0) then !sample at node
            CALL WELLBORE_Q_K_UPDATE3(STEPDIS,IENUM,ISTEP,IITER)
        elseif(solver_control.wellmethod==1) then !sample at centroid
            CALL WELLBORE_Q_K_UPDATE(STEPDIS,IENUM,ISTEP,IITER)
        else
            CALL WELLBORE_Q_K_UPDATE_SPMETHOD(STEPDIS,IENUM,ISTEP,IITER)
        endif
        BLOAD=ELEMENT(IENUM).FLUX 
    CASE(SPHFLOW,SEMI_SPHFLOW)
        SELECT CASE(solver_control.wellmethod) 
        CASE(0,2)
            CALL SPHFLOW_Q_K_UPDATE3(STEPDIS,IENUM,ISTEP,IITER)
        CASE(1)
            CALL SPHFLOW_Q_K_UPDATE(STEPDIS,IENUM,ISTEP,IITER)
        CASE(3)
            CALL SPHFLOW_Q_K_SPMETHOD(STEPDIS,IENUM,ISTEP,IITER)
        CASE DEFAULT
            CALL SPHFLOW_Q_K_UPDATE_ANALYTICAL(STEPDIS,IENUM,ISTEP,IITER)
        END SELECT
        
        BLOAD=ELEMENT(IENUM).FLUX
    CASE(ZT4_SPG,ZT6_SPG)
        CALL ZTE_Q_K_UPDATE(STEPDIS,IENUM,ISTEP,IITER)
        BLOAD=ELEMENT(IENUM).FLUX
    CASE(POREFLOW) !
        CALL POREFLOW_Q_K_UPDATE2(STEPDIS,IENUM,ISTEP,IITER)
        BLOAD=ELEMENT(IENUM).FLUX 

    CASE DEFAULT

	    ND1=ELEMENT(IENUM).ND
        N1=element(ienum).NGP
	    do j=1,N1
                   
		    hj=dot_product(ecp(element(ienum).et).Lshape(1:nbload,J),HHEAD(1:nbload))
        
		    R1=1.D0
		    if(element(ienum).ec==cax_spg) R1=ABS(element(ienum).xygp(1,j))
					
		    !gradient
		    element(ienum).igrad(1:nd1,j)=matmul(element(ienum).B(:,:,j),HHEAD(1:NBLOAD))
		    !velocity
		    CALL SPG_KT_UPDATE(KT1,HHEAD,HJ,nbload,IENUM,J,ISTEP,IITER)
            
		    element(ienum).velocity(1:nd1,j)=-matmul(KT1,element(ienum).igrad(1:nd1,j))
		    !flux
		    if(.not.stepinfo(istep).issteady) then
			    hj_ini=dot_product(ecp(element(ienum).et).Lshape(1:nbload,j),INIHEAD)
			    call slope_SWWC_spg(ienum,hj,element(ienum).xygp(ndimension,j),slope1,sita1,element(ienum).sita_ini(j),hj_ini,ISTEP)
                        
                element(ienum).mw(j)=slope1
                element(ienum).sita_fin(j)=sita1
						
			    Qstored=Qstored+sita1*ecp(element(ienum).et).weight(j)*element(ienum).detjac(j)
						
						
			    if(j==1) element(ienum).cmm=0.d0  !clear 
                if(slope1/=0.0d0) then
				    element(ienum).cmm=element(ienum).cmm+csproduct(ecp(element(ienum).et).Lshape(:,j),ecp(element(ienum).et).Lshape(:,j))* &
								    (R1*slope1*ecp(element(ienum).et).weight(j)*element(ienum).detjac(j)*MATERIAL(ELEMENT(ienum).MAT).GET(13,ISTEP)/DT)
                end if
                if((j==n1).and.any(element(ienum).cmm/=0.d0)) then
				    bload(1:NBLOAD)=bload(1:NBLOAD)+matmul(element(ienum).cmm,(HHEAD-INIHEAD))
			    end if
		    end if
										
		    bload(1:NBLOAD)=bload(1:NBLOAD)+ &
		    matmul(-element(ienum).velocity(1:nd1,j),element(ienum).b(:,:,j))* &
		    element(ienum).detjac(j)*ecp(element(ienum).et).weight(j)*R1
       
            if(solver_control.bfgm==continuum.or.solver_control.bfgm==consistent.or.solver_control.bfgm==inistress) then
                isref_spg=1
                if(j==1) element(ienum).km=0.0D0
			    element(ienum).km=element(ienum).km+ &
				    matmul(MATMUL(TRANSPOSE(element(ienum).b(:,:,j)),KT1),element(ienum).b(:,:,j))* &
				    (element(ienum).detjac(j)*ecp(element(ienum).et).weight(j)*R1)
       !         if(j==n1) then
			    !	bload(1:NBLOAD)=bload(1:NBLOAD)+matmul(element(ienum).km,HHEAD(1:NBLOAD))
			    !end if							
		    end if
					
	    end do
        
        element(ienum).flux(1:NBLOAD)=bload(1:NBLOAD)
    
    END SELECT
    
	    
    
END SUBROUTINE
    
!subroutine q_updated_control_volume_method()
!    use solverds
!    
!    if(.not.isIniSEdge) CALL Model_MESHTOPO_INI()
!    
!    
!    
!end subroutine
    
    
subroutine lamda_spg(ienum,hj,z,lamda,ISTEP)
	
	use solverds
	implicit none
	integer,intent(in)::ienum,ISTEP
	real(8),intent(in)::hj,z
	real(8),intent(out)::lamda
	real(8)::t1,epsilon1,epsilon2,krsml
	real(8)::alpha1,fn1,fm1,seff1,scale1=1.0
	real(8)::h1(50)
    
	scale1=1.0d0
	t1=hj-z
	krsml=1.e-3
    
   

!	epsilon1=max(eps1,0.1)
!	epsilon2=max(eps2,0.1)
	select case(material(element(ienum).mat).type)
		case(step_spg)
			if(t1>=0.d0) then
				lamda=1.d0
			else
				lamda=krsml
			end if
			
		case(vg_spg)

			if (t1 .ge. 0.d0) then
				lamda = 1.d0
			else
			   alpha1 = material(element(ienum).mat).GET(7,ISTEP)
			   fn1 = material(element(ienum).mat).GET(8,ISTEP)
				fm1 = 1.d0 - 1.d0 / fn1
				seff1 = (1.d0 + (-alpha1*t1)**fn1) ** (- fm1)
				lamda = (1.d0 - (1.d0 - seff1 ** (1.d0 / fm1)) ** fm1) ** 2 * dsqrt (seff1)
            endif
        case(lr_spg)
            if(t1>=0.d0) then
               lamda = 1.d0     
            else
               alpha1 = material(element(ienum).mat).GET(7,ISTEP)
			   fn1 = material(element(ienum).mat).GET(8,ISTEP)
               fm1= material(element(ienum).mat).GET(9,ISTEP)
               lamda=(dlog(dexp(1.d0)+(-t1/alpha1)**fn1))**(-fm1)
            end if
        case(exp_spg)
            if(t1>=0.d0) then
               lamda = 1.d0     
            else
               alpha1 = material(element(ienum).mat).GET(7,ISTEP)            
               lamda=dexp(alpha1*t1)
            end if
!		case(slope_spg)
			
		case default
!			if(mod(iiter,1000)==0) scale1=scale1*2 
			epsilon1=element(ienum).PROPERTY(3)
			epsilon2=element(ienum).PROPERTY(2)
			
			if(t1>=epsilon2) then
				lamda=1.0
			else
				if(t1<=-epsilon1) then
					lamda=krsml
				else
					lamda=(1-krsml)/(epsilon1+epsilon2)*(t1-epsilon2)+1
				end if
			end if
			
			!if(lamda<1e-3) lamda=1.e-3
			
!			if(iiter>30) then
!		         if(lamda>10*element(ienum).lamda(igp)) then
!                       lamda=10*element(ienum).lamda(igp)
!		         else
!		              if(lamda<0.1*element(ienum).lamda(igp)) lamda=0.1*element(ienum).lamda(igp)
!		          end if
!			end if
!			
!			element(ienum).lamda(igp)=lamda
				
	end select

	
end subroutine

subroutine slope_SWWC_spg(ienum,hj,z,slope,sita,sita_ini,hj_ini,ISTEP)
    use solverds
    implicit none
    integer,intent(in)::ienum,ISTEP
    Real(kind=DPN),intent(in)::hj,z,sita_ini,hj_ini
    Real(kind=DPN),intent(out)::slope,sita
    real(kind=DPN)::t1,alpha1,fn1,fm1,seff1,sita_s,sita_r,rw1
    
    t1=hj-z
    sita_s = material(element(ienum).mat).GET(11,ISTEP)!饱和体积含水量
    
    if(t1>=0.d0) then
        slope=material(element(ienum).mat).GET(10,ISTEP)
		sita=sita_s
        return
    end if
    
    !sita
    alpha1 = material(element(ienum).mat).GET(7,ISTEP)
    fn1 = material(element(ienum).mat).GET(8,ISTEP)
    sita_r = material(element(ienum).mat).GET(12,ISTEP) !!残余体积含水量(默认为0)
    rw1=material(element(ienum).mat).GET(13,ISTEP) !水的重度    
    
   select case(material(element(ienum).mat).type)
					
        case(vg_spg)

            fm1 = 1.d0 - 1.d0 / fn1
            seff1 = (1.d0 + (-alpha1*t1)**fn1) ** (-fm1)
            sita=sita_r+(sita_s-sita_r)*seff1
			
			!slope:mw2
			if(abs(hj-hj_ini)<1e-7.or.solver_control.mur==0) then
				slope=fm1*fn1*(-t1*alpha1)**fn1*(sita_s-sita_r)/(-t1*rw1*(1.d0+(-t1*alpha1)**fn1)**(fm1+1))
			else
				slope=(sita-sita_ini)/((hj-hj_ini)*rw1)
			end if
			
        case(lr_spg)
            fm1= material(element(ienum).mat).GET(9,ISTEP)
            sita=sita_r+(sita_s-sita_r)*(dlog(dexp(1.d0)+(-t1/alpha1)**fn1))**(-fm1)
            !slope=mw2
			if(abs(hj-hj_ini)<1e-7.or.solver_control.mur==0) then
				slope=dexp(1.d0)+(-t1/alpha1)**fn1
				slope=fm1*fn1*(-t1/alpha1)**fn1*(sita_s-sita_r)/(-t1*rw1*slope*(dlog(slope))**(fm1+1))
			else
				slope=(sita-sita_ini)/((hj-hj_ini)*rw1)
			end if
           
        case(exp_spg)
            sita=sita_r+(sita_s-sita_r)*dexp(alpha1*t1)
			if(abs(hj-hj_ini)<1e-7.or.solver_control.mur==0) then
				slope=(sita_s-sita_r)*dexp(alpha1*t1)*alpha1/rw1
			else
				slope=(sita-sita_ini)/((hj-hj_ini)*rw1)
			end if
    
        case default            
				sita=sita_ini
				slope=0
		!stop "No such an SWCC function."
    end select 
   
end subroutine

SUBROUTINE SPG_KT_UPDATE(KT,HHEAD,HJ,NHH,IENUM,IGP,ISTEP,IITER)

	USE SOLVERDS
	IMPLICIT NONE
	INTEGER,INTENT(IN)::NHH,IENUM,IGP,ISTEP,IITER
    REAL(DPN),INTENT(IN)::HHEAD(NHH),HJ
	REAL(DPN),INTENT(OUT)::KT(NDIMENSION,NDIMENSION)
	REAL(DPN)::ROT1(2,2),COS1,SIN1,T1,FAC1,SITA1,LAMDA,C1,PHI1,SFR1(8),SIGMA1(6),TS1,PI1,MU1
    INTEGER::I,IENUMC1,IE1,TSFAIL1=0	
    
    
	
	IF(MINVAL(HHEAD-NODE(ELEMENT(IENUM).NODE(1:NHH)).COORD(NDIMENSION))>-1.D-4) THEN
        LAMDA=1.0D0
    ELSE
	    call lamda_spg(ienum,hj,element(ienum).xygp(ndimension,IGP),lamda,ISTEP)                   						
    ENDIF                
    if(iiter>1) lamda=(element(ienum).kr(IGP)+lamda)/2.0D0
                    
    element(ienum).kr(IGP)=lamda 
    
    KT=ELEMENT(IENUM).D*LAMDA*element(ienum).fqw
    
        
	IF(ESET(ELEMENT(IENUM).SET).COUPLESET<0.OR.ESET(ELEMENT(IENUM).SET).COUPLESET==ELEMENT(IENUM).SET) RETURN
    
    PI1=3.14159265358979
    IE1=IENUM-(ESET(ELEMENT(IENUM).SET).ENUMS-1)
    IENUMC1=ESET(ESET(ELEMENT(IENUM).SET).COUPLESET).ENUMS-1+IE1
    
    !!计算弹性（第一次迭代计算时为弹性计算）ko状态下的sfr，只计算一次.
    !IF(.not.allocated(element(ienumc1).sfrko)) then
    !    allocate(element(ienumc1).sfrko(element(ienumc1).ngp))
    !    mu1=material(ELEMENT(IENUMC1).MAT).GET(2,ISTEP)
    !    C1=material(ELEMENT(IENUMC1).MAT).GET(3,ISTEP)
	   ! Phi1=material(ELEMENT(IENUMC1).MAT).GET(4,ISTEP)
    !    TS1=material(ELEMENT(IENUMC1).MAT).GET(5,ISTEP)
    !    IF(material(ELEMENT(IENUMC1).MAT).TYPE==MC)THEN
    !        TS1=C1/MAX(tan(phi1/180*PI1),1E-7) 
    !    ENDIF
    !    !!假定ko应力，ko=v/(1-v),sxx=k0*syy,szz=sxx,txy=0
    !    DO I=1,element(ienumc1).ngp
    !        SIGMA1=ELEMENT(IENUMC1).STRESS(:,I)+ELEMENT(IENUMC1).DSTRESS(:,I) !!迭代时应力还没更新
    !        sigma1(1)=mu1/(1-mu1)*sigma1(2);sigma1(3)=sigma1(1);sigma1(4:6)=0
	   !     call stress_in_failure_surface(sfr1,SIGMA1,2,C1,Phi1,solver_control.slidedirection,ELEMENT(IENUMC1).XYGP(1:NDIMENSION,I),TS1)
    !        element(ienumc1).sfrko(I)=SFR1(1)
    !    ENDDO
    !ENDIF    
        
    SFR1=ELEMENT(IENUMC1).SFR(:,IGP)

    !应力场与渗流场分开计算时，maxsfr被重新reset为负大数。
    IF(ISTEP>1.AND.MAXSFR<-1.D10) THEN
        MAXSFR=MAXSFR_LAST
        IF(MAXSFR<-1.D10) STOP "ERROR IN MAXSFR VALUE."
    ENDIF
    
    TSFAIL1=0
    IF(ABS(SFR1(1)+SFR1(8)+1.0D0)<1E-6) TSFAIL1=1 !令受拉破坏的SFR1=maxsfr    
    IF(TSFAIL1==1) SFR1(1)=MAXSFR
    !SFR1(1)=ABS(SFR1(1)-element(ienumc1).sfrko(IGP))
    
    
    FAC1=MAX((ABS(SFR1(1)))**solver_control.slope_kscale,1.D-9)
    SITA1=-SFR1(2)/180.*PI()
    
    !element(ienum).kr(IGP)=element(ienum).kr(IGP)*fac1
    
	!KT=0.D0
	!T1=MATERIAL(ELEMENT(IENUM).MAT).GET(1,ISTEP)
	!KT(2,2)=MATERIAL(ELEMENT(IENUM).MAT).GET(2,ISTEP)*FAC1
    IF(solver_control.slope_kbase>=0) THEN
        lamda=(solver_control.slope_kbase+FAC1)    
    ELSE
        !solver_control.slope_kbase=0.1
	    lamda=(ABS(solver_control.slope_kbase)*MAXSFR+FAC1)
    ENDIF
    IF(ABS(SFR1(1))/MAXSFR>=solver_control.slope_sfrpeak) THEN
        lamda=lamda*exp(2*abs(sfr1(1)))
    ENDIF
    !element(ienum).kr(IGP)=lamda
    ELEMENT(IENUMC1).SFR(7,IGP)=LAMDA
    if(solver_control.slope_kratio>0) then
        kt(2,2)=kt(1,1)/solver_control.slope_kratio
    endif
    KT=KT*lamda

    IF(ABS(KT(1,1)-KT(2,2))>1E-7) THEN        
	    COS1=COS(SITA1);SIN1=SIN(SITA1)
	    ROT1(1,1)=COS1;ROT1(2,2)=COS1;ROT1(1,2)=SIN1;ROT1(2,1)=-SIN1
	    KT=MATMUL(MATMUL(ROT1,KT),TRANSPOSE(ROT1))
    ENDIF
    
    !if(any(isnan(kt))) then
    !    pause
    !endif	

END SUBROUTINE

	!update Signorini boundary condition
SUBROUTINE SPG_Signorini_BC_UPDATE(DISQ,STEPDIS,ISTEP)
	USE SOLVERDS
	IMPLICIT NONE
	INTEGER,INTENT(IN)::ISTEP
	REAL(8),INTENT(IN)::DISQ(NDOF),STEPDIS(NDOF)
	INTEGER::I,N1
	REAL(8)::T1,T2
	

	isref_spg=0
	do i=1,numNseep
        if(Nseep(i).node<1) cycle
		if(sf(NSeep(i).sf).factor(istep)==-999.D0) cycle
		
		if(Nseep(i).isdual>0) then
			if(bc_disp(Nseep(i).isdual).isdead==0) cycle            
		end if
		
		n1=node(Nseep(i).node).dof(Nseep(i).dof)
		t1=DISQ(n1)
		
		if(stepinfo(istep).issteady) then
			t2=stepdis(n1)-node(Nseep(i).node).coord(ndimension)
		else
			t2=stepdis(n1)-node(Nseep(i).node).coord(ndimension) !for a transient problem, tdisp is passed to stepdisp in the form of an initial value. 
		end if		
		
				
		if(Nseep(i).isdead==0) then
			if(t1>1E-7) then
				Nseep(i).isdead=1
				isref_spg=1
			end if
		else
			if(t2>1E-3) then
				Nseep(i).isdead=0
				isref_spg=1
			end if
		end if
!		write(99,20) iiter,i,Nseep(i).isdead,n1,t1,t2
	end do
	
	
ENDSUBROUTINE

SUBROUTINE ZT_SPG_INI2(IELT)
    USE solverds
    USE SolverMath
    use MESHADJ,only:Setup_Solver_MESHTOPO
    IMPLICIT NONE
    INTEGER,INTENT(IN)::IELT
    !INTEGER,INTENT(INOUT)::ISPF(:,:)
    integer::i,j,K
	integer::n1,n2,n3,n4,AELT1(2),IN1(10)
	real(kind=DPN)::t1=0,vcos=0,vsin=0,rpi,coord1(3,4)=0,cent1(3,2)=0,c1(3,3)=0,VEC1(3),t2
    
    
    if(.not.isIniSEdge) CALL Setup_Solver_MESHTOPO()
    !找出zt单元的由face1指向face2的向量VEC1
    IF(ELEMENT(IELT).ET==ZT4_SPG.OR.ELEMENT(IELT).ET==ZT4_SPG2) THEN
        AELT1=ELEMENT(IELT).ADJELT([1,3])
    ELSE
        AELT1=ELEMENT(IELT).ADJELT([1,2])
    ENDIF
    where(AELT1<1)   AELT1=IELT
    
    !把防渗墙临空面上的节点设为出溢边界
    if(.not.allocated(ispf)) then
        allocate(ispf(nnum,2))
        ispf=0
    endif    
    N2=ELEMENT(IELT).NNUM/2
    DO I=1,2
        N1=(I-1)*N2
        IF(AELT1(I)==IELT) THEN
            ISPF(ELEMENT(IELT).NODE(N1+1:N1+N2),1)=ELEMENT(IELT).NODE(N1+1:N1+N2)
            ISPF(ELEMENT(IELT).NODE(N1+1:N1+N2),2)=ELEMENT(IELT).SF
        ENDIF
    ENDDO

    DO I=1,2
        DO J=1,NDIMENSION
            CENT1(J,I)=SUM(NODE(ELEMENT(AELT1(I)).NODE).COORD(J))/ELEMENT(AELT1(I)).NNUM            
        ENDDO
    ENDDO
    VEC1=(CENT1(:,2)-CENT1(:,1))
    t1=NORM2(VEC1)
    
    IF(ABS(t1)>1.E-6) THEN
    !    PRINT *, 'ZT单元I的指向向量无法确定.I=',ielt
    !ELSE
        vec1=vec1/t1
    ENDIF
    
    T2=1.0D0
    T1=NORM2(NODE(ELEMENT(IELT).NODE(1)).COORD-NODE(ELEMENT(IELT).NODE(2)).COORD)
    C1=0.d0
    IF(ELEMENT(IELT).ET==ZT4_SPG.OR.ELEMENT(IELT).ET==ZT6_SPG) THEN
        allocate(element(IELT).km(ELEMENT(IELT).NNUM,ELEMENT(IELT).NNUM))
        ELEMENT(IELT).KM=0.D0
        DO I=1,ELEMENT(IELT).NNUM
            ELEMENT(IELT).KM(I,I)=1.D0
        ENDDO
    
        ELEMENT(IELT).KM(1,4)=-1.D0
        ELEMENT(IELT).KM(4,1)=-1.D0
    ENDIF
    
	IF(ELEMENT(IELT).ET==ZT4_SPG.OR.ELEMENT(IELT).ET==ZT4_SPG2) THEN

		!LOCAL TO GLOBAL 
		VCOS=(NODE(ELEMENT(IELT).NODE(2)).COORD(1)-NODE(ELEMENT(IELT).NODE(1)).COORD(1))/T1
		VSIN=(NODE(ELEMENT(IELT).NODE(2)).COORD(2)-NODE(ELEMENT(IELT).NODE(1)).COORD(2))/T1
		C1(1,1)=VCOS
		C1(1,2)=VSIN
		C1(2,2)=VCOS
		C1(2,1)=-VSIN
        C1(3,3)=1.D0
        
        ELEMENT(IELT).SIGN=INT(sign(1.0,dot_product(vec1,c1(2,:))))
        
        !使薄元的局部坐标正向(厚度方向)与两相邻实体单元的方向一致(相邻单元1(与薄元单元面1相邻的单元)指向相邻单元2.)。
        IF(ELEMENT(IELT).SIGN<0) THEN
            IN1(1:4)=ELEMENT(IELT).NODE([2,1,4,3])
            ELEMENT(IELT).NODE=IN1(1:4)
            IN1(1:4)=-ELEMENT(IELT).EDGE([1,4,3,2])
            ELEMENT(IELT).EDGE=IN1(1:4)
            IN1(1:4)=-ELEMENT(IELT).ADJELT([1,4,3,2])
            ELEMENT(IELT).ADJELT=IN1(1:4)
            
            ELEMENT(IELT).SIGN=1
		    C1(1:2,1:2)=-C1(1:2,1:2)            
        ENDIF        
        
        ELEMENT(IELT).PROPERTY(1)=T1
        IF(ELEMENT(IELT).ET==ZT4_SPG) THEN
            ELEMENT(IELT).KM(2,3)=-1.D0
            ELEMENT(IELT).KM(3,2)=-1.D0        
        ENDIF    
	ELSE
	
	    !LOCAL TO GLOBAL
                    
        do j=1,3
            C1(:,j)=NODE(ELEMENT(IELT).NODE(j)).COORD;
        enddo
        C1(3,:)=NORMAL_TRIFACE(C1)
        ELEMENT(IELT).PROPERTY(1)=0.5*norm2(C1(3,:))
        C1(3,:)=C1(3,:)/ELEMENT(IELT).PROPERTY(1)/2.0
        C1(1,:)=(NODE(ELEMENT(IELT).NODE(2)).COORD-NODE(ELEMENT(IELT).NODE(1)).COORD)/T1
        C1(2,:)=0.d0
        C1(2,:)=NORMAL_TRIFACE(RESHAPE([C1(2,:),C1(1,:),C1(3,:)],([3,3])))
        C1(2,:)=C1(2,:)/norm2(C1(2,:))
        
        ELEMENT(IELT).SIGN=INT(sign(1.0,dot_product(vec1,c1(3,:))))
        !使薄元的局部坐标正向(厚度方向)与两相邻实体单元的方向一致(相邻单元1(与薄元单元面1相邻的单元)指向相邻单元2.)。
        IF(ELEMENT(IELT).SIGN<0) THEN
           
            IN1(1:6)=ELEMENT(IELT).NODE([1,3,2,4,6,5])
            ELEMENT(IELT).NODE=IN1(1:6)
            IN1(1:9)=ELEMENT(IELT).EDGE([3,2,1,6,5,4,7,8,9])
            ELEMENT(IELT).EDGE=[-IN1(1:6),IN1(7:9)]
            IN1(1:5)=ELEMENT(IELT).FACE([1,2,5,4,3])
            ELEMENT(IELT).FACE=-IN1(1:5)
            IN1(1:5)=ELEMENT(IELT).ADJELT([1,2,5,4,3])
            ELEMENT(IELT).ADJELT=IN1(1:5)
            
            ELEMENT(IELT).SIGN=1
            !UPDATE SYS
            C1(3,:)=-C1(3,:)
            C1(1,:)=-C1(1,:)
            
        ENDIF        
        IF(ELEMENT(IELT).ET==ZT6_SPG) THEN
            ELEMENT(IELT).KM(2,5)=-1.D0
            ELEMENT(IELT).KM(5,2)=-1.D0        
            ELEMENT(IELT).KM(3,6)=-1.D0
            ELEMENT(IELT).KM(6,3)=-1.D0          
        ENDIF
    ENDIF
    ELEMENT(IELT).G2L=C1
    IF(ABS(MATERIAL(ELEMENT(IELT).MAT).PROPERTY(14))<1E-10) MATERIAL(ELEMENT(IELT).MAT).PROPERTY(14)=1.D0
    
    !generate ghost nodes
    N2=NDIMENSION
	IF(.NOT.ALLOCATED(GNODE)) ALLOCATE(GNODE(3,1000))
	IF(NGNODE+2*N2>SIZE(GNODE,DIM=2)) CALL ENLARGE_GNODE(1000)
    COORD1=0.D0;COORD1(N2,1)=MATERIAL(ELEMENT(IELT).MAT).PROPERTY(14)

    
	COORD1(1:N2,1)=MATMUL(TRANSPOSE(C1(1:N2,1:N2)),COORD1(1:N2,1))
   
    DO I=1,N2    
	    !NGNODE=NGNODE+1;
        GNODE(:,NGNODE+I)=NODE(ELEMENT(IELT).NODE(I)).COORD
        IF(ELEMENT(IELT).ET==ZT6_SPG.OR.ELEMENT(IELT).ET==ZT6_SPG2) THEN
            N3=I
        ELSE
            N3=MOD(I,N2)+1
        ENDIF
        GNODE(:,NGNODE+N2+N3)=COORD1(:,1)+NODE(ELEMENT(IELT).NODE(I)).COORD
    ENDDO

	ALLOCATE(ELEMENT(IELT).NODE2(2*N2))
	ELEMENT(IELT).NODE2=[NGNODE+1:NGNODE+2*N2]
    
    NGNODE=NGNODE+2*N2    
    
    IF(ELEMENT(IELT).ET==ZT4_SPG.OR.ELEMENT(IELT).ET==ZT6_SPG) THEN
        T1=MATERIAL(ELEMENT(IELT).MAT).PROPERTY(NDIMENSION)*ELEMENT(IELT).PROPERTY(1)/MATERIAL(ELEMENT(IELT).MAT).PROPERTY(14)/ELEMENT(IELT).NNUM*2 !KA/L
     
        ELEMENT(IELT).KM=ELEMENT(IELT).KM*T1
    
        IF(.NOT.ALLOCATED(element(IELT).igrad)) then
            ALLOCATE(element(IELT).igrad(ELEMENT(IELT).ND,ELEMENT(IELT).NGP+ELEMENT(IELT).NNUM))
            ALLOCATE(element(IELT).VELOCITY(ELEMENT(IELT).ND,ELEMENT(IELT).NGP+ELEMENT(IELT).NNUM))
            element(IELT).igrad=0.d0;element(IELT).VELOCITY=0.d0        
        endif
    
        DO I=1,ELEMENT(IELT).NNUM
            IF(.NOT.ALLOCATED(NODE(ELEMENT(IELT).NODE(I)).IGRAD)) THEN
                ALLOCATE(NODE(ELEMENT(IELT).NODE(I)).IGRAD(NDIMENSION),NODE(ELEMENT(IELT).NODE(I)).VELOCITY(NDIMENSION))
                NODE(ELEMENT(IELT).NODE(I)).IGRAD=0.D0;NODE(ELEMENT(IELT).NODE(I)).VELOCITY=0.D0
            ENDIF
        ENDDO
    ENDIF

CONTAINS
        
    
ENDSUBROUTINE
    
SUBROUTINE ZTE_Q_K_UPDATE(STEPDIS,IELT,ISTEP,IITER)  
    USE solverds    
    IMPLICIT NONE
    INTEGER,INTENT(IN)::IELT,ISTEP,IITER
    REAL(8),INTENT(IN)::STEPDIS(NDOF)
    
    REAL(8),ALLOCATABLE::NH1(:)
    INTEGER::CNODE1(6),ND1,NGP1,N1,J
    REAL(8)::K1(4),IGRAD1,L1,KR1,KP1,Z1,V1,A1,T1,NQ1(8)
   
        
 
    NH1=STEPDIS(ELEMENT(IELT).G)
    
    NGP1=element(IELT).NGP
	N1=element(IELT).NNUM/2
    IF(ELEMENT(IELT).ET==ZT4_SPG) THEN
        CNODE1(1:2)=[4,3]
        ND1=2
    ELSE
        CNODE1(1:3)=[4,5,6]
        ND1=3
    ENDIF
    L1=MATERIAL(ELEMENT(IELT).MAT).PROPERTY(14)
    KP1=MATERIAL(ELEMENT(IELT).MAT).PROPERTY(1)
    A1=ELEMENT(IELT).PROPERTY(1)/N1
	do j=1,N1
                   
		      
					
		!gradient ELEMENT(IELT).SIGN
        IGRAD1=(NH1(CNODE1(J))-NH1(J))/L1


        !节点梯度
        IF(ELEMENT(IELT).ET==ZT4_SPG) THEN
		    element(IELT).igrad(1:nd1,NGP1+J)=ELEMENT(IELT).G2L(2,:)*IGRAD1
        ELSE
            element(IELT).igrad(1:nd1,NGP1+J)=ELEMENT(IELT).G2L(3,:)*IGRAD1
        ENDIF
        
        element(IELT).igrad(1:nd1,NGP1+CNODE1(J))=element(IELT).igrad(1:nd1,NGP1+J)
        
		!velocity
        Z1=NODE(ELEMENT(IELT).NODE(J)).COORD(NDIMENSION)
        IF(.NOT.ALLOCATED(element(IELT).kr)) THEN
            allocate(element(IELT).Kr(NGP1+N1*2),element(IELT).Mw(NGP1+N1*2))
			element(IELT).kr=1.0D0
			element(IELT).Mw=0.0D0
        ENDIF
        IF(NH1(J)-Z1<0.D0.OR.NH1(CNODE1(J))-Z1<0.D0) THEN
            KR1=1E-3
        ELSE
            KR1=1.D0
        ENDIF
        KR1=(KR1+element(IELT).kr(NGP1+J))/2.0
        element(IELT).kr(NGP1+J)=KR1
        element(IELT).kr(NGP1+CNODE1(J))=KR1
        V1=-KR1*KP1*IGRAD1        

        !节点速度
        IF(ELEMENT(IELT).ET==ZT4_SPG) THEN
		    element(IELT).velocity(1:nd1,NGP1+J)=ELEMENT(IELT).G2L(2,:)*V1
        ELSE
            element(IELT).velocity(1:nd1,NGP1+J)=ELEMENT(IELT).G2L(3,:)*V1
        ENDIF
        
        element(IELT).velocity(1:nd1,NGP1+CNODE1(J))=element(IELT).velocity(1:nd1,NGP1+J)
        
		!节点流量
        NQ1(J)=V1*A1;NQ1(CNODE1(J))=-NQ1(J)
        
        T1=KR1*KP1*A1/L1
        if(solver_control.bfgm==continuum.or.solver_control.bfgm==consistent.or.solver_control.bfgm==inistress) then
            isref_spg=1
            if(j==1) element(IELT).km=0.0D0
			ELEMENT(IELT).KM(J,J)=T1
            ELEMENT(IELT).KM(CNODE1(J),CNODE1(J))=T1
            ELEMENT(IELT).KM(J,CNODE1(J))=-T1
            ELEMENT(IELT).KM(CNODE1(J),J)=-T1
		end if
					
    end do
    
    !NQ1=MATMUL(ELEMENT(IELT).KM,NH1)
    IF(.NOT.ALLOCATED(element(IELT).flux)) ALLOCATE(element(IELT).flux(N1*2))    
    element(IELT).flux=NQ1(1:N1*2)
    
ENDSUBROUTINE
    
SUBROUTINE ZT_SPFACE_BC()
    use solverds,only:bc_tydef,nseep,NumNSeep,NODE,NNUM,NDIMENSION,ISPF
    implicit none
    !integer,ALLOCATABLE::ISPF(:,:) 
    TYPE(bc_tydef),ALLOCATABLE::NSEEP1(:)
    INTEGER::N1,I,N2,N3
    
    IF(.not.ALLOCATED(ISPF)) RETURN
    
    N1=COUNT(ISPF(:,1)>0)
    ALLOCATE(NSEEP1(N1))
    N2=0
    DO I=1,NNUM
        N3=ISPF(I,1)
        IF(N3>0)THEN
            N2=N2+1
            NSEEP1(N2).NODE=N3
            NSEEP1(N2).SF=ISPF(I,2)
            NSEEP1(N2).DOF=4
            NSEEP1(N2).ISDEAD=0
            NSEEP1(N2).VALUE=NODE(N3).COORD(NDIMENSION)
        ENDIF
    ENDDO
    NSEEP=[NSEEP,NSEEP1]
    NumNSeep=NumNSeep+N1
    DEALLOCATE(NSEEP1,ISPF)
ENDSUBROUTINE
    
SUBROUTINE POREFLOW_Q_K_UPDATE(STEPDIS,IEL,ISTEP,IITER)
    USE solverds

    IMPLICIT NONE
    INTEGER,INTENT(IN)::IEL,ISTEP,IITER
    REAL(8),INTENT(IN)::STEPDIS(NDOF)
    INTEGER::NDOF1,ET1,MODEL1
    REAL(8)::NHEAD1(2),PI1,X2(3),L1,D1,RW1,FD1,REW1,AREA1,VA1,QA1,G1,RE1,DIS1,LAMDA1
    
    
    NDOF1=ELEMENT(IEL).NDOF;ET1=ELEMENT(IEL).ET
    NHEAD1(1:NDOF1)=STEPDIS(ELEMENT(IEL).G)
    PI1=PI(); 
    
    X2=NODE(ELEMENT(IEL).NODE(1)).COORD-NODE(ELEMENT(IEL).NODE(2)).COORD
    DIS1=NORM2(X2)/2.0 !half length of the wellbore
    rw1=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(1)
    if(isporeflow) rw1=element(iel).property(2)/2.0
    FD1=ELEMENT(IEL).PROPERTY(1)
    REW1=0.D0;
        
    MODEL1=INT(MATERIAL(ELEMENT(IEL).MAT).PROPERTY(6))
    
    L1=2*DIS1
    D1=2*rw1
    AREA1=PI1*D1**2/4.0
    QA1=ABS((NHEAD1(1)-NHEAD1(2))*ELEMENT(IEL).KM(1,1))
    VA1=QA1/AREA1
    
       
    RE1=Re_W(VA1,D1,MATERIAL(ELEMENT(IEL).MAT).PROPERTY(2))
    ELEMENT(IEL).PROPERTY(3)=RE1 
    
    IF(model1/=4) THEN

        g1=solver_control.get_g() 
        
        IF(MODEL1==3) THEN
            FD1=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(8)
        ELSE
            FD1=FD_PF(RE1,MATERIAL(ELEMENT(IEL).MAT).PROPERTY(3),REW=0.0D0,POROSITY=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(9),MODEL=MODEL1) 
        ENDIF
        
        FD1=FD1*L1*QA1/(2*D1*G1*AREA1**2)
        IF(RE1<2260.AND.MODEL1==0)  FD1=ELEMENT(IEL).FD !PIPE-FLOW
          
    END IF 
    
    ELEMENT(IEL).PROPERTY(1)=FD1;
    
    IF(ABS(FD1)<1E-15) THEN
        FD1=1.D-15
    ENDIF
    ELEMENT(IEL).KM=RESHAPE([1.D0,-1.D0,-1.D0,1.D0],([2,2]))    
    
    
	IF(MINVAL(NHEAD1-NODE(ELEMENT(IEL).NODE(1:NDOF1)).COORD(NDIMENSION))>-1.D-6) THEN
        LAMDA1=1.0D0 !承压
    ELSE
	    LAMDA1=1.D-3 !某一节点在非饱和区                 						
    ENDIF 
    ELEMENT(IEL).KM=ELEMENT(IEL).KM/FD1*(ELEMENT(IEL).KR(1)+LAMDA1)/2.0
    ELEMENT(IEL).KR(1)=LAMDA1
    
    IF(.NOT.ALLOCATED(ELEMENT(IEL).FLUX)) ALLOCATE(ELEMENT(IEL).FLUX(ELEMENT(IEL).NNUM))
    
    ELEMENT(IEL).FLUX=MATMUL(ELEMENT(IEL).KM,NHEAD1(1:ELEMENT(IEL).NDOF))
    ELEMENT(IEL).CDT=ABS(ELEMENT(IEL).FLUX(1)/AREA1)


    
END SUBROUTINE    
    
SUBROUTINE POREFLOW_Q_K_UPDATE2(STEPDIS,IEL,ISTEP,IITER)
    USE solverds
    USE PoreNetWork
    IMPLICIT NONE
    INTEGER,INTENT(IN)::IEL,ISTEP,IITER
    REAL(8),INTENT(IN)::STEPDIS(NDOF)
    INTEGER::NDOF1,ET1,MODEL1
    REAL(8)::NHEAD1(2),PI1,X2(3),L1,D1,RW1,FD1,REW1,AREA1,VA1,QA1,G1,RE1,DIS1,LAMDA1,FD0
    
    
    NDOF1=ELEMENT(IEL).NDOF;ET1=ELEMENT(IEL).ET
    NHEAD1(1:NDOF1)=STEPDIS(ELEMENT(IEL).G)
    AREA1=PI()*ELEMENT(IEL).PROPERTY(2)**2/4.0
    
    FD0=ELEMENT(IEL).FD
    CALL ELEMENT(IEL).RIJ()
    FD1=ELEMENT(IEL).FD
    IF(ABS(FD1)<1E-15) THEN
        FD1=1.D-15
    ENDIF
    ELEMENT(IEL).KM=RESHAPE([1.D0,-1.D0,-1.D0,1.D0],([2,2]))    
    
    !假定承压
	!IF(MINVAL(NHEAD1-NODE(ELEMENT(IEL).NODE(1:NDOF1)).COORD(NDIMENSION))>-1.D-6) THEN
 !       LAMDA1=1.0D0 !承压
 !   ELSE
	!    LAMDA1=1.D-3 !某一节点在非饱和区                 						
 !   ENDIF 
    ELEMENT(IEL).KM=2.0*ELEMENT(IEL).KM/(FD0+FD1)   !*(ELEMENT(IEL).KR(1)+LAMDA1)/2.0
    !ELEMENT(IEL).KR(1)=LAMDA1
    
    IF(.NOT.ALLOCATED(ELEMENT(IEL).FLUX)) ALLOCATE(ELEMENT(IEL).FLUX(ELEMENT(IEL).NNUM))
    
    ELEMENT(IEL).FLUX=MATMUL(ELEMENT(IEL).KM,NHEAD1(1:ELEMENT(IEL).NDOF))
    !ELEMENT(IEL).CDT=ELEMENT(IEL).PROPERTY(4)/ABS(ELEMENT(IEL).FLUX(1)/AREA1)


    
END SUBROUTINE    