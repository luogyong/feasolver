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
        if(solver_control.wellmethod==0.or.solver_control.wellmethod==2) then 
            CALL SPHFLOW_Q_K_UPDATE3(STEPDIS,IENUM,ISTEP,IITER)
        elseIF(solver_control.wellmethod==1) THEN
            CALL SPHFLOW_Q_K_UPDATE(STEPDIS,IENUM,ISTEP,IITER)
        ELSE
            CALL SPHFLOW_Q_K_SPMETHOD(STEPDIS,IENUM,ISTEP,IITER)
        endif        
        
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


SUBROUTINE WELLBORE_Q_K_UPDATE(STEPDIS,IEL,ISTEP,IITER)
    USE solverds
    USE MESHGEO,ONLY:SNADJL
    IMPLICIT NONE
    INTEGER,INTENT(IN)::IEL,ISTEP,IITER
    REAL(8),INTENT(IN)::STEPDIS(NDOF)
    INTEGER::I,J,K,I0,IN1,IN2,IN3,N1,N2,IEL1,IP1,MODEL1,II1,NDOF1,ET1
    REAL(8)::H1,X1(3),R1,Q1(4),PI1,LAMDA1,K1,DIS1,DH1,A1(5),NHEAD1(5),KM1(5,5),HZ1,RA1,X2(3),&
    W1,TW1,RE1,VA1,QR1,FD1,L1,D1,QA1,VR1,AREA1,g1,FACC1,REW1,cf1,QN1(5),KX1,KY1,KZ1,SITA1(3),LP1,&
    KR1,XC1(3),VW1(3),T1,C1,B1,PO1,VS1,CF2
    LOGICAL::ISO1=.FALSE.,ISIN1=.FALSE.
    
    NDOF1=ELEMENT(IEL).NDOF;ET1=ELEMENT(IEL).ET
    NHEAD1(1:NDOF1)=STEPDIS(ELEMENT(IEL).G)
    PI1=PI(); 
    
    X2=NODE(ELEMENT(IEL).NODE(1)).COORD-NODE(ELEMENT(IEL).NODE(2)).COORD
    DIS1=NORM2(X2)/2.0 !half length of the wellbore
    VW1=X2/DIS1/2.0 !wellbore unit vection    
    
    FD1=ELEMENT(IEL).PROPERTY(1)
    FACC1=0.0D0;REW1=0.D0
        
    MODEL1=INT(MATERIAL(ELEMENT(IEL).MAT).PROPERTY(6))
 
    
    IF(ET1/=WELLBORE_SPGFACE.and.model1/=4) THEN
        L1=2*DIS1
        D1=2*MATERIAL(ELEMENT(IEL).MAT).PROPERTY(1)
        if(isporeflow>0) D1=element(iel).property(2)
        AREA1=PI1*D1**2/4.0
        QA1=ABS((NHEAD1(1)-NHEAD1(2))*ELEMENT(IEL).KM(1,1))
        VA1=QA1/AREA1
        !if(ABS(MATERIAL(ELEMENT(IEL).MAT).PROPERTY(5))<1.E-7) THEN
        !    g1=73156608000.00 
        !else
        !    g1=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(5)
        !endif
        g1=solver_control.get_g() 
        RE1=Re_W(VA1,D1,MATERIAL(ELEMENT(IEL).MAT).PROPERTY(2))
        if(isporeflow>0) ELEMENT(IEL).PROPERTY(3)=RE1
        !ACCELERATION LOSS.
        
        IF(ELEMENT(IEL).ET==WELLBORE) THEN
            QR1=(NHEAD1(3)-NHEAD1(2))*ELEMENT(IEL).KM(3,3)+(NHEAD1(4)-NHEAD1(1))*ELEMENT(IEL).KM(4,4)
            FACC1=QR1/(G1*AREA1**2)
            IF(MODEL1==1.AND.ABS(QR1)>1.E-7) THEN
            ! Siwoń, Z. (1987), Solutions for Lateral Inflow in Perforated Conduits, Journal of Hydraulic Engineering, 113(9), 1117-1132.
                PO1=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(9)
                B1=10/(1000*PO1)**4.2+4/10**7                
                C1=1.05+1./(1.235+B1*(QA1/QR1*4*L1*PO1/D1)**2)
            ELSE
                C1=2.0
            ENDIF
            FACC1=C1*FACC1
        ELSE
            QR1=0.D0
        ENDIF
        
        IF(MODEL1==3) THEN
            FD1=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(8)
        ELSE

            
            VR1=QR1/(L1*PI1*D1)
            REW1=Re_W(VR1,D1,MATERIAL(ELEMENT(IEL).MAT).PROPERTY(2))           
            FD1=FD_PF(RE1,MATERIAL(ELEMENT(IEL).MAT).PROPERTY(3),REW=REW1,POROSITY=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(9),MODEL=MODEL1) 
        ENDIF
        
        FD1=FD1*L1*QA1/(2*D1*G1*AREA1**2)
        IF(RE1<2260.AND.MODEL1==0)  FD1=ELEMENT(IEL).FD !PIPE-FLOW
        
        !FRICTION PART, CONSIDERING THE TURBULENT AND POROUS SIDE EFFECT
        !IF(MODEL1==2) THEN
        !    !ASSUME LENGHTH UNIT IS M
        !    IF(MATERIAL(ELEMENT(IEL).MAT).PROPERTY(8)>0.D0) THEN
        !        FD1=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(8)
        !    ELSE
        !        !IF(ABS(MATERIAL(ELEMENT(IEL).MAT).PROPERTY(4))<1.D-7) THEN
        !        !    CF1=1./24/3600 !to m/s
        !        !ELSE
        !        !    CF1=1.0D0
        !        !ENDIF
        !        !RE1=Re_W(VA1*CF1,D1,MATERIAL(ELEMENT(IEL).MAT).PROPERTY(2))
        !        VR1=QR1/(L1*PI1*D1)
        !        REW1=Re_W(VR1*CF1,D1,MATERIAL(ELEMENT(IEL).MAT).PROPERTY(2))
        !        FD1=FD_PF(RE1,MATERIAL(ELEMENT(IEL).MAT).PROPERTY(3),REW1) 
        !    ENDIF
        !    FD1=FD1*L1*QA1/(2*D1*G1*AREA1**2)
        !ENDIF            
    END IF 
    
    ELEMENT(IEL).PROPERTY(1)=FD1;ELEMENT(IEL).PROPERTY(4)=FACC1; 
    
    IF(ABS(FD1)<1E-15) THEN
        FD1=1.D-15
    ENDIF
    
    
    A1(1)=1./(FACC1+FD1)
    IF(A1(1)<0.D0) A1(1)=1./FD1    
    
    
    

    
    IF(ET1==PIPE2) THEN
        KM1(1,1)=A1(1);KM1(2,2)=A1(1);KM1(1,2)=-A1(1);KM1(2,1)=-A1(1)
    ELSE
        !WRITE(99,10) 
        
        Q1=0.D0
        
        QN1(1:NDOF1)=MATMUL(ELEMENT(IEL).KM,NHEAD1(1:NDOF1))
        QWAN(2,ELEMENT(IEL).NODE(1:NDOF1))=QWAN(2,ELEMENT(IEL).NODE(1:NDOF1))+QN1(1:NDOF1)        
        
        DO I=1,2
    
            !IF(NHEAD1(2+I)-NODE(ELEMENT(IEL).NODE(2+I)).COORD(NDIMENSION)<0.D0) THEN        
            !    A1(I+1)=1.D-10
            !    A1(1)=1.D-10 !!!!!
            !    CYCLE
            !ENDIF
        
            IN1=ELEMENT(IEL).NODE(2+I)
            IF(ELEMENT(IEL).ISTOPO>0) IN1=ELEMENT(IEL).NODE(5+I)
            !IN2=ELEMENT(IEL).NODE(2+MOD(I,2)+1)
            !DIS1=NORM2(NODE(IN1).COORD-NODE(IN2).COORD)/2.
            IF(I==1) THEN
                N2=2 !2-3
                IF(ELEMENT(IEL).ISTOPO>0) IN1=ELEMENT(IEL).NODE(6)
            ELSE
                N2=1 !1-4
                IF(ELEMENT(IEL).ISTOPO>0) IN1=ELEMENT(IEL).NODE(5)
            ENDIF
            !IF(ET1==WELLBORE) THEN

            DH1=STEPDIS(NODE(IN1).DOF(4))-STEPDIS(NODE(ELEMENT(IEL).NODE(N2)).DOF(4))
            !ELSE
            !    !WELLBORE_SPGFACE
            !    
            !    DH1=STEPDIS(NODE(IN1).DOF(4))-NODE(ELEMENT(IEL).NODE(N2)).COORD(NDIMENSION)
            !ENDIF
            
            DO II1=1,2
                IP1=0;RA1=0.D0;W1=0;TW1=0;     
                IF(II1==1) THEN
                    ISIN1=.TRUE.
                ELSE
                    ISIN1=.FALSE.
                ENDIF
                DO J=1,SNADJL(IN1).ENUM
                    IEL1=SNADJL(IN1).ELEMENT(J)
                    IF(ELEMENT(IEL1).ET==WELLBORE.OR.ELEMENT(IEL1).ET==PIPE2 &
                    .OR.ELEMENT(IEL1).ET==WELLBORE_SPGFACE) CYCLE !ITSELF EXCLUDED 
                    IF(ELEMENT(IEL1).ET==SPHFLOW.OR.ELEMENT(IEL1).ET==SEMI_SPHFLOW) CYCLE !ITSELF EXCLUDED
                    IF(ELEMENT(IEL1).EC/=SPG) CYCLE
                
                    KR1=0.D0
                    IF(NORM2(MATERIAL(ELEMENT(IEL1).MAT).PROPERTY(1:NDIMENSION)-MATERIAL(ELEMENT(IEL1).MAT).PROPERTY(1))<1E-10) THEN
                        K=5
                        ISO1=.TRUE.
                    ELSE
                        K=1
                        ISO1=.FALSE.
                    ENDIF
                
                    DO K=K,5
                        X1(1)=DOT_PRODUCT(NODE(ELEMENT(IEL1).NODE(1:4)).COORD(1),GQ_TET4_O3(1:4,K))
                        X1(2)=DOT_PRODUCT(NODE(ELEMENT(IEL1).NODE(1:4)).COORD(2),GQ_TET4_O3(1:4,K))
                        X1(3)=DOT_PRODUCT(NODE(ELEMENT(IEL1).NODE(1:4)).COORD(3),GQ_TET4_O3(1:4,K))
                        IF(K==5) THEN !K==5 IS CENTER OF THE ELEMENT
                            H1=DOT_PRODUCT(STEPDIS(ELEMENT(IEL1).G(1:4)),GQ_TET4_O3(1:4,K))
                            XC1=X1
                        ENDIF
                    
                        X1=X1-NODE(ELEMENT(IEL).NODE(3)).COORD !ORIGIN IS THE THIRD NODE. 

                        LP1=DOT_PRODUCT(X1,VW1) 
                        X2=LP1*VW1 !X1在井线上的投影 
                        X2=X1-X2 !投影线向量，为此单元的流量的方向矢量
                        R1=NORM2(X2)                   
                
                        !DIRECTIONAL K
                        SITA1=X2/R1            
                        K1=0.D0
                        DO I0=1,3
                            K1=K1+1.0/MATERIAL(ELEMENT(IEL1).MAT).PROPERTY(I0)*(SITA1(I0))**2
                        ENDDO
                        IF(ISO1) THEN
                            KR1=KR1+1/K1
                        ELSE
                            KR1=KR1+1/K1*GQ_TET4_O3(5,K)*6.D0 
                        ENDIF
                    ENDDO
                
                    CALL lamda_spg(IEL1,H1,XC1(NDIMENSION),lamda1,ISTEP)
                    KR1=LAMDA1*KR1
                
                    IF(R1<=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(1)) THEN
                        PRINT *,'THE ELEMENT SIZE AROUND THE NODE I SEEMS TO BE TOO SMALL. I=', NODE(ELEMENT(IEL).NODE(IN1)).COORD
                        CYCLE 
                        
                    ENDIF
                    
                    IF((LP1<0.D0.OR.LP1>2*DIS1).AND.ISIN1) CYCLE !跳过不在单元范围内的节点
               
                    HZ1=NHEAD1(2)+(NHEAD1(1)-NHEAD1(2))/(2*DIS1)*LP1
                    !要考虑井损
                    RA1=(H1-HZ1)/(LOG(R1/MATERIAL(ELEMENT(IEL).MAT).PROPERTY(1))/(ELEMENT(IEL).PROPERTY(6)*KR1*DIS1)+ELEMENT(IEL).PROPERTY(5))
                    W1=ELEMENT(IEL1).ANGLE(SNADJL(IN1).SUBID(J))
                    Q1(I+2)=Q1(I+2)+RA1*W1
                    TW1=TW1+W1;IP1=IP1+1
                    !WRITE(99,20) ISTEP,IITER,IEL,I+2,IP1,X2,X1,R1,H1,HZ1,RA1,W1
                ENDDO
                
                IF(IP1>0) EXIT
                
                IF(II1==1) THEN
                   PRINT *, "井单元I周边单元的形心的投影均不在单元I内，尝试区域外的单元.I=",iel                 
                ELSE
                   PRINT *, "井单元I周边没有相邻的单元.请检查.I=",iel 
                   STOP
                ENDIF
            ENDDO
            Q1(I+2)=Q1(I+2)/TW1
            !A1(I+1)=ABS(Q1(I+2)/DH1)
            
            QWAN(1,ELEMENT(IEL).NODE(I+2))=QWAN(1,ELEMENT(IEL).NODE(I+2))+Q1(I+2)
            
            IF(ABS(Q1(I+2))<1E-5) THEN
                Q1(I+2)=1.D-5
            ENDIF
            
            IF(ABS(DH1)<1.D-5) DH1=1.D-5
            
            !ELEMENT(IEL).PROPERTY(I+1)=ABS(2*DH1/(Q1(I+2)+QN1(I+2))) !!!包含了井损。
            ELEMENT(IEL).PROPERTY(I+1)=ABS(DH1/Q1(I+2))!!!包含了井损。
            
          
            A1(I+1)=1.0/(ELEMENT(IEL).PROPERTY(I+1))
        ENDDO
        
        IF(ANY(ISNAN(A1(1:3)))) THEN
            PRINT *, 'NAN'
            STOP "ERRORS IN WELLBORE_Q_K_UPDATE"
        ENDIF
        
        KM1(1:NDOF1,1:NDOF1)=KM_WELLBORE(A1(1),A1(2),A1(3))
        KM1(1:NDOF1,1:NDOF1)=(ELEMENT(IEL).KM+KM1(1:NDOF1,1:NDOF1))*0.5
    ENDIF
    
    
    IF(.NOT.ALLOCATED(ELEMENT(IEL).FLUX)) ALLOCATE(ELEMENT(IEL).FLUX(ELEMENT(IEL).NNUM))
    
    ELEMENT(IEL).FLUX=MATMUL(KM1(1:ELEMENT(IEL).NDOF,1:ELEMENT(IEL).NDOF),NHEAD1(1:ELEMENT(IEL).NDOF))
    
    ELEMENT(IEL).KM=KM1(1:ELEMENT(IEL).NDOF,1:ELEMENT(IEL).NDOF)
    !if(solver_control.bfgm==continuum.or.solver_control.bfgm==consistent) then
    !          							
    !end if

10  FORMAT("ISTEP,IITER,IEL,NODE1,P1,XG,YG,ZG,XL,YL,ZL,RA,HA,HW,Qunit,WGT   //WELLBORE INFO")
20  FORMAT(5I7,11F15.7)
    
END SUBROUTINE

SUBROUTINE WELLBORE_Q_K_UPDATE_SPMETHOD(STEPDIS,IEL,ISTEP,IITER)
    USE solverds
    USE MESHGEO,ONLY:SNADJL,GETVAL_SOLVER
    USE SolverMath,ONLY:Nint_gauss
    IMPLICIT NONE
    INTEGER,INTENT(IN)::IEL,ISTEP,IITER
    REAL(8),INTENT(IN)::STEPDIS(NDOF)
    INTEGER::I,J,K,I0,IN1,IN2,IN3,N1,N2,IEL1,IP1,MODEL1,II1,NDOF1,ET1,ANODE1,isok
    INTEGER,ALLOCATABLE::NSP1(:),NSP2(:)
    REAL(8)::H1,X1(3),R1,Q1(4),PI1,LAMDA1,K1,DIS1,DH1,A1(5),NHEAD1(5),KM1(5,5),HZ1,RA1,X2(3),&
    W1,TW1,RE1,VA1,QR1,FD1,L1,D1,QA1,VR1,AREA1,g1,FACC1,REW1,cf1,QN1(5),KX1,KY1,KZ1,SITA1(3),LP1,&
    KR1,XC1(3),VW1(3),T1,C1,B1,PO1,VS1,DIS2,H2,CF2
    REAL(8)::SPH1(200),NHEAD2(10,1),f1,scale1,Ru1,rwu1,rw1,val1
    LOGICAL::ISO1=.FALSE.,ISIN1=.FALSE.
    COMPLEX(8)::Z1,U1
    
    
    NDOF1=ELEMENT(IEL).NDOF;ET1=ELEMENT(IEL).ET
    NHEAD1(1:NDOF1)=STEPDIS(ELEMENT(IEL).G)
    PI1=PI(); 
    
    X2=NODE(ELEMENT(IEL).NODE(1)).COORD-NODE(ELEMENT(IEL).NODE(2)).COORD
    DIS1=NORM2(X2)/2.0 !half length of the wellbore
    VW1=X2/DIS1/2.0 !wellbore unit vection    
    rw1=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(1)
    if(isporeflow) rw1=element(iel).property(2)/2.0
    FD1=ELEMENT(IEL).PROPERTY(1)
    FACC1=0.0D0;REW1=0.D0
        
    MODEL1=INT(MATERIAL(ELEMENT(IEL).MAT).PROPERTY(6))
 
    
    IF(ET1/=WELLBORE_SPGFACE.and.model1/=4) THEN
        L1=2*DIS1
        D1=2*rw1
        AREA1=PI1*D1**2/4.0
        QA1=ABS((NHEAD1(1)-NHEAD1(2))*ELEMENT(IEL).KM(1,1))
        VA1=QA1/AREA1
        !if(ABS(MATERIAL(ELEMENT(IEL).MAT).PROPERTY(5))<1.E-7) THEN
        !    g1=73156608000.00 
        !else
        !    g1=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(5)
        !endif
        g1=solver_control.get_g()    
        !ACCELERATION LOSS.
        RE1=Re_W(VA1,D1,MATERIAL(ELEMENT(IEL).MAT).PROPERTY(2))
        if(isporeflow>0) ELEMENT(IEL).PROPERTY(3)=RE1
        IF(ELEMENT(IEL).ET==WELLBORE) THEN
            QR1=(NHEAD1(3)-NHEAD1(2))*ELEMENT(IEL).KM(3,3)+(NHEAD1(4)-NHEAD1(1))*ELEMENT(IEL).KM(4,4)
            FACC1=QR1/(G1*AREA1**2)
            IF(MODEL1==1.AND.ABS(QR1)>1.E-7) THEN
            ! Siwoń, Z. (1987), Solutions for Lateral Inflow in Perforated Conduits, Journal of Hydraulic Engineering, 113(9), 1117-1132.
                PO1=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(9)
                B1=10/(1000*PO1)**4.2+4/10**7                
                C1=1.05+1./(1.235+B1*(QA1/QR1*4*L1*PO1/D1)**2)
            ELSE
                C1=2.0
            ENDIF
            FACC1=C1*FACC1
        ELSE
            QR1=0.D0
        ENDIF
        
        IF(MODEL1==3) THEN
            FD1=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(8)
        ELSE
            !IF(ABS(MATERIAL(ELEMENT(IEL).MAT).PROPERTY(4))<1.D-7) THEN
            !    CF1=1./24./3600. !to m/s
            !ELSE
            !    CF1=1.0D0
            !ENDIF
            !time to s 
            VR1=QR1/(L1*PI1*D1)
            REW1=Re_W(VR1,D1,MATERIAL(ELEMENT(IEL).MAT).PROPERTY(2))
            FD1=FD_PF(RE1,MATERIAL(ELEMENT(IEL).MAT).PROPERTY(3),REW=REW1,POROSITY=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(9),MODEL=MODEL1) 
        ENDIF
        
        FD1=FD1*L1*QA1/(2*D1*G1*AREA1**2)
        IF(RE1<2260.AND.MODEL1==0)  FD1=ELEMENT(IEL).FD !PIPE-FLOW
        
          
    END IF 
    
    ELEMENT(IEL).PROPERTY(1)=FD1;ELEMENT(IEL).PROPERTY(4)=FACC1; 
    
    IF(ABS(FD1)<1E-15) THEN
        FD1=1.D-15
    ENDIF
    
    
    A1(1)=1./(FACC1+FD1)
    IF(A1(1)<0.D0) A1(1)=1./FD1    

    
    IF(ET1==PIPE2) THEN
        KM1(1,1)=A1(1);KM1(2,2)=A1(1);KM1(1,2)=-A1(1);KM1(2,1)=-A1(1)
    ELSE
        !WRITE(99,10) 
        QN1(1:NDOF1)=MATMUL(ELEMENT(IEL).KM,NHEAD1(1:NDOF1))
        QWAN(2,ELEMENT(IEL).NODE(1:NDOF1))=QWAN(2,ELEMENT(IEL).NODE(1:NDOF1))+QN1(1:NDOF1)
        
        Q1=0.D0
        
        DO I=1,SIZE(ELEMENT(IEL).NSPLOC)
            N1=ELEMENT(IEL).NSPLOC(I)
            IF(N1<1) cycle
            N2=ELEMENT(N1).NDOF
            NHEAD2(1:N2,1)=STEPDIS(ELEMENT(N1).G)
            CALL GETVAL_SOLVER(ELEMENT(IEL).A12(:,I),N1,SPH1(I:I),NHEAD2(1:N2,:))
            
        ENDDO
        
        N1=SIZE(ELEMENT(IEL).NSPLOC)/3 
        !同一单元各采样点到井轴的距离都相等。
!        N2=MAXLOC(ELEMENT(IEL).NSPLOC(1:N1),dim=1) !排除不在范围内的采样点。
!        R1=NORM2(ELEMENT(IEL).A12(:,N2)-NODE(ELEMENT(IEL).NODE(3)).COORD)
        
        !IF(R1<=rw1) THEN
        !    STOP "SAMPLE POINT IS TOO CLOSE TO THE WELLAXIS. ERROR IN WELLBORE_Q_K_UPDATE_SPMETHOD."            
        !ENDIF
        

        
        DO I=1,2
    
            
            IN1=ELEMENT(IEL).NODE(2+I)
            !IN2=ELEMENT(IEL).NODE(2+MOD(I,2)+1)
            !DIS1=NORM2(NODE(IN1).COORD-NODE(IN2).COORD)/2.
            IF(I==1) THEN
                N2=2 !2-3                
                NSP1=[1:N1,2*N1+1:3*N1] !内圈
                !NSP2=[1:N1,2*N1+1:3*N1]+3*N1 !外圈
            ELSE
                N2=1 !1-4
                NSP1=[N1+1:3*N1]
                !NSP2=[N1+1:3*N1]+3*N1 !外圈
            ENDIF
            !IF(ET1==WELLBORE) THEN

            DH1=STEPDIS(NODE(IN1).DOF(4))-STEPDIS(NODE(ELEMENT(IEL).NODE(N2)).DOF(4))
  
            

            IP1=0;RA1=0.D0;W1=0;TW1=0;     
            
            !本端采样节点+中间采样节点           
            
            DO J=1,2*N1
                IF(ELEMENT(IEL).NSPLOC(NSP1(J))<=0) CYCLE
                !采样点
                X1=ELEMENT(IEL).A12(:,NSP1(J))
                H1=SPH1(NSP1(J))
                !H2=SPH1(NSP2(J))
                !IF(H1<X1(NDIMENSION)) CYCLE
                
                !采样点投影点
                IF(J<=N1) THEN
                    XC1=NODE(ELEMENT(IEL).NODE(N2)).COORD
                ELSE
                    XC1=(NODE(ELEMENT(IEL).NODE(1)).COORD+NODE(ELEMENT(IEL).NODE(2)).COORD)/2.0
                ENDIF
                
                LP1=NORM2(XC1-NODE(ELEMENT(IEL).NODE(3)).COORD)
                
                X2=X1-XC1
                R1=NORM2(X2)
                
                !DIRECTIONAL K
                SITA1=X2/R1            
                K1=0.D0
                DO I0=1,3
                    !假定中间采样点所在的单元材料能代表为井周材料(两端节点可能是材料分界面，材料可能不确定)。
                    IF(J<=N1) THEN
                        IEL1=ELEMENT(IEL).NSPLOC(NSP1(N1+J))
                        IF(IEL1==0) IEL1=ELEMENT(IEL).NSPLOC(NSP1(J))
                    ELSE
                        IEL1=ELEMENT(IEL).NSPLOC(NSP1(J))
                    ENDIF
                    K1=K1+1.0/MATERIAL(ELEMENT(IEL1).MAT).PROPERTY(I0)*(SITA1(I0))**2
                ENDDO
                KR1=1/K1
                RU1=R1;
                RWU1=rw1                
                val1=log(ru1/rwu1)
                !VAL1=LOG((RU1+0.1)/(RU1))
                
                kx1=MATERIAL(ELEMENT(IEL1).MAT).PROPERTY(1)
                ky1=MATERIAL(ELEMENT(IEL1).MAT).PROPERTY(2)
                IF(KX1/=KY1.and.solver_control.wellaniso>0) THEN
                    KR1=(kx1*ky1)**0.5
                    if(solver_control.wellaniso==1) then 
                        !REf:Fitts, C.R., Exact Solution for Two-Dimensional Flow to a Well in an Anisotropic Domain. Groundwater, 2006. 44(1): p. 99-101.
                        !各向异性，转化为各向同性进行处理。
                        scale1=(max(kx1,ky1)/min(kx1,ky1))**0.5
                        f1=rw1*(scale1**2-1)**0.5
                        IF(KX1<KY1) THEN
                            Z1=CMPLX(X2(1)*SCALE1,X2(2))
                        ELSE
                            Z1=CMPLX(X2(2)*SCALE1,X2(1))
                        ENDIF
                        U1=0.5*(Z1+(Z1-F1)**0.5*(Z1+F1)**0.5)
                        Ru1=ABS(U1)
                        IF(KX1<KY1) THEN
                            Z1=CMPLX(rw1*SCALE1,0)
                        ELSE
                            Z1=CMPLX(0,rw1)
                        ENDIF
                        U1=0.5*(Z1+(Z1-F1)**0.5*(Z1+F1)**0.5)
                        Rwu1=ABS(U1)
                        val1=log(RU1/RWU1)
                    else
                        !王建荣, 各向异性介质中大口径潜水完整井流公式. 勘察科学技术, 1991(02): p. 10-13.
                        call Nint_gauss(fun_lnkr,0.,PI1/2,val1,1e-7,isok)
                       
                        val1=val1/PI1-log(rw1)                
                        val1=val1+0.5*log(x2(1)**2/kx1+x2(2)**2/ky1)
                    endif
                ENDIF

        

                
                !KR1=(MATERIAL(ELEMENT(IEL1).MAT).PROPERTY(1)*MATERIAL(ELEMENT(IEL1).MAT).PROPERTY(2)*MATERIAL(ELEMENT(IEL1).MAT).PROPERTY(3))**(1/3)
                IF(H1>X1(NDIMENSION)) THEN
                    LAMDA1=1.0
                ELSE
                    LAMDA1=1.D-3
                ENDIF
                KR1=LAMDA1*KR1

             
                HZ1=NHEAD1(2)+(NHEAD1(1)-NHEAD1(2))/(2*DIS1)*LP1
                
                !要考虑井损
                !DIS2=DIS1/(MATERIAL(ELEMENT(IEL1).MAT).PROPERTY(3))**0.5
                !RA1=(H1-HZ1)/(LOG(Ru1/rwu1)/(ELEMENT(IEL).PROPERTY(6)*KR1*DIS1)+ELEMENT(IEL).PROPERTY(5))
                RA1=(H1-HZ1)/(val1/(ELEMENT(IEL).PROPERTY(6)*KR1*DIS1)+ELEMENT(IEL).PROPERTY(5))
                !RA1=(H2-H1)/(val1/(ELEMENT(IEL).PROPERTY(6)*KR1*DIS1)+ELEMENT(IEL).PROPERTY(5))
                !W1=ELEMENT(IEL1).ANGLE(SNADJL(IN1).SUBID(J))
                w1=1
                !IF(J<=N1) w1=0.5  !按所辖长度为权，两端点的权比中间的权小一半。 
                Q1(I+2)=Q1(I+2)+RA1*W1
                TW1=TW1+W1;IP1=IP1+1
                !WRITE(99,20) ISTEP,IITER,IEL,I+2,IP1,X2,X1,R1,H1,HZ1,RA1,W1                    

            ENDDO
                
           
            
            Q1(I+2)=Q1(I+2)/TW1
            !A1(I+1)=ABS(Q1(I+2)/DH1)
            QWAN(1,ELEMENT(IEL).NODE(I+2))=QWAN(1,ELEMENT(IEL).NODE(I+2))+Q1(I+2)
            
            IF(ABS(Q1(I+2))<1E-5) THEN
                Q1(I+2)=1.D-5
            ENDIF            
            

            IF(ABS(DH1)<1.D-5) DH1=1.D-5
            !!
            !IF(abs(Qn1(i+2))<abs(Q1(I+2))) then                
            !    print *, Qn1(i+2)
            !ENDIF
            
            
            !ELEMENT(IEL).PROPERTY(I+1)=ABS(2*DH1/(Q1(I+2)+QN1(I+2))) !!!包含了井损。
            ELEMENT(IEL).PROPERTY(I+1)=ABS(DH1/Q1(I+2))!!!包含了井损。
            
          
            A1(I+1)=1.0/(ELEMENT(IEL).PROPERTY(I+1))
        ENDDO
        
        IF(ANY(ISNAN(A1(1:3)))) THEN
            PRINT *, 'NAN'
            STOP "ERRORS IN WELLBORE_Q_K_UPDATE_SPMETHOD."
        ENDIF
        
        KM1(1:NDOF1,1:NDOF1)=KM_WELLBORE(A1(1),A1(2),A1(3))
        KM1(1:NDOF1,1:NDOF1)=(ELEMENT(IEL).KM+KM1(1:NDOF1,1:NDOF1))*0.5
    ENDIF
    
    
    IF(.NOT.ALLOCATED(ELEMENT(IEL).FLUX)) ALLOCATE(ELEMENT(IEL).FLUX(ELEMENT(IEL).NNUM))
    
    ELEMENT(IEL).FLUX=MATMUL(KM1(1:ELEMENT(IEL).NDOF,1:ELEMENT(IEL).NDOF),NHEAD1(1:ELEMENT(IEL).NDOF))
    
    ELEMENT(IEL).KM=KM1(1:ELEMENT(IEL).NDOF,1:ELEMENT(IEL).NDOF)
    !if(solver_control.bfgm==continuum.or.solver_control.bfgm==consistent) then
    !          							
    !end if

10  FORMAT("ISTEP,IITER,IEL,NODE1,P1,XG,YG,ZG,XL,YL,ZL,RA,HA,HW,Qunit,WGT   //WELLBORE INFO")
20  FORMAT(5I7,11F15.7)
    
    contains

    real(8) function fun_lnkr(sita)

    real(8),intent(in)::sita
    !real(8)::kx,ky

    !kx=fparas(1);ky=fparas(2)

    fun_lnkr=log(kx1*ky1/(kx1*(sin(sita))**2+ky1*(cos(sita))**2))
    
    !IF(ISNAN(fun_lnkr)) THEN
    !    PRINT *,'NAN'
    !ENDIF
    !fun_lnkr=exp(sita)


    endfunction
    
END SUBROUTINE

SUBROUTINE WELLBORE_Q_K_UPDATE3(STEPDIS,IEL,ISTEP,IITER)
    USE solverds
    USE MESHGEO,ONLY:SNADJL
    IMPLICIT NONE
    INTEGER,INTENT(IN)::IEL,ISTEP,IITER
    REAL(8),INTENT(IN)::STEPDIS(NDOF)
    INTEGER::I,J,K,I0,IN1,IN2,IN3,N1,N2,IEL1,IP1,MODEL1,II1,NDOF1,ET1,ANODE1
    REAL(8)::H1,X1(3),R1,Q1(4),PI1,LAMDA1,K1,DIS1,DH1,A1(5),NHEAD1(5),KM1(5,5),HZ1,RA1,X2(3),&
    W1,TW1,RE1,VA1,QR1,FD1,L1,D1,QA1,VR1,AREA1,g1,FACC1,REW1,cf1,QN1(5),KX1,KY1,KZ1,SITA1(3),LP1,&
    KR1,XC1(3),VW1(3),T1,C1,B1,PO1,VS1,DIS2,CF2
    LOGICAL::ISO1=.FALSE.,ISIN1=.FALSE.
    
    NDOF1=ELEMENT(IEL).NDOF;ET1=ELEMENT(IEL).ET
    NHEAD1(1:NDOF1)=STEPDIS(ELEMENT(IEL).G)
    PI1=PI(); 
    
    X2=NODE(ELEMENT(IEL).NODE(1)).COORD-NODE(ELEMENT(IEL).NODE(2)).COORD
    DIS1=NORM2(X2)/2.0 !half length of the wellbore
    VW1=X2/DIS1/2.0 !wellbore unit vection    
    
    FD1=ELEMENT(IEL).PROPERTY(1)
    FACC1=0.0D0;REW1=0.D0
        
    MODEL1=INT(MATERIAL(ELEMENT(IEL).MAT).PROPERTY(6))
 
    
    IF(ET1/=WELLBORE_SPGFACE.AND.model1/=4) THEN
        L1=2*DIS1
        D1=2*MATERIAL(ELEMENT(IEL).MAT).PROPERTY(1)
        if(isporeflow>0) D1=element(iel).property(2)
        AREA1=PI1*D1**2/4.0
        QA1=ABS((NHEAD1(1)-NHEAD1(2))*ELEMENT(IEL).KM(1,1))
        VA1=QA1/AREA1
        !if(ABS(MATERIAL(ELEMENT(IEL).MAT).PROPERTY(5))<1.E-7) THEN
        !    g1=73156608000.00 
        !else
        !    g1=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(5)
        !endif
        g1=solver_control.get_g()
        RE1=Re_W(VA1,D1,MATERIAL(ELEMENT(IEL).MAT).PROPERTY(2)) 
        if(isporeflow>0) ELEMENT(IEL).PROPERTY(3)=RE1
        !ACCELERATION LOSS.
        
        IF(ELEMENT(IEL).ET==WELLBORE) THEN
            QR1=(NHEAD1(3)-NHEAD1(2))*ELEMENT(IEL).KM(3,3)+(NHEAD1(4)-NHEAD1(1))*ELEMENT(IEL).KM(4,4)
            FACC1=QR1/(G1*AREA1**2)
            IF(MODEL1==1.AND.ABS(QR1)>1.E-7) THEN
            ! Siwoń, Z. (1987), Solutions for Lateral Inflow in Perforated Conduits, Journal of Hydraulic Engineering, 113(9), 1117-1132.
                PO1=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(9)
                B1=10/(1000*PO1)**4.2+4/10**7                
                C1=1.05+1./(1.235+B1*(QA1/QR1*4*L1*PO1/D1)**2)
            ELSE
                C1=2.0
            ENDIF
            FACC1=C1*FACC1
        ELSE
            QR1=0.D0
        ENDIF
        
        IF(MODEL1==3) THEN
            FD1=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(8)
        ELSE
            !IF(ABS(MATERIAL(ELEMENT(IEL).MAT).PROPERTY(4))<1.D-7) THEN
            !    CF1=1./24/3600 !to m/s
            !ELSE
            !    CF1=1.0D0
            !ENDIF
            !
            
            VR1=QR1/(L1*PI1*D1)
            REW1=Re_W(VR1,D1,MATERIAL(ELEMENT(IEL).MAT).PROPERTY(2))
            FD1=FD_PF(RE1,MATERIAL(ELEMENT(IEL).MAT).PROPERTY(3),REW=REW1,POROSITY=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(9),MODEL=MODEL1) 
        ENDIF
        
        FD1=FD1*L1*QA1/(2*D1*G1*AREA1**2)
        IF(RE1<2260.AND.MODEL1==0)  FD1=ELEMENT(IEL).FD !PIPE-FLOW
        
        !FRICTION PART, CONSIDERING THE TURBULENT AND POROUS SIDE EFFECT
        !IF(MODEL1==2) THEN
        !    !ASSUME LENGHTH UNIT IS M
        !    IF(MATERIAL(ELEMENT(IEL).MAT).PROPERTY(8)>0.D0) THEN
        !        FD1=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(8)
        !    ELSE
        !        !IF(ABS(MATERIAL(ELEMENT(IEL).MAT).PROPERTY(4))<1.D-7) THEN
        !        !    CF1=1./24/3600 !to m/s
        !        !ELSE
        !        !    CF1=1.0D0
        !        !ENDIF
        !        !RE1=Re_W(VA1*CF1,D1,MATERIAL(ELEMENT(IEL).MAT).PROPERTY(2))
        !        VR1=QR1/(L1*PI1*D1)
        !        REW1=Re_W(VR1*CF1,D1,MATERIAL(ELEMENT(IEL).MAT).PROPERTY(2))
        !        FD1=FD_PF(RE1,MATERIAL(ELEMENT(IEL).MAT).PROPERTY(3),REW1) 
        !    ENDIF
        !    FD1=FD1*L1*QA1/(2*D1*G1*AREA1**2)
        !ENDIF            
    END IF 
    
    ELEMENT(IEL).PROPERTY(1)=FD1;ELEMENT(IEL).PROPERTY(4)=FACC1; 
    
    IF(ABS(FD1)<1E-15) THEN
        FD1=1.D-15
    ENDIF
    
    
    A1(1)=1./(FACC1+FD1)
    IF(A1(1)<0.D0) A1(1)=1./FD1    

    
    IF(ET1==PIPE2) THEN
        KM1(1,1)=A1(1);KM1(2,2)=A1(1);KM1(1,2)=-A1(1);KM1(2,1)=-A1(1)
    ELSE
        !WRITE(99,10) 
        
        Q1=0.D0
        
        QN1(1:NDOF1)=MATMUL(ELEMENT(IEL).KM,NHEAD1(1:NDOF1))
        QWAN(2,ELEMENT(IEL).NODE(1:NDOF1))=QWAN(2,ELEMENT(IEL).NODE(1:NDOF1))+QN1(1:NDOF1)
        
        DO I=1,2
    
        
            IN1=ELEMENT(IEL).NODE(2+I)
            !IN2=ELEMENT(IEL).NODE(2+MOD(I,2)+1)
            !DIS1=NORM2(NODE(IN1).COORD-NODE(IN2).COORD)/2.
            IF(I==1) THEN
                N2=2 !2-3
                IF(ELEMENT(IEL).ISTOPO>0) IN1=ELEMENT(IEL).NODE(6)
            ELSE
                N2=1 !1-4
                IF(ELEMENT(IEL).ISTOPO>0) IN1=ELEMENT(IEL).NODE(5)
            ENDIF
            !IF(ET1==WELLBORE) THEN

            DH1=STEPDIS(NODE(IN1).DOF(4))-STEPDIS(NODE(ELEMENT(IEL).NODE(N2)).DOF(4))
            !ELSE
            !    !WELLBORE_SPGFACE
            !    
            !    DH1=STEPDIS(NODE(IN1).DOF(4))-NODE(ELEMENT(IEL).NODE(N2)).COORD(NDIMENSION)
            !ENDIF
            

            IP1=0;RA1=0.D0;W1=0;TW1=0;     
               
            DO J=1,SNADJL(IN1).ENUM
                IEL1=SNADJL(IN1).ELEMENT(J)
                IF(ELEMENT(IEL1).ET==WELLBORE.OR.ELEMENT(IEL1).ET==PIPE2 &
                .OR.ELEMENT(IEL1).ET==WELLBORE_SPGFACE) CYCLE !ITSELF EXCLUDED 
                IF(ELEMENT(IEL1).ET==SPHFLOW.OR.ELEMENT(IEL1).ET==SEMI_SPHFLOW) CYCLE !ITSELF EXCLUDED
                IF(ELEMENT(IEL1).EC/=SPG) CYCLE
                
                    
                DO K=1,ELEMENT(IEL1).NNUM
                    ANODE1=ELEMENT(IEL1).NODE(K)
                    IF(ANY(ELEMENT(IEL).NODE==ANODE1)) CYCLE
                    KR1=0.D0
                    X1=NODE(ANODE1).COORD
                    XC1=X1
                    H1=STEPDIS(NODE(ANODE1).DOF(4))
                    X1=X1-NODE(ELEMENT(IEL).NODE(3)).COORD !ORIGIN IS THE THIRD NODE.
                    IF(NORM2(X1)<MATERIAL(ELEMENT(IEL).MAT).PROPERTY(1)) CYCLE
                    
                    LP1=DOT_PRODUCT(X1,VW1)
                
                    IF(LP1<0.D0.OR.LP1>2*DIS1) CYCLE !跳过不在单元范围内的节点
                
                    X2=LP1*VW1 !X1在井线上的投影 
                    X2=X1-X2 !投影线向量，为此单元的流量的方向矢量
                    R1=NORM2(X2) 
                
                    IF(R1<=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(1)) CYCLE
                    
                    !DIRECTIONAL K
                    SITA1=X2/R1            
                    K1=0.D0
                    DO I0=1,3
                        K1=K1+1.0/MATERIAL(ELEMENT(IEL1).MAT).PROPERTY(I0)*(SITA1(I0))**2
                    ENDDO
                    KR1=1/K1
                    !KR1=(MATERIAL(ELEMENT(IEL1).MAT).PROPERTY(1)*MATERIAL(ELEMENT(IEL1).MAT).PROPERTY(2)*MATERIAL(ELEMENT(IEL1).MAT).PROPERTY(3))**(1/3)
                    IF(H1>XC1(NDIMENSION)) THEN
                        LAMDA1=1.0
                    ELSE
                        LAMDA1=1.D-3
                    ENDIF
                    KR1=LAMDA1*KR1

                    
               
                    HZ1=NHEAD1(2)+(NHEAD1(1)-NHEAD1(2))/(2*DIS1)*LP1
                    !要考虑井损
                    !DIS2=DIS1/(MATERIAL(ELEMENT(IEL1).MAT).PROPERTY(3))**0.5
                    RA1=(H1-HZ1)/(LOG(R1/MATERIAL(ELEMENT(IEL).MAT).PROPERTY(1))/(ELEMENT(IEL).PROPERTY(6)*KR1*DIS1)+ELEMENT(IEL).PROPERTY(5))
                    !W1=ELEMENT(IEL1).ANGLE(SNADJL(IN1).SUBID(J))
                    w1=1
                    Q1(I+2)=Q1(I+2)+RA1*W1
                    TW1=TW1+W1;IP1=IP1+1
                    !WRITE(99,20) ISTEP,IITER,IEL,I+2,IP1,X2,X1,R1,H1,HZ1,RA1,W1                    
                        
                ENDDO

            ENDDO
                
            IF(IP1<1) THEN
                PRINT *, "井单元I周边与其第j节点相连节点的投影点均不在单元内.I=,j=",iel,2+i
                STOP
            ENDIF
            
            
            Q1(I+2)=Q1(I+2)/TW1
            !A1(I+1)=ABS(Q1(I+2)/DH1)
            
            QWAN(1,ELEMENT(IEL).NODE(I+2))=QWAN(1,ELEMENT(IEL).NODE(I+2))+Q1(I+2)
            
            IF(ABS(Q1(I+2))<1E-5) THEN
                Q1(I+2)=1.D-5
            ENDIF
            
            IF(ABS(DH1)<1.D-5) DH1=1.D-5
            
            !ELEMENT(IEL).PROPERTY(I+1)=ABS(2*DH1/(Q1(I+2)+QN1(I+2))) !!!包含了井损。
            ELEMENT(IEL).PROPERTY(I+1)=ABS(DH1/Q1(I+2))!!!包含了井损。
            
          
            A1(I+1)=1.0/(ELEMENT(IEL).PROPERTY(I+1))
        ENDDO
        
        IF(ANY(ISNAN(A1(1:3)))) THEN
            PRINT *, 'NAN'
            STOP "ERRORS IN WELLBORE_Q_K_UPDATE"
        ENDIF
        
         KM1(1:NDOF1,1:NDOF1)=KM_WELLBORE(A1(1),A1(2),A1(3))
         KM1(1:NDOF1,1:NDOF1)=(ELEMENT(IEL).KM+KM1(1:NDOF1,1:NDOF1))*0.5
    ENDIF
    
    
    IF(.NOT.ALLOCATED(ELEMENT(IEL).FLUX)) ALLOCATE(ELEMENT(IEL).FLUX(ELEMENT(IEL).NNUM))
    
    ELEMENT(IEL).FLUX=MATMUL(KM1(1:ELEMENT(IEL).NDOF,1:ELEMENT(IEL).NDOF),NHEAD1(1:ELEMENT(IEL).NDOF))
    
    ELEMENT(IEL).KM=KM1(1:ELEMENT(IEL).NDOF,1:ELEMENT(IEL).NDOF)
    !if(solver_control.bfgm==continuum.or.solver_control.bfgm==consistent) then
    !          							
    !end if

10  FORMAT("ISTEP,IITER,IEL,NODE1,P1,XG,YG,ZG,XL,YL,ZL,RA,HA,HW,Qunit,WGT   //WELLBORE INFO")
20  FORMAT(5I7,11F15.7)
    
END SUBROUTINE

SUBROUTINE Model_MESHTOPO_INI()
    USE MESHADJ,ONLY:SETUP_EDGE_ADJL,SEDGE,NSEDGE,SNADJL,GETGMSHET,ELTTYPE, &
                        SETUP_FACE_ADJL,SFACE,NSFACE, &
                        SETUP_ADJACENT_ELEMENT_SOLVER,SETUP_SUBZONE_SOLVER
    USE solverds,ONLY:isIniSEdge,NODE,ELEMENT,NDIMENSION,NQWNODE,QWELLNODE,WELLBORE,WELLBORE_SPGFACE,PIPE2,POREFLOW,ENUM,QWAN,NNUM
    IMPLICIT NONE
    INTEGER I,J,K,N1,N2,N3,N4,A1(2),A2(200)
    INTEGER,ALLOCATABLE::IELT1(:),ielt2(:)
    REAL(8)::XYLMT(2,3)
    
    
    if(.not.isIniSEdge) then
        CALL SETUP_EDGE_ADJL(SEDGE,NSEDGE,SNADJL)
        CALL SETUP_FACE_ADJL(SFACE,NSFACE,SEDGE,NSEDGE)
        CALL SETUP_ADJACENT_ELEMENT_SOLVER(SEDGE,SFACE,ELEMENT,NDIMENSION)
        !IF(NDIMENSION==2) CALL SETUP_EDGE_BC_SOLVER()
        DO I=1,3
            XYLMT(1,I)=MAXVAL(NODE.COORD(I))
            XYLMT(2,I)=MINVAL(NODE.COORD(I))
        ENDDO
        CALL SETUP_SUBZONE_SOLVER(XYLMT(1,1),XYLMT(2,1),XYLMT(1,2),XYLMT(2,2),XYLMT(1,3),XYLMT(2,3),ELEMENT)
        A1=[4,3]
        
        IF(.NOT.ALLOCATED(IELT1)) ALLOCATE(IELT1(ENUM))
        IELT1=0
        ielt2=ielt1
        DO I=1,NQWNODE
            N3=0
            DO K=1,SNADJL(QWELLNODE(I).NODE(1)).ENUM !如果井口节点位于中间，上下各搜索一次
                N1=QWELLNODE(I).NODE(1)
10              DO J=1,SNADJL(N1).ENUM
                    N2=SNADJL(N1).ELEMENT(J)                    
                    IF(IELT1(N2)>0) CYCLE !假定[WELLBORE,WELLBORE_SPGFACE,PIPE2]单元类中的每个单元只属于唯一井线
                    !与井口相连的必然是这三种单元,且假定井线不相交
                    IF(ANY([WELLBORE,WELLBORE_SPGFACE,PIPE2]-ELEMENT(N2).ET==0)) THEN
                        IELT1(N2)=I
                        IF(ELEMENT(N2).ET/=PIPE2) THEN
                            IF(N3<1) THEN
                                N3=N3+1
                                A2(N3)=ELEMENT(N2).NODE(A1(SNADJL(N1).SUBID(J)))
                                WHERE(ELEMENT(SNADJL(A2(N3)).ELEMENT).ESHAPE>300) IELT2(SNADJL(A2(N3)).ELEMENT)=I
                            ENDIF
                            N4=MOD(SNADJL(N1).SUBID(J),2)+1
                            N3=N3+1
                            A2(N3)=ELEMENT(N2).NODE(A1(N4))
                            WHERE(ELEMENT(SNADJL(A2(N3)).ELEMENT).ESHAPE>300) IELT2(SNADJL(A2(N3)).ELEMENT)=I
                            N1=ELEMENT(N2).NODE(N4)
                            GOTO 10
                        ENDIF
                    ENDIF  
                ENDDO
            ENDDO
            IF(N3<1) THEN
                PRINT *, "NO ELEMENT CONNECTING TO WELLHEAD I,I=",N1
            ELSE
                QWELLNODE(I).NNODE2=N3
                QWELLNODE(I).NODE2=A2(1:N3)
                QWELLNODE(I).ELEMENT=PACK([1:enum],IELT2==I)
                ALLOCATE(QWELLNODE(I).QAN(2,N3))                
            ENDIF
        ENDDO
        ALLOCATE(QWAN(2,NNUM))        
        ISINISEDGE=.TRUE.
    endif

ENDSUBROUTINE

SUBROUTINE INI_WELLBORE(IELT)
    
    USE MESHADJ,ONLY:SEDGE,NSEDGE,SNADJL,GETGMSHET,ELTTYPE,SFACE,NSFACE,&
                    POINTlOC_BC_SOLVER,getval_solver
    USE forlab,ONLY:ISOUTLIER
    USE solverds
    USE SolverMath
    USE IFPORT
    IMPLICIT NONE
    INTEGER,INTENT(IN)::IELT
    INTEGER::I,J,K,IELT1,IEDGE1,GMET1,SUBID1,AF1(2),N1,NADJL1(200),N2,N3,K1,N4
    REAL(8)::AV1(3,2),AR1(3,3),ZV1(3),YV1(3),XV1(3),D1,T1,T2,RDIS1(200),HK1,PHI1,TPHI1,WR1,L1,&
        ORG1(3),AREA1,XYLMT(2,3),SR1,try1,try2,TRY0
    REAL(8),ALLOCATABLE::RDIS2(:),SPT1(:,:)
    
    if(.not.isIniSEdge) CALL Model_MESHTOPO_INI()
    
    
    
    IEDGE1=ABS(ELEMENT(IELT).EDGE(3))
    L1=SEDGE(IEDGE1).DIS/2.0
    IF(ELEMENT(IELT).ISTOPO>0) IEDGE1=ABS(ELEMENT(IELT).EDGE(5))
    N1=0
    
    IF(SEDGE(IEDGE1).ENUM<2) THEN
        PRINT *, "NO ELEMENTS SURROUNDING THE WELLBORE ELEMENT I. I=",IELT
        PRINT *, "TRY TO EMBED THE WELL LINE INTO THE CORRESPOINDING VOLUME."
        STOP
    ENDIF
    !统计以井流单元为边的单元数
    DO I=1,SEDGE(IEDGE1).ENUM
        IELT1=SEDGE(IEDGE1).ELEMENT(I)
        GMET1=GETGMSHET(ELEMENT(IELT1).ET)
        IF(ELTTYPE(GMET1).DIM==3.AND.ELEMENT(IELT1).EC==SPG) THEN
            N1=N1+1
            IF(.NOT.ALLOCATED(ELEMENT(IELT1).ANGLE)) CALL calangle(IELT1) 
        ENDIF
    ENDDO
    
       
    ALLOCATE(ELEMENT(IELT).NODE2(N1),ELEMENT(IELT).ANGLE(N1))
    
    DO I=1,2
        DO J=1,SNADJL(SEDGE(IEDGE1).V(I)).ENUM
            IELT1=SNADJL(SEDGE(IEDGE1).V(I)).ELEMENT(J)
            IF(IELT1<1)  CYCLE
            IF(ELEMENT(IELT1).ET==WELLBORE.OR.ELEMENT(IELT1).ET==PIPE2 &
            .OR.ELEMENT(IELT1).ET==WELLBORE_SPGFACE) CYCLE !ITSELF EXCLUDED
            IF(ELEMENT(IELT1).ET==SPHFLOW.OR.ELEMENT(IELT1).ET==SEMI_SPHFLOW) CYCLE !ITSELF EXCLUDED
            IF(ELEMENT(IELT1).EC/=SPG) CYCLE
            IF(ALLOCATED(ELEMENT(IELT1).ANGLE)) CYCLE
            CALL calangle(IELT1)
        ENDDO    
    ENDDO
    
    N1=0;HK1=0.D0
    DO I=1,SEDGE(IEDGE1).ENUM
        IELT1=SEDGE(IEDGE1).ELEMENT(I)
        GMET1=GETGMSHET(ELEMENT(IELT1).ET)
        IF(ELTTYPE(GMET1).DIM<3) CYCLE
        SUBID1=SEDGE(IEDGE1).SUBID(I)
        N1=N1+1
        ELEMENT(IELT).NODE2(N1)=IELT1
        AF1=PACK([1:ELTTYPE(GMET1).NFACE],MASK=ANY(ABS(ELTTYPE(GMET1).FACEEDGE(1:4,:))-SUBID1==0,DIM=1))
        DO J=1,2
            AR1(:,1)=NODE(ELEMENT(IELT1).NODE(ELTTYPE(GMET1).FACE(1,AF1(J)))).COORD
            AR1(:,2)=NODE(ELEMENT(IELT1).NODE(ELTTYPE(GMET1).FACE(2,AF1(J)))).COORD
            AR1(:,3)=NODE(ELEMENT(IELT1).NODE(ELTTYPE(GMET1).FACE(3,AF1(J)))).COORD 
 
            AV1(:,J)=NORMAL_TRIFACE(AR1)            
        ENDDO
        ELEMENT(IELT).ANGLE(N1)=DihedralAngle(AV1(:,1),AV1(:,2))
        !井周边单元的的平均渗透系数，
        HK1=HK1+ELEMENT(IELT).ANGLE(N1)*(MATERIAL(ELEMENT(IELT1).MAT).PROPERTY(1)*MATERIAL(ELEMENT(IELT1).MAT).PROPERTY(2))**0.5        
    ENDDO
    
    TPHI1=SUM(ELEMENT(IELT).ANGLE)
    HK1=HK1/TPHI1
    ELEMENT(IELT).PROPERTY(6)=TPHI1
    IF(ABS(TPHI1)<1E-7) THEN
        PRINT *, "NO ELEMENTS ADJACENT TO THE WELLBORE ELEMENT I. IELT=",IELT
        STOP
    ENDIF
    PHI1=TPHI1/SIZE(ELEMENT(IELT).ANGLE)
    
    !LOCAL COORDINATE SYSTEM
    !以wellbore单元，3节点为原点，3-4为z轴的局部坐标系。
    ALLOCATE(ELEMENT(IELT).G2L(3,3))
    ORG1=NODE(ELEMENT(IELT).NODE(3)).COORD
    !local z
    ZV1=NODE(ELEMENT(IELT).NODE(4)).COORD-ORG1
    ZV1=ZV1/NORM2(ZV1)
    !local x
    !过ogr1与zv1垂直的平面方程：Ax+By+Cz+D=0,先求出D,然后假定（x,y,z）任意两个，求第三个。
    D1=-DOT_PRODUCT(ZV1,ORG1) !
    !
    IF(ABS(ZV1(3))>1E-7) THEN
        XV1=[ORG1(1)+1.0,0.D0,-(D1+ZV1(1)*(ORG1(1)+1.0))/ZV1(3)]
    ELSEIF(ABS(ZV1(2))>1E-7) THEN
        XV1=[ORG1(1)+1.0,-(D1+ZV1(1)*(ORG1(1)+1.0))/ZV1(2),0.D0]
    ELSE
        XV1=[-(D1+ZV1(2)*(ORG1(2)+1.0))/ZV1(1),ORG1(2)+1.0,0.D0]
    ENDIF
    XV1=XV1-ORG1
    XV1=XV1/NORM2(XV1)
    !local  y
    YV1=NORMAL_TRIFACE(RESHAPE([0.D0,0.D0,0.D0,ZV1,XV1],([3,3])))       
    YV1=YV1/NORM2(YV1)
    ELEMENT(IELT).G2L(1,:)=XV1;ELEMENT(IELT).G2L(2,:)=YV1;ELEMENT(IELT).G2L(3,:)=ZV1;
    
    WR1=MATERIAL(ELEMENT(IELT).MAT).PROPERTY(1)
    T2=PHI1*SIN(PHI1)/(1-COS(PHI1))
    SR1=0;N2=0
    DO J=1,2        
        N3=SEDGE(IEDGE1).V(J)
        DO I=1,SNADJL(N3).NNUM
            N1=SNADJL(N3).NODE(I)
            XV1=NODE(N1).COORD-NODE(ELEMENT(IELT).NODE(3)).COORD
            XV1=MATMUL(ELEMENT(IELT).G2L(:,:),XV1)
            T1=NORM2(XV1(1:2))
            IF(T1>WR1) THEN
                N2=N2+1
                !NADJL1(N2)=N1;
                RDIS1(N2)=T1
            ENDIF
           
        ENDDO
        !T1=SUM(RDIS1(1:N2))/N2
        
        !生成采样点半径取相邻最大半径，里面受井节点的影响，可能不准。
        !SR1=MAXVAL([SR1,RDIS1(1:N2)])
        !IF(.NOT.ALLOCATED(RDIS2)) THEN
        !    RDIS2=RDIS1(1:N2)
        !ELSE
        !    RDIS2=[RDIS2,RDIS1(1:N2)]
        !ENDIF
    ENDDO
    
    SR1=MAXVAL(RDIS1(1:N2),MASK=isoutlier(RDIS1(1:N2))==.FALSE.) !取异常值（大于三倍标准差）之外的最大值
    ELEMENT(IELT).PROPERTY(2:3)=1./(TPHI1*HK1*L1/(LOG(SR1/WR1)-T2))
    !ELEMENT(IELT).PROPERTY(2:3)=1.d-7
        !IF(J==3) THEN
        !    !ALLOCATE(ELEMENT(IELT).WELL_SP3(N2),ELEMENT(IELT).WELL_SP3_R(N2))
        !    !ELEMENT(IELT).NWSP3=N2
        !    !ELEMENT(IELT).WELL_SP3=NADJL1(1:N2)
        !    !ELEMENT(IELT).WELL_SP3_R=RDIS1(1:N2)
        !    ELEMENT(IELT).PROPERTY(2)=1./(TPHI1*HK1*L1/(LOG(T1/WR1)-T2))
        !ELSE
        !    !ALLOCATE(ELEMENT(IELT).WELL_SP4(N2),ELEMENT(IELT).WELL_SP4_R(N2))
        !    !ELEMENT(IELT).NWSP4=N2
        !    !ELEMENT(IELT).WELL_SP4=NADJL1(1:N2)
        !    !ELEMENT(IELT).WELL_SP4_R=RDIS1(1:N2)
        !    ELEMENT(IELT).PROPERTY(3)=1./(TPHI1*HK1*L1/(LOG(T1/WR1)-T2))
        !ENDIF
        
   
    
    
    
    !生成采样点
    
    IF(solver_control.WELLMETHOD>1) THEN
        N1=solver_control.nspwell
        T1=2*PI()/N1
        ALLOCATE(SPT1(3,N1))
        DO J=1,N1
            SPT1(1,J)=SR1*COS(T1*(J-1));SPT1(2,J)=SR1*SIN(T1*(J-1));SPT1(3,J)=0;
            SPT1(:,J)=MATMUL(TRANSPOSE(ELEMENT(IELT).G2L(:,:)),SPT1(:,J))
        ENDDO
        !SPT1(:,N1+1:2*N1)=(1+0.1/SR1)*SPT1(:,1:N1)
        
        
    
        ALLOCATE(ELEMENT(IELT).A12(3,3*N1),ELEMENT(IELT).NSPLOC(3*N1)) 
        ELEMENT(IELT).NSPLOC=0
        N2=0;N4=0
        DO I=1,1
            
            DO J=3,5
                IF(J<5) THEN
                    ZV1=NODE(ELEMENT(IELT).NODE(J)).COORD
                ELSE
                    ZV1=(NODE(ELEMENT(IELT).NODE(3)).COORD+NODE(ELEMENT(IELT).NODE(4)).COORD)/2
                ENDIF
                IF(N2<1) N2=ELEMENT(IELT).NODE2(MAX(INT(rand(1)*SIZE(ELEMENT(IELT).NODE2)),1))
                
                DO K=(I-1)*N1+1,I*N1
                    !!因为sr1取最大，所以有可能一些点不在区域内,采样二分法找到最大的区域内的点。
                    try1=1.;try2=WR1/NORM2(SPT1(:,K));TRY0=1;
                    N3=0
                    DO K1=10,1,-1
                        YV1=try0*SPT1(:,K)+ZV1                    
                        N2=POINTlOC_BC_SOLVER(YV1,N2)
                        IF(N2==0) THEN
                            IF(TRY0<TRY1) TRY1=TRY0
                        ELSE
                            IF(TRY0>TRY2) THEN
                                TRY2=TRY0
                                N3=N2
                            ENDIF
                        ENDIF
                        IF(ABS(TRY1-TRY2)<0.1) THEN
                            N2=N3
                            EXIT
                        ENDIF
                    
                        IF(K1==10) THEN
                            !检查最里面的点是否在区域内
                            YV1=try2*SPT1(:,K)+ZV1                    
                            N2=POINTlOC_BC_SOLVER(YV1,N2)
                            IF(N2==0) THEN
                                EXIT
                            ELSE
                                N3=N2
                            ENDIF
                        ENDIF
                    
                        try0=(try1+try2)/2                    
                    ENDDO
                    N4=N4+1
                    IF(N2==0) THEN
                    
                        PRINT *, 'FAILED TO LOCATE SAMPLE POINT (J,K),IN SUB INI_WELLBORE. (IELT,J,K)=',IELT,J,K
                        !PRINT *,'建议控制井周单元尺寸>5倍井半径.'
                        !STOP
                    ELSE
                        
                        ELEMENT(IELT).A12(:,N4)=YV1
                        ELEMENT(IELT).NSPLOC(N4)=N2
 
                        !CALL getval_solver(ELEMENT(IELT).A12(:,N1*(J-3)+K),N2,YV1,reshape([node(element(n2).NODE).COORD(1),&
                        !node(element(n2).NODE).COORD(2),node(element(n2).NODE).COORD(3)],[element(n2).NNUM,3]))
                        !T1=NORM2(YV1-ELEMENT(IELT).A12(:,N1*(J-3)+K))
                        !IF(T1>1E-6) THEN
                        !    PRINT *,"FAILED IN GETVAL.ERROR=",T1
                        !ELSE
                        !    PRINT *,"SUCCESSED IN GETVAL.ERROR=", T1
                        !ENDIF               
                
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
    ENDIF
    
    
    
    !井损
    IF(ABS(MATERIAL(ELEMENT(IELT).MAT).PROPERTY(7))>1E-12) THEN
        AREA1=2.0*PI()*MATERIAL(ELEMENT(IELT).MAT).PROPERTY(1)*L1
        ELEMENT(IELT).PROPERTY(5)=MATERIAL(ELEMENT(IELT).MAT).PROPERTY(7)/AREA1
        ELEMENT(IELT).PROPERTY(2:3)=ELEMENT(IELT).PROPERTY(2:3)+ELEMENT(IELT).PROPERTY(5)
        !ELEMENT(IELT).PROPERTY(2:3)=ELEMENT(IELT).PROPERTY(2:3)
    ELSE
        ELEMENT(IELT).PROPERTY(5)=0.D0
    ENDIF
    !ALLOCATE(ELEMENT(IELT).KM(4,4))
    
    ELEMENT(IELT).KM=KM_WELLBORE(1./ELEMENT(IELT).PROPERTY(1),1./ELEMENT(IELT).PROPERTY(2),1./ELEMENT(IELT).PROPERTY(3))


    
ENDSUBROUTINE

SUBROUTINE INI_SPHFLOW(IELT)
    
    USE MESHADJ
    USE solverds
    USE SolverMath
    USE MESHGEO,ONLY:getval_solver
    USE forlab ,ONLY:ISOUTLIER
    IMPLICIT NONE
    INTEGER,INTENT(IN)::IELT
    INTEGER::I,J,IELT1,N1,N2,N3,K
    REAL(8)::XV1(3),T1,HK1,PHI1,TPHI1,WR1,L1,DV1(3),PI1,ZV1(3),YV1(3),D1,ANGLE1(7),SPT1(3,100),SR1,ORG1(3),try0,try1,try2
    REAL(8)::RDIS2(400)
    
    if(.not.isIniSEdge) CALL Model_MESHTOPO_INI()
    
        
    N1=ELEMENT(IELT).NODE(2)
    TPHI1=0.D0;HK1=0.D0
    DV1(1:NDIMENSION)=ELEMENT(IELT).PROPERTY(1:NDIMENSION)
    WR1=MATERIAL(ELEMENT(IELT).MAT).PROPERTY(1)
    !DO J=1,SNADJL(N1).ENUM        
    !    IELT1=SNADJL(N1).ELEMENT(J)
    !    IF(ELEMENT(IELT1).ET==WELLBORE &
    !    .OR.ELEMENT(IELT1).ET==PIPE2 &
    !    .OR.ELEMENT(IELT1).ET==SPHFLOW &
    !    .OR.ELEMENT(IELT1).ET==SEMI_SPHFLOW &
    !    .OR.ELEMENT(IELT1).ET==WELLBORE_SPGFACE) CYCLE !ITSELF EXCLUDED
    !    
    !    IF(ELEMENT(IELT1).EC/=SPG) CYCLE
    !    IF(ALLOCATED(ELEMENT(IELT1).ANGLE)) CYCLE
    !    CALL calangle(IELT1)
    !    !AVERAGED K, ASSUME ISOTROPIC        
    !    IF(ELEMENT(IELT).ET==SEMI_SPHFLOW) THEN
    !        !centroid
    !        DO I=1,3
    !            XV1(I)=SUM(NODE(ELEMENT(IELT1).NODE).COORD(I))
    !        ENDDO
    !        
    !        XV1=XV1/ELEMENT(IELT1).NNUM
    !        print *, j,xv1
    !        XV1=XV1-NODE(N1).COORD
    !        T1=DOT_PRODUCT(XV1(1:NDIMENSION),DV1(1:NDIMENSION))/ &
    !        (NORM2(XV1(1:NDIMENSION))*NORM2(DV1(1:NDIMENSION)))
    !        IF(T1<0.D0) CYCLE
    !    ENDIF
    !    PHI1=ELEMENT(IELT1).ANGLE(SNADJL(N1).SUBID(J)) 
    !    TPHI1=TPHI1+PHI1
    !    HK1=HK1+PHI1*(MATERIAL(ELEMENT(IELT1).MAT).PROPERTY(1)*MATERIAL(ELEMENT(IELT1).MAT).PROPERTY(2))**0.5        
    !ENDDO
    !HK1=HK1/TPHI1 !平均的K
    !统计边长
    HK1=1
    N2=0;L1=0.D0;SR1=0.D0    
    DO J=1,SNADJL(N1).NNUM        
        N3=SNADJL(N1).NODE(J)
        XV1=NODE(N3).COORD-NODE(N1).COORD
        T1=NORM2(XV1)
        IF(T1<=WR1) CYCLE
        IF(ELEMENT(IELT).ET==SEMI_SPHFLOW) THEN
            IF(DOT_PRODUCT(XV1(1:NDIMENSION),DV1(1:NDIMENSION))<0.D0) CYCLE
        ENDIF            
        
        !SR1=MAX(T1,SR1)

        N2=N2+1
        RDIS2(N2)=T1
        L1=L1+T1
                   
    ENDDO 
    L1=L1/N2
    
    !生成采样点    
    
    IF(solver_control.WELLMETHOD>2) THEN
        SR1=MAXVAL(RDIS2(1:N2),MASK=.NOT.ISOUTLIER(RDIS2(1:N2))) !取剔除异常值后的最大值
        !LOCAL COORDINATE SYSTEM
        !以wellbore单元，3节点为原点，3-4为z轴的局部坐标系。
        ALLOCATE(ELEMENT(IELT).G2L(3,3))
        ORG1=NODE(ELEMENT(IELT).NODE(2)).COORD
        !local z
        ZV1=DV1
        ZV1=ZV1/NORM2(ZV1)
        !local x
        !过ogr1与zv1垂直的平面方程：Ax+By+Cz+D=0,先求出D,然后假定（x,y,z）任意两个，求第三个。
        D1=-DOT_PRODUCT(ZV1,ORG1) !
        !
        IF(ABS(ZV1(3))>1E-7) THEN
            XV1=[ORG1(1)+1.0,0.D0,-(D1+ZV1(1)*(ORG1(1)+1.0))/ZV1(3)]
        ELSEIF(ABS(ZV1(2))>1E-7) THEN
            XV1=[ORG1(1)+1.0,-(D1+ZV1(1)*(ORG1(1)+1.0))/ZV1(2),0.D0]
        ELSE
            XV1=[-(D1+ZV1(2)*(ORG1(2)+1.0))/ZV1(1),ORG1(2)+1.0,0.D0]
        ENDIF
        XV1=XV1-ORG1
        XV1=XV1/NORM2(XV1)
        !local  y
        YV1=NORMAL_TRIFACE(RESHAPE([0.D0,0.D0,0.D0,ZV1,XV1],([3,3])))       
        YV1=YV1/NORM2(YV1)
        ELEMENT(IELT).G2L(1,:)=XV1;ELEMENT(IELT).G2L(2,:)=YV1;ELEMENT(IELT).G2L(3,:)=ZV1; 
        
        
        N1=0
        angle1=[360,60,36,30,36,60,360]
        DO I=1,7
            T1=SR1*SIND(30*(REAL(I)-1))
            DO J=0,359,angle1(I) !359 防止生成与0重复的点
                N1=N1+1
                SPT1(1,N1)=T1*COSD(REAL(J));SPT1(2,N1)=T1*SIND(REAL(J));SPT1(3,N1)=SR1*COSD(30*(REAL(I)-1));
                SPT1(:,N1)=MATMUL(TRANSPOSE(ELEMENT(IELT).G2L(:,:)),SPT1(:,N1))
                IF(ELEMENT(IELT).ET==SEMI_SPHFLOW) THEN
                    D1=DOT_PRODUCT(SPT1(1:NDIMENSION,N1),DV1(1:NDIMENSION))                    
                    IF(D1<0.D0) THEN
                        N1=N1-1
                        CYCLE
                    ENDIF
                ENDIF  
                SPT1(:,N1)=SPT1(:,N1)+ORG1
            ENDDO
        ENDDO
    
        ALLOCATE(ELEMENT(IELT).A12(3,N1),ELEMENT(IELT).NSPLOC(N1))
        ELEMENT(IELT).NSPLOC=0
        N2=0;N3=ELEMENT(IELT).NODE(2);J=0
        DO I=1,N1
 
            IF(N2<1) THEN
                N2=MINLOC(ABS(ELEMENT(SNADJL(N3).ELEMENT).ESHAPE-304),DIM=1)
                
                IF(N2<1) THEN
                    N2=0
                ELSE
                    N2=SNADJL(N3).ELEMENT(N2)
                ENDIF
            ENDIF
            
            !DO K=10,1,-1                
            !    YV1=(SPT1(:,I)-ORG1)*REAL(K)/10.+ORG1
            !    N2=POINTlOC_BC_SOLVER(YV1,N2)
            !    IF(N2>0.OR.NORM2(YV1-ORG1)<=WR1) EXIT
            !ENDDO
            !!因为sr1取最大，所以有可能一些点不在区域内,采样二分法找到最大的区域内的点。
            try1=1.;try2=WR1/NORM2(SPT1(:,I)-ORG1);TRY0=1;
            N3=0
            DO K=10,1,-1
                YV1=(SPT1(:,I)-ORG1)*try1+ORG1                    
                N2=POINTlOC_BC_SOLVER(YV1,N2)
                IF(N2==0) THEN
                    IF(TRY0<TRY1) TRY1=TRY0
                ELSE
                    IF(TRY0>TRY2) THEN
                        TRY2=TRY0
                        N3=N2
                    ENDIF
                ENDIF
                IF(ABS(TRY1-TRY2)<0.1) THEN
                    N2=N3
                    EXIT
                ENDIF 
                IF(K==10) THEN
                    !检查最里面的点是否在区域内
                    YV1=(SPT1(:,I)-ORG1)*try2+ORG1                     
                    N2=POINTlOC_BC_SOLVER(YV1,N2)
                    IF(N2==0) THEN
                        EXIT
                    ELSE
                        N3=N2
                    ENDIF
                ENDIF                
                try0=(try1+try2)/2                    
            ENDDO
            
            
            
            !YV1=SPT1(:,I)
            !N2=POINTlOC_BC_SOLVER(YV1,N2)
            !因为sr1取最大，所以有可能一些点不在区域内。
            IF(N2>0) THEN
                J=J+1
                ELEMENT(IELT).A12(:,J)=YV1
                ELEMENT(IELT).NSPLOC(J)=N2
                !CALL getval_solver(ELEMENT(IELT).A12(:,I),N2,YV1,reshape([node(element(n2).NODE).COORD(1),&
                !node(element(n2).NODE).COORD(2),node(element(n2).NODE).COORD(3)],[element(n2).NNUM,3]))
                !T1=NORM2(YV1-ELEMENT(IELT).A12(:,I))
                !IF(T1>1E-6) THEN
                !    PRINT *,"FAILED IN GETVAL.ERROR=",T1
                !    
                !ELSE
                !    PRINT *,"SUCCESSED IN GETVAL.", T1
                !ENDIF               
            ENDIF
            
        ENDDO
        
        IF(J<1) THEN
            PRINT *,"NO SAMPLE POINTS AROUND SPHFLOW ELEMENT I. I= ",IELT
            STOP
        ELSE
            ELEMENT(IELT).A12=ELEMENT(IELT).A12(:,:J)
            ELEMENT(IELT).NSPLOC=ELEMENT(IELT).NSPLOC(:J)
        ENDIF
    ENDIF
    
    PI1=PI()
    T1=2*PI1*HK1/(1/WR1-(4*PI1+3)/(3*L1))
    !T1=1.D10
    IF(ELEMENT(IELT).ET==SPHFLOW) T1=2*T1
    ELEMENT(IELT).KM(1,1)=1.D0;ELEMENT(IELT).KM(2,2)=1.D0
    ELEMENT(IELT).KM(2,1)=-1.D0;ELEMENT(IELT).KM(1,2)=-1.D0
    ELEMENT(IELT).KM=ELEMENT(IELT).KM*T1
    ELEMENT(IELT).PROPERTY(1)=1/T1
    
ENDSUBROUTINE

SUBROUTINE SPHFLOW_Q_K_SPMETHOD(STEPDIS,IELT,ISTEP,IITER)
    
    USE solverds
    USE MESHGEO,ONLY:SNADJL,GETVAL_SOLVER
    IMPLICIT NONE
    INTEGER,INTENT(IN)::IELT,ISTEP,IITER
    REAL(8),INTENT(IN)::STEPDIS(NDOF)
    
    INTEGER::I,J,K,IELT1,N1,N2,N3,ANODE1
    REAL(8)::XV1(3),T1,HK1,PHI1,TPHI1,WR1,L1,DV1(3),PI1,RO1,Q1,X1(3),H1,LAMDA1,K1,KM1(2,2),NHEAD1(2),NHEAD2(10,1)
    REAL(8)::SITA1(3),KX1,KY1,KZ1,KR1,SPH1(10),QN1(2),T2
    LOGICAL::ISISOTROPIC=.FALSE.
   
        
    N1=ELEMENT(IELT).NODE(2)
    TPHI1=0.D0;HK1=0.D0
    DV1(1:NDIMENSION)=ELEMENT(IELT).PROPERTY(1:NDIMENSION)
    WR1=MATERIAL(ELEMENT(IELT).MAT).PROPERTY(1)
    PI1=PI();Q1=0.D0
    NHEAD1=STEPDIS(ELEMENT(IELT).G)
    QN1(1:2)=MATMUL(ELEMENT(IELT).KM,NHEAD1(1:2))
    QWAN(2,ELEMENT(IELT).NODE(1:2))=QWAN(2,ELEMENT(IELT).NODE(1:2))+QN1(1:2)
    
    IF(ANY(NHEAD1-NODE(ELEMENT(IELT).NODE).COORD(NDIMENSION)<0.D0)) THEN
        T1=1.D-7
    ELSE        
        
        DO J=1,SIZE(ELEMENT(IELT).NSPLOC)
            IF(ELEMENT(IELT).NSPLOC(J)<1) CYCLE 
            IELT1=ELEMENT(IELT).NSPLOC(J)

            X1=ELEMENT(IELT).A12(:,J)
            N2=ELEMENT(IELT1).NDOF
            NHEAD2(1:N2,1)=STEPDIS(ELEMENT(IELT1).G)
            
            
            CALL GETVAL_SOLVER(X1,IELT1,SPH1(1:1),NHEAD2(1:N2,:))            
 
            H1=SPH1(1)
            
            XV1=X1-NODE(N1).COORD !ORIGIN .
            RO1=NORM2(XV1)
            IF(RO1<WR1) CYCLE
            KR1=0.D0             
 
            SITA1=XV1/RO1            
            K1=0.D0
            DO I=1,3
                K1=K1+1.0/MATERIAL(ELEMENT(IELT1).MAT).PROPERTY(I)*(SITA1(I))**2
            ENDDO
            KR1=1/K1;
                
            IF(H1>X1(NDIMENSION)) THEN
                LAMDA1=1.0
            ELSE
                LAMDA1=1.D-3
            ENDIF
            KR1=LAMDA1*KR1

            IF(ELEMENT(IELT).ET==SEMI_SPHFLOW) THEN
                T1=DOT_PRODUCT(XV1(1:NDIMENSION),DV1(1:NDIMENSION))                
                IF(T1<0.D0) CYCLE
            ENDIF
            !H1=SUM(STEPDIS(ELEMENT(IELT1).G))/ELEMENT(IELT1).NNUM        
            !PHI1=ELEMENT(IELT1).ANGLE(SNADJL(N1).SUBID(J)) 
            PHI1=1.D0
            TPHI1=TPHI1+PHI1
            IF(SOLVER_CONTROL.well_bottom_type==0) THEN
                T2=4.0
            ELSE
                T2=2.0*PI1
            ENDIF

            Q1=Q1+(H1-NHEAD1(1))*T2*KR1/(1/WR1-1/RO1)*PHI1    
                        
                   

       
        ENDDO
        Q1=Q1/TPHI1
        IF(ELEMENT(IELT).ET==SPHFLOW) Q1=2.D0*Q1
        QWAN(1,ELEMENT(IELT).NODE(2))=QWAN(1,ELEMENT(IELT).NODE(2))+Q1
        
        IF(ABS(Q1)<1.D-5) Q1=1.D-5
        T1=(NHEAD1(2)-NHEAD1(1))
        IF(ABS(T1)<1.D-5) T1=1.D-5
        T1=ABS(Q1/T1) !!!!abs
    ENDIF
    
    KM1(1,1)=1.D0;KM1(2,2)=1.D0
    KM1(2,1)=-1.D0;KM1(1,2)=-1.D0
    T1=(T1+1.0/element(ielt).property(1))/2.0 !!!!
    KM1=KM1*T1
    element(ielt).property(1)=1/T1

    IF(.NOT.ALLOCATED(ELEMENT(IELT).FLUX)) ALLOCATE(ELEMENT(IELT).FLUX(ELEMENT(IELT).NNUM))
    
    ELEMENT(IELT).FLUX=MATMUL(KM1,NHEAD1)
    
    !if(solver_control.bfgm==continuum.or.solver_control.bfgm==consistent) then
    ELEMENT(IELT).KM=KM1    							
    !end if
    
ENDSUBROUTINE


SUBROUTINE SPHFLOW_Q_K_UPDATE3(STEPDIS,IELT,ISTEP,IITER)
    
    USE solverds
    USE MESHGEO,ONLY:SNADJL
    IMPLICIT NONE
    INTEGER,INTENT(IN)::IELT,ISTEP,IITER
    REAL(8),INTENT(IN)::STEPDIS(NDOF)
    
    INTEGER::I,J,K,IELT1,N1,N2,N3,ANODE1
    REAL(8)::XV1(3),T1,HK1,PHI1,TPHI1,WR1,L1,DV1(3),PI1,RO1,Q1,X1(3),H1,LAMDA1,K1,KM1(2,2),NHEAD1(2)
    REAL(8)::SITA1(3),KX1,KY1,KZ1,KR1,QN1(2),T2
    LOGICAL::ISISOTROPIC=.FALSE.
   
        
    N1=ELEMENT(IELT).NODE(2)
    TPHI1=0.D0;HK1=0.D0
    DV1(1:NDIMENSION)=ELEMENT(IELT).PROPERTY(1:NDIMENSION)
    WR1=MATERIAL(ELEMENT(IELT).MAT).PROPERTY(1)
    PI1=PI();Q1=0.D0
    NHEAD1=STEPDIS(ELEMENT(IELT).G)
    QN1(1:2)=MATMUL(ELEMENT(IELT).KM,NHEAD1(1:2))
    QWAN(2,ELEMENT(IELT).NODE(1:2))=QWAN(2,ELEMENT(IELT).NODE(1:2))+QN1(1:2)    
    
    
    IF(ANY(NHEAD1-NODE(ELEMENT(IELT).NODE).COORD(NDIMENSION)<0.D0)) THEN
        T1=1.D-7
    ELSE        
    
        DO J=1,SNADJL(N1).ENUM        
            IELT1=SNADJL(N1).ELEMENT(J)
            IF(ELEMENT(IELT1).ET==WELLBORE &
            .OR.ELEMENT(IELT1).ET==PIPE2 &
            .OR.ELEMENT(IELT1).ET==SPHFLOW &
            .OR.ELEMENT(IELT1).ET==SEMI_SPHFLOW &
            .OR.ELEMENT(IELT1).ET==WELLBORE_SPGFACE) CYCLE !ITSELF EXCLUDED
        
            IF(ELEMENT(IELT1).EC/=SPG) CYCLE
            
            DO K=1,ELEMENT(IELT1).NNUM
                ANODE1=ELEMENT(IELT1).NODE(K)
                IF(ANY(ELEMENT(IELT).NODE==ANODE1)) CYCLE

                X1=NODE(ANODE1).COORD                
                H1=STEPDIS(NODE(ANODE1).DOF(4))
                XV1=X1-NODE(N1).COORD !ORIGIN .
                RO1=NORM2(XV1)
                IF(RO1<WR1) CYCLE
                KR1=0.D0             
 
                SITA1=XV1/RO1            
                K1=0.D0
                DO I=1,3
                    K1=K1+1.0/MATERIAL(ELEMENT(IELT1).MAT).PROPERTY(I)*(SITA1(I))**2
                ENDDO
                KR1=1/K1;
                
                IF(H1>X1(NDIMENSION)) THEN
                    LAMDA1=1.0
                ELSE
                    LAMDA1=1.D-3
                ENDIF
                KR1=LAMDA1*KR1

                IF(ELEMENT(IELT).ET==SEMI_SPHFLOW) THEN
                    T1=DOT_PRODUCT(XV1(1:NDIMENSION),DV1(1:NDIMENSION))                
                    IF(T1<0.D0) CYCLE
                ENDIF
                !H1=SUM(STEPDIS(ELEMENT(IELT1).G))/ELEMENT(IELT1).NNUM        
                !PHI1=ELEMENT(IELT1).ANGLE(SNADJL(N1).SUBID(J)) 
                PHI1=1.D0
                TPHI1=TPHI1+PHI1
                IF(SOLVER_CONTROL.well_bottom_type==0) THEN
                    T2=4.0
                ELSE
                    T2=2.0*PI1
                ENDIF
                Q1=Q1+(H1-NHEAD1(1))*T2*KR1/(1/WR1-1/RO1)*PHI1    
                        
            ENDDO            

       
        ENDDO
        Q1=Q1/TPHI1
        IF(ELEMENT(IELT).ET==SPHFLOW) Q1=2.D0*Q1
        QWAN(1,ELEMENT(IELT).NODE(2))=QWAN(1,ELEMENT(IELT).NODE(2))+Q1
        IF(ABS(Q1)<1.D-5) Q1=1.D-5
        T1=(NHEAD1(2)-NHEAD1(1))
        IF(ABS(T1)<1.D-5) T1=1.D-5
        T1=ABS(Q1/T1) !!!!abs
    ENDIF
    
    KM1(1,1)=1.D0;KM1(2,2)=1.D0
    KM1(2,1)=-1.D0;KM1(1,2)=-1.D0
    T1=(T1+1.0/element(ielt).property(1))/2.0 !!!!
    KM1=KM1*T1
    element(ielt).property(1)=1/T1

    IF(.NOT.ALLOCATED(ELEMENT(IELT).FLUX)) ALLOCATE(ELEMENT(IELT).FLUX(ELEMENT(IELT).NNUM))
    
    ELEMENT(IELT).FLUX=MATMUL(KM1,NHEAD1)
    
    !if(solver_control.bfgm==continuum.or.solver_control.bfgm==consistent) then
    ELEMENT(IELT).KM=KM1    							
    !end if
    
ENDSUBROUTINE

SUBROUTINE SPHFLOW_Q_K_UPDATE(STEPDIS,IELT,ISTEP,IITER)
    
    USE solverds
    USE MESHGEO,ONLY:SNADJL
    IMPLICIT NONE
    INTEGER,INTENT(IN)::IELT,ISTEP,IITER
    REAL(8),INTENT(IN)::STEPDIS(NDOF)
    
    INTEGER::I,J,K,IELT1,N1,N2,N3
    REAL(8)::XV1(3),T1,HK1,PHI1,TPHI1,WR1,L1,DV1(3),PI1,RO1,Q1,X1(3),H1,LAMDA1,K1,KM1(2,2),NHEAD1(2)
    REAL(8)::SITA1(3),KX1,KY1,KZ1,KR1,QN1(2),T2
    LOGICAL::ISISOTROPIC=.FALSE.
   
        
    N1=ELEMENT(IELT).NODE(2)
    TPHI1=0.D0;HK1=0.D0
    DV1(1:NDIMENSION)=ELEMENT(IELT).PROPERTY(1:NDIMENSION)
    WR1=MATERIAL(ELEMENT(IELT).MAT).PROPERTY(1)
    PI1=PI();Q1=0.D0
    NHEAD1=STEPDIS(ELEMENT(IELT).G)
    QN1(1:2)=MATMUL(ELEMENT(IELT).KM,NHEAD1(1:2))
    QWAN(2,ELEMENT(IELT).NODE(1:2))=QWAN(2,ELEMENT(IELT).NODE(1:2))+QN1(1:2)    
    
    IF(ANY(NHEAD1-NODE(ELEMENT(IELT).NODE).COORD(NDIMENSION)<0.D0)) THEN
        T1=1.D-7
    ELSE        
    
        DO J=1,SNADJL(N1).ENUM        
            IELT1=SNADJL(N1).ELEMENT(J)
            IF(ELEMENT(IELT1).ET==WELLBORE &
            .OR.ELEMENT(IELT1).ET==PIPE2 &
            .OR.ELEMENT(IELT1).ET==SPHFLOW &
            .OR.ELEMENT(IELT1).ET==SEMI_SPHFLOW &
            .OR.ELEMENT(IELT1).ET==WELLBORE_SPGFACE) CYCLE !ITSELF EXCLUDED
        
            IF(ELEMENT(IELT1).EC/=SPG) CYCLE
            
            
            IF(NORM2(MATERIAL(ELEMENT(IELT1).MAT).PROPERTY(1:NDIMENSION)-MATERIAL(ELEMENT(IELT1).MAT).PROPERTY(1))<1E-10) THEN
                K=5 !isotropic
                ISISOTROPIC=.TRUE.
            ELSE
                K=1 !anisotropic
                ISISOTROPIC=.FALSE.
            ENDIF
            KR1=0.D0             
            DO K=K,5
                X1(1)=DOT_PRODUCT(NODE(ELEMENT(IELT1).NODE(1:4)).COORD(1),GQ_TET4_O3(1:4,K))
                X1(2)=DOT_PRODUCT(NODE(ELEMENT(IELT1).NODE(1:4)).COORD(2),GQ_TET4_O3(1:4,K))
                X1(3)=DOT_PRODUCT(NODE(ELEMENT(IELT1).NODE(1:4)).COORD(3),GQ_TET4_O3(1:4,K))
                IF(K==5) THEN !K==5 IS CENTER OF THE ELEMENT
                    H1=DOT_PRODUCT(STEPDIS(ELEMENT(IELT1).G(1:4)),GQ_TET4_O3(1:4,K))                    
                ENDIF  
                XV1=X1-NODE(N1).COORD !ORIGIN .                
                RO1=NORM2(XV1)                   
                
                !DIRECTIONAL K
                SITA1=XV1/RO1            
                K1=0.D0
                DO I=1,3
                    K1=K1+1.0/MATERIAL(ELEMENT(IELT1).MAT).PROPERTY(I)*(SITA1(I))**2
                ENDDO
                IF(ISISOTROPIC) THEN
                    KR1=KR1+1/K1
                ELSE
                    KR1=KR1+1/K1*GQ_TET4_O3(5,K)*6.D0
                ENDIF
            ENDDO
      
            !DO I=1,3
            !    X1(I)=SUM(NODE(ELEMENT(IELT1).NODE).COORD(I))
            !ENDDO
            !X1=X1/ELEMENT(IELT1).NNUM
            !XV1=X1-NODE(N1).COORD
            !RO1=NORM2(XV1)
            
            IF(ELEMENT(IELT).ET==SEMI_SPHFLOW) THEN
                T1=DOT_PRODUCT(XV1(1:NDIMENSION),DV1(1:NDIMENSION))                
                IF(T1<0.D0) CYCLE
            ENDIF
            !H1=SUM(STEPDIS(ELEMENT(IELT1).G))/ELEMENT(IELT1).NNUM        
            PHI1=ELEMENT(IELT1).ANGLE(SNADJL(N1).SUBID(J))         
            TPHI1=TPHI1+PHI1
            
            CALL lamda_spg(IELT1,H1,X1(NDIMENSION),lamda1,ISTEP)
            
            !SITA1=XV1/RO1            
            !K1=0.D0
            !DO I=1,3
            !    K1=K1+1.0/MATERIAL(ELEMENT(IELt1).MAT).PROPERTY(I)*(SITA1(I))**2
            !ENDDO
            !K1=1/K1
            
            
            KR1=LAMDA1*KR1 
            IF(SOLVER_CONTROL.well_bottom_type==0) THEN
                T2=4.0
            ELSE
                T2=2.0*PI1
            ENDIF
            Q1=Q1+(H1-NHEAD1(1))*T2*KR1/(1/WR1-1/RO1)*PHI1        
        ENDDO
        Q1=Q1/TPHI1
        IF(ELEMENT(IELT).ET==SPHFLOW) Q1=2.D0*Q1
        QWAN(1,ELEMENT(IELT).NODE(2))=QWAN(1,ELEMENT(IELT).NODE(2))+Q1
        IF(ABS(Q1)<1.D-5) Q1=1.D-5
        T1=(NHEAD1(2)-NHEAD1(1))
        IF(ABS(T1)<1.D-5) T1=1.D-5
        T1=ABS(Q1/T1) !!!!abs
        
    ENDIF
    
    KM1(1,1)=1.D0;KM1(2,2)=1.D0
    KM1(2,1)=-1.D0;KM1(1,2)=-1.D0
    T1=(T1+1.0/element(ielt).property(1))/2.0 !!!!
    KM1=KM1*T1
    element(ielt).property(1)=1/T1

    IF(.NOT.ALLOCATED(ELEMENT(IELT).FLUX)) ALLOCATE(ELEMENT(IELT).FLUX(ELEMENT(IELT).NNUM))
    
    ELEMENT(IELT).FLUX=MATMUL(KM1,NHEAD1)
    
    !if(solver_control.bfgm==continuum.or.solver_control.bfgm==consistent) then
    ELEMENT(IELT).KM=KM1    							
    !end if
    
ENDSUBROUTINE

SUBROUTINE WELL_ELEMENT_OUT(ISTEP,ISUBTS,IITER)
    USE solverds
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ISTEP,ISUBTS,IITER
    INTEGER::I,J,ISET1,NDOF1,FILE_UNIT
    REAL(8)::Q1(4),H1(4),t1
    CHARACTER(16)::ETNAME1(401:405)=["PIPE2","WELLBORE","SPHFLOW","SEMI_SPHFLOW","WELLBORE_SPGFACE"]
    
    IF(.NOT.ISOUT_WELL_FILE) RETURN
    
    FILE_UNIT=30    
    IF(ISTEP==1.AND.ISUBTS==1) THEN	   
        OPEN(FILE_UNIT,FILE=well_file,STATUS='REPLACE')        
    ELSE
        OPEN(FILE_UNIT,FILE=well_file,STATUS='UNKNOWN',ACCESS='APPEND')	
    END IF
    
    DO i=1,neset
	    iset1=esetid(i)
	    if(sf(eset(iset1).sf).factor(ISTEP)==0) cycle
        IF(.NOT.(ANY([WELLBORE,WELLBORE_SPGFACE,PIPE2,SPHFLOW,SEMI_SPHFLOW]-ESET(ISET1).ET==0))) CYCLE
        
        IF(ESET(ISET1).ET==WELLBORE.OR.ESET(ISET1).ET==WELLBORE_SPGFACE) THEN
            WRITE(FILE_UNIT,10) 
        ELSE
            WRITE(FILE_UNIT,20) 
        ENDIF
        
	    DO j=eset(iset1).enums,eset(iset1).enume
            H1(1:ELEMENT(J).NDOF)=TDISP(ELEMENT(J).G)
            Q1(1:ELEMENT(J).NDOF)=MATMUL(ELEMENT(J).KM,H1(1:ELEMENT(J).NDOF))
            IF(ESET(ISET1).ET==WELLBORE.OR.ESET(ISET1).ET==WELLBORE_SPGFACE) THEN
                WRITE(FILE_UNIT,11) ISTEP,ISUBTS,J,ESET(ISET1).ET,NODE(ELEMENT(J).NODE(1)).COORD,&
                NODE(ELEMENT(J).NODE(2)).COORD,H1(1:ELEMENT(J).NDOF),Q1(1:ELEMENT(J).NDOF),&
                MATERIAL(ELEMENT(J).MAT).PROPERTY(1),ELEMENT(J).PROPERTY([1,4,2,3,5]),&
                ABS(ELEMENT(J).KM(1,2)),ABS(ELEMENT(J).KM(2,3)),ABS(ELEMENT(J).KM(1,4))
            ELSE
                t1=MATERIAL(ELEMENT(J).MAT).PROPERTY(1);
                if(isporeflow) t1=element(j).property(2)
                WRITE(FILE_UNIT,21) ISTEP,ISUBTS,J,ETNAME1(ESET(ISET1).ET),NODE(ELEMENT(J).NODE(1)).COORD,&
                NODE(ELEMENT(J).NODE(2)).COORD,H1(1:ELEMENT(J).NDOF),Q1(1),&
                t1,ELEMENT(J).PROPERTY(1),&
                ELEMENT(J).KM(1,1)
            ENDIF
        ENDDO
        
    ENDDO
    
    CLOSE(FILE_UNIT)
    
10 FORMAT(2X,'ISTEP',2X,'ISUBTS',3X,'ELTNO',11X,'ET',9X,'X1',14X,'Y1',14X,'Z1',14X,'X2',14X,'Y2',14X,'Z2',14X,'H1',14X,'H2',14X,'H3',14X,'H4',14X,'Q1',14X,'Q2',14X,'Q3',14X,'Q4',14X,'RW',7X,'Rfriction',12X,'Racc',11X,'Rgeo3',11X,'Rgeo4',4X,'Rskin(N3=N4)',13X,'|K12|',13X,'|K23|',13X,'|K14|')    
11 FORMAT(3(I7,1X),A16,1X,23(G15.8,1X))  
20 FORMAT(2X,'ISTEP',2X,'ISUBTS',3X,'ELTNO',11X,'ET',9X,'X1',14X,'Y1',14X,'Z1',14X,'X2',14X,'Y2',14X,'Z2',14X,'H1',14X,'H2',14X,'Q',14X,'RW',11X,'Rgeo',11X,'|K11|')    
21 FORMAT(3(I7,1X),(A12,1X),12(G15.8,1X)) 
   
ENDSUBROUTINE

SUBROUTINE DIRECTION_K(KR,IEL,IWN,Vec)
!求含井节点IWN的单元elt的方向渗透系数，利用3阶四面体单元的采用点
!如果VEC出现，则方向为垂直于VEC，否则，方向指向井节点IWN
    USE solverds
    IMPLICIT NONE
    
    
    INTEGER,INTENT(IN)::IEL,IWN
    REAL(8),INTENT(IN),OPTIONAL::VEC(3)
    REAL(8),INTENT(OUT)::KR

    INTEGER::I,J,K
    REAL(8)::X1(3),LP1,SITA1(3),R1,K1
    LOGICAL::ISISOTROPIC=.FALSE.

    IF(NORM2(MATERIAL(ELEMENT(IEL).MAT).PROPERTY(1:NDIMENSION)-MATERIAL(ELEMENT(IEL).MAT).PROPERTY(1))<1E-10) THEN
        K=5 !isotropic
        ISISOTROPIC=.TRUE.
    ELSE
        K=1 !anisotropic
        ISISOTROPIC=.FALSE.
    ENDIF
    KR=0.D0             
    DO K=K,5
        X1(1)=DOT_PRODUCT(NODE(ELEMENT(IEL).NODE(1:4)).COORD(1),GQ_TET4_O3(1:4,K))
        X1(2)=DOT_PRODUCT(NODE(ELEMENT(IEL).NODE(1:4)).COORD(2),GQ_TET4_O3(1:4,K))
        X1(3)=DOT_PRODUCT(NODE(ELEMENT(IEL).NODE(1:4)).COORD(3),GQ_TET4_O3(1:4,K))

                    
        X1=X1-NODE(ELEMENT(IEL).NODE(IWN)).COORD !ORIGIN IS THE THIRD NODE.                
        IF(PRESENT(VEC)) THEN
            LP1=DOT_PRODUCT(X1,VEC) 
            X1=X1-LP1*VEC   !投影线向量，为此单元的流量方向矢量            
        ENDIF
        R1=NORM2(X1)                   
                
        !DIRECTIONAL K
        SITA1=X1/R1            
        K1=0.D0
        DO I=1,3
            K1=K1+1.0/MATERIAL(ELEMENT(IEL).MAT).PROPERTY(I)*(SITA1(I))**2
        ENDDO
        IF(ISISOTROPIC) THEN
            KR=KR+1/K1
        ELSE
            KR=KR+1/K1*GQ_TET4_O3(5,K)*6.D0
        ENDIF
    ENDDO
ENDSUBROUTINE
    
SUBROUTINE ZT_SPG_INI2(IELT)
    USE solverds
    USE SolverMath
    IMPLICIT NONE
    INTEGER,INTENT(IN)::IELT
    !INTEGER,INTENT(INOUT)::ISPF(:,:)
    integer::i,j,K
	integer::n1,n2,n3,n4,AELT1(2),IN1(10)
	real(kind=DPN)::t1=0,vcos=0,vsin=0,rpi,coord1(3,4)=0,cent1(3,2)=0,c1(3,3)=0,VEC1(3),t2
    
    
    if(.not.isIniSEdge) CALL Model_MESHTOPO_INI()
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