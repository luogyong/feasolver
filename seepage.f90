SUBROUTINE SPG_Q_UPDATE(STEPDIS,bload,HHEAD,INIHEAD,DT,nbload,ienum,iiter,istep,iscon)
    USE solverds
    IMPLICIT NONE
    logical,intent(in)::iscon
	integer,intent(in)::iiter,ienum,istep,nbload
	real(kind=DPN),intent(in)::STEPDIS(NDOF),HHEAD(nbload),INIHEAD(nbload),DT
	real(kind=DPN),intent(out)::bload(nbload)
    INTEGER::J,N1,ND1
    REAL(DPN)::KT1(NDIMENSION,NDIMENSION),HJ,hj_ini,slope1,sita1,R1,A1,A2,A3
    
    
    SELECT CASE(ELEMENT(IENUM).ET)
        
    CASE(PIPE2,WELLBORE,WELLBORE_SPGFACE)
        CALL WELLBORE_Q_K_UPDATE(STEPDIS,IENUM,ISTEP,IITER)
        BLOAD=ELEMENT(IENUM).FLUX 
    CASE(SPHFLOW,SEMI_SPHFLOW)
        CALL SPHFLOW_Q_K_UPDATE(STEPDIS,IENUM,ISTEP,IITER)
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

            if(solver_control.bfgm==continuum.or.solver_control.bfgm==consistent) then
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
    
    
subroutine lamda_spg(ienum,hj,z,lamda,ISTEP)
	
	use solverds
	implicit none
	integer,intent(in)::ienum,ISTEP
	real(8),intent(in)::hj,z
	real(8),intent(out)::lamda
	real(8)::t1,epsilon1,epsilon2,krsml
	real(8)::alpha1,fn1,fm1,seff1,scale1=1.0
	
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
	REAL(DPN)::ROT1(2,2),COS1,SIN1,T1,FAC1,SITA1,LAMDA,C1,PHI1,SFR1(6),SIGMA1(6),TS1,PI1
    INTEGER::IENUMC1,IE1	
    
    
	
					
	call lamda_spg(ienum,hj,element(ienum).xygp(ndimension,IGP),lamda,ISTEP)                   						
                    
    if(iiter>1) lamda=(element(ienum).kr(IGP)+lamda)/2.0D0
                    
    element(ienum).kr(IGP)=lamda 
    
    KT=ELEMENT(IENUM).D*LAMDA
    
	IF(ESET(ELEMENT(IENUM).SET).COUPLESET<0.OR.ESET(ELEMENT(IENUM).SET).COUPLESET==ELEMENT(IENUM).SET) RETURN
    
    PI1=3.14159265358979
    IE1=IENUM-(ESET(ELEMENT(IENUM).SET).ENUMS-1)
    IENUMC1=ESET(ESET(ELEMENT(IENUM).SET).COUPLESET).ENUMS-1+IE1
    C1=material(ELEMENT(IENUMC1).MAT).GET(3,ISTEP)
	Phi1=material(ELEMENT(IENUMC1).MAT).GET(4,ISTEP)
    TS1=material(ELEMENT(IENUMC1).MAT).GET(5,ISTEP)
    IF(material(ELEMENT(IENUMC1).MAT).TYPE==MC)THEN
        TS1=C1/MAX(tan(phi1/180*PI1),1E-7) 
    ENDIF
    !!!!HERE,THE STRESS IS STILL NOT UPDATED.
    SIGMA1=ELEMENT(IENUMC1).STRESS(:,IGP)+ELEMENT(IENUMC1).DSTRESS(:,IGP)
	call stress_in_failure_surface(sfr1,SIGMA1,2,C1,Phi1,solver_control.slidedirection,ELEMENT(IENUMC1).XYGP(1:NDIMENSION,IGP),TS1)
        
        

    FAC1=MAX((ABS(SFR1(1)))**solver_control.slope_kscale,1.D-9)
    SITA1=-SFR1(2)/180.*PI()
    !element(ienum).kr(IGP)=element(ienum).kr(IGP)*fac1
    
	!KT=0.D0
	!T1=MATERIAL(ELEMENT(IENUM).MAT).GET(1,ISTEP)
	!KT(2,2)=MATERIAL(ELEMENT(IENUM).MAT).GET(2,ISTEP)*FAC1
	KT=KT*(solver_control.slope_kbase+FAC1)
    IF(ABS(KT(1,1)-KT(2,2))>1E-7) THEN        
	    COS1=COS(SITA1);SIN1=SIN(SITA1)
	    ROT1(1,1)=COS1;ROT1(2,2)=COS1;ROT1(1,2)=SIN1;ROT1(2,1)=-SIN1
	    KT=MATMUL(MATMUL(ROT1,KT),TRANSPOSE(ROT1))
    ENDIF
	

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
    KR1,XC1(3),VW1(3),T1,C1,B1,PO1,VS1
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
 
    
    IF(ET1/=WELLBORE_SPGFACE) THEN
        L1=2*DIS1
        D1=2*MATERIAL(ELEMENT(IEL).MAT).PROPERTY(1)
        AREA1=PI1*D1**2/4.0
        QA1=ABS((NHEAD1(1)-NHEAD1(2))*ELEMENT(IEL).KM(1,1))
        VA1=QA1/AREA1
        if(ABS(MATERIAL(ELEMENT(IEL).MAT).PROPERTY(5))<1.E-7) THEN
            g1=73156608000.00 
        else
            g1=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(5)
        endif
            
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
            IF(ABS(MATERIAL(ELEMENT(IEL).MAT).PROPERTY(4))<1.D-7) THEN
                CF1=1./24/3600 !to m/s
            ELSE
                CF1=1.0D0
            ENDIF
            RE1=Re_W(VA1*CF1,D1,MATERIAL(ELEMENT(IEL).MAT).PROPERTY(2))
            VR1=QR1/(L1*PI1*D1)
            REW1=Re_W(VR1*CF1,D1,MATERIAL(ELEMENT(IEL).MAT).PROPERTY(2))
            FD1=FD_PF(RE1,MATERIAL(ELEMENT(IEL).MAT).PROPERTY(3),REW=REW1,POROSITY=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(9),MODEL=MODEL1) 
        ENDIF
        
        FD1=FD1*L1*QA1/(2*D1*G1*AREA1**2)
       
        
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
        
        DO I=1,2
    
            IF(NHEAD1(2+I)-NODE(ELEMENT(IEL).NODE(2+I)).COORD(NDIMENSION)<0.D0) THEN        
                A1(I+1)=1.D-7
                A1(1)=1.D-7 !!!!!
                CYCLE
            ENDIF
        
            IN1=ELEMENT(IEL).NODE(2+I)
            IN2=ELEMENT(IEL).NODE(2+MOD(I,2)+1)
            !DIS1=NORM2(NODE(IN1).COORD-NODE(IN2).COORD)/2.
            IF(I==1) THEN
                N2=2
            ELSE
                N2=1
            ENDIF
            IF(ET1==WELLBORE) THEN
                DH1=STEPDIS(NODE(IN1).DOF(4))-STEPDIS(NODE(ELEMENT(IEL).NODE(N2)).DOF(4))
            ELSE
                !WELLBORE_SPGFACE
                DH1=STEPDIS(NODE(IN1).DOF(4))-NODE(ELEMENT(IEL).NODE(N2)).COORD(NDIMENSION)
            ENDIF
            
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
                
                    HZ1=HZ1+QN1(I+2)*ELEMENT(IEL).PROPERTY(5)
                
                    RA1=ELEMENT(IEL).PROPERTY(6)*KR1*(H1-HZ1)/ &
                    LOG(R1/MATERIAL(ELEMENT(IEL).MAT).PROPERTY(1))
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
            Q1(I+2)=Q1(I+2)/TW1*DIS1
            !A1(I+1)=ABS(Q1(I+2)/DH1)
            
            IF(ABS(Q1(I+2))<1E-10) THEN
                Q1(I+2)=1.D-10
            ENDIF
            
            IF(ABS(DH1)<1.D-10) DH1=1.D-10
            
            ELEMENT(IEL).PROPERTY(I+1)=ABS(DH1/Q1(I+2)) !!!

            
            !井损
            A1(I+1)=1.0/(ELEMENT(IEL).PROPERTY(I+1)+ELEMENT(IEL).PROPERTY(5))
        ENDDO

        IF(ANY(ISNAN(A1(1:3)))) THEN
            PRINT *, 'NAN'
        ENDIF
        
         KM1(1:NDOF1,1:NDOF1)=KM_WELLBORE(A1(1),A1(2),A1(3))

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

SUBROUTINE INI_WELLBORE(IELT)
    
    USE MESHADJ,ONLY:SETUP_EDGE_ADJL,SEDGE,NSEDGE,SNADJL,GETGMSHET,ELTTYPE
    USE solverds
    USE SolverMath
    IMPLICIT NONE
    INTEGER,INTENT(IN)::IELT
    INTEGER::I,J,IELT1,IEDGE1,GMET1,SUBID1,AF1(2),N1,NADJL1(200),N2,N3
    REAL(8)::AV1(3,2),AR1(3,3),ZV1(3),YV1(3),XV1(3),D1,T1,T2,RDIS1(200),HK1,PHI1,TPHI1,WR1,L1,&
        ORG1(3),AREA1
    
    if(.not.isIniSEdge) then
        CALL SETUP_EDGE_ADJL(SEDGE,NSEDGE,SNADJL)
        ISINISEDGE=.TRUE.
    endif
    
    IEDGE1=ELEMENT(IELT).EDGE(3)
    N1=0
    L1=SEDGE(IEDGE1).DIS/2.0
    IF(SEDGE(IEDGE1).ENUM<2) THEN
        PRINT *, "NO ELEMENTS SURROUNDING THE WELLBORE ELEMENT I. I=",IELT
        PRINT *, "TRY TO EMBED THE WELL LINE INTO THE CORRESPOINDING VOLUME."
        STOP
    ENDIF
    DO I=1,SEDGE(IEDGE1).ENUM
        IELT1=SEDGE(IEDGE1).ELEMENT(I)
        GMET1=GETGMSHET(ELEMENT(IELT1).ET)
        IF(ELTTYPE(GMET1).DIM==3.AND.ELEMENT(IELT1).EC==SPG) THEN
            N1=N1+1
            IF(.NOT.ALLOCATED(ELEMENT(IELT1).ANGLE)) CALL calangle(IELT1) 
        ENDIF
    ENDDO
    
       
    ALLOCATE(ELEMENT(IELT).NODE2(N1),ELEMENT(IELT).ANGLE(N1))
    
    DO I=3,4
        DO J=1,SNADJL(ELEMENT(IELT).NODE(I)).ENUM
            IELT1=SNADJL(ELEMENT(IELT).NODE(I)).ELEMENT(J)
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
    ZV1=NODE(ELEMENT(IELT).NODE(4)).COORD-ORG1
    ZV1=ZV1/NORM2(ZV1)
    D1=-DOT_PRODUCT(ZV1,ORG1)
    
    IF(ABS(ZV1(3))>1E-7) THEN
        XV1=[ORG1(1)+1.0,0.D0,-(D1+ZV1(1)*(ORG1(1)+1.0))/ZV1(3)]
    ELSEIF(ABS(ZV1(2))>1E-7) THEN
        XV1=[ORG1(1)+1.0,-(D1+ZV1(1)*(ORG1(1)+1.0))/ZV1(2),0.D0]
    ELSE
        XV1=[-(D1+ZV1(2)*(ORG1(2)+1.0))/ZV1(1),ORG1(2)+1.0,0.D0]
    ENDIF
    XV1=XV1-ORG1
    XV1=XV1/NORM2(XV1)
    
    YV1=NORMAL_TRIFACE(RESHAPE([0.D0,0.D0,0.D0,ZV1,XV1],([3,3])))       
    YV1=YV1/NORM2(YV1)
    ELEMENT(IELT).G2L(1,:)=XV1;ELEMENT(IELT).G2L(2,:)=YV1;ELEMENT(IELT).G2L(3,:)=ZV1;
    
    WR1=MATERIAL(ELEMENT(IELT).MAT).PROPERTY(1)
    T2=PHI1*SIN(PHI1)/(1-COS(PHI1))
    DO J=3,4
        N2=0
        N3=ELEMENT(IELT).NODE(J)
        DO I=1,SNADJL(N3).NNUM
            N1=SNADJL(N3).NODE(I)
            XV1=NODE(N1).COORD-NODE(ELEMENT(IELT).NODE(3)).COORD
            XV1=MATMUL(ELEMENT(IELT).G2L(:,:),XV1)
            T1=NORM2(XV1(1:2))
            IF(T1>WR1) THEN
                N2=N2+1
                NADJL1(N2)=N1;
                RDIS1(N2)=T1
            ENDIF
           
        ENDDO
        T1=SUM(RDIS1(1:N2))/N2

        IF(J==3) THEN
            !ALLOCATE(ELEMENT(IELT).WELL_SP3(N2),ELEMENT(IELT).WELL_SP3_R(N2))
            !ELEMENT(IELT).NWSP3=N2
            !ELEMENT(IELT).WELL_SP3=NADJL1(1:N2)
            !ELEMENT(IELT).WELL_SP3_R=RDIS1(1:N2)
            ELEMENT(IELT).PROPERTY(2)=1./(TPHI1*HK1*L1/(LOG(T1/WR1)-T2))
        ELSE
            !ALLOCATE(ELEMENT(IELT).WELL_SP4(N2),ELEMENT(IELT).WELL_SP4_R(N2))
            !ELEMENT(IELT).NWSP4=N2
            !ELEMENT(IELT).WELL_SP4=NADJL1(1:N2)
            !ELEMENT(IELT).WELL_SP4_R=RDIS1(1:N2)
            ELEMENT(IELT).PROPERTY(3)=1./(TPHI1*HK1*L1/(LOG(T1/WR1)-T2))
        ENDIF
        
    ENDDO
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
    IMPLICIT NONE
    INTEGER,INTENT(IN)::IELT
    INTEGER::I,J,IELT1,N1,N2,N3
    REAL(8)::XV1(3),T1,HK1,PHI1,TPHI1,WR1,L1,DV1(3),PI1
    
    if(.not.isIniSEdge) then
        CALL SETUP_EDGE_ADJL(SEDGE,NSEDGE,SNADJL)
        ISINISEDGE=.TRUE.
    endif
    
        
    N1=ELEMENT(IELT).NODE(2)
    TPHI1=0.D0;HK1=0.D0
    DV1(1:NDIMENSION)=ELEMENT(IELT).PROPERTY(1:NDIMENSION)
    WR1=MATERIAL(ELEMENT(IELT).MAT).PROPERTY(1)
    DO J=1,SNADJL(N1).ENUM        
        IELT1=SNADJL(N1).ELEMENT(J)
        IF(ELEMENT(IELT1).ET==WELLBORE &
        .OR.ELEMENT(IELT1).ET==PIPE2 &
        .OR.ELEMENT(IELT1).ET==SPHFLOW &
        .OR.ELEMENT(IELT1).ET==SEMI_SPHFLOW &
        .OR.ELEMENT(IELT1).ET==WELLBORE_SPGFACE) CYCLE !ITSELF EXCLUDED
        
        IF(ELEMENT(IELT1).EC/=SPG) CYCLE
        IF(ALLOCATED(ELEMENT(IELT1).ANGLE)) CYCLE
        CALL calangle(IELT1)
        !AVERAGED K, ASSUME ISOTROPIC        
        IF(ELEMENT(IELT).ET==SEMI_SPHFLOW) THEN
            DO I=1,3
                XV1(I)=SUM(NODE(ELEMENT(IELT1).NODE).COORD(I))
            ENDDO
            XV1=XV1/ELEMENT(IELT1).NNUM
            XV1=XV1-NODE(N1).COORD
            T1=DOT_PRODUCT(XV1(1:NDIMENSION),DV1(1:NDIMENSION))/ &
            (NORM2(XV1(1:NDIMENSION))*NORM2(DV1(1:NDIMENSION)))
            IF(T1<0.D0) CYCLE
        ENDIF
        PHI1=ELEMENT(IELT1).ANGLE(SNADJL(N1).SUBID(J)) 
        TPHI1=TPHI1+PHI1
        HK1=HK1+PHI1*(MATERIAL(ELEMENT(IELT1).MAT).PROPERTY(1)*MATERIAL(ELEMENT(IELT1).MAT).PROPERTY(2))**0.5        
    ENDDO
    HK1=HK1/TPHI1
    N2=0;L1=0.D0
    DO J=1,SNADJL(N1).NNUM        
        N3=SNADJL(N1).NODE(J)
        XV1=NODE(N3).COORD-NODE(N1).COORD
        IF(ELEMENT(IELT).ET==SEMI_SPHFLOW) THEN
            T1=DOT_PRODUCT(XV1(1:NDIMENSION),DV1(1:NDIMENSION))/ &
            (NORM2(XV1(1:NDIMENSION))*NORM2(DV1(1:NDIMENSION)))
            IF(T1<0.D0) CYCLE
        ENDIF            
        T1=NORM2(XV1)
        IF(T1>WR1) THEN
            N2=N2+1                
            L1=L1+T1
        ENDIF           
    ENDDO 
    L1=L1/N2 
    
    PI1=PI()
    T1=2*PI1*HK1/(1/WR1-(4*PI1+3)/(3*L1))
    IF(ELEMENT(IELT).ET==SPHFLOW) T1=2*T1
    ELEMENT(IELT).KM(1,1)=1.D0;ELEMENT(IELT).KM(2,2)=1.D0
    ELEMENT(IELT).KM(2,1)=-1.D0;ELEMENT(IELT).KM(1,2)=-1.D0
    ELEMENT(IELT).KM=ELEMENT(IELT).KM*T1
    ELEMENT(IELT).PROPERTY(1)=1/T1
    
ENDSUBROUTINE

SUBROUTINE SPHFLOW_Q_K_UPDATE(STEPDIS,IELT,ISTEP,IITER)
    
    USE solverds
    USE MESHGEO,ONLY:SNADJL
    IMPLICIT NONE
    INTEGER,INTENT(IN)::IELT,ISTEP,IITER
    REAL(8),INTENT(IN)::STEPDIS(NDOF)
    
    INTEGER::I,J,K,IELT1,N1,N2,N3
    REAL(8)::XV1(3),T1,HK1,PHI1,TPHI1,WR1,L1,DV1(3),PI1,RO1,Q1,X1(3),H1,LAMDA1,K1,KM1(2,2),NHEAD1(2)
    REAL(8)::SITA1(3),KX1,KY1,KZ1,KR1
    LOGICAL::ISISOTROPIC=.FALSE.
   
        
    N1=ELEMENT(IELT).NODE(2)
    TPHI1=0.D0;HK1=0.D0
    DV1(1:NDIMENSION)=ELEMENT(IELT).PROPERTY(1:NDIMENSION)
    WR1=MATERIAL(ELEMENT(IELT).MAT).PROPERTY(1)
    PI1=PI();Q1=0.D0
    NHEAD1=STEPDIS(ELEMENT(IELT).G)
    
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
            
            Q1=Q1+(H1-NHEAD1(1))*2.D0*PI1*KR1/(1/WR1-1/RO1)*PHI1        
        ENDDO
        Q1=Q1/TPHI1
        IF(ELEMENT(IELT).ET==SPHFLOW) Q1=2.D0*Q1
    
        T1=ABS(Q1/(NHEAD1(2)-NHEAD1(1))) !!!!abs
    ENDIF
    
    KM1(1,1)=1.D0;KM1(2,2)=1.D0
    KM1(2,1)=-1.D0;KM1(1,2)=-1.D0
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
    REAL(8)::Q1(4),H1(4)
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
        IF(.NOT.(ANY([WELLBORE,PIPE2,SPHFLOW,SEMI_SPHFLOW]-ESET(ISET1).ET==0))) CYCLE
        
        IF(ESET(ISET1).ET==WELLBORE.OR.ESET(ISET1).ET==WELLBORE_SPGFACE) THEN
            WRITE(FILE_UNIT,10) 
        ELSE
            WRITE(FILE_UNIT,20) 
        ENDIF
        
	    DO j=eset(iset1).enums,eset(iset1).enume
            H1(1:ELEMENT(J).NDOF)=TDISP(ELEMENT(J).G)
            Q1(1:ELEMENT(J).NDOF)=MATMUL(ELEMENT(J).KM,H1(1:ELEMENT(J).NDOF))
            IF(ESET(ISET1).ET==WELLBORE) THEN
                WRITE(FILE_UNIT,11) ISTEP,ISUBTS,J,NODE(ELEMENT(J).NODE(1)).COORD,&
                NODE(ELEMENT(J).NODE(2)).COORD,H1(1:ELEMENT(J).NDOF),Q1(1:ELEMENT(J).NDOF),&
                MATERIAL(ELEMENT(J).MAT).PROPERTY(1),ELEMENT(J).PROPERTY([1,4,2,3,5]),&
                ABS(ELEMENT(J).KM(1,2)),ABS(ELEMENT(J).KM(2,3)),ABS(ELEMENT(J).KM(1,4))
            ELSE
                WRITE(FILE_UNIT,21) ISTEP,ISUBTS,J,ETNAME1(ESET(ISET1).ET),NODE(ELEMENT(J).NODE(1)).COORD,&
                NODE(ELEMENT(J).NODE(2)).COORD,H1(1:ELEMENT(J).NDOF),Q1(1),&
                MATERIAL(ELEMENT(J).MAT).PROPERTY(1),ELEMENT(J).PROPERTY(1),&
                ELEMENT(J).KM(1,1)
            ENDIF
        ENDDO
        
    ENDDO
    
    CLOSE(FILE_UNIT)
    
10 FORMAT(2X,'ISTEP',2X,'ISUBTS',3X,'ELTNO',11X,'ET',9X,'X1',14X,'Y1',14X,'Z1',14X,'X2',14X,'Y2',14X,'Z2',14X,'H1',14X,'H2',14X,'H3',14X,'H4',14X,'Q1',14X,'Q2',14X,'Q3',14X,'Q4',14X,'RW',7X,'Rfriction',12X,'Racc',11X,'Rgeo3',11X,'Rgeo4',4X,'Rskin(N3=N4)',13X,'|K12|',13X,'|K23|',13X,'|K14|')    
11 FORMAT(3(I7,1X),'    WELLBORE ',23(G15.8,1X))  
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