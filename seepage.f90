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
        
    CASE(PIPE2,WELLBORE)
        CALL WELLBORE_Q_K_UPDATE(STEPDIS,IENUM,ISTEP,IITER)
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
    INTEGER::I,J,K,IN1,IN2,IN3,N1,N2,IEL1,IP1
    REAL(8)::H1,X1(3),R1,Q1(4),PI1,LAMDA1,K1,DIS1,DH1,A1(3),NHEAD1(4),KM1(4,4),HZ1,RA1,X2(3),&
    W1,TW1
    
    
    NHEAD1=STEPDIS(ELEMENT(IEL).G)
    
    IF(ANY(NHEAD1(1:2)-NODE(ELEMENT(IEL).NODE(1:2)).COORD(NDIMENSION)<0.D0)) THEN
        A1(1)=1.D-7
    ELSE
        A1(1)=ELEMENT(IEL).PROPERTY(1)
    ENDIF
    
    IF(ELEMENT(IEL).ET==PIPE2) THEN
        KM1(1,1)=A1(1);KM1(2,2)=A1(1);KM1(1,2)=-A1(1);KM1(2,1)=-A1(1)
     
    ELSE
        !WRITE(99,10) 
        
        PI1=PI();Q1=0.D0
        DO I=1,2
    
            IF(NHEAD1(2+I)-NODE(ELEMENT(IEL).NODE(2+I)).COORD(NDIMENSION)<0.D0) THEN        
                A1(I+1)=1.D-7
                CYCLE
            ENDIF
        
            IN1=ELEMENT(IEL).NODE(2+I)
            IN2=ELEMENT(IEL).NODE(2+MOD(I,2)+1)
            DIS1=NORM2(NODE(IN1).COORD-NODE(IN2).COORD)/2.
            IF(I==1) THEN
                N2=2
            ELSE
                N2=1
            ENDIF
            DH1=STEPDIS(NODE(IN1).DOF(4))-STEPDIS(NODE(ELEMENT(IEL).NODE(N2)).DOF(4))
            IP1=0;RA1=0.D0;W1=0;TW1=0;
            DO J=1,SNADJL(IN1).ENUM
                IEL1=SNADJL(IN1).ELEMENT(J)
                IF(ELEMENT(IEL1).ET==WELLBORE.OR.ELEMENT(IEL1).ET==PIPE2) CYCLE !ITSELF EXCLUDED
                IF(ELEMENT(IEL1).EC/=SPG) CYCLE
                
                H1=0;X1=0.D0;N1=0
                DO K=1,ELEMENT(IEL1).NNUM
                    IN3=ELEMENT(IEL1).NODE(K)
                    !IF(IN3==IN1) CYCLE
                    !IF(IN3==IN2) CYCLE
                    N1=N1+1
                    H1=H1+STEPDIS(NODE(IN3).DOF(4))
                    X1=X1+NODE(IN3).COORD
                ENDDO
                
                IF(N1<1) CYCLE
                
                H1=H1/N1;X1=X1/N1;X2=X1
                
                CALL lamda_spg(IEL1,H1,X1(NDIMENSION),lamda1,ISTEP)
                K1=MATERIAL(ELEMENT(IEL1).MAT).PROPERTY(1)*LAMDA1 !ASSUME ISOTROPIC
                X1=X1-NODE(ELEMENT(IEL).NODE(3)).COORD !ORIGIN IS THE THIRD NODE.
                X1=MATMUL(ELEMENT(IEL).G2L,X1)
                R1=NORM2(X1(1:2))
                
                IF(R1<=MATERIAL(ELEMENT(IEL).MAT).PROPERTY(1)) CYCLE
                
                HZ1=NHEAD1(2)+(NHEAD1(1)-NHEAD1(2))/(NODE(ELEMENT(IEL).NODE(1)).COORD(NDIMENSION)-NODE(ELEMENT(IEL).NODE(2)).COORD(NDIMENSION))*X1(NDIMENSION)  !X1已经局部坐标
                RA1=ELEMENT(IEL).PROPERTY(4)*K1*(H1-HZ1)/ &
                LOG(R1/MATERIAL(ELEMENT(IEL).MAT).PROPERTY(1))
                W1=ELEMENT(IEL1).ANGLE(SNADJL(IN1).SUBID(J))
                Q1(I+2)=Q1(I+2)+RA1
                TW1=TW1+W1;IP1=IP1+1
                !WRITE(99,20) ISTEP,IITER,IEL,I+2,IP1,X2,X1,R1,H1,HZ1,RA1,W1
            ENDDO  
            Q1(I+2)=Q1(I+2)/IP1*DIS1
            A1(I+1)=ABS(Q1(I+2)/DH1)
        ENDDO
        KM1=KM_WELLBORE(A1(1),A1(2),A1(3))
    ENDIF
    
    
    IF(.NOT.ALLOCATED(ELEMENT(IEL).FLUX)) ALLOCATE(ELEMENT(IEL).FLUX(ELEMENT(IEL).NNUM))
    
    ELEMENT(IEL).FLUX=MATMUL(KM1(1:ELEMENT(IEL).NDOF,1:ELEMENT(IEL).NDOF),NHEAD1(1:ELEMENT(IEL).NDOF))
    
    if(solver_control.bfgm==continuum.or.solver_control.bfgm==consistent) then
        ELEMENT(IEL).KM=KM1(1:ELEMENT(IEL).NDOF,1:ELEMENT(IEL).NDOF)      							
    end if

10  FORMAT("ISTEP,IITER,IEL,NODE1,P1,XG,YG,ZG,XL,YL,ZL,RA,HA,HW,Qunit,WGT   //WELLBORE INFO")
20  FORMAT(5I7,11F15.7)
    
END SUBROUTINE