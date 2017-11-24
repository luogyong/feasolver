subroutine stree_failure_ratio_cal(ienum,ISTEP)
	use solverds
	implicit none
	integer,intent(in)::ienum,ISTEP
	integer::i,et1,nsh1,ndim1,dof1(3)
	real(8)::ss1(6),pss1(4),t1,t2,C1,phi1,R1,PI1,dis1(3,30),disg1(3)
	
	PI1=Pi()
	et1=element(ienum).et
	nsh1=ecp(et1).nshape
	ndim1=ecp(et1).ndim
	!do concurrent(i=1:element(ienum).nnum)	
	!	dis1(1:ndim1,i)=Tdisp(node(element(ienum).node(i)).dof(1:ndim1))		
	!enddo
		!MC material		
		phi1=material(element(ienum).mat).GET(4,ISTEP)
		c1=material(element(ienum).mat).GET(3,ISTEP)
		
	do i=1,element(ienum).ngp
		ss1=element(ienum).stress(:,i)
		!积分点的位移
		!disg1(1:ndim1)=matmul(dis1(1:ndim1,1:nsh1),ecp(et1).Lshape(:,i))		

		call stress_in_failure_surface(element(ienum).sfr(:,i),ss1,ecp(element(ienum).et).ndim,c1,Phi1,solver_control.slidedirection,&
                                    element(ienum).xygp(ndimension,i))
		
        
	enddo
	
end subroutine

subroutine stress_in_failure_surface(sfr,stress,ndim,cohesion,PhiD,slidedirection,Y)
    USE DS_SlopeStability
	implicit none
	integer,intent(in)::ndim,slidedirection
	real(8),intent(in)::stress(6),cohesion,PhiD,Y
	real(8),intent(inout)::sfr(6)
	
	integer::i
	real(8)::ss1(6),pss1(4),t1,t2,C1,phi1,PI1,sin1,cos1,t3,t4,t5,ta1(4),A,B
	real(8)::sigmaC1,R1,sita1
    
	PI1=dATAN(1.0)*4.0
	C1=cohesion
	phi1=PhiD/180*PI1
	sfr=0.d0
   
	call principal_stress_cal(pss1,stress,ndim)
	
	t2=(pss1(1)-pss1(2))/2.0
    
    sigmaC1=(stress(1)+stress(2))/2.0
    R1=t2
    if(R1>0) then
        if(pss1(1)<=0) then
        !if(sigmaC1<=0) then
            
            B=-tan(phi1)
            A=(C1+sigmaC1*B)/R1
            !最大破坏面与主应力面的夹角
            if(r1<=C1*cos(phi1)-sigmaC1*sin(phi1)) then
                sita1=asin((1.d0-(B/A)**2)**0.5)
            else
                sita1=pi1/2.0-phi1
            endif
            sfr(1)=sin(sita1)/(A+B*cos(sita1)) !sfr for stress failure ratio
            if(sfr(1)>1.d0) sfr(1)=1.0d0
        else
            !tension
            sfr(1)=-1.            
            sita1=pi1/2.0
        endif
        !if(isnan(sfr(1))) then
        !    pause
        !endif
    else
        sfr(1)=0.d0
        sita1=0.d0
    endif      
    
    !if((pss1(1)+pss1(2))<0.d0) then 
    !
	   ! !if(phi1>0.d0) then
		  ! ! t1=(-(pss1(1)+pss1(2))/2.0+c1/dtan(phi1))*dsin(phi1) !受拉为正			
	   ! !else
		  ! ! t1=C1
	   ! !endif
	   ! !IF(ABS(T1)>1E-14) THEN
    ! !       sfr(1)=t2/t1
    ! !   ELSE
    ! !       SFR(1)=1
    ! !   ENDIF
    !  
    !    
    !    
    !else !>0.拉坏
    !    sfr(1)=-1. !表拉坏
    !endif
    !
    t1=slidedirection
    !一般在坡顶和坡脚，会出现sxx<syy(水平方向的压力大于竖向)的情况，按下面的方法计算，其破坏面的夹角与x轴的夹角大于0
    !这时，当slidedirection=right时，坡顶会出现不相容的破坏面方向(sfr(2)>0),这时令其取另一个方向。
    IF(slidedirection==1.and.SLOPEPARAMETER.ISYDWZ.AND.Y>SLOPEPARAMETER.YDOWNWARDZONE.AND.stress(2)>stress(1)) THEN
        T1=-T1
    ENDIF

    IF(sfr(1)>=0.d0) then
        if(stress(2)>stress(1)) then !sigmax<sigmay         
            sfr(2)=pss1(4)+(sita1)/2.0*T1
        else
            sfr(2)=pss1(4)-(PI1-sita1)/2.0*T1
        endif
        !sfr(3)=(pss1(1)+pss1(2))/2+t2*cos(PI1/2-phi1)
        !sfr(4)=-t2*sin(PI1/2-phi1)*t1 
        sfr(3)=(pss1(1)+pss1(2))/2+t2*cos(sita1)
        sfr(4)=-t2*sin(sita1)*t1 
        sfr(5)=sign(sfr(1),-sfr(4))*dcos(sfr(2)) !输出归一化的剪应力
        sfr(6)=sign(sfr(1),-sfr(4))*dsin(sfr(2))
        !sfr(5)=abs(sfr(1))*dcos(sfr(2)) !输出归一化的剪应力
        !sfr(6)=abs(sfr(1))*dsin(sfr(2))
    ELSE
        !TENSION FAILURE,THE FAILURE SURFACE IS PLANE THE MAJOR PRINCIPLE STRESS.
        IF(stress(2)>stress(1)) THEN
            SFR(2)=pss1(4)
        ELSE
            SFR(2)=pss1(4)+PI1/2.0
        ENDIF
        SFR(5)=0.D0
        SFR(6)=-1
    ENDIF
    
    sfr(2)=sfr(2)/PI1*180.0 !to degree
    
    
    !t1=0.5*(stress(1)+stress(2))
    !t2=0.5*(stress(1)-stress(2))	
    !t3=2*PI1
    !do i=1,2
    !    !找出破坏面上的剪应力与位移矢量较小的破坏面
    !    if(i==1) then
    !        t5=dasin(dsin(pss1(4)-pi1/2.-(pi1/4.+phi1/2.)))
    !    else
    !        t5=dasin(dsin(pss1(4)-pi1/2.+(pi1/4.+phi1/2.)))
    !    endif
    !
    !    sin1=dsin(2*(t5+PI1/2.));cos1=dcos(2.*(t5+PI1/2.))
    !    ta1(1)=t1+t2*cos1+stress(4)*sin1
	   ! ta1(2)=-t2*sin1+stress(4)*cos1
    !    ta1(3)=-ta1(2)*dcos(t5)
	   ! ta1(4)=-ta1(2)*dsin(t5)    
    !    t4=dacos(dot_product(dis(1:2),ta1(3:4))/(norm2(dis(1:2))*norm2(ta1(3:4))))
    !    if(t4<t3) then
    !        sfr(2)=t5
    !        sfr(3:6)=ta1
    !        t3=t4
    !    endif
    !
    !enddo 
    
    
    
	
	

endsubroutine

subroutine principal_stress_cal(pstress,stress,ndim)
	implicit none
	integer,intent(in)::ndim
	real(8),intent(in)::stress(6)
	real(8),intent(inout)::pstress(4) 
	real(8)::t1,t2,sita,Pi
    
    Pi=datan(1.d0)*4.0
	
	pstress=0.0d0
	if(ndim==3) then 
		Print *, 'Not available yet. sub=principal_stress_cal'
		stop
	else
		t1=(stress(1)+stress(2))/2.0d0
		t2=(stress(1)-stress(2))/2.0d0
		pstress(1)=t1+(t2**2+stress(4)**2)**0.5		
		pstress(2)=t1-(t2**2+stress(4)**2)**0.5
		Pstress(3)=stress(3)
		!默认z方向为应力为中应力。
        IF(ABS(stress(1)-stress(2))>1E-14) THEN 
		    sita=0.5*atan(2*stress(4)/(stress(1)-stress(2)))
            !sita为莫尔圆主应力平面与x轴的夹角
		    Pstress(4)=sita
            !if(stress(1)>stress(2)) Pstress(4)=sita+Pi/2.0
        ELSE
            Pstress(4)=0 
        ENDIF
			
		
		
		!ref::https://en.wikipedia.org/wiki/Mohr%27s_circle 
		!Mohr Circel sign cenverstion: tension is +, shear stress, clockwise is +
		
		
	endif
	
end subroutine

!calculate stress and strain at gauss points and element nodes for each element needed  
subroutine average_stress_strain_cal(ienum)
	use solverds
	implicit none
	integer,intent(in)::ienum
	integer::i,j,k,n1,n2
	real(8)::un(50)=0.0D0,t1=0.0

	!average at nodes
	!do i=1,enum

		!yield the average variables for each element,weighting with gaussian coefficients.
	n1=element(ienum).ngp+1
	n2=element(ienum).ngp+element(ienum).nnum
	t1=0.0
	do j=1,element(ienum).ngp
		element(ienum).stress(:,n1)=element(ienum).stress(:,n1)+element(ienum).stress(:,j)* &
					ecp(element(ienum).et).weight(j)
		element(ienum).strain(:,n1)=element(ienum).strain(:,n1)+element(ienum).strain(:,j)* &
					ecp(element(ienum).et).weight(j)
		element(ienum).pstrain(:,n1)=element(ienum).pstrain(:,n1)+element(ienum).pstrain(:,j)* &
					ecp(element(ienum).et).weight(j)						
		t1=t1+ecp(element(ienum).et).weight(j)
	end do
	element(ienum).stress(:,n1)=element(ienum).stress(:,n1)/t1		
	element(ienum).strain(:,n1)=element(ienum).strain(:,n1)/t1	
	element(ienum).pstrain(:,n1)=element(ienum).pstrain(:,n1)/t1
	do j=n1+1,n2
		element(ienum).stress(:,j)=element(ienum).stress(:,n1)	
		element(ienum).strain(:,j)=element(ienum).strain(:,n1)
		element(ienum).pstrain(:,j)=element(ienum).pstrain(:,n1)
	end do
	!end do
	
	!call E2N_stress_strain()
	
end subroutine


!expolating from the integration point to nodes
subroutine extrapolation_stress_strain_cal(ienum)
	use solverds
	implicit none
	integer,intent(in)::ienum
	integer::i,j,k=0,n1,n2

	!do i=1,enum
	n1=element(ienum).ngp+1
	n2=element(ienum).ngp+element(ienum).nnum
	!SELECT CASE(ELEMENT(IENUM).EC) 
	
		!CASE(SPG,SPG2D,cax_spg)
	SELECT CASE(ELEMENT(IENUM).ET)
		CASE(CAX3_SPG,CPE3_SPG,CPE4R_SPG,TET4_SPG,CAX4R_SPG,&
			 CAX3,CPE3,CPE4R,TET4,CAX4R,&
			 CAX3_CPL,CPE3_CPL,CPE4R_CPL,TET4_CPL,CAX4R_CPL)
			do concurrent (i=n1:n2)
				IF(ALLOCATED(ELEMENT(IENUM).IGRAD)) ELEMENT(IENUM).IGRAD(:,i)=ELEMENT(IENUM).IGRAD(:,1)
				IF(ALLOCATED(ELEMENT(IENUM).VELOCITY)) ELEMENT(IENUM).VELOCITY(:,i)=ELEMENT(IENUM).VELOCITY(:,1)
				IF(ALLOCATED(ELEMENT(IENUM).KR)) ELEMENT(IENUM).kr(i)=ELEMENT(IENUM).kr(1)
				IF(ALLOCATED(ELEMENT(IENUM).MW)) ELEMENT(IENUM).mw(i)=ELEMENT(IENUM).mw(1)
				IF(ALLOCATED(ELEMENT(IENUM).SFR)) ELEMENT(IENUM).SFR(:,i)=ELEMENT(IENUM).SFR(:,1)
				IF(ALLOCATED(ELEMENT(IENUM).STRESS)) ELEMENT(IENUM).STRESS(:,i)=ELEMENT(IENUM).STRESS(:,1)
				IF(ALLOCATED(ELEMENT(IENUM).STRAIN)) ELEMENT(IENUM).STRAIN(:,i)=ELEMENT(IENUM).STRAIN(:,1)
				IF(ALLOCATED(ELEMENT(IENUM).PSTRAIN)) ELEMENT(IENUM).PSTRAIN(:,i)=ELEMENT(IENUM).PSTRAIN(:,1)
			end do
		CASE(CPE6_SPG,CPE8R_SPG,CPE4_SPG,PRM15_SPG,TET10_SPG,CAX6_SPG,CAX8R_SPG,CAX4_SPG,&
			 CPE6_CPL,CPE8R_CPL,CPE4_CPL,PRM15_CPL,TET10_CPL,CAX6_CPL,CAX8R_CPL,CAX4_CPL,&
			 CPE6,CPE8R,CPE4,PRM15,TET10,CAX6,CAX8R,CAX4,ZT4_SPG)
			do concurrent (i=n1:n2)
				
				do CONCURRENT (j=1:ecp(element(ienum).et).ndim)
					IF(ALLOCATED(ELEMENT(IENUM).IGRAD)) ELEMENT(IENUM).IGRAD(j,i)=dot_product( &
						ELEMENT(IENUM).IGRAD(j,1:element(ienum).ngp), &
						ecp(element(ienum).et).expolating_Lshape(1:element(ienum).ngp,i-n1+1))
					IF(ALLOCATED(ELEMENT(IENUM).VELOCITY)) ELEMENT(IENUM).VELOCITY(j,i)=dot_product( &
						ELEMENT(IENUM).VELOCITY(j,1:element(ienum).ngp), &
						ecp(element(ienum).et).expolating_Lshape(1:element(ienum).ngp,i-n1+1))						
				end do
				
				IF(ALLOCATED(ELEMENT(IENUM).KR)) ELEMENT(IENUM).kr(i)=dot_product( &
					ELEMENT(IENUM).kr(1:element(ienum).ngp), &
					ecp(element(ienum).et).expolating_Lshape(1:element(ienum).ngp,i-n1+1))
				IF(ALLOCATED(ELEMENT(IENUM).MW)) ELEMENT(IENUM).mw(i)=dot_product( &
					ELEMENT(IENUM).mw(1:element(ienum).ngp), &
					ecp(element(ienum).et).expolating_Lshape(1:element(ienum).ngp,i-n1+1))
				do CONCURRENT (j=1:4)
					IF(ALLOCATED(ELEMENT(IENUM).SFR)) ELEMENT(IENUM).SFR(j,i)=dot_product( &
						ELEMENT(IENUM).SFR(j,1:element(ienum).ngp), &
						ecp(element(ienum).et).expolating_Lshape(1:element(ienum).ngp,i-n1+1))	
				enddo
				do CONCURRENT (j=1:ELEMENT(IENUM).ND)
					IF(ALLOCATED(ELEMENT(IENUM).STRESS)) ELEMENT(IENUM).STRESS(j,i)=dot_product( &
						ELEMENT(IENUM).STRESS(j,1:element(ienum).ngp), &
						ecp(element(ienum).et).expolating_Lshape(1:element(ienum).ngp,i-n1+1))
					IF(ALLOCATED(ELEMENT(IENUM).STRAIN)) ELEMENT(IENUM).STRAIN(j,i)=dot_product( &
						ELEMENT(IENUM).STRAIN(j,1:element(ienum).ngp), &
						ecp(element(ienum).et).expolating_Lshape(1:element(ienum).ngp,i-n1+1))
					IF(ALLOCATED(ELEMENT(IENUM).PSTRAIN)) ELEMENT(IENUM).PSTRAIN(j,i)=dot_product( &
						ELEMENT(IENUM).PSTRAIN(j,1:element(ienum).ngp), &
						ecp(element(ienum).et).expolating_Lshape(1:element(ienum).ngp,i-n1+1))					
				end do
				
			end do				
		CASE(CPE15_SPG,CAX15_SPG,&
			 CPE15_CPL,CAX15_CPL,&
			 CPE15,CAX15)
			DO J=1,ecp(element(ienum).et).ndim
				IF(ALLOCATED(ELEMENT(IENUM).IGRAD)) THEN
					CALL I2N_TRI15(element(ienum).IGRAD(J,N1:N2),element(ienum).IGRAD(J,1:ELEMENT(IENUM).NGP),ELEMENT(IENUM).ET)
				ENDIF
				IF(ALLOCATED(ELEMENT(IENUM).VELOCITY)) THEN
					CALL I2N_TRI15(element(ienum).VELOCITY(J,N1:N2),element(ienum).VELOCITY(J,1:ELEMENT(IENUM).NGP),ELEMENT(IENUM).ET)
				ENDIF
			ENDDO		
			IF(ALLOCATED(ELEMENT(IENUM).KR)) THEN
				CALL I2N_TRI15(element(ienum).KR(N1:N2),element(ienum).KR(1:ELEMENT(IENUM).NGP),ELEMENT(IENUM).ET)
			ENDIF	
			IF(ALLOCATED(ELEMENT(IENUM).MW)) THEN
				CALL I2N_TRI15(element(ienum).MW(N1:N2),element(ienum).MW(1:ELEMENT(IENUM).NGP),ELEMENT(IENUM).ET)
			ENDIF
			do j=1,4
				IF(ALLOCATED(ELEMENT(IENUM).SFR)) THEN
					CALL I2N_TRI15(element(ienum).SFR(j,N1:N2),element(ienum).SFR(j,1:ELEMENT(IENUM).NGP),ELEMENT(IENUM).ET)
				ENDIF
			enddo
			do j=1,ELEMENT(IENUM).ND
				IF(ALLOCATED(ELEMENT(IENUM).STRESS)) THEN
					CALL I2N_TRI15(element(ienum).STRESS(J,N1:N2),element(ienum).STRESS(J,1:ELEMENT(IENUM).NGP),ELEMENT(IENUM).ET)
				ENDIF			
				IF(ALLOCATED(ELEMENT(IENUM).STRAIN)) THEN
					CALL I2N_TRI15(element(ienum).STRAIN(J,N1:N2),element(ienum).STRAIN(J,1:ELEMENT(IENUM).NGP),ELEMENT(IENUM).ET)
				ENDIF				
				IF(ALLOCATED(ELEMENT(IENUM).PSTRAIN)) THEN
					CALL I2N_TRI15(element(ienum).PSTRAIN(J,N1:N2),element(ienum).PSTRAIN(J,1:ELEMENT(IENUM).NGP),ELEMENT(IENUM).ET)
				ENDIF				
			end do			
			
			!do CONCURRENT (i=n1:n2)
			!	!k=k+1
			!	K=I-N1+1
			!	do j=1,ecp(element(ienum).et).ndim
			!		IF(ALLOCATED(ELEMENT(IENUM).IGRAD)) element(ienum).IGRAD(J,I)=dot_product( &
			!		matmul(ecp(element(ienum).et).expolating_Lshape,&
			!				 element(ienum).igrad(j,1:element(ienum).ngp)),&
			!						ecp(element(ienum).et).termval(:,k))
			!		IF(ALLOCATED(ELEMENT(IENUM).VELOCITY)) element(ienum).velocity(J,I)=dot_product( &
			!		matmul(ecp(element(ienum).et).expolating_Lshape,&
			!				 element(ienum).velocity(j,1:element(ienum).ngp)),&
			!						ecp(element(ienum).et).termval(:,k))							
			!	end do
			!	IF(ALLOCATED(ELEMENT(IENUM).KR)) element(ienum).kr(I)=dot_product( &
			!	matmul(ecp(element(ienum).et).expolating_Lshape,&
			!			 element(ienum).kr(1:element(ienum).ngp)),&
			!					ecp(element(ienum).et).termval(:,k))
			!	IF(ALLOCATED(ELEMENT(IENUM).MW)) element(ienum).mw(I)=dot_product( &
			!	matmul(ecp(element(ienum).et).expolating_Lshape,&
			!			 element(ienum).mw(1:element(ienum).ngp)),&
			!					ecp(element(ienum).et).termval(:,k))	
			!	
			!	do CONCURRENT (j=1:6)
			!		IF(ALLOCATED(ELEMENT(IENUM).STRESS)) element(ienum).STRESS(J,I)=dot_product( &
			!		matmul(ecp(element(ienum).et).expolating_Lshape,&
			!				 element(ienum).STRESS(j,1:element(ienum).ngp)),&
			!						ecp(element(ienum).et).termval(:,k))
			!		IF(ALLOCATED(ELEMENT(IENUM).STRAIN)) element(ienum).STRAIN(J,I)=dot_product( &
			!		matmul(ecp(element(ienum).et).expolating_Lshape,&
			!				 element(ienum).STRAIN(j,1:element(ienum).ngp)),&
			!						ecp(element(ienum).et).termval(:,k))
			!		IF(ALLOCATED(ELEMENT(IENUM).PSTRAIN)) element(ienum).PSTRAIN(J,I)=dot_product( &
			!		matmul(ecp(element(ienum).et).expolating_Lshape,&
			!				 element(ienum).PSTRAIN(j,1:element(ienum).ngp)),&
			!						ecp(element(ienum).et).termval(:,k))					
			!	end do
			!	
			!	
			!end do
			
		CASE(prm6_spg,PRM6,PRM6_CPL)
			DO CONCURRENT (I=1:2)
				N1=3*I
				N2=3*I+2
				DO CONCURRENT (J=1:ecp(element(ienum).et).ndim)
					IF(ALLOCATED(ELEMENT(IENUM).IGRAD)) element(ienum).IGRAD(J,N1:N2)=dot_product(ELEMENT(IENUM).IGRAD(j,1:element(ienum).ngp), &
								ecp(element(ienum).et).expolating_Lshape(1:element(ienum).ngp,I))
					IF(ALLOCATED(ELEMENT(IENUM).VELOCITY)) element(ienum).VELOCITY(J,N1:N2)=dot_product(ELEMENT(IENUM).VELOCITY(j,1:element(ienum).ngp), &
								ecp(element(ienum).et).expolating_Lshape(1:element(ienum).ngp,I))
				END DO
				IF(ALLOCATED(ELEMENT(IENUM).KR)) element(ienum).kr(N1:N2)=dot_product(ELEMENT(IENUM).kr(1:element(ienum).ngp), &
							ecp(element(ienum).et).expolating_Lshape(1:element(ienum).ngp,I))
				IF(ALLOCATED(ELEMENT(IENUM).MW)) element(ienum).mw(N1:N2)=dot_product(ELEMENT(IENUM).mw(1:element(ienum).ngp), &
							ecp(element(ienum).et).expolating_Lshape(1:element(ienum).ngp,I))
				DO CONCURRENT (J=1:4) 
					IF(ALLOCATED(ELEMENT(IENUM).SFR)) element(ienum).SFR(j,N1:N2)=dot_product(ELEMENT(IENUM).SFR(j,1:element(ienum).ngp), &
							ecp(element(ienum).et).expolating_Lshape(1:element(ienum).ngp,I))
				enddo
				DO CONCURRENT (J=1:ELEMENT(IENUM).ND)
					IF(ALLOCATED(ELEMENT(IENUM).STRESS)) element(ienum).STRESS(J,N1:N2)=dot_product(ELEMENT(IENUM).STRESS(j,1:element(ienum).ngp), &
								ecp(element(ienum).et).expolating_Lshape(1:element(ienum).ngp,I))
					IF(ALLOCATED(ELEMENT(IENUM).STRAIN)) element(ienum).STRAIN(J,N1:N2)=dot_product(ELEMENT(IENUM).STRAIN(j,1:element(ienum).ngp), &
								ecp(element(ienum).et).expolating_Lshape(1:element(ienum).ngp,I))
					IF(ALLOCATED(ELEMENT(IENUM).PSTRAIN)) element(ienum).STRAIN(J,N1:N2)=dot_product(ELEMENT(IENUM).PSTRAIN(j,1:element(ienum).ngp), &
								ecp(element(ienum).et).expolating_Lshape(1:element(ienum).ngp,I))
				END DO				
			END DO

		CASE DEFAULT
			PRINT *, 'NO SUCH A ELEMENT TYPE. SUB extrapolation_stress_strain_cal'
			STOP					
	END SELECT
			
		!CASE DEFAULT
			
		!	PRINT *, 'POSTPROCEDURE IN SOLID ELEMENT CLASS IS NOT COMPLETED.'
			!STOP
	!END SELECT
	
	
end subroutine

SUBROUTINE I2N_TRI15(VIN,VIGP,ET) !直接平移
	USE SOLVERLIB
	IMPLICIT NONE
	INTEGER,INTENT(IN)::ET
	REAL(KIND=DPN),INTENT(IN)::VIGP(ECP(ET).NGP)
	REAL(KIND=DPN),INTENT(OUT)::VIN(ECP(ET).NNUM)
	
	
    VIN=(/VIGP(1:3),(VIGP(7)+VIGP(8))/2.,(VIGP(9)+VIGP(10))/2.,(VIGP(11)+VIGP(12))/2., &
        (VIGP(1)+VIGP(7))/2.,(VIGP(8)+VIGP(2))/2.,(VIGP(2)+VIGP(9))/2., &
        (VIGP(10)+VIGP(3))/2.,(VIGP(3)+VIGP(11))/2.,(VIGP(1)+VIGP(12))/2., &
        VIGP(5:6),VIGP(4)/)
    
	!VIN(1:ECP(ET).NGP)=MATMUL(ECP(ET).expolating_Lshape,& 
	!					VIGP-MATMUL(TRANSPOSE(ecp(et).Lshape(13:15,1:ECP(ET).NGP)),VIN(13:15)))
	
ENDSUBROUTINE
    


!according the location of integration point, shift it to the closest element node
subroutine shift_stress_strain_cal(ienum)
	use solverds
	implicit none
	integer,intent(in)::ienum
	integer::i,j,n1,n2

	!do i=1,enum
	n1=element(ienum).ngp+1
	n2=element(ienum).ngp+element(ienum).nnum
	select case(element(ienum).et)
		case(cpe3,cps3,cax3,cpe4R,cps4R,cax4R)
			do j=n1,n2
				element(ienum).stress(:,j)=element(ienum).stress(:,1)
				element(ienum).strain(:,j)=element(ienum).strain(:,1)
				element(ienum).pstrain(:,j)=element(ienum).pstrain(:,1)
			end do
		case(cpe6,cps6,cax6)
			element(ienum).stress(:,7:9)=element(ienum).stress(:,1:3)
			element(ienum).strain(:,7:9)=element(ienum).strain(:,1:3)
			element(ienum).pstrain(:,7:9)=element(ienum).pstrain(:,1:3)
			element(ienum).stress(:,1)=element(ienum).stress(:,7)+element(ienum).stress(:,9) &
				-element(ienum).stress(:,8)
			element(ienum).stress(:,2)=element(ienum).stress(:,7)+element(ienum).stress(:,8) &
				-element(ienum).stress(:,9)
			element(ienum).stress(:,3)=element(ienum).stress(:,8)+element(ienum).stress(:,9) &
				-element(ienum).stress(:,7)
				
			element(ienum).strain(:,4)=element(ienum).strain(:,7)+element(ienum).strain(:,9) &
				-element(ienum).strain(:,8)
			element(ienum).strain(:,5)=element(ienum).strain(:,7)+element(ienum).strain(:,8) &
				-element(ienum).strain(:,9)
			element(ienum).strain(:,6)=element(ienum).strain(:,8)+element(ienum).strain(:,9) &
				-element(ienum).strain(:,7)
				
			element(ienum).pstrain(:,4)=element(ienum).pstrain(:,7)+element(ienum).pstrain(:,9) &
				-element(ienum).pstrain(:,8)
			element(ienum).pstrain(:,5)=element(ienum).pstrain(:,7)+element(ienum).pstrain(:,8) &
				-element(ienum).pstrain(:,9)
			element(ienum).pstrain(:,6)=element(ienum).pstrain(:,8)+element(ienum).pstrain(:,9) &
				-element(ienum).pstrain(:,7)
				
		case(cpe15,cps15,cax15)
			element(ienum).stress(:,13:15)=element(ienum).stress(:,1:3)
			element(ienum).strain(:,13:15)=element(ienum).strain(:,1:3)
			element(ienum).pstrain(:,13:15)=element(ienum).pstrain(:,1:3)
			
			element(ienum).stress(:,12+7:12+12)=element(ienum).stress(:,7:12)
			element(ienum).strain(:,12+7:12+12)=element(ienum).strain(:,7:12)
			element(ienum).pstrain(:,12+7:12+12)=element(ienum).pstrain(:,7:12)
			
			element(ienum).stress(:,12+13:12+14)=element(ienum).stress(:,5:6)
			element(ienum).strain(:,12+13:12+14)=element(ienum).strain(:,5:6)
			element(ienum).pstrain(:,12+13:12+14)=element(ienum).pstrain(:,5:6)
            
			element(ienum).stress(:,12+15)=element(ienum).stress(:,4)
			element(ienum).strain(:,12+15)=element(ienum).strain(:,4)
			element(ienum).pstrain(:,12+15)=element(ienum).pstrain(:,4)
            
			element(ienum).stress(:,12+4)=(element(ienum).stress(:,7)+element(ienum).stress(:,8))/2
			element(ienum).stress(:,12+5)=(element(ienum).stress(:,9)+element(ienum).stress(:,10))/2
			element(ienum).stress(:,12+6)=(element(ienum).stress(:,11)+element(ienum).stress(:,12))/2
			
			element(ienum).strain(:,12+4)=(element(ienum).strain(:,7)+element(ienum).strain(:,8))/2
			element(ienum).strain(:,12+5)=(element(ienum).strain(:,9)+element(ienum).strain(:,10))/2
			element(ienum).strain(:,12+6)=(element(ienum).strain(:,11)+element(ienum).strain(:,12))/2
			
			element(ienum).pstrain(:,12+4)=(element(ienum).pstrain(:,7)+element(ienum).pstrain(:,8))/2
			element(ienum).pstrain(:,12+5)=(element(ienum).pstrain(:,9)+element(ienum).pstrain(:,10))/2
			element(ienum).pstrain(:,12+6)=(element(ienum).pstrain(:,11)+element(ienum).pstrain(:,12))/2
			
		case(cpe8,cps8,cax8)
			element(ienum).stress(:,9+1:9+8)=element(ienum).stress(:,1:8)
			element(ienum).strain(:,9+1:9+8)=element(ienum).strain(:,1:8)
			element(ienum).pstrain(:,9+1:9+8)=element(ienum).pstrain(:,1:8)				
		case(cpe8r,cps8r,cax8r)
			element(ienum).stress(:,4+1:4+4)=element(ienum).stress(:,1:4)
			element(ienum).strain(:,4+1:4+4)=element(ienum).strain(:,1:4)
			element(ienum).pstrain(:,4+1:4+4)=element(ienum).pstrain(:,1:4)
	
			element(ienum).stress(:,4+5)=(element(ienum).stress(:,1)+element(ienum).stress(:,2))/2
			element(ienum).stress(:,4+6)=(element(ienum).stress(:,2)+element(ienum).stress(:,3))/2
			element(ienum).stress(:,4+7)=(element(ienum).stress(:,3)+element(ienum).stress(:,4))/2
			element(ienum).stress(:,4+8)=(element(ienum).stress(:,4)+element(ienum).stress(:,1))/2
			
			element(ienum).strain(:,4+5)=(element(ienum).strain(:,1)+element(ienum).strain(:,2))/2
			element(ienum).strain(:,4+6)=(element(ienum).strain(:,2)+element(ienum).strain(:,3))/2
			element(ienum).strain(:,4+7)=(element(ienum).strain(:,3)+element(ienum).strain(:,4))/2
			element(ienum).strain(:,4+8)=(element(ienum).strain(:,4)+element(ienum).strain(:,1))/2
	
			element(ienum).pstrain(:,4+5)=(element(ienum).pstrain(:,1)+element(ienum).pstrain(:,2))/2
			element(ienum).pstrain(:,4+6)=(element(ienum).pstrain(:,2)+element(ienum).pstrain(:,3))/2
			element(ienum).pstrain(:,4+7)=(element(ienum).pstrain(:,3)+element(ienum).pstrain(:,4))/2
			element(ienum).pstrain(:,4+8)=(element(ienum).pstrain(:,4)+element(ienum).pstrain(:,1))/2				
	end select
	!end do
	
	!call E2N_stress_strain()
	
end subroutine

!subroutine spr_stress_strain_cal()
!	include "link_fnl_shared.h"
!	use solverds
!	!use dfimsl
!	implicit none	
!	integer::i,j,k,k1,iel,nitem,nc,ixy,jxy,ivp,jvp,is
!	real(8)::b(30,20),xy(2,1),vp(30,1),x(30,20),res(30)
!	
!	do i=1,nnum
!		
!		if(SPRLIST(i).ispatch) then
!			nitem=sprlist(i).item
!			nc=0
!			b=0.0
!			do j=1,sprlist(i).nelist
!				iel=sprlist(i).elist(j)
!				do k=1,element(iel).ngp
!					nc=nc+1
!					do k1=1,6
!						b(1:nitem,k1)=b(1:nitem,k1)+sprlist(i).polynomial(:,nc)*element(iel).stress(k1,k)
!						b(1:nitem,6+k1)=b(1:nitem,6+k1)+sprlist(i).polynomial(:,nc)*element(iel).strain(k1,k)
!						b(1:nitem,12+k1)=b(1:nitem,12+k1)+sprlist(i).polynomial(:,nc)*element(iel).pstrain(k1,k)
!					end do
!				end do
!			end do
!!			write(10,'(i)') i
!!			write(10,999) (b(j,2),j=1,nitem)
!!			write(10,'(a)') ''
!!			write(10,999) ((sprlist(i).invA(j,k),j=1,nitem),k=1,nitem)
!!			x=0.0			
!			do j=1,18
!				res=0.0
!				!b(1:nitem,j)=matmul(sprlist(i).invA(1:nitem,1:nitem),b(1:nitem,j))
!				CALL DLFSSF (nitem, sprlist(i).A, nitem,sprlist(i).ipvt, b(1:nitem,j), b(1:nitem,j))
!!				CALL DLFISF (nitem, sprlist(i).A, nitem, sprlist(i).invA, nitem, &
!!					sprlist(i).ipvt, b(1:nitem,j), X(1:nitem,j),RES(1:nitem))
!				!CALL DLFSDS (nitem, sprlist(i).invA, nitem, b(1:nitem,j), b(1:nitem,j))
!			end do
!!			b=x
!			!write(10,999) (b(j,2),j=1,nitem)
!			
!			do j=1,sprlist(i).nelist
!				do k=1,element(sprlist(i).elist(j)).nnum
!					if((.not.sprlist(element(sprlist(i).elist(j)).node(k)).ispatch) &
!							.or.(element(sprlist(i).elist(j)).node(k)==i)) then
!						xy(1,1)=(node(element(sprlist(i).elist(j)).node(k)).coord(1)- &
!							node(i).coord(1))/sprlist(i).dxmax
!						xy(2,1)=(node(element(sprlist(i).elist(j)).node(k)).coord(2)- &
!							node(i).coord(2))/sprlist(i).dymax						
!						ixy=2
!						jxy=1
!						ivp=nitem
!						jvp=1
!						vp=0.0
!						call vpolynomial(xy(1:ixy,1:jxy),ixy,jxy,sprlist(i).p,vp(1:ivp,1:jvp),ivp,jvp)
!						is=element(sprlist(i).elist(j)).ngp+k
!						!Attentioin. the order should be consistent.
!						do k1=1,6
!							if(.not.sprlist(element(sprlist(i).elist(j)).node(k)).ispatch) then
!								element(sprlist(i).elist(j)).stress(k1,is)=element(sprlist(i).elist(j)).stress(k1,is)+ &
!									dot_product(b(1:nitem,k1),vp(1:nitem,1))/element(sprlist(i).elist(j)).nspr
!								element(sprlist(i).elist(j)).strain(k1,is)=element(sprlist(i).elist(j)).strain(k1,is)+ &
!									dot_product(b(1:nitem,k1+6),vp(1:nitem,1))/element(sprlist(i).elist(j)).nspr
!								element(sprlist(i).elist(j)).pstrain(k1,is)=element(sprlist(i).elist(j)).pstrain(k1,is)+ &
!									dot_product(b(1:nitem,k1+12),vp(1:nitem,1))/element(sprlist(i).elist(j)).nspr
!							else
!								element(sprlist(i).elist(j)).stress(k1,is)=dot_product(b(1:nitem,k1),vp(1:nitem,1))
!								element(sprlist(i).elist(j)).strain(k1,is)=dot_product(b(1:nitem,k1+6),vp(1:nitem,1))
!								element(sprlist(i).elist(j)).pstrain(k1,is)=dot_product(b(1:nitem,k1+12),vp(1:nitem,1))					
!							end if
!						end do
!					end if
!				end do
!			end do			
!		end if
!	end do	
!
!	call E2N_stress_strain()
!999 format(<sprlist(i).item>e15.7)	
!end subroutine

!according average the nodal stress in elements.
subroutine E2N_stress_strain(ISTEP)
	use solverds
	implicit none
	integer::i,j,k,n1,n2,ISTEP
	real(8)::un(50)=0.0D0,vangle(15)=0.0,C1,PHI1,dis1(3)

	!clear zero
	do i=1,nnum
		if(allocated(node(i).stress))	node(i).stress=0.0d0
		if(allocated(node(i).strain)) node(i).strain=0.0d0
		if(allocated(node(i).pstrain))	node(i).pstrain=0.0d0
		if(allocated(node(i).igrad)) node(i).igrad=0.0d0
		if(allocated(node(i).velocity)) node(i).velocity=0.0d0
		if(allocated(node(i).sfr)) then
            node(i).sfr=0.0d0
            node(i).sfr(1)=-1.0d6
        endif
		node(i).q=0.0d0
		node(i).kr=0.0d0
		node(i).mw=0.0d0
        !node(i).sfr(1)=-1.0d6
	end do
	
	!averaged simplily at nodes
	do i=1,enum
			
		!vangle=pi()
		!select case(element(i).et)
		!	case(cps3,cpe3,cax3,cpe6,cps6,cax6,cpe15,cps15,cax15, &
		!		  cpe3_spg,cpe6_spg,cax6_spg,cpe15_spg,cax15_spg)
		!		vangle(1:3)=element(i).angle(1:3)
		!		vangle(13:15)=2*vangle(13:15)

		!	case(cpe4,cps4,cax4,cpe4r,cps4r,cax4r,cpe8,cps8,cax8,cpe8r,cps8r,cax8r, &
		!		  cpe4_spg,cax4_spg,cpe4r_spg,cax4r_spg,cpe8_spg,cpe8r_spg,cax8r_spg)
		!		vangle(1:4)=element(i).angle(1:4)
		!end select
		n1=element(i).ngp
		
		select case(element(i).ec)
			case(C3D,CPE,CPS,CAX,CPL,CAX_CPL)
				do j=1,element(i).nnum
					n2=node(element(i).node(j)).nelist-node(element(i).node(j)).nelist_SPG
					node(element(i).node(j)).stress=node(element(i).node(j)).stress &
					+element(i).stress(:,n1+j)/n2
					node(element(i).node(j)).strain=node(element(i).node(j)).strain &
					+element(i).strain(:,n1+j)/n2
					node(element(i).node(j)).pstrain=node(element(i).node(j)).pstrain &
					+element(i).pstrain(:,n1+j)/n2
					!get the mat for later nodal sfr cal. 
                    
					if(element(i).sfr(1,n1+j)>node(element(i).node(j)).sfr(1)) then					
						node(element(i).node(j)).sfr(1)=element(i).sfr(1,n1+j)
						node(element(i).node(j)).sfr(2)=element(i).mat !borrow 
					endif
				
				end do
			case(spg2d,spg,cax_spg)
				do j=1,element(i).nnum
					n2=node(element(i).node(j)).nelist_SPG
					node(element(i).node(j)).igrad=node(element(i).node(j)).igrad &
					+element(i).igrad(:,n1+j)/n2
					node(element(i).node(j)).velocity=node(element(i).node(j)).velocity &
					+element(i).velocity(:,n1+j)/n2
					node(element(i).node(j)).kr=node(element(i).node(j)).kr &
					+element(i).kr(n1+j)/n2	
					node(element(i).node(j)).mw=node(element(i).node(j)).mw &
					+element(i).mw(n1+j)/n2						
					node(element(i).node(j)).q=node(element(i).node(j)).q+element(i).flux(j)
                end do
            case(stru,spring,soilspring)
                
			case default
				print *, 'Not Completed in E2N_stress_strain. EC=', element(i).ec
		end select
		
	end do

	! generalized shear stress and strain
	do i=1,nnum
		!shear strain
		if(allocated(node(i).stress)) then
			node(i).mises=0.50*((node(i).stress(1)-node(i).stress(2))**2 &
									+(node(i).stress(2)-node(i).stress(3))**2 &
									+(node(i).stress(1)-node(i).stress(3))**2 &
									+6*(node(i).stress(4)**2+node(i).stress(5)**2 &
									+node(i).stress(6)**2))
			node(i).mises=node(i).mises**0.5
			!节点的材料假定为破坏比最大的单元材料，以模拟成层土中的软弱夹层
			C1=material(node(i).sfr(2)).GET(3,ISTEP)
			Phi1=material(node(i).sfr(2)).GET(4,ISTEP)
            !dis1(1:ndimension)=Tdisp(node(i).dof(1:ndimension))
			NODE(I).SFR(7)=C1;NODE(I).SFR(8)=PHI1;
			call stress_in_failure_surface(node(i).sfr,node(i).stress,2,C1,Phi1,solver_control.slidedirection,node(i).coord(ndimension))
			
		end if
		!total generalized shear stain
		if(allocated(node(i).strain)) then
			node(i).eeq=1.0/6.0*((node(i).strain(1)-node(i).strain(2))**2 &
									+(node(i).strain(2)-node(i).strain(3))**2 &
									+(node(i).strain(1)-node(i).strain(3))**2 &
									+6*(node(i).strain(4)**2+node(i).strain(5)**2 &
									+node(i).strain(6)**2))
			node(i).eeq=2.0/3**0.5*node(i).eeq**0.5
		end if
		!total generalized plastic shear stain
		if(allocated(node(i).pstrain)) then
			node(i).peeq=1.0/6.0*((node(i).pstrain(1)-node(i).pstrain(2))**2 &
									+(node(i).pstrain(2)-node(i).pstrain(3))**2 &
									+(node(i).pstrain(1)-node(i).pstrain(3))**2 &
									+6*(node(i).pstrain(4)**2+node(i).pstrain(5)**2 &
									+node(i).pstrain(6)**2))
			node(i).peeq=2.0/3**0.5*node(i).peeq**0.5	
		end if
		

		
		
	end do
	
end subroutine

! form the element list sharing the same vertex(mid-nodes are excluded)
subroutine spr_elist()
	use solverds
	implicit none
	integer::i,j
	integer::nnum1,p1,nsc1
	real(8)::t1
	
	allocate(sprlist(nnum)) !assume that there  are 10 elements at most sharing the same vertex.
	do i=1,enum
		select case(element(i).et)
			case(cps3,cpe3,cax3)
				nnum1=3
				nsc1=1
				p1=1
			case(cps4r,cpe4r,cax4r)
				nnum1=4
				nsc1=1
				p1=1
			case(cps4,cpe4,cax4)
				nnum1=4
				nsc1=4
				p1=1
			case(cps6,cpe6,cax6)
				nnum1=3
				nsc1=3
				p1=2
			case(cpe15,cps15,cax15)
				nnum1=3
				nsc1=12
				p1=4	
			case(cps8r,cpe8r,cax8r)
				nnum1=4
				nsc1=4
				p1=2
			case(cps8,cpe8,cax8)
				nnum1=4
				nsc1=9
				p1=2
		end select
		do j=1,nnum1
			sprlist(element(i).node(j)).nelist=sprlist(element(i).node(j)).nelist+1
			SPRLIST(element(i).node(j)).elist(sprlist(element(i).node(j)).nelist)=i
			sprlist(element(i).node(j)).nsc=sprlist(element(i).node(j)).nsc+nsc1
			if(p1>sprlist(element(i).node(j)).p) sprlist(element(i).node(j)).p=p1
			sprlist(element(i).node(j)).ispatch=.true.
		end do		
	end do
	
	do i=1,nnum
		if(sprlist(i).ispatch) then
			select case(sprlist(i).p) 
				case(1)
					if(sprlist(i).nsc<3) sprlist(i).ispatch=.false.
					if(sprlist(i).nsc>6) sprlist(i).p=2
				case(2)
					if(sprlist(i).nsc<6) sprlist(i).ispatch=.false.
					if(sprlist(i).nsc>=12) sprlist(i).p=3
				case(3)
					if(sprlist(i).nsc<10) sprlist(i).ispatch=.false.
					if(sprlist(i).nsc>=20) sprlist(i).p=4
				case(4)
					if(sprlist(i).nsc<15) sprlist(i).ispatch=.false.
					if(sprlist(i).nsc>=30) sprlist(i).p=5
			end select
			if(sprlist(i).ispatch) then
				do j=1,sprlist(i).nelist
					element(sprlist(i).elist(j)).nspr=element(sprlist(i).elist(j)).nspr+1
					t1=maxval(node(element(sprlist(i).elist(j)).node).coord(1))- &
					minval(node(element(sprlist(i).elist(j)).node).coord(1))
					if(t1>sprlist(i).dxmax) sprlist(i).dxmax=t1
					t1=maxval(node(element(sprlist(i).elist(j)).node).coord(2))- &
					minval(node(element(sprlist(i).elist(j)).node).coord(2))
					if(t1>sprlist(i).dymax) sprlist(i).dymax=t1
				end do				
			end if 
		end if
	end do
end subroutine



!calculate the interpolating polynomials
subroutine vpolynomial(xy,ixy,jxy,order,p,ip,jp)
	implicit none
	integer::ixy,jxy,ip,jp,order,i,j
	real(8)::xy(ixy,jxy),p(ip,jp),dx,dy
	
	if(order>5) then
		print *, "the order of interpolating polynomials is assumed to be less or equel to 5." 
		stop
	end if
	p=0.0	
	do i=1,jxy
		dx=xy(1,i)
		dy=xy(2,i)
		if(order>=1) then
			p(1,i)=1.0
			p(2,i)=dx
			p(3,i)=dy
		end if
		if(order>=2) then
			p(4,i)=dx**2
			p(5,i)=dx*dy
			p(6,i)=dy**2
		end if
		if(order>=3) then
			p(7,i)=dx**3
			p(8,i)=dx**2*dy
			p(9,i)=dx*dy**2
			p(10,i)=dy**3
		end if
		if(order>=4) then
			p(11,i)=dx**4
			p(12,i)=dx**3*dy
			p(13,i)=dx**2*dy**2
			p(14,i)=dx*dy**3
			p(15,i)=dy**4
		end if
		if(order>=5) then
			p(16,i)=dx**5
			p(17,i)=dx**4*dy
			p(18,i)=dx**3*dy**2
			p(19,i)=dx**2*dy**3
			p(20,i)=dx*dy**4
			p(21,i)=dy**5
		end if		
	end do
	
end subroutine



!use SPR method to calculate the polynomials for each patch.
!subroutine spr_polynomial_cal()
!	use solverds
!	use operation_i
!	implicit none
!	integer::i,j,k,k1,nc,iel,item,ixy
!	real(8)::xy(3,120),w(100),rcond
!	
!	do i=1,nnum
! 		xy=0.0
!		if(SPRLIST(i).ispatch) then
!			nc=0
!			do j=1,sprlist(i).nelist
!				iel=sprlist(i).elist(j)
!				do k=1,element(iel).ngp
!					nc=nc+1
!					xy(1,nc)=(element(iel).xygp(1,k)-node(i).coord(1))/SPRLIST(i).dxmax
!					xy(2,nc)=(element(iel).xygp(2,k)-node(i).coord(2))/SPRLIST(i).dymax
!				end do
!			end do
!			select case(sprlist(i).p)
!				case(1)
!					sprlist(i).item=3
!				case(2)
!					sprlist(i).item=6
!				case(3)
!					sprlist(i).item=10
!				case(4)
!					sprlist(i).item=15
!				case(5)
!					sprlist(i).item=21
!			end select
!			allocate(sprlist(i).polynomial(sprlist(i).item,sprlist(i).nsc),&
!				sprlist(i).ipvt(sprlist(i).item),&
!				sprlist(i).A(sprlist(i).item,sprlist(i).item))
!			ixy=2
!			sprlist(i).polynomial=0.0
!			sprlist(i).A=0.0
!			sprlist(i).ipvt=0.0
!			call vpolynomial(xy(1:ixy,1:nc),ixy,nc,sprlist(i).p,sprlist(i).polynomial(1:sprlist(i).item,& 
!							1:sprlist(i).nsc),sprlist(i).item,sprlist(i).nsc)
!			do j=1,sprlist(i).nsc
!				do k=1,sprlist(i).item
!					do k1=1,sprlist(i).item
!						sprlist(i).A(k,k1)=sprlist(i).A(k,k1)+ &
!							sprlist(i).polynomial(k,j)*sprlist(i).polynomial(k1,j)
!					end do
!				end do
!			end do
!			rcond=0.0
!			CALL DLFCSF(sprlist(i).item, sprlist(i).A, sprlist(i).item, sprlist(i).A, sprlist(i).item, &
!				sprlist(i).ipvt,rcond) 
!			if(abs(rcond)<1e-15) then
!				sprlist(i).ispatch=.false.
!				do j=1,sprlist(i).nelist
!					element(sprlist(i).elist(j)).nspr=element(sprlist(i).elist(j)).nspr-1
!				end do
!			end if
!!			CALL DLFTDS (sprlist(i).item, sprlist(i).invA, sprlist(i).item, sprlist(i).invA, sprlist(i).item)
!!			CALL DLINDS (sprlist(i).item, sprlist(i).invA, sprlist(i).item, sprlist(i).invA, sprlist(i).item)
!!			w=0.0
!!			CALL MB01CD(sprlist(i).invA,sprlist(i).item,sprlist(i).item,sprlist(i).item,W(1:sprlist(i).item))
!!			sprlist(i).invA=.i.(sprlist(i).invA(1:sprlist(i).item,1:sprlist(i).item))
!!			write(10,999 ) (sprlist(i).invA(:,j),j=1,sprlist(i).item)
!!			write(10,'(a)') ''
!		end if
!	end do
!	
!999 format(<sprlist(i).item>e15.7)
!	
!end subroutine




