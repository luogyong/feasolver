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
				do CONCURRENT (j=1:6)
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
			do j=1,6
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
				DO CONCURRENT (J=1:6) 
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

subroutine spr_stress_strain_cal(ISTEP,ISUBST)
	
	use SPR_RECOVERY
	USE solverds
	implicit none
    INTEGER,INTENT(IN)::ISTEP,ISUBST
    INTEGER::I,J,K,E1,N1,N2,NPATCH1
    REAL(DPN)::VAL(100)
    
    IF(ISUBST==1) CALL spr_initialize()
    
        
    
    DO I=1,NNUM
        IF(ALLOCATED(NODE(I).STRESS)) NODE(I).STRESS=0.D0
        IF(ALLOCATED(NODE(I).STRAIN)) NODE(I).STRAIN=0.D0
        IF(ALLOCATED(NODE(I).PSTRAIN)) NODE(I).PSTRAIN=0.D0
        IF(ALLOCATED(NODE(I).IGRAD)) NODE(I).IGRAD=0.D0
        IF(ALLOCATED(NODE(I).VELOCITY)) NODE(I).VELOCITY=0.D0
        NODE(I).KR=0.D0;NODE(I).MW=0.D0
    ENDDO
    
    DO I=1,NNUM
        IF(.NOT.SPRLIST(I).ISPATCH) CYCLE
        
        DO K=1,SPRLIST(I).NNODE
            N2=SPRLIST(I).NODE(K)
            NPATCH1=SPRLIST(N2).NPATCH
            IF(ALLOCATED(NODE(N2).STRESS)) THEN           
                DO J=1,2*NDIMENSION            
                    NODE(N2).STRESS(J)=NODE(N2).STRESS(J)+SPRLIST(I).GETNODALDERIVATIVE(NODE(N2).COORD(1:NDIMENSION),SPRLIST(I).GETGPVALUE(J+SXX-1))/NPATCH1            
                ENDDO
            ENDIF
            IF(ALLOCATED(NODE(N2).STRAIN)) THEN           
                DO J=1,2*NDIMENSION            
                    NODE(N2).STRAIN(J)=NODE(N2).STRAIN(J)+SPRLIST(I).GETNODALDERIVATIVE(NODE(N2).COORD(1:NDIMENSION),SPRLIST(I).GETGPVALUE(J+EXX-1))/NPATCH1            
                ENDDO
            ENDIF        
            IF(ALLOCATED(NODE(N2).PSTRAIN)) THEN           
                DO J=1,2*NDIMENSION            
                    NODE(N2).PSTRAIN(J)=NODE(N2).PSTRAIN(J)+SPRLIST(I).GETNODALDERIVATIVE(NODE(N2).COORD(1:NDIMENSION),SPRLIST(I).GETGPVALUE(J+PEXX-1))/NPATCH1            
                ENDDO
            ENDIF        
            IF(ALLOCATED(NODE(N2).IGRAD)) THEN           
                DO J=1,NDIMENSION            
                    NODE(N2).IGRAD(J)=NODE(N2).IGRAD(J)+SPRLIST(I).GETNODALDERIVATIVE(NODE(N2).COORD(1:NDIMENSION),SPRLIST(I).GETGPVALUE(J+GRADX-1))/NPATCH1            
                ENDDO
            ENDIF 
            IF(ALLOCATED(NODE(N2).VELOCITY)) THEN           
                DO J=1,NDIMENSION            
                    NODE(N2).VELOCITY(J)=NODE(N2).VELOCITY(J)+SPRLIST(I).GETNODALDERIVATIVE(NODE(N2).COORD(1:NDIMENSION),SPRLIST(I).GETGPVALUE(J+VX-1))/NPATCH1            
                ENDDO
            ENDIF
            IF(OUTVAR(KR_SPG).VALUE>0) THEN
                NODE(N2).KR=NODE(N2).KR+SPRLIST(I).GETNODALDERIVATIVE(NODE(N2).COORD(1:NDIMENSION),SPRLIST(I).GETGPVALUE(KR_SPG))/NPATCH1        
            ENDIF
            IF(OUTVAR(MW_SPG).VALUE>0) THEN
                NODE(N2).MW=NODE(N2).MW+SPRLIST(I).GETNODALDERIVATIVE(NODE(N2).COORD(1:NDIMENSION),SPRLIST(I).GETGPVALUE(MW_SPG))/NPATCH1        
            ENDIF
        
        ENDDO
    ENDDO
			

!
!	call E2N_stress_strain()
!999 format(<sprlist(i).item>e15.7)	
end subroutine

!according average the nodal stress in elements.
subroutine E2N_stress_strain(ISTEP)
	use solverds
	implicit none
	integer::i,j,k,n1,n2,ISTEP
	real(8)::un(50)=0.0D0,vangle(15)=0.0,C1,PHI1,dis1(3)

	!clear zero
	do i=1,nnum
        IF(NODE(I).ISACTIVE==0) CYCLE
        IF(solver_control.i2ncal/=SPR) THEN
		    if(allocated(node(i).stress))	node(i).stress=0.0d0
		    if(allocated(node(i).strain)) node(i).strain=0.0d0
		    if(allocated(node(i).pstrain))	node(i).pstrain=0.0d0
		    if(allocated(node(i).igrad)) node(i).igrad=0.0d0
		    if(allocated(node(i).velocity)) node(i).velocity=0.0d0
		    node(i).kr=0.0d0
		    node(i).mw=0.0d0
        ENDIF
        if(allocated(node(i).sfr)) then
            node(i).sfr=0.0d0
            node(i).sfr(1)=-1.0d6
        endif
		node(i).q=0.0d0
        !node(i).sfr(1)=-1.0d6
	end do
	
	!averaged simplily at nodes
	do i=1,enum
		IF(ELEMENT(I).ISACTIVE==0) CYCLE	
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
                    IF(solver_control.i2ncal/=SPR) THEN
					    n2=node(element(i).node(j)).nelist-node(element(i).node(j)).nelist_SPG
					    node(element(i).node(j)).stress=node(element(i).node(j)).stress &
					    +element(i).stress(:,n1+j)/n2
					    node(element(i).node(j)).strain=node(element(i).node(j)).strain &
					    +element(i).strain(:,n1+j)/n2
					    node(element(i).node(j)).pstrain=node(element(i).node(j)).pstrain &
					    +element(i).pstrain(:,n1+j)/n2
                    ENDIF
                    
                    
                    
					!get the mat for later nodal sfr cal. 
                    IF(OUTVAR(SFR).VALUE>0) THEN
                        if(solver_control.i2ncal==SPR) call sfr_extrapolation_stress_strain_cal(i) !for spr method, no nodal sfr was calculated in an element way. 
                        
					    if(element(i).sfr(1,n1+j)>node(element(i).node(j)).sfr(1)) then					
						    node(element(i).node(j)).sfr(1)=element(i).sfr(1,n1+j)
						    node(element(i).node(j)).sfr(2)=element(i).mat !borrow 
					    endif
				    ENDIF
				end do
			case(spg2d,spg,cax_spg)
				do j=1,element(i).nnum
                    IF(solver_control.i2ncal/=SPR) THEN
					    n2=node(element(i).node(j)).nelist_SPG
					    node(element(i).node(j)).igrad=node(element(i).node(j)).igrad &
					    +element(i).igrad(:,n1+j)/n2
					    node(element(i).node(j)).velocity=node(element(i).node(j)).velocity &
					    +element(i).velocity(:,n1+j)/n2
					    node(element(i).node(j)).kr=node(element(i).node(j)).kr &
					    +element(i).kr(n1+j)/n2	
					    node(element(i).node(j)).mw=node(element(i).node(j)).mw &
					    +element(i).mw(n1+j)/n2	
                    ENDIF
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
        if(node(i).isactive==0) cycle
		if(allocated(node(i).stress)) then
			node(i).mises=0.50*((node(i).stress(1)-node(i).stress(2))**2 &
									+(node(i).stress(2)-node(i).stress(3))**2 &
									+(node(i).stress(1)-node(i).stress(3))**2 &
									+6*(node(i).stress(4)**2+node(i).stress(5)**2 &
									+node(i).stress(6)**2))
			node(i).mises=node(i).mises**0.5
            
            IF(OUTVAR(SFR).VALUE>0) THEN
			    !节点的材料假定为破坏比最大的单元材料，以模拟成层土中的软弱夹层
			    C1=material(node(i).sfr(2)).GET(3,ISTEP)
			    Phi1=material(node(i).sfr(2)).GET(4,ISTEP)
                !dis1(1:ndimension)=Tdisp(node(i).dof(1:ndimension))
			    NODE(I).SFR(7)=C1;NODE(I).SFR(8)=PHI1;
			    call stress_in_failure_surface(node(i).sfr,node(i).stress,2,C1,Phi1,solver_control.slidedirection,node(i).coord(ndimension))
			ENDIF
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


subroutine sfr_extrapolation_stress_strain_cal(ienum)
!it is an identical sub with extrapolation_stress_strain_cal except that is only cal the sfr.
	use solverds
	implicit none
	integer,intent(in)::ienum
	integer::i,j,k=0,n1,n2

	!do i=1,enum
	n1=element(ienum).ngp+1
	n2=element(ienum).ngp+element(ienum).nnum

	SELECT CASE(ELEMENT(IENUM).ET)
		CASE(CAX3_SPG,CPE3_SPG,CPE4R_SPG,TET4_SPG,CAX4R_SPG,&
			 CAX3,CPE3,CPE4R,TET4,CAX4R,&
			 CAX3_CPL,CPE3_CPL,CPE4R_CPL,TET4_CPL,CAX4R_CPL)
			do concurrent (i=n1:n2)
				IF(ALLOCATED(ELEMENT(IENUM).SFR)) ELEMENT(IENUM).SFR(:,i)=ELEMENT(IENUM).SFR(:,1)
			end do
		CASE(CPE6_SPG,CPE8R_SPG,CPE4_SPG,PRM15_SPG,TET10_SPG,CAX6_SPG,CAX8R_SPG,CAX4_SPG,&
			 CPE6_CPL,CPE8R_CPL,CPE4_CPL,PRM15_CPL,TET10_CPL,CAX6_CPL,CAX8R_CPL,CAX4_CPL,&
			 CPE6,CPE8R,CPE4,PRM15,TET10,CAX6,CAX8R,CAX4,ZT4_SPG)
			do concurrent (i=n1:n2)
				

				do CONCURRENT (j=1:6)
					IF(ALLOCATED(ELEMENT(IENUM).SFR)) ELEMENT(IENUM).SFR(j,i)=dot_product( &
						ELEMENT(IENUM).SFR(j,1:element(ienum).ngp), &
						ecp(element(ienum).et).expolating_Lshape(1:element(ienum).ngp,i-n1+1))	
				enddo			
				
			end do				
		CASE(CPE15_SPG,CAX15_SPG,&
			 CPE15_CPL,CAX15_CPL,&
			 CPE15,CAX15)
			
			do j=1,6
				IF(ALLOCATED(ELEMENT(IENUM).SFR)) THEN
					CALL I2N_TRI15(element(ienum).SFR(j,N1:N2),element(ienum).SFR(j,1:ELEMENT(IENUM).NGP),ELEMENT(IENUM).ET)
				ENDIF
			enddo

			
			
		CASE(prm6_spg,PRM6,PRM6_CPL)
			DO CONCURRENT (I=1:2)
				N1=3*I
				N2=3*I+2
				DO CONCURRENT (J=1:6) 
					IF(ALLOCATED(ELEMENT(IENUM).SFR)) element(ienum).SFR(j,N1:N2)=dot_product(ELEMENT(IENUM).SFR(j,1:element(ienum).ngp), &
							ecp(element(ienum).et).expolating_Lshape(1:element(ienum).ngp,I))
				enddo			
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


