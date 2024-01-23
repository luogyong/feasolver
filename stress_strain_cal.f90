

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
				IF(ALLOCATED(ELEMENT(IENUM).SFR)) ELEMENT(IENUM).SFR(1:9,i)=ELEMENT(IENUM).SFR(1:9,1)
				IF(ALLOCATED(ELEMENT(IENUM).STRESS)) ELEMENT(IENUM).STRESS(:,i)=ELEMENT(IENUM).STRESS(:,1)
				IF(ALLOCATED(ELEMENT(IENUM).STRAIN)) ELEMENT(IENUM).STRAIN(:,i)=ELEMENT(IENUM).STRAIN(:,1)
				IF(ALLOCATED(ELEMENT(IENUM).PSTRAIN)) ELEMENT(IENUM).PSTRAIN(:,i)=ELEMENT(IENUM).PSTRAIN(:,1)
			end do
		CASE(CPE6_SPG,CPE8R_SPG,CPE4_SPG,PRM15_SPG,TET10_SPG,CAX6_SPG,CAX8R_SPG,CAX4_SPG,&
			 CPE6_CPL,CPE8R_CPL,CPE4_CPL,PRM15_CPL,TET10_CPL,CAX6_CPL,CAX8R_CPL,CAX4_CPL,&
			 CPE6,CPE8R,CPE4,PRM15,TET10,CAX6,CAX8R,CAX4,ZT4_SPG2)
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
				do CONCURRENT (j=1:9)
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
			do j=1,9
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
			
			
		CASE(prm6_spg,PRM6,PRM6_CPL,ZT6_SPG2)
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
				DO CONCURRENT (J=1:9) 
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
        CASE(WELLBORE,WELLBORE_SPGFACE,PIPE2,SPHFLOW,SEMI_SPHFLOW,ZT4_SPG,ZT6_SPG,poreflow)
        
		CASE DEFAULT
			PRINT *, 'NO SUCH AN ELEMENT TYPE. SUB extrapolation_stress_strain_cal'
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
    
    IF(ISUBST==1) CALL spr_initialize(ISTEP)
    
        
    
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
subroutine E2N_stress_strain(ISTEP,isubts)
	use solverds
    use stress_failure_ratio
	implicit none
	integer::i,j,k,n1,n2,ISTEP,isubts
	real(8)::un(50)=0.0D0,vangle(15)=0.0,C1,PHI1,dis1(3),T2,ts1,SFR_MAX1=-1.D20,SIGMA1(6),MU1,SFR1(9),DCOS1(3,3),t1

    if(isubts==1) call NodalWeight(ISTEP) !同一步中单元的生死不发生改变
    
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
			node(i).sfr(12)=SOLVER_CONTROL.slidedirection
        endif
		IF(ALLOCATED(NODE(I).PSIGMA)) NODE(I).PSIGMA=0.D0
		node(i).q=0.0d0
        !node(i).sfr(1)=-1.0d6
	end do
	
	!averaged simplily at nodes
	do i=1,enum
		IF(ELEMENT(I).ISACTIVE==0) CYCLE
        


		n1=element(i).ngp
		
		select case(element(i).ec)
			case(C3D,CPE,CPS,CAX,CPL,CAX_CPL)
				do j=1,element(i).nnum
                    
					IF(SOLVER_CONTROL.I2NWEIGHT==WEIGHT_ANGLE) THEN
                        T2=element(i).ANGLE(J)/node(element(i).node(j)).ANGLE
                    ELSE
                        T2=1./(node(element(i).node(j)).nelist-node(element(i).node(j)).nelist_SPG) 
                    ENDIF
                    
                    IF(solver_control.i2ncal/=SPR) THEN
                        

                        
					    node(element(i).node(j)).stress=node(element(i).node(j)).stress &
					    +element(i).stress(:,n1+j)*T2
					    node(element(i).node(j)).strain=node(element(i).node(j)).strain &
					    +element(i).strain(:,n1+j)*T2
					    node(element(i).node(j)).pstrain=node(element(i).node(j)).pstrain &
					    +element(i).pstrain(:,n1+j)*T2
                    ENDIF
                    
                    
                    
					!get the mat for later nodal sfr cal. 
                    IF(OUTVAR(SFR).VALUE>0) THEN
                        if(solver_control.i2ncal==SPR) call sfr_extrapolation_stress_strain_cal(i) !for spr method, no nodal sfr was calculated in an element way. 
                        
					    if(element(i).sfr(1,n1+j)>node(element(i).node(j)).sfr(1)) then					
						    node(element(i).node(j)).sfr(1)=element(i).sfr(1,n1+j)
						    node(element(i).node(j)).sfr(2)=element(i).mat !borrow 
                        endif
                        !SFR_KR
                        node(element(i).node(j)).sfr(7)=node(element(i).node(j)).sfr(7)+element(i).sfr(7,n1+j)*T2
				    ENDIF
				end do
			case(spg2d,spg,cax_spg)
            

				do j=1,element(i).nnum
                    node(element(i).node(j)).q=node(element(i).node(j)).q+element(i).flux(j)
                    
                    IF(ELEMENT(I).ET==PIPE2.OR.ELEMENT(I).ET==POREFLOW) CYCLE
                    !recover the head value at the nodes along the wellbore
                    IF(ELEMENT(I).ET==WELLBORE.OR.ELEMENT(I).ET==WELLBORE_SPGFACE) THEN
                        TDISP(NODE(ELEMENT(I).NODE(4)).DOF(4))=TDISP(NODE(ELEMENT(I).NODE(1)).DOF(4))
                        TDISP(NODE(ELEMENT(I).NODE(3)).DOF(4))=TDISP(NODE(ELEMENT(I).NODE(2)).DOF(4))
                        CYCLE
                    ENDIF
                        
                    IF(ELEMENT(I).ET==SPHFLOW.OR.ELEMENT(I).ET==SEMI_SPHFLOW) THEN
                        TDISP(ELEMENT(I).G(2))=TDISP(ELEMENT(I).G(1))
                        CYCLE
                    ENDIF
                    
                    !IF(ELEMENT(I).ET==ZT4_SPG.OR.ELEMENT(I).ET==ZT6_SPG) CYCLE
                    
                    IF(solver_control.i2ncal/=SPR) THEN
					    IF(SOLVER_CONTROL.I2NWEIGHT==WEIGHT_ANGLE) THEN
                            T2=element(i).ANGLE(J)/node(element(i).node(j)).ANGLE
                        ELSE
                            T2=1./node(element(i).node(j)).nelist_SPG
                        ENDIF                    
                    
					    
                        !T2=node(element(i).node(j)).ANGLE/VANGLE(J)

					    node(element(i).node(j)).igrad=node(element(i).node(j)).igrad &
					    +element(i).igrad(:,n1+j)*T2
					    node(element(i).node(j)).velocity=node(element(i).node(j)).velocity &
					    +element(i).velocity(:,n1+j)*T2
					    node(element(i).node(j)).kr=node(element(i).node(j)).kr &
					    +element(i).kr(n1+j)*T2	
					    node(element(i).node(j)).mw=node(element(i).node(j)).mw &
					    +element(i).mw(n1+j)*T2
                        
                        !N2=element(i).node(j)
                        !IF(N2==1) THEN
                        !    N2=1
                        !ENDIF
                    ENDIF
					
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
			    !节点的材料假定为破坏比最大的单元材料，以模拟成层土中的软弱夹屄1?7
                MU1=material(node(i).sfr(2)).GET(2,ISTEP)
			    C1=material(node(i).sfr(2)).GET(3,ISTEP)
			    Phi1=material(node(i).sfr(2)).GET(4,ISTEP)
                TS1=material(node(i).sfr(2)).GET(5,ISTEP)
                IF(material(node(i).sfr(2)).TYPE==MC) THEN
                    TS1=1.D6 !NO TENSION CUT OFF
                ENDIF
                !dis1(1:ndimension)=Tdisp(node(i).dof(1:ndimension))
			    NODE(I).SFR(10)=C1;NODE(I).SFR(11)=PHI1;
                !call stress_in_failure_surface(node(i).sfr,node(i).stress,ndimension,C1,Phi1,solver_control.slidedirection,node(i).coord(1:ndimension),TS1)
          
			    call stress_in_failure_surface(node(i).sfr(1:9),node(i).stress,ndimension,C1,Phi1,solver_control.slidedirection,node(i).coord(1:ndimension),TS1,DCOS1)
                !假定ko应力，ko=v/(1-v),sxx=k0*syy,szz=sxx,txy=0
                SIGMA1(ndimension)=node(i).stress(ndimension)
                if(ndimension==2) then                    
            	    SIGMA1(1)=mu1/(1-mu1)*SIGMA1(ndimension);SIGMA1(3)=SIGMA1(1);SIGMA1(4:6)=0 
			    else
				    SIGMA1(1)=mu1/(1-mu1)*SIGMA1(ndimension);SIGMA1(2)=SIGMA1(1);SIGMA1(4:6)=0 
			    endif                    
	            call stress_in_failure_surface(sfr1,SIGMA1,ndimension,C1,Phi1,solver_control.slidedirection,node(i).coord(1:ndimension),TS1,DCOS1)
                NODE(I).SFR(8)=SFR1(1)
                if(solver_control.slope_mko>0.and.NODE(I).SFR(1)>=0.d0) then
                    
                    NODE(I).SFR(1)=NODE(I).SFR(1)-SFR1(1)
                    
                endif 
                IF(NODE(I).SFR(1)>SFR_MAX1) SFR_MAX1=NODE(I).SFR(1)
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
    
	DO i=1,nnum
		
        if(node(i).isactive==0) cycle
        !IF(OUTVAR(SFR).VALUE>0) THEN
        !    IF(NODE(I).SFR(1)<0.D0) NODE(I).SFR(1)=-SFR_MAX1
        !ENDIF
		
		IF(OUTVAR(PSIGMA).VALUE>0) THEN
            CALL principal_stress_cal(node(i).stress,node(i).psigma(1:3),DCOS1)
            DO J=1,3
				if(dcos1(j,j)<0.d0) then
					t1=-1.d0
				else
					t1=1.d0
				endif
                node(i).psigma(3*J+1:3*J+3)=node(i).psigma(J)*DCOS1(J,:)*t1
            ENDDO
			IF(NDIMENSION==2) THEN
			!call principal_stress_cal(node(i).psigma,node(i).stress,ndimension)
                node(i).psigma(13)=SIGN(ACOS(DCOS1(1,1)),DCOS1(1,2))
            ELSE
                node(i).psigma(13)=SIGN(ACOS(DCOS1(1,1)),DCOS1(1,3))
			    !IF(node(i).stress(1)>node(i).stress(2)) THEN
				   ! node(i).psigma(13)=node(i).psigma(13)+1.570796326794897/2.0D0
			    !ENDIF
			ENDIF
		ENDIF
	ENDDO
	
	
	
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
				IF(ALLOCATED(ELEMENT(IENUM).SFR)) ELEMENT(IENUM).SFR(1:9,i)=ELEMENT(IENUM).SFR(1:9,1)
			end do
		CASE(CPE6_SPG,CPE8R_SPG,CPE4_SPG,PRM15_SPG,TET10_SPG,CAX6_SPG,CAX8R_SPG,CAX4_SPG,&
			 CPE6_CPL,CPE8R_CPL,CPE4_CPL,PRM15_CPL,TET10_CPL,CAX6_CPL,CAX8R_CPL,CAX4_CPL,&
			 CPE6,CPE8R,CPE4,PRM15,TET10,CAX6,CAX8R,CAX4,ZT4_SPG2)
			do concurrent (i=n1:n2)
				

				do CONCURRENT (j=1:9)
					IF(ALLOCATED(ELEMENT(IENUM).SFR)) ELEMENT(IENUM).SFR(j,i)=dot_product( &
						ELEMENT(IENUM).SFR(j,1:element(ienum).ngp), &
						ecp(element(ienum).et).expolating_Lshape(1:element(ienum).ngp,i-n1+1))	
				enddo			
				
			end do				
		CASE(CPE15_SPG,CAX15_SPG,&
			 CPE15_CPL,CAX15_CPL,&
			 CPE15,CAX15)
			
			do j=1,9
				IF(ALLOCATED(ELEMENT(IENUM).SFR)) THEN
					CALL I2N_TRI15(element(ienum).SFR(j,N1:N2),element(ienum).SFR(j,1:ELEMENT(IENUM).NGP),ELEMENT(IENUM).ET)
				ENDIF
			enddo

			
			
		CASE(prm6_spg,PRM6,PRM6_CPL,ZT6_SPG2)
			DO CONCURRENT (I=1:2)
				N1=3*I
				N2=3*I+2
				DO CONCURRENT (J=1:9) 
					IF(ALLOCATED(ELEMENT(IENUM).SFR)) element(ienum).SFR(j,N1:N2)=dot_product(ELEMENT(IENUM).SFR(j,1:element(ienum).ngp), &
							ecp(element(ienum).et).expolating_Lshape(1:element(ienum).ngp,I))
				enddo			
			END DO
        CASE(WELLBORE,WELLBORE_SPGFACE,PIPE2,SPHFLOW,SEMI_SPHFLOW,ZT4_SPG,ZT6_SPG,POREFLOW)
        
		CASE DEFAULT
			PRINT *, 'NO SUCH A ELEMENT TYPE. SUB extrapolation_stress_strain_cal'
			STOP					
	END SELECT
			
		!CASE DEFAULT
			
		!	PRINT *, 'POSTPROCEDURE IN SOLID ELEMENT CLASS IS NOT COMPLETED.'
			!STOP
	!END SELECT
	
	
end subroutine


subroutine NodalWeight(ISTEP)
	use solverds
	implicit none
    integer,intent(in)::istep
	integer::i,j,k,n1,CESET1=0
	LOGICAL::ISSPG1
    
	
	
	node.nelist=0;node.nelist_SPG=0;node.angle=0
	do i=1,enum
        if(solver_control.i2nweight==weight_angle.OR.solver_control.I2NCAL==SPR) then
            if(.not.allocated(element(i).angle))    then
                  call calangle(i)
            endif
        endif
        
        if(element(i).isactive==0) cycle
        !line element,such as wellbore and pipe2 are excluded.
        if(element(i).et==wellbore.or.element(i).et==pipe2.OR.element(i).et==WELLBORE_SPGFACE.OR.element(i).et==POREFLOW) cycle
        if(element(i).et==SPHFLOW.or.element(i).et==SEMI_SPHFLOW) cycle
        
		ISSPG1=ELEMENT(I).EC==SPG.OR.ELEMENT(I).EC==SPG2D.OR.ELEMENT(I).EC==CAX_SPG
		do j=1,element(i).nnum
			n1=element(i).node(j)
			node(n1).nelist=node(n1).nelist+1
			IF(ISSPG1) THEN
				node(n1).nelist_SPG=node(n1).nelist_SPG+1
			ENDIF
            if(allocated(element(i).angle)) THEN
                CESET1=ESET(ELEMENT(I).SET).COUPLESET
                IF(.NOT.(CESET1>0.AND.CESET1<ELEMENT(I).SET)) THEN
                    node(n1).angle=node(n1).angle+element(i).angle(j)
                ENDIF
            ENDIF
		end do
	end do
    
	!do i=1,nnum
	!	allocate(node(i).elist(node(i).nelist),node(i).elist_SPG(node(i).nelist_SPG))
	!end do
	!node.nelist=0;node.nelist_SPG=0
	!
	!do i=1,enum
	!	ISSPG1=ELEMENT(I).EC==SPG.OR.ELEMENT(I).EC==SPG2D.OR.ELEMENT(I).EC==CAX_SPG
	!	do j=1,element(i).nnum
	!		n1=element(i).node(j)
	!		node(n1).nelist=node(n1).nelist+1
	!		node(n1).elist(node(n1).nelist)=i
	!		IF(ISSPG1) THEN
	!			node(n1).nelist_SPG=node(n1).nelist_SPG+1
	!			node(n1).elist_SPG(node(n1).nelist_SPG)=i
	!		ENDIF
	!	end do
	!end do	

end subroutine


! calculate the internal angles for each node in a element, for the post-procedure.
subroutine calangle(ienum)
    use SolverMath
	use solverds
	implicit none
	integer::i,j,IESET1
	integer::ienum,idim=2,jdim=2,v1,v2,IA1(4),IA2D1(2,10),IA2D2(4,6)
	real(8)::vec1(2,2)=0.0,angle=0.0,vangle(15)=0.0,ar1(3,10),VN1(3,15)=0
	
	!assume that all element edges are straight lines. such that the angle
	!of the nodes inside a line equel to pi.
	! the angle for inner nodes equel to 2pi	
	IF(ALLOCATED(ELEMENT(IENUM).ANGLE)) RETURN
    
	allocate(element(ienum).angle(element(ienum).nnum))
    IESET1=ELEMENT(IENUM).SET
    IF(ESET(IESET1).COUPLESET>0.AND.ESET(IESET1).COUPLESET<IESET1) THEN !GHOST
        I=ESET(ESET(IESET1).COUPLESET).ENUMS+(IENUM-ESET(IESET1).ENUMS)
        ELEMENT(IENUM).ANGLE=ELEMENT(I).ANGLE
        RETURN
    ENDIF
    
	select case(ecp(element(ienum).et).shtype)
		case(tri3,tri6,tri15)
            element(ienum).angle=pi()
			do i=1,2
				v1=i-1
				if(v1<1) v1=3
				v2=i+1
				vec1(:,1)=node(element(ienum).node(v1)).coord(1:2)-node(element(ienum).node(i)).coord(1:2)
				vec1(:,2)=node(element(ienum).node(v2)).coord(1:2)-node(element(ienum).node(i)).coord(1:2)
				call vecangle(vec1,idim,jdim,angle)
				element(ienum).angle(i)=angle
				element(ienum).angle(3)=element(ienum).angle(3)-angle
			end do
            IF(ELEMENT(IENUM).NNUM==15) then
			    element(ienum).angle(13:15)=2.*pi()
            endif
            
		case(qua4,qua8)
		
			element(ienum).angle=pi()
			element(ienum).angle(4)=2*PI()
			do i=1,3
				v1=i-1
				if(v1<1) v1=4
				v2=i+1
				vec1(:,1)=node(element(ienum).node(v1)).coord(1:2)-node(element(ienum).node(i)).coord(1:2)
				vec1(:,2)=node(element(ienum).node(v2)).coord(1:2)-node(element(ienum).node(i)).coord(1:2)
				IF(ELEMENT(IENUM).ET==ZT4_SPG2) THEN
					vec1(:,1)=GNODE(1:2,element(ienum).node2(v1))-Gnode(1:2,element(ienum).node2(i))
					vec1(:,2)=GNODE(1:2,element(ienum).node2(v2))-Gnode(1:2,element(ienum).node2(i))			
				ENDIF
				call vecangle(vec1,idim,jdim,angle)
				element(ienum).angle(i)=angle
				element(ienum).angle(4)=element(ienum).angle(4)-angle				
			end do
			

        case(tet4,tet10)
            
            !do i=1,4
            !    IA1(1)=I;IA1(2)=MOD(IA1(1),4)+1;IA1(3)=MOD(IA1(2),4)+1;IA1(4)=MOD(IA1(3),4)+1
            !    do j=1,4
            !        ar1(:,j)=node(element(ienum).node(IA1(j))).coord
            !    enddo
            !    element(ienum).angle(i)=solidangle(AR1(:,1:4))            
            !enddo
            do j=1,4
				ar1(:,j)=node(element(ienum).node(j)).coord
            enddo
            call tetrahedron_solid_angles_3d(ar1(:,1:4),element(ienum).angle(1:4))
            
            IF(ELEMENT(IENUM).NNUM>4) THEN
                do i=1,4
                    ar1(:,i)=node(element(ienum).node(i)).coord 
                enddo
                !NOMRAL VECTOR OF THE FACE
                VN1(:,1)=NORMAL_TRIFACE(AR1(:,[2,1,3]))
                VN1(:,2)=NORMAL_TRIFACE(AR1(:,[1,2,4]))
                VN1(:,3)=NORMAL_TRIFACE(AR1(:,[2:4]))
                VN1(:,4)=NORMAL_TRIFACE(AR1(:,[3,1,4]))
                IA2D1(:,1:6)=RESHAPE([1,2,1,3,1,4,2,4,2,3,3,4],([2,6]))
                DO I=5,ELEMENT(IENUM).NNUM
                    ELEMENT(IENUM).ANGLE(I)=DihedralAngle(VN1(:,IA2D1(1,I-4)),VN1(:,IA2D1(2,I-4)))        
                ENDDO
            
            ENDIF
            
            
        CASE(PRM6,PRM15)
        
            
            IA2D2=RESHAPE([1,2,3,4,&
                           2,3,1,5,&
                           3,1,2,6,&
                           4,6,5,1,&
                           5,4,6,2,&
                           6,5,4,3],([4,6]))
                           
            do i=1,6
                do j=1,4
                    IF(ELEMENT(IENUM).ET==ZT6_SPG2) THEN
                        ar1(:,j)=Gnode(:,element(ienum).node(IA2D2(j,i)))
                    ELSE
                        ar1(:,j)=node(element(ienum).node(IA2D2(j,i))).coord
                    ENDIF
                enddo
                call tetrahedron_solid_angles_3d(ar1(:,1:4),vangle(1:4))
                element(ienum).angle(i)=vangle(1)            
            enddo
            
            IF(ELEMENT(IENUM).NNUM>6) THEN
                do i=1,6
                    IF(ELEMENT(IENUM).ET==ZT6_SPG2) THEN
                        ar1(:,i)=Gnode(:,element(ienum).node(i))
                    ELSE
                        ar1(:,i)=node(element(ienum).node(i)).coord 
                    ENDIF
                        
                enddo
                !NOMRAL VECTOR OF THE FACE
                VN1(:,1)=NORMAL_TRIFACE(AR1(:,[2,1,3]))
                VN1(:,2)=NORMAL_TRIFACE(AR1(:,[1,2,4]))
                VN1(:,3)=NORMAL_TRIFACE(AR1(:,[2,3,5]))
                VN1(:,4)=NORMAL_TRIFACE(AR1(:,[3,1,4]))
                VN1(:,5)=NORMAL_TRIFACE(AR1(:,[4,5,6]))
                
                IA2D1(:,1:9)=RESHAPE([1,2,1,3,1,4,2,5,3,5,4,5,2,4,2,3,3,4],([2,9]))
                DO I=7,ELEMENT(IENUM).NNUM
                    ELEMENT(IENUM).ANGLE(I)=DihedralAngle(VN1(:,IA2D1(1,I-6)),VN1(:,IA2D1(2,I-6)))        
                ENDDO
            
            ENDIF
            
    case default 
        
        stop 'no sub element shapetype. To be improved in cal element angle.'
        
        
	end select
	

	
end subroutine

! given vec(2,2),calculate the angle
! vec(:,1) for the first vector, and vec(:,2) for the other
subroutine vecangle(vec,irow,jcol,angle)
	implicit none
	integer::irow,jcol
	real(8)::vec(irow,jcol),angle
	real(8)::r1=0.0,r2=0.0,r12=0.0
	
	r12=vec(1,1)*vec(1,2)+vec(2,1)*vec(2,2)
	r1=(vec(1,1)**2+vec(2,1)**2)**0.5
	r2=(vec(1,2)**2+vec(2,2)**2)**0.5
	angle=r12/r1/r2
	angle=dacos(angle)
end subroutine
