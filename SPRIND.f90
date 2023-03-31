module stress_failure_ratio
    use solverds,ONLY:NDIMENSION,SOLVER_CONTROL
    use SolverMath
    USE DS_SlopeStability
    implicit none
    PUBLIC::stress_in_failure_surface,principal_stress_cal,stress_on_plane,stress_in_inclined_plane

    PRIVATE
    
    INTERFACE stress_in_failure_surface
        MODULE PROCEDURE stress_in_failure_surface_3D,stress_in_failure_surface_pe
    END INTERFACE
    INTERFACE principal_stress_cal
        MODULE PROCEDURE principal_stress_cal_3d,principal_stress_cal_pe
    END INTERFACE    
    
    
    CONTAINS


    subroutine stress_in_failure_surface_3D(sfr,stress,ndim,cohesion,PhiD,slidedirection,Y,TensileS,trans,direction,iflag)
    
        implicit none
        integer,intent(in)::ndim,slidedirection
        real(8),intent(in)::stress(6),cohesion,PhiD,Y(NDIM),TensileS
        real(8),intent(inout)::sfr(:),trans(3,3)
        real(8),optional,intent(in)::direction(3)
        integer,optional,intent(in)::iflag !if iflag<>0, only calculate sfr(1) is neended.
        
        integer::i,n1,n2,iflag1=0
        real(8)::ss1(6),pss1(3),t1,t2,C1,phi1,PI1,sin1,cos1,t3,t4,t5,ta1(4),A,B,an(3,3),TRANS1(3,3),D1(3)
        real(8)::sigmaC1,R1,sita1,YF1,hstress1,BETA1
        logical::isazone=.true.
        
        PI1=dATAN(1.0)*4.0
        C1=cohesion
        phi1=PhiD/180*PI1
        sfr=0.d0
        trans=0.d0

        if(present(iflag)) then
            iflag1=iflag
        else
            iflag1=0
        endif

        CALL principal_stress_cal(STRESS, PSS1, AN)

            
        t2=(pss1(1)-pss1(3))/2.0
        
        sigmaC1=(pss1(1)+pss1(3))/2.0
        R1=t2
        if(R1>0) then
            T1=1E-7
            IF(PHI1>1E-7) T1=TAN(PHI1)
        
            if(pss1(1)<MIN(TensileS,C1/T1).and.SIGMAC1<0.0d0) then
                !CHECK YIELD
                !YF1=SIGMAC1*SIN(PHI1)+R1-C1*COS(PHI1)
                IF(PHI1>0.D0) THEN
                    YF1=R1/((C1/TAN(PHI1)-SIGMAC1)*SIN(PHI1))
                ELSE
                    YF1=R1/MAX(C1,1E-6)
                ENDIF
                !IF(YF1>0.D0.and.pss1(1)<MIN(TensileS,C1/T1)) THEN
                IF(YF1<1.D0.AND.YF1>0.D0) THEN
            
                    B=-tan(phi1)
                    A=(C1+sigmaC1*B)/R1
                    !W破坏面与主应力面的夹角
                    !if(r1<=C1*cos(phi1)-sigmaC1*sin(phi1)) then
                    if(abs(A)>=abs(B)) then
                        sita1=asin((1.d0-(B/A)**2)**0.5)
                    else
                        sita1=pi1/2.0-phi1
                    endif
                    sfr(1)=min(abs(sin(sita1)/(A+B*cos(sita1))),20.d0) !sfr for stress failure ratio,max is less than 20
                
                ELSE

                    sita1=pi1/2.0-phi1
                    
                    SFR(1)=min(YF1,20.0)

                ENDIF
                !if(sfr(1)>1.d0) sfr(1)=1.0d0
            else
                !tension
                sfr(1)=-1.            
                sita1=pi1/2.0
                if(Solver_control.slope_isTensionCrack<1) sita1=pi1/2.0-phi1
            endif
        else
            sfr(1)=0.001
            sita1=pi1/2.0-phi1
        endif

        if(iflag1/=0) return
        !
        t1=slidedirection

        !一般在坡顶和坡脚，会出现sxx<syy(水平方向的压力大于竖各压力的情况，按下面的方法计算，其破坏面的夹角与x轴的夹角大于0
        !这时，当slidedirection=right时，坡顶会出现不相容的破坏面方向(sfr(2)>0),这时令其取另一个方向的破坏面
       
        if(ndim==2) then
            hstress1=(stress(1)+stress(3))/2.0d0            
        else
            hstress1=(stress(1)+stress(2))/2.0d0 
        endif
        if(stress(ndim)>hstress1) then
            isazone=.false. 
        else
            isazone=.true.
        endif

        IF(slidedirection==1.and.SLOPEPARAMETER.ISYDWZ.AND.Y(NDIM)>SLOPEPARAMETER.YDOWNWARDZONE.AND.(.not.isazone)) THEN
            T1=-T1
        ENDIF
        IF(PRESENT(DIRECTION)) THEN
            N1=2
            T2=NORM2(DIRECTION)
            IF(T2>1.E-10) THEN
                D1=DIRECTION/T2
            ELSE
                D1=0.0d0
                D1(NDIM)=-1.D0 !downward
            ENDIF  
        ELSE
            N1=1
        ENDIF
        T2=1.D20
        DO I=1,N1
            !BETA1为破坏面与大主应力方向的夹角
            IF(sfr(1)>=0.d0.or.Solver_control.slope_isTensionCrack<1) then
                IF(PRESENT(DIRECTION)) THEN
                    BETA1=0.5*(PI1+((-1)**I)*SITA1)
                ELSE                    
                    BETA1=0.5*(PI1+slidedirection*SITA1*sign(1.0,stress(2*ndim))*(-1)**NDIM)
                ENDIF
                if(isazone) then          
                    !调整坡脚主动区的顶部最大滑动面方向
                    !IF(ABS(SLOPEPARAMETER.TOEZONE(2))/=0.D0) THEN
                    !    T3=-(SLOPEPARAMETER.TOEZONE(1)*Y(1)+SLOPEPARAMETER.TOEZONE(3))/SLOPEPARAMETER.TOEZONE(2)
                    !    T3=Y(ndim)-T3
                    !    IF(T3<=0.D0) THEN
                    !        sfr(2)=(PI1-sita1)/2.0*T1
                    !    ENDIF
                    !ELSEIF(ABS(SLOPEPARAMETER.TOEZONE(1))/=0.D0) THEN
                    !    T3=-SLOPEPARAMETER.TOEZONE(3)/SLOPEPARAMETER.TOEZONE(1)
                    !    T3=Y(1)-T3
                    !    IF(T3*slidedirection>=0) sfr(2)=(PI1-sita1)/2.0*T1
                    !ENDIF
                endif

            ELSE
                !TENSION FAILURE,THE FAILURE SURFACE IS ACTING PLANE OF THE MAJOR PRINCIPLE STRESS.
                BETA1=0.5*PI1
                N1=1
            ENDIF
        
            !failure plane in principal stress system
            trans(1,1)=cos(beta1);trans(1,3)=sin(beta1);trans(1,2)=0.d0
            !in globle system
            trans(1,:)=matmul(TRANSPOSE(an),trans(1,:)) !kmax axis
            !this trans(1,:) and the middle principal forms the failure plane. the Kmax is along the trans(1,:)
            trans(2,:)=an(2,:)
            trans(3,:)=cs_vector (trans(1,:),trans(2,:)) !direction cosine of the failure plane,kmin axis
            IF(N1==2) THEN
                T3=ABS(DOT_PRODUCT(TRANS(3,:),D1))
                IF(T3<T2) THEN
                    T2=T3;
                    TRANS1=TRANS 
                ENDIF
                IF(I==2) TRANS=TRANS1
            ENDIF   
        ENDDO

        sfr(2)=sign(acos(trans(1,1)),trans(1,NDIM)) !破坏面与x轴的夹角

        IF(sfr(1)>=0.d0.or.Solver_control.slope_isTensionCrack<1) then
            !破坏面上的正应力及剪应力
            call stress_on_plane(stress,trans(3,:),sfr(3:4)) !与原来二维计算相比，仅平面上的应力已经没有正负之分.
            
            sfr(5)=sign(sfr(1),-sfr(4))*TRANS(1,1) !输出归一化的剪应力
            sfr(6)=sign(sfr(1),-sfr(4))*TRANS(1,2)
            sfr(9)=sign(sfr(1),-sfr(4))*TRANS(1,3) 
        else
            SFR(5)=0.D0
            SFR(6)=0.D0
            sfr(9)=-1
            IF(NDIM==2) THEN
                SFR(6)=-1;SFR(9)=0.D0
            ENDIF            
        endif

        sfr(2)=sfr(2)/PI1*180.0 !to degree
        
        
    
    endsubroutine   
    
   
    subroutine stress_on_plane(s,pcos,sop,sdir,tau_sdir)
    !given:stress(s(6)) and plane direction cosine(pcos(3))
    !output:normal and tangent stress on the plane. (sop(2))  
        implicit none
        real(8),intent(in)::s(6),pcos(3)
        real(8),intent(out)::sop(2)
        real(8),intent(in),optional::sdir(3) !方向矢量sdir
        real(8),intent(out),optional::tau_sdir !该平面沿方向sdir上的应力
        real(8)::t(3),T1

        sop(1)=s(1)*pcos(1)**2+s(2)*pcos(2)**2+s(3)*pcos(3)**2 + &
            2*(s(4)*pcos(1)*pcos(2)+s(5)*pcos(2)*pcos(3)+s(6)*pcos(1)*pcos(3))
        
        t(1)=s(1)*pcos(1)+s(4)*pcos(2)+s(6)*pcos(3)
        t(2)=s(4)*pcos(1)+s(2)*pcos(2)+s(5)*pcos(3)
        t(3)=s(6)*pcos(1)+s(5)*pcos(2)+s(3)*pcos(3)
        T1=t(1)**2+t(2)**2+t(3)**2-sop(1)**2
        
        IF(T1>0.D0) THEN
            sop(2)=(T1)**0.5
        ELSE
            sop(2)=0.d0
        ENDIF

        if(present(sdir).and.present(tau_sdir)) then
            tau_sdir=dot_product(t,sdir)
        endif
    endsubroutine
    !THE 2D VERSION OF SUB stress_on_plane 
    SUBROUTINE stress_in_inclined_plane(ss,RAD,SNT)
    !rad，stress plane (not its direction ) angle with x axis(),in rad.
    !ss:sx,sy,sxy
    !ss1:sn,st
        implicit none
        real(8),intent(in)::ss(3),RAD
        real(8),INTENT(OUT)::Snt(2)
    
        !Snt(1)=0.5*(ss(1)+ss(2))+0.5*(ss(1)-ss(2))*cos(2*RAD)+ss(3)*sin(2*RAD)
        !Snt(2)=-0.5*(ss(1)-ss(2))*sin(2*RAD)+ss(3)*cos(2*RAD) 
        Snt(1)=0.5*(ss(1)+ss(2))-0.5*(ss(1)-ss(2))*cos(2*RAD)-ss(3)*sin(2*RAD)
        Snt(2)=0.5*(ss(1)-ss(2))*sin(2*RAD)-ss(3)*cos(2*RAD)
    
    endSUBROUTINE    
    
    subroutine stress_in_failure_surface_pe(sfr,stress,ndim,cohesion,PhiD,slidedirection,Y,TensileS)

	    implicit none
	    integer,intent(in)::ndim,slidedirection
	    real(8),intent(in)::stress(6),cohesion,PhiD,Y(NDIM),TensileS
	    real(8),intent(inout)::sfr(:)
	
	    integer::i
	    real(8)::ss1(6),pss1(4),t1,t2,C1,phi1,PI1,sin1,cos1,t3,t4,t5,ta1(4),A,B
	    real(8)::sigmaC1,R1,sita1,YF1
    
	    PI1=dATAN(1.0)*4.0
	    C1=cohesion
	    phi1=PhiD/180*PI1
	    sfr=0.d0
   
	    call principal_stress_cal(stress,pss1)
	
	    t2=(pss1(1)-pss1(3))/2.0
    
        sigmaC1=(stress(1)+stress(2))/2.0
        R1=t2
        if(R1>0) then
            T1=1E-7
            IF(PHI1>1E-7) T1=TAN(PHI1)
        
            !if(pss1(1)<MIN(TensileS,C1/T1).or.Solver_control.slope_isTensionCrack/=1) then  !c=0, sand, no tension crack
            if(pss1(1)<MIN(TensileS,C1/T1).and.SIGMAC1<0.0d0) then
                !CHECK YIELD
                !YF1=SIGMAC1*SIN(PHI1)+R1-C1*COS(PHI1)
                IF(PHI1>0.D0) THEN
                    YF1=R1/((C1/TAN(PHI1)-SIGMAC1)*SIN(PHI1))
                ELSE
                    YF1=R1/MAX(C1,1E-6)
                ENDIF
                !IF(YF1>0.D0.and.pss1(1)<MIN(TensileS,C1/T1)) THEN
                IF(YF1<1.D0.AND.YF1>0.D0) THEN
            
                    B=-tan(phi1)
                    A=(C1+sigmaC1*B)/R1
                    !朄1?7大破坏面与主应力面的夹角的1?7?1?7
                    !if(r1<=C1*cos(phi1)-sigmaC1*sin(phi1)) then
                    if(abs(A)>=abs(B)) then
                        sita1=asin((1.d0-(B/A)**2)**0.5)
                    else
                        sita1=pi1/2.0-phi1
                    endif
                    sfr(1)=min(abs(sin(sita1)/(A+B*cos(sita1))),20.d0) !sfr for stress failure ratio,max is less than 20
                
                ELSE
                    !if slope_isTensionCrack==0, 不进行这个检查sfr(1)会产生大敄1?7原因时莫尔圆圆心与破坏线与x轴的交点非常迄1?7
                    !IF(YF1<0.D0.OR.pss1(1)>=MIN(TensileS,C1/T1)) THEN 
                    !    !tension
                    !    sfr(1)=-1.            
                    !    sita1=pi1/2.0
                    !    if(Solver_control.slope_isTensionCrack<1) sita1=pi1/2.0-phi1
                    !ELSE
                    sita1=pi1/2.0-phi1
                    
                    SFR(1)=min(YF1,20.0)
                    !ENDIF
                ENDIF
                !if(sfr(1)>1.d0) sfr(1)=1.0d0
            else
                !tension
                sfr(1)=-1.            
                sita1=pi1/2.0
                if(Solver_control.slope_isTensionCrack<1) sita1=pi1/2.0-phi1
            endif
        !if(isnan(sfr(1))) then
        !    pause
        !endif
        else
            sfr(1)=0.001
            sita1=pi1/2.0-phi1
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
        !    sfr(1)=-1. !表拉坄1?7
        !endif
        !
        t1=slidedirection
        !丄1?7般在坡顶和坡脚，会出现sxx<syy(水平方向的压力大于竖各1?7的情况，按下面的方法计算，其破坏面的夹角与x轴的夹角大于0
        !这时，当slidedirection=right时，坡顶会出现不相容的破坏面方向(sfr(2)>0),这时令其取另丄1?7个方向?1?7?1?7
        IF(slidedirection==1.and.SLOPEPARAMETER.ISYDWZ.AND.Y(NDIM)>SLOPEPARAMETER.YDOWNWARDZONE.AND.stress(2)>stress(1)) THEN
            T1=-T1
        ENDIF

        IF(sfr(1)>=0.d0.or.Solver_control.slope_isTensionCrack<1) then
            if(stress(2)>stress(1)) then !sigmax<sigmay,passive zone         
                sfr(2)=pss1(4)+(sita1)/2.0*T1
            else
                sfr(2)=pss1(4)+(sita1)/2.0*T1
                !sfr(2)=pss1(4)-(PI1-sita1)/2.0*T1 !active zone
                !调整坡脚主动区的屄1?7部最大滑动面方向
                IF(ABS(SLOPEPARAMETER.TOEZONE(2))/=0.D0) THEN
                    T3=-(SLOPEPARAMETER.TOEZONE(1)*Y(1)+SLOPEPARAMETER.TOEZONE(3))/SLOPEPARAMETER.TOEZONE(2)
                    T3=Y(2)-T3
                    IF(T3<=0.D0) THEN
                        sfr(2)=pss1(4)+(PI1-sita1)/2.0*T1
                    ENDIF
                ELSEIF(ABS(SLOPEPARAMETER.TOEZONE(1))/=0.D0) THEN
                    T3=-SLOPEPARAMETER.TOEZONE(3)/SLOPEPARAMETER.TOEZONE(1)
                    T3=Y(1)-T3
                    IF(T3*slidedirection>=0) sfr(2)=pss1(4)+(PI1-sita1)/2.0*T1
                ENDIF
            endif
            !sfr(3)=(pss1(1)+pss1(2))/2+t2*cos(PI1/2-phi1)
            !sfr(4)=-t2*sin(PI1/2-phi1)*t1 
            sfr(3)=(pss1(1)+pss1(3))/2+t2*cos(sita1)
            sfr(4)=-t2*sin(sita1)*t1 
            sfr(5)=sign(sfr(1),-sfr(4))*dcos(sfr(2)) !输出归一化的剪应劄1?7
            sfr(6)=sign(sfr(1),-sfr(4))*dsin(sfr(2))
            !sfr(5)=abs(sfr(1))*dcos(sfr(2)) !输出归一化的剪应劄1?7
            !sfr(6)=abs(sfr(1))*dsin(sfr(2))
        ELSE
            !TENSION FAILURE,THE FAILURE SURFACE IS PLANE THE MAJOR PRINCIPLE STRESS.
            IF(stress(2)>stress(1)) THEN
                SFR(2)=pss1(4)
            ELSE
                SFR(2)=pss1(4)
            ENDIF
            SFR(5)=0.D0
            SFR(6)=-1
        ENDIF
    
        sfr(2)=sfr(2)/PI1*180.0 !to degree
    
 
    endsubroutine

    subroutine principal_stress_cal_pe(stress,pstress)
	    implicit none
	    real(8),intent(in)::stress(6)
	    real(8),intent(inout)::pstress(4) 
	    real(8)::t1,t2,sita,Pi
    
        Pi=datan(1.d0)*4.0
	
	    pstress=0.0d0

		t1=(stress(1)+stress(2))/2.0d0
		t2=(stress(1)-stress(2))/2.0d0
		pstress(1)=t1+(t2**2+stress(4)**2)**0.5		
		pstress(3)=t1-(t2**2+stress(4)**2)**0.5
		Pstress(2)=stress(3)
		!默认z方向为应力为中应力?1?7?1?7
        IF(ABS(stress(1)-stress(2))>1E-14) THEN 
		    sita=0.5*atan(2*stress(4)/(stress(1)-stress(2)))
            !sita为莫尔圆主应力平面与x轴的夹角
		    Pstress(4)=sita
            if(stress(1)>stress(2)) Pstress(4)=sita+Pi/2.0
        ELSE
            Pstress(4)=0 
        ENDIF
			
		
		
		!ref::https://en.wikipedia.org/wiki/Mohr%27s_circle 
		!Mohr Circel sign cenverstion: tension is +, shear stress, clockwise is +
		
		

	
    end subroutine
    
    
    ! ****************************************************************************************
    ! **  UTILITY SUBROUTINE FOR INTERFACE: SPRINC                                          **
    ! **  CALULATE PRINCIPAL VALUES AND DIRECTION COSINES                                   **
    ! ****************************************************************************************
    ! ****************************************************************************************

    SUBROUTINE principal_stress_cal_3d(S, PS, AN, LSTR )
     ! C
    ! C======================================================================+
    ! C-----------
    ! C  INPUT :
    ! C-----------
    ! C  S     	: STRESS OR STRAIN TENSOR
    ! C  LSTR 	: FLAG DETERMINING STRESS OR STRAIN CALCULATION

    ! C-----------
    ! C  OUTPUT :
    ! C-----------
    ! C  PS(I), I=1,2,3       : THE THREE PRINCIPAL VALUES
    ! C  AN(K1,I), I=1,2,3	: THE DIRECTION COSINES OF THE PRINCIPAL DIRECTIONS CORRESPONDING TO PS(K1)
    ! C----------------------------------------------------------------------+
    ! C=======================================================================
    ! C
    ! C     Calculate stress or strain invariants based on LSTR value
 
        IMPLICIT NONE
        INTEGER,OPTIONAL:: LSTR 
        REAL(8):: S(6), PS(3), A(3), B(3), C(3), V(3), AN(3,3)
        REAL(8) I1, I2, I3, J1, J2, J3, R, T, Q, ALP, SCALC, &
            PRINC1, PRINC2, PRINC3, ARG,T1,T2,T3,RA1(3)
        REAL(8),PARAMETER::ONE=1.D0,TWO=2.D0,THREE=3.D0, TWOSEVEN = 27.D0, &
            THIRD = ONE/THREE, TWENTYSEVENTH= ONE/TWOSEVEN, &
            TWOTWENTYSEVENTH = TWO/TWOSEVEN,TTPI=2.0943951023932D0, &
            FTPI=4.1887902047864D0
        INTEGER::K,IA1(3),N1,N2,N3,L,M,N,LSTR1,N4,K1

        IF(PRESENT(LSTR)) THEN
            LSTR1=LSTR
        ELSE
            LSTR1=1
        ENDIF
        !IF(PRESENT(NDIM)) THEN
        !    NDIM1=NDIM
        !ELSE
        !    NDIM1=3
        !ENDIF 
        
        IF (LSTR1.EQ.1) THEN
            I1 = S(1)+S(2)+S(3)
            I2 = (S(1)*S(2))+(S(2)*S(3))+(S(1)*S(3)) &
                -(S(4)**2)-(S(5)**2)-(S(6)**2)
            I3 = (S(1)*S(2)*S(3))+(2*S(4)*S(5)*S(6))-(S(1)*S(5)**2) &
                -(S(2)*S(6)**2)-(S(3)*S(4)**2)

            R = (THIRD*I1**2)-I2
            T = SQRT(TWENTYSEVENTH*R**3)
            Q = (THIRD*I1*I2)-I3-(TWOTWENTYSEVENTH*I1**3)
            IF(ABS(T)<1.D-10) T=1.D-10
            ARG = -Q/(TWO*T)
            IF (ABS(ARG).GT.1.D0) THEN
            ARG = SIGN(1.D0,ARG)
            END IF
            ALP = ACOS(ARG)
            SCALC = SQRT(THIRD*R)

            PRINC1 = (2*SCALC*COS(ALP/THREE))+(THIRD*I1)
            PRINC2 = (2*SCALC*COS((ALP/THREE)+FTPI))+(THIRD*I1)
            PRINC3 = (2*SCALC*COS((ALP/THREE)+TTPI))+(THIRD*I1)

        ELSE
            J1 = S(1)+S(2)+S(3)
            J2 = S(1)*S(2)+S(2)*S(3)+S(1)*S(3) &
                -(S(4)**2)-(S(5)**2)-(S(6)**2)
            J3 = S(1)*S(2)*S(3)+(2*S(4)*S(5)*S(6))-(S(1)*S(5)**2) &
                -(S(2)*S(6)**2)-(S(3)*S(4)**2)
            R = (THIRD*J1**2)-I2
            T = SQRT(TWENTYSEVENTH*R**3)
            Q = (THIRD*J1*J2)-J3-(TWOTWENTYSEVENTH*J1**3)
            ARG = -Q/(TWO*T)
            IF (ARG.GT.1.D0) THEN
            ARG = ARG - 1.E-10
            ELSEIF(ARG.LT.-1.D0) THEN
            ARG = ARG + 1.E-10 
            END IF
            ALP = ACOS(ARG)
            SCALC = SQRT(THIRD*R)

            PRINC1 = (2*SCALC*COS(ALP/THREE))+(THIRD*J1)
            PRINC2 = (2*SCALC*COS((ALP/THREE)+FTPI))+(THIRD*J1)
            PRINC3 = (2*SCALC*COS((ALP/THREE)+TTPI))+(THIRD*J1)

        END IF

        ! C     Assign Principal Stress/Strains values to array
        ! C
        PS(1) = PRINC1
        PS(2) = PRINC2
        PS(3) = PRINC3
        DO K=1,2
            DO K1=K+1,3
                IF(PS(K1)>PS(K)) THEN
                    T1=PS(K);PS(K)=PS(K1);PS(K1)=T1
                ENDIF
            ENDDO
        ENDDO
        
        ! C
        ! C     Calculate cofactors and factor
        ! C


        !S(4:6)=TXY,TYZ,TXZ
        !参考：王凯.主应力方向余弦的计算公式[J].力学与实践,2015,37(03):378-380.

        IA1=0;AN=0.D0
        N1=COUNT(ABS(S(4:6))<1.E-10)
        IF(N1==3) THEN
            DO K=1,3
                RA1=ABS(PS-S(K))                
                N2=MINLOC(RA1,MASK=IA1==0,DIM=1)
                IA1(N2)=1
                AN(K,N2)=1.D0
                     
                !IF(ABS(PS(K)-S(1))<1.E-5.AND.IA1(1)==0) THEN
                !    AN(K,:)=[1.D0,0.D0,0.D0]
                !    IA1(1)=1                        
                !ELSEIF(ABS(PS(K)-S(2))<1.E-5.AND.IA1(2)==0) THEN
                !    AN(K,:)=[0.D0,1.D0,0.D0]
                !    IA1(2)=1
                !ELSEIF(ABS(PS(K)-S(3))<1.E-5.AND.IA1(3)==0) THEN   
                !    AN(K,:)=[0.D0,0.D0,1.D0] 
                !    IA1(3)=1                    
                !ENDIF                  
            ENDDO
        ELSEIF(N1==2) THEN
            N2=MAXLOC(ABS(S(4:6)),DIM=1)            
            SELECT CASE(N2)
            CASE(1) !TXY/=0
                N3=3;L=2;N=1;M=2
            CASE(2) !TYZ/=0
                N3=1;L=2;N=3;M=2
            CASE(3) !TXZ/=0
                N3=2;L=1;N=3;M=1                  
            END SELECT
            
            N4=MINLOC(ABS(PS-S(N3)),DIM=1)
            DO K=1,3
                T1=(PS(K)-S(M))/S(3+N2)
                T2=1.0/(1+T1**2)**0.5
                T3=T1*T2
                
                IF(K==N4) THEN
                    AN(K,:)=0.D0
                    AN(K,N3)=1.D0
                ELSE
                    AN(K,N3)=0.0D0
                    AN(K,L)=T2
                    AN(K,N)=T3
                ENDIF
            ENDDO
        ELSE
        ! **  EQUATIONS FROM: ADVANCED STRENGTH AND APPLIED ELASTICITY: FOURTH EDITION          **
        ! **  ANSEL C. UGURAL, SAUL K. FENSTER PG. 505 APPENDIX B1 AND B2                       **          
            DO K=1, 3
                A(K)=((S(2)-PS(K))*(S(3)-PS(K)))-(S(5)**2)
                B(K)=-(S(4)*(S(3)-PS(K))-S(5)*S(6))
                C(K)=S(4)*S(5)-(S(2)-PS(K))*S(6)
                T1=SQRT(A(K)**2+B(K)**2+C(K)**2)
                IF(ABS(T1)<1.D-10) T1=SIGN(1.D-10,T1)
                V(K)=1./T1
                AN(K,1)=A(K)*V(K)
                AN(K,2)=B(K)*V(K)
                AN(K,3)=C(K)*V(K)         
            END DO
        ENDIF     
        
        
        !对于平面应变问题,通常假定sigma1和sigma3为平面内应力,sigma2为平面外应力。如不是，强制之
        if(NDIMENSION==2) then
            t1=1.e10            
            do K=1,3
                t2=abs((abs(dot_product([0.d0,0.d0,1.d0],an(K,:))))**0.5-1.0d0)
                if(t2<t1) then
                    t1=t2;n1=K
                endif
            enddo
            if(n1/=2) then
                t2=PS(2);V(1:3)=an(2,:)
                PS(2)=PS(n1);an(2,:)=an(n1,:)
                PS(n1)=t2;an(n1,:)=V(1:3)
            endif
        endif        

        RETURN
    END
   
     
END MODULE