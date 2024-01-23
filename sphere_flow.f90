
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
    
    if(.not.IS_WELL_GEO_INI) CALL WELL_GEO_INI(.false.)
    
        
    N1=ELEMENT(IELT).NODE(2)
    TPHI1=0.D0;HK1=0.D0
    DV1(1:NDIMENSION)=ELEMENT(IELT).PROPERTY(4:6)
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
            !centroid
            DO I=1,3
                XV1(I)=SUM(NODE(ELEMENT(IELT1).NODE).COORD(I))
            ENDDO
            
            XV1=XV1/ELEMENT(IELT1).NNUM
            !print *, j,xv1
            XV1=XV1-NODE(N1).COORD
            T1=DOT_PRODUCT(XV1(1:NDIMENSION),DV1(1:NDIMENSION))/ &
            (NORM2(XV1(1:NDIMENSION))*NORM2(DV1(1:NDIMENSION)))
            IF(T1<0.D0) CYCLE
        ENDIF
        PHI1=ELEMENT(IELT1).ANGLE(SNADJL(N1).SUBID(J)) 
        TPHI1=TPHI1+PHI1
        HK1=HK1+PHI1*(MATERIAL(ELEMENT(IELT1).MAT).PROPERTY(1)*MATERIAL(ELEMENT(IELT1).MAT).PROPERTY(2))**0.5 
              
    ENDDO
    HK1=HK1/TPHI1 !平均的K
    !HK1=1
    !统计边长
    
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
    !参考:毛昶熙,渗流计算分析与控制,第二版.P340
    T1=2*PI1*HK1/(1/WR1-(4*PI1+3)/(3*L1))
    IF(SOLVER_CONTROL.WELL_BOTTOM_TYPE==0) T1=T1/PI1*2
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
    DV1(1:NDIMENSION)=ELEMENT(IELT).PROPERTY(4:6)
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
    DV1(1:NDIMENSION)=ELEMENT(IELT).PROPERTY(4:6)
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
    DV1(1:NDIMENSION)=ELEMENT(IELT).PROPERTY(4:6)
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

SUBROUTINE sphere_flow_element(ielt)
    
    USE MESHADJ,only:SNADJL
    USE solverds

    IMPLICIT NONE
    INTEGER,INTENT(IN)::IELT
    INTEGER::I,J,IELT1,N1,N2,N3,K
    REAL(8)::XV1(3),T1,HK1,PHI1,TPHI1,WR1,L1,DV1(3),PI1,SK1
    REAL(8)::RDIS2(400)
    
    if(.not.IS_WELL_GEO_INI) CALL WELL_GEO_INI(.true.)
    
        
    N1=ELEMENT(IELT).NODE(2)
    TPHI1=0.D0;HK1=0.D0
    DV1(1:NDIMENSION)=ELEMENT(IELT).PROPERTY(4:6)
    WR1=MATERIAL(ELEMENT(IELT).MAT).PROPERTY(1)
    SK1=0
    DO J=1,SNADJL(N1).ENUM        
        IELT1=SNADJL(N1).ELEMENT(J)
        IF(ELEMENT(IELT1).ET==WELLBORE &
        .OR.ELEMENT(IELT1).ET==PIPE2 &
        .OR.ELEMENT(IELT1).ET==SPHFLOW &
        .OR.ELEMENT(IELT1).ET==SEMI_SPHFLOW &
        .OR.ELEMENT(IELT1).ET==WELLBORE_SPGFACE) CYCLE !ITSELF EXCLUDED
        
        IF(ELEMENT(IELT1).EC/=SPG) CYCLE
        IF(.NOT.ALLOCATED(ELEMENT(IELT1).ANGLE)) CALL calangle(IELT1)

        !AVERAGED K, ASSUME ISOTROPIC        
        IF(ELEMENT(IELT).ET==SEMI_SPHFLOW) THEN
            !centroid
            DO I=1,3
                XV1(I)=SUM(NODE(ELEMENT(IELT1).NODE).COORD(I))
            ENDDO
            
            XV1=XV1/ELEMENT(IELT1).NNUM
            !print *, j,xv1
            XV1=XV1-NODE(N1).COORD
            T1=DOT_PRODUCT(XV1(1:NDIMENSION),DV1(1:NDIMENSION))/ &
            (NORM2(XV1(1:NDIMENSION))*NORM2(DV1(1:NDIMENSION)))
            IF(T1<0.D0) CYCLE
        ENDIF
        PHI1=ELEMENT(IELT1).ANGLE(SNADJL(N1).SUBID(J)) 
        TPHI1=TPHI1+PHI1
        HK1=HK1+PHI1*(MATERIAL(ELEMENT(IELT1).MAT).PROPERTY(1)*MATERIAL(ELEMENT(IELT1).MAT).PROPERTY(2))**0.5 
        call cal_element_k(IELT1,SNADJL(N1).SUBID(J),WR1)
        SK1=SK1+ELEMENT(IELT1).FD     
    ENDDO
    HK1=HK1/TPHI1 !平均的K
    !HK1=1
    !统计边长
    SK1=(TPHI1)**2*WR1/SK1
    !SK1=(TPHI1)*WR1/SK1 
    !IF(solver_control.well_bottom_type==0)  SK1=SK1*0.95 
    
    !参考:毛昶熙,渗流计算分析与控制,第二版.P340
    !PI1=4*ATAN(1.D0);L1=1
    !T1=2*PI1*HK1/(1/WR1-(4*PI1+3)/(3*L1))
    !SK1=T1
    
    ELEMENT(IELT).KM(1,1)=1.D0;ELEMENT(IELT).KM(2,2)=1.D0
    ELEMENT(IELT).KM(2,1)=-1.D0;ELEMENT(IELT).KM(1,2)=-1.D0
    ELEMENT(IELT).KM=ELEMENT(IELT).KM*SK1
    ELEMENT(IELT).PROPERTY(1)=1/SK1
    ELEMENT(IELT).FD=SK1 !FD此处为导水系数
    
ENDSUBROUTINE

subroutine cal_element_k(ielt,inwell,rw)
!目前仅适用于各向同性的土体
    use solverds
    use SolverMath
    implicit none
    integer,intent(in)::ielt,inwell
    real(8),intent(in)::rw
    integer::node1(4),i,j
    real(8)::xy1(3,4),cofactor1(4,4),vol1,beta1(4),beta,alpha1,Sb,Si,K1(3),Ki,XSCALE1,t1,r1, &
                        XC1(3),cos1,sm,sp,grad1

    select case(inwell)
    case(1)
        node1=[1,2,3,4]        
    case(2)
        node1=[2,3,1,4]
    case(3)
        node1=[3,4,1,2]
    case(4)
        node1=[4,1,3,2]
    end select
    K1=MATERIAL(ELEMENT(IELT).MAT).PROPERTY(1:3)
    XSCALE1=(K1(3)/K1(1))**0.5 !!ASSUME KX=KY.
    !水平各向异性的处理,kx=ky/=kz,按沙金煊调整x和y轴的方法进行处理(沙金煊. 各向异性土渗流的转化问题[J]. 水利水运科学研究, 1987, (01): 15-28.)
    Ki=(K1(1)*K1(2))**0.5

    XY1(:,1)=0.D0
    do i=2,4
        xy1(:,i)=node(element(ielt).node(node1(i))).coord-node(element(ielt).node(node1(1))).coord
        beta1(i)=((xy1(1,i)*xscale1)**2 + &
            (xy1(2,i)*xscale1)**2 + &
            xy1(3,i)**2)**0.5
        t1=beta1(i)/norm2(xy1(:,i))
        !beta1(i)=(beta1(i)-t1*rw)/beta1(i)
        ! if(solver_control.well_bottom_type==0) then
        !     beta1(i)=pi()/2.0*beta1(i)
        ! endif 
        if(solver_control.well_bottom_type/=0) then        
            beta1(i)=(beta1(i)-t1*rw)/beta1(i) 
        ELSE
            r1=((xy1(1,i)*xscale1)**2 + (xy1(2,i)*xscale1)**2)**0.5
            if(r1>1.d-7) then
                beta1(i)=Pi()/2.0-asin(((xy1(3,i)**2+(r1+rw)**2)**0.5-(xy1(3,i)**2+(r1-rw)**2)**0.5)/(2*r1))
            else
                beta1(i)=Pi()/2.0-asin(rw/(rw**2+xy1(3,i)**2)**0.5)
            endif
            !if(i==4) beta1(2:4)=sum(beta1(2:4))/3.0
        ENDIF
    enddo
     
    do j=1,4        
        do i=1,4
            cofactor1(i,j)=tet_shape_factor(xy1,j*10+i)
        enddo
    enddo
    vol1=tet_shape_factor(xy1,0)
    Sb= &
    beta1(2)*(K1(1)*cofactor1(1,2)*cofactor1(2,2) + K1(2)*cofactor1(1,3)*cofactor1(2,3) + K1(3)*cofactor1(1,4)*cofactor1(2,4)) + &
    beta1(3)*(K1(1)*cofactor1(1,2)*cofactor1(3,2) + K1(2)*cofactor1(1,3)*cofactor1(3,3) + K1(3)*cofactor1(1,4)*cofactor1(3,4)) + &
    beta1(4)*(K1(1)*cofactor1(1,2)*cofactor1(4,2) + K1(2)*cofactor1(1,3)*cofactor1(4,3) + K1(3)*cofactor1(1,4)*cofactor1(4,4)) 
    Sb= -Sb
    Si = K1(1)*cofactor1(1,2)**2 + K1(2)*cofactor1(1,3)**2 + K1(3)*cofactor1(1,4)**2 
    alpha1=element(ielt).angle(inwell)
    element(ielt).fd=alpha1*(Sb-36*Ki*Vol1*rw*alpha1)/(Ki*Si) 
    
    !if(solver_control.well_bottom_type==0) then
    !    DO I=1,3
    !        XC1(I)=SUM(XY1(I,2:4))/3.0            
    !    ENDDO
    !    COS1=ABS(XC1(3))/NORM2(XC1(1:3))
    !    Sp=sqrt(2.*(1. + cos1));Sm=sqrt(2.*(1. - cos1));
    !    grad1= ((Sm + Sp)*cos1 + Sm - Sp)/(sqrt(2.*(cos1**2-1.) + Sp*Sm )*Sm*Sp)
    !    alpha1=alpha1*abs(grad1)/(2*pi())
    !else
    !    alpha1=alpha1/(2*pi())
    !endif  
    !element(ielt).fd=alpha1*(Sb-72*Ki*PI()*Vol1*rw*alpha1)/(Ki*Si)
    
    
endsubroutine 
    
SUBROUTINE SPHFLOW_Q_K_UPDATE_ANALYTICAL(STEPDIS,IELT,ISTEP,IITER)
    
    USE solverds
    USE MESHGEO,ONLY:SNADJL
    IMPLICIT NONE
    INTEGER,INTENT(IN)::IELT,ISTEP,IITER
    REAL(8),INTENT(IN)::STEPDIS(NDOF)
    REAL(8)::T1,KM1(2,2),NHEAD1(2)
    REAL(8)::QN1(2)
   
        
    NHEAD1=STEPDIS(ELEMENT(IELT).G)
    QN1(1:2)=MATMUL(ELEMENT(IELT).KM,NHEAD1(1:2))
    QWAN(2,ELEMENT(IELT).NODE(1:2))=QWAN(2,ELEMENT(IELT).NODE(1:2))+QN1(1:2)    
    
    
    IF(ANY(NHEAD1-NODE(ELEMENT(IELT).NODE).COORD(NDIMENSION)<0.D0)) THEN
        T1=1.D-7
    ELSE        
        T1=ELEMENT(IELT).FD
       
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