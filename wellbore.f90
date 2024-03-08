SUBROUTINE WELL_GEO_INI(IS_ONLY_SNADJL)
    USE MESHADJ,ONLY:SEDGE,NSEDGE,SNADJL,GETGMSHET,ELTTYPE, &                        
                        Setup_Solver_MESHTOPO,SETUP_NODE2ELEMENT_ADJL_SOLVER
    USE solverds,ONLY:isIniSEdge,NODE,ELEMENT,NDIMENSION,NQWNODE,QWELLNODE,WELLBORE,WELLBORE_SPGFACE,PIPE2,POREFLOW,ENUM,QWAN,NNUM, &
                IS_WELL_GEO_INI,spg,material
    IMPLICIT NONE
    LOGICAL,INTENT(IN)::IS_ONLY_SNADJL
    INTEGER I,J,K,N1,N2,N3,N4,A1(2),A2(200),v1,v2
    INTEGER,ALLOCATABLE::IELT1(:),ielt2(:)
    REAL(8)::XYLMT(2,3),wr1,area1
    logical::isp1
    
    IF(.NOT.isIniSEdge) THEN
        if(.not.IS_ONLY_SNADJL) then
            call Setup_Solver_MESHTOPO()
        else
            call SETUP_NODE2ELEMENT_ADJL_SOLVER(SNADJL)
        endif
    ENDIF

    IF(IS_WELL_GEO_INI) RETURN

    A1=[4,3]    
    IF(.NOT.ALLOCATED(IELT1)) ALLOCATE(IELT1(ENUM))
    IELT1=0
    ielt2=ielt1
    DO I=1,NQWNODE
        N3=0
        DO K=1,SNADJL(QWELLNODE(I).NODE(1)).ENUM !如果井口节点位于中间，上下各搜索一次
            N1=QWELLNODE(I).NODE(1)
10          DO J=1,SNADJL(N1).ENUM
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
    
    !借用sign
    element.sign=0
    !initialize weights
    !!假定井足够远，每个单元只属于一个井
    !do i=1,enum
    !    if(element(i).et==wellbore.or.element(i).et==WELLBORE_SPGFACE) then
    !        n1=abs(element(i).edge(3)) !命名iedge为第三边3-4
    !        wr1=material(element(i).mat).property(1)
    !        !统计以井流单元为边的单元数
    !        do j=1,sedge(n1).enum
    !            n2=sedge(n1).element(j)
    !            n3=getgmshet(element(n2).et)
    !            if(elttype(n3).dim==3.and.element(n2).ec==spg) then
    !                !v1=minloc(abs(element(n2).node-sedge(n1).v(1)),dim=1)
    !                !v2=minloc(abs(element(n2).node-element(i).node(4)),dim=1)
    !                element(n2).property(6)=1.0d0 !weights
    !                !call wellbore_cal_element_k(n2,v1,wr1,sedge(n1).v(1),sedge(n1).v(2),area1,isp1)
    !                
    !            
    !                element(n2).sign=-i !井轴附属单元(2个节点在井轴上的四面体单元)
    !                
    !                !if(.not.allocated(element(n2).angle)) call calangle(n2)
    !            endif
    !        enddo
    !        
    !        do j=1,2
    !            do k=1,snadjl(sedge(n1).v(j)).enum
    !                n2=snadjl(sedge(n1).v(j)).element(k)
    !                if(element(n2).sign<0) cycle !跳过井轴附属单元
    !                n3=getgmshet(element(n2).et)
    !                if(elttype(n3).dim==3.and.element(n2).ec==spg) then
    !                    element(n2).sign=element(n2).sign+1
    !                    element(n2).property(6)=1.0d0/element(n2).sign
    !                endif
    !            enddo
    !        enddo
    !        
    !        
    !        
    !    endif
    !enddo    
    !!total weights
    !do i=1,enum
    !    if(element(i).et==wellbore.or.element(i).et==WELLBORE_SPGFACE) then
    !        n1=abs(element(i).edge(3)) !命名iedge为第三边3-4
    !        element(i).ww=0.0d0
    !        do j=1,2
    !            do k=1,snadjl(sedge(n1).v(j)).enum
    !                n2=snadjl(sedge(n1).v(j)).element(k)
    !                if(element(n2).sign==0) cycle
    !                if(element(n2).sign<0.and.element(n2).sign/=-i) cycle 
    !                if(j==2.and.element(n2).sign==-i) cycle !只计算一次
    !                element(i).ww=element(i).ww+element(n2).property(6)
    !            enddo
    !        enddo
    !    endif
    !enddo     
    
    
    do i=1,enum
        if(element(i).et==wellbore.or.element(i).et==WELLBORE_SPGFACE) then
            n1=abs(element(i).edge(3)) !命名iedge为第三边3-4
            wr1=material(element(i).mat).property(1)
            !统计以井流单元为边的单元数
            do j=1,sedge(n1).enum
                n2=sedge(n1).element(j)
                n3=getgmshet(element(n2).et)
                if(elttype(n3).dim==3.and.element(n2).ec==spg) then
                    v1=minloc(abs(element(n2).node-sedge(n1).v(1)),dim=1)
                    !v2=minloc(abs(element(n2).node-element(i).node(4)),dim=1)
                    !element(n2).property(6)=1.0 !weights
                    call wellbore_cal_element_k(n2,v1,wr1,sedge(n1).v(1),sedge(n1).v(2),area1,isp1)
                    element(n2).property(6)=area1
                
                    element(n2).sign=-i !井轴附属单元(2个节点在井轴上的四面体单元)
                    
                    if(.not.allocated(element(n2).angle)) call calangle(n2)
                endif
            enddo
        endif    
    enddo
           
    IS_WELL_GEO_INI=.TRUE.

    IF(ALLOCATED(IELT1)) DEALLOCATE(IELT1)
    IF(ALLOCATED(IELT2)) DEALLOCATE(IELT2)

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
    
    if(.not.IS_WELL_GEO_INI) CALL WELL_GEO_INI(.false.)
    
    
    
    IEDGE1=ABS(ELEMENT(IELT).EDGE(3))
    !L1=SEDGE(IEDGE1).DIS/2.0
    L1=NORM2(NODE(SEDGE(IEDGE1).V(1)).COORD-NODE(SEDGE(IEDGE1).V(2)).COORD)/2.0
    IF(ELEMENT(IELT).ISTOPO>0) IEDGE1=ABS(ELEMENT(IELT).EDGE(5))
    N1=0
    WR1=MATERIAL(ELEMENT(IELT).MAT).PROPERTY(1)

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
            !call wellbore_area(IELT1,SEDGE(IEDGE1).V(I),SEDGE(IEDGE1).V(mod(i,2)+1),wr1)
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
    
   
 subroutine wellbore_element(ielt)
    
    use meshadj,only:sedge,nsedge,snadjl,getgmshet,elttype,sface,nsface
    use solverds
    use solvermath
    use ifport
    implicit none
    integer,intent(in)::ielt
    integer::i,j,k,ielt1,iedge1,gmet1,subid1,af1(2),n1,nadjl1(200),n2,n3,k1,n4
    real(8)::d1,t1,t2,hk(2),tphi,wr1,area,sk,L1
    logical::isp1
    
    if(.not.IS_WELL_GEO_INI) CALL WELL_GEO_INI(.false.)
    iedge1=abs(element(ielt).edge(3)) !命名iedge为第三边3-4
    if(element(ielt).istopo>0) iedge1=abs(element(ielt).edge(5))
    n1=0
    wr1=material(element(ielt).mat).property(1)
    L1=NORM2(NODE(SEDGE(IEDGE1).V(1)).COORD-NODE(SEDGE(IEDGE1).V(2)).COORD)
    if(sedge(iedge1).enum<2) then
        print *, "no elements surrounding the wellbore element i. i=",ielt
        print *, "try to embed the well line into the correspoinding volume."
        stop
    endif

    tphi=0
    sk=0
   
    !open(66,file='test.txt',status='replace')
    !借用sign
    !where(element.sign>0)  element.sign=0 
    n1=0
    do i=1,2
        do j=1,snadjl(sedge(iedge1).v(i)).enum
            ielt1=snadjl(sedge(iedge1).v(i)).element(j)!与第三边3-4邻接的单元
            if(ielt1<1)  cycle
            if(element(ielt1).et==wellbore.or.element(ielt1).et==pipe2 &
            .or.element(ielt1).et==wellbore_spgface) cycle !itself excluded
            if(element(ielt1).et==sphflow.or.element(ielt1).et==semi_sphflow) cycle !itself excluded
            if(element(ielt1).ec/=spg) cycle
            
            !if(element(ielt1).sign>0) cycle
            
            !if(element(ielt1).sign/=0.and.element(ielt1).sign/=ielt) cycle !井轴单元(2个节点在井轴上的单元)
            !if(i==2.and.element(ielt1).sign==ielt) cycle !每个单元只计算一次            
            
            
            if(element(ielt1).sign==0) then
                call wellbore_cal_element_k(ielt1,snadjl(sedge(iedge1).v(i)).subid(j),wr1,sedge(iedge1).v(i),sedge(iedge1).v(mod(i,2)+1),area,isp1)
                element(ielt1).property(6)=area
                
                if(isp1) then 
                    !should not be here
                    print *, 'unexpected error in sub=wellbore_element' !在sub=WELL_GEO_INI已经计算
                    stop
                    !element(ielt1).sign=2
                else
                    element(ielt1).sign=1
                endif
                
            elseif(element(ielt1).sign==-ielt) then
                if(i==2) cycle !只考虑一次
            elseif(element(ielt1).sign<0.and.element(ielt1).sign/=-ielt) then
                cycle !跳过其他井轴上的单元
            endif
            
            
            
            n1=n1+1
            
            !tphi=tphi+area
            tphi=tphi+element(ielt1).property(6)
            sk=sk+element(ielt1).fd   
        
            
                !write(66,*) i,'tphi=', tphi,'area=', area,'sk=', sk,'fd=',element(ielt1).fd
                !write(66,*)'-------------------------------------------------------------------------------' 
                
        enddo    
    enddo
 
    
    !sk=tphi**2.d0/sk!单刚
    !sk(2)=tphi(2)**2/sk(2)
    !sk(2)=sk(1)
    element(ielt).property(2:3)=2.0*sk/(tphi**2)  !阻力,并联阻力加倍
    !element(ielt).property(2:3)=2.0*sk/(2*pi()*wr1*L1*n1) 
    !井损
    IF(ABS(MATERIAL(ELEMENT(IELT).MAT).PROPERTY(7))>1E-12) THEN
        t1=2.0*PI()*WR1*L1/2.0
        ELEMENT(IELT).PROPERTY(5)=MATERIAL(ELEMENT(IELT).MAT).PROPERTY(7)/t1
        ELEMENT(IELT).PROPERTY(2:3)=ELEMENT(IELT).PROPERTY(2:3)+ELEMENT(IELT).PROPERTY(5)
        !ELEMENT(IELT).PROPERTY(2:3)=ELEMENT(IELT).PROPERTY(2:3)
    ELSE
        ELEMENT(IELT).PROPERTY(5)=0.D0
    ENDIF
    
    !element(ielt).property(3)=1/sk(2)
    element(ielt).km=km_wellbore(1./element(ielt).property(1),1./element(ielt).property(2),1./element(ielt).property(3))!单刚
    endsubroutine   

        
!    subroutine wellbore_cal_element_k2(ielt,inwell,rw,v1,v2,wq,isp1)
!!目前仅适用于各向同性的土体
!    use solverds
!    use solvermath
!    implicit none
!    integer,intent(in)::ielt,inwell,v1,v2
!    real(8),intent(in)::rw
!    integer::node1(4),i,j
!    real(8)::cofactor1(4,4),vol1,r(4),beta,sb,si,k1(3),ki,xscale1,t1,r1,mi,wq,xy1(3,4)
!    logical::isp1,iscalbyK1
!    real(8)::km1(element(ielt).ndof,element(ielt).ndof),ar1(1:10),t2,t3,alpha1(4)
!    integer::n1,ia1(2)
!    
!    !call wellbore_area(ielt,v1,v2,rw,area,xy1,isp1)
!    !select case(inwell)
!    !case(1)
!    !    node1=[1,2,3,4]        
!    !case(2)
!    !    node1=[2,3,1,4]
!    !case(3)
!    !    node1=[3,4,1,2]
!    !case(4)
!    !    node1=[4,1,3,2]
!    !end select
!    k1=material(element(ielt).mat).property(1:3)
!    xscale1=(k1(3)/k1(1))**0.5 !!assume kx=ky.
!    !水平各向异性的处理,kx=ky/=kz,按沙金煊调整x和y轴的方法进行处理(沙金煊. 各向异性土渗流的转化问题[j]. 水利水运科学研究, 1987, (01): 15-28.)
!    ki=(k1(1)*k1(2))**0.5
!
!    
!    do i=2,4!计算每个边的模长r
!        r(i)=((xy1(1,i)*xscale1)**2 + (xy1(2,i)*xscale1)**2)**0.5
!        if(r(i)==0) then
!            r(i)=rw
!        end if
!        
!    enddo
!   
!    km1=element(ielt).km
!    ia1=[v1,v2]   
!    alpha1(2:4)=-rw/(ki)*log(r(2:4)/rw) !!! -,解析解的流量正负与数值解相反????
!    do i=1,2
!        n1=minloc(abs(element(ielt).node(1:element(ielt).nnum)-ia1(i)),dim=1)
!        if(n1/=i) then
!            !交换列            
!            ar1(:element(ielt).ndof)=km1(:,i)
!            km1(:,i)=km1(:,n1)
!            km1(:,n1)=ar1(:element(ielt).ndof)
!            !交换行
!            ar1(:element(ielt).ndof)=km1(i,:)
!            km1(i,:)=km1(n1,:)
!            km1(n1,:)=ar1(:element(ielt).ndof)
!        endif
!        if(.not.isp1) exit
!    enddo
!        
!    if(isp1) then !2个节点都在井轴上
!        !area=area/2.0
!        !t1=(area-dot_product(km1(1,3:4),r(3:4))) 
!        !t2=km1(1,1)+Km1(1,2) 
!        !if(abs(t2)<1.d-8)  t2=km1(1,1)*1.d-4
!        t1=area-(km1(1,3)+km1(2,3))*alpha1(3)-(km1(1,4)+km1(2,4))*alpha1(4)
!        t2=km1(1,1)+Km1(1,2)+km1(2,1)+Km1(2,2)
!        !t3=1.0d0
!        !area=area/2.0
!    else
!        !if(v1==3) area=area/2.0
!        t1=(area-dot_product(km1(1,2:4),alpha1(2:4)))
!        t2=km1(1,1)
!        !t3=1.d0
!    endif
!    if(abs(t2)>1.d-8) then 
!        !element(ielt).fd=-t3*area*(t1/t2)
!        element(ielt).fd=-(t1/t2)
!    else
!        element(ielt).fd=0.d0
!    endif
!
!    !write(66,*)'r(1-4)=',r,'v=', vol1,'bi=',cofactor1(1,2),'sb=',sb,'si=',si
!endsubroutine     
    
    subroutine wellbore_cal_element_k(ielt,inwell,rw,v1,v2,area,isp1)
!目前仅适用于各向同性的土体
    use solverds
    use solvermath
    implicit none
    integer,intent(in)::ielt,inwell,v1,v2
    real(8),intent(in)::rw
    integer::node1(4),i,j
    real(8)::cofactor1(4,4),vol1,r(4),beta,sb,si,k1(3),ki,xscale1,t1,r1,mi,area,xy1(3,4)
    logical::isp1,iscalbyK1
    real(8)::km1(element(ielt).ndof,element(ielt).ndof),ar1(1:10),t2,t3,alpha1(4)
    integer::n1,ia1(2)
    
    call wellbore_area(ielt,v1,v2,rw,area,xy1,isp1)
    !select case(inwell)
    !case(1)
    !    node1=[1,2,3,4]        
    !case(2)
    !    node1=[2,3,1,4]
    !case(3)
    !    node1=[3,4,1,2]
    !case(4)
    !    node1=[4,1,3,2]
    !end select
    k1=material(element(ielt).mat).property(1:3)
    xscale1=(k1(3)/k1(1))**0.5 !!assume kx=ky.
    !水平各向异性的处理,kx=ky/=kz,按沙金煊调整x和y轴的方法进行处理(沙金煊. 各向异性土渗流的转化问题[j]. 水利水运科学研究, 1987, (01): 15-28.)
    ki=(k1(1)*k1(2))**0.5

    
    do i=2,4!计算每个边的模长r
        r(i)=((xy1(1,i)*xscale1)**2 + (xy1(2,i)*xscale1)**2)**0.5
        if(r(i)==0) then
            r(i)=rw
        end if
        
    enddo
    iscalbyK1=.true.
    if(iscalbyK1) then
        km1=element(ielt).km
        ia1=[v1,v2]   
        alpha1(2:4)=-rw/(ki)*log(r(2:4)/rw) !!! -,解析解的流量正负与数值解相反????
        do i=1,2
            n1=minloc(abs(element(ielt).node(1:element(ielt).nnum)-ia1(i)),dim=1)
            if(n1/=i) then
                !交换列            
                ar1(:element(ielt).ndof)=km1(:,i)
                km1(:,i)=km1(:,n1)
                km1(:,n1)=ar1(:element(ielt).ndof)
                !交换行
                ar1(:element(ielt).ndof)=km1(i,:)
                km1(i,:)=km1(n1,:)
                km1(n1,:)=ar1(:element(ielt).ndof)
            endif
            if(.not.isp1) exit
        enddo
        
        if(isp1) then !2个节点都在井轴上
            !area=area/2.0
            !t1=(area-dot_product(km1(1,3:4),r(3:4))) 
            !t2=km1(1,1)+Km1(1,2) 
            !if(abs(t2)<1.d-8)  t2=km1(1,1)*1.d-4
            t1=(area-(km1(1,3)+km1(2,3))*alpha1(3)-(km1(1,4)+km1(2,4))*alpha1(4))*area
            t2=km1(1,1)+Km1(1,2)+km1(2,1)+Km1(2,2)
            !t3=1.0d0
            !area=area/2.0
        else
            !if(v1==3) area=area/2.0
            t1=(area-dot_product(km1(1,2:4),alpha1(2:4)))*area
            t2=km1(1,1)
            !t3=1.d0
        endif
        if(abs(t2)>1.d-8) then 
            !element(ielt).fd=-t3*area*(t1/t2)
            element(ielt).fd=-(t1/t2)
        else
            element(ielt).fd=0.d0
        endif
    else
        do j=1,4        
            do i=1,4
                cofactor1(i,j)=tet_shape_factor(xy1,j*10+i)
            enddo
        enddo
        vol1=tet_shape_factor(xy1,0)
        if (isp1==.true.) then 
            sb= &
            log(r(3)/rw)*(k1(1)*(cofactor1(1,2)+cofactor1(2,2))*cofactor1(3,2) + k1(2)*(cofactor1(1,3)+cofactor1(2,3))*cofactor1(3,3) + k1(3)*(cofactor1(1,4)+cofactor1(2,4))*cofactor1(3,4)) + &
            log(r(4)/rw)*(k1(1)*(cofactor1(1,2)+cofactor1(2,2))*cofactor1(4,2) + k1(2)*(cofactor1(1,3)+cofactor1(2,3))*cofactor1(4,3) + k1(3)*(cofactor1(1,4)+cofactor1(2,4))*cofactor1(4,4)) 
            sb= -sb
            si = k1(1)*(cofactor1(1,2)+cofactor1(2,2))**2 + k1(2)*(cofactor1(1,3)+cofactor1(2,3))**2 + k1(3)*(cofactor1(1,4)+cofactor1(2,4))**2
            mi=norm2(node(v1).coord-node(v2).coord)            
            element(ielt).fd=area*(sb*rw-36*area*ki*vol1)/(ki*si) 
            !area=area/2.0
        else
            sb= &
            log(r(2)/rw)*(k1(1)*cofactor1(1,2)*cofactor1(2,2) + k1(2)*cofactor1(1,3)*cofactor1(2,3) + k1(3)*cofactor1(1,4)*cofactor1(2,4)) + &
            log(r(3)/rw)*(k1(1)*cofactor1(1,2)*cofactor1(3,2) + k1(2)*cofactor1(1,3)*cofactor1(3,3) + k1(3)*cofactor1(1,4)*cofactor1(3,4)) + &
            log(r(4)/rw)*(k1(1)*cofactor1(1,2)*cofactor1(4,2) + k1(2)*cofactor1(1,3)*cofactor1(4,3) + k1(3)*cofactor1(1,4)*cofactor1(4,4)) 
            sb= -sb
            si = k1(1)*cofactor1(1,2)**2 + k1(2)*cofactor1(1,3)**2 + k1(3)*cofactor1(1,4)**2 
            mi=norm2(node(v1).coord-node(v2).coord)
            element(ielt).fd=area*(sb*rw-36*area*ki*vol1)/(ki*si) 
        end if
    endif
    
    !write(66,*)'r(1-4)=',r,'v=', vol1,'bi=',cofactor1(1,2),'sb=',sb,'si=',si
endsubroutine     
   
    subroutine wellbore_area(ielt,v1,v2,wr,area,xy1,isp1)
    !求以v1为其中一个节点的四面体单元ielt的面与半径为wr,轴线为v1-v2的圆柱面的交线所围成的柱面面积(v1-v2范围内的面积)。
    !假定：1)单元ielt的其中一个(v1)或2个(v1,v2)节点落在圆柱轴线v1-v2上.
    !算法：1）先求出各边与圆柱面的交点，然后数值积分，数值积分区域的边界由三角面与圆柱面的交线确定
        use solverds,only:node,element
        use geometricalgorithm
        implicit none
        integer,intent(in)::ielt,v1,v2
        real(8),intent(in)::wr
        real(8),intent(out)::area,xy1(3,4)
        real(8)::axis(3),ip1(3,4),g2l1(3,3),org1(3),t1,t2,t3,s0,s1,xy(3,4),p1(3),p2(3),p3(3),p4(3),p5(3),d1,d2,r(3,3),a(3,3),b(3,3),x1,y1,z1
       integer::i,j,k,node1(4),i1
        logical::isp1
      
      

        node1=pack([1:4],element(ielt).node(1:4)==v1,[0,0,0,0])!出来的是[v1,v1,v1,v1]
        node1(2:4)=pack([1:4],element(ielt).node(1:4)/=v1)!然后再把后三个覆盖掉
        isp1=.false.
        do i=1,4
            xy(:,i)=node(element(ielt).node(i)).coord
        if(all( xy(:,i)==node(v2).coord)) then!判断是否有列与v2相同
            isp1=.true.!棱在轴线上判断isp1真
             node1(2:4)=i!把v2放第二位
             node1(3:4)=pack([1:4],element(ielt).node(1:4)/=v1.and.element(ielt).node(1:4)/=v2)!把其余两点放第三四位，此时1是v1位置，i是v2位置
        endif
        enddo
       
!write(66,*)'xy=',xy
!write(66,*)'v1=',node(v1).coord
!write(66,*) 'v2=',node(v2).coord
!write(66,*) 'node1=',node1
!write(66,*)' ' 
        axis=node(v2).coord-node(v1).coord!轴线向量
        do i=1,4
            xy1(:,i)=xy(:,node1(i))-node(v1).coord!以v1为原点用三行四列矩阵存储各个点相对值,排第一是v1，第二是v2
        enddo
!write(66,*)'xy1=',xy1
        !!set up local system
        !org1=node(v1).coord               
        !call gen_cordinate_system(org1,g2l1,zv=axis)
        !axis=matmul(g2l1,axis)
        !do i=2,4
        !   xy1(:,i)=matmul(g2l1,xy1(:,i)) !将坐标标准化
        !   if(isp1.and.i==2) cycle
        !   !若还有点在轴线上跳过该点
        !   t1=norm2(xy1(1:2,i))!2次范数求与原点直线距离
        !   if(t1<wr) then
        !       print *, 'node(i) is inside the wellbore,i=',element(ielt).node(node1(i))
        !       !error stop
        !   endif
        !end do  !以上循环只是将坐标标准化顺便判断是否有点在圆柱外（准备工作已完成）
        !s1=  xy1(1,1)
        if (axis(2)==0.and.axis(3)==0) then!向量就在x轴上，绕y轴旋转90°
        r = reshape((/0.0, 0.0, 1.0, 0.0, 1.0, 0.0,-1.0,0.0, 0.0/),(/3, 3/))
    else
        x1=axis(1)
        y1=axis(2)
        z1=axis(3)
   a = reshape((/1.0_8, 0.0_8, 0.0_8, 0.0_8, z1/sqrt(y1**2+z1**2), -y1/sqrt(y1**2+z1**2),0.0_8, y1/sqrt(y1**2+z1**2), z1/sqrt(y1**2+z1**2)/),(/3, 3/))
    b=reshape((/sqrt(y1**2+z1**2)/sqrt(x1**2+y1**2+z1**2), -x1/sqrt(x1**2+y1**2+z1**2), 0.0_8, x1/sqrt(x1**2+y1**2+z1**2), sqrt(y1**2+z1**2)/sqrt(x1**2+y1**2+z1**2), 0.0_8, 0.0_8, 0.0_8, 1.0_8/),(/3, 3/))
    r=MATMUL(a, b)!如果点不在x轴上则绕xz轴旋转至z轴
        endif
        axis=matmul(r,axis)
        do i=2,4
           xy1(:,i)=matmul(r,xy1(:,i)) !将坐标标准化
           if(isp1.and.i==2) cycle
           !若还有点在轴线上跳过该点
           t1=norm2(xy1(1:2,i))!2次范数求与原点直线距离
           if(t1<wr) then
               print *, 'node(i) is inside the wellbore,i=',element(ielt).node(node1(i))
               !error stop
           endif
        end do  !以上循环只是将坐标标准化顺便判断是否有点在圆柱外（准备工作已完成）
 
  !write(66,*)'xy1=',xy1
  !write(66,*)' ' 
       
  if(isp1==.false.)then!开始计算点在轴线上情况的面积
        do i=2,4
            ip1(:,i-1)=xy1(:,i)*wr/norm2(xy1(1:2,i))
            ip1(2,i-1)=atan2(ip1(2,i-1),ip1(1,i-1))
            ip1(1,i-1)=wr
             !坐标系变换为柱坐标只需要除二次范数乘半径即可保证点落在圆柱面上（相似三角形），第一行r第二行角度（要替换）第三行高度（不用替换）
        end do
            i=minloc(ip1(2,1:3),dim=1) !找角最小的位置赋值给i
            k=maxloc(ip1(2,1:3),dim=1) !找角最大的位置赋值给k
            do i1=1,3
                if(i1/=i.and.i1/=k) then
                    j=i1!赋值j为角度中间值位置
                    exit
                endif
            enddo
            !o-i-j plane equation
            p1=normal_triface(reshape([0.d0,0.d0,0.d0,ip1(1,i)*cos(ip1(2,i)),ip1(1,i)*sin(ip1(2,i)),ip1(3,i),&
            ip1(1,j)*cos(ip1(2,j)),ip1(1,j)*sin(ip1(2,j)),ip1(3,j)],([3,3])),.true.)
            !o-j-k plane equation
            p2=normal_triface(reshape([0.d0,0.d0,0.d0,ip1(1,j)*cos(ip1(2,j)),ip1(1,j)*sin(ip1(2,j)),ip1(3,j),&
            ip1(1,k)*cos(ip1(2,k)),ip1(1,k)*sin(ip1(2,k)),ip1(3,k)],([3,3])),.true.)
            !o-i-k plane equation
            p3=normal_triface(reshape([0.d0,0.d0,0.d0,ip1(1,i)*cos(ip1(2,i)),ip1(1,i)*sin(ip1(2,i)),ip1(3,i),&
            ip1(1,k)*cos(ip1(2,k)),ip1(1,k)*sin(ip1(2,k)),ip1(3,k)],([3,3])),.true.)!三角曲面矩阵0,0,0，最小角xyz，最大角xyz对应的法向量
            !p1=normal_triface(reshape([node(v1).coord,xy1(:,2),xy1(:,3)],([3,3])),.true.)
            !p2=normal_triface(reshape([node(v1).coord,xy1(:,3),xy1(:,4)],([3,3])),.true.)
            !p3=normal_triface(reshape([node(v1).coord,xy1(:,2),xy1(:,4)],([3,3])),.true.)最好不用直角坐标系求法向量因为没法判断极坐标角度大小
        !确定积分上下限sita        
        t1=ip1(2,i)!最小值
        t2=ip1(2,j)!中间值
        t3=ip1(2,k)!最大值    
        if(t2-t1<PI) then
            s0=t1;s1=t2
        else
            s0=t2;s1=t1+2*PI
        endif
        if (p1(3)*p3(3)==0) then
        t1=0
        else
        t1=((p1(1)*(sin(s1)-sin(s0))+p1(2)*(cos(s0)-cos(s1)))/-p1(3))-((p3(1)*(sin(s1)-sin(s0))+p3(2)*(cos(s0)-cos(s1)))/-p3(3))
        end if
        if(t3-t2<PI) then
            s0=t2;s1=t3
        else
            s0=t3;s1=t2+2*PI
        endif
        if (p2(3)*p3(3)==0) then
        t2=0
        else
        t2=((p2(1)*(sin(s1)-sin(s0))+p2(2)*(cos(s0)-cos(s1)))/-p2(3))-((p3(1)*(sin(s1)-sin(s0))+p3(2)*(cos(s0)-cos(s1)))/-p3(3))
        end if
        area=wr**2*(abs(t1)+abs(t2))
        
        
    else!开始计算棱在轴线上情况的面积
        do i=3,4
            ip1(:,i-2)=xy1(:,i)*wr/norm2(xy1(1:2,i))
            ip1(2,i-2)=atan2(ip1(2,i-2),ip1(1,i-2))
            ip1(1,i-2)=wr
             !坐标系变换为柱坐标，第一行r第二行角度（要替换）第三行高度（不用替换,但高度发生了比例变化）
        end do
         i=minloc(ip1(2,1:2),dim=1) !找角最小的位置赋值给i
         k=maxloc(ip1(2,1:2),dim=1) !找角最大的位置赋值给k
         !o-i-k plane equation
         !p4=normal_triface(reshape([0.d0,0.d0,0.d0,ip1(1,i)*cos(ip1(2,i)),ip1(1,i)*sin(ip1(2,i)),ip1(3,i),&
         !ip1(1,k)*cos(ip1(2,k)),ip1(1,k)*sin(ip1(2,k)),ip1(3,k)],([3,3])),.true.)
        p4=normal_triface(reshape([xy1(:,1),xy1(:,3),xy1(:,4)],([3,3])),.true.)
         !v2-i-k plane equation
         p5=normal_triface(reshape([xy1(:,2),xy1(:,3),xy1(:,4)],([3,3])),.true.)
        !确定积分上下限sita        
        t1=max(ip1(2,1),ip1(2,2))
        t2=min(ip1(2,1),ip1(2,2))        
        if(t1-t2<PI) then
            s0=t2;s1=t1
        else
            s0=t1;s1=t2+2*PI
        endif
       d1=-dot_product(p4,xy1(:,3))
       d2=-dot_product(p5,xy1(:,3))
        area=wr*abs(((wr*p4(1)*(sin(s1)-sin(s0))+wr*p4(2)*(cos(s0)-cos(s1))+d1*(s1-s0))/-p4(3))-((wr*p5(1)*(sin(s1)-sin(s0))+wr*p5(2)*(cos(s0)-cos(s1))+d2*(s1-s0))/-p5(3)))
         
 endif
!write(66,*)'p1=', p1
!write(66,*)'p2=', p2
!write(66,*)'p3=',p3
!write(66,*)'p4=',p4
!write(66,*)'p5=',p5
!write(66,*)'isp1=',isp1
!write(66,*)'S=',area
!write(66,*)' ' 
    endsubroutine  

SUBROUTINE WELLBORE_Q_K_UPDATE_ANALYTICAL(STEPDIS,IEL,ISTEP,IITER)
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
        A1(2:3)=1.0/(ELEMENT(IEL).PROPERTY(2:3))
        
        
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
    ! SUBROUTINE wellbore_area(IELT,V1,V2,WR,AREA)
    ! !求以V1为其中一个节点的四面体单元IELT的面与半径为WR,轴线为V1-V2的圆柱面的交线所围成的柱面面积(v1-v2范围内的面积)。
    ! !假定：1)单元ielt的其中一个(v1)或2个(v1,v2)节点落在圆柱轴线V1-V2上.
    ! !算法：1）先求出各边与圆柱面的交点，然后数值积分，数值积分区域的边界由三角面与圆柱面的交线确定
    !     USE solverds,ONLY:NODE,ELEMENT
    !     USE GeoMetricAlgorithm
    !     IMPLICIT NONE
    !     INTEGER,INTENT(IN)::IELT,V1,V2
    !     REAL(8),INTENT(IN)::WR
    !     REAL(8),INTENT(OUT)::AREA
    !     REAL(8)::AXIS(3),IP1(3,4),g2l1(3,3),org1(3),T1,XY1(3,4),P1(3)
    !     INTEGER::I,J,K,NODE1(4),N1,V11,V21
    !     LOGICAL::isint,ISP1
        
        
        
    !     NODE1=PACK([1:4],ELEMENT(IELT).NODE(1:4)==V1,[0,0,0,0])
    !     NODE1(2:4)=PACK([1:4],ELEMENT(IELT).NODE(1:4)/=V1)
        
    !     ISP1=.FALSE.
    !     IF(ANY(ELEMENT(IELT).NODE(1:4)==V2)) THEN
    !         ISP1=.TRUE.
    !         V21=MINLOC(ABS(ELEMENT(IELT).NODE(1:4)-V2),DIM=1)
    !     ENDIF
    !     DO I=1,4
    !         XY1(:,I)=NODE(ELEMENT(IELT).NODE(NODE1(I))).COORD-NODE(V1).COORD
    !     ENDDO
    !     AXIS=node(V2).coord-node(v1).coord
        
    !     !SET UP LOCAL SYSTEM
    !     ORG1=NODE(V1).COORD               
    !     CALL GEN_CORDINATE_SYSTEM(ORG1,G2L1,ZV=AXIS)
    !     AXIS=MATMUL(G2L1,AXIS)
    !     N1=0
    !     DO I=2,4
    !         XY1(:,I)=MATMUL(G2L1,XY1(:,I))
    !         if(ISP1.AND.I==V21) CYCLE
    !         !交点
    !         T1=norm2(xy1(1:2,i))
    !         IF(T1<WR) THEN
    !             PRINT *, 'NODE(I) IS INSIDE THE WELLBORE,I=',ELEMENT(IELT).NODE(NODE1(I))
    !             !ERROR STOP
    !         ENDIF
    !         N1=N1+1
    !         IP1(:,N1)=XY1(:,I)*WR/T1
    !         !to cylinder system
    !         IP1(2,N1)=atan2(ip1(2,N1),ip1(1,N1))
    !         IP1(1,N1)=WR;
        
    !         !当四面体的一边平行axis时，另一交点可根据相似三角形得到。
    !         !IP13//IP24
    !         IF(ISP1) THEN
    !             IP1(1:2,2+N1)=IP1(1:2,N1)
    !             IP1(3,2+N1)=IP1(3,N1)+AXIS(3)*(T1-WR)/T1
    !         ENDIF
            
            
    !     END DO
        
    !     if(N1==3) THEN
    !         I=MINLOC(IP1(2,1:3),DIM=1)
    !         K=MAXLOC(IP1(2,1:3),DIM=1)
    !         DO I1=1,3
    !             IF(I1/=I.AND.I1/=K) THEN
    !                 J=I1
    !                 EXIT
    !             ENDIF
    !         ENDDO
    !         !O-I-K PLANE EQUATION
    !         P1=NORMAL_TRIFACE(RESHAPE([0.D0,0.D0,0.D0,IP1(1,I)*COS(IP1(2,I)),IP1(1,I)*SIN(IP1(2,I)),IP1(3,I),&
    !         IP1(1,K)*COS(IP1(2,K)),IP1(1,K)*SIN(IP1(2,K)),IP1(3,K)],([3,3])),.TRUE.)
            
    !         CALL AREA_ON_CYLINDER_SURFACE_CUT_BY_TWO_TRIANGLES(IP1(:,I),IP1(:,J),IP1(:,I),IP1(:,K),AREA)
    !     ELSE

    !     ENDIF
    
             
       
    
    
    ! ENDSUBROUTINE


