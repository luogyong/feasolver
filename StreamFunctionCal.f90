MODULE STREAMFUNCTION
!假定：
!1)渗透系数的主轴与坐标轴平行
!2)承压
!3)单联通域，无奇异点
    
    USE solverds, ONLY:ELEMENT,ENUM,NODE,NNUM,ECP,SPG2D,CAX_SPG,NDOF,TDISP,ADOF,KM,SOLVER_CONTROL,&
        ROWINDEX,JCOL,diaglkmloc,bc_disp,BD_NUM,UM,isIniSEdge
    USE MESHGEO,ONLY:BLOOPS,NBLOOPS
    IMPLICIT NONE
    
    PRIVATE
    PUBLIC::STREAMFUNCTIONCAL
    
    integer, parameter :: RPRE = 8
    
    TYPE SF_KQ
        LOGICAL::ISINI=.FALSE.
        REAL(KIND=RPRE),ALLOCATABLE::K(:,:),QC(:,:)
    CONTAINS
        PROCEDURE,PASS::INITIALIZE=>KQ_INI
    END TYPE
    TYPE(SF_KQ),ALLOCATABLE::SFKQ(:)
    REAL(KIND=RPRE),ALLOCATABLE::Q1(:)
    
    CONTAINS
    

    
    FUNCTION STREAMFUNCTIONCAL() RESULT(VSF)

        INTEGER::I,J,K
        REAL(KIND=RPRE)::VSF(NNUM)
       
        
        IF(.NOT.ALLOCATED(SFKQ)) ALLOCATE(SFKQ(ENUM))
        IF(.NOT.ALLOCATED(Q1)) ALLOCATE(Q1(NDOF))
        
        Q1=0.D0
        DO I=1,ENUM
            if(element(I).isactive==0) CYCLE 
            CALL SFKQ(I).INITIALIZE(I)
            if(.NOT.SFKQ(I).ISINI) CYCLE 
            Q1(ELEMENT(I).G)=Q1(ELEMENT(I).G)+MATMUL(SFKQ(I).QC,TDISP(ELEMENT(I).G))
        ENDDO
        
        if(.not.isIniSEdge) CALL Model_MESHTOPO_INI()
        IF(NBLOOPS>1) then
            VSF=0.D0
            PRINT *, "simple connected domain is assumed."
            RETURN
        endif    
        
        CALL assemble_km_SF()
        CALL BC_SF()
        CALL SOLVE_SF()
        
        WHERE(NODE.DOF(4)>0)  VSF=Q1(NODE.DOF(4))
        
    END FUNCTION
    
    
    
    
    SUBROUTINE KQ_INI(SELF,IELT)
        CLASS(SF_KQ)::SELF
        INTEGER,INTENT(IN)::IELT
        REAL(KIND=RPRE)::KM1(2,2),r1,DET1,t1
        REAL(KIND=RPRE),ALLOCATABLE::B1(:,:)
        INTEGER::I,J,K,NDOF1
    

        if(SELF.ISINI.OR.element(IELT).isactive==0) RETURN 
        if(element(iELT).ec==spg2d.or.element(iELT).ec==cax_spg)  then
            NDOF1=ELEMENT(IELT).NDOF
            allocate(SELF.K(NDOF1,NDOF1),SELF.QC(NDOF1,NDOF1))
            allocate(B1(2,NDOF1))

		    r1=1.0
		    if(element(ielt).ec==CAX_SPG) then
			    r1=abs(element(ielt).xygp(1,i))
			    if(abs(r1)<1e-7) r1=1.0e-7
		    end if            
            
            
        !    call BTDB(element(IELT).b,element(IELT).nd,element(IELT).ndof,element(IELT).ngp, &
				    !DM1,element(IELT).nd,element(IELT).nd, &
				    !SELF.K,element(IELT).ndof,element(IELT).ndof, &
				    !ecp(element(IELT).et).weight,ecp(element(IELT).et).ngp, &
				    !element(IELT).detjac,element(IELT).ngp,ielt)
            
            self.k=0.0D0
            self.qc=0.0d0
	        do i=1,element(IELT).ngp
		        self.k=self.k+ecp(element(IELT).et).weight(I)*element(IELT).DETJAC(I)*R1* MATMUL( &
				        TRANSPOSE(element(IELT).b(:,:,I)),element(IELT).b(:,:,I))
                
                KM1=ELEMENT(ielt).D*element(ielt).kr(i)
                !DET1=KM1(1,1)*KM1(2,2)-KM1(1,2)*KM1(2,1)
                !KM1=KM1/DET1
                T1=km1(1,1);km1(1,1)=km1(2,2);km1(2,2)=t1;
                B1(1,:)=element(IELT).b(2,:,I)
                B1(2,:)=-element(IELT).b(1,:,I)                
		        self.QC=self.QC+ecp(element(IELT).et).weight(I)*element(IELT).DETJAC(I)*R1* MATMUL( &
				        MATMUL(TRANSPOSE(element(IELT).b(:,:,I)),KM1),B1)	                
            ENDDO
            
        endif
        
        SELF.ISINI=.TRUE.
    
    ENDSUBROUTINE
    
    subroutine assemble_km_SF()
	   
	        
	    integer::i,j,k,nj,n1=0
	    integer::loc,rdof,cdof,loc1
	    real(8)::var1,t1,t2
	
	    km=0.0D0
	    n1=0
	    do i=1,enum
            if(element(i).isactive==0) cycle 
		    !renew elements stiffness
		    t1=1.d0
			
		    !assemble the total matrix
		    do j=1,element(i).ndof !row
			
			    rdof=element(i).g(j)

			    do k=1,element(i).ndof  !colume
				    cdof=element(i).g(k)
				    if(solver_control.issym) then
					    if(solver_control.ismkl) then
						    !mkl solver, stored upper part
						    if(rdof>cdof) cycle						
					    else
						    !default solver, stored lower part
						    if(rdof<cdof) cycle  !if colume>row, skip for it is in upper part
					    end if
				    end if
				    if(solver_control.ismkl) then
				
                        do loc=rowindex(rdof),rowindex(rdof+1)-1
                            if(jcol(loc)==cdof) exit
                        enddo
                    else
                        PRINT *, 'KML IS REQUIRED.'
                        !loc=bw(rdof)-(Lmre(rdof)-cdof)
                    endif

				
				    km(loc)=km(loc)+t1*SFKQ(I).K(j,k)
                    if(isnan(km(loc))) then
                        print *, 'NAN VALUE IN ELEMENT(I).KM(J,K),I,J,K',I,J,K
                        STOP
                    endif
				
			    end do
		    end do
	    end do
	
	    !为避免对角线元素为零，保持正定，令没激活的自由度的对角线元素等于1，不会影响计算结果。
	    do i=1,ndof
		    if(adof(i)==0) then
			    km(diaglkmloc(i))=1.d0
		    endif	
	    enddo
	

    end subroutine    
    
    SUBROUTINE BC_SF()
        integer:: j,dof1,NA1(NNUM),N1,n2,N3,outbc1
        real,allocatable::ra1(:)
        
      !  na1=0
      !  do j=1,bd_num
		    !if(bc_disp(j).isdead==1) cycle
      !      if(bc_disp(j).dof==4) then
      !          na1(bc_disp(j).node)=1
      !      endif
      !  enddo
      !  
      !  !假定为单连通域
      !  !假定每个水头边界段的至少包含2个节点，否则为奇异点
      !  !找到水头边界段的第一个节点
      !  do j=2,bloops(1).nnode !首尾节点相同
      !      if(na1(bloops(1).node(j))==1.and.na1(bloops(1).node(j-1))==0) THEN
      !          N1=J
      !          exit
      !      ENDIF
      !  enddo
      !  
      !  allocate(ra1(bloops(1).nnode))
      !  ra1=0.d0
      !  
      !  ra1(n1)=node(bloops(1).node(n1)).q
      !  do j=1,bloops(1).nnode-2
      !      n2=mod(n1+j-1,bloops(1).nnode-1)+1
      !      IF(N2==1) THEN 
      !          N3=bloops(1).nnode-1
      !      ELSE
      !          N3=N2-1
      !      ENDIF
      !      ra1(n2)=ra1(N3)+node(bloops(1).node(n2)).q            
      !  enddo
      !  !修正各水头段的第一个节点的流量        
      !  do j=2,bloops(1).nnode
      !      if(na1(bloops(1).node(j))==1.and.na1(bloops(1).node(j-1))==0) THEN
      !          ra1(J)=ra1(J-1)
      !      endif
      !  enddo
      !  !修正各水头段的首尾节点为定流函数边界        
      !  do j=1,bloops(1).nnode-1
      !      if(j<2) then
      !          n1=bloops(1).nnode-1
      !      else
      !          n1=j-1
      !      endif
      !      if((na1(bloops(1).node(j))==1.and.na1(bloops(1).node(n1))==0).or. &
      !          (na1(bloops(1).node(j))==1.and.na1(bloops(1).node(j+1))==0)) THEN
      !          na1(bloops(1).node(j))=0
      !      endif
      !  enddo        
        
        outbc1=maxloc(bloops.isoutbc,dim=1)
        do j=1,bloops(outbc1).nnode-1
		    !if(bc_disp(j).isdead==1) cycle
            !if(na1(bloops(1).node(j))/=0) cycle
		    dof1=node(bloops(outbc1).node(j)).dof(4)
            if(adof(dof1)==0) cycle
            
		    km(diaglkmloc(dof1))=UM
            !令原第一个水头边界节点的流函数的值为0
            Q1(dof1)=0.d0
      !      if(solver_control.solver==N_R.OR.solver_control.solver==INISTIFF) then
			   ! load1(dof1)=(0.0d0-stepdis(dof1))*UM	!		
		    !else
			   ! load1(dof1)=0.0d0*UM			
		    !end if
            
            exit !一个即可,发现只施加一个节点的边界要比施加所有定流函数边界所得流网的质量要好。
	    end do
    ENDSUBROUTINE
    
    
    SUBROUTINE SOLVE_SF()
        USE MKLDS
        INTEGER::error
        REAL(KIND=RPRE)::SOL1(NDOF)
    	! Factor the matrix.！MKL_DSS_POSITIVE_DEFINITE
        error = DSS_FACTOR_REAL( handle, MKL_DSS_POSITIVE_DEFINITE,km)
        IF (error /= MKL_DSS_SUCCESS) THEN
            !PRINT *, 'MKL_DSS_POSITIVE_DEFINITE FAILED.TRY MKL_DSS_INDEFINITE OPTION.'
			error = DSS_FACTOR_REAL( handle, MKL_DSS_INDEFINITE,km)
        ENDIF
                        
		IF (error /= MKL_DSS_SUCCESS) call mkl_solver_error(error)	
    
        
		error = DSS_SOLVE_REAL(handle, MKL_DSS_DEFAULTS, Q1, 1, SOL1 )
		IF (error /= MKL_DSS_SUCCESS) call mkl_solver_error(error)
		Q1=SOL1        
    
    
    ENDSUBROUTINE
    
    
    
END MODULE
    