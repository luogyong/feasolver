
subroutine solve()
	use solverds
	implicit none
	
!	select case(solver_control.type)
!		case(SPG)
!		!	call solve_SPG()
!		case default
!			call solve_SLD()
!	end select
	
	call solve_SLD()


end subroutine



subroutine solver_initialization(kref,iincs,iiter)
	use solverds
	implicit none
	integer::kref !kref=1, need to reformulat the stiffness matrix.
			     !kref=2, no need to reformulate the stiffness matrix
	integer::iincs,iiter

	
	if(iiter==1) then
		kref=1
		return
	end if
	
	kref=2
	
	select case(solver_control.solver)
		case(N_R)
			if(solver_control.bfgm==consistent.or.solver_control.bfgm==continuum) KREF=1
			IF(solver_control.bfgm==iniflux) then
				if(isref_spg.and.(.not.solver_control.isfu)) KREF=1
				if(iiter==1) kref=1
			end if
			
			if(iiter==1) kref=1
			
		case(DIRECTI)
			kref=1
		case(CINIATAN)
			IF(iiter==1) kref=1		
		case(CINIATAN2)
			IF(iiter==2)  kref=1
		!case(isexca)
	
	end select	
end subroutine

subroutine assemble_km(istep)
	use solverds
	implicit none
	integer,intent(in)::istep
	integer::i,j,k,nj,n1=0
	integer::loc,rdof,cdof,loc1
	real(8)::var1,t1,t2
	
	km=0.0D0
	IF(.NOT.ALLOCATED(DIAGLKM)) ALLOCATE(DIAGLKM(NDOF))
	DIAGLKM=0.0D0
		

	n1=0
	do i=1,enum
        if(element(i).isactive==0) cycle 
		!renew elements stiffness
		t1=1.d0
		select case(element(i).ec)
			case(CND)
				var1=((tdisp(element(i).g(1))+load(element(i).g(1)))+(tdisp(element(i).g(2))+load(element(i).g(2))))/2.0
				element(i).property(1)=material(element(i).mat).GET(1,ISTEP)+material(element(i).mat).GET(2,ISTEP)*var1				
				t1=element(i).property(1)
			case(cps,cpe,cax)
				element(i).property(1)=1.0
				t1=element(i).property(1)
			case(stru)
				t1=element(i).property(1)
		end select

				
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
                    loc=bw(rdof)-(Lmre(rdof)-cdof)
                endif
                !if(loc1/=loc) then
                !    stop "error in assemble_km"
                !endif
                
                !abs(jcol(:1)-cdof)
				!loc=rowindex(rdof)+
				!if(solver_control.ismkl) loc=adrn(loc) !非对称默认的求解为mkl solver.
				
				km(loc)=km(loc)+t1*element(i).km(j,k)
                if(isnan(km(loc))) then
                    print *, 'NAN VALUE IN ELEMENT(I).KM(J,K),I,J,K',I,J,K
                    STOP
                endif
                if(.not.stepinfo(istep).issteady) then !目前仅考虑为渗流问题，每个节点只有一个自由度。
					if(element(i).ec==spg.or.element(i).ec==spg2d.or.element(i).ec==cax_spg) then
						
						km(loc)=km(loc)+element(i).cmm(j,k)
					else
						pause 'TO BE IMPROVED.SUB=assemble_km(istep)'
					end if
					
                end if
                
				IF(RDOF==CDOF) DIAGLKM(RDOF)=KM(LOC)
				
			end do
		end do
	end do
	
	!为避免对角线元素为零，保持正定，令没激活的自由度的对角线元素等于1，应该不会影响计算结果。
	do i=1,ndof
		if(adof(i)==0) then
			km(rowindex(i))=1.d0
		endif	
	enddo
	

end subroutine

subroutine dof_activate(istep)
	use solverds
	implicit none
	integer,intent(in)::istep
    integer::i
	
	adof=0
	do i=1,enum
		if(element(i).isactive==0) cycle 
		adof(element(i).g)=1
	enddo

endsubroutine

!转换成坐标格式及mkl格式进行总刚的存储
!subroutine irowjcol()
!	use solverds
!	implicit none
!	INTEGER*8::i,j,k,nj,n1=0,NNZ1
!	INTEGER*8::loc,rdof,cdof !RDOF FOR ROW DOF INDEX; CDOF FOR COLUME DOF INDEX
!    integer::err
!	real(8)::var1,t1
!	integer,allocatable::irow1(:),jcol1(:)
!	
!    NNZ1=BW(NDOF)
!	allocate(irow(nnz1),STAT=ERR)
!    allocate(jcol(nnz1),STAT=ERR)
!    allocate(adrn(nnz1),STAT=ERR)
!    !allocate(adrn(nnz1),jcol(nnz1),adrn(nnz1))
!	irow=0
!	jcol=0
!	adrn=0
!	
!	do i=1,enum
!		!assemble the total matrix
!		do j=1,element(i).ndof !row
!			rdof=element(i).g(j)
!			do k=1,element(i).ndof  !colume
!				cdof=element(i).g(k)
!				if(solver_control.issym) then
!					if(solver_control.ismkl) then						
!						if(rdof>cdof) cycle  !mkl solver, 存上三角（包含对角线元素）,下三角部分跳过
!					else
!						if(rdof<cdof) cycle  !default solver,存下三角（包含对角线元素）,上三角部分跳过				
!					end if
!				end if
!				
!				loc=bw(rdof)-(Lmre(rdof)-cdof)	!BW:THE MOST RIGHT ENTRY INDEX IN TOTAL MATRIX
!															!Lmre(rdof):the column of the most right entry in the row of rdof
!				irow(loc)=rdof
!				jcol(loc)=cdof	
!			end do
!		end do
!	end do
!	
!	n1=0
!	!adrn STORED THE NONZERO ENTRIES(IROW(I)/=0) 
!	do i=1,NNZ1
!		if(irow(i)/=0) then
!			n1=n1+1
!			irow(n1)=irow(i) !n1<i
!			jcol(n1)=jcol(i)
!			adrn(i)=n1  ! 
!			if(irow(i)==jcol(i)) diaglkmloc(irow(i))=n1 !FOR SETTING THE BOUNDARY CONDITIONS 
!			
!		end if
!	end do
!    
!	nnz=n1
!	IF(NNZ<0) THEN
!        PRINT *, "THE NUMBER NONZERO ENTRIES IN KM IS TOO LARGE TO BE HANDEL BY INTEGER*4."
!        STOP
!    ENDIF
!	!rowindex for mkl format
!	allocate(irow1(nnz),jcol1(nnz),rowindex(ndof+1))
!	irow1=irow(1:nnz)
!	jcol1=jcol(1:nnz)
!	deallocate(irow,jcol)
!	allocate(irow(nnz),jcol(nnz))
!	irow=irow1
!	jcol=jcol1
!	deallocate(irow1,jcol1)
!	if(.not.allocated(km)) allocate(km(nnz))
!	km=0.D0
!	
!	rowindex(1)=1
!	j=1
!	i=1
!	do while(j<=nnz)
!		if(irow(j)/=i) then
!			i=i+1
!			rowindex(i)=j			
!		end if
!		j=j+1
!	end do
!	rowindex(i+1)=j !最后一元素，i=ndof
!	
!end subroutine

subroutine irowjcol_SEC()
	use solverds
	implicit none
	INTEGER*8::i,n1=0	
    integer::err

    !diaglkmloc(irow(i))=n1 !FOR SETTING THE BOUNDARY CONDITIONS 
			
	    
	nnz=SUM(DOFADJL.NDOF)
	IF(NNZ<0) THEN
        PRINT *, "THE NUMBER NONZERO ENTRIES IN KM IS TOO LARGE TO BE HANDEL BY INTEGER*4."
        STOP
    ENDIF
    !rowindex for mkl format


	allocate(jcol(nnz),rowindex(ndof+1),STAT=ERR)
    n1=0;    
    DO I=1,NDOF
        rowindex(I)=N1+1
        JCOL(N1+1:N1+DOFADJL(I).NDOF)=DOFADJL(I).DOF
        N1=N1+DOFADJL(I).NDOF        
    ENDDO
	rowindex(NDOF+1)=rowindex(NDOF)+DOFADJL(NDOF).NDOF
	
end subroutine




subroutine incremental_load(iincs,iiter,method,isubts)
	use solverds
	implicit none
	integer,intent(in)::iincs,iiter,method,isubts
	integer::i,j
	integer::dof1
	real(kind=dpn)::t1=0,t2=0
    real(kind=dpn),external::tsfactor
	
	
	if(iiter==1) then
		load=0.0D0
		if(iincs>0)	then		
			do i=1,bl_num
                dof1=node(bc_load(i).node).dof(bc_load(i).dof)
                if(adof(dof1)==0) cycle
				if(sf(bc_load(i).sf).factor(iincs)==-999.D0) cycle
				
				t2=0
				if(bc_load(i).dof/=4) t2=sf(bc_load(i).sf).factor(max(iincs-1,0))
				if(t2==-999.d0 .OR. bc_load(i).ISINCREMENT==1) t2=0
				t1=(sf(bc_load(i).sf).factor(iincs)-t2)*tsfactor(iincs,isubts,stepinfo(iincs).loadtype)
				load(dof1)=load(dof1)+bc_load(i).value*t1
			end do		
		else !generate initial load for initial stress field
			call bf_initialstress()
		endif
		
		stepload=load
		Tload=Tload+load !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		TLOAD(DOFHEAD)=LOAD(DOFHEAD) !对于seepage，由于荷载是流量，子步之间流量这里认为不具有可加性（注意是单位时间，相当于流速）。
		
		!if(solver_control.type==spg) then
		!	Tload=load 
		!else
		!	Tload=Tload+load
		!end if
        
        if(.not.stepinfo(iincs).issteady) then
            do i=1,enum
				if(element(i).isactive==0) cycle
                if(any(element(i).cmm/=0.0d0)) then
                    load(element(i).g)=load(element(i).g)+matmul(element(i).cmm,inivaluedof(element(i).g))
                end if
            end do
        end if
		
!		write(10,*) load(2:ndof:2)
	else
		select case(method)
			case(DIRECTI,LELASTIC)
				load=stepload
			case DEFAULT
				load=0.0D0 !the load for N_r in every interation is the incrementals. 
		end select
	end if
	

	
end subroutine

!对于全量边界时间处理。
function tsfactor_total(istepfun,istep,isubts)
	use solverds
	implicit none
    integer,intent(in)::istep,isubts,istepfun
    real(kind=DPN)::tsfactor_total,tsfactor,tsfactor1
    external::tsfactor
    
	
    tsfactor1=tsfactor(istep,isubts,stepinfo(istep).bctype) !to be improved.
	
	if(istep>1) then
		!Q1=Q_i-1+(Qi-Q_i-1)*tsfactor
		tsfactor_total=sf(istepfun).factor(istep-1)+ &
							(sf(istepfun).factor(istep)-sf(istepfun).factor(istep-1))*tsfactor1
	else
		tsfactor_total=sf(istepfun).factor(istep)*tsfactor1
	end if
end function

function tsfactor(istep,isubts,itype)
	use solverds
	implicit none
    real(kind=DPN)::tsfactor
    integer,intent(in)::istep,isubts,itype
	
	if((stepinfo(istep).issteady).or.(itype==step)) then
		tsfactor=1.0D0
	else
		tsfactor=sum(timestep(istep).subts(1:isubts))/sum(timestep(istep).subts(1:timestep(istep).nsubts))
		if(itype==Reramp) tsfactor=1.0d0-tsfactor
	end if
	
end function

subroutine bc(iincs,iiter,load1,stepdis,isubts)
	use solverds
	use ds_hyjump
    use DS_SlopeStability, only:slopeparameter
    USE quicksort
	implicit none
	integer,INTENT(in)::iincs,iiter,isubts
	real(kind=dpn),intent(in)::stepdis(ndof)
	REAL(kind=dpn),INTENT(IN OUT)::LOAD1(NDOF)
	integer::i,j,dof1,N1
	real(kind=dpn)::t1=0,t2=0
	real(kind=dpn),external::tsfactor
    REAL(KIND=DPN),ALLOCATABLE::AR1(:)
    INTEGER,ALLOCATABLE::IAR1(:)
    type(bc_tydef),allocatable::BC1(:)
	
	if(iiter==1) then
		do i=1,nhjump
			call HJ_WaterSurfaceProfile_RC(i,iincs,isubts)
			do j=1,hjump(i).nnode
				if(HJump_CP==2) then !水跃区边界按不考虑水跃作用，急流缓流的水面线中两者的大值进行计算
					bc_disp(hjump(i).bc_node(j)).value=hjump(i).xy(2,j)+max(hjump(i).xy(3,j),hjump(i).xy(4,j))
				else !水跃区边界考虑水跃作用
					bc_disp(hjump(i).bc_node(j)).value=hjump(i).xy(9,j)
				endif
			end do
		end do
	end if
    
    !SORT BY X. just for slope bc pa
    IF(solver_control.isslopepa==4.AND.iiter==1.and.iincs==1) THEN
        IF(ALLOCATED(AR1)) DEALLOCATE(AR1,IAR1,BC1)
        N1=COUNT(BC_DISP.DOF==4,DIM=1)
        ALLOCATE(AR1(N1),IAR1(N1))
        ALLOCATE(BC1(N1))
        N1=0
        DO I=1,BD_NUM
            IF(BC_DISP(I).DOF/=4) CYCLE
            N1=N1+1
            AR1(N1)=NODE(BC_DISP(I).NODE).COORD(1)
            IAR1(N1)=I
        ENDDO
        !SORT BY X
        CALL QUICK_SORT(AR1(1:N1),IAR1(1:N1))
        DO I=1,N1
            BC1(I)=BC_DISP(IAR1(I))
        ENDDO
        !PUT DOF=4 TO THE TAIL
        N1=0 
        DO I=1,BD_NUM
            IF(BC_DISP(I).DOF/=4) THEN
                N1=N1+1
                BC_DISP(N1)=BC_DISP(I)            
            ENDIF  
        ENDDO
        BC_DISP(N1+1:)=BC1
        
        DEALLOCATE(BC1,AR1,IAR1)
    ENDIF
    
	
	do i=1,bd_num
        dof1=node(bc_disp(i).node).dof(bc_disp(i).dof)	
        if(adof(dof1)==0)  then
            bc_disp(i).isdead=1
            cycle
        endif
        
		t2=0
		if(bc_disp(i).dof/=4) t2=sf(bc_disp(i).sf).factor(max(iincs-1,0))
		if(t2==-999.d0 .OR. bc_disp(i).ISINCREMENT==1 ) t2=0.d0
		t1=(sf(bc_disp(i).sf).factor(iincs)-t2)*tsfactor(iincs,isubts,stepinfo(iincs).bctype)
		
        
		if(iiter==1) then			
			
			if(sf(bc_disp(i).sf).factor(iincs)==-999.D0) then
				bc_disp(i).isdead=1
				cycle			
			end if		

			!对于水头边界，如果边界值小于其相应的位置水头，则认为该边界无效，不起作用。			
			if(bc_disp(i).dof==4) then
				if(bc_disp(i).value*t1<node(bc_disp(i).node).coord(ndimension)) then
					bc_disp(i).isdead=1
					cycle
				end if
			end if
            
            !for slope bc PA only            
            if(slopeparameter.NBCPA*slopeparameter.IBCPA>0 &
            .AND.slopeparameter.NBCPA>=slopeparameter.IBCPA &            
            .AND.bc_disp(i).dof==4) then
                !TOE
                
                IF((NODE(BC_DISP(I).NODE).COORD(1)>=slopeparameter.BCEXIT(1,SLOPEPARAMETER.IBCPA) &
                   .AND. NODE(BC_DISP(I).NODE).COORD(1)<=slopeparameter.BCEXIT(2,SLOPEPARAMETER.IBCPA)) &
                   .OR.(NODE(BC_DISP(I).NODE).COORD(1)-slopeparameter.BCEXIT(1,SLOPEPARAMETER.IBCPA))* &  !前面已经排序，位于两节点之间，左节点也激活
                       (NODE(BC_DISP(min(I+1,bd_num)).NODE).COORD(1)-slopeparameter.BCEXIT(1,SLOPEPARAMETER.IBCPA))<0 &
                   .OR.(NODE(BC_DISP(I).NODE).COORD(1)-slopeparameter.BCEXIT(2,SLOPEPARAMETER.IBCPA))* &   !位于两节点之间，右节点也激活
                       (NODE(BC_DISP(MAX(I-1,1)).NODE).COORD(1)-slopeparameter.BCEXIT(2,SLOPEPARAMETER.IBCPA))<0 ) THEN
                   
                   IF(abs(slopeparameter.BCEXIT(3,SLOPEPARAMETER.IBCPA)+999.d0)>1e-7) then
                         IF(abs(slopeparameter.BCEXIT(3,SLOPEPARAMETER.IBCPA)+9999.d0)<1e-7) then
                            BC_DISP(I).VALUE=slopeparameter.BCEXIT(4,SLOPEPARAMETER.IBCPA)*NODE(BC_DISP(I).NODE).COORD(1) &
                            +slopeparameter.BCEXIT(5,SLOPEPARAMETER.IBCPA)*NODE(BC_DISP(I).NODE).COORD(2) &
                            +slopeparameter.BCEXIT(6,SLOPEPARAMETER.IBCPA)                            
                         ELSE
                            bc_disp(i).value=slopeparameter.BCEXIT(3,SLOPEPARAMETER.IBCPA)                             
                        ENDIF
                   endif
                !CREST   
                ELSEIF((NODE(BC_DISP(I).NODE).COORD(1)>=slopeparameter.BCENTRY(1,SLOPEPARAMETER.IBCPA) &
                   .AND. NODE(BC_DISP(I).NODE).COORD(1)<=slopeparameter.BCENTRY(2,SLOPEPARAMETER.IBCPA)) &
                   .OR.(NODE(BC_DISP(I).NODE).COORD(1)-slopeparameter.BCENTRY(1,SLOPEPARAMETER.IBCPA))* &  !位于两节点之间，左节点也激活
                       (NODE(BC_DISP(min(I+1,bd_num)).NODE).COORD(1)-slopeparameter.BCENTRY(1,SLOPEPARAMETER.IBCPA))<0 &
                   .OR.(NODE(BC_DISP(I).NODE).COORD(1)-slopeparameter.BCENTRY(2,SLOPEPARAMETER.IBCPA))* & !位于两节点之间，右节点也激活
                       (NODE(BC_DISP(MAX(I-1,1)).NODE).COORD(1)-slopeparameter.BCENTRY(2,SLOPEPARAMETER.IBCPA))<0 ) THEN
                   
                   IF(abs(slopeparameter.BCENTRY(3,SLOPEPARAMETER.IBCPA)+999.d0)>1e-7) then
                         IF(abs(slopeparameter.BCENTRY(3,SLOPEPARAMETER.IBCPA)+9999.d0)<1e-7) then
                            BC_DISP(I).VALUE=slopeparameter.BCENTRY(4,SLOPEPARAMETER.IBCPA)*NODE(BC_DISP(I).NODE).COORD(1) &
                            +slopeparameter.BCENTRY(5,SLOPEPARAMETER.IBCPA)*NODE(BC_DISP(I).NODE).COORD(2) &
                            +slopeparameter.BCENTRY(6,SLOPEPARAMETER.IBCPA)                            
                         ELSE
                            bc_disp(i).value=slopeparameter.BCENTRY(3,SLOPEPARAMETER.IBCPA)                             
                        ENDIF
                   endif
                   
                ELSE
					bc_disp(i).isdead=1
                    cycle                 
                ENDIF
            
            endif
			
			bc_disp(i).isdead=0
        else
            if(bc_disp(i).iswellhead) then
            !对于减压自流井，如其流量为正，则令其水头边界失效；反之，如果该点水头大于井口高程，则恢复激活
                
                if(bc_disp(i).isdead==0) then
                    if(NI_NodalForce(node(bc_disp(i).node).dof(4))>1.d-7) bc_disp(i).isdead=1
                else
                    if(stepdis(node(bc_disp(i).node).dof(4))>node(bc_disp(i).node).coord(ndimension)) then
                        bc_disp(i).isdead=0
                        bc_disp(i).value=node(bc_disp(i).node).coord(ndimension)
                    endif
                endif
            endif
            
            if(bc_disp(i).isdead==1) cycle
		end if
		
        
			
		!when Newton Raphson is used, the none zero displacement is applied at the first
		!iteration,a zero displacement is applied at the left iterations.
		if(solver_control.solver==N_R.OR.solver_control.solver==INISTIFF) then
			load1(dof1)=(bc_disp(i).value*t1-stepdis(dof1))*UM			
		else
			load1(dof1)=bc_disp(i).value*t1*UM			
		end if
		
		!km(bw(dof1))=UM
	end do

	do i=1,numNseep
        dof1=node(Nseep(i).node).dof(Nseep(i).dof)	
         if(adof(dof1)==0)  then
            Nseep(i).isdead=1
            cycle
        endif
        
		if(sf(Nseep(i).sf).factor(iincs)==-999.D0) cycle
		
		if(Nseep(i).isdual>0) then
            if(bc_disp(Nseep(i).isdual).isdead==0) cycle            
        end if
		
		if(Nseep(i).isdead==0) then
			dof1=node(Nseep(i).node).dof(Nseep(i).dof)														  
			if(solver_control.solver==N_R) then
				load1(dof1)=(Nseep(i).value-stepdis(dof1))*UM
			else
				load1(dof1)=Nseep(i).value*UM
			end if
		!load1(dof1)=(Nseep(j).value-stepdis(dof1))*sf(Nseep(j).sf).factor(iincs)*UM
		end if
	end do

!###lacy method	
	!IF(IITER>1) THEN
	!	DO I=1,NNUM
	!		dof1=node(I).dof(4)
	!		t1=minNPH+node(I).coord(ndimension)
	!		IF(STEPDIS(DOF1)<t1) then
	!			if(solver_control.solver==N_R) then
	!				load1(dof1)=(T1-stepdis(dof1))*UM
	!			else
	!				load1(dof1)=(T1)*UM
	!			end if
	!		end if		
	!	END DO
	!!	do i=1,numNseep
	!!		dof1=node(Nseep(i).node).dof(Nseep(i).dof)
	!!		t1=STEPDIS(DOF1)-node(Nseep(i).node).coord(ndimension)
	!!		if(t1>0) then
	!!																	  
	!!			if(solver_control.solver==N_R) then
	!!				load1(dof1)=(Nseep(i).value-stepdis(dof1))*sf(Nseep(i).sf).factor(iincs)*UM
	!!			else
	!!				load1(dof1)=Nseep(i).value*sf(Nseep(i).sf).factor(iincs)*UM
	!!			end if
	!!		!load1(dof1)=(Nseep(j).value-stepdis(dof1))*sf(Nseep(j).sf).factor(iincs)*UM
	!!		end if
	!!	end do		
	!END IF


	
end subroutine


subroutine checon(iscon,pvalue,load,ndof,tol)
	implicit none
	integer::ndof
	logical::iscon
	real(8)::pvalue,load(ndof),tol
	real(8)::cvalue=0

	integer::i
	
	cvalue=0.0D0

	do i=1, ndof
		cvalue=cvalue+load(i)**2
	end do
	cvalue=cvalue**0.5
	
	if(tol>abs(cvalue-pvalue)) then 
		iscon=.true.
	else
		iscon=.false.
	end if

	pvalue=cvalue		
	
end subroutine

SUBROUTINE checon_sec(iscon,PUBForce,UBForce,ndof,tol,resdis,convratio,niter)
!
! This subroutine sets converged to .FALSE. if relative change in loads
! and oldlds is greater than tol and updates oldlds.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,INTENT(IN)::NDOF,niter
 REAL(iwp),INTENT(IN)::UBForce(ndof),tol
 REAL(iwp),INTENT(IN OUT)::PUBForce(ndof),resdis,convratio
 real(iwp)::t1=0
 LOGICAL,INTENT(OUT)::iscon
 
 iscon=.false.
 t1=maxval(abs(UBForce))
 if(t1>1e-6) then  
	resdis=MAXVAL(ABS(UBForce-PUBForce))/t1
	iscon=(resdis<=tol)
	
	convratio=count(ABS(UBForce-PUBForce)/t1<=tol)/real(size(ubforce))*100
 else
	!if(niter>1) 
	iscon=.TRUE. !when niter=1,ubforce=0,however,
	!if niter>1 and unborce=0 indicates convergence.
	convratio=100.0	
 end if
 !convratio=count(ABS(UBForce-PUBForce)/t1<=tol)/size(ubforce)*100
 PUBForce=UBForce
 
RETURN
END SUBROUTINE

SUBROUTINE checon_thd(iscon,stepload,UBForce,ndof,tol,resdis,sumforce,convratio,ndofhead,dofhead,niter,ISTOCONV)
!
!
 IMPLICIT NONE
 INTEGER::I
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,INTENT(IN)::NDOF,NDOFHEAD,niter,DOFHEAD(NDOFHEAD)
 REAL(iwp),INTENT(IN)::UBForce(ndof),tol
 REAL(iwp),INTENT(IN OUT)::stepload(ndof),resdis,sumforce,convratio
 LOGICAL,INTENT(OUT)::iscon,ISTOCONV
 REAL(IWP)::RESDIS_SPG,SUMFORCE_SPG
 REAL(IWP),SAVE::LASTCONVRATIO=0.0
 
 iscon=.false.
 resdis_SPG=0.D0;SUMFORCE_SPG=0.D0
 
 
 
 if(ndofhead>0) resdis_SPG=(dot_product(UBFORCE(DOFHEAD),UBFORCE(DOFHEAD)))
 
 resdis=dsqrt(dot_product(UBFORCE,UBFORCE)-resdis_SPG)

 
 if(ndofhead>0) SumForce_SPG=(dot_product(STEPLOAD(DOFHEAD),STEPLOAD(DOFHEAD)))
 SumForce=dsqrt(dot_product(STEPLOAD,STEPLOAD)-SumForce_SPG)
 if(ndofhead>0) then
     resdis_SPG=DSQRT(resdis_SPG)
     SumForce_SPG=DSQRT(SumForce_SPG)
 endif
 
 if(abs(SumForce)<1.d-7) then
     if(abs(resdis)<1.d-7) then
         !iscon=.true.
		 convratio=TOL*0.1
         !return
	 
     end if
 ELSE
	convratio=resdis/SumForce
 end if
 if(ndofhead>0) then
     if(abs(SumForce_SPG)<1.d-7) then
         if(abs(resdis_SPG)<1.d-7) then
             !iscon_SPG=.true.
		     convratio=MAX(TOL*0.1,convratio)
             !return
         end if
     ELSE
	    convratio=MAX(convratio,resdis_SPG/SumForce_SPG)
     end if
 endif
 
 IF(NITER<2) THEN    
    ISTOCONV=.TRUE.
 ELSE
    ISTOCONV=(convratio/LASTconvratio<0.9999)    
 ENDIF
 LASTconvratio=convratio
 !convratio=resdis/SumForce
 
 ISCON=(convratio<=TOL)

 
 
RETURN
END SUBROUTINE

subroutine chodec(tm_t1,sbw,nnum_t) !cholesky decomposition
	implicit none
	integer::tmn,nnum_t
    INTEGER*8::sbw(nnum_t)
	real(8)::tm_t1(sbw(nnum_t))
	integer::i,j
	integer::K1,i1,j1
	real(8)::sum,art(25)
	real(8),allocatable::T(:)
	allocate(T(nnum_t))
	t=0

	tmn=sbw(nnum_t)

	do i=2,nnum_t
		i1=i-sbw(i)+sbw(i-1)+1  !第i行第一个非零元素的列号
		do j=i1,i-1 
			if(j==1) then
				j1=1
			else
				j1=j-sbw(j)+sbw(j-1)+1
			end if
			if(i1>j1) j1=i1
				sum=0
			do K1=j1,j-1
				sum=sum+T(K1)*tm_t1(sbw(j)-j+K1)
			end do
				T(j)=tm_t1(sbw(i)-i+j)-sum
				if(abs(tm_t1(sbw(j)))<1e-15) tm_t1(sbw(j))=1e-15
				tm_t1(sbw(i)-i+j)=T(j)/tm_t1(sbw(j))
				!if(isnan(tm_t1(sbw(i)-i+j))) pause
				tm_t1(sbw(i))=tm_t1(sbw(i))-T(j)*tm_t1(sbw(i)-i+j)
		end do     
	end do 
	deallocate(T)
end subroutine

subroutine chosol(tm_t1,sbw,value,nnum_t,maxbdw_t)

	implicit none
	integer::tmn,nnum_t,maxbdw_t
    INTEGER*8::sbw(nnum_t)
	real(8)::tm_t1(sbw(nnum_t)),value(nnum_t)
	integer::i,j
	integer::i1,j1

	tmn=sbw(nnum_t)

	do i=2,nnum_t
		i1=i-sbw(i)+sbw(i-1)+1
		do j=i1,i-1
			value(i)=value(i)-tm_t1(sbw(i)-i+j)*value(j)
		end do
	end do

	do i=1,nnum_t
		value(i)=value(i)/tm_t1(sbw(i))
!		if(isnan(value(i))) then
!			pause
!		end if
	end do
	do i=nnum_t-1,1,-1
		i1=i+Maxbdw_t
		!i1=i+ldlubw(i)-1
		if(i1>nnum_t) i1=nnum_t
		do j=i+1,i1
			j1=j-sbw(j)+sbw(j-1)+1
			if(i>=j1) value(i)=value(i)-tm_t1(sbw(j)-j+i)*value(j)
		end do  
	end do

end subroutine

subroutine LDLFACTORUPDATE(L,Y,ALPHA,BW,NDOF,MAXBW)
	IMPLICIT NONE
	INTEGER,INTENT(IN)::NDOF,MAXBW
	INTEGER*8,INTENT(IN)::BW(NDOF)
	REAL(8),INTENT(INOUT)::ALPHA
	REAL(8),INTENT(INOUT)::L(BW(NDOF)),Y(NDOF)
	
	INTEGER::I,J,IDIAGL,IJL,NI,I1
	REAL(8)::DELTA=0.,BETA=0.,P=0.
	
	DO J=1,NDOF
		IF(ABS(Y(J))<1.E-15) CYCLE
		IDIAGL=BW(J)		
		DELTA=L(IDIAGL) 
		P=Y(J)
		L(IDIAGL)=L(IDIAGL)+ALPHA*P**2
		BETA=ALPHA*P/L(IDIAGL)
		ALPHA=ALPHA*DELTA/L(IDIAGL)
		
		NI=MIN(NDOF,J+MAXBW)
		DO I=J+1,NI
			I1=I-BW(I)+BW(I-1)+1
			IF(J<I1) CYCLE
			IJL=BW(I)-I+J
			Y(I)=Y(I)-P*L(IJL)
			L(IJL)=L(IJL)+BETA*Y(I)
		END DO
	END DO	
	
END SUBROUTINE



SUBROUTINE LDLFACTORUPDATE_SPG(Y1,NC1,ISTEP)
	USE SOLVERDS 
	IMPLICIT NONE
	integer::i,IDOF1,NC1,ISTEP
	real(8)::ALPHA=0
	real(8)::Y1(NDOF)
	
	if(.not.allocated(isbcchange)) then
		allocate(isbcchange(Numnseep))
		isbcchange=0
	end if
	
	NC1=0
	do i=1,Numnseep
		if(sf(NSEEP(i).sf).factor(ISTEP)==-999.D0) cycle
		
		if(Nseep(i).isdual>0) then
			if(bc_disp(Nseep(i).isdual).isdead==0) cycle
		end if
		
		Y1=0.D0
		if(isbcchange(i)==Nseep(i).isdead) cycle
		NC1=NC1+1
		IDOF1=NODE(NSEEP(I).NODE).DOF(NSEEP(I).DOF)
		if(Nseep(i).isdead==0) then
			ALPHA=UM-DIAGLKM(IDOF1)
		else
			ALPHA=DIAGLKM(IDOF1)-UM
		end if
		Y1(IDOF1)=1.0d0
		CALL LDLFACTORUPDATE(KM,Y1,ALPHA,BW,NDOF,BWMAX)
		
		isbcchange(i)=Nseep(i).isdead
	end do
	
END SUBROUTINE

!subroutine LFTXG_IMSL(NFAC,NL,RFAC,IRFAC,JCFAC,IPVT,JPVT)
!   include "link_fnl_static.h"
!   !use dfimsl
!	use solverds
!	implicit none
!	integer::i,j,k,loc,dof1
!	integer::IPARAM(6),NFAC,NL,IRFAC(NFAC),JCFAC(NFAC),IPVT(NDOF),JPVT(NDOF)
!	real(8)::RPARAM(5),RFAC(NFAC)
!
!	
!	IPARAM(1)=0
!	
!
!	CALL DLFTXG(NDOF,NNZ,km,IROW,JCOL,IPARAM,RPARAM,NFAC,NL,&
!			& RFAC,IRFAC,JCFAC,IPVT,JPVT)
!
!end subroutine
!
!SUBROUTINE LFSXG_IMSL(NDOF,NFAC,NL,RFAC,IRFAC,JCFAC,IPVT,JPVT,LOAD,IPATH)
!	!use dfimsl
!	include "link_fnl_static.h" 
!	IMPLICIT NONE
!
!	INTEGER::NDOF,NFAC,NL,IPATH
!	INTEGER::IRFAC(NFAC),JCFAC(NFAC),IPVT(NDOF),JPVT(NDOF)
!	REAL(8)::LOAD(NDOF),RFAC(NFAC)
!	REAL(8)::X(NDOF)
!	
!	CALL DLFSXG(NDOF,NFAC,NL,RFAC,IRFAC,JCFAC,IPVT,JPVT,LOAD,IPATH,X)
!	load=x
!
!END SUBROUTINE

subroutine ko_initialstress()
	use solverds
	implicit none
	integer::i,j,k
	logical::islast=.false.
	REAL(8)::BLOAD(100),r1
	
	do i=1,enum
		if(element(i).isactive==0) cycle
		BLOAD=0.D0
		do j=1,element(i).ngp
			islast=.false.
			do k=1,geostatic.nsoil
				if(element(i).xygp(2,j)<=geostatic.height(k)) then
					element(i).stress(2,j)=element(i).stress(2,j)+ &
						geostatic.weight(k)*(geostatic.height(k-1)- &
						geostatic.height(k))
				else
					element(i).stress(2,j)=element(i).stress(2,j)+ &
						geostatic.weight(k)*(geostatic.height(k-1)- &
						element(i).xygp(2,j))					
					islast=.true.
					exit
				end if
			end do
			if(islast) k=k+1
			element(i).stress(1,j)=element(i).stress(2,j)*geostatic.ko(k-1)
			element(i).stress(3,j)=element(i).stress(1,j)
			r1=1.0d0 !for axis-sysmetrical element
			if(element(I).ec==cax) r1=abs(element(i).xygp(1,j))
			bload(1:element(I).ndof)=bload(1:element(I).ndof)+ &
				matmul(element(i).stress(1:ELEMENT(I).ND,J),element(I).b(:,:,j))* &
				element(I).detjac(j)*ecp(element(I).et).weight(j)*r1
		end do
		TLOAD(ELEMENT(I).G)=TLOAD(ELEMENT(I).G)+bload(1:element(I).ndof)
	end do
	
end subroutine

subroutine bf_initialstress() !for cal_geo
	use solverds
	implicit none
	integer::i,j,k,nnum1
	real(8)::bf1=0,t1=0,w1=0,r1

	t1=0	
	do i=1,enum
		w1=material(element(i).mat).weight
		if(abs(w1)<1e-7) cycle
		nnum1=element(i).nnum
		BF1=0.D0
		do j=1,nnum1
			bf1=0.0
			do k=1,element(i).ngp
				r1=1.0d0 !for axis-sysmetrical element
				if(element(i).ec==cax) r1=abs(element(i).xygp(1,k))
				bf1=bf1+ecp(element(i).et).Lshape(j,k)* &
					ecp(element(i).et).weight(k)* &
					element(i).detjac(k)*r1*w1
			end do
			LOAD(NODE(ELEMENT(I).NODE(J)).DOF(NDIMENSION))=LOAD(NODE(ELEMENT(I).NODE(J)).DOF(NDIMENSION))+BF1

		end do
		
	end do
	
	
end subroutine

  