
subroutine mkl_solver_error(error)
	implicit none
	integer::error
	 ! Print an error message and exit
	WRITE(*,*) "Solver returned error code ", error
	STOP 	

end subroutine

INCLUDE 'mkl_dss.f90' 

subroutine solve_SLD()
	use mkl_dss
	use solverds
    use ifqwin
    USE IFCORE
	implicit none
	integer::i,j,iincs,iiter,istep=1,isubts=0,isref_spgcount=0
	integer::kref=0,dof1,NC1
	integer::nnslope=0,npslope=0
	real(kind=dpn)::minrelax=0.1d0,maxrelax=1.0d0
	
	logical::iscon=.false.,isfirstcall=.true.,isfirstcall2=.true.,isoscilated=.false.
	real(kind=dpn)::NormBL=0.0,resdis=0,t1,relax=1.0,convratio=0.0, &
                    normres=0.0,sumforce=0.0,TTime1=0.d0,R0,R1
	real(kind=dpn),ALLOCATABLE::bdylds(:),exdis(:),stepdis(:),stepstress(:,:,:),stepstrain(:,:,:),YFACUPDATE(:)
	logical,allocatable::IsBCDis(:)
    character(1)::key
	
	
	DOUBLE PRECISION START_TIME, STOP_TIME, DCLOCK
	EXTERNAL DCLOCK
	
	!variables for mkl solver.
	
	INTEGER :: error,nRhs=1
	INTEGER, PARAMETER :: bufLen = 20
	TYPE(MKL_DSS_HANDLE) :: handle ! Allocate storage for the solver handle.
	REAL(KIND=DPN),ALLOCATABLE::statOUt( : ),solution(:)
	CHARACTER*15 statIn
	INTEGER perm(1)
	INTEGER buff(bufLen)
	
	if(solver_control.ismkl) then
		! Initialize the solver.
		error = DSS_CREATE( handle, MKL_DSS_DEFAULTS )
		IF (error /= MKL_DSS_SUCCESS) call mkl_solver_error(error)
		! Define the non-zero structure of the matrix.
		if(solver_control.issym) then
			error = DSS_DEFINE_STRUCTURE( handle, MKL_DSS_SYMMETRIC, rowIndex, ndof, ndof, jcol, nnz )
		else
			error = DSS_DEFINE_STRUCTURE( handle, MKL_DSS_NON_SYMMETRIC, rowIndex, ndof, ndof, jcol, nnz )
		end if
		IF (error /= MKL_DSS_SUCCESS) call mkl_solver_error(error)
		! Reorder the matrix.
		error = DSS_REORDER( handle, MKL_DSS_DEFAULTS, perm )
		IF (error /= MKL_DSS_SUCCESS) call mkl_solver_error(error)	
		
		allocate(solution(ndof))
		solution=0.D0
	end if

		
	allocate(bdylds(ndof),exdis(ndof),stepload(ndof),stepdis(ndof),isBCdis(ndof),NI_NodalForce(ndof),&
				YFACUPDATE(NDOF),Tload(ndof))

	bdylds=0.0D0
	stepload=0.0D0
	Tload=0.0D0
	exdis=tdisp !用于存储上一荷载步的总位移。
	
	NormRes=1.d20
	istep=1	
	ttime1=0.0D0
	
	if(geostatic.isgeo) istep=0
	do iincs=istep,nstep
        !active element
        call element_activate(iincs)
        
        if(isexca2d/=0) call excavation(iincs)
        
		if(solver_control.solver==LPSOLVER.OR.solver_control.solver==MOSEK) then
			call lpsolvefile()
			exit 
		end if
		if(iincs==0) then
			if(geostatic.method==ko_geo) then
				call ko_initialstress()
				call outdata(iincs,iiter,iscon,isfirstcall,isubts)
				cycle
			else
				call bf_initialstress()
			end if
		end if
		
		!For SPG problem only.
		if(stepinfo(iincs).issteady) then
			!initialize Tstepdis
			Tstepdis(:,iincs)=Tstepdis(:,Stepinfo(iincs).matherstep)
		else
			!initial value(default is 0), borrow
			if(.not.allocated(inivaluedof)) allocate(inivaluedof(ndof))
			inivaluedof=Tstepdis(:,stepinfo(iincs).matherstep)
			call cmm_spg_cal(iincs)		
		end if
		
		
		
		do isubts=1,timestep(iincs).nsubts
			START_TIME = DCLOCK()
			ttime1=ttime1+timestep(iincs).subts(isubts)	
			
			iiter=0
            nnslope=0
            npslope=0
            isref_spgcount=0
            relax=1.0d0
			stepdis=0.0d0 !the incremental displacement of the current step.
            !if(iincs==4) solver_control.niteration=4
			do while(iiter<=solver_control.niteration)
				
				iiter=iiter+1
				call solver_initialization(kref,iincs,iiter)

				call incremental_load(iincs,iiter,solver_control.solver,isubts)

				load=load+bdylds
	!			t1=sum(load)
				call bc(iincs,iiter,load,stepdis(1:ndof),isubts)					
				
				
				if(kref==1) then
				
					call assemble_km(iincs)
					
					do j=1,bd_num
						if(bc_disp(j).isdead==1) cycle
						dof1=node(bc_disp(j).node).dof(bc_disp(j).dof)
						km(diaglkmloc(dof1))=UM
					end do
					
					call KM_UPDATE_SPG2(stepdis,iiter,iincs)

					if(solver_control.issym.and.(.not.solver_control.ismkl)) then
						!START_TIME = DCLOCK()
						call  chodec(km,bw,ndof)
						!STOP_TIME= DCLOCK()
						!write(99,20) iincs,iiter,STOP_TIME-START_TIME
					else
						! Factor the matrix.
						error = DSS_FACTOR_REAL( handle, MKL_DSS_POSITIVE_DEFINITE,km)
						IF (error /= MKL_DSS_SUCCESS) call mkl_solver_error(error)			
							
					end if
				end if	
				
			

				
				
				if(solver_control.issym.and.(.not.solver_control.ismkl)) then
					call chosol(km,bw,load,ndof,bwmax)
				else
					error = DSS_SOLVE_REAL(handle, MKL_DSS_DEFAULTS, load, nRhs, solution )
					IF (error /= MKL_DSS_SUCCESS) call mkl_solver_error(error)
					load=solution	
				end if
			   
				select case(solver_control.solver) 
					case(N_R)
						
						!if(SOLVER_CONTROL.ISLS==1) then							
                        CALL Linesearch(iincs,isubts,iiter,iscon,relax,stepdis,LOAD,bdylds)
						!endif
						
						stepdis=stepdis+load*relax
						
						!call Cal_UnbalanceForce(iincs,iiter,iscon,stepdis,bdylds,isubts) 

						if(.not.solver_control.isfc) then
							!call checon_sec(iscon,exdis,stepdis,ndof,solver_control.disp_tol,resdis,convratio,iiter)
							call CHECON_THD(iscon,STEPDIS,load,ndof,solver_control.disp_tol,resdis,sumforce,convratio,iiter)
							
							t1=dsqrt(dot_product(bdylds,bdylds))	
                        else
							!call checon_sec(iscon,pload,bdylds,ndof,solver_control.force_tol,resdis,convratio,iiter)
							call CHECON_THD(iscon,NI_NodalForce,bdylds,ndof,solver_control.FORCE_tol,resdis,sumforce,convratio,iiter)
							t1=resdis
							

						end if
						
						if(isref_spg==1.and.iscon.and.isref_spgcount<=10) then
                            
                            iscon=.false.
                            isref_spgcount=isref_spgcount+1
						endif
						!if(iscon.or.iiter==solver_control.niteration) load=stepdis						
						
					
                        
						!if(iiter>10) then
      !                      if(iiter<=100) then
						!		if(t1>1.05*NormRes) relax=max(relax-0.1,0.1)
						!		if(t1<0.90*NormRes) relax=min(relax+0.1,1.0) 
      !                      else
      !                          
      !                         
						!		if(t1>1.05*NormRes) relax=max(relax/2.0,minrelax)
      !                          
      !                          if(mod(iiter,100)==0) minrelax=max(minrelax/2.,0.001)
      !                          
						!	end if
						!
						!end if
						

						
						NormRes=t1                    
						
						!call Cal_AcceleratingFactor(iiter,iincs,iscon,stepdis,bdylds,relax,handle)
					
					case default
					
						stepdis=load
						call Cal_UnbalanceForce(iincs,iiter,iscon,stepdis,bdylds,isubts)
						
						!call checon(iscon,pvalue,load,ndof,solver_control.tolerance)
						if(solver_control.solver==LELASTIC) then
							iscon=.true.
						else
							if(.not.solver_control.isfc) then
								!call checon_sec(iscon,exdis,load,ndof,solver_control.disp_tol,resdis,convratio,iiter)
								call CHECON_THD(iscon,STEPDIS,load,ndof,solver_control.disp_tol,resdis,sumforce,convratio,iiter)
							end if
						end if					
				end select
				
				write(*,10) solver_control.isfc,KREF==1,isref_spg,convratio,sumforce,resdis, &
						relax*maxval(abs(load)),maxloc(abs(load)),relax,iiter,isubts,ttime1,iincs
				write(99,10) solver_control.isfc,KREF==1,isref_spg,convratio,sumforce,resdis, &
						relax*maxval(abs(load)),maxloc(abs(load)),relax,iiter,isubts,ttime1,iincs

				

				if((.not.iscon).and.iiter==solver_control.niteration) then
					write(*,30) iiter,isubts,ttime1,iincs
					write(99,30) iiter,isubts,ttime1,iincs
				end if
				
				if(iscon.or.iiter==solver_control.niteration) then
					STOP_TIME= DCLOCK()
					!call Cal_UnbalanceForce(iincs,iiter,iscon,stepdis,bdylds,isubts)
					write(*,40) iscon,Qinput,(Qstored-Qstorted_ini),(Qstored-Qstorted_ini)/Qinput,iiter,isubts,ttime1,iincs,STOP_TIME-START_TIME
					write(99,40) iscon,Qinput,(Qstored-Qstorted_ini),(Qstored-Qstorted_ini)/Qinput,iiter,isubts,ttime1,iincs,STOP_TIME-START_TIME
					load=stepdis
					exit
				end if
				
			end do	
            
            !各时间子步之间的位移更新。
			IF(SOLVER_CONTROL.TYPE==SPG) THEN
				Qstorted_ini=Qstored
				do concurrent (i=1:enum)
					if(element(i).isactive==1.and.allocated(element(i).sita_ini)) element(i).sita_ini=element(i).sita_fin
                end do
            else
                do concurrent (i=1:enum)
					if(element(i).isactive==1) then
						if(allocated(element(i).gforce)) element(i).gforce=element(i).gforce+element(i).Dgforce
						if(allocated(element(i).stress)) element(i).stress=element(i).stress+element(i).Dstress
						if(allocated(element(i).strain)) element(i).strain=element(i).strain+element(i).Dstrain
						if(allocated(element(i).pstrain)) element(i).pstrain=element(i).pstrain+element(i).evp
					end if
				end do
				

            END IF
            
			if(stepinfo(iincs).issteady) then

				IF(SOLVER_CONTROL.TYPE==SPG) then 
					Tstepdis(:,iincs)=stepdis
				else
					Tstepdis(:,iincs)=Tstepdis(:,iincs)+stepdis
				end if
				 
			else
				Tstepdis(:,iincs)=stepdis  !For a transient problem, tdisp is passed to stepdisp in the form of an initial value. 			
			    
				inivaluedof=Tstepdis(:,iincs) !!!!!
				
				 
			end if
			
			
			
			Tdisp=Tstepdis(:,iincs) !!for outdata
			if(isexca2d/=0) call Beam_Result_EXCA(iincs)			
			call outdata(iincs,iiter,iscon,isfirstcall,isubts)
		

		end do
		
	end do	
	
	CALL pointout_barfamily_diagram()
	
	! Deallocate solver storage and various local arrays.
    if(solver_control.ismkl) then
	    error = DSS_DELETE( handle, MKL_DSS_DEFAULTS )
	    IF (error /= MKL_DSS_SUCCESS ) call mkl_solver_error(error)
    end if
    
    if(isexca2d/=0) then
        call SCIPLOT()
        do while(ichar(key)/=ichar('q').and.ichar(key)/=27)
            key=getcharqq()
        enddo
		nc1 = setexitqq(QWIN$EXITNOPERSIST)
    endif
    
	if(allocated(solution)) deallocate(solution)
	if(allocated(bdylds)) deallocate(bdylds)
	if(allocated(YFACUPDATE)) deallocate(YFACUPDATE)
	if(allocated(exdis)) deallocate(exdis)	
	if(allocated(stepdis)) deallocate(stepdis)	
	if(allocated(NI_NodalForce)) deallocate(NI_NodalForce)
	if(allocated(stepload)) deallocate(stepload)
	if(allocated(km)) deallocate(km)
	if(allocated(load)) deallocate(load)
	if(allocated(tdisp)) deallocate(tdisp)
	if(allocated(irow)) deallocate(irow)
	if(allocated(jcol)) deallocate(jcol)
	if(allocated(Lmre)) deallocate(Lmre)
	if(allocated(adrn)) deallocate(adrn)
	if(allocated(ROWINDEX)) deallocate(ROWINDEX)
	if(allocated(bw)) deallocate(bw)
    if(allocated(node)) deallocate(node)
	if(allocated(element)) deallocate(element)
    if(allocated(inivaluedof)) deallocate(inivaluedof)
    if(allocated(Tstepdis)) deallocate(Tstepdis)

10 format('ISFC=',L1,',ISFAC=',L1,',IsRef_SPG=',L1 ',ConvCoff.=',f6.3,',SumForce.=', E10.3,',SumRes=',E10.3 ',MaxDiff.=',f7.3,'(N=',I7,'),R.F.=',f8.5, &
				',NIter=',I4,',NSubTS(E.Time)=',I4,'(',F8.3,'),NIncr=',I2,'.')	
20 format('TOTAL REFACTORIZATION. NINCR=',I4,' NITE=',I4,' Duration=',f12.6) 
21 format('PARTIAL REFACTORIZATION. NINCR=',I4,' NITE=',I4,' Duration=',f12.6,' NUMBERS OF BCs UPDATED=',I7) 
30 format('Exit with the Max iteration without convergence. Niteration=',I4,',NSubTS(E.Time)=',I3,'(',F10.3,'), NIncr=',I3,'.')
40 format('Isconverged=',L1,'.Qinput=',E15.7,',Qstored=',E15.7,',Massbalance=',E15.7,',Niteration=',I4,',NSubTS(E.Time)=',I3,'(',F10.3,'), NIncr=',I3,'TIME TAKEN=',E15.7,'.')
end subroutine

subroutine Linesearch(iincs,isubts,iiter,iscon,relax,stepdis,Ddis,bdylds)
	use solverlib 	
	implicit none
	integer,intent(in)::iincs,isubts,iiter
    LOGICAL,INTENT(IN)::iscon
	real(kind=DPN),intent(in)::stepdis(ndof),Ddis(ndof)
	real(kind=DPN)::relax,bdylds(ndof)
	real(kinD=DPN)::stepdis1(ndof)
    INTEGER::ILS1=0,INEG1,IPOS1,N1,N2,I
    real(kind=DPN)::ETA1(SOLVER_CONTROL.NLS),RATIO1(SOLVER_CONTROL.NLS),T1,ETAALT1,S0,T2
	
   
	stepdis1=stepdis+RELAX*Ddis	
	call Cal_UnbalanceForce(iincs,iiter,iscon,stepdis1,bdylds,isubts)
    T1=dot_product(bdylds,Ddis)
    IF(IITER<30) THEN
        SOLVER_CONTROL.S0=T1
        !RELAX=1.0D0
        RETURN
    ELSEIF(IITER>=50) THEN
        IF(MOD(IITER,50)==0) RELAX=MAX(RELAX/2.0D0,0.01)
        RETURN
    ENDIF
    IF(SOLVER_CONTROL.ISLS==0) RETURN
    
    
    ETA1=0.D0
    RATIO1=0.D0
        
    RATIO1(1)=1.D0
    ETA1(2)=RELAX   
	RATIO1(2)=T1/SOLVER_CONTROL.S0
    ILS1=2

    DO WHILE(ABS(RATIO1(ILS1))>SOLVER_CONTROL.STOL.AND.ILS1<SOLVER_CONTROL.NLS)
        !INEG1=MINLOC(ETA1,MASK=.FALSE.,DIM=1)
        !INEG1=MINLOC(ETA1,MASK=(RATIO1<0.D0),DIM=1)
        INEG1=0
        IPOS1=0
        DO I=1,ILS1
            T1=1D20
            IF(RATIO1(I)<0.D0.AND.ETA1(I)<T1) THEN
                T1=ETA1(I)
                INEG1=I
            ENDIF
        ENDDO
        
        IF(INEG1/=0) THEN
            !IPOS1=MAXLOC(ETA1,MASK=(RATIO1>0.D0.AND.ETA1<ETA1(INEG1)),DIM=1)
            DO I=1,ILS1
                T1=-1D20
                IF(RATIO1(I)>0.D0.AND.ETA1(I)<ETA1(INEG1).AND.ETA1(1)>T1) THEN
                    T1=ETA1(I)
                    IPOS1=I
                ENDIF
            ENDDO            
        ENDIF
        IF(INEG1*IPOS1/=0) THEN
            RELAX=(RATIO1(IPOS1)*ETA1(INEG1)-RATIO1(INEG1)*ETA1(IPOS1))/(RATIO1(IPOS1)-RATIO1(INEG1))
            ETAALT1=ETA1(IPOS1)+0.2D0*(ETA1(INEG1)-ETA1(IPOS1))
            RELAX=MAXVAL((/RELAX,ETAALT1,0.1D0/))
            !IF(ETA1(ILS1)<0.1D0) ETA1(ILS1)=0.1
        ELSE
            N1=ILS1
            N2=ILS1-1
            RELAX=(RATIO1(N2)*ETA1(N1)-RATIO1(N1)*ETA1(N2))/(RATIO1(N2)-RATIO1(N1))            
        ENDIF
        
        IF(RELAX<0.1D0) RELAX=0.1
        IF(RELAX>10D0) RELAX=10D0
        
        ILS1=ILS1+1
        ETA1(ILS1)=RELAX
        stepdis1=stepdis+RELAX*Ddis	
        
	    call Cal_UnbalanceForce(iincs,iiter,iscon,stepdis1,bdylds,isubts)
        T1=dot_product(bdylds,Ddis)
        RATIO1(ILS1)=T1/SOLVER_CONTROL.S0
    ENDDO
    
	SOLVER_CONTROL.S0=T1
    
endsubroutine

subroutine KM_UPDATE_SPG2(stepdis,iiter,IINCS)
	USE SOLVERDS
	implicit none
	real(8),intent(in)::stepdis(ndof)
	integer,intent(in)::iiter,IINCS
	integer::J,DOF1
	real(8)::t1
	
		
	do j=1,numNseep
		if(sf(NSEEP(J).sf).factor(IINCS)==-999.D0) cycle
		
		if(Nseep(j).isdual>0) then
			if(bc_disp(Nseep(j).isdual).isdead==0) cycle
		end if
		
		if(Nseep(j).isdead==0) then
			dof1=node(Nseep(j).node).dof(Nseep(j).dof)														  
			km(diaglkmloc(dof1))=UM
			!load(dof1)=(Nseep(j).value-stepdis(dof1))*sf(Nseep(j).sf).factor(iincs)*UM
		end if
	end do
	
	!!!lacy method
	!IF(IITER>1) THEN
	!	DO J=1,NNUM
	!		dof1=node(J).dof(4)
	!		t1=minNPH+node(j).coord(ndimension)
	!		IF(STEPDIS(DOF1)<t1) then
	!			km(bw(dof1))=UM
	!		end if		
	!	END DO
	!!	
	!!	do j=1,numNseep
	!!		dof1=node(Nseep(j).node).dof(Nseep(j).dof)
	!!		t1=STEPDIS(DOF1)-node(Nseep(j).node).coord(ndimension)
	!!		if(t1>0) then
	!!			km(bw(dof1))=UM
	!!			 !load(dof1)=(Nseep(j).value-stepdis(dof1))*sf(Nseep(j).sf).factor(iincs)*UM
	!!		end if			
	!!	end do
	!END IF
END SUBROUTINE


subroutine Cal_UnbalanceForce(iincs,iiter,iscon,stepdis,bdylds,isubts)
	use solverds
	implicit none
	integer,intent(in)::iincs,iiter,isubts
	logical,intent(in)::iscon
	real(kind=DPN),intent(in)::stepdis(ndof)
	real(kind=DPN)::bdylds(ndof)
	
	select case(solver_control.bfgm)
		case(viscop)
			call bload_viscoplasticity(iiter,iscon,bdylds)
		case(inistress)
			call bload_inistress(iiter,iscon,bdylds)
		case(consistent,constant,continuum,iniflux)
			!print *, iincs,iiter,resdis
			call bload_consistent(iiter,iscon,bdylds,stepdis,iincs,isubts)
		!case(iniflux)
		!	call bload_iniflux(iiter,iscon,bdylds,stepdis,iincs)
		case default
			print *, 'No Such BODY FORCE GENERATION METHOD:',SOLVER_CONTROL.BFGM
	end select	
end subroutine


subroutine Cal_AcceleratingFactor(iiter,iincs,iscon,stepdis,bdylds,relax,handle,isubts)
!Sloan, S.W., D. Sheng, and A.J. Abbo, Accelerated initial stiffness schemes for elastoplasticity.
!International Journal for Numerical and Analytical Methods in Geomechanics, 2000. 24(6): p. 579-599.
	use solverds
	use mkl_dss
	implicit none
	integer,intent(in)::iincs,iiter,isubts
	logical,intent(in)::iscon
	real(KIND=DPN),intent(in)::stepdis(ndof)
	real(KIND=DPN),intent(in out)::bdylds(ndof),relax
	integer::j,dof1
	
	INTEGER :: error,nRhs=1
	TYPE(MKL_DSS_HANDLE):: handle ! Allocate storage for the solver handle.
	REAL(KIND=DPN)::solution(ndof)

	
	call Cal_UnbalanceForce(iincs,iiter,iscon,stepdis,bdylds,isubts)

	call bc(iincs,iiter,bdylds,stepdis,isubts)
	
	if(isref_spg==1) then
		call assemble_km(iincs)

		do j=1,bd_num
			if(bc_disp(j).isdead==1) cycle
			dof1=node(bc_disp(j).node).dof(bc_disp(j).dof)
			km(diaglkmloc(dof1))=UM
		end do
		
		call KM_UPDATE_SPG2(stepdis,iiter,iincs)
		
		
		if(solver_control.issym.and.(.not.solver_control.ismkl)) then
			call  chodec(km,bw,ndof)
			!write(99,20) iincs,iiter,STOP_TIME-START_TIME
		else
			! Factor the matrix.
			error = DSS_FACTOR_REAL(handle, MKL_DSS_POSITIVE_DEFINITE,km)
			IF (error /= MKL_DSS_SUCCESS) call mkl_solver_error(error)			
				
		end if

	end if	

	if(solver_control.issym.and.(.not.solver_control.ismkl)) then
		call chosol(km,bw,load,ndof,bwmax)
	else
		error = DSS_SOLVE_REAL(handle, MKL_DSS_DEFAULTS, load, nRhs, solution )
		IF (error /= MKL_DSS_SUCCESS) call mkl_solver_error(error)
		load=solution	
	end if
	
	relax=relax+dot_product(load,bdylds)/dot_product(load,load) 


end subroutine

subroutine Continuum_stress_update(iiter,iscon,istep,ienum,bload,Ddis,nbload)

	use solverds
	!use operation_ix
	implicit none
	logical,intent(in)::iscon
	integer,intent(in)::iiter,ienum,istep,nbload
	real(kind=DPN),intent(in)::Ddis(nbload)
	real(kind=DPN),intent(out)::bload(nbload)
    
	real(kind=DPN)::stress1(6)=0.0,strain1(6)=0.0, &
					dee(4,4)=0.0,sigma(6)=0.0,inv(3)=0.0,vyf=0.0, &
					devp(6)=0.0,m(6,3)=0.0,dywi(3)=0.0,dyf(7)=0.0, &
					evp(6),vyf_old=0.0,fac=0.0, &
					stress2(6)=0.0,dqwi(3)=0.0,dqf(7)=0.0,dp(6,6)=0.0,&
					de(6,6)=0.0,t1=0.0,dsigma(7)=0.0,lamda=0.0,&
					buf1(7)=0.0,buf2(7)=0.0,pstrain(6)=0.0,dlamda=0.0,&
					Rmat(7,7)=0.0,Rpstrain(7)=0,r1=1.0,hj=0.0,gforce1(100)=0.0
	real(kind=DPN)::vyf2=0,vyf3=0,dqf2(7)=0,dqf3(7)=0,dyf2(7)=0,dyf3(7)=0, &
		c(3,3)=0,lamda2=0,dlamda2=0,lamda3=0,dlamda3=0,stressb(6)=0, &
		invb(3)=0,dqwi2(3)=0,dqwi3(3)=0,ev=0,e=0,DPMWS(6)=0,DPMWEV=0, &
		PLM=0,DYWEV=0
	integer::nd1=0
    integer::j,k,n1,n2=0,ayf=0,ayf1=0,i1,j1,k1

	n1=element(ienum).ngp
	nd1=element(ienum).nd
	bload=0.0
 
	if(solver_control.bfgm==continuum) element(ienum).km=0
 
	do j=1,n1					
		strain1=0.0
		stress1=0.0
		pstrain=element(ienum).pstrain(:,j)
		ev=element(ienum).ev(j)
		!strain of this incrementals
		!UN(1:element(ienum).ndof)=DDIS(ELEMENT(IENUM).G)
		strain1(1:nd1)=matmul(element(ienum).b(1:nd1,:,j),DDIS(1:element(ienum).ndof))
		e=element(ienum).e(j)+sum(strain1(1:3))
		!stress incremental
		de(1:nd1,1:nd1)=element(ienum).d
		sigma(1:nd1)=matmul(de(1:nd1,1:nd1),strain1(1:nd1))
		stress1=element(ienum).stress(:,j)+sigma
		!check the stress resulted in the last interation whether it penetrate the yield surface.
		!ayf=0
		!call invanrant(stress1,inv)
		!call yieldfun(vyf,element(ienum).mat,inv,ayf,ev)
		!n2=0
		!lamda=0.0
		!dlamda=0.0
		!lamda2=0.0
		!dlamda2=0.0
		!lamda3=0.0
		!dlamda3=0.0
		!stressb=stress1
		!invb=inv
		!do while(vyf>=0.and.n2<=20)
		!	n2=n2+1
		!	call invanrant(stress1,inv)
		!	call yieldfun(vyf,element(ienum).mat,inv,ayf,ev)
 
		!	if(abs(vyf)<solver_control.ftol.or.n2+1>20.or.&
		!		(material(element(ienum).mat).type==MC.and.n2>1)) then
		!		
		!		if(solver_control.bfgm==constant) exit
		!		!algorithmic tangent modulus ,Depc
		!		select case(material(element(ienum).mat).type) 
		!			case(MC) !updated stress in trial position,ref to crisfield book.
		!				call MC_Depc(ayf,nd1,ienum,lamda,lamda2,lamda3,dqf,dqf2,dqf3,dyf,dyf2,dyf3, &
		!						dqwi,dqwi2,dqwi3,stressb,invb,m,Rmat,de,c)
		!			case default
		!				call FQderivatives(element(ienum).mat,ayf,ev,e,stress1,inv,m,dywi, &
		!								dyf,dqwi,dqf,PLM,DPMWS,DPMWEV,DYWEV)			
		!				call formRmat(Rmat,stress1,inv,m,element(ienum).mat,nd1,lamda, &
		!					lamda2,lamda3,element(ienum).d,dqwi,dqwi2,dqwi3,ayf,e,DPMWS, &
		!					DPMWEV)
		!				dyf(nd1+1)=DYWEV
		!				dqf(nd1+1)=PLM								
		!				buf1(1:nd1+1)=matmul(Rmat(1:nd1+1,1:nd1+1),dqf(1:nd1+1))
		!				buf2(1:nd1+1)=matmul(dyf(1:nd1+1),Rmat(1:nd1+1,1:nd1+1))
		!				de(1:nd1+1,1:nd1+1)=Rmat(1:nd1+1,1:nd1+1)-csproduct(buf1(1:nd1+1), &
		!					buf2(1:nd1+1))/dot_product(buf2(1:nd1+1),dqf(1:nd1+1))
		!		end select
		!		
		!		exit
		!	end if	
 
		!	call FQderivatives(element(ienum).mat,ayf,ev,e,stress1,inv,m,dywi, &
		!					dyf,dqwi,dqf,PLM,DPMWS,DPMWEV,DYWEV)				
		!	
		!	select case(material(element(ienum).mat).type)
		!		case(MC)
		!			call MC_UPDATE(ayf,nd1,ienum,vyf,ev,dlamda,dlamda2,dlamda3,lamda,lamda2,&
		!					lamda3,dqf,dqf2,dqf3,dyf,dyf2,dyf3,dqwi,dqwi2,dqwi3,dsigma,inv,m, &
		!					Rmat,de,c)					
		!		case default
		!			call formRmat(Rmat,stress1,inv,m,element(ienum).mat,nd1,lamda, &
		!			lamda2,lamda3,element(ienum).d,dqwi,dqwi2,dqwi3,ayf,e,DPMWS, &
		!			DPMWEV)
		!			Rpstrain(1:nd1)=lamda*dqf(1:nd1)-pstrain(1:nd1)+ &
		!				element(ienum).pstrain(1:nd1,j)
		!			Rpstrain(nd1+1)=lamda*PLM*-ev+element(ienum).ev(j)
		!			dyf(nd1+1)=DYWEV
		!			dqf(nd1+1)=PLM
		!			dyf(1:nd1+1)=matmul(dyf(1:nd1+1),Rmat(1:nd1+1,1:nd1+1))
		!			dlamda=(vyf-dot_product(dyf(1:nd1+1),Rpstrain(1:nd1+1)))/ &
		!				dot_product(dyf(1:nd1+1),dqf(1:nd1+1))
		!			dsigma(1:nd1+1)=-matmul(Rmat(1:nd1+1,1:nd1+1),Rpstrain(1:nd1+1))-dlamda* &
		!				matmul(Rmat(1:nd1+1,1:nd1+1),dqf(1:nd1+1))						
		!	end select
		!	!update 
		!	pstrain(1:nd1)=pstrain(1:nd1)-((element(ienum).D).ix.dsigma(1:nd1))
		!	stress1(1:nd1)=stress1(1:nd1)+dsigma(1:nd1)
		!	lamda=lamda+dlamda
		!	ev=ev+dsigma(nd1+1)	
		!			
		!end do
		
		if(solver_control.bfgm==continuum) then
			r1=1.0 !for axis-sysmetrical element
			if(element(ienum).ec==cax) r1=element(ienum).xygp(1,j)
			element(ienum).km=element(ienum).km+matmul(matmul( &
				transpose(element(ienum).b(:,:,j)),de(1:nd1,1:nd1)),element(ienum).b(:,:,j)) &
				*element(ienum).detjac(j)*ecp(element(ienum).et).weight(j)*r1
		end if
		!stress incremental of this incremental
		element(ienum).Dstress(1:nd1,j)=stress1(1:nd1)-element(ienum).stress(1:nd1,j)
		element(ienum).Dstrain(1:nd1,j)=strain1(1:nd1)
		element(ienum).evp(:,j)=pstrain-element(ienum).pstrain(:,j)
		element(ienum).ev(j)=ev
		element(ienum).e(j)=e
		bload(1:element(ienum).ndof)=bload(1:element(ienum).ndof)+ &
								matmul(stress1(1:nd1),element(ienum).b(:,:,j))* &
								element(ienum).detjac(j)*ecp(element(ienum).et).weight(j)
		!if(iscon.or.iiter==solver_control.niteration) then
		!	element(ienum).stress(:,j)=stress1
		!	element(ienum).strain(:,j)=element(ienum).strain(:,j)+strain1
		!	element(ienum).evp(:,j)=pstrain-element(ienum).pstrain(:,j)
		!	!element(ienum).pstrain(:,j)=element(ienum).pstrain(:,j)+ pstrain
		!	element(ienum).pstrain(:,j)=pstrain
		!	element(ienum).ev(j)=ev
		!	element(ienum).e(j)=e
		!end if			
	end do

end subroutine

subroutine FQderivatives(mat,ayf,ev,e,stress1,inv,m,dywi,dyf,dqwi,dqf, &
	PLM,DPMWS,DPMWEV,DYWEV)
	implicit none
	integer::mat,ayf
	real(8)::ev,e,stress1(6),inv(3),m(6,3),dywi(3),dqwi(3),dyf(7),dqf(7), &
		PLM,DPMWS(6),DPMWEV,DYWEV
	
	call deriv_sinv(m,stress1)
	call deriv_yf_with_inv(dywi,inv,mat,ayf,ev)
	call deriv_yf(dyf(1:6),dywi,m)
	call deriv_qf_with_inv(dqwi,inv,mat,ayf,ev)
	call deriv_yf(dqf(1:6),dqwi,m)
	call PlasticModuli(PLM,inv,mat,e,ev)
	call deri_PM_W_S(DPMWS,inv,mat,e,ev)
	call deri_PM_W_EV(DPMWEV,inv,mat,e,ev)
	call deri_yf_w_EV(DYWEV,inv,mat)	
	
end subroutine

subroutine formRmat(Rmat,sigma,inv,fstM,mat,nd,lamda,lamda2,lamda3,de,& 
	fdqwi,fdqwi2,fdqwi3,ayf,e,DPMWS,DPMWEV)
	use solverds
!	use operation_ix
!	use operation_i
	implicit none
	integer::mat,nd,i,ayf,ayf1
	real(8)::Rmat(7,7),secdqf(6,6)=0,inv(3),sigma(6),&
		lamda,de(nd,nd),buf1(6,6)=0,fdqwi(3),fstM(6,3),lamda2,lamda3, &
		fdqwi2(3),fdqwi3(3),SDQWEV(6)=0.0,DPMWS(6),DPMWEV,e
!	
!	Rmat=0.0	
!	select case(material(mat).type)
!		case(mises)
!			!call sec_deriv_qf(secdqf,sigma,inv,mat,fdqwi,fdinv,activeyf)
!			Rmat(1,1)=2
!			Rmat(2,2)=2
!			Rmat(3,3)=2
!			Rmat(4,4)=6
!			Rmat(5,5)=6
!			Rmat(6,6)=6	
!			Rmat(1,2)=-1
!			Rmat(2,3)=-1
!			Rmat(1,3)=-1
!			Rmat(2,1)=-1
!			Rmat(3,1)=-1
!			Rmat(3,2)=-1				
!			Rmat=0.5*Rmat/inv(2)
!			do i=1,6
!				buf1(i,:)=fstM(i,2)*fstM(:,2)
!			end do
!			buf1=buf1*2.25/inv(2)**3
!			secdqf=Rmat(1:6,1:6)-buf1
!			Rmat=0.0
!			Rmat(1:nd,1:nd)=lamda*matmul(de(1:nd,1:nd), &
!				secdqf(1:nd,1:nd))
!			do i=1,nd
!				Rmat(i,i)=Rmat(i,i)+1.0
!			end do
!			Rmat(1:nd,1:nd)=Rmat(1:nd,1:nd).ix.de(1:nd,1:nd)
!			return
!		case(MC) !in MC R is the T in the crisfield text. vol2, equ(14.83), P121 
!			select case(ayf)
!				case(12,23)
!					ayf1=0
!					call sec_deriv_qf(secdqf,sigma,inv,mat,fdqwi,fstM,ayf1)
!					Rmat(1:nd,1:nd)=-lamda*matmul(secdqf(1:nd,1:nd),de(1:nd,1:nd))
!					ayf1=ayf
!					call sec_deriv_qf(secdqf,sigma,inv,mat,fdqwi2,fstM,ayf1)
!					Rmat(1:nd,1:nd)=Rmat(1:nd,1:nd)-lamda2*matmul(secdqf(1:nd,1:nd),de(1:nd,1:nd))
!				case(123)
!					ayf1=0
!					call sec_deriv_qf(secdqf,sigma,inv,mat,fdqwi,fstM,ayf1)
!					Rmat(1:nd,1:nd)=-lamda*matmul(secdqf(1:nd,1:nd),de(1:nd,1:nd))
!					ayf1=12
!					call sec_deriv_qf(secdqf,sigma,inv,mat,fdqwi2,fstM,ayf1)
!					Rmat(1:nd,1:nd)=Rmat(1:nd,1:nd)-lamda2*matmul(secdqf(1:nd,1:nd),de(1:nd,1:nd))
!					ayf1=23
!					call sec_deriv_qf(secdqf,sigma,inv,mat,fdqwi3,fstM,ayf1)
!					Rmat(1:nd,1:nd)=Rmat(1:nd,1:nd)-lamda3*matmul(secdqf(1:nd,1:nd),de(1:nd,1:nd))
!				case default
!					call sec_deriv_qf(secdqf,sigma,inv,mat,fdqwi,fstM,ayf)
!					Rmat(1:nd,1:nd)=-lamda*matmul(secdqf(1:nd,1:nd),de(1:nd,1:nd))
!			end select
!			do i=1,nd
!				Rmat(i,i)=Rmat(i,i)+1.0
!			end do
!!			write(20,'(4e15.7)') (rmat(i,1:4),i=1,4)
!			return
!		case(camclay)
!			call sec_deriv_qf(secdqf,sigma,inv,mat,fdqwi,fstM,ayf)
!			call Sec_deri_qf_w_EV(SDQWEV,inv,mat)
!			Rmat(1:nd,1:nd)=(.i.de(1:nd,1:nd))+lamda*secdqf(1:nd,1:nd)
!			Rmat(nd+1,1:nd)=lamda*DPMWS(1:nd)
!			Rmat(1:nd,nd+1)=lamda*SDQWEV(1:nd)
!			Rmat(nd+1,nd+1)=-1+lamda*DPMWEV
!			Rmat(1:nd+1,1:nd+1)=(.i.Rmat(1:nd+1,1:nd+1))
			
!	end select

	
end subroutine

!yield the second derivative plastic potential function
subroutine sec_deriv_qf(Secdqf,sigma,inv,mat,fdqwi,fstM,activeyf)

	use solverds
	implicit none
	integer::i,j,mat,activeyf,ii,ji,ie,je
	real(8)::Secdqf(6,6),inv(3),sigma(6),secM(6,6,3)=0,fdqwi(3),fstM(6,3), &
	sdqwi(3,3)=0.0

	secdqf=0.0
!	cpm=0.0
	ii=1
	ji=1
	ie=3
	je=3
	select case(material(mat).type)
		case(mises)
			ii=2
			ji=2
			ie=2
			je=2			
		case(mc)
			ii=2
			ji=2
			ie=3
			je=3
		case(camclay)
			ii=1
			ji=1
			ie=1
			je=1		
	end select

	call second_deriv_sinv_2nd(secM,sigma,inv) 
	call sec_deriv_qfwi(sdqwi,inv,mat,activeyf)
	secdqf=fdqwi(2)*secM(:,:,2)+fdqwi(3)*secM(:,:,3)
	do i=ii,ie
		do j=ji,je
!			cpm(:,:,i,j)=csproduct(fstM(:,i),fstm(:,j))
			secdqf=secdqf+csproduct(fstM(:,i),fstm(:,j))*sdqwi(i,j)
!			m22s=csproduct(fstM(:,2),fstm(:,2))
!			m23s=csproduct(fstM(:,2),fstm(:,3))
!			m32s=csproduct(fstM(:,3),fstm(:,2))
!			m33s=csproduct(fstM(:,3),fstm(:,3))
		end do
	end do
!	secdqf=fdqwi(2)*secM(:,:,2)+fdqwi(3)*secM(:,:,3)+sdqwi(2,2)*m22s+ &
!		sdqwi(2,3)*(m23s+m32s)+sdqwi(3,3)*m33s
	
	return		
end subroutine

!according to the stress vector(stress)calculate the second derivertives of the p,J2,J3
!stress(1-6)=sx,sy,sz,sxy,syz,sxz.
subroutine second_deriv_sinv_2nd(secM,stress,inv)
	implicit none
	integer::i,j,p,q,m,n
	integer,external::delta
	real(8)::secm(6,6,3),stress(6),inv(3),devs(3,3)=0	
	
		
	!devs=stress
	devs(1,1)=2.0/3.0*stress(1)-1.0/3.0*(stress(2)+stress(3))
	devs(2,2)=2.0/3.0*stress(2)-1.0/3.0*(stress(1)+stress(3))
	devs(3,3)=2.0/3.0*stress(3)-1.0/3.0*(stress(2)+stress(1))
	devs(1,2)=stress(4)
	devs(2,1)=stress(4)
	devs(2,3)=stress(5)
	devs(3,2)=stress(5)
	devs(1,3)=stress(6)
	devs(3,1)=stress(6)

	secM=0
	secM(1,1,2)=2
	secM(2,2,2)=2
	secM(3,3,2)=2
	secM(4,4,2)=6
	secM(5,5,2)=6
	secM(6,6,2)=6	
	secM(1,2,2)=-1
	secM(2,3,2)=-1
	secM(1,3,2)=-1
	secM(2,1,2)=-1
	secM(3,1,2)=-1
	secM(3,2,2)=-1
	secM(:,:,2)=secM(:,:,2)*1.0/3.0

	secM(1,1,3)=devs(1,1)
	secM(2,2,3)=devs(2,2)
	secM(3,3,3)=devs(3,3)
	secM(4,4,3)=-3.0*devs(3,3)
	secM(5,5,3)=-3.0*devs(1,1)
	secM(6,6,3)=-3.0*devs(2,2)
	secM(2,1,3)=devs(3,3)
	secM(3,1,3)=devs(2,2)
	secM(3,2,3)=devs(1,1)		
	secM(4,1,3)=devs(1,2)
	secM(4,2,3)=devs(1,2)
	secM(4,3,3)=-2.0*devs(1,2)
	secM(5,1,3)=-2.0*devs(2,3)
	secM(5,2,3)=devs(2,3)
	secM(5,3,3)=devs(2,3)
	secM(5,4,3)=3*devs(1,3)
	secM(6,1,3)=devs(1,3)
	secM(6,2,3)=2*devs(1,3)
	secM(6,3,3)=devs(1,3)
	secM(6,4,3)=3*devs(2,3)
	secM(6,5,3)=3*devs(1,3)
	do i=1,6
		do j=i+1,6
			secM(i,j,3)=secM(j,i,3)
		end do
	end do
	secM(:,:,3)=2.0/3.0*secM(:,:,3)
	
end subroutine



!second derivative of the potential function with invarants p,J2,J3
subroutine sec_deriv_qfwi(secdqfwi,inv,mat,activeyf)
	use solverds
	implicit none
	integer::mat,activeyf
	logical::isqf=.true.
	real(8)::Secdqfwi(3,3),inv(3),vcos,vtan,vsin,J2,mcA(3),c4,c1
	
	secdqfwi=0.0
	
	select case(material(mat).type)
		case(mises)
			secdqfwi(2,2)=-2.25/inv(2)**3
		case(mc)
			J2=inv(2)**2/3.0
			vsin=dsin(material(mat).property(5)/180.0*pi())
			if(abs(dsin(inv(3)))<0.49) then
				call Mc_A(mcA,inv,mat,activeyf,isqf)
				vcos=dcos(3*inv(3))
				vtan=dtan(3*inv(3))
				c4=mcA(3)+3*vtan*mcA(2)
				secdqfwi(2,2)=-(mcA(1)-vtan**2*c4-3*vtan*mcA(2))/(4*J2**1.5)
				secdqfwi(3,3)=3*c4/(4*J2**2.5*vcos**2)
				secdqfwi(2,3)=(0.5*vtan*c4+mcA(2))*3**0.5/(2*J2**2*vcos)
				secdqfwi(3,2)=secdqfwi(2,3)
			else
				c1=1.0
				if(dsin(inv(3))<0) c1=-1.0
				secdqfwi(2,2)=(3.0-c1*vsin)/(-8*3**0.5*J2**1.5)				
			end if
		case(camclay)
			secdqfwi(1,1)=-2*material(mat).property(1)**2			
	end select
	
end subroutine

subroutine deri_yf_w_EV(dYWEV,inv,mat)
	use solverds
	implicit none
	real(8)::dYWEV,inv(3)
	integer::mat
	
	dYWEV=0
	select case(material(mat).type)
		case(camclay)
			dYWEV=material(mat).property(1)**2*inv(1)
	end select 
	
end subroutine

!calculate the derivative of the gradient of potential function with evolution variable.
subroutine Sec_deri_qf_w_EV(SDQWEV,inv,mat)
	use solverds
	implicit none
	real(8)::SDQWEV(6),inv(3)
	integer::mat
	
	SDQWEV=0
	select case(material(mat).type)
		case(camclay)
			SDQWEV(1:3)=1.0/3.0*material(mat).property(1)**2
	end select 
	
end subroutine

!calculate the derivative of the plastic modulus with stress
subroutine deri_PM_W_S(DPMWS,inv,mat,e,ev)
	use solverds
	implicit none
	real(8)::DPMWS(6),inv(3),e,ev
	integer::mat
	
	DPMWS=0
	select case(material(mat).type)
		case(camclay)
			DPMWS(1:3)=2.0/3.0*material(mat).property(1)**2*(1+e)*ev/ &
				(material(mat).property(3)-material(mat).property(4))
	end select 
	
end subroutine

!calculate the derivative of the plastic modulus with evolution variable
subroutine deri_PM_W_EV(DPMWEV,inv,mat,e,ev)
	use solverds
	implicit none
	real(8)::DPMWEV,inv(3),e,ev
	integer::mat
	
	DPMWEV=0
	select case(material(mat).type)
		case(camclay)
			DPMWEV=2.0*material(mat).property(1)**2*(1+e)*(inv(1)-ev)/ &
				(material(mat).property(3)-material(mat).property(4))
	end select 
	
end subroutine

!calculate the hardening moduli
subroutine PlasticModuli(PLM,inv,mat,e,ev)
	use solverds
	implicit none
	real(8)::PLM,inv(3),e,ev
	integer::mat
	
	PLM=0
	select case(material(mat).type)
		case(camclay)
			PLM=material(mat).property(1)**2*(1+e)*ev*(2*inv(1)-ev)/ &
				(material(mat).property(3)-material(mat).property(4))
	end select 	
end subroutine

subroutine MC_Depc(ayf,nd1,ienum,lamda,lamda2,lamda3,dqf,dqf2,dqf3,dyf,dyf2,dyf3, &
	dqwi,dqwi2,dqwi3,stressb,invb,m,Rmat,de,c)
	use solverds
	implicit none
	integer::i1,j1
	integer::nd1,ienum,ayf
	real(8)::lamda,lamda2,lamda3,dqf(7),dqf2(7),dqf3(7),dyf(7), &
		dyf2(7),dyf3(7),dqwi(3),dqwi2(3),dqwi3(3),Rmat(7,7), &
		de(6,6),stressb(6),invb(3),m(6,3),c(3,3)
	real(8)::MCQ(6,3)=0,MAC(6,3)=0,DPMWS(6)=0,e=0,DPMWEV=0, &
		buf1(6)=0,buf2(6)=0
			
	select case(ayf)
		case(12,23)
			MCQ(1:nd1,1)=matmul(element(ienum).d,dqf(1:nd1))
			MCQ(1:nd1,2)=matmul(element(ienum).d,dqf2(1:nd1))
			MAC(1:nd1,1)=matmul(dyf(1:nd1),element(ienum).d)
			MAC(1:nd1,2)=matmul(dyf2(1:nd1),element(ienum).d)
			de(1:nd1,1:nd1)=element(ienum).d
			do i1=1,2
				do j1=1,2
					de(1:nd1,1:nd1)=de(1:nd1,1:nd1)-csproduct(MCQ(1:nd1,i1), &
						MAC(1:nd1,j1))*c(i1,j1)
				end do
			end do										
		case(123)
			MCQ(1:nd1,1)=matmul(element(ienum).d,dqf(1:nd1))
			MCQ(1:nd1,2)=matmul(element(ienum).d,dqf2(1:nd1))
			MCQ(1:nd1,3)=matmul(element(ienum).d,dqf3(1:nd1))
			MAC(1:nd1,1)=matmul(dyf(1:nd1),element(ienum).d)
			MAC(1:nd1,2)=matmul(dyf2(1:nd1),element(ienum).d)
			MAC(1:nd1,3)=matmul(dyf3(1:nd1),element(ienum).d)
			de(1:nd1,1:nd1)=element(ienum).d
			do i1=1,3
				do j1=1,3
					de(1:nd1,1:nd1)=de(1:nd1,1:nd1)-csproduct(MCQ(1:nd1,i1), &
						MAC(1:nd1,j1))*c(i1,j1)
				end do
			end do									
		case default
			buf1(1:nd1)=matmul(element(ienum).d,dqf(1:nd1))
			buf2(1:nd1)=matmul(dyf(1:nd1),element(ienum).d)
			de(1:nd1,1:nd1)=element(ienum).d-csproduct(buf1(1:nd1), &
				buf2(1:nd1))/dot_product(buf2(1:nd1),dqf(1:nd1))
	end select
	
	call formRmat(Rmat,stressb,invb,m,element(ienum).mat,nd1,lamda, &
		lamda2,lamda3,element(ienum).d,dqwi,dqwi2,dqwi3,ayf,e,DPMWS, &
		DPMWEV)
	de(1:nd1,1:nd1)=matmul(de(1:nd1,1:nd1),Rmat(1:nd1,1:nd1))
	!write(20,'(4e15.7)') (de(i1,1:4),i1=1,4)
	!write(20,'(a)') ''	
	
end subroutine

subroutine MC_UPDATE(ayf,nd1,ienum,vyf,ev,dlamda,dlamda2,dlamda3,lamda,lamda2,&
	lamda3,dqf,dqf2,dqf3,dyf,dyf2,dyf3,dqwi,dqwi2,dqwi3,dsigma,inv,m,Rmat,de,c)
	use solverds
	implicit none
	integer::i1,j1
	integer::nd1,ienum,ayf,ayf1
	real(8)::vyf,ev,dlamda,dlamda2,dlamda3,lamda,lamda2,lamda3,dqf(7),dqf2(7), &
		dqf3(7),dyf(7),dyf2(7),dyf3(7),dqwi(3),dqwi2(3),dqwi3(3),Rmat(7,7), &
		de(6,6),dsigma(7),inv(3),m(6,3),c(3,3)
	real(8)::vyf2=0,vyf3=0,dywi(3)=0
	
	call active_yf(vyf,inv,element(ienum).mat,ayf)
	select case(ayf)
		case(0)
			dlamda=vyf/dot_product(matmul(dyf(1:nd1),de(1:nd1,1:nd1)),&
				dqf(1:nd1))
			dsigma(1:nd1)=-dlamda*matmul(de(1:nd1,1:nd1),dqf(1:nd1))
		case(12,23) !return to edge
			call yieldfun(vyf2,element(ienum).mat,inv,ayf,ev)
			call deriv_qf_with_inv(dqwi2,inv,element(ienum).mat,ayf,ev)
			call deriv_yf(dqf2(1:6),dqwi2,m)
			call deriv_yf_with_inv(dywi,inv,element(ienum).mat,ayf,ev)
			call deriv_yf(dyf2(1:6),dywi,m)								
			call TVR_Coef(c,de(1:nd1,1:nd1),nd1,dyf(1:6),dqf(1:6), &
				dyf2(1:6),dqf2(1:6),dyf3(1:6),dqf3(1:6),ayf)
			
			dlamda=c(1,1)*vyf+c(1,2)*vyf2
			dlamda2=c(2,1)*vyf+c(2,2)*vyf2
			dsigma(1:nd1)=-dlamda*matmul(de(1:nd1,1:nd1),dqf(1:nd1))- &
				dlamda2*matmul(de(1:nd1,1:nd1),dqf2(1:nd1))
			lamda2=lamda2+dlamda2

		case(123) !return to apex
			ayf1=12
			call yieldfun(vyf2,element(ienum).mat,inv,ayf1,ev)
			call deriv_qf_with_inv(dqwi2,inv,element(ienum).mat,ayf1,ev)
			call deriv_yf(dqf2(1:6),dqwi2,m)
			call deriv_yf_with_inv(dywi,inv,element(ienum).mat,ayf,ev)
			call deriv_yf(dyf2(1:6),dywi,m)								
			ayf1=23
			call yieldfun(vyf3,element(ienum).mat,inv,ayf1,ev)
			call deriv_qf_with_inv(dqwi3,inv,element(ienum).mat,ayf1,ev)
			call deriv_yf(dqf3(1:6),dqwi3,m)
			call deriv_yf_with_inv(dywi,inv,element(ienum).mat,ayf,ev)
			call deriv_yf(dyf3(1:6),dywi,m)								
			call TVR_Coef(c,de(1:nd1,1:nd1),nd1,dyf(1:6),dqf(1:6),dyf2(1:6), &
				dqf2(1:6),dyf3(1:6),dqf3(1:6),ayf)
			dlamda=c(1,1)*vyf+c(1,2)*vyf2+c(1,3)*vyf3
			dlamda2=c(2,1)*vyf+c(2,2)*vyf2+c(2,3)*vyf3
			dlamda3=c(3,1)*vyf+c(3,2)*vyf2+c(3,3)*vyf3
			dsigma(1:nd1)=-dlamda*matmul(de(1:nd1,1:nd1),dqf(1:nd1))- &
				dlamda2*matmul(de(1:nd1,1:nd1),dqf2(1:nd1))- &
				dlamda3*matmul(de(1:nd1,1:nd1),dqf3(1:nd1))
			lamda2=lamda2+dlamda2	
			lamda3=lamda3+dlamda3
	end select

end subroutine

!yield the zero, first and second derivative of the A function implemented in MC criteria.
!REf to book written by Crisfield. vol 2, p102,equ (14.10) and p119 equ (14.19)
subroutine MC_A(mcA,inv,mat,activeyf,isqf)
	use solverds
	implicit none
	logical::isqf
	integer::activeyf,mat
	real(8)::McA(3),inv(3),phi,sq3,vcos,vsin,vsin2
	
	if(isqf) then
		phi=material(mat).property(5)/180.0*pi()
	else
		phi=material(mat).property(4)/180.0*pi()
	end if
	sq3=3**0.5
	vcos=dcos(inv(3))
	vsin=dsin(inv(3))
	vsin2=dsin(phi)

	select case(activeyf)
		case(12) !f3b
			McA(1)=0.5*vcos*(1-vsin2)+vsin*(3+vsin2)/sq3/2
			McA(2)=-0.5*vsin*(1-vsin2)+vcos*(3+vsin2)/sq3/2
		case(23) !f2b
			McA(1)=0.5*vcos*(1+vsin2)+vsin*(-3+vsin2)/sq3/2
			McA(2)=-0.5*vsin*(1+vsin2)+vcos*(-3+vsin2)/sq3/2		
		case default
			McA(1)=vcos-vsin*vsin2/sq3
			McA(2)=-vsin-vcos*vsin2/sq3		
	end select
	
	McA(3)=-McA(1)
end subroutine

!judge which yield function is active
!ayf=0, one vector return
!ayf=12, two vector return, yield funciton s1>=s2>s3 and s2>=s1>=s3 are activated.
!ayf=23, two vector return, yield funciton s1>=s2>s3 and s21>=s3>=s2 are activated.
!afy=123, return to apex.
subroutine active_yf(vyf,inv,mat,ayf)
	use solverds
	implicit none
	integer::mat,ayf
	real(8)::vyf,inv(3),u12,u23,v,J2,s,q,c1
	

	v=material(mat).property(2)
	J2=inv(2)**2/3.0
	s=dsin(material(mat).property(4)/180*pi())
	q=dsin(material(mat).property(5)/180*pi())
	c1=2*(1-2*v+s*q)/((1-2*v)*(q+1))*J2**0.5
	
	u12=vyf-c1*dsin(pi()/6.0-inv(3))
	u23=vyf-c1*dsin(pi()/6.0+inv(3))
	
	ayf=123
	if(u12<0.and.u23<0) ayf=0
	if(u12>0.and.u23<0) ayf=12
	if(u12<0.and.u23>0) ayf=23

	if(dsin(inv(3))>=0.49) ayf=0
		
	
end subroutine

!calculate the coefficients matrix C(3,3) for two vector return method 
!implemented in MC model, return its inverse. 
subroutine TVR_Coef(c,de,nd,dyf1,dqf1,dyf2,dqf2,dyf3,dqf3,ayf)
	use solverds
	!use operation_i
	implicit none
	integer::nd,ayf,afy1
	real(8)::c(3,3),de(nd,nd),dyf1(6),dqf1(6),dyf2(6),dqf2(6), &
		m2(nd,1),dyf3(6),dqf3(6),c1(3,3)=0,t1=0
	
	!c=0
	!select case(ayf)
	!	case(12,23)
	!		m2(:,1)=dyf1(1:nd)
	!		c(1,1)=dot_product(matmul(m2(:,1),de),dqf1(1:nd))
	!		!c(1,1)=dot_product(m2(:,1),matmul(de,dqf1(1:nd)))
	!		c(1,2)=dot_product(matmul(m2(:,1),de),dqf2(1:nd))
	!		m2(:,1)=dyf2(1:nd)
	!		c(2,1)=dot_product(matmul(m2(:,1),de),dqf1(1:nd))
	!		c(2,2)=dot_product(matmul(m2(:,1),de),dqf2(1:nd))
	!		t1=c(1,1)*c(2,2)-c(1,2)*c(2,1)
	!		c1=c
	!		c(1,1)=c1(2,2)
	!		c(2,2)=c1(1,1)
	!		c(1,2)=-c1(1,2)
	!		c(2,1)=-c1(2,1)
	!		c=c/t1
	!	case(123)
	!		m2(:,1)=dyf1(1:nd)
	!		c(1,1)=dot_product(matmul(m2(:,1),de),dqf1(1:nd))
	!		c(1,2)=dot_product(matmul(m2(:,1),de),dqf2(1:nd))
	!		c(1,3)=dot_product(matmul(m2(:,1),de),dqf3(1:nd))			
	!		m2(:,1)=dyf2(1:nd)
	!		c(2,1)=dot_product(matmul(m2(:,1),de),dqf1(1:nd))
	!		c(2,2)=dot_product(matmul(m2(:,1),de),dqf2(1:nd))
	!		c(2,3)=dot_product(matmul(m2(:,1),de),dqf3(1:nd))
	!		m2(:,1)=dyf3(1:nd)
	!		c(3,1)=dot_product(matmul(m2(:,1),de),dqf1(1:nd))
	!		c(3,2)=dot_product(matmul(m2(:,1),de),dqf2(1:nd))
	!		c(3,3)=dot_product(matmul(m2(:,1),de),dqf3(1:nd))
	!		c(1:3,1:3)=.i.c(1:3,1:3)
	!end select

	
end subroutine

!according to the stress vector(stress)calculate the second derivertives of the p,q,theta 
!stress(1-6)=sx,sy,sz,sxy,syz,sxz.
subroutine second_deriv_sinv(secM,stress,inv)
	implicit none
	integer::i,j,p,q,m,n
	integer,external::delta
	real(8)::secm(6,6,3),stress(6),inv(3),devs(3,3)=0,t(3,3)=0,J2=0.0, &
		coef(6)=0,vsin=0,vcos=0,SecM_tensor(3,3,3,3,3)=0, &
		Wpqmn(3,3,3,3)=0,Ppqmn(3,3,3,3)=0
	
		
	!devs=stress
	devs(1,1)=2.0/3.0*stress(1)-1.0/3.0*(stress(2)+stress(3))
	devs(2,2)=2.0/3.0*stress(2)-1.0/3.0*(stress(1)+stress(3))
	devs(3,3)=2.0/3.0*stress(3)-1.0/3.0*(stress(2)+stress(1))
	devs(1,2)=stress(4)
	devs(2,1)=stress(4)
	devs(2,3)=stress(5)
	devs(3,2)=stress(5)
	devs(1,3)=stress(6)
	devs(3,1)=stress(6)
	
	vsin=dcos(3*inv(3)) !attention,theta is different with the original ones.
	vcos=-dsin(3*inv(3))
	
	coef(1)=-4.5/inv(2)**4*(vcos/vsin+1.5*vcos/vsin**3)
	coef(2)=81.0/4.0/inv(2)**5/vsin**3
	coef(3)=81.0/4.0/inv(2)**5*(1.0/vsin+vcos**2/vsin**3)
	coef(4)=-243.0/4.0/inv(2)**6/vsin**3*vcos
	coef(5)=1.5/inv(2)**2/vsin*vcos
	coef(6)=-4.5/inv(2)**3/vsin
	
	t=0
	J2=1/2.0*(devs(1,1)**2+devs(2,2)**2+devs(3,3)**2)+devs(1,2)**2+devs(2,3)**2+devs(1,3)**2
	do i=1,3
		do j=1,3
			t(i,j)=dot_product(devs(i,:),devs(:,j))-2.0/3.0*delta(i,j)*J2
		end do
	end do
	
	!yield the tensor form of the derivative of the  invariants with respective 
	!to stress.
	do p=1,3
		do q=1,3
			do m=1,3
				do n=1,3
					secM_tensor(p,q,m,n,2)=1.5/inv(2)*(delta(p,m)*delta(n,q)-1.0/3.0* &
						delta(p,q)*delta(n,m))-2.25/inv(2)**3*devs(m,n)*devs(p,q)
					Wpqmn(p,q,m,n)=devs(n,p)*delta(q,m)+devs(q,m)*delta(n,p)-2.0/3.0* &
						(devs(q,p)*delta(n,m)+delta(p,q)*devs(m,n))
					Ppqmn(p,q,m,n)=delta(m,p)*delta(n,q)-1.0/3.0*delta(p,q)*delta(m,n)
					SecM_tensor(p,q,m,n,3)=coef(1)*devs(p,q)*devs(m,n)+coef(2)* &
						devs(p,q)*t(m,n)+coef(3)*t(p,q)*devs(m,n)+coef(4)*t(p,q)*t(m,n) &
						+coef(5)*Ppqmn(p,q,m,n)+coef(6)*Wpqmn(p,q,m,n)
				end do
			end do
		end do
	end do
	
	!vecterization. followed with Peter Helnwein
	secM=0.0
	do i=1,6
		call voigt(i,p,q)
		do j=1,6
			call voigt(j,m,n)
			if(i<=3.and.j<=3)then
				secM(i,j,2:3)=secM_tensor(p,q,m,n,2:3)
			else
!				if(i>3.and.j>3) then
					secM(i,j,2:3)=2.0*secM_tensor(p,q,m,n,2:3)
!				else
!					secM(i,j,2:3)=2**0.5*secM_tensor(p,q,m,n,2:3)
!				end if
			end if
		end do
	end do

end subroutine

!define the voigt rule
!and 
subroutine voigt(i,m,n)
	implicit none
	integer::i,m,n,dmn
	
	m=0
	n=0
	dmn=0
	select case(i)
		case(1)
			m=1
			n=1
		case(2)
			m=2
			n=2
		case(3)
			m=3
			n=3
		case(4)
			m=1
			n=2			
		case(5)
			m=2
			n=3
		case(6)
			m=3
			n=1
	end select	
end subroutine

!the kronecker delta function
function delta(m,n)
	implicit none
	integer::delta,m,n
	
	delta=0
	if(m==n) delta=1
end function

!according to the stress vector(stress)calculate the derivertives of the p,J2,J3 
!stress(1-6)=sx,sy,sz,sxy,syz,sxz.
subroutine deriv_sinv(m,stress)
	implicit none
	real(8)::m(6,3),stress(6),devs(6),J2=0.0
	
	devs=stress
	devs(1)=2.0/3.0*stress(1)-1.0/3.0*(stress(2)+stress(3))
	devs(2)=2.0/3.0*stress(2)-1.0/3.0*(stress(1)+stress(3))
	devs(3)=2.0/3.0*stress(3)-1.0/3.0*(stress(2)+stress(1))
	
	m=0.0
	m(1:3,1)=1.0/3.0
	m(:,2)=devs
	m(4:6,2)=devs(4:6)*2.0
	J2=1/2.0*(devs(1)**2+devs(2)**2+devs(3)**2)+devs(4)**2+devs(5)**2+devs(6)**2
	m(1,3)=devs(2)*devs(3)-devs(5)**2
	m(2,3)=devs(1)*devs(3)-devs(6)**2
	m(3,3)=devs(2)*devs(1)-devs(4)**2
	m(1:3,3)=m(1:3,3)+J2/3.0
	m(5,3)=2*(devs(4)*devs(6)-devs(1)*devs(5))
	m(6,3)=2*(devs(4)*devs(5)-devs(2)*devs(6))
	m(4,3)=2*(devs(5)*devs(6)-devs(3)*devs(4))
end subroutine


!according to the stress, calculate the stress invarants(p,(3J2)**0.5=q,Lode angle)
subroutine invanrant(stress,inv)
	implicit none
	real(8)::stress(6),inv(3),devs(6),J3=0.0
	
	devs=stress
	devs(1)=2.0/3.0*stress(1)-1.0/3.0*(stress(2)+stress(3))
	devs(2)=2.0/3.0*stress(2)-1.0/3.0*(stress(1)+stress(3))
	devs(3)=2.0/3.0*stress(3)-1.0/3.0*(stress(2)+stress(1))	
	
	inv=0.0
	inv(1)=(stress(1)+stress(2)+stress(3))/3.0
	inv(2)=1.0/2.0*((stress(1)-stress(2))**2+(stress(3)-stress(2))**2+ &
				(stress(1)-stress(3))**2+6*(stress(4)**2+stress(5)**2+stress(6)**2))
	inv(2)=inv(2)**0.5
	
	if(inv(2)>1e-10) then
		J3=devs(1)*devs(2)*devs(3)+2*devs(4)*devs(5)*devs(6)- &
				devs(1)*devs(5)**2-devs(2)*devs(6)**2-devs(3)*devs(4)**2
		inv(3)=-27.0/2.0*J3/inv(2)**3
		if(inv(3)>1) inv(3)=1.0
		if(inv(3)<-1) inv(3)=-1.0
		inv(3)=dasin(inv(3))/3.0
	else
		inv(3)=0.0
	end if
	
end subroutine

!according to the material(mat).type and stress invarants,calculate the value of the yield functions(vyt)
subroutine yieldfun(vyf,mat,inv,ayf,ev)
	use solverds	
	implicit none
	real(8)::vyf,inv(3),t1,t2,mcA(3),t3,ev
	integer::mat,ayf
	logical::isqf=.false.
	
	select case(material(mat).type)
		case(Mises)
			vyf=inv(2)-material(mat).property(3)
		case(MC)
			call Mc_A(mcA,inv,mat,ayf,isqf)
			t1=dsin(material(mat).property(4)/180.0*pi())
			t2=dcos(material(mat).property(4)/180.0*pi())
			vyf=inv(1)*t1+ &
					inv(2)/3**0.5*mcA(1)-material(mat).property(3)*t2
!			t3=(dcos(inv(3))/3**0.5-dsin(inv(3))*t1/3.0)*3**0.5
!			vyf=inv(1)*t1+ &
!					inv(2)*(dcos(inv(3))/3**0.5-dsin(inv(3))*t1/3.0)- &
!					material(mat).property(3)*t2
		case(camclay)
			vyf=inv(2)**2-material(mat).property(1)**2*(inv(1)*(inv(1)-ev))
		case default
			vyf=-1.0
	end select 
end subroutine

!according to derivatives of yield(plastic potential) functions
!with respect to invarants((sigma_m,J2,J3),stored in "dywi" 
!and the derivatives of invarants w.r.t. sigma which stored in "m",
!calculate the derivative of yield functions or potential function w.r.t.
!stress, sigma.
subroutine deriv_yf(dyf,dywi,m)
	use solverds
	implicit none
	real(8)::dyf(6),dywi(3),m(6,3)
	
	dyf=dywi(1)*m(:,1)+dywi(2)*m(:,2)+dywi(3)*m(:,3)
	
end subroutine

!according to material(mat).type and stress invarants, calculate the 
!derivative of yield function with respect to the stress invarants(sigma_m,J2,J3).
subroutine deriv_yf_with_inv(dywi,inv,mat,ayf,ev)
	use solverds
	implicit none
	real(8)::dywi(3),inv(3),c1=1.0,mcA(3),J2=0,ev
	integer::mat,ayf
	logical::isqf=.false.
	
	dywi=0.0
	select case(material(mat).type)
		case(mises)
			dywi(2)=3.0/2.0/inv(2)			
		case(mc)
			dywi(1)=dsin(material(mat).property(4)/180.0*pi())
			J2=inv(2)**2/3.0
			if(abs(dsin(inv(3)))<0.49) then
				call Mc_A(mcA,inv,mat,ayf,isqf)
				dywi(2)=0.5*(mcA(1)-dtan(3*inv(3))*mcA(2))/J2**0.5
				dywi(3)=-3**0.5*McA(2)/(2*J2*dcos(3*inv(3)))
!				dywi(2)=dcos(inv(3))/((4.0/3.0)**0.5*inv(2))* &
!								(1+dtan(inv(3))*dtan(3*inv(3))+dywi(1)/3**0.5* &
!								(dtan(3*inv(3))-dtan(inv(3))))
!				dywi(3)=(3**0.5*dsin(inv(3))+dywi(1)*dcos(inv(3)))/ & 
!								(2.0/3.0*inv(2)**2*dcos(3*inv(3)))
			else
				c1=1.0
				if(dsin(inv(3))<0) c1=-1.0
				dywi(2)=(3.0-c1*dywi(1))*0.25/inv(2)				
			end if
		case(camclay)
			dywi(1)=-material(mat).property(1)**2*(2*inv(1)-ev)
			dywi(2)=3.0
			
	end select	
end subroutine

!according to material(mat).type and stress invarants, calculate the 
!derivative of potential function with respect to the stress invarants(sigma_m,J2,J3).
subroutine deriv_qf_with_inv(dywi,inv,mat,ayf,ev)
	use solverds
	implicit none
	real(8)::dywi(3),inv(3),c1=1.0,J2=0,mcA(3),ev
	integer::mat,ayf
	logical::isqf=.true.
	
	dywi=0.0
	select case(material(mat).type)
		case(mises)
			dywi(2)=3.0/2.0/inv(2)			
		case(mc)
			dywi(1)=dsin(material(mat).property(5)/180.0*pi())
			J2=inv(2)**2/3.0
			if(abs(dsin(inv(3)))<0.49) then
				call Mc_A(mcA,inv,mat,ayf,isqf)
				dywi(2)=0.5*(mcA(1)-dtan(3*inv(3))*mcA(2))/J2**0.5
				dywi(3)=-3**0.5*McA(2)/(2*J2*dcos(3*inv(3)))
!				dywi(2)=dcos(inv(3))/((4.0/3.0)**0.5*inv(2))* &
!								(1+dtan(inv(3))*dtan(3*inv(3))+dywi(1)/3**0.5* &
!								(dtan(3*inv(3))-dtan(inv(3))))
!				dywi(3)=(3**0.5*dsin(inv(3))+dywi(1)*dcos(inv(3)))/ & 
!								(2.0/3.0*inv(2)**2*dcos(3*inv(3)))

			else
				c1=1.0
				if(dsin(inv(3))<0) c1=-1.0
				dywi(2)=(3.0-c1*dywi(1))*0.25/inv(2)				
			end if
		case(camclay)
			dywi(1)=-material(mat).property(1)**2*(2*inv(1)-ev)
			dywi(2)=3.0			
	end select	
end subroutine

!generation plastic body force using viscoplasticity method.
subroutine bload_viscoplasticity(iiter,iscon,bdylds)
	use solverds
	implicit none
	logical::iscon
	integer::i,j,k,n1,iiter,ayf=0
	real(8)::bdylds(ndof),un(50)=0.0,stress1(6)=0.0,strain1(6)=0.0, &
					dee(4,4)=0.0,sigma(6)=0.0,inv(3)=0.0,vyf=0.0, &
					devp(6)=0.0,m(6,3)=0.0,dywi(3)=0.0,dyf(6)=0.0, &
					evp(6),bload(50)=0.0,dt=0.0,ev=0
	integer::nd1=0

		
	do i=1,enum
		n1=element(i).ngp
		nd1=element(i).nd
		strain1=0.0
		stress1=0.0
		bload=0.0
		call timestep_vp(dt,element(i).mat)
		do k=1,element(i).nnum
			un(2*k-1)=load(node(element(i).node(k)).dof(1))	
			un(2*k)=load(node(element(i).node(k)).dof(2))	
		end do
		do j=1,n1
			!strain of this incrementals
			strain1(1:nd1)=matmul(element(i).b(1:nd1,:,j),un(1:element(i).ndof))
			!get the elatic strain
			strain1=strain1-element(i).evp(:,j)
			!stress incremental
			stress1(1:nd1)=matmul(element(i).d,strain1(1:nd1))
			!total stress
			stress1=stress1+element(i).stress(:,j)
			!check whether penetrate the yield surface.
			call invanrant(stress1,inv)
 			call yieldfun(vyf,element(i).mat,inv,ayf,ev)
			
			if(vyf>0) then
				call deriv_sinv(m,stress1)
				call deriv_qf_with_inv(dywi,inv,element(i).mat,ayf,ev)
				call deriv_yf(dyf,dywi,m)
				evp=dt*vyf*dyf
				element(i).evp(:,j)=element(i).evp(:,j)+evp
				devp(1:nd1)=matmul(element(i).d,evp(1:nd1))
				bload(1:element(i).ndof)=bload(1:element(i).ndof)+ &
										matmul(devp(1:nd1),element(i).b(:,:,j))* &
										element(i).detjac(j)*ecp(element(i).et).weight(j)
			end if
			!update the stress and strain
			if(iscon.or.iiter==solver_control.niteration) then
				element(i).stress(:,j)=stress1
				element(i).strain(:,j)=element(i).strain(:,j)+strain1+ &
					element(i).evp(:,j)
				element(i).pstrain(:,j)=element(i).pstrain(:,j)+element(i).evp(:,j)
				bdylds=0.0 !reset				
			end if
		end do
		!get the globe self-equilibrating global body loads
		bdylds(element(i).g)=bdylds(element(i).g)+bload(1:element(i).ndof)
	end do
	
end subroutine

subroutine timestep_vp(dt,mat)
	use solverds
	implicit none
	integer::mat
	real(8)::dt,E,v,vsin
	
	dt=0.0	
	select case(material(mat).type)
		case(mises)
			E=material(mat).property(1)
			v=material(mat).property(2)
			dt=4*(v+1)/(3*E)
		case(MC)
			E=material(mat).property(1)
			v=material(mat).property(2)
			vsin=dsin(material(mat).property(4)/180.0*pi())
			dt=4*(v+1)*(1-2*v)/(E*(1-2*v+vsin**2))		
	end select
end subroutine

!generation plastic body force using initial stress method.
subroutine bload_inistress(iiter,iscon,bdylds)
	use solverds
	!use operation_ix
	implicit none
	logical::iscon
	integer::i,j,k,n1,iiter,ayf=0
	real(8)::bdylds(ndof),un(50)=0.0,stress1(6)=0.0,strain1(6)=0.0, &
					dee(4,4)=0.0,sigma(6)=0.0,inv(3)=0.0,vyf=0.0, &
					devp(6)=0.0,m(6,3)=0.0,dywi(3)=0.0,dyf(6)=0.0, &
					evp(6),bload(50)=0.0,dt=0.0,vyf_old=0.0,fac=0.0, &
					stress2(6)=0.0,dqwi(3)=0.0,dqf(6)=0.0,dp(6,6)=0.0,&
					de(6,6)=0.0,t1=0.0,pstress1(6)=0.0,lamda=0.0,&
					buf1(6)=0.0,buf2(6)=0.0,pstrain(6)=0.0,ev=0

	integer::nd1=0
	
	!bdylds=0.0
	!
	!do i=1,enum
	!	n1=element(i).ngp
	!	nd1=element(i).nd
	!	bload=0.0
	!	do k=1,element(i).nnum
	!		un(2*k-1)=load(node(element(i).node(k)).dof(1))	
	!		un(2*k)=load(node(element(i).node(k)).dof(2))	
	!	end do
	!	do j=1,n1
	!		strain1=0.0
	!		stress1=0.0
	!		pstress1=0.0
	!		!strain of this incrementals
	!		strain1(1:nd1)=matmul(element(i).b(1:nd1,:,j),un(1:element(i).ndof))
	!		!stress incremental
	!		de(1:nd1,1:nd1)=element(i).d
	!		stress2(1:nd1)=matmul(de(1:nd1,1:nd1),strain1(1:nd1))
	!		!total stress
	!		stress1=stress2+element(i).stress(:,j)
	!		!check whether penetrate the yield surface.
	!		call invanrant(stress1,inv)
 !			call yieldfun(vyf,element(i).mat,inv,ayf,ev)
	!		
	!		if(vyf>0) then
	!			call invanrant(element(i).stress(:,j),inv)
	!			call yieldfun(vyf_old,element(i).mat,inv,ayf,ev)
	!			fac=vyf/(vyf-vyf_old)
	!			stress1=stress2*(1-fac)+element(i).stress(:,j)
	!			call invanrant(stress1,inv)
	!			call deriv_sinv(m,stress1)
	!			call deriv_qf_with_inv(dqwi,inv,element(i).mat,ayf,ev)
	!			call deriv_yf_with_inv(dywi,inv,element(i).mat,ayf,ev)
	!			call deriv_yf(dyf,dywi,m)
	!			call deriv_yf(dqf,dqwi,m)
	!			dp=0.0			
	!			buf1(1:nd1)=matmul(dyf(1:nd1),de(1:nd1,1:nd1))
	!			buf2(1:nd1)=matmul(de(1:nd1,1:nd1),dqf(1:nd1))
	!			do k=1,nd1
	!				dp(k,1:nd1)=buf2(k)*buf1(1:nd1)
	!			end do				
	!			t1=dot_product(buf1(1:nd1),dqf(1:nd1))
	!			dp(1:nd1,1:nd1)=dp(1:nd1,1:nd1)/t1
	!			dp=dp*fac
	!			pstress1(1:nd1)=matmul(dp(1:nd1,1:nd1),strain1(1:nd1))
	!			bload(1:element(i).ndof)=bload(1:element(i).ndof)+ &
	!									matmul(pstress1(1:nd1),element(i).b(:,:,j))* &
	!									element(i).detjac(j)*ecp(element(i).et).weight(j)
	!		end if
	!		!update the stress and strain
	!		if(iscon.or.iiter==solver_control.niteration) then
	!			element(i).stress(:,j)=element(i).stress(:,j)+stress2-pstress1
	!			element(i).strain(:,j)=element(i).strain(:,j)+strain1
	!			if(vyf>0) then
	!				element(i).evp(1:nd1,j)=strain1(1:nd1)-(de(1:nd1,1:nd1).ix. &
	!					(stress2(1:nd1)-pstress1(1:nd1)))
	!				element(i).pstrain(:,j)=element(i).pstrain(:,j)+element(i).evp(:,j)
	!			end if				
	!		end if
	!	end do
	!	!get the globe self-equilibrating global body loads
	!	bdylds(element(i).g)=bdylds(element(i).g)+bload(1:element(i).ndof)
	!	
	!end do
	return
end subroutine


subroutine element_activate(istep)
    use SOLVERLIB
    implicit none
    integer,intent(in)::istep
    integer::i,j,k
    
    do concurrent (i=1:enum)
        if(element(i).ec/=soilspring) then
            if(abs(sf(element(i).sf).factor(istep))>1e-7) then
                element(i).isactive=1
            else
                element(i).isactive=0
				!if(allocated(element(i).gforce)) element(i).gforce=0
				!if(allocated(element(i).gforceILS)) element(i).gforceILS=0
            endif
        endif
    enddo
end
