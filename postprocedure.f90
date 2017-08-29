
subroutine outdata(iincs,iiter,iscon,isfirstcall,isubts)	
	use solverds
	implicit none
	integer::iincs,iiter,isubts,file_unit,ieset,file_diagram
	logical::iscon,isfirstcall,isbarfamily,anybarfamily
	integer::i,j,nc,n1,k,k1,iset1
	character(1024)::cstring=''
	real(8)::rar(MNDOF),t1

	!
	if(iincs==0.and.(geostatic.method==cal_geo)) then
		load=0.0
		do j=1,enum
			element(j).strain=0.0
			element(j).pstrain=0.0
			element(j).evp=0.0
		end do
	end if

	if(.not.allocated(tdispInLS)) then
		allocate(tdispInLS(ndof))
		tdispInLS=0.0
	end if
	
!	call GS2LS_Displacement()
	
	!map stress in gauss point to nodes.
	anybarfamily=.false.
	do i=1,enum
		if(element(i).isactive==0) cycle
        
		if(outvar(SFR).value>0.AND.ELEMENT(I).EC==CPE) CALL stree_failure_ratio_cal(I,IINCS)	
		select case(element(i).ec)

			case(stru,SPRING,SOILSPRING)
				!transfor the nodal force in the globle system to the local system. 		
				!call NodalForceInLocalSystem(i)
				anybarfamily=.true.
			case default
				select case(solver_control.i2ncal) 
					case(spr)
						print *, 'CAUSION! SPR METHOD IS OUT OF MANTAINED NOW.SHT METHOD IS USED INSTEAD'
						!call spr_stress_strain_cal()						
					case(sht)
						call shift_stress_strain_cal(i)				
					case(avg)
						call average_stress_strain_cal(i)
					case default
						
						call extrapolation_stress_strain_cal(i)
				end select
						
		end select  

	end do
	
	call E2N_stress_strain(IINCS)
	
	CALL BC_RHS_OUT(iincs,iiter,ISUBTS)
	
	if(anybarfamily) call barfamily_outfordiagram(iincs,isubts)
	
	file_unit=10
	file_diagram=20
	if(isfirstcall) then 
		open(unit=file_unit,file=resultfile2,status='replace')
		!if(anybarfamily) open(unit=file_diagram,file=resultfile21,status='replace')		
		cstring='TITLE = '//'"'//trim(title)//'"'
		write(file_unit,'(a1024)') cstring
		call tecplot_variables(cstring)
		write(file_unit,'(a1024)') cstring
		
		
		call tecplot_zonetitle(iincs,iiter,isfirstcall,isubts)
		
		do i=1,neset
            iset1=esetid(i)
			write(file_unit,'(a1024)') eset(iset1).zonetitle
			!if(anybarfamily) write(file_diagram,'(a1024)') eset(iset1).zonetitle
			
			isbarfamily=eset(iset1).et==bar.or.eset(iset1).et==bar2d.or.eset(iset1).et==beam.or.eset(iset1).et==beam2d.or.eset(iset1).et==ssp2d
			
			if(isbarfamily) then
				if(solver_control.datapaking)then
					call pointout_barfamily(file_unit,iset1)
					!call pointout_barfamily_diagram(file_diagram,i,iincs,isubts)
                else
                    print *, "To be improved. SUB=outdata"    
					!call blockout_barfamily(file_unit,ieset)
				end if				
			else
				if(i==1) then
					if(solver_control.datapaking)then
						call pointout(file_unit,IINCS,ISUBTS,IITER)					  
					else
						call blockout(file_unit)
					end if
				end if
			end if
			do j=eset(iset1).enums,eset(iset1).enume				
				select case(eset(iset1).et)
					case(cpe8,cps8,CAX8,cpe8r,cps8r,CAX8R, &
							 cpe8_spg,cps8_spg,CAX8_spg,cpe8r_spg,cps8r_spg,CAX8R_spg, &
							 cpe8_cpl,cps8_cpl,CAX8_cpl,cpe8r_cpl,cps8r_cpl,CAX8R_cpl) 
						nc=4
						write(file_unit,9999) element(j).node(8),element(j).node(1),&
											element(j).node(5),element(j).node(5)
						write(file_unit,9999) element(j).node(5),element(j).node(2),&
											element(j).node(6),element(j).node(6)
						write(file_unit,9999) element(j).node(6),element(j).node(3),&
											element(j).node(7),element(j).node(7)
						write(file_unit,9999) element(j).node(7),element(j).node(4),&
											element(j).node(8),element(j).node(8)
						write(file_unit,9999) element(j).node(5),element(j).node(6),&
											element(j).node(7),element(j).node(8)								
						case(cpe6,cps6,CAX6,&
								 cpe6_spg,cps6_spg,CAX6_spg, &
								 cpe6_cpl,cps6_cpl,cax6_cpl)
							nc=3
							write(file_unit,9999) element(j).node(1),element(j).node(4),&
												element(j).node(6)
							write(file_unit,9999) element(j).node(2),element(j).node(5),&
												element(j).node(4)
							write(file_unit,9999) element(j).node(3),element(j).node(6),&
												element(j).node(5)
							write(file_unit,9999) element(j).node(4),element(j).node(5),&
												element(j).node(6)
						case(cpe15,cps15,CAX15, &
								 cpe15_spg,cps15_spg,CAX15_spg, &
								 cpe15_cpl,cps15_cpl,CAX15_cpl)
							nc=3
							write(file_unit,9999) element(j).node(3),element(j).node(11),element(j).node(10)
							write(file_unit,9999) element(j).node(11),element(j).node(6),element(j).node(14)
							write(file_unit,9999) element(j).node(11),element(j).node(14),element(j).node(10)
							write(file_unit,9999) element(j).node(10),element(j).node(14),element(j).node(5)
							write(file_unit,9999) element(j).node(6),element(j).node(12),element(j).node(15)
							write(file_unit,9999) element(j).node(6),element(j).node(15),element(j).node(14)
							write(file_unit,9999) element(j).node(14),element(j).node(15),element(j).node(13)
							write(file_unit,9999) element(j).node(14),element(j).node(13),element(j).node(5)
							write(file_unit,9999) element(j).node(5),element(j).node(13),element(j).node(9)
							write(file_unit,9999) element(j).node(12),element(j).node(1),element(j).node(7)
							write(file_unit,9999) element(j).node(12),element(j).node(7),element(j).node(15)
							write(file_unit,9999) element(j).node(15),element(j).node(7),element(j).node(4)
							write(file_unit,9999) element(j).node(15),element(j).node(4),element(j).node(13)
							write(file_unit,9999) element(j).node(13),element(j).node(4),element(j).node(8)
							write(file_unit,9999) element(j).node(13),element(j).node(8),element(j).node(9)
							write(file_unit,9999) element(j).node(9),element(j).node(8),element(j).node(2)
					case(prm6,prm6_spg,prm6_cpl)
							nc=8	
							write(file_unit,9999) element(j).node(1:3),element(j).node(3), &
																		element(j).node(4:6),element(j).node(6)
					case(prm15,prm15_spg,prm15_cpl) !desolve into 14 tetrahedral element 
							nc=4
							write(file_unit,9999) element(j).node(1),element(j).node(7),element(j).node(9),element(j).node(13)
							write(file_unit,9999) element(j).node(7),element(j).node(2),element(j).node(8),element(j).node(14)
							write(file_unit,9999) element(j).node(8),element(j).node(3),element(j).node(9),element(j).node(15)
							write(file_unit,9999) element(j).node(8),element(j).node(9),element(j).node(13),element(j).node(15)
							write(file_unit,9999) element(j).node(9),element(j).node(8),element(j).node(13),element(j).node(7)
							write(file_unit,9999) element(j).node(8),element(j).node(14),element(j).node(13),element(j).node(7)
							write(file_unit,9999) element(j).node(8),element(j).node(13),element(j).node(14),element(j).node(15)
							write(file_unit,9999) element(j).node(4),element(j).node(10),element(j).node(12),element(j).node(13)
							write(file_unit,9999) element(j).node(10),element(j).node(5),element(j).node(11),element(j).node(14)
							write(file_unit,9999) element(j).node(11),element(j).node(6),element(j).node(12),element(j).node(15)
							write(file_unit,9999) element(j).node(12),element(j).node(11),element(j).node(13),element(j).node(15)
							write(file_unit,9999) element(j).node(11),element(j).node(12),element(j).node(13),element(j).node(10)
							write(file_unit,9999) element(j).node(11),element(j).node(13),element(j).node(14),element(j).node(10)
							write(file_unit,9999) element(j).node(11),element(j).node(14),element(j).node(13),element(j).node(15)
						case(tet10,tet10_spg,tet10_cpl)
							nc=4
							write(file_unit,9999) element(j).node(8),element(j).node(7),element(j).node(5),element(j).node(1)
							write(file_unit,9999) element(j).node(5),element(j).node(6),element(j).node(9),element(j).node(2)
							write(file_unit,9999) element(j).node(6),element(j).node(7),element(j).node(10),element(j).node(3)
							write(file_unit,9999) element(j).node(8),element(j).node(9),element(j).node(10),element(j).node(4)
							write(file_unit,9999) element(j).node(6),element(j).node(7),element(j).node(8),element(j).node(10)
							write(file_unit,9999) element(j).node(6),element(j).node(8),element(j).node(9),element(j).node(10)
							write(file_unit,9999) element(j).node(8),element(j).node(7),element(j).node(6),element(j).node(5)
							write(file_unit,9999) element(j).node(9),element(j).node(8),element(j).node(6),element(j).node(5)
						case(bar,bar2d,beam,beam2d,ssp2d)
							nc=8
							write(file_unit,9999) (((element(j).node2(k)-1)*4+k1,k1=1,4),k=1,element(j).nnum)
						case default
							nc=element(j).nnum
							write(file_unit,9999) element(j).node(1:nc)
				end select
				
			end do			
			
		end do
			
		close(file_unit)

!		open(unit=1,file=resultfile,status='replace')
		isfirstcall=.false.

	else
		open(unit=file_unit,file=resultfile2,status='old',access='append')
		call tecplot_zonetitle(iincs,iiter,isfirstcall,isubts)
		do i=1,neset
            iset1=esetid(i)
			write(file_unit,'(a1024)') eset(iset1).zonetitle
			isbarfamily=eset(iset1).et==bar.or.eset(iset1).et==bar2d.or.eset(iset1).et==beam.or.eset(iset1).et==beam2d.or.eset(iset1).et==ssp2d
			if(isbarfamily) then
				if(solver_control.datapaking)then
					call pointout_barfamily(file_unit,iset1)					  
                else
                    print *, "To be improved. SUB=outdata"    
					!call blockout_barfamily(file_unit,ieset)
				end if
			else
				if(i==1) then
					if(solver_control.datapaking)then
						call  pointout(file_unit,IINCS,ISUBTS,IITER)				  
					else
						call blockout(file_unit)
					end if
				end if
			end if
			
		end do
	end if 
    close(file_unit)
	999 format(<nc>E15.7)
	9999 format(<nc>I7)
end subroutine

SUBROUTINE NODAL_ACTIVE_DOF(INODE,DOFS,NDOFS)
	USE SOLVERDS
	IMPLICIT NONE
	INTEGER,INTENT(IN)::INODE,NDOFS
	INTEGER,INTENT(OUT)::DOFS(NDOFS)
	INTEGER::N1,I,IDOF1
	
	N1=0
	DOFS=-1
	DO I=1,MNDOF
		IDOF1=NODE(INODE).DOF(I)
		IF(IDOF1>0) THEN
			N1=N1+1
			DOFS(N1)=IDOF1
			IF(N1==NODE(INODE).NDOF) EXIT
		ENDIF
	ENDDO	
	
ENDSUBROUTINE

subroutine BC_RHS_OUT(inc,iter,ISUBTS) !输出节点荷载（力、流量）
	use solverds
	implicit none
	INTEGER,INTENT(IN)::INC,ITER,ISUBTS
	integer::i,J,NITEM1,N1,DOFS1(MNDOF)
	real(KIND=DPN)::t1,t2,t3,dt1,VAL1(MNDOF)
	
	WRITE(99,20)
	WRITE(99,10)
	if(stepinfo(INC).issteady) then
		dt1=1.D0
	else
		dt1=timestep(INC).subts(isubts)
	end if

	DO I=1, BL_NUM
		if(sf(bc_load(i).sf).factor(inc)==-999.D0) cycle
		NITEM1=NODE(bc_load(i).node).NDOF
		CALL NODAL_ACTIVE_DOF(bc_load(i).node,DOFS1,MNDOF)		
		WRITE(99,11) bc_load(i).node,bc_load(i).dof,'LOAD',INC,ITER,NI_NodalForce(DOFS1(1:NITEM1))*dt1,TDISP(DOFS1(1:NITEM1))*dt1
		!QT=QT+NODE(bc_load(i).node).Q*dt1
	END DO
	DO I=1, BD_NUM
		if(bc_disp(I).isdead==1) cycle
		NITEM1=NODE(bc_DISP(i).node).NDOF
		CALL NODAL_ACTIVE_DOF(bc_DISP(i).node,DOFS1,MNDOF)		
		WRITE(99,11) bc_DISP(i).node,bc_DISP(i).dof,'B.C.',INC,ITER,NI_NodalForce(DOFS1(1:NITEM1))*dt1,TDISP(DOFS1(1:NITEM1))*dt1
		!QT=QT+NODE(bc_DISP(i).node).Q*dt1
	END DO
	DO I=1, NUMNSEEP
		if(sf(NSEEP(i).sf).factor(inc)==-999.D0) cycle
		IF(NSEEP(I).ISDEAD==1) CYCLE
		NITEM1=NODE(NSEEP(i).node).NDOF
		CALL NODAL_ACTIVE_DOF(NSEEP(i).node,DOFS1,MNDOF)			
		WRITE(99,11) NSEEP(i).node,NSEEP(i).dof,'S.F.',INC,ITER,NI_NodalForce(DOFS1(1:NITEM1))*dt1,TDISP(DOFS1(1:NITEM1))*dt1
		!QT=QT+NODE(NSEEP(i).node).Q*dt1
	END DO
	
	!t1=0.0d0
	!do i=1,enum
	!    t1=t1+(sum(element(i).sita_fin)-sum(element(i).sita_ini))/element(i).ngp*element(i).detjac(1)*sum(ecp(element(i).et).weight)
	!end do	
	!WRITE(99,12) T1,Qt,T1/Qt
	
	WRITE(99,21)

10  format(3X,"NODE",5X,"DOF",3X,"TYPE",4X,"INC",3X,"ITER",11X,"GENERALIZED_LOADS",11X,"GENERALIZED_DISPLACEMENTS")
11	format(I7,1X,I7,3X,A4,4X,I3,4X,I4,1X,<NITEM1>F15.7,<NITEM1>F15.7)
12	format("WATER STORED:",F15.7,X,"WATER FLOWED IN:",F15.7,X,"RATIO:",F15.7)
20	format("\N******************OUTPUT THE FLUX ON BOUNDARIES AND MASS CONVERSATION RATIO******************"C)
21	format("******************END THE OUTPUT******************\N"C)

end subroutine


subroutine LIMITOUT(IINCS,IITER,ISCON,ISFIRSTCALL,C)
	USE SOLVERDS
	IMPLICIT NONE
	INTEGER::I,J,N2,N3
	integer::iincs,iiter,isubts
	logical::iscon,isfirstcall
	REAL(8)::C(*)
	!calculate the dissipation work for each element
	do i=1,enum
		if(element(i).isactive==0) cycle
		select case(element(i).et)			
			case(UB3)
				n2=element(i).ndof-int(material(element(i).mat).GET(1,IINCS))+1						
			case(UBZT4)
				n2=element(i).ndof-4+1
		end select 
		n3=element(i).ndof
		element(i).property(2)=0.0D0
		do j=n2,n3
			element(i).property(2)=element(i).property(2)+c(element(i).g(j))*load(element(i).g(j))
		end do
		!cal the shear strain 

	end do

	call outdata(iincs,iiter,iscon,isfirstcall,isubts)

end subroutine


subroutine tecplot_variables(cstring)
!*************************************************************************
!function: write the VARIABLES control line used by tecplot. 
!input parameter: NONE
!use module: solverds(outvar)
!**************************************************************************
!
	use solverds
	implicit none
	integer::i,j,k,n1
	character(1024)::cstring
	
	
	n1=ubound(outvar,1)
	nvo=0
	vo=0
	!coordinates X and Y are always output.
	outvar(locx).name='X'
	outvar(locx).value=locx
	outvar(locy).name='Y'
	outvar(locy).value=locy
	if(ndimension==3) then
		outvar(locz).name='Z'
		outvar(locz).value=locz		
	end if
	cstring='Variables='

	!default output for upper bound analysis
	if(solver_control.solver==la04.or.solver_control.solver==lpsolver) then
		outvar(disx).name='Vx'
		outvar(disx).value=disx
		outvar(disy).name='Vy'
		outvar(disy).value=disy
		outvar(pw).name='pw'
		outvar(pw).value=pw	
		return
	end if
	

	do i=1,n1
		if(outvar(i).value==0) cycle
		nvo=nvo+1
		vo(nvo)=i
		outvar(i).ivo=nvo
		nvo=nvo+outvar(i).nval-1 !一个变量包含多个数的情况
		cstring=trim(adjustL(cstring))//'"'//trim(adjustL(outvar(i).name))//'",'
	end do

end subroutine

!when the datapacking format is activated.the quantivities output are
!located at element nodes. 
!This subroutine is used to output these nodal quantivites.
subroutine pointout(FILE_UNIT,ISTEP,ISUBTS,ITER)
	use solverds
	implicit none
	INTEGER,INTENT(IN)::FILE_UNIT,ISTEP,ISUBTS,ITER
	LOGICAL,SAVE::ISFIRSTCALL1	
	integer::i,j,k,idof,IDISQ1=0,NVO1=0
	
	REAL(8)::SUMQ1=0
	
	if(.NOT.allocated(NODALQ)) allocate(NodalQ(nnum,nvo))
    if(.NOT.allocated(VEC)) allocate(VEC(3,NNUM))
	i=1
	do while(i<=nvo)
		select case(vo(i))
			case(locx)
				NodalQ(:,i)=node.coord(1)
			case(locy)
				NodalQ(:,i)=node.coord(2)
			case(locz)
				NodalQ(:,i)=node.coord(3)
			case(disx)
				idof=1
				if(outvar(disx).system/=0) then
					call GS2LS_Displacement(outvar(disx).system,idof,nodalQ(:,i),nnum)
					!NodalQ(:,i)=tdispINLS(node.dof(idof))
				else
					NodalQ(:,i)=tdisp(node.dof(idof))
				end if
			case(disy)
				idof=2
				if(outvar(disy).system/=0) then					
					call GS2LS_Displacement(outvar(disy).system,idof,nodalQ(:,i),nnum)
				else
					NodalQ(:,i)=tdisp(node.dof(idof))
				end if
			case(disz)
				idof=3
				if(outvar(disz).system/=0) then					
					call GS2LS_Displacement(outvar(disz).system,idof,nodalQ(:,i),nnum)
				else
					NodalQ(:,i)=tdisp(node.dof(idof))
				end if
			case(head)
				idof=4
				NodalQ(:,i)=tdisp(node.dof(idof))
			case(Phead)
				idof=4
				if(ndimension==2) NodalQ(:,i)=tdisp(node.dof(idof))-node.coord(2)
				if(ndimension==3) NodalQ(:,i)=tdisp(node.dof(idof))-node.coord(3)												
			case(Rx)
				idof=5
				if(outvar(Rx).system/=0) then					
					call GS2LS_Displacement(outvar(rx).system,idof,nodalQ(:,i),nnum)
				else
					NodalQ(:,i)=tdisp(node.dof(idof))
				end if
			case(Ry)
				idof=6
				if(outvar(Ry).system/=0) then					
					call GS2LS_Displacement(outvar(ry).system,idof,nodalQ(:,i),nnum)
				else
					NodalQ(:,i)=tdisp(node.dof(idof))
				end if
			case(Rz)
				idof=7
				if(outvar(Rz).system/=0) then					
					call GS2LS_Displacement(outvar(rz).system,idof,nodalQ(:,i),nnum)
				else
					NodalQ(:,i)=tdisp(node.dof(idof))
				end if
								
			case(sxx)
				do j=1,nnum
					NodalQ(j,i)=node(j).stress(1)
				end do
			case(syy)
				do j=1,nnum
					nodalQ(j,i)=node(j).stress(2)
				end do  
			case(szz)
				do j=1,nnum
					NodalQ(:,i)=node(j).stress(3)
				end do
			case(sxy)
				do j=1,nnum
					nodalQ(j,i)=node(j).stress(4)	
				end do
			case(syz)
				do j=1,nnum
					nodalQ(j,i)=node(j).stress(5)
				end do
			case(szx)
				do j=1,nnum
					nodalQ(j,i)=node(j).stress(6)
				end do
			case(exx)
				do j=1,nnum
					nodalQ(j,i)=node(j).strain(1)
				end do
			case(eyy)
				do j=1,nnum
					nodalQ(j,i)=node(j).strain(2)
				end do
			case(ezz)
				do j=1,nnum
					nodalQ(j,i)=node(j).strain(3)
				end do
			case(exy)
				do j=1,nnum
					nodalQ(j,i)=node(j).strain(4)
				end do
			case(eyz)
				do j=1,nnum
					nodalQ(j,i)=node(j).strain(5)
				end do
			case(ezx)
				do j=1,nnum
					nodalQ(j,i)=node(j).strain(6)
				end do
			case(pexx)
				do j=1,nnum
					nodalQ(j,i)=node(j).pstrain(1)
				end do
			case(peyy)
				do j=1,nnum	
					nodalQ(j,i)=node(j).pstrain(2)
				end do
			case(pezz)
				do j=1,nnum	
					nodalQ(j,i)=node(j).pstrain(3)
				end do	
			case(pexy)
				do j=1,nnum	
					nodalQ(j,i)=node(j).pstrain(4)
				end do
			case(peyz)
				do j=1,nnum
					nodalQ(j,i)=node(j).pstrain(5)
				end do
			case(pezx)
				do j=1,nnum
					nodalQ(j,i)=node(j).pstrain(6)
				end do
			case(sigma_mises)
				NodalQ(:,i)=node.mises
			case(eeq)
				NodalQ(:,i)=node.eeq
			case(peeq)
				NodalQ(:,i)=node.peeq
			case(xf_out)
				NodalQ(:,i)=node.coord(1)+tdisp(node.dof(1))
			case(yf_out)
				NodalQ(:,i)=node.coord(2)+tdisp(node.dof(2))
			case(zf_out)
				NodalQ(:,i)=node.coord(3)+tdisp(node.dof(3))
			case(gradx)
				do j=1,nnum
                    if(allocated(node(j).igrad)) then
					    NodalQ(j,i)=node(j).igrad(1)
                    else
                        NodalQ(j,i)=0.0d0
                    endif
				end do
			case(grady)
				do j=1,nnum
                     if(allocated(node(j).igrad)) then
					    NodalQ(j,i)=node(j).igrad(2)
                    else
                        NodalQ(j,i)=0.0d0
                    endif
					!NodalQ(j,i)=node(j).igrad(2)
				end do
			case(gradz)
				do j=1,nnum
                     if(allocated(node(j).igrad)) then
					    NodalQ(j,i)=node(j).igrad(3)
                    else
                        NodalQ(j,i)=0.0d0
                    endif
					!NodalQ(j,i)=node(j).igrad(3)
				end do
			case(vx)
				do j=1,nnum
                    if(allocated(node(j).velocity)) then
					    NodalQ(j,i)=node(j).velocity(1)
                    else
                        NodalQ(j,i)=0.0d0
                    endif
					!NodalQ(j,i)=node(j).velocity(1)
				end do
			case(vy)
				do j=1,nnum
                    if(allocated(node(j).velocity)) then
					    NodalQ(j,i)=node(j).velocity(2)
                    else
                        NodalQ(j,i)=0.0d0
                    endif
					!NodalQ(j,i)=node(j).velocity(2)
				end do
			case(vz)
				do j=1,nnum
                    if(allocated(node(j).velocity)) then
					    NodalQ(j,i)=node(j).velocity(3)
                    else
                        NodalQ(j,i)=0.0d0
                    endif                    
					!NodalQ(j,i)=node(j).velocity(3)
				end do
			case(discharge)
				IDISQ1=I
				NodalQ(:,i)=node.q
			case(kr_spg)
				NodalQ(:,i)=node.kr
			case(mw_spg)
				NodalQ(:,i)=node.mw
			case(SFR)
				do j=1,6
					NodalQ(:,i+j-1)=node.SFR(j)
				enddo
			case(NF)
				do j=1,NDIMENSION
					DO K=1,NNUM
						NodalQ(K,i+j-1)=NI_NodalForce(NODE(K).DOF(J))
					ENDDO
				enddo
		end select
		
		i=i+outvar(vo(i)).nval
		
	end do
	
	do i=1,nnum
		write(file_unit,999) (NodalQ(i,j),j=1,nvo)
	end do
	
	IF(NDATAPOINT>0) THEN
		IF(ISTEP==1.AND.ISUBTS==1) THEN
			IF(SOLVER_CONTROL.ISPARASYS<=1) THEN
                OPEN(DATAPOINT_UNIT,FILE=resultfile,STATUS='REPLACE')
                WRITE(DATAPOINT_UNIT,100) (OUTVAR(VO(I)).NAME(1:15),I=1,NVO)
            ELSE
                OPEN(DATAPOINT_UNIT,FILE=resultfile,STATUS='UNKNOWN',ACCESS='APPEND')
            ENDIF
			
		END IF
		DO I=1,NDATAPOINT		
			IF (DATAPOINT(I).ISSUMQ==0) THEN
				DO J=1,DATAPOINT(I).NNODE
					write(DATAPOINT_UNIT,110) SOLVER_CONTROL.ISPARASYS,SOLVER_CONTROL.CaseID,I,ISTEP,ISUBTS,ITER,J,(NodalQ(DATAPOINT(I).NODE(J),K),K=1,nvo)
				END DO
			ELSE
				SUMQ1=0.0				
				DO J=1,DATAPOINT(I).NNODE					
					SUMQ1=SUMQ1+NodalQ(DATAPOINT(I).NODE(J),IDISQ1)
				END DO
				NVO1=1
				WRITE(DATAPOINT_UNIT,120) SOLVER_CONTROL.ISPARASYS,SOLVER_CONTROL.CaseID,I,ISTEP,ISUBTS,ITER,J,SUMQ1,'SUMQ'
				
			ENDIF
		
		END DO		
		
        !ISFIRSTCALL1=.FALSE.
	END IF
	

	
	
999 format(<nvo>E15.7)
100	FORMAT("IPARASYS",7X,"CaseID",9X,"DATASET",8X,"ISTEP",10X,"ISUBTS",9X,"ITER",11X,"NO",13X,<NVO>A15)
110 FORMAT(7(I14,X),<nvo>(E14.7,X))
120 FORMAT(7(I14,X),<NVO1>(E14.7,X),"//",<NVO1>A15)

	!if(allocated(NodalQ)) deallocate(NodalQ)
	
end subroutine

subroutine pointout_barfamily(file_unit,ieset)
	use solverds
	implicit none
	integer,intent(in)::file_unit,ieset
	integer::i,j,k,idof,n1
	real(kind=dpn)::Vector1(3)=0.0d0,x1,y1,z1
	real(kind=dpn),allocatable::NODALQ1(:,:),DisILS(:,:),QMILS(:,:)
	
	n1=eset(ieset).noutorder
    !IF(.NOT.ALLOCATED(NODALQ)) allocate(NodalQ(n1,nvo))
	allocate(NodalQ1(n1,nvo),DisILS(n1,6),QMILS(n1,6))
	
	if(vo(3)/=locz) stop "'Z' must be output in a Barfamily element. Please add it by modifying the 'OUTDATA' command."
	
	!displacement in local system
	
	call Dis_Rot_ILS_barfamily(ieset,DisILS)
	call Shear_Moment_ILS_barfamily(ieset,QMILS)
	
	do i=1,nvo
		select case(vo(i))
			case(locx)
				do j=1,n1
					NodalQ(j,i)=node(eset(ieset).outorder(j)).coord(1)
				end do
			case(locy)
				do j=1,n1
					NodalQ(j,i)=node(eset(ieset).outorder(j)).coord(2)
				end do
			case(locz)
				do j=1,n1
					NodalQ(j,i)=node(eset(ieset).outorder(j)).coord(3)
				end do
				
			case(disx,disy,disz)
				idof=vo(i)-3
				do j=1,n1					
					if(outvar(vo(i)).system/=0) then
						NodalQ(j,i)=DisILS(j,idof)
					else
						NodalQ(j,i)=tdisp(node(eset(ieset).outorder(j)).dof(idof))
					end if
				end do
											
			case(Rx,Ry,Rz)
				idof=vo(i)-45
				do j=1,n1
					
					if((outvar(vo(i)).system/=0).and.(ndimension==3)) then
						NodalQ(j,i)=DisILS(j,idof)
					else
						NodalQ(j,i)=tdisp(node(eset(ieset).outorder(j)).dof(idof))
					end if
				end do

			case(Qx,Qy,Qz) !56,57,58
				!在局部坐标下，Qx=QMILS(:,1),Qy=QMILS(:,2),Qz=QMILS(:,3),Mx=QMILS(:,4),,My=QMILS(:,5),Mz=QMILS(:,6)
				
				idof=vo(i)-55
!				if(vo(i)==Qz) idof=vo(i)-54
				
				NodalQ(:,i)=QMILS(:,idof)
								
			case(Mx,My,Mz) !53,54,55
				idof=vo(i)-49
!				if(vo(i)==Mz) idof=vo(i)-52
				
				NodalQ(:,i)=QMILS(:,idof)				
			
			case(xf_out)
				do j=1,n1
					NodalQ(j,i)=node(eset(ieset).outorder(j)).coord(1)+tdisp(node(eset(ieset).outorder(j)).dof(1))
				end do
			case(yf_out)
				do j=1,n1
					NodalQ(j,i)=node(eset(ieset).outorder(j)).coord(2)+tdisp(node(eset(ieset).outorder(j)).dof(2))
				end do
			case(zf_out)
				do j=1,n1
					NodalQ(j,i)=node(eset(ieset).outorder(j)).coord(3)+tdisp(node(eset(ieset).outorder(j)).dof(3))
				end do
			
				
		end select
	end do
	
	do i=1,n1
		!1个节点对应4个节点，这4个节点只是x,y,z不同，其它节点量相同。
		do k=1,4
		
			x1=nodalQ(i,1)+eset(ieset).xyz_section(1,k)
			y1=nodalQ(i,2)+eset(ieset).xyz_section(2,k)
			z1=nodalQ(i,3)+eset(ieset).xyz_section(3,k) !对于杆系单元，在命令outdata中，必须含有Z.
			write(file_unit,999) x1,y1,z1,(NodalQ(i,j),j=4,nvo)
		
		end do
	end do
	
999 format(<nvo>E15.7)
	if(allocated(NodalQ1)) deallocate(NodalQ1)
	if(allocated(DisILS)) deallocate(DisILS)
	if(allocated(QMILS)) deallocate(QMILS)
	
end subroutine

subroutine Dis_Rot_ILS_barfamily(ieset,disILS)
	use solverds
	implicit none
	integer,intent(in)::ieset
	real(kind=DPN),intent(out)::disILS(eset(ieset).noutorder,6)
	integer::i,j,n1
	real(kind=DPN)::Vector1(3)
	
	n1=eset(ieset).noutorder
	disILS=0.0d0
	do j=1,n1
		Vector1(1:ndimension)=tdisp(node(eset(ieset).outorder(j)).dof(1:ndimension))
		disILS(j,1:ndimension)=matmul(element(eset(ieset).enums).g2l(1:ndimension,1:ndimension),Vector1(1:ndimension))		
	end do
	
	if(eset(ieset).et==beam) then	
		do j=1,n1
			Vector1(1:ndimension)=tdisp(node(eset(ieset).outorder(j)).dof(5:7))
			disILS(j,4:6)=matmul(element(eset(ieset).enums).g2l(1:ndimension,1:ndimension),Vector1(1:ndimension))		
		end do
	
	end if
	
end subroutine

subroutine Shear_Moment_ILS_barfamily(ieset,QM)
	use solverds
	implicit none
	integer,intent(in)::ieset
	real(kind=DPN),intent(out)::QM(eset(ieset).noutorder,6)
	integer::i,j,n1,n2,n3,n4,ndof1
	
	n1=eset(ieset).noutorder
	QM=0.0d0
	
	n4=1
	if(eset(ieset).et==beam2d.or.eset(ieset).et==ssp2d) n4=3
	if(eset(ieset).et==beam) n4=6
	
	do i=1,n1
		
		n2=eset(ieset).elist(1,i) !eset(ieset).elist(1) 肯定不为零；elist(2)为零或与elist(1)反号。
		n3=eset(ieset).elist(2,i)
		if(n3/=0) then
			if(n2>0) then
				n2=eset(ieset).elist(2,i)
				n3=eset(ieset).elist(1,i)
			end if
			
			!第(-n2)单元中，第2个节点为eset(ieset).noutorder(i)
			!第n3单元中，第1个节点为eset(ieset).noutorder(i)		
			do j=1,n4
				QM(i,j)=(element(-n2).gforceILS(j+n4)+element(n3).gforceILS(j))/2.
			end do
		else
			if(n2<0) then
				!第(-n2)单元中，第2个节点为eset(ieset).noutorder(i)
				do j=1,n4
					QM(i,j)=element(-n2).gforceILS(j+n4)
				end do				
			else 
				!第(n2)单元中，第1个节点为eset(ieset).noutorder(i)
				do j=1,n4
					QM(i,j)=element(n2).gforceILS(j)
				end do				
			end if
		end if
	end do
	
	if(eset(ieset).et==beam2d.or.eset(ieset).et==ssp2d) then
		QM(:,6)=QM(:,3)
		QM(:,3)=0.0D0
	end if

	
end subroutine


!when the datapacking format BlOCK is activated.the quantivities output are
!located at element nodes. 
!This subroutine is used to output these nodal quantivites.
subroutine BlOCKout(file_unit)
	use solverds
	implicit none
	integer::i,j,k,file_unit,idof
	real(kind=DPN),allocatable::dis(:)
	
	allocate(dis(nnum))
	i=1
	do while(i<=nvo)
		select case(vo(i))
			case(locx)
				write(file_unit,999)  node.coord(1)
			case(locy)
				write(file_unit,999)  node.coord(2)
			case(locz)
				write(file_unit,999)  node.coord(3)
			case(disx)
				idof=1
				if(outvar(disx).system/=0) then					
					call GS2LS_Displacement(outvar(disx).system,idof,dis,nnum)
					write(file_unit,999)  dis
				else
					write(file_unit,999)  tdisp(node.dof(idof))
				end if
				
			case(disy)
				idof=2
				if(outvar(disy).system/=0) then					
					call GS2LS_Displacement(outvar(disy).system,idof,dis,nnum)
					write(file_unit,999)  dis
				else
					write(file_unit,999)  tdisp(node.dof(idof))
				end if
			case(disz)
				idof=3
				if(outvar(disz).system/=0) then					
					call GS2LS_Displacement(outvar(disz).system,idof,dis,nnum)
					write(file_unit,999)  dis
				else
					write(file_unit,999)  tdisp(node.dof(idof))
				end if
			case(head)
				idof=4
				write(file_unit,999)  tdisp(node.dof(idof))
			case(Phead)
				idof=4
				if(ndimension==2) write(file_unit,999) tdisp(node.dof(idof))-node.coord(2)
				if(ndimension==3) write(file_unit,999) tdisp(node.dof(idof))-node.coord(3)	
			case(Rx)
				idof=5
				if(outvar(rx).system/=0) then					
					call GS2LS_Displacement(outvar(rx).system,idof,dis,nnum)
					write(file_unit,999)  dis
				else
					write(file_unit,999)  tdisp(node.dof(idof))
				end if
			case(Ry)
				idof=6
				if(outvar(ry).system/=0) then					
					call GS2LS_Displacement(outvar(ry).system,idof,dis,nnum)
					write(file_unit,999)  dis
				else
					write(file_unit,999)  tdisp(node.dof(idof))
				end if
			case(Rz)
				idof=7
				if(outvar(rz).system/=0) then					
					call GS2LS_Displacement(outvar(rz).system,idof,dis,nnum)
					write(file_unit,999)  dis
				else
					write(file_unit,999)  tdisp(node.dof(idof))
				end if				
			case(sxx)
			  write(file_unit,999)  (node(j).stress(1),j=1,20)
			case(syy)
			  write(file_unit,999)  (node(j).stress(2),j=1,20)
			case(szz)
			  write(file_unit,999)  (node(j).stress(3),j=1,20)
			case(sxy)
				write(file_unit,999)  (node(j).stress(4),j=1,20)	
			case(syz)
				write(file_unit,999)  (node(j).stress(5),j=1,20)
			case(szx)
				write(file_unit,999)  (node(j).stress(6),j=1,20)
			case(exx)
				write(file_unit,999)  (node(j).strain(1),j=1,20)
			case(eyy)
				write(file_unit,999)  (node(j).strain(2),j=1,20)
			case(ezz)
				write(file_unit,999)  (node(j).strain(3),j=1,20)
			case(exy)
				write(file_unit,999)  (node(j).strain(4),j=1,20)
			case(eyz)
				write(file_unit,999)  (node(j).strain(5),j=1,20)
			case(ezx)
				write(file_unit,999)  (node(j).strain(6),j=1,20)
			case(pexx)
				write(file_unit,999)  (node(j).pstrain(1),j=1,20)
			case(peyy)
				write(file_unit,999)  (node(j).pstrain(2),j=1,20)
			case(pezz)
				write(file_unit,999)  (node(j).pstrain(3),j=1,20)
			case(pexy)
				write(file_unit,999)  (node(j).pstrain(4),j=1,20)
			case(peyz)
				write(file_unit,999)  (node(j).pstrain(5),j=1,20)
			case(pezx)
				write(file_unit,999)  (node(j).pstrain(6),j=1,20)
			case(sigma_mises)
				write(file_unit,999)  node.mises
			case(eeq)
				write(file_unit,999)  node.eeq
			case(peeq)
				write(file_unit,999)  node.peeq
			case(xf_out)
				write(file_unit,999) node.coord(1)+tdisp(node.dof(1))
			case(yf_out)
				write(file_unit,999) node.coord(2)+tdisp(node.dof(2))
			case(zf_out)
				write(file_unit,999) node.coord(3)+tdisp(node.dof(3))	
			case(gradx)
				write(file_unit,999) (node(j).igrad(1),j=1,20)		
			case(grady)
				write(file_unit,999) (node(j).igrad(2),j=1,20)	
			case(gradz)
				write(file_unit,999) (node(j).igrad(3),j=1,20)
			case(vx)
				write(file_unit,999) (node(j).velocity(1),j=1,20)
			case(vy)
				write(file_unit,999) (node(j).velocity(2),j=1,20)
			case(vz)
				write(file_unit,999) (node(j).velocity(3),j=1,20)
			case(discharge)
				write(file_unit,999) node.q
			case(kr_spg)
				write(file_unit,999) node.kr
			case(mw_spg)
				write(file_unit,999) node.mw
			case(SFR)
				do j=1,6
					write(file_unit,999) node.SFR(j)
				enddo
			case(NF)
				do j=1,NDIMENSION
					write(file_unit,999) NI_NodalForce(NODE.DOF(J))
				enddo
		end select		
		i=i+outvar(vo(i)).nval
	end do

999 format(20E15.7)
		
end subroutine

subroutine tecplot_zonetitle(iincs,iiter,isfirstcall,isubts)
	use solverds
	implicit none
	integer::i,j,k,nc,n1,iincs,iiter,isubts,iset1
    real(kind=dpn)::t1=0.d0
	logical::isfirstcall,isbarfamily
	character(1024)::cstring='',cstring2='',cstring3=''
	character(48)::cword1='',cword2='',cword3='',cword4='',cword5='',cword6=''
	character(48)::cword7='',cword8='',cword9=''
	
	
	do i=1,neset
		nzone_tec=nzone_tec+1
        iset1=esetid(i)
		write(cword1,*) i
		
		
		isbarfamily=eset(iset1).et==bar.or.eset(iset1).et==bar2d.or.eset(iset1).et==beam.or.eset(iset1).et==beam2d.or.eset(iset1).et==ssp2d
		
		if(isbarfamily) then
			write(cword2,*) 4*eset(iset1).noutorder !将杆系单元以六面体实体单元输出，方便后处理。
		else
			write(cword2,*) nnum
		end if
		
		
		
		select case(eset(iset1).et)
			case(cpe6,cps6,CAX6,&
					 cpe6_spg,cps6_spg,CAX6_spg, &
					 cpe6_cpl,cps6_cpl,cax6_cpl)
				n1=(eset(iset1).enume-eset(iset1).enums+1)*4
			case(cpe15,cps15,CAX15, &
					 cpe15_spg,cps15_spg,CAX15_spg, &
					 cpe15_cpl,cps15_cpl,CAX15_cpl)
				n1=(eset(iset1).enume-eset(iset1).enums+1)*16
			case(cpe8,cps8,CAX8,cpe8r,cps8r,CAX8R, &
					 cpe8_spg,cps8_spg,CAX8_spg,cpe8r_spg,cps8r_spg,CAX8R_spg, &
					 cpe8_cpl,cps8_cpl,CAX8_cpl,cpe8r_cpl,cps8r_cpl,CAX8R_cpl) 
				n1=(eset(iset1).enume-eset(iset1).enums+1)*5
			case(prm15,prm15_spg,prm15_cpl)
				n1=(eset(iset1).enume-eset(iset1).enums+1)*14
			case(tet10,tet10_spg,tet10_cpl)
				n1=(eset(iset1).enume-eset(iset1).enums+1)*8
			case default
				n1=eset(iset1).enume-eset(iset1).enums+1
		end select
		write(cword3,*) n1
		if(solver_control.datapaking) then
			write(cword4,*) 'point'
		else
			n1=0
			do j=1,nvo
				if(outvar(vo(j)).iscentre)then
				n1=n1+1	
				write(cword6,*) vo(j)
				cword5=trim(adjustL(cword5))//','//trim(adjustL(cword6))
				end if
			end do
			if(n1>0) then
				write(cword4,*) 'block,varlocation=(['//trim(adjustL(cword5))//']=cellcentered)'
			else
				write(cword4,*) 'block'
			end if
		end if
		if(len_trim(eset(iset1).grouptitle)==0) then
			eset(iset1).zonetitle ='ZONE,T=ESET'//trim(adjustL(cword1))//',N='//trim(adjustL(cword2))//',E=' &
					//trim(adjustL(cword3))//',ZONETYPE=' &
					//trim(adjustL(eset(iset1).stype))//',DATAPACKING='//trim(cword4) 
		else
			eset(iset1).zonetitle ='ZONE,T='//trim(adjustL(eset(iset1).grouptitle))//',N='//trim(adjustL(cword2))//',E=' &
					//trim(adjustL(cword3))//',ZONETYPE=' &
					//trim(adjustL(eset(iset1).stype))//',DATAPACKING='//trim(cword4) 		
		endif
        t1=0.d0
        do j=1,iincs
            if(j<iincs) then
                n1=timestep(j).nsubts
            else
                n1=isubts
            end if
            do k=1,n1
                t1=t1+timestep(j).subts(k)
            end do
        end do
		write(cword7,'(E15.7)') t1
!		write(cword8,'(i8)') iiter
		eset(iset1).zonetitle=trim(eset(iset1).zonetitle)//',StrandID='//trim(adjustL(cword1))//',Solutiontime='//trim(adjustL(cword7))
		if(.not.isfirstcall) then
			eset(iset1).zonetitle=trim(eset(iset1).zonetitle)//',connectivitysharezone='//trim(adjustL(cword1))
		end if
		
		if(i==1) varsharezone_tec=nzone_tec
		if((i>1).and.(.not.isbarfamily)) then
			write(cword1,*) nvo
			write(cword2,*) varsharezone_tec
			eset(iset1).zonetitle=trim(eset(iset1).zonetitle)//',VARSHARELIST=([1-'//trim(adjustL(cword1))//']=' &
												//trim(adjustL(cword2))//')'
		
		end if
		
	end do
end subroutine

!calculate the transform matrix for a vector from a cartesia system to a cylinder system
function Cart2cylind(CoordInCart,center)
	implicit none
	real(8)::Cart2Cylind(3,3)
	real(8),intent(in)::CoordInCart(:),center(:)	
	!center() is the center coordinate of the syliner in globel system 
	real(8)::r1,cos_phi=0.,sin_phi=0.
	
	r1=((CoordInCart(1)-center(1))**2+(CoordInCart(2)-center(2))**2)**0.5
	cos_phi=(CoordInCart(1)-center(1))/r1
	sin_phi=(CoordInCart(2)-center(2))/r1
	Cart2cylind=0.0
	Cart2cylind(1,1)=cos_phi
	Cart2cylind(2,2)=cos_phi
	Cart2cylind(1,2)=sin_phi
	Cart2cylind(2,1)=-sin_phi
	Cart2cylind(3,3)=1.0		
end function


subroutine Vec_Cart2spherical(CoordInCart,VecInCart,VecInsperical,center)
	use solverds
	implicit none
	real(8),intent(in)::CoordInCart(:),VecInCart(:),center(:)
	!center() is the center coordinate of the sphere in globel system 
	real(8),intent(out)::VecInsperical(:)
	real(8)::R2,r1,sin_sita=0.,cos_sita=0.,cos_phi=0.,sin_phi=0.,Cart2spherical(3,3)=0.0
	
	r1=((CoordInCart(1)-center(1))**2+(CoordInCart(2)-center(2))**2)**0.5
	cos_phi=(CoordInCart(1)-center(1))/r1
	sin_phi=(CoordInCart(2)-center(2))/r1
	R2=(r1**2+(CoordInCart(3)-center(3))**2)**0.5 !R2=Big R
	sin_sita=r1/R2
	cos_sita=(CoordInCart(3)-center(3))/R2
	Cart2spherical=0.0
	Cart2spherical(1,1)=sin_sita*cos_phi
	Cart2spherical(2,2)=cos_sita*sin_phi
	Cart2spherical(1,2)=cos_sita*cos_phi
	Cart2spherical(1,3)=-sin_phi
	Cart2spherical(2,1)=sin_sita*sin_phi
	Cart2spherical(2,3)=cos_phi
	Cart2spherical(3,1)=cos_sita
	Cart2spherical(3,2)=-sin_sita	
	VecInsperical=matmul(Cart2spherical,VecInCart)
	
end subroutine

subroutine GS2LS_Displacement(isystem,idf,dis,ndis)
	use solverds
	implicit none
	integer,intent(in)::isystem,idf,ndis
	real(kind=DPN),intent(out)::dis(ndis)
	integer::i
	real(8)::VecInGS(3)=0,CoordInCart(3)=0,TransMatrix(3,3)=0
	
	interface
		function Cart2cylind(CoordInCart,center)
			real(8)::Cart2cylind(3,3)
			real(8),intent(in)::CoordInCart(:),center(:)
		end function
	end interface

	
	do i=1,ndis
		select case(isystem)
			case(sys_cylinder)
				CoordInCart=node(i).coord
				TransMatrix=Cart2cylind(CoordInCart,Origin_Sys_cylinder)				
			case(sys_sphere)
				print *, "to be updated. SUB GS2LS_Displacement"
				pause
			case default
				TransMatrix=coordinate(isystem).c				
		end select

		if(idf<=3) then 
			VecInGS(1:ndimension)=tdisp(node(i).dof(1:ndimension)) !displacement
			dis(i)=dot_product( &
			& TransMatrix(idf,1:ndimension),VecInGS(1:ndimension))
		else
			if(ndimension==3) then
				VecInGS(1:3)=tdisp(node(i).dof(5:7)) !Rotatioin
				dis(i)=dot_product( &
				& TransMatrix(idf-4,1:3),VecInGS(1:3))
			end if
		end if
	
	end do
end subroutine



subroutine generate_brickelement_bar_beam(ieset)
	use solverds
	implicit none
	integer,intent(in)::ieset
	integer::i,j,k,n1,n2,n3,n4,nc1=0
	integer,allocatable::auxnode1(:)
	
	real(kind=DPN)::hy1=0.0d0,hz1=0.d0,transm1(3,3)
	
	n1=eset(ieset).enums
	n2=eset(ieset).enume
	allocate(eset(ieset).xyz_section(3,4))
    eset(ieset).xyz_section=0.d0
	
	!!!!假定同一单元集的材料是一样的,局部坐标是一样。
	if(eset(ieset).et==bar.or.eset(ieset).et==bar2d) hy1=material(element(n1).mat).property(3)
	if(eset(ieset).et==beam.or.eset(ieset).et==beam2d.or.eset(ieset).et==ssp2d) hy1=material(element(n1).mat).property(7)
	if(eset(ieset).et==bar.or.eset(ieset).et==bar2d) hz1=material(element(n1).mat).property(4)
	if(eset(ieset).et==beam.or.eset(ieset).et==beam2d.or.eset(ieset).et==ssp2d) hz1=material(element(n1).mat).property(8)
	
	if(abs(hy1)<1e-7) hy1=0.1
	if(abs(hz1)<1e-7) hz1=0.1
	!!!!在局部坐标下的坐标
	eset(ieset).xyz_section(2,1)=-hy1/2.
	eset(ieset).xyz_section(3,1)=-hz1/2.
	eset(ieset).xyz_section(2,2)=hy1/2.
	eset(ieset).xyz_section(3,2)=-hz1/2.
	eset(ieset).xyz_section(2,3)=hy1/2.
	eset(ieset).xyz_section(3,3)=hz1/2.
	eset(ieset).xyz_section(2,4)=-hy1/2.
	eset(ieset).xyz_section(3,4)=hz1/2.
	
	!!!!假定同一单元集的局部坐标是一样。
	transm1=transpose(element(n1).g2l)
	
	do i=1,4
		eset(ieset).xyz_section(:,i)=matmul(transm1,eset(ieset).xyz_section(:,i))	
    end do	
	
    allocate(auxnode1(nnum))
    auxnode1=-1
    
	nc1=0
	do i=n1,n2
		allocate(element(i).node2(element(i).nnum))
		do j=1,element(i).nnum
			n3=element(i).node(j)
			if(auxnode1(n3)==-1) then
				nc1=nc1+1
				auxnode1(n3)=nc1
			end if
			element(i).node2(j)=auxnode1(n3)		
		end do
	end do
	
	eset(ieset).noutorder=count(auxnode1>0)
	allocate(eset(ieset).outorder(eset(ieset).noutorder), &
			 eset(ieset).elist(2,eset(ieset).noutorder))
	
	 
	 
	do i=1,nnum
		if(auxnode1(i)>0) then
			eset(ieset).outorder(auxnode1(i))=i
		end if
	end do
	
	eset(ieset).elist=0
	do i=n1,n2
		do j=1,element(i).nnum
			n3=element(i).node2(j)
			if(minval(abs(eset(ieset).elist(:,n3)))/=0) stop "Error in sub generate_brickelement_bar_beam"
			
			n4=minloc(abs(eset(ieset).elist(:,n3)),dim=1)
			
			if(j==1) eset(ieset).elist(n4,n3)=i
			if(j==2) eset(ieset).elist(n4,n3)=-i !表示第i个单元的第2个节点为n3.
		end do
	end do
	
	deallocate(auxnode1)

end subroutine

subroutine pointout_barfamily_diagram()
	use solverds
	implicit none
	integer::i,j,k,K1,NL1,NL2,N2,N3,zc1=0,enum1=0,nnum1=0,n1,LSTR,OUTV(12)=0, &
           & file_diagram,barfamily_res,EF=0,iset1
	PARAMETER(LSTR=1024)
	integer,save::zonenum=0
	real(kind=DPN)::t1,VecILS1(12,2),V1(3),V2(3),LCS(3,3)
	real(kind=DPN),ALLOCATABLE::VRES1(:,:)
	character(16)::cword1,color(6),CSD(3)
	logical::ext=.false.
	CHARACTER(LSTR)::CSTRING
	
    INTERFACE    
        SUBROUTINE SKIPCOMMENT(BARFAMILY_RES,EF)
            INTEGER,INTENT(IN)::BARFAMILY_RES
            INTEGER,OPTIONAL::EF
        END SUBROUTINE
    END INTERFACE
    
	!tecplot variables
	file_diagram=20
	barfamily_res=30
	
	INQUIRE (FILE=resultfile22, EXIST=EXT)
	if(ext) then
		open(unit=barfamily_res,file=resultfile22,status='OLD')
	else
		RETURN 
    end if
	
    CALL BARFAMILY_DIAGRAM_SCALE_CAL()
	open(unit=file_diagram,file=resultfile21,status='replace')
	
	if(ndimension==3) then
		write(file_diagram,11) TITLE
	else
		write(file_diagram,10) TITLE
    end if
	
	EF=0
	DO WHILE(EF==0)
	
		DO I=1,NESET
            iset1=esetid(i)
			IF(eset(iset1).EC/=STRU) CYCLE
			
			!同一单元集中的所有单元的局部坐标一样
			v1=(/0.d0,1.d0,0.d0/)
			v1=matmul(transpose(element(eset(iset1).enums).g2l),v1)
			
			zc1=0
			enum1=eset(iset1).enume-eset(iset1).enums+1
			nnum1=4*enum1
			
			!OUTV(1...12)=DISX,DISY,DISZ,RX,RY,RZ,QX,QY,QZ,MX,MY,MZ
			IF(ALLOCATED(VRES1)) DEALLOCATE(VRES1)
			SELECT CASE(eset(iset1).ET)
				CASE(BAR,BAR2D)
					N1=4
					OUTV=(/DISX,0,0,0,0,0,QX,0,0,0,0,0/)
				CASE(BEAM)
					N1=24
					OUTV=(/DISX,DISY,DISZ,RX,RY,RZ,QX,QY,QZ,MX,MY,MZ/)
				CASE(BEAM2D,SSP2D)
					N1=12
					OUTV=(/DISX,DISY,0,0,0,0,QX,QY,0,0,0,MZ/)
			END SELECT
			
			ALLOCATE(VRES1(N1,eset(iset1).enumS:eset(iset1).enume))
			
			!READ DATA
!			INQUIRE(BARFAMILY_RES,IOSTAT=EF) 
			
			CALL SKIPCOMMENT(BARFAMILY_RES,EF)
            IF(EF<0) EXIT
			READ(BARFAMILY_RES,'(A<LSTR>)') CSTRING
			CALL LOWCASE(CSTRING,LSTR)
			CALL TRANSLATETOPROPERTY(CSTRING)	
			DO J=1,PRO_NUM
				IF(TRIM(ADJUSTL(PROPERTY(J).NAME))=='solutiontime') THEN
					T1=PROPERTY(J).VALUE
					EXIT
				END IF
			END DO
			CALL SKIPCOMMENT(BARFAMILY_RES)
			READ(BARFAMILY_RES,*) ((VRES1(K,J),K=1,N1),J=eset(iset1).enumS,eset(iset1).enume)
			
			
			
			ZC1=0
			
			DO J=1,12
				IF(OUTV(J)==0) CYCLE
				IF(.NOT.ANY(VO(1:NVO)==OUTV(J),DIM=1)) CYCLE
				
				NL2=LEN_TRIM(ADJUSTL(OUTVAR(OUTV(J)).NAME))
				
				ZC1=ZC1+1
				!TECPLOT ZONETITLE
				if(zc1==1) then
					write(file_diagram,20) trim(adjustL(outvar(OUTV(J)).name)),nnum1,enum1,J,t1
				else
					write(file_diagram,21) trim(adjustL(outvar(OUTV(J)).name)),nnum1,enum1,J,t1,zonenum+1
				end if
				
				DO K=eset(iset1).ENUMS,eset(iset1).ENUME
					VECILS1=0.D0
					
					SELECT CASE(eset(iset1).ET)
						CASE(BAR,BAR2D)
							DO K1=1,2
								VECILS1(1,K1)=VRES1((K1-1)*N1/2+1,K)
								VECILS1(7,K1)=VRES1((K1-1)*N1/2+2,K)
							END DO
						CASE(BEAM2D,SSP2D)
							DO K1=1,2
								VECILS1(1:2,K1)=VRES1((K1-1)*N1/2+1:(K1-1)*N1/2+2,K)
								VECILS1(6,K1)=VRES1((K1-1)*N1/2+3,K)
								VECILS1(7:8,K1)=VRES1((K1-1)*N1/2+4:(K1-1)*N1/2+5,K)
								VECILS1(12,K1)=VRES1((K1-1)*N1/2+6,K)
							END DO
						
						CASE(BEAM)
							DO K1=1,2
								VECILS1(1:12,K1)=VRES1((K1-1)*N1/2+1:(K1-1)*N1/2+12,K)
							END DO
					END SELECT
					
					SELECT CASE(OUTV(J))
						CASE(DISX,DISY,DISZ)
							N2=OUTV(J)-3 !4,5,6
						CASE(RX,RY,RZ)
							N2=OUTV(J)-46 !50,51,52
						CASE(QX,QY,QZ)
							N2=OUTV(J)-49  !56,57,58
						CASE(MX,MY,MZ)
							N2=OUTV(J)-43  !53,54,55							
					END SELECT
					
					WRITE(file_diagram,30) NODE(ELEMENT(K).NODE(1)).COORD(1:NDIMENSION),VECILS1(N2,1)
					WRITE(file_diagram,30) NODE(ELEMENT(K).NODE(2)).COORD(1:NDIMENSION),VECILS1(N2,2)
					V2(1:NDIMENSION)=V1(1:NDIMENSION)*VECILS1(N2,2)*BARFAMILY_DIAGRAM_SCALE(N2)+NODE(ELEMENT(K).NODE(2)).COORD(1:NDIMENSION)
					WRITE(file_diagram,30) V2(1:NDIMENSION),VECILS1(N2,2)
					V2(1:NDIMENSION)=V1(1:NDIMENSION)*VECILS1(N2,1)*BARFAMILY_DIAGRAM_SCALE(N2)+NODE(ELEMENT(K).NODE(1)).COORD(1:NDIMENSION)
					WRITE(file_diagram,30) V2(1:NDIMENSION),VECILS1(N2,1)
					
				END DO
				
				IF(ZC1==1) THEN
					N3=0
					DO K=eset(iset1).ENUMS,eset(iset1).ENUME						
						WRITE(file_diagram,'(4(I7,X))') N3+1,N3+2,N3+3,N3+4
						N3=N3+4
					END DO
				END IF
			
			END DO
				
			zonenum=zonenum+zc1

			
		END DO

    END DO
	
    t1=0.d0
	do i=1,ndimension
		t1=t1+(barfamily_maxxyz(i)-barfamily_minxyz(i))**2
	end do
	t1=t1**0.5
    T1=t1/SOLVER_CONTROL.BARFAMILYSCALE
    color(1:6)=(/'RED','GREEN','BLUE','CYAN','YELLOW','PURPLE'/)
    CSD(1:3)=(/'X','Y','Z'/)
    DO I=1,NESET
        iset1=esetid(i)
        IF(eset(iset1).EC/=STRU) CYCLE
        LCS=0.D0

        !同一集的局部坐标系相同
                
        lcs(1:ndimension,1:ndimension)=transpose(element(eset(iset1).enums).g2l(1:ndimension,1:ndimension))*T1
                                      
        DO J=1,NDIMENSION
            IF(NDIMENSION==3) WRITE(file_diagram,40) color(mod(i,6))
            IF(NDIMENSION==2) WRITE(file_diagram,41) color(mod(i,6))
            WRITE(file_diagram,'(I3)') 1
            WRITE(file_diagram,'(I3)') 2
            WRITE(file_diagram,'(3E15.7)') node(element(eset(iset1).enums).node(1)).coord(1:ndimension) 
            WRITE(file_diagram,'(3E15.7)') (LCS(K,J)+node(element(eset(iset1).enums).node(1)).coord(K),K=1,NDIMENSION)
            IF(NDIMENSION==3) WRITE(file_diagram,50) (LCS(K,J)+node(element(eset(iset1).enums).node(1)).coord(K),K=1,NDIMENSION),TRIM(ADJUSTL(CSD(J)))
            IF(NDIMENSION==2) WRITE(file_diagram,51) (LCS(K,J)+node(element(eset(iset1).enums).node(1)).coord(K),K=1,NDIMENSION),TRIM(ADJUSTL(CSD(J)))
        END DO
        
    END DO
    
    close(file_diagram)
    close(barfamily_res)
    
10	FORMAT('TITLE="',A<LEN_TRIM(TITLE)>,'"','\NVARIABLES=X,Y,V_inLS'C)
11  FORMAT('TITLE="',A<LEN_TRIM(TITLE)>,'"','\NVARIABLES=X,Y,Z,V_inLS'C)
20  FORMAT("ZONE,T=,","VinLS_",A<NL2>,",N=",I7,",E=",I7,",ZONETYPE=FEQUADRILATERAL,DATAPACKING=POINT,STRANDID=",I3,",SOLUTIONTIME=",E15.7,".")
21  FORMAT("ZONE,T=,","VinLS_",A<NL2>,",N=",I7,",E=",I7,",ZONETYPE=FEQUADRILATERAL,DATAPACKING=POINT,STRANDID=",I3,",SOLUTIONTIME=",E15.7,",CONNECTIVITYSHAREZONE=",I3,".")
30  FORMAT(<NDIMENSION+1>(E15.7,X))
40  FORMAT('GEOMETRY T=LINE3D,CS=GRID3D,C=',A6,',AST=Hollow,AAT=End,ASZ=3,AAN=20,LT=0.3')    
41  FORMAT('GEOMETRY T=LINE,CS=GRID,C=',A6,',AST=Hollow,AAT=End,ASZ=3,AAN=20,LT=0.3')
50  FORMAT('TEXT CS=GRID3D,X=',E15.7,',Y=',E15.7,',Z=',E15.7,',T="',A,'"')
51  FORMAT('TEXT CS=GRID,X=',E15.7,',Y=',E15.7,',T="',A,'"')     
end subroutine


subroutine G2L_barfamily_diagram(ienum,VecILS)
	use solverds
	implicit none
	integer,intent(in)::ienum
	real(kind=DPN),intent(out)::VecILS(12)
	integer::i,n1
    REAL(KIND=DPN)::VECILS1(12)=0.D0
	
	n1=ndimension
	if(element(ienum).et==beam2D.OR.element(ienum).et==ssp2d) n1=ndimension+1
    
    VECILS1(1:ELEMENT(IENUM).NDOF)=TDISP(ELEMENT(IENUM).G)
    
	do i=1,element(ienum).ndof,n1
		VecILS1(i:i+n1-1)=matmul(element(ienum).g2l(1:n1,1:n1),VecILS1(i:i+n1-1))	
    end do
    
    VECILS=0.D0
    SELECT CASE(ELEMENT(IENUM).ET)
        CASE(BAR2D,BAR)
            VECILS(1:NDIMENSION)=VECILS1(1:NDIMENSION)
            VECILS(7:6+NDIMENSION)=VECILS1(NDIMENSION+1:2*NDIMENSION)
        CASE(BEAM2D,SSP2D)
            VECILS(1:2)=VECILS1(1:2)
            VECILS(6)=VECILS1(3)
            VECILS(7:8)=VECILS1(4:5)
            VECILS(12)=VECILS1(6)
        CASE DEFAULT
            VECILS=VECILS1
    END SELECT
end subroutine

subroutine barfamily_outfordiagram(iincs,isubts)
	use solverds
	implicit none
	integer,intent(in)::iincs,isubts
	LOGICAL::ISFIRST1=.TRUE.
    logical,save::exists=.false.
	integer::barfamily_res,i,j,K,K1,N1,N2,iset1
	REAL(KIND=DPN)::T1,VOUT1(NVO)
	REAL(KIND=DPN)::VecILS(12)=0.0D0,forceILS(12)=0.0D0
    CHARACTER(1024)::CSTRING1=''
	
	barfamily_res=20
	if(.not.exists) then
		open(unit=barfamily_res,file=resultfile22,status='replace')	
	else
        open(unit=barfamily_res,file=resultfile22,status='OLD',ACCESS='APPEND')	
	end if
	exists=.true.
    

	
	T1=0.D0
	DO J=1,IINCS
		IF(J<IINCS) THEN
			N1=TIMESTEP(J).NSUBTS
		ELSE
			N1=ISUBTS
		END IF
		DO K=1,N1
			T1=T1+TIMESTEP(J).SUBTS(K)
		END DO
    END DO	
	
    ISFIRST1=.TRUE.
	do i=1,neset
        iset1=esetid(i)
		if(eset(iset1).ec/=stru) cycle
		SELECT CASE(eset(iset1).ET) 
			CASE(BAR,BAR2D)
				N1=1
				WRITE(barfamily_res,100) iset1,iincs,isubts,T1
			CASE(BEAM)
				N1=6
				WRITE(barfamily_res,102) iset1,iincs,isubts,T1
			CASE(BEAM2D,SSP2D)
				N1=3
				WRITE(barfamily_res,101) iset1,iincs,isubts,T1
		END SELECT

		VECILS=0.0D0
		FORCEILS=0.0D0
		DO J=eset(iset1).ENUMS,eset(iset1).ENUME
		
			call G2L_barfamily_diagram(j,VecILS)
				
			SELECT CASE(ELEMENT(J).ET) 
				CASE(BAR,BAR2D)
					 WRITE(barfamily_res,'(4(E15.7,X))') VECILS(1),ELEMENT(J).GFORCEILS(1),VECILS(7),ELEMENT(J).GFORCEILS(2)
				CASE(BEAM)
					WRITE(barfamily_res,'(24(E15.7,X))') VECILS(1:6),ELEMENT(J).GFORCEILS(1:6),VECILS(7:12),ELEMENT(J).GFORCEILS(7:12)
				CASE(BEAM2D,SSP2D)
					WRITE(barfamily_res,'(12(E15.7,X))') VECILS(1:2),VECILS(6),ELEMENT(J).GFORCEILS(1:3),VECILS(7:8),VECILS(12),ELEMENT(J).GFORCEILS(4:6)
			END SELECT
			
			
			DO K=1,2
			    N2=(K-1)*6
				DO K1=1,6
					if(ABS(VecILS(N2+K1))>BARFAMILY_DIAGRAM_SCALE(K1)) & 
						BARFAMILY_DIAGRAM_SCALE(K1)=ABS(VecILS(N2+K1))				
				END DO 

				DO K1=1,N1
					IF((K1<3).OR.(ELEMENT(J).ET/=BEAM2D.and.ELEMENT(J).ET/=ssp2d)) THEN
						
						if(ABS(ELEMENT(J).GFORCEILS(N1*(K-1)+K1))>BARFAMILY_DIAGRAM_SCALE(6+K1)) & 
						BARFAMILY_DIAGRAM_SCALE(6+K1)=ABS(ELEMENT(J).GFORCEILS(N1*(K-1)+K1))
					ELSE
						
						!MZ
						if(ABS(ELEMENT(J).GFORCEILS(N1*(K-1)+K1))>BARFAMILY_DIAGRAM_SCALE(12)) & 
						BARFAMILY_DIAGRAM_SCALE(12)=ABS(ELEMENT(J).GFORCEILS(N1*(K-1)+K1))
					END IF
					
					
				END DO
			END DO
       END DO
	
    end do
    
	
    close(barfamily_res)

100	FORMAT('BARFAMILY_EFORCE,ESET=',I3,',STEP=',I3,',SUBSTEP=',I3,',SOLUTIONTIME=',E15.7,',ET=BAR(BAR2D)\N/DX1,FX1,DX2,FX2'C)	
101	FORMAT('BARFAMILY_EFORCE,ESET=',I3,',STEP=',I3,',SUBSTEP=',I3,',SOLUTIONTIME=',E15.7,',ET=BEAM2D\N/DX1,DY1,RZ1,FX1,FY1,MZ1,DX2,DY2,RZ2,FX2,FY2,MZ2'C)	
102	FORMAT('BARFAMILY_EFORCE,ESET=',I3,',STEP=',I3,',SUBSTEP=',I3,',SOLUTIONTIME=',E15.7,',ET=BEAM\N/DX1,DY1,DZ1,RX1,RY1,RZ1,FX1,FY1,FZ1,MX1,MY1,MZ1,DX2,DY2,DZ2,RX2,RY2,RZ2,FX2,FY2,FZ2,MX2,MY2,MZ2'C)
!10	FORMAT('ZONE,T=','"','ESET ',I2,'"',',N=',I7',E=',I7,',ZONETYPE=FELINESEG,DATAPACKING=POINT,StrandID=',I2,',SOLUTIONTIME=',E15.7,',SYSTEM=',I2)
!11	FORMAT('#',X,<NDIMENSION+2*N1>(7X,A2,7X),'NODEIGS','(VALUES IN ELEMENTAL SYSTEM.)')
!20	FORMAT(<NVO>(E15.7,X))
!30	FORMAT('VARIABLES=X,Y,Z,DX,DY,DZ,RX,RY,RZ,FX,FY,FZ,MX,MY,MZ','\N# ALL VARIABLES ARE IN ELEMENTAL SYSTEM EXCEPT FOR X,Y AND Z.'C)
!31	FORMAT('VARIABLES=X,Y,DX,DY,RZ,FX,FY,MZ','\N# ALL VARIABLES ARE IN ELEMENTAL SYSTEM EXCEPT FOR X AND Y.'C)
!32	FORMAT('VARIABLES=X,Y,Z,DX,FX','\N# ALL VARIABLES ARE IN ELEMENTAL SYSTEM EXCEPT FOR X,Y AND Z.'C)
!33	FORMAT('VARIABLES=X,Y,DX,FX','\N# ALL VARIABLES ARE IN ELEMENTAL SYSTEM EXCEPT FOR X AND Y.'C)

end subroutine

 subroutine BARFAMILY_DIAGRAM_SCALE_CAL()
	use solverds
	implicit none
	integer::i,BARFAMILY_RES
	real(kind=DPN)::t1
	
	t1=0.d0
	do i=1,ndimension
		t1=t1+(barfamily_maxxyz(i)-barfamily_minxyz(i))**2
	end do
	t1=t1**0.5
    T1=t1/SOLVER_CONTROL.BARFAMILYSCALE
	
	DO I=1,12
		IF(ABS(BARFAMILY_DIAGRAM_SCALE(I))>1.D-7) BARFAMILY_DIAGRAM_SCALE(I)=t1/BARFAMILY_DIAGRAM_SCALE(I)
	END DO
	
 end subroutine

 


