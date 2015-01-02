
subroutine lpsolvefile()
	use solverds
	use dfport
	use dflib
	implicit none
	integer::i,j,k,n1,n2,n3,ef,ncstring,ntrie=0
	parameter(ncstring=1024)
	character(ncstring)::cstring=''
	character(48)::cs1,cs2,cs3
	REAL(8),allocatable::C(:)
	logical,allocatable::a1(:)
	integer(4)::ip,errnum
	real(8)::t1=0,t2=0,t3=0,rpi,yc=0
	integer::iincs=1,iiter=1,XCOL,ATCOL,VALCOL,LMCOL
	logical::iscon,isfirstcall=.true.
	character(1024)::nme
	CHARACTER(3)        drive
	CHARACTER(1024)      dir
	CHARACTER(1024)      name
	CHARACTER(1024)      ext

	rpi=pi()
	allocate(c(ndof))
	c=0.0D0
	open(2,file=resultfile1,status='replace')
	do i=1, enum
		if(element(i).isactive==0) cycle
		select case(element(i).et)			
			case(UB3)
				ntrie=ntrie+1
				n1=3
				n2=element(i).ndof-int(material(element(i).mat).property(1))
				n3=int(material(element(i).mat).property(1))				
				!gravity effect
				if(.not.solver_control.ismg) then
					t1=material(element(i).mat).property(4)
					if(material(element(i).mat).FF1D(4)>0) then
						t1=materialfun(material(element(i).mat).FF1D(4),element(i).property(3))
						t1=material(element(i).mat).property(4)*t1
					end if
					do j=2,6,2
						if(abs(t1)>1e-7)  then
							c(element(i).g(j))=-element(i).property(1)/6*t1
						end if
					end do
				end if
				t1=material(element(i).mat).property(3)
				if(material(element(i).mat).FF1D(3)>0) then
					t1=materialfun(material(element(i).mat).FF1D(3),element(i).property(3))
					t1=material(element(i).mat).property(3)*t1
				end if
				t1= t1*dcos(material(element(i).mat).property(2)/180*rpi) !cohesion*cos(friction angle)
				!pore water pressure effect
				t2=material(element(i).mat).property(5)
				if(abs(t2)>1e-7) then
					if(material(element(i).mat).FF1D(5)>0) then
						t3=materialfun(material(element(i).mat).FF1D(5),element(i).property(3))
						t2=t2*t3
					end if
					t2=t2*dsin(material(element(i).mat).property(2)/180*rpi) !Pore Pressure*sin(friction angle)
				end if
				t1=element(i).property(1)*(t1-t2) !2A*(c*cos(phi)-p*sin(phi))
				do j=1,n3					
					c(element(i).g(n2+j))=t1				
				end do		
			case(UBZT4)
				n1=4
				n2=element(i).ndof-4
				n3=4
				t1=material(element(i).mat).property(3)
				if(material(element(i).mat).FF1D(3)>0) then
					t1=materialfun(material(element(i).mat).FF1D(3),element(i).property(3))
					t1=material(element(i).mat).property(3)*t1
				end if
				!pore water pressure effect
				t2=material(element(i).mat).property(5)
				if(abs(t2)>1e-7) then
					if(material(element(i).mat).FF1D(5)>0) then
						t3=materialfun(material(element(i).mat).FF1D(5),element(i).property(3))
						t2=t2*t3
					end if
					t2=t2*dtan(material(element(i).mat).property(2)/180*rpi) !Pore Pressure*tan(friction angle)
				end if

				t1=element(i).property(1)*(t1-t2)/2.0D0	!L/2*(C-P*tan(phi))	
				do j=1,n3
					c(element(i).g(n2+j))=t1	
				end do
		end select
	end do	
	
	write(2,'(a4)') 'MIN:'
	n1=0
	cstring=''
	do i=1,ndof,20
		do j=1,20
			n1=n1+1
			if(n1>ndof) exit
			if(abs(c(n1))>1e-7) then
				write(cs1,'(f15.7)') c(n1)
				write(cs2,'(i15)')  n1
				cstring=trim(adjustl(cstring))//'+'//trim(adjustl(cs1))//'x'//trim(adjustl(cs2))
			end if
		end do
		if(j<=20.or.n1==ndof) then
			cstring=trim(adjustL(cstring))//';'
		end if	
		if(len_trim(adjustl(cstring))/=0) then	
			write(2,10) cstring
			cstring=''
		end if				
	end do

	write(2,'(/,a)') ''

	do i=1, enum
		select case(element(i).et)			
			case(UB3)
				n1=3
				n2=element(i).ndof-int(material(element(i).mat).property(1))
				n3=int(material(element(i).mat).property(1))
			case(UBZT4)
				n1=4
				n2=element(i).ndof-4
				n3=4				
		end select
		
		do j=1,n1
			do k=1,n2
				if(abs(element(i).km(j,k))>1e-7) then
					write(cs1,'(f15.7)') element(i).km(j,k)
					write(cs2,'(i15)') element(i).g(k)
					cstring=trim(adjustl(cstring))//'+'//trim(adjustL(cs1))//'x'//trim(adjustL(cs2))
				end if						
			end do
			do k=1,n3
				if(abs(element(i).A12(j,k))>1e-7) then
					write(cs1,'(f15.7)') element(i).A12(j,k)
					write(cs2,'(i15)') element(i).g(n2+k)
					cstring=trim(adjustl(cstring))//'+'//trim(adjustL(cs1))//'x'//trim(adjustL(cs2))
				end if						
			end do
			cstring=trim(adjustl(cstring))//"=0;"
			write(2,10) cstring
			cstring=''
		end do
	end do
	!gravity effect
	if(solver_control.ismg) then
		n1=0
		n2=0
		do i=1,enum
			if(element(i).et==UB3) then
				t1=element(i).property(1)
				if(material(element(i).mat).FF1D(4)>0) then
					t1=materialfun(material(element(i).mat).FF1D(4),element(i).property(3))
					t1=element(i).property(1)*t1
				end if				
				write(cs1,'(f15.7)') t1/6.0D0
				do j=1,3					
						write(cs2,'(i15)') element(i).g(2*j)
						n1=n1+1
						n2=n2+1
						cstring=trim(adjustl(cstring))//'+'//trim(adjustL(cs1))//'x'//trim(adjustL(cs2))				
				end do
				if(n1>=20.or.n2==3*ntrie) then
					if(n2==3*ntrie) then
						cstring=trim(adjustl(cstring))//"=-100.0;"
					end if
					write(2,10) cstring
					cstring=''				
					n1=0
				end if
			end if
		end do
	end if


	allocate(a1(nuvar))
	a1=.false.
	cstring=''
	do i=1,bd_num
		write(cs1,'(f15.7)') bc_disp(i).value
		write(cs2,'(i15)') node(bc_disp(i).node).dof(bc_disp(i).dof)
		a1(node(bc_disp(i).node).dof(bc_disp(i).dof))=.true.
		cstring=trim(adjustl(cstring))//'1.0'//'x'//trim(adjustl(cs2))//"="//trim(adjustL(cs1))//';'
		write(2,10) cstring
		cstring=''			
	end do

	write(2,'(/,a4)') 'FREE'
	n1=0
	n2=0
	cstring=''
	do i=1,nuvar,20
		do j=1,20
			n1=n1+1
			if(n1>nuvar) exit
			if(a1(n1)) cycle
			n2=n2+1
			write(cs2,'(i15)') n1
			cstring=trim(adjustl(cstring))//'x'//trim(adjustl(cs2))//','
		end do
		if((n2==nuvar-bd_num)) then
			cstring=adjustr(cstring)
			cstring(ncstring:ncstring)=';'
			cstring=adjustL(cstring)
		end if
		if(len_trim(cstring)>1) write(2,10) cstring
		cstring=''
	end do

	10 format(a<len_trim(cstring)>)

	close(2)

	if(.not.solver_control.ispg) then
		open(3,file=trim(adjustL(resultfile1))//'.bat', status='replace')
		write(3, *) '@ECHO OFF'
		if(solver_control.solver==lpsolver) then
			cstring='F:\program\LP_solve\lp_solve -s -S -v4 -time -R -bfp bfp_LUSOL '//trim(adjustL(resultfile1))//' >'// &
								trim(adjustL(resultfile1))//'.dat'
			write(3, 10) cstring
			cstring='::F:\program\LP_solve\lp_solve -parse_only -wmps '//trim(adjustL(resultfile1))//'.mps'// &
							' '//trim(adjustL(resultfile1))
			write(3, 10) cstring
			cstring='::mosek '//trim(adjustL(resultfile1))//'.mps'
			write(3, 10) cstring
			write(3, *) 'PAUSE'	
			close(3)
		else
			cstring='::F:\program\LP_solve\lp_solve -s -S -v4 -time -R -bfp bfp_LUSOL '//trim(adjustL(resultfile1))//' >'// &
								trim(adjustL(resultfile1))//'.dat'	
			write(3, 10) cstring
			cstring='F:\program\LP_solve\lp_solve -parse_only -wmps '//trim(adjustL(resultfile1))//'.mps'// &
							' '//trim(adjustL(resultfile1))
			write(3, 10) cstring
			write(3, *) 'PAUSE Any Problem?'	
			cstring='mosek '//trim(adjustL(resultfile1))//'.mps'
			write(3, 10) cstring							
			write(3, *) 'PAUSE'	
			close(3)									
		end if	

		

		IP = MESSAGEBOXQQ('Input file for LP_SOLVE is ready! Press OK to solve it now and CANCEL to do it later.', &
								  'Wait for your Choice'C, MB$ICONASTERISK.OR.MB$OKCANCEL .OR.MB$DEFBUTTON1)
		if(IP==MB$IDOK) then
				IP = SYSTEMqq(trim(adjustL(resultfile1))//'.bat')
				If (IP .eq. -1) then
					errnum = ierrno( )
					print *, 'Error ', errnum
					stop
				end if	
		 end if	
			
	end if

		IP = MESSAGEBOXQQ('If the result data is ready, Press OK to Continue, Else CANCEL to exit.', &
					'Wait for your Choice'C, MB$ICONASTERISK.OR.MB$OKCANCEL .OR.MB$DEFBUTTON1)
		if(IP==MB$IDOK) then

			cstring="Data Files(*.dat),*.dat;Mosek Files(*.sol),*.sol;Data Files(*.bas),*.bas;CSV Files(*.csv),*.csv;  All Files(*.*),*.*;"
			cstring=trim(cstring)
			call setmessageqq(cstring,QWIN$MSG_FILEOPENDLG)
			cstring=''
			open(2, file=' ',status='old')
			inquire(2,name=nme)
			errnum = SPLITPATHQQ(nme, drive, dir, name, ext)
			if('.csv'==trim(ext)) then
				n1=0
				do while(n1==0)				
					read(2,'(a1024)',iostat=n1) cstring
					if(n1<0) exit	
					cstring=adjustl(cstring)
					n2=index(cstring,';')
					cs1=cstring(1:n2-1)
					cs3=cs1(2:len_trim(cs1))
					n3=0
					read(cs3(1:len_trim(cs3)),'(I16)',iostat=ef) n3
					if(ef>0) cycle
					if(n3>0) then
						cs2=cstring(n2+1:len_trim(cstring))
						read(cs2(1:len_trim(cs2)),*)  t1
						load(n3)=t1
					end if

				end do
			end if

			if('.dat'==trim(ext)) then
				n1=0
				iscon=.false.
				do while(n1==0)				
					read(2,'(a1024)',iostat=n1) cstring
					if(n1<0) exit	
					cstring=adjustl(cstring)
					
					if(.not.iscon) then
						n2=index(cstring,'Actual values of the variables:')						
						if(n2==1) iscon=.true.
						cycle
					end if
					n2=index(cstring,' ')	
					cs1=cstring(1:n2-1)
					cs3=cs1(2:len_trim(cs1))
					n3=0
					read(cs3(1:len_trim(cs3)),'(I16)',iostat=ef) n3
					if(ef>0) cycle
					if(n3>0) then
						cs2=cstring(n2+1:len_trim(adjustL(cstring)))
						read(cs2(1:len_trim(cs2)),*)  t1
						load(n3)=t1
					end if
				end do			
			end if

			if('.sol'==trim(ext).or.'.bas'==trim(ext)) then
				n1=0
				iscon=.false.
				do while(n1==0)				
					read(2,'(a1024)',iostat=n1) cstring
					if(n1<0) exit	
					cstring=adjustl(cstring)
					
					if(.not.iscon) then
						n2=index(cstring,'VARIABLES')						
						if(n2==1) then
							iscon=.true.
							read(2,'(a1024)') cstring
							XCOL=index(cstring,'NAME')
							ATCOL=index(cstring,'AT')
							VALCOL=index(cstring,'ACTIVITY')
							LMCOL=index(cstring,'LOWER LIMIT')
						end if
						cycle
					end if
					cs3=cstring(xcol+1:atcol-1)
					
					!cs3=cs1(2:len_trim(cs1))
					n3=0
					read(cs3(1:len_trim(cs3)),'(I16)',iostat=ef) n3
					if(ef>0) cycle
					if(n3>0) then
						cs2=cstring(valcol:lmcol-1)
						read(cs2(1:len_trim(cs2)),*)  t1
						load(n3)=t1
					end if
				end do			
			end if
									
			close(2)
			
			ISCON=.TRUE.
			if(solver_control.solver==mosek) solver_control.solver=lpsolver
			CALL LIMITOUT(IINCS,IITER,ISCON,ISFIRSTCALL,C)

		 else
			stop			
		 end if		
 	   

end subroutine

