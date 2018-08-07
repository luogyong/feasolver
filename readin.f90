 subroutine readin(itype)
	use SOLVERLIB
    use ds_hyjump
    use ExcaDS
	use dflib
    use ifport
	implicit none
	integer:: itype,unit,i,j,k
	LOGICAL(4)::tof
	character(1024) term,keyword
	character(1024)::nme
	CHARACTER(3)        drive
	CHARACTER(1024)      dir
	CHARACTER(1024)      name
	CHARACTER(1024)      ext
	type(qwinfo) winfo
	integer(4)::length,msg
    real(DPN)::AR1(100)
	type(bc_tydef),allocatable::bf1(:),bf2(:)
    EXTERNAL::SetLineColor,Marker,LineStyle,SETBGCOLOR,solvercommand
    
	term="Solver Files(*.sinp),*.sinp; &
          Plot File(*.plot),*.plot; &  
			  Data Files(*.dat),*.dat; &
			  Prof.Cao Program files(*.z7),*.z7; &
			  All Files(*.*),*.*;"
	term=trim(term)
	call setmessageqq(term,QWIN$MSG_FILEOPENDLG)
 
	term=' '
	title=''
	open(1,file=' ',status='old' )
	unit=1
	
    inquire(1,name=nme)
    length = SPLITPATHQQ(nme, drive, dir, name, ext)
    CALL LOWCASE(EXT)
    msg = CHDIR(trim(drive)//trim(dir))
    
    IF(TRIM(ADJUSTL(EXT))=='.plot') THEN
        CALL plot_func(nme)
        STOP
    ENDIF        
   

    
    print *, 'Begin to read in data...'
	call read_execute(unit,itype,keyword,solvercommand)

	if(itype==0) then
		
		resultfile=trim(drive)//trim(dir)//trim(name)//'_datapoint.dat'
		IF(solver_control.solver==LPSOLVER) THEN
			resultfile1=trim(drive)//trim(dir)//trim(name)//'_lpsolve.lp'
		END IF
		IF(solver_control.solver==MOSEK) THEN
			resultfile1=trim(drive)//trim(dir)//trim(name)//'_mosek.lp'
		END IF

		resultfile2=trim(drive)//trim(dir)//trim(name)//'_tec.plot'
		resultfile21=trim(drive)//trim(dir)//trim(name)//'_barfamilydiagram_tec.plot'
		resultfile22=trim(drive)//trim(dir)//trim(name)//'_barfamilydiagram_res.dat'
		resultfile3=trim(drive)//trim(dir)//trim(name)//'_msg.dat'
        hydraulicjumpfile=trim(drive)//trim(dir)//trim(name)//'_hjump.dat'
		EXCAMSGFILE=trim(drive)//trim(dir)//trim(name)//'_exca_msg.dat'
		EXCAB_BEAMRES_FILE=trim(drive)//trim(dir)//trim(name)//'_exca_res.dat'
		EXCAB_STRURES_FILE=trim(drive)//trim(dir)//trim(name)//'_exca_stru.dat'
		EXCAB_EXTREMEBEAMRES_FILE=trim(drive)//trim(dir)//trim(name)//'_exca_minmaxbeam.dat'
        Slope_file=trim(drive)//trim(dir)//trim(name)//'_slope_res.dat'
        helpfile=trim(drive)//trim(dir)//'IFSOLVER_HELP.TXT'
	end if
    
	open(99,file=resultfile3,status='replace')
	!the default value of Title=resultfile2.
	msg=len_trim(title)
	if(len_trim(title)==0) title=resultfile2

	close(1)
	print * ,'Read in data completed!' 
    if (solver_control.isParasys>0) then
        msg = setexitqq(QWIN$EXITNOPERSIST)
        winfo%TYPE = QWIN$MIN
    else
    	winfo%TYPE = QWIN$MAX
    endif
	tof=SETWSIZEQQ(QWIN$FRAMEWINDOW, winfo)
	tof=SETWSIZEQQ(0, winfo) 
    
    
    
	!INITIALIZATION 
    
    IF(ISEXCA2D/=0) THEN
		
		
        !CHECKDATA
        
         CALL   checksoilprofile()
      !   DO I=1,NSOILPROFILE
      !      
      !      ar1=-9.99999D20
		    !do j=1,soilprofile(i).nasoil
			   ! IF(KPOINT(NDIMENSION,soilprofile(i).asoil(j).z(1))<KPOINT(NDIMENSION,soilprofile(i).asoil(j).z(2))) THEN
      !              PRINT *, "Z1 IS SMALLER THAN Z2.PLEASE CKECK. SOILPROFILE=, ACTIVE SIDE SOILLAYER=",I,J
      !              STOP
      !          ENDIF
      !          DO K=1,NSTEP
      !              IF(ABS(SF(soilprofile(i).asoil(j).SF).FACTOR(K))>1E-6) THEN
      !                  IF(AR1(K)/=-9.99999D20)THEN
      !                      IF(ABS(AR1(K)-KPOINT(NDIMENSION,soilprofile(i).asoil(j).Z(1)))<1E-6) THEN
      !                          AR1(K)=KPOINT(NDIMENSION,soilprofile(i).asoil(j).Z(2))
      !                      ELSE
      !                          PRINT *, "THE ELEVATIONS ARE NOT CONTINUOUS BETWEEN ASOIL LAYER J-1 AND J IN STEP K,OF SOILPROFILE I. (J-1,J,K,I)=", J-1,J,K,I
      !                          STOP "INPUT ERROR IN SOILPROFILE1."
      !                      ENDIF                            
      !                  ELSE
      !                      AR1(K)=KPOINT(NDIMENSION,soilprofile(i).asoil(j).Z(2))
      !                  ENDIF                        
      !              ENDIF                    
      !          ENDDO
      !                         
      !      enddo
      !      
      !      ar1=-9.99999D20
		    !do j=1,soilprofile(i).npsoil
      !         IF(KPOINT(NDIMENSION,soilprofile(i).Psoil(j).z(1))<KPOINT(NDIMENSION,soilprofile(i).Psoil(j).z(2))) THEN
      !              PRINT *, "Z1 IS SMALLER THAN Z2.PLEASE CHECK. SOILPROFILE=, PASSIVE SIDE SOILLAYER=",I,J
      !              STOP
      !         ENDIF
      !         DO K=1,NSTEP
      !              IF(ABS(SF(soilprofile(i).Psoil(j).SF).FACTOR(K))>1E-6) THEN
      !                  IF(AR1(K)/=-9.99999D20)THEN
      !                      IF(ABS(AR1(K)-KPOINT(NDIMENSION,soilprofile(i).Psoil(j).Z(1)))<1E-6) THEN
      !                          AR1(K)=KPOINT(NDIMENSION,soilprofile(i).Psoil(j).Z(2))
      !                      ELSE
      !                          PRINT *, "THE ELEVATIONS ARE NOT CONTINUOUS BETWEEN PSOIL LAYER J-1 AND J IN STEP K,OF SOILPROFILE I. (J-1,J,K,I)=", J-1,J,K,I
      !                          STOP "INPUT ERROR IN SOILPROFILE2."
      !                      ENDIF                            
      !                  ELSE
      !                      AR1(K)=KPOINT(NDIMENSION,soilprofile(i).Psoil(j).Z(2))
      !                  ENDIF                        
      !              ENDIF                    
      !          ENDDO
      !         
      !      enddo 
      !      
      !  ENDDO
        
        CALL GenElement_EXCA2()
        solver_control.bfgm=continuum
        
        msg=INSERTMENUQQ (5, 0, $MENUENABLED, 'GraphSetting'c,NUL)
		msg=INSERTMENUQQ (5, 1, $MENUENABLED, 'Gray'c, SetLineColor)
		msg=INSERTMENUQQ (5, 2, $MENUENABLED, 'NoMarker'c, Marker)
		msg=INSERTMENUQQ (5, 3, $MENUENABLED, 'ThickLine'c, LineStyle)
        msg=INSERTMENUQQ (5, 4, $MENUENABLED, 'BGC_BLACK'c, SETBGCOLOR)
        !solver_control.ismkl=NO
    ENDIF    
    
    IF(ISSLOPE/=0) THEN
        CALL Gen_slope_model()
    ENDIF
    
    if(enum==0) then
        STOP 'NO ELEMENT WAS INPUT.'
    endif
    
	LF1D(0,1)=0.0D0
	LF1D(0,2)=1.0D0	
	if(nsf==0) then
		nsf=1
		nstep=max(1,nstep)
		allocate(sf(0:nsf))		
		allocate(sf(1).factor(0:nstep),sf(0).factor(0:nstep))
		sf(0).factor(0)=0.0d0
		sf(0).factor(1:nstep)=1.0d0
		sf(1).factor(1:nstep)=1.0d0
		sf(1).factor(0)=0.0d0
	end if
	
	if(.not.allocated(timestep)) then
		allocate(timestep(0:nstep))		
		do i=0,nstep
			timestep(I).nsubts=1
			allocate(timestep(i).subts(1)) 
			timestep(i).subts(1)=1.d0  !如为稳态分析，则此时间步长为虚步长。
		end do  
		timestep(0).nsubts=0
		Timestep(0).subts(1)=0.d0 
    end if
    
	if(.not.allocated(stepinfo)) then
		allocate(stepinfo(0:nstep))
        nstepinfo=nstep
        stepinfo(0).matherstep=0
        stepinfo(0).issteady=.true.
        do i=1,nstepinfo
            stepinfo(i).matherstep=i-1
            stepinfo(i).issteady=.true.
        end do
    end if
    
    
    if(nhjump/=0) open(unit_hj,file=hydraulicjumpfile,status='replace')
    if(HJump_CP==1) then
        do i=1,nstep
			do j=1,timestep(i).nsubts
				do k=1,nhjump
					call HJ_WaterSurfaceProfile_RC(k,i,j)
				end do
			end do
        end do
        stop "WaterSurfaceProfile Calculation Completed!"
	else
		do j=1,nhjump
			allocate(bf1(bd_num+hjump(j).nnode))
			do i=bd_num+1,bd_num+hjump(j).nnode
				bf1(i).node=hjump(j).node(i-bd_num)
				hjump(j).bc_node(i-bd_num)=i
				bf1(i).dof=4
				!bf1(i).value=ar(3)
				
				bf1(i).sf=0								

				bf1(i).isdual=1
                OUTVAR(90+BF1(I).DOF).VALUE=90+BF1(I).DOF   
			end do
			if(bd_num>0) then
				bf1(1:bd_num)=bc_disp(1:bd_num)
				deallocate(bc_disp)
			end if			
			allocate(bc_disp(bd_num+hjump(j).nnode))
			bc_disp=bf1
			bd_num=bd_num+hjump(j).nnode
			
			deallocate(bf1)	
						
			!deallocate(bf1)
		end do
    end if
    
	if(ncoord==0) then
		allocate(coordinate(-1:0))
		coordinate(0).c=0.0
		coordinate(0).c(1,1)=1.0
		coordinate(0).c(2,2)=1.0
		coordinate(0).c(3,3)=1.0
	end if
	if(nueset==1) then
		ueset(0).enum=enum
		ueset(0).name='all'
		allocate(ueset(0).element(ueset(0).enum))
		do i=1,ueset(0).enum
			ueset(0).element(i)=i
		end do
	end if
	if(nunset==1) then
		unset(0).nnum=nnum
		unset(0).name='all'
		allocate(unset(0).node(unset(0).nnum))
		do i=1,unset(0).nnum
			unset(0).node(i)=i
		end do
	end if
	
	!!!以水头为未知量的渗流模型，水头不具有叠加性，导致其多步求解时与以位移为未知量的其它模型不同，必须注意。
    !!!假定，如果为渗流模型，则模型中所有的单元均为渗流单元。
	!if(solver_control.type/=spg.and.(element(1).ec==spg2d.or.element(1).ec==spg.or.element(1).ec==cax_spg)) then
	!	solver_control.type=spg
	!end if
	
	do i=1, bd_num
		if(bc_disp(i).isdual>0) then
			do j=1,numNseep
				if((Nseep(j).node==bc_disp(i).node).and.(Nseep(j).dof==bc_disp(i).dof)) then
					bc_disp(i).isdual=j
					Nseep(j).isdual=i
					exit
				end if
			end do
			if(j>numNseep) bc_disp(i).isdual=0
		end if
	end do
	
	!intialize 
	do i=1,nSMNP
		do j=1,bd_num
			if(bc_disp(j).node==smnp(i).master.and.bc_disp(j).dof==smnp(i).mdof) then
				print *, "displacement condition is applied on the master node=",	i
				stop
			end if
		end do
		do j=1,bl_num
			if(bc_load(j).node==smnp(i).master.and.bc_load(j).dof==smnp(i).mdof) then
				smnp(i).nmbl=smnp(i).nmbl+1
				if(smnp(i).nmbl>10) stop "smnp(i).nmbl>10.作用在master节点mdof度的荷载个数最多10."
				smnp(i).mbl(smnp(i).nmbl)=j
			end if
		end do
    end do
	
    if(nfreedof>0) then
		call enlarge_node(node,nnum,nfreedof,k)
        do i=1,nfreedof
           do j=1,element(freedof(i).element).nnum
                if(element(freedof(i).element).node(j)==freedof(i).node) then
                    element(freedof(i).element).ifreedof=i
                    freedof(i).newnode=k
                    element(freedof(i).element).node(j)=k
                    node(k)=node(freedof(i).node)                    
                    k=k+1
                    exit
                endif                
           enddo                       
        enddo    
        
	endif

    if(solver_control.bfgm==inistress) solver_control.issym=.true.
    
    !OUTVAR(90).VALUE=90;OUTVAR(90).NAME='BC_TYPE'; 
    IF(OUTVAR(DISX_BC).VALUE>0) OUTVAR(DISX_BC).NAME='DISX_BC'
    IF(OUTVAR(DISY_BC).VALUE>0) OUTVAR(DISY_BC).NAME='DISY_BC'
    IF(OUTVAR(DISZ_BC).VALUE>0) OUTVAR(DISZ_BC).NAME='DISZ_BC'
    IF(OUTVAR(H_BC).VALUE>0) OUTVAR(H_BC).NAME='H_BC'
    IF(OUTVAR(RX_BC).VALUE>0) OUTVAR(RX_BC).NAME='RX_BC'
    IF(OUTVAR(RY_BC).VALUE>0) OUTVAR(RY_BC).NAME='RY_BC'
    IF(OUTVAR(RZ_BC).VALUE>0) OUTVAR(RZ_BC).NAME='RZ_BC'
    
	!if(solver_control.bfgm==lacy) then
	!	do i=1,bd_num
	!		if(bc_disp(i).dof==4) then
	!			!convert the hydraulic head boundaries to pressure head boundaries.
	!			bc_disp(i).value=bc_disp(i).value-node(bc_disp(i).node).coord(ndimension)
	!		end if
	!	end do
	!end if	
	return

 end subroutine



subroutine read_execute(unit,itype,keyword,COMMAND_PARSER)

!**************************************************************************************************************
!IF ITYPE=0, READ IN DATAS  from THE UNIT FILE TO ITS END.
!IF ITYPE>0, JUST READ IN DATA RELATED WITH THE KEYWORD BLOCK IN THE UNIT FILE
!INPUT VARIABLES:
!UNIT: FILE NUMBER, 
!ITYPE: DEFAULT VALUE=0, IF VALUE>0, IT WILL WORK WITH THE KEYWORD.  
!KEYWORD: DATA BLOCK KEYWORD
!COMMAND_PARSER: SUBROUTINE TO HANDLE COMMANDS FROM READING
!A LINE STARTED WITH '/' IS A COMMENT LINE, IT WILL BE SKIPPED DURING READING
!OUPUT VARIABLES:
!NO EXPLICIT OUTPUT VARIABLES. ALL THE READ IN DATA STORED IN THE VARIABLES DEFINED IN THE MODULE SOLVERDS
!SUBROUTINES CALLED: 
!COMMAND()
!Programer: LUO Guanyong
!Last update: 2008.03.16
!**************************************************************************************************************
	!use solverds
	implicit none
    INTEGER,INTENT(IN)::UNIT,ITYPE    
	integer::ef,iterm,i,strL,N1
	parameter(iterm=1024)
    character(*),INTENT(IN)::keyword
	character(iterm)::term,term2
	character(1)::ch
    EXTERNAL::COMMAND_PARSER
	
	ef=0
	
	do while(ef==0)
		
		term=''
		do while(.true.)
			read(unit,999,iostat=ef) term2
			if(ef<0) exit	
			term2=adjustL(term2)
			strL=len_trim(term2)
			if(strL==0.or.term2(1:2)=='//'.or.term2(1:1)=='#') cycle		

			!每行后面以'/'开始的后面的字符是无效的。
			if(index(term2,'//')/=0) then
				strL=index(term2,'//')-1
				term2=term2(1:strL)
				strL=len_trim(term2)
			end if			

			if(term2(strL:strL)/="&") then
				term=trim(adjustL(term))//trim(term2)
				exit
			else
				term=trim(adjustL(term))//term2(1:strL-1)			
			end if
		end do
		
		if(ef<0) exit
		
		term=adjustl(term)
		strL=len_trim(term)
		if(strL==0) cycle
		do i=1,strL !remove 'Tab'
			if(term(i:i)/=char(9)) exit
		end do
		term=term(i:strL)
		term=adjustl(term)
		strL=len_trim(term)
		if(strL==0) cycle		
		write(ch,'(a1)') term
		if(TERM(1:2)/='//'.and.TERM(1:1)/='#') then
			!backspace(unit)
			!read(unit,999) term
			call lowcase(term)
			call translatetoproperty(term)			
			term=adjustl(trim(term))
			call COMMAND_PARSER(term,unit)			 	
		end if
	end do


	
999	format(a<iterm>)

end subroutine


subroutine solvercommand(term,unit)

!**************************************************************************************************************
!Function:
!read in data block: "term" from the unit file.
!Input Variables:
!Term: data block keyword
!Unit: data file unit number
!Output varibibles:
!No explicit output variable. All the data read in is stored in variables defined in modulus SOLVERDS.
!Modulus Used:
!dflib, SOLVERDS
!Subroutines Called:
!strtoint()
!Programer: LUO Guanyong
!Last Update: 2008.03.16
!**************************************************************************************************************
	use dflib
	use solverlib
	use ExcaDSLIB
    use ds_hyjump
    USE DS_SlopeStability
	implicit none
    
	integer::unit
	character(1024) term
	integer::i,j,k
	integer::n1,n2,n3,n4,n5,n_toread,nmax
	real(8)::ar(MaxNumRead)=0,t1,T2
	integer(4)::msg
	type(mat_tydef),allocatable::mat1(:)
	type(element_tydef),allocatable::element1(:),element2(:)
	type(bc_tydef),allocatable::bf1(:),bf2(:)
	integer::enum1=0,nnum1=0,et1=0,set1=0,material1=0,&
		ndof1=0,nd1=0,ngp1=0,matid1=0,nset=0,sf1=0,system1=0,nbf1=0
	integer::ec1=0,id1=0,excelformat=0
	character(16)::stype
	character(128)::cstring
	character(64)::name1,set(50)
	logical::isset=.false.,TOF1=.FALSE.
	
	type(bc_tydef),allocatable::Nseep1(:)

    INTERFACE
        SUBROUTINE SKIPCOMMENT(BARFAMILY_RES,EF)
            INTEGER,INTENT(IN)::BARFAMILY_RES
            INTEGER,OPTIONAL::EF
        END SUBROUTINE
    END INTERFACE   
    
	nmax=MaxNumRead
	n_toread=MaxNumRead
	ar=0.0
	set=''

	term=trim(term)

	select case(term)
		case('title')
			print *, 'Reading TITLE data...'
			excelformat=0
			do i=1, pro_num
				select case(property(i).name)
					case('excelformat')
						excelformat=int(property(i).value)
					case default
						call Err_msg(property(i).name)
				end select
			end do
			call skipcomment(unit)			
			read(unit,'(a1024)') title
        CASE('helpfile')
            print *, 'WRITING A HELPFILE IN THE CURRENT DIR'
            call write_readme_FEASOLVER()
			do i=1, pro_num
				select case(property(i).name)
					case('exit')
						IF(INT(PROPERTY(I).VALUE)==1) STOP "DONE.A HELPFILE IS OUT IN THE CURRENT DIR."
					case default
						call Err_msg(property(i).name)
				end select
			end do            
		    
        
		case('node')
			print *, 'Reading NODE data...'
			do i=1, pro_num
				select case(property(i).name)
					case('num')
						nnum=int(property(i).value)
					case('datapacking','dp')
						datapacking=int(property(i).value)
					case('dimension','d')
						ndimension=int(property(i).value)
					case default
						call Err_msg(property(i).name)
				end select
			end do
			allocate(node(nnum))
			if (datapacking==1) then
				call skipcomment(unit)
				read(unit,*) ((node(i).coord(j),j=1,ndimension),i=1,nnum)
			else
				call skipcomment(unit)
				read(unit,*) ((node(i).coord(j),i=1,nnum),j=1,ndimension)
			end if
		case('kpoint','kp')
			print *, 'Reading KEYPOINT data...'
			do i=1, pro_num
				select case(property(i).name)
					case('num')
						nkp=int(property(i).value)
					!case('datapacking','dp')
					!	datapacking=int(property(i).value)
					!case('dimension','d')
					!	ndimension=int(property(i).value)
					case default
						call Err_msg(property(i).name)
				end select
			end do
			allocate(KPOINT(1:NDIMENSION+1,nkp),kpnode(nkp))
            KPOINT(NDIMENSION+1,:)=DEFAULTSIZE
            
			do i=1,nkp
				call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
				!read(unit,*) n1,kpoint(1:ndimension,n1)
                n2=int(ar(1))
                kpoint(1:ndimension,n2)=ar(2:ndimension+1)
                if(n1>ndimension+1) kpoint(ndimension+1,n2)=ar(ndimension+2)
				if(kpoint(ndimension+1,n2)<1e-6) kpoint(ndimension+1,n2)=DEFAULTSIZE
                if(kpoint(1,n2)<Minx) minx=kpoint(1,n2)
                if(kpoint(1,n2)>MaxX) maxX=kpoint(1,n2)
                if(kpoint(2,n2)<MinY) minY=kpoint(2,n2)
                if(kpoint(2,n2)>Minx) maxY=kpoint(2,n2)
                
            end do
            TOF1=.FALSE.
			do i=1,nkp-1
                DO J=I+1,NKP
                    T1=0.D0
                    DO K=1,NDIMENSION
                        T1=T1+(KPOINT(K,J)-KPOINT(K,I))**2
                    ENDDO
                    IF(ABS(T1)<1E-7) THEN
                        PRINT *, "THE POINTS OF I AND J ARE IDENTICAL. I,J=",I,J
                        TOF1=.TRUE.
                    ENDIF
                ENDDO
            end do            
            IF(TOF1) STOP "ERROR STOP IN KPOINT READ IN."
        case('geoline')
            print *, 'Reading geological LINE data...'
            do i=1,pro_num
                select case(property(i).name)
					case('num')
						ngeoline=int(property(i).value)                                            
					case default
						call Err_msg(property(i).name)
				end select    
            end do
            allocate(geoline(ngeoline))
            do i=1,ngeoline
                call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
                n2=int(ar(1))
                geoline(n2).mat=int(ar(2))
                n3=int(ar(3))
                geoline(n2).npoint=n3
                allocate(geoline(n2).point(n3))
                
                geoline(n2).point(1:n3)=int(ar(4:3+n3))
                do j=2,n3
                    if(kpoint(1,geoline(n2).point(j-1))>kpoint(1,geoline(n2).point(j))) then
                        print *, "x(j-1)>x(j) in geoline(i). j,i=",j,n2
                        stop
                    endif
                enddo
                if(nset==1) geoline(n2).title=set(1)
            enddo
		CASE('waterlevel')
			print *, 'Reading WALTERLEVEL LINE data...'
            N2=0
            do i=1,pro_num
                select case(property(i).name)
					case('num')
						WATERLEVEL.NPOINT=int(property(i).value)
                    case('fmt','format')
                        WATERLEVEL.FMT=int(property(i).value)
                    case('var')
                        WATERLEVEL.VAR=int(property(i).value)
					case default
						call Err_msg(property(i).name)
				end select    
            end do
            allocate(WATERLEVEL.POINT(WATERLEVEL.NPOINT),WATERLEVEL.H(2,WATERLEVEL.NPOINT))
            IF(WATERLEVEL.FMT==0) THEN
			    call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
			    WATERLEVEL.POINT=AR(1:n1)
                WATERLEVEL.H(1,:)=KPOINT(WATERLEVEL.VAR,WATERLEVEL.POINT)
                WATERLEVEL.H(2,:)=KPOINT(NDIMENSION,WATERLEVEL.POINT)
            ELSE
                DO J=1,WATERLEVEL.NPOINT
                    call skipcomment(unit)
                    READ(UNIT,*) WATERLEVEL.H(:,J)
                ENDDO
            ENDIF
			
        case('right turn point','rtpoint','rtp')
			print *, 'Reading RIGHT TURN POINTS ...'
            do i=1,pro_num
                select case(property(i).name)
					case('num')
						NRTPOINT=int(property(i).value)                                            
					case default
						call Err_msg(property(i).name)
				end select    
            end do
			IF(NRTPOINT>0) THEN
				ALLOCATE(RTPOINT(NRTPOINT))
				n2=0
				DO WHILE(n2<nrtpoint)
					call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
					rtpoint(n2+1:n2+n1)=int(ar(1:n1))
					n2=n2+n1
				ENDDO
			ENDIF
		case('hbeam','pile') !for retaining structure
			print *, 'Reading beam/pile(Retaining Structure) data...'
			do i=1, pro_num
				select case(property(i).name)
					case('num')
						npile=int(property(i).value)
					!case('datapacking','dp')
					!	datapacking=int(property(i).value)
					!case('dimension','d')
					!	ndimension=int(property(i).value)
					case('isplot','hbeam')
						ishbeam=int(property(i).value)
						isExca2D=2
					case default
						call Err_msg(property(i).name)
				end select
			end do
			allocate(pile(npile))
			do i=1,npile				
				call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
				pile(i).nseg=int(ar(1))
				if(n1>1) pile(i).system=int(ar(2))
				allocate(pile(i).kpoint(pile(i).nseg+1),pile(i).mat(pile(i).nseg))
				call skipcomment(unit)
				read(unit,*) pile(i).kpoint
				call skipcomment(unit)
				read(unit,*) pile(i).mat
			end do			
		case('strut') !for retaining structure
			print *, 'Reading STRUT(Retaining Structure) data...'
			n2=0
			n3=0
			n4=0
            do i=1, pro_num
				select case(property(i).name)
					case('num')
						n3=int(property(i).value)
					case('isbar')
						n2=int(property(i).value)
					case default
						call Err_msg(property(i).name)
				end select
			end do
			call enlarge_strut(strut,nstrut,n3,n4)
			!allocate(strut(nstrut))
			
			do i=1,n3
				call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
				IF(N1==4.AND.INT(AR(4))>0) then
					N2=1
				else
					n2=0
				endif
				if(n2==0) then
					strut(n4-1+i).z(1)=int(ar(1))
					strut(n4-1+i).mat=int(ar(2))
					strut(n4-1+i).sf=int(ar(3))
					!if(n1>3) strut(n4-1+i).preLoad=int(ar(4))
					!if(n1>4) strut(n4-1+i).preDis=int(ar(5))
				else
					if (excelformat==1) then
						strut(n4-1+i).z(1)=int(ar(1))
						strut(n4-1+i).z(2)=int(ar(4))
						strut(n4-1+i).mat=int(ar(2))
						strut(n4-1+i).sf=int(ar(3))
					else
						strut(n4-1+i).z(1)=int(ar(1))
						strut(n4-1+i).z(2)=int(ar(2))
						strut(n4-1+i).mat=int(ar(3))
						strut(n4-1+i).sf=int(ar(4))
					endif
					!if(n1>4) strut(n4-1+i).preLoad=int(ar(5))
					!if(n1>5) strut(n4-1+i).preDis=int(ar(6))				
				endif
				
				strut(n4-1+i).isbar=n2
				
			end do				

		case('action') !for retaining structure
			print *, 'Reading ACTION(Retaining Structure) data...'
			do i=1, pro_num
				select case(property(i).name)
					case('num')
						naction=int(property(i).value)
					!case('datapacking','dp')
					!	datapacking=int(property(i).value)
					!case('dimension','d')
					!	ndimension=int(property(i).value)
					case default
						call Err_msg(property(i).name)
				end select
			end do
			allocate(action(naction))
			do i=1,naction
				call skipcomment(unit)
				read(unit,'(A64)') action(i).title
				n2=0
				n3=0
				n4=0
				call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
				action(i).nkp=int(ar(1))
				action(i).type=int(ar(2))
				action(i).dof=int(ar(3))
				action(i).ndim=int(ar(4))
				if(n1>4) action(i).sf=int(ar(5))
				if(n1>5) n3=int(ar(6))
				if(n1>6) n4=int(ar(7))
				allocate(action(i).kpoint(action(i).nkp),action(i).value(action(i).nkp), &
						 action(i).vsf(action(i).nkp),action(i).exvalue(action(i).nkp,2),&
						 action(i).node(action(i).nkp))
				call skipcomment(unit)
				read(unit,*) action(i).kpoint
				call skipcomment(unit)
				read(unit,*) action(i).value
                MAXACTION=MAX(ABS(MAXVAL(ACTION(I).value)),ABS(MINVAL(ACTION(I).VALUE)),MAXACTION)                
				if(n3==1) then
					call skipcomment(unit)
					read(unit,*) action(i).vsf
				else
					action(i).vsf=0
				endif
				if(n4==1) then
					call skipcomment(unit)
					read(unit,*) action(i).exvalue
				else
					action(i).exvalue(:,1)=-1E20
					action(i).exvalue(:,2)=1E20
				endif				
				
			end do
			
			
		case('element')
			print *, 'Reading ELEMENT data...'
			matid1=0
            system1=0
			name1=""
			sf1=0
            n2=-1
			do i=1,pro_num
				select case(property(i).name)
					case('num')
						enum1=int(property(i).value)
					case('set')
						set1=int(property(i).value)
						isset=.true.
                    CASE('coupleset')
                        n2=int(property(i).value)
					case('et','type')
						et1=int(property(i).value)
!						if(et1>maxet) maxet=et1
!						if(et1<minet) minet=et1
						!according to the element type, return the element node number,and the dofs number
						call ettonnum(et1,nnum1,ndof1,ngp1,nd1,stype,ec1)
					case('material')
						material1=int(property(i).value)
					case('mat','matid','material id')
						matid1=int(property(i).value)
					case('system') !local coordinate
						system1=int(property(i).value)
					case('title','name')
						name1=property(i).cvalue
					CASE('sf','step function')
						sf1=int(property(i).value)
					case default
						call Err_msg(property(i).name)						
				end select
			end do
			
			!if(matid1==0) matid1=material1 !If there is only one set of such material of this type in this model 
			neset=neset+1
            
            if(isset) then
				isset=.false.
			else
				set1=neset
			end if
            
            SET1=ESET_GETFREEID(SET1)
            
            esetid(neset)=set1
            
            if(n2<0) n2=set1
			
			if(material(matid1).type==0) then
				select case(ec1)
					CASE(SPG,SPG2D,CAX_SPG)
						MATERIAL(MATID1).TYPE=linear_spg
					CASE DEFAULT
						MATERIAL(MATID1).TYPE=ELASTIC
				end select
			endif
            
			allocate(element1(enum1))
			
			
			
			do i=1,enum1
				element1(i).nnum=nnum1
				allocate(element1(i).node(nnum1))
				select case(et1)
					!case(beam) !beam element, a local system must be input.
					!	call skipcomment(unit)
					!	read(unit,*) element1(i).node,element1(i).system
					case(dkt3,shell3,shell3_KJB) !.h  is the thickness of the element.
						call skipcomment(unit)
						read(unit,*) element1(i).node,element1(i).PROPERTY(3)
					case default
						call skipcomment(unit)
						read(unit,*) element1(i).node
				end select
				element1(i).id=i
				element1(i).et=et1
				element1(i).set=set1                
				element1(i).mat=matid1
				element1(i).mattype=material1
				element1(i).ndof=ndof1
				element1(i).ngp=ngp1
				element1(i).nd=nd1
				element1(i).ec=ec1
				element1(i).sf=sf1
				if(et1==beam) element1(i).system=system1
			end do
			!eset(set1).num=set1
			eset(set1).stype=stype
			eset(set1).grouptitle=name1
			eset(set1).et=et1
			eset(set1).ec=ec1
            eset(set1).system=system1
			eset(set1).enums=enum+1
            eset(set1).coupleset=n2
            eset(set1).sf=sf1
			allocate(element2(enum+enum1))
			element2(1:enum)=element(1:enum)
			element2(enum+1:enum+enum1)=element1(1:enum1)
			if(allocated(element))	deallocate(element)
			deallocate(element1)
			enum=enum+enum1
			eset(set1).enume=enum
			allocate(element(enum))
			element=element2
			deallocate(element2)
        case('rcd')
            print *, 'Reading Rigid Connected Dof data...'
            do i=1,pro_num
				select case(property(i).name)
					case('num')
						enum1=int(property(i).value)
					case('dof')
						material1=int(property(i).value)
					case default
						call Err_msg(property(i).name)						
				end select
			end do
		case('material')
			print *, 'Reading MATERIAL data...'
			n1=0
			matid1=0
			j=0
			name1=""
            n3=0
			n2=1
			do i=1,pro_num
				select case(property(i).name)
					CASE('num')
						n2=int(property(i).value)
					case('type')
						j=int(property(i).value)
						!material(j).id=j
					case('isff')					
						if(int(property(i).value)==YES)  n1=1
					case('matid')
						matid1=int(property(i).value)
					case('issf')
						n3=int(property(i).value)
					case('name','title')
						name1=property(i).cvalue
					case default
						call Err_msg(property(i).name)
				end select
			end do
			
			
			do i=1,n2
				if(excelformat==1) then
					call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
					matid1=int(ar(1));j=int(ar(2));n1=int(ar(3));n3=int(ar(4))
					name1=set(1)
				endif
				
				if(matid1==0) matid1=j !If there is only one set of such material of this type in this model 
				
				if(j/=0) material(matid1).type=j
				material(matid1).name=name1
				if(n1==1) material(matid1).isff=.true.
				call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
				!n1=size(material(matid1).property)
				material(matid1).property(1:n1)=ar(1:n1)	
				select case(material(matid1).type)
					case(mises)
						material(matid1).weight=material(matid1).property(4)
					case(mc)
						material(matid1).weight=material(matid1).property(6)
                        if(abs(material(matid1).property(4)-material(matid1).property(5))>1e-7) solver_control.issym=.false.
						if(n1<=6) then
							if(abs(material(matid1).property(4))>1e-7) then
								material(matid1).property(23)=0.05*material(matid1).property(3) &
																/dtan(material(matid1).property(4)/180.*PI())
							endif
							if(abs(material(matid1).property(5))>1e-7) then
								material(matid1).property(24)=0.05*material(matid1).property(3) &
															/dtan(material(matid1).property(5)/180.*PI())
							endif							
							material(matid1).property(8)=28.d0
						elseif(n1<=7) then                            
							material(matid1).property(8)=28.d0
						endif
                        
                        !FOR CLAUSEN MC PARAMETERS
                        T2=DSIN(material(matid1).property(4)/180.*PI())
                        MATERIAL(MATID1).PROPERTY(19)=(1+T2)/(1-T2) !k
                        MATERIAL(MATID1).PROPERTY(20)=2*MATERIAL(MATID1).PROPERTY(3)*sqrt(MATERIAL(MATID1).PROPERTY(19))
                        T2=DSIN(material(matid1).property(5)/180.*PI()) 
                        MATERIAL(MATID1).PROPERTY(21)=(1+T2)/(1-T2) !M
						!INVERSE(d)
						ALLOCATE(MATERIAL(MATID1).DINV(NDIMENSION*2,NDIMENSION*2)) 
						MATERIAL(MATID1).DINV=DINV(MATERIAL(MATID1).PROPERTY(1),MATERIAL(MATID1).PROPERTY(2),NDIMENSION*2)
                        
                        IF (ABS(material(matid1).property(8)-30.D0)<1E-14.AND.ABS(material(matid1).property(7))<1E-14 ) THEN
                            material(matid1).property(22)=0.D0 !NO ROUNDED AND GRIFFITHS ALGORITHM IS USED
                        ELSE
						    !ROUNDOFF PARAMETERS
                            material(matid1).property(22)=1.D0
							if(n1>6) then
								material(matid1).property(23:24)=material(matid1).property(7)
							endif
						    T1=material(matid1).property(8)/180.*PI()
						    T2=DSIN(material(matid1).property(4)/180.*PI())
							!A1,A2*SIN(PHI),B1,B2*SIN(PHI)
						    material(matid1).property(29)=dcos(T1)/3.d0*(3.D0+DTAN(T1)*DTAN(3.*T1))
						    material(matid1).property(30)=dcos(T1)/3.d0/SQRT(3.0)*(DTAN(3.*T1)-3.*DTAN(T1))*T2
						    material(matid1).property(31)=DSIN(T1)/(3.d0*dcos(3.*T1))
						    material(matid1).property(32)=T2*DCOS(T1)/(3.d0*SQRT(3.0)*dcos(3.*T1))
						
						    T2=DSIN(material(matid1).property(5)/180.*PI())
						    material(matid1).property(25)=dcos(T1)/3.d0*(3.D0+DTAN(T1)*DTAN(3.*T1))
						    material(matid1).property(26)=dcos(T1)/3.d0/SQRT(3.0)*(DTAN(3.*T1)-3.*DTAN(T1))*T2
						    material(matid1).property(27)=DSIN(T1)/(3.d0*dcos(3.*T1))
						    material(matid1).property(28)=T2*DCOS(T1)/(3.d0*SQRT(3.0)*dcos(3.*T1))						
						ENDIF
					case(eip_bar)
						if(n1<=4) then
							material(matid1).property(5)=-1.0D20 !最大轴向压力
							material(matid1).property(6)=1.0D20	 !最大的轴向拉力
						end if
				end select
				
				if(material(matid1).isff) then
					call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
					material(matid1).ff1d(1:n1)=int(ar(1:n1))
				end if
				material(matid1).sf=0
				if(n3==1) then
					call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
					material(matid1).sf(1:n1)=int(ar(1:n1))
				end if
			
			enddo
		case('load')
			print *,'Reading LOAD data...'
			n2=0
			n3=0
			n4=0
			n5=0
			do i=1,pro_num
				select case(property(i).name)
					case('num')
						nbf1=int(property(i).value)
					case('ssp_onepile')
						n2=int(property(i).value)
					case('spg_isdual')
						n3=int(property(i).value)
					case('sf','stepfunction','stepfunc')					
						n4=int(property(i).value)
					case('isinc')
						n5=int(property(i).value)
					case default
						call Err_msg(property(i).name)
				end select
			end do
			allocate(bf1(bl_num+nbf1))
			do i=bl_num+1,bl_num+nbf1
				call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
				bf1(i).node=int(ar(1))
				bf1(i).dof=int(ar(2))
				bf1(i).value=ar(3)
				bf1(i).isincrement=n5
				bf1(i).sf=n4								
				if(n1>=4) bf1(i).sf=int(ar(4))
				bf1(i).isdual=n3
				if(n1>=5) bf1(i).isdual=int(ar(5)) !同是也可能是出溢边界
				bf1(i).ssp_onepile=n2
				if(n1>=6) bf1(i).ssp_onepile=int(ar(6))
                IF(N1>=7) bf1(i).isincrement=int(ar(7))
                OUTVAR(90+BF1(I).DOF).VALUE=90+BF1(I).DOF   
			end do
			if(bl_num>0) then
				bf1(1:bl_num)=bc_load(1:bl_num)
				deallocate(bc_load)
			end if			
			allocate(bc_load(bl_num+nbf1))
			bc_load=bf1
			bl_num=bl_num+nbf1
			deallocate(bf1)
		case('ncf','normal contact force')
			print *,'Reading normal contact force data...'
			n2=0
			n3=0
			n4=0
			do i=1,pro_num
				select case(property(i).name)
					case('num')
						nbf1=int(property(i).value)
					case('ssp_onepile')
						n2=int(property(i).value)
					case('spg_isdual')
						n3=int(property(i).value)
					case('sf','stepfunction','stepfunc')
						n4=int(property(i).value)
					case default
						call Err_msg(property(i).name)
				end select
			end do
			allocate(bf1(ncfn+nbf1))
			do i=ncfn+1,ncfn+nbf1
				call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
				bf1(i).node=int(ar(1))
				bf1(i).dof=int(ar(2))
				bf1(i).value=ar(3)
				
				bf1(i).sf=n4								
				if(n1>=4) bf1(i).sf=int(ar(4))
				bf1(i).isdual=n3
				if(n1>=5) bf1(i).isdual=int(ar(5)) !同是也可能是出溢边界
				bf1(i).ssp_onepile=n2
				if(n1>=6) bf1(i).ssp_onepile=int(ar(6))

			end do
			if(ncfn>0) then
				bf1(1:ncfn)=cfn(1:ncfn)
				deallocate(cfn)
			end if			
			allocate(cfn(ncfn+nbf1))
			cfn=bf1
			ncfn=ncfn+nbf1
			deallocate(bf1)	
			
		case('bf','body force','pressure','elt_load')
			print *, 'Reading BODY FORCE data...'
			do i=1,pro_num
				select case(property(i).name)
					case('num')
						nbf1=int(property(i).value)					
					case default
						call Err_msg(property(i).name)
				end select
			end do			
			allocate(bf1(nbf1))
			bf1.value=0.0			
			n1=0
			n2=0			
			do while(n2<nbf1)
				!the structure for each data line must be kept the same.
				call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
				do i=n2+1,n2+n1-3
					bf1(i).node=int(ar(i-n2))
					bf1(i).dof=int(ar(n1-2))
					bf1(i).value=bf1(i).value+ar(n1-1) 
					!Attention. the unit is force/volume if the element is an planar element, then unit should 
					!be Force/volume*(element thickness).
					bf1(i).sf=int(ar(n1))
				end do
				n2=n2+n1-3
				do i=1,nset
					do j=0,nueset
						if(index(ueset(j).name,set(i))>0)  exit
					end do
					if(j==nueset+1) then
						print *, 'No such element set. '//trim(set(i))
						stop
					end if
					bf1(n2+1:n2+ueset(j).enum).node= &
							ueset(j).element(1:ueset(j).enum)
					bf1(n2+1:n2+ueset(j).enum).dof=int(ar(n1-2))					
					bf1(n2+1:n2+ueset(j).enum).value= & 
					bf1(n2+1:n2+ueset(j).enum).value+ar(n1-1)
					bf1(n2+1:n2+ueset(j).enum).sf=int(ar(n1))
					n2=n2+ueset(j).enum						
				end do
								
			end do			
			if(bfnum>0) then
				allocate(bf2(bfnum))
				bf2(1:bfnum)=bf(1:bfnum)
				deallocate(bf)
			end if			
			allocate(bf(bfnum+nbf1))
			if(bfnum>0) bf(1:bfnum)=bf2(1:bfnum)
			bf(bfnum+1:bfnum+nbf1)=bf1
			bfnum=bfnum+nbf1
			deallocate(bf1)
			if(allocated(bf2)) deallocate(bf2)
		case('bc','boundary condition')
			print *,'Reading BOUNDARY CONDITION data...'
			n2=0
			n3=0
			n4=0
			n5=0
			do i=1,pro_num
				select case(property(i).name)
					case('num')
						nbf1=int(property(i).value)
					case('ssp_onepile')
						n2=int(property(i).value)
					case('spg_isdual')
						n3=int(property(i).value)
					case('sf','stepfunction','stepfunc')
						n4=int(property(i).value)
					case('isinc')
						n5=int(property(i).value)	
					case default
						call Err_msg(property(i).name)
				end select
			end do
			allocate(bf1(bd_num+nbf1))
			do i=bd_num+1,bd_num+nbf1
				call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
				bf1(i).node=int(ar(1))
				bf1(i).dof=int(ar(2))
				bf1(i).value=ar(3)
				bf1(i).isincrement=n5
				bf1(i).sf=n4								
				if(n1>=4) bf1(i).sf=int(ar(4))
				bf1(i).isdual=n3
				if(n1>=5) bf1(i).isdual=int(ar(5)) !同是也可能是出溢边界
				bf1(i).ssp_onepile=n2
				if(n1>=6) bf1(i).ssp_onepile=int(ar(6))
                if(n1>=7) bf1(i).isincrement=int(ar(7))
                OUTVAR(90+BF1(I).DOF).VALUE=90+BF1(I).DOF                
			end do
			if(bd_num>0) then
				bf1(1:bd_num)=bc_disp(1:bd_num)
				deallocate(bc_disp)
			end if			
			allocate(bc_disp(bd_num+nbf1))
			bc_disp=bf1
			bd_num=bd_num+nbf1
			deallocate(bf1)	
			
            
            
        case('hinge','freedof')
            print *, 'Reading HINGE/FREEDOF data...'
  			do i=1,pro_num
				select case(property(i).name)
					case('num')
						nfreedof=int(property(i).value)
					case default
						call Err_msg(property(i).name)
				end select
            end do 
            allocate(freedof(nfreedof))
            do i=1,nfreedof
                call skipcomment(unit)
				read(unit,*) freedof(i).element,freedof(i).node,freedof(i).dof
            enddo
		case('seepage face')
			print *,'Reading Nodes In SEEPAGEFACE data...'
			n3=0
			n4=0
			do i=1,pro_num
				select case(property(i).name)
					case('num')
						n4=int(property(i).value)
					case('step function','sf')
						n3=int(property(i).value)
					case default
						call Err_msg(property(i).name)
				end select
			end do
			allocate(Nseep1(n4+Numnseep))
			n2=0
			do while(n2<n4)
				call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
				Nseep1(n2+1:n2+n1).node=int(ar(1:n1))
				Nseep1(n2+1:n2+n1).sf=n3
				n2=n1+n2				
			end do
			if(Numnseep>0)	then
				Nseep1(n4+1:n4+Numnseep)=Nseep
				deallocate(Nseep)
			end if
			Numnseep=n4+Numnseep
			allocate(Nseep(Numnseep))
			Nseep=Nseep1
			deallocate(Nseep1)
			
			Nseep.dof=4
			Nseep.isdead=0

			Nseep.value=Node(Nseep.node).coord(ndimension)
!			node(Nseep.node).property=1
		case('datapoint')
			print *, 'Reading DATAPOINT data...'
			do i=1,pro_num
				select case(property(i).name)
					case('num')
						ndatapoint=int(property(i).value)
					case default
						call Err_msg(property(i).name)
				end select
			end do
			allocate(datapoint(ndatapoint))
			do i=1,ndatapoint
				call skipcomment(unit)
				call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
				datapoint(i).nnode=int(ar(1))
				if(n1>1) datapoint(i).issumq=1
				allocate(datapoint(i).node(datapoint(i).nnode))
				call skipcomment(unit)
				read(unit,*) datapoint(i).node
			end do
		case('soilprofile')
			print *, "Reading SOILPROFILE Data..."
			n3=0
			n4=0
			do i=1,pro_num
				select case(property(i).name)
					case('num')
						nsoilprofile=int(property(i).value)
					case('spm','soilpressuremethod','spmethod')
						n3=int(property(i).value)
					case('kmethod')
						n4=int(property(i).value)
                    case('rf_epp')
                        solver_control.rf_epp=int(property(i).value)
                    case('rf_app')
                        solver_control.rf_app=int(property(i).value)
                    case('iniepp')
						solver_control.iniepp=int(property(i).value)
					case('soilspringmodel')
						material(-2:-1).type=int(property(i).value)
						
					case default
						call Err_msg(property(i).name)
				end select
			end do
			allocate(soilprofile(nsoilprofile))
			soilprofile.spm=n3
			soilprofile.kmethod=n4
            ISEXCA2D=1
			
			soilprofile(1).soilspringmodel=material(-1).type !!!all is the same
			
			do i=1,nsoilprofile
				call skipcomment(unit)
				read(unit,'(A64)') soilprofile(i).title
				call skipcomment(unit)
				call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
				!read(unit,*) soilprofile(i).nasoil,soilprofile(i).npsoil,soilprofile(i).beam,soilprofile(i).naction,SOILPROFILE(I).NSTRUT
				soilprofile(i).nasoil=int(ar(1));soilprofile(i).npsoil=int(ar(2));soilprofile(i).beam=int(ar(3))
				if(n1>3) soilprofile(i).naction=int(ar(4))
				if(n1>4) soilprofile(i).NSTRUT=int(ar(5))
				if(soilprofile(i).nasoil<0) then 
					soilprofile(i).aside=-1
					soilprofile(i).nasoil=-soilprofile(i).nasoil
				endif
				allocate(soilprofile(i).asoil(soilprofile(i).nasoil),soilprofile(i).psoil(soilprofile(i).npsoil), &
						 soilprofile(i).iaction(soilprofile(i).naction),soilprofile(i).istrut(soilprofile(i).nstrut))
				do j=1,soilprofile(i).nasoil
					call skipcomment(unit)
					call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
					!read(unit,*) soilprofile(i).asoil(j).z,soilprofile(i).asoil(j).mat,soilprofile(i).asoil(j).wpflag,soilprofile(i).asoil(j).sf
					soilprofile(i).asoil(j).z=int(ar(1:2))
					soilprofile(i).asoil(j).mat=int(ar(3))
					soilprofile(i).asoil(j).wpflag=int(ar(4))
					soilprofile(i).asoil(j).sf=int(ar(5))
					if(n1>5) soilprofile(i).asoil(j).pv=ar(6)
                    if(n1>6) soilprofile(i).asoil(j).soiltype=int(ar(7))
				enddo
				
				do j=1,soilprofile(i).npsoil
					call skipcomment(unit)
					call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
					!read(unit,*) soilprofile(i).psoil(j).z,soilprofile(i).psoil(j).mat,soilprofile(i).psoil(j).wpflag,soilprofile(i).psoil(j).sf
					soilprofile(i).psoil(j).z=int(ar(1:2))
					soilprofile(i).psoil(j).mat=int(ar(3))
					soilprofile(i).psoil(j).wpflag=int(ar(4))
					soilprofile(i).psoil(j).sf=int(ar(5))
					if(n1>5) soilprofile(i).psoil(j).pv=ar(6)
                    if(n1>6) soilprofile(i).psoil(j).soiltype=int(ar(7))
                enddo
                
				call skipcomment(unit)
				read(unit,*) soilprofile(i).awL,soilprofile(i).sf_awL,soilprofile(i).pwL,soilprofile(i).sf_pwL
				call skipcomment(unit)
				read(unit,*) soilprofile(i).aLoad,soilprofile(i).sf_aLoad,soilprofile(i).pLoad,soilprofile(i).sf_pLoad
				if(soilprofile(i).naction>0) then
					call skipcomment(unit)
					read(unit,*) soilprofile(i).iaction
				endif
				if(soilprofile(i).NSTRUT>0) then
					call skipcomment(unit)
					read(unit,*) soilprofile(i).istrut
                endif
                
			end do
		
			
		case('solvercontrol','solver','solver_control')
			print *, 'Reading SOLVER_CONTROL data'
            n1=0
			do i=1,pro_num
				select case(property(i).name)
					case('type') 
						solver_control.type=int(property(i).value)
						if(solver_control.type==ssa) isslope=1
					case('solver')
						solver_control.solver=int(property(i).value)
!					case('nincrement','ninc')
!						solver_control.nincrement=int(property(i).value)
!						allocate(solver_control.factor(solver_control.nincrement))
					case('tolerance','tol')
						solver_control.tolerance=property(i).value
					case('dtol')
						solver_control.disp_tol=property(i).value
					case('ftol')
						solver_control.force_tol=property(i).value
					case('niteration','nite','maxiter')
						solver_control.niteration=int(property(i).value)
					case('output')
						solver_control.output=int(property(i).value)
					case('symmetric','sys')
						if(int(property(i).value)==0) solver_control.issym=.false.
					case('datapaking')
						if(int(property(i).value)==BLOCK) solver_control.datapaking=.false.
					case('ismg')
						if(int(property(i).value)==YES) solver_control.ismg=.true.
					case('islaverify')
						if(int(property(i).value)==YES) solver_control.islaverify=.true.
					case('ispg')
						if(int(property(i).value)==YES) solver_control.ispg=.true.
					case('i2ncal','i2n')
						solver_control.i2ncal=int(property(i).value)
					case('i2nweight','i2nw')
						solver_control.i2nweight=int(property(i).value)                        
					case('bfgm','sim')
						solver_control.bfgm=int(property(i).value)                        
                    case('bfgm_spg')
                        solver_control.bfgm_spg=int(property(i).value)
                        n1=1
					case('isfc','force_criteria')
						if(int(property(i).value)==YES) then
							solver_control.isfc=.true.
						else
							solver_control.isfc=.false.
						end if						
!					case('para_spg')
!						solver_control.para_spg=int(property(i).value)
					case('isfu')
						if(int(property(i).value)==YES) then
							solver_control.isfu=.true.
						else
							solver_control.isfu=.false.
						end if
!					case('steady')
!						if(int(property(i).value)==YES) then
!							solver_control.issteady=.true.
!						else
!							solver_control.issteady=.false.
!						end if
					case('mkl')
						if(int(property(i).value)==YES) then
							solver_control.ismkl=.true.
						else
							solver_control.ismkl=.false.
                        end if
                    case('ls','linesearch')
						if(int(property(i).value)==YES) then
							solver_control.isls=.true.
						else
							solver_control.isls=.false.
                        end if
                    CASE('acc') 
						solver_control.isacc=int(property(i).value)
                    case('mur')
                        solver_control.mur=int(property(i).value)
					case('barfamilyscale')
						solver_control.barfamilyscale=int(property(i).value)
                    case('nopopup')
                        solver_control.nopopup=int(property(i).value)
                    case('isparasys')
                        solver_control.isparasys=int(property(i).value)
                    case('caseid')
                        solver_control.caseid=int(property(i).value)
                    case('slidedirection')
                        solver_control.slidedirection=int(property(i).value)
                    case('slope_kscale')
                        solver_control.slope_kscale=property(i).value
                    case('slope_kbase')
                        solver_control.slope_kbase=property(i).value
                    case('slope_istensioncrack')
                        solver_control.slope_isTensionCrack=int(property(i).value)
					!case('ispostcal')
					!	solver_control.ispostcal=int(property(i).value)
					case default
						call Err_msg(property(i).name)
				end select
            end do

			if(solver_control.BFGM==INISTRESS) then
				solver_control.solver=INISTIFF
			endif	
			if(solver_control.BFGM==CONTINUUM.OR.solver_control.BFGM==CONSISTENT) then
				solver_control.solver=N_R
            endif	
            if(n1==0) solver_control.bfgm_spg=solver_control.bfgm
!			if(associated(solver_control.factor)) then
!				read(unit,*)   solver_control.factor
!			else
!				!only one increment, set the factor=1.0
!				allocate(solver_control.factor(1))
!				solver_control.factor(1)=1.0
!			end if
        case('bfgm_step')
            print *,'Reading BFGM INFO FOR EACH STEP...'
            do i=1,pro_num
				select case(property(i).name)
                case('step','num')
                    n1=int(property(i).value)
                case default
                    call Err_msg(property(i).name)                    
                endselect
            enddo
            allocate(bfgm_step(n1))
            call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
            bfgm_step=int(ar(1:n1))
        CASE('slopeparameter')
            print *,'Reading SLOPE PARAMETER data...'
            n1=0
            do i=1,pro_num
				select case(property(i).name)
                case('slopemthod','sm')
                    slopeparameter.slopemethod=int(property(i).value)
                case('optimizationmthod','om') 
                    slopeparameter.optmethod=int(property(i).value)
                case('slipshape','ss')
                    slopeparameter.slipshape=int(property(i).value)
                case('slicewdith','sw')
                    slopeparameter.slicewidth=property(i).value
                CASE('xmin_mc')
                    slopeparameter.xmin_mc=property(i).value
                CASE('xmax_mc')
                    slopeparameter.xmax_mc=property(i).value
                case('downwardzone')
                    slopeparameter.ISYDWZ=.TRUE.
                    slopeparameter.ydownwardzone=property(i).value
                case('toezone')
                    IF(INT(property(i).value)/=0) N1=1
                endselect
            enddo
            IF(N1/=0) THEN
                call skipcomment(unit)
                READ(UNIT,*) AR(1:4)
                T1=AR(3)-AR(1) !X2-X1
                IF(ABS(T1)>1.E-7) THEN
                    slopeparameter.TOEZONE(1)=(AR(4)-AR(2))/T1
                    slopeparameter.TOEZONE(2)=-1.D0
                    slopeparameter.TOEZONE(3)=AR(2)-slopeparameter.TOEZONE(1)*AR(1)
                ELSE
                    slopeparameter.TOEZONE(1)=1.D0
                    slopeparameter.TOEZONE(2)=0.D0
                    slopeparameter.TOEZONE(3)=AR(1)
                ENDIF

            ENDIF
		case('ueset')
			print *, 'Reading USER DEFINED ELEMENT SET data'
			if(nueset==1) then !intialize ueset(0)
				ueset(0).enum=enum
				ueset(0).name='all'
				allocate(ueset(0).element(ueset(0).enum))
				do i=1,ueset(0).enum
					ueset(0).element(i)=i
				end do
			end if
			do i=1,pro_num
				select case(property(i).name)
					case('num')
						enum1=int(property(i).value)
					case('name')
						name1=property(i).cvalue
					case default
						call Err_msg(property(i).name)
				end select
			end do
			nueset=nueset+1
			allocate(ueset(nueset).element(enum1))
			n1=0
			n2=0
			ueset(nueset).name=name1
			ueset(nueset).enum=enum1
			do while(n2<enum1)
				call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
				ueset(nueset).element(n2+1:n2+n1)=int(ar(1:n1))
				n2=n2+n1
				do i=1,nset
					do j=0,nueset-1
						if(index(ueset(j).name,set(i))>0)  exit
					end do
					if(j==nueset) then
						print *, 'No such element set. '//trim(set(i))
						stop
					end if
					ueset(nueset).element(n2+1:n2+ueset(j).enum)= &
							ueset(j).element(1:ueset(j).enum)
				end do
				n2=n2+ueset(j).enum
			end do
		case('unset')
			print *, 'Reading USER DEFINED NODE SET data'
			if(nunset==1) then !initialize the unset(0)
				unset(0).nnum=nnum
				unset(0).name='all'
				allocate(unset(0).node(unset(0).nnum))
				do i=1,unset(0).nnum
					unset(0).node(i)=i
				end do
			end if			
			do i=1,pro_num
				select case(property(i).name)
					case('num')
						nnum1=int(property(i).value)
					case('name')
						name1=property(i).cvalue					
					case default
						call Err_msg(property(i).name)
				end select
			end do
			nunset=nunset+1
			allocate(unset(nunset).node(nnum1))
			n2=0
			n1=0
			unset(nunset).name=name1
			unset(nunset).nnum=nnum1
			do while(n2<nnum1)
				call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
				unset(nunset).node(n2+1:n2+n1)=ar(1:n1)
				n2=n2+n1
				do i=1,nset
					do j=0,nunset-1
						if(index(unset(j).name,set(i))>0)  exit
					end do
					unset(nunset).node(n2+1:n2+unset(j).nnum)= &
							unset(j).node(1:unset(j).nnum)
				end do
				n2=n2+unset(j).nnum
			end do
		case('sf','step_function','step function')
			print *, 'Reading STEP FUNCTION data...'
            n2=1
			do i=1,pro_num
				select case(property(i).name)
					case('num')
						nsf=int(property(i).value)
					case('step','nstep')
						nstep=int(property(i).value)
                    case('base')
                        n2=int(property(i).value)
					case default
						call Err_msg(property(i).name)
				end select
			end do		
			allocate(sf(0:nsf))
            sf(1:nsf).base=n2
			do i=1,nsf
				allocate(sf(i).factor(0:nstep))
				sf(i).factor=0.0
				call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
				if(excelformat==0) then
					sf(i).factor(n2:nstep)=ar(1:n1)                    
				else
					sf(int(ar(1))).factor(n2:nstep)=ar(2:n1)
				endif
                IF(NSET==1) SF(I).TITLE=SET(1)
			end do
			nsf=nsf+1 !!
			allocate(sf(0).factor(0:nstep))
			sf(0).factor=1.0d0
			sf(0).factor(0)=0.0d0
		case('time step','ts')
			print *, 'Reading Incremental Time for each step...'
			do i=1,pro_num
				select case(property(i).name)
					case('num')
						nts=int(property(i).value)
					case default
						call Err_msg(property(i).name)
				end select
			end do		
			allocate(timestep(0:nts))
			!usually, nts=nstep,it may be nts=nstep+1
			do i=1,nts
				call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
				j=int(ar(1))
				timestep(j).nsubts=int(ar(2))
				allocate(timestep(j).subts(timestep(j).nsubts))
				timestep(j).subts(1:timestep(j).nsubts)=ar(3:n1)
			end do
			if(.not.allocated(timestep(0).subts)) then
				timestep(0).nsubts=0
				allocate(timestep(0).subts(1))
				timestep(0).subts(1)=0.d0
			end if
		case('stepinfo')
			print *, 'Reading Stepinfo data...'
			do i=1,pro_num
				select case(property(i).name)
					case('num')
						nstepinfo=int(property(i).value)
					case default
						call Err_msg(property(i).name)
				end select
			end do		
			allocate(stepinfo(0:nstepinfo))
			!usually, nstepinfo=nstep
			do i=1,nstepinfo
				call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
				j=int(ar(1))
				stepinfo(j).matherstep=int(ar(2))
				if(int(ar(3))==1) then
					stepinfo(j).issteady=.true.
				else
					stepinfo(j).issteady=.false.
				end if
				if(n1>=4) then
					if(int(ar(4))==1) stepinfo(j).bctype=ramp 
					if(int(ar(4))==-1) stepinfo(j).bctype=Reramp
					if(int(ar(4))==2) stepinfo(j).bctype=step
				end if
				if(n1>=5) then
					if(int(ar(5))==1) stepinfo(j).loadtype=ramp 
					if(int(ar(5))==-1) stepinfo(j).loadtype=Reramp
					if(int(ar(5))==2) stepinfo(j).loadtype=step				
				end if				
			end do
		
        case('wsp','watersurfaceprofile')
            print *, 'Reading WaterSurfaceProfile data...'
            do i=1,pro_num
				select case(property(i).name)
                    CASE('num')
                        nHJump=int(property(i).value)
					!case('caltype')
						!HJump.caltype=int(property(i).value)
                    case('cp')
                        HJump_CP=int(property(i).value)
					case('method')
						n1=int(property(i).value)
					case default
						call Err_msg(property(i).name)
				end select
            end do
			allocate(hjump(nhjump))
			hjump.method=n1
            do i=1,nhjump
				call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
				HJump(i).nseg=int(ar(1))
				HJump(i).nnode=int(ar(2))
				HJump(i).Q=ar(3)				
				HJump(i).B=ar(4)
				HJump(i).UYBC=ar(5)
				HJump(i).DYBC=ar(6)
				if(n1>6) HJump(i).Q_SF=int(ar(7))
				if(n1>7) HJump(i).UYBC_SF=int(ar(8))
				if(n1>8) HJump(i).DYBC_SF=int(ar(9))
				if(n1>9) HJump(i).Caltype=int(ar(10))				
				if(n1>10) HJump(i).g=ar(11)
				if(n1>11) HJump(i).kn=ar(12)
				if(HJump(i).UYBC<0) HJump(i).UYBC_Type=int(HJump(i).UYBC)
				if(HJump(i).DYBC<0) HJump(i).DYBC_Type=int(HJump(i).DYBC)
				allocate(HJump(i).segment(HJump(i).nseg),HJump(i).node(HJump(i).nnode),HJump(i).xy(11,HJump(i).nnode), &
						 HJump(i).HJump(11,2),HJump(i).JTinfo(5,HJump(i).nnode))
				if(HJUMP_CP/=1) allocate(HJump(i).BC_node(HJump(i).nnode))
				HJump(i).xy=0.0
				HJump(i).Hjump=0.d0
				do j=1,HJump(i).nseg
					 call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
					 HJump(i).segment(j).unode=int(ar(1))
					 HJump(i).segment(j).dnode=int(ar(2))
					 HJump(i).segment(j).n=ar(3)
					 HJump(i).segment(j).So=ar(4)
					 if(nset>0) HJump(i).segment(j).profileshape(1)=adjustL(set(1)(1:2))
					 if(nset>1) HJump(i).segment(j).profileshape(2)=adjustL(set(2)(1:2))
				end do    
						   
				call skipcomment(unit)
				read(unit,*) HJump(i).node
				do k=1,HJump(i).nnode
					HJump(i).xy(1:2,k)=node(HJump(i).node(k)).coord(1:2)                    
				end do
			end do

		case('initial value','iv')
			print *,'Reading Initial Value data...'
			do i=1,pro_num
				select case(property(i).name)
					case('num')
						NiniV=int(property(i).value)
					case default
						call Err_msg(property(i).name)
				end select
			end do
			allocate(inivalue(NiniV))
			do i=1,NiniV
				call strtoint(unit,ar,nmax,n1,n_toread,set,maxset,nset)
				inivalue(i).node=int(ar(1))
				inivalue(i).dof=int(ar(2))
				inivalue(i).value=ar(3)
				if(n1==4) inivalue(i).sf=int(ar(4))
			end do
		
		case('slave_master_node_pair','smnp')
			print *, 'Reading SLAVE-MASTER NODE PAIR'
			do i=1,pro_num
				select case(property(i).name)
					case('num')
						nSMNP=int(property(i).value)
					case default
						call Err_msg(property(i).name)
				end select
			end do
			allocate(SMNP(nSMNP))
			do i=1,nSMNP
				call skipcomment(unit)
				read(unit,*) smnp(i).slave,smnp(i).sdof,smnp(i).master,smnp(i).mdof
			end do
			
		case('geostatic')
			print *, 'Reading GEOSTATIC data...'
			do i=1,pro_num
				select case(property(i).name)
					case('method')
						geostatic.method=int(property(i).value)
					case('nsoil','soil')
						geostatic.nsoil=int(property(i).value)
					case default
						call Err_msg(property(i).name)
				end select
			end do
			geostatic.isgeo=.true.
			if(geostatic.method==ko_geo) then
				allocate(geostatic.ko(geostatic.nsoil),& 
					geostatic.weight(geostatic.nsoil), &
					geostatic.height(0:geostatic.nsoil))
				call skipcomment(unit)
				read(unit,*) geostatic.ko,geostatic.weight,geostatic.height
			end if
		case('feasible')
			print *, 'Reading LIMIT ANALYSIS FEASIBLE SOLUCTION data'
			do i=1,pro_num
				select case(property(i).name)
					case('num')
						n1=int(property(i).value)
					case('eset','element set')
						n2=int(property(i).value)
					case('col')
						n3=int(property(i).value)						
					case default
						call Err_msg(property(i).name)
				end select
			end do	
			if(n3-1/=element(eset(n2).enums).ndof) then
				cstring='the NDOF of the element set  is not consistent with what are read in'
				call Mistake_msg(cstring) 
			end if					
			do i=1,n1
				call skipcomment(unit)
				read(unit,*) n4, ar(1:n3-1)
				n5=eset(n2).enums+n4-1
				do j=1,n3-1
					load(element(n5).g(j))=ar(j)
				end do				
			end do
		case('lf1d')
			print *, 'Reading ONE DIMENSIONAL LINEAR FUNCTION data...'
			do i=1,pro_num
				select case(property(i).name)
					case('num')
						n1=int(property(i).value)
					case default
						call Err_msg(property(i).name)
				end select
				call skipcomment(unit)
				read(unit,*) (LF1D(J,:),J=1,n1)
			end do
		case('coordinate','system')
			print *, 'Reading coordinate system data...'	
				do i=1,pro_num
					select case(property(i).name)
						case('num')
							ncoord=int(property(i).value)
						case default
							call Err_msg(property(i).name)
					end select
				end do
                
				allocate(coordinate(-1:ncoord))
				
                !-1 IS FOR ZT4_SPG,ZT6_SPG, WHICH C IS DETERMINED FROM ELEMENT GEOMETRY.
				do i=1,ncoord
					call skipcomment(unit)
					read(unit,*) ((coordinate(i).c(j,k),k=1,3),j=1,3)
				end do
                coordinate(-1).C=0.0D0;coordinate(0).C=0.0D0
				coordinate(-1:0).c(1,1)=1.0d0
				coordinate(-1:0).c(2,2)=1.0d0
				coordinate(-1:0).c(3,3)=1.0d0
		case('origin_of_syscylinder')
			print *, 'Reading the origin of a cylinder system...'
			call skipcomment(unit)
			read(unit,*) Origin_Sys_cylinder
		case('output data','output variable','outdata','outvar')
			print *, 'Reading OUTPUT DATA control keyword'
			!**Attention**, to date, all the output variables keyword must be given in a line.
			!that is, it must be input like that: outvar,x,y,z,...(limited to 1024 character)
			do i=1,pro_num
				select case(property(i).name)
					case('x')
					  outvar(locx).name='X'
					  outvar(locx).value=locx
					  outvar(locx).system=property(i).value
					case('y')
					  outvar(locy).name='Y'  
					  outvar(locy).value=locy
					  outvar(locy).system=property(i).value
					case('z')
					  outvar(locz).name='Z'  
					  outvar(locz).value=locz
					  outvar(locz).system=property(i).value
					case('disx')
					  outvar(disx).name='disx'  
					  outvar(disx).value=disx
					  outvar(disx).system=property(i).value
					case('disy')
					  outvar(disy).name='disy'  
					  outvar(disy).value=disy
					  outvar(disy).system=property(i).value
					case('disz')
					  outvar(disz).name='disz'  
					  outvar(disz).value=disz
					  outvar(disz).system=property(i).value
					
					case('sxx')
					  outvar(sxx).name='sxx'  
					  outvar(sxx).value=sxx
					  outvar(sxx).system=property(i).value
					case('syy')
					  outvar(syy).name='syy'  
					  outvar(syy).value=syy
					  outvar(syy).system=property(i).value
					case('szz')
					  outvar(szz).name='szz'  
					  outvar(szz).value=szz
					  outvar(szz).system=property(i).value
					case('sxy')
					  outvar(sxy).name='sxy'  
					  outvar(sxy).value=sxy
					  outvar(sxy).system=property(i).value
					case('syz')
					  outvar(syz).name='syz'  
					  outvar(syz).value=syz
					  outvar(syz).system=property(i).value
					case('szx')
					  outvar(szx).name='szx'  
					  outvar(szx).value=szx
					  outvar(szx).system=property(i).value
					case('exx')
					  outvar(exx).name='exx'  
					  outvar(exx).value=exx
					  outvar(exx).system=property(i).value
					case('eyy')
					  outvar(eyy).name='eyy'  
					  outvar(eyy).value=eyy
					  outvar(eyy).system=property(i).value
					case('ezz')
					  outvar(ezz).name='ezz'  
					  outvar(ezz).value=ezz
					  outvar(ezz).system=property(i).value
					case('exy')
					  outvar(exy).name='exy'  
					  outvar(exy).value=exy
					  outvar(exy).system=property(i).value
					case('eyz')
					  outvar(eyz).name='eyz'  
					  outvar(eyz).value=eyz
					  outvar(eyz).system=property(i).value
					case('ezx')
					  outvar(ezx).name='ezx'  
					  outvar(ezx).value=ezx
					  outvar(ezx).system=property(i).value
					case('pexx')
					  outvar(pexx).name='pexx'  
					  outvar(pexx).value=pexx
					  outvar(pexx).system=property(i).value
					case('peyy')
					  outvar(peyy).name='peyy'  
					  outvar(peyy).value=peyy
					  outvar(peyy).system=property(i).value
					case('pezz')
					  outvar(pezz).name='pezz'  
					  outvar(pezz).value=pezz
					  outvar(pezz).system=property(i).value
					case('pexy')
					  outvar(pexy).name='pexy'  
					  outvar(pexy).value=pexy
					  outvar(pexy).system=property(i).value
					case('peyz')
					  outvar(peyz).name='peyz'  
					  outvar(peyz).value=peyz
					  outvar(peyz).system=property(i).value
					case('pezx')
					  outvar(pezx).name='pezx'  
					  outvar(pezx).value=pezx
					  outvar(pezx).system=property(i).value
					case('sxxg')
					  outvar(sxxg).name='sxxg'  
					  outvar(sxxg).value=sxxg
					  outvar(sxxg).system=property(i).value
					case('syyg')
					  outvar(syyg).name='syyg'  
					  outvar(syyg).value=syyg
					  outvar(syyg).system=property(i).value
					case('szzg')
					  outvar(szzg).name='szzg'  
					  outvar(szzg).value=szzg
					  outvar(szzg).system=property(i).value
					case('sxyg')
					  outvar(sxyg).name='sxyg'  
					  outvar(sxyg).value=sxyg
					  outvar(sxyg).system=property(i).value
					case('syzg')
					  outvar(syzg).name='syzg'  
					  outvar(syzg).value=syzg
					  outvar(syzg).system=property(i).value
					case('szxg')
					  outvar(szxg).name='szxg'  
					  outvar(szxg).value=szxg
					  outvar(szxg).system=property(i).value
					case('exxg')
					  outvar(exxg).name='exxg'  
					  outvar(exxg).value=exxg
					  outvar(exxg).system=property(i).value
					case('eyyg')
					  outvar(eyyg).name='eyyg'  
					  outvar(eyyg).value=eyyg
					  outvar(eyyg).system=property(i).value
					case('ezzg')
					  outvar(ezzg).name='ezzg'  
					  outvar(ezzg).value=ezzg
					  outvar(ezzg).system=property(i).value
					case('exyg')
					  outvar(exyg).name='exyg'  
					  outvar(exyg).value=exyg
					  outvar(exyg).system=property(i).value
					case('eyzg')
					  outvar(eyzg).name='eyzg'  
					  outvar(eyzg).value=eyzg
					  outvar(eyzg).system=property(i).value
					case('ezxg')
						outvar(ezxg).name='ezxg'  
						outvar(ezxg).value=ezxg
						outvar(ezxg).system=property(i).value
					case('pexxg')
						outvar(pexxg).name='pexxg'  
						outvar(pexxg).value=pexxg
						outvar(pexxg).system=property(i).value
					case('peyyg')
						outvar(peyyg).name='peyyg'  
						outvar(peyyg).value=peyyg
						outvar(peyyg).system=property(i).value
					case('pezzg')
						outvar(pezzg).name='pezzg'  
						outvar(pezzg).value=pezzg
						outvar(pezzg).system=property(i).value
					case('pexyg')
						outvar(pexyg).name='pexyg'  
						outvar(pexyg).value=pexyg
						outvar(pexyg).system=property(i).value
					case('peyzg')
						outvar(peyzg).name='peyzg'  
						outvar(peyzg).value=peyzg
						outvar(peyzg).system=property(i).value
					case('pezxg')
						outvar(pezxg).name='pezxg'  
						outvar(pezxg).value=pezxg
						outvar(pezxg).system=property(i).value
					case('pw')
						outvar(pw).name='pw'  
						outvar(pw).value=pw
						outvar(pw).iscentre=.true.
					case('mises')
						outvar(sigma_mises).name='mises'
						outvar(sigma_mises).value=sigma_mises
					case('eeq')
						outvar(eeq).name='eeq'
						outvar(eeq).value=eeq
					case('peeq')
						outvar(peeq).name='peeq'
						outvar(peeq).value=peeq
					case('xf')
						outvar(xf_out).name='xf'
						outvar(xf_out).value=xf_out
						outvar(xf_out).system=property(i).value
					case('yf')
						outvar(yf_out).name='yf'
						outvar(yf_out).value=yf_out
						outvar(yf_out).system=property(i).value
					case('zf')
						outvar(zf_out).name='zf'
						outvar(zf_out).value=zf_out
						outvar(zf_out).system=property(i).value
					case('rx')
						outvar(rx).name='rx'
						outvar(rx).value=rx
						outvar(rx).system=property(i).value
					case('ry')
						outvar(ry).name='ry'
						outvar(ry).value=ry
						outvar(ry).system=property(i).value
					case('rz')
						outvar(rz).name='rz'
						outvar(rz).value=rz
						outvar(rz).system=property(i).value
					case('qx')
						outvar(qx).name='qx'
						outvar(qx).value=qx
						outvar(qx).system=property(i).value
					case('qy')
						outvar(qy).name='qy'
						outvar(qy).value=qy
						outvar(qy).system=property(i).value
					case('qz')
						outvar(qz).name='qz'
						outvar(qz).value=qz
						outvar(qz).system=property(i).value	
					case('mx')
						outvar(mx).name='mx'
						outvar(mx).value=mx
						outvar(mx).system=property(i).value
					case('my')
						outvar(my).name='my'
						outvar(my).value=my
						outvar(my).system=property(i).value
					case('mz')
						outvar(mz).name='mz'
						outvar(mz).value=mz
						outvar(mz).system=property(i).value							
					case('gradx','ix')
						outvar(Gradx).name='Ix'
						outvar(Gradx).value=Gradx
					case('grady','iy')
						outvar(Grady).name='Iy'
						outvar(Grady).value=Grady
					case('gradz','iz')
						outvar(Gradz).name='Iz'
						outvar(Gradz).value=Gradz
					case('vx')
						outvar(Vx).name='Vx'
						outvar(Vx).value=Vx
					case('vy')
						outvar(vy).name='Vy'
						outvar(vy).value=vy
					case('vz')
						outvar(vz).name='Vz'
						outvar(vz).value=vz	
					case('head','h')
						outvar(head).name='H'
						outvar(head).value=head
					case('q','discharge')
						outvar(discharge).name='Q'
						outvar(discharge).value=discharge
					case('phead','ph','pressure head')
						outvar(phead).name='PH'
						outvar(phead).value=Phead
					case('kr')
						outvar(kr_spg).name='Kr'
						outvar(kr_spg).value=kr_spg
					case('mw')
						outvar(mw_spg).name='mw'
						outvar(mw_spg).value=mw_spg
					case('sfr')
						outvar(SFR).name='SFR'
						outvar(SFR).value=SFR
						outvar(sfr_sita).name='sfr_sita(deg(CCW+))'
						outvar(sfr_sita).value=sfr_sita
						outvar(SFR_SN).name='Sn(Tension+)'
						outvar(SFR_SN).value=SFR_SN
						outvar(SFR_TN).name='Tn(CCW+)'
						outvar(SFR_TN).value=SFR_TN
						outvar(SFR_SFRX).name='SFRX'
						outvar(SFR_SFRX).value=SFR_SFRX
						outvar(SFR_SFRY).name='SFRY'
						outvar(SFR_SFRY).value=SFR_SFRY
                        outvar(MC_C).name='MC_C'
                        outvar(MC_C).VALUE=MC_C
                        OUTVAR(MC_PHI).NAME='MC_PHI'
                        OUTVAR(MC_PHI).VALUE=MC_PHI
                        OUTVAR(SLOPE_SD).NAME='SlideDirection'
                        OUTVAR(SLOPE_SD).VALUE=SLOPE_SD
                        !OUTVAR(SFR_ALPHA).NAME='PSIGMASURFACE'
                        !OUTVAR(SFR_ALPHA).VALUE=SFR_ALPHA
                        !OUTVAR(SFR_PSITA).NAME='SFRMAX_ANGLE_WITH_PSS'
                        !OUTVAR(SFR_PSITA).VALUE=SFR_PSITA
                        !OUTVAR(SFR_MCSITA).NAME='MC_FAILURESURFACE_WITH_PSS'
                        !OUTVAR(SFR_MCSITA).VALUE=SFR_MCSITA                        
						!outvar(SFR).nval=6
					case('spg')
						outvar(Gradx).name='Ix'
						outvar(Gradx).value=Gradx
						outvar(Grady).name='Iy'
						outvar(Grady).value=Grady
						outvar(vx).name='Vx'
						outvar(vx).value=Vx
						outvar(Vy).name='Vy'
						outvar(Vy).value=Vy
						outvar(head).name='H'
						outvar(head).value=head	
						outvar(discharge).name='Q'
						outvar(discharge).value=discharge
						outvar(phead).name='PH'
						outvar(phead).value=Phead
						outvar(kr_spg).name='Kr'
						outvar(kr_spg).value=kr_spg
						outvar(mw_spg).name='mw'
						outvar(mw_spg).value=mw_spg
						
						IF(NDIMENSION>2) THEN
							outvar(Gradz).name='Iz'
							outvar(Gradz).value=Gradz
							outvar(Vz).name='Vz'
							outvar(Vz).value=Vz							
						ENDIF						
					case('dis')
					  outvar(disx).name='disx'  
					  outvar(disx).value=disx
					  outvar(disx).system=property(i).value						
					  outvar(disy).name='disy'  
					  outvar(disy).value=disy
					  outvar(disy).system=property(i).value
					  IF(NDIMENSION>2) THEN
						outvar(disz).name='disz'  
						outvar(disz).value=disz
						outvar(disz).system=property(i).value
					  ENDIF
					case('disf')
					  outvar(xf_out).name='xf'  
					  outvar(xf_out).value=xf_out
					  outvar(xf_out).system=property(i).value						
					  outvar(yf_out).name='yf'  
					  outvar(yf_out).value=xf_out
					  outvar(yf_out).system=property(i).value
					  IF(NDIMENSION>2) THEN
						outvar(zf_out).name='zf'  
						outvar(zf_out).value=zf_out
						outvar(zf_out).system=property(i).value
					  ENDIF					  
					case('nf','force')
						outvar(NFX).name='NFX'
						outvar(NFX).value=NFX
						outvar(NFY).name='NFY'
						outvar(NFY).value=NFY
						IF(NDIMENSION>2) THEN
							outvar(NFZ).name='NFZ'
							outvar(NFZ).value=NFZ
						ENDIF
						!IF(NDIMENSION==2) outvar(NF).name='NFX","NFY'
						
						!outvar(NF).nval=NDIMENSION	
					CASE('stress','s')
						outvar(sxx).name='sxx'  
						outvar(sxx).value=sxx
						outvar(syy).name='syy'  
						outvar(syy).value=syy
						outvar(szz).name='szz'  
						outvar(szz).value=szz
						outvar(sxy).name='sxy'  
						outvar(sxy).value=sxy
						outvar(SXX:SXY).system=property(i).value
						outvar(sigma_mises).name='mises'
						outvar(sigma_mises).value=sigma_mises
						OUTVAR(PSIGMA1).VALUE=PSIGMA1
						OUTVAR(PSIGMA1).NAME='PSIGMA1'
						OUTVAR(PSIGMA3).VALUE=PSIGMA3
						OUTVAR(PSIGMA3).NAME='PSIGMA3'						
						OUTVAR(PSIGMA2).VALUE=PSIGMA2
						OUTVAR(PSIGMA2).NAME='PSIGMA2'
						OUTVAR(APSIGMA1).VALUE=APSIGMA1
						OUTVAR(APSIGMA1).NAME='APSIGMA1'						
						IF(NDIMENSION>2) THEN
							outvar(syz).name='syz'  
							outvar(syz).value=syz
							outvar(szx).name='szx'  
							outvar(szx).value=szx
							outvar(SYZ:SZX).system=property(i).value
						ENDIF
					CASE('strain','e')
						outvar(exx).name='exx'  
						outvar(exx).value=exx						
						outvar(eyy).name='eyy'  
						outvar(eyy).value=eyy	
						outvar(ezz).name='ezz'  
						outvar(ezz).value=ezz
						outvar(exy).name='exy'  
						outvar(exy).value=exy
						outvar(eeq).name='eeq'
						outvar(eeq).value=eeq
						outvar(EXX:EXY).system=property(i).value
						IF(NDIMENSION>2) THEN
						  outvar(eyz).name='eyz'  
						  outvar(eyz).value=eyz	
						  outvar(ezx).name='ezx'  
						  outvar(ezx).value=ezx
						  outvar(EYZ:EZX).system=property(i).value
						ENDIF
					CASE('pstrain','pe')
						outvar(pexx).name='pexx'  
						outvar(pexx).value=pexx						
						outvar(peyy).name='peyy'  
						outvar(peyy).value=peyy	
						outvar(pezz).name='pezz'  
						outvar(pezz).value=pezz
						outvar(pexy).name='pexy'  
						outvar(pexy).value=pexy
						outvar(peeq).name='peeq'
						outvar(peeq).value=peeq
						outvar(PEXX:PEXY).system=property(i).value
						IF(NDIMENSION>2) THEN
						  outvar(peyz).name='peyz'  
						  outvar(peyz).value=peyz	
						  outvar(pezx).name='pezx'  
						  outvar(pezx).value=pezx
						  outvar(PEYZ:PEZX).system=property(i).value
						ENDIF
					CASE('m','moment')
						outvar(mz).name='mz'
						outvar(mz).value=mz
						outvar(mz).system=property(i).value
						if(ndimension>2) then
							outvar(mx).name='mx'
							outvar(mx).value=mx
							outvar(mx).system=property(i).value	
							outvar(my).name='my'
							outvar(my).value=my
							outvar(my).system=property(i).value
						endif
					CASE('r','rotate')
						if(ndimension>2) then
							outvar(rx).name='rx'
							outvar(rx).value=rx
							outvar(rx).system=property(i).value
							outvar(ry).name='ry'
							outvar(ry).value=ry
							outvar(ry).system=property(i).value
						endif
						outvar(rz).name='rz'
						outvar(rz).value=rz
						outvar(rz).system=property(i).value                    
                        
					case default
						call Err_msg(property(i).name)
				end select
			end do			
								
		case default
			call Err_msg(term)
	end select


end subroutine



   !把字符串中相当的数字字符(包括浮点型)转化为对应的数字
   !如 '123'转为123,'14-10'转为14,13,12,11,10
   !string中转化后的数字以数组ar(n1)返回，其中,n1为字符串中数字的个数:(注　1-3转化后为3个数字：1,2,3)
   !nmax为数组ar的大小,string默认字符长度为1024。
   !num_read为要读入数据的个数。
   !unit为文件号
   !每次只读入一个有效行（不以'/'开头的行）
   !每行后面以'/'开始的后面的字符是无效的。
   subroutine  strtoint(unit,ar,nmax,n1,num_read,set,maxset,nset)
	  implicit none
      INTEGER,INTENT(IN)::unit,nmax,num_read,maxset
      INTEGER,INTENT(INOUT)::N1,NSET
      REAL(8),INTENT(INOUT)::ar(nmax)
      character(*)::set(maxset)
	  logical::tof1,tof2
	  integer::i,j,k,strl,ns,ne,n2,n3,n4,step,& 
			ef,n5,nsubs
	  real(8)::t1	  
	  character(1024)::string
	  character(32)::substring(100)
	  character(16)::legalC,SC

		LegalC='0123456789.-+eE*'
		sc=',; '//char(9)
		n1=0
		nset=0
		ar=0
		!set(1:maxset)=''
	  do while(.true.)
		 read(unit,'(a1024)',iostat=ef) string
		 if(ef<0) then
			print *, 'file ended unexpected. sub strtoint()'
			stop
		 end if

		 string=adjustL(string)
		 strL=len_trim(string)
		 
		do i=1,strL !remove 'Tab'
			if(string(i:i)/=char(9)) exit
		end do
		string=string(i:strL)
		string=adjustl(string)
		strL=len_trim(string)
		if(strL==0) cycle

		 if(string(1:2)/='//'.and.string(1:1)/='#') then
			
			!每行后面以'/'开始的后面的字符是无效的。
			if(index(string,'//')/=0) then
				strL=index(string,'//')-1
				string=string(1:strL)
				strL=len_trim(string)
			end if

			nsubs=0
			n5=1
			do i=2,strL+1
				if(index(sc,string(i:i))/=0.and.index(sc,string(i-1:i-1))==0) then
					nsubs=nsubs+1					
					substring(nsubs)=string(n5:i-1)					
				end if
				if(index(sc,string(i:i))/=0) n5=i+1
			end do
			
			do i=1, nsubs
				substring(i)=adjustl(substring(i))				
				n2=len_trim(substring(i))
				!the first character should not be a number if the substring is a set.
				if(index('0123456789-+.', substring(i)(1:1))==0) then
					!set
					nset=nset+1
					set(nset)=substring(i)
					cycle
				end if
				n3=index(substring(i),'-')
				n4=index(substring(i),'*')
				tof1=.false.
				if(n3>1) then
				    tof1=(substring(i)(n3-1:n3-1)/='e'.and.substring(i)(n3-1:n3-1)/='E')
				end if
				if(tof1) then !处理类似于'1-5'这样的形式的读入数据
					read(substring(i)(1:n3-1),'(i8)') ns
					read(substring(i)(n3+1:n2),'(i8)') ne
					if(ns>ne) then
						step=-1
					else
						step=1
					end if
					do k=ns,ne,step
						n1=n1+1
						ar(n1)=k
					end do				     	
				else
				     tof2=.false.
				     if(n4>1) then
				             tof2=(substring(i)(n4-1:n4-1)/='e'.and.substring(i)(n4-1:n4-1)/='E')
				     end if
					if(tof2) then !处理类似于'1*5'(表示5个1)这样的形式的读入数据
						read(substring(i)(1:n4-1),*) t1
						read(substring(i)(n4+1:n2),'(i8)') ne
						ar((n1+1):(n1+ne))=t1
						n1=n1+ne
					else
						n1=n1+1
						read(substring(i),*) ar(n1)
					end if	
				end if			
			end do
		 else
			cycle
		 end if
		
		 if(n1<=num_read) then
		    exit
		 else
		    if(n1>num_read)  print *, 'error!nt2>num_read. i=',n1
		 end if
	
	  end do	

   end subroutine

subroutine translatetoproperty(term)

!**************************************************************************************************************
!Get a keyword and related property values from a control line (<1024)
!input variables:
!term, store control data line content.
!ouput variables:
!property,pro_num
!mudulus used:
!None
!Subroutine called:
!None
!Programmer:LUO Guanyong
!Last updated: 2008,03,20

!Example: 
!term='element, num=10,type=2,material=1,set=3'
!after processed,the following data will be returned:
!term='element'
!property(1).name=num
!property(1).value=1
!.....
!**************************************************************************************************************
	use solverds
	implicit none
	integer::i,strL
	character(1024)::term,keyword
	integer::ns,ne,nc
	character(128)::str(50)
	
	if(index(term,'//')/=0) then !每一行‘/’后面的内容是无效的。
		strL=index(term,'//')-1
		term=term(1:strL)
	end if
	
	term=adjustl(term)
	ns=1
	ne=0
	nc=0
	property.name=''
	property.value=0.0
	property.cvalue=''
	do while(len_trim(term)>0) 
		nc=nc+1
		if(nc>51) then
			print *, 'nc>51,subroutine translatetoproperty()'
			stop
		end if
		ne=index(term,',')
		if(ne>0.and.len_trim(term)>1) then
			str(nc)=term(ns:ne-1)
			str(nc)=adjustL(str(nc))
		else 
		!no commas in the end
			ne=min(len_trim(term),len(str(nc)))
			str(nc)=term(ns:ne)
			str(nc)=adjustL(str(nc))
		end if
		term=term(ne+1:len_trim(term))
		term=adjustL(term)		
	end do

	
    term=str(1)    
   
    
    !TRANSLATE "A=1" TO"A,A=1" 
    ne=index(str(1),'=')
    IF(NE>1) THEN
        TERM=STR(1)(1:NE-1)        
        STR(NC+1:3:-1)=STR(NC:2:-1)
        STR(2)=STR(1)(NE+1:LEN_TRIM(ADJUSTL(STR(1))))
        NC=NC+1
    ENDIF
    
	pro_num=nc-1
	do i=2,nc
		ne=index(str(i),'=')
		if(ne>0) then
			property(i-1).name=str(i)(1:ne-1)
            !ns=len_trim(str(i))
			ns=len_trim(adjustl(str(i)))-ne
            
			call inp_ch_c_to_int_c(str(i)(ne+1:len_trim(adjustl(str(i)))),ns,property(i-1).value,property(i-1).cvalue)
		else
			property(i-1).name=str(i)(1:len_trim(str(i)))
		end if
		!read(str(i)(ne+1:len_trim(str(i))),*) property(i-1).value
	end do

end subroutine


subroutine ettonnum(et1,nnum1,ndof1,ngp1,nd1,stype,EC1)
!according to the element type return how many nodes for each this type element, and how many dofs
!for each this type element
!initializing the element class property
	use solverds
	implicit none
	integer::et1,nnum1,ndof1,ngp1,nd1,EC1
	character(16)::stype

	ngp1=0
	nd1=0

	select case(et1)
		case(CONDUCT1D)
			nnum1=2
			ndof1=2
			stype='FELINESEG'
			ec1=CND
		case(UB3)
			nnum1=3
			ndof1=6
			stype='FETRIANGLE'
			ec1=LMT
		case(UBZT4)
			nnum1=4
			ndof1=8
			stype='FEQUADRILATERAL'
			ec1=LMT
		case(LB3)
			nnum1=3
			ndof1=9
			stype='FETRIANGLE'
			ec1=LMT
		case(LBZT4)
			nnum1=4
			ndof1=12
			stype='FEQUADRILATERAL'
			ec1=LMT
		case(CPE3,CPS3,CAX3)
			nnum1=3
			ndof1=6
			ngp1=1
			nd1=4
			stype='FETRIANGLE'
			if(et1==cpe3)then
				EC1=CPE
			end if
			if(et1==cps3)EC1=cps
			if(et1==CAX3)then
				EC1=CAX
			end if
			call EL_SFR2(ET1)
		case(CPE6,CPS6,CAX6)
			nnum1=6
			ndof1=12
			ngp1=3
			nd1=4
			stype='FETRIANGLE'
			if(et1==cpe6) then
				ec1=CPE
			end if
			if(et1==cpS6)	ec1=cps
			if(et1==CAX6)then
				EC1=CAX
			end if
			call EL_SFR2(ET1)			
		case(CPE4,CPS4,cax4)
			nnum1=4
			ndof1=8
			ngp1=4			
			nd1=4
			stype='FEQUADRILATERAL'
			if(et1==cpe4) then
				ec1=CPE
			end if
			if(et1==cpS4)	ec1=cps
			if(et1==CAX4)then
				EC1=CAX
			end if
			call EL_SFR2(ET1)			
		case(CPE4R,CPS4R,CAX4R)
			nnum1=4
			ndof1=8
			ngp1=1			
			nd1=4
			stype='FEQUADRILATERAL'	
			if(et1==cpe4r) then
				ec1=CPE
			end if
			if(et1==cpS4R)	ec1=cps
			if(et1==CAX4r)then
				EC1=CAX
			end if
			call EL_SFR2(ET1)			
		case(CPE8,CPS8,CAX8)
			nnum1=8
			ndof1=16
			ngp1=9
			nd1=4
			stype='FEQUADRILATERAL'
			if(et1==cpe8) then
				ec1=CPE
			end if
			if(et1==CPS8)	ec1=cps
			if(et1==CAX8)then
				EC1=CAX
			end if
			call EL_SFR2(ET1)			
		case(CPE8R,CPS8R,CAX8R)
			nnum1=8
			ndof1=16
			ngp1=4
			nd1=4
			stype='FEQUADRILATERAL'
			if(et1==cpe8r) then
				ec1=CPE
			end if
			if(et1==CPS8R)	ec1=cps
			if(et1==CAX8r)then
				EC1=CAX
			end if
			call EL_SFR2(ET1)			
		case(CPE15,CPS15,CAX15)
			nnum1=15
			ndof1=30
			ngp1=12
			nd1=4
			stype='FETRIANGLE'
			if(et1==cpe15) then
				ec1=CPE
			end if
			if(et1==CPS15)	ec1=cps
			if(et1==CAX15)then
				EC1=CAX
			end if
			call EL_SFR2(ET1)
		!case(PRM6)
		!	nnum1=6
		!	ndof1=18
		!	ngp1=2
		!	nd1=6
		!	stype='FEBRICK'
		!	EC1=C3D
		!	call EL_SFR2(ET1)
		!case(PRM15)
		!	nnum1=15
		!	ndof1=45
		!	ngp1=9
		!	nd1=6
		!	stype='FETETRAHEDRON' !
		!	EC1=C3D
		!	call EL_SFR2(ET1)
		!case(TET10)
		!	nnum1=10
		!	ndof1=30
		!	ngp1=4
		!	nd1=6
		!	stype='FETETRAHEDRON' !
		!	EC1=C3D
		!	call EL_SFR2(ET1)			
		case(CPE3_SPG,CAX3_SPG,CPE3_CPL,CAX3_CPL)
			nnum1=3
			!ndof1=3
			ngp1=1
			
			stype='FETRIANGLE'
			IF(ET1==CPE3_SPG) then
				EC1=SPG2D;NDOF1=3;nd1=2
			endif
			if(et1==CAX3_SPG) then
				EC1=CAX_SPG;NDOF1=3;nd1=2
			endif
			IF(ET1==CPE3_CPL) THEN
				EC1=CPL;NDOF1=9;nd1=4
			ENDIF
			if(et1==CAX3_CPL) THEN
				EC1=CAX_CPL;NDOF1=9;nd1=4
            ENDIF
			call EL_SFR2(ET1)

		case(CPE6_SPG,CAX6_SPG,CPE6_CPL,CAX6_CPL)
			nnum1=6
			!ndof1=6
			ngp1=3
			
			stype='FETRIANGLE'
            
			IF(ET1==CPE6_SPG) THEN
				EC1=SPG2D;NDOF1=6;nd1=2
			ENDIF
			if(et1==CAX6_SPG) THEN
				EC1=CAX_SPG;NDOF1=6;nd1=2
			ENDIF
			IF(ET1==CPE6_CPL) THEN
				EC1=CPL;NDOF1=18;nd1=4
			ENDIF
			if(et1==CAX6_CPL) THEN
				EC1=CAX_CPL;NDOF1=18;nd1=4
            ENDIF
			call EL_SFR2(ET1)			
		case(CPE4_SPG,cax4_SPG,CPE4_CPL,cax4_CPL,ZT4_SPG)
			nnum1=4
			!ndof1=4
			ngp1=4			
			
			stype='FEQUADRILATERAL'
			IF(ET1==CPE4_SPG.OR.ET1==ZT4_SPG) THEN
				EC1=SPG2D;NDOF1=4;nd1=2
			ENDIF
			if(et1==CAX4_SPG) THEN
			EC1=CAX_SPG;NDOF1=4;nd1=2
			ENDIF
			IF(ET1==CPE4_CPL) THEN
			EC1=CPL;NDOF1=12;nd1=4
			ENDIF
			if(et1==CAX4_CPL) THEN
			EC1=CAX_CPL;NDOF1=12;nd1=4
			ENDIF
			call EL_SFR2(ET1)
            
		case(CPE4R_SPG,CAX4R_SPG)
			nnum1=4
			!ndof1=4
			ngp1=1			
			
			stype='FEQUADRILATERAL'	
			IF(ET1==CPE4R_SPG) THEN
			EC1=SPG2D;NDOF1=4;nd1=2
			ENDIF
			if(et1==CAX4R_SPG) THEN
			EC1=CAX_SPG;NDOF1=4;nd1=2
			ENDIF
			IF(ET1==CPE4R_CPL) THEN
			EC1=CPL;NDOF1=12;nd1=4
			ENDIF
			if(et1==CAX4R_CPL) THEN
			EC1=CAX_CPL;NDOF1=12;nd1=4
			ENDIF
            
			call EL_SFR2(ET1)			
		case(CPE8_SPG,CAX8_SPG,CPE8_CPL,CAX8_CPL)
			nnum1=8
			!ndof1=8
			ngp1=9
			
			stype='FEQUADRILATERAL'
			IF(ET1==CPE8_SPG) THEN
			EC1=SPG2D;NDOF1=8;nd1=2
			ENDIF
			if(et1==CAX8_SPG) THEN
			EC1=CAX_SPG;NDOF1=8;nd1=2
			ENDIF
			IF(ET1==CPE8_CPL) THEN
			EC1=CPL;NDOF1=24;nd1=4
			ENDIF
			if(et1==CAX8_CPL) THEN
			EC1=CAX_CPL;NDOF1=24;nd1=4  
			ENDIF
            
			call EL_SFR2(ET1)			
		case(CPE8R_SPG,CAX8R_SPG,CPE8R_CPL,CAX8R_CPL)
			nnum1=8
			!ndof1=8
			ngp1=4
			
			stype='FEQUADRILATERAL'
			IF(ET1==CPE8R_SPG) THEN
			EC1=SPG2D;NDOF1=8;nd1=2
			ENDIF
			if(et1==CAX8R_SPG) THEN
			EC1=CAX_SPG;NDOF1=8;nd1=2
			ENDIF
			IF(ET1==CPE8R_CPL) THEN
			EC1=CPL;NDOF1=24;nd1=4
			ENDIF
			if(et1==CAX8R_CPL) THEN
			EC1=CAX_CPL;NDOF1=24;nd1=4 
  			ENDIF
			call EL_SFR2(ET1)	
            
		case(CPE15_SPG,CAX15_SPG,CPE15_CPL,CAX15_CPL)
			nnum1=15
			!ndof1=15
			ngp1=12
			
			stype='FETRIANGLE'
			IF(ET1==CPE15_SPG) THEN
			EC1=SPG2D;NDOF1=15;nd1=2
			ENDIF
			if(et1==CAX15_SPG) THEN
			EC1=CAX_SPG;NDOF1=15;nd1=2
			ENDIF
			IF(ET1==CPE15_CPL) THEN
			EC1=CPL;NDOF1=45;nd1=4
			ENDIF
			if(et1==CAX15_CPL) THEN
			EC1=CAX_CPL;NDOF1=45;nd1=4           
			ENDIF			
			call EL_SFR2(ET1)
		case(PRM6_SPG,PRM6,PRM6_CPL,ZT6_SPG)
			nnum1=6
			!ndof1=6
			ngp1=2
			!nd1=3
			stype='FETETRAHEDRON'
			IF(ET1==PRM6_SPG.OR.ET1==ZT6_SPG) THEN
			EC1=SPG;NDOF1=6;nd1=3
			ENDIF
            IF(ET1==PRM6) THEN
			EC1=C3D;NDOF1=18;nd1=6
			ENDIF
            IF(ET1==PRM6_CPL) THEN
			EC1=CPL;NDOF1=24;nd1=6
            ENDIF
			call EL_SFR2(ET1)
		case(PRM15_SPG,PRM15,PRM15_CPL)
			nnum1=15
			!ndof1=15
			ngp1=9
			!nd1=3
			stype='FETETRAHEDRON' !
			IF(ET1==PRM15_SPG) THEN
			EC1=SPG;NDOF1=15;nd1=3
			ENDIF
            IF(ET1==PRM15) THEN
			EC1=C3D;NDOF1=45;nd1=6
			ENDIF
            IF(ET1==PRM15_CPL) THEN
			EC1=CPL;NDOF1=60;nd1=6
            ENDIF
			call EL_SFR2(ET1)
		case(tet4_spg,TET4,TET4_CPL)
			nnum1=4
			!ndof1=4
			ngp1=1
			
			stype='FETETRAHEDRON'
			IF(ET1==TET4_SPG) THEN
			ec1=SPG;NDOF1=4;nd1=3
			ENDIF
			IF(ET1==TET4_CPL) THEN
			ec1=CPL;NDOF1=16;nd1=6
			ENDIF
			IF(ET1==TET4) THEN
			ec1=C3D;NDOF1=12;nd1=6
			ENDIF
			CALL EL_SFR2(ET1)
		case(tet10_spg,TET10,TET10_CPL)
			nnum1=10
			!ndof1=10
			ngp1=4
			
			stype='FETETRAHEDRON'
			
			IF(ET1==TET10_SPG) THEN
			ec1=SPG;NDOF1=10;nd1=3
			ENDIF
			IF(ET1==TET10_CPL) THEN
			ec1=CPL;NDOF1=40;;nd1=6
			ENDIF
			IF(ET1==TET10) THEN
			ec1=C3D;NDOF1=30;nd1=6
			ENDIF
			CALL EL_SFR2(ET1)	
		case(BAR) !3d BAR ELEMENT
			nnum1=2
			ndof1=6
			nd1=6
			stype='FEBRICK' !
			EC1=STRU
			!call EL_SFR2(ET1)
		case(BAR2D) !3d BAR ELEMENT
			nnum1=2
			ndof1=4
			nd1=4
			stype='FEBRICK' !
			EC1=STRU			
		case(BEAM) !3D Beam element
			nnum1=2
			ndof1=12
			nd1=12
			stype='FEBRICK'
			EC1=STRU
		case(BEAM2D,SSP2D) !2D Beam element
			nnum1=2
			ndof1=6
			nd1=6
			stype='FEBRICK'
			EC1=STRU
		case(PE_SSP2D)
			nnum1=2
			ndof1=2
			nd1=2
			stype='FELINESEG'
			EC1=PE
		case(SSP2D1) !2D Beam element
			nnum1=2
			ndof1=10
			nd1=10
			stype='FEBRICK'
			EC1=STRU

		case(SHELL3,SHELL3_KJB)
			NNUM1=3
			NDOF1=18
			ND1=3
			STYPE='FETRIANGLE'
			EC1=STRU
		case(DKT3)
			NNUM1=3
			NDOF1=9
			ND1=3
			STYPE='FETRIANGLE'
			EC1=STRU
		case(pipe2,ppipe2)
			NNUM1=2
			NDOF1=2
			ND1=2
			STYPE='FELINESEG'
			EC1=PIPE
		case(springx,springy,springz,springmx,springmy,springmz)
			nnum1=1
			ndof1=1
			nd1=1
			stype='FELINESEG'
			ec1=spring
		case(SOILSPRINGX,SOILSPRINGY,SOILSPRINGZ)
			nnum1=1
			ndof1=1
			nd1=1
			stype='FELINESEG'
			ec1=soilspring
        CASE DEFAULT
            PRINT *, "NO SUCH ELEMENT TYPE. IN SUB ettonnum. TO BE IMPROVED."
            STOP
	end select

end subroutine

subroutine Err_msg(cstring)
	use dflib
	implicit none
	character(*)::cstring
	character(64)::term
	integer(4)::msg

	term='No such Constant: '//trim(cstring)
	msg = MESSAGEBOXQQ(term,'Caution'C,MB$ICONASTERISK.OR.MB$OKCANCEL.OR.MB$DEFBUTTON1)
	if(msg==MB$IDCANCEL) then
		stop
	end if	
	
end subroutine

subroutine Mistake_msg(cstring)
	use dflib
	implicit none
	character(*)::cstring
	character(128)::term
	integer(4)::msg

	term="Mistake:  "//trim(cstring)
	msg = MESSAGEBOXQQ(term,'Mistake'C,MB$ICONSTOP.OR.MB$OK.OR.MB$DEFBUTTON1)
	if(msg==MB$IDOK) then
		stop
	end if	
	
end subroutine

subroutine skipcomment(unit,EF)
	implicit none
	integer,intent(in)::unit
    INTEGER,OPTIONAL::EF
	integer::i,strL,EF1
	character(1024) string
	
	do while(.true.)
		read(unit,'(a1024)',iostat=ef1) string
		if(ef1<0) then
            IF(PRESENT(EF)) THEN
                EF=EF1
                RETURN
            END IF
			print *, 'file ended unexpected. sub skipcomment()'
			stop
		end if

		string=adjustL(string)
		strL=len_trim(string)

		do i=1,strL !remove 'Tab'
			if(string(i:i)/=char(9)) exit
		end do
		string=string(i:strL)
		string=adjustl(string)
		strL=len_trim(string)
		if(strL==0) cycle

		if(string(1:2)/='//'.and.string(1:1)/='#') then
			backspace(unit)
			exit
		end if
	end do
	
end subroutine

!=============================================================  
subroutine StringSplit(InStr,delimiter,StrArray,nsize)  
!----------------------------------------------  
!---将字符串InStr进行分割,结果放入StrArray中  
!---delimiter::分隔符号,例如';,,' 使用;和,分割字符串  
!---nsize:分割数目  
!---吴徐平2011-04-29(wxp07@qq.com)  
!----------------------------------------------  
implicit none  
character(len = *) , Intent( IN ) :: InStr  
character(len = *)  , Intent( IN ) :: delimiter  
character(len = LEN(InStr)),dimension(LEN(InStr)),Intent( OUT ) :: StrArray  
integer, Intent( OUT ) :: nsize ! Effective Size of StrArray  
integer:: i,j ! loop variable  
integer:: istart ! split index for Start Position  
nsize=0  
istart=1  
do i=1,LEN(InStr)  
    do j=1,LEN(delimiter)  
        if (InStr(i:i) == delimiter(j:j)) then  
            if (istart == i) then  
            istart=i+1 ! ---可防止分隔符相连的情况  
            end if  
            if (istart<i) then  
                nsize=nsize+1  
                StrArray(nsize)=InStr(istart:i-1)  
                istart=i+1  
            end if  
        end if  
    end do  
end do  
! ---匹配最后一个子字符串  
if (nsize>0) then  
    if (istart<LEN(InStr)) then  
        nsize=nsize+1  
        StrArray(nsize)=InStr(istart:LEN(InStr))  
    end if  
end if  
! ---如果无可分割的子字符串,则包含整个字符串为数组的第一元素  
if ( (nsize<1) .AND. (LEN(TRIM(InStr)) > 0 )) then  
        nsize=1  
        StrArray(1)=InStr  
end if  
end subroutine StringSplit



