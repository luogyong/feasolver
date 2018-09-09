
program main
	use SOLVERLIB
	use dflib
	implicit none
	integer::itype=0,ef
	LOGICAL(4) status
	logical::tof1
	character*8 char_time
	character(1)::key
	type(qwinfo) winfo
	REAL time_begin, time_end,time_begin1, time_end1

	CALL CPU_TIME (time_begin )
	!winfo%TYPE = QWIN$MAX 
	!status = SETWSIZEQQ(QWIN$FRAMEWINDOW, winfo) 
	!status = SETWSIZEQQ(0, winfo)
	Print *, 'FEASOLVER. LGY WORK.'
    !ef = setexitqq(QWIN$EXITNOPERSIST)
	call TIME(char_time) 
    !IF(SOLVER_CONTROL.NOPOPUP==0) THEN
	   ! write(*, 10) 
	   ! key=getcharqq()
	   ! if(ichar(key)==ichar('h').or.ichar(key)==ichar('H')) then		
		  !  call write_readme_FEASOLVER()			
		  !  stop
    !    end if
    !ENDIF
	
	call TIME(char_time) 
	PRINT *, 'Reading data...',char_time
	call readin(itype)
	call TIME(char_time)

	PRINT *, 'Initializing process...',char_time
	call Initialization()
	call TIME(char_time)
	PRINT *, 'Begin to solve...',char_time
	call solve()
	call TIME(char_time)
	PRINT *, 'Solve Completed!',char_time
	
    if(isexca2d/=0) then
        call SCIPLOT()
        key=getcharqq()
		ef = setexitqq(QWIN$EXITPERSIST)
    endif
    
   ef = setexitqq(QWIN$EXITNOPERSIST)
10 format("Press 'H/h' to write a keyword help file in the current direction . Any other key to read an SINP file.")		

end program main
