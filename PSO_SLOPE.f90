MODULE SLOPE_PSO
    USE stochoptim,ONLY: Evolutionary,EA_ITER
    USE solverds,ONLY:INPUTFILE
    USE INPUT_PARSER
    IMPLICIT NONE
    
    PRIVATE
    
    PUBLIC::SLOPE_OPTIM,PSO_SLIP,POLYLINE_FOS_CAL,PSO_SLIP_PLOT_MENU
    
    character(512)::PARAFILE,PARAFILENAME,dir,ext
	CHARACTER(3)::drive
    INTEGER::IFLAG_OPTIM=0,NITER=0
    REAL(8),ALLOCATABLE,DIMENSION(:)::IFLAG_OPTIM_PARA
    
    INTEGER,PARAMETER::SEARCH_BY_EA=0,EA_SHOW_ALL=1,EA_SHOW_MIN=2,EA_SHOW_TOPTEN=3,EA_SHOW_XSTART=4,&
                        EA_SHOW_NONE=5,EA_SHOW_XSTART_NONE=6
    TYPE PSO_SLIP_INFO
        INTEGER::NSLIP=0,NX=0,NFOSCAL=0   !NSLIP=POPSIZE     
        INTEGER,ALLOCATABLE::A2Z(:),COLOR(:)
        REAL(8),ALLOCATABLE::FOS(:),FA(:),FM(:)       
        REAL(8),ALLOCATABLE::X(:,:,:),XSTART(:,:,:)
        LOGICAL,ALLOCATABLE::ISSHOW(:),ISSHOW_XSTART(:)
        CONTAINS
        PROCEDURE::INIT=>PSO_SLIP_INIT
        PROCEDURE::SORT=>PSO_SLIP_SORT
        PROCEDURE::PLOT=>PSO_SLIP_PLOT
        PROCEDURE,NOPASS::PLOT_SETTING=>PSO_SLIP_PLOT_SETTING
        PROCEDURE,NOPASS::PLOT_MENU=>PSO_SLIP_PLOT_MENU
        PROCEDURE::EXPORT=>PSO_SLIP_EXPORT
    ENDTYPE
    TYPE(PSO_SLIP_INFO)::PSO_SLIP
    
!    REAL(8),ALLOCATABLE::XSI(:),YSI(:) !XSI,YSI IS CURRENT SLIP SURFACE 
    
    TYPE SLOPE_PSO_PARAR
        !SLOPE PARAMETER
        INTEGER::NGX,NRX,NSLICE=30,SHAPE=0 !SHAPE=0,NON-CIRCULAR,=1,CIRCULAR
        REAL(8)::XTU,XTL,XCU,XCL !ENTRY AND EXIT LIMITS
        REAL(8),ALLOCATABLE::GX(:,:),RX(:,:)
        !EA PARAMETER
        integer(kind = 4) :: popsize = 30, max_iter = 1000,ndim=31
        REAL(8)::EPS1=1.D-4,EPS2=1.D-4,EPS3=1.D-3,w = 0.7298, c1 = 1.49618,c2 = 1.49618, gamma = 1., &
            F=0.5,CR=0.1,DL=0.5D0,mu_perc=0.5,sigma=0.5  !DL=SLICE BASE STRESS SAMPLE DISTANCE
        character(len = 5) :: solver='de', strategy='rand1'
        
        !trial slip surface generation method
        INTEGER::SSGM=0 !slip =0，chen,=1,circular,=2circular+softband
        INTEGER::NWSP=0
        REAL(8),ALLOCATABLE::WSP(:,:)
        
        !!FOR STREAMLINE ADMISSIBLITY REPAIR ONLY
        INTEGER::DIRECTION=-1 !=-1 TOE,=1,(DOWNHILL),CREST(UPHILL).
        REAL(8)::A=0,B=0,C=0 !REPAIRED RANGE LINE ,A*X+B*Y+C
        CONTAINS
        PROCEDURE::INIT=>FILEREADIN
        PROCEDURE::GETY=>GET_BC_Y
        
    ENDTYPE
    TYPE(SLOPE_PSO_PARAR)::SLOPE_PARA
    
    
    !interface SLOPEPARAMETER
    !    module procedure GET_SLOPE_PARA,SET_SLOPE_PARA
    !end interface SLOPEPARAMETER
    
    
    


CONTAINS


SUBROUTINE SLOPE_OPTIM(PARAFILE,IFLAG,OPTS)
    
    INTEGER,OPTIONAL,INTENT(IN)::IFLAG
    REAL(8),DIMENSION(:),OPTIONAL,INTENT(IN)::OPTS 
    CHARACTER(1024),OPTIONAL,INTENT(IN)::PARAFILE
    
    CHARACTER(1024)::PFILE1="",IFILE1=""
    
    real(kind = 8), parameter :: PI = 3.141592653589793238460d0    
    real(kind = 8), dimension(:), allocatable :: lower, upper
    REAL(8),DIMENSION(:,:),ALLOCATABLE::X,Y,RDF,XSLIP1
    REAL(8)::FA1,FM1
    LOGICAL::FILE_EXIST1
    INTEGER::ISERROR1,I,N1
    
    type(Evolutionary) :: ea
        
    IF(PRESENT(PARAFILE)) THEN
        PFILE1=PARAFILE
    ELSE
        IFILE1=TRIM(ADJUSTL(INPUTFILE))
        call lowcase(IFILE1)
        PFILE1=TRIM(ADJUSTL(IFILE1))//'_eap.dat'        
        FILE_EXIST1=.FALSE.
        inquire(file=TRIM(ADJUSTL(PFILE1)), exist=FILE_EXIST1) 
        IF(.NOT.FILE_EXIST1) THEN
            
            N1=INDEX(TRIM(ADJUSTL(IFILE1)),'_tec')
            IF(N1>0) THEN
                PFILE1=INPUTFILE(1:N1-1)//'_eap.dat'
            ELSE
                PFILE1=""            
            ENDIF            
        ENDIF         
    ENDIF

    IF(PRESENT(IFLAG)) THEN
        IFLAG_OPTIM=IFLAG
    ELSE
        IFLAG_OPTIM=0
    ENDIF

    IF(PRESENT(OPTS)) THEN
        IF(ALLOCATED(IFLAG_OPTIM_PARA)) DEALLOCATE(IFLAG_OPTIM_PARA)
        ALLOCATE(IFLAG_OPTIM_PARA,SOURCE=OPTS)
    ENDIF
    
    CALL SLOPE_PARA.INIT(PFILE1) 
    !IF(SLOPE_PARA.IS_SLOPEINIT==0) THEN        
    !    CALL SLOPE_PARA.INIT(PFILE1) 
    !    SLOPE_PARA.IS_SLOPEINIT=1        
    !ENDIF
    
    IF(IFLAG_OPTIM==1) THEN !OPTIMIZE THE STREAMLINE 
        SLOPE_PARA.XCL=IFLAG_OPTIM_PARA(1)
        SLOPE_PARA.XCU=IFLAG_OPTIM_PARA(1)
    ENDIF
    

    
  ! Define search boundaries
    allocate(lower(SLOPE_PARA.NDIM), upper(SLOPE_PARA.NDIM))
    allocate(X(SLOPE_PARA.NSLICE+1,SLOPE_PARA.POPSIZE),Y(SLOPE_PARA.NSLICE+1,SLOPE_PARA.POPSIZE),RDF(SLOPE_PARA.NDIM,SLOPE_PARA.POPSIZE))
    lower =0.D0
    upper =1.D0
    
    CALL PSO_SLIP.INIT(SLOPE_PARA.NSLICE+1,SLOPE_PARA.POPSIZE) !SAVE TOPTENSLIP
    
   CALL GEN_SLIP_SAMPLE(SLOPE_PARA.XTL,SLOPE_PARA.XTU,SLOPE_PARA.XCL,SLOPE_PARA.XCU, &
        SLOPE_PARA.NSLICE,SLOPE_PARA.popsize,SLOPE_PARA.GX,SLOPE_PARA.NGX,  & 
        SLOPE_PARA.RX,SLOPE_PARA.NRX,X,Y,RDF,SLOPE_PARA.SHAPE)
   PSO_SLIP.XSTART(1,:,:)=X;PSO_SLIP.XSTART(2,:,:)=Y
    
  ! Initialize evolutionary optimizer
  PSO_SLIP.NFOSCAL=0      
  ea = Evolutionary(SLOPE_FOS_CAL, lower, upper, eps1=slope_para.eps1,eps2=slope_para.eps2,&
                    eps3=slope_para.eps3,popsize = SLOPE_PARA.popsize, max_iter = SLOPE_PARA.max_iter)

  ! Optimize
  call ea % optimize(solver =SLOPE_PARA.SOLVER ,XSTART=TRANSPOSE(RDF),&
                    w=SLOPE_PARA.W, c1=SLOPE_PARA.C1, c2=SLOPE_PARA.C2, gamma=SLOPE_PARA.GAMMA, &
                    F=SLOPE_PARA.F, CR=SLOPE_PARA.CR, strategy=SLOPE_PARA.strategy)
 

  
  CALL PSO_SLIP.SORT()
  !FM1=0.;FA1=0.
  !CALL POLYLINE_FOS_CAL(PSO_SLIP.X(:,:,PSO_SLIP.A2Z(1)),SLOPE_PARA.NDIM,SLOPE_PARA.DL,FM1,FA1)
  !PRINT *,'FM1,FA1,FOS=',FM1,FA1,ABS(FA1/FM1)
  IF(IFLAG_OPTIM/=1) THEN
      
      CALL PSO_SLIP.PLOT()
      CALL PSO_SLIP.EXPORT(EA)
  
   ENDIF
   
  ! Display results
  call ea % print()
  
  !IF(.NOT.ALLOCATED(XSLIP1)) ALLOCATE(XSLIP1(2,SLOPE_PARA.NDIM))
  !CALL UNKNOWN2SLIPS(EA.XOPT,XSLIP1,ISERROR1)
  !IF(ISERROR1==0) THEN
  !  DO I=1,SLOPE_PARA.NDIM
  !      PRINT *, XSLIP1(:,I)
  !  ENDDO
  !ELSE
  !   PRINT *, 'CANNOT RECOVERED SLIP SURFACE FROM THE FINAL SOLUTION.'  
  !ENDIF
  
!  IF(ALLOCATED(SLOPE_PARA.GX)) DEALLOCATE(SLOPE_PARA.GX)
!  IF(ALLOCATED(SLOPE_PARA.RX)) DEALLOCATE(SLOPE_PARA.RX)
  !IF(ALLOCATED(XSI)) DEALLOCATE(XSI)
  !IF(ALLOCATED(YSI)) DEALLOCATE(YSI)
  IF(ALLOCATED(LOWER)) DEALLOCATE(LOWER)
  IF(ALLOCATED(UPPER)) DEALLOCATE(UPPER)
  IF(ALLOCATED(X)) DEALLOCATE(X,Y,RDF)
  IF(ALLOCATED(XSLIP1)) DEALLOCATE(XSLIP1)
ENDSUBROUTINE

REAL(8) FUNCTION SLOPE_FOS_CAL(UNS)
!GIVEN:
!THE VARIABLES: UNKNOWNS=UNS
    USE,INTRINSIC :: ISO_FORTRAN_ENV
    real(kind = 8), dimension(:), intent(in) :: UNS
    REAL(8)::XT1(2),XC1(2),FM1,FA1,T1
    REAL(8)::X1(2,SLOPE_PARA.NSLICE+1)
    
    INTEGER::ISERROR=0,N1,N2
    INTEGER::IT1=0
    
    PSO_SLIP.NFOSCAL=PSO_SLIP.NFOSCAL+1
    IT1=PSO_SLIP.NFOSCAL

       
    CALL UNKNOWN2SLIPS(UNS,X1,ISERROR)
    IF(ISERROR/=0) THEN
        PRINT *, "WARNING IN SLOPE_FOS_CAL.FAIL TO GET THE SLIP SURFACE WITH THE SOLUTION"
        SLOPE_FOS_CAL=1.D10
        FA1=-999;FM1=1
    ELSE
        CALL POLYLINE_FOS_CAL(X1,SLOPE_PARA.NSLICE+1,SLOPE_PARA.DL,FM1,FA1)
        SLOPE_FOS_CAL=ABS(FA1/FM1)
    ENDIF
    
    
    N2=MINLOC(PSO_SLIP.FOS,DIM=1)
    IF(PSO_SLIP.FOS(N2)>SLOPE_FOS_CAL) THEN
        T1=PSO_SLIP.FOS(N2)-SLOPE_FOS_CAL
        PRINT 110,EA_ITER,IT1,PSO_SLIP.FOS(N2),SLOPE_FOS_CAL,T1,T1/PSO_SLIP.FOS(N2)*100    
    ENDIF
    
    N1=MOD(IT1-1,PSO_SLIP.NSLIP)+1
    !IF(N1==0) N1=PSO_SLIP.NSLIP
    
    IF(PSO_SLIP.FOS(N1)>SLOPE_FOS_CAL) THEN
        PSO_SLIP.FOS(N1)=SLOPE_FOS_CAL
        PSO_SLIP.X(:,:,N1)=X1
        PSO_SLIP.FA(N1)=FA1;PSO_SLIP.FM(N1)=FM1;
    ENDIF
!    IF(MOD(IT1,SLOP_PARA.NDIM*10)==0) PRINT 100, MINVAL(PSO_SLIP.FOS,DIM=1)
    
    !write(*,100,ADVANCE='NO') achar(13),MINVAL(PSO_SLIP.FOS,DIM=1)
    !CALL FLUSH(6)
    !CALL SLEEP(0.001)
100 FORMAT('TILL NOW,THE MINIMAL FOS=',F7.4)
110 FORMAT('NITER=',I7,1X,'NSLIPTRIAL=',I7,1X,'SOLUTION WAS IMPROVED./OLD/NEW/DIF./%/=',4('/',F7.4),'%/')
END FUNCTION

SUBROUTINE UNKNOWN2SLIPS(UNS,X1,ISERROR)
!GIVEN:
!OPTIMIZATION VARIALBES:UNS
!RETURN:
!THE CORRESPONDING SLIP SURFACE:X1(2,:)
!ERROR STATUS:ISERROR. =0,NO ERROR. 

    real(kind = 8), dimension(:), intent(IN) :: UNS
    REAL(8),INTENT(OUT)::X1(:,:)
    INTEGER,INTENT(OUT)::ISERROR
    REAL(8)::XT1(2),XC1(2),UNS1(SIZE(UNS))
    INTEGER::IFLAG=1,N1
    
    UNS1=UNS
    IF(SLOPE_PARA.SHAPE==0) THEN
        N1=SLOPE_PARA.NDIM
    ELSE
        N1=2
    ENDIF
    CALL GET_ENTRY_EXIT([UNS(1),UNS(N1)],XT1,XC1)
    IFLAG=1;ISERROR=0
    CALL GEN_ADMISSIBLE_SLIP_SURFACE(XT1,XC1,SLOPE_PARA.NSLICE,SLOPE_PARA.GX,SLOPE_PARA.NGX,&
        SLOPE_PARA.RX,SLOPE_PARA.NRX,UNS1,IFLAG,X1(1,:),X1(2,:),ISERROR,ISHAPE=SLOPE_PARA.SHAPE)
    
ENDSUBROUTINE


SUBROUTINE GET_ENTRY_EXIT(RAD1,XT1,XC1)
    REAL(8),INTENT(IN)::RAD1(2)
    REAL(8),INTENT(OUT)::XT1(2),XC1(2)
    
    XT1(1)=SLOPE_PARA.XTL+(SLOPE_PARA.XTU-SLOPE_PARA.XTL)*RAD1(1)
    XT1(2)=INTERPOLATION(SLOPE_PARA.GX(1,:),SLOPE_PARA.GX(2,:),SLOPE_PARA.NGX,XT1(1))
    XC1(1)=SLOPE_PARA.XCL+(SLOPE_PARA.XCU-SLOPE_PARA.XCL)*RAD1(2)
    IF(IFLAG_OPTIM/=1) THEN
        XC1(2)=INTERPOLATION(SLOPE_PARA.GX(1,:),SLOPE_PARA.GX(2,:),SLOPE_PARA.NGX,XC1(1))
    ELSE
        XC1(2)=IFLAG_OPTIM_PARA(2) !STREAMLINE REPAIR
    ENDIF
END SUBROUTINE

SUBROUTINE FILEREADIN(SELF,FILE)
    USE IFPORT
    CLASS(SLOPE_PSO_PARAR)::SELF
    CHARACTER(1024)::FILE
    INTEGER::UNIT,LENGTH

    !EXTERNAL SLOPE_PARA_PARSER
    UNIT=11
    OPEN(UNIT,FILE=FILE,STATUS='OLD')
    
    inquire(UNIT,name=PARAFILE)
    length = SPLITPATHQQ(PARAFILE, drive, dir, PARAFILENAME, ext)
    
    CALL READ_FILE(SELF,UNIT,SLOPE_PARA_PARSER)
    CLOSE(UNIT)


END SUBROUTINE

subroutine READ_FILE(INDATA,UNIT,COMMAND_READ)
        
	implicit none
    TYPE(SLOPE_PSO_PARAR)::INDATA
    INTEGER,INTENT(IN)::UNIT   
	EXTERNAL::COMMAND_READ
    integer::ef,iterm,i,strL,N1
	parameter(iterm=1024)
	character(iterm)::term,term2
	character(1)::ch
                
	TYPE(COMMAND_TYDEF)::COMMAND
        
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
			if(index(term2,'/')/=0) then
				strL=index(term2,'/')-1
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

			call lowcase(term)
			call COMMAND.PARSER(term)			
			!term=adjustl(trim(term))
			call COMMAND_READ(INDATA,COMMAND,unit)			 	
		end if
	end do


	
999	format(a<iterm>)

end subroutine  

SUBROUTINE SLOPE_PARA_PARSER(INDATA,COMMAND,UNIT)
        
    TYPE(SLOPE_PSO_PARAR):: INDATA
    TYPE(COMMAND_TYDEF),INTENT(IN)::COMMAND
    INTEGER,INTENT(IN)::UNIT
    INTEGER::I,NC1,NC2,IZONE,N1
        
    !SELECT TYPE(INDATA)
    !CLASS IS(SLOPE_PSO_PARAR)
    
        SELECT CASE(TRIM(ADJUSTL(COMMAND.KEYWORD)))
	    case('ge','gx')
		    print *, 'Reading SLOPE GROUND SURFACE data...'
		    do i=1, COMMAND.NOPT
			    select case(COMMAND.OPTION(i).NAME)
				    case('num')
					    INDATA.NGX=COMMAND.OPTION(i).VALUE
				    case default
					    call Err_msg(COMMAND.OPTION(i).name)
			    end select
		    end do
            IF(ALLOCATED(INDATA.GX)) DEALLOCATE(INDATA.GX)
            ALLOCATE(INDATA.GX(2,INDATA.NGX))
            READ(UNIT,*)  INDATA.GX
        CASE('rx','re')
            print *, 'Reading SLOPE BEDROCK data...'
		    do i=1, COMMAND.NOPT
			    select case(COMMAND.OPTION(i).NAME)
				    case('num')
					    INDATA.NRX=COMMAND.OPTION(i).VALUE
				    case default
					    call Err_msg(COMMAND.OPTION(i).name)
			    end select
		    end do
            IF(ALLOCATED(INDATA.RX)) DEALLOCATE(INDATA.RX)
            ALLOCATE(INDATA.RX(2,INDATA.NRX))
            READ(UNIT,*)  INDATA.RX
        CASE('wsp','ws','wsb')
            print *, 'Reading Weak soil bottom data...'
		    do i=1, COMMAND.NOPT
			    select case(COMMAND.OPTION(i).NAME)
				    case('num')
					    INDATA.NWSP=COMMAND.OPTION(i).VALUE
				    case default
					    call Err_msg(COMMAND.OPTION(i).name)
			    end select
		    end do
            IF(ALLOCATED(INDATA.wsp)) DEALLOCATE(INDATA.wsp)
            ALLOCATE(INDATA.WSP(2,INDATA.nwsp))
            READ(UNIT,*)  INDATA.wsp            
        CASE('x_e&e','xlimit','x_entry_and_exit','x_search_limits')
            print *, 'Reading SLOPE ENTRY AND EXIT LIMITS data...'
		    do i=1, COMMAND.NOPT
			    select case(COMMAND.OPTION(i).NAME)
				    case('xtu','xtr')
					    INDATA.xtu=COMMAND.OPTION(i).VALUE
				    case('xtl')
					    INDATA.xtl=COMMAND.OPTION(i).VALUE
				    case('xcu','xcr')
					    INDATA.xcu=COMMAND.OPTION(i).VALUE
				    case('xcl')
					    INDATA.xcl=COMMAND.OPTION(i).VALUE 
                   
				    case default
					    call Err_msg(COMMAND.OPTION(i).name)
			    end select
		    end do
        CASE('repair','streamline_repair')
            print *, 'Reading STREAMLINE ADMISSIBLITY REPAIR data...'
            DO i=1, COMMAND.NOPT
			    select case(COMMAND.OPTION(i).NAME)
				    case('direction','zone')
					    INDATA.direction=int(COMMAND.OPTION(i).VALUE)
				    case('a')                    
					    INDATA.a=COMMAND.OPTION(i).VALUE
                    case('b')
                        INDATA.b=COMMAND.OPTION(i).VALUE
                    case('c')
                        INDATA.c=COMMAND.OPTION(i).VALUE
				    case default
					    call Err_msg(COMMAND.OPTION(i).name)                        
                end select
            ENDDO
        CASE('ssgm')
            print *, 'Reading initial slip surface generation method data...'
            DO i=1, COMMAND.NOPT
			    select case(COMMAND.OPTION(i).NAME)
				    case('method')
					    INDATA.ssgm=int(COMMAND.OPTION(i).VALUE)
				    case default
					    call Err_msg(COMMAND.OPTION(i).name)                        
                end select
            ENDDO 
            IF(INDATA.ssgm==2) THEN
                READ(UNIT,*)  INDATA.NWSP
                IF(ALLOCATED(INDATA.WSP)) DEALLOCATE(INDATA.WSP)
                ALLOCATE(INDATA.WSP(2,INDATA.NWSP))
                READ(UNIT,*)  INDATA.WSP
            ENDIF
        CASE('optims','optimpara')
            print *, 'Reading SLOPE OPTIMPARAMETERS data...'
            N1=0
		    do i=1, COMMAND.NOPT
			    select case(COMMAND.OPTION(i).NAME)
				    case('solver')
					    INDATA.solver=trim(adjustl(COMMAND.OPTION(i).CVALUE))
				    case('p_eps','eps1')
					    INDATA.eps1=COMMAND.OPTION(i).VALUE
				    case('eps2','f_eps')
					    INDATA.eps2=COMMAND.OPTION(i).VALUE
				    case('c1')
					    INDATA.c1=COMMAND.OPTION(i).VALUE 
				    case('c2')
					    INDATA.c2=COMMAND.OPTION(i).VALUE   
				    case('w')
					    INDATA.w=COMMAND.OPTION(i).VALUE
                    CASE('f')
                        INDATA.f=COMMAND.OPTION(i).VALUE
                    CASE('cr')
                        INDATA.cr=COMMAND.OPTION(i).VALUE
                    case('eps3')
                        INDATA.eps3=COMMAND.OPTION(i).VALUE
				    case('gamma')
					    INDATA.gamma=COMMAND.OPTION(i).VALUE
                    case('popsize')
                        INDATA.popsize=int(COMMAND.OPTION(i).VALUE)
                        N1=1
                    case('max_iter','maxiter')
                        INDATA.max_iter=int(COMMAND.OPTION(i).VALUE)
                    case('nslice')
                        INDATA.nslice=int(COMMAND.OPTION(i).VALUE)                        
                    case('dl')
                        INDATA.dl=COMMAND.OPTION(i).VALUE
                    case('strategy')
                        INDATA.strategy=trim(adjustl(COMMAND.OPTION(i).CVALUE))
                    CASE('sigma')
                        INDATA.sigma=COMMAND.OPTION(i).VALUE
                    CASE('mu_perc')
                        INDATA.mu_perc=COMMAND.OPTION(i).VALUE
                    CASE('shape')
                        INDATA.shape=int(COMMAND.OPTION(i).VALUE)
				    case default
					    call Err_msg(COMMAND.OPTION(i).name)
			    end select
		    end do
            
            IF(INDATA.shape==0) THEN
                INDATA.ndim=INDATA.nslice+1                
            ELSE
                INDATA.ndim=3 !CIRCULAR SLIP SURFACE
            ENDIF
            IF(N1==0) INDATA.popsize=max(INDATA.ndim*2,60)
        CASE DEFAULT
             call Err_msg(COMMAND.KEYWORD)
        END SELECT
    !END SELECT 
END SUBROUTINE


SUBROUTINE PSO_SLIP_INIT(SELF,NX,NSLIP)
    CLASS(PSO_SLIP_INFO)::SELF
    INTEGER,INTENT(IN)::NSLIP,NX
    
    SELF.NX=NX;SELF.NSLIP=NSLIP
    IF(ALLOCATED(SELF.FOS)) DEALLOCATE(SELF.FOS,SELF.X,SELF.A2Z,SELF.FA,SELF.FM,SELF.XSTART,SELF.ISSHOW,SELF.ISSHOW_XSTART,SELF.COLOR)
    ALLOCATE(SELF.FOS(NSLIP),SELF.X(2,NX,NSLIP),SELF.A2Z(NSLIP),SELF.FA(NSLIP),SELF.FM(NSLIP), &
            SELF.XSTART(2,NX,NSLIP),SELF.ISSHOW(NSLIP),SELF.ISSHOW_XSTART(NSLIP),SELF.COLOR(NSLIP))
    
    SELF.FOS=1.D30;SELF.X=0.D0;SELF.XSTART=0.D0;SELF.ISSHOW=.FALSE.;SELF.ISSHOW_XSTART=.FALSE.;
    
END SUBROUTINE


SUBROUTINE PSO_SLIP_SORT(SELF)
    USE quicksort
    CLASS(PSO_SLIP_INFO)::SELF
    REAL(8)::FOS1(SIZE(SELF.FOS))
    
    FOS1=SELF.FOS;SELF.A2Z=[1:SELF.NSLIP]
    CALL QUICK_SORT(FOS1,SELF.A2Z)
    
    SELF.ISSHOW(SELF.A2Z(1))=.TRUE.
    
END SUBROUTINE

SUBROUTINE PSO_SLIP_EXPORT(SELF,EA)
    USE IFPORT
    CLASS(PSO_SLIP_INFO)::SELF
    type(Evolutionary) :: ea
    INTEGER::I,N1,N2
    CHARACTER(8):: CHAR_TIME
    IF(SELF.NX*SELF.NSLIP==0) RETURN
    
    OPEN(11,FILE=TRIM(PARAFILENAME)//'_SLIP_OPTIM_BY_EA.TXT',STATUS='UNKNOWN',POSITION='APPEND')
    N2=LEN(TRIM(PARAFILENAME)//TRIM(EXT))
    
    N1=SELF.A2Z(1)
    CALL TIME(CHAR_TIME)
    WRITE(11,100) SELF.NX,SELF.FOS(N1),PSO_SLIP.NFOSCAL,SLOPE_PARA.SOLVER,TRIM(PARAFILENAME)//TRIM(EXT),DATE(),CHAR_TIME
    
    DO I=1,SELF.NX
        WRITE(11,110) SELF.X(:,I,N1)
    ENDDO
    
      ! Display results
    WRITE(11,120)
    call ea % print(11)
    WRITE(11,130)
    
    CLOSE(11) 
    
100 FORMAT(I5,1X,'FOS=',F7.4,1X,'NFOSCAL=',I7,1X,'SOLVER=',A5,1X,'INPUTFILE=',A<N2>,1X,'TIME=',A8,1X,'DATE=',A8) 
110 FORMAT(2(F8.3,1X))
120 FORMAT('//OPTIMIZATION PARAMETERS:BEGIN')
130 FORMAT('//OPTIMIZATION PARAMETERS:END')
END SUBROUTINE

SUBROUTINE PSO_SLIP_PLOT(SELF)
    use function_plotter
    implicit none
    CLASS(PSO_SLIP_INFO)::SELF
    
    integer :: i,j,k,n1,DIRECTION1=1,AI1(SELF.NSLIP)
    REAL(8)::DX1,DY1,DEG1,T1,SCALE1,PPM1,FS1,DEG2,MAX1
    CHARACTER(16)::STR1
    REAL(GLFLOAT)::COLOR1(4),COLOR2(4)
 
    call glDeleteLists(PSO_SLIP_PLOT_LIST, 1_glsizei)
    IF(SELF.NSLIP==0.OR.SELF.NX==0) RETURN
    
    call reset_view    
    call glNewList(PSO_SLIP_PLOT_LIST, gl_compile_and_execute)
 
    call glPolygonMode(gl_front_and_back, gl_fill)
	call gldisable(GL_CULL_FACE);  
    !MAX1=MAXVAL(STREAMLINE.SF_SLOPE,MASK=STREAMLINE.SF_SLOPE<10)
    
    !AI1(SELF.A2Z)=[1:SELF.NSLIP]
    !I=SELF.A2Z(1)
    PSO_SLIP.COLOR=BLUE
    PSO_SLIP.COLOR(PSO_SLIP.A2Z(1:10))=COLOR_TOP_TEN
    DO I=1,SELF.NSLIP
        IF(.NOT.SELF.ISSHOW(I)) CYCLE
        !CALL glEnable(GL_LINE_STIPPLE)
        !!
        !CALL glLineStipple(1,0xAAAA)
        
        CALL glLineWidth(3.0_glfloat)
 
        COLOR1=mycolor(:,SELF.COLOR(I))
       
	    call glBegin(gl_LINE_STRIP)
		DO J=1,SELF.NX
            COLOR2=COLOR1
            CALL glcolor4fv(COLOR2)
			call glvertex3dv([SELF.X(:,J,I),0.D0]) 
		ENDDO
	    CALL GLEND()
        
        
        DX1=SELF.X(1,SELF.NX,I)-SELF.X(1,SELF.NX-1,I) 
        DY1=SELF.X(2,SELF.NX,I)-SELF.X(2,SELF.NX-1,I) 
 
        IF(ABS(DX1)>1E-7) THEN
            DEG1=ATAN(DY1/DX1)/PI*180.
        ELSE
            DEG1=SIGN(PI/2.0,DY1)/PI*180.
        ENDIF
            
        IF(DEG1<0) DEG1=DEG1+180.
 
        DX1=SELF.X(1,1,I)-SELF.X(1,2,I) 
        DY1=SELF.X(2,1,I)-SELF.X(2,2,I) 
        !DEG1=ASIN(DY1/T1)/PI*180.0 
        IF(ABS(DX1)>1E-7) THEN
            DEG2=ATAN(DY1/DX1)/PI*180.
        ELSE
            DEG2=SIGN(PI/2.0,DY1)/PI*180.
        ENDIF
            
        IF(DEG2<0) DEG2=DEG2+180.
            
        T1=I
        N1=0
        DO WHILE(T1>1.D0)
            T1=T1/10
            N1=N1+1
        ENDDO
        
        WRITE(STR1,110) SELF.FOS(I),I
            
        !PPM1=glutget(GLUT_SCREEN_WIDTH)/REAL(glutget(GLUT_SCREEN_WIDTH_MM)) !PIXELS PER MM
        !10 pound 
        FS1=glutStrokeWidth(GLUT_STROKE_ROMAN,ICHAR("X")) 
        SCALE1=POSDATA.MODELR/80/FS1*stroke_fontsize
        !scale1=PPM1*3.527777778/119.05*0.02*stroke_fontsize
        call glLineWidth(1.0_glfloat)
            
        CALL drawStrokeText(DEG1,&
                            SELF.X(1,SELF.NX,I), &
                            SELF.X(2,SELF.NX,I),0.0, &
                            scale1,&
                            STR1)
        CALL drawStrokeText(DEG2,&
                            SELF.X(1,1,I), &
                            SELF.X(2,2,I),0.0, &
                            scale1,&
                            STR1)                                
    ENDDO
    
    !PLOT XSTART
    !PSO_SLIP.COLOR=[MOD(I,)]
    DO I=1,SELF.NSLIP
        IF(.NOT.SELF.ISSHOW_XSTART(I)) CYCLE
        
        !CALL glEnable(GL_LINE_STIPPLE)
        !!!
        !CALL glLineStipple(1,0xAAAA)
        
        CALL glLineWidth(3.0_glfloat)
 
        COLOR1=DISTINCT_COLOR(:,MOD(I-1,44)+1)
       
	    call glBegin(gl_LINE_STRIP)
		DO J=1,SELF.NX
            COLOR2=COLOR1
            CALL glcolor4fv(COLOR2)
			call glvertex3dv([SELF.XSTART(:,J,I),0.D0]) 
		ENDDO
	    CALL GLEND()
        
        
        DX1=SELF.XSTART(1,SELF.NX,I)-SELF.XSTART(1,SELF.NX-1,I) 
        DY1=SELF.XSTART(2,SELF.NX,I)-SELF.XSTART(2,SELF.NX-1,I) 
 
        IF(ABS(DX1)>1E-7) THEN
            DEG1=ATAN(DY1/DX1)/PI*180.
        ELSE
            DEG1=SIGN(PI/2.0,DY1)/PI*180.
        ENDIF
            
        IF(DEG1<0) DEG1=DEG1+180.
 
        DX1=SELF.XSTART(1,1,I)-SELF.XSTART(1,2,I) 
        DY1=SELF.XSTART(2,1,I)-SELF.XSTART(2,2,I) 
        !DEG1=ASIN(DY1/T1)/PI*180.0 
        IF(ABS(DX1)>1E-7) THEN
            DEG2=ATAN(DY1/DX1)/PI*180.
        ELSE
            DEG2=SIGN(PI/2.0,DY1)/PI*180.
        ENDIF
            
        IF(DEG2<0) DEG2=DEG2+180.
            
        T1=I
        N1=0
        DO WHILE(T1>1.D0)
            T1=T1/10
            N1=N1+1
        ENDDO
        
        WRITE(STR1,100) I
            
        !PPM1=glutget(GLUT_SCREEN_WIDTH)/REAL(glutget(GLUT_SCREEN_WIDTH_MM)) !PIXELS PER MM
        !10 pound 
        FS1=glutStrokeWidth(GLUT_STROKE_ROMAN,ICHAR("X")) 
        SCALE1=POSDATA.MODELR/80/FS1*stroke_fontsize
        !scale1=PPM1*3.527777778/119.05*0.02*stroke_fontsize
        call glLineWidth(1.0_glfloat)
            
        CALL drawStrokeText(DEG1,&
                            SELF.XSTART(1,SELF.NX,I), &
                            SELF.XSTART(2,SELF.NX,I),0.0, &
                            scale1,&
                            STR1)
        CALL drawStrokeText(DEG2,&
                            SELF.XSTART(1,1,I), &
                            SELF.XSTART(2,2,I),0.0, &
                            scale1,&
                            STR1)                                
    ENDDO    
    
    CALL glLineWidth(1.0_glfloat)
    call glEndList
 
    
    call glutPostRedisplay

100 FORMAT('P=',I<N1>)
110 FORMAT(F7.3,'/',I<N1>)

ENDSUBROUTINE

SUBROUTINE PSO_SLIP_PLOT_SETTING(FLAG_EA_SHOW)
    
    implicit none
    !CLASS(PSO_SLIP_INFO)::SELF
    integer,INTENT(INOUT)::FLAG_EA_SHOW
    

    
    SELECTCASE(FLAG_EA_SHOW)
    CASE(SEARCH_BY_EA)
        CALL SLOPE_OPTIM()
    CASE(EA_SHOW_ALL)
        PSO_SLIP.ISSHOW=.TRUE.
        
    CASE(EA_SHOW_MIN)
        PSO_SLIP.ISSHOW=.FALSE.
        PSO_SLIP.ISSHOW(PSO_SLIP.A2Z(1))=.TRUE.
        
    CASE(EA_SHOW_TOPTEN)
        PSO_SLIP.ISSHOW=.FALSE.
        PSO_SLIP.ISSHOW(PSO_SLIP.A2Z(1:10))=.TRUE.
    CASE(EA_SHOW_XSTART)
        PSO_SLIP.ISSHOW_XSTART=.TRUE.
    CASE(EA_SHOW_NONE)
        PSO_SLIP.ISSHOW=.FALSE. 
    CASE(EA_SHOW_XSTART_NONE)
        PSO_SLIP.ISSHOW_XSTART=.FALSE.        
    ENDSELECT
    
    CALL PSO_SLIP.PLOT()
    
ENDSUBROUTINE

INTEGER FUNCTION PSO_SLIP_PLOT_MENU()
USE OPENGL_GLUT
!CLASS(PSO_SLIP_INFO)::SELF

PSO_SLIP_PLOT_MENU=glutCreateMenu(PSO_SLIP_PLOT_SETTING)

CALL glutAddMenuEntry("SearchByEA...",SEARCH_BY_EA)
call glutAddMenuEntry("ShowTopTenParticle",EA_SHOW_TOPTEN)
call glutAddMenuEntry("ShowMinimalParticle",EA_SHOW_MIN)
call glutAddMenuEntry("ShowAllParticles",EA_SHOW_ALL)
call glutAddMenuEntry("ShowNoneParticles",EA_SHOW_NONE)
call glutAddMenuEntry("ShowAllInitialParticles",EA_SHOW_XSTART)
call glutAddMenuEntry("ShowNoneInitialParticles",EA_SHOW_XSTART_NONE)

ENDFUNCTION



SUBROUTINE GEN_SLIP_SAMPLE(XTL,XTU,XCL,XCU,NSEG,NSLIP,GX,NGX,RX,NRX,X,Y,RDF,ISHAPE)

!GIVEN:
!THE LIMITS OF THE SLOPE TOE: XTL,XTU
!THE LIMITS OF THE SLOPE CREST: XCL,XCU
!SAMPLE NUMBER: NSLIP
!GROUND SURFACE: GX(2,NGX)
!BEDROCK SURFACE: RX(2,NRX)
!NUMBER OF SUBDIVISION OF EACH SLIP LINE: NSEG>1
!SLIP SURFACE SHAPE: ISHAPE,=0,NONCIRCULAR(BYDEFAULT),=1,CIRCULAR
!RETURN:
!THE SLIP NODES: X(NSEG+1,NSLIP),Y(NSEG+1,NSLIP)
!(RADOM) FACTOR: RDF(NSEG+1)
    
    IMPLICIT NONE

    INTEGER,INTENT(IN)::NSEG,NSLIP,NGX,NRX
    INTEGER,OPTIONAL,INTENT(IN)::ISHAPE
    REAL(8),INTENT(IN)::XTL,XTU,XCL,XCU,GX(2,NGX),RX(2,NRX)
    REAL(8),INTENT(OUT)::X(NSEG+1,NSLIP),Y(NSEG+1,NSLIP),RDF(:,:)
    

    INTEGER::I,J,ISERROR=1,NITER=0,IFLAG=0,NITER2=0,ISHAPE1
    REAL(8)::R1,R2,XT(2),XC(2),YL(NSEG+1),YU(NSEG+1)

    IF(PRESENT(ISHAPE)) THEN
        ISHAPE1=ISHAPE
    ELSE
        ISHAPE1=0
    ENDIF
    
    CALL RANDOM_SEED()
    J=1;NITER2=0
    DO WHILE (J<=NSLIP)
        call random_number(R1)
        call random_number(R2)
        
        XT(1)=XTL+(XTU-XTL)*R1
        XC(1)=XCL+(XCU-XCL)*R2
        XT(2)=INTERPOLATION(GX(1,:),GX(2,:),NGX,XT(1))
        IF(IFLAG_OPTIM/=1) THEN
            XC(2)=INTERPOLATION(GX(1,:),GX(2,:),NGX,XC(1))
        ELSE
            XC(2)=IFLAG_OPTIM_PARA(2)
            R2=1.d0
        ENDIF
        ISERROR=1;NITER=0;IFLAG=0
        DO WHILE(ISERROR/=0.AND.NITER<10)
            CALL GEN_ADMISSIBLE_SLIP_SURFACE(XT,XC,NSEG,GX,NGX,RX,NRX,RDF(:,J),IFLAG,X(:,J),Y(:,J),ISERROR,YL,YU,ISHAPE1)
            RDF(1,J)=R1;
            IF(ISHAPE1==0) THEN
                RDF(NSEG+1,J)=R2
            ELSE
                RDF(2,J)=R2
            ENDIF
            IF(ISERROR/=0) THEN
                NITER=NITER+1
            ENDIF            
        END DO
        IF(ISERROR==0) THEN            
            J=J+1
            NITER2=0
        ELSE
            NITER2=NITER2+1
            IF(NITER2>100) THEN
                PRINT *,"ERROR IN GEN_SLIP_SAMPLE. FAIL TO GENERATE J SLIPS.J=",J
                STOP
            ENDIF
        ENDIF
        
    ENDDO

ENDSUBROUTINE
    
SUBROUTINE GEN_ADMISSIBLE_SLIP_SURFACE(XT,XC,NSEG,GX,NGX,RX,NRX,RDF,IFLAG,X,Y,ISERROR,YL,YU,ISHAPE)
!REFERRENC: Cheng YM, Li L, Chi SC. Performance studies on six heuristic global optimization methods in the
! location of critical slip surface. Computers and Geotechnics. 2007;34(6):462-484.
!GIVEN:
!THE TOE AND CREST POINT LOCATIONS OF THE SLIP SURFACE: XT(2),XC(2)
!NUMBER OF SUBDIVISION: NSEG>1
!GROUND SURFACE: GX(2,NGX)
!BEDROCK SURFACE: RX(2,NRX)
!(RADOM) FACTOR: RDF(NSEG+1),ONLY RDF(2:NSEG) IS USED.
!IFLAG:
!=1 USE RDF TO GENERATE A SLIP SURFACE
!/=1, USE RADOM FACTOR TO GENERATE A SLIP SURFACE, AND RETURN THE FACTOR IN RDF.
!ISHAPE: SHAPE OF SLIP SURFACE, =0(BYDEFAULT) NONCIRCULAR, =1, CIRCULAR
!RETURN: 
!THE X COMPONENT OF THE VETEX OF THE SLIP LINE:X(NSEG+1)
!THE Y COMPONENT : Y(NSEG+1)
!ISERROR=0,SUCCESSFULLY GENERATE A SLIP SURFACE.
!/=0, NO ADMISSIBLE SLIP SURFACE FOUND.
IMPLICIT NONE
INTEGER,INTENT(IN)::NGX,NRX,NSEG,IFLAG
REAL(8),INTENT(IN)::XT(2),XC(2),GX(2,NGX),RX(2,NRX)
REAL(8),INTENT(INOUT)::RDF(:)
REAL(8),INTENT(OUT)::X(NSEG+1),Y(NSEG+1)
REAL(8),OPTIONAL,INTENT(OUT)::YU(NSEG+1),YL(NSEG+1) 
INTEGER,INTENT(OUT)::ISERROR
INTEGER,OPTIONAL,INTENT(IN)::ISHAPE

INTEGER::I,J,NX,ISHAPE1
REAL(8)::DX1,YU1,YL1,K1,YCL1,T1
REAL(8)::RMIN,RMAX,CENTER(2,2),R1,CENT1(2),RA1(7)

!REAL(8),EXTERNAL::INTERPOLATION

ISERROR=0

IF(ABS(XT(1)-XC(1))<1.D-7) THEN
    ISERROR=1 
    RETURN
ENDIF

IF(NSEG==1) THEN
    X=[XT(1),XC(1)];
    Y=[XT(2),XC(2)];
    RETURN
ENDIF

IF(PRESENT(ISHAPE)) THEN
    ISHAPE1=ISHAPE
ELSE
    ISHAPE1=0
ENDIF



NX=NSEG+1
DX1=(XC(1)-XT(1))/NSEG

!ALLOCATE(X(NX),YU(NX),YL(NX))

DO I=1,NX
    X(I)=XT(1)+DX1*(I-1)
ENDDO
!add characteristic slope surface point
DO I=1,NGX
    IF(MIN(XT(1),XC(1))<GX(1,I).AND.MAX(XT(1),XC(1))>GX(1,I)) THEN
        X(MINLOC(ABS(X(2:NX-1)-GX(1,I)),DIM=1)+1)=GX(1,I)    
    ENDIF
ENDDO

Y(1)=XT(2);Y(NX)=XC(2)
IF(PRESENT(YU)) THEN
YU(1)=Y(1);YU(NX)=Y(NX)
ENDIF
IF(PRESENT(YL)) THEN
    YL(1)=Y(1);YL(NX)=Y(NX)
ENDIF

IF(IFLAG/=1) THEN
    call random_number(RDF)
ENDIF

IF(ISHAPE/=0) THEN
    CALL GEN_CIRCULAR_SLIP(XT(1),XC(1),RMIN,RMAX,CENTER)
    IF(RMIN>RMAX) THEN
        ISERROR=1
        PRINT *, "UNADMISSIBILITY OCCURS IN GEN_ADMISSIBLE_SLIP_SURFACE. RMAX<RMIN" 
        RETURN
    ENDIF
    !R1=RMIN+RDF(3)*(RMAX-RMIN)
    CENT1=CENTER(:,1)+RDF(3)*(CENTER(:,2)-CENTER(:,1))
    R1=NORM2(CENT1-XT)
    !IF(ABS(R1-NORM2(CENT1-XC))>0.01) THEN
    !    PRINT *, 'UNEXPECTED 1'
    !    PAUSE
    !ENDIF
ENDIF


DO I=2,NSEG
    
    
    YL1=INTERPOLATION(RX(1,:),RX(2,:),NRX,X(I))    
    YU1=INTERPOLATION(GX(1,:),GX(2,:),NGX,X(I))
    
    IF(ISHAPE1/=0) THEN
        RA1=CIRCLE_LINE_INTERSECTION(CENT1,R1,[X(I),YU1],[X(I),-1.D10])
        IF(NINT(RA1(1))>1.1) THEN
            IF(RA1(4)<=1.D0.AND.RA1(4)>=0.D0) THEN
                Y(I)=RA1(3)
            ELSEIF(RA1(7)<=1.D0.AND.RA1(7)>=0.D0) THEN
                Y(I)=RA1(6)
            ELSE
                PRINT *, 'SOMETHING WRONG IN GEN_ADMISSIBLE_SLIP_SURFACE.IT SHOULD NOT RUN INTO HERE.'
            ENDIF
        ELSE
            PRINT *, 'SOMETHING WRONG IN GEN_ADMISSIBLE_SLIP_SURFACE.IT SHOULD HAVE TWO INTERSECTION POINTS.'            
        ENDIF
        IF(Y(I)<YL1) THEN
            Y(I)=YL1
        ENDIF
        IF(SLOPE_PARA.NWSP>0) THEN
            T1=SLOPE_PARA.GETY(X(I),2)
            Y(I)=MAX(T1,Y(I))            
        ENDIF
        CYCLE
    ENDIF
    
    !IMPROVED BY LGY
    K1=(Y(I-1)-Y(NX))/(X(I-1)-X(NX))
    YU1=MIN(K1*(X(I)-X(I-1))+Y(I-1),YU1) 
    IF(SLOPE_PARA.XCL>X(I)) THEN
        YCL1=INTERPOLATION(GX(1,:),GX(2,:),NGX,SLOPE_PARA.XCL)    
        K1=(Y(I-1)-YCL1)/(X(I-1)-SLOPE_PARA.XCL)
        YU1=MIN(K1*(X(I)-X(I-1))+Y(I-1),YU1)
    ENDIF
    
    DO J=1,NGX
        IF(GX(1,J)>X(I).AND.GX(1,J)<X(NX)) THEN
            K1=(Y(I-1)-GX(2,J))/(X(I-1)-GX(1,J))
            YU1=MIN(K1*(X(I)-X(I-1))+Y(I-1),YU1) 
        ENDIF
    ENDDO
    !IMPROVED BY LGY
    
    IF(I>2) THEN 
        K1=(Y(I-1)-Y(I-2))/(X(I-1)-X(I-2))
        YL1=MAX(K1*(X(I)-X(I-1))+Y(I-1),YL1)
    ENDIF
    !ASSUMMING THE EXIT ANGLE <45 DEGREE
    IF(I==2) THEN 
        IF((Y(I-1)-YL1)>ABS(X(I)-X(I-1))) THEN
            YL1=Y(I-1)-ABS(X(I)-X(I-1))
        ENDIF
    ENDIF
    
    IF(YU1-YL1<-0.001d0) THEN
        ISERROR=1
        PRINT *, "UNADMISSIBILITY OCCURS IN GEN_ADMISSIBLE_SLIP_SURFACE. YU<YL"
    ENDIF
    
    Y(I)=YL1+(YU1-YL1)*RDF(I)
    IF(PRESENT(YU)) YU(I)=YU1;
    IF(PRESENT(YL)) YL(I)=YL1
ENDDO

END SUBROUTINE    



SUBROUTINE GEN_CIRCULAR_SLIP(XT,XC,RMIN,RMAX,CENTER)
!GIVEN:
!X ORDINATE OF THE EXIT ANE ENTRY POINT: XT AND XC
!RETURN:
!RADIUS LIMITS OF ADMISSIBLE CIRCLAR SLIPS: RMIN,RMAX
!CENTER(2,2):CENTER(:,1)=CENTER FOR RMIN,...
!DEPENDS:
!SLOPE_PARA:Gx,Rx

REAL(8),INTENT(IN)::XT,XC
REAL(8),INTENT(OUT)::RMIN,RMAX,CENTER(2,2)
REAL(8)::YT,YC,TCL1,RA1(7),P1(2),P2(2),P3(2),R1,RA2(11),K1,T1,DX,DY
INTEGER::I,J,K,SD1,ISACW1


YT=SLOPE_PARA.GETY(XT)
YC=SLOPE_PARA.GETY(XC)

CENTER(:,1)=([XT,YT]+[XC,YC])/2.0
RMIN=NORM2([XT,YT]-[XC,YC])/2.0


RMAX=100*RMIN

IF(ABS(YC-YT)<1.D-8) THEN
    K1=1.D8
ELSE
    K1=-(XC-XT)/(YC-YT)
ENDIF
T1=(RMAX**2-RMIN**2)**0.5
DY=T1/(1+(1/K1)**2)**0.5
DX=DY/K1

!RA1(1:4)=CIRCLE_3P([XT,YT,XC,YC,CENTER(1,1),CENTER(2,1)-0.001])
!RMAX=RA1(4)
CENTER(:,2)=CENTER(:,1)+[DX,DY]

!CHECK ADMISSIBILITY

IF(ABS(XC-XT)>1E-7) THEN
    SD1=INT(SIGN(1.,(YC-YT)/(XC-XT))) !LEFT SLIDE SLOPE=1,RIGHT SLIDE SLOPE=-1
ELSE
    SD1=1
ENDIF


!Rmax
!CHECK WHETHER THERE IS ANY POINT BETWEEN XT,XC ON THE GE OUTSIDE THE CIRCLE
DO I=1,SLOPE_PARA.NGX
    IF((SLOPE_PARA.GX(1,I)-XT)*(SLOPE_PARA.GX(1,I)-XC)<0.0D0) THEN
        ISACW1=ISACW(XT,YT,0.D0,XC,YC,0.D0,SLOPE_PARA.GX(1,I),SLOPE_PARA.GX(2,I),0.D0)
        IF(ISACW1>1) CYCLE !COLINEAR
        IF(ISACW1*SD1>0) CYCLE 
        RA1(1:4)=CIRCLE_3P([XT,YT,XC,YC,SLOPE_PARA.GX(:,I)])
        IF(ABS(RA1(1))<1.E-7.AND.RA1(4)<RMAX) THEN !NO ERROR
             CENTER(:,2)=RA1(2:3);RMAX=RA1(4)
        ENDIF
    ENDIF
ENDDO
!Rmin
!DO I=0,SLOPE_PARA.NRX
!    IF(I==0) THEN
!        P1(1)=SLOPE_PARA.RX(1,1)
!        P1(2)=SLOPE_PARA.GETY(P1(1))
!    ELSE
!        P1=SLOPE_PARA.RX(:,I)
!    ENDIF
!    IF(I==SLOPE_PARA.NRX) THEN
!        P2(1)=SLOPE_PARA.RX(1,I)
!        P2(2)=SLOPE_PARA.GETY(P2(1))
!    ELSE
!        P2=SLOPE_PARA.RX(:,I+1)
!    ENDIF
!    
!    RA1=CIRCLE_LINE_INTERSECTION(CENTER(:,1),RMIN,P1,P2)
!    
!    
!    
!    IF(RA1(1)>1.1) THEN 
!        IF((RA1(4)<0.D0.AND.RA1(7)<0.D0).OR.(RA1(4)>1.D0.AND.RA1(7)>1.D0)) CYCLE
!        !P3=[(RA1(2)+RA1(5))/2.0,(RA1(3)+RA1(6))/2.0]
!        !RA1(1:4)=CIRCLE_3P([XT,YT,XC,YC,P3])
!        RA2(:11)=CIRCLE_2P_TANGENTLINE(RESHAPE([XT,YT,XC,YC,P1,P2],([2,4])))
!        
!        IF(ABS(RA2(1))<1E-7) THEN
!            R1=RA2(6);P3=RA2(2:3)
!        ELSEIF(RA2(1)>0.5) THEN
!            IF(RA2(6)<=RA2(11)) THEN
!                R1=RA2(6);P3=RA2(2:3)
!            ELSE
!                R1=RA2(11);P3=RA2(7:8)
!            ENDIF
!        ENDIF
!        
!                
!        IF(R1>RMIN) THEN
!            CENTER(:,1)=P3;RMIN=R1               
!        ENDIF
!        
!        !IF(CENTER(2,1)-RMIN<19.999) THEN
!        !    PRINT *, 'UNEXPECTED ERROR 3'
!        !    PAUSE
!        !ENDIF
!       
!    ENDIF 
!    
!ENDDO

!MAX(X OF CIRCULAR SLIP) MUST <XC FOR LEFT SLOPE, AND >XC FOR RIGHT SLOPE.
P1(1)=XC
P1(2)=YC+0.01
RA1(1:4)=CIRCLE_3P([XT,YT,XC,YC,P1])
IF(RA1(4)>RMIN) THEN
    CENTER(:,1)=RA1(2:3);RMIN=RA1(4)                
ENDIF

IF(RMAX<RMIN) THEN
    WRITE(*,100)
    WRITE(*,200) XT,YT,XC,YC,CENTER(:,1),RMIN,CENTER(:,2),RMAX
    PAUSE
ENDIF

100 FORMAT('CANNOT FIND A COMPATABLE CIRCULAR SLIP.')
200 FORMAT('[XT,YT,XC,YC,XMIN,YMIM,RMIN,XMAX,YMAX,RMAX]=',10(F10.4,1X))

ENDSUBROUTINE

FUNCTION CIRCLE_LINE_INTERSECTION(CENTER,R,P1,P2) RESULT(INTERSECT)
!GIVEN:
!CIRCLE:CENTER(2),R
!LINE:P1,P2
!RETURN
!INTERSECTION:INTERSECT(0:4),INTERSECT(0)=-1,ERROR,INTERSECT(0)=0,NO INTERSECT,=1,ONE INTERSECT;=2,TWO INTERSECTS.
REAL(8),INTENT(IN)::CENTER(2),R,P1(2),P2(2)
REAL(8)::INTERSECT(0:6) 
!INTERSECT(1:6)=[X1,Y1,U1,X2,Y2,U2] 
!IF 0<=U<=1, THE INTERSECT POINT (IP) INSIDE [P1,P2],P1-IP-P2
!IF U<0 IP-P1-P2
!IF U>1 P1-P2-IP
REAL(8)::A,B,C,T1,U

!a = (x2 - x1)2 + (y2 - y1)2 + (z2 - z1)2
!b = 2[ (x2 - x1) (x1 - x3) + (y2 - y1) (y1 - y3) + (z2 - z1) (z1 - z3) ]
!c = x32 + y32 + z32 + x12 + y12 + z12 - 2[x3 x1 + y3 y1 + z3 z1] - r2

A=SUM((P2-P1)**2)
B=2*DOT_PRODUCT((P2-P1),(P1-CENTER))
C=SUM(CENTER**2+P1**2)-2*DOT_PRODUCT(CENTER,P1)-R**2

INTERSECT=-1.D0

IF(ABS(A)<1.D-10) THEN
    IF(IS_IN_CIRCLE(CENTER,R,P1)==0) THEN
        INTERSECT(0)=1.0
        INTERSECT(1:2)=P1
        INTERSECT(3)=0.D0
    ENDIF
    RETURN
ENDIF

T1=B**2-4*A*C

IF(T1>=0.D0) THEN
    T1=SQRT(T1)
    IF(T1<1.E-7) THEN
        U=-B/(2*A)
        INTERSECT(0)=1.0
        INTERSECT(1:2)=P1+U*(P2-P1)
        INTERSECT(3)=U
    ELSE
        INTERSECT(0)=2.0
        U=(-B+T1)/(2*A)
        INTERSECT(1:2)=P1+U*(P2-P1)
        INTERSECT(3)=U
        U=(-B-T1)/(2*A)
        INTERSECT(4:5)=P1+U*(P2-P1)
        INTERSECT(6)=U
    ENDIF
ENDIF



END FUNCTION

INTEGER FUNCTION IS_IN_CIRCLE(CENTER,R,P)
!GIVEN:
!CIRCLE:CENTER(2),R,
!POINT P(2)
!RETURN:
!1,INSIDE,-1=OUSIDE,0=ON THE CIRCLE LINE

REAL(8),INTENT(IN)::CENTER(2),R,P(2)
REAL(8)::T1

T1=R-NORM2(CENTER-P)

IF(ABS(T1)<1.E-7) THEN
    IS_IN_CIRCLE=0
ELSE
    IS_IN_CIRCLE=NINT(SIGN(1.D0,T1))
ENDIF


END FUNCTION

REAL(8) FUNCTION GET_BC_Y(SELF,X,FLAG)
!GIVEN:
!X
!FLAG, =0, GETY FROM GE,=1 GET Y FROM RE,=OTHERS,GET Y FROM WSP
!RETURN:
!Y

CLASS(SLOPE_PSO_PARAR)::SELF
REAL(8),INTENT(IN)::X
INTEGER,OPTIONAL,INTENT(IN)::FLAG
REAL(8),DIMENSION(:,:),ALLOCATABLE:: PLINE1
INTEGER::IFLAG1

IF(PRESENT(FLAG)) THEN
    IFLAG1=FLAG    
ELSE
    IFLAG1=0
ENDIF

IF(IFLAG1==0) THEN
    PLINE1=SELF.GX
ELSEIF(IFLAG1==1) THEN
    PLINE1=SELF.RX
ELSE
    PLINE1=SELF.WSP
ENDIF

GET_BC_Y=interpolation(PLINE1(1,:),PLINE1(2,:),SIZE(PLINE1,DIM=2),X)

IF(ALLOCATED(PLINE1)) DEALLOCATE(PLINE1)

END FUNCTION


FUNCTION CIRCLE_3P(PT) RESULT(CR)
!REF:https://www.qc.edu.hk/math/Advanced%20Level/circle%20given%203%20points.htm
    !GIVEN:
    !THREE POINTS(NOT COLINEAR),PT(2,3)
    !RETURN:
    !CENTER AND R IN CR(3)
    REAL(8),INTENT(IN)::PT(2,3)
    REAL(8)::CR(0:3) ![ERROR,CX,CY,R] ERROR=-1,COLINEAR
    REAL(8)::K1,K2

    CR=0.D0
    !TEST WHETHER COLINEAR
    IF(ISACW(PT(1,1),PT(2,1),0.D0,PT(1,2),PT(2,2),0.D0,PT(1,3),PT(2,3),0.D0)>1) THEN
        PRINT *, 'THREE POINTS ARE COLINEAR.'
        CR(0)=-1.    
        RETURN
    ENDIF

    IF(ABS(PT(2,2)-PT(2,1))>1E-8) THEN
        K1=-((PT(1,2)-PT(1,1)))/(PT(2,2)-PT(2,1))
    ELSE
        K1=1.D8
    ENDIF

    IF(ABS(PT(2,3)-PT(2,1))>1E-8) THEN
        K2=-((PT(1,3)-PT(1,1)))/(PT(2,3)-PT(2,1))
    ELSE
        K2=1.D8
    ENDIF
    CR(1)=(-K1*PT(1,2)+K2*PT(1,3)+PT(2,2)-PT(2,3))/(K2-K1)
    CR(2)=K1*(CR(1)-PT(1,2))+PT(2,2)
    CR(3)=NORM2(PT(:,1)-CR(1:2))/2.0
    CR(1:2)=(PT(:,1)+CR(1:2))/2.0

END FUNCTION

FUNCTION CIRCLE_2P_TANGENTLINE(PT) RESULT(CR)
!REF:https://www.cut-the-knot.org/Curriculum/Geometry/GeoGebra/PPL.shtml#solution
    !GIVEN:
    !TWO POINTS,PT(2,1:2) AND LINE DEFINED BY PT(2,3:4)
    !RETURN:
    !CENTER, TANGENT POINT AND  R IN CR(0:10), CR(0)=ERROR CODE, =0, ONE CIRCLE,=1,TWO CIRCE.
    REAL(8),INTENT(IN)::PT(2,4)
    REAL(8)::CR(0:10) ![ERROR,CX1,CY1,PT1X,PT1Y,R1,CX2,CY2,PT2X,PT2Y,R2] ERROR=-1,COLINEAR
    REAL(8)::K1,K2,PT5(2),DX1,DY1,PT6(2),PT7(2),XC1(2),T1

    
    CR(0)=0.D0
    !INTERSECTION OF P1P2 AND P3P4

    IF(ABS(PT(1,2)-PT(1,1))>1E-8) THEN
        K1=(PT(2,2)-PT(2,1))/(PT(1,2)-PT(1,1))
    ELSE
        K1=1.D8
    ENDIF

    IF(ABS(PT(1,3)-PT(1,4))>1E-8) THEN
        K2=((PT(2,3)-PT(2,4)))/(PT(1,3)-PT(1,4))
    ELSE
        K2=1.D8
    ENDIF
    XC1=(PT(:,1)+PT(:,2))/2.0

    IF(ABS(K1-K2)<1E-7)THEN
        IF(ABS(K1)<1E-7) K1=1E-7
        K1=-1/K1
        PT6(1)=(-K1*XC1(1)+K2*PT(1,3)+XC1(2)-PT(2,3))/(K2-K1)
        PT6(2)=K1*(PT6(1)-XC1(1))+XC1(2)
        CR(:3)=CIRCLE_3P([PT(:,1:2),PT6])
        CR(3:4)=PT6
        CR(5)=CR(3)
        RETURN
    ENDIF


    PT5(1)=(-K1*PT(1,2)+K2*PT(1,3)+PT(2,2)-PT(2,3))/(K2-K1)
    PT5(2)=K1*(PT5(1)-PT(1,2))+PT(2,2)

    !Power of a Point Theorem
    T1=NORM2(PT(:,1)-PT5)*NORM2(PT(:,2)-PT5)
    DX1=(T1/(1+K2**2))**0.5
    DY1=DX1*K2
    PT6=PT5+[DX1,DY1];PT7=PT5-[DX1,DY1]

    IF(ABS(K1)<1E-8) K1=1E-8
    IF(ABS(K2)<1E-8) K2=1E-8

    K1=-1/K1;K2=-1/K2

    CR(0)=1.D0
    CR(1)=(-K1*XC1(1)+K2*PT6(1)+XC1(2)-PT6(2))/(K2-K1)
    CR(2)=K1*(CR(1)-XC1(1))+XC1(2)
    CR(3:4)=PT6
    CR(5)=NORM2(CR(1:2)-PT6)

    CR(6)=(-K1*XC1(1)+K2*PT7(1)+XC1(2)-PT7(2))/(K2-K1)
    CR(7)=K1*(CR(4)-XC1(1))+XC1(2)
    CR(8:9)=PT7
    CR(10)=NORM2(CR(6:7)-PT7)

ENDFUNCTION

SUBROUTINE POLYLINE_FOS_CAL(X,NX,DL,FM,FA)
!GIVEN:
!POLYLINE VETEX:X(2,NX)
!SUBDIVISION LENGTH: DL
!RETURN:
!MOBOLIZED SHEAR FORCE ALONG THE POLYLINE: FM
!AVAILABLE SHEAR FORCE ALONG THEPOLYLINE:FA
!THE FOS CAN BE CALCULATED AS: FOS=FA/FM
    IMPLICIT NONE
    INTEGER,INTENT(IN)::NX
    REAL(8),INTENT(IN)::X(2,NX),DL
    REAL(8),INTENT(OUT)::FM,FA

    INTEGER::I
    REAL(8)::FM1,FA1
    !OPEN(10,FILE='POLYLINE_FOS_CAL_DEBUG.DAT',STATUS='REPLACE')
    !WRITE(10,100)
    FM=0.D0;FA=0.D0
    DO I=1,NX-1
        CALL SEGMENT_SHEAR_FORCE_CAL(X(:,I),X(:,I+1),DL,FM1,FA1)
        FM=FM+FM1;FA=FA+FA1
        !WRITE(10,110) I,FA1,FM1,X(1:2,I),X(1:2,I+1)
    ENDDO
!    CLOSE(10)
!100 FORMAT(4X,'ISEGM',1X,'TFORCE_A',1X,'TFORCE_M',7X,'X1',7X,'Y1',7X,'X2',7X,'Y2')
!110 FORMAT(I8,6(F8.3,1X))
    
END SUBROUTINE

SUBROUTINE SEGMENT_SHEAR_FORCE_CAL(XA,XB,DL,FM,FA)
!GIVEN A LINE SEGMENT DEFINED BY XA,XB,AND ITS DIVISION LENGTH DL
!CALCULATE TEH SUM OF THE MOBLIZED SHEAR FORCE FM AND IT SHEAR FORCE AVALABLE FA.

	IMPLICIT NONE
	REAL(8),INTENT(IN)::XA(2),XB(2),DL
	REAL(8),INTENT(OUT)::FM,FA
	REAL(8)::DX1(2),XI1(2),T1,TAU1,TAUF1,RAD1
	REAL(8),PARAMETER::PI1=3.141592653589793
	INTEGER::N1,K
    
	FM=0.D0;FA=0.D0
    
	IF(ABS(DL)<1.D-7) THEN
		STOP "ERROR IN  SEGMENT_SHEAR_FORCE_CAL. DIVISION LENGTH =0." 		
	ENDIF
	N1=NINT(NORM2(XA-XB)/DL)
	N1=MAX(N1,1)
	DX1=(XB-XA)/N1
    
	
	
	IF(ABS(DX1(1))>1.D-7) THEN
		RAD1=ATAN(DX1(2)/DX1(1))
	ELSE
		RAD1=SIGN(PI1/2.0,DX1(2))
	ENDIF
	T1=NORM2(DX1)
    !IF(N1>2) THEN
    !    PRINT *, N1
    !ENDIF
	DO K=1,N1
		XI1=XA+(real(K)-0.5)*DX1		
		CALL POINT_SFR_CAL(XI1,RAD1,TAU1,TAUF1)
		FM=FM+TAU1*T1;FA=FA+TAUF1*T1
	ENDDO
    
    !CALL POINT_SFR_CAL((XB+XA)/2,RAD1,TAU1,TAUF1)
    !TAU1=TAU1*NORM2(XB-XA); TAUF1=TAUF1*NORM2(XB-XA);
    !PRINT  *, FM-TAU1,FA-TAUF1
    
    RETURN
    
END SUBROUTINE

SUBROUTINE POINT_SFR_CAL(X,RAD,TAU,TAUF)
!GIVEN POSITION=X(3) NAD ANGLE=RAD,
!CALCULATE TAU AND TAUF
	USE POS_IO,ONLY : POSDATA 
	IMPLICIT NONE	
	REAL(8),INTENT(IN)::X(2),RAD
	REAL(8),INTENT(OUT)::TAU,TAUF
	INTEGER,SAVE::IEL1=0
	REAL(8)::SS(3),C1,PHI1,SFR1,SNT1(2),VAL1(100)
	REAL(8),PARAMETER::PI1=3.141592653589793
	
	!INTEGER,EXTERNAL::POINTlOC_BC
    
	VAL1=0.D0
    !IEL1=POINTlOC_BC(X,IEL1)
	!TRYIEL1=IEL1
	!IF(IEL1>0) THEN 
	!	call getval(X,iel1,VAL1(1:POSDATA.NVAR))
	!ELSE
	!	STOP "ERROR IN POINT_SFR_CAL. POINT CANNOT BE LOCATED. "
	!ENDIF
    
	CALL ProbeatPhyscialspace([X,0.D0],VAL1(1:POSDATA.NVAR),iel1)
    IF(IEL1==0) THEN
        TAU=0.D10;TAUF=1.D10
        PRINT *, 'POINT LOCATION FAILED.X=',X 
        RETURN
    ENDIF
			
	SS(1:3)=[VAL1(POSDATA.ISXX),VAL1(POSDATA.ISYY),VAL1(POSDATA.ISXY)]
	C1=VAL1(POSDATA.IMC_C)
	PHI1=VAL1(POSDATA.IMC_PHI)
	SFR1=VAL1(POSDATA.ISFR)
	CALL stress_in_inclined_plane(SS,RAD,SNT1)
    TAU=SNT1(2)
	IF(SNT1(1)<0.D0) THEN
		TAUF=(C1-SNT1(1)*dTAN(PHI1/180.0*PI1)) !TAUF		
	ELSE
		TAUF=0.D0;
	ENDIF	

END SUBROUTINE

integer function isacw(x1,y1,z1,x2,y2,z2,x3,y3,z3)
	implicit none
	real(8),intent(in)::x1,y1,z1,x2,y2,z2,x3,y3,z3
    real(8)::t1,yn2,xn2,zn2,yn3,xn3,zn3,norm1(3),t2
	
	isacw=0
    !isacw=1,=2,onedge,=3,coline;
    
	yn2=y2-y1
	xn2=x2-x1
    zn2=z2-z1
	yn3=y3-y1
	xn3=x3-x1    
    zn3=z3-z1
    norm1=[yn2*zn3-zn2*yn3,-(xn2*zn3-zn2*xn3),xn2*yn3-yn2*xn3]
    t2=norm2(norm1)
    if(t2<1e-10) then !共线
        ISACW=3
        if(x1<min(x3,x2)) return
        if(X1>max(x3,x2)) return
        if(Y1<min(y3,y2)) return
        if(Y1>max(y3,y2)) return 
        if(Z1<min(z3,z2)) return
        if(Z1>max(z3,z2)) return 
        isacw=2 !on the edge
    else
        !t1=(xn2*yn3-yn2*xn3)+(yn2*zn3-zn2*yn3)-(xn2*zn3-zn2*xn3)
        !if(abs(t1)<1.d-10) then 
            !与(1,1,1)垂直            
        if(abs(norm1(3))>1e-10) then
            isacw=sign(1.,norm1(3)) !从+z看
        elseif(abs(norm1(1))>1e-10) then
            isacw=sign(1.,norm1(1)) !从+x看
        elseif(abs(norm1(2))>1e-10) then
            isacw=sign(1.,norm1(2)) !从+y看
        endif
        !else
            !从(1,1,1)方向看
        !    isacw=sign(1.,t1)
        !endif
    endif
	

end function

function interpolation(x,y,nx,xi)
!x,y must be in order.
	implicit none
	INTEGER,PARAMETER::Double=KIND(1.0D0)
    integer,intent(in)::nx
	real(double),intent(in)::x(nx),y(nx),xi
	real(double)::interpolation,t1
	integer::i
    
    interpolation=0.D0
    
    if(nx==1.AND.ABS(XI-X(1))<1.D-6) then
       interpolation=y(1)
       return
    endif
    do i=1,nx-1
        if((xi<=x(i+1).and.xi>=x(i)).or.(xi<=x(i).and.xi>=x(i+1))) then
	        t1=x(i+1)-x(i)
	        if(abs(t1)<1e-7) then
		        print *, "Warning! 分母=0,function=Interpolation()"
		        interpolation=(y(i)+y(i+1))/2.0d0
	        else
		        interpolation=(y(i+1)-y(i))/(t1)*(xi-x(i))+y(i)
            endif
            return
        endif
    enddo
    if(i==nx) then
        interpolation=-HUGE(1.0D0)
        print *, "xi is out of the range.function=Interpolation()"
    endif
    
    
endfunction

END MODULE