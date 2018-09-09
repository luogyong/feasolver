MODULE SLOPE_PSO
    USE stochoptim,ONLY: Evolutionary,EA_ITER
    USE solverds,ONLY:INPUTFILE
    !USE IndexColor
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
        INTEGER::NGX,NRX,NSLICE=30,IS_SLOPEINIT=0
        REAL(8)::XTU,XTL,XCU,XCL !ENTRY AND EXIT LIMITS
        REAL(8),ALLOCATABLE::GX(:,:),RX(:,:)
        !EA PARAMETER
        integer(kind = 4) :: popsize = 30, max_iter = 1000,ndim=31
        REAL(8)::EPS1=1.D-4,EPS2=1.D-4,EPS3=1.D-3,w = 0.7298, c1 = 1.49618,c2 = 1.49618, gamma = 1., &
            F=0.5,CR=0.1,DL=0.5D0,mu_perc=0.5,sigma=0.5  !DL=SLICE BASE STRESS SAMPLE DISTANCE
        character(len = 5) :: solver='de', strategy='rand1'
        
        !!FOR STREAMLINE ADMISSIBLITY REPAIR ONLY
        INTEGER::DIRECTION=-1 !=-1 TOE,=1,(DOWNHILL),CREST(UPHILL).
        REAL(8)::A=0,B=0,C=0 !REPAIRED RANGE LINE ,A*X+B*Y+C
        CONTAINS
        PROCEDURE::INIT=>FILEREADIN
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
    allocate(X(SLOPE_PARA.NDIM,SLOPE_PARA.POPSIZE),Y(SLOPE_PARA.NDIM,SLOPE_PARA.POPSIZE),RDF(SLOPE_PARA.NDIM,SLOPE_PARA.POPSIZE))
    lower =0.D0
    upper =1.D0
    
    CALL PSO_SLIP.INIT(SLOPE_PARA.NDIM,SLOPE_PARA.POPSIZE) !SAVE TOPTENSLIP
    
   CALL GEN_SLIP_SAMPLE(SLOPE_PARA.XTL,SLOPE_PARA.XTU,SLOPE_PARA.XCL,SLOPE_PARA.XCU, &
        SLOPE_PARA.NSLICE,SLOPE_PARA.popsize,SLOPE_PARA.GX,SLOPE_PARA.NGX,  & 
        SLOPE_PARA.RX,SLOPE_PARA.NRX,X,Y,RDF)
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
    REAL(8)::XT1(2),XC1(2),X1(2,SIZE(UNS)),FM1,FA1,T1
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
        CALL POLYLINE_FOS_CAL(X1,SLOPE_PARA.NDIM,SLOPE_PARA.DL,FM1,FA1)
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
    REAL(8),INTENT(OUT)::X1(2,SIZE(UNS))
    INTEGER,INTENT(OUT)::ISERROR
    REAL(8)::XT1(2),XC1(2),UNS1(SIZE(UNS))
    INTEGER::IFLAG=1
    
    UNS1=UNS
    CALL GET_ENTRY_EXIT([UNS(1),UNS(SLOPE_PARA.NDIM)],XT1,XC1)
    IFLAG=1;ISERROR=0
    CALL GEN_ADMISSIBLE_SLIP_SURFACE(XT1,XC1,SLOPE_PARA.NSLICE,SLOPE_PARA.GX,SLOPE_PARA.NGX,&
        SLOPE_PARA.RX,SLOPE_PARA.NRX,UNS1,IFLAG,X1(1,:),X1(2,:),ISERROR)
    
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
        XC1(2)=IFLAG_OPTIM_PARA(2)
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
                    case('ndim')
                        INDATA.ndim=int(COMMAND.OPTION(i).VALUE)
                        INDATA.nslice=INDATA.ndim-1
                    case('dl')
                        INDATA.dl=COMMAND.OPTION(i).VALUE
                    case('strategy')
                        INDATA.strategy=trim(adjustl(COMMAND.OPTION(i).CVALUE))
                    CASE('sigma')
                        INDATA.sigma=COMMAND.OPTION(i).VALUE
                    CASE('mu_perc')
                        INDATA.mu_perc=COMMAND.OPTION(i).VALUE
				    case default
					    call Err_msg(COMMAND.OPTION(i).name)
			    end select
		    end do
            IF(N1==0) INDATA.popsize=INDATA.ndim*2
            
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



SUBROUTINE GEN_SLIP_SAMPLE(XTL,XTU,XCL,XCU,NSEG,NSLIP,GX,NGX,RX,NRX,X,Y,RDF)

!GIVEN:
!THE LIMITS OF THE SLOPE TOE: XTL,XTU
!THE LIMITS OF THE SLOPE CREST: XCL,XCU
!SAMPLE NUMBER: NSLIP
!GROUND SURFACE: GX(2,NGX)
!BEDROCK SURFACE: RX(2,NRX)
!NUMBER OF SUBDIVISION OF EACH SLIP LINE: NSEG>1
!RETURN:
!THE SLIP NODES: X(NSEG+1,NSLIP),Y(NSEG+1,NSLIP)
!(RADOM) FACTOR: RDF(NSEG+1)
    
    IMPLICIT NONE

    INTEGER,INTENT(IN)::NSEG,NSLIP,NGX,NRX
    REAL(8),INTENT(IN)::XTL,XTU,XCL,XCU,GX(2,NGX),RX(2,NRX)
    REAL(8),INTENT(OUT)::X(NSEG+1,NSLIP),Y(NSEG+1,NSLIP),RDF(NSEG+1,NSLIP)

    INTEGER::I,J,ISERROR=1,NITER=0,IFLAG=0,NITER2=0
    REAL(8)::R1,R2,XT(2),XC(2),YL(NSEG+1),YU(NSEG+1)

    
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
            CALL GEN_ADMISSIBLE_SLIP_SURFACE(XT,XC,NSEG,GX,NGX,RX,NRX,RDF(:,J),IFLAG,X(:,J),Y(:,J),ISERROR,YL,YU)
            RDF(1,J)=R1;RDF(NSEG+1,J)=R2
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
    
SUBROUTINE GEN_ADMISSIBLE_SLIP_SURFACE(XT,XC,NSEG,GX,NGX,RX,NRX,RDF,IFLAG,X,Y,ISERROR,YL,YU)
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
!RETURN: 
!THE X COMPONENT OF THE VETEX OF THE SLIP LINE:X(NSEG+1)
!THE Y COMPONENT : Y(NSEG+1)
!ISERROR=0,SUCCESSFULLY GENERATE A SLIP SURFACE.
!/=0, NO ADMISSIBLE SLIP SURFACE FOUND.
IMPLICIT NONE
INTEGER,INTENT(IN)::NGX,NRX,NSEG,IFLAG
REAL(8),INTENT(IN)::XT(2),XC(2),GX(2,NGX),RX(2,NRX)
REAL(8),INTENT(INOUT)::RDF(NSEG+1)
REAL(8),INTENT(OUT)::X(NSEG+1),Y(NSEG+1)
REAL(8),OPTIONAL,INTENT(OUT)::YU(NSEG+1),YL(NSEG+1) 
INTEGER,INTENT(OUT)::ISERROR

INTEGER::I,J,NX
REAL(8)::DX1,YU1,YL1,K1,YCL1

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

DO I=2,NSEG    
    YL1=INTERPOLATION(RX(1,:),RX(2,:),NRX,X(I))    
    YU1=INTERPOLATION(GX(1,:),GX(2,:),NGX,X(I))
    
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
        stop "xi is out of the range.function=Interpolation()"
    endif
    
    
endfunction



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
        PRINT *, 'POINT LOCATION FAILED.X=',X    
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


END MODULE