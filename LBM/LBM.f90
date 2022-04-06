
    
module lbm_model
    use INPUT_PARSER
    implicit none
    private
    public::lbm_readin
    
    integer,parameter::d2q9=1,d3q19=2,d1q3=0
    integer,parameter::dif=0,dif_con=1,iiflow=2
    integer::ltype=1; !lattice type
    integer::ptype=2; !problem type
    integer::dim=2,q=9,nx,ny,nz,mstep
    real(8)::dt=1.0d0,dx=1.0d0,xmin(3)=0.d0
    real(8)::rh=1.0d0 ! density
    real(8)::cs=1/3.**0.5 !
    
    real(8),allocatable::w(:)
    type lbm_bc_tydef
        integer::type=-1 
        !-1,free；0,fixed(Dielectric); 1,flow(numen);2,periodic 
        integer::dir=0 !direction
        !1,2,3=x,y,z;-1,-2,-3=-x,-y,-z
        real(8)::value=0
    endtype
        
    type lattice_tydef
        integer::mat=0
        real(8),allocatable::f(:),feq(:),s(:)
        !macro variable
        real(8)::rh
        real(8),allocatable::v(:)
        type(lbm_bc_tydef)::bc
    contains
        
    endtype
    type(lattice_tydef),allocatable::latdata(:,:,:)
    
    type lbm_mat_tydef
        integer::type        
        real(8)::property(10)        
    end type
    type(lbm_mat_tydef),allocatable::material(:)
    integer::nmat=0    
    
    type boxzone_tydef
        integer::box(2,3) !xmin,xmax,...
        integer::mat=1
    end type
    type(boxzone_tydef),allocatable::box(:)
    integer::nbox=0
    
    contains
    
    subroutine lbm_initialcondition()
        implicit none
        
    endsubroutine
    
    
    
    subroutine lbm_space_allocate()
        implicit none
        integer::ERR
        
        allocate(latdata(0:nx,0:ny,0:nz),STAT=ERR)
        allocate(latdata.f(Q),latdata.feq(Q),STAT=ERR)
        if(ptype>0) then
            allocate(latdata.v(dim),STAT=ERR)
        endif
        allocate(w(0:q),STAT=ERR)
        
        if(ltype==d1q3) then
            w(0)=2/3.
            w(1:2)=1./6.
        elseif(ltype==d2q9) then
            w(0)=4/9.
            w(1:4)=1./9.
            w(5:8)=1./36.
        elseif(ltype==d3q19) then
            w(0)=1./3.
            w(1:6)=1./18.
            w(7:18)=1./36.
        else
            print *,'no such element type=',ltype
            error stop
        endif
        
        
    endsubroutine
    
    subroutine lbm_readin(file)
        implicit none
        
        character(len=*),intent(in)::file
        integer::ifile=10
        !external::lbmreadin_parser
        
        open(ifile,file=file,status='old')
        
        call read_execute(ifile,lbm_readin_parser,lbm_string2value)
        
        
        close(ifile)
        
    endsubroutine
    
    subroutine lbm_readin_parser(COMMAND,unit)
        implicit none
        TYPE(COMMAND_TYDEF),INTENT(IN)::COMMAND
        integer,intent(in)::unit
        integer::i,ival1,n1,n2,n3,ix1,iy1,iz1
        real(8)::val1
        select case(trim(COMMAND.KEYWORD))
        case('lbm_paras')
            print *, 'reading lbm_paras data...'
            do i=1, COMMAND.NOPT
                val1=COMMAND.OPTION(i).value
                ival1=int(val1)
				select case(COMMAND.OPTION(i).NAME)
					case('ltype')
						ltype=ival1
                    case('nx')
                        nx=ival1
                    case('ny')
                        ny=ival1
                    case('nz')
                        nz=ival1
                    case('dt')
                        dt=val1
                    case('dx')
                        dx=val1
                    case('ptype')
                        ptype=ival1
                    case('rh')
                        rh=val1
					case default
						call Err_msg(COMMAND.OPTION(i).name)
				end select
            end do
            
            call lbm_space_allocate()
            
        case('lbm_bc')
            print *, 'reading lbm_bc data...'
            do i=1, COMMAND.NOPT
                val1=COMMAND.OPTION(i).value
                ival1=int(val1)
				select case(COMMAND.OPTION(i).NAME)
					case('num')
						n1=ival1
					case default
						call Err_msg(COMMAND.OPTION(i).name)
                end select                    
			end do            
            do i=1,n1
                read(unti,*) ix1,iy1,iz1,latdata(ix1,iy1,iz1).bc.type,latdata(ix1,iy1,iz1).bc.dir,latdata(ix1,iy1,iz1).bc.value
            enddo            
        case('lbm_initialcons')
            print *, 'reading lbm_initialcons data...'
            do i=1, COMMAND.NOPT
                val1=COMMAND.OPTION(i).value
                ival1=int(val1)
				select case(COMMAND.OPTION(i).NAME)
					case('rh')
						latdata.rh=val1
					case default
						call Err_msg(COMMAND.OPTION(i).name)
                end select                    
			end do            
        case('lbm_material')
            print *, 'reading lbm_material data...'
            do i=1, COMMAND.NOPT
                val1=COMMAND.OPTION(i).value
                ival1=int(val1)
				select case(COMMAND.OPTION(i).NAME)
				case('num')
					n1=val1
				case default
					call Err_msg(COMMAND.OPTION(i).name)
                end select                  
            end do 
            allocate(material(n1))
            do i=1,n1
                read(unit,*) n2,aterial(n2).type,n3,material(n2).property(1:n3) !matid,type,npro,property(npro)
            enddo
        case('lbm_zone')
            print *, 'reading lbm_zone data...'
            do i=1, COMMAND.NOPT
                val1=COMMAND.OPTION(i).value
                ival1=int(val1)
				select case(COMMAND.OPTION(i).NAME)
				case('num')
					nbox=val1
				case default
					call Err_msg(COMMAND.OPTION(i).name)
                end select                  
            end do 
            allocate(box(nbox))
            do i=1,nbox
                read(unit,*) box(i).box,box(i).mat
            enddo
            
        case default
            call Err_msg(COMMAND.KEYWORD)
        end select
        
        
    endsubroutine    
    
    subroutine lbm_string2value(str,value,cvalue)

	
	implicit none
	character(*)::str
	character(*)::cvalue
	real(8)::value
	integer(4)::msg
	
	!Justify the content of the STR whehter a string or a number
	!if the str(1:1) is a number, then, STR is regarded as a number. That means 
	!the of character of the character constant should not be a number or a '.,+,-'.
	str=adjustL(str)
	
	if(index('0123456789.+-',str(1:1))>0) then
		read(str,*) value
		return
	end if
 
	call lowcase(str)
    cvalue=str
	select case(trim(str))
    case('d1q3')
        value=d1q3
    case('d2q9')
        value=d2q9
    case('d3q19')
        value=d3q19
    case('dif')
        value=dif
    case('dif_con')
        value=dif_con
    case('iiflow')
        value=iiflow
       
    case default
    
	    print *, 'No such Constant:'//trim(str)//',It will be returned as a character variable.'
    
    end select
    
    end subroutine
    


    
end module
    
    
!    
!parameter (m=100) !m is the number of lattice nodes
!real fo(0:m),f1(0:m),f2(0:m),rho(0:m),feq(0:m),x(0:m)
!integer i
!open(2,ﬁle=’result’)
!dt=1.0
!dx=1.0
!x(0)=0.0
!do i=1,m
!x(i)=x(i-1)+dx
!end do
!csq=dx*dx/(dt*dt)
!alpha=0.25
!omega=1.0/(alpha/(dt*csq)+0.5)
!mstep=200 ! The total number of time steps
!twall=1.0 ! Left hand wall temperature
!! Initial condition
!do i=0,m
!rho(i)=0.0 ! Initial value of the domain temperature
!f1(i)=0.5*rho(i)
!f2(i)=0.5*rho(i)
!end do
!do kk=1,mstep
!! main loop
!! collision process:
!do i=0,m
!rho(i)=f1(i)+f2(i)
!feq(i)=0.5*rho(i)
!! since k1=k2=0.5, then feq1=feq2=feq
!f1(i)=(1.-omega)*f1(i)+omega*feq(i)
!f2(i)=(1.-omega)*f2(i)+omega*feq(i)
!end do
!!Streaming process:
!do i=1,m-1
!f1(m-i)=f1(m-i-1) ! f1 streaming
!f2(i-1)=f2(i) ! f2 streaming
!end do
!! Boundary condition
!f1(0)=twall-f2(0)! constant temperature boundary condition, x=0
!f1(m)=f1(m-1) ! adiabatic boundary condition, x=L
!f2(m)=f2(m-1)! adiabatic boundary condition, x=L
!end do
!! end of the main loop
!do i=0,m
!write(2,*)x(i), rho(i)
!end do
!stop
!end
