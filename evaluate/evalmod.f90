module evaluate

! The user can assign values to parameters that can be used in expressions with 
! the subroutine defparam. The calling syntax is
!
!     call defparam(symbol,value) or call defparam(symbol,expr)
!
! where symbol is the desired parameter name; value is a real, integer, or 
! complex variable (single or double precision); and expr is a string
! containing an expression to be evaluated. The value obtained by evaluating the
! expression expr is associated with the parameter symbol. Parameter names must
! begin with a letter (a-z, A-Z) and must not be longer than 24 characters.
! Parameter names are not case dependent.
!
! An expression can be evaluated with the subroutine evalexpr. The calling
! syntax is 
!
!          call evalexpr(expr,value)
!
! where expr is a string containing the expression to be evaluated; value is the
! result (single or double precision real, complex or integer). The 
! expression can contain the arithmetic operations +, -, *, /, or ^ as well as
! the functions sin, cos, tan, log, ln, abs, exp, sqrt, real, imag, conjg, and 
! ang (the function ang calculates the phase angle of its complex argument). The
! expression can also contain numerical values and previously defined parameters
! Grouping by nested levels of parentheses is also allowed. The parameters pi
! and i (imaginary unit) are predefined. Complex numbers can be entered as a+i*b 
! if the parameter i has not been redefined by the user. Complex numbers can also
! be entered using complex(a,b).
! Example expression: 
!
!          conjg(((cos(x) + sqrt(a+i*b))^2+complex(ln(1.6e-4),20))/2)
!
! An equation of the form <symbol> = <expression> can be evaluated using the 
! subroutine evaleqn. The calling syntax is
!
!          call evaleqn(eqn)       
!
! where eqn is a string containing the equation. The right-hand-side of the
! equation is evaluated and assigned to the symbol given by the left-hand-side.
!
! The value assigned to a symbol can be retrieved using the subroutine getparam.
! The calling syntax is
!
!          call getparam(sym,value)
!
! where sym is a symbol string; value is a numeric variable (any of the six 
! standard types).
!
! The symbols and their values in the symbol table can be listed using the
! subroutine listvar. The variable ierr is always available following a call
! to any of the above subroutines and is zero if there were no errors. The
! possible nonzero values for ierr are
!
! 1       Expression empty
! 2       Parentheses don't match
! 3       Number string does not correspond to a valid number
! 4       Undefined symbol
! 5       Less than two operands for binary operation
! 6       No operand for unary plus or minus operators
! 7       No argument(s) for function
! 8       Zero or negative real argument for logarithm
! 9       Negative real argument for square root
! 10      Division by zero
! 11      Improper symbol format
! 12      Missing operator
! 13      Undefined function
! 14      Argument of tangent function a multiple of pi/2
!
use precision
use strings

save
private
public :: valuep,evalexpr,defparam,evaleqn,getparam,listvar,ierr

type item                         
  character(len=24):: char
  character :: type
end type item

type param                        
  character (len=24):: symbol
  complex(kc8):: value
end type param

interface defparam                
  module procedure strdef       ! value given by expression      
  module procedure valdef_dc    ! Double precision complex value
  module procedure valdef_sc    ! Single precision complex value
  module procedure valdef_dr    ! Double precision real value
  module procedure valdef_sr    ! Single precision real value
  module procedure valdef_di    ! Double precision integer value
  module procedure valdef_si    ! Single precision integer value
end interface

interface evalexpr
  module procedure evalexpr_dc  ! Double precision complex result
  module procedure evalexpr_sc  ! Single precision complex result
  module procedure evalexpr_dr  ! Double precision real result
  module procedure evalexpr_sr  ! Single precision real result
  module procedure evalexpr_di  ! Double precision integer result
  module procedure evalexpr_si  ! Single precision integer result
end interface

interface getparam
  module procedure getparam_dc  ! Double precision complex result
  module procedure getparam_sc  ! Single precision complex result
  module procedure getparam_dr  ! Double precision real result
  module procedure getparam_sr  ! Single precision real result
  module procedure getparam_di  ! Double precision integer result
  module procedure getparam_si  ! Single precision integer result
end interface

integer,parameter :: numtok=100  ! Maximum number of tokens 
type(param) :: params(100)       ! Symbol table
integer :: nparams=0,itop,ibin
complex(kc8) :: valstack(numtok) ! Stack used in evaluation of expression
type(item):: opstack(numtok)     ! Operator stack used in conversion to postfix
integer :: ierr                  ! Error flag

!**********************************************************************

contains

!**********************************************************************


SUBROUTINE EVALEXPR_DC(expr,val)    ! Evaluate expression expr for
                                    ! val double precision complex

character (len=*),intent(in) :: expr
complex(kc8) :: val
character (len=len(expr)+1) :: tempstr
character :: cop
integer :: isp(numtok)          ! On stack priority of operators in opstack
integer :: lstr
complex(kc8) :: cval,oper1,oper2
real(kr8) :: valr,vali
type(item):: token(numtok)      ! List of tokens ( a token is an operator or 
                                ! operand) in postfix order
type(item) :: x,junk,tok

ierr=0
token(1:)%char=' '

if(nparams == 0) then                  ! Initialize symbol table
  params(1)%symbol='PI'
  params(1)%value=(3.14159265358979_kr8,0.0_kr8)
  params(2)%symbol='I'
  params(2)%value=(0.0_kr8,1.0_kr8)
  nparams=2
end if

if(len_trim(expr) == 0) then           ! Expression empty
  ierr=1
  write(*,*) 'Error: expression being evaluated is empty'
  return
end if

tempstr=adjustl(expr)
call removesp(tempstr)   ! Removes spaces, tabs, and control characters

! ****************************************************************************
! STEP 1:  Convert string to token array. Each token is either an operator or
!          an operand. Token array will be in postfix (reverse Polish) order.
!*****************************************************************************

ntok=0
ibin=0
itop=0
do
  lstr=len_trim(tempstr)
  call get_next_token(tempstr(1:lstr),tok,icp,insp)
  select case(tok%type)
  case('S')
    ntok=ntok+1
    token(ntok)=tok
  case('E')
    do 
      if(itop < 1)exit
      call popop(x)        ! Output remaining operators on stack
      ntok=ntok+1
      token(ntok)=x
    end do
    ntok=ntok+1
    token(ntok)=tok
    exit
  case('R')  ! Token is right parenenthesis
    do 
      if(opstack(itop)%type == 'L') exit  ! Output operators on stack down
      call popop(x)                       ! to left parenthesis
      ntok=ntok+1
      token(ntok)=x
    end do                             
    call popop(junk)                      ! Remove left parenthesis from stack
    if(opstack(itop)%type == 'F') then    ! Output function name if present
      call popop(x)
      ntok=ntok+1
      token(ntok)=x
    end if                               
  case('D')  ! Token is comma
    do 
      if(opstack(itop)%type == 'L') exit  ! Output operators on stack down
      call popop(x)                       ! to left parenthesis
      ntok=ntok+1
      token(ntok)=x
    end do                              
  case('U','B','L','F') ! Token is operator, left parenthesis or function name
    do 
      if(isp(itop) < icp) exit            ! Output operators on stack having
      call popop(x)                       ! an instack priority that is
      ntok=ntok+1                         ! greater than or equal to the
      token(ntok)=x                       ! priority of the incoming operator  
    end do                                
    call pushop(tok)     ! Put incoming operator on stack                       
    isp(itop)=insp
  end select
end do

isum=0                                 ! Error check for matching parentheses
do i=1,ntok
  if(token(i)%type == 'L' ) isum=isum+1
  if(token(i)%type == 'R' ) isum=isum-1
end do
if(isum /= 0) then
  ierr=2
  write(*,*) 'Error in the evaluation of the expression ',trim(expr)
  write(*,*) "Parentheses don't match"
  write(*,*)
  return
end if


!*****************************************************************************
! STEP 2: Evaluate token string in postfix order
!*****************************************************************************

itop=0
do i=1,ntok
  x=token(i)
  select case(x%type)
  case('E')  ! Token is end token
    if(itop>1) then                
      ierr=12
      write(*,*) 'Error: missing operator in expression ',trim(expr)
      write(*,*)
      return
    end if
    call popval(val)               ! Final result left on stack of values
    exit
  case('S')  ! Token is operand
    call valuep(x%char,cval)       ! Evaluate operand
    if(ierr/=0) return
    call pushval(cval)             ! Put value of operand on stack
  case('B')  ! Token is a binary operator
    if(itop < 2) then
      ierr=5
      write(*,*) 'Error in evaluation of expression ',trim(expr)
      write(*,*) 'Less than two operands for binary operator  '&
                 ,trim(x%char)
      write(*,*)
      return
    end if                         
    call popval(oper1)             ! Pull off top two values from stack
    call popval(oper2)
    select case(trim(x%char))      ! Perform operation on values
    case('^')
      cval=oper2**oper1
    case('*')
      cval=oper2*oper1
    case('/')
      if(oper1 == (0._kr8,0._kr8)) then
        ierr=10
        write(*,*) 'Error in expression ',trim(expr)
        write(*,*) 'Division by zero'
        write(*,*)
        return
      end if
      cval=oper2/oper1
    case('+')
      cval=oper2+oper1
    case('-')
      cval=oper2-oper1
    end select
    call pushval(cval)             ! Put result back on stack
  case('U')  ! Token is unary operator
    if(itop == 0) then
      ierr=6
      write(*,*) 'Error in expression ',trim(expr)
      write(*,*) 'No operand for unary operator ',trim(x%char)
      write(*,*)
      return
    else
      call popval(oper1)           ! Pull top value off stack
    end if
    select case(trim(x%char))      ! Operate on value
    case('+')
      cval=oper1
    case('-')
      cval=-oper1
    end select
    call pushval(cval)             ! Put result back on stack
  case('F')  ! Token is a function name
    if(itop == 0) then
      ierr=7
      write(*,*) 'Error in expression ',trim(expr)
      write(*,*) 'Missing argument(s) for function ',trim(x%char)
      write(*,*)
      return
    else  
      call popval(oper1)           ! Pull top value off stack
    end if 
    tempstr=uppercase(x%char)
    select case(trim(tempstr))      ! Evaluate function
    case('SIN')
      cval=sin(oper1)
    case('COS')
      cval=cos(oper1)
    case('TAN')
      oper2=cos(oper1)
      if(abs(oper2) == 0.0_kr8) then
        ierr=14
        write(*,*) 'Error: argument of tan function a multiple',&
        ' of pi/2 in expression ',trim(expr)
        write(*,*)
        return
      else 
        cval=sin(oper1)/oper2
      endif
    case('SQRT')
      if(real(oper1,kr8) < 0. .and. aimag(oper1)==0.) then
        ierr=9
        write(*,*) 'Warning: square root of negative real number',&
                   ' in expression ',trim(expr)
        write(*,*)
      end if
      cval=sqrt(oper1)
    case('ABS')
      cval=abs(oper1)
    case('LN')
      if(real(oper1,kr8) <= 0. .and. aimag(oper1)==0.) then
        ierr=8
        write(*,*) 'Error: negative real or zero argument for',&
                   ' natural logarithm in expression ',trim(expr)
        write(*,*)
        return
      end if
      cval=log(oper1)
    case('LOG')
      if(real(oper1,kr8) <= 0. .and. aimag(oper1)==0.) then
        ierr=8
        write(*,*) 'Error: negative real or zero argument for base',&
                   '10 logarithm in expression ',trim(expr)
        write(*,*)
        return
      end if
      cval=log(oper1)/2.30258509299405_kr8
    case('EXP')
      cval=exp(oper1)
    case('COMPLEX')
      if(itop == 0) then
        ierr=7
        write(*,*) 'Error in expression ',trim(expr)
        write(*,*) 'Missing argument(s) for function ',trim(x%char)
        write(*,*)
        return
      else  
        call popval(oper2)  ! Pull second argument off stack
      end if 
      valr=real(oper2,kr8)
      vali=real(oper1,kr8)
      cval=cmplx(valr,vali,kc8)
    case('CONJG')
      cval=conjg(oper1)
    case('ANG')
      cval=atan2(aimag(oper1),real(oper1,kr8))
    case('REAL')
      cval=real(oper1,kr8)
    case('IMAG')
      cval=aimag(oper1)
    case default ! Undefined function
      ierr=13
      write(*,*) 'Error: the function ',trim(x%char), ' is undefined',&
                 ' in the expression ',trim(expr)
      write(*,*)
      return
    end select
    call pushval(cval)    ! Put result back on stack
  end select
end do

end subroutine evalexpr_dc

!**********************************************************************

SUBROUTINE GET_NEXT_TOKEN(str,tok,icp,isp)

character(len=*) :: str
character :: cop,chtemp
type(item) :: tok
integer :: icp

lstr=len_trim(str)
if(lstr == 0) then
  tok%char='#'             ! Output end token
  tok%type='E'
  return
end if
ipos=scan(str,'+-*/^(),')  ! Look for an arithmetic operator 
                           ! + - * / ^ ( ) or ,
cop=str(ipos:ipos)                 
select case (ipos)              
case(0)    ! Operators not present
  ntok=ntok+1
  tok%char=str
  tok%type='S'
  str=''
  icp=0
  isp=0
case(1) 
  tok%char=cop
  select case(cop)
  case('+','-')
    if(ibin==0) then
      tok%type='U'
      icp=4
      isp=3
    else
      tok%type='B'
      icp=1
      isp=1
    end if
    ibin=0
  case('*','/')
    tok%type='B'
    icp=2
    isp=2
    ibin=0
  case('^')
    tok%type='B'
    icp=4
    isp=3
    ibin=0
  case('(')
    tok%type='L'
    icp=4
    isp=0
    ibin=0
  case(')')
    tok%type='R'
    icp=0
    isp=0
    ibin=1
  case(',')
    tok%type='D'
    icp=0
    isp=0
    ibin=0
  end select
  str=str(2:)
case(2:)
  select case(cop)
  case('(')
    tok%char=str(1:ipos-1)
    tok%type='F'
    icp=4
    isp=0
    ibin=0
    str=str(ipos:)
  case('+','-')
    chtemp=uppercase(str(ipos-1:ipos-1))
    if(is_letter(str(1:1))==.true. .or. chtemp/='E') then
      tok%char=str(1:ipos-1)
      tok%type='S'
      icp=0
      isp=0
      ibin=1
      str=str(ipos:)
    else
      inext=scan(str(ipos+1:),'+-*/^(),')
      if(inext==0) then
        tok%char=str
        tok%type='S'
        icp=0
        isp=0
        ibin=0
        str=''
      else
        tok%char=str(1:ipos+inext-1)
        tok%type='S'
        icp=0
        isp=0
        ibin=1
        str=str(ipos+inext:)
      end if
    end if
  case default
    tok%char=str(1:ipos-1)
    tok%type='S'
    icp=0
    isp=0
    ibin=1
    str=str(ipos:)
  end select
end select

end subroutine get_next_token


!**********************************************************************

SUBROUTINE EVALEXPR_SC(expr,val)    ! Evaluate expression expr for
                                    ! val single precision complex
character(len=*) :: expr
complex(kc4) :: val
complex(kc8) :: vald

call evalexpr_dc(expr,vald)
val=vald

end subroutine evalexpr_sc

!**********************************************************************

SUBROUTINE EVALEXPR_SR(expr,val)    ! Evaluate expression expr for
                                    ! val single precision real
character(len=*) :: expr
real(kr4) :: val
complex(kc8) :: vald

call evalexpr_dc(expr,vald)
val=real(vald)

end subroutine evalexpr_sr

!**********************************************************************

SUBROUTINE EVALEXPR_DR(expr,val)    ! Evaluate expression expr for 
                                    ! val double precision real
character(len=*) :: expr
real(kr8) :: val
complex(kc8) :: vald

call evalexpr_dc(expr,vald)
val=real(vald,kr8)

end subroutine evalexpr_dr

!**********************************************************************

SUBROUTINE EVALEXPR_SI(expr,ival)   ! Evaluate expression expr for 
                                    ! ival single precision integer
character(len=*) :: expr
integer(ki4) :: ival
complex(kc8) :: vald

call evalexpr_dc(expr,vald)
ival=nint(real(vald,kr8),ki4)

end subroutine evalexpr_si

!**********************************************************************

SUBROUTINE EVALEXPR_DI(expr,ival)   ! Evaluate expression expr for 
                                    ! ival double precision integer
character(len=*) :: expr
integer(ki8) :: ival
complex(kc8) :: vald

call evalexpr_dc(expr,vald)
ival=nint(real(vald,kr8),ki8)

end subroutine evalexpr_di


!**********************************************************************
SUBROUTINE VALDEF_DC(sym,val)    ! Associates sym with val in symbol table,
                                 ! val double precision complex
character(len=*) :: sym
character(len=len_trim(sym)) :: usym
complex(kc8) :: val

ierr=0
if(nparams == 0) then               ! Initialize symbol table
  params(1)%symbol='PI'
  params(1)%value=(3.14159265358979_kr8,0.0_kr8)
  params(2)%symbol='I'
  params(2)%value=(0.0_kr8,1.0_kr8)
  nparams=2
end if

! Assign val to sym if sym is already in symbol table
usym=uppercase(sym)
if(is_letter(sym(1:1))==.false. .or. len_trim(sym)>24) then
  ierr=11
  write(*,*) 'Error: symbol ',trim(sym),' has improper format'
  write(*,*)
  return
end if
do i=1,nparams
  if(trim(usym)==trim(params(i)%symbol)) then
    params(i)%value=val
    return
  end if
end do

nparams=nparams+1    ! Otherwise assign val to new symbol sym
params(nparams)%symbol=usym
params(nparams)%value=val

end subroutine valdef_dc


!**********************************************************************

SUBROUTINE VALDEF_SC(sym,val)     ! Associates sym with val in symbol table,
                                  ! val single precision complex
character(len=*) :: sym
complex(kc4) :: val
complex(kc8) :: vald

vald=val
call valdef_dc(sym,vald)

end subroutine valdef_sc


!**********************************************************************

SUBROUTINE VALDEF_DR(sym,val)    ! Associates sym with val in symbol table,
                                 ! val double precision real
character(len=*) :: sym
real(kr8) :: val
complex(kc8) :: vald

vald=cmplx(val,0.0_kr8,kc8)
call valdef_dc(sym,vald)

end subroutine valdef_dr


!**********************************************************************

SUBROUTINE VALDEF_SR(sym,val)    ! Associates sym with val in symbol table,
                                 ! val single precision real
character(len=*) :: sym
real(kr4) :: val
complex(kc8) :: vald

vald=cmplx(val,0.0,kc8)
call valdef_dc(sym,vald)

end subroutine valdef_sr


!**********************************************************************

SUBROUTINE VALDEF_DI(sym,ival)   ! Associates sym with ival in symbol table,
                                 ! ival double precision integer 
character(len=*) :: sym
integer(ki8) :: ival
complex(kc8) :: vald

vald=cmplx(real(ival,kr8),0.0_kr8,kc8)
call valdef_dc(sym,vald)

end subroutine valdef_di


!**********************************************************************

SUBROUTINE VALDEF_SI(sym,ival)   ! Associates sym with ival in symbol table,
                                 ! ival single precision integer
character(len=*) :: sym
integer(ki4) :: ival
complex(kc8) :: vald

vald=cmplx(real(ival,kr8),0.0,kc8)
call valdef_dc(sym,vald)

end subroutine valdef_si


!**********************************************************************

SUBROUTINE STRDEF(sym,expr)      ! Associates sym with the value of the
                                 ! expression expr

character(len=*) :: sym,expr
complex(kc8) :: val

if(nparams == 0) then            ! Initialize symbol table
  params(1)%symbol='PI'
  params(1)%value=(3.14159265358979_kr8,0.0_kr8)
  params(2)%symbol='I'
  params(2)%value=(0.0_kr8,1.0_kr8)
  nparams=2
end if

call evalexpr_dc(expr,val)       ! val is value of expression expr
if(ierr==0 .or. ierr==9) then
  call valdef_dc(sym,val)          ! Assign val to symbol sym
end if

end subroutine strdef


!**********************************************************************

SUBROUTINE VALUEP(xinchar,cval)  ! Finds double precision complex value 
                                 ! corresponding to number string xinchar 
                                 ! or value in symbol table corresponding 
                                 ! to symbol name xinchar.

character (len=*):: xinchar
complex(kc8) :: cval
real(kr8) :: rval

ierr=0

if(is_letter(xinchar(1:1))==.true.) then   ! xinchar is a symbol
  call getparam(xinchar,cval)
else                               ! xinchar is a number string
  call value(xinchar,rval,ios)     ! rval is the value of xinchar
  if(ios > 0) then
    ierr=3
    write(*,*) 'Error: number string ',trim(xinchar),' does not correspond to a valid number' 
    write(*,*)
  end if
  cval=cmplx(rval,0.0_kr8,kc8)
  return
end if

end subroutine valuep


!**********************************************************************


SUBROUTINE PUSHOP(op)  ! Puts an operator on operator stack

type(item):: op

itop=itop+1
if(itop > numtok) then
  write(*,*) 'Error: operator stack overflow in evaluation of expression'
  write(*,*)
  return
end if
opstack(itop)=op

end subroutine pushop

SUBROUTINE POPOP(op) ! Takes top operator of operator stack and assigns it to op

type(item):: op

op=opstack(itop)
itop=itop-1

end subroutine popop

SUBROUTINE PUSHVAL(val) ! Puts value on value stack

complex(kc8) :: val

itop=itop+1
if(itop > numtok) then
  write(*,*) 'Error: value stack overflow in evaluation of expression'
  write(*,*)
  return
end if
valstack(itop)=val

end subroutine pushval

SUBROUTINE POPVAL(val) ! Takes top value off value stack and assigns it to val

complex(kc8) :: val

val=valstack(itop)
itop=itop-1

end subroutine popval

!**********************************************************************

SUBROUTINE GETPARAM_DC(sym,var)  ! Find double precision complex value var
                                 ! corresponding to symbol sym

character(len=*) :: sym
character(len=len_trim(sym)) :: usym
complex(kc8) :: var

ierr=0
sym=adjustl(sym)
if(is_letter(sym(1:1))==.false. .or. len_trim(sym)>24) then
  ierr=11
  write(*,*) 'Error: symbol ',trim(sym),' has incorrect format'
  write(*,*)
  return
end if
ifind=0
usym=uppercase(sym)
do j=1,nparams
  if(trim(usym) == trim(params(j)%symbol)) then
    var=params(j)%value
    ifind=j
    exit
  end if
end do
if(ifind == 0) then          
  ierr=4
  write(*,*) 'Error: symbol ',trim(sym), ' not in symbol table'
  write(*,*) 
  return
end if

end subroutine getparam_dc

!**********************************************************************

SUBROUTINE GETPARAM_SC(sym,var)  ! Find single precision complex value var
                                 ! corresponding to symbol sym


character(len=*) :: sym
complex(kc4) :: var
complex(kc8) :: vard

call getparam_dc(sym,vard)
var=vard

end subroutine getparam_sc

!**********************************************************************

SUBROUTINE GETPARAM_DR(sym,var)  ! Find double precision real value var
                                 ! corresponding to symbol sym


character(len=*) :: sym
real(kr8) :: var
complex(kc8) :: vard

call getparam_dc(sym,vard)
var=real(vard,kr8)

end subroutine getparam_dr

!**********************************************************************

SUBROUTINE GETPARAM_SR(sym,var)  ! Find single precision real value var
                                 ! corresponding to symbol sym


character(len=*) :: sym
real(kr4) :: var
complex(kc8) :: vard

call getparam_dc(sym,vard)
var=real(vard)

end subroutine getparam_sr

!**********************************************************************

SUBROUTINE GETPARAM_DI(sym,ivar)  ! Find double precision integer value ivar
                                  ! corresponding to symbol sym


character(len=*) :: sym
integer(ki8) :: ivar
complex(kc8) :: vard

call getparam_dc(sym,vard)
ivar=nint(real(vard,kr8),ki8)

end subroutine getparam_di

!**********************************************************************

SUBROUTINE GETPARAM_SI(sym,ivar)  ! Find single precision integer value ivar
                                  ! corresponding to symbol sym


character(len=*) :: sym
integer(ki4) :: ivar
complex(kc8) :: vard

call getparam_dc(sym,vard)
ivar=nint(real(vard,kr8),ki4)

end subroutine getparam_si

!**********************************************************************

SUBROUTINE EVALEQN(eqn)  ! Evaluate an equation

character(len=*) :: eqn
character(len=len(eqn)) :: args(2)
complex(kc8) :: val

call parse(eqn,'=',args,nargs)   ! Seperate right- and left-hand-sides
call defparam(adjustl(args(1)),args(2)) ! Evaluate right-hand-side and
                                        ! assign to symbol on the
                                        ! left-hand-side.
end subroutine evaleqn

!**********************************************************************

SUBROUTINE LISTVAR      ! List all variables and their values

write(*,'(/a)') ' VARIABLE LIST:'
if(nparams == 0) then            ! Initialize symbol table
  params(1)%symbol='PI'
  params(1)%value=(3.14159265358979_kr8,0.0_kr8)
  params(2)%symbol='I'
  params(2)%value=(0.0_kr8,1.0_kr8)
  nparams=2
end if
do i=1,nparams
  write(*,*) trim(params(i)%symbol),' = ',params(i)%value
end do

end subroutine listvar

!**********************************************************************

end module evaluate
