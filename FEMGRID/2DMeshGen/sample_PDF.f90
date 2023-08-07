module sample_pdf
    use meshDS,only:write_help,Err_msg,property,pro_num,path_name,strtoint
    implicit none
    
    private
    
    public sample_tydef,example_pdf
    
    !INTERFACE sample
    !    MODULE PROCEDURE uniform_pdf,nomral_pdf,general_pdf
    !END INTERFACE 
    
    type sample_tydef
        integer::isout=0
        integer::nval=0,nprob=0,ipdf=0
        real(8)::minv,maxv,mean,std
        real(8),allocatable::rangediv(:),prob(:)
        character(1024):: helpstring= &
            &"pdf_sample的功能是实现按照某分布进行取样. \n &
            & pdf_sample的输入格式为: \n &
            & pdf_sample,pdf=I,nsample=I[,nprob=I,isout=I] \n &
            & <parameters> \n &
            & NOTE: \n &
            & 0)pdf=0,均匀分布,对应的参数行数据parameters=[minv,maxv],共1行. \n &
            &   pdf=1,正态分布,对应的参数行数据parameters=[minv,maxv,mean,std],共1行. \n &
            &   pdf=2,一般的分布,对应的参数行数据parameters=[rangediv(nprob+1),prob(nprob)],分2行输入,rangediv及prob各1行. \n &
            & 2)nsample--要采用的样本数;\n & 
            & 3)nprob--概率分布区间的个数，仅当pdf=2时输入;\n & 
            & 4)minv--生成样本的最小值;\n & 
            & 5)maxv--生成样本的最大值;\n & 
            & 6)mean,std--正态分布的期望和标准差;\n & 
            & 7)rangediv(nprob+1)--一般分布的各区间的样本的限值,共nprob+1个;\n & 
            & 8)prob(nprob)--一般分布的各区间的出现的概率,共nprob个;\n &
            & 9)isout=0--不输出,默认为0;\n &
            &   isout=1--输出,且输出后退出程序;\n &
            & "C                
        
    contains
        procedure,nopass::help=>write_help
        procedure::readin=>sample_read
        procedure,nopass::init=>rng
        procedure::uniform_pdf,nomral_pdf,general_pdf
        procedure::output=>sample_output
        procedure::sample
    endtype
    
    type(sample_tydef)::example_pdf
    
    contains

    function sample(self) result(x)
      implicit none
      class(sample_tydef)::self
      real(8)::x(self.nval)

      select case(self.ipdf)
      case(0)
        x=self.uniform_pdf()
      case(1)
        x=self.nomral_pdf()
      case(2)
        x=self.general_pdf()        
      endselect

    end
    subroutine sample_output(self)
        implicit none
        class(sample_tydef)::self
        integer::i
        character(1024)::file
        real(8)::area1=0.0,t1
        
        file=trim(path_name)//'_pdf_sample.txt'
        
        open(unit=10,file=file,status='replace')
        select case(self.ipdf)
        case(0)
          write(10,20) self.minv,self.maxv,self.nval       
        case(1)
          write(10,30) self.minv,self.maxv,self.nval,self.mean,self.std
        case(2)
          write(10,42) self.nval
          write(10,40) self.rangediv
          write(10,41) self.prob
        ENDSELECT
        call self.init()
        write(10,50) self.sample()
        
        close(10)
        
        stop
        
    20  format('uniform pdf,[min,max,nsample]=[',2(f7.4,','),i7,']')
    30  format('normal pdf,[min,max,nsample,mean,std]=[',2(f7.4,','),i7,f7.4,',',f7.4,']')
    40  format('general pdf,rangediv=',<self.nprob+1>(f7.4,X))
    41  format('general pdf,prob=',<self.nprob>(f7.4,X))
    42  format('general pdf,nsample=',i7)
    50  format(<10>(f7.4,x)) 
    end
    subroutine sample_read(self,unit)
        implicit none
        class(sample_tydef)::self
        integer,intent(in)::unit
        integer::i,n1
        INTEGER::DN=0
        INTEGER,PARAMETER::DNMAX=1000
        REAL(8)::AR(DNMAX) 
        
        print *,'Reading sample_pdf data...'
        call self.help(self.helpstring)        
        do i=1, pro_num
            select case(property(i).name)
            case('pdf')
                self.ipdf=int(property(i).value)
            case('nsample')
	             self.nval=int(property(i).value)
            case('nprob')
	             self.nprob=int(property(i).value)
            case('isout')
	             self.isout=int(property(i).value)     
            case default
	            call Err_msg(property(i).name)
            end select
        end do 
        call strtoint(unit,ar,dnmax,dn,dnmax)
        select case(self.ipdf)
        case(0)
          self.minv=ar(1)
          self.maxv=ar(2)
        case(1)
          self.minv=ar(1)
          self.maxv=ar(2)
          self.mean=ar(3)
          self.std=ar(4)
        case(2)
          if(dn/=self.nprob+1) THEN
            error stop "wrong data number in readin general_pdf parameters1."
          ENDIF
          self.rangediv=ar(1:self.nprob+1)
          call strtoint(unit,ar,dnmax,dn,dnmax)
          if(dn/=self.nprob) THEN
            error stop "wrong data number in readin general_pdf parameters2."
          ENDIF          
          self.prob=ar(1:self.nprob) 
          if(abs(sum(self.prob)-1.0d0)>1e-6) then
            error stop "the sum probability is not equal 1."
          endif
          self.minv=self.rangediv(1);self.maxv=self.rangediv(self.nprob+1) 
        endselect

    endsubroutine
    
    function uniform_pdf(self)  result(x)
    !返回nval个在区间[minv,max]服从均匀分布的数。
        implicit none
        class(sample_tydef)::self
        real(8)::x(self.nval)
        
        call random_number(x)
        x=self.minv+x*(self.maxv-self.minv)        
    
    endfunction

    function nomral_pdf(self) result(X)
    !返回nval个在区间[minv,max]服从期望为mean，方差为std的正态分布的数。
        implicit none
        class(sample_tydef)::self
        real(8)::x(self.nval),x1
        integer::i

        do i=1,self.nval
            call truncated_normal_ab_sample ( self.mean, self.std, self.minv, self.maxv, x1)
            x(i)=x1
        enddo

    endfunction

    function general_pdf(self) result(X)
    !返回服从广义分布的nval个数。
    !rangediv:颗粒分布区间:v1,v2,...,vn,共n+1个
    !prob:各区间的概率：v1-v2,prob1;v2-v3,prob2;...,vn-1:vn=probn-1.共n个sum(probi)==1
    !
        implicit none
        class(sample_tydef)::self
        real(8)::x(self.nval)
        integer::i,j
        real(8)::r1,t1,prob1(0:self.nprob)

        
        
        prob1(0)=0.d0
        do i=1,self.nprob
            prob1(i)=prob1(i-1)+self.prob(i)
        enddo
        do j=1,self.nval
            call random_number(r1)
            do i=1,self.nprob            
                if(self.prob(i)>0.d0.and.r1<=prob1(i).and.r1>prob1(i-1)) then                    
                    call random_number(t1)
                    x(j:j)=self.rangediv(i)+t1*(self.rangediv(i+1)-self.rangediv(i))    
                    exit
                endif
            enddo
        enddo

    endfunction

    
subroutine truncated_normal_ab_sample ( mu, sigma, a, b, x )

!*****************************************************************************80
!
!! TRUNCATED_NORMAL_AB_SAMPLE samples the truncated Normal PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) MU, SIGMA, the mean and standard deviation of the
!    parent Normal distribution.
!
!    Input, real ( kind = 8 ) A, B, the lower and upper truncation limits.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, real ( kind = 8 ) X, a sample of the PDF.
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) alpha
  real ( kind = 8 ) alpha_cdf
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) beta_cdf
  real ( kind = 8 ) mu
  real ( kind = 8 ) r8_uniform_01
  real ( kind = 8 ) sigma
  !integer ( kind = 4 ) seed
  real ( kind = 8 ) u
  real ( kind = 8 ) x
  real ( kind = 8 ) xi
  real ( kind = 8 ) xi_cdf

  alpha = ( a - mu ) / sigma
  beta = ( b - mu ) / sigma

  call normal_01_cdf ( alpha, alpha_cdf )
  call normal_01_cdf ( beta, beta_cdf )

  !u = r8_uniform_01 ( seed ) //lgy
  call random_number(u)
  xi_cdf = alpha_cdf + u * ( beta_cdf - alpha_cdf )
  call normal_01_cdf_inv ( xi_cdf, xi )

  x = mu + sigma * xi

  return
end

subroutine normal_01_cdf ( x, cdf )

!*****************************************************************************80
!
!! NORMAL_01_CDF evaluates the Normal 01 CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    AG Adams,
!    Algorithm 39,
!    Areas Under the Normal Curve,
!    Computer Journal,
!    Volume 12, pages 197-198, 1969.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the CDF.
!
!    Output, real ( kind = 8 ) CDF, the value of the CDF.
!
  implicit none

  real ( kind = 8 ), parameter :: a1 = 0.398942280444D+00
  real ( kind = 8 ), parameter :: a2 = 0.399903438504D+00
  real ( kind = 8 ), parameter :: a3 = 5.75885480458D+00
  real ( kind = 8 ), parameter :: a4 = 29.8213557808D+00
  real ( kind = 8 ), parameter :: a5 = 2.62433121679D+00
  real ( kind = 8 ), parameter :: a6 = 48.6959930692D+00
  real ( kind = 8 ), parameter :: a7 = 5.92885724438D+00
  real ( kind = 8 ), parameter :: b0 = 0.398942280385D+00
  real ( kind = 8 ), parameter :: b1 = 3.8052D-08
  real ( kind = 8 ), parameter :: b2 = 1.00000615302D+00
  real ( kind = 8 ), parameter :: b3 = 3.98064794D-04
  real ( kind = 8 ), parameter :: b4 = 1.98615381364D+00
  real ( kind = 8 ), parameter :: b5 = 0.151679116635D+00
  real ( kind = 8 ), parameter :: b6 = 5.29330324926D+00
  real ( kind = 8 ), parameter :: b7 = 4.8385912808D+00
  real ( kind = 8 ), parameter :: b8 = 15.1508972451D+00
  real ( kind = 8 ), parameter :: b9 = 0.742380924027D+00
  real ( kind = 8 ), parameter :: b10 = 30.789933034D+00
  real ( kind = 8 ), parameter :: b11 = 3.99019417011D+00
  real ( kind = 8 ) cdf
  real ( kind = 8 ) q
  real ( kind = 8 ) x
  real ( kind = 8 ) y
!
!  |X| <= 1.28.
!
  if ( abs ( x ) <= 1.28D+00 ) then

    y = 0.5D+00 * x * x

    q = 0.5D+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
      + a6 / ( y + a7 ) ) ) )
!
!  1.28 < |X| <= 12.7
!
  else if ( abs ( x ) <= 12.7D+00 ) then

    y = 0.5D+00 * x * x

    q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
      + b2 / ( abs ( x ) + b3 &
      + b4 / ( abs ( x ) - b5 &
      + b6 / ( abs ( x ) + b7 &
      - b8 / ( abs ( x ) + b9 &
      + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
!
!  12.7 < |X|
!
  else

    q = 0.0D+00

  end if
!
!  Take account of negative X.
!
  if ( x < 0.0D+00 ) then
    cdf = q
  else
    cdf = 1.0D+00 - q
  end if

  return
end
subroutine normal_01_cdf_inv ( p, x )

!*****************************************************************************80
!
!! NORMAL_01_CDF_INV inverts the standard normal CDF.
!
!  Discussion:
!
!    The result is accurate to about 1 part in 10^16.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 February 2015
!
!  Author:
!
!    Original FORTRAN77 version by Michael Wichura.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Michael Wichura,
!    Algorithm AS241:
!    The Percentage Points of the Normal Distribution,
!    Applied Statistics,
!    Volume 37, Number 3, pages 477-484, 1988.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P, the value of the cumulative probability
!    densitity function.  0 < P < 1.  If P is outside this range, an
!    "infinite" value will be returned.
!
!    Output, real ( kind = 8 ) X, the normal deviate value
!    with the property that the probability of a standard normal deviate being
!    less than or equal to the value is P.
!
  implicit none

  real ( kind = 8 ), parameter, dimension ( 8 ) :: a = (/ &
    3.3871328727963666080D+00, &
    1.3314166789178437745D+02, &
    1.9715909503065514427D+03, &
    1.3731693765509461125D+04, &
    4.5921953931549871457D+04, &
    6.7265770927008700853D+04, &
    3.3430575583588128105D+04, &
    2.5090809287301226727D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: b = (/ &
    1.0D+00, &
    4.2313330701600911252D+01, &
    6.8718700749205790830D+02, &
    5.3941960214247511077D+03, &
    2.1213794301586595867D+04, &
    3.9307895800092710610D+04, &
    2.8729085735721942674D+04, &
    5.2264952788528545610D+03 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: c = (/ &
    1.42343711074968357734D+00, &
    4.63033784615654529590D+00, &
    5.76949722146069140550D+00, &
    3.64784832476320460504D+00, &
    1.27045825245236838258D+00, &
    2.41780725177450611770D-01, &
    2.27238449892691845833D-02, &
    7.74545014278341407640D-04 /)
  real ( kind = 8 ), parameter :: const1 = 0.180625D+00
  real ( kind = 8 ), parameter :: const2 = 1.6D+00
  real ( kind = 8 ), parameter, dimension ( 8 ) :: d = (/ &
    1.0D+00, &
    2.05319162663775882187D+00, &
    1.67638483018380384940D+00, &
    6.89767334985100004550D-01, &
    1.48103976427480074590D-01, &
    1.51986665636164571966D-02, &
    5.47593808499534494600D-04, &
    1.05075007164441684324D-09 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: e = (/ &
    6.65790464350110377720D+00, &
    5.46378491116411436990D+00, &
    1.78482653991729133580D+00, &
    2.96560571828504891230D-01, &
    2.65321895265761230930D-02, &
    1.24266094738807843860D-03, &
    2.71155556874348757815D-05, &
    2.01033439929228813265D-07 /)
  real ( kind = 8 ), parameter, dimension ( 8 ) :: f = (/ &
    1.0D+00, &
    5.99832206555887937690D-01, &
    1.36929880922735805310D-01, &
    1.48753612908506148525D-02, &
    7.86869131145613259100D-04, &
    1.84631831751005468180D-05, &
    1.42151175831644588870D-07, &
    2.04426310338993978564D-15 /)
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) r
  !real ( kind = 8 ) r8poly_value_horner
  real ( kind = 8 ), parameter :: split1 = 0.425D+00
  real ( kind = 8 ), parameter :: split2 = 5.0D+00
  real ( kind = 8 ) x

  if ( p <= 0.0D+00 ) then
    x = - huge ( x )
    return
  end if

  if ( 1.0D+00 <= p ) then
    x = huge ( x )
    return
  end if

  q = p - 0.5D+00

  if ( abs ( q ) <= split1 ) then

    r = const1 - q * q
    x = q * r8poly_value_horner ( 7, a, r ) &
          / r8poly_value_horner ( 7, b, r )

  else

    if ( q < 0.0D+00 ) then
      r = p
    else
      r = 1.0D+00 - p
    end if

    if ( r <= 0.0D+00 ) then

      x = huge ( x )

    else

      r = sqrt ( - log ( r ) )

      if ( r <= split2 ) then

        r = r - const2
        x = r8poly_value_horner ( 7, c, r ) &
          / r8poly_value_horner ( 7, d, r )

      else

        r = r - split2
        x = r8poly_value_horner ( 7, e, r ) &
          / r8poly_value_horner ( 7, f, r )

      end if

    end if

    if ( q < 0.0D+00 ) then
      x = -x
    end if

  end if

  return
end
function r8poly_value_horner ( m, c, x )

!*****************************************************************************80
!
!! R8POLY_VALUE_HORNER evaluates a polynomial using Horner's method.
!
!  Discussion:
!
!    The polynomial 
!
!      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m
!
!    is to be evaluated at the value X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the degree.
!
!    Input, real ( kind = 8 ) C(0:M), the polynomial coefficients.  
!    C(I) is the coefficient of X^I.
!
!    Input, real ( kind = 8 ) X, the evaluation point.
!
!    Output, real ( kind = 8 ) R8POLY_VALUE_HORNER, the polynomial value.
!
  implicit none

  integer ( kind = 4 ) m

  real ( kind = 8 ) c(0:m)
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8poly_value_horner
  real ( kind = 8 ) value
  real ( kind = 8 ) x

  value = c(m)
  do i = m - 1, 0, -1
    value = value * x + c(i)
  end do

  r8poly_value_horner = value

  return
end

!=======================================================================
! rng
!-----------------------------------------------------------------------
! rng controls random number generation.
!
! Syntax
!-----------------------------------------------------------------------
! call rng()
! call rng(seed)
!
! Description
!-----------------------------------------------------------------------
! call rng() uses the current date and time as seed for random number
! generation.
!
! call rng(seed) sets the input seed for random number generation.
!
! Notes
!-----------------------------------------------------------------------
! It is advised to call rng at the beginning of a program so that each
! run of the program produces different sequences of random numbers.
!=======================================================================

  subroutine rng(seed)
    integer, parameter :: IPRE = 4
    integer(kind = IPRE), intent(in), optional :: seed
    integer(kind = 4) :: seed_size, values(8)
    integer(kind = 4), dimension(:), allocatable :: seed_put

    call random_seed(size = seed_size)
    allocate(seed_put(seed_size))
    if (present(seed)) then
      seed_put = seed
    else
      call date_and_time(values = values)
      seed_put = values(8) * values(7) * values(6)
    end if
    call random_seed(put = seed_put)
    return
  end subroutine rng
    
end module
    