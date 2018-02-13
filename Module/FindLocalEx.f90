module FindLocalEx1D
    
contains
  !PEAKDET Detect peaks in a vector
  !
  !        call PEAKDET(MAXTAB, MINTAB, N, V, DELTA) finds the local
  !        maxima and minima ("peaks") in the vector V of size N.
  !        MAXTAB and MINTAB consists of two columns. Column 1
  !        contains indices in V, and column 2 the found values.
  !
  !        call PEAKDET(MAXTAB, MINTAB, N, V, DELTA, X) replaces the 
  !        indices in MAXTAB and MINTAB with the corresponding X-values.
  !
  !        A point is considered a maximum peak if it has the maximal
  !        value, and was preceded (to the left) by a value lower by
  !        DELTA.
  !
  ! Eli Billauer, 3.4.05 (http://billauer.co.il)
  ! Translated into Fortran by Brian McNoldy (http://andrew.rsmas.miami.edu/bmcnoldy)
  ! This function is released to the public domain; Any use is allowed.

  subroutine peakdet(maxtab,mintab,n,v,delta,x)

    use, intrinsic :: ieee_arithmetic
    implicit none

    integer, intent(in)            :: n
    real, intent(in)               :: v(n), delta
    real, intent(in), optional     :: x(n)
    real, intent(out), allocatable :: maxtab(:,:), mintab(:,:)
    integer                        :: nargin, lookformax, i, j, c, d
    real                           :: a, NaN, Pinf, Minf, &
                                      mn, mx, mnpos, mxpos, this, &
                                      x2(n), maxtab_tmp(2,n), mintab_tmp(2,n)

    !nargin=command_argument_count()
    if (.not.present(x) ) then
      forall(j=1:n) x2(j)=dble(j)
    else
      x2=x
      if (size(v) /= size(x)) then
        print*,'Input vectors v and x must have same length'
      end if
    end if

    if (size((/ delta /)) > 1) then
      print*,'Input argument DELTA must be a scalar'
    end if

    if (delta <= 0) then
      print*,'Input argument DELTA must be positive'
    end if

    NaN=ieee_value(a, ieee_quiet_nan)
    Pinf=ieee_value(a, ieee_positive_inf)
    Minf=ieee_value(a, ieee_negative_inf)

    mn=Pinf
    mx=Minf
    mnpos=NaN
    mxpos=NaN
    lookformax=1

    c=0
    d=0
    maxtab_tmp(:,:)=NaN
    mintab_tmp(:,:)=NaN
    do i=1,n
      this = v(i)
      if (this > mx) then
        mx = this
        mxpos = x2(i)
      end if
      if (this < mn) then
        mn = this
        mnpos = x2(i)
      end if
      if (lookformax==1) then
        if (this < mx-delta) then
          c=c+1
          maxtab_tmp(:,C)=(/ mxpos, mx /)
          mn = this
          mnpos = x2(i)
          lookformax = 0
        end if
      else
        if (this > mn+delta) then
          d=d+1
          mintab_tmp(:,d)=(/ mnpos, mn /)
          mx = this
          mxpos = x2(i)
          lookformax = 1
        end if
      end if
    end do

    allocate(maxtab(2,c))
    allocate(mintab(2,d))
    where (.not.ieee_is_nan(maxtab_tmp))
      maxtab=maxtab_tmp
    end where
    where (.not.ieee_is_nan(mintab_tmp))
      mintab=mintab_tmp
    end where
  end subroutine peakdet

end module FindLocalEx1D