module mod_interfaces
use iso_fortran_env, only:i1 => int8, i2 => int16, i4 => int32, i8 => int64, &
                         r4 => real32, r8 => real64, r16 => real128
implicit none

private

public :: f_i1, f_i2, f_i4, f_i8
public :: f_r4, f_r8, f_r16
public :: f_c4, f_c8, f_c16
public :: f_array_i1, f_array_i2, f_array_i4, f_array_i8
public :: f_array_r4, f_array_r8, f_array_r16
public :: f_array_c4, f_array_c8, f_array_c16
public :: f2_i1, f2_i2, f2_i4, f2_i8
public :: f2_r4, f2_r8, f2_r16
public :: f2_c4, f2_c8, f2_c16
public :: f_i1_logical, f_i2_logical, f_i4_logical, f_i8_logical
public :: f_r4_logical, f_r8_logical, f_r16_logical
public :: f_c4_logical, f_c8_logical, f_c16_logical
public :: f_char_logical

interface

  pure integer(i1) function f_i1(x)
    !! f :: i1 -> i1
    import :: i1
    integer(i1), intent(in) :: x
  end function f_i1

  pure integer(i2) function f_i2(x)
    !! f :: i2 -> i2
    import :: i2
    integer(i2), intent(in) :: x
  end function f_i2

  pure integer(i4) function f_i4(x)
    !! f :: i4 -> i4
    import :: i4
    integer(i4), intent(in) :: x
  end function f_i4

  pure integer(i8) function f_i8(x)
    !! f :: i8 -> i8
    import :: i8
    integer(i8), intent(in) :: x
  end function f_i8

  pure real(r4) function f_r4(x)
    !! f :: r4 -> r4
    import :: r4
    real(r4), intent(in) :: x
  end function f_r4

  pure real(r8) function f_r8(x)
    !! f :: r8 -> r8
    import :: r8
    real(r8), intent(in) :: x
  end function f_r8

  pure real(r16) function f_r16(x)
    !! f :: r16 -> r16
    import :: r16
    real(r16), intent(in) :: x
  end function f_r16

  pure complex(r4) function f_c4(x)
    !! f :: c4 -> c4
    import :: r4
    complex(r4), intent(in) :: x
  end function f_c4

  pure complex(r8) function f_c8(x)
    !! f :: c8 -> c8
    import :: r8
    complex(r8), intent(in) :: x
  end function f_c8

  pure complex(r16) function f_c16(x)
    !! f :: c16 -> c16
    import :: r16
    complex(r16), intent(in) :: x
  end function f_c16

  pure function f_array_i1(x) result(f)
    !! f :: [i1] -> [i1]
    import :: i1
    integer(i1), dimension(:), intent(in) :: x
    integer(i1), dimension(:), allocatable :: f
  end function f_array_i1

  pure function f_array_i2(x) result(f)
    !! f :: [i2] -> [i2]
    import :: i2
    integer(i2), dimension(:), intent(in) :: x
    integer(i2), dimension(:), allocatable :: f
  end function f_array_i2

  pure function f_array_i4(x) result(f)
    !! f :: [i4] -> [i4]
    import :: i4
    integer(i4), dimension(:), intent(in) :: x
    integer(i4), dimension(:), allocatable :: f
  end function f_array_i4

  pure function f_array_i8(x) result(f)
    !! f :: [i8] -> [i8]
    import :: i8
    integer(i8), dimension(:), intent(in) :: x
    integer(i8), dimension(:), allocatable :: f
  end function f_array_i8

  pure function f_array_r4(x) result(f)
    !! f :: [r4] -> [r4]
    import :: r4
    real(r4), dimension(:), intent(in) :: x
    real(r4), dimension(:), allocatable :: f
  end function f_array_r4

  pure function f_array_r8(x) result(f)
    !! f :: [r8] -> [r8]
    import :: r8
    real(r8), dimension(:), intent(in) :: x
    real(r8), dimension(:), allocatable :: f
  end function f_array_r8

  pure function f_array_r16(x) result(f)
    !! f :: [r16] -> [r16]
    import :: r16
    real(r16), dimension(:), intent(in) :: x
    real(r16), dimension(:), allocatable :: f
  end function f_array_r16

  pure function f_array_c4(x) result(f)
    !! f :: [c4] -> [c4]
    import :: r4
    complex(r4), dimension(:), intent(in) :: x
    complex(r4), dimension(:), allocatable :: f
  end function f_array_c4

  pure function f_array_c8(x) result(f)
    !! f :: [c8] -> [c8]
    import :: r8
    complex(r8), dimension(:), intent(in) :: x
    complex(r8), dimension(:), allocatable :: f
  end function f_array_c8

  pure function f_array_c16(x) result(f)
    !! f :: [c16] -> [c16]
    import :: r16
    complex(r16), dimension(:), intent(in) :: x
    complex(r16), dimension(:), allocatable :: f
  end function f_array_c16

  pure integer(i1) function f2_i1(x, y)
    !! f :: i1 i1 -> i1
    import :: i1
    integer(i1), intent(in) :: x, y
  end function f2_i1

  pure integer(i2) function f2_i2(x, y)
    !! f :: i2 i2 -> i2
    import :: i2
    integer(i2), intent(in) :: x, y
  end function f2_i2

  pure integer(i4) function f2_i4(x, y)
    !! f :: i4 i4 -> i4
    import :: i4
    integer(i4), intent(in) :: x, y
  end function f2_i4

  pure integer(i8) function f2_i8(x, y)
    !! f :: i8 i8 -> i8
    import :: i8
    integer(i8), intent(in) :: x, y
  end function f2_i8

  pure real(r4) function f2_r4(x, y)
    !! f :: r4 r4 -> r4
    import :: r4
    real(r4), intent(in) :: x, y
  end function f2_r4

  pure real(r8) function f2_r8(x, y)
    !! f :: r8 r8 -> r8
    import :: r8
    real(r8), intent(in) :: x, y
  end function f2_r8

  pure real(r16) function f2_r16(x, y)
    !! f :: r16 r16 -> r16
    import :: r16
    real(r16), intent(in) :: x, y
  end function f2_r16

  pure complex(r4) function f2_c4(x, y)
    !! f :: c4 c4 -> c4
    import :: r4
    complex(r4), intent(in) :: x, y
  end function f2_c4

  pure complex(r8) function f2_c8(x, y)
    !! f :: c8 c8 -> c8
    import :: r8
    complex(r8), intent(in) :: x, y
  end function f2_c8

  pure complex(r16) function f2_c16(x, y)
    !! f :: c16 c16 -> c16
    import :: r16
    complex(r16), intent(in) :: x, y
  end function f2_c16

  pure logical function f_i1_logical(x)
    !! f :: i1 -> logical
    import :: i1
    integer(i1), intent(in) :: x
  end function f_i1_logical

  pure logical function f_i2_logical(x)
    !! f :: i2 -> logical
    import :: i2
    integer(i2), intent(in) :: x
  end function f_i2_logical

  pure logical function f_i4_logical(x)
    !! f :: i4 -> logical
    import :: i4
    integer(i4), intent(in) :: x
  end function f_i4_logical

  pure logical function f_i8_logical(x)
    !! f :: i8 -> logical
    import :: i8
    integer(i8), intent(in) :: x
  end function f_i8_logical

  pure logical function f_r4_logical(x)
    !! f :: r4 -> logical
    import :: r4
    real(r4), intent(in) :: x
  end function f_r4_logical

  pure logical function f_r8_logical(x)
    !! f :: r8 -> logical
    import :: r8
    real(r8), intent(in) :: x
  end function f_r8_logical

  pure logical function f_r16_logical(x)
    !! f :: r16 -> logical
    import :: r16
    real(r16), intent(in) :: x
  end function f_r16_logical

  pure logical function f_c4_logical(x)
    !! f :: c4 -> logical
    import :: r4
    complex(r4), intent(in) :: x
  end function f_c4_logical

  pure logical function f_c8_logical(x)
    !! f :: c8 -> logical
    import :: r8
    complex(r8), intent(in) :: x
  end function f_c8_logical

  pure logical function f_c16_logical(x)
    !! f :: c16 -> logical
    import :: r16
    complex(r16), intent(in) :: x
  end function f_c16_logical

  pure logical function f_char_logical(x)
    !! f :: character -> logical
    character(len=1), intent(in) :: x
  end function f_char_logical

end interface

end module mod_interfaces
