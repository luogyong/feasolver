module quicksort

! sort routine to arrange array elements from smallest to largest
!
! grabbed from A millers web site http://users.bigpond.net.au/amiller/
! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives
! the positions of the elements in the original order. a(order(1))=smallest... order=order2node
! pjr added module declaration
! mvr modified integer array to intent inout - may now be any integer 
!     array that gets sorted along with associated REAL(8) array
!added by lgy: ulist sorted unique list  and n2ulist index array
implicit none
save
private
public quick_sort

INTERFACE quick_sort
    MODULE PROCEDURE quick_sort_INT,quick_sort_REAL,unique_int,unique_real
END INTERFACE


contains

subroutine unique_real(list,ulist,n2ulist,o2n)    
    implicit none

    REAL(8), DIMENSION (:), INTENT(INOUT)  :: list
    INTEGER, OPTIONAL,DIMENSION (:), INTENT(INOUT)  :: o2n
    INTEGER,DIMENSION (:), INTENT(OUT)  :: n2ulist
    REAL(8),allocatable,intent(out)::ulist(:)
    
    integer,allocatable::o2n2(:)
    
    integer::n1,n2,j
    if(present(o2n)) then
        call quick_sort_REAL(list, o2n)
    else
        call quick_sort_REAL(list)
    endif
    
    ulist=list
    n1=size(ulist)
    n2=1
    O2N2=[1:n1]
    do j=2,n1
        if(ulist(j)>ulist(j-1)) then
            n2=n2+1
            if(j>n2) ulist(n2)=ulist(j)
        endif
        O2N2(j)=n2
    enddo
    
    n2ulist(o2n)=O2N2    
    ulist=ulist(1:n2)
    
endsubroutine

subroutine unique_int(list,ulist,n2ulist,o2n)    
    implicit none

    INTEGER, DIMENSION (:), INTENT(INOUT)  :: list
    INTEGER, OPTIONAL,DIMENSION (:), INTENT(INOUT)  :: o2n
    INTEGER,DIMENSION (:), INTENT(OUT)  :: n2ulist
    INTEGER,allocatable,intent(out)::ulist(:)
    
    integer,allocatable::o2n2(:)
    
    integer::n1,n2,j
    if(present(o2n)) then
        call quick_sort_int(list, o2n)
    else
        call quick_sort_int(list)
    endif
    
    ulist=list
    n1=size(ulist)
    n2=1
    O2N2=[1:n1]
    do j=2,n1
        if(ulist(j)>ulist(j-1)) then
            n2=n2+1
            if(j>n2) ulist(n2)=ulist(j)
        endif
        O2N2(j)=n2
    enddo
    
    n2ulist(o2n)=O2N2    
    ulist=ulist(1:n2)
    
endsubroutine

RECURSIVE SUBROUTINE quick_sort_REAL(list, order)

implicit none

REAL(8), DIMENSION (:), INTENT(INOUT)  :: list
INTEGER, OPTIONAL,DIMENSION (:), INTENT(INOUT)  :: order

! Local variable
INTEGER :: i

CALL quick_sort_1(1, SIZE(list))

CONTAINS


RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

implicit none
INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
REAL(8)                :: reference, temp
INTEGER, PARAMETER  :: max_simple_sort_size = 6

IF (right_end < left_end + max_simple_sort_size) THEN
  ! Use interchange sort for small lists
  CALL interchange_sort(left_end, right_end)

ELSE
  ! Use partition ("quick") sort
  reference = list((left_end + right_end)/2)
  i = left_end - 1; j = right_end + 1

  DO
    ! Scan list from left end until element >= reference is found
    DO
      i = i + 1
      IF (list(i) >= reference) EXIT
    END DO
    ! Scan list from right end until element <= reference is found
    DO
      j = j - 1
      IF (list(j) <= reference) EXIT
    END DO


    IF (i < j) THEN
      ! Swap two out-of-order elements
      temp = list(i); list(i) = list(j); list(j) = temp
      IF(PRESENT(ORDER)) THEN
        itemp = order(i); order(i) = order(j); order(j) = itemp
      ENDIF
    ELSE IF (i == j) THEN
      i = i + 1
      EXIT
    ELSE
      EXIT
    END IF
  END DO

  IF (left_end < j) CALL quick_sort_1(left_end, j)
  IF (i < right_end) CALL quick_sort_1(i, right_end)
END IF

END SUBROUTINE quick_sort_1



SUBROUTINE interchange_sort(left_end, right_end)

implicit none
INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
REAL(8)                :: temp

DO i = left_end, right_end - 1
  DO j = i+1, right_end
    IF (list(i) > list(j)) THEN
      temp = list(i); list(i) = list(j); list(j) = temp
      IF(PRESENT(ORDER)) THEN
        itemp = order(i); order(i) = order(j); order(j) = itemp
      ENDIF
    END IF
  END DO
END DO

END SUBROUTINE interchange_sort

END SUBROUTINE quick_sort_REAL

RECURSIVE SUBROUTINE quick_sort_INT(list, order)

implicit none

INTEGER, DIMENSION (:), INTENT(INOUT)  :: list
INTEGER, OPTIONAL,DIMENSION (:), INTENT(INOUT)  :: order

! Local variable
INTEGER :: i

CALL quick_sort_1(1, SIZE(list))

CONTAINS


RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

implicit none
INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
INTEGER                :: reference, temp
INTEGER, PARAMETER  :: max_simple_sort_size = 6

IF (right_end < left_end + max_simple_sort_size) THEN
  ! Use interchange sort for small lists
  CALL interchange_sort(left_end, right_end)

ELSE
  ! Use partition ("quick") sort
  reference = list((left_end + right_end)/2)
  i = left_end - 1; j = right_end + 1

  DO
    ! Scan list from left end until element >= reference is found
    DO
      i = i + 1
      IF (list(i) >= reference) EXIT
    END DO
    ! Scan list from right end until element <= reference is found
    DO
      j = j - 1
      IF (list(j) <= reference) EXIT
    END DO


    IF (i < j) THEN
      ! Swap two out-of-order elements
      temp = list(i); list(i) = list(j); list(j) = temp
      IF(PRESENT(ORDER)) THEN
        itemp = order(i); order(i) = order(j); order(j) = itemp
      ENDIF
    ELSE IF (i == j) THEN
      i = i + 1
      EXIT
    ELSE
      EXIT
    END IF
  END DO

  IF (left_end < j) CALL quick_sort_1(left_end, j)
  IF (i < right_end) CALL quick_sort_1(i, right_end)
END IF

END SUBROUTINE quick_sort_1



SUBROUTINE interchange_sort(left_end, right_end)

implicit none
INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
INTEGER                :: temp

DO i = left_end, right_end - 1
  DO j = i+1, right_end
    IF (list(i) > list(j)) THEN
      temp = list(i); list(i) = list(j); list(j) = temp
      IF(PRESENT(ORDER)) THEN
        itemp = order(i); order(i) = order(j); order(j) = itemp
      ENDIF
    END IF
  END DO
END DO

END SUBROUTINE interchange_sort

END SUBROUTINE quick_sort_INT


end module quicksort