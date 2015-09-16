! Module implementing an OO hash table (dictionary) in Fortran 2003.
! Compiles and runs with accompanying test program under the Intel 
! Fortran Compiler, version 11.1.046

! Copyright (c) Izaak Beekman 2010

    ! This program is free software: you can redistribute it and/or modify
    ! it under the terms of the GNU Lesser General Public License as published by
    ! the Free Software Foundation, either version 3 of the License, or
    ! (at your option) any later version.

    ! This program is distributed in the hope that it will be useful,
    ! but WITHOUT ANY WARRANTY; without even the implied warranty of
    ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    ! GNU Lesser General Public License for more details.

    ! You should have received a copy of the GNU Lesser General Public License
    ! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE DS_MYDATA
	IMPLICIT NONE
	TYPE MYDATA
		INTEGER::IEL=0 !ELEMENT NO
		INTEGER::ISE=0 !THE SUB IDEX OF THE COMMON EDGE/FACE
	ENDTYPE
END MODULE

MODULE hashtbl
  USE DS_MYDATA, DICT_DATA=>MYDATA
  IMPLICIT NONE ! Use strong typing
  INTEGER:: tbl_size = 10000,NLISTVAL=10
	

  
  TYPE sllist
     TYPE(sllist), POINTER :: child => NULL()
     CHARACTER(len=:), ALLOCATABLE :: key 
     integer::nval=0
     TYPE(DICT_DATA),allocatable::val(:)
   CONTAINS
     PROCEDURE :: put  => put_sll
     PROCEDURE :: get  => get_sll
     PROCEDURE :: free => free_sll
     PROCEDURE,NOPASS :: ENLARGEVAL => ENLARGE_VAL 
  END TYPE sllist

  TYPE hash_tbl_sll
     TYPE(sllist), DIMENSION(:), ALLOCATABLE :: vec
     INTEGER                                 :: vec_len = 0
     LOGICAL                                 :: is_init = .FALSE.
  CONTAINS
     PROCEDURE::KEY=>I2C_KEY_hash_tbl_sll
     PROCEDURE :: init => init_hash_tbl_sll
     PROCEDURE :: put  => put_hash_tbl_sll
     PROCEDURE :: get  => get_hash_tbl_sll
     PROCEDURE :: free => free_hash_tbl_sll
  END TYPE hash_tbl_sll

  PUBLIC :: hash_tbl_sll
  
  TYPE(hash_tbl_sll)::EDGE_TBL,FACE_TBL
  
  CONTAINS

  RECURSIVE SUBROUTINE put_sll(list,key,val)
    CLASS(sllist),    INTENT(inout) :: list
    CHARACTER(len=*), INTENT(in)    :: key  !, val
    INTEGER                         :: keylen !, vallen
	TYPE(DICT_DATA)::VAL
	
    keylen = LEN(key)
    !vallen = LEN(val)
    IF (ALLOCATED(list%key)) THEN
       IF (list%key /= key) THEN
          IF ( .NOT. ASSOCIATED(list%child) ) ALLOCATE(list%child)
          CALL put_sll(list%child,key,val)
       ELSE
          LIST.NVAL=LIST.NVAL+1
          IF(LIST.NVAL>NLISTVAL) CALL LIST.ENLARGEVAL(LIST.VAL)
          list%val(LIST.NVAL) = val
       END IF
    ELSE
       IF (.NOT. ALLOCATED(list%key)) &
            ALLOCATE(CHARACTER(len=keylen) :: list%key)
       list%key = key
       IF (ALLOCATED(list%val)) DEALLOCATE(list%val)
       !ALLOCATE(CHARACTER(len=vallen) :: list%val)       
       ALLOCATE(LIST.VAL(NLISTVAL))       
       LIST.NVAL=1
       !LIST.VAL(0)=LIST.VAL(0)+1
       list%val(LIST.NVAL) = val
    END IF
  END SUBROUTINE put_sll


  RECURSIVE SUBROUTINE get_sll(list,key,val)
    CLASS(sllist),                 INTENT(in)    :: list
    CHARACTER(len=*),              INTENT(in)    :: key
    TYPE(DICT_DATA), ALLOCATABLE, INTENT(out)   :: val(:)
    !INTEGER                                      :: vallen

    !vallen = 0
    IF (ALLOCATED(list%key) .AND. (list%key == key)) THEN
       !vallen = SIZE(list%val,DIM=1)
       IF (ALLOCATED(val)) DEALLOCATE(val)
       ALLOCATE(val,SOURCE=LIST.VAL(1:LIST.NVAL))
       !val = list%val
    ELSE IF(ASSOCIATED(list%child)) THEN ! keep going
       CALL get_sll(list%child,key,val)
    ELSE ! At the end of the list, no key found
       IF (ALLOCATED(val)) DEALLOCATE(val) ! Exit indication
       RETURN
    END IF
  END SUBROUTINE get_sll


  RECURSIVE SUBROUTINE free_sll(list)
    CLASS(sllist), INTENT(inout) :: list
    IF (ASSOCIATED(list%child)) THEN
       CALL free_sll(list%child)
       DEALLOCATE(list%child)
    END IF
    list%child => NULL()
    IF (ALLOCATED(list%key)) DEALLOCATE(list%key)
    IF (ALLOCATED(list%val)) DEALLOCATE(list%val)
  END SUBROUTINE free_sll
  
  SUBROUTINE ENLARGE_VAL(AVAL)
    TYPE(DICT_DATA),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
    TYPE(DICT_DATA),ALLOCATABLE::VAL1(:)
    INTEGER::LB1=0,UB1=0
    
    LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
    ALLOCATE(VAL1,SOURCE=AVAL)
    DEALLOCATE(AVAL)
    ALLOCATE(AVAL(LB1:UB1+10))
    AVAL(LB1:UB1)=VAL1
    !AVAL(UB1+1:UB1+10)=0
    DEALLOCATE(VAL1)
  END SUBROUTINE
  
  SUBROUTINE init_hash_tbl_sll(tbl,tbl_len)
    CLASS(hash_tbl_sll),   INTENT(inout) :: tbl
    INTEGER,     OPTIONAL, INTENT(in)    :: tbl_len

    IF (ALLOCATED(tbl%vec)) DEALLOCATE(tbl%vec)
    IF (PRESENT(tbl_len)) THEN
       ALLOCATE(tbl%vec(0:tbl_len-1))
       tbl%vec_len = tbl_len
    ELSE
       ALLOCATE(tbl%vec(0:tbl_size-1))
       tbl%vec_len = tbl_size
    END IF
    tbl%is_init = .TRUE.
  END SUBROUTINE init_hash_tbl_sll
  
  SUBROUTINE I2C_KEY_hash_tbl_sll(tbl,KEY,AA,NAA)
    CLASS(hash_tbl_sll),   INTENT(inout) :: tbl
    INTEGER,INTENT(in)    ::NAA,AA(NAA)
    CHARACTER(LEN=:),ALLOCATABLE,INTENT(OUT)::KEY
    INTEGER::I,J,N1,AB1(10)
    CHARACTER(64)::KEY1,KEY2
    
    AB1(1:NAA)=AA(1:NAA)
    DO I=1,NAA
        DO J=I+1,NAA
            IF(AB1(J)<AB1(I)) THEN
                N1=AB1(I);AB1(I)=AB1(J);AB1(J)=N1
            ENDIF
        ENDDO
        WRITE(KEY1,'(I8)') AB1(I)
        KEY2=TRIM(KEY2)//'+'//TRIM(ADJUSTL(KEY1))
    ENDDO
    ALLOCATE(CHARACTER(LEN=LEN_TRIM(KEY2))::KEY) 
    KEY=TRIM(KEY2)

  END SUBROUTINE I2C_KEY_hash_tbl_sll  
 
  SUBROUTINE put_hash_tbl_sll(tbl,key,val)
    CLASS(hash_tbl_sll), INTENT(inout) :: tbl
    CHARACTER(len=*),    INTENT(in)    :: key
    INTEGER                            :: hash
    TYPE(DICT_DATA)::VAL
     
    hash = MOD(ABS(HASH_DJB(KEY)),tbl%vec_len)
    CALL tbl%vec(hash)%put(key=key,val=val)
  END SUBROUTINE put_hash_tbl_sll


  SUBROUTINE get_hash_tbl_sll(tbl,key,val)
    CLASS(hash_tbl_sll),           INTENT(in)    :: tbl
    CHARACTER(len=*),              INTENT(in)    :: key
    TYPE(DICT_DATA), ALLOCATABLE, INTENT(out)   :: val(:)
    INTEGER                                      :: hash

    hash = MOD(ABS(HASH_DJB(KEY)),tbl%vec_len)
    CALL tbl%vec(hash)%get(key=key,val=val)
  END SUBROUTINE get_hash_tbl_sll


  SUBROUTINE free_hash_tbl_sll(tbl)
    CLASS(hash_tbl_sll), INTENT(inout) :: tbl    
    INTEGER     :: i, low, high

    low  = LBOUND(tbl%vec,dim=1)
    high = UBOUND(tbl%vec,dim=1) 
    IF (ALLOCATED(tbl%vec)) THEN
       DO i=low,high
          CALL tbl%vec(i)%free()
       END DO
       DEALLOCATE(tbl%vec)
    END IF
    tbl%is_init = .FALSE.
  END SUBROUTINE free_hash_tbl_sll

INTEGER FUNCTION HASH_DJB(KEY)
    CHARACTER(LEN=*), INTENT(IN) :: KEY
    INTEGER :: I, KEYLEN, UV
    !INTEGER (KIND=4) :: HASHVAL

    HASH_DJB = 5381
    KEYLEN = LEN_TRIM(KEY)
    DO I = 1, KEYLEN
        UV = MOD(HASH_DJB * 33, 65536_4)
        HASH_DJB = MOD(UV + IACHAR(KEY(I:I)), 65536_4)
    END DO
    !HASH = MOD (HASHVAL, HASHSIZE) + 1
END FUNCTION


END MODULE hashtbl
    
    
    
