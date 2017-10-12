module SolverMath
    
contains    
    
SUBROUTINE invert(matrix)
!
! This subroutine inverts a small square matrix onto itself.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),INTENT(IN OUT)::matrix(:,:)
 REAL(iwp)::det,j11,j12,j13,j21,j22,j23,j31,j32,j33,con
 INTEGER::ndim,i,k
 ndim=UBOUND(matrix,1)

 IF(ndim==2)THEN
   det=matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)
   j11=matrix(1,1)
   matrix(1,1)=matrix(2,2)
   matrix(2,2)=j11
   matrix(1,2)=-matrix(1,2)
   matrix(2,1)=-matrix(2,1)
   matrix=matrix/det
 ELSE IF(ndim==3)THEN
   det=matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3))
   det=det-matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(3,1)*matrix(2,3))
   det=det+matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2))
   j11=matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3)
   j21=-matrix(2,1)*matrix(3,3)+matrix(3,1)*matrix(2,3)
   j31=matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2)
   j12=-matrix(1,2)*matrix(3,3)+matrix(3,2)*matrix(1,3)
   j22=matrix(1,1)*matrix(3,3)-matrix(3,1)*matrix(1,3)
   j32=-matrix(1,1)*matrix(3,2)+matrix(3,1)*matrix(1,2)
   j13=matrix(1,2)*matrix(2,3)-matrix(2,2)*matrix(1,3)
   j23=-matrix(1,1)*matrix(2,3)+matrix(2,1)*matrix(1,3)
   j33=matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2)
   matrix(1,1)=j11
   matrix(1,2)=j12
   matrix(1,3)=j13
   matrix(2,1)=j21
   matrix(2,2)=j22
   matrix(2,3)=j23
   matrix(3,1)=j31
   matrix(3,2)=j32
   matrix(3,3)=j33
   matrix=matrix/det
 ELSE
   DO k=1,ndim
     con=matrix(k,k)
     matrix(k,k)=1.0_iwp
     matrix(k,:)=matrix(k,:)/con
     DO i=1,ndim
       IF(i/=k)THEN
         con=matrix(i,k)
         matrix(i,k)=0.0_iwp
         matrix(i,:)=matrix(i,:)-matrix(k,:)*con
       END IF
     END DO
   END DO
 END IF
RETURN
END SUBROUTINE invert   

FUNCTION determinant(jac) RESULT(det)
!
! This function returns the determinant of a 1x1, 2x2 or 3x3
! Jacobian matrix.
!
 IMPLICIT NONE    
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 REAL(iwp),DIMENSION(:,:),INTENT(IN)::jac
 REAL(iwp)::det
 INTEGER::it 
 
 it=ubound(jac,1)
 !print *,it
 !pause
 
 SELECT CASE(it)
 CASE(1)
   det=1.0_iwp
 CASE(2)
   det=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
 CASE(3)
   det=jac(1,1)*(jac(2,2)*jac(3,3)-jac(3,2)*jac(2,3))
   det=det-jac(1,2)*(jac(2,1)*jac(3,3)-jac(3,1)*jac(2,3))
   det=det+jac(1,3)*(jac(2,1)*jac(3,2)-jac(3,1)*jac(2,2))
 CASE DEFAULT
   WRITE(*,*)' wrong dimension for Jacobian matrix',IT
   PAUSE 
 END SELECT

RETURN
END FUNCTION determinant 

end module