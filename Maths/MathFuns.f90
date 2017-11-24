module SolverMath
IMPLICIT NONE

INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)

contains    
    
SUBROUTINE invert(matrix)
!
! This subroutine inverts a small square matrix onto itself.
!
 IMPLICIT NONE
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


SUBROUTINE EQUDIVIDE(X1,X2,DT,NODE,NNODE) 
!IF ABS(X1-X2)>0 AND DT>0,DIVIDE SEGMENT(X1,X2) INTO MAX(INT((X2-X1)/DT),1) SUBSEGS AND RETURN ALL NODES
!IF ABS(X1-X2)>0 AND DT=0, NODE=[X1,X2]
!IF  ABS(X1-X2)=0, RETURN NODE=[X1]
    IMPLICIT NONE
    REAL(IWP),INTENT(IN)::X1(:),X2(:),DT
    INTEGER,INTENT(OUT)::NNODE
    REAL(IWP),ALLOCATABLE::NODE(:,:)
    INTEGER::N1,I
    REAL(IWP)::T1,DT1(3)
    
    T1=NORM2(X1-X2)
    N1=SIZE(X1)
    IF(ALLOCATED(NODE)) DEALLOCATE(NODE)
    IF(T1<1.D-10) THEN
        NNODE=1
        ALLOCATE(NODE(N1,1))
        NODE(:,1)=X1
    ELSE
        IF(DT<1.D-10) THEN
            NNODE=2
            ALLOCATE(NODE(N1,2))
            NODE(:,1)=X1;NODE(:,2)=X2;        
        ELSE
            NNODE=MAX(INT(T1/DT),1)+1
            ALLOCATE(NODE(N1,NNODE))
            DT1(1:N1)=(X2-X1)/(NNODE-1)
            DO I=1,NNODE
                NODE(:,I)=X1+DT1*(I-1)
            ENDDO
        ENDIF
        
    ENDIF
    
ENDSUBROUTINE

end module