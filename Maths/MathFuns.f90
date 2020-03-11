module SolverMath
IMPLICIT NONE

INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)

contains
    
    
subroutine Nint_gauss(fun,a,b,val,tol,isok)
!ref:http://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature#Fortran
!the weight and x is updated and calculated by the FASTGL.https://people.sc.fsu.edu/~jburkardt/f_src/fastgl/fastgl.html
  implicit none
  integer, parameter :: p = 8 ! quadruple precision
  integer,intent(out) :: isok  ! =0,OK. =1.Not OK.
  real(kind=p),intent(in)::a, b,tol !integration limits, precision
  real(kind=p),intent(out)::val !integrated value.
  real(kind=p),external           :: fun !integrand function
  integer            :: i,n=5,k,miter=50,iter=0
  real(kind=p), allocatable :: r(:,:)
  real(kind=p)       :: last,theta, weight, x
  
  iter=0
  isok=0
  val=0.D0
  n=4
  miter=10
  iter=0
  do while(abs(last-val)>tol.or.iter<2)
    last=val
    
    !r = gaussquad(n)
    
    val=0.d0
    do i=1,n
        call glpair ( n, i, theta, weight, x )
        !val = val+r(2,i)*fun((a+b)/2+r(1,i)*(b-a)/2)
        val = val+weight*fun((a+b)/2+x*(b-a)/2)
    enddo
    val = (b-a)/2*val
    
    n=n*2
    iter=iter+1
    
    if(miter<iter) then
        isok=1
        exit 
    endif
    
  enddo
  
  contains 
 
  function gaussquad(n) result(r)
  integer                 :: n
  real(kind=p), parameter :: pi = 4*atan(1._p)
  real(kind=p)            :: r(2, n), x, f, df, dx
  integer                 :: i,  iter
  real(kind = p), allocatable :: p0(:), p1(:), tmp(:)
 
  p0 = [1._p]
  p1 = [1._p, 0._p]
 
  do k = 2, n
     tmp = ((2*k-1)*[p1,0._p]-(k-1)*[0._p, 0._p,p0])/k
     p0 = p1; p1 = tmp
  end do
  do i = 1, n
    x = cos(pi*(i-0.25_p)/(n+0.5_p))
    do iter = 1, 10
      f = p1(1); df = 0._p
      do k = 2, size(p1)
        df = f + x*df
        f  = p1(k) + x * f
      end do
      dx =  f / df
      x = x - dx
      if (abs(dx)<10*epsilon(dx)) exit
    end do
    r(1,i) = x
    r(2,i) = 2/((1-x**2)*df**2)
  end do
  end function
  
end subroutine    
    
    
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
   !DO k=1,ndim
   !  con=matrix(k,k)
   !  matrix(k,k)=1.0_iwp
   !  matrix(k,:)=matrix(k,:)/con
   !  DO i=1,ndim
   !    IF(i/=k)THEN
   !      con=matrix(i,k)
   !      matrix(i,k)=0.0_iwp
   !      matrix(i,:)=matrix(i,:)-matrix(k,:)*con
   !    END IF
   !  END DO
   !END DO
    call InvByGauss(matrix,ndim)
 END IF
RETURN
END SUBROUTINE invert   


! --------------------------------------------------------------------
SUBROUTINE InvByGauss (a,n)       ! Invert matrix by Gauss method
! --------------------------------------------------------------------
IMPLICIT NONE

INTEGER,intent(in):: n
REAL(8),intent(inout) :: a(n,n)

! - - - Local Variables - - -
REAL(8) :: b(n,n), c, d, temp(n)
INTEGER :: i, j, k, m, imax(1), ipvt(n)
! - - - - - - - - - - - - - -

b = a
ipvt = (/ (i, i = 1, n) /)

DO k = 1,n
   imax = MAXLOC(ABS(b(k:n,k)))
   m = k-1+imax(1)

   IF (m /= k) THEN
      ipvt( (/m,k/) ) = ipvt( (/k,m/) )
      b((/m,k/),:) = b((/k,m/),:)
   END IF
   d = 1/b(k,k)

   temp = b(:,k)
   DO j = 1, n
      c = b(k,j)*d
      b(:,j) = b(:,j)-temp*c
      b(k,j) = c
   END DO
   b(:,k) = temp*(-d)
   b(k,k) = d
END DO

a(:,ipvt) = b

END SUBROUTINE InvByGauss


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
            NNODE=MAX(NINT(T1/DT),1)+1
            ALLOCATE(NODE(N1,NNODE))
            DT1(1:N1)=(X2-X1)/(NNODE-1)
            DO I=1,NNODE
                NODE(:,I)=X1+DT1(1:N1)*(I-1)
            ENDDO
        ENDIF
        
    ENDIF
    
ENDSUBROUTINE



   !计算四面体的其中的一个顶点（第一个）ar(:,1)的立体角，以angle返回。    
	real(8) function solidangle(ar)
	   implicit none
	   
	   real(8),intent(in)::ar(3,4)
       integer::i
       real(8)::br(3,3),t1,v1
	   real(8)::cosA,cosB,cosC,A,B,C,p
	   !br存储由ar形成三个向量，同时化为单位向量
	   do i=1,3
	      br(:,i)=ar(:,i+1)-ar(:,1) 
		  t1=(br(1,i)**2+br(2,i)**2+br(3,i)**2)**0.5
		  if(t1<1e-10) then
		     print *, 'sub solid angle,the distance between two vertex is 0.'
		     stop
		  end if
		  br(:,i)=br(:,i)/t1
       end do
	   
       !求br中第一、二向量所组成的夹角
       cosA=br(1,1)*br(1,2)+br(2,1)*br(2,2)+br(3,1)*br(3,2)
       cosB=br(1,1)*br(1,3)+br(2,1)*br(2,3)+br(3,1)*br(3,3)
	   cosC=br(1,3)*br(1,2)+br(2,3)*br(2,2)+br(3,3)*br(3,2)
	   !因为三个向量都化为了单位向量。
       A=dacos(cosA)
       B=dacos(cosB)
	   C=dacos(cosC)
       p=(A+B+C)/2
	   t1=(sin(p)*sin(p-A)*sin(p-B)*sin(p-C))**0.5/(2*cos(A/2)*cos(B/2)*cos(C/2))
       solidangle=2*dasin(t1)
       
!	   !求第一、二向量所组成三角形的面积
!       s=sin(A)/2
!	   !求四面体的体积:abs(v1)/6
!	   call dt(br,v1)
!	   v1=abs(v1)/6
!	   !求高h
!       h=3*v1/s
!	   !三角形球面的面积s
!	   s=1*h*A/2
       !s=A+B+C-3.1415926536
	   !立体角，球面度
       !angle=s/4/3.1415926536
       !solidangle=s
	end function

    function NORMAL_TRIFACE(V) result (Normal)

    !*****************************************************************************80
    !V, XY OF THE FACET.
    !CALCULATE THE NORMAL VECTOR OF A TRI-FACET. 
    !

    !
      implicit none

      REAL(8),INTENT(IN)::V(:,:) !3*3

      real ( kind = 8 ) v1(SIZE(V,DIM=1))
      real ( kind = 8 ) v2(SIZE(V,DIM=1))
      real ( kind = 8 ) normal(SIZE(V,DIM=1))
      
      V1=V(:,2)-V(:,1);V2=V(:,3)-V(:,1);

      normal(1) = v1(2) * v2(3) - v1(3) * v2(2)
      normal(2) = v1(3) * v2(1) - v1(1) * v2(3)
      normal(3) = v1(1) * v2(2) - v1(2) * v2(1)
      !normal=normal/norm2(normal)
      return
  
    end function
    
    !二面角，单位弧度,N1,N2为面的方向矢量
    real(8) function DihedralAngle(N1,N2) 
        implicit none
        integer ( kind = 4 ), parameter :: dim_num = 3
        real ( kind = 8 ),INTENT(IN):: N1(dim_num)
        real ( kind = 8 ),INTENT(IN):: N2(dim_num)
        REAL(8),PARAMETER::PI=3.141592653589793
        
        DihedralAngle=PI-DACOS(DOT_PRODUCT(N1,N2)/(NORM2(N1)*NORM2(N2)))
        
        
    
    end function

end module