module cubic_root

!https://blog.csdn.net/he_nan/article/details/78069950
!https://zhuanlan.zhihu.com/p/137077558
!https://fortran-lang.discourse.group/t/numerical-resolution-using-cardan-method/344/5    
implicit none
!Cardano's solutions:
! x1 = s + t - b/(3a) ,  real root
! x2 = -(s+t)/2 - b/(3a) + i(sqrt(3)/2)(s-t) , first complex root
! x3 = -(s+t)/2 - b/(3a) - i(sqrt(3)/2)(s-t) , second complex root
!where:
! s = ( r + sqrt(q**3 + r**2) )**(1/3)
! t = ( r - sqrt(q**3 + r**2) )**(1/3)
! q = (3ac - b**2)/(9a**2)
! r = (9abc - 27da**2 - 2b**3)/(54a**3)
private 
public t_cubic_solution, cubic_solve
real(8),parameter::pi=3.141592653589793d0

type t_cubic_solution
    integer::rtype=1 
    !1: one real roots and two conjugated complex roots
    !3: three real roots
    !2: three real roots,but at least one root is duplicated.
    real(8) :: x1
    complex :: x2
    complex :: x3
end type t_cubic_solution


contains
    pure type(t_cubic_solution) function cubic_solve(a, b, c, d) result(res)
        real(8), intent(in)    :: a, b, c, d
        real(8)                :: s, t, q, r, rpart, ipart, temp,delta,sita1,t1
        
        q       = (3.0*a*c - b**2)/(9.0*a**2) 
        r       = (9.0*a*b*c - 27.0*d*a**2 - 2.0*b**3)/(54.0*a**3)
        delta=q**3+r**2
        if(delta>0.d0) then 
            res.rtype=1 !one real roots
            temp    = r + sqrt(delta)
            s       = sign(1.0, temp) * abs(temp)**(1.0/3.0)

            temp    =  r - sqrt(delta)
            t       = sign(1.0, temp) * abs(temp)**(1.0/3.0)


            res%x1  = s + t - b/(3.0*a)

            rpart   = -(s+t)/2.0 - b/(3.0*a)
            ipart   = (sqrt(3.0)/2.0)*(s-t) 

            res%x2  = cmplx( rpart, ipart )
            res%x3  = cmplx( rpart, -ipart )
        elseif(delta<0) then            
            res.rtype=3 !three real roots
            t1=r/q**2*sqrt(-q)
            sita1=1/3.0d0*acos(t1)
            res.x1=2*sqrt(-q)*cos(sita1)-b/(3.0*a)
            res.x2=cmplx(2*sqrt(-q)*cos(sita1-2*pi/3.)-b/(3.0*a),0.d0)
            res.x3=cmplx(2*sqrt(-q)*cos(sita1-4*pi/3.)-b/(3.0*a),0.d0) 
        else
            res.rtype=2
            t1=(-r)**(1/3.0)
            res.x1=-2*t1
            res.x2=cmplx(t1,0.d0)
            res.x3=res.x2
        endif
    end function cubic_solve
end module cubic_root