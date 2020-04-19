!Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracy
!and  adjust  stepsize.  Input  are  the  dependent  variable  vector y(1:n) and  its  derivative
!dydx(1:n) at the starting value of the independent variable x.  Also input are the stepsize
!to be attempted htry, the required accuracy eps, and the vector yscal(1:n) against
!which the error is scaled.  On output, y and x are replaced by their new values, hdid is the
!stepsize that was actually accomplished, and hnext is the estimated next stepsize. derivs
!is the user-supplied subroutine that computes the right-hand side derivatives.
      
      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      USE ODE_SOLVER
      INTEGER n,NMAX,BCITER1
      REAL eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
CU    USES derivs,rkck
      INTEGER i
      REAL errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,
     *PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)
      h=htry;BCITER1=1
1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs) !x is time
      !IF(RKINFO.ISOUTOFRANGE) THEN
      !    BCITER1=BCITER1+1
      !    H=(RKINFO.HF_SUC+RKINFO.HF_FAIL)/2.0*HTRY
      !    IF(BCITER1<=10.AND.H>0.01*HTRY) THEN
      !        GOTO 1
      !    ELSE
      !        RETURN
      !    ENDIF
      !ELSE
      !    BCITER1=1
      !ENDIF
      errmax=0.
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
11    continue
      errmax=errmax/eps
      if(errmax.gt.1.)then
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(abs(htemp),0.1*abs(h)),h)
        xnew=x+h
        if(xnew.eq.x) pause 'stepsize underflow in rkqs'
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      END
