      SUBROUTINE rk4(y,dydx,n,x,h,yout,derivs)
      INTEGER n,NMAX
      REAL h,x,dydx(n),y(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
      INTEGER i
      REAL h6,hh,xh,dym(NMAX),dyt(NMAX),yt(NMAX)
      hh=h*0.5
      h6=h/6.
      xh=x+hh
      do 11 i=1,n
        yt(i)=y(i)+hh*dydx(i)
11    continue
      call derivs(xh,yt,dyt)
      do 12 i=1,n
        yt(i)=y(i)+hh*dyt(i)
12    continue
      call derivs(xh,yt,dym)
      do 13 i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
13    continue
      call derivs(x+h,yt,dyt)
      do 14 i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
14    continue
      return
      END
