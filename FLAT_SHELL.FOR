C=DECK SM3SHELL SM3SHELL FORTRAN
C=PURPOSE Form  stiffness of 3-node, 18-dof flat triangle shell
C=AUTHOR C. A. Felippa, October 1996
C=VERSION October 1996
C=EQUIPMENT Machine independent
C=KEYWORDS Kirchhoff thin shell isotropic homogeneous
C=KEYWORDS finite element triangle stiffness matrix
C=BLOCK ABSTRACT
C 
C     SM3SHELL forms the material element stiffness matrix
C     18-dof flat triangular shell element constructed with a
C     comination of two high performance components.
C     Isotropic material, homogeneous wall construction.
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C     call SM3SHELL (x, y, z, em, nu, h, sm, m, status)
C
C     The input arguments are:
C
C     X       (3 x 1) array of x coordinates of triangle nodes
C     Y       (3 x 1) array of y coordinates of triangle nodes
C     Z       (3 x 1) array of z coordinates of triangle nodes
C     EM      Elastic modulus
C     NU      Poisson's ratio
C     H       Thicknes of element at nodes (3 x 1) array
C     M       First dimension of SM in calling program.
C
C     The outputs are:
C
C     SM      Output stiffness array with basic stiffness
C             coefficients added in.  the (I,J)-th entry of the
C             (9 x 9) element bending stiffness is added to
C             SM(K,L), where K=LS(I) and L=LS(J).
C     STATUS  Status character variable. Blank if no error detected.
C		AREA		TRIANGLE ELEMENT AREA
C=END USAGE
      subroutine SM3SHELL
     &      (xg, yg, zg, h, em, nu, esm, m, dcm, status,area)
	implicit none
C
C                   A R G U M E N T S
C
      character*(*),intent(out):: status
      integer,intent(in)::m
      double precision,intent(in)::xg(3), yg(3), zg(3) 
      double precision,intent(in)::em, nu, h
      double precision,intent(out)::esm(m,m),dcm(3,3),area
C
C                   L O C A L   V A R I A B L E S
C
      double precision  sm(18,18), db(3,3), dm(3,3)
      double precision  xlp(3), ylp(3), zlp(3) 
      double precision  h0, clr, cqr, betab
      integer           lb(9), le(9)
      integer           i, j, nd
C
C                   D A T A
C
      data       lb    /3,4,5, 9,10,11, 15,16,17/
      data       le    /1,2, 7,8, 13,14, 6,12,18/
C
C                   L O G I C
C
      status =  ' '
C
C               Establish local system
C
      call SM3SHLOCALSYS (xg,yg,zg, xlp,ylp,zlp, dcm,area, status)
      if (status(1:1) .ne. ' ')       return
      if (area .le. 0.0)      then
        status = 'SM3SHELLST: Negative area'
        if (area .eq. 0.0)   status = 'SM3SHELLST: Zero area'
        return
      end if
C
C               Form constitutive matrices
C
C      h0 =    (h(1) + h(2) + h(3))/3.d0
	h0=h
      call   SM3SHISODM (em, nu, h0, dm)
      call   SM3SHISODB (em, nu, h0, db)
C
C               Form local stiffnesses
C
      do 1800  i = 1,18
        do 1800  j = 1,18
          sm(i,j)  = 0.0
 1800   continue
c     betab =  max(0.01,0.5*(1.-4.*nu**2))
      betab =  0.32d+00
      clr =    0.0d0
      cqr =    1.0d0
      nd =     18
      call   SM3SHMEMBB 
     &       (xlp,ylp,dm,1.5d0,1.0d0,le,sm,nd,status)
      call   SM3SHMEMBH (xlp,ylp,dm,betab,le,sm,nd,status)
      call   SM3SHBENDB (xlp,ylp,db,1.0d0,clr,cqr,lb,sm,nd,status)
      call   SM3SHBENDH (xlp,ylp,db,1.0d+0,lb,sm,nd,status)

C
C              Transform to global coordinates
C
      call   SM3SHTRANSFORM (sm, dcm,dcm,dcm,dcm,dcm,dcm)
      do 3000  i = 1,18
        do 2500  j = 1,18
          esm(i,j)  = sm(i,j)
 2500   continue
 3000 continue
      return
      end
C=END FORTRAN
C=DECK SM3SHLOCALSYTEM
C=PURPOSE Define local system of 3-node triangle in 3D space
C=AUTHOR C. A. Felippa
C=VERSION September 1996
C=EQUIPMENT Machine independent
C=KEYWORDS finite element shell triangle local direction cosines
C=BLOCK ABSTRACT
C
C     SM3SHLOCALSYS computes the local corner coordinates
C     of a flat triangular element in 3D space and the direction
C     cosines of the local system
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       call  SM3SHLOCALSYS (xg,yg,zg, xl,yl,zl, dcm, area,status)
C
C     where the input arguments are:
C
C       XG,YG,ZG  Corner coordinates of triangle in global system
C
C     The outputs are:
C
C       XL,YL,ZL   Corner coordinates of triangle in local system 
C       DCM        Matrix of direction cosines of local system
C       AREA       Signed area of triangle
C       STATUS     Blank if no error detected, else error message
C
C     The local system is defined as follows:
C        x' is directed parallel to the 2-1 side
C        z' is the external normal (counterclockwise).
C        y' computed as z' x x'
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine  SM3SHLOCALSYS 
     &            (xg,yg,zg, xl,yl,zl, dcm, area, status)
C
C                   A R G U M E N T S
C
      character*(*)     status
      double precision  xg(3), yg(3), zg(3)
      double precision  xl(3), yl(3), zl(3)
      double precision  dcm(3,3), area
C
C                   L O C A L   V A R I A B L E S
C
      double precision  x21, y21, z21, x32, y32, z32
      double precision  xlr, ylr, zlr, x0, y0, z0
      double precision  dx(3), dy(3), dz(3)
      integer           i
C
C                   L O G I C
C
      status = ' '
      x21 =   xg(2) - xg(1)
      y21 =   yg(2) - yg(1)
      z21 =   zg(2) - zg(1)
      x32 =   xg(3) - xg(2)
      y32 =   yg(3) - yg(2)
      z32 =   zg(3) - zg(2)
      xlr =   sqrt( x21**2+y21**2+z21**2 )
      if (xlr .eq. 0.0)               then
         status =  'SM3SHLOCALSYTEM: nodes 1-2 coincide'
         return
      end if
      dx(1) =  x21/xlr
      dx(2) =  y21/xlr
      dx(3) =  z21/xlr
      dz(1) =  y21*z32 - z21*y32
      dz(2) =  z21*x32 - x21*z32
      dz(3) =  x21*y32 - y21*x32
      zlr  =   sqrt( dz(1)**2 + dz(2)**2 + dz(3)**2 )
      if (zlr .eq. 0.0)               then
         status =  'SM3SHLOCALSYTEM: nodes 1-2-3 are colinear'
         return
      end if
      dz(1) =  dz(1)/zlr
      dz(2) =  dz(2)/zlr
      dz(3) =  dz(3)/zlr
      dy(1) =  dz(2) * dx(3) - dz(3) * dx(2)
      dy(2) =  dz(3) * dx(1) - dz(1) * dx(3)
      dy(3) =  dz(1) * dx(2) - dz(2) * dx(1)
      ylr  =   sqrt( dy(1)**2 + dy(2)**2 + dy(3)**2 )
      dy(1) =  dy(1)/ylr
      dy(2) =  dy(2)/ylr
      dy(3) =  dy(3)/ylr
      x0 =     (xg(1) + xg(2) + xg(3))/3.0d0
      y0 =     (yg(1) + yg(2) + yg(3))/3.0d0
      z0 =     (zg(1) + zg(2) + zg(3))/3.0d0 
      do 2000  i = 1,3
        xl(i) = dx(1)*(xg(i)-x0) + dx(2)*(yg(i)-y0) + dx(3)*(zg(i)-z0)
        yl(i) = dy(1)*(xg(i)-x0) + dy(2)*(yg(i)-y0) + dy(3)*(zg(i)-z0)
        zl(i) = dz(1)*(xg(i)-x0) + dz(2)*(yg(i)-y0) + dz(3)*(zg(i)-z0)
        dcm(1,i) =  dx(i)
        dcm(2,i) =  dy(i)
        dcm(3,i) =  dz(i)
 2000   continue
      area =  0.5d0*( (xl(2)-xl(1))*(yl(3)-yl(1)) 
     &              - (xl(3)-xl(1))*(yl(1)-yl(2)) )
      return
      end
C=END FORTRAN
C=DECK SM3SHISODB SM3SHISODB FORTRAN
C=PURPOSE Form isotropic bending constitutive matrix DB
C=AUTHOR C. A. Felippa
C=BLOCK FORTRAN
      subroutine  SM3SHISODB (em, nu, h, db)
      double precision       em, nu, h, db(3,3)
      double precision       c
      c =          em*h**3/(12.d0*(1.-nu**2))
      db(1,1) =    c
      db(1,2) =    nu*c
      db(2,1) =    nu*c
      db(2,2) =    c
      db(3,3) =    0.5d0*(1.-nu)*c
      db(1,3) =    0.0
      db(2,3) =    0.0
      db(3,1) =    0.0
      db(3,2) =    0.0
      return
      end
C=END FORTRAN
C=DECK SM3SHISODM SM3SHISODM FORTRAN
C=PURPOSE Form isotropic membrane constitutive matrix DM
C=AUTHOR C. A. Felippa
C=BLOCK FORTRAN
      subroutine  SM3SHISODM (em, nu, h, dm)
      double precision       em, nu, h, dm(3,3)
      double precision       c
      c =          em*h/(1.-nu**2)
      dm(1,1) =    c
      dm(1,2) =    nu*c
      dm(2,1) =    nu*c
      dm(2,2) =    c
      dm(3,3) =    0.5d0*(1.-nu)*c
      dm(1,3) =    0.0
      dm(2,3) =    0.0
      dm(3,1) =    0.0
      dm(3,2) =    0.0
      return
      end
C=END FORTRAN
C=DECK SM3SHTRANSFORM
C=PURPOSE Transform stiffness of 3-node shell element to global axes
C=AUTHOR C. A. Felippa
C=VERSION September 1996
C=EQUIPMENT Machine independent
C=KEYWORDS finite element shell triangle local global transform
C=BLOCK ABSTRACT
C
C     SM3SHTRANSFORM transforms the local element stiffness to
C     global coordinates. 
C
C    The stiffness transformation is assumed to have the block 
C    diagonal form
C
C [T1' 0  0  0  0  0 ] [S11 S12 S13 S14 S15 S16] [T1  0  0  0  0  0]
C [ 0 T2' 0  0  0  0 ] [    S22 S23 S24 S25 S26] [ 0 T2  0  0  0  0]
C [ 0  0 T3' 0  0  0 ] [        S33 S34 S35 S36] [ 0  0 T3  0  0  0]
C [ 0  0  0 T4' 0  0 ] [            S44 S45 S46] [ 0  0  0 T4  0  0] 
C [ 0  0  0  0 T5' 0 ] [                S55 S56] [ 0  0  0  0 T5  0]
C [ 0  0  0  0  0 T6'] [   symm             S66] [ 0  0  0  0  0 T6]
C
C where Ti are direction cosines of local directions wrt global.
C Each block is 3 x 3. Indices 1, 3 and 5 pertains to the three 
C translations at corners 1, 2 and 3, respectively. Indices 2, 4
C and 6 pertain to the three rotations at corners 1, 2 and 3.
C    
C The subroutine accounts for the possibility that T1, T2, ... T6 
C may be different, although they will often be the same.
C
C The transformed block Sij is evidently Ti'.Sij.Tj, which can be
C computed in 54 multiply-add operations.  For all blocks this
C amounts to about 1000 operations. The implementation below uses
C 50% more operations: 1539, to streamline loops.  A brute force   
C version ignoring the block-diagonal form would require 2 x18^3/2 
C ~ 6000  operations.  Tests on the Sun show that the present
C implementation is 8 times faster than the brute force, because
C of the careful use of local variables to reduce indexing overhead. 
C     
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       call  SM3SHTRANSFORM (sm, T1,T2,T3,T4,T5,T6)
C
C     where the input arguments are:
C
C      SM               Incoming 18 x 18 matrix in local system
C      T1 through T6    Direction cosine matrices - see abstract
C     
C
C     The outputs are:
C
C      SM               Output 18 x 18 matrix in global system
C     
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine  SM3SHTRANSFORM 
     &            (sm, t1, t2, t3, t4, t5, t6)
C
C                   A R G U M E N T S
C
      double precision  t1(3,3), t2(3,3), t3(3,3)
      double precision  t4(3,3), t5(3,3), t6(3,3)
      double precision  sm(18,18)
C
C                   L O C A L   V A R I A B L E S
C
      double precision  t(3,3,6), st1(18), st2(18), st3(18)
      double precision  t11, t12, t13, t21, t22, t23, t31, t32, t33
      double precision  tst11, tst12, tst13, tst21, tst22, tst23  
      double precision  tst31, tst32, tst33 
      integer           i, j, ib, jb
C
C                   L O G I C
C
      do 1400  j = 1,3
        do 1200  i = 1,3
          t(i,j,1) = t1(i,j)
          t(i,j,2) = t2(i,j)
          t(i,j,3) = t3(i,j)
          t(i,j,4) = t4(i,j)
          t(i,j,5) = t5(i,j)
          t(i,j,6) = t6(i,j)
 1200   continue
 1400 continue
      do 3000  jb = 1,6
        j =     3*(jb-1)
        t11 =   t(1,1,jb)
        t21 =   t(2,1,jb)
        t31 =   t(3,1,jb)
        t12 =   t(1,2,jb)
        t22 =   t(2,2,jb)
        t32 =   t(3,2,jb)
        t13 =   t(1,3,jb)
        t23 =   t(2,3,jb)
        t33 =   t(3,3,jb)
        do 2200  i = 1,18
          st1(i) =  sm(i,j+1)*t11 + sm(i,j+2)*t21 + sm(i,j+3)*t31
          st2(i) =  sm(i,j+1)*t12 + sm(i,j+2)*t22 + sm(i,j+3)*t32  
          st3(i) =  sm(i,j+1)*t13 + sm(i,j+2)*t23 + sm(i,j+3)*t33
 2200   continue
        do 2800  ib = 1,jb
          i =     3*(ib-1)
          t11 =   t(1,1,ib)
          t21 =   t(2,1,ib)
          t31 =   t(3,1,ib)
          t12 =   t(1,2,ib)
          t22 =   t(2,2,ib)
          t32 =   t(3,2,ib)
          t13 =   t(1,3,ib)
          t23 =   t(2,3,ib)
          t33 =   t(3,3,ib)
          tst11 = t11*st1(i+1) + t21*st1(i+2) + t31*st1(i+3)
          tst21 = t12*st1(i+1) + t22*st1(i+2) + t32*st1(i+3)  
          tst31 = t13*st1(i+1) + t23*st1(i+2) + t33*st1(i+3)
          tst12 = t11*st2(i+1) + t21*st2(i+2) + t31*st2(i+3)
          tst22 = t12*st2(i+1) + t22*st2(i+2) + t32*st2(i+3)  
          tst32 = t13*st2(i+1) + t23*st2(i+2) + t33*st2(i+3)
          tst13 = t11*st3(i+1) + t21*st3(i+2) + t31*st3(i+3)
          tst23 = t12*st3(i+1) + t22*st3(i+2) + t32*st3(i+3)  
          tst33 = t13*st3(i+1) + t23*st3(i+2) + t33*st3(i+3)
          sm(i+1,j+1) =   tst11
          sm(j+1,i+1) =   tst11
          sm(i+1,j+2) =   tst12
          sm(j+2,i+1) =   tst12
          sm(i+1,j+3) =   tst13
          sm(j+3,i+1) =   tst13
          sm(i+2,j+1) =   tst21
          sm(j+1,i+2) =   tst21
          sm(i+2,j+2) =   tst22
          sm(j+2,i+2) =   tst22
          sm(i+2,j+3) =   tst23
          sm(j+3,i+2) =   tst23
          sm(i+3,j+1) =   tst31
          sm(j+1,i+3) =   tst31
          sm(i+3,j+2) =   tst32
          sm(j+2,i+3) =   tst32
          sm(i+3,j+3) =   tst33
          sm(j+3,i+3) =   tst33
 2800   continue
 3000 continue
      return
      end
C=END FORTRAN
C=DECK SM3SHBENDB SM3SHBENDB FORTRAN
C=PURPOSE Form bending stiffness of 3-node, 9dof Kirchhoff triangle
C=AUTHOR C. A. Felippa, May 1984
C=VERSION September 1989
C=EQUIPMENT Machine independent
C=KEYWORDS Kirchhoff thin plate bending
C=KEYWORDS finite element triangle basic stiffness matrix
C=BLOCK ABSTRACT
C
C     SM3SHBENDB forms the basic material element stiffness matrix
C     of a Kirchhoff thin plate bending triangle constructed with the 
C     generalized FF/ANS formulation
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C     call  SM3SHBENDB (x,y,db, f,clr,cqr, ls,sm,m,status)
C
C     where the input arguments are:
C
C     X        (3 x 1) array of x coordinates of triangle nodes
C     Y        (3 x 1) array of y coordinates of triangle nodes
C     DB       (3 x 3) moment-curvature constitutive matrix
C     F        Factor by which stiffness entries will be multiplied.
C     CLR,CQR  Use CLR*LLR+CQR*LQR for lumping matrix L
C              thus CLR+CQR  must add up to 1.
C              LLR=linear rotation, LLR=quadratic rotation (Kirchhoff)
C     LS       (9 x 1) array of stiffness location pointers
C              (see output SM).
C     SM       Incoming material stiffness array.
C     M        First dimension of SM in calling program.
C
C     The outputs are:
C
C     SM       Output stiffness array with basic stiffness
C              coefficients added in.  the (I,J)-th entry of the
C              (9 x 9) element bending stiffness is added to
C              SM(K,L), where K=LS(I) and L=LS(J).
C     STATUS   Status character variable.  blank if no error detected.
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    SM3SHBENDB
     &             ( x, y, db, f, clr, cqr, ls, sm, m, status)
C
C                   A R G U M E N T S
C
      double precision   x(3), y(3), db(3,3), f, sm(m,m)
      double precision   clr, cqr
      integer            m, ls(9)
      character*(*)      status
C
C                   L O C A L   V A R I A B L E S
C
      double precision   llr(9,3), lqr(9,3), l(9,3)
      double precision   db11, db12, db13, db22, db23, db33
      double precision   x0, y0, cab, a1, a2, a3, b1, b2, b3
      double precision   x21, x32, x13, y21, y32, y13
      double precision   x12, x23, x31, y12, y23, y31
      double precision   xl12, xl23, xl31
      double precision   c12, c23, c31, s12, s23, s31
      double precision   cc12, cc23, cc31, ss12, ss23, ss31
      double precision   cs12, cs23, cs31, s1, s2, s3
      double precision   area2, c
      integer            i, j, ii, jj
C
C                   L O G I C
C
      status =   ' '
      x21 =      x(2) - x(1)
      x12 =     -x21
      x32 =      x(3) - x(2)
      x23 =     -x32
      x13 =      x(1) - x(3)
      x31 =     -x13
      y21 =      y(2) - y(1)
      y12 =     -y21
      y32 =      y(3) - y(2)
      y23 =     -y32
      y13 =      y(1) - y(3)
      y31 =     -y13
      area2 =    y21*x13 - x21*y13
      if (area2 .le. 0.0)      then
        status = 'SM3SHBENDB: Negative area'
        if (area2 .eq. 0.0)   status = 'SM3SHBENDB: Zero area'
        return
      end if
      x0 =      (x(1)+x(2)+x(3))/3.0d0
      y0 =      (y(1)+y(2)+y(3))/3.0d0
      cab =      3.0d0/area2
      a1 =      -cab*(y(3)-y0)
      a2 =      -cab*(y(1)-y0)
      a3 =      -cab*(y(2)-y0)
      b1 =       cab*(x(3)-x0)
      b2 =       cab*(x(1)-x0)
      b3 =       cab*(x(2)-x0)
      xl12 =     sqrt(x12**2+y12**2)
      xl23 =     sqrt(x23**2+y23**2)
      xl31 =     sqrt(x31**2+y31**2)
      do 1200  j = 1,3
        do 1200  i = 1,9
          llr(i,j) =  0.0
          lqr(i,j) =  0.0
 1200   continue
      if (clr .ne. 0.0)         then
        llr(3,1) =   y32*0.5d0
        llr(6,1) =   y13*0.5d0
        llr(9,1) =   y21*0.5d0
        llr(2,2) =   x32*0.5d0
        llr(5,2) =   x13*0.5d0
        llr(8,2) =   x21*0.5d0
        llr(2,3) =  -y32*0.5d0
        llr(3,3) =  -x32*0.5d0
        llr(5,3) =  -y13*0.5d0
        llr(6,3) =  -x13*0.5d0
        llr(8,3) =  -y21*0.5d0
        llr(9,3) =  -x21*0.5d0
      end if
c
      if (cqr .ne. 0.0)         then
        c12 =       y21/xl12
        s12 =       x12/xl12
        c23 =       y32/xl23
        s23 =       x23/xl23
        c31 =       y13/xl31
        s31 =       x31/xl31
        cc12 =      c12*c12
        cc23 =      c23*c23
        cc31 =      c31*c31
        ss12 =      s12*s12
        ss23 =      s23*s23
        ss31 =      s31*s31
        cs12 =      c12*s12
        cs23 =      c23*s23
        cs31 =      c31*s31
        lqr(1,1) =   cs12 - cs31
        lqr(1,2) =  -lqr(1,1)
        lqr(1,3) =  (cc31-ss31) - (cc12-ss12)
        lqr(2,1) =  (cc12*x12 + cc31*x31)*0.5d0
        lqr(2,2) =  (ss12*x12 + ss31*x31)*0.5d0
        lqr(2,3) =   ss12*y21 + ss31*y13
        lqr(3,1) = -(cc12*y21 + cc31*y13)*0.5d0
        lqr(3,2) = -0.5d0*lqr(2,3)
        lqr(3,3) =  -2.d0*lqr(2,1)
        lqr(4,1) =  cs23 - cs12 
        lqr(4,2) =  -lqr(4,1)
        lqr(4,3) =  (cc12-ss12) - (cc23-ss23)
        lqr(5,1) =  (cc12*x12 + cc23*x23)*0.5d0
        lqr(5,2) =  (ss12*x12 + ss23*x23)*0.5d0
        lqr(5,3) =   ss12*y21 + ss23*y32
        lqr(6,1) = -(cc12*y21 + cc23*y32)*0.5d0
        lqr(6,2) = -0.5d0*lqr(5,3)
        lqr(6,3) =  -2.d0*lqr(5,1)
        lqr(7,1) =  cs31 - cs23
        lqr(7,2) =  -lqr(7,1)
        lqr(7,3) =  (cc23-ss23) - (cc31-ss31)
        lqr(8,1) =  (cc23*x23 + cc31*x31)*0.5d0
        lqr(8,2) =  (ss23*x23 + ss31*x31)*0.5d0
        lqr(8,3) =   ss23*y32 + ss31*y13
        lqr(9,1) = -(cc23*y32 + cc31*y13)*0.5d0
        lqr(9,2) = -0.5d0*lqr(8,3)
        lqr(9,3) =  -2.d0*lqr(8,1)
      end if
C
      do 1600  j = 1,9
        l(j,1) =   clr*llr(j,1) + cqr*lqr(j,1)
        l(j,2) =   clr*llr(j,2) + cqr*lqr(j,2)
        l(j,3) =   clr*llr(j,3) + cqr*lqr(j,3)
 1600   continue
c
      c =        2.d0*f/area2
      db11 =     c*db(1,1)
      db22 =     c*db(2,2)
      db33 =     c*db(3,3)
      db12 =     c*db(1,2)
      db13 =     c*db(1,3)
      db23 =     c*db(2,3)
      do 4000  j = 1,9
        jj =     ls(j)
        s1 =     db11*l(j,1) + db12*l(j,2) + db13*l(j,3)
        s2 =     db12*l(j,1) + db22*l(j,2) + db23*l(j,3)
        s3 =     db13*l(j,1) + db23*l(j,2) + db33*l(j,3)
        do 3500  i = 1,j
          ii =      ls(i)
          sm(jj,ii) =  sm(jj,ii) + (s1*l(i,1)+s2*l(i,2)+s3*l(i,3))
          sm(ii,jj) =  sm(jj,ii)
 3500     continue
 4000   continue
      return
      end
C=END FORTRAN
C=DECK SM3SHMEMBB
C=PURPOSE Form basic membrane stiffness of 9-dof triangle
C=AUTHOR C. A. Felippa, June 1984
C=VERSION June 1984
C=EQUIPMENT Machine independent
C=KEYWORDS finite element membrane plane stress
C=KEYWORDS basic material stiffness matrix
C=BLOCK ABSTRACT
C
C     SM3SHMEMBB forms the material element stiffness matrix 
C     associated with the basic displacement modes 
C     (rigid modes + constant strain modes) of a 9-dof 
C     plane-stress triangle based on the free formulation 
C     of Bergan and Nygard.
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       CALL      SM3SHMEMBB (X, Y, DM, ALPHA, F, LS, SM, M, STATUS)
C
C     where the input arguments are
C
C     X      (3 x 1) array of x coordinates of triangle nodes.
C     Y      (3 x 1) array of y coordinates of triangle nodes.
C     DM     (3 x 3) matrix relating in-plane forces to strains.
C     ALPHA  Rotational lumping factor; if zero form CST.
C     F      Factor by which stiffness entries will be multiplied.
C     LS     (9 x 1) array of stiffness location pointers, for SM
C     SM     Incoming material stiffness array.
C     M      First dimension of SM in calling program.
C
C     The outputs are:
C
C     SM     Output stiffness array with basic stiffness
C            coefficients added in.  The (i,j)-th entry of the
C            basic element stiffness is added to SM(K,L),
C            where K=LS(I) and L=LS(J).
C
C     STATUS  Status character variable.  Blank if no error detected.
C
C     Note: the internal ordering of degrees of freedom is
C
C          ux1,uy1, ux2,uy2, ux3,uy3, thetaz1,thetaz2,thetaz3
C
C     This simplifies the stiffness formation loop if alpha=0.
C     
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    SM3SHMEMBB
     &             (x, y, dm, alpha, f, ls, sm, m, status)
C
C                   A R G U M E N T S
C
      character*(*)      status
      integer            m, ls(9)
      double precision   x(3), y(3), dm(3,3), alpha, f, sm(m,m)
C
C                   L O C A L   V A R I A B L E S
C
      double precision   area2, c, p(9,3)
      double precision   dm11, dm12, dm13, dm22, dm23, dm33
      double precision   x21, x32, x13, y21, y32, y13
      double precision   x12, x23, x31, y12, y23, y31
      double precision   s1, s2, s3
      integer            i, j, k, l, n
C
C                   L O G I C
C
      status =   ' '
      x21 =      x(2) - x(1)
      x12 =     -x21
      x32 =      x(3) - x(2)
      x23 =     -x32
      x13 =      x(1) - x(3)
      x31 =     -x13
      y21 =      y(2) - y(1)
      y12 =     -y21
      y32 =      y(3) - y(2)
      y23 =     -y32
      y13 =      y(1) - y(3)
      y31 =     -y13
      area2 =    y21*x13 - x21*y13
      if (area2 .le. 0.0)      then
        status = 'SM3SHMEMBB: Negative area'
        if (area2 .eq. 0.0)   status = 'SM3SHMEMBB: Zero area'
        return
      end if
      p(1,1) =   y23
      p(2,1) =   0.0
      p(3,1) =   y31
      p(4,1) =   0.0
      p(5,1) =   y12
      p(6,1) =   0.0
      p(1,2) =   0.0
      p(2,2) =   x32
      p(3,2) =   0.0
      p(4,2) =   x13
      p(5,2) =   0.0
      p(6,2) =   x21
      p(1,3) =   x32
      p(2,3) =   y23
      p(3,3) =   x13
      p(4,3) =   y31
      p(5,3) =   x21
      p(6,3) =   y12
      n =        6
      if (alpha .ne. 0.0)         then
        c =      alpha/6.0d0
        p(7,1) =   y23*(y13-y21)    *c
        p(7,2) =   x32*(x31-x12)    *c
        p(7,3) =  (x31*y13-x12*y21) *c*2.0d0
        p(8,1) =   y31*(y21-y32)    *c
        p(8,2) =   x13*(x12-x23)    *c
        p(8,3) =  (x12*y21-x23*y32) *c*2.0d0
        p(9,1) =   y12*(y32-y13)    *c
        p(9,2) =   x21*(x23-x31)    *c
        p(9,3) =  (x23*y32-x31*y13) *c*2.0d0
        n =        9
      end if
      c =        0.5d0*f/area2
      dm11 =     c * dm(1,1)
      dm22 =     c * dm(2,2)
      dm33 =     c * dm(3,3)
      dm12 =     c * dm(1,2)
      dm13 =     c * dm(1,3)
      dm23 =     c * dm(2,3)
      do 3000  j = 1,n
        l =      ls(j)
        s1 =     dm11*p(j,1) + dm12*p(j,2) + dm13*p(j,3)
        s2 =     dm12*p(j,1) + dm22*p(j,2) + dm23*p(j,3)
        s3 =     dm13*p(j,1) + dm23*p(j,2) + dm33*p(j,3)
        do 2500  i = 1,j
          k =      ls(i)
          sm(k,l) =  sm(k,l) + (s1*p(i,1) + s2*p(i,2) + s3*p(i,3))
          sm(l,k) =  sm(k,l)
 2500     continue
 3000   continue
      return
      end
C=END FORTRAN


C=DECK SM3SHMEMBH
C=PURPOSE Form H.O. material stiffness of 9-dof ANDES membrane tria.
C=AUTHOR C. A. Felippa, June 1991
C=VERSION July 1991
C=EQUIPMENT Machine independent
C=KEYWORDS finite element
C=KEYWORDS material stiffness matrix high-order
C=KEYWORDS triangle membrane assumed natural deviatoric strain
C=BLOCK ABSTRACT
C
C     SM3MANDES forms the higher order element stiffness matrix
C     of a 9-dof membrane triangle based on the ANDES formulation.
C     Implementation moderately optimized for speed.
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       call SM3SHMEMBH (X, Y, DM, F, LS, SM, M, STATUS)
C
C     The inputs are:
C
C     X    (3 x 1) array of x coordinates of triangle nodes
C     Y    (3 x 1) array of y coordinates of triangle nodes
C     DM   (3 x 3) membrane constitutive matrix already
C          integrated through the thickness
C     F    Factor by which all stiffness entries will
C          be multiplied.
C     SM   Incoming material stiffness array.
C     LS   (9 x 1) array of stiffness location pointers
C          (see examples in SM3SHMEMBB)
C     M    First dimension of SM in calling program.
C
C     The outputs are:
C
C     SM      Output stiffness array with higher order stiffness
C             coefficients added in.
C             The (i,j)-th entry of the basic element stiffness is 
C             added to SM(K,L), where K=LS(I) and L=LS(J).
C     STATUS  Status character variable.  Blank if no error
C             detected.
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    SM3SHMEMBH
     &              (x, y, dm, f, ls, sm, m, status)
C
C                   A R G U M E N T S
C
      integer           ls(9), m
      double precision  x(3),y(3), dm(3,3), f, sm(m,m)
      character         status*(*)
C
C                   L O C A L   V A R I A B L E S
C
      double precision  x12, x21, x23, x32, x31, x13
      double precision  y12, y21, y23, y32, y31, y13
      double precision  l21,l32,l13
      double precision  chi213,chi321,chi132
      double precision  area, area2, area43
      double precision  c(3,3), e(3,3), et(3), d(3), qm(3,3,3)
      double precision  t(3,3), tfac, kth(3,3)
      double precision  s(3), xyij(6), sum, w(3), wfac
      integer           i, j, k, l
C
C                   L O G I C
C
      status =   ' '
      if (f .eq. 0.0)         return
      x12 =      x(1) - x(2)
      x21 =     -x12
      x23 =      x(2) - x(3)
      x32 =     -x23
      x31 =      x(3) - x(1)
      x13 =     -x31
      y12 =      y(1) - y(2)
      y21 =     -y12
      y23 =      y(2) - y(3)
      y32 =     -y23
      y31 =      y(3) - y(1)
      y13 =     -y31
      area2 =    x21*y31-x31*y21
      if (area2 .le. 0.0)      then
        status = 'SM3SHMEMBH: Negative area'
        if (area2 .eq. 0.0)   status = 'SM3SHMEMBH: Zero area'
        return
      end if
      area  =     0.5d0*area2
      l21 =       sqrt(x21**2+y21**2)
      l32 =       sqrt(x32**2+y32**2)
      l13 =       sqrt(x13**2+y13**2)
      tfac =      0.25d0/area**2
      t(1,1) =    tfac*y23*y13*l21**2
      t(1,2) =    tfac*y31*y21*l32**2
      t(1,3) =    tfac*y12*y32*l13**2
      t(2,1) =    tfac*x23*x13*l21**2
      t(2,2) =    tfac*x31*x21*l32**2
      t(2,3) =    tfac*x12*x32*l13**2
      t(3,1) =    tfac*(y23*x31+x32*y13)*l21**2
      t(3,2) =    tfac*(y31*x12+x13*y21)*l32**2
      t(3,3) =    tfac*(y12*x23+x21*y32)*l13**2
      wfac =      0.75d0*f*area
      e(1,1) =    wfac*dm(1,1)
      e(1,2) =    wfac*dm(1,2)
      e(1,3) =    wfac*dm(1,3)
      e(2,1) =    wfac*dm(2,1)
      e(2,2) =    wfac*dm(2,2)
      e(2,3) =    wfac*dm(2,3)
      e(3,1) =    wfac*dm(3,1)
      e(3,2) =    wfac*dm(3,2)
      e(3,3) =    wfac*dm(3,3)
      do 1600  j = 1,3
        do 1400  i = 1,3
          et(i) =  e(i,1)*t(1,j)+e(i,2)*t(2,j)+e(i,3)*t(3,j)
 1400     continue
        do 1500  i = 1,3
          c(i,j) = t(1,i)*et(1)+t(2,i)*et(2)+t(3,i)*et(3)
 1500     continue
 1600   continue
      area43 =     (2.d0/3.d0)*area2
      chi213 =     area43/l21**2
      chi321 =     area43/l32**2
      chi132 =     area43/l13**2
      qm(1,1,1) =  -0.25d0*chi213
      qm(1,2,1) =  -qm(1,1,1)
      qm(1,3,1) =   0.0
      qm(2,1,1) =   0.25d0*chi321
      qm(2,2,1) =   0.50d0*chi321
      qm(2,3,1) =   qm(2,1,1)
      qm(3,1,1) =  -0.50d0*chi132
      qm(3,2,1) =  -0.25d0*chi132
      qm(3,3,1) =   qm(3,2,1)
      qm(1,1,2) =  -0.25d0*chi213
      qm(1,2,2) =  -0.50d0*chi213
      qm(1,3,2) =   qm(1,1,2)
      qm(2,1,2) =   0.0
      qm(2,2,2) =  -0.25d0*chi321
      qm(2,3,2) =  -qm(2,2,2)
      qm(3,1,2) =   0.25d0*chi132
      qm(3,2,2) =   qm(3,1,2)
      qm(3,3,2) =   0.50d0*chi132
      qm(1,1,3) =   0.50d0*chi213
      qm(1,2,3) =   0.25d0*chi213
      qm(1,3,3) =   qm(1,2,3)
      qm(2,1,3) =  -0.25d0*chi321
      qm(2,2,3) =   qm(2,1,3)
      qm(2,3,3) =  -0.50d0*chi321
      qm(3,1,3) =   0.25d0*chi132
      qm(3,2,3) =   0.0
      qm(3,3,3) =  -qm(3,1,3)
      kth(1,1) =   0.0
      kth(1,2) =   0.0
      kth(1,3) =   0.0
      kth(2,2) =   0.0
      kth(2,3) =   0.0
      kth(3,3) =   0.0
      do 2800  k = 1,3
        do 2600  j = 1,3
          d(1) = c(1,1)*qm(1,j,k)+c(1,2)*qm(2,j,k)+c(1,3)*qm(3,j,k)
          d(2) = c(2,1)*qm(1,j,k)+c(2,2)*qm(2,j,k)+c(2,3)*qm(3,j,k)
          d(3) = c(3,1)*qm(1,j,k)+c(3,2)*qm(2,j,k)+c(3,3)*qm(3,j,k)
          do 2500  i = 1,j
            kth(i,j) = kth(i,j) +
     $                 qm(1,i,k)*d(1)+qm(2,i,k)*d(2)+qm(3,i,k)*d(3)
            kth(j,i) = kth(i,j)
 2500       continue
 2600     continue
 2800   continue
      s(1) =   kth(1,1) + kth(1,2) + kth(1,3)
      s(2) =   kth(2,1) + kth(2,2) + kth(2,3)
      s(3) =   kth(3,1) + kth(3,2) + kth(3,3)
      xyij(1) =  0.25d0*x32/area
      xyij(2) =  0.25d0*y32/area
      xyij(3) =  0.25d0*x13/area
      xyij(4) =  0.25d0*y13/area
      xyij(5) =  0.25d0*x21/area
      xyij(6) =  0.25d0*y21/area
      do 4000  j = 1,9
        l =    ls(j)
        do 3600  i = 1,3
          if (j .le. 6)      then
            w(i) =  s(i)*xyij(j)
          else
            w(i) =  kth(i,j-6)
          end if
 3600     continue
        sum =    w(1) + w(2) + w(3)
        do 3700  i = 1,j
          k =      ls(i)
          if (i .le. 6)      then
             sm(k,l) =  sm(k,l) + sum*xyij(i)
          else
             sm(k,l) =  sm(k,l) + w(i-6)
          end if
          sm(l,k) =  sm(k,l)
 3700     continue
 4000   continue
      return
      end
C=END FORTRAN
C=DECK SM3SHBENDH
C=PURPOSE Form Kh for 3-node Kirchhoff ANDES plate bending elem
C=AUTHOR C. Militello and C. Felippa, May 1989
C=VERSION July 1989
C=EQUIPMENT Machine independent
C=KEYWORDS thin plate bending
C=KEYWORDS finite element triangle higher order stiffness matrix
C=BLOCK ABSTRACT
C
C     SM3SHBENDH forms the higher order material stiffness matrix of a
C     9-dof thin-plate-bending triangle obtained by using linear
C     curvatures over the sides.  Implementation optimized for speed
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       CALL      SM3SHBENDH (X, Y, DB, F, LS, SM, M, STATUS)
C
C     where the input arguments are
C
C       X     (3 x 1) array of x coordinates of triangle nodes
C       Y     (3 x 1) array of y coordinates of triangle nodes
C       DB    (3 x 3) moment-curvature matrix.
C       F     Factor by which stiffness entries will be multiplied.
C       LS    (9 x 1) array of stiffness location pointers
C             (see Output SM).
C       SM    Incoming material stiffness array.
C       M     First dimension of SM in calling program.
C
C     The outputs are:
C
C       SM      Output stiffness array with higher order stiffness
C               coefficients added in.  The (i,j)-th entry of the
C               (9 by 9) element bending stiffness is added to
C               SM(K,L), where K=LS(I) and L=LS(J).
C       STATUS  Status character variable.  Blank if no error
C               detected.
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine  SM3SHBENDH (x, y, db, f, ls, sm, m, status)
C
C                   A R G U M E N T S
C
      character*(*)     status
      integer           ls(9), m
      double precision  x(3), y(3), db(3,3), f
      double precision  sm(m,m)
C
C                   L O C A L   V A R I A B L E S
C
      double precision  x21, x32, x13, y21, y32, y13
      double precision  x12, x23, x31, y12, y23, y31
      double precision  area, area2, tfac, cfac
      double precision  lam12, lam23, lam31
      double precision  dbt11, dbt12, dbt13, dbt21, dbt22, dbt23
      double precision  dbt31, dbt32, dbt33
      double precision  t(3,3), sh(9,9)
      double precision  mu11, mu12, mu13, mu22, mu23, mu33
      double precision  c11, c12, c13, c22, c23, c33
      integer           i, j, k, l
C
C                   L O G I C
C
      status =   ' '
      x21 =  x(2) - x(1)
      x12 =  -x21
      x32 =  x(3) - x(2)
      x23 =  -x32
      x13 =  x(1) - x(3)
      x31 =  -x13
      y21 =  y(2) - y(1)
      y12 =  -y21
      y32 =  y(3) - y(2)
      y23 =  -y32
      y13 =  y(1) - y(3)
      y31 =  -y13
      area2 =    y21*x13 - x21*y13
      if (area2 .le. 0.0)      then
        status = 'SM3SHBENDH: Negative area'
        if (area2 .eq. 0.0)   status = 'SM3SHBENDH: Zero area'
        return
      end if
      lam12 =  (x12*x13+y12*y13)/(x21*x21+y21*y21)
      lam23 =  (x23*x21+y23*y21)/(x32*x32+y32*y32)
      lam31 =  (x31*x32+y31*y32)/(x13*x13+y13*y13)
      mu11 =   2.d0*(lam12**2+1-lam12)
      mu22 =   2.d0*(lam23**2+1-lam23)
      mu33 =   2.d0*(lam31**2+1-lam31)
      mu12 =   2.d0*lam12-1-(1+lam12)*lam23
      mu23 =   2.d0*lam23-1-(1+lam23)*lam31
      mu13 =   2.d0*lam31-1-(1+lam31)*lam12
      tfac =     1.d0/area2**2
      t(1,1) =   tfac* y23*y13 
      t(1,2) =   tfac* y31*y21 
      t(1,3) =   tfac* y12*y32
      t(2,1) =   tfac* x23*x13 
      t(2,2) =   tfac* x31*x21 
      t(2,3) =   tfac* x12*x32 
      t(3,1) =   tfac* (x31*y23+x32*y13)
      t(3,2) =   tfac* (x12*y31+x13*y21) 
      t(3,3) =   tfac* (x23*y12+x21*y32) 
      dbt11 =   db(1,1)*t(1,1)+db(1,2)*t(2,1)+db(1,3)*t(3,1)
      dbt12 =   db(1,1)*t(1,2)+db(1,2)*t(2,2)+db(1,3)*t(3,2)
      dbt13 =   db(1,1)*t(1,3)+db(1,2)*t(2,3)+db(1,3)*t(3,3)
      dbt21 =   db(2,1)*t(1,1)+db(2,2)*t(2,1)+db(2,3)*t(3,1)
      dbt22 =   db(2,1)*t(1,2)+db(2,2)*t(2,2)+db(2,3)*t(3,2)
      dbt23 =   db(2,1)*t(1,3)+db(2,2)*t(2,3)+db(2,3)*t(3,3)
      dbt31 =   db(3,1)*t(1,1)+db(3,2)*t(2,1)+db(3,3)*t(3,1)
      dbt32 =   db(3,1)*t(1,2)+db(3,2)*t(2,2)+db(3,3)*t(3,2)
      dbt33 =   db(3,1)*t(1,3)+db(3,2)*t(2,3)+db(3,3)*t(3,3)
      area =    0.5d0*area2
      cfac =    area*f
      c11 =    cfac*mu11*( t(1,1)*dbt11+t(2,1)*dbt21+t(3,1)*dbt31 )
      c12 =    cfac*mu12*( t(1,2)*dbt11+t(2,2)*dbt21+t(3,2)*dbt31 )
      c13 =    cfac*mu13*( t(1,3)*dbt11+t(2,3)*dbt21+t(3,3)*dbt31 )
      c22 =    cfac*mu22*( t(1,2)*dbt12+t(2,2)*dbt22+t(3,2)*dbt32 )
      c23 =    cfac*mu23*( t(1,3)*dbt12+t(2,3)*dbt22+t(3,3)*dbt32 )
      c33 =    cfac*mu33*( t(1,3)*dbt13+t(2,3)*dbt23+t(3,3)*dbt33 )
      sh(1,1) =  4.d0*(c33-c13-c13+c11)
      sh(1,2) =  2.d0*((c11-c13)*y21+(c13-c33)*y13)
      sh(1,3) =  2.d0*((c13-c11)*x21+(c33-c13)*x13)
      sh(1,4) =  4.d0*(-c23+c13+c12-c11)
      sh(1,5) =  2.d0*((c12-c23)*y32+(c11-c13)*y21)
      sh(1,6) =  2.d0*((c23-c12)*x32+(c13-c11)*x21)
      sh(1,7) =  4.d0*(-c33+c23+c13-c12)
      sh(1,8) =  2.d0*((c12-c23)*y32+(c13-c33)*y13)
      sh(1,9) =  2.d0*((c23-c12)*x32+(c33-c13)*x13)
      sh(2,2) =  c11*y21**2+2.d0*c13*y13*y21+c33*y13**2
      sh(2,3) =  (-c11*x21-c13*x13)*y21+(-c13*x21-c33*x13)*y13
      sh(2,4) =  2.d0*((c12-c11)*y21+(c23-c13)*y13)
      sh(2,5) =  (c12*y21+c23*y13)*y32+c11*y21**2+c13*y13*y21
      sh(2,6) =  (-c12*x32-c11*x21)*y21+(-c23*x32-c13*x21)*y13
      sh(2,7) =  2.d0*((c13-c12)*y21+(c33-c23)*y13)
      sh(2,8) =  (c12*y21+c23*y13)*y32+c13*y13*y21+c33*y13**2
      sh(2,9) =  (-c12*x32-c13*x13)*y21+(-c23*x32-c33*x13)*y13
      sh(3,3) =  c11*x21**2+2.d0*c13*x13*x21+c33*x13**2
      sh(3,4) =  2.d0*((c11-c12)*x21+(c13-c23)*x13)
      sh(3,5) =  (-c12*x21-c23*x13)*y32+(-c11*x21-c13*x13)*y21
      sh(3,6) =  (c12*x21+c23*x13)*x32+c11*x21**2+c13*x13*x21
      sh(3,7) =  2.d0*((c12-c13)*x21+(c23-c33)*x13)
      sh(3,8) =  (-c12*x21-c23*x13)*y32+(-c13*x21-c33*x13)*y13
      sh(3,9) =  (c12*x21+c23*x13)*x32+c13*x13*x21+c33*x13**2
      sh(4,4) =  4.d0*(c22-c12-c12+c11)
      sh(4,5) =  2.d0*((c22-c12)*y32+(c12-c11)*y21)
      sh(4,6) =  2.d0*((c12-c22)*x32+(c11-c12)*x21)
      sh(4,7) =  4.d0*(c23-c22-c13+c12)
      sh(4,8) =  2.d0*((c22-c12)*y32+(c23-c13)*y13)
      sh(4,9) =  2.d0*((c12-c22)*x32+(c13-c23)*x13)
      sh(5,5) =  c22*y32**2+2.*c12*y21*y32+c11*y21**2
      sh(5,6) =  (-c22*x32-c12*x21)*y32+(-c12*x32-c11*x21)*y21
      sh(5,7) =  2.d0*((c23-c22)*y32+(c13-c12)*y21)
      sh(5,8) =  c22*y32**2+(c12*y21+c23*y13)*y32+c13*y13*y21
      sh(5,9) =  (-c22*x32-c23*x13)*y32+(-c12*x32-c13*x13)*y21
      sh(6,6) =  c22*x32**2+2.*c12*x21*x32+c11*x21**2
      sh(6,7) =  2.d0*((c22-c23)*x32+(c12-c13)*x21)
      sh(6,8) =  (-c22*x32-c12*x21)*y32+(-c23*x32-c13*x21)*y13
      sh(6,9) =  c22*x32**2+(c12*x21+c23*x13)*x32+c13*x13*x21
      sh(7,7) =  4.d0*(c33-c23-c23+c22)
      sh(7,8) =  2.d0*((c23-c22)*y32+(c33-c23)*y13)
      sh(7,9) =  2.d0*((c22-c23)*x32+(c23-c33)*x13)
      sh(8,8) =  c22*y32**2+2.*c23*y13*y32+c33*y13**2
      sh(8,9) =  (-c22*x32-c23*x13)*y32+(-c23*x32-c33*x13)*y13
      sh(9,9) =  c22*x32**2+2.*c23*x13*x32+c33*x13**2
      do 3500  i = 1,9
        k =   ls(i)
        do 3400  j = i,9
          l =   ls(j)
            sm(k,l) =   sm(k,l) + sh(i,j)
            sm(l,k) =   sm(k,l)
 3400      continue
 3500    continue
      return
      end
C=END FORTRAN
c

