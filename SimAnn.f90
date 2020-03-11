! ABSTRACT:
!   Simulated annealing is a global optimization method that distinguishes
!   between different local optima. Starting from an initial point, the
!   algorithm takes a step and the function is evaluated. When minimizing a
!   function, any downhill step is accepted and the process repeats from this
!   new point. An uphill step may be accepted. Thus, it can escape from local
!   optima. This uphill decision is made by the Metropolis criteria. As the
!   optimization process proceeds, the length of the steps decline and the
!   algorithm closes in on the global optimum. Since the algorithm makes very
!   few assumptions regarding the function to be optimized, it is quite
!   robust with respect to non-quadratic surfaces. The degree of robustness
!   can be adjusted by the user. In fact, simulated annealing can be used as
!   a local optimizer for difficult functions.
!
!   This implementation of simulated annealing was used in "Global Optimizatio
!   of Statistical Functions with Simulated Annealing," Goffe, Ferrier and
!   Rogers, Journal of Econometrics, vol. 60, no. 1/2, Jan./Feb. 1994, pp.
!   65-100. Briefly, we found it competitive, if not superior, to multiple
!   restarts of conventional optimization routines for difficult optimization
!   problems.
!
!   For more information on this routine, contact its author:
!   Bill Goffe, bgoffe@whale.st.usm.edu

    
MODULE SIMANN

IMPLICIT NONE

PRIVATE
PUBLIC:: SA,GSA,PRTVEC,optifunc

abstract interface
    REAL(8) function optifunc(x) 
        real(8), DIMENSION(:),intent(in) :: x
    end function optifunc
end interface


CONTAINS    
    

!            Generalized Simulated Annealing - Code
!    Program developed by Members of Lab. of Molecular Modeling
!    Kleber C. Mundim (1995)
!            ------------------------------------------------------
!       Global optimization method using Generalized Simulated Annealing
!    First Version   Jan./1994   (Kleber Mundim and Constantino Tsallis)
!    Second Version  Jun./1995   (Marcelo Moret)
!    Third Version   Sep./1995   (Thierry Lemaire and Amin Bassrei)
!
!    Some Basic References and literature citation:
!
!1-  Title   : Geometry Optimization and Conformational Analysis Through
!              Generalized Simulated Annealing
!    Authors : Kleber C. Mundim and Constantino Tsallis
!    Journal : Int.Journal of Quantum Chemistry, 58 (1996),373-381
!
!2-  Title   : Stochastic Molecular Optimization using Generalized
!              Simulated Annealing
!    Authors : Marcelo Moret, Pedro G. Pascutti, Paulo M. Bish
!              and Kleber C. Mundim
!    Journal : Journal of Computational Chemistry, 19 (1998) 647-657
!
!3-  Title   : Modeling Gravity Anomalies Through Generalized
!              Simulated Annealing
!    Authors : Kleber C. Mundim, Thierry Lemaire and Amin Bassrei
!    Journal : Physica A, 252 (1998) 405-416
!            ------------------------------------------------------
!     Most important vectors and parameters used in the GSA routines.
!     NDimension   -->   Number of parameters to be optimized (problem dependent).
!     X_t(i)       -->   This vector contains the parameters.
!     X_0(i)       -->   X_0(i) contains X_t at the time t-1.
!     X_Min        -->   Minimal parameters obtained.
!     To           -->   Initial Temperature.
!     qV			 -->   qV caracterize the Visiting Probability Function.
!     qA			 -->   qA Acceptance Parameter.
!	qT			 -->   qT Temperature parameter
!	Maxevl     -->   Maximum number of GSA loops (problem dependent).
 
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * |
    SUBROUTINE GSA(FCN,TO,Fopt,Xopt,QAI,QTI,QVI,EPSI,MaxevlI,XSTART,LbI,UbI,IsSGSAI,PITERI,NITERI)
    !  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * |
    !     This routine start the GSA-loop
 
    IMPLICIT NONE
    procedure(optifunc) :: FCN
    REAL(8),DIMENSION(:),INTENT(IN),OPTIONAL:: XSTART,LBI,UBI
    REAL(8),DIMENSION(:),INTENT(OUT)::Xopt
    REAL(8),INTENT(in),OPTIONAL:: EpsI,QAI,QTI,QVI 
    REAL(8),INTENT(inout):: TO
    REAL(8),INTENT(OUT)::Fopt      
    INTEGER,TARGET,INTENT(IN),OPTIONAL::MaxevlI,PITERI,NITERI     
    LOGICAL,INTENT(IN),OPTIONAL::ISSGSAI
    
    REAL(8)::QA,QT,QV,EPS
    INTEGER,TARGET::Maxevl,PITER,NITER
    LOGICAL::ISSGSA=.TRUE.
!  Type all internal variables.
    REAL(8)::X_T(SIZE(XOPT)),X_0(SIZE(XOPT)),LB(SIZE(XOPT)),UB(SIZE(XOPT)) 
    REAL(8),ALLOCATABLE::FOPT1(:)
      

    REAL(8):: COEf , D , deltae , EXP1 , EXP2 , func_0 ,&
                    & func_t , one , oneqa1 , pqa ,    &
                    & QA1 , QT1 , QV1 , QV2 
    REAL(8):: rand , t , time , TOScale , TQT , tup ,FOPT_PITER,T1
                     
    INTEGER i ,  NDImension,NC1,NC2 
    INTEGER,POINTER::NCYCLE
    
    DATA one/1.0D+00/
 
    CALL GSAINI()    ! Read initial GSA parameters in the file "GSA.in"
 
! IMPORTANT!  Is it necessary to initialize the vectors X_0 e X_Min.
 

    XOPT = X_0
    X_T = X_0
    FOPT1=1.D10
    func_0 = FCN(X_0)
    FOPT = func_0
    func_t = func_0
    FOPT1(1)=FOPT
    FOPT_PITER=FOPT
    NC1=0
    oneqa1 = one/QA1
 
    time = 0.0D0
    ncycle = 0
 
! Starting minimization loop
 
      DO WHILE ( ncycle<=Maxevl )
 
         time = time + one
         ncycle = ncycle + 1
 
         t = TQT/((one+time)**QT1-one)   ! Temperature T(t)
         IF ( D.EQ.0.0D0 ) THEN
            tup = one
         ELSE
            tup = t**(D/(QT-3.0D0))
         ENDIF
 
!--Create the new vector  X_t  using the g(qV,t) function.
         CALL GFUNC(t,tup)
 
!--Evaluate the Potential function in terms of X_t.
         func_t = FCN(X_T)
 
!--Verify if the Potential(X_n+1) is going down. In this case
!  change  X(n) ---> X(n+1). In the other case retains X(n).
 
        IF ( func_t<func_0 ) THEN
        X_0 = X_T
        func_0 = func_t
        IF ( func_t<=FOPT ) THEN
            FOPT = func_t
            XOPT = X_T
            NC1=NC1+1
            NC2=MIN(NC1,10)
            FOPT1(NC2:2:-1)=FOPT1(NC2-1:1:-1)
            FOPT1(1)=FOPT
            IF(NC2==10) THEN
                T1=ABS(FOPT1(10)-FOPT)/MAX(FOPT,1.D-8)
                IF(T1<EPS) THEN
                    WRITE(*,100) T1,EPS,CHAR(10),NC1,NCYCLE,T  
                    RETURN
                ENDIF
                
            ENDIF
        ENDIF
        ELSE
        !            Evaluate the Acceptance Probability Function [0,1].
        deltae = func_t - func_0
        !---         Xiang (Phys. Let. A 233, (1997) 216-220)
        !		   qA = -3 -0.85*Time
        QA1 = QA - one
        !---
        !            PqA = One/(One+(One+qA1*DeltaE/T)**OneqA1)
        pqa = one/((one+QA1*deltae/t)**oneqa1)   ! Corre玢o
        CALL RANDOM_NUMBER(rand)
        !            If  Rand > Acc_Prob  retains   X_0.
        !            If  Rand < Acc_Prob  replace   X_0 by the new X_t.
        IF ( rand.LT.pqa ) THEN
                X_0 = X_T;func_0 = func_t
        ENDIF
        ENDIF
         
        !IF(MOD(NCYCLE,PITER)==0) THEN
        !    T1=ABS(FOPT-FOPT_PITER)/MAX(FOPT,1.D-8)
        !    IF(T1<EPS) THEN                    
        !        WRITE(*,100) T1,EPS,PITER,CHAR(10),NCYCLE,T    
        !        RETURN
        !    ELSE                
        !        FOPT_PITER=FOPT
        !    ENDIF
        !                               
        !ENDIF
       
        IF(NCYCLE==MAXEVL) THEN
            WRITE(*,200) NCYCLE,T    
            RETURN
        ENDIF

         
      ENDDO
      
      IF(ALLOCATED(FOPT1)) DEALLOCATE(FOPT1)
      
100 FORMAT('EXIT WITH OPTIMIZATIONE RATIO ABS(FOPT(I)-FOPT(I-10))/MAX(FOPT(I),1.D-8)=',E12.5,'<',E12.5,' PER 10 SUCCESSIVE UPDATED.',A1,'IUPDATED=',I10,',NITER=',I10,' ,TEMPERATURE=',E12.5) 
200 FORMAT('EXIT WITH MAXEVL REACHED.,MAXEVL=',I10,',TEMPERATURE=',E12.5)

CONTAINS

!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * |
      SUBROUTINE GSAINI()
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * |
!     This routine read and initialize the GSA parameters
!     The number of parameters to be optimized (NDIMENSION) is problem
!     dependent.
 
      IMPLICIT NONE
      REAL(8):: coef1, gamadown, gamaup,pi,tmp 

      !CALL RANDOM_SEED()
      
      NDImension=SIZE(XOPT)
      IF(PRESENT(PITERI)) THEN
        PITER=PITERI
      ELSE
        PITER=100
      ENDIF 
      IF(ALLOCATED(FOPT1)) DEALLOCATE(FOPT1)
      ALLOCATE(FOPT1(100))
      
      IF(PRESENT(QAI)) THEN
        QA=QAI
      ELSE
        QA=1.5
        !QA=-5
      ENDIF
      IF(PRESENT(QVI)) THEN
        QV=QVI
      ELSE
        QV=1.5
        !QV=2.62
      ENDIF 
      IF(PRESENT(QTI)) THEN
        QT=QTI
      ELSE
        QT=1.5
      ENDIF   
      IF(PRESENT(MaxevlI)) THEN
        Maxevl=MaxevlI
      ELSE
        Maxevl=1000000
      ENDIF 
      IF(PRESENT(EPSI)) THEN
        EPS=EPSI
      ELSE
        EPS=1.0D-4
      ENDIF       
    IF(PRESENT(LBI)) THEN
        LB=LBI
    ELSE
        LB=-1.D10
    ENDIF    
    IF(PRESENT(UBI)) THEN
        UB=UBI
    ELSE
        UB=1.D10
    ENDIF
 
    IF(PRESENT(XSTART)) THEN
        X_0=XSTART
    ELSE
        CALL RANDOM_NUMBER(X_0)
        X_0=LB+X_0*(UB-LB)
    ENDIF      
    IF(PRESENT(ISSGSAI)) THEN
        ISSGSA=ISSGSAI
    ELSE
        ISSGSA=.TRUE.
        
    ENDIF
    
    IF(ISSGSA) THEN
      
        D = 0    ! Dimension problem 
    ELSE
        D=NDImension
        !QV=2.62;QA=-5
    ENDIF
    
    IF(PRESENT(NITERI)) THEN
        NCYCLE=>NITERI
    ELSE
        NCYCLE=>NITER    
    ENDIF
    
      pi = 3.14159265359D0 
!    Acceptance probability
      QA1 = QA - one 
!    Temperature
      QT1 = QT - one
      TQT = TO*(2.0D0**QT1-one) 
!    Visiting probability
      QV1 = QV - one
      QV2 = 2.0D0**QV1 - one
      tmp = one/QV1 - 0.5D0
      gamadown = DGAMMA(tmp)
      EXP1 = 2.0D0/(3.0D0-QV)
 
      IF ( D.EQ.0.0D0 ) THEN
         coef1 = one
         EXP2 = one/QV1 - 0.5D0
         gamaup = gamadown
      ELSE
         coef1 = (QV1/pi)**(D*0.5D0)
         EXP2 = one/QV1 + 0.5D0*D - 0.5D0
         gamaup = DGAMMA(EXP2)
      ENDIF
 
      COEf = coef1*gamaup/gamadown
 
 
      END
      

!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * |
      SUBROUTINE GFUNC(T,Tup)
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * |
!     This routine evaluate the new set o parameters using the
!     visiting probability function g(q,T).
 
      IMPLICIT NONE
      REAL(8),INTENT(IN)::T,Tup 
      REAL(8)::deltax,r , s,RND1   

 
      DO i = 1 , NDImension
         CALL RANDOM_NUMBER(R)
         CALL RANDOM_NUMBER(S)
         deltax = COEf*Tup/(one+QV1*r*r/T**EXP1)**EXP2
         IF ( s.LE.0.5 ) deltax = -deltax
         X_T(i) = X_0(i) + deltax
        !  If X_T is out of bounds, select a point in bounds for the trial.
        IF ( (X_T(i).LT.Lb(i)) .OR. (X_T(i).GT.Ub(i)) ) THEN
            CALL RANDOM_NUMBER(RND1)
            X_T(i) = Lb(i) + (Ub(i)-Lb(i))*RND1
        ENDIF
      ENDDO
 
 
      END

 
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * |
      FUNCTION DGAMMA(R)
!  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * |
!     This routine evaluate usual Gamma function.
 
      IMPLICIT NONE


      REAL*8 DGAMMA , P0 , P1 , P10 , P11 , P12 , P13 , P2 , P3 , P4 ,  &
           & P5 , P6 , P7 , P8 , P9 , R , w , y
      INTEGER k , n

 
      PARAMETER (P0=0.999999999999999990D+00,                           &
               & P1=-0.422784335098466784D+00,                          &
               & P2=-0.233093736421782878D+00,                          &
               & P3=0.191091101387638410D+00,                           &
               & P4=-0.024552490005641278D+00,                          &
               & P5=-0.017645244547851414D+00,                          &
               & P6=0.008023273027855346D+00)
      PARAMETER (P7=-0.000804329819255744D+00,                          &
               & P8=-0.000360837876648255D+00,                          &
               & P9=0.000145596568617526D+00,                           &
               & P10=-0.000017545539395205D+00,                         &
               & P11=-0.000002591225267689D+00,                         &
               & P12=0.000001337767384067D+00,                          &
               & P13=-0.000000199542863674D+00)
 
      n = NINT(R-2)
      w = R - (n+2)
      y = ((((((((((((P13*w+P12)*w+P11)*w+P10)*w+P9)*w+P8)*w+P7)*w+P6)*w&
        & +P5)*w+P4)*w+P3)*w+P2)*w+P1)*w + P0
      IF ( n.GT.0 ) THEN
         w = R - 1
         DO k = 2 , n
            w = w*(R-k)
         ENDDO
      ELSE
         w = 1
         DO k = 0 , -n - 1
            y = y*(R+k)
         ENDDO
      ENDIF
      DGAMMA = w/y
 
      END
    
 
      
END SUBROUTINE GSA

 


 
 


SUBROUTINE SA(FCN,T,Maxevl,Nfcnev,Ier,Xopt,Fopt,Vm,&
            & XSTART,Lb1,Ub1,STEPFACTOR,Rt1,Eps1,Ns1,Nt1,Neps1,Nacc1,Nobds1,Iprint1,IsMax)
IMPLICIT NONE

 
!  Version: 3.2
!  Date: 1/22/94.
!  Differences compared to Version 2.0:
!     1. If a trial is out of bounds, a point is randomly selected
!        from LB(i) to UB(i). Unlike in version 2.0, this trial is
!        evaluated and is counted in acceptances and rejections.
!        All corresponding documentation was changed as well.
!  Differences compared to Version 3.0:
!     1. If VM(i) > (UB(i) - LB(i)), VM is set to UB(i) - LB(i).
!        The idea is that if T is high relative to LB & UB, most
!        points will be accepted, causing VM to rise. But, in this
!        situation, VM has little meaning; particularly if VM is
!        larger than the acceptable region. Setting VM to this size
!        still allows all parts of the allowable region to be selected.
!  Differences compared to Version 3.1:
!     1. Test made to see if the initial temperature is positive.
!     2. WRITE statements prettied up.
!     3. References to paper updated.
!
!  Synopsis:
!  This routine implements the continuous simulated annealing global
!  optimization algorithm described in Corana et al.'s article
!  "Minimizing Multimodal Functions of Continuous Variables with the
!  "Simulated Annealing" Algorithm" in the September 1987 (vol. 13,
!  no. 3, pp. 262-280) issue of the ACM Transactions on Mathematical
!  Software.
!
!  A very quick (perhaps too quick) overview of SA:
!     SA tries to find the global optimum of an N dimensional function.
!  It moves both up and downhill and as the optimization process
!  proceeds, it focuses on the most promising area.
!     To start, it randomly chooses a trial point within the step length
!  VM (a vector of length N) of the user selected starting point. The
!  function is evaluated at this trial point and its value is compared
!  to its value at the initial point.
!     In a maximization problem, all uphill moves are accepted and the
!  algorithm continues from that trial point. Downhill moves may be
!  accepted; the decision is made by the Metropolis criteria. It uses T
!  (temperature) and the size of the downhill move in a probabilistic
!  manner. The smaller T and the size of the downhill move are, the more
!  likely that move will be accepted. If the trial is accepted, the
!  algorithm moves on from that point. If it is rejected, another point
!  is chosen instead for a trial evaluation.
!     Each element of VM periodically adjusted so that half of all
!  function evaluations in that direction are accepted.
!     A fall in T is imposed upon the system with the RT variable by
!  T(i+1) = RT*T(i) where i is the ith iteration. Thus, as T declines,
!  downhill moves are less likely to be accepted and the percentage of
!  rejections rise. Given the scheme for the selection for VM, VM falls.
!  Thus, as T declines, VM falls and SA focuses upon the most promising
!  area for optimization.
!
!  The importance of the parameter T:
!     The parameter T is crucial in using SA successfully. It influences
!  VM, the step length over which the algorithm searches for optima. For
!  a small intial T, the step length may be too small; thus not enough
!  of the function might be evaluated to find the global optima. The user
!  should carefully examine VM in the intermediate output (set IPRINT =
!  1) to make sure that VM is appropriate. The relationship between the
!  initial temperature and the resulting step length is function
!  dependent.
!     To determine the starting temperature that is consistent with
!  optimizing a function, it is worthwhile to run a trial run first. Set
!  RT = 1.5 and T = 1.0. With RT > 1.0, the temperature increases and VM
!  rises as well. Then select the T that produces a large enough VM.
!
!  For modifications to the algorithm and many details on its use,
!  (particularly for econometric applications) see Goffe, Ferrier
!  and Rogers, "Global Optimization of Statistical Functions with
!  Simulated Annealing," Journal of Econometrics, vol. 60, no. 1/2,
!  Jan./Feb. 1994, pp. 65-100.
!  For more information, contact
!              Bill Goffe
!              Department of Economics and International Business
!              University of Southern Mississippi
!              Hattiesburg, MS  39506-5072
!              (601) 266-4484 (office)
!              (601) 266-4920 (fax)
!              bgoffe@whale.st.usm.edu (Internet)
!
!  As far as possible, the parameters here have the same name as in
!  the description of the algorithm on pp. 266-8 of Corana et al.
!
!  In this description, SP is single precision, DP is REAL(8),
!  INT is integer, L is logical and (N) denotes an array of length n.
!  Thus, DP(N) denotes a REAL(8) array of length n.
!
!  Input Parameters:
!    Note: The suggested values generally come from Corana et al. To
!          drastically reduce runtime, see Goffe et al., pp. 90-1 for
!          suggestions on choosing the appropriate RT and NT.
!    N - Number of variables in the function to be optimized. (INT)
!    X - The starting values for the variables of the function to be
!        optimized. (DP(N))
!    IsMax - Denotes whether the function should be maximized or
!          minimized. A true value denotes maximization while a false
!          value denotes minimization. Intermediate output (see IPRINT)
!          takes this into account. (L)
!    RT - The temperature reduction factor. The value suggested by
!         Corana et al. is .85. See Goffe et al. for more advice. (DP)
!    EPS - Error tolerance for termination. If the final function
!          values from the last neps temperatures differ from the
!          corresponding value at the current temperature by less than
!          EPS and the final function value at the current temperature
!          differs from the current optimal function value by less than
!          EPS, execution terminates and IER = 0 is returned. (EP)
!    NS - Number of cycles. After NS*N function evaluations, each
!         element of VM is adjusted so that approximately half of
!         all function evaluations are accepted. The suggested value
!         is 20. (INT)
!    NT - Number of iterations before temperature reduction. After
!         NT*NS*N function evaluations, temperature (T) is changed
!         by the factor RT. Value suggested by Corana et al. is
!         MAX(100, 5*N). See Goffe et al. for further advice. (INT)
!    NEPS - Number of final function values used to decide upon termi-
!           nation. See EPS. Suggested value is 4. (INT)
!    MAXEVL - The maximum number of function evaluations. If it is
!             exceeded, IER = 1. (INT)
!    LB - The lower bound for the allowable solution variables. (DP(N))
!    UB - The upper bound for the allowable solution variables. (DP(N))
!         If the algorithm chooses X(I) .LT. LB(I) or X(I) .GT. UB(I),
!         I = 1, N, a point is from inside is randomly selected. This
!         This focuses the algorithm on the region inside UB and LB.
!         Unless the user wishes to concentrate the search to a par-
!         ticular region, UB and LB should be set to very large positive
!         and negative values, respectively. Note that the starting
!         vector X should be inside this region. Also note that LB and
!         UB are fixed in position, while VM is centered on the last
!         accepted trial set of variables that optimizes the function.
!    STEPFACTOR - Vector that controls the step length adjustment. The suggested
!        value for all elements is 2.0. (DP(N))
!    IPRINT - controls printing inside SA. (INT)
!             Values: 0 - Nothing printed.
!                     1 - Function value for the starting value and
!                         summary results before each temperature
!                         reduction. This includes the optimal
!                         function value found so far, the total
!                         number of moves (broken up into uphill,
!                         downhill, accepted and rejected), the
!                         number of out of bounds trials, the
!                         number of new optima found at this
!                         temperature, the current optimal X and
!                         the step length VM. Note that there are
!                         N*NS*NT function evalutations before each
!                         temperature reduction. Finally, notice is
!                         is also given upon achieveing the termination
!                         criteria.
!                     2 - Each new step length (VM), the current optimal
!                         X (XOPT) and the current trial X (X). This
!                         gives the user some idea about how far X
!                         strays from XOPT as well as how VM is adapting
!                         to the function.
!                     3 - Each function evaluation, its acceptance or
!                         rejection and new optima. For many problems,
!                         this option will likely require a small tree
!                         if hard copy is used. This option is best
!                         used to learn about the algorithm. A small
!                         value for MAXEVL is thus recommended when
!                         using IPRINT = 3.
!             Suggested value: 1
!             Note: For a given value of IPRINT, the lower valued
!                   options (other than 0) are utilized.
!  Input/Output Parameters:
!    T - On input, the initial temperature. See Goffe et al. for advice.
!        On output, the final temperature. (DP)
!    VM - The step length vector. On input it should encompass the
!         region of interest given the starting value X. For point
!         X(I), the next trial point is selected is from X(I) - VM(I)
!         to  X(I) + VM(I). Since VM is adjusted so that about half
!         of all points are accepted, the input value is not very
!         important (i.e. is the value is off, SA adjusts VM to the
!         correct value). (DP(N))
!
!  Output Parameters:
!    XOPT - The variables that optimize the function. (DP(N))
!    FOPT - The optimal value of the function. (DP)
!    NACC - The number of accepted function evaluations. (INT)
!    NFCNEV - The total number of function evaluations. In a minor
!             point, note that the first evaluation is not used in the
!             core of the algorithm; it simply initializes the
!             algorithm. (INT).
!    NOBDS - The total number of trial function evaluations that
!            would have been out of bounds of LB and UB. Note that
!            a trial point is randomly selected between LB and UB.
!            (INT)
!    IER - The error return number. (INT)
!          Values: 0 - Normal return; termination criteria achieved.
!                  1 - Number of function evaluations (NFCNEV) is
!                      greater than the maximum number (MAXEVL).
!                  2 - The starting value (X) is not inside the
!                      bounds (LB and UB).
!                  3 - The initial temperature is not positive.
!                  99 - Should not be seen; only used internally.
!

 
!  Type all external variables.
REAL(8),DIMENSION(:),INTENT(IN),OPTIONAL:: XSTART,STEPFACTOR,LB1,UB1
REAL(8),DIMENSION(:),INTENT(OUT)::Xopt
REAL(8),DIMENSION(:),INTENT(INOUT)::Vm  
REAL(8),INTENT(in),OPTIONAL:: Eps1,Rt1 
REAL(8),INTENT(inout):: T
REAL(8),INTENT(OUT)::Fopt      
INTEGER,INTENT(IN),OPTIONAL::Ns1,Nt1,Neps1
INTEGER,INTENT(IN),TARGET,OPTIONAL::Nobds1,Nacc1
INTEGER,INTENT(IN)::Maxevl
INTEGER,INTENT(OUT)::Ier,Nfcnev
INTEGER,INTENT(IN),OPTIONAL::Iprint1 
LOGICAL,OPTIONAL::IsMax 
      
procedure(optifunc) :: FCN
   
!  Type all internal variables.
REAL(8)::Xp(SIZE(XOPT)),X(SIZE(XOPT)),C(SIZE(XOPT)),LB(SIZE(XOPT)),UB(SIZE(XOPT))
REAL(8),DIMENSION(:),ALLOCATABLE::Fstar
INTEGER::Nacp(SIZE(XOPT))
      
REAL(8)::EPS,RND1,RT,f , fp , p , pp , ratio
INTEGER,TARGET::N,IPRINT,Ns,Nt,Neps,nup , ndown , nrej , nnew , lnobds , h , i , j , m  
INTEGER,TARGET::NOBDS2,NACC2
INTEGER,POINTER::Nobds,Nacc
LOGICAL quit,MAX1
 
      
CALL SA_INIT()
 

!  If the initial temperature is not positive, notify the user and
!  return to the calling routine.
IF ( T.LE.0.0 ) THEN
    WRITE (*,                                                      &
&'(/,''  THE INITIAL TEMPERATURE IS NOT POSITIVE. ''               &
&  /,''  RESET THE VARIABLE T. ''/)')
    Ier = 3
    RETURN
ENDIF
 
!  If the initial value is out of bounds, notify the user and return
!  to the calling routine.
DO i = 1 , N
    IF ( (X(i).GT.Ub(i)) .OR. (X(i).LT.Lb(i)) ) THEN
    CALL PRT1
    Ier = 2
    RETURN
    ENDIF
ENDDO
 
!  Evaluate the function with input X and return value as F.
!CALL FCN(N,X,f)
F=FCN(X)
!  If the function is to be minimized, switch the sign of the function.
!  Note that all intermediate and final output switches the sign back
!  to eliminate any possible confusion for the user.
IF ( .NOT.MAX1 ) f = -f
Nfcnev = Nfcnev + 1
Fopt = f
Fstar(1) = f
IF ( Iprint.GE.1 ) CALL PRT2(MAX1,N,X,f)
 
!  Start the main loop. Note that it terminates if (i) the algorithm
!  succesfully optimizes the function or (ii) there are too many
!  function evaluations (more than MAXEVL).
DO WHILE (Nfcnev<=Maxevl)
      
    nup = 0
    nrej = 0
    nnew = 0
    ndown = 0
    lnobds = 0
 
    DO m = 1 , Nt
        DO j = 1 , Ns
            DO h = 1 , N
                !Generate XP, the trial value of X. Note use of VM to choose XP.
                DO i = 1 , N
                    IF ( i.EQ.h ) THEN
                        CALL RANDOM_NUMBER(RND1)
                        Xp(i) = X(i) + (RND1*2.-1.)*Vm(i)
                    ELSE
                        Xp(i) = X(i)
                    ENDIF
 
                    !  If XP is out of bounds, select a point in bounds for the trial.
                    IF ( (Xp(i).LT.Lb(i)) .OR. (Xp(i).GT.Ub(i)) ) THEN
                        CALL RANDOM_NUMBER(RND1)
                        Xp(i) = Lb(i) + (Ub(i)-Lb(i))*RND1
                        lnobds = lnobds + 1
                        Nobds = Nobds + 1
                        IF ( Iprint.GE.3 ) CALL PRT3(Max1,N,Xp,X,fp,f)
                    ENDIF
                ENDDO
 
                !  Evaluate the function with the trial point XP and return as FP.
        
                FP=FCN(XP)
                IF ( .NOT.Max1 ) fp = -fp
                Nfcnev = Nfcnev + 1
                IF ( Iprint.GE.3 ) CALL PRT4(MAX1,N,Xp,X,fp,f)
 
                !  If too many function evaluations occur, terminate the algorithm.
                IF ( Nfcnev.GE.Maxevl ) THEN
                    CALL PRT5
                    IF ( .NOT.MAX1 ) Fopt = -Fopt
                    Ier = 1
                    RETURN
                ENDIF
 
                !  Accept the new point if the function value increases.
                IF ( fp.GE.f ) THEN
                    IF ( Iprint.GE.3 ) WRITE (*,'(''  POINT ACCEPTED'')')

                    X = Xp
                    f = fp
                    Nacc = Nacc + 1
                    Nacp(h) = Nacp(h) + 1
                    nup = nup + 1
 
                    !  If greater than any other point, record as new optimum.
                    IF ( fp.GT.Fopt ) THEN
                        IF ( Iprint.GE.3 ) WRITE (*,'(''  NEW OPTIMUM'')')
                        Xopt = Xp
                        Fopt = fp
                        nnew = nnew + 1
                    ENDIF
 
                !  If the point is lower, use the Metropolis criteria to decide on
                !  acceptance or rejection.
                ELSE
                    p = EXP((fp-f)/T)
                    CALL RANDOM_NUMBER(PP)
                    !pp = RANMAR()
                    IF ( pp.LT.p ) THEN
                        IF ( Iprint.GE.3 ) CALL PRT6(MAX1)

                        X = Xp
                        f = fp
                        Nacc = Nacc + 1
                        Nacp(h) = Nacp(h) + 1
                        ndown = ndown + 1
                    ELSE
                        nrej = nrej + 1
                        IF ( Iprint.GE.3 ) CALL PRT7(MAX1)
                    ENDIF
                ENDIF
 
            ENDDO
        ENDDO
 
        !  Adjust VM so that approximately half of all evaluations are accepted.
        DO i = 1 , N
            ratio = DFLOAT(Nacp(i))/DFLOAT(Ns)
            IF ( ratio.GT..6 ) THEN
                Vm(i) = Vm(i)*(1.+C(i)*(ratio-.6)/.4)
            ELSEIF ( ratio.LT..4 ) THEN
                Vm(i) = Vm(i)/(1.+C(i)*((.4-ratio)/.4))
            ENDIF
            IF ( Vm(i).GT.(Ub(i)-Lb(i)) ) Vm(i) = Ub(i) - Lb(i)
        ENDDO
 
        IF ( Iprint.GE.2 ) CALL PRT8(N,Vm,Xopt,X)
 
         
        Nacp = 0
         
 
    ENDDO
 
    IF ( Iprint.GE.1 ) CALL PRT9(MAX1,N,T,Xopt,Vm,Fopt,nup,ndown,nrej, &
                            & lnobds,nnew)
 
    !  Check termination criteria.
    quit = .FALSE.
    Fstar(1) = f
    IF ( (Fopt-Fstar(1)).LE.Eps ) quit = .TRUE.
    DO i = 1 , Neps
        IF ( ABS(f-Fstar(i)).GT.Eps ) quit = .FALSE.
    ENDDO
 
    !  Terminate SA if appropriate.
    IF ( quit ) THEN
        X = Xopt
        Ier = 0
        IF ( .NOT.MAX1 ) Fopt = -Fopt
        IF ( Iprint.GE.1 ) CALL PRT10
        RETURN
    ENDIF
 
    !  If termination criteria is not met, prepare for another loop.
    T = Rt*T
    DO i = Neps , 2 , -1
        Fstar(i) = Fstar(i-1)
    ENDDO
    
    f = Fopt      
    X = Xopt
      
 

END DO
      
CONTAINS


SUBROUTINE  SA_INIT()
    
    CALL RANDOM_SEED()
        
    N=SIZE(X,DIM=1)
    
    IF(PRESENT(ISMAX)) THEN
        MAX1=ISMAX
    ELSE
        MAX1=.FALSE.
    ENDIF
    IF(PRESENT(NEPS1)) THEN
        NEPS=NEPS1
    ELSE
        NEPS=4
    ENDIF         
    IF(ALLOCATED(FSTAR)) DEALLOCATE(FSTAR)
    ALLOCATE(FSTAR(NEPS))
    Fstar = 1.0D+20
        
    IF(PRESENT(NS1)) THEN
        NS=NS1
    ELSE
        NS=20
    ENDIF 
        
    IF(PRESENT(NT1)) THEN
        NT=NT1
    ELSE
        NT=MAX(100,5*N)
    ENDIF 
    IF(PRESENT(EPS1)) THEN
        EPS=EPS1
    ELSE
        EPS=1.D-4
    ENDIF 

    IF(PRESENT(RT1)) THEN
        RT=RT1
    ELSE
        RT=0.9D0
    ENDIF 

    IF(PRESENT(NOBDS1)) THEN
        NOBDS=>NOBDS1
    ELSE
        NOBDS=>NOBDS2
    ENDIF 

    IF(PRESENT(NACC1)) THEN
        NACC=>NACC1
    ELSE
        NACC=>NACC2
    ENDIF 
 
    IF(PRESENT(LB1)) THEN
        LB=LB1
    ELSE
        LB=-1.D10
    ENDIF    
    IF(PRESENT(UB1)) THEN
        UB=UB1
    ELSE
        UB=1.D10
    ENDIF
        
    IF(PRESENT(IPRINT1)) THEN
        IPRINT=IPRINT1
    ELSE
        IPRINT=0
    ENDIF
    
    IF(PRESENT(STEPFACTOR)) THEN
        C=STEPFACTOR
    ELSE
        C=2.0D0
    ENDIF  
 
    IF(PRESENT(XSTART)) THEN
        X=XSTART
    ELSE
        CALL RANDOM_NUMBER(X)
        X=LB+X*(UB-LB)
    ENDIF
    
    !  Set initial values.
    Nacc = 0
    Nobds = 0
    Nfcnev = 0
    Ier = 99
 
    Xopt = X
    Nacp = 0
        
        
    
END SUBROUTINE
 
END SUBROUTINE

 
SUBROUTINE PRT1
IMPLICIT NONE
!*--PRT1105
!  This subroutine prints intermediate output, as does PRT2 through
!  PRT10. Note that if SA is minimizing the function, the sign of the
!  function value and the directions (up/down) are reversed in all
!  output to correspond with the actual function optimization. This
!  correction is because SA was written to maximize functions and
!  it minimizes by maximizing the negative a function.
 
WRITE (*,                                                         &
&'(/,''  THE STARTING VALUE (X) IS OUTSIDE THE BOUNDS ''           &
&  /,''  (LB AND UB). EXECUTION TERMINATED WITHOUT ANY''           &
&  /,''  OPTIMIZATION. RESPECIFY X, UB OR LB SO THAT  ''           &
&  /,''  LB(I) .LT. X(I) .LT. UB(I), I = 1, N. ''/)')
 
END
!*==PRT2.spg  processed by SPAG 6.72Dc at 05:07 on 29 Sep 2018
 
SUBROUTINE PRT2(MAX,N,X,F)
IMPLICIT NONE
!*--PRT2124
 
REAL(8) X(*) , F
INTEGER N
LOGICAL Max
 
WRITE (*,'(''  '')')
CALL PRTVEC(X,N,'INITIAL X')
IF ( Max ) THEN
    WRITE (*,'(''  INITIAL F: '',/, G25.18)') F
ELSE
    WRITE (*,'(''  INITIAL F: '',/, G25.18)') -F
ENDIF
 
END
!*==PRT3.spg  processed by SPAG 6.72Dc at 05:07 on 29 Sep 2018
 
SUBROUTINE PRT3(Max,N,Xp,X,Fp,F)
IMPLICIT NONE
!*--PRT3143
 
REAL(8) Xp(*) , X(*) , Fp , F
INTEGER N
LOGICAL Max
 
WRITE (*,'(''  '')')
CALL PRTVEC(X,N,'CURRENT X')
IF ( Max ) THEN
    WRITE (*,'(''  CURRENT F: '',G25.18)') F
ELSE
    WRITE (*,'(''  CURRENT F: '',G25.18)') -F
ENDIF
CALL PRTVEC(Xp,N,'TRIAL X')
WRITE (*,'(''  POINT REJECTED SINCE OUT OF BOUNDS'')')
 
END
!*==PRT4.spg  processed by SPAG 6.72Dc at 05:07 on 29 Sep 2018
 
SUBROUTINE PRT4(Max,N,Xp,X,Fp,F)
IMPLICIT NONE
!*--PRT4164
 
REAL(8) Xp(*) , X(*) , Fp , F
INTEGER N
LOGICAL Max
 
WRITE (*,'(''  '')')
CALL PRTVEC(X,N,'CURRENT X')
IF ( Max ) THEN
    WRITE (*,'(''  CURRENT F: '',G25.18)') F
    CALL PRTVEC(Xp,N,'TRIAL X')
    WRITE (*,'(''  RESULTING F: '',G25.18)') Fp
ELSE
    WRITE (*,'(''  CURRENT F: '',G25.18)') -F
    CALL PRTVEC(Xp,N,'TRIAL X')
    WRITE (*,'(''  RESULTING F: '',G25.18)') -Fp
ENDIF
 
END
!*==PRT5.spg  processed by SPAG 6.72Dc at 05:07 on 29 Sep 2018
 
SUBROUTINE PRT5
IMPLICIT NONE
!*--PRT5187
 
WRITE (*,                                                         &
&'(/,''  TOO MANY FUNCTION EVALUATIONS; CONSIDER ''                &
&  /,''  INCREASING MAXEVL OR EPS, OR DECREASING ''                &
&  /,''  NT OR RT. THESE RESULTS ARE LIKELY TO BE ''               &
&  /,''  POOR.'',/)')
 
END
!*==PRT6.spg  processed by SPAG 6.72Dc at 05:07 on 29 Sep 2018
 
SUBROUTINE PRT6(Max)
IMPLICIT NONE
!*--PRT6200
 
LOGICAL Max
 
IF ( Max ) THEN
    WRITE (*,'(''  THOUGH LOWER, POINT ACCEPTED'')')
ELSE
    WRITE (*,'(''  THOUGH HIGHER, POINT ACCEPTED'')')
ENDIF
 
END
!*==PRT7.spg  processed by SPAG 6.72Dc at 05:07 on 29 Sep 2018
 
SUBROUTINE PRT7(Max)
IMPLICIT NONE
!*--PRT7215
 
LOGICAL Max
 
IF ( Max ) THEN
    WRITE (*,'(''  LOWER POINT REJECTED'')')
ELSE
    WRITE (*,'(''  HIGHER POINT REJECTED'')')
ENDIF
 
END
!*==PRT8.spg  processed by SPAG 6.72Dc at 05:07 on 29 Sep 2018
 
SUBROUTINE PRT8(N,Vm,Xopt,X)
IMPLICIT NONE
!*--PRT8230
 
REAL(8) Vm(*) , Xopt(*) , X(*)
INTEGER N
 
WRITE (*,                                                         &
&'(/,                                                        '' INT&
&ERMEDIATE RESULTS AFTER STEP LENGTH ADJUSTMENT'',/)')
CALL PRTVEC(Vm,N,'NEW STEP LENGTH (VM)')
CALL PRTVEC(Xopt,N,'CURRENT OPTIMAL X')
CALL PRTVEC(X,N,'CURRENT X')
WRITE (*,'('' '')')
 
END
!*==PRT9.spg  processed by SPAG 6.72Dc at 05:07 on 29 Sep 2018
 
SUBROUTINE PRT9(Max,N,T,Xopt,Vm,Fopt,Nup,Ndown,Nrej,Lnobds,Nnew)
IMPLICIT NONE
!*--PRT9248
 
REAL(8) Xopt(*) , Vm(*) , T , Fopt
INTEGER N , Nup , Ndown , Nrej , Lnobds , Nnew , totmov
LOGICAL Max
 
totmov = Nup + Ndown + Nrej
 
WRITE (*,                                                         &
&'(/,                                                        '' INT&
&ERMEDIATE RESULTS BEFORE NEXT TEMPERATURE REDUCTION'',/)')
WRITE (*,'(''  CURRENT TEMPERATURE:            '',G12.5)') T
IF ( Max ) THEN
    WRITE (*,'(''  MAX FUNCTION VALUE SO FAR:  '',G25.18)') Fopt
    WRITE (*,'(''  TOTAL MOVES:                '',I8)') totmov
    WRITE (*,'(''     UPHILL:                  '',I8)') Nup
    WRITE (*,'(''     ACCEPTED DOWNHILL:       '',I8)') Ndown
    WRITE (*,'(''     REJECTED DOWNHILL:       '',I8)') Nrej
    WRITE (*,'(''  OUT OF BOUNDS TRIALS:       '',I8)') Lnobds
    WRITE (*,'(''  NEW MAXIMA THIS TEMPERATURE:'',I8)') Nnew
ELSE
    WRITE (*,'(''  MIN FUNCTION VALUE SO FAR:  '',G25.18)') -Fopt
    WRITE (*,'(''  TOTAL MOVES:                '',I8)') totmov
    WRITE (*,'(''     DOWNHILL:                '',I8)') Nup
    WRITE (*,'(''     ACCEPTED UPHILL:         '',I8)') Ndown
    WRITE (*,'(''     REJECTED UPHILL:         '',I8)') Nrej
    WRITE (*,'(''  TRIALS OUT OF BOUNDS:       '',I8)') Lnobds
    WRITE (*,'(''  NEW MINIMA THIS TEMPERATURE:'',I8)') Nnew
ENDIF
CALL PRTVEC(Xopt,N,'CURRENT OPTIMAL X')
CALL PRTVEC(Vm,N,'STEP LENGTH (VM)')
WRITE (*,'('' '')')
 
END
!*==PRT10.spg  processed by SPAG 6.72Dc at 05:07 on 29 Sep 2018
 
SUBROUTINE PRT10
IMPLICIT NONE
!*--PRT10286
 
WRITE (*,'(/,''  SA ACHIEVED TERMINATION CRITERIA. IER = 0. '',/)'&
    & )
 
END
!*==PRTVEC.spg  processed by SPAG 6.72Dc at 05:07 on 29 Sep 2018
 
SUBROUTINE PRTVEC(Vector,Ncols,Name)
IMPLICIT NONE
!*--PRTVEC296
!*** Start of declarations inserted by SPAG
INTEGER i , j , lines , ll
!*** End of declarations inserted by SPAG
!  This subroutine prints the REAL(8) vector named VECTOR.
!  Elements 1 thru NCOLS will be printed. NAME is a character variable
!  that describes VECTOR. Note that if NAME is given in the call to
!  PRTVEC, it must be enclosed in quotes. If there are more than 10
!  elements in VECTOR, 10 elements will be printed on each line.
 
INTEGER Ncols
REAL(8) Vector(Ncols)
CHARACTER*(*) Name
 
WRITE (*,99001) Name
99001 FORMAT (/,25X,A)
 
IF ( Ncols.GT.10 ) THEN
    lines = INT(Ncols/10.)
 
    DO i = 1 , lines
    ll = 10*(i-1)
    WRITE (*,99002) (Vector(j),j=1+ll,10+ll)
    ENDDO
 
    WRITE (*,99002) (Vector(j),j=11+ll,Ncols)
ELSE
    WRITE (*,99002) (Vector(j),j=1,Ncols)
ENDIF
 
99002 FORMAT (10(G12.5,1X))
 
END
      
END MODULE



!PROGRAM SIMANN_TRY
!    USE SIMANN
!      IMPLICIT NONE
!!*--SIMANN28
!!*** Start of declarations inserted by SPAG
!      INTEGER i , N , NEPS
!!*** End of declarations inserted by SPAG
!!  This file is an example of the Corana et al. simulated annealing
!!  algorithm for multimodal and robust optimization as implemented
!!  and modified by Goffe, Ferrier and Rogers. Counting the above line
!!  ABSTRACT as 1, the routine itself (SA), with its supplementary
!!  routines, is on lines 232-990. A multimodal example from Judge et al.
!!  (FCN) is on lines 150-231. The rest of this file (lines 1-149) is a
!!  driver routine with values appropriate for the Judge example. Thus, this
!!  example is ready to run.
!!
!!  To understand the algorithm, the documentation for SA on lines 236-
!!  484 should be read along with the parts of the paper that describe
!!  simulated annealing. Then the following lines will then aid the user
!!  in becomming proficient with this implementation of simulated
!!  annealing.
!!
!!  Learning to use SA:
!!      Use the sample function from Judge with the following suggestions
!!  to get a feel for how SA works. When you've done this, you should be
!!  ready to use it on most any function with a fair amount of expertise.
!!    1. Run the program as is to make sure it runs okay. Take a look at
!!       the intermediate output and see how it optimizes as temperature
!!       (T) falls. Notice how the optimal point is reached and how
!!       falling T reduces VM.
!!    2. Look through the documentation to SA so the following makes a
!!       bit of sense. In line with the paper, it shouldn't be that hard
!!       to figure out. The core of the algorithm is described on pp. 68-70
!!       and on pp. 94-95. Also see Corana et al. pp. 264-9.
!!    3. To see how it selects points and makes decisions about uphill
!!       and downhill moves, set IPRINT = 3 (very detailed intermediate
!!       output) and MAXEVL = 100 (only 100 function evaluations to limit
!!       output).
!!    4. To see the importance of different temperatures, try starting
!!       with a very low one (say T = 10E-5). You'll see (i) it never
!!       escapes from the local optima (in annealing terminology, it
!!       quenches) & (ii) the step length (VM) will be quite small. This
!!       is a key part of the algorithm: as temperature (T) falls, step
!!       length does too. In a minor point here, note how VM is quickly
!!       reset from its initial value. Thus, the input VM is not very
!!       important. This is all the more reason to examine VM once the
!!       algorithm is underway.
!!    5. To see the effect of different parameters and their effect on
!!       the speed of the algorithm, try RT = .95 & RT = .1. Notice the
!!       vastly different speed for optimization. Also try NT = 20. Note
!!       that this sample function is quite easy to optimize, so it will
!!       tolerate big changes in these parameters. RT and NT are the
!!       parameters one should adjust to modify the runtime of the
!!       algorithm and its robustness.
!!    6. Try constraining the algorithm with either LB or UB.
! 
!      PARAMETER (N=2,NEPS=4)
! 
!      REAL(8) lb(N) , ub(N) , x(N) , xopt(N) , c(N) , vm(N) ,  &
!                     & fstar(NEPS) , xp(N) , t , eps , rt , fopt
! 
!      INTEGER nacp(N) , ns , nt , nfcnev , ier , iseed1 , iseed2 ,      &
!            & maxevl , iprint , nacc , nobds
! 
!      LOGICAL max
! 
!      PROCEDURE(optifunc)::FCN
! 
!!  Set underflows to zero on IBM mainframes.
!!     CALL XUFLOW(0)
! 
!!  Set input parameters.
!      max = .FALSE.
!      eps = 1.0D-6
!      rt = .5
!      iseed1 = 1
!      iseed2 = 2
!      ns = 20
!      nt = 5
!      maxevl = 1000000
!      iprint = 1
!      DO i = 1 , N
!         lb(i) = -1.0D10
!         ub(i) = 1.0D10
!         c(i) = 2.0
!      ENDDO
! 
!!  Note start at local, but not global, optima of the Judge function.
!      !x(1) = 2.354471
!      !x(2) = -0.319186
!      call RANDOM_SEED()
!      call RANDOM_NUMBER(x)
!      !x=lb+(ub-lb)*x
!      x=1.e6
!!  Set input values of the input/output parameters.
!      t = 5.0
!      DO i = 1 , N
!         vm(i) = 1.0
!      ENDDO
!      vm=1.d10
!      WRITE (*,99001) N , max , t , rt , eps , ns , nt , NEPS , maxevl ,&
!                    & iprint , iseed1 , iseed2
! 
!99001 FORMAT (/,' SIMULATED ANNEALING EXAMPLE',/,/,                     &
!             &' NUMBER OF PARAMETERS: ',I3,'   MAXIMAZATION: ',L5,/,    &
!             &' INITIAL TEMP: ',G8.2,'   RT: ',G8.2,'   EPS: ',G8.2,/,  &
!             &' NS: ',I3,'   NT: ',I2,'   NEPS: ',I2,/,' MAXEVL: ',I10, &
!             &'   IPRINT: ',I1,'   ISEED1: ',I4,'   ISEED2: ',I4)
! 
!      CALL PRTVEC(x,N,'STARTING VALUES')
!      CALL PRTVEC(vm,N,'INITIAL STEP LENGTH')
!      CALL PRTVEC(lb,N,'LOWER BOUND')
!      CALL PRTVEC(ub,N,'UPPER BOUND')
!      CALL PRTVEC(c,N,'C VECTOR')
!      WRITE (*,                                                         &
!     &'(/,''  ****   END OF DRIVER ROUTINE OUTPUT   ****''              &
!     &  /,''  ****   BEFORE CALL TO SA.             ****'')')
! 
!      !CALL SA(N,x,max,rt,eps,ns,nt,NEPS,maxevl,lb,ub,c,iprint,iseed1,   &
!      !      & iseed2,t,vm,xopt,fopt,nacc,nfcnev,nobds,ier,fstar,xp,nacp)
!       !CALL SA(FCN,T,Maxevl,Nfcnev,Ier,Xopt,Fopt,Vm , &
!       !     & X ,LB,Ub ,C,Rt,Eps,Ns,Nt,Neps,Nacc,Nobds,Iprint,MAX)
!            
!       CALL GSA(FCN,T,Fopt,Xopt,QAI=1.5D0,QTI=1.5D0,QVI=1.5D0,EPSI=EPS,MaxevlI=MAXEVL,XSTART=X,LbI=LB,UbI=UB,IsSGSAI=.FALSE.,PITERI=100)
! 
!      WRITE (*,'(/,''  ****   RESULTS AFTER SA   ****   '')')
!      CALL PRTVEC(xopt,N,'SOLUTION')
!      !CALL PRTVEC(vm,N,'FINAL STEP LENGTH')
!      WRITE (*,99002) fopt , nfcnev , nacc , nobds , t , ier
!99002 FORMAT (/,' OPTIMAL FUNCTION VALUE: ',G20.13/,                    &
!             &' NUMBER OF FUNCTION EVALUATIONS:     ',I10,/,            &
!             &' NUMBER OF ACCEPTED EVALUATIONS:     ',I10,/,            &
!             &' NUMBER OF OUT OF BOUND EVALUATIONS: ',I10,/,            &
!             &' FINAL TEMP: ',G20.13,'  IER: ',I3)
!      PAUSE
! 
!      END
!
! 
!REAl(8) FUNCTION FCN(Theta)
!      IMPLICIT NONE
!      REAL(8),DIMENSION(:),INTENT(IN)::Theta
!        
!      INTEGER i , N
!
!!  This subroutine is from the example in Judge et al., The Theory and
!!  Practice of Econometrics, 2nd ed., pp. 956-7. There are two optima:
!!  F(.864,1.23) = 16.0817 (the global minumum) and F(2.35,-.319) = 20.9805.
! 
!      
! 
!      REAL(8) y(20) , x2(20) , x3(20) 
!      
!      N=SIZE(Theta)
!      
!      y(1) = 4.284
!      y(2) = 4.149
!      y(3) = 3.877
!      y(4) = 0.533
!      y(5) = 2.211
!      y(6) = 2.389
!      y(7) = 2.145
!      y(8) = 3.231
!      y(9) = 1.998
!      y(10) = 1.379
!      y(11) = 2.106
!      y(12) = 1.428
!      y(13) = 1.011
!      y(14) = 2.179
!      y(15) = 2.858
!      y(16) = 1.388
!      y(17) = 1.651
!      y(18) = 1.593
!      y(19) = 1.046
!      y(20) = 2.152
! 
!      x2(1) = .286
!      x2(2) = .973
!      x2(3) = .384
!      x2(4) = .276
!      x2(5) = .973
!      x2(6) = .543
!      x2(7) = .957
!      x2(8) = .948
!      x2(9) = .543
!      x2(10) = .797
!      x2(11) = .936
!      x2(12) = .889
!      x2(13) = .006
!      x2(14) = .828
!      x2(15) = .399
!      x2(16) = .617
!      x2(17) = .939
!      x2(18) = .784
!      x2(19) = .072
!      x2(20) = .889
! 
!      x3(1) = .645
!      x3(2) = .585
!      x3(3) = .310
!      x3(4) = .058
!      x3(5) = .455
!      x3(6) = .779
!      x3(7) = .259
!      x3(8) = .202
!      x3(9) = .028
!      x3(10) = .099
!      x3(11) = .142
!      x3(12) = .296
!      x3(13) = .175
!      x3(14) = .180
!      x3(15) = .842
!      x3(16) = .039
!      x3(17) = .103
!      x3(18) = .620
!      x3(19) = .158
!      x3(20) = .704
! 
!      FCN = 0.0
!      DO i = 1 , 20
!         FCN = (Theta(1)+Theta(2)*x2(i)+(Theta(2)**2)*x3(i)-y(i))**2 + FCN
!      ENDDO
! 
!      END
