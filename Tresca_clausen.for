!      	Subroutine SIGINI(SIGMA,COORDS,NTENS,NCRDS,NOEL,NPT,LAYER,
!     1 KSPT,LREBAR,NAMES)
!!-----------------------------------------------------------------------------
!! Subroutine used to define the initial stress state. In this version
!! The initial stress is linearly dependent on the vertical coordinate z, which
!! is assumed to be zero at the surface and negative below. It is intended
!! for use on with a homogenous soil with horizontal surface. The vertical stress
!! is found from the weight of the soil, and the horizontal stresses are found
!! from the K0 factor, Specific mass, gravity acceleration and K0 must be defined
!! in the subroutine, i.e. they are not given as input parameters.
!!
!! INPUT (most explanations taken from ABAQUS help)
!!
!!	COORDS(ncrds)	An array containing the initial coordinates of this point.
!!	NTENS			Number of stresses to be defined, which depends on the element type.
!!	NCRDS			Number of coordinates.
!!	NOEL			Element number.
!!	NPT				Integration point number in the element.
!!	LAYER			Layer number (for composite shells and layered solids).
!!	KSPT			Section point number within the current layer.
!!	LREBAR			Rebar flag. If LREBAR=1, the current integration point is associated with element rebar. Otherwise, LREBAR=0.
!!	NAMES(1)		Name of the rebar to which the current integration point belongs, which is the name given in the rebar or rebar layer definition (defining reinforcement,?Section 2.2.3 of the ABAQUS Analysis User's Manual, or defining rebar as an element property,?Section 2.2.4 of the ABAQUS Analysis User's Manual). If no name was given in the rebar or rebar layer definition, this variable will be blank. This variable is relevant only when LREBAR=1.
!!	NAMES(2)		Element type name (see Section I.1, abAQUS/Standard Element Index,?of the ABAQUS Analysis User's Manual).
!!
!! OUTPUT
!!
!!	Sigma(ntens)	Stress vector. z is taken to be the vertical coordinate,
!!					  while x and y are horizontal.
!!					  Some examples of the ordering of the stress components are
!!					  3D: x, y, z, xy, xz, yz
!!
!!-----------------------------------------------------------------------------
!! Johan Clausen
!! Department of Civil Engineering
!! Aalborg University
!! January 2008
!!-----------------------------------------------------------------------------
!      INCLUDE 'ABA_PARAM.INC'
!	
!	integer(4) ntens, ncrds
!	DIMENSION SIGMA(NTENS),COORDS(NCRDS)
!      CHARACTER NAMES(2)*80
!	real(8), parameter :: g	= 9.82_8   ! [m/s^2] Gravitational accelaration
!      real(8), parameter :: rho = 2036.7_8 ! [kg/m^3] Mass density
!      real(8), parameter :: K0  = 1.0_8	   ! [-] Earth pressure coefficient at rest
!	
!!------------------------------------------------------------
!!		User coding
!!------------------------------------------------------------

!      Sigma = 0.0_8
!	
!	if (ntens == 6) then ! General 3D with z (or 3) as the vertical coordinate
!		Sigma(3) = g*rho*Coords(3)
!		Sigma(1) = K0*Sigma(3)
!		Sigma(2) = Sigma(1)
!	elseif (ntens == 4) then ! Plane situations with y (or 2) as the vertical coordinate
!		Sigma(2) = g*rho*Coords(2)
!		Sigma(1) = K0*Sigma(2)
!		Sigma(3) = Sigma(1)
!	end if
!!------------------------------------------------------------------------------      
!	end subroutine SigIni !--------------------------------------------------
!!------------------------------------------------------------------------------	
!	
!      subroutine UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
!     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
!     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,NDI,NSHR,NTENS,
!     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
!     4 DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

!!-------------------------------------------------------------------------
!! Plastic strain update for a linearly elastic - perfectly plastic
!! Tresca material with input parameters
!!     1) Young's modulus, E
!!     2) Poisson's ratio, nu
!!     3) Uniaxial yield stress, yield (yield = 2*cohesion)
!!
!! This version is intended for use with ABAQUS as a user programmable 
!! feature. Because of this, the dummy arguments have the the names of the 
!! example on user plasticity provided by ABAQUS.
!!
!! Input arguments:
!!       
!!   STRAN(NTENS)         -  total strain at beginning of increment
!!   DSTRAN(NTENS)        -  strain increments
!!   TIME(1)              -  step time at beginning of increment
!!   TIME(2)              -  total time at beginning of increment
!!   DTIME                -  time increment
!!   TEMP                 -  temperature at beginning of increment
!!   DTEMP                -  increment of temperature
!!   PREDEF               -  interp. values of predef. field variables
!!   DPRED                -  increment in predef. field variables
!!   CMNAME               -  name given on *material option.
!!   NDI                  -  no. of direct stress/strain components
!!   NSHR                 -  no. of shear stress/strain components
!!   NTENS                -  total no. of stress/strain components
!!   NSTATV               -  no. of state variables
!!   PROPS(NPROPS)        -  array with material properties
!!   NPROPS               -  no. of material properties
!!   COORDS(3)            -  coordinates at point
!!   DROT(3,3)            -  rotation increment vector
!!   CELENT               -  characteristic element length
!!   DFGRD0(3,3)          -  deformation increment at start of step
!!   DFGRD1(3,3)          -  deformation increment at end of step
!!   NOEL                 -  element number
!!   NPT                  -  integration point number
!!   LAYER                -  layer number (not used)
!!   KSPT                 -  section point number (not used)
!!   KSTEP                -  step number
!!   KINC                 -  increment number
!!       
!! Arguments that can be updated:
!!
!!   PNEWDT               -  ratio of suggested new time increment
!!       
!! Output arguments:
!!
!!   DDSDDE(NTENS,NTENS)  -  the elastoplastic constitutive tensor
!!   STRESS(NTENS)        -  array with stresses
!!   STATEV(NSTATEV)      -  statevariables**
!!   SSE                  -  specific elastic strain energy
!!   SPE                  -  specific plastic dissipation
!!   SCD                  -  specific "creep" dissipation
!!   RPL                  -  volumetric heat generation*
!!   DDSDDT(NTENS)        -  stress variation due to temperature*
!!   DRPLDE(NTENS)        -  variation of RPL due to strain*
!!   DRPLDT               -  variation of RPL due to temperature*
!!
!! *:  Only required in fully coupled temperature-displacement analysis
!! **: Four values are saved in STATEV, which means that Depvar = 4 in the edit
!!     material menu in Abaqus CAE.
!!	- STATEV(1): Region of the updated stress point
!!        0) Elastic
!!        1) Stress point on a plane of the failure envelope
!!        2) Stress point on the triaxial compression meridian
!!        3) Stress point on the triaxial tension meridian
!!	- STATEV(2): Increment in plastic shear strains
!!	- STATEV(3): Increment in equivalent plastic strain
!!	- STATEV(4): Total equivalent plastic strain
!!
!! REFERENCES:
!! [1] Johan Clausen, Lars Damkilde and Lars Andersen: "Efficient return
!!	  algorithms for associated plasticity with multiple yield planes", 
!!	  International Journal for Numerical Methods in Engineering, 2006,
!!	  volume 66, pages 1036-1059.
!! [2] Johan Clausen, Lars Damkilde and Lars Andersen: An efficient return
!!	  algorithm for non-associated plasticity with linear yield criteria
!!	  in principal stress space. Computers & Structures, 2007, volume 85,
!!	  pages 1795-1807.
!! [3] Johan Clausen: Efficient Non-Linear Finite Element Implementation of
!!	  Elasto-Plasticity for Geotechnical Problems. Ph.D. Thesis. 2006.
!!	  Esbjerg Technical Institute, Aalborg University. Can be downloaded
!!	  from http://vbn.aau.dk/fbspretrieve/14058639/JCthesis.pdf
!! [4] Johan Clausen, Lars Damkilde and Lars Andersen: "Robust and efficient
!!       handling of yield surface discontinuities in elasto-plastic finite
!!       element calculations", Engineering Computations, 2015, volume 32,
!!       No. 6, pages 1722-1752.
!!-----------------------------------------------------------------------
!! Johan Clausen                          Lars Andersen
!! Department of Civil Engineering    &   Department of Civil Engineering
!! Aalborg University                     Aalborg University
!! Denmark                                Denmark
!!-----------------------------------------------------------------------
!! Original code: May 2008
!! Latest modification: May 2016
!!-----------------------------------------------------------------------
!	implicit none
!      !INCLUDE 'ABA_PARAM.INC'

!!-----------------------------------------------------------------------
!!     Variables defined by ABAQUS input/output
!!-----------------------------------------------------------------------

!      CHARACTER*80 CMNAME
!      integer(4) NTENS,NDI,NSHR,NSTATV,NPROPS,NOEL,NPT,LAYER,KSPT,
!     &	KSTEP,KINC
!      real(8) SPD,SCD,RPL,DRPLDT,DTIME,DTEMP,PNEWDT,CELENT
!      real(8) STRESS(NTENS),STATEV(NSTATV),
!     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
!     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
!     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

!!-----------------------------------------------------------------------
!!     Variables defined by this subroutine
!!-----------------------------------------------------------------------

!      integer(4) region, i
!      real(8) SigB(NTENS), SigC(NTENS), 
!     1 D(NTENS,NTENS), Dinv(NTENS,NTENS), Depc(NTENS,NTENS),
!     2 StrainElas(NTENS), DelStrainPlas(NTENS),
!     3 yield, E, nu, ElastPar(2),
!	4 DSTRES(NTENS), SSE, SPE, DEE, TDE, DPE_eq, temp, DPE12

!      real(8), parameter :: pi = 3.1415926535897932_8

!!---- Explicit interfaces ----------------------------------
!	interface
!		pure subroutine TrescaStressReturn(Sigma,nsigma,yield,D,Dinv,
!     &											Sigma_up,Depc,region)
!			integer(4), intent(in)  :: nsigma
!			real(8), intent(in)     :: Sigma(nsigma), yield
!			real(8), intent(in)     :: D(nsigma,nsigma), Dinv(nsigma,nsigma)

!			real(8), intent(out)	:: Sigma_up(nsigma), Depc(nsigma,nsigma)
!			integer(4), intent(out) :: region
!		end subroutine TrescaStressReturn
!		
!		pure subroutine DlinElas(E,nu,nsigma,D,Dinv)
!			integer(4), intent(in) :: nsigma
!			real(8), intent(in) :: E, nu

!			real(8), intent(out) :: D(nsigma,nsigma)
!			real(8), intent(out), optional :: Dinv(nsigma,nsigma)
!		end subroutine DlinElas
!	end interface
!!-----------------------------------------------------------

!!-----------------------------------------------------------------------
!!	Some info about output to ABAQUS data and message files
!!-----------------------------------------------------------------------
!!
!!	write(6,*) 'FILE ID# 6' ! .DAT-file
!!	write(7,*) 'FILE ID# 7' ! .MSG-file
!!
!!-----------------------------------------------------------------------

!!-----------------------------------------------------------------------
!!     Material Properties
!!-----------------------------------------------------------------------

!      E      =  PROPS(1)  !  E    - Young's modulus
!      nu     =  PROPS(2)  !  nu   - Poisson's ratio
!      yield  =  PROPS(3)  !  yield - uniaxial yield stress (yield = 2*cohesion)

!!-----------------------------------------------------------------------
!!     Predictor stress
!!-----------------------------------------------------------------------

!      call DlinElas(E,nu,NTENS,D,Dinv)

!	DSTRES = matmul(D,DSTRAN) ! stress increment
!	SigB   = STRESS + DSTRES  ! elastic predictor stress

!!-----------------------------------------------------------------------
!!     Stress update and creation of consistent constitutive matrix
!!-----------------------------------------------------------------------

!	if ((KINC.eq.1).and.(NOEL.eq.1).and.(NPT.eq.1)) then
!		write(6,'(''      UMAT Material Model: '',
!	1		        ''Tresca'')')
!	endif

!	call TrescaStressReturn(SigB,NTENS,yield,D,Dinv,SigC,Depc,region)
!      
!!-----------------------------------------------------------------------
!!	Plastic shear strain increment
!!-----------------------------------------------------------------------

!	DPE12 = DSTRAN(NDI+1) - 2*(1D0+nu)*(SigC(NDI+1)-STRESS(NDI+1))/E

!!-----------------------------------------------------------------------
!!	Plastic strain increment
!!-----------------------------------------------------------------------
!	
!	if ((KSTEP.EQ.1).and.(KINC.eq.1)) then
!		DPE_eq = 0.0_8
!		if (NSTATV.GE.4) STATEV(4) = 0.0_8
!	else
!		DelStrainPlas = matmul(Dinv,(SigB-SigC)) ! Plastic strain increment
!		DPE_eq = sqrt(dot_product(DelStrainPlas,DelStrainPlas)
!     &									*0.6666666666666667_8) ! Equivalent plastic strain increment
!	end if

!!-----------------------------------------------------------------------
!!     Postconditioner for ABAQUS
!!-----------------------------------------------------------------------

!      DDSDDE    = Depc          ! consistent constitutive matrix
!      STRESS    = SigC          ! updated stress at the Gauss point

!!-----------------------------------------------------------------------
!!	Type of failure and Lode angle
!!-----------------------------------------------------------------------

!      STATEV(1) = real(region,8)	! type of stress return
!      if (NSTATV >= 1) STATEV(2) = DPE12			! Plastic shear strain increment
!	if (NSTATV >= 2) STATEV(3) = DPE_eq			! Equivalent plastic strain increment
!	if (NSTATV >= 3) STATEV(4) = STATEV(4) + DPE_eq ! Total equivalent plastic strain

!!-----------------------------------------------------------------------
!      end subroutine UMAT
!!-----------------------------------------------------------------------


      pure subroutine TrescaStressReturn(Sigma,nsigma,yield,D,Dinv,
     &											Sigma_up,Depc,region)
!------------------------------------------------------------------------
! Wrapping subroutine for return mapping with a linearly elastic perfectly
! plastic Tresca material with an associated flow rule. Performed in
! principal coordinates.
!
! The Tresca criterion is defined by
!		f = sigP(1) - sigP(3) - yield = 0
! 
! INPUT
!  - Name -		-type and description --
!	Sigma		real array(nsigma). Stress vector. Size is dependent on the type of stress.
!				  The ordering of the yieldonents must be as follows:
!				    Plane situations (plane stress, plane strain and axisymmetry):
!					Sigma = [sig_x sig_y sig_z tau_xy] (i.e. sig_z is out-of-plane)
!					General 3D:  Sigma = [sig_x sig_y sig_z tau_xy tau_xz tau_yz]
!	nsigma		integer,scalar. Number of stress components, i.e. the size of Sigma
!	yield		real, scalar. Uniaxial yield stress. Is defined by yield = 2*c,
!							  where c is the cohesion strength
!	D			real, array(nsigma,nsigma) Elastic constitutive matrix.
!	Dinv		real, array(nsigma,nsigma) Inverted elastic constitutive matrix (compliance matrix).
!
! OUTPUT
!	Sigma_up	real, array(nsigma) Stress vector containg the updated stresses.
!	Depc		real, array(nsigma,nsigma) Consistent constitutive elasto-plastic matrix
!				  If the material is not yielding Depc = D
!	region		integer, scalar. Number of the region of the returned stress
!					region = 0: the elastic region, i.e. no yielding
!					region = 1: Return to the Tresca plane
!					region = 2: Return to the triaxial compressive line, sigp_1 = sigp_2
!					region = 3: Return to the triaxial tensile line, sigp_2 = sigp_3
! REFERENCES:
! [1] Johan Clausen, Lars Damkilde and Lars Andersen: "Efficient return
!	  algorithms for associated plasticity with multiple yield planes", 
!	  International Journal for Numerical Methods in Engineering, 2006,
!	  volume 66, pages 1036-1059.
! [2] Johan Clausen, Lars Damkilde and Lars Andersen: An efficient return
!	  algorithm for non-associated plasticity with linear yield criteria
!	  in principal stress space. Computers & Structures, 2007, volume 85,
!	  pages 1795-1807.
! [4] Lars Andersen and Johan Clausen: "Comments on "Flow rule effects in 
!       the TresTresca model?by H.A. Taiebat and J.P. Carter [Computers
!       and Geotechnics 35 (2008) 500?03]", Computers and Geotechnics, 2009,
!       volume 36, pages 911-913. 
! [3] Johan Clausen, Lars Damkilde and Lars Andersen: "Robust and efficient
!       handling of yield surface discontinuities in elasto-plastic finite
!       element calculations", Engineering Computations, 2015, volume 32,
!       No. 6, pages 1722-1752.
!------------------------------------------------------------------------
! Johan Clausen
! Department of Civil Engineering
! Aalborg University
! Original code: May 2008
! Latest revision: January 2011
!------------------------------------------------------------------------
	implicit none

	!-----------------------
	integer(4), intent(in)  :: nsigma
	real(8), intent(in)     :: Sigma(nsigma), yield
	real(8), intent(in)     :: D(nsigma,nsigma), Dinv(nsigma,nsigma)

	real(8), intent(out)	:: Sigma_up(nsigma), Depc(nsigma,nsigma)
	integer(4), intent(out) :: region
	!-----------------------

	integer(4) hoop,ouplP, s1, s2, nshear
	real(8) psi(3,3), f, SigP(3), SigP_up(3)
	real(8), dimension(nsigma,nsigma) :: DepP, DepcP, A, Atrans
	real(8), dimension(nsigma,nsigma) :: MatProd, T
	real(8) Fnorm(3), Lfdir(3), SiPla(3)
	real(8) Tshear(nsigma-3,nsigma-3),Tshearinv(nsigma-3,nsigma-3)
	
	!--- Settings regarding calculation of the constitutive matrix, DepP -------------
	integer(4), parameter :: DepFormLine = 1 ! Form of the consistent constitutive matrix on a line
					! DepFormLine = 0: Double singular constitutive matrix
					! DepFormLine = 1: Modified "almost double singular" matrix on a line ("Tiln鎋b=axr" in the MatLab code)
					! DepFormLine = 2: A single singular matrix in the Koiter direction, i.e. the direction of the plastic strain.
	real(8), parameter :: beta = 40.0_8 ! The value used in the modification of DepP when DepFormLine = 1
	!--------------------------------------------------------------------------------

	!---- Explicit interfaces ----------------------
	interface
		pure subroutine PrinStressAna(Sigma,nsigma,SigP)
			integer(4), intent(in) :: nsigma
			real(8), intent(in) :: Sigma(nsigma)
			real(8), intent(out) :: SigP(3)
		end subroutine PrinStressAna
		
		pure subroutine PrinRetTresca(SigP,f,yield,Dfull,nsigma,SigP_up,
     &														region)
			integer(4), intent(in) :: nsigma
			real(8), intent(in)	   :: SigP(3), f, yield
			real(8), intent(in)	   :: Dfull(nsigma,nsigma)

			real(8), intent(out)	:: SigP_up(3)
			integer(4), intent(out) :: region
		end subroutine PrinRetTresca
		
		pure subroutine TshearPrinPerfect(SigP,SigP_up,nshear,s1,s2,
     1												Tshear,Tshearinv)
			real(8), intent(in) :: SigP(3), SigP_up(3)
			integer(4), intent(in) :: nshear, s1, s2
			real(8), dimension(nshear,nshear), intent(out) :: Tshear, Tshearinv
		end subroutine TshearPrinPerfect
		
		pure subroutine FormDepPerfect(D,Norm,Edir,nsigma,Dep)
			integer(4), intent(in) :: nsigma
			real(8), intent(in) :: D(nsigma,nsigma), Norm(nsigma)
			real(8), intent(in) :: Edir(nsigma)
			real(8), intent(out) :: Dep(nsigma,nsigma)
		end subroutine FormDepPerfect
		
		pure subroutine FormModDepLine(mod_type,beta,SiPla,Fline,Gline,Dninv,
     1											Dc,Dcinv,nsigma,Depc)
			integer(4), intent(in) :: mod_type, nsigma
			real(8), intent(in) :: Dninv(3,3),Dc(nsigma,nsigma),Gline(3)
			real(8), intent(in) :: beta,Sipla(3),Fline(3),Dcinv(nsigma,nsigma)
			real(8), intent(out) :: Depc(nsigma,nsigma)
		end subroutine FormModDepLine
		
		pure subroutine PrinDirect(Sigma,nsigma,SigP,psi)
			integer(4), intent(in) :: nsigma
			real(8), intent(in) :: Sigma(nsigma), SigP(3)
			real(8), intent(out)  :: psi(3,3)
		end subroutine PrinDirect
		
		pure subroutine TransMatrix(psi,nsigma,ouplP,A)
			integer(4), intent(in) :: nsigma, ouplP
			real(8), intent(in)    :: psi(3,3)
			real(8), intent(out)   :: A(nsigma,nsigma)
		end subroutine TransMatrix
      
	end interface
	!------------------------------------------
!------------------------------------------------------------------------
	call PrinStressAna(Sigma,nsigma,SigP) ! Principal stresses SigP
!
!	--- Value of Yield function f ---
      f = SigP(1) - SigP(3) - yield
!     ---------------------------------
!
      if (f > 0.0_8) then
!		--- Position of out-of plane or hoop stress in principal stress vector SigP ----
		if (nsigma == 4) then
			if (Sigma(3) == SigP(1)) then
				ouplP = 1 ; s1 = 2 ; s2 = 3
			elseif (Sigma(3) == SigP(2)) then
				ouplP = 2 ; s1 = 1 ; s2 = 3
			elseif (Sigma(3) == SigP(3)) then
				ouplP = 3 ; s1 = 1 ; s2 = 2
			end if
		end if
!		---------------------------------------------------------------------------------
!
!		--- Stress return in principal coordinates --
		call PrinRetTresca(SigP,f,yield,D,nsigma, SigP_up,region)
		SiPla = SigP - SigP_up ! Plastic corrector stress in principal stress space
!		---------------------------------------------
!
!		--- Modification matrix, T ----------------------
		nshear = nsigma - 3 ! Number of shear components
		call TshearPrinPerfect(SigP,SigP_up,nshear,s1,s2,
     &										Tshear,Tshearinv) ! Forms the shear part modification matrix T
		T = 0.0_8
		T(1:3,1:3) = reshape((/1.0_8, 0.0_8, 0.0_8,	
     &						   0.0_8, 1.0_8, 0.0_8,	
     &						   0.0_8, 0.0_8, 1.0_8/), (/3,3/)) ! Unit Matrix
		T(4:nsigma,4:nsigma) = Tshear
!		-------------------------------------------------


!		--- Infinitesimal constitutive matrix, DepP -----------------	
		DepP = 0.0_8
		if (region == 1) then ! Return to the yield surface 
			DepP(4:nsigma,4:nsigma) = D(4:nsigma,4:nsigma)	  ! Shear components of consistent constitutive matrix in principal stress space

			Fnorm = [1.0_8, 0.0_8, -1.0_8] ! Yield plane normal			
			call FormDepPerfect(D(1:3,1:3),Fnorm,Fnorm,3,DepP(1:3,1:3)) ! Normal components of consistent constitutive matrix in principal stress space

		elseif (region == 2 .or. region == 3) then ! return to a line 1 or line 2 (triaxial compression, sigp1 = sigp2, or triaxial extension, sigp2 = sigp3)
			Lfdir = (/1.0_8, 1.0_8, 1.0_8/) ! Edge line direction
			
			call FormModDepLine(DepFormLine,beta,SiPla,Lfdir,Lfdir,
     &								Dinv(1:3,1:3),D,Dinv,nsigma,DepP)
		end if ! region = ...
!		------------------------------------------------------------

!		--- Consistent constitutive matrix, DepcP -------	
		DepcP = matmul(T,DepP)
!		-------------------------------------------------

!		--- Tranformation matrix A ------------------
		call PrinDirect(Sigma,nsigma,SigP,psi)
		call TransMatrix(psi,nsigma,ouplP,A)
!		---------------------------------------------

!		--- Coordinate transformation ---------------
		Atrans = transpose(A)
		Sigma_up = 0.0_8
      
		Sigma_up = matmul(Atrans(:,1:3),SigP_up)
		MatProd = matmul(Atrans,DepcP)
		Depc = matmul(MatProd,A)
!		---------------------------------------------
!
      else ! no yielding
		Sigma_up = Sigma
		Depc = D
		region = 0
      end if ! f > 0
!
!--------------------------------------------------------------------------
      end subroutine TrescaStressReturn !-----------------------------
!--------------------------------------------------------------------------

		
      pure subroutine PrinRetTresca(SigP,f,yield,Dfull,nsigma,SigP_up,
     &														region)
!-----------------------------------------------------------------------
! Stress return for a linearly elastic - perfectly plastic Tresca 
! material in principal stress space. The flow rule is associated.
! The Tresca criterion is written as f = sigP_1 - sigP_3 - yield = 0
!
! INPUT
!  - Name -	-type,size -	-- Description --
!	SigP	(real,ar(3))	Principal predictor stresses in descending order
!							  SigP = [sigP_1 sigP_2 sigP_3]
!	f		(real,sc)		Value of the yield function given by
!							  f = k*sigP_1 - sigP_3 - yield
!	yield	(real,ar(3))	Uniaxial yield stress. Is defined by yield = 2*c,
!							  where c is the cohesion strength
!	Dfull	(real,ar(nsigma,nsigma)) Elastic constitutive matrix
!	nsigma	(integer,sc)	Number of stress components, Dfull
!
! OUTPUT
!	SigP_up	(real,ar(3))	Updated principal stresses, i.e. they obey the
!							  Tresca yield criterion
!	region	(integer,sc)	Type of return, i.e. the region of the
!							  predictor stress
! REFERENCES:
! [1] Johan Clausen, Lars Damkilde and Lars Andersen: "Efficient return
!	  algorithms for associated plasticity with multiple yield planes", 
!	  International Journal for Numerical Methods in Engineering, 2006,
!	  volume 66, pages 1036-1059.
! [2] Johan Clausen, Lars Damkilde and Lars Andersen: An efficient return
!	  algorithm for non-associated plasticity with linear yield criteria
!	  in principal stress space. Computers & Structures, 2007, volume 85,
!	  pages 1795-1807.
!----------------------------------------------------------------------
! Johan Clausen
! Department of Civil Engineering
! Aalborg University
! May 2008
!----------------------------------------------------------------------
	implicit none
	
	!-------------------
	integer(4), intent(in) :: nsigma
	real(8), intent(in)	 :: SigP(3), f, yield
	real(8), intent(in)	 :: Dfull(nsigma,nsigma)

	real(8), intent(out)	:: SigP_up(3)
	integer(4), intent(out) :: region
	!-------------------

	real(8) D(3,3), P1(3), P2(3)
	real(8) den, Rp(3), den_t
	real(8) Rp2(3), Rp3(3), N2(3), N3(3), t1, t2
	real(8) N(3), pI_II, pI_III

	!----- Stress return ------------------------------------

	!--- Preliminary parameters ----------------
	D = Dfull(1:3,1:3)! Relates to normal stresses

	P1 = (/1.0_8,  1.0_8, -2.0_8/)*yield/3 ! Stress point on line 1 (triaxial compression)
	P2 = (/2.0_8, -1.0_8, -1.0_8/)*yield/3 ! Stress point on line 2 (triaxial extension)
	
	den = D(1,1) - D(1,3) - D(3,1) + D(3,3) ! denominator a'*D*a
	Rp(1) = (D(1,1)-D(1,3))/den ! Rp is the scaled direction of
	Rp(2) = 0.0_8				  ! the update direction,
	Rp(3) = (D(3,1)-D(3,3))/den ! Rp = D*b/(a'*D*a) a is the gradient of the
							      ! yield surface, a = [1 0 -1]'
	!-------------------------------------------

	!--- Boundary planes -----------------------

	! Boundary plane normal
	N(1) = Rp(2) - Rp(3) ! NI_II = cross(Rp,R1)
	N(2) = Rp(3) - Rp(1) ! R1 = [1 1 1], direction
	N(3) = Rp(1) - Rp(2) ! vector of line 1
 
	! Boundary plane between regions I and II
	pI_II  = dot_product(N,SigP-P1)
 
	! Boundary plane between regions I and III
	pI_III = dot_product(N,SigP-P2)
	!-------------------------------------------


	!--- Region determination and update -------
	if (pI_II < 0.0_8) then
		region = 2
	
		!--- t-paramter for line 1 ------
		! secondary surface in region II a = [0 1 -1]
		den = D(2,2) - D(2,3) - D(3,2) + D(3,3) ! denominator a'*D*a
		Rp2(1) = 0.0_8				 ! Rp is the scaled direction of
		Rp2(2) = (D(2,2)-D(2,3))/den ! the update direction,
		Rp2(3) = (D(3,2)-D(3,3))/den ! Rp = D*a/(a'*D*a) a is the gradient the secondary
									 ! yield surface in Region II
	
		N2(1) = Rp(2)*Rp2(3) - Rp(3)*Rp2(2) ! N2 = cross(Rp,Rp2)
		N2(2) = Rp(3)*Rp2(1) - Rp(1)*Rp2(3)
		N2(3) = Rp(1)*Rp2(2) - Rp(2)*Rp2(1)
	
		den_t = N2(1) + N2(2) + N2(3) ! N2'*R1, R1 = [1 1 1] Direction of line 1
		! t-parameter of line 1
		t1 = dot_product(N2,SigP-P1)/den_t
		!--------------------------------

		!-- Updated stress ----------------------
		SigP_up(1) = t1 + P1(1) ! SigP_up = t1*R1 + P1
		SigP_up(2) = t1 + P1(2) ! R1 = [1 1 1], direction
		SigP_up(3) = t1 + P1(3) ! vector of line 1
		!----------------------------------------

	elseif (pI_III <= 0.0_8) then
		region = 1

		SigP_up = SigP - f*Rp ! SigP_up = SigP - SiPla
							  ! SiPla is the plastic corrector
							  ! given by f*Rp

	else
		region = 3	
	
		!--- t-paramter for line 2 ------
		 ! secondary surface in region III a = [1 -1 0]
		den = D(1,1) - D(1,2) - D(2,1) + D(2,2) ! denominator a'*D*a
		Rp3(1) = (D(1,1)-D(1,2))/den ! Rp is the scaled direction of
		Rp3(2) = (D(2,1)-D(2,2))/den ! the update direction,
		Rp3(3) = 0.0_8				 ! Rp = D*a/(a'*D*a) a is the gradient of the secondary
				  				     ! yield surface in region III
	
		N3(1) = Rp(2)*Rp3(3) - Rp(3)*Rp3(2) ! N3 = cross(Rp,Rp3)
		N3(2) = Rp(3)*Rp3(1) - Rp(1)*Rp3(3)
		N3(3) = Rp(1)*Rp3(2) - Rp(2)*Rp3(1)
	
		den_t = N3(1) + N3(2) + N3(3) ! N3'*R2, R2 = [1 1 1] Direction of line 2
		! t-parameter of line 2
		t2 = dot_product(N3,SigP-P2)/den_t
		!--------------------------------
	
		!-- Updated stress ----------------------
		SigP_up(1) = t2 + P2(1) ! SigP_up = t2*R2 + P2
		SigP_up(2) = t2 + P2(2) ! R2 = [1 1 1], direction
		SigP_up(3) = t2 + P2(3) ! vector of line 2
	!----------------------------------------

	end if ! Regions and update
	!----------------------------------------------

	!------------------------------------------------------------------

!--------------------------------------------------------------------------------
	end subroutine PrinRetTresca !---------------------------------------------
!--------------------------------------------------------------------------------


!	pure subroutine FormModDepLine(mod_type,beta,SiPla,Fline,Gline,Dninv,
!     &											Dc,Dcinv,nsigma,Depc)
!!--------------------------------------------------------------------------
!! Calculates the modified elasto-plastic constitutive matrix at on a line
!! for a perfectly plastic material, It is modified in the sence that is not
!! the double singular matrix that theory prescribes. Only one type of
!! modification is possible. The calculation are carried out in principal
!! stress space.
!!
!! INPUT
!! - Name -		-type and description --
!!	mod_type	(integer,scalar) Detemines the modified type
!!				  mod_type= 0: A Doubly singular matrix, as given in the
!!							   References [1] and [2].
!!				  mod_type= 1: A single singular matrix in the plastic
!!							   strain direction, but almost singular in
!!							   any direction perpendicular to the potential
!!							   line. See Reference [3] for a detailed explanation.
!!	beta		(real,scalar) The factor by which to divide the stiffness in
!!				  the directions perpendicular to the line.
!!	SiPla		(real,array(3)) The plastic corrector principal stresses.
!!	Fline		(real,array(3)) Direction of the yield surface line in principal
!!				  stress space.
!!	Gline		(real,array(3)) Direction of the plastic potential line in principal
!!				  stress space.
!!	Dninv		(real,arrey(3,3)) Normal components of the elastic
!!				  compliance matrix.
!!	Dc			(real,array(nsigma,nsigma)) Modified elastic stiffness.
!!				  If the yield criterion is linear Dc = D.
!!	Dcinv		(real,array(nsigma,nsigma)) Inverse of the modified
!!				  elastic stiffness. If the yield criterion is linear,
!!				  Dcinv = Dinv.
!!	nsigma		(integer, scalar) Number of stress components
!!
!! OUTPUT
!!	Depc	(real,array(nsigma,nsigma) Modified Elasto-plastic constitutive
!!			  matrix on a line. If the criterion is linear Depc is the
!!			  infinitesimal version, otherwise it is the consistent. In
!!			  both cases it is expressed in principal stres space.
!!
!! REFERENCES:
!! REFERENCES:
!! [1] Johan Clausen, Lars Damkilde and Lars Andersen: "Efficient return
!!	  algorithms for associated plasticity with multiple yield planes", 
!!	  International Journal for Numerical Methods in Engineering, 2006,
!!	  volume 66, pages 1036-1059.
!! [2] Johan Clausen, Lars Damkilde and Lars Andersen: "An efficient return
!!	  algorithm for non-associated plasticity with linear yield criteria
!!	  in principal stress space". Computers & Structures, 2007, volume 85,
!!	  pages 1795-1807.
!! [3] Johan Clausen: Efficient Non-Linear Finite Element Implementation of
!!	  Elasto-Plasticity for Geotechnical Problems. Ph.D. Thesis. 2006.
!!	  Esbjerg Technical Institute, Aalborg University. Can be downloaded
!!	  from http://vbn.aau.dk/fbspretrieve/14058639/JCthesis.pdf
!! [4] Johan Clausen, Lars Damkilde and Lars Andersen: "Robust and efficient
!!       handling of yield surface discontinuities in elasto-plastic finite
!!       element calculations", Engineering Computations, 2015, volume 32,
!!       No. 6, pages 1722-1752.
!!------------------------------------------------------------------------
!! Johan Clausen
!! Department of Civil Engineering
!! Aalborg University
!! April 2008
!!------------------------------------------------------------------------
!	implicit none
!	!----- Input/Output ------
!	integer(4), intent(in) :: mod_type, nsigma
!	real(8), intent(in) :: Dninv(3,3),Dc(nsigma,nsigma),Gline(3)
!	real(8), intent(in) :: beta,Sipla(3),Fline(3),Dcinv(nsigma,nsigma)

!	real(8), intent(out) :: Depc(nsigma,nsigma)
!	!-------------------------

!	real(8) KoitDir(3),KoitPerDir(3), Dper(3,3)
!	
!	!---- Explicit interfaces ----------------------
!	interface
!		pure subroutine FormDepPerfect(D,Norm,Edir,nsigma,Dep)
!			integer(4), intent(in) :: nsigma
!			real(8), intent(in) :: D(nsigma,nsigma), Norm(nsigma)
!			real(8), intent(in) :: Edir(nsigma)
!			real(8), intent(out) :: Dep(nsigma,nsigma)
!		end subroutine FormDepPerfect
!		!------
!		pure function Cross(A,B)
!		  real(8), intent(in)	::	A(3), B(3)
!		  real(8)				::  Cross(3)
!		end function Cross
!		!-----
!		pure subroutine FormDepLinePerfect(D,Dinv,Ra,Rb,nsigma,DepLine)
!			integer(4), intent(in) :: nsigma
!			real(8), intent(in), dimension(nsigma,nsigma) :: D, Dinv
!			real(8), intent(in), dimension(3) :: Ra, Rb
!			real(8), intent(out) :: DepLine(nsigma,nsigma)
!		end subroutine FormDepLinePerfect
!	end interface
!	!-----------------------------------------------
!!------------------------------------------------------------------------
!	Depc = 0.0_8

!	call FormDepLinePerfect(Dc,Dcinv,Fline,Gline,nsigma,Depc)
!	
!	if (mod_type == 1 .and. dot_product(SiPla,SiPla) > 0.0_8) then
!		KoitDir = matmul(Dninv,SiPla) ! the plastic strain direction
!		KoitPerDir = Cross(KoitDir,Gline) ! Direction perpendicular to the strain direction
!											 ! and the plastic potential line
!		call FormDepLinePerfect(Dc(1:3,1:3),Dcinv(1:3,1:3),
!     &									KoitPerDir,KoitPerDir,3,Dper)
!		Depc(1:3,1:3) = Depc(1:3,1:3) + Dper/beta
!	elseif (mod_type == 2 .and. dot_product(SiPla,SiPla) > 0.0_8) then
!		KoitDir = matmul(Dninv,SiPla) ! The plastic strain direction in principal stress space
!		call FormDepPerfect(Dc(1:3,1:3),KoitDir,KoitDir,3,Depc(1:3,1:3))
!		Depc(4:nsigma,4:nsigma) = Dc(4:nsigma,4:nsigma)
!	end if

!!--------------------------------------------------------------------------------
!	end subroutine FormModDepLine !--------------------------------------------
!!--------------------------------------------------------------------------------

	
!	pure subroutine TshearPrinPerfect(SigP,SigP_up,nshear,s1,s2,
!     &												Tshear,Tshearinv)
!!--------------------------------------------------------------------------
!! Forms the shear part of the modification matrix T in principal stress space
!! for at perfectly plastic material. T is used when calculation the modified
!! elastic stiffness matrix Dc which is, in turn, used to form the consistent
!! constitutive matrix. The shear part of T is unaffected whether the yield
!! criterion is linear or not.
!! 
!!
!! INPUT
!!  - Name -	-type,size -	-- Description --
!!	SigP	 real,ar(3)		Principal predictor stresses in descending order
!!							  SigP = [sigP_1 sigP_2 sigP_3]
!!	SigP_up	 real,ar(3)		Updated principal stresses in descending order
!!	nshear	 int,sc			Number of shear stress components. 1 in plane problems
!!							  and 3 in full 3D
!!	s1		 int,sc			Position of largest principal in-plane stress in SigP.
!!							  Only used in plane situations
!!	s2		 int,sc			Position of smallest principal in-plane stress in SigP.
!!							  Only used in plane situations
!!
!! OUTPUT
!!	Tshear		real,ar(nshear,nshear) Shear part of the modification matrix T
!!	Tshearinv	real,ar(nshear,nshear) Inverse of Tshear
!! REFERENCES:
!! [1] Johan Clausen, Lars Damkilde and Lars Andersen: "Efficient return
!!	  algorithms for associated plasticity with multiple yield planes", 
!!	  International Journal for Numerical Methods in Engineering, 2006,
!!	  volume 66, pages 1036-1059.
!! [2] Johan Clausen: Efficient Non-Linear Finite Element Implementation of
!!	  Elasto-Plasticity for Geotechnical Problems. Ph.D. Thesis. 2006.
!!	  Esbjerg Technical Institute, Aalborg University. Can be downloaded
!!	  from http://vbn.aau.dk/fbspretrieve/14058639/JCthesis.pdf
!!------------------------------------------------------------------------
!! Johan Clausen
!! Department of Civil Engineering
!! Aalborg University
!! April 2008
!!------------------------------------------------------------------------
!	implicit none
!	real(8), intent(in) :: SigP(3), SigP_up(3)
!	integer(4), intent(in) :: nshear, s1, s2

!	real(8), dimension(nshear,nshear), intent(out) :: Tshear, Tshearinv
!!------------------------------------------------------------------------

!	Tshear = 0.0_8
!	Tshearinv = 0.0_8

!	if (nshear == 1) then
!		if (SigP_up(s1) - SigP_up(s2) > 0.0_8) then
!			Tshear = (SigP_up(s1)-SigP_up(s2)) / (SigP(s1)-SigP(s2))
!			Tshearinv = 1.0_8/Tshear
!		end if
!	elseif (nshear == 3) then
!		if (SigP_up(1) - SigP_up(2) > 0.0_8  .and.
!     &			SigP(1) - SigP(2) > 0.0_8) then
!			Tshear(1,1) = (SigP_up(1)-SigP_up(2)) / (SigP(1)-SigP(2))
!			Tshearinv(1,1) = 1.0_8/Tshear(1,1)
!		end if
!		if (SigP_up(1) - SigP_up(3) > 0.0_8  .and.
!     &			SigP(1) - SigP(3) > 0.0_8) then
!			Tshear(2,2) = (SigP_up(1)-SigP_up(3)) / (SigP(1)-SigP(3))
!			Tshearinv(2,2) = 1.0_8/Tshear(2,2)
!		end if
!		if (SigP_up(2) - SigP_up(3) > 0.0_8  .and.
!     &			SigP(2) - SigP(3) > 0.0_8) then
!			Tshear(3,3) = (SigP_up(2)-SigP_up(3)) / (SigP(2)-SigP(3))
!			Tshearinv(3,3) = 1.0_8/Tshear(3,3)
!		end if
!	end if

!!--------------------------------------------------------------------------------
!	end subroutine TshearPrinPerfect !-----------------------------------------
!!--------------------------------------------------------------------------------


!	pure subroutine PrinStressAna(Sigma,nsigma,SigP)
!!--------------------------------------------------------------------------
!! Calculates the principal stress of a stress vector Sigma. Analytical 
!! expressions are used in the calculations.
!!
!! INPUT
!!  - Name - -type,size -  -- description --
!!     Sigma (real,ar(nsigma)Stress vector. Size is dependent on the type of stress
!!                   nsigma = 4 in plane situations and nsigma = 6 in 3D
!!                   The ordering of the components must be as follows:
!!                   Plane situations (plane stress, plane strain and axisymmetry):
!!                   Sigma = [sig_x sig_y sig_z tau_xy] (i.e. sig_z is out-of-plane)
!!                   General 3D:  Sigma = [sig_x sig_y sig_z tau_xy tau_xz tau_yz]
!!     nsigma  (integer,sc)  Number of stress components, i.e. the size of Sigma
!!
!! OUTPUT
!!     SigP  (real,ar(3))  Principal stresses in descending order
!!
!!---------------------------------------------------------------------------
!! Johan Clausen
!! Esbjerg Institute of Technology
!! Aalborg University
!! November 2005
!!--------------------------------------------------------------------------
!	implicit none
!      !----- Input/Output -----------
!      integer(4), intent(in) :: nsigma
!      real(8), intent(in) :: Sigma(nsigma)

!      real(8), intent(out) :: SigP(3)
!	!------------------------------

!      real(8) sig_av, sig_hj, I1, J2, J3, lode, sigm
!      real(8) sqJ2, sin3lode
!      real(8) S(6)

!      !---- Explicit interfaces ----------------------
!	interface
!		pure subroutine Invariants(Sigma,nsigma,S,I1,J2,J3,lode,sin3lode)
!			integer(4), intent(in) :: nsigma
!			real(8), intent(in) :: Sigma(nsigma)
!			real(8), intent(out) :: S(nsigma), I1, J2
!			real(8), intent(out) :: J3, lode, sin3lode
!		end subroutine Invariants
!	end interface
!	!-----------------------------------------------

!!------------------------------------------------------------------------
!      SigP = 0.0_8

!!----- Plane situations including axisymmetry -----------------------------

!      if (nsigma == 4) then
!        
!        sig_av = 0.5_8*(Sigma(1) + Sigma(2)) ;
!        sig_hj = sqrt((0.5_8*(Sigma(1) - Sigma(2)))**2 + Sigma(4)**2)
!    
!        SigP(1) = sig_av + sig_hj
!        SigP(2) = sig_av - sig_hj
!        SigP(3) = Sigma(3) ! The out-of-plane stress

!      !-- Sorting: SigP(1) >= SigP(2) >= SigP(3) --
!        if (Sigma(3) > SigP(1)) then
!          SigP(3) = SigP(2)  ;  SigP(2) = SigP(1)  ; SigP(1) = Sigma(3)
!        elseif (Sigma(3) > SigP(2) ) then
!          SigP(3) = SigP(2)  ;  SigP(2) = Sigma(3) 
!        end if
!      !--------------------------------------------
!  
!!------------------------------------------------------------------------

!!----- General stress state ---------------------------------------------
!      elseif (nsigma == 6) then
!      
!!   ----- Invariants ----------------------------------------------------
!		call Invariants(Sigma,nsigma,S,I1,J2,J3,lode,sin3lode)
!		sigm = I1/3.0_8 ! Hydrostatic stress

!		if (sqrt(J2) < 1.0D-12*abs(I1)) then
!			SigP(1) = Sigma(1)
!			SigP(2) = Sigma(2)
!			SigP(3) = Sigma(3)
!		else
!			
!!   ----- Principal stresses -----------------------------------------
!			sqJ2 = 1.154700538379252_8*sqrt(J2)
!      !				2/sqrt(3)
!      !                                    2*pi/3
!			SigP(1) = sqJ2 * sin(lode + 2.094395102393195_8) + sigm
!			SigP(2) = sqJ2 * sin(lode)                       + sigm
!			SigP(3) = sqJ2 * sin(lode - 2.094395102393195_8) + sigm
!		endif

!      end if ! nsigma == ...

!!----------------------------------------------------------------------------------
!      end subroutine PrinStressAna !-----------------------------------------------
!!----------------------------------------------------------------------------------


!      pure subroutine FormDepPerfect(D,Norm,Edir,nsigma,Dep)
!!--------------------------------------------------------------------------
!! Calculates elasto-plastic constitutive matrix of a material with 
!! non-associated and perfect plasticity. If Norm = Edir the flow rule
!! is associated
!!
!! INPUT
!!     D   (real) elastic constitutive matrix. Size = nsigma x nsigma
!!     Norm  (real) Gradient (normal) of the yield surface. Size = nsigma
!!     Edir  (real) Gradient (normal) of the plastic potential. Size = nsigma
!!     nsigma  (integer) Number of stress components
!!
!! OUTPUT
!!     Dep   (real) Elasto-plastic infinitesimal constitutive matrix
!!
!!--------------------------------------------------------------------------
!! Johan Clausen
!! Esbjerg Department of Engineering
!! Aalborg University
!! July 2005
!!--------------------------------------------------------------------------
!	implicit none
!      integer(4), intent(in) :: nsigma
!	real(8), intent(in) :: D(nsigma,nsigma), Norm(nsigma)
!	real(8), intent(in) :: Edir(nsigma)
!      real(8), intent(out) :: Dep(nsigma,nsigma)
!!--------------------------------------------------------------------------
!      real(8), dimension(nsigma) :: Num1, Num2
!      real(8) den, Num(nsigma,nsigma)
!      integer(4) i,j
!!--------------------------------------------------------------------------
!      Num1 = matmul(D,Edir)
!      Num2 = matmul(Norm,D)
!      forall(i = 1:nsigma, j = 1:nsigma) Num(i,j) = Num1(i)*Num2(j)
!      den  = dot_product(Norm,Num1)
!      Dep = D - Num/den
!!--------------------------------------------------------------------------
!      end subroutine FormDepPerfect !--------------------------------------
!!--------------------------------------------------------------------------

!	pure subroutine FormDepLinePerfect(D,Dinv,Ra,Rb,nsigma,DepLine)
!!--------------------------------------------------------------------------
!! Calculates double singular elasto-plastic constitutive matrix of a
!! material with non-associated and perfect plasticity. If Ra = Rb the flow
!! rule is associated. The formula is only valid in principal stress space
!!
!! INPUT
!!     D   (real) elastic constitutive matrix. Size = nsigma x nsigma
!!     Dinv  (real) Inverted elastic constitutive matrix. Size = nsigma x nsigma
!!     Ra    (real) Direction of the line on the yield surface in principal
!!           stress space. If the line is curved Ra is the tangent
!!     Rb    (real) Direction of the line on the plastic potential in principal
!!           stress space. If the line is curved Rb is the tangent
!!
!! OUTPUT
!!     Depline (real) Double singular elasto-plastic infinitesimal
!!           constitutive matrix on a line. Size = nsigma x nsigma
!! REFERENCES:
!! [1] Johan Clausen, Lars Damkilde and Lars Andersen: "Efficient return
!!	  algorithms for associated plasticity with multiple yield planes", 
!!	  International Journal for Numerical Methods in Engineering, 2006,
!!	  volume 66, pages 1036-1059.
!! [2] Johan Clausen, Lars Damkilde and Lars Andersen: An efficient return
!!	  algorithm for non-associated plasticity with linear yield criteria
!!	  in principal stress space. Computers & Structures, 2007, volume 85,
!!	  pages 1795-1807.
!!------------------------------------------------------------------------
!! Johan Clausen
!! Esbjerg Department of Engineering
!! Aalborg University
!! July 2005
!!------------------------------------------------------------------------

!      integer(4), intent(in) :: nsigma
!      real(8), intent(in), dimension(nsigma,nsigma) :: D, Dinv
!      real(8), intent(in), dimension(3) :: Ra, Rb

!      real(8), intent(out) :: DepLine(nsigma,nsigma)
!      !-------------------------

!      real(8) den, Den1(3), Num(3,3), Dprin(3,3)
!      integer(4) i,j

!      !---------------------------------------------

!      forall(i = 1:3, j = 1:3) Num(i,j) = Ra(i)*Rb(j)

!      Den1 = matmul(Dinv(1:3,1:3),Rb)
!      den = dot_product(Ra,Den1) ! denomimator, Ra'*Dinv*Rb
!	
!      Dprin = Num/den

!      DepLine = 0.0_8
!      DepLine(1:3,1:3) = Dprin
!      if (nsigma > 3) then
!		DepLine(4,4) = D(4,4) ! In plane situations nsigma == 4
!		if (nsigma == 6) then ! General three dimensional stress state 
!			DepLine(5,5) = D(5,5)
!			DepLine(6,6) = D(6,6)
!		end if
!      end if
!      

!!-------------------------------------------------------------------------------
!      end subroutine FormDepLinePerfect !---------------------------------------
!!-------------------------------------------------------------------------------

				
!	pure subroutine PrinDirect(Sigma,nsigma,SigP,psi)
!!--------------------------------------------------------------------------
!! Calculates the principal directions of a principal stress state
!! Analytical expressions are used in the calculations.
!!
!! INPUT
!!  - Name -	-type,size -	-- description --
!!	Sigma	(real,ar(nsigma)Stress vector. Size is dependent on the type of stress
!!							  nsigma = 4 in plane situations and nsigma = 6 in 3D
!!							  The ordering of the components must be as follows:
!!							  Plane situations (plane stress, plane strain and axisymmetry):
!!							  Sigma = [sig_x sig_y sig_z tau_xy] (i.e. sig_z is out-of-plane)
!!							  General 3D:  Sigma = [sig_x sig_y sig_z tau_xy tau_xz tau_yz]
!!	nsigma	(integer,sc)	Number of stress components, i.e. the size of Sigma
!!	SigP	(real,ar(3))	Principal stresses in descending order
!!
!! OUTPUT
!!	psi		(real,ar(3,3))	Principal angle if nsigma = 4 or
!!							  vector of normalized eigenvectors (principal directions)
!!							  if nsigma = 6. When nsigma = 4 the angle is stored
!!							  in the upper left corner, i.e. psi(1,1).
!!
!!---------------------------------------------------------------------------
!! Johan Clausen
!! Esbjerg Institute of Technology
!! Aalborg University
!! November 2005
!!--------------------------------------------------------------------------
!	implicit none
!	integer(4), intent(in) :: nsigma
!	real(8), intent(in) :: Sigma(nsigma), SigP(3)

!	real(8), intent(out)  :: psi(3,3)

!!----- Internal variables -------------------------------------------------
!	real(8), parameter :: pi = 3.14159265358979323_8
!	real(8) tol, tol1, tol2, nl(3), normal(3,3), leng
!	real(8) denom1, denom2, denom3
!	integer(4) i, hj1, hj2, nr

!	interface
!		pure function Cross(A,B)
!			real(8), intent(in)	::	A(3), B(3)
!			real(8)				::  Cross(3)
!		end function Cross
!	end interface
!!----- Plane situations including axisymmetry -----------------------------

!	if (nsigma == 4) then
!        
!		if (Sigma(1) > Sigma(2)  .and.  Sigma(4) >= 0.0_8) then
!			psi = 0.5_8*atan(2*Sigma(4)/(Sigma(1)-Sigma(2)))
!		elseif (Sigma(1) < Sigma(2)  .and.  Sigma(4) >= 0.0_8) then
!			psi = 0.5_8*(pi - atan(2*Sigma(4)/(Sigma(2)-Sigma(1))))
!		elseif (Sigma(1) < Sigma(2)  .and.  Sigma(4) < 0.0_8) then
!			psi = 0.5_8*(atan(-2*Sigma(4)/(Sigma(2)-Sigma(1))) + pi)
!		elseif (Sigma(1) > Sigma(2)  .and.  Sigma(4) < 0.0_8) then
!			psi = 0.5_8*(2*pi - atan(-2*Sigma(4)/(Sigma(1)-Sigma(2))))
!		elseif (Sigma(1) == Sigma(2) .and. Sigma(4) > 0.0_8) then
!			psi = 0.25_8*pi
!		elseif (Sigma(1) == Sigma(2) .and. Sigma(4) < 0.0_8) then
!			psi = 0.75_8*pi
!		elseif (Sigma(1) == Sigma(2) .and. Sigma(4) == 0.0_8) then
!			psi = 0.0_8 ;
!		end if

!!------------------------------------------------------------------------


!!----- General 3D ----------------------------------------------------
!	elseif (nsigma == 6) then
!		tol1 = abs(1.0e-10_8*(SigP(1) - SigP(3))) ! tolerance value 1
!		tol2 = abs(1.0e-12_8*(SigP(1) + SigP(2) + SigP(3))) ! tolerance value 2
!		tol = maxval((/tol1,tol2,1.0e-13_8/))
!		
!		nl = 0.0_8
!		normal = 0.0_8
!!	-- Sigma is already expressed in principal coordinates ---
!		if (abs(Sigma(4)) < tol  .and.
!     &		abs(Sigma(5)) < tol  .and.
!     &		abs(Sigma(6)) < tol)		then
!			if (Sigma(1) >= Sigma(2) .and. Sigma(1) >= Sigma(3)) then
!				normal(1,1) = 1.0_8
!				if (Sigma(2) >= Sigma(3)) then
!					normal(2,2) = 1.0_8 ; normal(3,3) = 1.0_8
!				else
!					normal(3,2) = 1.0_8 ; normal(2,3) = 1.0_8
!				end if
!			elseif (Sigma(2) >= Sigma(1) .and. Sigma(2) >= Sigma(3))
!     &															 then
!				normal(2,1) = 1.0_8
!				if (Sigma(1) >= Sigma(3)) then
!					normal(1,2) = 1.0_8 ; normal(3,3) = 1.0_8
!				else
!					normal(3,2) = 1.0_8 ; normal(1,3) = 1.0_8
!				end if
!			elseif (Sigma(3) >= Sigma(1) .and. Sigma(3) >= Sigma(2))
!     &															 then
!				normal(3,1) = 1.0_8
!				if (Sigma(1) >= Sigma(2)) then
!					normal(1,2) = 1.0_8 ; normal(2,3) = 1.0_8
!				else
!					normal(2,2) = 1.0_8 ; normal(1,3) = 1.0_8
!				end if
!			end if
!!	-- Three distinct eigenvalues -----
!		elseif (abs(SigP(1)-SigP(2)) > tol  .and.
!     &			abs(SigP(2)-SigP(3)) > tol  .and.
!     &			abs(SigP(1)-SigP(3)) > tol)			then
!			do i = 1,2
!				denom1 = (Sigma(2)-SigP(i))*(Sigma(3)-SigP(i)) - Sigma(6)**2 ! Denominator if first eigenvector component is set to 1
!				denom2 = (Sigma(1)-SigP(i))*(Sigma(3)-SigP(i)) - Sigma(5)**2 ! Denominator if second eigenvector component is set to 1
!				denom3 = (Sigma(1)-SigP(i))*(Sigma(2)-SigP(i)) - Sigma(4)**2 ! Denominator if third eigenvector component is set to 1
!			
!				if (abs(denom1)  >=  max(abs(denom2),abs(denom3))) then
!				    nl(1) = 1.0_8
!				    nl(2) = (Sigma(4)*(SigP(i)-Sigma(3))+Sigma(6)*Sigma(5))/denom1
!				    nl(3) = (Sigma(5)*(SigP(i)-Sigma(2))+Sigma(6)*Sigma(4))/denom1
!			    elseif (abs(denom2)  >=  max(abs(denom1),abs(denom3))) then
!				    nl(1) = (Sigma(4)*(SigP(i)-Sigma(3))+Sigma(6)*Sigma(5))/denom2
!				    nl(2) = 1.0_8
!				    nl(3) = (Sigma(6)*(SigP(i)-Sigma(1))+Sigma(4)*Sigma(5))/denom2
!			    else
!				    nl(1) = (Sigma(5)*(SigP(i)-Sigma(2))+Sigma(6)*Sigma(4))/denom3
!				    nl(2) = (Sigma(6)*(SigP(i)-Sigma(1))+Sigma(4)*Sigma(5))/denom3
!				    nl(3) = 1.0_8
!                  end if
!				leng = sqrt(nl(1)**2 + nl(2)**2 + nl(3)**2)
!				normal(:,i) = nl/leng
!			end do
!			normal(:,3) = cross(Normal(:,1),Normal(:,2))
!			
!!	-- Three equal eigenvalues ----
!		elseif (abs(SigP(1)-SigP(2)) < tol  .and.
!     &			abs(SigP(2)-SigP(3)) < tol  .and.
!     &			abs(SigP(1)-SigP(3)) < tol)			then
!			
!			normal(1,1) = 1.0_8
!			normal(2,2) = 1.0_8
!			normal(3,3) = 1.0_8

!!	-- two equal eigenvalues -----
!		else
!			if (abs(SigP(1) - SigP(2)) <= tol) then
!				hj1 = 1 ; hj2 = 2 ; nr = 3
!			elseif (abs(SigP(1) - SigP(3)) <= tol) then
!				hj1 = 1 ; nr = 2 ; hj2 = 3
!			else 
!				nr = 1 ; hj1 = 2 ; hj2 = 3
!			end if
!            
!!       --- First principal direction -------------------------            
!			nl(1) = (Sigma(2)-SigP(nr)) * (Sigma(3)-SigP(nr))
!     &				- Sigma(6)**2
!			nl(2) = -Sigma(4) * (Sigma(3)-SigP(nr))
!     &				+ Sigma(6) * Sigma(5)
!			nl(3) = Sigma(4) * Sigma(6)
!     &				- (Sigma(2)-SigP(nr)) * Sigma(5)
!			leng = sqrt(nl(1)**2 + nl(2)**2 + nl(3)**2)
!			normal(:,nr) = nl/leng ;

!!       --- Second principal direction ------------------------            
!			if (normal(3,nr) > 0.0_8  .or. normal(3,nr) < 0.0_8) then
!				nl(1) = 1.0_8 ; nl(2) = 1.0_8 ;
!				nl(3) = -(normal(1,nr)+normal(2,nr))/normal(3,nr) ;
!			elseif (normal(2,nr) > 0.0_8  .or.
!     &								normal(2,nr) < 0.0_8) then
!				nl(1) = 1.0_8 ; nl(3) = 1.0_8 ;
!				nl(2) = -(normal(1,nr)+normal(3,nr))/normal(2,nr) ;
!			elseif (normal(1,nr) > 0.0_8  .or.
!     &								normal(1,nr) < 0.0_8) then
!				nl(2) = 1.0_8 ; nl(3) = 1.0_8
!				nl(1) = -(normal(2,nr)+normal(3,nr))/normal(1,nr) ;             
!			end if
!			leng = sqrt(nl(1)**2 + nl(2)**2 + nl(3)**2)
!			normal(:,hj1) = nl/leng ;
!!       --- Third principal direction -------------------------
!			normal(:,hj2) = cross(normal(:,nr),Normal(:,hj1))
!		end if
!		psi = normal
!	end if ! nsigma == ....

!!----------------------------------------------------------------------------------
!	end subroutine PrinDirect !--------------------------------------------------
!!----------------------------------------------------------------------------------

    
!      pure function Cross(A,B)
!!--------------------------------------------------------------------------
!! Calculates the cross product (vector product) between the two
!! three dimensional vectors A and B
!!
!! INPUT
!!	A, B	(real) Vectors with three components
!!
!! OUTPUT
!!	Cross	(real) A vector perpendicular to both A and B
!!
!!------------------------------------------------------------------------
!! Johan Clausen
!! Esbjerg Department of Engineering
!! Aalborg University
!! July 2005
!!------------------------------------------------------------------------
!	implicit none
!	real(8), intent(in) :: A(3), B(3)
!!-------------------------
!	real(8) Cross(3)
!!---------------------------------------------

!	Cross(1) = A(2)*B(3) - A(3)*B(2)
!	Cross(2) = A(3)*B(1) - A(1)*B(3)
!	Cross(3) = A(1)*B(2) - A(2)*B(1)
!!--------------------------------------------------------------------------------
!	end function Cross !-------------------------------------------------------
!!--------------------------------------------------------------------------------

!      pure subroutine TransMatrix(psi,nsigma,ouplP,A)
!!--------------------------------------------------------------------
!! Calculates transformation matrix, A, for stress transformation
!! based on either an angle (plane and axisymmetry) or direction
!! cosines (general 3D). The ordering of the shear stresses
!! must be specified, as it influences how A is built
!!
!! INPUT
!!  - Name- -type,size -   -- Description --
!!     psi   (real,ar(3,3))  Matrix of direction cosines in general
!!                   three-dimensional stress state or the rotation
!!                   angle otherwise. In 3D we have 
!!                   sig_tranformed = psi'*sig*psi, where sig
!!                   is a stress tensor
!!     nsigma  (integer,sc)  Number of stress components, i.e. size of A
!!     ouplP (integer,sc)  Only used in plan situations in which it is similar to
!!                   hoop but for the principal stress vector SigP, i.e.
!!                   it is the position of the out-of-plane principal stress.
!!
!! OUTPUT
!!     A Transformation matrix, with elements ordered according to hoop and ouplP
!!         A stress is transformed according to ("'" signifies matrix transpose)
!!         Sigma_transformed = inv(A')*Sigma and conversely
!!         Sigma = A'*Sigma_transformed
!!         A strain vector, Epsilon, is tranformed according to
!!         Epsilon_transformed = A*Epsilon and conversely
!!         Epsilon = inv(A)*Epsilon_transformed
!!         See e.g. [Cook, Malkus and Plesha, 1989] for further details
!!--------------------------------------------------------------------
!! Johan Clausen
!! Esbjerg Department of Engineering
!! Aalborg University
!! July 2005
!!--------------------------------------------------------------------

!      !-----------------
!      integer(4), intent(in) :: nsigma, ouplP
!      real(8), intent(in)    :: psi(3,3)

!      real(8), intent(out)   :: A(nsigma,nsigma)
!      !-----------------

!      real(8) sin_psi, cos_psi, sin_psi2, cos_psi2, sin_2psi
!      real(8), dimension(3,3) :: A1, A2, A3, A4, Asmall
!      real(8), dimension(3,4) :: Aint
!      integer(4), dimension(3):: Hj1, Hj2
!      integer(4) Ext(3,4), ExtP(4,3)
!      integer(4) i, j

!!---- Plane situations --------------------------------------------
!! Assumes that the stress components are ordered according to
!! Sigma = [sig_x sig_y sig_z tau_xy]
!!------------------------------------------------------------------ 
!      if (nsigma == 4) then
!        sin_psi  = sin(psi(1,1)) ; cos_psi =dcos(psi(1,1))
!        sin_psi2 = sin_psi**2    ; cos_psi2 = cos_psi**2
!        sin_2psi = sin(2.0_8*psi(1,1))
!      
!      
!        Asmall(1,1) = cos_psi2  ; Asmall(1,2) = sin_psi2
!        Asmall(2,1) = sin_psi2  ; Asmall(2,2) = cos_psi2
!        Asmall(3,1) = -sin_2psi ; Asmall(3,2) = sin_2psi
!        
!        Asmall(1,3) =  0.5_8*sin_2psi
!        Asmall(2,3) = -0.5_8*sin_2psi
!        Asmall(3,3) = cos_psi2 - sin_psi2

!        ExtP = 0.0_8
!        if (ouplP == 2) then
!          ExtP(1,1) = 1  ;  ExtP(3,2) = 1  ;  ExtP(4,3) = 1
!        elseif (ouplP == 1) then
!          ExtP(2,1) = 1  ;  ExtP(3,2) = 1  ;  ExtP(4,3) = 1
!        elseif (ouplP == 3) then
!          ExtP(1,1) = 1  ;  ExtP(2,2) = 1  ;  ExtP(4,3) = 1
!        end if
!				    ! IMPORTANT! This order of
!        Ext = 0.0_8 ! Ext assumes that the out-of-plane stress i in position 3!!
!        Ext(1,1) = 1 ; Ext(2,2) = 1 ; Ext(3,4) = 1

!        Aint = matmul(Asmall,Ext)
!        A = matmul(ExtP,Aint) ! Transformation matrix
!        A(ouplP,3) = 1.0_8

!!---- General three dimensional stress state ------------------------
!! Assumes that the stress components are ordered according to
!! Sigma = [sig_x sig_y sig_z tau_xy tau_xz tau_yz]
!!------------------------------------------------------------------ 
!      elseif (nsigma == 6) then
!      
!        Hj1 = (/1, 3, 2/)
!        Hj2 = (/2, 1, 3/)
!      
!        A1 = psi**2
!      
!        do i = 1,3
!          do j = 1,3
!            A2(i,j) = psi(i,Hj1(j))*psi(i,Hj2(j))
!            A3(i,j) = psi(Hj1(i),j)*psi(Hj2(i),j)
!            A4(i,j) = psi(Hj1(i),Hj1(j))*psi(Hj2(i),Hj2(j))
!     &              + psi(Hj2(i),Hj1(j))*psi(Hj1(i),Hj2(j))
!          end do
!        end do
!        A2 = 2*A2
!        
!        A(1:3,1:3) = A1 ; A(1:3,4:6) = A2
!        A(4:6,1:3) = A3 ; A(4:6,4:6) = A4
!		A = transpose(A)
!      end if ! nsigma == ....

!!----------------------------------------------------------------------------------
!      end subroutine TransMatrix !-------------------------------------------------
!!----------------------------------------------------------------------------------

!      pure subroutine DlinElas(E,nu,nsigma,D,Dinv)
!!--------------------------------------------------------------------------
!! Calculates linear elastic stiffness matrix D. If nsigma = 3, i.e.
!! only three stress components plane stress is assumed
!!
!! INPUT
!!     E   (real) Youngs module, module of elasticity
!!     nu    (real) Poisson's ratio
!!     nsigma  (int)  Number of stress components. 
!!             nsigma = 4 in plane problems in general
!!             if nsigma = 3 plane stress is assumed
!!             nsigma = 6 in general stress states
!! OUTPUT
!!     D   (real) Linear elastic stiffness matrix D
!!           Size is dependent on stress type
!!     Dinv  (real) Inverse of D, optional
!!------------------------------------------------------------------------
!! Johan Clausen
!! Esbjerg Department of Engineering
!! Aalborg University
!! July 2005
!!------------------------------------------------------------------------

!	implicit none

!      integer(4), intent(in) :: nsigma
!      real(8), intent(in) :: E, nu

!      real(8), intent(out) :: D(nsigma,nsigma)
!      real(8), intent(out), optional :: Dinv(nsigma,nsigma)

!!-------------------------

!      real(8) c, g

!      D = 0.0_8

!!----- Plane stress with three stress components -----------------------
!      if (nsigma == 3) then
!      
!        D(1,1) = 1.0_8 ; D(1,2) = nu
!        D(2,1) = nu    ; D(2,2) = 1.0_8
!                         D(3,3) = (1-nu)/2
!        D = E/(1-nu*nu) * D

!        if (present(Dinv)) then
!          Dinv = 0.0_8
!          Dinv(1,1) = 1.0_8 ; Dinv(1,2) = -nu
!          Dinv(2,1) = -nu   ; Dinv(2,2) = 1.0_8
!                              Dinv(3,3) = 2.0_8*(1+nu)
!          Dinv = Dinv/E
!        end if

!!----- Plane stress, strain or axisymmetry (four stress components) ---
!      elseif(nsigma == 4) then
!		D(1,1) = 1.0_8 - nu ; D(1,2) = nu         ; D(1,3) = nu
!		D(2,1) = nu         ; D(2,2) = 1.0_8 - nu ; D(2,3) = nu
!		D(3,1) = nu         ; D(3,2) = nu         ; D(3,3) = 1.0_8- nu
!		D(4,4) = (1-2.0_8*nu)/2
!        D = E/((1+nu)*(1-2.0_8*nu)) * D

!        if (present(Dinv)) then
!          Dinv = 0.0_8
!          Dinv(1,1) = 1.0_8 ; Dinv(1,2) = -nu     ; Dinv(1,3) = -nu
!          Dinv(2,1) = -nu   ; Dinv(2,2) = 1.0_8   ; Dinv(2,3) = -nu
!          Dinv(3,1) = -nu   ; Dinv(3,2) = -nu     ; Dinv(3,3) = 1.0_8
!          Dinv(4,4) = 2.0_8*(1+nu)
!          Dinv = Dinv/E
!        end if
!!----- General stress state ---------------------------------------------
!      elseif(nsigma == 6) then
!        c = E/((1+nu)*(1-2.0_8*nu))
!        g = E/(2.0_8*(1+nu))

!        D(1,1) = (1-nu)*c ; D(1,2) = nu*c   ; D(1,3) = D(1,2)
!        D(2,1) = D(1,2)   ; D(2,2) = D(1,1) ; D(2,3) = D(1,2)
!        D(3,1) = D(1,2)   ; D(3,2) = D(1,2) ; D(3,3) = D(1,1)
!        D(4,4) = g
!        D(5,5) = g
!        D(6,6) = g

!        if (present(Dinv)) then
!          Dinv = 0.0_8
!          Dinv(1,1) = 1.0_8 ; Dinv(1,2) = -nu   ; Dinv(1,3) = -nu
!          Dinv(2,1) = -nu   ; Dinv(2,2) = 1.0_8 ; Dinv(2,3) = -nu
!          Dinv(3,1) = -nu   ; Dinv(3,2) = -nu   ; Dinv(3,3) = 1.0_8
!          Dinv(4,4) = 2.0_8*(1+nu)
!          Dinv(5,5) = 2.0_8*(1+nu)
!          Dinv(6,6) = 2.0_8*(1+nu)
!          Dinv = Dinv/E
!        end if
!      end if ! nsigma = ...

!!--------------------------------------------------------------------------------
!      end subroutine  DlinElas !-------------------------------------------------
!!--------------------------------------------------------------------------------

      
!	pure subroutine Invariants(Sigma,nsigma,S,I1,J2,J3,lode,sin3lode)
!!--------------------------------------------------------------------------
!! Calculates the deviator stress vector of a stress vector Sigma, along
!! with several stress invariants.
!!
!! INPUT
!!  - Name -	-- description --
!!	Sigma	(real,ar(nsigma): Stress vector. Size is dependent on the type of stress
!!			  as given by nsigma.
!!			  nsigma = 3: It is assumed that Sigma is principal stresses,
!!						  Sigma = [SigP_1 SigP_2 SigP_3]
!!			  nsigma = 4: Plane situations, i.e. plane stress, plane strain
!!						  and axisymmetry. The ordering of the stresses must
!!						  follow Sigma = [sig_x sig_y sig_z tau_xy] , i.e.
!!						  sig_z is the out-of-plane stress.
!!			  nsigma = 6  General 3D. The stresses must be ordered according
!!						  to Sigma = [sig_x sig_y sig_z tau_xy tau_xz tau_yz]
!!	nsigma (integer,sc): Number of stress components, i.e. the size of Sigma.
!!
!! OUTPUT
!!	S		(real, ar(nsigma)): Deviator stress vector. Same size as Sigma
!!	I1	(real,scalar): First stress invariant.
!!	J2	(real,scalar): Second deviator stress invariant.
!!	J3	(real,scalar): Third deviator stress invariant.
!!	lode	(real,scalar): Lode angle in radians. In the 
!!			  range -pi/6 <= lode <= pi/6,    (pi/6 = 30 deg)
!!	sin3lode (real,scalar): sin(3*lode)
!!
!!---------------------------------------------------------------------------
!! Johan Clausen
!! Deparment of Civil Engineering
!! Aalborg University
!! May 2008
!!--------------------------------------------------------------------------
!	implicit none
!	integer(4), intent(in) :: nsigma
!	real(8), intent(in) :: Sigma(nsigma)

!	real(8), intent(out) :: S(nsigma), I1, J2
!	real(8), intent(out) :: J3, lode, sin3lode

!!--------------------------------------------------------------------------

!	I1 = sum(Sigma(1:3))
!	if (nsigma == 3) then ! Sigma is the principal stresses
!		S = Sigma - I1*0.333333333333333333_8 ! Deviator stress vector
!		J2 = 0.5_8*sum(S**2)
!		J3 = S(1)*S(2)*S(3)
!    
!	elseif (nsigma == 4) then ! Plane situation, including axisymmetry
!		S(1:3) = Sigma(1:3) - I1*0.333333333333333333_8 ! Deviator stress vector
!		S(4) = Sigma(4)
!		J2 = 0.5_8 * sum(S(1:3)**2) + S(4)**2
!    
!		J3 = (S(1)**3 + S(2)**3 + S(3)**3 +3.0_8*S(4)**2 *(S(1)+S(2)))
!     &								*0.3333333333333333

!	elseif (nsigma == 6) then ! General three-dimensional stress state
!		S(1:3) = Sigma(1:3) - I1*0.33333333333333333_8 ! Deviator stress vector
!		S(4:6) = Sigma(4:6)

!		J2 = 0.5_8 * sum(S(1:3)**2) + sum(S(4:6)**2)

!		J3 = ( sum(S(1:3)**3) + 6.0_8*S(4)*S(5)*S(6) +
!     &		 3.0_8*(S(1)*(S(4)**2 + S(5)**2) +
!     &			    S(2)*(S(4)**2 + S(6)**2) +
!     &				S(3)*(S(5)**2 + S(6)**2)) )*0.3333333333333333_8
!    
!	end if

!	if (J2 > 0.0_8) then
!		sin3lode = -2.598076211353316_8 * J3 / J2**(1.5_8)
!	else
!		sin3lode = 0.0_8
!	endif

!	if (sin3lode <= -1.0_8) then
!		lode = -0.5235987755982988_8 ! -pi/6
!	elseif (sin3lode >= 1.0_8) then
!		lode = +0.5235987755982988_8 ! +pi/6
!	else
!		lode = asin(sin3lode)/3.0_8
!	endif

!!----------------------------------------------------------------------------------
!	end subroutine Invariants !--------------------------------------------------
!!----------------------------------------------------------------------------------
