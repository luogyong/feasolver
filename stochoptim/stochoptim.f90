!=======================================================================
! StochOptim provides user friendly functions to solve optimization
! problems using stochastic algorithms.
!
! Created by
!     Keurfon Luu <keurfon.luu@mines-paristech.fr>
!     MINES ParisTech - Centre de Géosciences
!     PSL - Research University
!=======================================================================

module stochoptim
  USE, INTRINSIC :: ISO_FORTRAN_ENV
  use forlab, only: IPRE, RPRE, progress_perc, num2str, randu, randn, &
    randi, randperm, norm, argsort, zeros, ones, eye, prctile, mean, &
    median, diag, signum, argsort, triu, eig, savebin
!dec$ if defined( do_mpi )
  use forlab, only: mpi_rpre
  use mpi
!dec$endif

  implicit none
  
  private
  public :: optifunc, EvolutionaryKeywords, MonteCarloKeywords, &
    Evolutionary, MonteCarlo,EA_ITER

    
    INTEGER::EA_ITER=0
    
  abstract interface
    real(kind = RPRE) function optifunc(x)
      use forlab, only : RPRE
      import
      real(kind = RPRE), dimension(:), intent(in) :: x
    end function optifunc
  end interface

  type EvolutionaryKeywords
    integer(kind = IPRE) :: popsize = 30, max_iter = 1000, verbose = 0
    real(kind = RPRE) :: eps1 = 1d-4, eps2 = 1d-3, w = 0.7298, c1 = 1.49618, &
      c2 = 1.49618, gamma = 1., F = 0.5, CR = 0.1, sigma = 0.5, mu_perc = 0.5,&
      eps3=1.d-3
    logical :: constrain = .true., snap = .false.
    character(len = 5) :: solver, strategy
  contains
    procedure, private, pass :: print_EvolutionaryKeywords
    generic, public :: print => print_EvolutionaryKeywords
  end type EvolutionaryKeywords

  type MonteCarloKeywords
    integer(kind = IPRE) :: max_iter = 1000, verbose = 0
    real(kind = RPRE) :: perc
    logical :: constrain = .true.
    real(kind = RPRE), dimension(:), allocatable :: stepsize, xstart
    character(len = 4) :: sampler
  contains
    procedure, private, pass :: print_MonteCarloKeywords
    generic, public :: print => print_MonteCarloKeywords
  end type MonteCarloKeywords

  type Evolutionary
    procedure(optifunc), pointer, nopass :: func
    real(kind = RPRE), dimension(:), allocatable :: lower, upper
    integer(kind = IPRE) :: n_dim, popsize = 30, max_iter = 1000, verbose = 0
    real(kind = RPRE) :: eps1 = 1e-4, eps2 = 1e-4,eps3=1e-4
    logical :: constrain = .true., snap = .false.
    character(len = 5) :: solver, strategy

    real(kind = RPRE), dimension(:), allocatable :: xopt
    real(kind = RPRE) :: gfit
    integer(kind = IPRE) :: iflag, n_iter, n_eval, n_restart
    real(kind = RPRE), dimension(:,:,:), allocatable :: models
    real(kind = RPRE), dimension(:,:), allocatable :: energy, means !,bestmodel(:,:)
    real(kind = RPRE), dimension(:), allocatable :: fitness, mu_scale, std_scale
    real(kind = RPRE) :: w, c1, c2, gamma, F, CR, sigma, mu_perc

    integer(kind = IPRE) :: mpi_rank = 0, mpi_size = 1, mpi_ierr
  contains
    procedure, private, pass :: print_Evolutionary_parameters, print_Evolutionary_result, &
      print_Evolutionary, flag, standardize1_evolutionary, standardize2_evolutionary, &
      unstandardize1_evolutionary, unstandardize2_evolutionary, init_models_evolutionary, &
      eval_models, constrain_de, constrain_cpso, de_mutation, de, cpso, cmaes, &
      save_models_evolutionary
    generic, private :: init_models => init_models_evolutionary
    generic, private :: standardize => standardize1_evolutionary, standardize2_evolutionary
    generic, private :: unstandardize => unstandardize1_evolutionary, unstandardize2_evolutionary
    generic, public :: print_parameters => print_Evolutionary_parameters
    generic, public :: print_result => print_Evolutionary_result
    generic, public :: print => print_Evolutionary
    generic, public :: save_models => save_models_evolutionary
    procedure, public, pass :: optimize
  end type Evolutionary

  type MonteCarlo
    procedure(optifunc), pointer, nopass :: func
    real(kind = RPRE), dimension(:), allocatable :: lower, upper
    integer(kind = IPRE) :: n_dim, n_dim_per_iter, max_iter = 1000, verbose = 0
    logical :: constrain = .true.
    character(len = 4) :: sampler

    real(kind = RPRE), dimension(:), allocatable :: xopt
    real(kind = RPRE) :: gfit, acceptance_ratio, perc
    real(kind = RPRE), dimension(:,:), allocatable :: models
    real(kind = RPRE), dimension(:), allocatable :: energy
    real(kind = RPRE), dimension(:), allocatable :: mu_scale, std_scale
  contains
    procedure, private, pass :: print_MonteCarlo_parameters, print_MonteCarlo_result, &
      print_MonteCarlo, standardize_montecarlo, unstandardize1_montecarlo, unstandardize2_montecarlo, &
      init_models_montecarlo, in_search_space, mcmc, save_models_montecarlo
    generic, private :: init_models => init_models_montecarlo
    generic, private :: standardize => standardize_montecarlo
    generic, private :: unstandardize => unstandardize1_montecarlo, unstandardize2_montecarlo
    generic, public :: print_parameters => print_MonteCarlo_parameters
    generic, public :: print_result => print_MonteCarlo_result
    generic, public :: print => print_MonteCarlo
    generic, public :: save_models => save_models_montecarlo
    procedure, public, pass :: sample
  end type MonteCarlo

  interface EvolutionaryKeywords
    module procedure init_EvolutionaryKeywords
  end interface EvolutionaryKeywords

  interface MonteCarloKeywords
    module procedure init_MonteCarloKeywords
  end interface MonteCarloKeywords

  interface Evolutionary
    module procedure init_Evolutionary
  end interface Evolutionary

  interface MonteCarlo
    module procedure init_MonteCarlo
  end interface MonteCarlo

contains

  function init_EvolutionaryKeywords(popsize, max_iter, eps1, eps2, eps3,constrain, snap, solver, &
                                     w, c1, c2, gamma, F, CR, strategy, sigma, mu_perc, verbose) result(evo_kws)
    type(EvolutionaryKeywords) :: evo_kws
    integer(kind = IPRE), intent(in), optional :: popsize, max_iter, verbose
    real(kind = RPRE), intent(in), optional :: eps1, eps2, eps3,w, c1, c2, gamma, &
      F, CR, sigma, mu_perc
    logical, intent(in), optional :: constrain, snap
    character(len = *), intent(in), optional :: solver, strategy

    if ( present(popsize) ) evo_kws % popsize = popsize
    if ( present(max_iter) ) evo_kws % max_iter = max_iter
    if ( present(eps1) ) evo_kws % eps1 = eps1
    if ( present(eps2) ) evo_kws % eps2 = eps2
    if ( present(eps3) ) evo_kws % eps2 = eps3
    if ( present(constrain) ) evo_kws % constrain = constrain
    if ( present(snap) ) evo_kws % snap = snap
    if ( present(w) ) evo_kws % w = w
    if ( present(c1) ) evo_kws % c1 = c1
    if ( present(c2) ) evo_kws % c2 = c2
    if ( present(gamma) ) evo_kws % gamma = gamma
    if ( present(F) ) evo_kws % F = F
    if ( present(CR) ) evo_kws % CR = CR
    if ( present(sigma) ) evo_kws % sigma = sigma
    if ( present(mu_perc) ) evo_kws % mu_perc = mu_perc
    if ( present(verbose) ) evo_kws % verbose = verbose
    if ( present(solver) ) then
      evo_kws % solver = trim(solver)
    else
      evo_kws % solver = "cpso"
    end if
    if ( present(strategy) ) then
      evo_kws % strategy = trim(strategy)
    else
      evo_kws % strategy = "rand1"
    end if
    return
  end function init_EvolutionaryKeywords

  subroutine print_EvolutionaryKeywords(self,UNIT)
    class(EvolutionaryKeywords), intent(inout) :: self
    character(len = 18) :: attributes
    INTEGER,INTENT(IN),OPTIONAL::UNIT    
    INTEGER::UNIT1
    
    UNIT1=OUTPUT_UNIT
    IF(PRESENT(UNIT)) UNIT1=UNIT

    attributes = "solver:"
    WRITE(UNIT1,*)  adjustr(attributes) // " '" // trim(self % solver) // "'"
    attributes = "parameters:"
    select case(trim(self % solver))
    case("pso")
      WRITE(UNIT1,*)  adjustr(attributes) &
        // " w = " // num2str(self % w, "(F5.2)") &
        // ", c1 = " // num2str(self % c1, "(F5.2)") &
        // ", c2 = " // num2str(self % c2, "(F5.2)")
    case("cpso")
      WRITE(UNIT1,*)  adjustr(attributes) &
        // " w = " // num2str(self % w, "(F5.2)") &
        // ", c1 = " // num2str(self % c1, "(F5.2)") &
        // ", c2 = " // num2str(self % c2, "(F5.2)") &
        // ", gamma = " // num2str(self % gamma, "(F5.2)")
    case("de")
      WRITE(UNIT1,*)  adjustr(attributes) &
        // " F = " // num2str(self % F, "(F5.2)") &
        // ", CR = " // num2str(self % CR, "(F5.2)") &
        // ", strategy = '" // trim(self % strategy) // "'"
    case("cmaes")
      WRITE(UNIT1,*)  adjustr(attributes) &
        // " sigma = " // num2str(self % sigma, "(F5.2)") &
        // ", mu_perc = " // num2str(self % mu_perc, "(F5.2)")
    end select
    attributes = "popsize:"
    WRITE(UNIT1,*)  adjustr(attributes) // " " // num2str(self % popsize)
    attributes = "max_iter:"
    WRITE(UNIT1,*)  adjustr(attributes) // " " // num2str(self % max_iter)
    return
  end subroutine print_EvolutionaryKeywords

  function init_MonteCarloKeywords(max_iter, perc, constrain, sampler, &
                                   stepsize, xstart, verbose) result(mc_kws)
    type(MonteCarloKeywords) :: mc_kws

    integer(kind = IPRE), intent(in), optional :: max_iter, verbose
    real(kind = RPRE), intent(in), optional :: perc
    logical, intent(in), optional :: constrain
    character(len = *), intent(in), optional :: sampler
    real(kind = RPRE), dimension(:), intent(in), optional :: stepsize, xstart

    if ( present(max_iter) ) mc_kws % max_iter = max_iter
    if ( present(perc) ) mc_kws % perc = perc
    if ( present(constrain) ) mc_kws % constrain = constrain
    if ( present(stepsize) ) mc_kws % stepsize = stepsize
    if ( present(xstart) ) mc_kws % xstart = xstart
    if ( present(verbose) ) mc_kws % verbose = verbose
    if ( present(sampler) ) then
      mc_kws % sampler = trim(sampler)
    else
      mc_kws % sampler = "mcmc"
    end if
    return
  end function init_MonteCarloKeywords

  subroutine print_MonteCarloKeywords(self)
    class(MonteCarloKeywords), intent(inout) :: self
    character(len = 18) :: attributes

    attributes = "sampler:"
    print *, adjustr(attributes) // " '" // trim(self % sampler) // "'"
    attributes = "max_iter:"
    print *, adjustr(attributes) // " " // num2str(self % max_iter)
    return
  end subroutine print_MonteCarloKeywords

  function init_Evolutionary(func, lower, upper, popsize, max_iter, eps1, eps2,eps3, &
                             constrain, snap, verbose) result(ea)
    type(Evolutionary) :: ea
    procedure(optifunc) :: func
    real(kind = RPRE), dimension(:), intent(in) :: lower, upper
    integer(kind = IPRE), intent(in), optional :: popsize, max_iter, verbose
    real(kind = RPRE), intent(in), optional :: eps1, eps2,eps3
    logical, intent(in), optional :: constrain, snap

    ea % func => func
    if ( size(lower) .ne. size(upper) ) then
      print *, "Error: lower and upper must have the same length"
      stop
    else if ( any(lower .gt. upper) ) then
      print *, "Error: upper must be greater than lower"
      stop
    else
      ea % lower = lower
      ea % upper = upper
      ea % n_dim = size(lower)
    end if

    if ( present(popsize) ) then
      if ( popsize .lt. 2 ) then
        print *, "Error: popsize must be an integer > 1"
        stop
      else
        ea % popsize = popsize
      end if
    end if

    if ( present(max_iter) ) then
      if ( max_iter .lt. 1 ) then
        print *, "Error: max_iter must be a positive integer"
        stop
      else
        ea % max_iter = max_iter
      end if
    end if

    if ( present(verbose) ) then
      if ( verbose .ne. 0 .and. verbose .ne. 1 ) then
        print *, "Error: verbose must be either 0 or 1"
        stop
      else
        ea % verbose = verbose
      end if
    end if

    if ( present(eps1) ) then
      if ( eps1 .lt. 0. ) then
        print *, "Error: eps1 must be positive"
        stop
      else
        ea % eps1 = eps1
      end if
    end if

    if ( present(eps2) ) ea % eps2 = eps2
    if ( present(eps3) ) ea % eps3 = eps3
    if ( present(constrain) ) ea % constrain = constrain
    if ( present(snap) ) ea % snap = snap
    return
  end function init_Evolutionary

  subroutine print_Evolutionary_parameters(self,UNIT)
    class(Evolutionary), intent(inout) :: self
    character(len = 18) :: attributes
    type(EvolutionaryKeywords) :: evo_kws
    INTEGER,INTENT(IN),OPTIONAL::UNIT    
    INTEGER::UNIT1
    
    UNIT1=OUTPUT_UNIT
    IF(PRESENT(UNIT)) UNIT1=UNIT

    if ( self % mpi_rank .eq. 0 ) then
      if ( allocated(self % xopt) ) then
        evo_kws = EvolutionaryKeywords(&
          solver = self % solver, &
          popsize = self % popsize, &
          max_iter = self % max_iter, &
          w = self % w, &
          c1 = self % c1, &
          c2 = self % c2, &
          gamma = self % gamma, &
          F = self % F, &
          CR = self % CR, &
          strategy = self % strategy, &
          sigma = self % sigma, &
          mu_perc = self % mu_perc)
        call evo_kws % print(UNIT1)
        attributes = "n_dim:"
        print *, adjustr(attributes) // " " // num2str(self % n_dim)
      end if
    end if
    return
  end subroutine print_Evolutionary_parameters

  subroutine print_Evolutionary_result(self,UNIT)
    class(Evolutionary), intent(inout) :: self
    integer(kind = IPRE) :: i
    character(len = 18) :: attributes
    INTEGER,INTENT(IN),OPTIONAL::UNIT    
    INTEGER::UNIT1
    
    UNIT1=OUTPUT_UNIT
    IF(PRESENT(UNIT)) UNIT1=UNIT
    
    if ( self % mpi_rank .eq. 0 ) then
      if ( allocated(self % xopt) ) then
        attributes = "solution:"
        WRITE(UNIT1,*)  adjustr(attributes)
        do i = 1, self % n_dim
          WRITE(UNIT1,*)  char(9), char(9), char(9), self % xopt(i)
        end do
        attributes = "fitness:"
        WRITE(UNIT1,*)  adjustr(attributes) // " " // num2str(self % gfit)
        attributes = "n_iter:"
        WRITE(UNIT1,*)  adjustr(attributes) // " " // num2str(self % n_iter)
        attributes = "n_eval:"
        WRITE(UNIT1,*)  adjustr(attributes) // " " // num2str(self % n_eval)
        if ( trim(self % solver) .eq. "cpso" ) then
          attributes = "n_restart:"
          WRITE(UNIT1,*)  adjustr(attributes) // " " // num2str(self % n_restart)
        end if
        attributes = "flag:"
        WRITE(UNIT1,*)  adjustr(attributes) // " " // self % flag()
      end if
    end if
    return
  end subroutine print_Evolutionary_result

  subroutine print_Evolutionary(self,UNIT)
    class(Evolutionary), intent(inout) :: self
    INTEGER,INTENT(IN),OPTIONAL::UNIT    
    INTEGER::UNIT1
    
    UNIT1=OUTPUT_UNIT
    IF(PRESENT(UNIT)) UNIT1=UNIT
    
    call self % print_parameters(UNIT1)
    call self % print_result(UNIT1)
    return
  end subroutine print_Evolutionary

  function init_MonteCarlo(func, lower, upper, max_iter, perc, constrain, verbose) result (mc)
    type(MonteCarlo) :: mc
    procedure(optifunc) :: func
    real(kind = RPRE), dimension(:), intent(in) :: lower, upper
    integer(kind = IPRE), intent(in), optional :: max_iter, verbose
    real(kind = RPRE), intent(in), optional :: perc
    logical, intent(in), optional :: constrain

    mc % func => func
    if ( size(lower) .ne. size(upper) ) then
      print *, "Error: lower and upper must have the same length"
      stop
    else if ( any(lower .gt. upper) ) then
      print *, "Error: upper must be greater than lower"
      stop
    else
      mc % lower = lower
      mc % upper = upper
      mc % n_dim = size(lower)
    end if

    if ( present(max_iter) ) then
      if ( max_iter .lt. 1 ) then
        print *, "Error: max_iter must be a positive integer"
        stop
      else
        mc % max_iter = max_iter
      end if
    end if

    if ( present(perc) ) then
      if ( perc .lt. 0. .or. perc .gt. 1. ) then
        print *, "Error: perc must be a scalar in [ 0, 1 ]"
        stop
      else
        mc % n_dim_per_iter = max(1, int(perc * mc % n_dim))
      end if
    else
      mc % n_dim_per_iter = mc % n_dim
    end if

    if ( present(verbose) ) then
      if ( verbose .ne. 0 .and. verbose .ne. 1 ) then
        print *, "Error: verbose must be either 0 or 1"
        stop
      else
        mc % verbose = verbose
      end if
    end if

    if ( present(constrain) ) mc % constrain = constrain
    return
  end function init_MonteCarlo

  subroutine print_MonteCarlo_parameters(self)
    class(MonteCarlo), intent(inout) :: self
    character(len = 18) :: attributes
    type(MonteCarloKeywords) :: mc_kws

    if ( allocated(self % xopt) ) then
      mc_kws = MonteCarloKeywords(&
        sampler = self % sampler, &
        max_iter = self % max_iter)
      call mc_kws % print()
      attributes = "n_dim:"
      print *, adjustr(attributes) // " " // num2str(self % n_dim)
    end if
    return
  end subroutine print_MonteCarlo_parameters

  subroutine print_MonteCarlo_result(self)
    class(MonteCarlo), intent(inout) :: self
    integer(kind = IPRE) :: i
    character(len = 18) :: attributes

    if ( allocated(self % xopt) ) then
      attributes = "solution:"
      print *, adjustr(attributes)
      do i = 1, self % n_dim
        print *, char(9), char(9), char(9), self % xopt(i)
      end do
      attributes = "fitness:"
      print *, adjustr(attributes) // " " // num2str(self % gfit)
      attributes = "acceptance_ratio:"
      print *, adjustr(attributes) // " " // num2str(self % acceptance_ratio, "(F6.2)")
    end if
    return
  end subroutine print_MonteCarlo_result

  subroutine print_MonteCarlo(self)
    class(MonteCarlo), intent(inout) :: self

    call self % print_parameters()
    call self % print_result()
    return
  end subroutine print_MonteCarlo

  function flag(self) result(str)
    character(len = :), allocatable :: str
    class(Evolutionary), intent(inout) :: self

    select case(self % iflag)
    case(-1)
      str = "maximum number of iterations is reached"
    case(0)
      str = "best individual position changes less than eps1 (" // num2str(self % eps1) // ")"
    case(1)
      str = "fitness is lower than threshold eps2 (" // num2str(self % eps2) // ")"
    case(2)
      str = "NoEffectAxis"
    case(3)
      str = "NoEffectCoord"
    case(4)
      str = "ConditionCov"
    case(5)
      str = "EqualFunValues"
    case(6)
      str = "TolXUp"
    case(7)
      str = "TolFun"
    case(8)
      str = "TolX"
    case(101)
      str = "fitness is improved less than eps3("//num2str(self % eps3)//",[(gold-gnew)/gold]) every 100 iter." 
    case(201)
      str =  "fitness is improved less than eps3("//num2str(self % eps3)//",[(gold-gnew)/gold]) every 200 iter."
    end select
    return
  end function flag

  subroutine optimize(self, solver, xstart, w, c1, c2, gamma, &
                      F, CR, strategy, sigma, mu_perc)
    class(Evolutionary), intent(inout) :: self
    character(len = *), intent(in), optional :: solver, strategy
    real(kind = RPRE), intent(in), optional :: w, c1, c2, gamma, F, CR, sigma, mu_perc
    real(kind = RPRE), dimension(:,:), intent(in), optional :: xstart
    real(kind = RPRE) :: opt_w = 0.7298, opt_c1 = 1.49618, opt_c2 = 1.49618, opt_gamma = 1., &
      opt_F = 0.5, opt_CR = 0.1, opt_sigma = 0.5, opt_mu_perc = 0.5
    real(kind = RPRE), dimension(:,:), allocatable :: opt_xstart
    character(len = :), allocatable :: opt_solver, opt_strategy
    integer(kind = 4) :: master = 0, num_procs, rank, ierr

    ! Check inputs
    if ( present(solver) ) then
      if ( trim(solver) .ne. "pso" .and. trim(solver) .ne. "cpso" &
        .and. trim(solver) .ne. "de" .and. trim(solver) .ne. "cmaes" ) then
        print *, "Error: solver must be either 'cpso', 'pso', 'de' or 'cmaes'"
        stop
      else
        opt_solver = trim(solver)
      end if
    else
      opt_solver = "cpso"
    end if
    if ( present(xstart) ) then
      select case(opt_solver)
      case("cmaes")
        if ( size(xstart, 1) .ne. 1 &
          .or. size(xstart, 2) .ne. self % n_dim ) then
          print *, "Error: xstart must be of shape [ 1, " // num2str(self % n_dim) &
            // " ] for '" // opt_solver // "'"
          stop
        else
          opt_xstart = xstart
        end if
      case default
        if ( size(xstart, 1) .ne. self % popsize &
          .or. size(xstart, 2) .ne. self % n_dim ) then
          print *, "Error: xstart must be of shape [ " // num2str(self % popsize) &
            // ", " // num2str(self % n_dim) // " ] for '" // opt_solver // "'"
          stop
        else
          opt_xstart = xstart
        end if
      end select
    end if
    if ( present(w) ) then
      if ( w .lt. 0. .or. w .gt. 1. ) then
        print *, "Error: w must be a scalar in [ 0, 1 ]"
        stop
      else
        opt_w = w
      end if
    end if
    if ( present(c1) ) then
      if ( ( c1 .lt. 0. ) .or. ( c1 .gt. 4. ) ) then
        print *, "Error: c1 must be a scalar in [ 0, 4 ]"
        stop
      else
        opt_c1 = c1
      end if
    end if
    if ( present(c2) ) then
      if ( c2 .lt. 0. .or. c2 .gt. 4. ) then
        print *, "Error: c2 must be a scalar in [ 0, 4 ]"
        stop
      else
        opt_c2 = c2
      end if
    end if
    if ( present(gamma) ) then
      if ( gamma .lt. 0. .or. gamma .gt. 2. ) then
        print *, "Error: gamma must be a scalar in [ 0, 2 ]"
        stop
      else
        opt_gamma = gamma
      end if
    end if
    if ( present(CR) ) then
      if ( ( CR .lt. 0. ) .or. ( CR .gt. 1. ) ) then
        print *, "Error: CR must be a scalar in [ 0, 1 ]"
        stop
      else
        opt_CR = CR
      end if
    end if
    if ( present(F) ) then
      if ( F .lt. 0. .or. F .gt. 2. ) then
        print *, "Error: F must be a scalar in [ 0, 2 ]"
        stop
      else
        opt_F = F
      end if
    end if
    if ( present(strategy) ) then
      if ( trim(strategy) .ne. "rand1" .and. trim(strategy) .ne. "rand2" &
        .and. trim(strategy) .ne. "best1" .and. trim(strategy) .ne. "best2" ) then
        print *, "Error: solver must be either 'rand1', 'rand2', 'best1' or 'best2'"
        stop
      else
        opt_strategy = trim(strategy)
      end if
    else
      opt_strategy = "rand1"
    end if
    if ( present(sigma) ) then
      if ( sigma .le. 0. ) then
        print *, "Error: sigma must be positive"
        stop
      else
        opt_sigma = sigma
      end if
    end if
    if ( present(mu_perc) ) then
      if ( mu_perc .le. 1. / self % popsize .or. mu_perc .gt. 1. ) then
        print *, "Error: mu_perc must be a scalar in [ " // num2str(1. / self % popsize) // ", 1 ]"
        stop
      else
        opt_mu_perc = mu_perc
      end if
    end if

    ! Check population size
    select case(opt_solver)
    case("de")
      if ( self % popsize .le. 3 .and. ( opt_strategy .eq. "rand1" .or. opt_strategy .eq. "rand2" ) ) then
        print *, "Error: popsize cannot be lower than 4 for DE '" // opt_strategy // "' strategy"
        stop
      else if ( self % popsize .le. 4 .and. opt_strategy .eq. "best2" ) then
        print *, "Error: popsize cannot be lower than 5 for DE '" // opt_strategy // "' strategy"
        stop
      else if ( self % popsize .le. 5 .and. opt_strategy .eq. "rand2" ) then
        print *, "Error: popsize cannot be lower than 6 for DE '" // opt_strategy // "' strategy"
        stop
      end if
    case("cmaes")
      if ( self % popsize .le. 3 ) then
        print *, "Error: popsize cannot be lower than 4 for CMA-ES"
        stop
      end if
    end select

    ! Initialize
    self % solver = opt_solver
    self % n_iter = 0
    self % n_eval = 0
    self % n_restart = 0
    self % mu_scale = 0.5 * (self % upper + self % lower)
    self % std_scale = 0.5 * (self % upper - self % lower)
!dec$ if defined(do_mpi)
    call mpi_comm_size(mpi_comm_world, self % mpi_size, self % mpi_ierr)
    call mpi_comm_rank(mpi_comm_world, self % mpi_rank, self % mpi_ierr)
!dec$ endif

    ! Solve
    select case(opt_solver)
    case("pso")
      self % w = opt_w
      self % c1 = opt_c1
      self % c2 = opt_c2
      call self % cpso(w = opt_w, c1 = opt_c1, c2 = opt_c2, gamma = real(0., RPRE), xstart = opt_xstart)
    case("cpso")
      self % w = opt_w
      self % c1 = opt_c1
      self % c2 = opt_c2
      self % gamma = opt_gamma
      call self % cpso(w = opt_w, c1 = opt_c1, c2 = opt_c2, gamma = opt_gamma, xstart = opt_xstart)
    case("de")
      self % CR = opt_CR
      self % F = opt_F
      self % strategy = opt_strategy
      call self % de(CR = opt_CR, F = opt_F, strategy = opt_strategy, xstart = opt_xstart)
    case("cmaes")
      self % sigma = opt_sigma
      self % mu_perc = opt_mu_perc
      call self % cmaes(sigma = opt_sigma, mu_perc = opt_mu_perc, xstart = opt_xstart)
    end select
    return
  end subroutine optimize

  function standardize1_evolutionary(self, models) result(out)
    real(kind = RPRE), dimension(:), allocatable :: out
    class(Evolutionary), intent(inout) :: self
    real(kind = RPRE), dimension(:), intent(in) :: models
    real(kind = RPRE), dimension(:), allocatable :: std_scale
    logical, dimension(:), allocatable :: mask

    mask = self % std_scale .ne. 0.
    std_scale = merge(self % std_scale, real(1., RPRE), mask)
    out = ( models - self % mu_scale ) / std_scale
    return
  end function standardize1_evolutionary

  function standardize2_evolutionary(self, models) result(out)
    real(kind = RPRE), dimension(:,:), allocatable :: out
    class(Evolutionary), intent(inout) :: self
    real(kind = RPRE), dimension(:,:), intent(in) :: models
    real(kind = RPRE), dimension(:), allocatable :: std_scale
    logical, dimension(:), allocatable :: mask

    mask = self % std_scale .ne. 0.
    std_scale = merge(self % std_scale, real(1., RPRE), mask)
    out = ( models - spread(self % mu_scale, 1, self % popsize) ) &
          / spread(std_scale, 1, self % popsize)
    return
  end function standardize2_evolutionary

  function unstandardize1_evolutionary(self, models) result(out)
    real(kind = RPRE), dimension(:), allocatable :: out
    class(Evolutionary), intent(inout) :: self
    real(kind = RPRE), dimension(:), intent(in) :: models
    out = models * self % std_scale + self % mu_scale
    return
  end function unstandardize1_evolutionary

  function unstandardize2_evolutionary(self, models) result(out)
    real(kind = RPRE), dimension(:,:), allocatable :: out
    class(Evolutionary), intent(inout) :: self
    real(kind = RPRE), dimension(:,:), intent(in) :: models
    out = models * spread(self % std_scale, 1, self % popsize) &
          + spread(self % mu_scale, 1, self % popsize)
    return
  end function unstandardize2_evolutionary

  subroutine init_models_evolutionary(self)
    class(Evolutionary), intent(inout) :: self
    if ( .not. allocated(self % models) ) then
      allocate(self % models(self % popsize, self % n_dim, self % max_iter))
      allocate(self % energy(self % popsize, self % max_iter))
      allocate(self % fitness(self % max_iter))
    else
      self % models = 0.
      self % energy = 0.
      self % fitness = 0.
    end if
    return
  end subroutine init_models_evolutionary

  function eval_models(self, models) result(fit)
    real(kind = RPRE), dimension(:), allocatable :: fit
    class(Evolutionary), intent(inout) :: self
    real(kind = RPRE), dimension(:), allocatable :: fit_mpi
    real(kind = RPRE), dimension(:,:), intent(in) :: models
    integer(kind = IPRE) :: i, n

    n = size(models, 1)
!dec$ if defined( do_mpi )
    allocate(fit(n), fit_mpi(n))
    fit = 0.
    fit_mpi = 0.
    call mpi_barrier(mpi_comm_world, self % mpi_ierr)
    call mpi_bcast(models, n * self % n_dim, mpi_rpre(), 0, mpi_comm_world, self % mpi_ierr)
    do i = self % mpi_rank+1, n, self % mpi_size
      fit_mpi(i) = self % func(self % unstandardize(models(i,:)))
    end do
    call mpi_barrier(mpi_comm_world, self % mpi_ierr)
    call mpi_allreduce(fit_mpi, fit, n, mpi_rpre(), mpi_sum, mpi_comm_world, self % mpi_ierr)
!dec$ else
    fit = [ ( self % func(self % unstandardize(models(i,:))), i = 1, n ) ]
!dec$ endif
    return
  end function eval_models

  function constrain_de(self, models) result(out)
    real(kind = RPRE), dimension(:,:), allocatable :: out
    class(Evolutionary), intent(inout) :: self
    real(kind = RPRE), dimension(:,:), intent(in) :: models

    out = merge(randu(self % popsize, self % n_dim), models, models .lt. -1. .or. models .gt. 1.)
    return
  end function constrain_de

  function constrain_cpso(self, models, models_old) result(out)
    real(kind = RPRE), dimension(:,:), allocatable :: out
    class(Evolutionary), intent(inout) :: self
    real(kind = RPRE), dimension(:,:), intent(in) :: models, models_old
    integer(kind = IPRE) :: i
    real(kind = RPRE) :: beta_l, beta_u, beta
    logical, dimension(:), allocatable :: maskl, masku

    out = models
    do i = 1, self % popsize
      maskl = models(i,:) .lt. -1.
      masku = models(i,:) .gt. 1.
      if ( any(maskl) .and. any(masku) ) then
        beta_l = minval( ( models_old(i,:) + 1. ) / ( models_old(i,:) - models(i,:) ), mask = maskl )
        beta_u = minval( ( models_old(i,:) - 1. ) / ( models_old(i,:) - models(i,:) ), mask = masku )
        beta = min(beta_l, beta_u)
        out(i,:) = models_old(i,:) + beta * ( models(i,:) - models_old(i,:) )
      else if ( any(maskl) .and. .not. any(masku) ) then
        beta = minval( ( models_old(i,:) + 1. ) / ( models_old(i,:) - models(i,:) ), mask = maskl )
        out(i,:) = models_old(i,:) + beta * ( models(i,:) - models_old(i,:) )
      else if ( .not. any(maskl) .and. any(masku) ) then
        beta = minval( ( models_old(i,:) - 1. ) / ( models_old(i,:) - models(i,:) ), mask = masku )
        out(i,:) = models_old(i,:) + beta * ( models(i,:) - models_old(i,:) )
      end if
    end do
    return
  end function constrain_cpso

  function de_mutation(self, X, F, gbest, strategy) result(V)
    real(kind = RPRE), dimension(:,:), allocatable :: V
    class(Evolutionary), intent(inout) :: self
    real(kind = RPRE), dimension(:,:), intent(in) :: X
    real(kind = RPRE), intent(in) :: F
    real(kind = RPRE), dimension(:), intent(in) :: gbest
    character(len = *), intent(in) :: strategy
    integer(kind = IPRE), dimension(:,:), allocatable :: idx
    real(kind = RPRE), dimension(:,:), allocatable :: X1, X2, X3, X4, X5

    select case(trim(strategy))
    case("rand1")
      idx = make_idx(self % popsize, 3)
      X1 = X(idx(:,1),:)
      X2 = X(idx(:,2),:)
      X3 = X(idx(:,3),:)
      V = X1 + F * ( X2 - X3 )
    case("rand2")
      idx = make_idx(self % popsize, 5)
      X1 = X(idx(:,1),:)
      X2 = X(idx(:,2),:)
      X3 = X(idx(:,3),:)
      X4 = X(idx(:,4),:)
      X5 = X(idx(:,5),:)
      V = X1 + F * ( X2 + X3 - X4 - X5 )
    case("best1")
      idx = make_idx(self % popsize, 2)
      X1 = X(idx(:,1),:)
      X2 = X(idx(:,2),:)
      V = spread(gbest, 1, self % popsize) + F * ( X1 - X2 )
    case("best2")
      idx = make_idx(self % popsize, 4)
      X1 = X(idx(:,1),:)
      X2 = X(idx(:,2),:)
      X3 = X(idx(:,3),:)
      X4 = X(idx(:,4),:)
      V = spread(gbest, 1, self % popsize) + F * ( X1 + X2 - X3 - X4 )
    end select
    return
  contains
    function make_idx(popsize, n) result(idx)
      integer(kind = IPRE), dimension(:,:), allocatable :: idx
      integer(kind = IPRE), intent(in) :: popsize, n
      integer(kind = IPRE) :: i
      integer(kind = IPRE), dimension(:), allocatable :: tmp

      allocate(idx(popsize, n))
      do i = 1, popsize
        tmp = randperm(popsize, n+1)
        tmp = pack(tmp, tmp .ne. i)
        idx(i,:) = tmp(:n)
      end do
      return
    end function make_idx
  end function de_mutation

  subroutine de(self, CR, F, strategy, xstart)
    class(Evolutionary), intent(inout) :: self
    real(kind = RPRE), intent(in) :: CR, F
    character(len = *), intent(in) :: strategy
    real(kind = RPRE), dimension(:,:), allocatable, intent(in) :: xstart
    integer(kind = IPRE) :: i, it, gbidx,gbidx0,IT0
    integer(kind = IPRE), dimension(:), allocatable :: irand
    real(kind = RPRE) :: gfit
    real(kind = RPRE), dimension(:), allocatable :: pfit, pbestfit, gbest, xopt,gfitH1
    real(kind = RPRE), dimension(:,:), allocatable :: X, V, U, r1
    logical :: converge
    logical, dimension(:), allocatable :: mask1
    logical, dimension(:,:), allocatable :: mask2

    ! Verbosity
    if ( self % verbose .eq. 1 ) call progress_perc(0, self % max_iter, " Processing: ")

    ! Population initial positions
    if ( allocated(xstart) ) then
      X = self % standardize(xstart)
    else
      X = 2. * randu(self % popsize, self % n_dim) - 1.
    end if

    allocate(gfith1(self.max_iter))
    ! Compute fitness
    pfit = self % eval_models(X)
    pbestfit = pfit
    self % n_eval = self % popsize
    if ( self % verbose .eq. 1 ) call progress_perc(1, self % max_iter, " Processing: ")

    ! Initialize best individual
    gbidx = minloc(pbestfit, dim = 1)
    gfit = pbestfit(gbidx)
    gbest = X(gbidx,:)

    if ( self % snap ) then
      call self % init_models()
      self % models(:,:,1) = self % unstandardize(X)
      self % energy(:,1) = pbestfit
      self % fitness(1) = gfit
    end if
    !LGY
    !self % bestmodel(:,:) = self % unstandardize(X) !bestmodel[popsize,ndim]

    ! Iterate until one of the termination criterion is satisfied
    it = 1
    gfith1(it)=gfit
    converge = .false.
    allocate(mask2(self % popsize, self % n_dim))
    do while ( .not. converge )
      it = it + 1
      r1 = randu(self % popsize, self % n_dim)

      ! Mutation
      V = self % de_mutation(X, F, gbest, strategy)

      ! Recombination
      mask2 = .false.
      irand = randi(self % n_dim, self % popsize)
      do i = 1, self % popsize
        mask2(i,irand(i)) = .true.
      end do
      mask2 = mask2 .or. r1 .le. CR
      U = merge(V, X, mask2)
      if ( self % constrain ) U = self % constrain_de(U)

      ! Selection
      pfit = self % eval_models(U)
      self % n_eval = self % n_eval + self % popsize
      mask1 = pfit .lt. pbestfit
      pbestfit = merge(pfit, pbestfit, mask1)
      X = merge(U, X, spread(mask1, 2, self % n_dim))
      !LGY
      !self % bestmodel(:,:) = self % unstandardize(X) !bestmodel[popsize,ndim]
      
      
      ! Update best individual
      
      gbidx0 = MAXloc(pbestfit, dim = 1)

      gbidx = minloc(pbestfit, dim = 1)

      ! Stop if best individual position changes less than eps1
      if ( norm(gbest - X(gbidx,:)) .lt. self % eps1 &
        .and.  abs(pbestfit(gbidx)-pbestfit(gbidx0)) .lt. self % eps2 ) then
        converge = .true.
        xopt = self % unstandardize(X(gbidx,:))
        gfit = pbestfit(gbidx)
        self % iflag = 0

      ! Stop if maxfitness and minfitness is less than eps2
      else if (  abs(pbestfit(gbidx)-pbestfit(gbidx0)) .lt. self % eps2  ) then
        converge = .true.
        xopt = self % unstandardize(X(gbidx,:))
        gfit = pbestfit(gbidx)
        self % iflag = 1
     !stop if fitness is improved less than 0.001 for each 100 iteration
      elseif(mod(it,100)==1) then
        if(abs((pbestfit(gbidx)-gfith1(it-100))/pbestfit(gbidx)) .lt. self.eps3 ) then
            converge = .true.
            xopt = self % unstandardize(X(gbidx,:))
            gfit = pbestfit(gbidx)
            self % iflag = 101
        endif
      ! Stop if maximum iteration is reached
      else if ( it .ge. self % max_iter ) then
        converge = .true.
        xopt = self % unstandardize(X(gbidx,:))
        gfit = pbestfit(gbidx)
        self % iflag = -1

      ! Otherwise, update best individual
      else
        gbest = X(gbidx,:)
        gfit = pbestfit(gbidx)
      end if

      gfith1(it)=gfit
      self % n_iter = it
      EA_ITER=IT
      ! Save models and energy
      if ( self % snap ) then
        self % models(:,:,it) = self % unstandardize(X)
        self % energy(:,it) = pfit
        self % fitness(it) = gfit
      end if

      if ( self % verbose .eq. 1 ) call progress_perc(it, self % max_iter, " Processing: ")
    end do

    self % xopt = xopt
    self % gfit = gfit
    !self % n_iter = it
    if ( self % snap ) then
      self % models = self % models(:,:,:it)
      self % energy = self % energy(:,:it)
      self % fitness = self % fitness(:it)
    end if
    if ( self % verbose .eq. 1 ) then
      call progress_perc(self % max_iter, self % max_iter, " Processing: ")
      print *
    end if
    return
  end subroutine de

  subroutine cpso(self, w, c1, c2, gamma, xstart)
    class(Evolutionary), intent(inout) :: self
    real(kind = RPRE), intent(in) :: w, c1, c2, gamma
    real(kind = RPRE), dimension(:,:), allocatable, intent(in) :: xstart
    integer(kind = IPRE) :: i, it, gbidx, nw,gbidx0 
    integer(kind = IPRE), dimension(:), allocatable :: idx
    real(kind = RPRE) :: gfit, delta, inorm, swarm_radius
    real(kind = RPRE), dimension(:), allocatable :: pfit, pbestfit, gbest, xopt,GFITH1
    real(kind = RPRE), dimension(:,:), allocatable :: X, pbest, V, r1, r2
    logical :: converge
    logical, dimension(:), allocatable :: mask

    ! Verbosity
    if ( self % verbose .eq. 1 ) call progress_perc(0, self % max_iter, " Processing: ")

    ! Particles initial positions
    if ( allocated(xstart) ) then
      X = self % standardize(xstart)
    else
      X = 2. * randu(self % popsize, self % n_dim) - 1.
    end if
    pbest = X

    ! Initialize particle velocity
    allocate(V(self % popsize, self % n_dim))
    V = 0.

    ! Compute fitness
    pfit = self % eval_models(X)
    pbestfit = pfit
    self % n_eval = self % popsize
    if ( self % verbose .eq. 1 ) call progress_perc(1, self % max_iter, " Processing: ")

    ! Initialize best individual
    gbidx = minloc(pbestfit, dim = 1)
    gbidx0=gbidx
    gfit = pbestfit(gbidx)
    gbest = X(gbidx,:)

    if ( self % snap ) then
      call self % init_models()
      self % models(:,:,1) = self % unstandardize(X)
      self % energy(:,1) = pbestfit
      self % fitness(1) = gfit
    end if
    !LGY
    !self % bestmodel(:,:) = self % unstandardize(X) !bestmodel[popsize,ndim]

    ! Swarm maximum radius
    delta = log(1. + 0.003 * self % popsize) / max(0.2, log(0.01*real(self % max_iter)))

    ! Iterate until one of the termination criterion is satisfied
    allocate(gfith1(self.max_iter))
    it = 1
    converge = .false.
    gfith1(it)=gfit
    do while ( .not. converge )
      it = it + 1
      r1 = randu(self % popsize, self % n_dim)
      r2 = randu(self % popsize, self % n_dim)

      ! Mutation
      V = w * V + c1 * r1 * (pbest - X) + c2 * r2 * (spread(gbest, 1, self % popsize) - X)
      if ( self % constrain ) then
        X = self % constrain_cpso(X + V, X)
      else
        X = X + V
      end if

      ! Selection
      pfit = self % eval_models(X)
      self % n_eval = self % n_eval + self % popsize
      mask = pfit .lt. pbestfit
      pbestfit = merge(pfit, pbestfit, mask)
      pbest = merge(X, pbest, spread(mask, 2, self % n_dim))

      ! Update best individual
      gbidx0= maxloc(pbestfit, dim = 1)
      gbidx = minloc(pbestfit, dim = 1)

      
      ! Stop if best individual position changes less than eps1
      if ( norm(gbest - pbest(gbidx,:)) .lt. self % eps1 &
        .and. abs(pbestfit(gbidx)-pbestfit(gbidx0)) .lt. self % eps2   ) then
        converge = .true.
        xopt = self % unstandardize(pbest(gbidx,:))
        gfit = pbestfit(gbidx)
        self % iflag = 0

      ! Stop if fitness is less than eps2
      else if ( abs(pbestfit(gbidx)-pbestfit(gbidx0)) .lt. self % eps2 ) then
        converge = .true.
        xopt = self % unstandardize(pbest(gbidx,:))
        gfit = pbestfit(gbidx)
        self % iflag = 1
      !stop if fitness is improved less than 0.001 for each 200 iteration
      elseif(mod(it,200)==1) then
        if(abs((pbestfit(gbidx)-gfith1(it-200))/pbestfit(gbidx)) .lt. self.eps3 ) then
            converge = .true.
            xopt = self % unstandardize(X(gbidx,:))
            gfit = pbestfit(gbidx)
            self % iflag = 201
        endif  

      ! Stop if maximum iteration is reached
      else if ( it .ge. self % max_iter ) then
        converge = .true.
        xopt = self % unstandardize(pbest(gbidx,:))
        gfit = pbestfit(gbidx)
        self % iflag = -1

      ! Otherwise, update best individual
      else
        gbest = pbest(gbidx,:)
        gfit = pbestfit(gbidx)
      end if
      
      gfith1(it)=gfit
      self % n_iter = it
      EA_ITER=IT
      
      ! Save models and energy
      if ( self % snap ) then
        self % models(:,:,it) = self % unstandardize(X)
        self % energy(:,it) = pfit
        self % fitness(it) = gfit
      end if
      !LGY
      !self % bestmodel(:,:) = self % unstandardize(X) !bestmodel[popsize,ndim]
      
      ! Competitive PSO algorithm
      if ( .not. converge .and. gamma .gt. 0. ) then
        swarm_radius = maxval([ ( norm(X(i,:) - gbest), i = 1, self % popsize ) ])
        swarm_radius = swarm_radius / sqrt(4.*self % n_dim)

        ! Reset positions if swarm size if lower than threshold
        if ( swarm_radius .lt. delta ) then
          inorm = real(it, RPRE) / real(self % max_iter, RPRE)
          nw = int((self % popsize - 1.) / (1. + exp(1./0.09*(inorm-gamma+0.5))))

          ! Reset positions, velocities and personal bests
          if ( nw .gt. 0 ) then
            self % n_restart = self % n_restart + 1
            idx = argsort(pbestfit, 2)
            idx = idx(:nw)
            V(idx,:) = 0.
            do i = 1, nw
              X(idx(i),:) = merge(2.*randu(self % n_dim)-1., real(0., RPRE), self % std_scale .ne. 0.)
            end do
            pbest(idx,:) = X(idx,:)
            pbestfit(idx) = 1e30
          end if
        end if
      end if

      if ( self % verbose .eq. 1 ) call progress_perc(it, self % max_iter, " Processing: ")
    end do

    self % xopt = xopt
    self % gfit = gfit
    self % n_iter = it
    if ( self % snap ) then
      self % models = self % models(:,:,:it)
      self % energy = self % energy(:,:it)
      self % fitness = self % fitness(:it)
    end if
    if ( self % verbose .eq. 1 ) then
      call progress_perc(self % max_iter, self % max_iter, " Processing: ")
      print *
    end if
    return
  end subroutine cpso

  subroutine cmaes(self, sigma, mu_perc, xstart)
    class(Evolutionary), intent(inout) :: self
    real(kind = RPRE), intent(in) :: sigma, mu_perc
    real(kind = RPRE), dimension(:,:), allocatable, intent(in) :: xstart
    integer(kind = IPRE) :: i, k, it, mu, eigeneval, ilim
    real(kind = RPRE) :: mueff, cc, cs, c1, cmu, damps, chind, hsig, perc(2), &
      delta, insigma
    integer(kind = IPRE), dimension(:), allocatable :: arindex
    real(kind = RPRE), dimension(:), allocatable :: xmean, xold, weights, &
      pc, ps, D, arfitness, arbestfitness, bnd_weights, bnd_scale, dfithist, tx
    real(kind = RPRE), dimension(:,:), allocatable :: B, C, invsqrtC, arx, arxvalid, &
      artmp, pctmp
    logical :: validfitval, iniphase, converge
    logical, dimension(:), allocatable :: ti, idx

    ! Initialize saved outputs
    if ( self % snap ) then
      call self % init_models()
      allocate(self % means(self % max_iter, self % n_dim))
    end if

    ! Population initial positions
    if ( allocated(xstart) ) then
      xmean = self % standardize([ xstart ])
    else
      xmean = 2. * randu(self % n_dim) - 1.
    end if
    allocate(xold(self % n_dim))

    ! Number of parents
    mu = int(mu_perc * self % popsize)

    ! Strategy parameter setting: Selection
    weights = log(mu + 0.5) - log([ ( real(i, RPRE), i = 1, mu ) ])
    weights = weights / sum(weights)
    mueff = sum(weights)**2 / sum(weights**2)

    ! Strategy parameter setting: Adaptation
    cc = ( 4. + mueff / self % n_dim ) / ( self % n_dim + 4. + 2. * mueff / self % n_dim )
    cs = ( mueff + 2. ) / ( self % n_dim + mueff + 5. )
    c1 = 2. / ( ( self % n_dim + 1.3 )**2 + mueff )
    cmu = min(1. - c1, 2. * ( mueff - 2. + 1. / mueff ) / ( ( self % n_dim + 2. )**2 + mueff ) )
    damps = 1. + 2. * max(0., sqrt( ( mueff - 1. ) / ( self % n_dim + 1. ) ) - 1.) + cs

    ! Initialize dynamic (internal) strategy parameters and constants
    arx = zeros(self % popsize, self % n_dim)
    pc = zeros(self % n_dim)
    ps = zeros(self % n_dim)
    B = eye(self % n_dim, self % n_dim)
    D = ones(self % n_dim)
    C = eye(self % n_dim, self % n_dim)
    invsqrtC = eye(self % n_dim, self % n_dim)
    chind = sqrt(real(self % n_dim, RPRE)) * ( 1. - 1. / ( 4. * self % n_dim ) + 1. / ( 21. * self % n_dim**2 ) )

    ! Initialize boundaries weights
    bnd_weights = zeros(self % n_dim)
    bnd_scale = zeros(self % n_dim)
    dfithist = [ real(1., RPRE) ]

    ! (opt_mu, opt_lambda)-CMA-ES
    it = 0
    eigeneval = 0
    arbestfitness = zeros(self % max_iter)
    ilim = int(10 + 30 * self % n_dim / self % popsize)
    insigma = sigma
    validfitval = .false.
    iniphase = .true.
    converge = .false.

    do while ( .not. converge )
      it = it + 1

      ! Generate lambda offspring
      do k = 1, self % popsize
        arx(k,:) = xmean + insigma * matmul(B, D*randn(self % n_dim))
      end do
      arxvalid = arx

      ! Evaluate Fitness
      if ( self % constrain ) then
        ! Clip to boundaries
        arxvalid = merge(real(-1., RPRE), arxvalid, arxvalid .lt. -1.)
        arxvalid = merge(real(1., RPRE), arxvalid, arxvalid .gt. 1.)
        arfitness = self % eval_models(arxvalid)

        ! Get delta fitness values
        perc = prctile(arfitness, [ 25, 75 ])
        delta = ( perc(2) - perc(1) ) / real(self % n_dim, RPRE) / mean(diag(C)) / insigma**2

        ! Catch non-sensible values
        if ( delta .eq. 0. ) then
          delta = minval(dfithist, mask = dfithist .gt. 0.)
        else if ( .not. validfitval ) then
          dfithist = zeros(0)
          validfitval = .true.
        end if

        ! Store delta fitness values
        if ( size(dfithist) .lt. 20. + (3.*self % n_dim) / self % popsize ) then
          dfithist = [ dfithist, delta ]
        else
          dfithist = [ dfithist(2:), delta ]
        end if

        ! Corrected mean
        ti = xmean .lt. -1. .or. xmean .gt. 1.
        tx = merge(xmean, real(-1., RPRE), xmean .ge. -1.)
        tx = merge(xmean, real(1., RPRE), xmean .le. 1.)

        ! Set initial weights
        if ( iniphase ) then
          if ( any(ti) ) then
            bnd_weights = 2.0002 * median(dfithist)
            if ( validfitval .and. it .gt. 2 ) iniphase = .false.
          end if
        end if

        if ( any(ti) ) then
          tx = xmean - tx
          idx = ti .and. abs(tx) .gt. 3. * max(1., sqrt(real(self % n_dim, RPRE))/mueff) &
                * insigma * sqrt(diag(C))
          idx = idx .and. signum(tx) .eq. signum(xmean - xold)
          where ( idx )
            bnd_weights = bnd_weights * 1.2**(min(1., mueff/10./real(self % n_dim, RPRE)))
          end where
        end if

        ! Calculate scaling biased to unity, product is one
        bnd_scale = exp( 0.9 * ( log(diag(C)) - mean(log(diag(C))) ) )

        ! Assigned penalized fitness
        arfitness = self % eval_models(arxvalid) &
                    + matmul((arxvalid - arx)**2, bnd_weights / bnd_scale)
      else
        arfitness = self % eval_models(arxvalid)
      end if
      self % n_eval = self % n_eval + self % popsize
      if ( self % snap ) then
        self % models(:,:,it) = self % unstandardize(arx)
        self % energy(:,it) = arfitness
        self % fitness(it) = minval(arfitness)
        self % means(it,:) = self % unstandardize(xmean)
      end if

      ! Sort by fitness and compute weighted mean into xmean
      arindex = argsort(arfitness)
      xold = xmean
      xmean = matmul(weights, arx(arindex(:mu),:))

      ! Save best Fitness
      arbestfitness(it) = arfitness(arindex(1))

      ! Cumulation
      ps = ( 1. - cs ) * ps &
           + sqrt( cs * ( 2. - cs ) * mueff ) * matmul(invsqrtC, xmean - xold) / insigma
      if ( norm(ps) / sqrt( 1. - ( 1. - cs )**(2.*real(self % n_eval, RPRE) / real(self % popsize, RPRE))) / chind &
        .lt. 1.4 + 2. / ( real(self % n_dim) + 1. ) ) then
        hsig = 1.
        pc = ( 1. - cc ) * pc &
             + sqrt( cc * ( 2. - cc ) * mueff ) * (xmean - xold) / insigma
      else
        hsig = 0.
        pc = ( 1. - cc ) * pc
      end if

      ! Adapt covariance matrix C
      artmp = ( arx(arindex(:mu),:) - spread(xold, 1, mu) ) / insigma
      pctmp = reshape(pc, shape = [ size(pc), 1 ])
      if ( hsig .eq. 1. ) then
        C = ( 1. - c1 - cmu ) * C &
            + c1 * matmul(pctmp, transpose(pctmp)) &
            + cmu * matmul(matmul(transpose(artmp), diag(weights)), artmp)
      else
        C = ( 1. - c1 - cmu ) * C &
            + c1 * ( matmul(pctmp, transpose(pctmp)) + cc * ( 2. - cc ) * C ) &
            + cmu * matmul(matmul(transpose(artmp), diag(weights)), artmp)
      end if

      ! Adapt step size sigma
      insigma = insigma * exp( ( cs / damps ) * ( norm(ps) / chind - 1. ) )

      ! Diagonalization of C
      if ( self % n_eval - eigeneval .gt. real(self % popsize, RPRE) / ( c1 + cmu ) / self % n_dim / 10. ) then
        eigeneval =  self % n_eval
        C = triu(C) + transpose(triu(C, 1))
        call eig(C, B, D)
        D = sqrt(D)
        invsqrtC = matmul(matmul(B, diag(real(1./D, RPRE))), transpose(B))
      end if

      ! Stop if maximum number of iteration is reached
      if ( it .ge. self % max_iter ) then
        converge = .true.
        self % iflag = -1
      end if

      ! Stop if mean position changes less than eps1
      if ( .not. converge .and. norm(xold - xmean) .le. self % eps1 ) then
        converge = .true.
        self % iflag = 0
      end if

      ! Stop if fitness is less than eps2
      if ( .not. converge .and. arfitness(arindex(1)) .le. self % eps2 ) then
        converge = .true.
        self % iflag = 1
      end if

      ! NoEffectAxis: stop if numerical precision problem
      i = floor( mod( real(it), real(self % n_dim) ) ) + 1
      if ( .not. converge .and. all( 0.1 * insigma * B(:,i) * D(i) .lt. 1e-10 ) ) then
        converge = .true.
        self % iflag = 2
      end if

      ! NoEffectCoord: stop if too low coordinate axis deviations
      if ( .not. converge .and. any( 0.2 * insigma * sqrt(diag(C)) .lt. 1e-10 ) ) then
        converge = .true.
        self % iflag = 3
      end if

      ! ConditionCov: stop if the condition number exceeds 1e14
      if ( .not. converge .and. maxval(D) .gt. 1e7 * minval(D) ) then
        converge = .true.
        self % iflag = 4
      end if

      ! EqualFunValues: stop if the range of fitness values is zero
      if ( .not. converge .and. it .ge. ilim ) then
        if ( maxval(arbestfitness(it-ilim+1:it)) &
             - minval(arbestfitness(it-ilim+1:it)) .lt. 1e-10 ) then
          converge = .true.
          self % iflag = 5
        end if
      end if

      ! TolXUp: stop if x-changes larger than 1e3 times initial sigma
      if ( .not. converge .and. any( insigma * sqrt(diag(C)) .gt. 1e3 * sigma ) ) then
        converge = .true.
        self % iflag = 6
      end if

      ! TolFun: stop if fun-changes smaller than 1e-12
      if ( .not. converge .and.  it .gt. 2 &
           .and. maxval( [ arfitness, arbestfitness ] ) &
                 - minval( [ arfitness, arbestfitness ] ) .lt. 1e-12 ) then
        converge = .true.
        self % iflag = 7
      end if

      ! TolX: stop if x-changes smaller than 1e-11 times initial sigma
      if ( .not. converge .and. all( insigma * ( max( abs(pc), sqrt(diag(C)) ) ) .lt. 1e-11 * sigma ) ) then
        converge = .true.
        self % iflag = 8
      end if

      if ( self % verbose .eq. 1 ) call progress_perc(it, self % max_iter, " Processing: ")
    end do

    self % xopt = self % unstandardize(arx(arindex(1),:))
    self % gfit = arfitness(arindex(1))
    self % n_iter = it
    if ( self % snap ) then
      self % models = self % models(:,:,:it)
      self % energy = self % energy(:,:it)
      self % fitness = self % fitness(:it)
      self % means = self % means(:it,:)
    end if
    if ( self % verbose .eq. 1 ) then
      call progress_perc(self % max_iter, self % max_iter, " Processing: ")
      print *
    end if
    return
  end subroutine cmaes

  subroutine save_models_evolutionary(self, prefix)
    class(Evolutionary), intent(inout) :: self
    character(len = *), intent(in), optional :: prefix
    character(len = :), allocatable :: opt_prefix, str_popsize, str_n_dim, &
      str_n_iter, models_filename, energy_filename, fitness_filename, means_filename

    if ( self % mpi_rank .eq. 0 ) then
      if ( present(prefix) ) then
        opt_prefix = trim(prefix)
      else
        opt_prefix = "./"
      end if

      if ( self % snap ) then
        str_popsize = num2str(self % popsize)
        str_n_dim = num2str(self % n_dim)
        str_n_iter = num2str(self % n_iter)
        models_filename = "models_" // str_popsize // "_" // str_n_dim &
                          // "_" // str_n_iter // ".bin"
        energy_filename = "energy_" // str_popsize // "_" // str_n_iter // ".bin"
        fitness_filename = "fitness_" // str_n_iter // ".bin"
        call savebin(opt_prefix // models_filename, self % models)
        call savebin(opt_prefix // energy_filename, self % energy)
        call savebin(opt_prefix // fitness_filename, self % fitness)
        if ( self % solver .eq. "cmaes" ) then
          means_filename = "means_" // str_n_iter // "_" // str_n_dim // ".bin"
          call savebin(opt_prefix // means_filename, self % means)
        end if
      else
        print *, "Warning: cannot save models when snap is false"
      end if
    end if
    return
  end subroutine save_models_evolutionary

  subroutine sample(self, sampler, stepsize, xstart)
    class(MonteCarlo), intent(inout) :: self
    character(len = *), intent(in), optional :: sampler
    real(kind = RPRE), dimension(:), intent(in), optional :: stepsize, xstart
    real(kind = RPRE), dimension(:), allocatable :: opt_stepsize, opt_xstart
    character(len = :), allocatable :: opt_sampler
    integer(kind = 4) :: master = 0, num_procs, rank, ierr

    ! Check inputs
    if ( present(sampler) ) then
      if ( trim(sampler) .ne. "mcmc" ) then
        print *, "Error: solver must be 'mcmc'"
        stop
      else
        opt_sampler = trim(sampler)
      end if
    else
      opt_sampler = "mcmc"
    end if
    if ( present(stepsize) ) then
      if ( size(stepsize) .ne. self % n_dim ) then
        print *, "Error: stepsize must be a 1-D array of length " // num2str(self % n_dim)
        stop
      end if
      if ( any(stepsize .le. 0.) ) then
        print *, "Error: elements in stepsize must be positive"
        stop
      else
        opt_stepsize = stepsize
      end if
    else
      allocate(opt_stepsize(self % n_dim))
      opt_stepsize = 0.1
    end if
    if ( present(xstart) ) then
      if ( size(xstart) .ne. self % n_dim ) then
        print *, "Error: xstart must be a 1-D array of length " // num2str(self % n_dim)
        stop
      else
        opt_xstart = self % standardize(xstart)
      end if
    else
      opt_xstart = 2. * randu(self % n_dim) - 1.
    end if

    ! Initialize
    self % sampler = opt_sampler
    self % mu_scale = 0.5 * (self % upper + self % lower)
    self % std_scale = 0.5 * (self % upper - self % lower)

    ! Sample
    select case(opt_sampler)
    case("mcmc")
      call self % mcmc(stepsize = opt_stepsize, xstart = opt_xstart)
    end select
    return
  end subroutine sample

  function standardize_montecarlo(self, models) result(out)
    real(kind = RPRE), dimension(:), allocatable :: out
    class(MonteCarlo), intent(inout) :: self
    real(kind = RPRE), dimension(:), intent(in) :: models
    real(kind = RPRE), dimension(:), allocatable :: std_scale
    logical, dimension(:), allocatable :: mask

    out = ( models - self % mu_scale ) / std_scale
    return
  end function standardize_montecarlo

  function unstandardize1_montecarlo(self, models) result(out)
    real(kind = RPRE), dimension(:), allocatable :: out
    class(MonteCarlo), intent(inout) :: self
    real(kind = RPRE), dimension(:), intent(in) :: models
    out = models * self % std_scale + self % mu_scale
    return
  end function unstandardize1_montecarlo

  function unstandardize2_montecarlo(self, models) result(out)
    real(kind = RPRE), dimension(:,:), allocatable :: out
    class(MonteCarlo), intent(inout) :: self
    real(kind = RPRE), dimension(:,:), intent(in) :: models
    out = models * spread(self % std_scale, 1, self % max_iter) &
          + spread(self % mu_scale, 1, self % max_iter)
    return
  end function unstandardize2_montecarlo

  subroutine init_models_montecarlo(self)
    class(MonteCarlo), intent(inout) :: self
    if ( allocated(self % models) ) deallocate(self % models, self % energy)
    allocate(self % models(self % max_iter, self % n_dim))
    allocate(self % energy(self % max_iter))
    return
  end subroutine init_models_montecarlo

  function in_search_space(self, x) result(out)
    logical :: out
    class(MonteCarlo), intent(inout) :: self
    real(kind = RPRE), dimension(:), intent(in) :: x

    if ( self % constrain ) then
      out = all( x .le. 1. ) .and. all( x .ge. -1. )
    else
      out = .true.
    end if
    return
  end function in_search_space

  subroutine mcmc(self, stepsize, xstart)
    class(MonteCarlo), intent(inout) :: self
    real(kind = RPRE), dimension(:), intent(in) :: stepsize, xstart
    integer(kind = IPRE) :: i, j, jmax, rejected, idx
    real(kind = RPRE) :: log_alpha

    ! Verbosity
    if ( self % verbose .eq. 1 ) call progress_perc(0, self % max_iter, " Processing: ")

    ! Initialize models
    call self % init_models()
    self % models(1,:) = xstart
    self % energy(1) = self % func(self % unstandardize(self % models(1,:)))
    if ( self % verbose .eq. 1 ) call progress_perc(1, self % max_iter, " Processing: ")

    ! Metropolis-Hastings algorithm
    rejected = 0
    i = 1
    do while ( i .lt. self % max_iter )
      do j = 1, self % n_dim, self % n_dim_per_iter
        i = i + 1
        jmax = min(self % n_dim, j + self % n_dim_per_iter - 1)
        self % models(i,:) = self % models(i-1,:)
        self % models(i,j:jmax) = self % models(i,j:jmax) + randn(jmax-j+1) * stepsize(j:jmax)
        self % energy(i) = self % func(self % unstandardize(self % models(i,:)))

        log_alpha = min(real(0., RPRE), self % energy(i-1) - self % energy(i))
        if ( log_alpha .lt. log(randu()) &
          .or. .not. self % in_search_space(self % models(i,:)) ) then
          rejected = rejected + 1
          self % models(i,:) = self % models(i-1,:)
          self % energy(i) = self % energy(i-1)
        end if

        if ( i .eq. self % max_iter ) exit
      end do
      if ( self % verbose .eq. 1 ) call progress_perc(i, self % max_iter, " Processing: ")
    end do

    idx = minloc(self % energy, dim = 1)
    self % models = self % unstandardize(self % models)
    self % xopt = self % models(idx,:)
    self % gfit = self % energy(idx)
    self % acceptance_ratio = 1. - real(rejected) / real(self % max_iter)
    if ( self % verbose .eq. 1 ) print *
    return
  end subroutine mcmc

  subroutine save_models_montecarlo(self, prefix)
    class(MonteCarlo), intent(inout) :: self
    character(len = *), intent(in), optional :: prefix
    character(len = :), allocatable :: opt_prefix, str_max_iter, str_n_dim, &
      models_filename, energy_filename

    if ( present(prefix) ) then
      opt_prefix = trim(prefix)
    else
      opt_prefix = "./"
    end if

    str_max_iter = num2str(self % max_iter)
    str_n_dim = num2str(self % n_dim)
    models_filename = "models_" // str_max_iter // "_" // str_n_dim // ".bin"
    energy_filename = "energy_" // str_max_iter // ".bin"
    call savebin(opt_prefix // models_filename, self % models)
    call savebin(opt_prefix // energy_filename, self % energy)
    return
  end subroutine save_models_montecarlo

end module stochoptim
