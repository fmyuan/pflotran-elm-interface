module Solver_module
 
#include "petsc/finclude/petsc.h"
  use petsc
#include "petsc/finclude/petscts.h"
  use petscts
  use PFLOTRAN_Constants_module
  use CPR_Preconditioner_module
  use Solver_CPR_module

  implicit none

  private

  type, public :: solver_type
    !TODO(geh): remove itype in favor of setting prefix through call to
    !           nonlinear/linear solver read routine
    PetscInt :: itype            ! type: flow or transport
    PetscReal :: linear_atol       ! absolute tolerance
    PetscReal :: linear_rtol       ! relative tolerance
    PetscReal :: linear_dtol       ! divergence tolerance
    PetscInt :: linear_max_iterations     ! maximum number of iterations
    PetscReal :: linear_zero_pivot_tol  ! zero pivot tolerance for LU
    PetscBool :: linear_stop_on_failure ! flag determines whether the code is
                                        ! killed when the solver fails, as
                                        ! opposed ot cutting the step.
    PetscBool :: linear_shift      ! shift diagonal to alleviate zero pivots
    PetscReal :: newton_atol       ! absolute tolerance
    PetscReal :: newton_rtol       ! relative tolerance
    PetscReal :: newton_stol       ! relative tolerance (relative to previous iteration)
    PetscReal :: newton_dtol       ! divergence tolerance
    PetscReal :: newton_inf_res_tol    ! infinity tolerance for residual
    PetscReal :: newton_inf_upd_tol    ! infinity tolerance for update
    PetscReal :: newton_inf_rel_update_tol ! infinity norm on relative update (c(i)-c(i-1))/c(i-1)
    PetscReal :: newton_inf_scaled_res_tol ! infinity norm on scale residual (r(i)/accum(i))
    PetscReal :: newton_inf_res_tol_sec  ! infinity tolerance for secondary continuum residual
    PetscInt :: newton_max_iterations     ! maximum number of iterations
    PetscInt :: newton_min_iterations     ! minimum number of iterations
    PetscInt :: newton_maxf      ! maximum number of function evaluations
    PetscReal :: max_norm          ! maximum norm for divergence
    PetscBool :: use_galerkin_mg  ! If true, precondition linear systems with 
                                   ! Galerkin-type geometric multigrid.
    PetscInt :: galerkin_mg_levels  ! Number of discretization levels for 
                                    ! the Galerkin MG (includes finest level).
    PetscInt :: galerkin_mg_levels_x
    PetscInt :: galerkin_mg_levels_y
    PetscInt :: galerkin_mg_levels_z
    PetscBool :: verbose_logging
    PetscBool :: convergence_2r
    PetscBool :: convergence_2x
    PetscBool :: convergence_2u
    PetscBool :: convergence_ir
    PetscBool :: convergence_iu

    ! Jacobian matrix
    Mat :: J    ! Jacobian
    Mat :: Jpre ! Jacobian to be used in preconditioner
    MatType :: J_mat_type
    MatType :: Jpre_mat_type

    MatFDColoring :: matfdcoloring
      ! Coloring used for computing the Jacobian via finite differences.

    Mat, pointer :: interpolation(:)
      ! Hierarchy of interpolation operators for Galerkin multigrid.

    ! PETSc nonlinear solver context
    SNES :: snes
    KSPType :: ksp_type
    PCType :: pc_type
    KSP ::  ksp
    PC ::  pc
    TS :: ts
    
    PetscBool :: inexact_newton

    PetscBool :: print_convergence
    PetscBool :: print_detailed_convergence
    PetscBool :: print_linear_iterations
    PetscBool :: check_infinity_norm

    ! added for CPR option:
    type(cpr_pc_type), pointer :: cprstash
            
  end type solver_type
  
  public :: SolverCreate, &
            SolverDestroy, &
            SolverReadLinear, &
            SolverReadNewton, &
            SolverReadNewtonSelectCase, &
            SolverCreateKSP, &
            SolverSetKSPOptions, &
            SolverCreateSNES, &
            SolverSetSNESOptions, &
            SolverCreateTS, &
            SolverPrintNewtonInfo, &
            SolverPrintLinearInfo, &
            SolverCheckCommandLine, &
            SolverLinearPrintFailedReason, &
            SolverNewtonPrintFailedReason
  
contains

! ************************************************************************** !

function SolverCreate()
  ! 
  ! Allocates and initializes a new (empty) Solver object
  ! Note that this does not create the PETSc solver contexts associated
  ! with the Solver.  These contexts are created via a subsequent call to
  ! SolverCreateSNES().
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  implicit none
  
  type(solver_type), pointer :: SolverCreate
  
  type(solver_type), pointer :: solver
  
  allocate(solver)
  
  ! initialize to default values
  solver%itype = NULL_CLASS
  solver%linear_atol = PETSC_DEFAULT_REAL
  solver%linear_rtol = PETSC_DEFAULT_REAL
  solver%linear_dtol = PETSC_DEFAULT_REAL
  solver%linear_max_iterations = PETSC_DEFAULT_INTEGER
  solver%linear_zero_pivot_tol = UNINITIALIZED_DOUBLE
  solver%linear_stop_on_failure = PETSC_FALSE
  solver%linear_shift = PETSC_TRUE
  
  solver%newton_atol = PETSC_DEFAULT_REAL
  solver%newton_rtol = PETSC_DEFAULT_REAL
  solver%newton_stol = PETSC_DEFAULT_REAL
  solver%newton_dtol = PETSC_DEFAULT_REAL
  solver%max_norm = 1.d20     ! set to a large value
  solver%newton_inf_res_tol = UNINITIALIZED_DOUBLE
  solver%newton_inf_upd_tol = UNINITIALIZED_DOUBLE
  solver%newton_inf_rel_update_tol = UNINITIALIZED_DOUBLE
  solver%newton_inf_scaled_res_tol = UNINITIALIZED_DOUBLE
  solver%newton_inf_res_tol_sec = 1.d-10
  solver%newton_max_iterations = PETSC_DEFAULT_INTEGER
  solver%newton_min_iterations = 1
  solver%newton_maxf = PETSC_DEFAULT_INTEGER

  solver%use_galerkin_mg = PETSC_FALSE
  solver%galerkin_mg_levels = 1
  solver%galerkin_mg_levels_x = 1
  solver%galerkin_mg_levels_y = 1
  solver%galerkin_mg_levels_z = 1

  solver%verbose_logging = PETSC_FALSE
  solver%convergence_2r = PETSC_TRUE
  solver%convergence_2x = PETSC_TRUE
  solver%convergence_2u = PETSC_TRUE
  solver%convergence_ir = PETSC_TRUE
  solver%convergence_iu = PETSC_TRUE
  
  solver%J = PETSC_NULL_MAT
  solver%Jpre = PETSC_NULL_MAT
  solver%J_mat_type = PETSC_NULL_CHARACTER
  solver%Jpre_mat_type = PETSC_NULL_CHARACTER
!  solver%interpolation = 0
  nullify(solver%interpolation)
  solver%matfdcoloring = PETSC_NULL_MATFDCOLORING
  solver%snes = PETSC_NULL_SNES
  solver%ksp_type = KSPBCGS
  solver%pc_type = ""
  solver%ksp = PETSC_NULL_KSP
  solver%pc = PETSC_NULL_PC
  solver%ts = PETSC_NULL_TS
  
  solver%inexact_newton = PETSC_FALSE
  
  solver%print_convergence = PETSC_TRUE
  solver%print_detailed_convergence = PETSC_FALSE
  solver%print_linear_iterations = PETSC_FALSE
  solver%check_infinity_norm = PETSC_TRUE

  nullify(solver%cprstash)
    
  SolverCreate => solver
   
end function SolverCreate

! ************************************************************************** !

subroutine SolverCreateKSP(solver,comm)
  ! 
  ! Create PETSc KSP object
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  ! 

  implicit none
  
  type(solver_type) :: solver

  PetscMPIInt :: comm
  PetscErrorCode :: ierr
  
  call KSPCreate(comm,solver%ksp,ierr);CHKERRQ(ierr)
  call KSPSetFromOptions(solver%ksp,ierr);CHKERRQ(ierr)

  ! grab handle for pc
  call KSPGetPC(solver%ksp,solver%pc,ierr);CHKERRQ(ierr)

end subroutine SolverCreateKSP

! ************************************************************************** !

subroutine SolverSetKSPOptions(solver, option)
  ! 
  ! Sets options for KSP
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  ! 
  use Option_module

  implicit none
  
  type(solver_type) :: solver
  type(option_type) :: option

  PetscErrorCode :: ierr

  call SolverSetupCustomKSP(solver,option)
  call SolverSetupPCGalerkinMG(solver, option)
  ! KSPSetFromOptions must come after custom setup in order to override 
  ! from command line
  call KSPSetFromOptions(solver%ksp,ierr);CHKERRQ(ierr)
  call SolverSetupPCShiftAndPivoting(solver,option)
  call KSPGetTolerances(solver%ksp,solver%linear_rtol,solver%linear_atol, &
                        solver%linear_dtol,solver%linear_max_iterations, &
                        ierr);CHKERRQ(ierr)

end subroutine SolverSetKSPOptions

! ************************************************************************** !

subroutine SolverSetupCustomKSP(solver, option)
  ! 
  ! Sets options for KSP and PC when specified through input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  ! 
  use Option_module

  implicit none
  
  type(solver_type) :: solver
  type(option_type) :: option

  PetscErrorCode :: ierr

  ! if ksp_type or pc_type specified in input file, set them here
  if (len_trim(solver%ksp_type) > 1) then
    call KSPSetType(solver%ksp,solver%ksp_type,ierr);CHKERRQ(ierr)
  endif
  if (len_trim(solver%pc_type) > 1) then
    if (associated(solver%cprstash)) then
      call SolverCPRInit(solver%J, solver%cprstash, solver%pc, ierr, option)
    else
      call PCSetType(solver%pc,solver%pc_type,ierr);CHKERRQ(ierr)
    endif
  endif

  call KSPSetTolerances(solver%ksp,solver%linear_rtol,solver%linear_atol, &
                        solver%linear_dtol,solver%linear_max_iterations, &
                        ierr);CHKERRQ(ierr)
  ! as of PETSc 3.7, we need to turn on error reporting due to zero pivots
  ! as PETSc no longer reports zero pivots for very small concentrations
  !geh: this gets overwritten by ksp->errorifnotconverted
  if (solver%linear_stop_on_failure) then
    call KSPSetErrorIfNotConverged(solver%ksp,PETSC_TRUE,ierr);CHKERRQ(ierr)
  endif

end subroutine SolverSetupCustomKSP

! ************************************************************************** !

subroutine SolverSetupPCGalerkinMG(solver, option)
  ! 
  ! Sets up a Galerkin multigrid approach
  ! 
  ! Author: Richard Mills
  ! Date: 12/06/19
  ! 
  use Option_module

  implicit none
  
  type(solver_type) :: solver
  type(option_type) :: option

  PetscInt :: i
  PetscErrorCode :: ierr

  ! Setup for n-level Galerkin multigrid.
  if (solver%use_galerkin_mg) then
    call PCSetType(solver%pc, PCMG,ierr);CHKERRQ(ierr)
    call PCMGSetLevels(solver%pc, solver%galerkin_mg_levels, &
                       MPI_COMM_NULL,ierr);CHKERRQ(ierr)
    do i=1,solver%galerkin_mg_levels-1
      call PCMGSetInterpolation(solver%pc, i, solver%interpolation(i), &
                                ierr);CHKERRQ(ierr)
      !geh: not sure if this is the right type....
      call PCMGSetGalerkin(solver%pc,PC_MG_GALERKIN_MAT,ierr);CHKERRQ(ierr)
    enddo
  endif

end subroutine SolverSetupPCGalerkinMG

! ************************************************************************** !

subroutine SolverSetupPCShiftAndPivoting(solver, option)
  ! 
  ! Sets up a shift and pivoting in order to avoid near-zero values on 
  ! the matrix diagonal
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  ! 
  use Option_module

  implicit none
  
  type(solver_type) :: solver
  type(option_type) :: option

  KSP, pointer :: sub_ksps(:)
  PC :: pc
  PetscInt :: nsub_ksp
  PetscInt :: first_sub_ksp
  PetscInt :: i
  PetscErrorCode :: ierr

  if (solver%linear_shift) then
    ! the below must come after SNESSetFromOptions
    ! PETSc no longer performs a shift on matrix diagonals by default.  We 
    ! force the shift since it helps alleviate zero pivots.
    ! Note that if the preconditioner type does not support a shift, the shift 
    ! we've set is ignored; we don't need to check to see if the type supports 
    ! a shift before calling this.
    call PCFactorSetShiftType(solver%pc,MAT_SHIFT_INBLOCKS,ierr);CHKERRQ(ierr)
    if (solver%pc_type == PCBJACOBI .or. solver%pc_type == PCASM .or. &
        solver%pc_type == PCGASM) then
      call KSPSetup(solver%ksp,ierr);CHKERRQ(ierr)
      select case(solver%pc_type)
        case(PCBJACOBI)
          call PCBJacobiGetSubKSP(solver%pc,nsub_ksp,first_sub_ksp, &
                                  PETSC_NULL_KSP,ierr);CHKERRQ(ierr)
        case(PCASM)
          call PCASMGetSubKSP(solver%pc,nsub_ksp,first_sub_ksp, &
                              PETSC_NULL_KSP,ierr);CHKERRQ(ierr)
        case(PCGASM)
          call PCGASMGetSubKSP(solver%pc,nsub_ksp,first_sub_ksp, &
                               PETSC_NULL_KSP,ierr);CHKERRQ(ierr)
      end select
      allocate(sub_ksps(nsub_ksp))
      select case(solver%pc_type)
        case(PCBJACOBI)
          call PCBJacobiGetSubKSP(solver%pc,nsub_ksp,first_sub_ksp, &
                                  sub_ksps,ierr);CHKERRQ(ierr)
        case(PCASM)
          call PCASMGetSubKSP(solver%pc,nsub_ksp,first_sub_ksp, &
                              sub_ksps,ierr);CHKERRQ(ierr)
        case(PCGASM)
          call PCGASMGetSubKSP(solver%pc,nsub_ksp,first_sub_ksp, &
                               sub_ksps,ierr);CHKERRQ(ierr)
      end select
      do i = 1, nsub_ksp
        call KSPGetPC(sub_ksps(i),pc,ierr);CHKERRQ(ierr)
        call PCFactorSetShiftType(pc,MAT_SHIFT_INBLOCKS,ierr);CHKERRQ(ierr)
      enddo
      deallocate(sub_ksps)
      nullify(sub_ksps)
    endif
  endif
  
  if (Initialized(solver%linear_zero_pivot_tol)) then
    call PCFactorSetZeroPivot(solver%pc,solver%linear_zero_pivot_tol, &
                              ierr);CHKERRQ(ierr)
    if (solver%pc_type == PCBJACOBI) then
      call KSPSetup(solver%ksp,ierr);CHKERRQ(ierr)
      call PCBJacobiGetSubKSP(solver%pc,nsub_ksp,first_sub_ksp, &
                              PETSC_NULL_KSP,ierr);CHKERRQ(ierr)
      allocate(sub_ksps(nsub_ksp))
      call PCBJacobiGetSubKSP(solver%pc,nsub_ksp,first_sub_ksp, &
                              sub_ksps,ierr);CHKERRQ(ierr)
      do i = 1, nsub_ksp
        call KSPGetPC(sub_ksps(i),pc,ierr);CHKERRQ(ierr)
        call PCFactorSetZeroPivot(pc,solver%linear_zero_pivot_tol, &
                                  ierr);CHKERRQ(ierr)
      enddo
      deallocate(sub_ksps)
      nullify(sub_ksps)
    elseif (.not.(solver%pc_type == PCLU .or. solver%pc_type == PCILU)) then
      option%io_buffer = 'PCFactorSetZeroPivot for PC ' // &
        trim(solver%pc_type) // ' is not supported at this time.'
      call PrintErrMsg(option)
    endif
  endif

end subroutine SolverSetupPCShiftAndPivoting

! ************************************************************************** !

subroutine SolverCreateSNES(solver,comm)
  ! 
  ! Create PETSc SNES object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/12/08
  ! 

  implicit none
  
  type(solver_type) :: solver

  PetscMPIInt :: comm
  PetscErrorCode :: ierr
  
  call SNESCreate(comm,solver%snes,ierr);CHKERRQ(ierr)
  call SNESSetFromOptions(solver%snes,ierr);CHKERRQ(ierr)

  ! grab handles for ksp and pc
  call SNESGetKSP(solver%snes,solver%ksp,ierr);CHKERRQ(ierr)
  call KSPGetPC(solver%ksp,solver%pc,ierr);CHKERRQ(ierr)

end subroutine SolverCreateSNES

! ************************************************************************** !

subroutine SolverSetSNESOptions(solver, option)
  ! 
  ! Sets options for SNES
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/12/08
  ! 
  use Option_module

  implicit none
  
  type(solver_type) :: solver
  type(option_type) :: option

  SNESLineSearch :: linesearch
  KSP, pointer :: sub_ksps(:)
  PC :: pc
  PetscInt :: nsub_ksp
  PetscInt :: first_sub_ksp
  PetscErrorCode :: ierr
  PetscInt :: i
  
  ! if ksp_type or pc_type specified in input file, set them here
  call SolverSetupCustomKSP(solver,option)
  call SolverSetupPCGalerkinMG(solver,option)

  ! Set the tolerances for the Newton solver.
  call SNESSetTolerances(solver%snes, solver%newton_atol, solver%newton_rtol, &
                         solver%newton_stol,solver%newton_max_iterations, &
                         solver%newton_maxf,ierr);CHKERRQ(ierr)
  call SNESSetDivergenceTolerance(solver%snes,solver%newton_dtol, &
                                  ierr);CHKERRQ(ierr)

  ! set inexact newton, currently applies default settings
  if (solver%inexact_newton) then
    call SNESKSPSetUseEW(solver%snes,PETSC_TRUE,ierr);CHKERRQ(ierr)
  endif

  
  ! allow override from command line; for some reason must come before
  ! LineSearchParams, or they crash
  ! Note that SNESSetFromOptions() calls KSPSetFromOptions(), which calls
  ! PCSetFromOptions(), so these should not be called separately (doing so
  ! causes unintended results when PCCOMPOSITE is used).
  call SNESSetFromOptions(solver%snes,ierr);CHKERRQ(ierr)

  ! get the ksp_type and pc_type incase of command line override.
  call KSPGetType(solver%ksp,solver%ksp_type,ierr);CHKERRQ(ierr)
  call PCGetType(solver%pc,solver%pc_type,ierr);CHKERRQ(ierr)

  call SolverSetupPCShiftAndPivoting(solver,option)

  call SNESGetLineSearch(solver%snes, linesearch, ierr);CHKERRQ(ierr)
  call SNESLineSearchSetTolerances(linesearch, solver%newton_stol,       &
          PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL, &
          PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL, &
          PETSC_DEFAULT_INTEGER, ierr);CHKERRQ(ierr)

  call SNESGetTolerances(solver%snes,solver%newton_atol,solver%newton_rtol, &
                         solver%newton_stol,solver%newton_max_iterations, &
                         solver%newton_maxf,ierr);CHKERRQ(ierr)
  call SNESGetDivergenceTolerance(solver%snes,solver%newton_dtol, &
                                  ierr);CHKERRQ(ierr)
  call KSPGetTolerances(solver%ksp,solver%linear_rtol,solver%linear_atol, &
                        solver%linear_dtol,solver%linear_max_iterations, &
                        ierr);CHKERRQ(ierr)

end subroutine SolverSetSNESOptions

! ************************************************************************** !

subroutine SolverCreateTS(solver,comm)
  ! 
  ! This routine creates PETSc TS object.
  ! 
  ! Author: Gautam Bisht, LBL
  ! Date: 01/18/13
  ! 

  implicit none
  
  type(solver_type) :: solver

  PetscMPIInt :: comm
  PetscErrorCode :: ierr
  
  call TSCreate(comm,solver%ts,ierr);CHKERRQ(ierr)
  call TSSetFromOptions(solver%ts,ierr);CHKERRQ(ierr)

end subroutine SolverCreateTS

! ************************************************************************** !

subroutine SolverReadLinear(solver,input,option)
  ! 
  ! Reads parameters associated with linear solver
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/21/07
  ! 

  use Input_Aux_module
  use String_module
  use Option_module
  
  implicit none

  type(solver_type) :: solver
  type(input_type), pointer :: input
  type(option_type) :: option
  PetscErrorCode :: ierr
  
  character(len=MAXWORDLENGTH) :: keyword, word, word2, prefix
  character(len=MAXSTRINGLENGTH) :: string

  select case(solver%itype)
    case(FLOW_CLASS)
      prefix = '-flow_'
    case(TRANSPORT_CLASS)
      prefix = '-tran_'
  end select

  input%ierr = 0
  call InputPushBlock(input,option)
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword','LINEAR SOLVER')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('SOLVER','KSP_TYPE')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'ksp_type','LINEAR SOLVER')   
        call StringToUpper(word)
        select case(trim(word))
          case('NONE','PREONLY')
            solver%ksp_type = KSPPREONLY
          case('GMRES')
            solver%ksp_type = KSPGMRES
          case('FGMRES')
            solver%ksp_type = KSPFGMRES
          case('BCGS','BICGSTAB','BI-CGSTAB')
            solver%ksp_type = KSPBCGS
          case('IBCGS','IBICGSTAB','IBI-CGSTAB')
            solver%ksp_type = KSPIBCGS
          case('RICHARDSON')
            solver%ksp_type = KSPRICHARDSON
          case('CG')
            solver%ksp_type = KSPCG
          case('DIRECT')
            solver%ksp_type = KSPPREONLY
            solver%pc_type = PCLU
          case('ITERATIVE','KRYLOV')
            solver%ksp_type = KSPBCGS
            solver%pc_type = PCBJACOBI
          case default
            option%io_buffer  = 'Krylov solver type: ' // trim(word) // &
                                ' unknown.'
            call PrintErrMsg(option)
        end select

      case('PRECONDITIONER','PC_TYPE')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'pc_type','LINEAR SOLVER')   
        call StringToUpper(word)
        select case(trim(word))
          case('NONE','PCNONE')
            solver%pc_type = PCNONE
          case('ILU','PCILU')
            solver%pc_type = PCILU
          case('LU','PCLU')
            solver%pc_type = PCLU
          case('BJACOBI','BLOCK_JACOBI')
            solver%pc_type = PCBJACOBI
          case('JACOBI')
            solver%pc_type = PCJACOBI
          case('ASM','ADDITIVE_SCHWARZ')
            solver%pc_type = PCASM
          case('HYPRE')
            solver%pc_type = PCHYPRE
          case('SHELL')
            solver%pc_type = PCSHELL
          case('CPR')
            solver%pc_type = PCSHELL
            allocate(solver%cprstash)
            call SolverCPRInitializeStorage(solver%cprstash)
          case default
            option%io_buffer  = 'Preconditioner type: ' // trim(word) // &
                                ' unknown.'
            call PrintErrMsg(option)
        end select

      case('HYPRE_OPTIONS')
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          if (InputCheckExit(input,option)) exit  
          call InputReadCard(input,option,keyword)
          call InputErrorMsg(input,option,'keyword', &
                             'LINEAR SOLVER, HYPRE options')   
          call StringToUpper(keyword)
          select case(trim(keyword))
            case('TYPE')
              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'type', &
                                 'LINEAR SOLVER, HYPRE options')  
              call StringToUpper(word)
              select case(trim(word))
                case('PILUT','PARASAILS','BOOMERAMG','EUCLID')
                  string = trim(prefix) // 'pc_hypre_type'
                  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                            trim(string),trim(word), &
                                            ierr);CHKERRQ(ierr)
                case default
                  option%io_buffer  = 'HYPRE preconditioner type: ' // &
                                      trim(word) // ' unknown.'
                  call PrintErrMsg(option)
              end select
            case('BOOMERAMG_CYCLE_TYPE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG cycle type', &
                                 'LINEAR SOLVER, HYPRE options')  
              call StringToUpper(word)
              string = trim(prefix) // 'pc_hypre_boomeramg_cycle_type'
              select case(trim(word))
                case('V')
                  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                            trim(string),'1', &
                                            ierr);CHKERRQ(ierr)
                case('W')
                  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                            trim(string),'2', &
                                            ierr);CHKERRQ(ierr)
                case default
                  option%io_buffer  = 'HYPRE BoomerAMG cycle type: ' &
                                      // trim(word) // ' unknown.'
                  call PrintErrMsg(option)
              end select
            case('BOOMERAMG_MAX_LEVELS')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG maximum levels', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_max_levels'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_MAX_ITER')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG maximum iterations', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_max_iter'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_TOL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG convergence tolerance', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_tol'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_TRUNCFACTOR')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG interpolation truncation factor', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_truncfactor'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_AGG_NL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG # levels aggressive coarsening', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_agg_nl'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_AGG_NUM_PATHS')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                'BoomerAMG # paths for aggressive coarsening', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_agg_num_paths'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_STRONG_THRESHOLD')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                'BoomerAMG threshold for strong connectivity', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_strong_threshold'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_GRID_SWEEPS_ALL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                         'BoomerAMG number of grid sweeps up and down cycles', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_grid_sweeps_all'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_GRID_SWEEPS_DOWN')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                'BoomerAMG number of grid sweeps down cycles', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_grid_sweeps_down'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_GRID_SWEEPS_UP')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                  'BoomerAMG number of grid sweeps up cycles', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_grid_sweeps_up'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_GRID_SWEEPS_COARSE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                           'BoomerAMG number of grid sweeps for coarse level', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_grid_sweeps_coarse'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_RELAX_TYPE_ALL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                           'BoomerAMG relaxation type for up and down cycles', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_relax_type_all'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_RELAX_TYPE_DOWN')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                  'BoomerAMG relaxation type for down cycles', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_relax_type_down'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_RELAX_TYPE_UP')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG relaxation type for up cycles', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_relax_type_up'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_RELAX_TYPE_COARSE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG relaxation type for coarse grids', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_relax_type_coarse'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_RELAX_WEIGHT_ALL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG relaxation weight for all levels', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_relax_weight_all'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_RELAX_WEIGHT_LEVEL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputReadWord(input,option,word2,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG relaxation weight for a level', &
                                 'LINEAR SOLVER, HYPRE options')  
              word = trim(word) // ' ' // trim(word2)
              string = trim(prefix) // 'pc_hypre_boomeramg_relax_weight_level'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_OUTER_RELAX_WEIGHT_ALL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                           'BoomerAMG outer relaxation weight for all levels', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) //  &
                       'pc_hypre_boomeramg_outer_relax_weight_all'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_OUTER_RELAX_WEIGHT_LEVEL')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputReadWord(input,option,word2,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                              'BoomerAMG outer relaxation weight for a level', &
                                 'LINEAR SOLVER, HYPRE options')  
              word = trim(word) // ' ' // trim(word2)
              string = trim(prefix) // &
                       'pc_hypre_boomeramg_outer_relax_weight_level'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_NO_CF')
              string = trim(prefix) // 'pc_hypre_boomeramg_no_CF'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),'',ierr);CHKERRQ(ierr)
            case('BOOMERAMG_MEASURE_TYPE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG measure type', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_measure_type'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_COARSEN_TYPE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG coarsen type', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_coarsen_type'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_INTERPOLATION_TYPE','BOOMERAMG_INTERP_TYPE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'BoomerAMG interpolation type', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_interp_type'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),trim(word), &
                                        ierr);CHKERRQ(ierr)
            case('BOOMERAMG_NODAL_COARSEN')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG set nodal coarsening', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_nodal_coarsen'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),'',ierr);CHKERRQ(ierr)
            case('BOOMERAMG_NODAL_RELAXATION')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option, &
                                 'BoomerAMG nodal relaxation via Schwarz', &
                                 'LINEAR SOLVER, HYPRE options')  
              string = trim(prefix) // 'pc_hypre_boomeramg_nodal_relaxation'
              call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                        trim(string),'',ierr);CHKERRQ(ierr)
            case default
              option%io_buffer  = 'HYPRE option: ' // trim(keyword) // &
                                  ' unknown.'
              call PrintErrMsg(option)
          end select
        enddo
        call InputPopBlock(input,option)

      case('ATOL')
        call InputReadDouble(input,option,solver%linear_atol)
        call InputErrorMsg(input,option,'linear_atol','LINEAR_SOLVER')

      case('RTOL')
        call InputReadDouble(input,option,solver%linear_rtol)
        call InputErrorMsg(input,option,'linear_rtol','LINEAR_SOLVER')

      case('DTOL')
        call InputReadDouble(input,option,solver%linear_dtol)
        call InputErrorMsg(input,option,'linear_dtol','LINEAR_SOLVER')
   
      case('MAXIT')
        call InputKeywordDeprecated('MAXIT', &
                                    'MAXIMUM_NUMBER_OF_ITERATIONS',option)

      case('MAXIMUM_NUMBER_OF_ITERATIONS')
        call InputReadInt(input,option,solver%linear_max_iterations)
        call InputErrorMsg(input,option,'linear_max_iterations','LINEAR_SOLVER')

      case('LU_ZERO_PIVOT_TOL')
        call InputReadDouble(input,option,solver%linear_zero_pivot_tol)
        call InputErrorMsg(input,option,'linear_zero_pivot_tol', &
                           'LINEAR_SOLVER')

      case('DISABLE_SHIFT')
        solver%linear_shift = PETSC_FALSE

      case('STOP_ON_FAILURE')
        solver%linear_stop_on_failure = PETSC_TRUE

      case('MUMPS')
        string = trim(prefix) // 'pc_factor_mat_solver_type'
        word = 'mumps'
        call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                  trim(string),trim(word),ierr);CHKERRQ(ierr)

      case ('CPR_OPTIONS')
        call SolverCPRRead(solver%cprstash, input,option, ierr)

      case ('VERBOSE_LOGGING','VERBOSE_ERROR_MESSAGING')
        solver%verbose_logging = PETSC_TRUE 

      case('GMRES_RESTART')
        ! Equivalent to 
        ! -[prefix]_ksp_gmres_restart x
        ! on command line
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option, &
                           'GMRES restart','LINEAR SOLVER')  
        string = trim(prefix) // 'ksp_gmres_restart'
        call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                  trim(string),trim(word), &
                                  ierr);CHKERRQ(ierr)

      case('GMRES_MODIFIED_GS')
        ! Equivalent to 
        ! -[prefix]_ksp_gmres_modifiedgramschmidt
        ! on command line
        string = trim(prefix) // 'ksp_gmres_modifiedgramschmidt'
        word = ''
        call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                                  trim(string),trim(word),ierr);CHKERRQ(ierr)

      case default
        call InputKeywordUnrecognized(input,keyword,'LINEAR_SOLVER',option)
    end select 
  
  enddo 
  call InputPopBlock(input,option)

end subroutine SolverReadLinear

! ************************************************************************** !

subroutine SolverReadNewton(solver,input,option)
  ! 
  ! Reads parameters associated with linear solver
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/21/07, 03/16/20
  ! 

  use Input_Aux_module
  use String_module
  use Option_module
  
  implicit none

  type(solver_type) :: solver
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word, word2
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  error_string = 'SUBSURFACE,NEWTON_SOLVER'

  input%ierr = 0
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword','NEWTON SOLVER')
    call StringToUpper(keyword)   

    found = PETSC_TRUE
    call SolverReadNewtonSelectCase(solver,input,keyword,found, &
                                    error_string,option)
    if (.not.found) then
      call InputKeywordUnrecognized(input,keyword,error_string,option)
    endif
  
  enddo  
  call InputPopBlock(input,option)

end subroutine SolverReadNewton

! ************************************************************************** !

subroutine SolverReadNewtonSelectCase(solver,input,keyword,found, &
                                      error_string,option)
  ! 
  ! Reads keywords specific to the solver object and not process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/21/07, 03/16/20
  ! 

  use Input_Aux_module
  use String_module
  use Option_module
  
  implicit none

  type(solver_type) :: solver
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: word, word2
  PetscBool :: boolean

  found = PETSC_TRUE
  select case(trim(keyword))
  
    case ('INEXACT_NEWTON')
      solver%inexact_newton = PETSC_TRUE

    case ('NO_PRINT_CONVERGENCE')
      solver%print_convergence = PETSC_FALSE

    case ('NO_INF_NORM','NO_INFINITY_NORM')
      solver%check_infinity_norm = PETSC_FALSE

    case('MAXIT')
        call InputKeywordDeprecated('MAXIT', &
                                    'MAXIMUM_NUMBER_OF_ITERATIONS',option)

    case('MAXIMUM_NUMBER_OF_ITERATIONS')
      call InputReadInt(input,option,solver%newton_max_iterations)
      call InputErrorMsg(input,option,'maximum newton iterations',error_string)

    case('MINIMUM_NEWTON_ITERATION','MINIMUM_NEWTON_ITERATIONS')
      call InputReadInt(input,option,solver%newton_min_iterations)
      call InputErrorMsg(input,option,'minimum newton iterations',error_string)

    case ('PRINT_DETAILED_CONVERGENCE')
      solver%print_detailed_convergence = PETSC_TRUE

    case ('PRINT_LINEAR_ITERATIONS')
      solver%print_linear_iterations = PETSC_TRUE

    case('ATOL')
      call InputReadDouble(input,option,solver%newton_atol)
      call InputErrorMsg(input,option,'newton_atol',error_string)

    case('RTOL')
      call InputReadDouble(input,option,solver%newton_rtol)
      call InputErrorMsg(input,option,'newton_rtol',error_string)

    case('STOL')
      call InputReadDouble(input,option,solver%newton_stol)
      call InputErrorMsg(input,option,'newton_stol',error_string)
    
    case('DTOL')
      call InputReadDouble(input,option,solver%newton_dtol)
      call InputErrorMsg(input,option,'newton_dtol',error_string)

    case('MAX_NORM')
      call InputReadDouble(input,option,solver%max_norm)
      call InputErrorMsg(input,option,'max_norm',error_string)
 
    case('ITOL', 'INF_TOL', 'ITOL_RES', 'INF_TOL_RES')
      call InputReadDouble(input,option,solver%newton_inf_res_tol)
      call InputErrorMsg(input,option,'newton_inf_res_tol',error_string)
 
    case('ITOL_UPDATE', 'INF_TOL_UPDATE')
      call InputReadDouble(input,option,solver%newton_inf_upd_tol)
      call InputErrorMsg(input,option,'newton_inf_upd_tol',error_string)

    case('ITOL_SEC','ITOL_RES_SEC','INF_TOL_SEC')
      !TODO(geh): move to PM
      if (.not.option%use_mc) then
        option%io_buffer = 'NEWTON ITOL_SEC not supported without ' // &
          'MULTIPLE_CONTINUUM keyword.'
        call PrintErrMsg(option)
      endif
      if (.not.solver%itype == TRANSPORT_CLASS) then
        option%io_buffer = 'NEWTON ITOL_SEC supported in ' // &
          'TRANSPORT only.'
        call PrintErrMsg(option)
      endif         
      call InputReadDouble(input,option,solver%newton_inf_res_tol_sec)
      call InputErrorMsg(input,option,'newton_inf_res_tol_sec', &
                         error_string)
 
    case('MAXF')
      call InputReadInt(input,option,solver%newton_maxf)
      call InputErrorMsg(input,option,'newton_maxf',error_string)

    case('MATRIX_TYPE')
      call InputReadCard(input,option,word)
      call InputErrorMsg(input,option,'mat_type','NEWTON SOLVER')   
      call StringToUpper(word)
      select case(trim(word))
        case('BAIJ')
          solver%J_mat_type = MATBAIJ
        case('AIJ')
          solver%J_mat_type = MATBAIJ
          solver%J_mat_type = MATAIJ
        case('MFFD','MATRIX_FREE')
          solver%J_mat_type = MATMFFD
        case('HYPRESTRUCT')
          solver%J_mat_type = MATHYPRESTRUCT
        case('SELL')
          solver%J_mat_type = MATSELL
        case default
          option%io_buffer = 'Matrix type: ' // trim(word) // ' unknown.'
          call PrintErrMsg(option)
      end select
      
    case('PRECONDITIONER_MATRIX_TYPE')
      call InputReadCard(input,option,word)
      call InputErrorMsg(input,option,'mat_type','NEWTON SOLVER')   
      call StringToUpper(word)
      select case(trim(word))
        case('BAIJ')
          solver%Jpre_mat_type = MATBAIJ
        case('AIJ')
          solver%Jpre_mat_type = MATBAIJ
          solver%Jpre_mat_type = MATAIJ
        case('MFFD','MATRIX_FREE')
          solver%Jpre_mat_type = MATMFFD
        case('HYPRESTRUCT')
           solver%Jpre_mat_type = MATHYPRESTRUCT
        case('SELL')
          solver%J_mat_type = MATSELL
        case('SHELL')
           solver%Jpre_mat_type = MATSHELL
        case default
          option%io_buffer  = 'Preconditioner Matrix type: ' // &
            trim(word) // ' unknown.'
          call PrintErrMsg(option)
      end select
      
    case ('VERBOSE_LOGGING','VERBOSE_ERROR_MESSAGING')
      solver%verbose_logging = PETSC_TRUE 

    case ('CONVERGENCE_INFO')
      error_string = 'NEWTON_SOLVER,CONVERGENCE_INFO'
      call InputPushBlock(input,option)
      do
        call InputReadPflotranString(input,option)
        if (InputCheckExit(input,option)) exit  
        call InputReadCard(input,option,keyword)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(keyword)
        call InputReadCard(input,option,word)
        call StringToUpper(word)
        select case(StringYesNoOther(word))
          case(STRING_YES)
            boolean = PETSC_TRUE
          case(STRING_NO)
            boolean = PETSC_FALSE
          case(STRING_OTHER)
            error_string = trim(error_string) // ',' // keyword
            call InputKeywordUnrecognized(input,word,error_string,option)
        end select
        select case(trim(keyword))
          case('2R','FNORM','2NORMR')
            solver%convergence_2r = boolean
          case('2X','XNORM','2NORMX')
            solver%convergence_2x = boolean
          case('2U','UNORM','2NORMU')
            solver%convergence_2u = boolean
          case('IR','INORMR')
            solver%convergence_ir = boolean
          case('IU','INORMU')
            solver%convergence_iu = boolean
          case default
            call InputKeywordUnrecognized(input,keyword,error_string,option)
        end select
      enddo
    case default
      found = PETSC_FALSE
  end select 

end subroutine SolverReadNewtonSelectCase

! ************************************************************************** !

subroutine SolverPrintLinearInfo(solver,header,option)
  ! 
  ! Prints information about linear solver
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/23/08
  ! 

  use Option_module
  use String_module
  
  implicit none
  
  type(solver_type) :: solver
  character(len=*) :: header  
  type(option_type) :: option

  PetscInt :: fids(2)
  character(len=MAXSTRINGLENGTH) :: strings(20)
  PetscErrorCode :: ierr

#if !defined(PETSC_HAVE_MUMPS)
  if (option%mycommsize > 1) then
    if (solver%ksp_type == KSPPREONLY .and. solver%pc_type == PCLU) then
      option%io_buffer = 'Direct solver (KSPPREONLY + PCLU) not ' // &
        ' supported when running in parallel.  Switch to SOLVER ITERATIVE.'
      call PrintErrMsg(option)
    endif
  endif
#endif

  fids = OptionGetFIDs(option)
  call StringWriteToUnits(fids,'')
  call StringWriteToUnits(fids,trim(header) // ' Linear Solver')
  strings(:) = ''
  strings(1) = 'solver: '// StringWrite(solver%ksp_type)
  strings(2) = 'preconditioner: '// StringWrite(solver%pc_type)
  strings(3) = 'atol: '// StringWrite(solver%linear_atol)
  strings(4) = 'rtol: '// StringWrite(solver%linear_rtol)
  strings(5) = 'dtol: '// StringWrite(solver%linear_dtol)
  strings(6) = 'maximum iteration: ' // &
                                   StringWrite(solver%linear_max_iterations)
  if (Initialized(solver%linear_zero_pivot_tol)) then
    strings(7) = 'zero pivot tolerance: ' // &
                                   StringWrite(solver%linear_zero_pivot_tol)
  endif
  call StringsCenter(strings,30,':')
  call StringWriteToUnits(fids,strings)

  if (solver%verbose_logging) then
    call KSPView(solver%ksp,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
    call PCView(solver%pc,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
  endif
  

end subroutine SolverPrintLinearInfo

! ************************************************************************** !

subroutine SolverPrintNewtonInfo(solver,header,option)    
  ! 
  ! Prints information about Newton solver
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/23/08
  ! 
  use Option_module
  use String_module

  implicit none
  
  type(solver_type) :: solver
  character(len=*) :: header  
  type(option_type) :: option

  PetscInt :: fids(2)
  character(len=MAXSTRINGLENGTH) :: strings(20)
  PetscErrorCode :: ierr

  fids = OptionGetFIDs(option)
  call StringWriteToUnits(fids,'')
  call StringWriteToUnits(fids,trim(header) // ' Newton Solver')
  strings(:) = ''
  strings(1) = 'atol: '// StringWrite(solver%newton_atol)
  strings(2) = 'rtol: '// StringWrite(solver%newton_rtol)
  strings(3) = 'stol: '// StringWrite(solver%newton_stol)
  strings(4) = 'dtol: '// StringWrite(solver%newton_dtol)
  strings(5) = 'maxnorm: '// StringWrite(solver%max_norm)
  strings(6) = 'inftolres: ' // StringWrite(solver%newton_inf_res_tol)
  strings(7) = 'inftolupd: ' // StringWrite(solver%newton_inf_upd_tol)
  strings(8) = 'inftolrelupd: ' // StringWrite(solver%newton_inf_rel_update_tol)
  strings(9) = 'inftolsclres: ' // StringWrite(solver%newton_inf_scaled_res_tol)
  strings(10) = 'max iter: ' // StringWrite(solver%newton_max_iterations)
  strings(11) = 'min iter: ' // StringWrite(solver%newton_min_iterations)
  strings(12) = 'maxf: ' // StringWrite(solver%newton_maxf)
  strings(13) = 'matrix type: ' // StringWrite(solver%J_mat_type)
  strings(14) = 'precond. matrix type: ' // StringWrite(solver%Jpre_mat_type)
  strings(15) = 'inexact newton: ' // &
                   StringWrite(String1Or2(solver%inexact_newton,'on','off'))
  strings(16) = 'print convergence: ' // &
                StringWrite(String1Or2(solver%print_convergence,'on','off'))
  strings(17) = 'print detailed convergence: ' // &
       StringWrite(String1Or2(solver%print_detailed_convergence,'on','off'))
  strings(18) = 'check infinity norm: ' // &
              StringWrite(String1Or2(solver%check_infinity_norm,'on','off'))

  call StringsCenter(strings,30,':')
  call StringWriteToUnits(fids,strings)

  if (solver%verbose_logging) then
    call SNESView(solver%snes,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
  endif
  
end subroutine SolverPrintNewtonInfo

! ************************************************************************** !

subroutine SolverCheckCommandLine(solver)
  ! 
  ! Parses the command line for various solver
  ! options.
  ! Note: In order to use the PETSc OptionsPrefix associated with
  ! solver%snes in parsing the options, the call to SolverCheckCommandLine()
  ! should come after the SNESSetOptionsPrefix(solver%snes,...) call.
  ! 
  ! Author: Richard Tran Mills
  ! Date: 05/09/2008
  ! 

  implicit none
  
  type(solver_type) :: solver

  PetscErrorCode :: ierr
  character(len=MAXSTRINGLENGTH) :: prefix
  character(len=MAXSTRINGLENGTH) :: mat_type
  PetscBool :: is_present

  if (solver%snes /= PETSC_NULL_SNES) then
    call SNESGetOptionsPrefix(solver%snes, prefix, ierr);CHKERRQ(ierr)
  else
    prefix = PETSC_NULL_CHARACTER
  endif

  ! Parse the options to determine if the matrix type has been specified.
  call PetscOptionsGetString(PETSC_NULL_OPTIONS,prefix, '-mat_type', mat_type, &
                             is_present,ierr);CHKERRQ(ierr)
  if (is_present) solver%J_mat_type = trim(mat_type)
  
  call PetscOptionsGetString(PETSC_NULL_OPTIONS,prefix, '-pre_mat_type', &
                             mat_type, is_present,ierr);CHKERRQ(ierr)
  if (is_present) solver%Jpre_mat_type = trim(mat_type)

  ! Parse the options for the Galerkin multigrid solver.
  ! Users can specify the number of levels of coarsening via the
  ! 'galerkin_mg N' option, which will set the number of levels in the 
  ! x, y, and z directions all to N.  For semi-coarsening, however, 
  ! it is possible to set the number of levels in each direction 
  ! individually via options such as '-galerkin_mg_x N', which would 
  ! override the number of levels in the x direction set by '-galerkin_mg'.
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS,prefix, '-galerkin_mg', &
                          solver%galerkin_mg_levels, solver%use_galerkin_mg, &
                          ierr);CHKERRQ(ierr)
  if (solver%use_galerkin_mg) then
    solver%galerkin_mg_levels_x = solver%galerkin_mg_levels
    solver%galerkin_mg_levels_y = solver%galerkin_mg_levels
    solver%galerkin_mg_levels_z = solver%galerkin_mg_levels
  endif

  call PetscOptionsGetInt(PETSC_NULL_OPTIONS,prefix, '-galerkin_mg_x', &
                          solver%galerkin_mg_levels_x, is_present, &
                          ierr);CHKERRQ(ierr)
  if (is_present) solver%use_galerkin_mg = PETSC_TRUE
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS,prefix, '-galerkin_mg_y', &
                          solver%galerkin_mg_levels_y, is_present, &
                          ierr);CHKERRQ(ierr)
  if (is_present) solver%use_galerkin_mg = PETSC_TRUE
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS,prefix, '-galerkin_mg_z', &
                          solver%galerkin_mg_levels_z, is_present, &
                          ierr);CHKERRQ(ierr)
  if (is_present) solver%use_galerkin_mg = PETSC_TRUE

  if (solver%use_galerkin_mg) then
    solver%J_mat_type = MATAIJ
      ! Must use AIJ above, as BAIJ is not supported for Galerkin MG solver.
    solver%galerkin_mg_levels = max(solver%galerkin_mg_levels_x, &
                                    solver%galerkin_mg_levels_y, &
                                    solver%galerkin_mg_levels_z)
  endif
                             

end subroutine SolverCheckCommandLine

! ************************************************************************** !

subroutine SolverNewtonPrintFailedReason(solver,option)    
  ! 
  ! Prints the reason for the solver failing
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/02/16
  ! 
  use Option_module

  implicit none
  
  type(solver_type) :: solver
  type(option_type) :: option

  SNESConvergedReason :: snes_reason
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscReal :: abstol, rtol, stol, divtol
  PetscInt :: maxit, maxf
  PetscErrorCode :: ierr

  call SNESGetConvergedReason(solver%snes,snes_reason,ierr);CHKERRQ(ierr)
  call SNESGetTolerances(solver%snes,abstol,rtol,stol,maxit,maxf, &
                         ierr);CHKERRQ(ierr)
  select case(snes_reason)
    case(SNES_DIVERGED_FUNCTION_DOMAIN)
      if (solver%verbose_logging) then
        error_string = 'The new solution location passed to the function &
          &is not in the domain of F.'
      else
        error_string = 'SNES_DIVERGED_FUNCTION_DOMAIN'
      endif
    case(SNES_DIVERGED_FUNCTION_COUNT)
      if (solver%verbose_logging) then
        write(word,*) maxf
        error_string = 'The user provided function has been called &
          &more times than the final argument (' // trim(adjustl(word)) // &
          ') in SNESSetTolerances().'
      else
        error_string = 'SNES_DIVERGED_FUNCTION_COUNT'
      endif
    case(SNES_DIVERGED_LINEAR_SOLVE)
      if (solver%verbose_logging) then
        error_string = 'The linear solver failed.'
      else
        error_string = 'SNES_DIVERGED_LINEAR_SOLVE'
      endif
    case(SNES_DIVERGED_FNORM_NAN)
      if (solver%verbose_logging) then
        error_string = 'The norm of the residual is &
          &not a number (NaN). It is likely that the residual has NaNs &
          &in it.  This could be caused by errors in boundary conditions, &
          &using a constitutive relation evaluated outside prescribed &
          &bounds or a bug.'
      else
        error_string = 'SNES_DIVERGED_FNORM_NAN'
      endif
    case(SNES_DIVERGED_MAX_IT)
      if (solver%verbose_logging) then
        write(word,*) maxit
        error_string = 'The maximum number of Newton iterations (' // &
        trim(adjustl(word)) // ') was reached.'
      else
        error_string = 'SNES_DIVERGED_MAX_IT'
      endif
    case(SNES_DIVERGED_LINE_SEARCH)
      if (solver%verbose_logging) then
        error_string = 'The line search failed.'
      else
        error_string = 'SNES_DIVERGED_LINE_SEARCH'
      endif
    case(SNES_DIVERGED_INNER)
      if (solver%verbose_logging) then
        error_string = 'The inner solver failed.'
      else
        error_string = 'SNES_DIVERGED_INNER'
      endif
    case(SNES_DIVERGED_LOCAL_MIN)
      if (solver%verbose_logging) then
        error_string = '|| J^T b || is small, implying convergence to a &
          &local minimum of F().' 
      else
        error_string = 'SNES_DIVERGED_LOCAL_MIN'
      endif
    case(SNES_DIVERGED_DTOL)
      if (solver%verbose_logging) then
        call SNESGetDivergenceTolerance(solver%snes,divtol,ierr);CHKERRQ(ierr)
        write(word,'(es12.4)') divtol
        error_string = 'The nonlinear residual has diverged based on &
          &||F|| > divtol*||F_initial||, where divtol = ' // &
          trim(adjustl(word)) // '.'
      else
        error_string = 'SNES_DIVERGED_DTOL'
      endif
    case default
      write(word,*) snes_reason 
      error_string = 'Unknown(' // &
        trim(adjustl(word)) // ').'
  end select
  option%io_buffer = 'Newton solver reason: ' // trim(error_string)
  call PrintMsg(option)

  ! print out subsequent information specific to each case
  select case(snes_reason)
    case(SNES_DIVERGED_LINEAR_SOLVE)
      call SolverLinearPrintFailedReason(solver,option)    
  end select

end subroutine SolverNewtonPrintFailedReason

! ************************************************************************** !

subroutine SolverLinearPrintFailedReason(solver,option)    
  ! 
  ! Prints the reason for the solver failing
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/02/16
  ! 
  use Option_module

  implicit none
  
  type(solver_type) :: solver
  type(option_type) :: option

  KSP, pointer :: sub_ksps(:)
  PC :: pc
  Mat :: mat
  PCType :: pc_type
  PetscInt :: i
  PetscInt :: nsub_ksp
  PetscInt :: first_sub_ksp
  KSPConvergedReason :: ksp_reason
  PCFailedReason :: pc_failed_reason, global_pc_failed_reason
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: zero_pivot_tol, zero_pivot
  character(len=MAXWORDLENGTH) :: word, word2
  PetscInt :: irow, temp_int
  PetscReal :: rtol, abstol, dtol
  PetscReal :: rnorm
  PetscInt :: maxits
  PetscErrorCode :: ierr

  call KSPGetConvergedReason(solver%ksp,ksp_reason,ierr);CHKERRQ(ierr)
  call KSPGetTolerances(solver%ksp,rtol,abstol,dtol,maxits,ierr);CHKERRQ(ierr)
  select case(ksp_reason)
    case(KSP_DIVERGED_ITS)
      if (solver%verbose_logging) then
        write(word,*) maxits
        error_string = 'The linear solver took too many iterations, beyond &
          &the allowable number set by maxits ('// trim(adjustl(word)) // ').'
      else
        error_string = 'KSP_DIVERGED_ITS'
      endif
    case(KSP_DIVERGED_DTOL)
      if (solver%verbose_logging) then
        call KSPGetResidualNorm(solver%ksp,rnorm,ierr);CHKERRQ(ierr)
        write(word,'(es12.4)') rnorm
        error_string = 'The linear solver diverged due to dtol based on &
          &the equation: norm(r) >= dtol*norm(b) where r = b-Ax for the &
          &linear system Ax=b within the Krylov solver. norm(r) = ' // &
          trim(adjustl(word)) // ').'
      else
        error_string = 'KSP_DIVERGED_DTOL'
      endif
    case(KSP_DIVERGED_BREAKDOWN)
      if (solver%verbose_logging) then
        error_string = 'A breakdown in the Krylov method was detected &
          &so the method could not continue to enlarge the Krylov space. &
          &Could be due to a singlular matrix or preconditioner.'
      else
        error_string = 'KSP_DIVERGED_BREAKDOWN'
      endif
    case(KSP_DIVERGED_BREAKDOWN_BICG)
      if (solver%verbose_logging) then
        error_string = 'A breakdown in the KSPBICG method was detected &
          &so the method could not continue to enlarge the Krylov space.'
      else
        error_string = 'KSP_DIVERGED_BREAKDOWN_BICG'
      endif
    case(KSP_DIVERGED_NONSYMMETRIC)
      if (solver%verbose_logging) then
        ! must use '"' instead of "'" due to parentheses
        error_string = "It appears the operator or preconditioner is not &
          &symmetric and this Krylov method (KSPCG, KSPMINRES, KSPCR)&
          &requires symmetry."
      else
        error_string = 'KSP_DIVERGED_NONSYMMETRIC'
      endif
    case(KSP_DIVERGED_INDEFINITE_PC)
      if (solver%verbose_logging) then
        ! must use '"' instead of "'" due to parentheses
        error_string = "It appears the preconditioner is indefinite (has &
          &both positive and negative eigenvalues) and this Krylov method &
          &(KSPCG) requires it to be positive definite. This can happen &
          &with the PCICC preconditioner, use -pc_factor_shift_positive_&
          &definite to force the PCICC preconditioner to generate a &
          &positive definite preconditioner."
      else
        error_string = 'KSP_DIVERGED_INDEFINITE_PC'
      endif
    case(KSP_DIVERGED_NANORINF)
      if (solver%verbose_logging) then
        error_string = 'The linear solver produced a NaN (not a number) &
          for Inf (infinite number) likely due to a divide by zero (0/0).'
      else
        error_string = 'KSP_DIVERGED_NANORINF'
      endif
    case(KSP_DIVERGED_INDEFINITE_MAT)
      if (solver%verbose_logging) then
        error_string = 'The linear solver failed due to an indefinite matrix.'
      else
        error_string = 'KSP_DIVERGED_INDEFINITE_MAT'
      endif
    case(KSP_DIVERGED_PC_FAILED)
      if (solver%verbose_logging) then
        error_string = 'Preconditioner setup failed.'
        pc = solver%pc
        call PCGetType(pc,pc_type,ierr);CHKERRQ(ierr)
        call PCGetFailedReason(pc,pc_failed_reason, &
                               ierr);CHKERRQ(ierr)
        ! have to perform global reduction on pc_failed_reason
        temp_int = pc_failed_reason
        call MPI_Allreduce(MPI_IN_PLACE,temp_int,ONE_INTEGER_MPI, &
                           MPI_INTEGER,MPI_MAX,option%mycomm,ierr)
        global_pc_failed_reason = temp_int
        if (global_pc_failed_reason == PC_SUBPC_ERROR) then
          if (pc_type == PCBJACOBI) then
            call PCBJacobiGetSubKSP(pc,nsub_ksp,first_sub_ksp, &
                                    PETSC_NULL_KSP,ierr);CHKERRQ(ierr)
            allocate(sub_ksps(nsub_ksp))
            call PCBJacobiGetSubKSP(pc,nsub_ksp,first_sub_ksp, &
                                    sub_ksps,ierr);CHKERRQ(ierr)
            if (nsub_ksp > 1) then
              option%io_buffer = 'NSUB_KSP > 1.  What to do?'
              call PrintErrMsgToDev(option,'')
            endif
            do i = 1, nsub_ksp
              call KSPGetPC(sub_ksps(i),pc,ierr);CHKERRQ(ierr)
              call PCGetFailedReason(pc,pc_failed_reason, &
                                     ierr);CHKERRQ(ierr)
            enddo
            deallocate(sub_ksps)
            nullify(sub_ksps)
          else
            option%io_buffer = 'Error in SUB PC of unknown type "' // &
              trim(pc_type) // '".'
            call PrintErrMsg(option)
          endif
        endif
        ! have to perform global reduction (again) on pc_failed_reason
        temp_int = pc_failed_reason
        call MPI_Allreduce(MPI_IN_PLACE,temp_int,ONE_INTEGER_MPI, &
                           MPI_INTEGER,MPI_MAX,option%mycomm,ierr)
        global_pc_failed_reason = temp_int
        select case(global_pc_failed_reason)
          case(PC_FACTOR_STRUCT_ZEROPIVOT,PC_FACTOR_NUMERIC_ZEROPIVOT)
            select case(solver%itype)
              case(FLOW_CLASS)
                string = 'Flow'
              case(TRANSPORT_CLASS)
                string = 'Transport'
            end select
            call PCFactorGetZeroPivot(pc,zero_pivot_tol, &
                                      ierr);CHKERRQ(ierr)
            write(word,*) zero_pivot_tol
            ! In parallel, some processes will not have a zero pivot and
            ! will report zero as the error.  We must skip these processes.
            zero_pivot = 1.d20
            ! note that this is not the global pc reason
            select case(pc_failed_reason)
              case(PC_FACTOR_STRUCT_ZEROPIVOT,PC_FACTOR_NUMERIC_ZEROPIVOT)
              call PCFactorGetMatrix(pc,mat,ierr);CHKERRQ(ierr)
              call MatFactorGetErrorZeroPivot(mat,zero_pivot,irow, &
                                              ierr);CHKERRQ(ierr)
            end select
            call MPI_Allreduce(MPI_IN_PLACE,zero_pivot,ONE_INTEGER_MPI, &
                               MPI_DOUBLE_PRECISION,MPI_MIN,option%mycomm,ierr)
            write(word2,*) zero_pivot
            option%io_buffer = 'PC Setup failed for ' // trim(string) // &
              '. The ' // trim(string) // ' preconditioner zero pivot &
              &tolerance (' // trim(adjustl(word)) // &
              ') is too large due to a zero pivot of ' // &
              trim(adjustl(word2)) // '. Please set a ZERO_PIVOT_TOL smaller &
              &than that value.'
            call PrintErrMsgToDev(option,'if you need further help')
        end select
      else
        error_string = 'KSP_DIVERGED_PC_FAILED'
      endif
    case default
      write(word,*) ksp_reason 
      option%io_buffer = 'Unknown(' // &
        trim(adjustl(word)) // ').'
  end select
  option%io_buffer = 'Linear solver reason: ' // trim(error_string)
  call PrintMsg(option)

end subroutine SolverLinearPrintFailedReason

! ************************************************************************** !

subroutine SolverDestroy(solver)
  ! 
  ! Deallocates a solver
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 
  
  use String_module
   
  implicit none
  
  type(solver_type), pointer :: solver
  
  PetscErrorCode :: ierr
  PetscInt :: i

  if (.not.associated(solver)) return

  if (solver%Jpre == solver%J) then
    solver%Jpre = PETSC_NULL_MAT
  else if (solver%Jpre /= PETSC_NULL_MAT) then
    call MatDestroy(solver%Jpre,ierr);CHKERRQ(ierr)
  endif
  if (solver%J /= PETSC_NULL_MAT) then
    call MatDestroy(solver%J,ierr);CHKERRQ(ierr)
  endif
  if (associated(solver%interpolation)) then
    do i=1,solver%galerkin_mg_levels-1
      call MatDestroy(solver%interpolation(i),ierr);CHKERRQ(ierr)
    enddo
    deallocate(solver%interpolation)
  endif
  if (solver%matfdcoloring /= PETSC_NULL_MATFDCOLORING) then
    call MatFDColoringDestroy(solver%matfdcoloring,ierr);CHKERRQ(ierr)
  endif

  ! the highest level object frees everything within
  if (solver%ts /= PETSC_NULL_TS) then
    call TSDestroy(solver%ts,ierr);CHKERRQ(ierr)
  else if (solver%snes /= PETSC_NULL_SNES) then
    call SNESDestroy(solver%snes,ierr);CHKERRQ(ierr)
  else if (solver%ksp /= PETSC_NULL_KSP) then
    call KSPDestroy(solver%ksp,ierr);CHKERRQ(ierr)
  endif

  solver%ts = PETSC_NULL_TS
  solver%snes = PETSC_NULL_SNES
  solver%ksp = PETSC_NULL_KSP
  solver%pc = PETSC_NULL_PC

  if (associated(solver%cprstash)) then

      call DeallocateWorkersInCPRStash(solver%cprstash)

      if (solver%cprstash%T1_KSP /= PETSC_NULL_KSP) then
        call KSPDestroy(solver%cprstash%T1_KSP, ierr);CHKERRQ(ierr)
      endif
      if (solver%cprstash%T3_KSP /= PETSC_NULL_KSP .and. &
          solver%cprstash%CPR_type == "ADDITIVE") then
        call KSPDestroy(solver%cprstash%T3_KSP, ierr);CHKERRQ(ierr)
      endif
      if (solver%cprstash%T1_PC /= PETSC_NULL_PC) then
        call PCDestroy(solver%cprstash%T1_PC, ierr);CHKERRQ(ierr)
      endif
      if (solver%cprstash%T2_PC /= PETSC_NULL_PC) then
        call PCDestroy(solver%cprstash%T2_PC, ierr);CHKERRQ(ierr)
      endif
      if (solver%cprstash%T3_PC /= PETSC_NULL_PC .and. &
          solver%cprstash%CPR_type == "ADDITIVE") then
        call PCDestroy(solver%cprstash%T3_PC, ierr);CHKERRQ(ierr)
      endif
      if (solver%cprstash%Ap /= PETSC_NULL_MAT) then
        call MatDestroy(solver%cprstash%Ap, ierr);CHKERRQ(ierr)
      endif
      if (solver%cprstash%As /= PETSC_NULL_MAT) then
        call MatDestroy(solver%cprstash%As, ierr);CHKERRQ(ierr)
      endif
      if (solver%cprstash%T1r /= PETSC_NULL_VEC) then
        call VecDestroy(solver%cprstash%T1r, ierr);CHKERRQ(ierr)
      endif
      if (solver%cprstash%T3r /= PETSC_NULL_VEC) then
        call VecDestroy(solver%cprstash%T3r, ierr);CHKERRQ(ierr)
      endif
      if (solver%cprstash%r2 /= PETSC_NULL_VEC) then
        call VecDestroy(solver%cprstash%r2, ierr);CHKERRQ(ierr)
      endif
      if (solver%cprstash%s /= PETSC_NULL_VEC) then
        call VecDestroy(solver%cprstash%s, ierr);CHKERRQ(ierr)
      endif
      if (solver%cprstash%z /= PETSC_NULL_VEC) then
        call VecDestroy(solver%cprstash%z, ierr);CHKERRQ(ierr)
      endif
      if (solver%cprstash%factors1vec /= PETSC_NULL_VEC) then
        call VecDestroy(solver%cprstash%factors1vec, ierr);CHKERRQ(ierr)
      endif
      if (solver%cprstash%factors2vec /= PETSC_NULL_VEC) then
        call VecDestroy(solver%cprstash%factors2vec, ierr);CHKERRQ(ierr)
      endif
      if (solver%cprstash%factors3vec /= PETSC_NULL_VEC) then
        call VecDestroy(solver%cprstash%factors3vec, ierr);CHKERRQ(ierr)
      endif
      deallocate(solver%cprstash)
      nullify(solver%cprstash)
    endif
    
  deallocate(solver)
  nullify(solver)
  
end subroutine SolverDestroy
  
end module Solver_module
