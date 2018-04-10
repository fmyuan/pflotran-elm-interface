module PM_WIPP_Flow_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class
  use PM_WIPP_SrcSink_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscInt, parameter :: FORCE_ITERATION = 1
  PetscInt, parameter :: OUTSIDE_BOUNDS = 2
  PetscInt, parameter :: MAX_NORMAL_RES_LIQ = 3
  PetscInt, parameter :: MAX_NORMAL_RES_GAS = 4
  PetscInt, parameter :: MAX_RES_LIQ = 5
  PetscInt, parameter :: MAX_RES_GAS = 6
  PetscInt, parameter :: MAX_REL_CHANGE_LIQ_PRES_NI = 7
  PetscInt, parameter :: MAX_CHANGE_LIQ_PRES_NI = 8
  PetscInt, parameter :: MAX_CHANGE_GAS_SAT_NI = 9
  PetscInt, parameter :: MAX_CHANGE_GAS_SAT_NI_TRACK = 10
  PetscInt, parameter :: MAX_CHANGE_GAS_SAT_TS = 11
  PetscInt, parameter :: MAX_CHANGE_LIQ_PRES_TS = 12
  PetscInt, parameter :: MAX_REL_CHANGE_LIQ_PRES_TS = 13
  ! these must be the last two due to the need to calculate the minimum
  PetscInt, parameter :: MIN_LIQ_PRES = 14
  PetscInt, parameter :: MIN_GAS_PRES = 15

  type, public, extends(pm_subsurface_flow_type) :: pm_wippflo_type
    PetscInt, pointer :: max_change_ivar(:)
    PetscReal :: liquid_equation_tolerance
    PetscReal :: gas_equation_tolerance
    PetscReal :: liquid_pressure_tolerance
    PetscReal :: gas_saturation_tolerance
    PetscReal :: dsatlim
    PetscReal :: dprelim
    PetscReal :: satlimit
    PetscReal :: dsat_max  ! can this be this%saturation_change_limit?
    PetscReal :: dpres_max ! can this be this%pressure_change_limit?
    PetscReal :: eps_sat
    PetscReal :: eps_pres
    PetscReal :: satnorm
    PetscReal :: presnorm
    PetscReal :: tswitch
    PetscReal :: dtimemax
    PetscReal :: deltmin
    PetscInt :: iconvtest
    class(pm_wipp_srcsink_type), pointer :: pmwss_ptr
    Vec :: stored_residual_vec
    PetscInt :: convergence_flags(MIN_GAS_PRES)
    ! store maximum quantities for the above
    PetscReal :: convergence_reals(MIN_GAS_PRES)
  contains
    procedure, public :: Read => PMWIPPFloRead
    procedure, public :: InitializeRun => PMWIPPFloInitializeRun
    procedure, public :: InitializeTimestep => PMWIPPFloInitializeTimestep
    procedure, public :: Residual => PMWIPPFloResidual
    procedure, public :: Jacobian => PMWIPPFloJacobian
    procedure, public :: UpdateTimestep => PMWIPPFloUpdateTimestep
    procedure, public :: FinalizeTimestep => PMWIPPFloFinalizeTimestep
    procedure, public :: PreSolve => PMWIPPFloPreSolve
    procedure, public :: PostSolve => PMWIPPFloPostSolve
    procedure, public :: CheckUpdatePre => PMWIPPFloCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMWIPPFloCheckUpdatePost
    procedure, public :: CheckConvergence => PMWIPPFloConvergence
    procedure, public :: TimeCut => PMWIPPFloTimeCut
    procedure, public :: UpdateSolution => PMWIPPFloUpdateSolution
    procedure, public :: UpdateAuxVars => PMWIPPFloUpdateAuxVars
    procedure, public :: MaxChange => PMWIPPFloMaxChange
    procedure, public :: ComputeMassBalance => PMWIPPFloComputeMassBalance
    procedure, public :: InputRecord => PMWIPPFloInputRecord
    procedure, public :: CheckpointBinary => PMWIPPFloCheckpointBinary
    procedure, public :: RestartBinary => PMWIPPFloRestartBinary
    procedure, public :: Destroy => PMWIPPFloDestroy
  end type pm_wippflo_type
  
  public :: PMWIPPFloCreate, &
            PMWIPPFloInitObject, &
            PMWIPPFloReadSelectCase, &
            PMWIPPFloInitializeRun, &
            PMWIPPFloFinalizeTimestep, &
            PMWIPPFloCheckUpdatePre, &
            PMWIPPFloDestroy
  
contains

! ************************************************************************** !

function PMWIPPFloCreate()
  ! 
  ! Creates WIPPFlo process models shell
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
  use Variables_module, only : LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
                               LIQUID_MOLE_FRACTION, TEMPERATURE, &
                               GAS_SATURATION
  implicit none
  
  class(pm_wippflo_type), pointer :: PMWIPPFloCreate

  class(pm_wippflo_type), pointer :: wippflo_pm
  
  allocate(wippflo_pm)
  call PMWIPPFloInitObject(wippflo_pm)

  PMWIPPFloCreate => wippflo_pm
  
end function PMWIPPFloCreate

! ************************************************************************** !

subroutine PMWIPPFloInitObject(this)
  ! 
  ! Creates WIPPFlo process models shell
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/17
  ! 
  use Variables_module, only : LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
                               LIQUID_MOLE_FRACTION, TEMPERATURE, &
                               GAS_SATURATION
  implicit none
  
  class(pm_wippflo_type) :: this
  
  allocate(this%max_change_ivar(3))
  this%max_change_ivar = [LIQUID_PRESSURE, GAS_PRESSURE, GAS_SATURATION]
  nullify(this%pmwss_ptr)
  
  call PMSubsurfaceFlowCreate(this)
  this%name = 'WIPP Immiscible Multiphase Flow'

  this%check_post_convergence = PETSC_TRUE

  ! defaults from BRAGFLO input deck or recommended values from user manual
  this%liquid_equation_tolerance = 1.d-2
  this%gas_equation_tolerance = 1.d-2
  this%liquid_pressure_tolerance = 1.d-2
  this%gas_saturation_tolerance = 1.d-3
  this%dsatlim = 0.20d0   ! [-]
  this%dprelim = -1.0d8   ! [Pa]
  this%satlimit = 1.0d-3  ! [-]
  this%dsat_max = 1.d0    ! [-]
  this%dpres_max = 1.d7   ! [Pa]
  this%eps_sat = 3.0d0    ! [-]
  this%eps_pres = 1.0d-3  ! [-]
  this%satnorm = 3.d-1    ! [-]
  this%presnorm = 5.d5    ! [Pa]
  this%tswitch = 0.01d0   ! [-]
  this%dtimemax = 1.25    ! [-]
  this%deltmin = 8.64d-4  ! [sec]
  this%stored_residual_vec = PETSC_NULL_VEC
  this%iconvtest = 1      ! 0 = either, 1 = both
  this%convergence_flags = 0
  this%convergence_reals = 0.d0

end subroutine PMWIPPFloInitObject

! ************************************************************************** !

subroutine PMWIPPFloRead(this,input)
  ! 
  ! Sets up SNES solvers.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  !
  use WIPP_Flow_module
  use Input_Aux_module
  use String_module
  use Option_module

  implicit none
  
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: keyword, word
  class(pm_wippflo_type) :: this
  type(option_type), pointer :: option
  PetscReal :: tempreal
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  option => this%option

  error_string = 'WIPP Flow Options'
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)
    
    found = PETSC_FALSE
    call PMWIPPFloReadSelectCase(this,input,keyword,found, &
                                  error_string,option)    
    if (found) cycle
    
    select case(keyword)
      case default
        call InputKeywordUnrecognized(keyword,'WIPP Flow Mode',option)
    end select
    
  enddo  
  
  ! Check that SATLIMIT is smaller than DSATLIM
  if (this%satlimit > this%dsatlim) then
    option%io_buffer = 'The value of DSATLIM must be larger than SATLIMIT.'
    call printErrMsg(option)
  endif
  ! Check the sign of given variables
  if (this%dsatlim < 0.d0) then
    option%io_buffer = 'The value of DSATLIM must be positive.'
    call printErrMsg(option)
  endif
  if (this%dprelim > 0.d0) then
    option%io_buffer = 'The value of DPRELIM must be negative.'
    call printErrMsg(option)
  endif
  if (this%satlimit < 0.d0) then
    option%io_buffer = 'The value of SATLIMIT must be positive.'
    call printErrMsg(option)
  endif
  ! Assign tightest tolerence to EPS_SAT
  ! This code should be removed when GAS_SATURATION_TOLERANCE is removed because
  ! this%gas_saturation_tolerance will no longer exist
  this%eps_sat = max(this%eps_sat,(-1.d0*log10(this%gas_saturation_tolerance)))
   
  
  
end subroutine PMWIPPFloRead

! ************************************************************************** !

subroutine PMWIPPFloReadSelectCase(this,input,keyword,found, &
                                   error_string,option)
  ! 
  ! Reads input file parameters associated with the subsurface flow process 
  !       model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/05/16
  !
  use WIPP_Flow_Aux_module
  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  class(pm_wippflo_type) :: this
  type(input_type) :: input

  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option

  found = PETSC_FALSE
  call PMSubsurfaceFlowReadSelectCase(this,input,keyword,found, &
                                      error_string,option)
  if (found) return

  found = PETSC_TRUE
  select case(trim(keyword))
    case('LIQUID_EQUATION_TOLERANCE')
      call InputReadDouble(input,option,this%liquid_equation_tolerance)
      call InputDefaultMsg(input,option,'LIQUID_EQUATION_TOLERANCE')
    case('GAS_EQUATION_TOLERANCE')
      call InputReadDouble(input,option,this%gas_equation_tolerance)
      call InputDefaultMsg(input,option,'GAS_EQUATION_TOLERANCE')
    case('LIQUID_PRESSURE_TOLERANCE')
      call InputReadDouble(input,option,this%liquid_pressure_tolerance)
      call InputDefaultMsg(input,option,'LIQUID_PRESSURE_TOLERANCE')
    case('GAS_SATURATION_TOLERANCE')
      call InputReadDouble(input,option,this%gas_saturation_tolerance)
      call InputDefaultMsg(input,option,'GAS_SATURATION_TOLERANCE')
    case('GAS_COMPONENT_FORMULA_WEIGHT')
      call InputReadDouble(input,option,fmw_comp(2))
      call InputErrorMsg(input,option,'gas component formula wt.', &
                         error_string)
    case('MAXIMUM_PRESSURE_CHANGE')
      call InputReadDouble(input,option,wippflo_max_pressure_change)
      call InputErrorMsg(input,option,'maximum pressure change', &
                         error_string)
    case('MAX_ITERATION_BEFORE_DAMPING')
      call InputReadInt(input,option,wippflo_max_it_before_damping)
      call InputErrorMsg(input,option,'maximum iteration before damping', &
                         error_string)
    case('DAMPING_FACTOR')
      call InputReadDouble(input,option,wippflo_damping_factor)
      call InputErrorMsg(input,option,'damping factor',error_string)
    case('FIX_UPWIND_DIRECTION')
      wippflo_fix_upwind_direction = PETSC_TRUE
    case('UNFIX_UPWIND_DIRECTION')
      wippflo_fix_upwind_direction = PETSC_FALSE
    case('COUNT_UPWIND_DIRECTION_FLIP')
      wippflo_count_upwind_dir_flip = PETSC_TRUE
    case('UPWIND_DIR_UPDATE_FREQUENCY')
      call InputReadInt(input,option,wippflo_upwind_dir_update_freq)
      call InputErrorMsg(input,option,'upwind direction update frequency', &
                         error_string)
    case('NO_FRACTURE')
      wippflo_use_fracture = PETSC_FALSE
    case('NO_CREEP_CLOSURE')
      wippflo_use_creep_closure = PETSC_FALSE
    case('NO_GAS_GENERATION')
      wippflo_use_gas_generation = PETSC_FALSE
    case('BRAGFLO_RESIDUAL_UNITS')
      wippflo_use_bragflo_units = PETSC_TRUE
    case('DEBUG')
      wippflo_debug = PETSC_TRUE
    case('DEBUG_GAS_GENERATION')
      wippflo_debug_gas_generation = PETSC_TRUE
    case('DEBUG_FIRST_ITERATION')
      wippflo_debug = PETSC_TRUE
      wippflo_debug_first_iteration = PETSC_TRUE
    case('DEBUG_OSCILLATORY_BEHAVIOR')
      wippflo_check_oscillatory_behavior = PETSC_TRUE
    case('DEBUG_TS_UPDATE')
      wippflo_debug_ts_update = PETSC_TRUE
    case('MATCH_BRAGFLO_OUTPUT')
      wippflo_match_bragflo_output = PETSC_TRUE
    case('USE_LEGACY_PERTURBATION')
      wippflo_use_legacy_perturbation = PETSC_TRUE
    case('USE_BRAGFLO_CC')
      wippflo_use_bragflo_cc = PETSC_TRUE
    case('PRESSURE_REL_PERTURBATION')
      call InputReadDouble(input,option,wippflo_pres_rel_pert)
      call InputErrorMsg(input,option,'pressure relative perturbation', &
                         error_string)
    case('PRESSURE_MIN_PERTURBATION')
      call InputReadDouble(input,option,wippflo_pres_min_pert)
      call InputErrorMsg(input,option,'pressure minimum perturbation', &
                         error_string)
    case('SATURATION_REL_PERTURBATION')
      call InputReadDouble(input,option,wippflo_sat_rel_pert)
      call InputErrorMsg(input,option,'saturation relative perturbation', &
                         error_string)
    case('SATURATION_MIN_PERTURBATION')
      call InputReadDouble(input,option,wippflo_sat_min_pert)
      call InputErrorMsg(input,option,'saturation minimum perturbation', &
                         error_string)
    case('DSATLIM')
      call InputReadDouble(input,option,this%dsatlim)
      call InputDefaultMsg(input,option,'DSATLIM')
    case('DPRELIM')
      call InputReadDouble(input,option,this%dprelim)
      call InputDefaultMsg(input,option,'DPRELIM')
    case('SATLIMIT')
      call InputReadDouble(input,option,this%satlimit)
      call InputDefaultMsg(input,option,'SATLIMIT')
    case('DSAT_MAX')
      call InputReadDouble(input,option,this%dsat_max)
      call InputDefaultMsg(input,option,'DSAT_MAX')
    case('DPRES_MAX')
      call InputReadDouble(input,option,this%dpres_max)
      call InputDefaultMsg(input,option,'DPRES_MAX')
    case('EPS_SAT')
      call InputReadDouble(input,option,this%eps_sat)
      call InputDefaultMsg(input,option,'EPS_SAT')
    case('EPS_PRES')
      call InputReadDouble(input,option,this%eps_pres)
      call InputDefaultMsg(input,option,'EPS_PRES')
    case('SATNORM')
      call InputReadDouble(input,option,this%satnorm)
      call InputDefaultMsg(input,option,'SATNORM')
    case('PRESNORM')
      call InputReadDouble(input,option,this%presnorm)
      call InputDefaultMsg(input,option,'PRESNORM')
    case('TSWITCH')
      call InputReadDouble(input,option,this%tswitch)
      call InputDefaultMsg(input,option,'TSWITCH')
    case('DTIMEMAX')
      call InputReadDouble(input,option,this%dtimemax)
      call InputDefaultMsg(input,option,'DTIMEMAX')
    case('DELTMIN')
      call InputReadDouble(input,option,this%deltmin)
      call InputDefaultMsg(input,option,'DELTMIN')
    case('ICONVTEST')
      call InputReadInt(input,option,this%iconvtest)
      call InputDefaultMsg(input,option,'ICONVTEST')
    case('RESIDUAL_TEST')
      wippflo_residual_test = PETSC_TRUE
    case('RESIDUAL_TEST_CELL')
      call InputReadInt(input,option,wippflo_residual_test_cell)
      call InputErrorMsg(input,option,'residual test dof', error_string)
    case('JACOBIAN_TEST')
      wippflo_jacobian_test = PETSC_TRUE
    case('JACOBIAN_TEST_RDOF')
      call InputReadInt(input,option,wippflo_jacobian_test_rdof)
      call InputErrorMsg(input,option,'jacobian test rdof', error_string)
    case('JACOBIAN_TEST_XDOF')
      call InputReadInt(input,option,wippflo_jacobian_test_xdof)
      call InputErrorMsg(input,option,'jacobian test xdof', error_string)
    case('NO_ACCUMULATION')
      wippflo_calc_accum = PETSC_FALSE
    case('NO_FLUX')
      wippflo_calc_flux = PETSC_FALSE
    case('NO_BCFLUX')
      wippflo_calc_bcflux = PETSC_FALSE
    case('NO_CHEMISTRY')
      wippflo_calc_chem = PETSC_FALSE
    case('PRINT_RESIDUAL')
      wippflo_print_residual = PETSC_TRUE
    case('PRINT_SOLUTION')
      wippflo_print_solution = PETSC_TRUE
    case('PRINT_UPDATE')
      wippflo_print_update = PETSC_TRUE
    case('ALLOW_NEGATIVE_GAS_PRESSURE')
      wippflo_allow_neg_gas_pressure = PETSC_TRUE
    case default
      found = PETSC_FALSE
  end select

end subroutine PMWIPPFloReadSelectCase

! ************************************************************************** !

recursive subroutine PMWIPPFloInitializeRun(this)
  ! 
  ! Initializes the WIPP_FLOW mode run.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17

  use Realization_Base_class
  use WIPP_Flow_Aux_module
  use Input_Aux_module
  
  implicit none
  
  class(pm_wippflo_type) :: this
  
  PetscInt :: i
  PetscErrorCode :: ierr
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: block_string

  ! need to allocate vectors for max change
  call VecDuplicateVecsF90(this%realization%field%work,SIX_INTEGER, &
                           this%realization%field%max_change_vecs, &
                           ierr);CHKERRQ(ierr)
  ! set initial values
  do i = 1, 3
    call RealizationGetVariable(this%realization, &
                                this%realization%field%max_change_vecs(i), &
                                this%max_change_ivar(i),ZERO_INTEGER)
  enddo

  ! call parent implementation
  call PMSubsurfaceFlowInitializeRun(this)
  
  ! look for WIPP_SOURCE_SINK block 
  input => InputCreate(IN_UNIT,this%option%input_filename,this%option)
  block_string = 'WIPP_SOURCE_SINK'
  call InputFindStringInFile(input,this%option,block_string)
  if (input%ierr == 0 .and. wippflo_use_gas_generation) then
    this%pmwss_ptr => PMWSSCreate()
    this%pmwss_ptr%option => this%option
    call this%pmwss_ptr%Read(input)
  endif
  ! call setup/initialization of all WIPP process models
  if (associated(this%pmwss_ptr)) then
    call PMWSSSetRealization(this%pmwss_ptr,this%realization)
    call this%pmwss_ptr%Setup()
    call this%pmwss_ptr%InitializeRun()
  endif
  
end subroutine PMWIPPFloInitializeRun

! ************************************************************************** !

subroutine PMWIPPFloInitializeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use WIPP_Flow_module, only : WIPPFloInitializeTimestep
  use WIPP_Flow_Aux_module
  use Global_module
  use Variables_module, only : TORTUOSITY
  use Material_module, only : MaterialAuxVarCommunicate
  
  implicit none
  
  class(pm_wippflo_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)                                 
  if (this%option%print_screen_flag) then
    if (wippflo_use_bragflo_flux) then
      write(*,'(/,2("=")," BRAGFLO MODE ",64("="))')
    else
      write(*,'(/,2("=")," WIPP FLOW MODE ",62("="))')
    endif
  endif
  
  call WIPPFloInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)  
  
  ! initialize timestep of all WIPP process models
  if (associated(this%pmwss_ptr)) then
    call this%pmwss_ptr%InitializeTimestep()
  endif

  this%convergence_flags = 0
  this%convergence_reals = 0.d0
  wippflo_prev_liq_res_cell = 0
  wippflo_print_oscillatory_behavior = PETSC_FALSE
  
end subroutine PMWIPPFloInitializeTimestep

! ************************************************************************** !

subroutine PMWIPPFloFinalizeTimestep(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/21/17
  ! 
  implicit none
  
  class(pm_wippflo_type) :: this

  if (associated(this%pmwss_ptr)) then
    call this%pmwss_ptr%FinalizeTimestep()
  endif
  call PMSubsurfaceFlowFinalizeTimestep(this)

end subroutine PMWIPPFloFinalizeTimestep

! ************************************************************************** !

subroutine PMWIPPFloPreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17

  implicit none

  class(pm_wippflo_type) :: this

  call PMSubsurfaceFlowPreSolve(this)

end subroutine PMWIPPFloPreSolve

! ************************************************************************** !

subroutine PMWIPPFloPostSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  
  use WIPP_Flow_Common_module
  use WIPP_Flow_Aux_module
  use Option_module

  implicit none

  class(pm_wippflo_type) :: this

  PetscInt, save :: lr = 0, gr = 0, lbr = 0, gbr = 0
  PetscInt, save :: lj = 0, gj = 0, lbj = 0, gbj = 0

  if (wippflo_fix_upwind_direction .and. &
      wippflo_count_upwind_dir_flip .and. &
      OptionPrintToScreen(this%realization%option)) then
    write(*,'(6x,"Res: ",4i5," : ",4i7)') &
      liq_upwind_flip_count_by_res-lr, &
      gas_upwind_flip_count_by_res-gr, &
      liq_bc_upwind_flip_count_by_res-lbr, &
      gas_bc_upwind_flip_count_by_res-gbr, &
      liq_upwind_flip_count_by_res, &
      gas_upwind_flip_count_by_res, &
      liq_bc_upwind_flip_count_by_res, &
      gas_bc_upwind_flip_count_by_res
    write(*,'(6x,"Jac: ",4i5," : ",4i7)') &
      liq_upwind_flip_count_by_jac-lj, &
      gas_upwind_flip_count_by_jac-gj, &
      liq_bc_upwind_flip_count_by_jac-lbj, &
      gas_bc_upwind_flip_count_by_jac-gbj, &
      liq_upwind_flip_count_by_jac, &
      gas_upwind_flip_count_by_jac, &
      liq_bc_upwind_flip_count_by_jac, &
      gas_bc_upwind_flip_count_by_jac
  endif

  lr = liq_upwind_flip_count_by_res
  gr = gas_upwind_flip_count_by_res
  lbr = liq_bc_upwind_flip_count_by_res
  gbr = gas_bc_upwind_flip_count_by_res
  lj = liq_upwind_flip_count_by_jac
  gj = gas_upwind_flip_count_by_jac
  lbj = liq_bc_upwind_flip_count_by_jac
  gbj = gas_bc_upwind_flip_count_by_jac

end subroutine PMWIPPFloPostSolve

! ************************************************************************** !

subroutine PMWIPPFloUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                                   num_newton_iterations,tfac, &
                                   time_step_max_growth_factor)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Realization_Base_class, only : RealizationGetVariable
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL
  use Field_module
  use Global_module, only : GlobalSetAuxVarVecLoc
  use Variables_module, only : LIQUID_SATURATION, GAS_SATURATION
  use WIPP_Flow_Aux_module, only : wippflo_debug_ts_update

  implicit none
  
  class(pm_wippflo_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min ! DO NOT USE
  PetscReal :: dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor
  
  PetscReal :: dtime(2)
  type(field_type), pointer :: field

  PetscReal :: dt_prev

  dt_prev = dt
  
  ! calculate the time step ramping factor
  dtime(1) = (2.d0*this%satnorm)/(this%satnorm+this%max_saturation_change)
  dtime(2) = (2.d0*this%presnorm)/(this%presnorm+this%max_pressure_change)
  ! pick minimum time step from calc'd ramping factor or maximum ramping factor
  !TODO(geh) %dtimemax should be replace by time_step_max_growth_factor
  dt = min(min(dtime(1),dtime(2))*dt,this%dtimemax*dt)
  ! make sure time step is within bounds given in the input deck
  dt = min(dt,dt_max)
  ! do not use the PFLOTRAN dt_min as it will shut down the simulation from
  ! within timestepper_BE. use %deltmin, which is specific to bragflo.
  dt = max(dt,this%deltmin)

  if (wippflo_debug_ts_update) then
    if (minval(dtime(:)) < this%dtimemax .and. dt < dt_max) then
      write(*,'(" scaled dt: ",2es13.5)') dtime(:)
    endif
  endif

  if (Initialized(this%cfl_governor)) then
    ! Since saturations are not stored in global_auxvar for wipp flow mode, we
    ! must copy them over for the CFL check
    ! liquid saturation
    field => this%realization%field
    call RealizationGetVariable(this%realization,field%work, &
                                LIQUID_SATURATION,ZERO_INTEGER)
    call this%realization%comm1%GlobalToLocal(field%work,field%work_loc)
    call GlobalSetAuxVarVecLoc(this%realization,field%work_loc, &
                               LIQUID_SATURATION,TIME_NULL)
    call RealizationGetVariable(this%realization,field%work, &
                                GAS_SATURATION,ZERO_INTEGER)
    call this%realization%comm1%GlobalToLocal(field%work,field%work_loc)
    call GlobalSetAuxVarVecLoc(this%realization,field%work_loc, &
                               GAS_SATURATION,TIME_NULL)
    call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt)
  endif

end subroutine PMWIPPFloUpdateTimestep

! ************************************************************************** !

subroutine PMWIPPFloResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
  use WIPP_Flow_module, only : WIPPFloResidual
  use Debug_module

  implicit none
  
  class(pm_wippflo_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr

  PetscViewer :: viewer
  character(len=MAXSTRINGLENGTH) :: string
  
  call PMSubsurfaceFlowUpdatePropertiesNI(this)

  ! calculate residual
  call WIPPFloResidual(snes,xx,r,this%realization,this%pmwss_ptr,ierr)

  if (this%realization%debug%vecview_residual) then
    string = 'WFresidual'
    call DebugCreateViewer(this%realization%debug,string,this%option,viewer)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (this%realization%debug%vecview_solution) then
    string = 'WFxx'
    call DebugCreateViewer(this%realization%debug,string,this%option,viewer)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  call this%PostSolve()

end subroutine PMWIPPFloResidual

! ************************************************************************** !

subroutine PMWIPPFloJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use WIPP_Flow_module, only : WIPPFloJacobian
  use Debug_module
  use Option_module

  implicit none
  
  class(pm_wippflo_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr

  PetscViewer :: viewer
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: norm
  
  call WIPPFloJacobian(snes,xx,A,B,this%realization,this%pmwss_ptr,ierr)

  if (this%realization%debug%matview_Jacobian) then
    string = 'WFjacobian'
    call DebugCreateViewer(this%realization%debug,string,this%option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (this%realization%debug%norm_Jacobian) then
    call MatNorm(A,NORM_1,norm,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer,'("1 norm: ",es11.4)') norm
    call printMsg(this%option)
    call MatNorm(A,NORM_FROBENIUS,norm,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer,'("2 norm: ",es11.4)') norm
    call printMsg(this%option)
    call MatNorm(A,NORM_INFINITY,norm,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer,'("inf norm: ",es11.4)') norm
    call printMsg(this%option)
  endif

end subroutine PMWIPPFloJacobian

! ************************************************************************** !

subroutine PMWIPPFloCheckUpdatePre(this,line_search,X,dX,changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Option_module
  use Saturation_Function_module
  use Patch_module
  use WIPP_Flow_Aux_module
  use Global_Aux_module
  
  implicit none
  
  class(pm_wippflo_type) :: this
  SNESLineSearch :: line_search
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  PetscReal, pointer :: X_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: offset
  PetscInt :: saturation_index 
  PetscInt :: pressure_index
  PetscBool :: cut_timestep
  PetscBool :: force_another_iteration
  PetscBool :: outside_limits
  PetscReal :: saturation0, saturation1, del_saturation
  PetscReal :: pressure0, pressure1, del_pressure
  PetscReal :: max_gas_sat_outside_lim
  PetscInt :: max_gas_sat_outside_lim_cell

  SNES :: snes

  this%convergence_flags = 0
  this%convergence_reals = 0.d0
  changed = PETSC_FALSE
  
end subroutine PMWIPPFloCheckUpdatePre

! ************************************************************************** !

subroutine PMWIPPFloCheckUpdatePost(this,line_search,X0,dX,X1,dX_changed, &
                                    X1_changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  use Material_Aux_class  
  use WIPP_Flow_Aux_module
  
  implicit none
  
  class(pm_wippflo_type) :: this
  SNESLineSearch :: line_search
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch

  PetscReal, pointer :: X0_p(:)
  PetscReal, pointer :: X1_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscReal, pointer :: press_ptr(:)
  PetscReal, pointer :: sat_ptr(:)

  PetscInt :: local_id, ghosted_id
  PetscInt :: offset
  PetscInt :: saturation_index
  PetscInt :: pressure_index

  PetscReal :: abs_X
  PetscReal :: abs_dX
  PetscReal :: abs_dX_TS
  PetscReal :: abs_rel_dX_TS
  PetscBool :: converged_liquid_pressure
  PetscBool :: converged_gas_saturation
  PetscReal :: max_liq_pres_rel_change
  PetscReal :: max_gas_sat_change_NI
  PetscReal :: max_gas_sat_change_TS
  PetscReal :: max_abs_pressure_change_NI
  PetscReal :: max_abs_pressure_change_TS
  PetscReal :: max_rel_pressure_change_TS
  PetscReal :: min_liq_pressure
  PetscInt :: max_liq_pres_rel_change_cell
  PetscInt :: max_gas_sat_change_NI_cell
  PetscInt :: max_gas_sat_change_TS_cell
  PetscInt :: max_abs_pressure_change_NI_cell
  PetscInt :: max_abs_pressure_change_TS_cell
  PetscInt :: max_rel_pressure_change_TS_cell
  PetscInt :: min_liq_pressure_cell
  PetscReal :: abs_dX_over_absX

  PetscBool :: cut_timestep
  PetscBool :: force_another_iteration
  PetscReal :: pressure_outside_limits
  PetscReal :: saturation_outside_limits
  PetscReal :: max_gas_sat_outside_lim
  PetscInt :: max_gas_sat_outside_lim_cell
  PetscInt :: i
  ! dX_p is subtracted to update the solution.  The max values need to be 
  ! scaled by this delta_scale for proper screen output.
  PetscReal, parameter :: delta_scale = -1.d0
  
  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch

  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE
  
  call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  if (wippflo_print_update) then
    open(IUNIT_TEMP,file='pf_update.txt')
    do i = 1, grid%nlmax*2
      write(IUNIT_TEMP,'(1i5,es16.8)') i, dX_p(i)
    enddo
    close(IUNIT_TEMP)
  endif
  call VecGetArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)
  if (wippflo_print_solution) then
    open(IUNIT_TEMP,file='pf_solution.txt')
    do i = 1, grid%nlmax*2
      write(IUNIT_TEMP,'(1i5,es16.8)') i, X0_p(i)
    enddo
    close(IUNIT_TEMP)
  endif
  call VecGetArrayF90(X1,X1_p,ierr);CHKERRQ(ierr)
  ! max change variables: [LIQUID_PRESSURE, GAS_PRESSURE, GAS_SATURATION]
  call VecGetArrayReadF90(field%max_change_vecs(1),press_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%max_change_vecs(3),sat_ptr,ierr);CHKERRQ(ierr)
  converged_liquid_pressure = PETSC_TRUE
  converged_gas_saturation = PETSC_TRUE
  cut_timestep = PETSC_FALSE
  force_another_iteration = PETSC_FALSE
  max_liq_pres_rel_change = 0.d0
  max_gas_sat_change_NI = 0.d0
  max_gas_sat_change_TS = 0.d0
  max_abs_pressure_change_NI = 0.d0
  max_abs_pressure_change_TS = 0.d0
  max_rel_pressure_change_TS = 0.d0
  min_liq_pressure = 1.d20
  max_liq_pres_rel_change_cell = 0
  max_gas_sat_change_NI_cell = 0
  max_gas_sat_change_TS_cell = 0
  max_abs_pressure_change_NI_cell = 0
  max_abs_pressure_change_TS_cell = 0
  max_rel_pressure_change_TS_cell = 0
  min_liq_pressure_cell = 0
  max_gas_sat_outside_lim = 0.d0
  max_gas_sat_outside_lim_cell = 0
  do local_id = 1, grid%nlmax
    offset = (local_id-1)*option%nflowdof
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    pressure_index = offset + WIPPFLO_LIQUID_PRESSURE_DOF
    saturation_index = offset + WIPPFLO_GAS_SATURATION_DOF

    !NOTE: store the actual value, not the absolute value, better enabling the 
    !      pinpointing of oscillatory behavior.

    !TODO(geh): switch to flow_yy as it is cleaner and more precise.
    ! maximum relative change in liquid pressure
    abs_X = dabs(X1_p(pressure_index))
    if (abs_X > 0.d0) then
      abs_dX_over_absX = dabs(dX_p(pressure_index))/abs_X
      if (dabs(max_liq_pres_rel_change) < abs_dX_over_absX) then
        max_liq_pres_rel_change_cell = local_id
        max_liq_pres_rel_change = delta_scale*dX_p(pressure_index)/abs_X
      endif
      if (abs_dX_over_absX >= this%liquid_pressure_tolerance) then
        converged_liquid_pressure = PETSC_FALSE
      endif
    endif

    ! maximum absolute change in liquid pressure
    if (dabs(dX_p(pressure_index)) > dabs(max_abs_pressure_change_NI)) then
      max_abs_pressure_change_NI_cell = local_id
      max_abs_pressure_change_NI = delta_scale*dX_p(pressure_index)
    endif
    
    ! EPS_SAT maximum gas saturation change "digits of accuracy"
    abs_dX = dabs(dX_p(saturation_index))
    if (abs_dX > 0.d0) then
      if (dabs(max_gas_sat_change_NI) < abs_dX) then
        max_gas_sat_change_NI_cell = local_id
        max_gas_sat_change_NI = delta_scale*dX_p(saturation_index)
      endif
      !TODO(geh): change '<' to '<=' according to bragflo
      if ((-1.d0*log10(abs_dX)) < this%eps_sat) then
        converged_gas_saturation = PETSC_FALSE
      endif
    endif

    !TODO(geh): remove storage of signed max change over time step
    ! DSAT_MAX maximum absolute gas saturation change over time step
    abs_dX_TS = dabs(sat_ptr(local_id)-X1_p(saturation_index))
    if (abs_dX_TS > 0.d0) then
      if (dabs(max_gas_sat_change_TS) < abs_dX_TS) then
        max_gas_sat_change_TS_cell = local_id
        max_gas_sat_change_TS = delta_scale * & 
          (sat_ptr(local_id)-X1_p(saturation_index))
      endif
    endif
    
    ! DPRE_MAX maximum absolute liquid pressure change over time step
    abs_dX_TS = dabs(press_ptr(local_id)-X1_p(pressure_index))
    if (abs_dX_TS > 0.d0) then
      if (dabs(max_abs_pressure_change_TS) < abs_dX_TS) then
        max_abs_pressure_change_TS_cell = local_id
        max_abs_pressure_change_TS = delta_scale * &
          (press_ptr(local_id)-X1_p(pressure_index))
      endif
    endif
    
    ! EPS_PRES maximum relative liquid pressure change over time step
    !geh: BRAGFLO divides by DEPOUT(L), which is the updated solution (X1_p)
    abs_rel_dX_TS = dabs((press_ptr(local_id)-X1_p(pressure_index))/ &
                         X1_p(pressure_index))
    if (abs_rel_dX_TS > 0.d0) then
      if (dabs(max_rel_pressure_change_TS) < abs_rel_dX_TS) then
        max_rel_pressure_change_TS_cell = local_id
        max_rel_pressure_change_TS = delta_scale * &
          (press_ptr(local_id)-X1_p(pressure_index))/X1_p(pressure_index)
      endif
    endif

    ! liquid pressure
    if (X1_p(pressure_index) < min_liq_pressure) then
      min_liq_pressure = X1_p(pressure_index)
      min_liq_pressure_cell = local_id
    endif

    ! limits are checked after the calculations above since they can truncate
    pressure_outside_limits = 0.d0
    saturation_outside_limits = 0.d0
    if (X1_p(saturation_index) < 0.d0) then
      if (dabs(max_gas_sat_outside_lim) < dabs(X1_p(saturation_index))) then
        max_gas_sat_outside_lim = X1_p(saturation_index)
        max_gas_sat_outside_lim_cell = local_id
      endif
      if (X1_p(saturation_index) < (-1.d0*this%dsatlim)) then  ! DEPLIMIT(1)
        saturation_outside_limits = X1_p(saturation_index)
      else 
        if (X1_p(saturation_index) < (-1.d0*this%satlimit)) then  ! SATLIMIT
          force_another_iteration = PETSC_TRUE
        endif
        ! set saturation to zero
        X1_p(saturation_index) = 0.d0
        dX_p(saturation_index) = X0_p(saturation_index)
        dX_changed = PETSC_TRUE
        X1_changed = PETSC_TRUE
      endif
    else if (X1_p(saturation_index) > 1.d0) then
      if (abs(max_gas_sat_outside_lim) < X1_p(saturation_index) - 1.d0) then
        max_gas_sat_outside_lim = X1_p(saturation_index) - 1.d0
        max_gas_sat_outside_lim_cell = local_id
      endif
      if (X1_p(saturation_index) > 1.d0 + this%dsatlim) then  ! DEPLIMIT(1)
        saturation_outside_limits = X1_p(saturation_index)
      else 
        if (X1_p(saturation_index) > 1.d0 + this%satlimit) then  ! SATLIMIT
          force_another_iteration = PETSC_TRUE
        endif
        ! set saturation to one
        X1_p(saturation_index) = 1.d0
        dX_p(saturation_index) = X0_p(saturation_index) - 1.d0
        dX_changed = PETSC_TRUE
        X1_changed = PETSC_TRUE
      endif
    endif
    ! DPRELIM is designed to catch large negative values in liquid pressure
    ! and cut the timestep if this occurs
    if (X1_p(pressure_index) <= (this%dprelim)) then  ! DEPLIMIT(2)
      pressure_outside_limits = X1_p(pressure_index)
    endif

    if (dabs(pressure_outside_limits) > 0.d0 .or. &
        dabs(saturation_outside_limits) > 0.d0) then
      cut_timestep = PETSC_TRUE
      write(*,'(4x,"Outside Limits (PL,SG): ",i8,2es10.2)') local_id, &
        pressure_outside_limits, saturation_outside_limits
    endif

  enddo

  if (wippflo_debug_first_iteration) stop

  if (.not.converged_liquid_pressure) then 
    this%convergence_flags(MAX_REL_CHANGE_LIQ_PRES_NI) = &
      max_liq_pres_rel_change_cell
  endif
  if (.not.converged_gas_saturation) then 
    this%convergence_flags(MAX_CHANGE_GAS_SAT_NI) = max_gas_sat_change_NI_cell
  endif
  ! this following flags can always be set
  this%convergence_flags(MAX_CHANGE_GAS_SAT_NI_TRACK) = &
    max_gas_sat_change_NI_cell
  this%convergence_flags(MAX_CHANGE_LIQ_PRES_NI) = &
    max_abs_pressure_change_NI_cell
  this%convergence_flags(MIN_LIQ_PRES) = min_liq_pressure_cell
  this%convergence_reals(MAX_REL_CHANGE_LIQ_PRES_NI) = max_liq_pres_rel_change
  this%convergence_reals(MAX_REL_CHANGE_LIQ_PRES_TS) = &
    max_rel_pressure_change_TS
  this%convergence_reals(MAX_CHANGE_LIQ_PRES_NI) = max_abs_pressure_change_NI
  this%convergence_reals(MAX_CHANGE_LIQ_PRES_TS) = max_abs_pressure_change_TS
  this%convergence_reals(MAX_CHANGE_GAS_SAT_NI) = max_gas_sat_change_NI
  this%convergence_reals(MAX_CHANGE_GAS_SAT_NI_TRACK) = max_gas_sat_change_NI
  this%convergence_reals(MAX_CHANGE_GAS_SAT_TS) = max_gas_sat_change_TS
  this%convergence_reals(MIN_LIQ_PRES) = min_liq_pressure
  this%convergence_reals(FORCE_ITERATION) = max_gas_sat_outside_lim
  if (force_another_iteration) this%convergence_flags(FORCE_ITERATION) = &
    max_gas_sat_outside_lim_cell
  if (cut_timestep) this%convergence_flags(OUTSIDE_BOUNDS) = 1

  call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(X1,X1_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%max_change_vecs(1),press_ptr, &
                              ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%max_change_vecs(3),sat_ptr, &
                              ierr);CHKERRQ(ierr)
                               
end subroutine PMWIPPFloCheckUpdatePost

! ************************************************************************** !

subroutine PMWIPPFloConvergence(this,snes,it,xnorm,unorm, &
                                            fnorm,reason,ierr)
  ! Author: Glenn Hammond
  ! Date: 11/15/17
  ! 
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  use Material_Aux_class  
  use WIPP_Flow_Aux_module
  use Convergence_module

  implicit none

  class(pm_wippflo_type) :: this
  SNES :: snes
  PetscInt :: it
  PetscReal :: xnorm
  PetscReal :: unorm
  PetscReal :: fnorm
  SNESConvergedReason :: reason
  PetscErrorCode :: ierr

  Vec :: residual_vec
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum2_p(:)
  PetscReal, pointer :: X1_p(:)
  PetscReal :: residual
  PetscReal :: accumulation
  PetscReal :: abs_residual_over_accumulation
  character(len=10) :: reason_string

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  class(material_auxvar_type), pointer :: material_auxvars(:)  

  PetscInt :: local_id, ghosted_id
  PetscInt :: offset
  PetscInt :: liquid_equation_index
  PetscInt :: gas_equation_index
  PetscInt :: converged_flag

  PetscReal, parameter :: zero_saturation = 1.d-15
  PetscReal, parameter :: zero_accumulation = 1.d-15

  PetscBool :: converged_liquid_equation
  PetscBool :: converged_gas_equation
  PetscReal :: max_res_liq_
  PetscReal :: max_res_gas_
  PetscReal :: max_normal_res_liq_
  PetscReal :: max_normal_res_gas_
  PetscReal :: min_gas_pressure
  PetscInt :: max_res_liq_cell
  PetscInt :: max_res_gas_cell
  PetscInt :: max_normal_res_liq_cell
  PetscInt :: max_normal_res_gas_cell
  PetscInt :: min_gas_pressure_cell
  PetscReal :: pflotran_to_bragflo(2)
  PetscReal :: bragflo_residual(2)
  PetscReal :: bragflo_accum(2)
  PetscInt :: istart, iend
  PetscReal :: tempreal3
  PetscInt :: tempint3
  PetscReal :: tempreal4
  PetscInt :: tempint4
  PetscInt :: i
  PetscMPIInt :: int_mpi
  PetscBool :: cell_id_match
  
  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch
  wippflo_auxvars => patch%aux%WIPPFlo%auxvars
  material_auxvars => patch%aux%Material%auxvars

  residual_vec = tVec(0)
  ! check residual terms
  if (this%stored_residual_vec == PETSC_NULL_VEC) then
    residual_vec = field%flow_r
  else
    ! this vector is to be used if linear system scaling is employed in 
    ! BRAGFLO mode as the residual will be altered by that scaling.  This
    ! vector stores the original residual. 
    residual_vec = this%stored_residual_vec
  endif
  call VecGetArrayReadF90(residual_vec,r_p,ierr);CHKERRQ(ierr)
  if (wippflo_print_residual) then
    open(IUNIT_TEMP,file='pf_residual.txt')
    do i = 1, grid%nlmax*2
      write(IUNIT_TEMP,'(1i5,es16.8)') i, r_p(i)
    enddo
    close(IUNIT_TEMP)
  endif
  call VecGetArrayReadF90(field%flow_accum2,accum2_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%flow_xx,X1_p,ierr);CHKERRQ(ierr)
  converged_liquid_equation = PETSC_TRUE
  converged_gas_equation = PETSC_TRUE
  max_normal_res_liq_ = 0.d0
  max_normal_res_gas_ = 0.d0
  min_gas_pressure = 1.d20
  max_normal_res_liq_cell = 0
  max_normal_res_gas_cell = 0
  min_gas_pressure_cell = 0
  do local_id = 1, grid%nlmax
    offset = (local_id-1)*option%nflowdof
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    liquid_equation_index = offset + WIPPFLO_LIQUID_EQUATION_INDEX
    gas_equation_index = offset + WIPPFLO_GAS_EQUATION_INDEX

    
    bragflo_residual = r_p(liquid_equation_index:gas_equation_index)
    bragflo_accum = accum2_p(liquid_equation_index:gas_equation_index)
    if (.not.wippflo_use_bragflo_units) then
      pflotran_to_bragflo = fmw_comp * option%flow_dt / &
                            material_auxvars(ghosted_id)%volume
      bragflo_residual = bragflo_residual * pflotran_to_bragflo
      bragflo_accum = bragflo_accum * pflotran_to_bragflo
    else
      
    endif

    if (wippflo_debug) then
      ! in bragflo gas is first.
      print *, local_id, bragflo_residual(2), bragflo_residual(1)
    endif
    
    ! liquid component equation
    residual = bragflo_residual(WIPPFLO_LIQUID_EQUATION_INDEX)
    accumulation = bragflo_accum(WIPPFLO_LIQUID_EQUATION_INDEX)

    ! residual
    if (dabs(residual) > dabs(max_res_liq_)) then
      max_res_liq_cell = local_id
      max_res_liq_ = residual
    endif

    ! normalized residual
    if (dabs(accumulation) > zero_accumulation) then 
      abs_residual_over_accumulation = dabs(residual / accumulation)
      if (dabs(residual) > this%liquid_equation_tolerance) then
        if (abs_residual_over_accumulation > &
            this%liquid_equation_tolerance) then
          converged_liquid_equation = PETSC_FALSE
          if (dabs(max_normal_res_liq_) < dabs(residual)) then
            max_normal_res_liq_cell = local_id
          endif
        endif
      endif
      ! update outside to always record the maximum residual
      if (dabs(max_normal_res_liq_) < dabs(residual)) then
        if (wippflo_match_bragflo_output) then
          if (dabs(residual) > this%liquid_equation_tolerance) then
            max_normal_res_liq_ = abs_residual_over_accumulation
          else
            max_normal_res_liq_ = residual
          endif
        else
          max_normal_res_liq_ = residual
        endif
      endif
    endif

    ! gas component equation
    residual = bragflo_residual(WIPPFLO_GAS_EQUATION_INDEX)
    accumulation = bragflo_accum(WIPPFLO_GAS_EQUATION_INDEX)

    ! residual
    if (dabs(residual) > dabs(max_res_gas_)) then
      max_res_gas_cell = local_id
      max_res_gas_ = residual
    endif

    ! normalized residual
    if (dabs(accumulation) > zero_accumulation .and. &
        X1_p(gas_equation_index) > zero_saturation) then
      abs_residual_over_accumulation = abs(residual / accumulation)
      if (dabs(residual) > this%gas_equation_tolerance) then 
        if (abs_residual_over_accumulation > &
            this%gas_equation_tolerance) then
          converged_gas_equation = PETSC_FALSE
          if (dabs(max_normal_res_gas_) < dabs(residual)) then
            max_normal_res_gas_cell = local_id
          endif
        endif
      endif
      ! update outside to always record the maximum residual
      if (dabs(max_normal_res_gas_) < dabs(residual)) then
        if (wippflo_match_bragflo_output) then
          if (dabs(residual) > this%gas_equation_tolerance) then
            max_normal_res_gas_ = abs_residual_over_accumulation
          else
            max_normal_res_gas_ = residual
          endif
        else
          max_normal_res_gas_ = residual
        endif
      endif
    endif

    ! gas pressure
    if (wippflo_auxvars(0,ghosted_id)%pres(2) < min_gas_pressure) then
      min_gas_pressure = wippflo_auxvars(0,ghosted_id)%pres(2)
      min_gas_pressure_cell = local_id
    endif
  enddo

  if (.not.converged_liquid_equation) then
    this%convergence_flags(MAX_NORMAL_RES_LIQ) = max_normal_res_liq_cell
  endif
  if (.not.converged_gas_equation) then
    this%convergence_flags(MAX_NORMAL_RES_GAS) = max_normal_res_gas_cell
  endif
  if (min_gas_pressure < 0.d0 .and. .not.wippflo_allow_neg_gas_pressure) then
    this%convergence_flags(MIN_GAS_PRES) = min_gas_pressure_cell
  endif
  ! the following flags are not used for convergence purposes, and thus can
  ! always be set
  this%convergence_flags(MAX_RES_LIQ) = max_res_liq_cell
  this%convergence_flags(MAX_RES_GAS) = max_res_gas_cell

  this%convergence_reals(MAX_RES_LIQ) = max_res_liq_
  this%convergence_reals(MAX_RES_GAS) = max_res_gas_
  this%convergence_reals(MAX_NORMAL_RES_LIQ) = max_normal_res_liq_
  this%convergence_reals(MAX_NORMAL_RES_GAS) = max_normal_res_gas_
  this%convergence_reals(MIN_GAS_PRES) = min_gas_pressure

  int_mpi = size(this%convergence_flags)
  call MPI_Allreduce(MPI_IN_PLACE,this%convergence_flags,int_mpi, &
                     MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
  ! if running in parallel, we can no longer report the sign on the maximum
  ! change variables as the sign may differ across processes.
  if (option%mycommsize > 1) then
    this%convergence_reals(1:MIN_LIQ_PRES-1) = &
      dabs(this%convergence_reals(1:MIN_LIQ_PRES-1))
  endif
  ! flip sign on min pressure in order to use MPI_MAX to get minimum value
  this%convergence_reals(MIN_LIQ_PRES) = &
    -1.d0 * this%convergence_reals(MIN_LIQ_PRES)
  this%convergence_reals(MIN_GAS_PRES) = &
    -1.d0 * this%convergence_reals(MIN_GAS_PRES)
  int_mpi = size(this%convergence_reals)
  call MPI_Allreduce(MPI_IN_PLACE,this%convergence_reals,int_mpi, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
  ! flip sign back
  this%convergence_reals(MIN_LIQ_PRES) = &
    -1.d0 * this%convergence_reals(MIN_LIQ_PRES)
  this%convergence_reals(MIN_GAS_PRES) = &
    -1.d0 * this%convergence_reals(MIN_GAS_PRES)

  ! to catch oscillatory behavior
  if (wippflo_check_oscillatory_behavior) then
    do i = size(wippflo_prev_liq_res_cell), 2, -1
      wippflo_prev_liq_res_cell(i) = wippflo_prev_liq_res_cell(i-1)
    enddo
    wippflo_prev_liq_res_cell(1) = this%convergence_flags(MAX_NORMAL_RES_LIQ)
    if (wippflo_prev_liq_res_cell(size(wippflo_prev_liq_res_cell)) > 0) then
      cell_id_match = PETSC_TRUE
      do i = 1, size(wippflo_prev_liq_res_cell)-2
        if (wippflo_prev_liq_res_cell(i) /= wippflo_prev_liq_res_cell(i+2)) then
          cell_id_match = PETSC_FALSE
        endif
      enddo
      if (cell_id_match) then
        wippflo_print_oscillatory_behavior = PETSC_TRUE
      endif
    endif
  endif

  ! these conditionals cannot change order
  reason_string = '-------|--'
  converged_flag = CONVERGENCE_CONVERGED
  if (this%convergence_flags(OUTSIDE_BOUNDS) > 0) then
    converged_flag = CONVERGENCE_CUT_TIMESTEP
    reason_string(1:1) = 'B'
  endif
  if (this%convergence_flags(MIN_GAS_PRES) > 0) then
    converged_flag = CONVERGENCE_CUT_TIMESTEP
    reason_string(2:2) = 'N'
  endif
  if (this%convergence_flags(FORCE_ITERATION) > 0) then
    if (converged_flag == CONVERGENCE_CONVERGED) then
      converged_flag = CONVERGENCE_FORCE_ITERATION
    endif
    reason_string(3:3) = '!'
  endif
  converged_liquid_equation = PETSC_TRUE
  converged_gas_equation = PETSC_TRUE
  if (this%convergence_flags(MAX_NORMAL_RES_LIQ) > 0) then
    if (converged_flag /= CONVERGENCE_CUT_TIMESTEP) then
      converged_flag = CONVERGENCE_KEEP_ITERATING ! cannot override cut
    endif
    converged_liquid_equation = PETSC_FALSE
    reason_string(4:4) = 'L'
  endif
  if (this%convergence_flags(MAX_NORMAL_RES_GAS) > 0) then
    if (converged_flag /= CONVERGENCE_CUT_TIMESTEP) then
      converged_flag = CONVERGENCE_KEEP_ITERATING ! cannot override cut
    endif
    converged_gas_equation = PETSC_FALSE
    reason_string(5:5) = 'G'
  endif
  if (this%convergence_flags(MAX_REL_CHANGE_LIQ_PRES_NI) > 0) then
    if (this%iconvtest == 1 .or. .not.converged_liquid_equation) then
      if (converged_flag /= CONVERGENCE_CUT_TIMESTEP) then
        converged_flag = CONVERGENCE_KEEP_ITERATING ! cannot override cut
      endif
      reason_string(6:6) = 'P'
    endif
  endif
  if (this%convergence_flags(MAX_CHANGE_GAS_SAT_NI) > 0) then
    if (converged_flag /= CONVERGENCE_CUT_TIMESTEP) then
      if (this%iconvtest == 1 .or. .not.converged_gas_equation) then
        converged_flag = CONVERGENCE_KEEP_ITERATING ! cannot override cut
      endif
      reason_string(7:7) = 'S'
    endif
  endif
  if (converged_flag == CONVERGENCE_CONVERGED) then
    ! converged based on NI criteria, but need to check TS
    if (this%convergence_flags(MAX_CHANGE_LIQ_PRES_TS) > 0) then
      converged_flag = CONVERGENCE_CUT_TIMESTEP
      reason_string(9:9) = 'P'
    endif
    if (this%convergence_flags(MAX_CHANGE_GAS_SAT_TS) > 0) then
      converged_flag = CONVERGENCE_CUT_TIMESTEP
      reason_string(10:10) = 'S'
    endif
  endif
  if (OptionPrintToScreen(option)) then
    !TODO(geh): add the option to report only violated tolerances, zeroing
    !           the others.
    tempreal3 = this%convergence_reals(MAX_REL_CHANGE_LIQ_PRES_NI)
    tempint3 = this%convergence_flags(MAX_REL_CHANGE_LIQ_PRES_NI)
    tempreal4 = this%convergence_reals(MAX_CHANGE_GAS_SAT_NI)
    tempint4 = this%convergence_flags(MAX_CHANGE_GAS_SAT_NI)
    if (this%convergence_flags(MIN_GAS_PRES) > 0) then
      tempreal3 = this%convergence_reals(MIN_GAS_PRES)
      tempint3 = this%convergence_flags(MIN_GAS_PRES)
      reason_string(6:6) = 'N'
    endif
    if (this%convergence_flags(FORCE_ITERATION) > 0) then
      tempreal4 = this%convergence_reals(FORCE_ITERATION)
      tempint4 = this%convergence_flags(FORCE_ITERATION)
      reason_string(7:7) = '!'
    endif
    if (this%convergence_flags(OUTSIDE_BOUNDS) > 0) then
      ! just overwrite the character, the flag/real matches FORCE_ITERATION
      reason_string(7:7) = 'B'
    endif
    if (option%mycommsize > 1) then
      write(*,'(4x,"Rsn: ",a10,4es10.2)') reason_string, &
        this%convergence_reals(MAX_NORMAL_RES_LIQ), &
        this%convergence_reals(MAX_NORMAL_RES_GAS), &
        tempreal3, tempreal4
    else
      write(*,'(4x,"Rsn: ",a10,4(i5,es10.2))') reason_string, &
        this%convergence_flags(MAX_NORMAL_RES_LIQ), &
        this%convergence_reals(MAX_NORMAL_RES_LIQ), &
        this%convergence_flags(MAX_NORMAL_RES_GAS), &
        this%convergence_reals(MAX_NORMAL_RES_GAS), &
        tempint3, tempreal3, &
        tempint4, tempreal4
! for debugging minumum gas pressure cell
!      local_id = min_gas_pressure_cell
!      offset = (local_id-1)*option%nflowdof
!      istart = offset + 1
!      iend = offset + option%nflowdof
!      ghosted_id = grid%nL2G(local_id)
!      write(*,'(4x,i5,3es11.3)') local_id, &
!        wippflo_auxvars(0,ghosted_id)%pres(1:2), &
!        wippflo_auxvars(0,ghosted_id)%pres(option%capillary_pressure_id)
! for monitoring block
      if (wippflo_residual_test_cell > 0) then
        local_id = wippflo_residual_test_cell
        offset = (local_id-1)*option%nflowdof
        istart = offset + 1
        iend = offset + option%nflowdof
        ghosted_id = grid%nL2G(local_id)
        write(*,'(2x,"GEH: ",i5,3es12.4,2es12.4)') local_id, &
          wippflo_auxvars(0,ghosted_id)%pres(1:2), &
          wippflo_auxvars(0,ghosted_id)%pres(option%capillary_pressure_id), &
          wippflo_auxvars(0,ghosted_id)%sat(:)
      endif
    endif
    if (wippflo_match_bragflo_output) then
      write(*,'(x,"GEHMAX(SPGL): ",4(i4,es11.3))') &
        this%convergence_flags(MAX_CHANGE_GAS_SAT_NI_TRACK), &
        this%convergence_reals(MAX_CHANGE_GAS_SAT_NI_TRACK), &
        this%convergence_flags(MAX_CHANGE_LIQ_PRES_NI), &
        this%convergence_reals(MAX_CHANGE_LIQ_PRES_NI), &
        this%convergence_flags(MAX_RES_GAS), &
        this%convergence_reals(MAX_RES_GAS), &
        this%convergence_flags(MAX_RES_LIQ), &
        this%convergence_reals(MAX_RES_LIQ)
    endif
  endif
  option%convergence = converged_flag

  call VecRestoreArrayReadF90(residual_vec,r_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%flow_accum2,accum2_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%flow_xx,X1_p,ierr);CHKERRQ(ierr)

  call ConvergenceTest(snes,it,xnorm,unorm,fnorm,reason, &
                       this%realization%patch%grid, &
                       this%option,this%solver,ierr)

end subroutine PMWIPPFloConvergence

! ************************************************************************** !

subroutine PMWIPPFloTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use WIPP_Flow_module, only : WIPPFloTimeCut
  use WIPP_Flow_Aux_module, only : wippflo_prev_liq_res_cell, &
                                   wippflo_print_oscillatory_behavior, &
                                   wippflo_match_bragflo_output

  implicit none
  
  class(pm_wippflo_type) :: this

  if (wippflo_match_bragflo_output) then
    write(*,'(5x,"SPGL: ",4(i5,es11.3))') &
      this%convergence_flags(MAX_CHANGE_GAS_SAT_NI), &
      dabs(this%convergence_reals(MAX_CHANGE_GAS_SAT_NI)*10.d0**this%eps_sat), &
      this%convergence_flags(MAX_REL_CHANGE_LIQ_PRES_NI), &
      dabs(this%convergence_reals(MAX_REL_CHANGE_LIQ_PRES_NI)/ &
           this%liquid_pressure_tolerance), &
      this%convergence_flags(MAX_NORMAL_RES_GAS), &
      dabs(this%convergence_reals(MAX_NORMAL_RES_GAS)/ &
           this%gas_equation_tolerance), &
      this%convergence_flags(MAX_NORMAL_RES_LIQ), &
      dabs(this%convergence_reals(MAX_NORMAL_RES_LIQ)/ &
           this%liquid_equation_tolerance)
  endif
  
  call PMSubsurfaceFlowTimeCut(this)
  call WIPPFloTimeCut(this%realization)

  this%convergence_flags = 0
  this%convergence_reals = 0.d0
  wippflo_prev_liq_res_cell = 0
  wippflo_print_oscillatory_behavior = PETSC_FALSE

end subroutine PMWIPPFloTimeCut

! ************************************************************************** !

subroutine PMWIPPFloUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use WIPP_Flow_module, only : WIPPFloUpdateSolution, &
                               WIPPFloMapBCAuxVarsToGlobal
  use WIPP_Flow_Aux_module, only : wippflo_debug

  implicit none
  
  class(pm_wippflo_type) :: this

  if (wippflo_debug) then
    write(*,'("Final: ",10es14.6)') &
      this%realization%patch%aux%WIPPFlo%auxvars(0,1)%pres(1:2), &
      this%realization%patch%aux%WIPPFlo%auxvars(0,1)%effective_porosity, &
      this%realization%patch%aux%WIPPFlo%auxvars(0,1)%sat(1)
  endif
  
  call PMSubsurfaceFlowUpdateSolution(this)
  call WIPPFloUpdateSolution(this%realization)
  call WIPPFloMapBCAuxVarsToGlobal(this%realization)

end subroutine PMWIPPFloUpdateSolution     

! ************************************************************************** !

subroutine PMWIPPFloUpdateAuxVars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  use WIPP_Flow_module, only : WIPPFloUpdateAuxVars

  implicit none
  
  class(pm_wippflo_type) :: this

  call WIPPFloUpdateAuxVars(this%realization)

end subroutine PMWIPPFloUpdateAuxVars   

! ************************************************************************** !

subroutine PMWIPPFloMaxChange(this)
  ! 
  ! Not needed given WIPPFloMaxChange is called in PostSolve
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Realization_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use WIPP_Flow_Aux_module
  use Variables_module, only : LIQUID_PRESSURE, LIQUID_MOLE_FRACTION, &
                               TEMPERATURE, GAS_PRESSURE, AIR_PRESSURE, &
                               GAS_SATURATION
  implicit none
  
  class(pm_wippflo_type) :: this
  
  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_old_ptr(:), vec_new_ptr(:)
  PetscReal :: max_change_local(6)
  PetscReal :: max_change_global(6)
  PetscReal :: max_change, change
  PetscInt :: i, j
  
  PetscErrorCode :: ierr

  realization => this%realization
  option => realization%option
  field => realization%field
  grid => realization%patch%grid

  max_change_global = 0.d0
  max_change_local = 0.d0
  
  ! max change variables: [LIQUID_PRESSURE, GAS_PRESSURE, GAS_SATURATION]
  ! these are values from the previous time step
  do i = 1, 3
    call RealizationGetVariable(realization,field%work, &
                                this%max_change_ivar(i),ZERO_INTEGER)
    ! yes, we could use VecWAXPY and a norm here, but we need the ability
    ! to customize
    call VecGetArrayF90(field%work,vec_new_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%max_change_vecs(i),vec_old_ptr,ierr);CHKERRQ(ierr)
    max_change = 0.d0
    do j = 1, grid%nlmax
      ! have to weed out cells that changed state
      if (dabs(vec_old_ptr(j)) > 1.d-40 .and. &
          dabs(vec_new_ptr(j)) > 1.d-40) then
        ! this is an absolute change over a time step
        change = dabs(vec_new_ptr(j)-vec_old_ptr(j))
        if (i == 3) then ! (gas saturation)
          if (dabs(vec_new_ptr(j)) > this%tswitch) then
            ! use relative change only if gas saturation is higher than TSWITCH
            change = dabs(change/vec_new_ptr(j))
          endif
        endif
        max_change = max(max_change,change)  
      endif
    enddo
    max_change_local(i) = max_change
    call VecRestoreArrayF90(field%work,vec_new_ptr,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%max_change_vecs(i),vec_old_ptr, &
                            ierr);CHKERRQ(ierr)
    call VecCopy(field%work,field%max_change_vecs(i),ierr);CHKERRQ(ierr)
  enddo
  call MPI_Allreduce(max_change_local,max_change_global,SIX_INTEGER, &
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
  ! print them out
  if (OptionPrintToScreen(option)) then
    write(*,'("  --> max chng: dpl= ",1pe12.4, " dpg= ",1pe12.4, &
             &" dsg= ",1pe12.4)') &
      max_change_global(1:3)
  endif
  if (OptionPrintToFile(option)) then
    write(option%fid_out,'("  --> max chng: dpl= ",1pe12.4, " dpg= ", &
                          &1pe12.4, " dsg= ",1pe12.4)') &
      max_change_global(1:3)
  endif

  ! max change variables: [LIQUID_PRESSURE, GAS_PRESSURE, GAS_SATURATION]
  this%max_pressure_change = max_change_global(1)
  this%max_saturation_change = max_change_global(3)
  
end subroutine PMWIPPFloMaxChange

! ************************************************************************** !

subroutine PMWIPPFloComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use WIPP_Flow_module, only : WIPPFloComputeMassBalance

  implicit none
  
  class(pm_wippflo_type) :: this
  PetscReal :: mass_balance_array(:)
  
  call WIPPFloComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMWIPPFloComputeMassBalance

! ************************************************************************** !

subroutine PMWIPPFloInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 07/11/17
  ! 
  
  implicit none
  
  class(pm_wippflo_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'wipp flow'

end subroutine PMWIPPFloInputRecord

! ************************************************************************** !

subroutine PMWIPPFloCheckpointBinary(this,viewer)
  ! 
  ! Checkpoints data associated with WIPPFlo PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17

  use Checkpoint_module
  use Global_module
  use Variables_module, only : STATE

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_wippflo_type) :: this
  PetscViewer :: viewer
  
  call PMSubsurfaceFlowCheckpointBinary(this,viewer)
  
end subroutine PMWIPPFloCheckpointBinary

! ************************************************************************** !

subroutine PMWIPPFloRestartBinary(this,viewer)
  ! 
  ! Restarts data associated with WIPPFlo PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17

  use Checkpoint_module
  use Global_module
  use Variables_module, only : STATE

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_wippflo_type) :: this
  PetscViewer :: viewer
  
  call PMSubsurfaceFlowRestartBinary(this,viewer)
  
end subroutine PMWIPPFloRestartBinary

! ************************************************************************** !

subroutine PMWIPPFloDestroy(this)
  ! 
  ! Destroys WIPPFlo process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use WIPP_Flow_module, only : WIPPFloDestroy

  implicit none
  
  class(pm_wippflo_type) :: this

  PetscErrorCode :: ierr
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  deallocate(this%max_change_ivar)
  nullify(this%max_change_ivar)
  if (this%stored_residual_vec /= PETSC_NULL_VEC) then
    call VecDestroy(this%stored_residual_vec,ierr);CHKERRQ(ierr)
  endif

  ! preserve this ordering
  if (associated(this%pmwss_ptr)) then
    call this%pmwss_ptr%Destroy()
    deallocate(this%pmwss_ptr)
    nullify(this%pmwss_ptr)
  endif
  call WIPPFloDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)
  
end subroutine PMWIPPFloDestroy
  
end module PM_WIPP_Flow_class
