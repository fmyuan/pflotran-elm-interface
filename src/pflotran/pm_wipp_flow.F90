module PM_WIPP_Flow_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class
  use PM_WIPP_SrcSink_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

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
    class(pm_wipp_srcsink_type), pointer :: pmwss_ptr
    Vec :: stored_residual_vec
  contains
    procedure, public :: Read => PMWIPPFloRead
    procedure, public :: InitializeRun => PMWIPPFloInitializeRun
    procedure, public :: InitializeTimestep => PMWIPPFloInitializeTimestep
    procedure, public :: Residual => PMWIPPFloResidual
    procedure, public :: Jacobian => PMWIPPFloJacobian
    procedure, public :: UpdateTimestep => PMWIPPFloUpdateTimestep
    procedure, public :: PreSolve => PMWIPPFloPreSolve
    procedure, public :: PostSolve => PMWIPPFloPostSolve
    procedure, public :: CheckUpdatePre => PMWIPPFloCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMWIPPFloCheckUpdatePost
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
  this%stored_residual_vec = PETSC_NULL_VEC

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
    case('DEBUG_FIRST_ITERATION')
      wippflo_debug_first_iteration = PETSC_TRUE
    case('USE_LEGACY_PERTURBATION')
      wippflo_use_legacy_perturbation = PETSC_TRUE
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
  if (input%ierr == 0) then
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
  ! Should not need this as it is called in PreSolve.
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
  
end subroutine PMWIPPFloInitializeTimestep

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

  implicit none
  
  class(pm_wippflo_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor
  
  PetscReal :: dtime(2)
  type(field_type), pointer :: field
  
  ! calculate the time step ramping factor
  dtime(1) = (2.d0*this%satnorm)/(this%satnorm+this%max_saturation_change)
  dtime(2) = (2.d0*this%presnorm)/(this%presnorm+this%max_pressure_change)
  ! pick minimum time step from calc'd ramping factor or maximum ramping factor
  dt = min(min(dtime(1),dtime(2))*dt,this%dtimemax*dt, &
           time_step_max_growth_factor*dt)
  ! make sure time step is within bounds given in the input deck
  dt = min(dt,dt_max)
  dt = max(dt,dt_min)

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

  implicit none
  
  class(pm_wippflo_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  call PMSubsurfaceFlowUpdatePropertiesNI(this)

  ! calculate residual
  call WIPPFloResidual(snes,xx,r,this%realization,this%pmwss_ptr,ierr)

  call this%PostSolve()

end subroutine PMWIPPFloResidual

! ************************************************************************** !

subroutine PMWIPPFloJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use WIPP_Flow_module, only : WIPPFloJacobian

  implicit none
  
  class(pm_wippflo_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  call WIPPFloJacobian(snes,xx,A,B,this%realization,this%pmwss_ptr,ierr)

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
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  PetscReal, pointer :: X_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: offset
  PetscInt :: saturation_index 
  PetscInt :: pressure_index
  PetscInt :: temp_int(2)
  PetscBool :: cut_timestep
  PetscBool :: force_another_iteration
  PetscReal :: saturation0, saturation1, del_saturation
  PetscReal :: pressure0, pressure1, del_pressure

  SNES :: snes
  
  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  wippflo_auxvars => this%realization%patch%aux%WIPPFlo%auxvars
  global_auxvars => this%realization%patch%aux%Global%auxvars

  patch => this%realization%patch

  call SNESLineSearchGetSNES(line_search,snes,ierr)

  call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

  changed = PETSC_TRUE
  cut_timestep = PETSC_FALSE
  force_another_iteration = PETSC_FALSE

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    offset = (local_id-1)*option%nflowdof
    saturation_index = offset + WIPPFLO_GAS_SATURATION_DOF
    pressure_index = offset + WIPPFLO_LIQUID_PRESSURE_DOF
    del_saturation = dX_p(saturation_index)
    del_pressure = dX_p(pressure_index)
    saturation0 = X_p(saturation_index)
    pressure0 = X_p(pressure_index)
    saturation1 = saturation0 - del_saturation
    pressure1 = pressure0 - del_pressure
    ! DSATLIM is designed to catch saturations well outside the range of the 
    ! physically realistic range (0 - 1) and cut the time step if this occurs
    ! SATLIMIT is designed to force another newton iteration if the gas 
    ! saturation is "slightly" outside of the range of realistic values
    if (saturation1 < 0.d0) then
      if (saturation1 < (-1.d0*this%dsatlim)) then  ! DEPLIMIT(1)
        cut_timestep = PETSC_TRUE
      else 
        if (saturation1 < (-1.d0*this%satlimit)) then  ! SATLIMIT
          force_another_iteration = PETSC_TRUE
        endif
        ! set saturation to zero
        dX_p(saturation_index) = saturation0
      endif
    else if (saturation1 > 1.d0) then
      if (saturation1 > 1.d0 + this%dsatlim) then  ! DEPLIMIT(1)
        cut_timestep = PETSC_TRUE
      else 
        if (saturation1 > 1.d0 + this%satlimit) then  ! SATLIMIT
          force_another_iteration = PETSC_TRUE
        endif
        ! set saturation to one
        dX_p(saturation_index) = saturation0 - 1.d0
      endif
    endif
    ! DPRELIM is designed to catch large negative values in liquid pressure
    ! and cut the timestep if this occurs
    if (pressure1 <= (this%dprelim)) then  ! DEPLIMIT(2)
      cut_timestep = PETSC_TRUE
    endif
  enddo

  call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

  temp_int(:) = 0
  if (force_another_iteration) temp_int(1) = 1
  if (cut_timestep) temp_int(2) = 1
  call MPI_Allreduce(MPI_IN_PLACE,temp_int,TWO_INTEGER_MPI, &
                     MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
  option%convergence = CONVERGENCE_KEEP_ITERATING
  if (temp_int(1) > 0) option%convergence = CONVERGENCE_FORCE_ITERATION
  if (temp_int(2) > 0) option%convergence = CONVERGENCE_CUT_TIMESTEP

end subroutine PMWIPPFloCheckUpdatePre

! ************************************************************************** !

subroutine PMWIPPFloCheckUpdatePost(this,line_search,X0,dX,X1,dX_changed, &
                                    X1_changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
  use WIPP_Flow_Aux_module
  use Global_Aux_module
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  use Material_Aux_class  
  use Output_EKG_module
  
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
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  class(material_auxvar_type), pointer :: material_auxvars(:)  
  type(material_parameter_type), pointer :: material_parameter

  PetscReal, pointer :: X0_p(:)
  PetscReal, pointer :: X1_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum_p(:), accum_p2(:)
  PetscReal, pointer :: press_ptr(:)
  PetscReal, pointer :: sat_ptr(:)

  PetscInt :: local_id, ghosted_id
  PetscInt :: offset , idof
  PetscInt :: liquid_equation_index
  PetscInt :: gas_equation_index
  PetscInt :: saturation_index
  PetscInt :: pressure_index
  PetscInt :: temp_int(8)
  PetscReal :: temp_real(7)

  PetscReal, parameter :: zero_saturation = 1.d-15
  PetscReal, parameter :: zero_accumulation = 1.d-15

  PetscReal :: residual
  PetscReal :: accumulation
  PetscReal :: residual_over_accumulation
  PetscReal :: abs_X
  PetscReal :: abs_dX
  PetscReal :: abs_dX_TS
  PetscReal :: rel_dX_TS
  PetscBool :: converged_liquid_equation
  PetscBool :: converged_gas_equation
  PetscBool :: converged_liquid_pressure
  PetscBool :: converged_gas_saturation
  PetscReal :: max_liq_eq
  PetscReal :: max_gas_eq
  PetscReal :: max_liq_pres_rel_change
  PetscReal :: max_gas_sat_change_NI
  PetscReal :: max_gas_sat_change_TS
  PetscReal :: max_abs_pressure_change_TS
  PetscReal :: max_rel_pressure_change_TS
  PetscReal :: min_liq_pressure
  PetscReal :: min_gas_pressure
  PetscInt :: max_liq_eq_cell
  PetscInt :: max_gas_eq_cell
  PetscInt :: max_liq_pres_rel_change_cell
  PetscInt :: max_gas_sat_change_NI_cell
  PetscInt :: max_gas_sat_change_TS_cell
  PetscInt :: max_abs_pressure_change_TS_cell
  PetscInt :: max_rel_pressure_change_TS_cell
  PetscInt :: min_liq_pressure_cell
  PetscInt :: min_gas_pressure_cell
  PetscReal :: abs_dX_over_absX
  PetscReal :: pflotran_to_bragflo(2)
  character(len=4) :: reason
  PetscInt :: istart, iend
  
  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch
  wippflo_auxvars => patch%aux%WIPPFlo%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  material_parameter => patch%aux%Material%material_parameter

  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE
  
  call VecGetArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X1,X1_p,ierr);CHKERRQ(ierr)
  if (this%stored_residual_vec == PETSC_NULL_VEC) then
    call VecGetArrayReadF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
  else
    ! this vector is to be used if linear system scaling is employed in 
    ! BRAGFLO mode as the residual will be altered by that scaling.  This
    ! vector stores the original residual. 
    call VecGetArrayReadF90(this%stored_residual_vec,r_p,ierr);CHKERRQ(ierr)
  endif
  call VecGetArrayReadF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%flow_accum2,accum_p2,ierr);CHKERRQ(ierr)
  ! max change variables: [LIQUID_PRESSURE, GAS_PRESSURE, GAS_SATURATION]
  call VecGetArrayF90(field%max_change_vecs(1),press_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%max_change_vecs(3),sat_ptr,ierr);CHKERRQ(ierr)
  converged_liquid_equation = PETSC_TRUE
  converged_gas_equation = PETSC_TRUE
  converged_liquid_pressure = PETSC_TRUE
  converged_gas_saturation = PETSC_TRUE
  max_liq_eq = 0.d0
  max_gas_eq = 0.d0
  max_liq_pres_rel_change = 0.d0
  max_gas_sat_change_NI = 0.d0
  max_gas_sat_change_TS = 0.d0
  max_abs_pressure_change_TS = 0.d0
  max_rel_pressure_change_TS = 0.d0
  min_liq_pressure = 1.d20
  min_gas_pressure = 1.d20
  max_liq_eq_cell = 0
  max_gas_eq_cell = 0
  max_liq_pres_rel_change_cell = 0
  max_gas_sat_change_NI_cell = 0
  max_gas_sat_change_TS_cell = 0
  max_abs_pressure_change_TS_cell = 0
  max_rel_pressure_change_TS_cell = 0
  min_liq_pressure_cell = 0
  min_gas_pressure_cell = 0
  do local_id = 1, grid%nlmax
    offset = (local_id-1)*option%nflowdof
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    liquid_equation_index = offset + WIPPFLO_LIQUID_EQUATION_INDEX
    gas_equation_index = offset + WIPPFLO_GAS_EQUATION_INDEX
    pressure_index = offset + WIPPFLO_LIQUID_PRESSURE_DOF
    saturation_index = offset + WIPPFLO_GAS_SATURATION_DOF

    if (wippflo_debug_first_iteration) then
      pflotran_to_bragflo = fmw_comp * option%flow_dt / &
                            material_auxvars(ghosted_id)%volume
      print *, local_id, r_p(1:2) * pflotran_to_bragflo
    endif
    
    ! liquid component equation
    if (accum_p2(liquid_equation_index) > zero_accumulation) then 
      residual = r_p(liquid_equation_index)
      accumulation = accum_p2(liquid_equation_index)
      residual_over_accumulation = residual / accumulation
      if (max_liq_eq < dabs(residual_over_accumulation)) then
        max_liq_eq_cell = local_id
        max_liq_eq = residual_over_accumulation
      endif
      if (dabs(residual) > this%liquid_equation_tolerance) then
        if (dabs(residual_over_accumulation) > &
            this%liquid_equation_tolerance) then
          converged_liquid_equation = PETSC_FALSE
        endif
      endif
    endif

    ! gas component equation
    if (accum_p2(gas_equation_index) > zero_accumulation .and. &
        X1_p(gas_equation_index) > zero_saturation) then
      residual = r_p(gas_equation_index)
      accumulation = accum_p2(gas_equation_index)
      residual_over_accumulation = residual / accumulation
      if (max_gas_eq < dabs(residual_over_accumulation)) then
        max_gas_eq_cell = local_id
        max_gas_eq = residual_over_accumulation
      endif
      if (dabs(residual) > this%gas_equation_tolerance) then 
        if (dabs(residual_over_accumulation) > &
            this%gas_equation_tolerance) then
          converged_gas_equation = PETSC_FALSE
        endif
      endif
    endif

    ! maximum relative change in liquid pressure
    abs_X = dabs(X1_p(pressure_index))
    if (abs_X > 0.d0) then
      abs_dX_over_absX = dabs(dX_p(pressure_index))/abs_X
      if (dabs(max_liq_pres_rel_change) < abs_dX_over_absX) then
        max_liq_pres_rel_change_cell = local_id
        max_liq_pres_rel_change = dX_p(pressure_index)/abs_X
      endif
      if (abs_dX_over_absX >= this%liquid_pressure_tolerance) then
        converged_liquid_pressure = PETSC_FALSE
      endif
    endif
    
    ! EPS_SAT maximum relative gas saturation change "digits of accuracy"
    abs_dX = dabs(dX_p(saturation_index))
    if ((-1.d0*log10(abs_dX)) < this%eps_sat) then
      converged_gas_saturation = PETSC_FALSE
      if (dabs(max_gas_sat_change_NI) < abs_dX) then
        max_gas_sat_change_NI_cell = local_id
        max_gas_sat_change_NI = dX_p(saturation_index)
      endif
    endif
    
    ! DSAT_MAX maximum absolute gas saturation change over time step
    abs_dX_TS = dabs(sat_ptr(local_id)-X1_p(saturation_index))
    if (abs_dX_TS > 0.d0) then
      if (dabs(max_gas_sat_change_TS) < abs_dX_TS) then
        max_gas_sat_change_TS_cell = local_id
        max_gas_sat_change_TS = abs_dX_TS
      endif
    endif
    
    ! DPRE_MAX maximum absolute liquid pressure change over time step
    abs_dX_TS = dabs(press_ptr(local_id)-X1_p(pressure_index))
    if (abs_dX_TS > 0.d0) then
      if (dabs(max_abs_pressure_change_TS) < abs_dX_TS) then
        max_abs_pressure_change_TS_cell = local_id
        max_abs_pressure_change_TS = abs_dX_TS
      endif
    endif
    
    ! EPS_PRES maximum relative liquid pressure change over time step
    rel_dX_TS = dabs((press_ptr(local_id)-X1_p(pressure_index))/ &
                      press_ptr(local_id))
    if (rel_dX_TS > 0.d0) then
      if (dabs(max_rel_pressure_change_TS) < rel_dX_TS) then
        max_rel_pressure_change_TS_cell = local_id
        max_rel_pressure_change_TS = rel_dX_TS
      endif
    endif

    ! liquid pressure
    if (X1_p(pressure_index) < min_liq_pressure) then
      min_liq_pressure = X1_p(pressure_index)
      min_liq_pressure_cell = local_id
    endif

    ! gas pressure
    if (wippflo_auxvars(0,ghosted_id)%pres(2) < min_gas_pressure) then
      min_gas_pressure = wippflo_auxvars(0,ghosted_id)%pres(2)
      min_gas_pressure_cell = local_id
    endif
  enddo

  if (wippflo_debug_first_iteration) stop

  temp_int = 0
  if (.not.converged_liquid_equation) temp_int(1) = 1
  if (.not.converged_gas_equation) temp_int(2) = 1
  if (.not.converged_liquid_pressure) temp_int(3) = 1
  if (.not.converged_gas_saturation) temp_int(4) = 1
  if (min_gas_pressure < 0.d0) temp_int(5) = 1
  if (converged_gas_saturation .and. &
      (max_gas_sat_change_TS > this%dsat_max)) temp_int(6) = 1
  if (converged_liquid_pressure .and. &
      (max_abs_pressure_change_TS > this%dpres_max)) temp_int(7) = 1
  if (converged_liquid_pressure .and. &
      (max_rel_pressure_change_TS > this%eps_pres)) temp_int(8) = 1
  call MPI_Allreduce(MPI_IN_PLACE,temp_int,8, &
                     MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
  temp_real(1) = max_liq_eq
  temp_real(2) = max_gas_eq
  temp_real(3) = max_liq_pres_rel_change
  temp_real(4) = max_gas_sat_change_NI
  temp_real(5) = max_gas_sat_change_TS
  temp_real(6) = max_abs_pressure_change_TS
  temp_real(7) = max_rel_pressure_change_TS
  call MPI_Allreduce(MPI_IN_PLACE,temp_real,SEVEN_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
  max_liq_eq = temp_real(1)
  max_gas_eq = temp_real(2)
  max_liq_pres_rel_change = temp_real(3)
  max_gas_sat_change_NI = temp_real(4)
  max_gas_sat_change_TS = temp_real(5)
  max_abs_pressure_change_TS = temp_real(6)
  max_rel_pressure_change_TS = temp_real(7)

  if (option%convergence == CONVERGENCE_CUT_TIMESTEP) then
    reason = '!cut'
  else if (temp_int(8) > 0) then
    reason = 'PrTS'
    option%convergence = CONVERGENCE_CUT_TIMESTEP
  else if (temp_int(7) > 0) then
    reason = 'PaTS'
    option%convergence = CONVERGENCE_CUT_TIMESTEP
  else if (temp_int(6) > 0) then
    reason = 'SgTS'
    option%convergence = CONVERGENCE_CUT_TIMESTEP
  else if (temp_int(5) > 0) then
    reason = 'negP'
    option%convergence = CONVERGENCE_CUT_TIMESTEP
  else if (option%convergence == CONVERGENCE_FORCE_ITERATION) then
    reason = '!it '
  else if (maxval(temp_int) > 0) then
    option%convergence = CONVERGENCE_KEEP_ITERATING
    reason = '    '
    if (temp_int(1) > 0) reason(1:1) = 'L'
    if (temp_int(2) > 0) reason(2:2) = 'G'
    if (temp_int(3) > 0) reason(3:3) = 'P'
    if (temp_int(4) > 0) reason(4:4) = 'S'
  else
    reason = '----'
    option%convergence = CONVERGENCE_CONVERGED
  endif
  if (OptionPrintToScreen(option)) then
    if (option%mycommsize > 1) then
      write(*,'(4x,"Reason: ",a4,4es10.2)') reason, &
        max_liq_eq, max_gas_eq, &
        max_liq_pres_rel_change, max_gas_sat_change_NI
    else
      write(*,'(4x,"Reason: ",a4,4(i5,es10.2))') reason, &
        max_liq_eq_cell, max_liq_eq, &
        max_gas_eq_cell, max_gas_eq, &
        max_liq_pres_rel_change_cell, max_liq_pres_rel_change, &
        max_gas_sat_change_NI_cell, max_gas_sat_change_NI
! for debugging
#if 0
      local_id = min_gas_pressure_cell
      offset = (local_id-1)*option%nflowdof
      istart = offset + 1
      iend = offset + option%nflowdof
      ghosted_id = grid%nL2G(local_id)
      write(*,'(4x,i5,3es11.3)') local_id, &
        wippflo_auxvars(0,ghosted_id)%pres(1:2), &
        wippflo_auxvars(0,ghosted_id)%pres(option%capillary_pressure_id)
#endif
    endif
  endif

  call VecRestoreArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X1,X1_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%flow_accum2,accum_p2,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%max_change_vecs(1),sat_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%max_change_vecs(3),sat_ptr,ierr);CHKERRQ(ierr)
                               
end subroutine PMWIPPFloCheckUpdatePost

! ************************************************************************** !

subroutine PMWIPPFloTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use WIPP_Flow_module, only : WIPPFloTimeCut

  implicit none
  
  class(pm_wippflo_type) :: this
  
  call PMSubsurfaceFlowTimeCut(this)
  call WIPPFloTimeCut(this%realization)

end subroutine PMWIPPFloTimeCut

! ************************************************************************** !

subroutine PMWIPPFloUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use WIPP_Flow_module, only : WIPPFloUpdateSolution, &
                             WIPPFloMapBCAuxVarsToGlobal

  implicit none
  
  class(pm_wippflo_type) :: this
  
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
