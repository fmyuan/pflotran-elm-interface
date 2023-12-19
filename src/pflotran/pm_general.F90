module PM_General_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class

  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscInt, parameter :: ABS_UPDATE_INDEX = 1
  PetscInt, parameter :: REL_UPDATE_INDEX = 2
  PetscInt, parameter :: RESIDUAL_INDEX = 3
  PetscInt, parameter :: SCALED_RESIDUAL_INDEX = 4
  PetscInt, parameter :: MAX_INDEX = SCALED_RESIDUAL_INDEX
  PetscInt, pointer :: general_max_states
  PetscInt, pointer :: max_change_index
  PetscBool, public :: general_g_state_air_mass_dof = PETSC_FALSE
  PetscBool, public :: general_set_solute = PETSC_FALSE

  PetscBool, public :: transient_porosity

  type, public, extends(pm_subsurface_flow_type) :: pm_general_type
    PetscInt, pointer :: max_change_ivar(:)
    PetscInt, pointer :: max_change_isubvar(:)
    PetscBool, pointer :: converged_flag(:,:,:)
    PetscInt, pointer :: converged_cell(:,:,:)
    PetscReal, pointer :: converged_real(:,:,:)
    PetscReal, pointer :: residual_abs_inf_tol(:)
    PetscReal, pointer :: residual_scaled_inf_tol(:)
    PetscReal, pointer :: abs_update_inf_tol(:,:)
    PetscReal, pointer :: rel_update_inf_tol(:,:)
    PetscReal :: damping_factor
  contains
    procedure, public :: ReadSimulationOptionsBlock => &
                           PMGeneralReadSimOptionsBlock
    procedure, public :: ReadNewtonBlock => PMGeneralReadNewtonSelectCase
    procedure, public :: InitializeSolver => PMGeneralInitializeSolver
    procedure, public :: InitializeRun => PMGeneralInitializeRun
    procedure, public :: InitializeTimestep => PMGeneralInitializeTimestep
    procedure, public :: Residual => PMGeneralResidual
    procedure, public :: Jacobian => PMGeneralJacobian
    procedure, public :: UpdateTimestep => PMGeneralUpdateTimestep
    procedure, public :: PreSolve => PMGeneralPreSolve
    procedure, public :: PostSolve => PMGeneralPostSolve
    procedure, public :: CheckUpdatePre => PMGeneralCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMGeneralCheckUpdatePost
    procedure, public :: CheckConvergence => PMGeneralCheckConvergence
    procedure, public :: TimeCut => PMGeneralTimeCut
    procedure, public :: UpdateSolution => PMGeneralUpdateSolution
    procedure, public :: UpdateAuxVars => PMGeneralUpdateAuxVars
    procedure, public :: MaxChange => PMGeneralMaxChange
    procedure, public :: ComputeMassBalance => PMGeneralComputeMassBalance
    procedure, public :: InputRecord => PMGeneralInputRecord
    procedure, public :: CheckpointBinary => PMGeneralCheckpointBinary
    procedure, public :: RestartBinary => PMGeneralRestartBinary
    procedure, public :: Destroy => PMGeneralDestroy
  end type pm_general_type

  public :: PMGeneralCreate, &
            PMGeneralSetFlowMode

contains

! ************************************************************************** !

function PMGeneralCreate()
  !
  ! Creates General process models shell
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  use Upwind_Direction_module
  use Option_module

  implicit none

  class(pm_general_type), pointer :: PMGeneralCreate

  class(pm_general_type), pointer :: this

#ifdef PM_GENERAL_DEBUG
  print *, 'PMGeneralCreate()'
#endif

  allocate(this)
  call PMSubsurfaceFlowInit(this)
  this%name = 'General Multiphase Flow'
  this%header = 'GENERAL MULTIPHASE FLOW'

  ! turn off default upwinding which is set to PETSC_TRUE in
  !  upwind_direction.F90
  fix_upwind_direction = PETSC_FALSE
  ! set default to infinity norm convergence
  this%check_post_convergence = PETSC_TRUE

  PMGeneralCreate => this

end function PMGeneralCreate

! ************************************************************************** !

subroutine PMGeneralSetFlowMode(pm,option)
  !
  ! Sets the flow mode for general mode
  !
  ! Author: David Fukuyama
  ! Date: 06/10/21
  !

  use Option_module
  use Variables_module, only : LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
                               LIQUID_MOLE_FRACTION, TEMPERATURE, &
                               GAS_SATURATION, PRECIPITATE_SATURATION, &
                               LIQUID_SATURATION, POROSITY

  implicit none

  type(option_type) :: option
  class(pm_general_type) :: pm!, pointer, intent(inout) :: pm
  !class(pm_general_type), pointer, intent(inout) :: pm

  PetscReal, parameter :: ref_temp = 20.d0 !degrees C
  PetscReal, parameter :: ref_pres = 101325.d0 !Pa
  PetscReal, parameter :: ref_sat = 0.5
  PetscReal, parameter :: ref_xmol = 1.d-6

  !MAN optimized:
  PetscReal, parameter :: pres_abs_inf_tol = 1.d0 ! Reference tolerance [Pa]
  PetscReal, parameter :: temp_abs_inf_tol = 1.d-5
  PetscReal, parameter :: sat_abs_inf_tol = 1.d-5
  PetscReal, parameter :: xmol_abs_inf_tol = 1.d-9
  PetscReal, parameter :: por_abs_inf_tol = 1.d-5

  PetscReal, parameter :: pres_rel_inf_tol = 1.d-3
  PetscReal, parameter :: temp_rel_inf_tol = 1.d-3
  PetscReal, parameter :: sat_rel_inf_tol = 1.d-3
  PetscReal, parameter :: xmol_rel_inf_tol = 1.d-3

  !MAN optimized:
  PetscReal, parameter :: w_mass_abs_inf_tol = 1.d-5 !1.d-7 !kmol_water/sec
  PetscReal, parameter :: a_mass_abs_inf_tol = 1.d-5 !1.d-7
  PetscReal, parameter :: u_abs_inf_tol = 1.d-5 !1.d-7
  PetscReal, parameter :: s_mass_abs_inf_tol = 1.d-5

  PetscReal, parameter :: ref_density_w = 55.058 !kmol_water/m^3
  PetscReal, parameter :: ref_density_a = 0.0423 !kmol_air/m^3
  PetscReal, parameter :: ref_u = 83.8 !MJ/m^3

!   PetscReal, pointer :: residual_abs_inf_tol(:)
!   PetscReal, pointer :: residual_scaled_inf_tol(:)
!   PetscReal, pointer :: abs_update_inf_tol(:,:)
!   PetscReal, pointer :: rel_update_inf_tol(:,:)

   PetscReal, allocatable :: residual_abs_inf_tol(:)
   PetscReal, allocatable :: residual_scaled_inf_tol(:)
   PetscReal, allocatable :: abs_update_inf_tol(:,:)
   PetscReal, allocatable :: rel_update_inf_tol(:,:)
  option%iflowmode = G_MODE
  allocate(general_max_states)
  allocate(max_change_index)


  option%use_isothermal = PETSC_FALSE

  option%liquid_phase = 1  ! liquid_pressure
  option%gas_phase = 2     ! gas_pressure

  option%air_pressure_id = 3
  option%capillary_pressure_id = 4
  option%vapor_pressure_id = 5
  option%saturation_pressure_id = 6

  option%water_id = 1
  option%air_id = 2

  if (option%nflowdof == 0 .or. option%nflowdof == 3) then
    option%nphase = 2
    option%nflowdof = 3
    option%nflowspec = 2
    general_max_states = 3
    max_change_index = 6
    option%energy_id = 3
  elseif (option%nflowdof == 4) then
    option%nflowspec = 3
    option%salt_id = 3
    option%energy_id = 4
    option%nphase = 3
    option%precipitate_phase = 3
    general_max_states = 7
    max_change_index = 8
    if (.not.general_set_solute) then
      option%io_buffer = 'Solute must be acknowledged in the OPTIONS block &
                          &of GENERAL MODE.'
      call PrintErrMsg(option)
    endif
  else
    write(option%io_buffer,*) option%nflowdof
    option%io_buffer = trim(adjustl(option%io_buffer)) // &
     ' equations not currently supported in GENERAL mode.'
    call PrintErrMsgByRank(option)
  endif
  allocate(abs_update_inf_tol(option%nflowdof,general_max_states))
  allocate(rel_update_inf_tol(option%nflowdof,general_max_states))
  allocate(residual_abs_inf_tol(option%nflowdof))
  allocate(residual_scaled_inf_tol(option%nflowdof))

  allocate(pm%abs_update_inf_tol(option%nflowdof,general_max_states))
  allocate(pm%rel_update_inf_tol(option%nflowdof,general_max_states))
  allocate(pm%residual_abs_inf_tol(option%nflowdof))
  allocate(pm%residual_scaled_inf_tol(option%nflowdof))

  allocate(pm%converged_flag(option%nflowdof,general_max_states,MAX_INDEX))
  allocate(pm%converged_cell(option%nflowdof,general_max_states,MAX_INDEX))
  allocate(pm%converged_real(option%nflowdof,general_max_states,MAX_INDEX))

  if (option%nflowdof == 3) then
    rel_update_inf_tol = &
      reshape([pres_rel_inf_tol,xmol_rel_inf_tol,temp_rel_inf_tol, &
               pres_rel_inf_tol,pres_rel_inf_tol,temp_rel_inf_tol, &
               pres_rel_inf_tol,sat_rel_inf_tol,temp_rel_inf_tol], &
               shape(rel_update_inf_tol))*1.d0 !change to 0 to zero tolerances
    abs_update_inf_tol = &
      reshape([pres_abs_inf_tol,xmol_abs_inf_tol,temp_abs_inf_tol, &
               pres_abs_inf_tol,pres_abs_inf_tol,temp_abs_inf_tol, &
               pres_abs_inf_tol,sat_abs_inf_tol,temp_abs_inf_tol], &
               shape(abs_update_inf_tol))*1.d0 !change to 0 to zero tolerances
    residual_abs_inf_tol(:) = [w_mass_abs_inf_tol,a_mass_abs_inf_tol,&
                               u_abs_inf_tol]
    residual_scaled_inf_tol(:) = 1.d-6
    allocate(pm%max_change_ivar(6))
    pm%max_change_ivar = [LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
                          LIQUID_MOLE_FRACTION, TEMPERATURE, &
                          GAS_SATURATION]
    allocate(pm%max_change_isubvar(6))
                                     ! UNINITIALIZED_INTEGER avoids zeroing of
                                     ! pressures not represented in phase
                                         ! 2 = air in xmol(air,liquid)
    pm%max_change_isubvar = [0,0,0,2,0,0]
  elseif (option%nflowdof == 4) then
    rel_update_inf_tol = &
      reshape([pres_rel_inf_tol,xmol_rel_inf_tol,temp_rel_inf_tol, &
               xmol_rel_inf_tol,&
               pres_rel_inf_tol,pres_rel_inf_tol,temp_rel_inf_tol,999.d0,&
               pres_rel_inf_tol,sat_rel_inf_tol,temp_rel_inf_tol, &
               xmol_rel_inf_tol,&
               pres_rel_inf_tol,999.d0,temp_rel_inf_tol,999.d0,&
               pres_rel_inf_tol,xmol_rel_inf_tol,temp_rel_inf_tol, &
               sat_rel_inf_tol,&
               pres_rel_inf_tol,pres_rel_inf_tol,temp_rel_inf_tol, &
               sat_rel_inf_tol,&
               pres_rel_inf_tol,sat_rel_inf_tol,temp_rel_inf_tol, &
               sat_rel_inf_tol], &
               shape(rel_update_inf_tol)) * &
               1.d0 ! change to 0.d0 to zero tolerances
    abs_update_inf_tol = &
      reshape([pres_abs_inf_tol,xmol_abs_inf_tol,temp_abs_inf_tol,&
               xmol_abs_inf_tol,&
               pres_abs_inf_tol,pres_abs_inf_tol,temp_abs_inf_tol,999.d0,&
               pres_abs_inf_tol,sat_abs_inf_tol,temp_abs_inf_tol, &
               xmol_abs_inf_tol,&
               pres_abs_inf_tol,999.d0,temp_abs_inf_tol,999.d0,&
               pres_abs_inf_tol,xmol_abs_inf_tol,temp_abs_inf_tol, &
               sat_abs_inf_tol,&
               pres_abs_inf_tol,pres_abs_inf_tol,temp_abs_inf_tol, &
               sat_abs_inf_tol,&
               pres_abs_inf_tol,sat_abs_inf_tol,temp_abs_inf_tol, &
               sat_abs_inf_tol], &
               shape(abs_update_inf_tol)) * &
               1.d0 ! change to 0.d0 to zero tolerances
    residual_abs_inf_tol = [w_mass_abs_inf_tol,a_mass_abs_inf_tol,&
                            s_mass_abs_inf_tol,u_abs_inf_tol]
    residual_scaled_inf_tol(:) = 1.d-6
    allocate(pm%max_change_ivar(8))
    pm%max_change_ivar = [LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
                          LIQUID_MOLE_FRACTION, LIQUID_MOLE_FRACTION, &
                          TEMPERATURE, GAS_SATURATION, POROSITY]
    allocate(pm%max_change_isubvar(8))
    pm%max_change_isubvar = [0,0,0,2,3,0,0,0]
  endif
  pm%damping_factor = -1.d0
  pm%residual_abs_inf_tol = residual_abs_inf_tol
  pm%residual_scaled_inf_tol = residual_scaled_inf_tol
  pm%abs_update_inf_tol = abs_update_inf_tol
  pm%rel_update_inf_tol = rel_update_inf_tol

  pm%converged_flag(:,:,:) = PETSC_TRUE
  pm%converged_real(:,:,:) = 0.d0
  pm%converged_cell(:,:,:) = 0

  if (general_g_state_air_mass_dof) then
    pm%abs_update_inf_tol(2,2)=pm%abs_update_inf_tol(2,1)
    pm%rel_update_inf_tol(2,2)=pm%rel_update_inf_tol(2,1)
  endif

end subroutine PMGeneralSetFlowMode

! ************************************************************************** !
subroutine PMGeneralReadSimOptionsBlock(this,input)
  !
  ! Read simulation options for GENERAL mode
  !
  ! Author: Glenn Hammond
  ! Date: 03/04/15
  !
  use General_module
  use General_Aux_module
  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: keyword, word
  class(pm_general_type) :: this
  type(option_type), pointer :: option
  PetscReal :: tempreal
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
  PetscInt :: lid, gid, eid, sid

  option => this%option

  lid = option%liquid_phase
  gid = option%gas_phase
  eid = option%energy_id
  sid = option%salt_id

  error_string = 'General Options'

  input%ierr = 0
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    found = PETSC_FALSE
    call PMSubsurfFlowReadSimOptionsSC(this,input,keyword,found, &
                                       error_string,option)
    if (found) cycle

    select case(trim(keyword))
      case('DIFFUSE_XMASS')
        general_diffuse_xmol = PETSC_FALSE
      case('ARITHMETIC_GAS_DIFFUSIVE_DENSITY')
        general_harmonic_diff_density = PETSC_FALSE
      case('CHECK_MAX_DPL_LIQ_STATE_ONLY')
        gen_chk_max_dpl_liq_state_only = PETSC_TRUE
      case('DEBUG_CELL')
        call InputReadInt(input,option,general_debug_cell_id)
        call InputErrorMsg(input,option,keyword,error_string)
      case('GAS_COMPONENT_FORMULA_WEIGHT')
        !geh: assuming gas component is index 2
        call InputReadDouble(input,option,fmw_comp(2))
        call InputErrorMsg(input,option,keyword,error_string)
      case('GAS_STATE_AIR_MASS_DOF')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,keyword,error_string)
        call GeneralAuxSetAirMassDOF(word,option)
        general_g_state_air_mass_dof = PETSC_TRUE
      case('HARMONIC_GAS_DIFFUSIVE_DENSITY')
        general_harmonic_diff_density = PETSC_TRUE
      case('NEWTONTRDC_HOLD_INNER_ITERATIONS',&
           'HOLD_INNER_ITERATIONS','NEWTONTRDC_HOLD_INNER')
        !heeho: only used when using newtontrd-c
        general_newtontrdc_hold_inner = PETSC_TRUE
      case('IMMISCIBLE')
        general_immiscible = PETSC_TRUE
      case('ISOTHERMAL')
        general_isothermal = PETSC_TRUE
      case('NON_DARCY_FLOW')
        general_non_darcy_flow = PETSC_TRUE
      case('NON_DARCY_FLOW_A')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,error_string)
        non_darcy_A = tempreal
      case('NON_DARCY_FLOW_B')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,error_string)
        non_darcy_B = tempreal
      case('LIQUID_COMPONENT_FORMULA_WEIGHT')
        !heeho: assuming liquid component is index 1
        call InputReadDouble(input,option,fmw_comp(1))
        call InputErrorMsg(input,option,keyword,error_string)
      case('NO_AIR')
        general_no_air = PETSC_TRUE
      case('NO_STATE_TRANSITION_OUTPUT')
        general_print_state_transition = PETSC_FALSE
      case('NO_TEMP_DEPENDENT_DIFFUSION')
        general_temp_dep_gas_air_diff = PETSC_FALSE
      case('PHASE_CHANGE_EPSILON')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,error_string)
        general_phase_chng_epsilon = tempreal
      case('RESTRICT_STATE_CHANGE')
        general_restrict_state_chng = PETSC_TRUE
      case('TWO_PHASE_ENERGY_DOF')
        call InputKeywordDeprecated('TWO_PHASE_ENERGY_DOF', &
                                    'TWO_PHASE_STATE_ENERGY_DOF',option)
      case('TWO_PHASE_STATE_ENERGY_DOF')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,keyword,error_string)
        call GeneralAuxSetEnergyDOF(word,option)
      case('WINDOW_EPSILON')
        call InputReadDouble(input,option,window_epsilon)
        call InputErrorMsg(input,option,keyword,error_string)
     case('CALCULATE_SURFACE_TENSION')
        general_compute_surface_tension = PETSC_TRUE
      case('VAPOR_PRESSURE_KELVIN')
        general_kelvin_equation = PETSC_TRUE
      case('SOLUBLE_MATRIX')
        option%io_buffer = "Matrix solubility is now entered in material properties as SOLUBLE."
        call PrintErrMsg(option)
      case('UPDATE_PERMEABILITY')
        general_update_permeability = PETSC_TRUE
        call InputReadDouble(input,option,permeability_func_porosity_exp)
        call InputErrorMsg(input,option,keyword,error_string)
      case('SOLUTE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'solute',error_string)
        call GeneralAuxSetSolute(word,option)
        general_salt = PETSC_TRUE
        general_set_solute = PETSC_TRUE
        option%nflowdof = FOUR_INTEGER
      case default
        call InputKeywordUnrecognized(input,keyword,'GENERAL Mode',option)
    end select

  enddo
  call InputPopBlock(input,option)

  if (general_isothermal .and. &
      general_2ph_energy_dof == GENERAL_AIR_PRESSURE_INDEX) then
    option%io_buffer = 'Isothermal GENERAL mode may only be run with ' // &
                       'temperature as the two phase energy dof.'
    call PrintErrMsg(option)
  endif

end subroutine PMGeneralReadSimOptionsBlock

! ************************************************************************** !

subroutine PMGeneralReadNewtonSelectCase(this,input,keyword,found, &
                                         error_string,option)
  !
  ! Reads input file parameters associated with the GENERAL process model
  ! Newton solver convergence
  !
  ! Author: Glenn Hammond
  ! Date: 03/23/20

  use Input_Aux_module
  use String_module
  use Utility_module
  use Option_module
  use General_Aux_module

  implicit none

  class(pm_general_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option

  PetscBool :: found
  PetscReal :: tempreal
  PetscInt :: lid, gid, eid, sid

  option => this%option

  lid = option%liquid_phase
  gid = option%gas_phase
  eid = option%energy_id
  sid = option%salt_id

  error_string = 'GENERAL Newton Solver'

  found = PETSC_FALSE
  call PMSubsurfaceFlowReadNewtonSelectCase(this,input,keyword,found, &
                                            error_string,option)
  if (found) return

  found = PETSC_TRUE
  select case(trim(keyword))
    case('MAX_NEWTON_ITERATIONS')
      call InputKeywordDeprecated('MAX_NEWTON_ITERATIONS', &
                                  'MAXIMUM_NUMBER_OF_ITERATIONS.',option)
    ! Tolerances

    ! All Residual
    case('RESIDUAL_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%residual_abs_inf_tol(:) = tempreal
      this%residual_scaled_inf_tol(:) = tempreal

    ! Absolute Residual
    case('RESIDUAL_ABS_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%residual_abs_inf_tol(:) = tempreal
    case('LIQUID_RESIDUAL_ABS_INF_TOL')
      call InputReadDouble(input,option,this%residual_abs_inf_tol(lid))
      call InputErrorMsg(input,option,keyword,error_string)
    case('GAS_RESIDUAL_ABS_INF_TOL')
      call InputReadDouble(input,option,this%residual_abs_inf_tol(gid))
      call InputErrorMsg(input,option,keyword,error_string)
    case('ENERGY_RESIDUAL_ABS_INF_TOL')
      call InputReadDouble(input,option,this%residual_abs_inf_tol(eid))
      call InputErrorMsg(input,option,keyword,error_string)

    ! Scaled Residual
    case('ITOL_SCALED_RESIDUAL')
      call InputKeywordDeprecated('ITOL_SCALED_RESIDUAL', &
                                  'RESIDUAL_SCALED_INF_TOL',option)
    case('RESIDUAL_SCALED_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%residual_scaled_inf_tol(:) = tempreal
    case('LIQUID_RESIDUAL_SCALED_INF_TOL')
      call InputReadDouble(input,option,this%residual_scaled_inf_tol(lid))
      call InputErrorMsg(input,option,keyword,error_string)
    case('GAS_RESIDUAL_SCALED_INF_TOL')
      call InputReadDouble(input,option,this%residual_scaled_inf_tol(gid))
      call InputErrorMsg(input,option,keyword,error_string)
    case('ENERGY_RESIDUAL_SCALED_INF_TOL')
      call InputReadDouble(input,option,this%residual_scaled_inf_tol(eid))
      call InputErrorMsg(input,option,keyword,error_string)

    ! All Updates
    case('UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%abs_update_inf_tol(:,:) = tempreal
      this%rel_update_inf_tol(:,:) = tempreal

    ! Absolute Updates
    case('ABS_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%abs_update_inf_tol(:,:) = tempreal
    case('PRES_ABS_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%abs_update_inf_tol(1,:) = tempreal
      this%abs_update_inf_tol(2,2) = tempreal
    case('TEMP_ABS_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%abs_update_inf_tol(3,:) = tempreal
    case('SAT_ABS_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%abs_update_inf_tol(2,3) = tempreal
    case('XMOL_ABS_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%abs_update_inf_tol(2,1) = tempreal
    case('LIQUID_PRES_ABS_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%abs_update_inf_tol(1,1) = tempreal
    case('GAS_PRES_ABS_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%abs_update_inf_tol(1,2:3) = tempreal
    case('AIR_PRES_ABS_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%abs_update_inf_tol(2,2) = tempreal

    ! Relative Updates
    case('ITOL_RELATIVE_UPDATE')
      call InputKeywordDeprecated('ITOL_RELATIVE_UPDATE', &
                                  'REL_UPDATE_INF_TOL',option)
    case('REL_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%rel_update_inf_tol(:,:) = tempreal
    case('PRES_REL_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%rel_update_inf_tol(1,:) = tempreal
      this%rel_update_inf_tol(2,2) = tempreal
    case('TEMP_REL_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%rel_update_inf_tol(3,:) = tempreal
    case('SAT_REL_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%rel_update_inf_tol(2,3) = tempreal
    case('XMOL_REL_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%rel_update_inf_tol(2,1) = tempreal
    case('LIQUID_PRES_REL_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%rel_update_inf_tol(1,1) = tempreal
    case('GAS_PRES_REL_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%rel_update_inf_tol(1,2:3) = tempreal
    case('AIR_PRES_REL_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%rel_update_inf_tol(2,2) = tempreal
    case('MAXIMUM_PRESSURE_CHANGE')
      call InputReadDouble(input,option,general_max_pressure_change)
      call InputErrorMsg(input,option,keyword,error_string)
    case('MAX_ITERATION_BEFORE_DAMPING')
      call InputReadInt(input,option,general_max_it_before_damping)
      call InputErrorMsg(input,option,keyword,error_string)
    case('DAMPING_FACTOR')
      call InputReadDouble(input,option,general_damping_factor)
      call InputErrorMsg(input,option,keyword,error_string)
      this%damping_factor = general_damping_factor
    case default
      found = PETSC_FALSE

  end select

end subroutine PMGeneralReadNewtonSelectCase

! ************************************************************************** !

subroutine PMGeneralInitializeSolver(this)
  !
  ! Author: Glenn Hammond
  ! Date: 04/06/20

  use Solver_module

  implicit none

  class(pm_general_type) :: this

  call PMBaseInitializeSolver(this)

  ! helps accommodate rise in residual due to change in state
  this%solver%newton_dtol = 1.d9
  this%solver%newton_max_iterations = 8

end subroutine PMGeneralInitializeSolver

! ************************************************************************** !

recursive subroutine PMGeneralInitializeRun(this)
  !
  ! Initializes the time stepping
  !
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use Realization_Base_class

  implicit none

  class(pm_general_type) :: this

  PetscInt :: i
  PetscErrorCode :: ierr

  ! need to allocate vectors for max change
  call VecDuplicateVecsF90(this%realization%field%work,max_change_index, &
                           this%realization%field%max_change_vecs, &
                           ierr);CHKERRQ(ierr)
  ! set initial values
  do i = 1, max_change_index
    call RealizationGetVariable(this%realization, &
                                this%realization%field%max_change_vecs(i), &
                                this%max_change_ivar(i), &
                                this%max_change_isubvar(i))
  enddo

  ! call parent implementation
  call PMSubsurfaceFlowInitializeRun(this)

end subroutine PMGeneralInitializeRun

! ************************************************************************** !

subroutine PMGeneralInitializeTimestep(this)
  !
  ! Should not need this as it is called in PreSolve.
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  use General_module, only : GeneralInitializeTimestep
  use Global_module
  use Variables_module, only : TORTUOSITY
  use Material_module, only : MaterialAuxVarCommunicate
  use Option_module

  implicit none

  class(pm_general_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)
!geh:remove   everywhere
  call MaterialAuxVarCommunicate(this%comm1, &
                                 this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc,TORTUOSITY, &
                                 ZERO_INTEGER)

  call GeneralInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)

end subroutine PMGeneralInitializeTimestep

! ************************************************************************** !

subroutine PMGeneralPreSolve(this)
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13

  implicit none

  class(pm_general_type) :: this

  call PMSubsurfaceFlowPreSolve(this)

end subroutine PMGeneralPreSolve

! ************************************************************************** !

subroutine PMGeneralPostSolve(this)
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13

  implicit none

  class(pm_general_type) :: this

end subroutine PMGeneralPostSolve

! ************************************************************************** !

subroutine PMGeneralUpdateTimestep(this,update_dt, &
                                   dt,dt_min,dt_max,iacceleration, &
                                   num_newton_iterations,tfac, &
                                   time_step_max_growth_factor)
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  use Realization_Base_class, only : RealizationGetVariable
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL
  use Field_module
  use Global_module, only : GlobalSetAuxVarVecLoc
  use Variables_module, only : LIQUID_SATURATION, GAS_SATURATION
  use Utility_module, only : Equal
  use Option_module
  use String_module

  implicit none

  class(pm_general_type) :: this
  PetscBool :: update_dt
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor

  PetscReal :: fac
  PetscInt :: ifac
  PetscReal :: up, ut, ux, us, umin
  PetscReal :: dtt
  PetscReal :: governed_dt
  PetscReal :: umin_scale
  PetscReal :: value
  PetscReal :: governor_value
  character(MAXSTRINGLENGTH) :: string
  type(field_type), pointer :: field

#ifdef PM_GENERAL_DEBUG
  call PrintMsg(this%option,'PMGeneral%UpdateTimestep()')
#endif

  if (update_dt .and. iacceleration /= 0) then
    fac = 0.5d0
    if (num_newton_iterations >= iacceleration) then
      fac = 0.33d0
      umin = 0.d0
    else
      up = this%pressure_change_governor/(this%max_pressure_change+0.1)
      ut = this%temperature_change_governor/(this%max_temperature_change+1.d-5)
      ux = this%xmol_change_governor/(this%max_xmol_change+1.d-5)
      us = this%saturation_change_governor/(this%max_saturation_change+1.d-5)
      umin = min(up,ut,ux,us)
    endif
    ifac = max(min(num_newton_iterations,size(tfac)),1)
    umin_scale = fac * (1.d0 + umin)
    governed_dt = umin_scale * dt
    dtt = min(time_step_max_growth_factor*dt,governed_dt)
    dt = min(dtt,tfac(ifac)*dt,dt_max)
    dt = max(dt,dt_min)

    ! Inform user that time step is being limited by a state variable.
    if (Equal(dt,governed_dt)) then
      umin = umin * (1.d0 + 1.d-8)
      if (up < umin) then
        string = 'Pressure'
        value = this%max_pressure_change
        governor_value = this%pressure_change_governor
      else if (ut < umin) then
        string = 'Temperature'
        value = this%max_temperature_change
        governor_value = this%temperature_change_governor
      else if (ux < umin) then
        string = 'Mole Fraction'
        value = this%max_xmol_change
        governor_value = this%xmol_change_governor
      else if (us < umin) then
        string = 'Saturation'
        value = this%max_saturation_change
        governor_value = this%saturation_change_governor
      else
        string = 'Unknown'
        value = -999.d0
        governor_value = -999.d0
      endif
      string = ' Dt limited by ' // trim(string) // ': Val=' // &
        trim(StringWriteF('(es10.3)',value)) // ', Gov=' // &
        trim(StringWriteF('(es10.3)',governor_value)) // ', Scale=' // &
        trim(StringWriteF('(f4.2)',umin_scale))
      if (OptionPrintToScreen(this%option)) then
        write(*,'(a,/)') trim(string)
      endif
      if (OptionPrintToFile(this%option)) then
        write(this%option%fid_out,'(a,/)') trim(string)
      endif
    endif
  endif

  if (Initialized(this%cfl_governor)) then
    ! Since saturations are not stored in global_auxvar for general mode, we
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
    call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt,dt_max)
  endif

end subroutine PMGeneralUpdateTimestep

! ************************************************************************** !

subroutine PMGeneralResidual(this,snes,xx,r,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  use General_module, only : GeneralResidual

  implicit none

  class(pm_general_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr

  call PMSubsurfaceFlowUpdatePropertiesNI(this)
  call GeneralResidual(snes,xx,r,this%realization,ierr)

end subroutine PMGeneralResidual

! ************************************************************************** !

subroutine PMGeneralJacobian(this,snes,xx,A,B,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  use General_module, only : GeneralJacobian

  implicit none

  class(pm_general_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr

  call GeneralJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMGeneralJacobian

! ************************************************************************** !

subroutine PMGeneralCheckUpdatePre(this,snes,X,dX,changed,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 11/21/18
  !

  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Option_module
  use Saturation_Function_module
  use Patch_module
  use General_Aux_module
  use Global_Aux_module

  implicit none

  class(pm_general_type) :: this
  SNES :: snes
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr

  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(general_auxvar_type), pointer :: gen_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: offset
  PetscInt :: pgas_index, xmol_index, smol_index, pw_index
  PetscInt :: saturation_index
  PetscReal :: temp_real

  PetscReal, parameter :: ALMOST_ZERO = 1.d-10
  PetscReal, parameter :: ALMOST_ONE = 1.d0-ALMOST_ZERO

  PetscReal, pointer :: X_p(:),dX_p(:)

  type(field_type), pointer :: field
  PetscInt :: liquid_pressure_index, gas_pressure_index, air_pressure_index
  PetscInt :: temperature_index
  PetscInt :: apid, spid
  PetscReal :: liquid_pressure0, liquid_pressure1, del_liquid_pressure
  PetscReal :: gas_pressure0, gas_pressure1, del_gas_pressure
  PetscReal :: air_pressure0, air_pressure1, del_air_pressure
  PetscReal :: temperature0, temperature1, del_temperature
  PetscReal :: saturation0, saturation1, del_saturation
  PetscReal :: max_saturation_change = 0.125d0
  PetscReal :: scale, temp_scale
  PetscReal, parameter :: initial_scale = 1.d0
  PetscInt :: newton_iteration

  call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

  grid => this%realization%patch%grid
  patch => this%realization%patch
  option => this%realization%option
  gen_auxvars => this%realization%patch%aux%General%auxvars
  global_auxvars => this%realization%patch%aux%Global%auxvars

  changed = PETSC_TRUE

  ! MAN: OLD
  field => this%realization%field

  spid = option%saturation_pressure_id
  apid = option%air_pressure_id

  call SNESGetIterationNumber(snes,newton_iteration,ierr);CHKERRQ(ierr)

  ! MAN: END OLD
  if (this%check_post_convergence) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      offset = (local_id-1)*option%nflowdof
      select case(global_auxvars(ghosted_id)%istate)
        case(LIQUID_STATE)
          xmol_index = offset + GENERAL_LIQUID_STATE_X_MOLE_DOF
          pw_index = offset + GENERAL_LIQUID_PRESSURE_DOF
          if (X_p(xmol_index) - dX_p(xmol_index) < 0.d0) then
            dX_p(xmol_index) = X_p(xmol_index)
            changed = PETSC_TRUE
          endif
          if (X_p(pw_index)- dX_p(pw_index) <= 0.d0) then
           dX_p(pw_index) = X_p(pw_index) - ALMOST_ZERO
           changed = PETSC_TRUE
          endif
          smol_index = offset + GENERAL_LIQUID_STATE_S_MOLE_DOF !salt mole fraction
          if (X_p(smol_index) - dX_p(smol_index) < 0.d0) then
            dX_p(smol_index) = X_p(smol_index)
            changed = PETSC_TRUE
          endif
        case(GAS_STATE)
         pgas_index = offset + GENERAL_GAS_PRESSURE_DOF
         if (X_p(pgas_index)- dX_p(pgas_index) <= 0.d0) then
           dX_p(pgas_index) = X_p(pgas_index) - ALMOST_ZERO
           changed = PETSC_TRUE
         endif
        case(TWO_PHASE_STATE)
          pgas_index = offset + GENERAL_GAS_PRESSURE_DOF
          if (X_p(pgas_index) - dX_p(pgas_index) < &
                  gen_auxvars(ZERO_INTEGER,ghosted_id)% &
                  pres(option%saturation_pressure_id)) then
            dX_p(pgas_index) = X_p(pgas_index) - &
                    gen_auxvars(ZERO_INTEGER,ghosted_id)% &
                    pres(option%saturation_pressure_id)
            changed = PETSC_TRUE
          endif
          if (general_immiscible) then
            saturation_index = offset + GENERAL_GAS_SATURATION_DOF
            temp_real = X_p(saturation_index) - dX_p(saturation_index)
            if (temp_real > ALMOST_ONE) then
              dX_p(saturation_index) = X_p(saturation_index) - ALMOST_ONE
              changed = PETSC_TRUE
            else if (temp_real < ALMOST_ZERO) then
              dX_p(saturation_index) = X_p(saturation_index) - ALMOST_ZERO
              changed = PETSC_TRUE
            endif
          endif
        case(LGP_STATE)
          pgas_index = offset + GENERAL_GAS_PRESSURE_DOF
          if (X_p(pgas_index) - dX_p(pgas_index) < &
               gen_auxvars(ZERO_INTEGER,ghosted_id)% &
               pres(option%saturation_pressure_id)) then
             dX_p(pgas_index) = X_p(pgas_index) - &
                  gen_auxvars(ZERO_INTEGER,ghosted_id)% &
                  pres(option%saturation_pressure_id)
             changed = PETSC_TRUE
          endif
        case(LP_STATE)
          xmol_index = offset + GENERAL_LIQUID_STATE_X_MOLE_DOF
          pw_index = offset + GENERAL_LIQUID_PRESSURE_DOF
          if (X_p(xmol_index) - dX_p(xmol_index) < 0.d0) then
            dX_p(xmol_index) = X_p(xmol_index)
            changed = PETSC_TRUE
          endif
          if (X_p(pw_index)- dX_p(pw_index) <= 0.d0) then
           dX_p(pw_index) = X_p(pw_index) - ALMOST_ZERO
           changed = PETSC_TRUE
          endif
        case(GP_STATE)  
          pgas_index = offset + GENERAL_GAS_PRESSURE_DOF
          if (X_p(pgas_index)- dX_p(pgas_index) <= 0.d0) then
            dX_p(pgas_index) = X_p(pgas_index) - ALMOST_ZERO
            changed = PETSC_TRUE
          endif
      end select
    enddo

    if (this%damping_factor > 0.d0) then
      dX_p = dX_p*this%damping_factor
      changed = PETSC_TRUE
    endif

! MAN OLD
  else

    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      offset = (local_id-1)*option%nflowdof
      select case(global_auxvars(ghosted_id)%istate)
        case(LIQUID_STATE)
          xmol_index = offset + GENERAL_LIQUID_STATE_X_MOLE_DOF
          pw_index = offset + GENERAL_LIQUID_PRESSURE_DOF
          if (X_p(xmol_index) - dX_p(xmol_index) < 0.d0) then
            dX_p(xmol_index) = X_p(xmol_index)
            changed = PETSC_TRUE
          endif
        case(TWO_PHASE_STATE)
          pgas_index = offset + GENERAL_GAS_PRESSURE_DOF
          if (X_p(pgas_index) - dX_p(pgas_index) < &
                  gen_auxvars(ZERO_INTEGER,ghosted_id)% &
                  pres(option%saturation_pressure_id)) then
            dX_p(pgas_index) = X_p(pgas_index) - &
                    gen_auxvars(ZERO_INTEGER,ghosted_id)% &
                    pres(option%saturation_pressure_id)
            changed = PETSC_TRUE
          endif
          if (general_immiscible) then
            saturation_index = offset + GENERAL_GAS_SATURATION_DOF
            temp_real = X_p(saturation_index) - dX_p(saturation_index)
            if (temp_real > ALMOST_ONE) then
              dX_p(saturation_index) = X_p(saturation_index) - ALMOST_ONE
              changed = PETSC_TRUE
            else if (temp_real < ALMOST_ZERO) then
              dX_p(saturation_index) = X_p(saturation_index) - ALMOST_ZERO
              changed = PETSC_TRUE
            endif
          endif
      end select
    enddo

    scale = initial_scale
    if (general_max_it_before_damping > 0 .and. &
        newton_iteration > general_max_it_before_damping) then
      scale = general_damping_factor
    endif

#define LIMIT_MAX_PRESSURE_CHANGE
#define LIMIT_MAX_SATURATION_CHANGE
    ! scaling
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      offset = (local_id-1)*option%nflowdof
      temp_scale = 1.d0
      select case(global_auxvars(ghosted_id)%istate)
        case(LIQUID_STATE)
          liquid_pressure_index  = offset + GENERAL_LIQUID_PRESSURE_DOF
          temperature_index  = offset + GENERAL_ENERGY_DOF
          dX_p(liquid_pressure_index) = dX_p(liquid_pressure_index) * &
                                        GENERAL_PRESSURE_SCALE
          temp_scale = 1.d0
          del_liquid_pressure = dX_p(liquid_pressure_index)
          liquid_pressure0 = X_p(liquid_pressure_index)
          liquid_pressure1 = liquid_pressure0 - del_liquid_pressure
          del_temperature = dX_p(temperature_index)
          temperature0 = X_p(temperature_index)
          temperature1 = temperature0 - del_temperature
#ifdef LIMIT_MAX_PRESSURE_CHANGE
          if (dabs(del_liquid_pressure) > general_max_pressure_change) then
            temp_real = dabs(general_max_pressure_change/del_liquid_pressure)
            temp_scale = min(temp_scale,temp_real)
          endif
#endif
!LIMIT_MAX_PRESSURE_CHANGE
        case(TWO_PHASE_STATE)
          gas_pressure_index = offset + GENERAL_GAS_PRESSURE_DOF
!        air_pressure_index = offset + 2
          saturation_index = offset + GENERAL_GAS_SATURATION_DOF
          temperature_index  = offset + GENERAL_ENERGY_DOF
          dX_p(gas_pressure_index) = dX_p(gas_pressure_index) * &
                                     GENERAL_PRESSURE_SCALE
          if (general_2ph_energy_dof == GENERAL_AIR_PRESSURE_INDEX) then
            air_pressure_index = offset + GENERAL_ENERGY_DOF
            dX_p(air_pressure_index) = dX_p(air_pressure_index) * &
                                       GENERAL_PRESSURE_SCALE
            del_air_pressure = dX_p(air_pressure_index)
            air_pressure0 = X_p(air_pressure_index)
            air_pressure1 = air_pressure0 - del_air_pressure
          endif
          temp_scale = 1.d0
          del_gas_pressure = dX_p(gas_pressure_index)
          gas_pressure0 = X_p(gas_pressure_index)
          gas_pressure1 = gas_pressure0 - del_gas_pressure
          del_saturation = dX_p(saturation_index)
          saturation0 = X_p(saturation_index)
          saturation1 = saturation0 - del_saturation
#ifdef LIMIT_MAX_PRESSURE_CHANGE
          if (dabs(del_gas_pressure) > general_max_pressure_change) then
            temp_real = dabs(general_max_pressure_change/del_gas_pressure)
            temp_scale = min(temp_scale,temp_real)
          endif
#endif
#ifdef LIMIT_MAX_SATURATION_CHANGE
          if (dabs(del_saturation) > max_saturation_change) then
            temp_real = dabs(max_saturation_change/del_saturation)
            temp_scale = min(temp_scale,temp_real)
          endif
#endif
!LIMIT_MAX_SATURATION_CHANGE
        case(GAS_STATE)
          gas_pressure_index = offset + GENERAL_GAS_PRESSURE_DOF
          air_pressure_index = offset + GENERAL_GAS_STATE_AIR_PRESSURE_DOF
          dX_p(gas_pressure_index) = dX_p(gas_pressure_index) * &
                                     GENERAL_PRESSURE_SCALE
          if (general_gas_air_mass_dof == GENERAL_AIR_PRESSURE_INDEX) then
            dX_p(air_pressure_index) = dX_p(air_pressure_index) * &
                                      GENERAL_PRESSURE_SCALE
          endif
      end select
      scale = min(scale,temp_scale)
    enddo

    temp_scale = scale
    call MPI_Allreduce(temp_scale,scale,ONE_INTEGER_MPI,MPI_DOUBLE_PRECISION, &
                       MPI_MIN,option%mycomm,ierr);CHKERRQ(ierr)

    if (scale < 0.9999d0) then
      dX_p = scale*dX_p
    endif
  endif

  call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

end subroutine PMGeneralCheckUpdatePre

! ************************************************************************** !

subroutine PMGeneralCheckUpdatePost(this,snes,X0,dX,X1,dX_changed, &
                                    X1_changed,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 11/20/18
  !
  use General_Aux_module
  use Global_Aux_module
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module

  implicit none

  class(pm_general_type) :: this
  SNES :: snes
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr

  PetscReal, pointer :: X0_p(:)
  PetscReal, pointer :: dX_p(:)
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars(:)
  PetscInt :: local_id, ghosted_id, natural_id
  PetscInt :: offset, ival, idof
  PetscReal :: dX_abs, dX_X0
  PetscBool, pointer :: converged_abs_update_flag(:,:)
  PetscBool, pointer :: converged_rel_update_flag(:,:)
  PetscInt, pointer :: converged_abs_update_cell(:,:)
  PetscInt, pointer :: converged_rel_update_cell(:,:)
  PetscReal, pointer :: converged_abs_update_real(:,:)
  PetscReal, pointer :: converged_rel_update_real(:,:)
  PetscInt :: istate
  PetscBool :: converged_absolute
  PetscBool :: converged_relative

  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch
  global_auxvars => patch%aux%Global%auxvars

  allocate(converged_abs_update_flag(option%nflowdof,general_max_states))
  allocate(converged_rel_update_flag(option%nflowdof,general_max_states))
  allocate(converged_abs_update_cell(option%nflowdof,general_max_states))
  allocate(converged_rel_update_cell(option%nflowdof,general_max_states))
  allocate(converged_abs_update_real(option%nflowdof,general_max_states))
  allocate(converged_rel_update_real(option%nflowdof,general_max_states))

  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE

  call VecGetArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)
  converged_abs_update_flag = PETSC_TRUE
  converged_rel_update_flag = PETSC_TRUE
  converged_abs_update_cell = ZERO_INTEGER
  converged_rel_update_cell = ZERO_INTEGER
  converged_abs_update_real = 0.d0
  converged_rel_update_real = 0.d0

  do local_id = 1, grid%nlmax
    offset = (local_id-1)*option%nflowdof
    ghosted_id = grid%nL2G(local_id)
    natural_id = grid%nG2A(ghosted_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    istate = global_auxvars(ghosted_id)%istate

    do idof = 1, option%nflowdof

      ival = offset+idof

      ! infinity norms on update
      converged_absolute = PETSC_TRUE
      converged_relative = PETSC_TRUE
      dX_abs = dabs(dX_p(ival))
      if (X0_p(ival) > 0.d0) then
        dX_X0 = dabs(dX_abs/(X0_p(ival)))
      else
        dX_X0 = dabs(dX_abs/1.d-40)
      endif

      if (dX_abs > this%abs_update_inf_tol(idof,istate)) then
        converged_absolute = PETSC_FALSE
      endif
      if (converged_abs_update_real(idof,istate) < dX_abs) then
        converged_abs_update_real(idof,istate) = dX_abs
        converged_abs_update_cell(idof,istate) = natural_id
      endif
      if (dX_X0 > this%rel_update_inf_tol(idof,istate)) then
        converged_relative = PETSC_FALSE
      endif
      if (converged_rel_update_real(idof,istate) < dX_X0) then
        converged_rel_update_real(idof,istate) = dX_X0
        converged_rel_update_cell(idof,istate) = natural_id
      endif

      ! only enter this condition if both are not converged
      if (.not.(converged_absolute .or. converged_relative)) then
        converged_abs_update_flag(idof,istate) = PETSC_FALSE
        converged_rel_update_flag(idof,istate) = PETSC_FALSE
      endif
    enddo
  enddo
  call VecRestoreArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)

  this%converged_flag(:,:,ABS_UPDATE_INDEX) = converged_abs_update_flag(:,:)
  this%converged_flag(:,:,REL_UPDATE_INDEX) = converged_rel_update_flag(:,:)
  this%converged_real(:,:,ABS_UPDATE_INDEX) = converged_abs_update_real(:,:)
  this%converged_real(:,:,REL_UPDATE_INDEX) = converged_rel_update_real(:,:)
  this%converged_cell(:,:,ABS_UPDATE_INDEX) = converged_abs_update_cell(:,:)
  this%converged_cell(:,:,REL_UPDATE_INDEX) = converged_rel_update_cell(:,:)

end subroutine PMGeneralCheckUpdatePost


! ************************************************************************** !

subroutine PMGeneralCheckConvergence(this,snes,it,xnorm,unorm,fnorm, &
                                     reason,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 11/20/18
  !
  use Convergence_module
  use General_Aux_module
  use Global_Aux_module
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  use String_module

  implicit none

  class(pm_general_type) :: this
  SNES :: snes
  PetscInt :: it
  PetscReal :: xnorm
  PetscReal :: unorm
  PetscReal :: fnorm
  SNESConvergedReason :: reason
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars(:)
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum2_p(:)
  PetscInt :: local_id, ghosted_id, natural_id
  PetscInt :: offset, ival, idof, itol
  PetscReal :: R, A, R_A
  PetscReal, parameter :: A_zero = 1.d-15
  PetscBool, allocatable :: converged_abs_residual_flag(:,:)
  PetscReal, allocatable :: converged_abs_residual_real(:,:)
  PetscInt, allocatable :: converged_abs_residual_cell(:,:)
  PetscBool, allocatable :: converged_scaled_residual_flag(:,:)
  PetscReal, allocatable :: converged_scaled_residual_real(:,:)
  PetscInt, allocatable :: converged_scaled_residual_cell(:,:)
  PetscInt :: istate
  PetscBool :: converged_absolute
  PetscBool :: converged_scaled
  PetscMPIInt :: mpi_int
  PetscBool, allocatable :: flags(:)
  PetscBool :: rho_flag
  character(len=MAXSTRINGLENGTH) :: string
  character(len=12), allocatable :: state_string(:)
  character(len=17), allocatable :: dof_string(:,:)
  character(len=15), parameter :: tol_string(MAX_INDEX) = &
    ['Absolute Update','Relative Update','Residual       ','Scaled Residual']

  patch => this%realization%patch
  option => this%realization%option
  field => this%realization%field
  grid => patch%grid
  global_auxvars => patch%aux%Global%auxvars

  allocate(flags(MAX_INDEX*option%nflowdof*general_max_states+1))
  allocate(state_string(general_max_states))
  allocate(dof_string(option%nflowdof,general_max_states))

  allocate(converged_abs_residual_flag(option%nflowdof,general_max_states))
  allocate(converged_abs_residual_real(option%nflowdof,general_max_states))
  allocate(converged_abs_residual_cell(option%nflowdof,general_max_states))
  allocate(converged_scaled_residual_flag(option%nflowdof,general_max_states))
  allocate(converged_scaled_residual_real(option%nflowdof,general_max_states))
  allocate(converged_scaled_residual_cell(option%nflowdof,general_max_states))

  if (option%nflowdof == 3) then
    state_string = &
       ['Liquid State','Gas State   ','2Phase State']
    dof_string = &
      reshape(['Liquid Pressure  ','Air Mole Fraction','Temperature      ', &
               'Gas Pressure     ','Air Pressure     ','Temperature      ', &
               'Gas Pressure     ','Gas Saturation   ','Temperature      '], &
               shape(dof_string))
  elseif (option%nflowdof == 4) then
    state_string = &
      ['Liquid State','Gas State   ','LG State    ',&
       'Precip State','LP State    ','GP State    ',&
       'LGP State   ']
    dof_string = &
    reshape(['Liquid Pressure  ','Air Mole Fraction','Temperature      ',&
             'Salt Fraction    ',&
             'Gas Pressure     ','Air Mole Fraction','Temperature      ',&
             'Salt Fraction    ',&
             'Gas Pressure     ','Gas Saturation   ','Temperature      ',&
             'Salt Fraction    ',&
             'Liquid Pressure  ','Air Mole Fraction','Temperature      ',&
             'Salt Fraction    ',&
             'Liquid Pressure  ','Air Mole Fraction','Temperature      ',&
             'Salt Fraction    ',&
             'Gas Pressure     ','Air Pressure     ','Temperature      ',&
             'Precipitate Sat. ',&
             'Gas Pressure     ','Gas Saturation   ','Temperature      ',&
             'Precipitate Sat. '],shape(dof_string))
  endif
  call SNESNewtonTRDCGetRhoFlag(snes,rho_flag,ierr);CHKERRQ(ierr);

  if (this%option%flow%using_newtontrdc) then
    if (general_newtontrdc_prev_iter_num == it) then
      general_sub_newton_iter_num = general_sub_newton_iter_num + 1
    endif
    general_newtontrdc_prev_iter_num = it
  endif
  if (this%check_post_convergence) then
    call VecGetArrayReadF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(field%flow_accum2,accum2_p,ierr);CHKERRQ(ierr)
    converged_abs_residual_flag = PETSC_TRUE
    converged_abs_residual_real = 0.d0
    converged_abs_residual_cell = ZERO_INTEGER
    converged_scaled_residual_flag = PETSC_TRUE
    converged_scaled_residual_real = 0.d0
    converged_scaled_residual_cell = ZERO_INTEGER
    do local_id = 1, grid%nlmax
      offset = (local_id-1)*option%nflowdof
      ghosted_id = grid%nL2G(local_id)
      natural_id = grid%nG2A(ghosted_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      istate = global_auxvars(ghosted_id)%istate
      do idof = 1, option%nflowdof
        ival = offset+idof
        converged_absolute = PETSC_TRUE
        converged_scaled = PETSC_TRUE
        ! infinity norms on residual
        R = dabs(r_p(ival))
        A = dabs(accum2_p(ival))
!         R_A = R/A

        !TOUGH3 way:
        if (A > 1.d0) then
          R_A = R/A
        else
          R_A = R
        endif

        if (R > this%residual_abs_inf_tol(idof)) then
          converged_absolute = PETSC_FALSE
        endif

        ! find max value regardless of convergence
        if (converged_abs_residual_real(idof,istate) < R) then
          converged_abs_residual_real(idof,istate) = R
          converged_abs_residual_cell(idof,istate) = natural_id
        endif
        if (A > A_zero) then
          if (R_A > this%residual_scaled_inf_tol(idof)) then
            converged_scaled = PETSC_FALSE
          endif
          ! find max value regardless of convergence
          if (converged_scaled_residual_real(idof,istate) < R_A) then
            converged_scaled_residual_real(idof,istate) = R_A
            converged_scaled_residual_cell(idof,istate) = natural_id
          endif
        endif
        ! only enter this condition if both are not converged
        if (.not.(converged_absolute .or. converged_scaled)) then
          converged_abs_residual_flag(idof,istate) = PETSC_FALSE
          converged_scaled_residual_flag(idof,istate) = PETSC_FALSE
        endif
      enddo
    enddo
    call VecRestoreArrayReadF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(field%flow_accum2,accum2_p, &
                                ierr);CHKERRQ(ierr)

    this%converged_flag(:,:,RESIDUAL_INDEX) = converged_abs_residual_flag(:,:)
    this%converged_real(:,:,RESIDUAL_INDEX) = converged_abs_residual_real(:,:)
    this%converged_cell(:,:,RESIDUAL_INDEX) = converged_abs_residual_cell(:,:)
    this%converged_flag(:,:,SCALED_RESIDUAL_INDEX) = &
                                       converged_scaled_residual_flag(:,:)
    this%converged_real(:,:,SCALED_RESIDUAL_INDEX) = &
                                       converged_scaled_residual_real(:,:)
    this%converged_cell(:,:,SCALED_RESIDUAL_INDEX) = &
                                       converged_scaled_residual_cell(:,:)
    ! do not perform an all reduce on cell id as this info is not printed
    ! in parallel

    ! geh: since we need to pack other flags into this global reduction,
    !      convert to 1D array
    flags(1:option%nflowdof*general_max_states*MAX_INDEX) = &
      reshape(this%converged_flag,(/option%nflowdof*general_max_states*MAX_INDEX/))
    ! due to the 'and' operation, must invert the boolean using .not.
    flags(option%nflowdof*general_max_states*MAX_INDEX+1) =&
      .not.general_high_temp_ts_cut
    mpi_int = option%nflowdof*general_max_states*MAX_INDEX+1
    call MPI_Allreduce(MPI_IN_PLACE,flags,mpi_int, &
                       MPI_LOGICAL,MPI_LAND,option%mycomm,ierr);CHKERRQ(ierr)

    this%converged_flag = reshape(flags(1:option%nflowdof*general_max_states*MAX_INDEX),&
                            (/option%nflowdof,general_max_states,MAX_INDEX/))
    ! due to the 'and' operation, must invert the boolean using .not.
    general_high_temp_ts_cut = .not.flags(option%nflowdof*general_max_states*MAX_INDEX+1)

    mpi_int = option%nflowdof*general_max_states*MAX_INDEX
    call MPI_Allreduce(MPI_IN_PLACE,this%converged_real,mpi_int, &
                       MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)

    option%convergence = CONVERGENCE_CONVERGED
    do itol = 1, MAX_INDEX
      do istate = 1, general_max_states
        do idof = 1, option%nflowdof
          if (.not.this%converged_flag(idof,istate,itol)) then
            option%convergence = CONVERGENCE_KEEP_ITERATING
            if (this%logging_verbosity > 0) then
              if (trim(tol_string(itol)) == 'Residual' .or. &
                  trim(tol_string(itol)) == 'Scaled Residual') then
                if (idof == 1) then
                  string = '   ' // trim(tol_string(itol)) // ', ' // &
                   trim(state_string(istate)) // ', Water Mass'
                elseif (idof == 2) then
                  string = '   ' // trim(tol_string(itol)) // ', ' // &
                   trim(state_string(istate)) // ', Air Mass'
                elseif (idof == 3) then
                  string = '   ' // trim(tol_string(itol)) // ', ' // &
                   trim(state_string(istate)) // ', Energy'
                elseif (idof == 4) then
                  string = '   ' // trim(tol_string(itol)) // ', ' // &
                   trim(state_string(istate)) // ', Salt Mass'
                endif
              else
                string = '   ' // trim(tol_string(itol)) // ', ' // &
                 trim(state_string(istate)) // ', ' // dof_string(idof,istate)
              endif
              if (option%comm%size == 1) then
                string = trim(string) // ' (' // &
                  trim(StringFormatInt(this%converged_cell(idof,istate,itol))) &
                  // ')'
              endif
              string = trim(string) // ' : ' // &
                StringFormatDouble(this%converged_real(idof,istate,itol))
              call PrintMsg(option,string)
            endif
          endif
        enddo
      enddo
    enddo

    if (option%flow%using_newtontrdc .and. &
        general_state_changed .and. &
        .not.rho_flag) then
      if (general_newtontrdc_hold_inner) then
        ! if we hold inner iterations, we must not change state in
        ! the inner iteration. If we reach convergence in an inner
        ! newtontrdc iteration, then we must force an outer iteration
        ! to allow state change in case the solutions are
        ! out-of-bounds of the states -hdp
        general_force_iteration = PETSC_TRUE
        general_state_changed = PETSC_FALSE
      else
        ! if we have state changes, we exit out of inner iteration
        ! and go to the next newton iteration. the tr inner iteration
        !  should only be used when there is no state changes
        ! if rho is satisfied in inner iteration, the algorithm already
        ! exited the inner iteration. -heeho
        general_force_iteration = PETSC_TRUE
        general_state_changed = PETSC_FALSE
      endif
    endif

    call MPI_Allreduce(MPI_IN_PLACE,general_force_iteration,ONE_INTEGER, &
                       MPI_LOGICAL,MPI_LOR,option%mycomm,ierr)
    if (general_force_iteration) then
      if (.not.general_newtontrdc_hold_inner) then
        option%convergence = CONVERGENCE_BREAKOUT_INNER_ITER
        general_force_iteration = PETSC_FALSE
      else if (general_newtontrdc_hold_inner .and. &
               option%convergence == CONVERGENCE_CONVERGED) then
        option%convergence = CONVERGENCE_BREAKOUT_INNER_ITER
        general_force_iteration = PETSC_FALSE
      endif
    endif

    if (this%logging_verbosity > 0 .and. it > 0 .and. &
        option%convergence == CONVERGENCE_CONVERGED) then
      string = '   Converged'
      call PrintMsg(option,string)
      write(string,'(4x," R:",9es8.1)') this%converged_real(:,:,RESIDUAL_INDEX)
      call PrintMsg(option,string)
      write(string,'(4x,"SR:",9es8.1)') &
        this%converged_real(:,:,SCALED_RESIDUAL_INDEX)
      call PrintMsg(option,string)
      write(string,'(4x,"AU:",9es8.1)') &
        this%converged_real(:,:,ABS_UPDATE_INDEX)
      call PrintMsg(option,string)
      write(string,'(4x,"RU:",9es8.1)') &
        this%converged_real(:,:,REL_UPDATE_INDEX)
      call PrintMsg(option,string)
    endif

    if (it >= this%solver%newton_max_iterations) then
      option%convergence = CONVERGENCE_CUT_TIMESTEP

      if (this%logging_verbosity > 0) then
        string = '    Exceeded General Mode Max Newton Iterations'
        call PrintMsg(option,string)
      endif
    endif
    if (general_high_temp_ts_cut) then
      general_high_temp_ts_cut = PETSC_FALSE
      string = '    Exceeded General Mode EOS max temperature'
      call PrintMsg(option,string)
      option%convergence = CONVERGENCE_CUT_TIMESTEP
    endif
    if (general_sub_newton_iter_num > 20) then
      ! cut time step in case PETSC solvers are missing inner iterations
      option%convergence = CONVERGENCE_CUT_TIMESTEP
    endif

  endif

  call PMSubsurfaceFlowCheckConvergence(this,snes,it,xnorm,unorm,fnorm, &
                                        reason,ierr)

  if (this%check_post_convergence) then
    !RESET ALL FLAGS AND ZERO ALL THE CELL ID'S
    this%converged_flag(:,:,:) = PETSC_TRUE
    this%converged_real(:,:,:) = 0.d0
    this%converged_cell(:,:,:) = 0
  endif

end subroutine PMGeneralCheckConvergence

! ************************************************************************** !

subroutine PMGeneralTimeCut(this)
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  use General_module, only : GeneralTimeCut

  implicit none

  class(pm_general_type) :: this

  call PMSubsurfaceFlowTimeCut(this)
  call GeneralTimeCut(this%realization)

end subroutine PMGeneralTimeCut

! ************************************************************************** !

subroutine PMGeneralUpdateSolution(this)
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  use General_module, only : GeneralUpdateSolution, &
                             GeneralMapBCAuxVarsToGlobal

  implicit none

  class(pm_general_type) :: this

  call PMSubsurfaceFlowUpdateSolution(this)
  call GeneralUpdateSolution(this%realization)
  call GeneralMapBCAuxVarsToGlobal(this%realization)

end subroutine PMGeneralUpdateSolution

! ************************************************************************** !

subroutine PMGeneralUpdateAuxVars(this)
  !
  ! Author: Glenn Hammond
  ! Date: 04/21/14
  use General_module, only : GeneralUpdateAuxVars

  implicit none

  class(pm_general_type) :: this

                                                !update BCs (second PETSC_TRUE)
  call GeneralUpdateAuxVars(this%realization,PETSC_FALSE,PETSC_TRUE)

end subroutine PMGeneralUpdateAuxVars

! ************************************************************************** !

subroutine PMGeneralMaxChange(this)
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  use Realization_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use Global_Aux_module
  use General_Aux_module

  implicit none

  class(pm_general_type) :: this

  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscReal, pointer :: max_change_local(:)
  PetscReal, pointer :: max_change_global(:)
  PetscReal :: max_change
  PetscInt :: i, j, max_change_index
  PetscInt :: ghosted_id



  PetscErrorCode :: ierr

  realization => this%realization
  option => realization%option
  field => realization%field
  grid => realization%patch%grid

  global_auxvars => realization%patch%aux%global%auxvars

  if (option%nflowdof == 3) then
    max_change_index = SIX_INTEGER
  elseif (option%nflowdof == 4) then
    max_change_index = EIGHT_INTEGER
  endif
  allocate(max_change_local(max_change_index))
  allocate(max_change_global(max_change_index))
  max_change_global = 0.d0
  max_change_local = 0.d0

  do i = 1, max_change_index
    call RealizationGetVariable(realization,field%work, &
                                this%max_change_ivar(i), &
                                this%max_change_isubvar(i))
    ! yes, we could use VecWAXPY and a norm here, but we need the ability
    ! to customize
    call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%max_change_vecs(i),vec_ptr2,ierr);CHKERRQ(ierr)
    max_change = 0.d0
    if (i==1 .and. gen_chk_max_dpl_liq_state_only) then
      do j = 1,grid%nlmax
        ghosted_id = grid%nL2G(j)
        if (global_auxvars(ghosted_id)%istate /= LIQUID_STATE .or. &
            global_auxvars(ghosted_id)%istate /= LP_STATE) then
          vec_ptr(j) = 0.d0
        endif
      enddo
    endif

    do j = 1, grid%nlmax
      ! have to weed out cells that changed state
      if (dabs(vec_ptr(j)) > 1.d-40 .and. dabs(vec_ptr2(j)) > 1.d-40) then
        max_change = max(max_change,dabs(vec_ptr(j)-vec_ptr2(j)))
      endif
    enddo

    max_change_local(i) = max_change
    call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%max_change_vecs(i),vec_ptr2, &
                            ierr);CHKERRQ(ierr)
    call VecCopy(field%work,field%max_change_vecs(i),ierr);CHKERRQ(ierr)
  enddo
  call MPI_Allreduce(max_change_local,max_change_global,max_change_index, &
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr);CHKERRQ(ierr)
  ! print them out
  if (OptionPrintToScreen(option)) then
    if (option%nflowdof == 3) then
     write(*,'("  --> max chng: dpl= ",1pe12.4, " dpg= ",1pe12.4,&
        & " dpa= ",1pe12.4,/,15x," dxa= ",1pe12.4,"  dt= ",1pe12.4,&
        & " dsg= ",1pe12.4)') &
        max_change_global(1:max_change_index)
    elseif (option%nflowdof == 4) then
      write(*,'("  --> max chng: dpl= ",1pe12.4, " dpg= ",1pe12.4,&
        & " dpa= ",1pe12.4,/,15x," dxa= ",1pe12.4," dxs= ",1pe12.4,&
        & "  dt= ",1pe12.4,/,15x," dsg= ",1pe12.4," dpo= ",1pe12.4)') &
        max_change_global(1:max_change_index)
    endif
  endif
  if (OptionPrintToFile(option)) then
    if (option%nflowdof == 3) then
      write(option%fid_out,'("  --> max chng: dpl= ",1pe12.4, " dpg= ",1pe12.4,&
        & " dpa= ",1pe12.4,/,15x," dxa= ",1pe12.4,"  dt= ",1pe12.4, &
        & " dsg= ",1pe12.4)') &
        max_change_global(1:max_change_index)
    elseif (option%nflowdof == 4) then
      write(option%fid_out,'("  --> max chng: dpl= ",1pe12.4, " dpg= ",1pe12.4,&
        & " dpa= ",1pe12.4,/,15x," dxa= ",1pe12.4,"  dt= ",1pe12.4,&
        & " dsg= ",1pe12.4,/,15x," dpo= ",1pe12.4)') &
        max_change_global(1:max_change_index)
    endif
  endif

  this%max_pressure_change = maxval(max_change_global(1:3))
  this%max_temperature_change = max_change_global(5)
  if (option%nflowdof == 3) then
    this%max_saturation_change = max_change_global(6)
    this%max_xmol_change = max_change_global(4)
  elseif (option%nflowdof == 4) then
    this%max_xmol_change = maxval(max_change_global(4:5))
    this%max_saturation_change = maxval(max_change_global(6:7))
  endif

end subroutine PMGeneralMaxChange
! ************************************************************************** !

subroutine PMGeneralComputeMassBalance(this,mass_balance_array)
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  use General_module, only : GeneralComputeMassBalance

  implicit none

  class(pm_general_type) :: this
  PetscReal :: mass_balance_array(:)

  call GeneralComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMGeneralComputeMassBalance

! ************************************************************************** !

subroutine PMGeneralInputRecord(this)
  !
  ! Writes ingested information to the input record file.
  !
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  !

  implicit none

  class(pm_general_type) :: this

  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'general'
  if (this%check_post_convergence) then
    write(id,'(a29)',advance='no') 'ITOL_SCALED_RESIDUAL: '
    write(id,'(a)') 'ON'
    write(id,'(a29)',advance='no') 'ITOL_RELATIVE_RESIDUAL: '
    write(id,'(a)') 'ON'
  endif

end subroutine PMGeneralInputRecord

! ************************************************************************** !

subroutine PMGeneralCheckpointBinary(this,viewer)
  !
  ! Checkpoints data associated with General PM
  !
  ! Author: Glenn Hammond
  ! Date: 02/18/15

  use Checkpoint_module
  use Global_module

  implicit none
#include "petsc/finclude/petscviewer.h"

  class(pm_general_type) :: this
  PetscViewer :: viewer

  call PMSubsurfaceFlowCheckpointBinary(this,viewer)

end subroutine PMGeneralCheckpointBinary

! ************************************************************************** !

subroutine PMGeneralRestartBinary(this,viewer)
  !
  ! Restarts data associated with General PM
  !
  ! Author: Glenn Hammond
  ! Date: 02/18/15

  use Checkpoint_module
  use Global_module

  implicit none
#include "petsc/finclude/petscviewer.h"

  class(pm_general_type) :: this
  PetscViewer :: viewer

  call PMSubsurfaceFlowRestartBinary(this,viewer)

end subroutine PMGeneralRestartBinary

! ************************************************************************** !

subroutine PMGeneralDestroy(this)
  !
  ! Destroys General process model
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  use General_module, only : GeneralDestroy

  implicit none

  class(pm_general_type) :: this

  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  deallocate(this%max_change_ivar)
  nullify(this%max_change_ivar)
  deallocate(this%max_change_isubvar)
  nullify(this%max_change_isubvar)

  ! preserve this ordering
  call GeneralDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)

end subroutine PMGeneralDestroy

end module PM_General_class
