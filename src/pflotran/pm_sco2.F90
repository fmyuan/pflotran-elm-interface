module PM_SCO2_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class

  use PFLOTRAN_Constants_module

  private

  PetscBool :: sco2_use_governors = PETSC_FALSE
  PetscBool :: sco2_check_updates = PETSC_FALSE
  PetscBool :: sco2_stomp_convergence = PETSC_TRUE

  PetscInt, parameter :: ABS_UPDATE_INDEX = 1
  PetscInt, parameter :: REL_UPDATE_INDEX = 2
  PetscInt, parameter :: RESIDUAL_INDEX = 3
  PetscInt, parameter :: SCALED_RESIDUAL_INDEX = 4
  PetscInt, parameter :: MAX_INDEX = SCALED_RESIDUAL_INDEX
  PetscInt, parameter :: sco2_max_states = 4
  PetscInt, parameter :: max_change_index = 6 ! & with energy (temperature)

  ! DOF's: water mass, CO2 mass, salt mass fraction
  PetscInt, parameter :: MAX_DOF = 3
  ! DOF's: water mass, CO2 mass, salt mass fraction, energy
  ! PetscInt, parameter :: MAX_DOF = 4
  ! States: Liquid, Gas, Trapped Gas, Liquid-Gas
  ! MAN: is this redundant? Already have MAX_STATE in aux
  PetscInt, parameter :: MAX_STATE = 4 ! 7

  type, public, extends(pm_subsurface_flow_type) :: pm_sco2_type
    PetscInt, pointer :: max_change_ivar(:)
    PetscInt, pointer :: max_change_isubvar(:)
    PetscBool :: converged_flag(MAX_DOF,MAX_STATE,MAX_INDEX)
    PetscInt :: converged_cell(MAX_DOF,MAX_STATE,MAX_INDEX)
    PetscReal :: converged_real(MAX_DOF,MAX_STATE,MAX_INDEX)
    PetscReal :: residual_abs_inf_tol(MAX_DOF)
    PetscReal :: residual_scaled_inf_tol(MAX_DOF)
    PetscReal :: abs_update_inf_tol(MAX_DOF,MAX_STATE)
    PetscReal :: rel_update_inf_tol(MAX_DOF,MAX_STATE)
    PetscReal :: damping_factor
  contains
    procedure, public :: ReadSimulationOptionsBlock => &
                           PMSCO2ReadSimOptionsBlock
    procedure, public :: ReadNewtonBlock => PMSCO2ReadNewtonSelectCase
    procedure, public :: InitializeSolver => PMSCO2InitializeSolver
    procedure, public :: InitializeRun => PMSCO2InitializeRun
    procedure, public :: InitializeTimestep => PMSCO2InitializeTimestep
    procedure, public :: Residual => PMSCO2Residual
    procedure, public :: Jacobian => PMSCO2Jacobian
    procedure, public :: UpdateTimestep => PMSCO2UpdateTimestep
    procedure, public :: PreSolve => PMSCO2PreSolve
    procedure, public :: PostSolve => PMSCO2PostSolve
    procedure, public :: CheckUpdatePre => PMSCO2CheckUpdatePre
    procedure, public :: CheckUpdatePost => PMSCO2CheckUpdatePost
    procedure, public :: CheckConvergence => PMSCO2CheckConvergence
    procedure, public :: TimeCut => PMSCO2TimeCut
    procedure, public :: UpdateSolution => PMSCO2UpdateSolution
    procedure, public :: UpdateAuxVars => PMSCO2UpdateAuxVars
    procedure, public :: MaxChange => PMSCO2MaxChange
    procedure, public :: ComputeMassBalance => PMSCO2ComputeMassBalance
    procedure, public :: InputRecord => PMSCO2InputRecord
    procedure, public :: CheckpointBinary => PMSCO2CheckpointBinary
    procedure, public :: RestartBinary => PMSCO2RestartBinary
    procedure, public :: Destroy => PMSCO2Destroy

  end type pm_sco2_type

  public :: PMSCO2Create, &
            PMSCO2SetFlowMode

contains

! ************************************************************************** !

function PMSCO2Create()
  !
  ! Creates SCO2 mode process models
  !
  ! Author: Michael Nole
  ! Date: 11/20/2023
  !

  use Upwind_Direction_module
  use Option_module

  implicit none

  class(pm_sco2_type), pointer :: PMSCO2Create

  class(pm_sco2_type), pointer :: this

#ifdef PM_SCO2_DEBUG
  print *, 'PMSCO2Create()'
#endif

  allocate(this)
  call PMSubsurfaceFlowInit(this)
  this%name = 'SCO2 Flow'
  this%header = 'SCO2 FLOW'

  ! turn off default upwinding which is set to PETSC_TRUE in
  !  upwind_direction.F90
  fix_upwind_direction = PETSC_FALSE

  PMSCO2Create => this

end function PMSCO2Create

! ************************************************************************** !

subroutine PMSCO2SetFlowMode(pm,option)
  !
  ! Sets the flow mode to SCO2 mode
  !
  ! Author: Michael Nole
  ! Date: 11/17/2023
  !
  
  use Option_module
  use Variables_module, only: LIQUID_PRESSURE, GAS_PRESSURE, CO2_PRESSURE, &
                              LIQUID_MASS_FRACTION, LIQUID_SALT_MASS_FRAC, &
                              TEMPERATURE, LIQUID_SATURATION, GAS_SATURATION, &
                              PRECIPITATE_SATURATION, POROSITY

  implicit none

  class(pm_sco2_type) :: pm
  type(option_type) :: option

  PetscReal, parameter :: pres_abs_inf_tol = 1.d0 ! Reference tolerance [Pa]
  ! PetscReal, parameter :: temp_abs_inf_tol = 1.d-5 ! [C]
  PetscReal, parameter :: sat_abs_inf_tol = 1.d-5 ! [-]
  PetscReal, parameter :: xmass_abs_inf_tol = 1.d-5 ! [mass fraction]

  PetscReal, parameter :: pres_rel_inf_tol = 1.d-3
  ! PetscReal, parameter :: temp_rel_inf_tol = 1.d-3
  PetscReal, parameter :: sat_rel_inf_tol = 1.d-3
  PetscReal, parameter :: xmass_rel_inf_tol = 1.d-3

  PetscReal, parameter :: w_mass_abs_inf_tol = 1.d-5 !1.d-7 !kg_water/sec
  PetscReal, parameter :: co2_mass_abs_inf_tol = 1.d-5 !1.d-7 !kg_co2/sec
  ! PetscReal, parameter :: u_abs_inf_tol = 1.d-5 !1.d-7 !MW
  PetscReal, parameter :: s_mass_abs_inf_tol = 1.d-5 !1.d-7 !kg_salt/sec

  !With Energy
  ! PetscReal, parameter :: residual_abs_inf_tol(MAX_DOF) = &
  !                            (/w_mass_abs_inf_tol, co2_mass_abs_inf_tol, &
  !                              u_abs_inf_tol, s_mass_abs_inf_tol/)
  !Without Energy
  PetscReal, parameter :: residual_abs_inf_tol(MAX_DOF) = &
                             (/w_mass_abs_inf_tol, co2_mass_abs_inf_tol, &
                               s_mass_abs_inf_tol/)
  PetscReal, parameter :: residual_scaled_inf_tol(MAX_DOF) = 1.d-6

  ! With Energy
  ! PetscReal, parameter :: abs_update_inf_tol(MAX_DOF,MAX_STATE) = &
  !                         ! Liquid State: Pl, Xw_co2, T, X_salt
  !                         reshape([pres_abs_inf_tol,xmass_abs_inf_tol, &
  !                                  temp_abs_inf_tol,xmass_abs_inf_tol, &
  !                         ! Gas: Pg, Pco2, T, X_salt
  !                                  pres_abs_inf_tol,pres_abs_inf_tol, &
  !                                  temp_abs_inf_tol,xmass_abs_inf_tol,&
  !                         ! Trapped Gas: Pl, Sg, T, X_salt
  !                                  pres_abs_inf_tol,sat_abs_inf_tol, &
  !                                  temp_abs_inf_tol,xmass_abs_inf_tol, &
  !                         ! Liquid-Gas State: Pl, Pg, T, X_salt
  !                                  pres_abs_inf_tol,pres_abs_inf_tol, &
  !                                  temp_abs_inf_tol,xmass_abs_inf_tol], &
  !                                  shape(abs_update_inf_tol)) * &
  !                                  1.d0
  ! PetscReal, parameter :: rel_update_inf_tol(MAX_DOF,MAX_STATE) = &
  !                         ! Liquid State: Pl, Xw_co2, T, X_salt
  !                         reshape([pres_rel_inf_tol,xmass_rel_inf_tol, &
  !                                  temp_rel_inf_tol,xmass_rel_inf_tol,&
  !                         ! Gas State: Pg, Pco2, T, X_salt
  !                                  pres_rel_inf_tol,pres_rel_inf_tol, &
  !                                  temp_rel_inf_tol,xmass_rel_inf_tol,&
  !                         ! Trapped Gas State: Pl, Sg, T, X_salt
  !                                  pres_rel_inf_tol,sat_rel_inf_tol, &
  !                                  temp_rel_inf_tol,xmass_rel_inf_tol, &
  !                         ! Liquid-Gas State: Pl, Pg, T, X_salt
  !                                  pres_rel_inf_tol,pres_rel_inf_tol, &
  !                                  temp_rel_inf_tol,xmass_rel_inf_tol],&
  !                                  shape(rel_update_inf_tol)) * &
  !                                  1.d0
  
  !Without Energy
  PetscReal, parameter :: abs_update_inf_tol(MAX_DOF,MAX_STATE) = &
                          ! Liquid State: Pl, Xw_co2, T, X_salt
                          reshape([pres_abs_inf_tol,xmass_abs_inf_tol, &
                                   xmass_abs_inf_tol, &
                          ! Gas: Pg, Pco2, T, X_salt
                                   pres_abs_inf_tol,pres_abs_inf_tol, &
                                   xmass_abs_inf_tol,&
                          ! Trapped Gas: Pl, Sg, T, X_salt
                                   pres_abs_inf_tol,sat_abs_inf_tol, &
                                   xmass_abs_inf_tol, &
                          ! Liquid-Gas State: Pl, Pg, T, X_salt
                                   pres_abs_inf_tol,pres_abs_inf_tol, &
                                   xmass_abs_inf_tol], &
                                   shape(abs_update_inf_tol)) * &
                                   1.d0
  PetscReal, parameter :: rel_update_inf_tol(MAX_DOF,MAX_STATE) = &
                          ! Liquid State: Pl, Xw_co2, T, X_salt
                          reshape([pres_rel_inf_tol,xmass_rel_inf_tol, &
                                   xmass_rel_inf_tol,&
                          ! Gas State: Pg, Pco2, T, X_salt
                                   pres_rel_inf_tol,pres_rel_inf_tol, &
                                   xmass_rel_inf_tol,&
                          ! Trapped Gas State: Pl, Sg, T, X_salt
                                   pres_rel_inf_tol,sat_rel_inf_tol, &
                                   xmass_rel_inf_tol, &
                          ! Liquid-Gas State: Pl, Pg, T, X_salt
                                   pres_rel_inf_tol,pres_rel_inf_tol, &
                                   xmass_rel_inf_tol],&
                                   shape(rel_update_inf_tol)) * &
                                   1.d0

  option%iflowmode = SCO2_MODE

  ! liquid, gas, precipitate composite properties
  option%nphase = 3
  option%liquid_phase = 1
  option%gas_phase = 2
  option%precipitate_phase = 3

  option%capillary_pressure_id = 3
  option%co2_pressure_id = 4
  option%vapor_pressure_id = 5
  option%saturation_pressure_id = 6
  option%reduced_vapor_pressure_id = 7

  ! Extras: individual component properties, trapped gas phase
  option%pure_water_phase = 5
  option%pure_brine_phase = 6
  
  option%trapped_gas_phase = 4

  ! Water balance, CO2 balance, Salt balance
  option%nflowdof = 3
  ! Water balance, CO2 balance, Salt balance, Energy balance
  ! option%nflowdof = 4
  ! Components: water, co2, salt
  option%nflowspec = 3
  option%water_id = 1
  option%co2_id = 2
  option%air_id = 2
  option%salt_id = 3
  ! option%eid = 3
  
  

  allocate(pm%max_change_ivar(7))
  ! MAN: need to check this
  allocate(pm%max_change_isubvar(6))
  pm%max_change_isubvar = [0,0,0,2,0,0]
  pm%max_change_ivar = [LIQUID_PRESSURE, GAS_PRESSURE, CO2_PRESSURE, &
                                LIQUID_MASS_FRACTION, LIQUID_SALT_MASS_FRAC, &
                                GAS_SATURATION]
                                ! TEMPERATURE, GAS_SATURATION]
  pm%damping_factor = -1.d0

  pm%residual_abs_inf_tol = residual_abs_inf_tol
  pm%residual_scaled_inf_tol = residual_scaled_inf_tol
  pm%abs_update_inf_tol = abs_update_inf_tol
  pm%rel_update_inf_tol = rel_update_inf_tol

  pm%converged_flag(:,:,:) = PETSC_TRUE
  pm%converged_real(:,:,:) = 0.d0
  pm%converged_cell(:,:,:) = 0

end subroutine PMSCO2SetFlowMode

! ************************************************************************** !
subroutine PMSCO2ReadSimOptionsBlock(this,input)
  !
  ! Read simulation options for SCO2 mode
  !
  ! Author: Michael Nole
  ! Date: 11/21/2023
  !
  use SCO2_module
  use SCO2_Aux_module
  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: keyword
  class(pm_sco2_type) :: this
  type(option_type), pointer :: option
  PetscReal :: tempreal
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  option => this%option

  error_string = 'SCO2 Options'

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
      case('CHECK_MAX_DPL_LIQ_STATE_ONLY')
        sco2_chk_max_dpl_liq_state_only = PETSC_TRUE
      case('DEBUG_CELL')
        call InputReadInt(input,option,sco2_debug_cell_id)
        call InputErrorMsg(input,option,keyword,error_string)
      case('NO_STATE_TRANSITION_OUTPUT')
        sco2_print_state_transition = PETSC_FALSE
      case('PHASE_CHANGE_EPSILON')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,error_string)
        sco2_phase_chng_epsilon = tempreal
      case('RESTRICT_STATE_CHANGE')
        sco2_restrict_state_chng = PETSC_TRUE
      case('WINDOW_EPSILON')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,error_string)
        sco2_window_epsilon = tempreal
      case('TEMPERATURE')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,error_string)
        sco2_isothermal_temperature = tempreal
      case default
        call InputKeywordUnrecognized(input,keyword,'SCO2 Mode',option)
    end select

  enddo
  call InputPopBlock(input,option)

end subroutine PMSCO2ReadSimOptionsBlock

! ************************************************************************** !

subroutine PMSCO2ReadNewtonSelectCase(this,input,keyword,found, &
                                      error_string,option)
  !
  ! Reads input file parameters associated with the SCO2 process model for
  ! Newton solver convergence
  !
  ! Author: Michael Nole
  ! Date: 11/21/23

  use Input_Aux_module
  use String_module
  use Utility_module
  use Option_module
  use SCO2_Aux_module

  implicit none

  class(pm_sco2_type) :: this
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

  error_string = 'SCO2 Newton Solver'

  found = PETSC_FALSE
  call PMSubsurfaceFlowReadNewtonSelectCase(this,input,keyword,found, &
                                            error_string,option)
  if (found) return

  found = PETSC_TRUE
  select case(trim(keyword))
    case('USE_GOVERNORS')
      sco2_use_governors = PETSC_TRUE
    case('CHECK_SOLUTION_UPDATES')
      sco2_check_updates = PETSC_TRUE
    case('USE_FULL_CONVERGENCE_CRITERIA')
      sco2_stomp_convergence = PETSC_FALSE
    case('MAX_NEWTON_ITERATIONS')
      call InputKeywordDeprecated('MAX_NEWTON_ITERATIONS', &
                                  'MAXIMUM_NUMBER_OF_ITERATIONS.',option)
    ! Tolerances

    ! All Residual
    case('CENTRAL_DIFFERENCE_JACOBIAN')
      sco2_central_diff_jacobian = PETSC_TRUE 
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
    ! case('ENERGY_RESIDUAL_ABS_INF_TOL')
    !   call InputReadDouble(input,option,this%residual_abs_inf_tol(eid))
    !   call InputErrorMsg(input,option,keyword,error_string)
    case('SALT_RESIDUAL_ABS_INF_TOL')
      call InputReadDouble(input,option,this%residual_abs_inf_tol(sid))
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
    case('SALT_RESIDUAL_SCALED_INF_TOL')
      call InputReadDouble(input,option,this%residual_scaled_inf_tol(sid))
      call InputErrorMsg(input,option,keyword,error_string)
    ! case('ENERGY_RESIDUAL_SCALED_INF_TOL')
    !   call InputReadDouble(input,option,this%residual_scaled_inf_tol(eid))
    !   call InputErrorMsg(input,option,keyword,error_string)

    ! All Updates
    case('UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%abs_update_inf_tol(:,:) = tempreal
      this%rel_update_inf_tol(:,:) = tempreal

    ! Update Map:
    !  State   Primary1 Primary2 Primary3 Primary4
    !    L        Pl     Xw_CO2      T     X_salt
    !    G        Pg     Pco2        T     X_salt
    !    TG       Pl     Sg          T     X_salt
    !    LG       Pg     Pl          T     X_salt


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
      this%abs_update_inf_tol(2,4) = tempreal
    ! case('TEMP_ABS_UPDATE_INF_TOL')
    !   call InputReadDouble(input,option,tempreal)
    !   call InputErrorMsg(input,option,keyword,error_string)
    !   this%abs_update_inf_tol(3,:) = tempreal
    case('SAT_ABS_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%abs_update_inf_tol(2,3) = tempreal
    case('XMASS_ABS_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%abs_update_inf_tol(2,1) = tempreal
    case('XMASS_SALT_ABS_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%abs_update_inf_tol(3,1:4) = tempreal

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
      this%rel_update_inf_tol(2,4) = tempreal
    ! case('TEMP_REL_UPDATE_INF_TOL')
    !   call InputReadDouble(input,option,tempreal)
    !   call InputErrorMsg(input,option,keyword,error_string)
    !   this%rel_update_inf_tol(3,:) = tempreal
    case('SAT_REL_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%rel_update_inf_tol(2,3) = tempreal
    case('XMASS_REL_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%rel_update_inf_tol(2,1) = tempreal
    case('XMASS_SALT_REL_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%rel_update_inf_tol(3,1:4) = tempreal

    ! Other Controls
    case('MAXIMUM_PRESSURE_CHANGE')
      call InputReadDouble(input,option,sco2_max_pressure_change)
      call InputErrorMsg(input,option,keyword,error_string)
    case('MAX_ITERATION_BEFORE_DAMPING')
      call InputReadInt(input,option,sco2_max_it_before_damping)
      call InputErrorMsg(input,option,keyword,error_string)
    case('DAMPING_FACTOR')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%damping_factor = tempreal
    case default
      found = PETSC_FALSE

  end select

end subroutine PMSCO2ReadNewtonSelectCase

! ************************************************************************** !

subroutine PMSCO2InitializeSolver(this)
  !
  ! Author: Michael Nole
  ! Date: 11/21/23

  use Solver_module

  implicit none

  class(pm_sco2_type) :: this

  call PMBaseInitializeSolver(this)

  ! MAN: Do we need dtol?
  this%solver%newton_dtol = 1.d9
  this%solver%newton_max_iterations = 16

end subroutine PMSCO2InitializeSolver

! ************************************************************************** !

recursive subroutine PMSCO2InitializeRun(this)
  !
  ! Initializes the time stepping
  !
  ! Author: Michael Nole
  ! Date: 11/21/23

  use Realization_Base_class

  implicit none

  class(pm_sco2_type) :: this

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

end subroutine PMSCO2InitializeRun

! ************************************************************************** !

subroutine PMSCO2InitializeTimestep(this)
  !
  ! Should not need this as it is called in PreSolve, but
  ! the base routine must be extended.
  !
  ! Author: Michael Nole
  ! Date: 12/18/23
  !

  use SCO2_module, only : SCO2InitializeTimestep
  use Global_module
  use Variables_module, only : TORTUOSITY
  use Material_module, only : MaterialAuxVarCommunicate
  use Option_module

  implicit none

  class(pm_sco2_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)
!geh:remove   everywhere
  call MaterialAuxVarCommunicate(this%comm1, &
                                 this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc,TORTUOSITY, &
                                 ZERO_INTEGER)

  call SCO2InitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)

end subroutine PMSCO2InitializeTimestep

! ************************************************************************** !

subroutine PMSCO2PreSolve(this)
  !
  ! Author: Michael Nole
  ! Date: 11/21/23

  implicit none

  class(pm_sco2_type) :: this

  call PMSubsurfaceFlowPreSolve(this)

end subroutine PMSCO2PreSolve

! ************************************************************************** !

subroutine PMSCO2PostSolve(this)
  !
  ! MAN: Do we need this?
  !
  ! Author: Michael Nole
  ! Date: 11/21/23

  implicit none

  class(pm_sco2_type) :: this

end subroutine PMSCO2PostSolve

! ************************************************************************** !

subroutine PMSCO2UpdateTimestep(this,update_dt, &
                                dt,dt_min,dt_max,iacceleration, &
                                num_newton_iterations,tfac, &
                                time_step_max_growth_factor)
  !
  ! Author: Michael Nole
  ! Date: 11/21/23
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

  class(pm_sco2_type) :: this
  PetscBool :: update_dt
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor

  PetscReal :: fac
  PetscInt :: ifac
  PetscReal :: up, ux, us, um, umin
  ! PetscReal :: ut
  PetscReal :: dtt
  PetscReal :: governed_dt
  PetscReal :: umin_scale
  PetscReal :: value
  PetscReal :: governor_value
  character(MAXSTRINGLENGTH) :: string
  type(field_type), pointer :: field

#ifdef PM_SCO2_DEBUG
  call PrintMsg(this%option,'PMSCO2%UpdateTimestep()')
#endif

  if (update_dt .and. iacceleration /= 0) then
    fac = 0.5d0
    if (num_newton_iterations >= iacceleration) then
      fac = 0.33d0
      umin = 0.d0
    else
      up = this%pressure_change_governor/(this%max_pressure_change+0.1)
      ! ut = this%temperature_change_governor/(this%max_temperature_change+1.d-5)
      ux = this%xmol_change_governor/(this%max_xmol_change+1.d-5)
      us = this%saturation_change_governor/(this%max_saturation_change+1.d-5)
      um = this%salt_mass_change_governor/(this%max_salt_mass_change+1.d-5)
      umin = min(up,ux,us)
      ! umin = min(up,ut,ux,us)
    endif
    ifac = max(min(num_newton_iterations,size(tfac)),1)
    umin_scale = fac * (1.d0 + umin)
    if (sco2_use_governors) then
      governed_dt = umin_scale * dt
      dtt = min(time_step_max_growth_factor*dt,governed_dt)
      dt = min(dtt,tfac(ifac)*dt,dt_max)
      dt = max(dt,dt_min)
    else
      dtt = time_step_max_growth_factor*dt
      dt = min(dtt,dt_max)
      dt = max(dt,dt_min)
    endif
    dt = min(dtt,dt_max)
    dt = max(dt,dt_min)

    ! Inform user that time step is being limited by a state variable.
    if (Equal(dt,governed_dt)) then
      umin = umin * (1.d0 + 1.d-8)
      if (up < umin) then
        string = 'Pressure'
        value = this%max_pressure_change
        governor_value = this%pressure_change_governor
      ! else if (ut < umin) then
      !   string = 'Temperature'
      !   value = this%max_temperature_change
      !   governor_value = this%temperature_change_governor
      else if (ux < umin) then
        string = 'CO2 Mass Fraction'
        value = this%max_xmol_change
        governor_value = this%xmol_change_governor
      else if (us < umin) then
        string = 'Saturation'
        value = this%max_saturation_change
        governor_value = this%saturation_change_governor
      else if (um < umin) then
        string = 'Salt Mass Fraction'
        value = this%max_salt_mass_change
        governor_value = this%salt_mass_change_governor
      else
        string = 'Newton Iterations'
        value = num_newton_iterations
        governor_value = iacceleration + 0.d0
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
    ! Since saturations are not stored in global_auxvar for SCO2 mode, we
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

end subroutine PMSCO2UpdateTimestep

! ************************************************************************** !

subroutine PMSCO2Residual(this,snes,xx,r,ierr)
  !
  ! Author: Michael Nole
  ! Date: 11/21/23
  !

  use SCO2_module, only : SCO2Residual

  implicit none

  class(pm_sco2_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr

  call PMSubsurfaceFlowUpdatePropertiesNI(this)
  call SCO2Residual(snes,xx,r,this%realization,ierr)

end subroutine PMSCO2Residual

! ************************************************************************** !

subroutine PMSCO2Jacobian(this,snes,xx,A,B,ierr)
  !
  ! Author: Michael Nole
  ! Date: 11/21/23
  !

  use SCO2_module, only : SCO2Jacobian

  implicit none

  class(pm_sco2_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr

  call SCO2Jacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMSCO2Jacobian

! ************************************************************************** !

subroutine PMSCO2CheckUpdatePre(this,snes,X,dX,changed,ierr)
  !
  ! Author: Michael Nole
  ! Date: 11/21/23
  !

  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Option_module
  use Saturation_Function_module
  use Patch_module
  use SCO2_Aux_module
  use Global_Aux_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Common_module

  implicit none

  class(pm_sco2_type) :: this
  SNES :: snes
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr

  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(sco2_auxvar_type), pointer :: sco2_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)

  type(sco2_auxvar_type) :: sco2_auxvar
  class(characteristic_curves_type), pointer :: characteristic_curves
  PetscReal, pointer :: X_p(:),dX_p(:),dX_p2(:)
  PetscInt :: liq_pressure_index, gas_pressure_index, co2_frac_index, &
              gas_sat_index, co2_pressure_index, salt_index
  PetscInt :: local_id, ghosted_id, offset
  PetscInt :: lid
  PetscReal :: dP, dsg, Pc_max, Psb, Pvb, rho_b, Pc, Pc_entry
  PetscReal :: xsl
  PetscReal, parameter :: epsilon = 1.d-14
  PetscReal, parameter :: gravity = EARTH_GRAVITY

  field => this%realization%field
  grid => this%realization%patch%grid
  patch => this%realization%patch
  option => this%realization%option
  sco2_auxvars => this%realization%patch%aux%SCO2%auxvars
  global_auxvars => this%realization%patch%aux%Global%auxvars

  call VecCopy(dX,field%flow_dxx,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%flow_dxx,dX_p2,ierr);CHKERRQ(ierr)

  lid = option%liquid_phase

  dX_p = -1.d0 * dX_p
  dX_p2 = -1.d0 * dX_p2

  changed = PETSC_TRUE

  if (this%check_post_convergence) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      offset = (local_id-1)*option%nflowdof
      sco2_auxvar = sco2_auxvars(ZERO_INTEGER,ghosted_id)
      characteristic_curves => patch%characteristic_curves_array( &
                               patch%cc_id(ghosted_id))%ptr
      select case(global_auxvars(ghosted_id)%istate)

        case(SCO2_LIQUID_STATE)
          ! Limit pressure change
          liq_pressure_index = offset + ONE_INTEGER
          dP = 1.d6
          dX_p(liq_pressure_index) = sign( min(dabs(dP), &
          dabs(dX_p(liq_pressure_index))),dX_p(liq_pressure_index) )
          if ((X_p(liq_pressure_index) + dX_p(liq_pressure_index)) > 5.d8) &
            dX_p(liq_pressure_index) = 5.d8 - X_p(liq_pressure_index)

          ! Zero negative corrections for zero aqueous CO2
          co2_frac_index = offset + TWO_INTEGER
          if (X_p(co2_frac_index) / epsilon < epsilon .and. &
              dX_p(co2_frac_index) / epsilon < epsilon ) then
            dX_p(co2_frac_index) = 0.d0
            dX_p2(co2_frac_index) = 0.d0
          endif
          if ((X_p(co2_frac_index) + dX_p(co2_frac_index)) < 0.d0) &
             dX_p(co2_frac_index) = - X_p(co2_frac_index)

          ! Limit salt mass fraction changes to 0.25 of max
          salt_index = offset + THREE_INTEGER
          call SCO2ComputeSaltSolubility(sco2_auxvar%temp,xsl)
          if (X_p(salt_index) < xsl ) THEN
            dX_p(salt_index) = sign(min(dabs(2.5d-1*xsl), &
                               dabs(dX_p(salt_index))), dX_p(salt_index))
          endif
          ! Zero negative corrections without salt present
          if (X_p(salt_index) / epsilon < epsilon .and. &
             dX_p(salt_index) / epsilon < epsilon) then
            dX_p(salt_index) = 0.d0
            dX_p2(salt_index) = 0.d0
          endif
          if ((X_p(salt_index) + dX_p(salt_index)) < epsilon) &
             dX_p(salt_index) = - X_p(salt_index)

        case(SCO2_GAS_STATE)
          gas_pressure_index = offset + ONE_INTEGER
          co2_pressure_index = offset + TWO_INTEGER
          salt_index = offset + THREE_INTEGER

          ! Limit changes in gas pressure  ---

          dP = 2.5d-1*max(X_p(gas_pressure_index),1.d6)
          dX_p(gas_pressure_index) = sign(min(dabs(dP), &
                                     dabs(dX_p(gas_pressure_index))), &
                                     dabs(dX_p(gas_pressure_index)))

          ! Limit changes in CO2 pressure  ---

          if (X_p(co2_pressure_index) / epsilon < epsilon .and. &
              dX_p(co2_pressure_index) / epsilon < epsilon) then
            dX_p(co2_pressure_index) = 0.d0
            dX_p2(co2_pressure_index) = 0.d0
          endif
          call SCO2ComputeSaltSolubility(sco2_auxvar%temp,xsl)
          xsl = min(X_p(salt_index),xsl)
          call SCO2BrineSaturationPressure(sco2_auxvar%temp, &
                                         xsl,Psb)
          dP = X_p(gas_pressure_index) -  2.5d-1*Psb
          dX_p(co2_pressure_index) = sign(min(dabs(dP),&
                                     dabs(dX_p(co2_pressure_index))), &
                                     dX_p(co2_pressure_index))
          if ((X_p(co2_pressure_index) + dX_p(co2_pressure_index)) > &
            (1.d0 - 1.d-6) * X_p(gas_pressure_index)) &
            dX_p(co2_pressure_index) = X_p(gas_pressure_index) - &
            X_p(co2_pressure_index)

          ! Zero negative corrections for salt volumetric conc.  ---
            if (X_p(salt_index)/epsilon < epsilon .and. &
                dX_p(salt_index)/epsilon < epsilon ) then
              dX_p(salt_index) = 0.d0
              dX_p2(salt_index) = 0.d0
            endif
            if ((X_p(salt_index) + dX_p(salt_index)) < epsilon) &
              dX_p(salt_index) = - X_p(salt_index)

          call SCO2BrineDensity(sco2_auxvar%temp,Psb,xsl,rho_b,option)
          Pc = max(Psb - sco2_auxvar%pres(lid),0.d0 )
          call SCO2VaporPressureBrine(sco2_auxvar%temp, Psb, &
                                      Pc,rho_b,xsl,Pvb)
          if ((X_p(gas_pressure_index) + dX_p(gas_pressure_index)) < Pvb) &
            dX_p(gas_pressure_index) = Pvb - X_p(gas_pressure_index)

        case(SCO2_TRAPPED_GAS_STATE)
          liq_pressure_index = offset + ONE_INTEGER
          gas_sat_index = offset + TWO_INTEGER
          salt_index = offset + THREE_INTEGER
          ! Limit changes in pressure  ---

          dP = 2.5d-1*max(X_p(liq_pressure_index),1.d6)
          dX_p(liq_pressure_index) = sign(min(dabs(dP), &
                                     dabs(dX_p(liq_pressure_index))), &
                                     dX_p(liq_pressure_index))
          if ((X_p(liq_pressure_index) + dX_p(liq_pressure_index)) > &
            5.d8) dX_p(liq_pressure_index) = 5.d8 - X_p(liq_pressure_index)

          ! Limit changes in trapped gas  ---

          dsg = 1.d-1*characteristic_curves%saturation_function%Sgt_max
          dX_p(gas_sat_index) = sign(min(dabs(dsg), &
                                dabs(dX_p(gas_sat_index))),&
                                dX_p(gas_sat_index))
          if (X_p(gas_sat_index) + dX_p(gas_sat_index) < epsilon) &
              dX_p(gas_sat_index) = - X_p(gas_sat_index)

          ! Limit salt mass fraction changes to 0.25 of max
          salt_index = offset + THREE_INTEGER
          call SCO2ComputeSaltSolubility(sco2_auxvar%temp,xsl)
          if (X_p(salt_index) < xsl ) THEN
            dX_p(salt_index) = sign(min(dabs(2.5d-1*xsl), &
                               dabs(dX_p(salt_index))), dX_p(salt_index))
          endif
          ! Zero negative corrections without salt present
          if (X_p(salt_index) / epsilon < epsilon .and. &
             dX_p(salt_index) / epsilon < epsilon) then
            dX_p(salt_index) = 0.d0
            dX_p2(salt_index) = 0.d0
          endif
          if ((X_p(salt_index) + dX_p(salt_index)) < epsilon) &
             dX_p(salt_index) = - X_p(salt_index)

        case(SCO2_LIQUID_GAS_STATE)
          liq_pressure_index = offset + ONE_INTEGER
          gas_pressure_index = offset + TWO_INTEGER
          salt_index = offset + THREE_INTEGER

          Pc_entry = 0.d0
          !select type(sf => characteristic_curves%saturation_function)
          !  class is (sat_func_vg_type)
          !    Pc_entry = (1.d0 / characteristic_curves% &
          !                saturation_function%GetAlpha_())
          !  class is (sat_func_VG_STOMP_type)
          !    Pc_entry = characteristic_curves% &
          !               saturation_function%GetAlpha_() * &
          !               LIQUID_REFERENCE_DENSITY * gravity
          !  class default
          !end select

          !Limit changes in pressure
          dP = max(1.d6,2.5d-1*(X_p(gas_pressure_index) - &
              X_p(liq_pressure_index)))
          dX_p(liq_pressure_index) = sign(min(dabs(dP), &
                                     dabs(dX_p(liq_pressure_index))), &
                                     dX_p(liq_pressure_index))
          dX_p(gas_pressure_index) = sign(min(dabs(dP), &
                                     dabs(dX_p(gas_pressure_index))), & 
                                     dX_p(gas_pressure_index))
          ! Relax pressure updates when transitioning to unsaturated
          ! conditions
          if((X_p(gas_pressure_index) + dX_p(gas_pressure_index)) - &
             (X_p(liq_pressure_index) + dX_p(liq_pressure_index)) < &
              Pc_entry ) then
            dX_p(gas_pressure_index) = 6.d-1*dX_p(gas_pressure_index)
            dX_p(liq_pressure_index) = 6.d-1*dX_p(liq_pressure_index)
          endif
          if ((X_p(gas_pressure_index) + dX_p(gas_pressure_index)) > &
             5.d8) dX_p(gas_pressure_index) = 5.d8 - X_p(gas_pressure_index)
          select type(sf => characteristic_curves%saturation_function)
            class is (sat_func_vg_type)
              Pc_max = characteristic_curves%saturation_function%Pcmax
            class is (sat_func_VG_STOMP_type)
              Pc_max = characteristic_curves%saturation_function%Pcmax * &
                         LIQUID_REFERENCE_DENSITY * gravity
            class default
          end select
          if ((X_p(liq_pressure_index) + dX_p(liq_pressure_index)) < & 
            ((X_p(gas_pressure_index) + dX_p(gas_pressure_index)) - &
             Pc_max)) &
            dX_p(liq_pressure_index) = ((X_p(gas_pressure_index) + &
                 dX_p(gas_pressure_index)) - Pc_max) - X_p(liq_pressure_index)

          ! Limit salt mass fraction changes to 0.25 of max
          call SCO2ComputeSaltSolubility(sco2_auxvar%temp,xsl)
          if (X_p(salt_index) < xsl ) THEN
            dX_p(salt_index) = sign(min(dabs(2.5d-1*xsl), &
                               dabs(dX_p(salt_index))), dX_p(salt_index))
          endif
          ! Zero negative corrections without salt present
          if (X_p(salt_index) / epsilon < epsilon .and. &
             dX_p(salt_index) / epsilon < epsilon) then
            dX_p(salt_index) = 0.d0
            dX_p2(salt_index) = 0.d0
          endif
          if ((X_p(salt_index) + dX_p(salt_index)) < epsilon) &
             dX_p(salt_index) = - X_p(salt_index)

          !Maintain the gas pressure above or at the water vapor
          !pressure  
          xsl = min(X_p(salt_index),xsl)
          call SCO2BrineSaturationPressure(sco2_auxvar%temp, &
                                         xsl,Psb)
          call SCO2BrineDensity(sco2_auxvar%temp,Psb,xsl,rho_b,option)
          Pc = max(Psb-X_p(liq_pressure_index),0.d0)
          call SCO2VaporPressureBrine(sco2_auxvar%temp, Psb, &
                                      Pc,rho_b,xsl,Pvb)
          if ((X_p(gas_pressure_index) + dX_p(gas_pressure_index)) < &
              Pvb) dX_p(gas_pressure_index) = Pvb - X_p(gas_pressure_index)

      end select

      ! Throw an error or cut timestep
      !if ((X_p(gas_pressure_index) + dX_p(gas_pressure_index)) > 8.d8 .or. &
      !    (X_p(gas_pressure_index) + dX_p(gas_pressure_index)) < 0.d0)
      !if ((X_p(liq_pressure_index) + dX_p(liq_pressure_index)) > 8.d8)
      !if (sco2_auxvar%temp + 273.15 > H2O_CRITICAL_TEMPERATURE .or. &
      !    sco2_auxvar%temp < 0.d0 )

    enddo

    if (this%damping_factor > 0.d0) then
      dX_p = dX_p*this%damping_factor
      changed = PETSC_TRUE
    endif
  endif

  dX_p = -1.d0 * dX_p

  call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%flow_dxx,dX_p2,ierr);CHKERRQ(ierr)

end subroutine PMSCO2CheckUpdatePre

! ************************************************************************** !

subroutine PMSCO2CheckUpdatePost(this,snes,X0,dX,X1,dX_changed, &
                                 X1_changed,ierr)
  !
  ! Author: Michael Nole
  ! Date: 11/21/23
  !
  use SCO2_Aux_module
  use Global_Aux_module
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module

  implicit none

  class(pm_sco2_type) :: this
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

  allocate(converged_abs_update_flag(option%nflowdof,sco2_max_states))
  allocate(converged_rel_update_flag(option%nflowdof,sco2_max_states))
  allocate(converged_abs_update_cell(option%nflowdof,sco2_max_states))
  allocate(converged_rel_update_cell(option%nflowdof,sco2_max_states))
  allocate(converged_abs_update_real(option%nflowdof,sco2_max_states))
  allocate(converged_rel_update_real(option%nflowdof,sco2_max_states))

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
 
  if (sco2_check_updates) then 
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

        if (.not.(converged_absolute .or. converged_relative)) then
          converged_abs_update_flag(idof,istate) = PETSC_FALSE
          converged_rel_update_flag(idof,istate) = PETSC_FALSE
        endif
      enddo
    enddo
  endif

  call VecRestoreArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)

  this%converged_flag(:,:,ABS_UPDATE_INDEX) = converged_abs_update_flag(:,:)
  this%converged_flag(:,:,REL_UPDATE_INDEX) = converged_rel_update_flag(:,:)
  this%converged_real(:,:,ABS_UPDATE_INDEX) = converged_abs_update_real(:,:)
  this%converged_real(:,:,REL_UPDATE_INDEX) = converged_rel_update_real(:,:)
  this%converged_cell(:,:,ABS_UPDATE_INDEX) = converged_abs_update_cell(:,:)
  this%converged_cell(:,:,REL_UPDATE_INDEX) = converged_rel_update_cell(:,:)

end subroutine PMSCO2CheckUpdatePost

! ************************************************************************** !

subroutine PMSCO2CheckConvergence(this,snes,it,xnorm,unorm,fnorm, &
                                  reason,ierr)
  !
  ! Author: Michael Nole
  ! Date: 11/21/23
  !
  use Convergence_module
  use SCO2_Aux_module
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

  class(pm_sco2_type) :: this
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
  type(sco2_auxvar_type), pointer :: sco2_auxvars(:,:)
  type(sco2_auxvar_type) :: sco2_auxvar
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum2_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscInt :: local_id, ghosted_id, natural_id
  PetscInt :: offset, ival, idof, itol
  PetscInt :: lid, gid, co2_id, vpid
  PetscReal :: R, A, R_A
  PetscReal :: res_scaled, residual, accumulation, update
  PetscReal :: Psat, Pv, Prvap, Pco2
  PetscReal :: xco2g, xwg, xco2l, xsl, xwl, xmolco2g, xmolwg, xmolco2l, &
               xmolsl, xmolwl, x_salt_dissolved
  PetscReal :: Hc, den_salt
  PetscReal, parameter :: A_zero = 1.d-15
  PetscReal, parameter :: epsilon = 1.d-20
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
  character(len=20), allocatable :: state_string(:)
  character(len=17), allocatable :: dof_string(:,:)
  character(len=15), parameter :: tol_string(MAX_INDEX) = &
    ['Absolute Update','Relative Update','Residual       ','Scaled Residual']

  patch => this%realization%patch
  option => this%realization%option
  field => this%realization%field
  grid => patch%grid
  global_auxvars => patch%aux%Global%auxvars
  sco2_auxvars => this%realization%patch%aux%SCO2%auxvars

  allocate(flags(MAX_INDEX*option%nflowdof*sco2_max_states))
  ! allocate(flags(MAX_INDEX*option%nflowdof*sco2_max_states+1))
  allocate(state_string(sco2_max_states))
  allocate(dof_string(option%nflowdof,sco2_max_states))

  allocate(converged_abs_residual_flag(option%nflowdof,sco2_max_states))
  allocate(converged_abs_residual_real(option%nflowdof,sco2_max_states))
  allocate(converged_abs_residual_cell(option%nflowdof,sco2_max_states))
  allocate(converged_scaled_residual_flag(option%nflowdof,sco2_max_states))
  allocate(converged_scaled_residual_real(option%nflowdof,sco2_max_states))
  allocate(converged_scaled_residual_cell(option%nflowdof,sco2_max_states))

  lid = option%liquid_phase
  gid = option%gas_phase
  co2_id = option%co2_id
  vpid = option%vapor_pressure_id

  state_string = &
    ['Liquid State        ','Gas State           ', &
     'Trapped Gas State   ','Liquid-Gas State    ']
     !Without Energy
     dof_string = &
     reshape(['Liquid Pressure  ','CO2 Fraction     ', &
              'Salt Mass Frac   ', &
              'Gas Pressure     ','CO2 Partial Pres ', &
              'Salt Mass Frac   ', &
              'Gas Pressure     ','Gas Saturation   ', &
              'Salt Mass Frac   ', &
              'Gas Pressure     ','Liquid Pressure  ', &
              'Salt Mass Frac   '], &
             shape(dof_string))
  ! With Energy
  ! dof_string = &
  !   reshape(['Liquid Pressure  ','CO2 Fraction     ','Temperature      ', &
  !            'Salt Mass Frac   ', &
  !            'Gas Pressure     ','CO2 Partial Pres ','Temperature      ', &
  !            'Salt Mass Frac   ', &
  !            'Gas Pressure     ','Gas Saturation   ','Temperature      ', &
  !            'Salt Mass Frac   ', &
  !            'Gas Pressure     ','Liquid Pressure  ','Temperature      ', &
  !            'Salt Mass Frac   '], &
  !           shape(dof_string))

  call SNESNewtonTRDCGetRhoFlag(snes,rho_flag,ierr);CHKERRQ(ierr);

  if (this%option%flow%using_newtontrdc) then
    if (sco2_newtontrdc_prev_iter_num == it) then
      sco2_sub_newton_iter_num = sco2_sub_newton_iter_num + 1
    endif
    sco2_newtontrdc_prev_iter_num = it
  endif
  if (this%check_post_convergence) then
    call VecGetArrayReadF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(field%flow_accum2,accum2_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(field%flow_dxx,dX_p,ierr);CHKERRQ(ierr)
    converged_abs_residual_flag = PETSC_TRUE
    converged_abs_residual_real = 0.d0
    converged_abs_residual_cell = ZERO_INTEGER
    converged_scaled_residual_flag = PETSC_TRUE
    converged_scaled_residual_real = 0.d0
    converged_scaled_residual_cell = ZERO_INTEGER

    if (sco2_stomp_convergence) then
      do local_id = 1, grid%nlmax
        offset = (local_id-1)*option%nflowdof
        ghosted_id = grid%nL2G(local_id)
        natural_id = grid%nG2A(ghosted_id)

        sco2_auxvar = sco2_auxvars(ZERO_INTEGER,ghosted_id)

        if (patch%imat(ghosted_id) <= 0) cycle
        istate = global_auxvars(ghosted_id)%istate
        do idof = 1, option%nflowdof
          res_scaled = 0.d0
          ival = offset+idof
          !converged_absolute = PETSC_TRUE
          converged_absolute = PETSC_FALSE
          converged_scaled = PETSC_TRUE

          residual = r_p(ival)
          accumulation = accum2_p(ival)
          update = dX_p(ival)

          ! STOMP convergence criteria
          select case (istate)
            case(SCO2_LIQUID_STATE)

              if (idof == ONE_INTEGER) then
              !---      Water mass equation  ---

                res_scaled = min(dabs(update) / (dabs(sco2_auxvar%pres(lid))), &
                             dabs(residual/(accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif

              elseif (idof == TWO_INTEGER) then

              !---      CO2 mass equation, ignore residual for small aqueous-CO2
                call SCO2ComputeSaltSolubility(sco2_auxvar%temp,xsl)
                xsl = min(sco2_auxvar%m_salt(1),xsl)
                call SCO2BrineSaturationPressure(sco2_auxvar%temp,xsl,Psat) 
                Pv = sco2_auxvar%pres(vpid)
                Prvap = Psat
                x_salt_dissolved = xsl
                call SCO2Equilibrate(sco2_auxvar%temp,sco2_auxvar%pres(lid), &
                                     Pco2, Pv, Psat, Prvap, &
                                     xco2g, xwg, xco2l, xsl, xwl, &
                                     xmolco2g, xmolwg, xmolco2l, xmolsl, &
                                     xmolwl, option)
                x_salt_dissolved = x_salt_dissolved + &
                                   (xsl - x_salt_dissolved) * &
                                   (sco2_auxvar%xmass(co2_id,lid) / xco2l)
                if (sco2_auxvar%xmass(co2_id,lid) > (1.d-6 * xco2l)) then
                  Hc = SCO2Henry(sco2_auxvar%temp,x_salt_dissolved)
                  res_scaled = min(dabs(update) / &
                               max(sco2_auxvar%pres(gid)/Hc, &
                               SCO2_REFERENCE_PRESSURE/Hc), &
                               dabs(residual / (accumulation + epsilon)))
                  ! find max value regardless of convergence
                  if (converged_scaled_residual_real(idof,istate) < &
                      res_scaled) then
                   converged_scaled_residual_real(idof,istate) = res_scaled
                   converged_scaled_residual_cell(idof,istate) = natural_id
                  endif
                endif

              elseif (idof == THREE_INTEGER) then

              !---      Salt mass equation  ---
                call SCO2ComputeSaltSolubility(sco2_auxvar%temp,xsl)
                res_scaled = min(dabs(update)/xsl, &
                             dabs(residual / (accumulation + epsilon)))
                res_scaled = 1.d-1 * res_scaled

              ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                  converged_scaled_residual_real(idof,istate) = res_scaled
                  converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              endif
            case(SCO2_GAS_STATE)

              if (idof == ONE_INTEGER) then
              !---      Water mass equation  ---
                call SCO2ComputeSaltSolubility(sco2_auxvar%temp,xsl)
                if (sco2_auxvar%m_salt(2) <= 0.d0) xsl = 0.d0
                call SCO2BrineSaturationPressure(sco2_auxvar%temp,xsl,Psat)

                res_scaled = min(dabs(update) / Psat, &
                             dabs(residual / (accumulation + epsilon)))

                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                  converged_scaled_residual_real(idof,istate) = res_scaled
                  converged_scaled_residual_cell(idof,istate) = natural_id
                endif

              elseif (idof == TWO_INTEGER) then
              !---      Gas mass equation  ---
                res_scaled = min(dabs(update)/dabs(sco2_auxvar%pres(gid)), &
                             dabs(residual / (accumulation + epsilon)))

                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                  converged_scaled_residual_real(idof,istate) = res_scaled
                  converged_scaled_residual_cell(idof,istate) = natural_id
                endif

              elseif (idof == THREE_INTEGER) then
              !---      Salt mass equation, isobrine option  ---

                call SCO2ComputeSaltDensity(sco2_auxvar%temp, &
                                            sco2_auxvar%pres(gid), &
                                            den_salt)
                res_scaled = min(dabs(update) / den_salt, &
                             dabs(residual / (accumulation + epsilon)))
              endif
  
            case(SCO2_TRAPPED_GAS_STATE)
              if (idof == ONE_INTEGER) then
              !---      Water mass equation  ---

                res_scaled = min(dabs(update)/dabs(sco2_auxvar%pres(lid)), &
                           dabs(residual / (accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                    converged_scaled_residual_real(idof,istate) = res_scaled
                    converged_scaled_residual_cell(idof,istate) = natural_id
                endif

              elseif (idof == TWO_INTEGER) then

              !---      CO2 mass equation  ---
                if (sco2_auxvar%sat(gid) > 1.d-5) then
                  res_scaled = min(dabs(update/1.d1), &
                               dabs(residual / (accumulation + epsilon)))
                  ! find max value regardless of convergence
                  if (converged_scaled_residual_real(idof,istate) < &
                      res_scaled) then
                      converged_scaled_residual_real(idof,istate) = res_scaled
                      converged_scaled_residual_cell(idof,istate) = natural_id
                  endif

                endif

              elseif (idof == THREE_INTEGER) then
              !---      Salt mass equation  ---
                call SCO2ComputeSaltSolubility(sco2_auxvar%temp,xsl)
                res_scaled = min(dabs(update) / xsl, &
                             dabs(residual / (accumulation + epsilon)))
                res_scaled = 1.d-1 * res_scaled
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                  converged_scaled_residual_real(idof,istate) = res_scaled
                  converged_scaled_residual_cell(idof,istate) = natural_id
                endif

              endif

            case(SCO2_LIQUID_GAS_STATE)

              if (idof == ONE_INTEGER) then

              !---      Water mass equation  ---
                res_scaled = min(dabs(update) / dabs(sco2_auxvar%pres(lid)), &
                             dabs(residual / (accumulation + epsilon)))

                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                  converged_scaled_residual_real(idof,istate) = res_scaled
                  converged_scaled_residual_cell(idof,istate) = natural_id
                endif

              elseif (idof == TWO_INTEGER) then
              !---      CO2 mass equation  ---
                if (sco2_auxvar%sat(gid) > 1.d-3) then
                  res_scaled = min(dabs(update)/dabs(sco2_auxvar%pres(gid)), &
                               dabs(residual / (accumulation + epsilon)))

                  if (converged_scaled_residual_real(idof,istate) < &
                      res_scaled) then
                    converged_scaled_residual_real(idof,istate) = res_scaled
                    converged_scaled_residual_cell(idof,istate) = natural_id
                  endif

                endif
              elseif (idof == THREE_INTEGER) then
              !---      Salt mass equation  ---
                call SCO2ComputeSaltSolubility(sco2_auxvar%temp,xsl)
                res_scaled = min(dabs(update) / xsl, &
                             dabs(residual / (accumulation + epsilon)))
                res_scaled = 1.d-1 * res_scaled
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                  converged_scaled_residual_real(idof,istate) = res_scaled
                  converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              endif
          end select

          if (res_scaled > this%residual_scaled_inf_tol(idof)) then
            converged_scaled = PETSC_FALSE
          endif

          if (.not.(converged_absolute .or. converged_scaled)) then
           converged_abs_residual_flag(idof,istate) = PETSC_FALSE
           converged_scaled_residual_flag(idof,istate) = PETSC_FALSE
          endif
        enddo
      enddo
    else
      do local_id = 1, grid%nlmax
        offset = (local_id-1)*option%nflowdof
        ghosted_id = grid%nL2G(local_id)
        natural_id = grid%nG2A(ghosted_id)

        sco2_auxvar = sco2_auxvars(ZERO_INTEGER,ghosted_id)

        if (patch%imat(ghosted_id) <= 0) cycle
        istate = global_auxvars(ghosted_id)%istate
        do idof = 1, option%nflowdof
          res_scaled = 0.d0
          ival = offset+idof
          converged_absolute = PETSC_TRUE
          converged_scaled = PETSC_TRUE

        ! ! infinity norms on residual
         R = dabs(r_p(ival)) 
         A = dabs(accum2_p(ival))

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
    endif
    call VecRestoreArrayReadF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(field%flow_accum2,accum2_p, &
                                ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(field%flow_dxx,dX_p,ierr);CHKERRQ(ierr)

    this%converged_flag(:,:,RESIDUAL_INDEX) = converged_abs_residual_flag(:,:)
    this%converged_real(:,:,RESIDUAL_INDEX) = converged_abs_residual_real(:,:)
    this%converged_cell(:,:,RESIDUAL_INDEX) = converged_abs_residual_cell(:,:)
    this%converged_flag(:,:,SCALED_RESIDUAL_INDEX) = &
                                       converged_scaled_residual_flag(:,:)
    this%converged_real(:,:,SCALED_RESIDUAL_INDEX) = &
                                       converged_scaled_residual_real(:,:)
    this%converged_cell(:,:,SCALED_RESIDUAL_INDEX) = &
                                       converged_scaled_residual_cell(:,:)

    flags(1:option%nflowdof*sco2_max_states*MAX_INDEX) = &
      reshape(this%converged_flag,(/option%nflowdof*sco2_max_states* &
              MAX_INDEX/))

    ! flags(option%nflowdof*sco2_max_states*MAX_INDEX+1) =&
    !   .not.sco2_high_temp_ts_cut

    mpi_int = option%nflowdof*sco2_max_states*MAX_INDEX+1
    call MPI_Allreduce(MPI_IN_PLACE,flags,mpi_int, &
                       MPI_LOGICAL,MPI_LAND,option%mycomm,ierr);CHKERRQ(ierr)

    this%converged_flag = reshape(flags(1:option%nflowdof*sco2_max_states* &
                                  MAX_INDEX),(/option%nflowdof, &
                                  sco2_max_states,MAX_INDEX/))

    ! sco2_high_temp_ts_cut = .not.flags(option%nflowdof*sco2_max_states* &
    !                                    MAX_INDEX+1)

    mpi_int = option%nflowdof*sco2_max_states*MAX_INDEX
    call MPI_Allreduce(MPI_IN_PLACE,this%converged_real,mpi_int, &
                       MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr); &
                       CHKERRQ(ierr)

    option%convergence = CONVERGENCE_CONVERGED
    do itol = 1, MAX_INDEX
      do istate = 1, sco2_max_states
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
                   trim(state_string(istate)) // ', CO2 Mass'
                ! elseif (idof == 3) then
                !   string = '   ' // trim(tol_string(itol)) // ', ' // &
                !    trim(state_string(istate)) // ', Energy'
                ! elseif (idof == 4) then
                !   string = '   ' // trim(tol_string(itol)) // ', ' // &
                !    trim(state_string(istate)) // ', Salt Mass Frac'
                ! endif
                elseif (idof == 3) then
                  string = '   ' // trim(tol_string(itol)) // ', ' // &
                   trim(state_string(istate)) // ', Salt Mass Frac'
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
        sco2_state_changed .and. &
        .not.rho_flag) then
      if (sco2_newtontrdc_hold_inner) then
        sco2_force_iteration = PETSC_TRUE
        sco2_state_changed = PETSC_FALSE
      else
        sco2_force_iteration = PETSC_TRUE
        sco2_state_changed = PETSC_FALSE
      endif
    endif

    call MPI_Allreduce(MPI_IN_PLACE,sco2_force_iteration,ONE_INTEGER, &
                       MPI_LOGICAL,MPI_LOR,option%mycomm,ierr)
    if (sco2_force_iteration) then
      if (.not.sco2_newtontrdc_hold_inner) then
        option%convergence = CONVERGENCE_BREAKOUT_INNER_ITER
        sco2_force_iteration = PETSC_FALSE
      else if (sco2_newtontrdc_hold_inner .and. &
               option%convergence == CONVERGENCE_CONVERGED) then
        option%convergence = CONVERGENCE_BREAKOUT_INNER_ITER
        sco2_force_iteration = PETSC_FALSE
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
        string = '    Exceeded SCO2 Mode Max Newton Iterations'
        call PrintMsg(option,string)
      endif
    endif
    ! if (sco2_high_temp_ts_cut) then
    !   sco2_high_temp_ts_cut = PETSC_FALSE
    !   string = '    Exceeded SCO2 Mode EOS max temperature'
    !   call PrintMsg(option,string)
    !   option%convergence = CONVERGENCE_CUT_TIMESTEP
    ! endif
    if (sco2_sub_newton_iter_num > 20) then
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

end subroutine PMSCO2CheckConvergence

! ************************************************************************** !

subroutine PMSCO2TimeCut(this)
  !
  ! Author: Michael Nole
  ! Date: 11/21/23
  !

  use SCO2_module, only : SCO2TimeCut

  implicit none

  class(pm_sco2_type) :: this

  call PMSubsurfaceFlowTimeCut(this)
  call SCO2TimeCut(this%realization)

end subroutine PMSCO2TimeCut

! ************************************************************************** !

subroutine PMSCO2UpdateSolution(this)
  !
  ! Author: Michael Nole
  ! Date: 11/21/23
  !

  use SCO2_module, only : SCO2UpdateSolution, &
                          SCO2MapBCAuxVarsToGlobal

  implicit none

  class(pm_sco2_type) :: this

  call PMSubsurfaceFlowUpdateSolution(this)
  call SCO2UpdateSolution(this%realization)
  call SCO2MapBCAuxVarsToGlobal(this%realization)

end subroutine PMSCO2UpdateSolution

! ************************************************************************** !

subroutine PMSCO2UpdateAuxVars(this)
  !
  ! Author: Michael Nole
  ! Date: 11/21/23
  use SCO2_module, only : SCO2UpdateAuxVars

  implicit none

  class(pm_sco2_type) :: this

  call SCO2UpdateAuxVars(this%realization,PETSC_FALSE,PETSC_TRUE)

end subroutine PMSCO2UpdateAuxVars

! ************************************************************************** !

subroutine PMSCO2MaxChange(this)
  !
  ! Author: Michael Nole
  ! Date: 11/21/23
  !

  use Realization_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use Global_Aux_module
  use SCO2_Aux_module

  implicit none

  class(pm_sco2_type) :: this

  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscReal, pointer :: max_change_local(:)
  PetscReal, pointer :: max_change_global(:)
  PetscReal :: max_change
  PetscInt :: i, j
  PetscInt :: ghosted_id

  PetscErrorCode :: ierr

  realization => this%realization
  option => realization%option
  field => realization%field
  grid => realization%patch%grid

  global_auxvars => realization%patch%aux%global%auxvars

  allocate(max_change_local(max_change_index))
  allocate(max_change_global(max_change_index))
  max_change_global = 0.d0
  max_change_local = 0.d0

  do i = 1, max_change_index
    call RealizationGetVariable(realization,field%work, &
                                this%max_change_ivar(i), &
                                this%max_change_isubvar(i))
    call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%max_change_vecs(i),vec_ptr2,ierr);CHKERRQ(ierr)
    max_change = 0.d0
    if (i==1 .and. sco2_chk_max_dpl_liq_state_only) then
      do j = 1,grid%nlmax
        ghosted_id = grid%nL2G(j)
        if (global_auxvars(ghosted_id)%istate /= SCO2_LIQUID_STATE) then
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
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr);&
                      CHKERRQ(ierr)
  ! print them out
  if (OptionPrintToScreen(option)) then
    write(*,'("  --> max chng: dpl= ",1pe12.4, " dpg= ",1pe12.4,&
      & " dpco2= ",1pe12.4,/,15x," dxco2= ",1pe12.4," dxs= ",1pe12.4,&
      & " dsg= ",1pe12.4)') &
      max_change_global(1:max_change_index)
  endif

  if (OptionPrintToFile(option)) then
    write(option%fid_out,'("  --> max chng: dpl= ",1pe12.4," dpg= ",1pe12.4,&
    & " dpco2= ",1pe12.4,/,15x," dxco2= ",1pe12.4," dxs= ",1pe12.4, &
    & " dsg= ",1pe12.4)') &
      max_change_global(1:max_change_index)
  endif
  ! if (OptionPrintToScreen(option)) then
  !   write(*,'("  --> max chng: dpl= ",1pe12.4, " dpg= ",1pe12.4,&
  !     & " dpco2= ",1pe12.4,/,15x," dxco2= ",1pe12.4," dxs= ",1pe12.4,&
  !     & "dt= ",1pe12.4," dsg= ",1pe12.4)') &
  !     max_change_global(1:max_change_index)
  ! endif

  ! if (OptionPrintToFile(option)) then
  !   write(option%fid_out,'("  --> max chng: dpl= ",1pe12.4," dpg= ",1pe12.4,&
  !   & " dpco2= ",1pe12.4,/,15x," dxco2= ",1pe12.4," dxs= ",1pe12.4, &
  !   & "dt= ",1pe12.4," dsg= ",1pe12.4)') &
  !     max_change_global(1:max_change_index)
  ! endif

  ! MAN: check these
  this%max_pressure_change = maxval(max_change_global(1:3))
  this%max_xmol_change = max_change_global(4)
  this%max_salt_mass_change = max_change_global(5)
  this%max_saturation_change = max_change_global(6)
  ! this%max_temperature_change = max_change_global(6)
  ! this%max_saturation_change = max_change_global(7)

end subroutine PMSCO2MaxChange

! ************************************************************************** !

subroutine PMSCO2ComputeMassBalance(this,mass_balance_array)
  !
  ! Author: Michael Nole
  ! Date: 11/21/23
  !

  use SCO2_module, only : SCO2ComputeMassBalance

  implicit none

  class(pm_sco2_type) :: this
  PetscReal :: mass_balance_array(:)
  PetscReal :: mass_trapped(this%realization%option%nphase)

  call SCO2ComputeMassBalance(this%realization,mass_balance_array,mass_trapped)

end subroutine PMSCO2ComputeMassBalance

! ************************************************************************** !

subroutine PMSCO2InputRecord(this)
  !
  ! Author: Michael Nole
  ! Date: 11/21/23
  !

  implicit none

  class(pm_sco2_type) :: this

end subroutine PMSCO2InputRecord

! ************************************************************************** !

subroutine PMSCO2CheckpointBinary(this,viewer)
  !
  ! Checkpoints data associated with SCO2 PM
  !
  ! Author: Michael Nole
  ! Date: 11/21/23

  use Checkpoint_module
  use Global_module

  implicit none

#include "petsc/finclude/petscviewer.h"

  class(pm_sco2_type) :: this
  PetscViewer :: viewer

  call PMSubsurfaceFlowCheckpointBinary(this,viewer)

end subroutine PMSCO2CheckpointBinary

! ************************************************************************** !

subroutine PMSCO2RestartBinary(this,viewer)
  !
  ! Restarts data associated with SCO2 PM
  !
  ! Author: Michael Nole
  ! Date: 11/21/23

  use Checkpoint_module
  use Global_module

  implicit none

#include "petsc/finclude/petscviewer.h"

  class(pm_sco2_type) :: this
  PetscViewer :: viewer

  call PMSubsurfaceFlowRestartBinary(this,viewer)

end subroutine PMSCO2RestartBinary

! ************************************************************************** !

subroutine PMSCO2Destroy(this)
  !
  ! Destroys SCO2 process model
  !
  ! Author: Michael Nole
  ! Date: 11/21/23
  !

  use SCO2_module, only : SCO2Destroy

  implicit none

  class(pm_sco2_type) :: this

  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  deallocate(this%max_change_ivar)
  nullify(this%max_change_ivar)
  deallocate(this%max_change_isubvar)
  nullify(this%max_change_isubvar)

  ! preserve this ordering
  call SCO2Destroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)

end subroutine PMSCO2Destroy

end module PM_SCO2_class
