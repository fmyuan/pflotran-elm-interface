module PM_Hydrate_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class
  use Hydrate_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscInt, parameter :: ABS_UPDATE_INDEX = 1
  PetscInt, parameter :: REL_UPDATE_INDEX = 2
  PetscInt, parameter :: RESIDUAL_INDEX = 3
  PetscInt, parameter :: SCALED_RESIDUAL_INDEX = 4
  PetscInt, parameter :: MAX_INDEX = SCALED_RESIDUAL_INDEX

  type, public, extends(pm_subsurface_flow_type) :: pm_hydrate_type
    PetscInt, pointer :: max_change_ivar(:)
    PetscInt, pointer :: max_change_isubvar(:)
    type(hydrate_parameter_type), pointer :: hydrate_parameters
    PetscBool :: converged_flag(3,15,MAX_INDEX)
    PetscInt :: converged_cell(3,15,MAX_INDEX)
    PetscReal :: converged_real(3,15,MAX_INDEX)
    PetscReal :: residual_abs_inf_tol(3)
    PetscReal :: residual_scaled_inf_tol(3)
    PetscReal :: abs_update_inf_tol(3,15)
    PetscReal :: rel_update_inf_tol(3,15)
    PetscReal :: damping_factor
  contains
    procedure, public :: ReadSimulationOptionsBlock => &
                           PMHydrateReadSimOptionsBlock
    procedure, public :: ReadNewtonBlock => PMHydrateReadNewtonSelectCase
    procedure, public :: InitializeSolver => PMHydrateInitializeSolver
    procedure, public :: Setup => PMHydrateSetup
    procedure, public :: InitializeRun => PMHydrateInitializeRun
    procedure, public :: InitializeTimestep => PMHydrateInitializeTimestep
    procedure, public :: Residual => PMHydrateResidual
    procedure, public :: Jacobian => PMHydrateJacobian
    procedure, public :: UpdateTimestep => PMHydrateUpdateTimestep
    procedure, public :: PreSolve => PMHydratePreSolve
    procedure, public :: PostSolve => PMHydratePostSolve
    procedure, public :: CheckUpdatePre => PMHydrateCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMHydrateCheckUpdatePost
    procedure, public :: CheckConvergence => PMHydrateCheckConvergence
    procedure, public :: TimeCut => PMHydrateTimeCut
    procedure, public :: UpdateSolution => PMHydrateUpdateSolution
    procedure, public :: UpdateAuxVars => PMHydrateUpdateAuxVars
    procedure, public :: MaxChange => PMHydrateMaxChange
    procedure, public :: ComputeMassBalance => PMHydrateComputeMassBalance
    procedure, public :: InputRecord => PMHydrateInputRecord
    procedure, public :: CheckpointBinary => PMHydrateCheckpointBinary
    procedure, public :: RestartBinary => PMHydrateRestartBinary
    procedure, public :: Destroy => PMHydrateDestroy
  end type pm_hydrate_type

  public :: PMHydrateCreate, &
            PMHydrateSetFlowMode, &
            PMHydrateReadParameters, &
            PMHydrateAssignParameters

contains

! ************************************************************************** !

function PMHydrateCreate()
  !
  ! Creates Hydrate process models shell
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !
  use Variables_module, only : LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
                               LIQUID_MOLE_FRACTION, TEMPERATURE, &
                               GAS_SATURATION, HYDRATE_SATURATION, &
                               LIQUID_SATURATION, ICE_SATURATION
  use Upwind_Direction_module

  implicit none

  class(pm_hydrate_type), pointer :: PMHydrateCreate

  class(pm_hydrate_type), pointer :: this

  !MAN optimized:
  PetscReal, parameter :: pres_abs_inf_tol = 1.d0 ! Reference tolerance [Pa]
  PetscReal, parameter :: temp_abs_inf_tol = 1.d-5
  PetscReal, parameter :: sat_abs_inf_tol = 1.d-5
  PetscReal, parameter :: xmol_abs_inf_tol = 1.d-9

  PetscReal, parameter :: pres_rel_inf_tol = 1.d-3
  PetscReal, parameter :: temp_rel_inf_tol = 1.d-3
  PetscReal, parameter :: sat_rel_inf_tol = 1.d-3
  PetscReal, parameter :: xmol_rel_inf_tol = 1.d-3

  !MAN optimized:
  PetscReal, parameter :: w_mass_abs_inf_tol = 1.d-5 !1.d-7 !kmol_water/sec
  PetscReal, parameter :: a_mass_abs_inf_tol = 1.d-5 !1.d-7
  PetscReal, parameter :: u_abs_inf_tol = 1.d-5 !1.d-7

  PetscReal, parameter :: residual_abs_inf_tol(3) = (/w_mass_abs_inf_tol, &
                             a_mass_abs_inf_tol, u_abs_inf_tol/)
  PetscReal, parameter :: residual_scaled_inf_tol(3) = 1.d-6

  !For convergence using hydrate and ice formation capability
  PetscReal, parameter :: abs_update_inf_tol(3,15) = &
             !L_STATE
    reshape([pres_abs_inf_tol,xmol_abs_inf_tol,temp_abs_inf_tol, &
             !G_STATE
             pres_abs_inf_tol,pres_abs_inf_tol,temp_abs_inf_tol, &
             !H_STATE
             pres_abs_inf_tol,999.d0,temp_abs_inf_tol,  &
             !I_STATE
             pres_abs_inf_tol,999.d0,temp_abs_inf_tol, &
             !GA_STATE
             pres_abs_inf_tol,sat_abs_inf_tol,temp_abs_inf_tol, &
             !HG_STATE
             pres_abs_inf_tol,sat_abs_inf_tol,temp_abs_inf_tol, &
             !HA_STATE
             pres_abs_inf_tol,sat_abs_inf_tol,temp_abs_inf_tol, &
             !HI_STATE
             pres_abs_inf_tol,sat_abs_inf_tol,temp_abs_inf_tol, &
             !GI_STATE
             pres_abs_inf_tol,sat_abs_inf_tol,temp_abs_inf_tol, &
             !AI_STATE
             pres_abs_inf_tol,xmol_abs_inf_tol,sat_abs_inf_tol, &
             !HGA_STATE
             sat_abs_inf_tol,sat_abs_inf_tol,temp_abs_inf_tol, &
             !HAI_STATE
             pres_abs_inf_tol,sat_abs_inf_tol,sat_abs_inf_tol, &
             !HGI_STATE
             sat_abs_inf_tol,sat_abs_inf_tol,temp_abs_inf_tol, &
             !GAI_STATE
             pres_abs_inf_tol,sat_abs_inf_tol,temp_abs_inf_tol, &
             !HGAI_STATE
             sat_abs_inf_tol,sat_abs_inf_tol,sat_abs_inf_tol], &
            shape(abs_update_inf_tol)) * &
            1.d0 ! change to 0.d0 to zero tolerances
  PetscReal, parameter :: rel_update_inf_tol(3,15) = &
             !L_STATE
    reshape([pres_rel_inf_tol,xmol_rel_inf_tol,temp_rel_inf_tol, &
             !G_STATE
             pres_rel_inf_tol,pres_rel_inf_tol,temp_rel_inf_tol, &
             !H_STATE
             pres_rel_inf_tol,999.d0,temp_rel_inf_tol, &
             !I_STATE
             pres_rel_inf_tol,999.d0,temp_rel_inf_tol, &
             !GA_STATE
             pres_rel_inf_tol,sat_rel_inf_tol,temp_rel_inf_tol, &
             !HG_STATE
             pres_rel_inf_tol,sat_rel_inf_tol,temp_rel_inf_tol, &
             !HA_STATE
             pres_rel_inf_tol,sat_rel_inf_tol,temp_rel_inf_tol, &
             !HI_STATE
             pres_rel_inf_tol,sat_rel_inf_tol,temp_rel_inf_tol, &
             !GI_STATE
             pres_rel_inf_tol,sat_rel_inf_tol,temp_rel_inf_tol, &
             !AI_STATE
             pres_rel_inf_tol,xmol_rel_inf_tol,sat_rel_inf_tol, &
             !HGA_STATE
             sat_rel_inf_tol,sat_rel_inf_tol,temp_rel_inf_tol, &
             !HAI_STATE
             pres_rel_inf_tol,sat_rel_inf_tol,sat_rel_inf_tol, &
             !HGI_STATE
             sat_rel_inf_tol,sat_rel_inf_tol,temp_rel_inf_tol, &
             !GAI_STATE
             pres_rel_inf_tol,sat_rel_inf_tol,temp_rel_inf_tol, &
             !HGAI_STATE
             sat_rel_inf_tol,sat_rel_inf_tol,sat_rel_inf_tol], &
            shape(rel_update_inf_tol)) * &
            1.d0 ! change to 0.d0 to zero tolerances
  allocate(this)
  allocate(this%max_change_ivar(9))
  this%max_change_ivar = [LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
                                LIQUID_MOLE_FRACTION, TEMPERATURE, &
                                GAS_SATURATION, HYDRATE_SATURATION, &
                                LIQUID_SATURATION, ICE_SATURATION]
  allocate(this%max_change_isubvar(9))
                                   ! UNINITIALIZED_INTEGER avoids zeroing of
                                   ! pressures not represented in phase
                                       ! 2 = air in xmol(air,liquid)
  this%max_change_isubvar = [0,0,0,2,0,0,0,0,0]
  this%damping_factor = -1.d0

  call PMSubsurfaceFlowInit(this)
  this%name = 'Hydrate Multiphase Flow'
  this%header = 'GAS HYDRATE + MULTIPHASE FLOW'

  ! turn off default upwinding which is set to PETSC_TRUE in
  !  upwind_direction.F90
  fix_upwind_direction = PETSC_FALSE

! this should be set explicitly in input file using
! USE_INFINITY_NORM_CONVERGENCE specified in input block
!  this%check_post_convergence = PETSC_TRUE
  this%residual_abs_inf_tol = residual_abs_inf_tol
  this%residual_scaled_inf_tol = residual_scaled_inf_tol
  this%abs_update_inf_tol = abs_update_inf_tol
  this%rel_update_inf_tol = rel_update_inf_tol

  this%converged_flag(:,:,:) = PETSC_TRUE
  this%converged_real(:,:,:) = 0.d0
  this%converged_cell(:,:,:) = 0

  allocate(this%hydrate_parameters)
  nullify(this%hydrate_parameters%methanogenesis)

  PMHydrateCreate => this

end function PMHydrateCreate

! ************************************************************************** !

subroutine PMHydrateSetFlowMode(option)
!
! Sets the flow mode for equilibrium hydrate formation and dissociation
!
! Author: Michael Nole
! Date: 01/02/19
!

  use Option_module

  implicit none

  type(option_type) :: option

  option%iflowmode = H_MODE
  option%nphase = 5
  option%liquid_phase = 1  ! liquid_pressure
  option%gas_phase = 2     ! gas_pressure
  option%hydrate_phase = 3
  option%ice_phase = 4
  option%precipitate_phase = 5
  option%trapped_gas_phase = 6

  option%air_pressure_id = 4
  option%capillary_pressure_id = 5
  option%vapor_pressure_id = 6
  option%saturation_pressure_id = 7
  option%reduced_vapor_pressure_id = 8
  
  option%pure_water_phase = 5
  option%pure_brine_phase = 6
  

  option%water_id = 1
  option%air_id = 2
  option%energy_id = 3
  option%salt_id = 3 ! Salt component

  option%nflowdof = 3
  option%nflowspec = 3
  option%use_isothermal = PETSC_FALSE

end subroutine PMHydrateSetFlowMode

! ************************************************************************** !

subroutine PMHydrateReadParameters(input,pm_hydrate,option)

  use Input_Aux_module
  use Option_module
  use String_module

  implicit none

  type(input_type), pointer :: input
  class(pm_hydrate_type) :: pm_hydrate
  type(option_type), pointer :: option

  type(methanogenesis_type), pointer :: methanogenesis
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: temp_int

  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)
    if (input%ierr /= 0) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword','HYDRATE')
    call StringToUpper(word)

    select case(trim(word))
      case('NO_SOLID_SATURATION_PERM_SCALING')
        ! This turns of scaling of the intrinsic permeability
        ! as a function of (ice + hydrate) saturations
        hydrate_perm_scaling = PETSC_FALSE
      case('HYDRATE_PHASE_BOUNDARY')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'keyword','hydrate phase boundary')
        call StringToUpper(word)
        select case(word)
          case('MORIDIS')
            hydrate_phase_boundary = 2
          case('MORIDIS_SIMPLE')
            hydrate_phase_boundary = 3
          case default
            call InputKeywordUnrecognized(input,word,&
                 'HYDRATE_PHASE_BOUNDARY',option)
        end select
      case('HENRYS_CONSTANT')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'keyword','hydrate henrys constant')
        call StringToUpper(word)
        select case(word)
          case('CRAMER')
            hydrate_henrys_constant = 2
          case('CO2')
            hydrate_use_henry_co2 = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,word,&
                 'HYDRATE_HENRYS_CONSTANT',option)
        end select
      case('GAS')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'keyword','hydrate mode gas')
        call StringToUpper(word)
        select case(word)
          case('AIR')
            hydrate_former = HYDRATE_FORMER_NULL
            hydrate_fmw_comp(2) = FMWAIR
          case('METHANE','CH4')
            hydrate_former = HYDRATE_FORMER_CH4
            hydrate_fmw_comp(2) = FMWCH4
          case('CO2')
            hydrate_former = HYDRATE_FORMER_CO2
            hydrate_fmw_comp(2) = FMWCO2
          case default
            call InputKeywordUnrecognized(input,word,&
                 'HYDRATE_GAS',option)
        end select
      case('NO_ICE_VOLUME_CHANGE')
        ! This sets the ice density to the water density (no volume change)
        hydrate_no_ice_density_change = PETSC_TRUE
      case('NO_EFFECTIVE_SATURATION_SCALING')
        ! This turns off normalizing the liquid and gas saturations by the
        ! sum of mobile phases when computing relative permeabilities.
        hydrate_eff_sat_scaling = PETSC_FALSE
      case('WITH_GIBBS_THOMSON')
        ! Scales methane solubility as a function of pore size.
        hydrate_with_gibbs_thomson = PETSC_TRUE
      case('GT_3PHASE')
        hydrate_gt_3phase = PETSC_TRUE
      case('ADJUST_SOLUBILITY_WITHIN_GHSZ')
        hydrate_adjust_ghsz_solubility = PETSC_TRUE
      case('WITH_SEDIMENTATION')
        hydrate_with_sedimentation = PETSC_TRUE
      case('NO_PC')
        hydrate_no_pc = PETSC_TRUE
      case('METHANOGENESIS')
        hydrate_with_methanogenesis = PETSC_TRUE
        if (.not. associated(pm_hydrate%hydrate_parameters%methanogenesis)) then
          pm_hydrate%hydrate_parameters%methanogenesis => &
                     HydrateMethanogenesisCreate()
        endif
        methanogenesis => pm_hydrate%hydrate_parameters%methanogenesis
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          if (input%ierr /= 0) exit
          if (InputCheckExit(input,option)) exit

          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword','HYDRATE')
          call StringToUpper(word)
          select case(trim(word))
            case('NAME')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'methanogenesis source name', &
                                 error_string)
              call StringToUpper(word)
              methanogenesis%source_name = trim(word)
            case('ALPHA')
              call InputReadDouble(input,option,methanogenesis%alpha)
              call InputErrorMsg(input,option,'alpha',error_string)
            case('LAMBDA')
              call InputReadDouble(input,option,methanogenesis%lambda)
              call InputErrorMsg(input,option,'lambda',error_string)
            case('V_SED')
              call InputReadDouble(input,option,methanogenesis%omega)
              call InputErrorMsg(input,option,'v_sed',error_string)
            case('SMT_DEPTH')
              call InputReadDouble(input,option,methanogenesis%z_smt)
              call InputErrorMsg(input,option,'smt_depth',error_string)
            case('K_ALPHA')
              call InputReadDouble(input,option,methanogenesis%k_alpha)
              call InputErrorMsg(input,option,'k_alpha',error_string)
          end select
        enddo
        call InputPopBlock(input,option)
      case('PERM_SCALING_FUNCTION')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'keyword','hyd_perm_scaling_function')
        call StringToUpper(word)
        select case(trim(word))
          case('DAI_AND_SEOL')
            temp_int = 1
            hydrate_perm_scaling_function = temp_int
          case default
            call InputKeywordUnrecognized(input,word,&
                 'PERM_SCALING_FUNCTION',option)
        end select
      case('SALINITY')
        call InputReadDouble(input,option,hydrate_xmass_nacl)
        call InputErrorMsg(input,option,'SALINITY',error_string)
      case('THERMAL_CONDUCTIVITY')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'keyword','thermal_conductivity_func')
        call StringToUpper(word)
        select case(trim(word))
          case('IGHCC2')
            hydrate_tcond = 1
          case default
            call InputKeywordUnrecognized(input,word,&
                 'HYDRATE COMPOSITE THERMAL CONDUCTIVITY MODEL',option)
        end select
      case('BC_REFERENCE_PRESSURE')
        call InputReadDouble(input,option,hydrate_bc_reference_pressure)
        call InputErrorMsg(input,option,'BC Reference Pressure',error_string)
      case default
        call InputKeywordUnrecognized(input,word,'HYDRATE block',option)
    end select

  enddo
  call InputPopBlock(input,option)

end subroutine PMHydrateReadParameters

! ************************************************************************** !

subroutine PMHydrateAssignParameters(realization, pm)

  ! Points hydrate parameters to those read into the PM.
  !
  ! Author: Michael Nole
  ! Date: 11/20/19
  !

  use Realization_Subsurface_class
  use Hydrate_Aux_module
  use Fluid_module
  use Option_module

  implicit none

  class(realization_subsurface_type), pointer :: realization
  class(pm_hydrate_type) :: pm

  type(fluid_property_type), pointer :: cur_fluid_property
  type(option_type), pointer :: option

  option => realization%option

  ! initialize parameters
  allocate(pm%hydrate_parameters%diffusion_coefficient(option%nflowspec, &
                                                       option%nphase))
  cur_fluid_property => realization%fluid_properties
  do
    if (.not.associated(cur_fluid_property)) exit
    pm%hydrate_parameters% &
      diffusion_coefficient(:,cur_fluid_property%phase_id) = &
        cur_fluid_property%diffusion_coefficient
    cur_fluid_property => cur_fluid_property%next
  enddo

  realization%patch%aux%hydrate%hydrate_parameter => pm%hydrate_parameters

end subroutine PMHydrateAssignParameters

! ************************************************************************** !

subroutine PMHydrateReadSimOptionsBlock(this,input)
  !
  ! Sets up SNES solvers.
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !
  use Hydrate_Aux_module
  use Input_Aux_module
  use String_module
  use Option_module


  implicit none

  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: keyword, word
  class(pm_hydrate_type) :: this
  type(option_type), pointer :: option
  PetscReal :: tempreal
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
  PetscInt :: lid, gid, eid

  option => this%option

  lid = 1 !option%liquid_phase
  gid = 2 !option%gas_phase
  eid = 3 !option%energy_id

  error_string = 'Hydrate Options'

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
      case('ARITHMETIC_GAS_DIFFUSIVE_DENSITY')
        hydrate_harmonic_diff_density = PETSC_FALSE
      case('CHECK_MAX_DPL_LIQ_STATE_ONLY')
        hyd_chk_max_dpl_liq_state_only = PETSC_TRUE
      case('DEBUG_CELL')
        call InputReadInt(input,option,hydrate_debug_cell_id)
        call InputErrorMsg(input,option,keyword,error_string)
      case('DIFFUSE_XMASS')
        hydrate_diffuse_xmol = PETSC_FALSE
      case('GAS_COMPONENT_FORMULA_WEIGHT')
        !geh: assuming gas component is index 2
        call InputReadDouble(input,option,hydrate_fmw_comp(2))
        call InputErrorMsg(input,option,keyword,error_string)
      case('HARMONIC_GAS_DIFFUSIVE_DENSITY')
        hydrate_harmonic_diff_density = PETSC_TRUE
      case('NEWTONTRDC_HOLD_INNER_ITERATIONS',&
           'HOLD_INNER_ITERATIONS','NEWTONTRDC_HOLD_INNER')
        !heeho: only used when using newtontrd-c
        hydrate_newtontrdc_hold_inner = PETSC_TRUE
      case('IMMISCIBLE')
        hydrate_immiscible = PETSC_TRUE
      case('LIQUID_COMPONENT_FORMULA_WEIGHT')
        !heeho: assuming liquid component is index 1
        call InputReadDouble(input,option,hydrate_fmw_comp(1))
        call InputErrorMsg(input,option,keyword,error_string)
      case('NO_STATE_TRANSITION_OUTPUT')
        hydrate_print_state_transition = PETSC_FALSE
      case('NO_TEMP_DEPENDENT_DIFFUSION')
        hydrate_temp_dep_gas_air_diff = PETSC_FALSE
      case('PHASE_CHANGE_EPSILON')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,error_string)
        hydrate_phase_chng_epsilon = tempreal
      case('RESTRICT_STATE_CHANGE')
        hydrate_restrict_state_chng = PETSC_TRUE
      case('TWO_PHASE_ENERGY_DOF')
        call InputKeywordDeprecated('TWO_PHASE_ENERGY_DOF', &
                                    'TWO_PHASE_STATE_ENERGY_DOF',option)
      case('TWO_PHASE_STATE_ENERGY_DOF')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,keyword,error_string)
        call HydrateAuxSetEnergyDOF(word,option)
      case('WINDOW_EPSILON')
        call InputReadDouble(input,option,window_epsilon)
        call InputErrorMsg(input,option,keyword,error_string)
      case('CALCULATE_SURFACE_TENSION')
        hydrate_compute_surface_tension = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(input,keyword,'HYDRATE Mode',option)
    end select

  enddo
  call InputPopBlock(input,option)

end subroutine PMHydrateReadSimOptionsBlock

! ************************************************************************** !

subroutine PMHydrateReadNewtonSelectCase(this,input,keyword,found, &
                                         error_string,option)
  !
  ! Reads input file parameters associated with the HYDRATE process model
  ! Newton solver convergence
  !
  ! Author: Glenn Hammond
  ! Date: 03/23/20

  use Input_Aux_module
  use String_module
  use Utility_module
  use Option_module
  use Hydrate_Aux_module

  implicit none

  class(pm_hydrate_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option

  PetscBool :: found
  PetscReal :: tempreal
  PetscInt :: lid, gid, eid

  option => this%option

  lid = 1 !option%liquid_phase
  gid = 2 !option%gas_phase
  eid = 3 !option%energy_id

  error_string = 'HYDRATE Newton Solver'

  found = PETSC_FALSE
  call PMSubsurfaceFlowReadNewtonSelectCase(this,input,keyword,found, &
                                            error_string,option)
  if (found) return

  found = PETSC_TRUE
  select case(trim(keyword))
      case('USE_GOVERNORS')
        hydrate_use_governors = PETSC_TRUE
      case('CHECK_SOLUTION_UPDATES')
        hydrate_check_updates = PETSC_TRUE
      case('USE_FULL_CONVERGENCE_CRITERIA')
        hydrate_full_convergence = PETSC_TRUE
      case('NO_UPDATE_TRUNCATION')
        hydrate_truncate_updates = PETSC_FALSE
      case('CENTRAL_DIFFERENCE_JACOBIAN')
        hydrate_central_diff_jacobian = PETSC_TRUE
      case('HYDRATE_UPDATE_INF_TOL')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,error_string)
        this%abs_update_inf_tol(2,3) = tempreal
        this%abs_update_inf_tol(2,6) = tempreal
        this%abs_update_inf_tol(2,7) = tempreal
        this%abs_update_inf_tol(2,8) = tempreal
        this%abs_update_inf_tol(2,9) = tempreal
        this%abs_update_inf_tol(3,10) = tempreal
        this%abs_update_inf_tol(1:2,11) = tempreal
        this%abs_update_inf_tol(2:3,12) = tempreal
        this%abs_update_inf_tol(1:2,13) = tempreal
        this%abs_update_inf_tol(2:3,14) = tempreal
        this%abs_update_inf_tol(:,15) = tempreal
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
        this%abs_update_inf_tol(1,1:10) = tempreal
        this%abs_update_inf_tol(2,2) = tempreal
        this%abs_update_inf_tol(1,12) = tempreal
        this%abs_update_inf_tol(1,14) = tempreal
      case('TEMP_ABS_UPDATE_INF_TOL')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,error_string)
        this%abs_update_inf_tol(3,1:9) = tempreal
        this%abs_update_inf_tol(3,11) = tempreal
        this%abs_update_inf_tol(3,13) = tempreal
      case('SAT_ABS_UPDATE_INF_TOL')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,error_string)
        this%abs_update_inf_tol(2,3) = tempreal
        this%abs_update_inf_tol(2,5:9) = tempreal
        this%abs_update_inf_tol(2,11:15) = tempreal
        this%abs_update_inf_tol(3,10) = tempreal
        this%abs_update_inf_tol(3,12) = tempreal
        this%abs_update_inf_tol(3,14:15) = tempreal
        this%abs_update_inf_tol(1,11) = tempreal
        this%abs_update_inf_tol(1,13) = tempreal
        this%abs_update_inf_tol(1,15) = tempreal
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
      case('REL_UPDATE_INF_TOL','ITOL_RELATIVE_UPDATE')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,error_string)
        this%rel_update_inf_tol(:,:) = tempreal
      case('PRES_REL_UPDATE_INF_TOL')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,error_string)
        this%rel_update_inf_tol(1,:) = tempreal
        this%rel_update_inf_tol(2,2) = tempreal
        this%rel_update_inf_tol(1,1:10) = tempreal
        this%rel_update_inf_tol(1,12) = tempreal
        this%rel_update_inf_tol(1,14) = tempreal
        this%rel_update_inf_tol(2,2) = tempreal
      case('TEMP_REL_UPDATE_INF_TOL')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,error_string)
        this%rel_update_inf_tol(3,:) = tempreal
        this%rel_update_inf_tol(3,1:9) = tempreal
        this%rel_update_inf_tol(3,11) = tempreal
        this%rel_update_inf_tol(3,13) = tempreal
      case('SAT_REL_UPDATE_INF_TOL')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,error_string)
        this%rel_update_inf_tol(2,3) = tempreal
        this%rel_update_inf_tol(2,3) = tempreal
        this%rel_update_inf_tol(2,6:9) = tempreal
        this%rel_update_inf_tol(3,10) = tempreal
        this%rel_update_inf_tol(2,11:15) = tempreal
        this%rel_update_inf_tol(3,12) = tempreal
        this%rel_update_inf_tol(3,14:15) = tempreal
        this%rel_update_inf_tol(1,11) = tempreal
        this%rel_update_inf_tol(1,13) = tempreal
        this%rel_update_inf_tol(1,15) = tempreal
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
        call InputReadDouble(input,option,hydrate_max_pressure_change)
        call InputErrorMsg(input,option,keyword,error_string)
      case('MAX_ITERATION_BEFORE_DAMPING')
        call InputReadInt(input,option,hydrate_max_it_before_damping)
        call InputErrorMsg(input,option,keyword,error_string)
      case('DAMPING_FACTOR')
        call InputReadDouble(input,option,hydrate_damping_factor)
        call InputErrorMsg(input,option,keyword,error_string)
        this%damping_factor = hydrate_damping_factor
    case default
      found = PETSC_FALSE

  end select

end subroutine PMHydrateReadNewtonSelectCase

! ************************************************************************** !

subroutine PMHydrateInitializeSolver(this)
  !
  ! Author: Glenn Hammond
  ! Date: 04/06/20

  use Solver_module

  implicit none

  class(pm_hydrate_type) :: this

  call PMBaseInitializeSolver(this)

  ! helps accommodate rise in residual due to change in state
  this%solver%newton_dtol = 1.d9
  this%solver%newton_max_iterations = 16

end subroutine PMHydrateInitializeSolver

! ************************************************************************** !

subroutine PMHydrateSetup(this)
  !
  ! Sets up auxvars and parameters
  !
  ! Author: Glenn Hammond
  ! Date: 04/11/24

  use Hydrate_module
  use Material_module

  implicit none

  class(pm_hydrate_type) :: this

  call this%SetRealization()
  call MaterialSetupThermal( &
         this%realization%patch%aux%Material%material_parameter, &
         this%realization%patch%material_property_array, &
         this%realization%option)
  call HydrateSetup(this%realization)
  call PMHydrateAssignParameters(this%realization,this)
  call PMSubsurfaceFlowSetup(this)

end subroutine PMHydrateSetup

! ************************************************************************** !

recursive subroutine PMHydrateInitializeRun(this)
  !
  ! Initializes the time stepping
  !
  ! Author: Michael Nole
  ! Date: 07/23/19

  use Realization_Base_class

  implicit none

  class(pm_hydrate_type) :: this

  PetscInt :: i
  PetscErrorCode :: ierr

  ! need to allocate vectors for max change
  call VecDuplicateVecsF90(this%realization%field%work,NINE_INTEGER, &
                           this%realization%field%max_change_vecs, &
                           ierr);CHKERRQ(ierr)
  ! set initial values
  do i = 1, 9
    call RealizationGetVariable(this%realization, &
                                this%realization%field%max_change_vecs(i), &
                                this%max_change_ivar(i), &
                                this%max_change_isubvar(i))
  enddo

  ! call parent implementation
  call PMSubsurfaceFlowInitializeRun(this)

end subroutine PMHydrateInitializeRun

! ************************************************************************** !

subroutine PMHydrateInitializeTimestep(this)
  !
  ! Should not need this as it is called in PreSolve.
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !

  use Hydrate_module, only : HydrateInitializeTimestep
  use Global_module
  use Variables_module, only : TORTUOSITY
  use Material_module, only : MaterialAuxVarCommunicate
  use Option_module

  implicit none

  class(pm_hydrate_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)
!geh:remove   everywhere
  call MaterialAuxVarCommunicate(this%comm1, &
                                 this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc,TORTUOSITY, &
                                 ZERO_INTEGER)

  call HydrateInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)

end subroutine PMHydrateInitializeTimestep

! ************************************************************************** !

subroutine PMHydratePreSolve(this)
  !
  ! Author: Michael Nole
  ! Date: 07/23/19

  implicit none

  class(pm_hydrate_type) :: this

  call PMSubsurfaceFlowPreSolve(this)

end subroutine PMHydratePreSolve

! ************************************************************************** !

subroutine PMHydratePostSolve(this)
  !
  ! Author: Michael Nole
  ! Date: 07/23/19

  implicit none

  class(pm_hydrate_type) :: this

end subroutine PMHydratePostSolve

! ************************************************************************** !

subroutine PMHydrateUpdateTimestep(this,update_dt, &
                                   dt,dt_min,dt_max,iacceleration, &
                                   num_newton_iterations,tfac, &
                                   time_step_max_growth_factor)
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
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

  class(pm_hydrate_type) :: this
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
    if (hydrate_use_governors) then
      governed_dt = umin_scale * dt
      dtt = min(time_step_max_growth_factor*dt,governed_dt)
      dt = min(dtt,tfac(ifac)*dt,dt_max)
      dt = max(dt,dt_min)
    else
      dtt = time_step_max_growth_factor*dt
      dt = min(dtt,dt_max)
      dt = max(dt,dt_min)
    endif

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
    ! Since saturations are not stored in global_auxvar for hydrate mode, we
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

end subroutine PMHydrateUpdateTimestep

! ************************************************************************** !

subroutine PMHydrateResidual(this,snes,xx,r,ierr)
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !

  use Hydrate_module, only : HydrateResidual

  implicit none

  class(pm_hydrate_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr

  call PMSubsurfaceFlowUpdatePropertiesNI(this)
  call HydrateResidual(snes,xx,r,this%realization,ierr)

end subroutine PMHydrateResidual

! ************************************************************************** !

subroutine PMHydrateJacobian(this,snes,xx,A,B,ierr)
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !

  use Hydrate_module, only : HydrateJacobian

  implicit none

  class(pm_hydrate_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr

  call HydrateJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMHydrateJacobian

! ************************************************************************** !

subroutine PMHydrateCheckUpdatePre(this,snes,X,dX,changed,ierr)
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !

  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Option_module
  use Saturation_Function_module
  use Patch_module
  use Hydrate_Aux_module
  use Global_Aux_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Common_module

  implicit none

  class(pm_hydrate_type) :: this
  SNES :: snes
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr

  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(hydrate_auxvar_type), pointer :: hyd_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(hydrate_auxvar_type) :: hyd_auxvar
  class(characteristic_curves_type), pointer :: characteristic_curves
  PetscInt :: local_id, ghosted_id, lid, gid, hid, iid, apid
  PetscInt :: offset
  PetscInt :: liq_pressure_index, gas_pressure_index, air_pressure_index, &
              air_frac_index, temp_index, liq_sat_index, gas_sat_index, &
              hyd_sat_index, ice_sat_index
  PetscReal :: s_extra
  PetscReal :: Pc_entry, dP

  PetscReal, parameter :: ALMOST_ZERO = 1.d-10
  PetscReal, parameter :: ALMOST_ONE = 1.d0-ALMOST_ZERO
  PetscReal, parameter :: eps_sat = 1.d-14
  PetscReal, parameter :: epsilon = 1.d-14
  PetscReal, parameter :: gravity = EARTH_GRAVITY

  PetscReal, pointer :: X_p(:),dX_p(:),dX_p2(:)

  field => this%realization%field
  grid => this%realization%patch%grid
  patch => this%realization%patch
  option => this%realization%option
  hyd_auxvars => this%realization%patch%aux%Hydrate%auxvars
  global_auxvars => this%realization%patch%aux%Global%auxvars

  lid = option%liquid_phase
  gid = option%gas_phase
  hid = option%hydrate_phase
  iid = option%ice_phase
  apid = option%air_pressure_id

  call VecCopy(dX,field%flow_dxx,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%flow_dxx,dX_p2,ierr);CHKERRQ(ierr)

  dX_p = -1.d0 * dX_p
  dX_p2 = -1.d0 * dX_p2

  changed = PETSC_TRUE

  if (this%check_post_convergence .and. hydrate_truncate_updates) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      offset = (local_id-1)*option%nflowdof
      hyd_auxvar = hyd_auxvars(ZERO_INTEGER,ghosted_id)
      characteristic_curves => patch%characteristic_curves_array( &
                               patch%cc_id(ghosted_id))%ptr

      ! Compute surface tension
      ! Compute entry pressure
      Pc_entry = 0.d0
      select type(sf => characteristic_curves%saturation_function)
        class is (sat_func_vg_type)
          Pc_entry = (1.d0 / characteristic_curves% &
                      saturation_function%GetAlpha_())
        !class is (sat_func_VG_STOMP_type)
        !  Pc_entry = characteristic_curves% &
        !             saturation_function%GetAlpha_() * &
        !             LIQUID_REFERENCE_DENSITY * gravity
        class default
      end select
      Pc_entry = 0.d0
      
      select case(global_auxvars(ghosted_id)%istate)
        case(L_STATE)
          liq_pressure_index = offset + ONE_INTEGER
          air_frac_index = offset + TWO_INTEGER
          temp_index = offset + THREE_INTEGER

          !Limit changes in liquid pressure
          dP = 5.d-2 * hyd_auxvar%pres(gid)
          dX_p(liq_pressure_index) = sign( min(dabs(dP), &
          dabs(dX_p(liq_pressure_index))),dX_p(liq_pressure_index))

          ! Zero negative corrections for zero aqueous CH4 or CO2
          if (X_p(air_frac_index) / epsilon < epsilon .and. &
              dX_p(air_frac_index) / epsilon < epsilon ) then
            dX_p(air_frac_index) = 0.d0
            dX_p2(air_frac_index) = 0.d0
          endif
          if ((X_p(air_frac_index) + dX_p(air_frac_index)) < 0.d0) &
             dX_p(air_frac_index) = - X_p(air_frac_index)

          ! Limit changes in temperature
          if (hyd_auxvar%sat(hid) > epsilon) then
            dX_p(temp_index) = sign(min(2.5d-1,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          else
            dX_p(temp_index) = sign(min(1.d0,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          endif

        case(G_STATE)
          gas_pressure_index = offset + ONE_INTEGER
          air_pressure_index = offset + TWO_INTEGER
          temp_index = offset + THREE_INTEGER

          !Limit changes in gas phase pressure
          dP = max(1.d6,2.5d-1*(X_p(gas_pressure_index)))
          dX_p(gas_pressure_index) = sign(min(dabs(dP), &
                                     dabs(dX_p(gas_pressure_index))), &
                                     dX_p(gas_pressure_index))

          !Limit changes in gas partial pressure
          if ((hyd_auxvar%pres(apid) / epsilon < epsilon) .and. &
              (dX_p(air_pressure_index)/epsilon < epsilon)) then
            dX_p(air_pressure_index) = 0.d0
          endif
          dP = max(1.d-1*hyd_auxvar%pres(apid),1.d4) 
          dX_p(air_pressure_index) = sign(min(dabs(dP), &
                                     dabs(dX_p(air_pressure_index))), &
                                     dX_p(air_pressure_index))
          if (hyd_auxvar%pres(apid) + dX_p(air_pressure_index) < 1.d-6) &
             dX_p(air_pressure_index) = - hyd_auxvar%pres(apid)

          ! Limit changes in temperature
          if (hyd_auxvar%sat(hid) > epsilon) then
            dX_p(temp_index) = sign(min(2.5d-1,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          else
            dX_p(temp_index) = sign(min(1.d0,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          endif

        case(H_STATE)
          gas_pressure_index = offset + ONE_INTEGER
          temp_index = offset + THREE_INTEGER

          !Limit changes in gas phase pressure
          dP = max(1.d6,2.5d-1*(X_p(gas_pressure_index)))
          dX_p(gas_pressure_index) = sign(min(dabs(dP), &
                                     dabs(dX_p(gas_pressure_index))), &
                                     dX_p(gas_pressure_index))
          ! Limit changes in temperature
          if (hyd_auxvar%sat(hid) > epsilon) then
            dX_p(temp_index) = sign(min(2.5d-1,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          else
            dX_p(temp_index) = sign(min(1.d0,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          endif

        case(I_STATE)
          gas_pressure_index = offset + ONE_INTEGER
          temp_index = offset + THREE_INTEGER

          !Limit changes in gas phase pressure
          dP = max(1.d6,2.5d-1*(X_p(gas_pressure_index)))
          dX_p(gas_pressure_index) = sign(min(dabs(dP), &
                                     dabs(dX_p(gas_pressure_index))), &
                                     dX_p(gas_pressure_index))

          ! Limit changes in temperature
          if (hyd_auxvar%sat(hid) > epsilon) then
            dX_p(temp_index) = sign(min(2.5d-1,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          else
            dX_p(temp_index) = sign(min(1.d0,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          endif

        case(GA_STATE)
          gas_pressure_index = offset + ONE_INTEGER
          gas_sat_index = offset + TWO_INTEGER
          temp_index = offset + THREE_INTEGER

          !Limit changes in gas phase pressure
          dP = max(1.d6,2.5d-1*(X_p(gas_pressure_index)-hyd_auxvar%pres(lid)))
          dX_p(gas_pressure_index) = sign(min(dabs(dP), &
                                     dabs(dX_p(gas_pressure_index))), &
                                     dX_p(gas_pressure_index))
          ! Relax pressure updates when transitioning to unsaturated
          ! conditions
          if((X_p(gas_pressure_index) + dX_p(gas_pressure_index)) - &
             (hyd_auxvar%pres(lid)) < Pc_entry ) then
            dX_p(gas_pressure_index) = 6.d-1*dX_p(gas_pressure_index)
          endif

          !Limit changes in gas saturation
          dP = 1.d-1
          dX_p(gas_sat_index) = sign(min(dabs(dP),dabs(dX_p(gas_sat_index))), &
                                dX_p(gas_sat_index))
          if (X_p(gas_sat_index) + dX_p(gas_sat_index) > 1.d0) &
             dX_p(gas_sat_index) = 1.d0 - X_p(gas_sat_index)
          if (X_p(gas_sat_index) + dX_p(gas_sat_index) < 0.d0) &
             dX_p(gas_sat_index) = - X_p(gas_sat_index)


          ! Limit changes in temperature
          if (hyd_auxvar%sat(hid) > epsilon) then
            dX_p(temp_index) = sign(min(2.5d-1,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          else
            dX_p(temp_index) = sign(min(1.d0,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          endif

        case(HG_STATE)
          gas_pressure_index = offset + ONE_INTEGER
          gas_sat_index = offset + TWO_INTEGER
          temp_index = offset + THREE_INTEGER

          !Limit changes in gas phase pressure
          dP = max(1.d6,2.5d-1*(X_p(gas_pressure_index)))
          dX_p(gas_pressure_index) = sign(min(dabs(dP), &
                                     dabs(dX_p(gas_pressure_index))), &
                                     dX_p(gas_pressure_index))

          !Limit changes in gas saturation
          dP = 1.d-1
          dX_p(gas_sat_index) = sign(min(dabs(dP),dabs(dX_p(gas_sat_index))), &
                                dX_p(gas_sat_index))
          if (X_p(gas_sat_index) + dX_p(gas_sat_index) > 1.d0) &
             dX_p(gas_sat_index) = 1.d0 - X_p(gas_sat_index)
          if (X_p(gas_sat_index) + dX_p(gas_sat_index) < 0.d0) &
             dX_p(gas_sat_index) = - X_p(gas_sat_index)


          ! Limit changes in temperature
          if (hyd_auxvar%sat(hid) > epsilon) then
            dX_p(temp_index) = sign(min(2.5d-1,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          else
            dX_p(temp_index) = sign(min(1.d0,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          endif

        case(HA_STATE)
          gas_pressure_index = offset + ONE_INTEGER
          hyd_sat_index = offset + TWO_INTEGER
          temp_index = offset + THREE_INTEGER

          !Limit changes in gas phase pressure
          dP = max(1.d6,2.5d-1*(X_p(gas_pressure_index)))
          dX_p(gas_pressure_index) = sign(min(dabs(dP), &
                                     dabs(dX_p(gas_pressure_index))), &
                                     dX_p(gas_pressure_index))
          ! Relax pressure updates when transitioning to unsaturated
          ! conditions
          if((X_p(gas_pressure_index) + dX_p(gas_pressure_index)) - &
             (hyd_auxvar%pres(lid)) < Pc_entry ) then
            dX_p(gas_pressure_index) = 6.d-1*dX_p(gas_pressure_index)
          endif

          !Limit changes in hydrate saturation
          dP = 1.d-1
          dX_p(hyd_sat_index) = sign(min(dabs(dP),dabs(dX_p(hyd_sat_index))), &
                                dX_p(hyd_sat_index))
          if(hyd_auxvar%sat(lid) < eps_sat) dX_p(hyd_sat_index) = &
                                            min(dX_p(hyd_sat_index),0.d0)
          if (X_p(hyd_sat_index) + dX_p(hyd_sat_index) > 1.d0) &
             dX_p(hyd_sat_index) = 1.d0 - X_p(hyd_sat_index)
          if (X_p(hyd_sat_index) + dX_p(hyd_sat_index) < 0.d0) &
             dX_p(hyd_sat_index) = - X_p(hyd_sat_index)

          ! Limit changes in temperature
          if (hyd_auxvar%sat(hid) > epsilon) then
            dX_p(temp_index) = sign(min(2.5d-1,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          else
            dX_p(temp_index) = sign(min(1.d0,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          endif

        case(HI_STATE)
          gas_pressure_index = offset + ONE_INTEGER
          hyd_sat_index = offset + TWO_INTEGER
          temp_index = offset + THREE_INTEGER

          !Limit changes in gas phase pressure
          dP = max(1.d6,2.5d-1*(X_p(gas_pressure_index)))
          dX_p(gas_pressure_index) = sign(min(dabs(dP), &
                                     dabs(dX_p(gas_pressure_index))), &
                                     dX_p(gas_pressure_index))

          !Limit changes in hydrate saturation
          dP = 1.d-1
          dX_p(hyd_sat_index) = sign(min(dabs(dP),dabs(dX_p(hyd_sat_index))), &
                                dX_p(hyd_sat_index))
          if(hyd_auxvar%sat(iid) < eps_sat) dX_p(hyd_sat_index) = &
                                            min(dX_p(hyd_sat_index),0.d0)
          if (X_p(hyd_sat_index) + dX_p(hyd_sat_index) > 1.d0) &
             dX_p(hyd_sat_index) = 1.d0 - X_p(hyd_sat_index)
          if (X_p(hyd_sat_index) + dX_p(hyd_sat_index) < 0.d0) &
             dX_p(hyd_sat_index) = - X_p(hyd_sat_index)

          ! Limit changes in temperature
          if (hyd_auxvar%sat(hid) > epsilon) then
            dX_p(temp_index) = sign(min(2.5d-1,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          else
            dX_p(temp_index) = sign(min(1.d0,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          endif

        case(GI_STATE)
          gas_pressure_index = offset + ONE_INTEGER
          ice_sat_index = offset + TWO_INTEGER
          temp_index = offset + THREE_INTEGER

          !Limit changes in gas phase pressure
          dP = max(1.d6,2.5d-1*(X_p(gas_pressure_index)))
          dX_p(gas_pressure_index) = sign(min(dabs(dP), &
                                     dabs(dX_p(gas_pressure_index))), &
                                     dX_p(gas_pressure_index))

          !Limit changes in ice saturation
          if(hyd_auxvar%sat(gid) < eps_sat) dX_p(ice_sat_index) = &
                                            min(dX_p(ice_sat_index),0.d0)
          if (X_p(ice_sat_index) + dX_p(ice_sat_index) > 1.d0) &
             dX_p(ice_sat_index) = 1.d0 - X_p(ice_sat_index)
          if (X_p(ice_sat_index) + dX_p(ice_sat_index) < 0.d0) &
             dX_p(ice_sat_index) = - X_p(ice_sat_index)

          ! Limit changes in temperature
          if (hyd_auxvar%sat(hid) > epsilon) then
            dX_p(temp_index) = sign(min(2.5d-1,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          else
            dX_p(temp_index) = sign(min(1.d0,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          endif

        case(AI_STATE)
          liq_pressure_index = offset + ONE_INTEGER
          air_frac_index = offset + TWO_INTEGER
          liq_sat_index = offset + THREE_INTEGER

          !Limit changes in liquid pressure
          dP = 5.d-2 * hyd_auxvar%pres(gid)
          dX_p(liq_pressure_index) = sign( min(dabs(dP), &
          dabs(dX_p(liq_pressure_index))),dX_p(liq_pressure_index))

          ! Zero negative corrections for zero aqueous CH4 or CO2
          if (X_p(air_frac_index) / epsilon < epsilon .and. &
              dX_p(air_frac_index) / epsilon < epsilon ) then
            dX_p(air_frac_index) = 0.d0
            dX_p2(air_frac_index) = 0.d0
          endif
          if ((X_p(air_frac_index) + dX_p(air_frac_index)) < 0.d0) &
             dX_p(air_frac_index) = - X_p(air_frac_index)

          !Limit changes in liquid saturation
          dP = 1.d-1
          dX_p(liq_sat_index) = sign(min(dabs(dP),dabs(dX_p(liq_sat_index))), &
                                dX_p(liq_sat_index))
          if (X_p(liq_sat_index) + dX_p(liq_sat_index) > 1.d0) &
             dX_p(liq_sat_index) = 1.d0 - X_p(liq_sat_index)
          if (X_p(liq_sat_index) + dX_p(liq_sat_index) < 0.d0) &
             dX_p(liq_sat_index) = - X_p(liq_sat_index)

        case(HGA_STATE)
          liq_sat_index = offset + ONE_INTEGER
          hyd_sat_index = offset + TWO_INTEGER
          temp_index = offset + THREE_INTEGER

          !Limit changes in liquid saturation
          dP = 1.d-1
          dX_p(liq_sat_index) = sign(min(dabs(dP),dabs(dX_p(liq_sat_index))), &
                                dX_p(liq_sat_index))
          if (X_p(liq_sat_index) + dX_p(liq_sat_index) > 1.d0) &
             dX_p(liq_sat_index) = 1.d0 - X_p(liq_sat_index)
          if (X_p(liq_sat_index) + dX_p(liq_sat_index) < 0.d0) &
             dX_p(liq_sat_index) = - X_p(liq_sat_index)

          !Limit changes in hydrate saturation
          dX_p(hyd_sat_index) = sign(min(dabs(dP),dabs(dX_p(hyd_sat_index))), &
             dX_p(hyd_sat_index))

          if (X_p(hyd_sat_index) + dX_p(hyd_sat_index) > 1.d0) &
             dX_p(hyd_sat_index) = 1.d0 - X_p(hyd_sat_index)
          if (X_p(hyd_sat_index) + dX_p(hyd_sat_index) < 0.d0) &
             dX_p(hyd_sat_index) = - X_p(hyd_sat_index)

          if ((X_p(hyd_sat_index) + dX_p(hyd_sat_index) + &
              X_p(liq_sat_index) + dX_p(liq_sat_index)) > 1.d0) then
            s_extra = 1.d0 - ((X_p(hyd_sat_index) + dX_p(hyd_sat_index) + &
                      X_p(liq_sat_index) + dX_p(liq_sat_index)))
            dX_p(hyd_sat_index) = dX_p(hyd_sat_index) + s_extra / 2.d0
            dX_p(liq_sat_index) = dX_p(liq_sat_index) + s_extra / 2.d0
          endif

          ! Limit changes in temperature
          if (hyd_auxvar%sat(hid) > epsilon) then
            dX_p(temp_index) = sign(min(2.5d-1,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          else
            dX_p(temp_index) = sign(min(1.d0,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          endif

        case(HAI_STATE)
          gas_pressure_index = offset + ONE_INTEGER
          liq_sat_index = offset + TWO_INTEGER
          ice_sat_index = offset + THREE_INTEGER

          !Limit changes in gas phase pressure
          dP = max(1.d6,2.5d-1*(X_p(gas_pressure_index)))
          dX_p(gas_pressure_index) = sign(min(dabs(dP), &
                                     dabs(dX_p(gas_pressure_index))), &
                                     dX_p(gas_pressure_index))

          !Limit changes in liquid saturation
          dP = 1.d-1
          dX_p(liq_sat_index) = sign(min(dabs(dP),dabs(dX_p(liq_sat_index))), &
                                dX_p(liq_sat_index))
          if (X_p(liq_sat_index) + dX_p(liq_sat_index) > 1.d0) &
             dX_p(liq_sat_index) = 1.d0 - X_p(liq_sat_index)
          if (X_p(liq_sat_index) + dX_p(liq_sat_index) < 0.d0) &
             dX_p(liq_sat_index) = - X_p(liq_sat_index)


          !Limit changes in ice saturation
          dP = 1.d-1
          dX_p(ice_sat_index) = sign(min(dabs(dP),dabs(dX_p(ice_sat_index))), &
                                dX_p(ice_sat_index))
          if (X_p(ice_sat_index) + dX_p(ice_sat_index) > 1.d0) &
             dX_p(ice_sat_index) = 1.d0 - X_p(ice_sat_index)
          if (X_p(ice_sat_index) + dX_p(ice_sat_index) < 0.d0) &
             dX_p(ice_sat_index) = - X_p(ice_sat_index)

        case(HGI_STATE)
          ice_sat_index = offset + ONE_INTEGER
          hyd_sat_index = offset + TWO_INTEGER
          temp_index = offset + THREE_INTEGER

          !Limit changes in ice saturation
          dP = 1.d-1
          dX_p(ice_sat_index) = sign(min(dabs(dP),dabs(dX_p(ice_sat_index))), &
                                dX_p(ice_sat_index))
          if (X_p(ice_sat_index) + dX_p(ice_sat_index) > 1.d0) &
             dX_p(ice_sat_index) = 1.d0 - X_p(ice_sat_index)
          if (X_p(ice_sat_index) + dX_p(ice_sat_index) < 0.d0) &
             dX_p(ice_sat_index) = - X_p(ice_sat_index)

          !Limit changes in hydrate saturation
          dP = 1.d-1
          dX_p(hyd_sat_index) = sign(min(dabs(dP),dabs(dX_p(hyd_sat_index))), &
                                dX_p(hyd_sat_index))
          if (X_p(hyd_sat_index) + dX_p(hyd_sat_index) > 1.d0) &
             dX_p(hyd_sat_index) = 1.d0 - X_p(hyd_sat_index)
          if (X_p(hyd_sat_index) + dX_p(hyd_sat_index) < 0.d0) &
             dX_p(hyd_sat_index) = - X_p(hyd_sat_index)

          ! Limit changes in temperature
          if (hyd_auxvar%sat(hid) > epsilon) then
            dX_p(temp_index) = sign(min(2.5d-1,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          else
            dX_p(temp_index) = sign(min(1.d0,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          endif

        case(GAI_STATE)
          gas_pressure_index = offset + ONE_INTEGER
          gas_sat_index = offset + TWO_INTEGER
          temp_index = offset + THREE_INTEGER

          !Limit changes in gas phase pressure
          dP = max(1.d6,2.5d-1*(X_p(gas_pressure_index)))
          dX_p(gas_pressure_index) = sign(min(dabs(dP), &
                                     dabs(dX_p(gas_pressure_index))), &
                                     dX_p(gas_pressure_index))

          !Limit changes in gas saturation
          dP = 1.d-1
          dX_p(gas_sat_index) = sign(min(dabs(dP),dabs(dX_p(gas_sat_index))), &
                                dX_p(gas_sat_index))
          if (X_p(gas_sat_index) + dX_p(gas_sat_index) > 1.d0) &
             dX_p(gas_sat_index) = 1.d0 - X_p(gas_sat_index)
          if (X_p(gas_sat_index) + dX_p(gas_sat_index) < 0.d0) &
             dX_p(gas_sat_index) = - X_p(gas_sat_index)


          ! Limit changes in temperature
          if (hyd_auxvar%sat(hid) > epsilon) then
            dX_p(temp_index) = sign(min(2.5d-1,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          else
            dX_p(temp_index) = sign(min(1.d0,dabs(dX_p(temp_index))), &
                               dX_p(temp_index))
          endif

        case(HGAI_STATE)
          liq_sat_index = offset + ONE_INTEGER
          gas_sat_index = offset + TWO_INTEGER
          ice_sat_index = offset + THREE_INTEGER

          !Limit changes in liquid saturation
          dP = 1.d-1
          dX_p(liq_sat_index) = sign(min(dabs(dP),dabs(dX_p(liq_sat_index))), &
                                dX_p(liq_sat_index))
          if (X_p(liq_sat_index) + dX_p(liq_sat_index) > 1.d0) &
             dX_p(liq_sat_index) = 1.d0 - X_p(liq_sat_index)
          if (X_p(liq_sat_index) + dX_p(liq_sat_index) < 0.d0) &
             dX_p(liq_sat_index) = - X_p(liq_sat_index)

          !Limit changes in gas saturation
          dP = 1.d-1
          dX_p(gas_sat_index) = sign(min(dabs(dP),dabs(dX_p(gas_sat_index))), &
                                dX_p(gas_sat_index))
          if (X_p(gas_sat_index) + dX_p(gas_sat_index) > 1.d0) &
             dX_p(gas_sat_index) = 1.d0 - X_p(gas_sat_index)
          if (X_p(gas_sat_index) + dX_p(gas_sat_index) < 0.d0) &
             dX_p(gas_sat_index) = - X_p(gas_sat_index)

          !Limit changes in ice saturation
          dP = 1.d-1
          dX_p(ice_sat_index) = sign(min(dabs(dP),dabs(dX_p(ice_sat_index))), &
                                dX_p(ice_sat_index))
          if (X_p(ice_sat_index) + dX_p(ice_sat_index) > 1.d0) &
             dX_p(ice_sat_index) = 1.d0 - X_p(ice_sat_index)
          if (X_p(ice_sat_index) + dX_p(ice_sat_index) < 0.d0) &
             dX_p(ice_sat_index) = - X_p(ice_sat_index)

      end select
    enddo
  endif

  if (this%damping_factor > 0.d0) then
    dX_p = dX_p*this%damping_factor
    changed = PETSC_TRUE
  endif

  dX_p = -1.d0 * dX_p

  call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%flow_dxx,dX_p2,ierr);CHKERRQ(ierr)

end subroutine PMHydrateCheckUpdatePre

! ************************************************************************** !

subroutine PMHydrateCheckUpdatePost(this,snes,X0,dX,X1,dX_changed, &
                                    X1_changed,ierr)
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !
  use Hydrate_Aux_module
  use Global_Aux_module
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module

  implicit none

  class(pm_hydrate_type) :: this
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
  PetscBool :: converged_abs_update_flag(3,15)
  PetscBool :: converged_rel_update_flag(3,15)
  PetscInt :: converged_abs_update_cell(3,15)
  PetscInt :: converged_rel_update_cell(3,15)
  PetscReal :: converged_abs_update_real(3,15)
  PetscReal :: converged_rel_update_real(3,15)
  PetscInt :: istate
  PetscBool :: converged_absolute
  PetscBool :: converged_relative

  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch
  global_auxvars => patch%aux%Global%auxvars

  hydrate_allow_state_change = PETSC_TRUE

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

  if (hydrate_check_updates) then
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
        if (dX_X0 > this%rel_update_inf_tol(idof,istate)) then
          converged_relative = PETSC_FALSE
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
  endif
  call VecRestoreArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)

  this%converged_flag(:,:,ABS_UPDATE_INDEX) = converged_abs_update_flag(:,:)
  this%converged_flag(:,:,REL_UPDATE_INDEX) = converged_rel_update_flag(:,:)
  this%converged_real(:,:,ABS_UPDATE_INDEX) = converged_abs_update_real(:,:)
  this%converged_real(:,:,REL_UPDATE_INDEX) = converged_rel_update_real(:,:)
  this%converged_cell(:,:,ABS_UPDATE_INDEX) = converged_abs_update_cell(:,:)
  this%converged_cell(:,:,REL_UPDATE_INDEX) = converged_rel_update_cell(:,:)

end subroutine PMHydrateCheckUpdatePost

! ************************************************************************** !

subroutine PMHydrateCheckConvergence(this,snes,it,xnorm,unorm,fnorm, &
                                     reason,ierr)
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !
  use Convergence_module
  use Hydrate_Aux_module
  use Global_Aux_module
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  use String_module
  use EOS_Gas_module
  use Material_Aux_module

  implicit none

  class(pm_hydrate_type) :: this
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
  type(hydrate_auxvar_type), pointer :: hyd_auxvars(:,:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(hydrate_auxvar_type) :: hyd_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum2_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscInt :: local_id, ghosted_id, natural_id
  PetscInt :: offset, ival, idof, itol
  PetscInt :: gid, lid, acid, wid, eid, hid, iid, spid, apid, vpid
  PetscReal :: Psat, Pv, Prvap, Pa
  PetscReal :: xag, xwg, xal, xsl, xwl, xmolag, xmolwg, xmolal, &
               xmolsl, xmolwl, x_salt_dissolved
  PetscReal :: R, A, R_A
  PetscReal :: res_scaled, residual, accumulation, update
  PetscReal :: Hc
  PetscReal, parameter :: A_zero = 1.d-15
  PetscBool :: converged_abs_residual_flag(3,15)
  PetscReal :: converged_abs_residual_real(3,15)
  PetscInt :: converged_abs_residual_cell(3,15)
  PetscBool :: converged_scaled_residual_flag(3,15)
  PetscReal :: converged_scaled_residual_real(3,15)
  PetscInt :: converged_scaled_residual_cell(3,15)
  PetscInt :: istate
  PetscBool :: converged_absolute
  PetscBool :: converged_scaled
  PetscMPIInt :: mpi_int
  PetscBool :: flags(181)
  character(len=MAXSTRINGLENGTH) :: string

  PetscBool :: rho_flag
  PetscReal, parameter :: T_ref = 273.15d0
  PetscReal, parameter :: epsilon = 1.d-14

  character(len=14), parameter :: state_string(15) = &
    ['Liquid State  ','Gas State     ','Hydrate State ','Ice State     ', &
     'GA State      ','HG State      ','HA State      ','HI State      ', &
     'GI State      ','AI State      ','HGA State     ','HAI State     ', &
     'HGI State     ','GAI State     ','Quad State    ']
  character(len=17), parameter :: dof_string(3,15) = &
    reshape(['Liquid Pressure  ','Air Mole Fraction','Temperature      ', &
             'Gas Pressure     ','Air Pressure     ','Temperature      ', &
             'Gas Pressure     ','Air Mole Frac Hyd','Temperature      ', &
             'Gas Pressure     ','Air Mole Frac Ice','Temperature      ', &
             'Gas Pressure     ','Gas Saturation   ','Temperature      ', &
             'Gas Pressure     ','Gas Saturation   ','Temperature      ', &
             'Gas Pressure     ','Hydrate Sat      ','Temperature      ', &
             'Gas Pressure     ','Hydrate Sat      ','Temperature      ', &
             'Gas Pressure     ','Ice Saturation   ','Temperature      ', &
             'Gas Pressure     ','Air Mole Fraction','Liquid Saturation', &
             'Liquid Saturation','Hydrate Sat      ','Temperature      ', &
             'Gas Pressure     ','Liquid Saturation','Ice Saturation   ', &
             'Hydrate Sat      ','Ice Saturation   ','Temperature      ', &
             'Gas Pressure     ','Liquid Saturation','Ice Saturation   ', &
             'Gas Saturation   ','Liquid Saturation','Ice Saturation   '], &
             shape(dof_string))
  character(len=15), parameter :: tol_string(MAX_INDEX) = &
    ['Absolute Update','Relative Update','Residual       ','Scaled Residual']

  patch => this%realization%patch
  option => this%realization%option
  field => this%realization%field
  grid => patch%grid
  global_auxvars => patch%aux%Global%auxvars
  hyd_auxvars => this%realization%patch%aux%Hydrate%auxvars
  material_auxvars => patch%aux%Material%auxvars

  lid = option%liquid_phase
  gid = option%gas_phase
  hid = option%hydrate_phase
  iid = option%ice_phase
  acid = option%air_id
  apid = option%air_pressure_id
  spid = option%saturation_pressure_id
  vpid = option%vapor_pressure_id

  call SNESNewtonTRDCGetRhoFlag(snes,rho_flag,ierr);CHKERRQ(ierr);

  if (this%option%flow%using_newtontrdc) then
    if (hydrate_newtontrdc_prev_iter_num == it) then
      hydrate_sub_newton_iter_num = hydrate_sub_newton_iter_num + 1
    endif
    hydrate_newtontrdc_prev_iter_num = it
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

    if (hydrate_full_convergence) then
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

    else
      do local_id = 1, grid%nlmax
        offset = (local_id-1)*option%nflowdof
        ghosted_id = grid%nL2G(local_id)
        natural_id = grid%nG2A(ghosted_id)

        hyd_auxvar = hyd_auxvars(ZERO_INTEGER,ghosted_id)
        material_auxvar = material_auxvars(ghosted_id)

        if (patch%imat(ghosted_id) <= 0) cycle
        istate = global_auxvars(ghosted_id)%istate
        do idof = 1, option%nflowdof
          res_scaled = 0.d0
          ival = offset+idof
          converged_absolute = PETSC_FALSE
          converged_scaled = PETSC_TRUE

          residual = r_p(ival)
          accumulation = accum2_p(ival)
          update = dX_p(ival)

          select case (istate)
            case(L_STATE)
              if (idof == ONE_INTEGER) then
              ! Water equation
                res_scaled = min(dabs(update)/ (dabs(hyd_auxvar%pres(lid))), &
                             dabs(residual / (accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == TWO_INTEGER) then
              ! Air equation
                call HydrateComputeSaltSolubility(hyd_auxvar%temp,xsl)
                xsl = min(hyd_auxvar%m_salt(1),xsl)
                call HydrateBrineSaturationPressure(hyd_auxvar%temp,xsl,Psat)
                Pv = hyd_auxvar%pres(vpid)
                Prvap = Psat
                x_salt_dissolved = xsl
                call HydrateEquilibrate(hyd_auxvar%temp,hyd_auxvar%pres(lid), &
                                     istate, &
                                     hyd_auxvar%sat(hid), &
                                     Pa, Pv, Psat, Prvap, &
                                     xag, xwg, xal, xsl, xwl, &
                                     xmolag, xmolwg, xmolal, xmolsl, &
                                     xmolwl, &
                                     patch%characteristic_curves_array( &
                                     patch%cc_id(ghosted_id))%ptr, &
                                     material_auxvar, option)
                x_salt_dissolved = x_salt_dissolved + &
                                   (xsl - x_salt_dissolved) * &
                                   (hyd_auxvar%xmass(acid,lid) / xal)
                if (hyd_auxvar%xmol(acid,lid) > (1.d-6 * xmolal)) then
                  Hc = HydrateHenryCO2(hyd_auxvar%temp,x_salt_dissolved)
                  res_scaled = min(dabs(update) / &
                               max(hyd_auxvar%pres(gid)/Hc, &
                               HYDRATE_REFERENCE_PRESSURE/Hc), &
                               dabs(residual / (accumulation + epsilon)))
                  ! find max value regardless of convergence
                  if (converged_scaled_residual_real(idof,istate) < &
                      res_scaled) then
                   converged_scaled_residual_real(idof,istate) = res_scaled
                   converged_scaled_residual_cell(idof,istate) = natural_id
                  endif
                endif
              elseif (idof == THREE_INTEGER) then
              ! Energy equation
                res_scaled = dabs(update)/T_ref
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              endif
            case(G_STATE)
              if (idof == ONE_INTEGER) then
              ! Water equation
                call HydrateComputeSaltSolubility(hyd_auxvar%temp,xsl)
                if (hyd_auxvar%m_salt(2) <= 0.d0) xsl = 0.d0
                call HydrateBrineSaturationPressure(hyd_auxvar%temp,xsl,Psat)

                res_scaled = min(dabs(update) / Psat, &
                             dabs(residual / (accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == TWO_INTEGER) then
              ! Air equation
                res_scaled = min(dabs(update)/dabs(hyd_auxvar%pres(apid)), &
                             dabs(residual / (accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == THREE_INTEGER) then
              ! Energy equation
                res_scaled = dabs(update)/T_ref
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              endif
            case(H_STATE)
              if (idof == ONE_INTEGER) then
              ! Water equation
                res_scaled = min(dabs(update)/dabs(hyd_auxvar%pres(gid)), &
                             dabs(residual / (accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == TWO_INTEGER) then
              ! Air equation
              elseif (idof == THREE_INTEGER) then
              ! Energy equation
                res_scaled = dabs(update)/T_ref
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              endif
            case(I_STATE)
              if (idof == ONE_INTEGER) then
              ! Water equation
                res_scaled = min(dabs(update)/dabs(hyd_auxvar%pres(gid)), &
                             dabs(residual / (accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == TWO_INTEGER) then
              ! Air equation
              elseif (idof == THREE_INTEGER) then
              ! Energy equation
                res_scaled = dabs(update)/T_ref
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              endif
            case(GA_STATE)
              if (idof == ONE_INTEGER) then
              ! Water equation
                res_scaled = min(dabs(update)/dabs(hyd_auxvar%pres(lid)), &
                             dabs(residual / (accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == TWO_INTEGER) then
              ! Air equation
                if (hyd_auxvar%sat(gid) > 1.d-3) then
                  res_scaled = min(dabs(update)/dabs(hyd_auxvar%pres(gid)), &
                               dabs(residual/(accumulation + epsilon)))
                  ! find max value regardless of convergence
                  if (converged_scaled_residual_real(idof,istate) < &
                      res_scaled) then
                   converged_scaled_residual_real(idof,istate) = res_scaled
                   converged_scaled_residual_cell(idof,istate) = natural_id
                  endif
                endif
              elseif (idof == THREE_INTEGER) then
              ! Energy equation
                res_scaled = dabs(update)/T_ref

                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              endif
            case(HG_STATE)
              if (idof == ONE_INTEGER) then
              ! Water equation
                res_scaled = min(dabs(update)/dabs(hyd_auxvar%pres(gid)), &
                             dabs(residual / (accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == TWO_INTEGER) then
              ! Air equation
                res_scaled = min(dabs(update), &
                             dabs(residual/(accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == THREE_INTEGER) then
              ! Energy equation
                res_scaled = dabs(update)/T_ref

                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              endif
            case(HA_STATE)
              if (idof == ONE_INTEGER) then
              ! Water equation
                res_scaled = min(dabs(update)/dabs(hyd_auxvar%pres(gid)), &
                             dabs(residual / (accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == TWO_INTEGER) then
              ! Air equation
                res_scaled = min(dabs(update), &
                             dabs(residual/(accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == THREE_INTEGER) then
              ! Energy equation
                res_scaled = dabs(update)/T_ref

                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              endif
            case(HI_STATE)
              if (idof == ONE_INTEGER) then
              ! Water equation
                res_scaled = min(dabs(update)/dabs(hyd_auxvar%pres(gid)), &
                             dabs(residual / (accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == TWO_INTEGER) then
              ! Air equation
                res_scaled = min(dabs(update), &
                             dabs(residual/(accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == THREE_INTEGER) then
              ! Energy equation
                res_scaled = dabs(update)/T_ref

                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              endif
            case(GI_STATE)
              if (idof == ONE_INTEGER) then
              ! Water equation
                res_scaled = min(dabs(update)/dabs(hyd_auxvar%pres(gid)), &
                             dabs(residual / (accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == TWO_INTEGER) then
              ! Air equation
                res_scaled = min(dabs(update), &
                             dabs(residual/(accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == THREE_INTEGER) then
              ! Energy equation
                res_scaled = dabs(update)/T_ref

                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              endif
            case(AI_STATE)
              if (idof == ONE_INTEGER) then
              ! Water equation
                res_scaled = min(dabs(update)/dabs(hyd_auxvar%pres(lid)), &
                             dabs(residual / (accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == TWO_INTEGER) then
              ! Air equation
                if (hyd_auxvar%xmass(acid,lid) > 1.d-10) then
                  call EOSGasHenry(hyd_auxvar%temp,hyd_auxvar%pres(spid), &
                                   Hc,ierr)
                  res_scaled = min(dabs(update) / &
                               max(hyd_auxvar%pres(gid)/Hc, &
                               HYD_REFERENCE_PRESSURE/Hc), &
                               dabs(residual / (accumulation + epsilon)))
                  ! find max value regardless of convergence
                  if (converged_scaled_residual_real(idof,istate) < &
                      res_scaled) then
                   converged_scaled_residual_real(idof,istate) = res_scaled
                   converged_scaled_residual_cell(idof,istate) = natural_id
                  endif
                endif
              elseif (idof == THREE_INTEGER) then
              ! Energy equation
                res_scaled = residual/(dabs(accumulation + epsilon))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              endif
            case(HGA_STATE)
              if (idof == ONE_INTEGER) then
              ! Water equation
                res_scaled = min(dabs(update), &
                             dabs(residual/(accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == TWO_INTEGER) then
              ! Air equation
                res_scaled = min(dabs(update), &
                             dabs(residual/(accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == THREE_INTEGER) then
              ! Energy equation
                res_scaled = dabs(update)/T_ref

                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              endif
            case(HAI_STATE)
              if (idof == ONE_INTEGER) then
              ! Water equation
                res_scaled = min(dabs(update)/dabs(hyd_auxvar%pres(gid)), &
                             dabs(residual / (accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == TWO_INTEGER) then
              ! Air equation
                res_scaled = min(dabs(update), &
                             dabs(residual/(accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == THREE_INTEGER) then
              ! Energy equation
                res_scaled = residual/(dabs(accumulation + epsilon))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              endif
            case(HGI_STATE)
              if (idof == ONE_INTEGER) then
              ! Water equation
                res_scaled = min(dabs(update), &
                             dabs(residual/(accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == TWO_INTEGER) then
              ! Air equation
                res_scaled = min(dabs(update), &
                             dabs(residual/(accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == THREE_INTEGER) then
              ! Energy equation
                res_scaled = dabs(update)/T_ref

                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              endif
            case(GAI_STATE)
              if (idof == ONE_INTEGER) then
              ! Water equation
                res_scaled = min(dabs(update)/dabs(hyd_auxvar%pres(gid)), &
                             dabs(residual / (accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == TWO_INTEGER) then
              ! Air equation
                res_scaled = min(dabs(update), &
                             dabs(residual/(accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == THREE_INTEGER) then
              ! Energy equation
                res_scaled = dabs(update)/T_ref

                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              endif
            case(HGAI_STATE)
              if (idof == ONE_INTEGER) then
              ! Water equation
                res_scaled = min(dabs(update), &
                             dabs(residual/(accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == TWO_INTEGER) then
              ! Air equation
                res_scaled = min(dabs(update), &
                             dabs(residual/(accumulation + epsilon)))
                ! find max value regardless of convergence
                if (converged_scaled_residual_real(idof,istate) < &
                    res_scaled) then
                 converged_scaled_residual_real(idof,istate) = res_scaled
                 converged_scaled_residual_cell(idof,istate) = natural_id
                endif
              elseif (idof == THREE_INTEGER) then
              ! Energy equation
                res_scaled = residual/(dabs(accumulation + epsilon))
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
    !mpi_int = 9*MAX_INDEX
    flags(1:45*MAX_INDEX) = reshape(this%converged_flag,(/45*MAX_INDEX/))
    flags(181) = .not.hydrate_high_temp_ts_cut
    mpi_int = 181
    ! do not perform an all reduce on cell id as this info is not printed
    ! in parallel
    call MPI_Allreduce(MPI_IN_PLACE,flags,mpi_int,MPI_LOGICAL,MPI_LAND, &
                       option%mycomm,ierr);CHKERRQ(ierr)
    this%converged_flag = reshape(flags(1:45*MAX_INDEX),(/3,15,MAX_INDEX/))
    hydrate_high_temp_ts_cut = .not.flags(181)

    mpi_int = 45*MAX_INDEX
    call MPI_Allreduce(MPI_IN_PLACE,this%converged_real,mpi_int, &
                       MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm, &
                       ierr);CHKERRQ(ierr)

    option%convergence = CONVERGENCE_CONVERGED

    do itol = 1, MAX_INDEX
      do istate = 1, 15
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
                else
                  string = '   ' // trim(tol_string(itol)) // ', ' // &
                   trim(state_string(istate)) // ', Energy'
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
        hydrate_state_changed .and. &
        .not.rho_flag) then
      if (hydrate_newtontrdc_hold_inner) then
        ! if we hold inner iterations, we must not change state in
        ! the inner iteration. If we reach convergence in an inner
        ! newtontrdc iteration, then we must force an outer iteration
        ! to allow state change in case the solutions are
        ! out-of-bounds of the states -hdp
        hydrate_force_iteration = PETSC_TRUE
        hydrate_state_changed = PETSC_FALSE
      else
        ! if we have state changes, we exit out of inner iteration
        ! and go to the next newton iteration. the tr inner iteration
        !  should only be used when there is no state changes
        ! if rho is satisfied in inner iteration, the algorithm already
        ! exited the inner iteration. -heeho
        hydrate_force_iteration = PETSC_TRUE
        hydrate_state_changed = PETSC_FALSE
      endif
    endif

    call MPI_Allreduce(MPI_IN_PLACE,hydrate_force_iteration,ONE_INTEGER, &
                       MPI_LOGICAL,MPI_LOR,option%mycomm,ierr)
    if (hydrate_force_iteration) then
      if (.not.hydrate_newtontrdc_hold_inner) then
        option%convergence = CONVERGENCE_BREAKOUT_INNER_ITER
        hydrate_force_iteration = PETSC_FALSE
      elseif (hydrate_newtontrdc_hold_inner .and. &
               option%convergence == CONVERGENCE_CONVERGED) then
        option%convergence = CONVERGENCE_BREAKOUT_INNER_ITER
        hydrate_force_iteration = PETSC_FALSE
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
        string = '    Exceeded Hydrate Mode Max Newton Iterations'
        call PrintMsg(option,string)
      endif
    endif

    if (hydrate_high_temp_ts_cut) then
      hydrate_high_temp_ts_cut = PETSC_FALSE
      string = '    Exceeded Hydrate Mode EOS max temperature'
      call PrintMsg(option,string)
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

end subroutine PMHydrateCheckConvergence

! ************************************************************************** !

subroutine PMHydrateTimeCut(this)
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !

  use Hydrate_module, only : HydrateTimeCut

  implicit none

  class(pm_hydrate_type) :: this

  call PMSubsurfaceFlowTimeCut(this)
  call HydrateTimeCut(this%realization)

end subroutine PMHydrateTimeCut

! ************************************************************************** !

subroutine PMHydrateUpdateSolution(this)
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !

  use Hydrate_module, only : HydrateUpdateSolution, &
                             HydrateMapBCAuxVarsToGlobal

  implicit none

  class(pm_hydrate_type) :: this

  call PMSubsurfaceFlowUpdateSolution(this)
  call HydrateUpdateSolution(this%realization)
  call HydrateMapBCAuxVarsToGlobal(this%realization)

end subroutine PMHydrateUpdateSolution

! ************************************************************************** !

subroutine PMHydrateUpdateAuxVars(this)
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  use Hydrate_module, only : HydrateUpdateAuxVars

  implicit none

  class(pm_hydrate_type) :: this

  call HydrateUpdateAuxVars(this%realization,PETSC_FALSE)

end subroutine PMHydrateUpdateAuxVars

! ************************************************************************** !

subroutine PMHydrateMaxChange(this)
  !
  ! Not needed given HydrateMaxChange is called in PostSolve
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !

  use Realization_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use Global_Aux_module
  use Hydrate_Aux_module

  implicit none

  class(pm_hydrate_type) :: this

  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscReal :: max_change_local(9)
  PetscReal :: max_change_global(9)
  PetscReal :: max_change
  PetscInt :: i, j
  PetscInt :: ghosted_id


  PetscErrorCode :: ierr

  realization => this%realization
  option => realization%option
  field => realization%field
  grid => realization%patch%grid

  global_auxvars => realization%patch%aux%global%auxvars

  max_change_global = 0.d0
  max_change_local = 0.d0

  !max change variables:[LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
  !                      LIQUID_MOLE_FRACTION, TEMPERATURE, GAS_SATURATION,
  !                      HYDRATE_SATURATION, LIQUID_SATUARTION, ICE_SATUARTION]
  do i = 1, 9
    call RealizationGetVariable(realization,field%work, &
                                this%max_change_ivar(i), &
                                this%max_change_isubvar(i))
    ! yes, we could use VecWAXPY and a norm here, but we need the ability
    ! to customize
    call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%max_change_vecs(i),vec_ptr2,ierr);CHKERRQ(ierr)
    max_change = 0.d0
    if (i==1 .and. hyd_chk_max_dpl_liq_state_only) then
      do j = 1,grid%nlmax
        ghosted_id = grid%nL2G(j)
        if (global_auxvars(ghosted_id)%istate /= L_STATE) then
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
  call MPI_Allreduce(max_change_local,max_change_global,NINE_INTEGER, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm, &
                     ierr);CHKERRQ(ierr)
  ! print them out
  write(option%io_buffer,'("  --> max change: dpl= ",1pe12.4, " dpg= ",&
                         &1pe12.4," dpa= ",1pe12.4)') max_change_global(1:3)
  call PrintMsg(option)
  write(option%io_buffer,'(17x," dxa= ",1pe12.4,"  dt= ",1pe12.4,&
                         &" dsg= ",1pe12.4)') max_change_global(4:6)
  call PrintMsg(option)
  write(option%io_buffer,'(17x," dsh= ",1pe12.4," dsl= ",1pe12.4,&
                         &" dsi= ",1pe12.4)') max_change_global(7:9)
  call PrintMsg(option)

  ! max change variables:[LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
  !                       LIQUID_MOLE_FRACTION, TEMPERATURE, GAS_SATURATION,
  !                       HYDRATE_SATURATION, LIQUID_SATURATION, ICE_SATURATION]
  ! ignore air pressure as it jumps during phase change
  this%max_pressure_change = maxval(max_change_global(1:2))
  this%max_xmol_change = max_change_global(4)
  this%max_temperature_change = max_change_global(5)
  this%max_saturation_change = maxval(max_change_global(6:9))

end subroutine PMHydrateMaxChange

! ************************************************************************** !

subroutine PMHydrateComputeMassBalance(this,mass_balance_array)
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !

  use Hydrate_module, only : HydrateComputeMassBalance

  implicit none

  class(pm_hydrate_type) :: this
  PetscReal :: mass_balance_array(:)

  call HydrateComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMHydrateComputeMassBalance

! ************************************************************************** !

subroutine PMHydrateInputRecord(this)
  !
  ! Writes ingested information to the input record file.
  !
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  !

  implicit none

  class(pm_hydrate_type) :: this

  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'hydrate'
  if (this%check_post_convergence) then
    write(id,'(a29)',advance='no') 'ITOL_SCALED_RESIDUAL: '
    write(id,'(a)') 'ON'
    write(id,'(a29)',advance='no') 'ITOL_RELATIVE_RESIDUAL: '
    write(id,'(a)') 'ON'
  endif

end subroutine PMHydrateInputRecord

! ************************************************************************** !

subroutine PMHydrateCheckpointBinary(this,viewer)
  !
  ! Checkpoints data associated with Hydrate PM
  !
  ! Author: Michael Nole
  ! Date: 07/23/19

  use Checkpoint_module
  use Global_module

  implicit none
#include "petsc/finclude/petscviewer.h"

  class(pm_hydrate_type) :: this
  PetscViewer :: viewer

  call PMSubsurfaceFlowCheckpointBinary(this,viewer)

end subroutine PMHydrateCheckpointBinary

! ************************************************************************** !

subroutine PMHydrateRestartBinary(this,viewer)
  !
  ! Restarts data associated with Hydrate PM
  !
  ! Author: Michael Nole
  ! Date: 07/23/19

  use Checkpoint_module
  use Global_module

  implicit none
#include "petsc/finclude/petscviewer.h"

  class(pm_hydrate_type) :: this
  PetscViewer :: viewer

  call PMSubsurfaceFlowRestartBinary(this,viewer)

end subroutine PMHydrateRestartBinary
! ************************************************************************** !

subroutine PMHydrateDestroy(this)
  !
  ! Destroys Hydrate process model
  !
  ! Author: Michael Nole
  ! Date: 07/23/19
  !

  use Hydrate_module, only : HydrateDestroy

  implicit none

  class(pm_hydrate_type) :: this

  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  deallocate(this%max_change_ivar)
  nullify(this%max_change_ivar)
  deallocate(this%max_change_isubvar)
  nullify(this%max_change_isubvar)

  nullify(this%hydrate_parameters)

  ! preserve this ordering
  call HydrateDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)

end subroutine PMHydrateDestroy

end module PM_Hydrate_class
