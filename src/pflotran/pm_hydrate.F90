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
    PetscInt :: hydrate_newton_max_iter
  contains
    procedure, public :: ReadSimulationBlock => PMHydrateRead
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
                               GAS_SATURATION
  use Upwind_Direction_module

  implicit none
  
  class(pm_hydrate_type), pointer :: PMHydrateCreate

  class(pm_hydrate_type), pointer :: this
  
  PetscReal, parameter :: ref_temp = 20.d0 !degrees C
  PetscReal, parameter :: ref_pres = 101325.d0 !Pa
  PetscReal, parameter :: ref_sat = 0.5
  PetscReal, parameter :: ref_xmol = 1.d-6
                                             
  !MAN optimized:
  PetscReal, parameter :: pres_abs_inf_tol = 1.d0 ! Reference tolerance [Pa]
  PetscReal, parameter :: temp_abs_inf_tol = 1.d-5
  PetscReal, parameter :: sat_abs_inf_tol = 1.d-5
  PetscReal, parameter :: xmol_abs_inf_tol = 1.d-9
  
  PetscReal, parameter :: pres_rel_inf_tol = 1.d-3
  PetscReal, parameter :: temp_rel_inf_tol = 1.d-3
  PetscReal, parameter :: sat_rel_inf_tol = 1.d-3
  PetscReal, parameter :: xmol_rel_inf_tol = 1.d-3
  
  PetscReal, parameter :: ref_density_w = 55.058 !kmol_water/m^3
  PetscReal, parameter :: ref_density_a = 0.0423 !kmol_air/m^3
  PetscReal, parameter :: ref_u = 83.8 !MJ/m^3
  
  !MAN optimized:
  PetscReal, parameter :: w_mass_abs_inf_tol = 1.d-5 !1.d-7 !kmol_water/sec
  PetscReal, parameter :: a_mass_abs_inf_tol = 1.d-5 !1.d-7
  PetscReal, parameter :: u_abs_inf_tol = 1.d-5 !1.d-7
                                          
  PetscReal, parameter :: residual_abs_inf_tol(3) = (/w_mass_abs_inf_tol, &
                             a_mass_abs_inf_tol, u_abs_inf_tol/)
  PetscReal, parameter :: residual_scaled_inf_tol(3) = 1.d-6

  PetscReal, parameter :: hyd_sat_abs_inf_tol = 1.d-5 !1.d-10
  !For convergence using hydrate and ice formation capability
  PetscReal, parameter :: abs_update_inf_tol(3,15) = &
    reshape([pres_abs_inf_tol,xmol_abs_inf_tol,temp_abs_inf_tol, &
             pres_abs_inf_tol,pres_abs_inf_tol,temp_abs_inf_tol, &
             pres_abs_inf_tol,hyd_sat_abs_inf_tol,temp_abs_inf_tol,  &
             pres_abs_inf_tol,999.d0,temp_abs_inf_tol, &
             pres_abs_inf_tol,999.d0,temp_abs_inf_tol, &
             pres_abs_inf_tol,hyd_sat_abs_inf_tol,temp_abs_inf_tol, &
             pres_abs_inf_tol,hyd_sat_abs_inf_tol,temp_abs_inf_tol, &
             pres_abs_inf_tol,hyd_sat_abs_inf_tol,temp_abs_inf_tol, &
             pres_abs_inf_tol,hyd_sat_abs_inf_tol,temp_abs_inf_tol, &
             pres_abs_inf_tol,xmol_abs_inf_tol,hyd_sat_abs_inf_tol, &
             hyd_sat_abs_inf_tol,hyd_sat_abs_inf_tol,temp_abs_inf_tol, &
             pres_abs_inf_tol,hyd_sat_abs_inf_tol,hyd_sat_abs_inf_tol, &
             hyd_sat_abs_inf_tol,hyd_sat_abs_inf_tol,temp_abs_inf_tol, &
             pres_abs_inf_tol,hyd_sat_abs_inf_tol,hyd_sat_abs_inf_tol, &
             hyd_sat_abs_inf_tol,hyd_sat_abs_inf_tol,hyd_sat_abs_inf_tol], &
            shape(abs_update_inf_tol)) * &
            1.d0 ! change to 0.d0 to zero tolerances
  PetscReal, parameter :: rel_update_inf_tol(3,15) = &
    reshape([pres_rel_inf_tol,xmol_rel_inf_tol,temp_rel_inf_tol, &
             pres_rel_inf_tol,pres_rel_inf_tol,temp_rel_inf_tol, &
             pres_rel_inf_tol,sat_rel_inf_tol,temp_rel_inf_tol, &
             pres_rel_inf_tol,999.d0,temp_rel_inf_tol, &
             pres_rel_inf_tol,999.d0,temp_rel_inf_tol, &
             pres_rel_inf_tol,sat_rel_inf_tol,temp_rel_inf_tol, &
             pres_rel_inf_tol,sat_rel_inf_tol,temp_rel_inf_tol, &
             pres_rel_inf_tol,sat_rel_inf_tol,temp_rel_inf_tol, &
             pres_rel_inf_tol,sat_rel_inf_tol,temp_rel_inf_tol, &
             pres_rel_inf_tol,xmol_rel_inf_tol,sat_rel_inf_tol, &
             sat_rel_inf_tol,sat_rel_inf_tol,temp_rel_inf_tol, &
             pres_rel_inf_tol,sat_rel_inf_tol,sat_rel_inf_tol, &
             sat_rel_inf_tol,sat_rel_inf_tol,temp_rel_inf_tol, &
             pres_rel_inf_tol,sat_rel_inf_tol,sat_rel_inf_tol, &
             sat_rel_inf_tol,sat_rel_inf_tol,sat_rel_inf_tol], &
            shape(rel_update_inf_tol)) * &
            1.d0 ! change to 0.d0 to zero tolerances
  allocate(this)
  allocate(this%max_change_ivar(6))
  this%max_change_ivar = [LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
                                LIQUID_MOLE_FRACTION, TEMPERATURE, &
                                GAS_SATURATION]
  allocate(this%max_change_isubvar(6))
                                   ! UNINITIALIZED_INTEGER avoids zeroing of 
                                   ! pressures not represented in phase
                                       ! 2 = air in xmol(air,liquid)
  this%max_change_isubvar = [0,0,0,2,0,0]
  this%damping_factor = -1.d0
  this%hydrate_newton_max_iter = 8
  
  call PMSubsurfaceFlowCreate(this)
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
  option%nphase = 4
  option%liquid_phase = 1  ! liquid_pressure
  option%gas_phase = 2     ! gas_pressure
  option%hydrate_phase = 3
  option%ice_phase = 4

  option%air_pressure_id = 3
  option%capillary_pressure_id = 4
  option%vapor_pressure_id = 5
  option%saturation_pressure_id = 6

  option%water_id = 1
  option%air_id = 2
  option%energy_id = 3

  option%nflowdof = 3
  option%nflowspec = 2
  option%use_isothermal = PETSC_FALSE

end subroutine PMHydrateSetFlowMode

! ************************************************************************** !

subroutine PMHydrateReadParameters(input,pm_hydrate,option)

  use Input_Aux_module
  use Option_module
  use String_module

  implicit none

  type(input_type), pointer :: input
  type(pm_hydrate_type), pointer :: pm_hydrate
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
      case('SCALE_PERM_BY_HYD_SAT')
        HYDRATE_PERM_SCALING = PETSC_TRUE
      case('EFFECTIVE_SAT_SCALING')
        HYDRATE_EFF_SAT_SCALING = PETSC_TRUE
      case('WITH_GIBBS_THOMSON')
        HYDRATE_WITH_GIBBS_THOMSON = PETSC_TRUE
      case('ADJUST_SOLUBILITY_WITHIN_GHSZ')
        HYDRATE_ADJUST_GHSZ_SOLUBILITY = PETSC_TRUE
      case('WITH_SEDIMENTATION')
        HYDRATE_WITH_SEDIMENTATION = PETSC_TRUE
      case('NO_PC')
        HYDRATE_NO_PC = PETSC_TRUE
      case('METHANOGENESIS')
        HYDRATE_WITH_METHANOGENESIS = PETSC_TRUE
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
            HYDRATE_PERM_SCALING_FUNCTION = temp_int
        end select
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
  class(pm_hydrate_type), pointer ::pm

  type(fluid_property_type), pointer :: cur_fluid_property
  type(option_type), pointer :: option

  option => realization%option

  ! initialize parameters
  allocate(pm%hydrate_parameters%diffusion_coefficient(option%nphase))
  cur_fluid_property => realization%fluid_properties
  do
    if (.not.associated(cur_fluid_property)) exit
    pm%hydrate_parameters% &
      diffusion_coefficient(cur_fluid_property%phase_id) = &
        cur_fluid_property%diffusion_coefficient
    cur_fluid_property => cur_fluid_property%next
  enddo
  ! check whether diffusion coefficients are initialized.
  if (Uninitialized(pm%hydrate_parameters% &
      diffusion_coefficient(LIQUID_PHASE))) then
    option%io_buffer = &
      UninitializedMessage('Liquid phase diffusion coefficient','')
    call PrintErrMsg(option)
  endif
  if (Uninitialized(pm%hydrate_parameters% &
      diffusion_coefficient(GAS_PHASE))) then
    option%io_buffer = &
      UninitializedMessage('Gas phase diffusion coefficient','')
    call PrintErrMsg(option)
  endif

  realization%patch%aux%hydrate%hydrate_parameter => pm%hydrate_parameters 

end subroutine PMHydrateAssignParameters

! ************************************************************************** !
subroutine PMHydrateRead(this,input)
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
    call PMSubsurfaceFlowReadSelectCase(this,input,keyword,found, &
                                        error_string,option)    
    if (found) cycle
    
    select case(trim(keyword))
      case('CENTRAL_DIFFERENCE_JACOBIAN')
        hydrate_central_diff_jacobian = PETSC_TRUE
      case('NO_STATE_TRANSITION_PRINTING')    
        hydrate_print_state_transition = PETSC_FALSE
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
        
      !man: phase change
      case('MAX_NEWTON_ITERATIONS')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,error_string)
        this%hydrate_newton_max_iter = tempreal
      case('PHASE_CHANGE_EPSILON')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,error_string)
        hydrate_phase_chng_epsilon = tempreal
      
      case('RESTRICT_STATE_CHANGE')
        hydrate_restrict_state_chng = PETSC_TRUE
      ! Tolerances
      
      case('NO_STATE_TRANSITION_OUTPUT')
        hydrate_print_state_transition = PETSC_FALSE
      
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
      case('RESIDUAL_SCALED_INF_TOL','ITOL_SCALED_RESIDUAL')
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
        this%abs_update_inf_tol(1,1:10) = tempreal
        this%abs_update_inf_tol(1,12) = tempreal
        this%abs_update_inf_tol(1,14) = tempreal
        this%abs_update_inf_tol(2,2) = tempreal
      case('TEMP_ABS_UPDATE_INF_TOL')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,error_string)
        this%abs_update_inf_tol(3,:) = tempreal
        this%abs_update_inf_tol(3,1:9) = tempreal
        this%abs_update_inf_tol(3,11) = tempreal
        this%abs_update_inf_tol(3,13) = tempreal
      case('SAT_ABS_UPDATE_INF_TOL')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,error_string)
        this%abs_update_inf_tol(2,3) = tempreal
        this%abs_update_inf_tol(2,3) = tempreal
        this%abs_update_inf_tol(2,6:9) = tempreal
        this%abs_update_inf_tol(3,10) = tempreal
        this%abs_update_inf_tol(2,11:15) = tempreal
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

      case('WINDOW_EPSILON') 
        call InputReadDouble(input,option,window_epsilon)
        call InputErrorMsg(input,option,'window epsilon',error_string)
      case('GAS_COMPONENT_FORMULA_WEIGHT')
        !geh: assuming gas component is index 2
        call InputReadDouble(input,option,fmw_comp(2))
        call InputErrorMsg(input,option,'gas component formula wt.', &
             error_string)
      case('LIQUID_COMPONENT_FORMULA_WEIGHT')
         !heeho: assuming liquid component is index 1
         call InputReadDouble(input,option,fmw_comp(1))
         call InputErrorMsg(input,option,'liquid component formula wt.', &
             error_string)
      case('TWO_PHASE_ENERGY_DOF')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'two_phase_energy_dof',error_string)
        call HydrateAuxSetEnergyDOF(word,option)
      case('MAXIMUM_PRESSURE_CHANGE')
        call InputReadDouble(input,option,hydrate_max_pressure_change)
        call InputErrorMsg(input,option,'maximum pressure change', &
                           error_string)
      case('MAX_ITERATION_BEFORE_DAMPING')
        call InputReadInt(input,option,hydrate_max_it_before_damping)
        call InputErrorMsg(input,option,'maximum iteration before damping', &
                           error_string)
      case('DAMPING_FACTOR')
        call InputReadDouble(input,option,hydrate_damping_factor)
        call InputErrorMsg(input,option,'damping factor',error_string)
        this%damping_factor = hydrate_damping_factor
#if 0        
      case('GOVERN_MAXIMUM_PRESSURE_CHANGE')
        call InputReadDouble(input,option,this%dPmax_allowable)
        call InputErrorMsg(input,option,'maximum allowable pressure change', &
                           error_string)
      case('GOVERN_MAXIMUM_TEMPERATURE_CHANGE')
        call InputReadDouble(input,option,this%dTmax_allowable)
        call InputErrorMsg(input,option, &
                           'maximum allowable temperature change', &
                           error_string)
      case('GOVERN_MAXIMUM_SATURATION_CHANGE')
        call InputReadDouble(input,option,this%dSmax_allowable)
        call InputErrorMsg(input,option,'maximum allowable saturation change', &
                           error_string)
      case('GOVERN_MAXIMUM_MOLE_FRACTION_CHANGE')
        call InputReadDouble(input,option,this%dXmax_allowable)
        call InputErrorMsg(input,option, &
                           'maximum allowable mole fraction change', &
                           error_string)
#endif
      case('DEBUG_CELL')
        call InputReadInt(input,option,hydrate_debug_cell_id)
        call InputErrorMsg(input,option,'debug cell id',error_string)
      case('NO_TEMP_DEPENDENT_DIFFUSION')
        hydrate_temp_dep_gas_air_diff = PETSC_FALSE
      case('DIFFUSE_XMASS')
        hydrate_diffuse_xmol = PETSC_FALSE
      case('HARMONIC_GAS_DIFFUSIVE_DENSITY')
        hydrate_harmonic_diff_density = PETSC_TRUE
      case('ARITHMETIC_GAS_DIFFUSIVE_DENSITY')
        hydrate_harmonic_diff_density = PETSC_FALSE
      case('IMMISCIBLE')
        hydrate_immiscible = PETSC_TRUE
      case('CHECK_MAX_DPL_LIQ_STATE_ONLY')
        hyd_chk_max_dpl_liq_state_only = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(input,keyword,'HYDRATE Mode',option)
    end select
    
  enddo  
  call InputPopBlock(input,option)

end subroutine PMHydrateRead

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
  call VecDuplicateVecsF90(this%realization%field%work,SIX_INTEGER, &
                           this%realization%field%max_change_vecs, &
                           ierr);CHKERRQ(ierr)
  ! set initial values
  do i = 1, 6
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

subroutine PMHydrateUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
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

subroutine PMHydrateCheckUpdatePre(this,line_search,X,dX,changed,ierr)
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

  implicit none
  
  class(pm_hydrate_type) :: this
  SNESLineSearch :: line_search
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  type(patch_type), pointer :: patch  
  type(grid_type), pointer :: grid 
  type(option_type), pointer :: option 
  type(hydrate_auxvar_type), pointer :: hyd_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: offset
  PetscInt :: pgas_index, xmol_index, pw_index
  PetscInt :: saturation_index
  PetscReal :: temp_real

  PetscReal, parameter :: ALMOST_ZERO = 1.d-10
  PetscReal, parameter :: ALMOST_ONE = 1.d0-ALMOST_ZERO

  PetscReal, pointer :: X_p(:),dX_p(:)

  ! MAN: OLD
  PetscReal, pointer :: r_p(:)
  type(field_type), pointer :: field
  PetscInt :: liquid_pressure_index, gas_pressure_index, air_pressure_index
  PetscInt :: temperature_index
  PetscInt :: lid, gid, apid, cpid, vpid, spid
  PetscReal :: liquid_pressure0, liquid_pressure1, del_liquid_pressure
  PetscReal :: gas_pressure0, gas_pressure1, del_gas_pressure
  PetscReal :: air_pressure0, air_pressure1, del_air_pressure
  PetscReal :: temperature0, temperature1, del_temperature
  PetscReal :: saturation0, saturation1, del_saturation
  PetscReal :: xmol0, xmol1, del_xmol
  PetscReal :: max_saturation_change = 0.125d0
  PetscReal :: max_temperature_change = 10.d0
  PetscReal :: min_pressure
  PetscReal :: scale, temp_scale
  PetscReal, parameter :: tolerance = 0.99d0
  PetscReal, parameter :: initial_scale = 1.d0
  SNES :: snes
  PetscInt :: newton_iteration
  ! MAN: END OLD

  call VecGetArrayF90(dX,dX_p,ierr); CHKERRQ(ierr)
  call VecGetArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

  grid => this%realization%patch%grid
  patch => this%realization%patch
  option => this%realization%option
  hyd_auxvars => this%realization%patch%aux%Hydrate%auxvars
  global_auxvars => this%realization%patch%aux%Global%auxvars
  
  changed = PETSC_TRUE 

  ! MAN: OLD
  field => this%realization%field

  spid = option%saturation_pressure_id
  apid = option%air_pressure_id

  call SNESLineSearchGetSNES(line_search,snes,ierr)
  call SNESGetIterationNumber(snes,newton_iteration,ierr)

  hydrate_allow_state_change = PETSC_FALSE

  ! MAN: END OLD
  if (this%check_post_convergence) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      offset = (local_id-1)*option%nflowdof
      select case(global_auxvars(ghosted_id)%istate)
        case(L_STATE)
          xmol_index = offset + HYDRATE_L_STATE_X_MOLE_DOF
          pw_index = offset + HYDRATE_LIQUID_PRESSURE_DOF
          if (X_p(xmol_index) - dX_p(xmol_index) < 0.d0) then
            dX_p(xmol_index) = X_p(xmol_index)
            changed = PETSC_TRUE
          endif
          if (X_p(pw_index)- dX_p(pw_index) <= 0.d0) then
           dX_p(pw_index) = X_p(pw_index) - ALMOST_ZERO
           changed = PETSC_TRUE
          endif
        case(G_STATE)
         pgas_index = offset + HYDRATE_GAS_PRESSURE_DOF
         if (X_p(pgas_index)- dX_p(pgas_index) <= 0.d0) then
           dX_p(pgas_index) = X_p(pgas_index) - ALMOST_ZERO
           changed = PETSC_TRUE
         endif
        case(GA_STATE)
          pgas_index = offset + HYDRATE_GAS_PRESSURE_DOF
          if (X_p(pgas_index) - dX_p(pgas_index) < &
                  hyd_auxvars(ZERO_INTEGER,ghosted_id)% &
                  pres(option%saturation_pressure_id)) then
            dX_p(pgas_index) = X_p(pgas_index) - &
                    hyd_auxvars(ZERO_INTEGER,ghosted_id)% &
                    pres(option%saturation_pressure_id)
            changed = PETSC_TRUE
          endif
          if (hydrate_immiscible) then
            saturation_index = offset + HYDRATE_GAS_SATURATION_DOF
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
        case(L_STATE)
          xmol_index = offset + HYDRATE_L_STATE_X_MOLE_DOF
          pw_index = offset + HYDRATE_LIQUID_PRESSURE_DOF
          if (X_p(xmol_index) - dX_p(xmol_index) < 0.d0) then
            dX_p(xmol_index) = X_p(xmol_index)
            changed = PETSC_TRUE
          endif
        case(GA_STATE)
          pgas_index = offset + HYDRATE_GAS_PRESSURE_DOF
          if (X_p(pgas_index) - dX_p(pgas_index) < &
                  hyd_auxvars(ZERO_INTEGER,ghosted_id)% &
                  pres(option%saturation_pressure_id)) then
            dX_p(pgas_index) = X_p(pgas_index) - &
                    hyd_auxvars(ZERO_INTEGER,ghosted_id)% &
                    pres(option%saturation_pressure_id)
            changed = PETSC_TRUE
          endif
          if (hydrate_immiscible) then
            saturation_index = offset + HYDRATE_GAS_SATURATION_DOF
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
    if (hydrate_max_it_before_damping > 0 .and. &
        newton_iteration > hydrate_max_it_before_damping) then
      scale = hydrate_damping_factor
    endif

#define LIMIT_MAX_PRESSURE_CHANGE
#define LIMIT_MAX_SATURATION_CHANGE
    ! scaling
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      offset = (local_id-1)*option%nflowdof
      temp_scale = 1.d0
      select case(global_auxvars(ghosted_id)%istate)
        case(L_STATE)
          liquid_pressure_index  = offset + HYDRATE_LIQUID_PRESSURE_DOF
          temperature_index  = offset + HYDRATE_ENERGY_DOF
          dX_p(liquid_pressure_index) = dX_p(liquid_pressure_index) * &
                                        HYDRATE_PRESSURE_SCALE
          temp_scale = 1.d0
          del_liquid_pressure = dX_p(liquid_pressure_index)
          liquid_pressure0 = X_p(liquid_pressure_index)
          liquid_pressure1 = liquid_pressure0 - del_liquid_pressure
          del_temperature = dX_p(temperature_index)
          temperature0 = X_p(temperature_index)
          temperature1 = temperature0 - del_temperature
#ifdef LIMIT_MAX_PRESSURE_CHANGE
          if (dabs(del_liquid_pressure) > hydrate_max_pressure_change) then
            temp_real = dabs(hydrate_max_pressure_change/del_liquid_pressure)
            temp_scale = min(temp_scale,temp_real)
          endif
#endif
!LIMIT_MAX_PRESSURE_CHANGE
        case(GA_STATE)
          gas_pressure_index = offset + HYDRATE_GAS_PRESSURE_DOF
!        air_pressure_index = offset + 2
          saturation_index = offset + HYDRATE_GAS_SATURATION_DOF
          temperature_index  = offset + HYDRATE_ENERGY_DOF
          dX_p(gas_pressure_index) = dX_p(gas_pressure_index) * &
                                     HYDRATE_PRESSURE_SCALE
          if (hydrate_2ph_energy_dof == HYDRATE_AIR_PRESSURE_INDEX) then
            air_pressure_index = offset + HYDRATE_ENERGY_DOF
            dX_p(air_pressure_index) = dX_p(air_pressure_index) * &
                                       HYDRATE_PRESSURE_SCALE
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
          if (dabs(del_gas_pressure) > hydrate_max_pressure_change) then
            temp_real = dabs(hydrate_max_pressure_change/del_gas_pressure)
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
        case(G_STATE)
          gas_pressure_index = offset + HYDRATE_GAS_PRESSURE_DOF
          air_pressure_index = offset + HYDRATE_G_STATE_AIR_PRESSURE_DOF
          dX_p(gas_pressure_index) = dX_p(gas_pressure_index) * &
                                     HYDRATE_PRESSURE_SCALE
          dX_p(air_pressure_index) = dX_p(air_pressure_index) * &
                                     HYDRATE_PRESSURE_SCALE
      end select
      scale = min(scale,temp_scale)
    enddo

    temp_scale = scale
    call MPI_Allreduce(temp_scale,scale,ONE_INTEGER_MPI, &
                       MPI_DOUBLE_PRECISION, &
                       MPI_MIN,option%mycomm,ierr)

    if (scale < 0.9999d0) then
      dX_p = scale*dX_p
    endif
  endif

  call VecRestoreArrayF90(dX,dX_p,ierr); CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

end subroutine PMHydrateCheckUpdatePre

! ************************************************************************** !

subroutine PMHydrateCheckUpdatePost(this,line_search,X0,dX,X1,dX_changed, &
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
  SNESLineSearch :: line_search
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
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum2_p(:)
  PetscInt :: local_id, ghosted_id, natural_id
  PetscInt :: offset, ival, idof, itol
  PetscReal :: R, A, R_A
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
  character(len=MAXSTRINGLENGTH) :: string
  character(len=14), parameter :: state_string(15) = &
    ['Liquid State  ','Gas State     ','Hydrate State ','Ice State     ', &
     'GA State      ','HG State      ','HA State      ','HI State      ', &
     'GI State      ','AI State      ','HGA State     ','HAI State     ', &
     'HGI State     ','GAI State     ','Quad State    ']
  character(len=17), parameter :: dof_string(3,15) = &
    reshape(['Liquid Pressure  ','Air Mole Fraction','Temperature      ', &
             'Gas Pressure     ','Air Pressure     ','Temperature      ', &
             'Gas Pressure     ','Gas Saturation   ','Temperature      ', &
             'Placeholder      ','Placeholder      ','Placeholder      ', &
             'Placeholder      ','Placeholder      ','Placeholder      ', &
             'Placeholder      ','Placeholder      ','Placeholder      ', &
             'Placeholder      ','Placeholder      ','Placeholder      ', &
             'Placeholder      ','Placeholder      ','Placeholder      ', &
             'Placeholder      ','Placeholder      ','Placeholder      ', &
             'Placeholder      ','Placeholder      ','Placeholder      ', &
             'Placeholder      ','Placeholder      ','Placeholder      ', &
             'Placeholder      ','Placeholder      ','Placeholder      ', &
             'Placeholder      ','Placeholder      ','Placeholder      ', &
             'Placeholder      ','Placeholder      ','Placeholder      ', &
             'Placeholder      ','Placeholder      ','Placeholder      '], &
             shape(dof_string))
  character(len=15), parameter :: tol_string(MAX_INDEX) = &
    ['Absolute Update','Relative Update','Residual       ','Scaled Residual']
  
  patch => this%realization%patch
  option => this%realization%option
  field => this%realization%field
  grid => patch%grid
  global_auxvars => patch%aux%Global%auxvars

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
    call VecRestoreArrayReadF90(field%flow_accum2,accum2_p,ierr);CHKERRQ(ierr)
  
    this%converged_flag(:,:,RESIDUAL_INDEX) = converged_abs_residual_flag(:,:)
    this%converged_real(:,:,RESIDUAL_INDEX) = converged_abs_residual_real(:,:)
    this%converged_cell(:,:,RESIDUAL_INDEX) = converged_abs_residual_cell(:,:)
    this%converged_flag(:,:,SCALED_RESIDUAL_INDEX) = &
                                       converged_scaled_residual_flag(:,:)
    this%converged_real(:,:,SCALED_RESIDUAL_INDEX) = &
                                       converged_scaled_residual_real(:,:)
    this%converged_cell(:,:,SCALED_RESIDUAL_INDEX) = &
                                       converged_scaled_residual_cell(:,:)
    mpi_int = 9*MAX_INDEX
    ! do not perform an all reduce on cell id as this info is not printed 
    ! in parallel
    call MPI_Allreduce(MPI_IN_PLACE,this%converged_flag,mpi_int, &
                       MPI_LOGICAL,MPI_LAND,option%mycomm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,this%converged_real,mpi_int, &
                       MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
  
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
              if (option%mycommsize == 1) then
                string = trim(string) // ' (' // &
                  trim(StringFormatInt(this%converged_cell(idof,istate,itol))) &
                  // ')' 
              endif
              string = trim(string) // ' : ' // &
                StringFormatDouble(this%converged_real(idof,istate,itol))
              call OptionPrint(string,option)
            endif
          endif
        enddo
      enddo
    enddo
    if (this%logging_verbosity > 0 .and. it > 0 .and. &
        option%convergence == CONVERGENCE_CONVERGED) then
      string = '   Converged'
      call OptionPrint(string,option)
      write(string,'(4x," R:",9es8.1)') this%converged_real(:,:,RESIDUAL_INDEX)
      call OptionPrint(string,option)
      write(string,'(4x,"SR:",9es8.1)') &
        this%converged_real(:,:,SCALED_RESIDUAL_INDEX)
      call OptionPrint(string,option)
      write(string,'(4x,"AU:",9es8.1)') &
        this%converged_real(:,:,ABS_UPDATE_INDEX)
      call OptionPrint(string,option)
      write(string,'(4x,"RU:",9es8.1)') &
        this%converged_real(:,:,REL_UPDATE_INDEX)
      call OptionPrint(string,option)
    endif
  
    if (it >= this%hydrate_newton_max_iter) then
      option%convergence = CONVERGENCE_CUT_TIMESTEP
    
      if (this%logging_verbosity > 0) then
        string = '    Exceeded Hydrate Mode Max Newton Iterations'
        call OptionPrint(string,option)
      endif
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
  use Variables_module, only : LIQUID_PRESSURE, LIQUID_MOLE_FRACTION, &
                               TEMPERATURE, GAS_PRESSURE, AIR_PRESSURE, &
                               GAS_SATURATION
  implicit none
  
  class(pm_hydrate_type) :: this
  
  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscReal :: max_change_local(6)
  PetscReal :: max_change_global(6)
  PetscReal :: max_change
  PetscInt :: i, j
  PetscInt :: local_id, ghosted_id

  
  PetscErrorCode :: ierr
  
  realization => this%realization
  option => realization%option
  field => realization%field
  grid => realization%patch%grid

  global_auxvars => realization%patch%aux%global%auxvars

  max_change_global = 0.d0
  max_change_local = 0.d0
  
  ! max change variables: [LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
  !                        LIQUID_MOLE_FRACTION, TEMPERATURE, GAS_SATURATION]
  do i = 1, 6
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
  call MPI_Allreduce(max_change_local,max_change_global,SIX_INTEGER, &
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
  ! print them out
  if (OptionPrintToScreen(option)) then
    write(*,'("  --> max chng: dpl= ",1pe12.4, " dpg= ",1pe12.4,&
      & " dpa= ",1pe12.4,/,15x," dxa= ",1pe12.4,"  dt= ",1pe12.4,&
      & " dsg= ",1pe12.4)') &
      max_change_global(1:6)
  endif
  if (OptionPrintToFile(option)) then
    write(option%fid_out,'("  --> max chng: dpl= ",1pe12.4, " dpg= ",1pe12.4,&
      & " dpa= ",1pe12.4,/,15x," dxa= ",1pe12.4,"  dt= ",1pe12.4, &
      & " dsg= ",1pe12.4)') &
      max_change_global(1:6)
  endif

  ! max change variables: [LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
  !                        LIQUID_MOLE_FRACTION, TEMPERATURE, GAS_SATURATION]
  ! ignore air pressure as it jumps during phase change
  this%max_pressure_change = maxval(max_change_global(1:2))
  this%max_xmol_change = max_change_global(4)
  this%max_temperature_change = max_change_global(5)
  this%max_saturation_change = max_change_global(6)
  
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

  character(len=MAXWORDLENGTH) :: word
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
  use Variables_module, only : STATE

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
  use Variables_module, only : STATE

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
