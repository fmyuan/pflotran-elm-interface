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
    PetscReal :: liquid_residual_infinity_tol
    PetscReal :: gas_equation_infinity_tol
    PetscReal :: max_allow_rel_liq_pres_chang_ni
    PetscReal :: max_allow_rel_gas_sat_change_ni
    ! the below is set automatically to -log10(max_allow_rel_gas_sat_change_ni)
    PetscReal :: neg_log10_rel_gas_sat_change_ni 
    PetscReal :: gas_sat_thresh_force_ts_cut
    PetscReal :: min_liq_pres_force_ts_cut
    PetscReal :: gas_sat_thresh_force_extra_ni
    PetscReal :: max_allow_gas_sat_change_ts
    PetscReal :: max_allow_liq_pres_change_ts
    PetscReal :: gas_sat_change_ts_governor
    PetscReal :: liq_pres_change_ts_governor
    PetscReal :: gas_sat_gov_switch_abs_to_rel
    PetscReal :: minimum_timestep_size
    PetscBool :: convergence_test_both
    class(pm_wipp_srcsink_type), pointer :: pmwss_ptr
    Vec :: stored_residual_vec
    PetscInt :: convergence_flags(MIN_GAS_PRES)
    ! store maximum quantities for the above
    PetscReal :: convergence_reals(MIN_GAS_PRES)
    character(len=MAXWORDLENGTH) :: alpha_dataset_name
    character(len=MAXWORDLENGTH) :: elevation_dataset_name
    PetscReal :: rotation_angle ! radians
    PetscReal :: rotation_origin(3)
    PetscReal :: rotation_ceiling  ! elevation in z
    PetscReal :: rotation_basement ! elevation in z
    character(len=MAXWORDLENGTH), pointer :: rotation_region_names(:)
    PetscInt, pointer :: auto_pressure_material_ids(:)
    PetscReal :: auto_pressure_rho_b0
    PetscReal :: auto_pressure_c_b
    PetscReal :: auto_pressure_Pb_ref
    PetscReal :: auto_pressure_Pb_0
    PetscReal :: auto_press_shallow_origin(3)
    PetscReal :: linear_system_scaling_factor
    PetscBool :: scale_linear_system
    Vec :: scaling_vec
    ! When reading Dirichlet 2D Flared BC
    PetscInt, pointer :: dirichlet_dofs_ghosted(:) ! this array is zero-based indexing
    ! int_array has a natural_id and 1,2, or 3 that indicates
    ! pressure, satruation, or both to be zerod in the residual
    PetscInt, pointer :: dirichlet_dofs_ints(:,:)
    PetscInt, pointer :: dirichlet_dofs_local(:) ! this array is zero-based indexing
 
  contains
    procedure, public :: ReadSimulationBlock => PMWIPPFloRead
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
    procedure, public :: CheckConvergence => PMWIPPFloCheckConvergence
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
  this%header = 'WIPP IMMISCIBLE MULTIPHASE FLOW'

  this%check_post_convergence = PETSC_TRUE

  ! defaults from BRAGFLO input deck or recommended values from user manual
  this%liquid_residual_infinity_tol = 1.d-2
  this%gas_equation_infinity_tol = 1.d-2
  this%max_allow_rel_liq_pres_chang_ni = 1.d-2
  this%max_allow_rel_gas_sat_change_ni = 1.d-3
  ! the below is set automatically to -log10(max_allow_rel_gas_sat_change_ni)
  this%neg_log10_rel_gas_sat_change_ni = UNINITIALIZED_DOUBLE
  this%gas_sat_thresh_force_ts_cut = 0.20d0   ! [-]
  this%min_liq_pres_force_ts_cut = -1.0d8   ! [Pa]
  this%gas_sat_thresh_force_extra_ni = 1.0d-3  ! [-]
  this%max_allow_gas_sat_change_ts = 1.d0    ! [-]
  this%max_allow_liq_pres_change_ts = 1.d7   ! [Pa]
  this%gas_sat_change_ts_governor = 3.d-1    ! [-]
  this%liq_pres_change_ts_governor = 5.d5    ! [Pa]
  this%gas_sat_gov_switch_abs_to_rel = 0.01d0   ! [-]
  this%minimum_timestep_size = 8.64d-4  ! [sec]
  this%stored_residual_vec = PETSC_NULL_VEC
  this%alpha_dataset_name = ''
  this%elevation_dataset_name = ''
  this%rotation_angle = UNINITIALIZED_DOUBLE
  this%rotation_origin = UNINITIALIZED_DOUBLE
  this%rotation_ceiling = UNINITIALIZED_DOUBLE
  this%rotation_basement = UNINITIALIZED_DOUBLE
  nullify(this%rotation_region_names)
  nullify(this%auto_pressure_material_ids)
  this%auto_pressure_rho_b0 = 1220.d0
  this%auto_pressure_c_b = 3.1d-10
  this%auto_pressure_Pb_ref = 101325.d0
  this%auto_pressure_Pb_0 = UNINITIALIZED_DOUBLE !make user put this in
  this%auto_press_shallow_origin = UNINITIALIZED_DOUBLE !this will default to dip rotation origin later
  this%linear_system_scaling_factor = 1.d7
  this%scale_linear_system = PETSC_TRUE
  this%scaling_vec = PETSC_NULL_VEC
  nullify(this%dirichlet_dofs_ghosted)
  nullify(this%dirichlet_dofs_ints)
  nullify(this%dirichlet_dofs_local)
  this%convergence_test_both = PETSC_TRUE
  this%convergence_flags = 0
  this%convergence_reals = 0.d0

end subroutine PMWIPPFloInitObject

! ************************************************************************** !

subroutine PMWIPPFloRead(this,input)
  ! 
  ! Read WIPP FLOW input block
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  !
  use WIPP_Flow_module
  use WIPP_Flow_Aux_module
  use Input_Aux_module
  use String_module
  use Option_module
  use Utility_module

  implicit none
  
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: keyword, word, word2
  class(pm_wippflo_type) :: this
  type(option_type), pointer :: option
  PetscReal :: tempreal
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=MAXSTRINGLENGTH), pointer :: strings(:)
  PetscBool :: found
  PetscInt :: icount
  PetscInt :: temp_int
  ! temp_int_array has a natural_id and 1,2, or 3 that indicates
  ! pressure, satruation, or both to be zerod in the residual
  PetscInt, parameter :: max_dirichlet_bc = 1000 ! capped at 1000
  PetscInt :: temp_int_array(2,max_dirichlet_bc) 
  
  option => this%option

  error_string = 'WIPP Flow Options'
  
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
      case('LIQUID_RESIDUAL_INFINITY_TOL')
        call InputReadDouble(input,option,this%liquid_residual_infinity_tol)
        call InputErrorMsg(input,option,keyword,error_string)
      case('GAS_RESIDUAL_INFINITY_TOL')
        call InputReadDouble(input,option,this%gas_equation_infinity_tol)
        call InputErrorMsg(input,option,keyword,error_string)
      case('MAX_ALLOW_REL_LIQ_PRES_CHANG_NI')
        call InputReadDouble(input,option,this%max_allow_rel_liq_pres_chang_ni)
        call InputErrorMsg(input,option,keyword,error_string)
        ! no units conversion since it is relative
      case('MAX_ALLOW_REL_GAS_SAT_CHANGE_NI')
        call InputReadDouble(input,option,this%max_allow_rel_gas_sat_change_ni)
        call InputErrorMsg(input,option,keyword,error_string)
      case('GAS_COMPONENT_FORMULA_WEIGHT')
        call InputReadDouble(input,option,fmw_comp(2))
        call InputErrorMsg(input,option,keyword,error_string)
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
      case('REL_LIQ_PRESSURE_PERTURBATION')
        call InputReadDouble(input,option,wippflo_pres_rel_pert)
        call InputErrorMsg(input,option,keyword,error_string)
        ! no units conversion since it is relative
      case('MIN_LIQ_PRESSURE_PERTURBATION')
        call InputReadDouble(input,option,wippflo_pres_min_pert)
        call InputErrorMsg(input,option,keyword,error_string)
        call InputReadAndConvertUnits(input,wippflo_pres_min_pert, &
                                      'Pa',keyword,option)
      case('REL_GAS_SATURATION_PERTURBATION')
        call InputReadDouble(input,option,wippflo_sat_rel_pert)
        call InputErrorMsg(input,option,keyword,error_string)
      case('MIN_GAS_SATURATION_PERTURBATION')
        call InputReadDouble(input,option,wippflo_sat_min_pert)
        call InputErrorMsg(input,option,keyword,error_string)
      case('GAS_SAT_THRESH_FORCE_TS_CUT')
        call InputReadDouble(input,option,this%gas_sat_thresh_force_ts_cut)
        call InputErrorMsg(input,option,keyword,error_string)
      case('GAS_SAT_THRESH_FORCE_EXTRA_NI')
        call InputReadDouble(input,option,this%gas_sat_thresh_force_extra_ni)
        call InputErrorMsg(input,option,keyword,error_string)
      case('MIN_LIQ_PRES_FORCE_TS_CUT')
        call InputReadDouble(input,option,this%min_liq_pres_force_ts_cut)
        call InputErrorMsg(input,option,keyword,error_string)
        call InputReadAndConvertUnits(input,this%min_liq_pres_force_ts_cut, &
                                      'Pa',keyword,option)
      case('MAX_ALLOW_GAS_SAT_CHANGE_TS')
        call InputReadDouble(input,option,this%max_allow_gas_sat_change_ts)
        call InputErrorMsg(input,option,keyword,error_string)
      case('MAX_ALLOW_LIQ_PRES_CHANGE_TS')
        call InputReadDouble(input,option,this%max_allow_liq_pres_change_ts)
        call InputErrorMsg(input,option,keyword,error_string)
        ! units conversion since it is absolute
        call InputReadAndConvertUnits(input,this%max_allow_liq_pres_change_ts, &
                                      'Pa',keyword,option)
      case('GAS_SAT_CHANGE_TS_GOVERNOR')
        call InputReadDouble(input,option,this%gas_sat_change_ts_governor)
        call InputErrorMsg(input,option,keyword,error_string)
      case('LIQ_PRES_CHANGE_TS_GOVERNOR')
        call InputReadDouble(input,option,this%liq_pres_change_ts_governor)
        call InputErrorMsg(input,option,keyword,error_string)
        ! units conversion since it is absolute
        call InputReadAndConvertUnits(input,this%liq_pres_change_ts_governor, &
                                      'Pa',keyword,option)
      case('GAS_SAT_GOV_SWITCH_ABS_TO_REL')
        call InputReadDouble(input,option,this%gas_sat_gov_switch_abs_to_rel)
        call InputErrorMsg(input,option,keyword,error_string)
      case('MINIMUM_TIMESTEP_SIZE')
        call InputReadDouble(input,option,this%minimum_timestep_size)
        call InputErrorMsg(input,option,keyword,error_string)
        call InputReadAndConvertUnits(input,this%minimum_timestep_size, &
                                      'sec',keyword,option)
      case('CONVERGENCE_TEST')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,keyword,error_string)
        call StringToUpper(word)
        select case(word)
          case('BOTH')
            this%convergence_test_both = PETSC_TRUE
          case('EITHER')
            this%convergence_test_both = PETSC_FALSE
          case default
            call InputKeywordUnrecognized(input,keyword, &
                           trim(error_string)//','//keyword,option)
        end select
      case('RESIDUAL_TEST')
        wippflo_residual_test = PETSC_TRUE
      case('RESIDUAL_TEST_CELL')
        call InputReadInt(input,option,wippflo_residual_test_cell)
        call InputErrorMsg(input,option,keyword,error_string)
      case('JACOBIAN_TEST')
        wippflo_jacobian_test = PETSC_TRUE
      case('JACOBIAN_TEST_RDOF')
        call InputReadInt(input,option,wippflo_jacobian_test_rdof)
        call InputErrorMsg(input,option,keyword,error_string)
      case('JACOBIAN_TEST_XDOF')
        call InputReadInt(input,option,wippflo_jacobian_test_xdof)
        call InputErrorMsg(input,option,keyword,error_string)
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
      case('HARMONIC_PERMEABILITY_ONLY')
        wippflo_use_lumped_harm_flux = PETSC_FALSE
      case('DEFAULT_ALPHA')
        wippflo_default_alpha = PETSC_TRUE
      case('ALPHA_DATASET')
        call InputReadNChars(input,option,this%alpha_dataset_name,&
                             MAXWORDLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_string)
      case('ELEVATION_DATASET')
        call InputReadNChars(input,option,this%elevation_dataset_name,&
                             MAXWORDLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_string)
      case('DIP_ROTATION_ANGLE')
        call InputReadDouble(input,option,this%rotation_angle)
        call InputErrorMsg(input,option,keyword,error_string)
        ! convert to radians
        ! acos(-1) = pi
        this%rotation_angle = this%rotation_angle * acos(-1.d0) / 180.d0
      case('DIP_ROTATION_ORIGIN')
        call InputReadNDoubles(input,option,this%rotation_origin,THREE_INTEGER)
        call InputErrorMsg(input,option,keyword,error_string)
      case('DIP_ROTATION_CEILING')
        call InputReadDouble(input,option,this%rotation_ceiling)
        call InputErrorMsg(input,option,keyword,error_string)
      case('DIP_ROTATION_BASEMENT')
        call InputReadDouble(input,option,this%rotation_basement)
        call InputErrorMsg(input,option,keyword,error_string)
      case('DIP_ROTATION_REGIONS')
        strings => StringSplit(adjustl(input%buf),' ')
        allocate(this%rotation_region_names(size(strings)))
        this%rotation_region_names(:) = strings(:)
        deallocate(strings)
        nullify(strings)
      case('AUTO_PRESSURE_MATERIAL_IDS')
        strings => StringSplit(adjustl(input%buf),' ')
        allocate(this%auto_pressure_material_ids(size(strings)))
        do temp_int = 1, size(strings)
          call InputReadInt(strings(temp_int),option, &
                            this%auto_pressure_material_ids(temp_int), &
                            input%ierr)
          call InputErrorMsg(input,option,keyword,error_string)
        enddo
        deallocate(strings)
        nullify(strings)
      case('AUTO_PRESSURE_RHO_B0')
        call InputReadDouble(input,option,this%auto_pressure_rho_b0)
        call InputErrorMsg(input,option,keyword,error_string)
      case('AUTO_PRESSURE_C_B')
        call InputReadDouble(input,option,this%auto_pressure_c_b)
        call InputErrorMsg(input,option,keyword,error_string)
      case('AUTO_PRESSURE_PB_REF')
        call InputReadDouble(input,option,this%auto_pressure_Pb_ref)
        call InputErrorMsg(input,option,keyword,error_string)
      case('AUTO_PRESSURE_PB_0')
        call InputReadDouble(input,option,this%auto_pressure_Pb_0)
        call InputErrorMsg(input,option,keyword,error_string)
      case('AUTO_PRESS_SHALLOW_ORIGIN')
        call InputReadNDoubles(input,option,this%auto_press_shallow_origin,THREE_INTEGER)
        call InputErrorMsg(input,option,keyword,error_string)
      case('JACOBIAN_PRESSURE_DERIV_SCALE')
        call InputReadDouble(input,option,this%linear_system_scaling_factor)
        call InputErrorMsg(input,option,keyword,error_string)
      case('SCALE_JACOBIAN')
        this%scale_linear_system = PETSC_TRUE
      case('DO_NOT_SCALE_JACOBIAN')
        this%scale_linear_system = PETSC_FALSE
      case('2D_FLARED_DIRICHLET_BCS')
        icount = 0
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,keyword)
          if (InputCheckExit(input,option)) exit
          if (icount+1 > max_dirichlet_bc) then
            option%io_buffer = 'Must increase size of "max_dirichlet_bc" & 
              &in PMWIPPFloRead'
            call PrintErrMsg(option)
          endif
          call InputReadInt(input,option,temp_int)
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'pressure', &
                             '2D_FLARED_DIRICHLET_BCS')
          call InputReadWord(input,option,word2,PETSC_TRUE)
          call InputErrorMsg(input,option,'saturation', &
                             '2D_FLARED_DIRICHLET_BCS')
                             
          if (StringYesNoOther(word) == STRING_YES .and. &
              StringYesNoOther(word2) == STRING_YES) then
            icount = icount + 1
            temp_int_array(1,icount) = temp_int
            temp_int_array(2,icount) = 3
          else if (StringYesNoOther(word) == STRING_YES .and. &
                   StringYesNoOther(word2) == STRING_NO) then
            icount = icount + 1
            temp_int_array(1,icount) = temp_int
            temp_int_array(2,icount) = 1
          else if (StringYesNoOther(word) == STRING_NO .and. &
                   StringYesNoOther(word2) == STRING_YES) then
            icount = icount + 1
            temp_int_array(1,icount) = temp_int
            temp_int_array(2,icount) = 2
          endif
        enddo
        allocate(this%dirichlet_dofs_ints(2,icount))
        this%dirichlet_dofs_ints = temp_int_array(1:2,1:icount)
      case default
        call InputKeywordUnrecognized(input,keyword,'WIPP Flow Mode',option)
    end select
  enddo  
  call InputPopBlock(input,option)
  
  ! Check that gas_sat_thresh_force_extra_ni is smaller than 
  ! gas_sat_thresh_force_ts_cut
  if (this%gas_sat_thresh_force_extra_ni > &
      this%gas_sat_thresh_force_ts_cut) then
    option%io_buffer = 'The value of GAS_SAT_THRESH_FORCE_TS_CUT must &
                       &be larger than GAS_SAT_THRESH_FORCE_EXTRA_NI.'
    call PrintErrMsg(option)
  endif
  ! Check the sign of given variables
  if (this%gas_sat_thresh_force_ts_cut < 0.d0) then
    option%io_buffer = 'The value of GAS_SAT_THRESH_FORCE_TS_CUT &
                       &must be positive.'
    call PrintErrMsg(option)
  endif
  if (this%min_liq_pres_force_ts_cut > 0.d0) then
    option%io_buffer = 'The value of MIN_LIQ_PRES_FORCE_TS_CUT &
                       &must be negative.'
    call PrintErrMsg(option)
  endif
  if (this%gas_sat_thresh_force_extra_ni < 0.d0) then
    option%io_buffer = 'The value of GAS_SAT_THRESH_FORCE_EXTRA_NI &
                       &must be positive.'
    call PrintErrMsg(option)
  endif
  ! always calculate neg_log10_rel_gas_sat_change_ni automatically
  this%neg_log10_rel_gas_sat_change_ni = &
    -1.d0*log10(this%max_allow_rel_gas_sat_change_ni)
   
end subroutine PMWIPPFloRead

! ************************************************************************** !

recursive subroutine PMWIPPFloInitializeRun(this)
  ! 
  ! Initializes the WIPP_FLOW mode run.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17

  use Realization_Base_class
  use Patch_module
  use WIPP_Flow_module, only : WIPPFloUpdateAuxVars
  use WIPP_Flow_Aux_module
  use Input_Aux_module
  use Dataset_Base_class
  use Dataset_Common_HDF5_class
  use Dataset_module
  use Field_module
  use Grid_module
  use HDF5_module
  use Option_module
  use Discretization_module
  use Region_module
  
  implicit none
  
  class(pm_wippflo_type) :: this
  
  PetscInt :: i
  PetscErrorCode :: ierr
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: block_string
  class(dataset_base_type), pointer :: dataset
  type(field_type), pointer :: field
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  type(grid_type), pointer :: grid
  type(region_type), pointer :: region
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscReal, pointer :: work_p(:)
  PetscReal, pointer :: work_loc_p(:)
  PetscReal, pointer :: flow_xx_p(:)
  PetscReal :: h
  PetscReal :: x, z
  PetscInt :: idof, icell, iregion, icount, jcount
  PetscInt :: ghosted_id, local_id, natural_id
  PetscInt :: imat, nmat_id, ndof
  ! for auto pressure
  PetscReal :: rhob0
  PetscReal :: cb
  PetscReal :: Pbref
  PetscReal :: zref
  PetscReal :: zref2
  PetscReal :: Pb0
  PetscReal :: ze
  PetscReal :: rhobref
  PetscReal :: Phiref
  PetscReal :: Phiref2
  PetscReal :: rhob
  PetscReal :: Pb
  PetscBool :: found
  PetscReal, parameter :: gravity = 9.80665d0

  patch => this%realization%patch
  grid => patch%grid
  field => this%realization%field
  option => this%option

  ! need to allocate vectors for max change
  call VecDuplicateVecsF90(field%work,SIX_INTEGER,field%max_change_vecs, &
                           ierr);CHKERRQ(ierr)
  ! set initial values
  do i = 1, 3
    call RealizationGetVariable(this%realization,field%max_change_vecs(i), &
                                this%max_change_ivar(i),ZERO_INTEGER)
  enddo

  ! call parent implementation
  call PMSubsurfaceFlowInitializeRun(this)
  
  ! look for WIPP_SOURCE_SINK block 
  input => InputCreate(IN_UNIT,option%input_filename,option)
  block_string = 'WIPP_SOURCE_SINK'
  call InputFindStringInFile(input,option,block_string)
  if (input%ierr == 0 .and. wippflo_use_gas_generation) then
    this%pmwss_ptr => PMWSSCreate()
    this%pmwss_ptr%option => option
    call this%pmwss_ptr%ReadPMBlock(input)
  endif
  ! call setup/initialization of all WIPP process models
  if (associated(this%pmwss_ptr)) then
    call PMWSSSetRealization(this%pmwss_ptr,this%realization)
    call this%pmwss_ptr%Setup()
    call this%pmwss_ptr%InitializeRun()
  endif

  ! read in alphas
  if (len_trim(this%alpha_dataset_name) > 0) then
    string = 'BRAGFLO ALPHA Dataset'
    dataset => DatasetBaseGetPointer(this%realization%datasets, &
                                     this%alpha_dataset_name, &
                                     string,option)
    select type(d => dataset)
      class is(dataset_common_hdf5_type)
        string2 = ''
        string = d%hdf5_dataset_name
        call HDF5ReadCellIndexedRealArray(this%realization, &
                                          field%work, &
                                          d%filename,&
                                          string2, &
                                          string, &
                                          d%realization_dependent)
        wippflo_auxvars => patch%aux%WIPPFlo%auxvars
        call this%realization%comm1%GlobalToLocal(field%work,field%work_loc)
        call VecGetArrayReadF90(field%work_loc,work_loc_p,ierr);CHKERRQ(ierr)
        do ghosted_id = 1, grid%ngmax
          do idof = 0, option%nflowdof
            wippflo_auxvars(idof,ghosted_id)%alpha = work_loc_p(ghosted_id)
          enddo
        enddo
        call VecRestoreArrayReadF90(field%work_loc, &
                                    work_loc_p,ierr);CHKERRQ(ierr)
      class default
        option%io_buffer = 'Unsupported dataset type for BRAGFLO ALPHA.'
        call PrintErrMsg(option)
    end select
  else
    if (.not.wippflo_default_alpha .and. wippflo_use_lumped_harm_flux) then
      option%io_buffer = 'ALPHA should have been read from a dataset.'
      call PrintErrMsg(option)
    endif
  endif

  ! read in elevations
  wippflo_auxvars => patch%aux%WIPPFlo%auxvars
  if (len_trim(this%elevation_dataset_name) > 0) then
    string = 'BRAGFLO Elevation Dataset'
    dataset => DatasetBaseGetPointer(this%realization%datasets, &
                                     this%elevation_dataset_name, &
                                     string,option)
    select type(d => dataset)
      class is(dataset_common_hdf5_type)
        string2 = ''
        string = d%hdf5_dataset_name
        call HDF5ReadCellIndexedRealArray(this%realization, &
                                          field%work, &
                                          d%filename,&
                                          string2, &
                                          string, &
                                          d%realization_dependent)
        wippflo_auxvars => patch%aux%WIPPFlo%auxvars
        call this%realization%comm1%GlobalToLocal(field%work,field%work_loc)
        call VecGetArrayReadF90(field%work_loc,work_loc_p,ierr);CHKERRQ(ierr)
        do ghosted_id = 1, patch%grid%ngmax
          do idof = 0, option%nflowdof
            wippflo_auxvars(idof,ghosted_id)%elevation = work_loc_p(ghosted_id)
          enddo
          !geh: remove after 9/30/19
          !print *, ghosted_id, wippflo_auxvars(0,ghosted_id)%elevation
        enddo
        call VecRestoreArrayReadF90(field%work_loc, &
                                    work_loc_p,ierr);CHKERRQ(ierr)
      class default
        option%io_buffer = 'Unsupported dataset type for WIPP FLOW &
          &Elevation.'
        call PrintErrMsg(option)
    end select
  else if (Initialized(this%rotation_angle)) then
    if (.not.Initialized(this%rotation_origin(3))) then
      option%io_buffer = 'An origin must be defined for dip rotation.'
      call PrintErrMsg(option)
    endif
    call VecGetArrayF90(field%work,work_p,ierr);CHKERRQ(ierr)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      z = grid%z(ghosted_id)
      if (z > this%rotation_ceiling .or. z < this%rotation_basement) then
        work_p(local_id) = z
      else
        x = grid%x(ghosted_id)
        ! from PA.33 in CRA2014 appendix
        h = (x-this%rotation_origin(1))* &
            sin(this%rotation_angle) + &
            (z-this%rotation_origin(3))* &
            cos(this%rotation_angle)
        work_p(local_id) = h + this%rotation_origin(3)
      endif
    enddo
    if (associated(this%rotation_region_names)) then
      do iregion = 1, size(this%rotation_region_names)
        region => RegionGetPtrFromList(this%rotation_region_names(iregion), &
                                       patch%region_list)
        if (.not.associated(region)) then
          option%io_buffer = 'Region "' // &
               trim(this%rotation_region_names(iregion)) // &
               '" in WIPP FLOW Elevation definition not found in region list'
          call PrintErrMsg(option)
        endif
        do icell = 1, region%num_cells
          local_id = region%cell_ids(icell)
          ghosted_id = grid%nL2G(local_id)
          x = grid%x(ghosted_id)
          z = grid%z(ghosted_id)
          ! from PA.33 in CRA2014 appendix
          h = (x-this%rotation_origin(1))* &
              sin(this%rotation_angle) + &
              (z-this%rotation_origin(3))* &
              cos(this%rotation_angle)
          work_p(local_id) = h + this%rotation_origin(3)
        enddo
      enddo
    endif
    call VecRestoreArrayF90(field%work,work_p,ierr);CHKERRQ(ierr)
    call this%realization%comm1%GlobalToLocal(field%work,field%work_loc)
    call VecGetArrayReadF90(field%work_loc,work_loc_p,ierr);CHKERRQ(ierr)
    do ghosted_id = 1, grid%ngmax 
      do idof = 0, option%nflowdof
        wippflo_auxvars(idof,ghosted_id)%elevation = work_loc_p(ghosted_id)
      enddo
      !geh: remove after 9/30/19
      !print *, ghosted_id, wippflo_auxvars(0,ghosted_id)%elevation
    enddo
    call VecRestoreArrayReadF90(field%work_loc,work_loc_p,ierr);CHKERRQ(ierr)
  else ! or set them baesd on grid cell elevation
    do ghosted_id = 1, grid%ngmax
      do idof = 0, option%nflowdof
        wippflo_auxvars(idof,ghosted_id)%elevation = grid%z(ghosted_id)
      enddo
    enddo
  endif
  ! auto pressures by material id
  if (associated(this%auto_pressure_material_ids)) then
    if (Uninitialized(this%rotation_origin(3)) .or. &
        Uninitialized(this%rotation_angle) .or. &
        Uninitialized(this%rotation_ceiling) .or. &
        Uninitialized(this%rotation_basement) .or. &
        Uninitialized(this%auto_pressure_Pb_0)) then
      option%io_buffer = 'DIP_ROTATION_ORIGIN, DIP_ROTATION_CELINING, &
         & DIP_ROTATION_BASEMENT, AUTO_PRESSURE_PB_0 and DIP_ROTATION_ANGLE &
        & must be initialized under WIPP_FLOW,OPTIONS.'
      call PrintErrMsg(option)
   endif
   if (.not.Initialized(this%auto_press_shallow_origin(3))) then
      this%auto_press_shallow_origin = this%rotation_origin
   endif
    ndof = option%nflowdof
    rhob0 = this%auto_pressure_rho_b0
    cb = this%auto_pressure_c_b
    Pbref = this%auto_pressure_Pb_ref
    Pb0 = this%auto_pressure_Pb_0
    zref = this%rotation_origin(3)
    rhobref = rhob0*exp(-cb*(Pbref-Pb0))                         ! PA.56
    Phiref = zref + 1.d0/(gravity*cb)*(1.d0/rhob0-1.d0/rhobref)  ! PA.55
    zref2 = this%auto_press_shallow_origin(3)
    Phiref2 =  zref2 + 0  ! the second term is always zero because Pbref = Pb02
    call VecGetArrayF90(field%flow_xx, flow_xx_p, ierr);CHKERRQ(ierr)
    nmat_id = size(this%auto_pressure_material_ids)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      imat = patch%imat(ghosted_id)
      do i = 1, nmat_id
        if (imat == this%auto_pressure_material_ids(i)) exit
      enddo
      ! if found, i <= nmatid; otherwise i = nmat+1
      if (i <= nmat_id) then
        x = grid%x(ghosted_id)
        z = grid%z(ghosted_id)
        if (z > this%rotation_ceiling) then
          h = z - zref2  ! PA.33 without rotation using shallow origin
          ze = zref2 + h                                                ! PA.57
          rhob = 1.d0/(gravity*cb*(ze-Phiref2+1.d0/(gravity*cb*rhob0))) ! PA.54
        elseif (z < this%rotation_basement) then
          h = z - zref  ! PA.33 without rotation using rotation origin
          ze = zref + h                                                ! PA.57
          rhob = 1.d0/(gravity*cb*(ze-Phiref+1.d0/(gravity*cb*rhob0))) ! PA.54
        else  
          h = (x-this%rotation_origin(1))* &  ! PA.33
              sin(this%rotation_angle) + &
              (z-this%rotation_origin(3))* &
              cos(this%rotation_angle)
          ze = zref + h                                                ! PA.57
          rhob = 1.d0/(gravity*cb*(ze-Phiref+1.d0/(gravity*cb*rhob0))) ! PA.54
        endif
        Pb = Pbref + 1.d0/cb*log(rhob/rhob0)                         ! PA.53
        flow_xx_p((local_id-1)*ndof+WIPPFLO_LIQUID_PRESSURE_DOF) = Pb
        flow_xx_p((local_id-1)*ndof+WIPPFLO_GAS_SATURATION_DOF) = 0.d0
      endif
    enddo
    call VecRestoreArrayF90(field%flow_xx, flow_xx_p, ierr);CHKERRQ(ierr)
    ! have to ensure the auxvars are updated for initial condition output
    call DiscretizationGlobalToLocal(this%realization%discretization, &
                                     field%flow_xx,field%flow_xx_loc,NFLOWDOF)
    call WIPPFloUpdateAuxVars(this%realization)
  endif

  call DiscretizationCreateVector(this%realization%discretization,NFLOWDOF, &
                                  this%scaling_vec,GLOBAL,option)
  call VecDuplicate(this%scaling_vec,this%stored_residual_vec, &
                    ierr);CHKERRQ(ierr)

  ! if using Dirichlet cell-centered BC, set option for matrix
  if (associated(this%dirichlet_dofs_ints)) then
    jcount = 0
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      natural_id = grid%nG2A(ghosted_id)
      do icount = 1, size(this%dirichlet_dofs_ints)
        if (natural_id == this%dirichlet_dofs_ints(1,icount)) then
           if (this%dirichlet_dofs_ints(2,icount) == 1 .or. &
               this%dirichlet_dofs_ints(2,icount) == 3) then
            jcount = jcount + 1
          endif
          if (this%dirichlet_dofs_ints(2,icount) == 2 .or. &
              this%dirichlet_dofs_ints(2,icount) == 3) then
            jcount = jcount + 1
          endif
          exit
        endif   
      enddo 
    enddo
    allocate(this%dirichlet_dofs_ghosted(jcount))
    allocate(this%dirichlet_dofs_local(jcount))

    jcount = 0
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      natural_id = grid%nG2A(ghosted_id)
      do icount = 1, size(this%dirichlet_dofs_ints)
        if (natural_id == this%dirichlet_dofs_ints(1,icount)) then
           if (this%dirichlet_dofs_ints(2,icount) == 1 .or. &
               this%dirichlet_dofs_ints(2,icount) == 3) then
            jcount = jcount + 1
            this%dirichlet_dofs_local(jcount) = (local_id-1)*2+1
            ! zero based indexing
            this%dirichlet_dofs_ghosted(jcount) = (ghosted_id-1)*2
          endif
          if (this%dirichlet_dofs_ints(2,icount) == 2 .or. &
              this%dirichlet_dofs_ints(2,icount) == 3) then
            jcount = jcount + 1
            this%dirichlet_dofs_local(jcount) = (local_id-1)*2+2
            ! zero based indexing
            this%dirichlet_dofs_ghosted(jcount) = (ghosted_id-1)*2+1
          endif
          EXIT
        endif   
      enddo 
    enddo
    call MatSetOption(this%solver%J,MAT_NEW_NONZERO_ALLOCATION_ERR, &
         PETSC_FALSE,ierr);CHKERRQ(ierr)
    deallocate(this%dirichlet_dofs_ints)
  endif

  ! prevent use of block Jacobi preconditioning in parallel
  if (this%solver%pc_type == PCILU .or. &
      this%solver%pc_type == PCBJACOBI) then
    call PetscOptionsHasName(PETSC_NULL_OPTIONS, &
                             PETSC_NULL_CHARACTER,"-bypass_wipp_pc_check", &
                             found,ierr);CHKERRQ(ierr)
    if (.not.found) then
      option%io_buffer = 'Block Jacobi or ILU preconditioning is not allowed &
        &with WIPP_FLOW due to excessive error in the solvers. Please use &
        &a direct solver or FGMRES-CPR.'
      call PrintErrMsg(option)
    endif
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
  use Option_module
  
  implicit none
  
  class(pm_wippflo_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)
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
  
  use Upwind_Direction_module
  use Option_module

  implicit none

  class(pm_wippflo_type) :: this

  PetscInt, save :: lr = 0, gr = 0, lbr = 0, gbr = 0
  PetscInt, save :: lj = 0, gj = 0, lbj = 0, gbj = 0

  if (fix_upwind_direction .and. &
      count_upwind_direction_flip .and. &
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
  use Option_module
  use Realization_Base_class, only : RealizationGetVariable
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL
  use Field_module
  use Global_module, only : GlobalSetAuxVarVecLoc
  use Variables_module, only : LIQUID_SATURATION, GAS_SATURATION
  use WIPP_Flow_Aux_module, only : wippflo_debug_ts_update
  use Utility_module, only : Equal

  implicit none
  
  class(pm_wippflo_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min ! DO NOT USE (see comment below)
  PetscReal :: dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor
  character(len=MAXSTRINGLENGTH) :: string
  
  PetscReal :: dtime(2)
  type(field_type), pointer :: field

  PetscReal :: dt_prev

  dt_prev = dt
  
  ! calculate the time step ramping factor
  dtime(1) = (2.d0*this%gas_sat_change_ts_governor)/ &
             (this%gas_sat_change_ts_governor+this%max_saturation_change)
  dtime(2) = (2.d0*this%liq_pres_change_ts_governor)/ &
             (this%liq_pres_change_ts_governor+this%max_pressure_change)
  ! pick minimum time step from calc'd ramping factor or maximum ramping factor
  dt = min(min(dtime(1),dtime(2))*dt,time_step_max_growth_factor*dt)
  ! make sure time step is within bounds given in the input deck
  dt = min(dt,dt_max)
  if (this%logging_verbosity > 0) then
    if (Equal(dt,dt_max)) then
      string = 'maximum time step size'
    else if (minval(dtime) > time_step_max_growth_factor) then
      string = 'maximum time step growth factor'
    else if (dtime(1) < dtime(2)) then
      string = 'gas saturation governor'
    else
      string = 'liquid pressure governor'
    endif
    string = 'TS update: ' // trim(string)
    call OptionPrint(string,this%option)
  endif
  ! do not use the PFLOTRAN dt_min as it will shut down the simulation from
  ! within timestepper_BE. use %minimum_timestep_size, which is specific to 
  ! wipp_flow.
  dt = max(dt,this%minimum_timestep_size)

  if (wippflo_debug_ts_update) then
    if (minval(dtime(:)) < time_step_max_growth_factor .and. dt < dt_max) then
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
    call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt,dt_max)
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
  use Grid_module

  implicit none
  
  class(pm_wippflo_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr

  PetscViewer :: viewer
  character(len=MAXSTRINGLENGTH) :: string
  type(grid_type), pointer :: grid
  PetscReal, pointer :: r_p(:)
  PetscInt :: i, idof
  
  grid => this%realization%patch%grid

  call PMSubsurfaceFlowUpdatePropertiesNI(this)

  ! calculate residual
  call WIPPFloResidual(snes,xx,r,this%realization,this%pmwss_ptr,ierr)

  ! cell-centered dirichlet BCs
  if (associated(this%dirichlet_dofs_local)) then
    call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
    do i = 1, size(this%dirichlet_dofs_local)
      r_p(this%dirichlet_dofs_local(i)) = 0.d0
    enddo
    call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  endif
  call VecCopy(r,this%stored_residual_vec,ierr);CHKERRQ(ierr)

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
  Vec :: residual_vec
  Vec :: diagonal_vec
  PetscReal :: norm
  PetscInt :: matsize
  PetscInt :: i, irow
  PetscReal, allocatable :: diagonal_values(:)
  PetscReal :: array(1,1)
  PetscReal, pointer :: vec_p(:)


  call WIPPFloJacobian(snes,xx,A,B,this%realization,this%pmwss_ptr,ierr)

  ! cell-centered dirichlet BCs
  if (associated(this%dirichlet_dofs_ghosted)) then
    allocate(diagonal_values(size(this%dirichlet_dofs_local)))
    call VecDuplicate(this%stored_residual_vec,diagonal_vec,ierr);CHKERRQ(ierr)
    call MatGetDiagonal(A,diagonal_vec,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(diagonal_vec,vec_p,ierr);CHKERRQ(ierr)
    do i = 1, size(this%dirichlet_dofs_local)
      diagonal_values(i) = vec_p(this%dirichlet_dofs_local(i)) * 1.d8 + 1.d8
    enddo
    call VecRestoreArrayReadF90(diagonal_vec,vec_p,ierr);CHKERRQ(ierr)
    call VecDestroy(diagonal_vec,ierr);CHKERRQ(ierr)
    i = size(this%dirichlet_dofs_ghosted)
    norm = 1.d0
    ! replace all the rows with zero and the diagonals to 1
    ! on the location of dirchlet_dofs_ghosted indices
    call MatZeroRowsLocal(A,i,this%dirichlet_dofs_ghosted,norm, &
                          PETSC_NULL_VEC,PETSC_NULL_VEC, &
                          ierr);CHKERRQ(ierr)
    do i = 1, size(this%dirichlet_dofs_ghosted)
      irow = this%dirichlet_dofs_ghosted(i)
      array(1,1) = diagonal_values(i)
      call MatSetValuesLocal(A,1,irow,1,irow,array,INSERT_VALUES, &
                             ierr);CHKERRQ(ierr)
    enddo
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    deallocate(diagonal_values)
  endif

  if (this%realization%debug%matview_Jacobian) then
    string = 'WFjacobian'
    call DebugCreateViewer(this%realization%debug,string,this%option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  call SNESGetFunction(snes,residual_vec,PETSC_NULL_FUNCTION, &
                       PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
  if (this%scale_linear_system) then
!    if (this%option%mycommsize > 1) then
!      this%option%io_buffer = 'WIPP FLOW matrix scaling not allowed in &
!        &parallel.'
!      call PrintErrMsg(this%option)
!    endif
    call VecGetLocalSize(this%scaling_vec,matsize,ierr);CHKERRQ(ierr)
    call VecSet(this%scaling_vec,1.d0,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(this%scaling_vec,vec_p,ierr);CHKERRQ(ierr)
    do irow = 1, matsize, 2
      vec_p(irow) = this%linear_system_scaling_factor
    enddo
    call VecRestoreArrayF90(this%scaling_vec,vec_p,ierr);CHKERRQ(ierr)
    call MatDiagonalScale(A,PETSC_NULL_VEC,this%scaling_vec,ierr);CHKERRQ(ierr)
    call MatGetRowMaxAbs(A,this%scaling_vec,PETSC_NULL_INTEGER, &
                         ierr);CHKERRQ(ierr)

    call VecReciprocal(this%scaling_vec,ierr);CHKERRQ(ierr)

    call MatDiagonalScale(A,this%scaling_vec,PETSC_NULL_VEC, &
                          ierr);CHKERRQ(ierr)
    call VecPointwiseMult(residual_vec,residual_vec, &
                          this%scaling_vec,ierr);CHKERRQ(ierr)

    if (this%realization%debug%matview_Jacobian) then
      string = 'WFscale_vec'
      call DebugCreateViewer(this%realization%debug,string,this%option,viewer)
      call VecView(this%scaling_vec,viewer,ierr);CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
      string = 'WFjacobian_scaled'
      call DebugCreateViewer(this%realization%debug,string,this%option,viewer)
      call MatView(A,viewer,ierr);CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    endif
    if (this%realization%debug%vecview_residual) then
      string = 'WFresidual_scaled'
      call DebugCreateViewer(this%realization%debug,string,this%option,viewer)
      call VecView(residual_vec,viewer,ierr);CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    endif
  endif

  if (this%realization%debug%norm_Jacobian) then
    call MatNorm(A,NORM_1,norm,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer,'("1 norm: ",es11.4)') norm
    call PrintMsg(this%option)
    call MatNorm(A,NORM_FROBENIUS,norm,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer,'("2 norm: ",es11.4)') norm
    call PrintMsg(this%option)
    call MatNorm(A,NORM_INFINITY,norm,ierr);CHKERRQ(ierr)
    write(this%option%io_buffer,'("inf norm: ",es11.4)') norm
    call PrintMsg(this%option)
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
  
  this%convergence_flags = 0
  this%convergence_reals = 0.d0
  changed = PETSC_FALSE

  if (this%scale_linear_system) then
    changed = PETSC_TRUE
    call VecStrideScale(dX,ZERO_INTEGER,this%linear_system_scaling_factor, &
                        ierr);CHKERRQ(ierr)
  endif
  
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

  ! If these are changed from true, we must add a global reduction on both
  ! variables to ensure that their values match across all processes. Otherwise
  ! PETSc will throw an error in debug mode or ignore the error in optimized.
  dX_changed = PETSC_TRUE
  X1_changed = PETSC_TRUE
  
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
      if (abs_dX_over_absX >= this%max_allow_rel_liq_pres_chang_ni) then
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
      if ((-1.d0*log10(abs_dX)) < this%neg_log10_rel_gas_sat_change_ni) then
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
      if (X1_p(saturation_index) < &
          (-1.d0*this%gas_sat_thresh_force_ts_cut)) then  ! DEPLIMIT(1)
        saturation_outside_limits = X1_p(saturation_index)
      else 
        if (X1_p(saturation_index) < &
            (-1.d0*this%gas_sat_thresh_force_extra_ni)) then  ! SATLIMIT
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
      if (X1_p(saturation_index) > &
          1.d0 + this%gas_sat_thresh_force_ts_cut) then  ! DEPLIMIT(1)
        saturation_outside_limits = X1_p(saturation_index)
      else 
        if (X1_p(saturation_index) > &
            1.d0 + this%gas_sat_thresh_force_extra_ni) then  ! SATLIMIT
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
    if (X1_p(pressure_index) <= (this%min_liq_pres_force_ts_cut)) then  ! DEPLIMIT(2)
      pressure_outside_limits = X1_p(pressure_index)
    endif

    if (dabs(pressure_outside_limits) > 0.d0 .or. &
        dabs(saturation_outside_limits) > 0.d0) then
      cut_timestep = PETSC_TRUE
      write(*,'(4x,"Outside Limits (PL,SG): ",i8,2es10.2)') &
        grid%nG2A(ghosted_id), &
        pressure_outside_limits, saturation_outside_limits
    endif

  enddo

  if (wippflo_debug_first_iteration) stop

  ! the following flags are used in detemining convergence
  if (.not.converged_liquid_pressure) then 
    this%convergence_flags(MAX_REL_CHANGE_LIQ_PRES_NI) = &
      max_liq_pres_rel_change_cell
  endif
  if (.not.converged_gas_saturation) then 
    this%convergence_flags(MAX_CHANGE_GAS_SAT_NI) = max_gas_sat_change_NI_cell
  endif
  if (max_abs_pressure_change_TS > this%max_allow_liq_pres_change_ts) then
    this%convergence_flags(MAX_CHANGE_LIQ_PRES_TS) = &
                                             max_abs_pressure_change_TS_cell
  endif
  if (max_gas_sat_change_TS > this%max_allow_gas_sat_change_ts) then
    this%convergence_flags(MAX_CHANGE_GAS_SAT_TS) = max_gas_sat_change_TS_cell
  endif
  if (force_another_iteration) this%convergence_flags(FORCE_ITERATION) = &
    max_gas_sat_outside_lim_cell
  if (cut_timestep) this%convergence_flags(OUTSIDE_BOUNDS) = 1

  ! the following flags are for REPORTING purposes only
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

  call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(X1,X1_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%max_change_vecs(1),press_ptr, &
                              ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%max_change_vecs(3),sat_ptr, &
                              ierr);CHKERRQ(ierr)
                               
end subroutine PMWIPPFloCheckUpdatePost

! ************************************************************************** !

subroutine PMWIPPFloCheckConvergence(this,snes,it,xnorm,unorm, &
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

  ! check residual terms
  if (this%stored_residual_vec == PETSC_NULL_VEC) then
    residual_vec = field%flow_r
  else
    ! this vector is to be used if linear system scaling is employed when 
    ! the residual is altered by scaling.  This vector stores the original 
    ! residual. 
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
      if (dabs(residual) > this%liquid_residual_infinity_tol) then
        if (abs_residual_over_accumulation > &
            this%liquid_residual_infinity_tol) then
          converged_liquid_equation = PETSC_FALSE
          if (dabs(max_normal_res_liq_) < dabs(residual)) then
            max_normal_res_liq_cell = local_id
          endif
        endif
      endif
      ! update outside to always record the maximum residual
      if (dabs(max_normal_res_liq_) < dabs(residual)) then
        if (wippflo_match_bragflo_output) then
          if (dabs(residual) > this%liquid_residual_infinity_tol) then
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
      if (dabs(residual) > this%gas_equation_infinity_tol) then 
        if (abs_residual_over_accumulation > &
            this%gas_equation_infinity_tol) then
          converged_gas_equation = PETSC_FALSE
          if (dabs(max_normal_res_gas_) < dabs(residual)) then
            max_normal_res_gas_cell = local_id
          endif
        endif
      endif
      ! update outside to always record the maximum residual
      if (dabs(max_normal_res_gas_) < dabs(residual)) then
        if (wippflo_match_bragflo_output) then
          if (dabs(residual) > this%gas_equation_infinity_tol) then
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

  ! the following flags are used in detemining convergence
  if (.not.converged_liquid_equation) then
    this%convergence_flags(MAX_NORMAL_RES_LIQ) = max_normal_res_liq_cell
  endif
  if (.not.converged_gas_equation) then
    this%convergence_flags(MAX_NORMAL_RES_GAS) = max_normal_res_gas_cell
  endif
  if (min_gas_pressure < 0.d0 .and. .not.wippflo_allow_neg_gas_pressure) then
    this%convergence_flags(MIN_GAS_PRES) = min_gas_pressure_cell
  endif
  ! the following flags are for REPORTING purposes only
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
    if (this%convergence_test_both .or. .not.converged_liquid_equation) then
      if (converged_flag /= CONVERGENCE_CUT_TIMESTEP) then
        converged_flag = CONVERGENCE_KEEP_ITERATING ! cannot override cut
      endif
      reason_string(6:6) = 'P'
    endif
  endif
  if (this%convergence_flags(MAX_CHANGE_GAS_SAT_NI) > 0) then
    if (converged_flag /= CONVERGENCE_CUT_TIMESTEP) then
      if (this%convergence_test_both .or. .not.converged_gas_equation) then
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
    if (this%convergence_flags(MAX_CHANGE_LIQ_PRES_TS) > 0) then
      tempreal3 = this%convergence_reals(MAX_CHANGE_LIQ_PRES_TS)
      tempint3 = this%convergence_flags(MAX_CHANGE_LIQ_PRES_TS)
      reason_string(6:6) = 'p'
    else
      tempreal3 = this%convergence_reals(MAX_REL_CHANGE_LIQ_PRES_NI)
      tempint3 = this%convergence_flags(MAX_REL_CHANGE_LIQ_PRES_NI)
    endif
    if (this%convergence_flags(MAX_CHANGE_GAS_SAT_TS) > 0) then
      tempreal4 = this%convergence_reals(MAX_CHANGE_GAS_SAT_TS)
      tempint4 = this%convergence_flags(MAX_CHANGE_GAS_SAT_TS)
    else
      tempreal4 = this%convergence_reals(MAX_CHANGE_GAS_SAT_NI)
      tempint4 = this%convergence_flags(MAX_CHANGE_GAS_SAT_NI)
    endif
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
    if (option%mycommsize > 1 .or. grid%nmax > 9999) then
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

  call PMSubsurfaceFlowCheckConvergence(this,snes,it,xnorm,unorm,fnorm, &
                                        reason,ierr)

end subroutine PMWIPPFloCheckConvergence

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
      dabs(this%convergence_reals(MAX_CHANGE_GAS_SAT_NI)/ &
           this%max_allow_rel_gas_sat_change_ni), &
      this%convergence_flags(MAX_REL_CHANGE_LIQ_PRES_NI), &
      dabs(this%convergence_reals(MAX_REL_CHANGE_LIQ_PRES_NI)/ &
           this%max_allow_rel_liq_pres_chang_ni), &
      this%convergence_flags(MAX_NORMAL_RES_GAS), &
      dabs(this%convergence_reals(MAX_NORMAL_RES_GAS)/ &
           this%gas_equation_infinity_tol), &
      this%convergence_flags(MAX_NORMAL_RES_LIQ), &
      dabs(this%convergence_reals(MAX_NORMAL_RES_LIQ)/ &
           this%liquid_residual_infinity_tol)
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
          if (dabs(vec_new_ptr(j)) > this%gas_sat_gov_switch_abs_to_rel) then
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
  use Utility_module, only : DeallocateArray

  implicit none
  
  class(pm_wippflo_type) :: this

  PetscErrorCode :: ierr
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  call DeallocateArray(this%max_change_ivar)
  call DeallocateArray(this%rotation_region_names)
  call DeallocateArray(this%auto_pressure_material_ids)
  if (this%stored_residual_vec /= PETSC_NULL_VEC) then
    call VecDestroy(this%stored_residual_vec,ierr);CHKERRQ(ierr)
  endif
  if (this%scaling_vec /= PETSC_NULL_VEC) then
    call VecDestroy(this%scaling_vec,ierr);CHKERRQ(ierr)
  endif
  if (associated(this%dirichlet_dofs_ghosted)) then
    deallocate(this%dirichlet_dofs_ghosted)
    nullify(this%dirichlet_dofs_ghosted)
  endif
  if (associated(this%dirichlet_dofs_local)) then
    deallocate(this%dirichlet_dofs_local)
    nullify(this%dirichlet_dofs_local)
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
