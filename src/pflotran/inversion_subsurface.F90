module Inversion_Subsurface_class

#include "petsc/finclude/petscksp.h"
  use petscksp

  use PFLOTRAN_Constants_module
  use Inversion_Aux_module
  use Inversion_Coupled_Aux_module
  use Inversion_TS_Aux_module
  use Inversion_Measurement_Aux_module
  use Inversion_Parameter_module
  use Inversion_Base_class
  use Realization_Subsurface_class
  use Simulation_Subsurface_class
  use Option_Inversion_module

  implicit none

  private

  PetscInt, parameter, public :: GET_MATERIAL_VALUE = 0
  PetscInt, parameter :: OVERWRITE_MATERIAL_VALUE = 1
  PetscInt, parameter :: COPY_TO_VEC = 3
  PetscInt, parameter :: COPY_FROM_VEC = 4

  PetscInt, parameter, public :: INVSUBSCATFORWARD = 0
  PetscInt, parameter, public :: INVSUBSCATREVERSE = 1

  type, public, extends(inversion_base_type) :: inversion_subsurface_type
    character(len=MAXSTRINGLENGTH) :: forward_simulation_filename
    class(simulation_subsurface_type), pointer :: forward_simulation
    class(realization_subsurface_type), pointer :: realization
    type(inversion_aux_type), pointer :: inversion_aux
    type(inversion_measurement_aux_type), pointer :: measurements(:)
    type(inversion_coupled_aux_type), pointer :: inversion_coupled_aux
    type(inversion_parameter_type), pointer :: parameters(:)
    type(perturbation_type), pointer :: perturbation
    PetscInt :: dist_measurement_offset
    PetscInt :: parameter_offset  ! needed?
    PetscInt :: num_parameters_local
    PetscInt :: n_qoi_per_cell
    Vec :: quantity_of_interest       ! reserved for inversion_ert
    Vec :: ref_quantity_of_interest   ! reserved for inversion_ert
    Vec :: measurement_vec
    Vec :: dist_measurement_vec
    Vec :: parameter_vec
    Vec :: dist_parameter_vec
    Vec :: ref_parameter_vec    ! needed?
    character(len=MAXWORDLENGTH) :: ref_qoi_dataset_name
    ! flags reference in objects below the inversion objects should be
    ! stored in option_inversion_type
    PetscBool :: print_sensitivity_jacobian
    PetscBool :: qoi_is_full_vector
    PetscBool :: first_inversion_interation
    PetscBool :: annotate_output
    PetscBool :: perturbation_risk_acknowledged
    PetscBool :: debug_adjoint
    PetscInt :: debug_verbosity
    VecScatter :: scatter_measure_to_dist_measure
    VecScatter :: scatter_param_to_dist_param
    VecScatter :: scatter_global_to_dist_param
    PetscReal, pointer :: local_measurement_values(:)
    PetscReal, pointer :: local_dobs_dunknown_values(:)
    PetscReal, pointer :: local_dobs_dparam_values(:)
    PetscInt, pointer :: local_measurement_map(:)
  contains
    procedure, public :: Init => InversionSubsurfaceInit
    procedure, public :: ReadBlock => InversionSubsurfReadBlock
    procedure, public :: InitializeForwardRun => InvSubsurfInitForwardRun
    procedure, public :: SetupForwardRunLinkage => &
                           InvSubsurfSetupForwardRunLinkage
    procedure, public :: ConnectToForwardRun => InvSubsurfConnectToForwardRun
    procedure, public :: ExecuteForwardRun => InvSubsurfExecuteForwardRun
    procedure, public :: DestroyForwardRun => InvSubsurfDestroyForwardRun
    procedure, public :: EvaluateCostFunction => InvSuburfSkipThisOnly
    procedure, public :: CheckConvergence => InvSuburfSkipThisOnly
    procedure, public :: WriteIterationInfo => InvSuburfSkipThisOnly
    procedure, public :: CalculateSensitivity => InvSubsurfCalculateSensitivity
    procedure, public :: ScaleSensitivity => InvSuburfSkipThisOnly
    procedure, public :: CalculateUpdate => InvSuburfSkipThisOnly
    procedure, public :: UpdateRegularizationParameters => &
                           InvSuburfSkipThisOnly
    procedure, public :: Strip => InversionSubsurfaceStrip
  end type inversion_subsurface_type

  type perturbation_type
    Vec :: base_parameter_vec
    Vec :: base_measurement_vec
    PetscInt :: ndof
    PetscInt :: idof_pert
    PetscReal :: pert
    PetscReal :: base_value
    PetscReal :: tolerance
    PetscInt, pointer :: select_cells(:)
  end type perturbation_type

  public :: InversionSubsurfaceCreate, &
            InversionSubsurfaceInit, &
            InversionSubsurfReadSelectCase, &
            InvSubsurfSetupForwardRunLinkage, &
            InvSubsurfConnectToForwardRun, &
            InvSubsurfOutputSensitivity, &
            InversionSubsurfaceStrip

  public :: InvSubsurfScatMeasToDistMeas, &
            InvSubsurfScatParamToDistParam, &
            InvSubsurfScatGlobalToDistParam, &
            InvSubsurfGetParamValueByCell, &
            InvSubsurfGetSetParamValueByMat

contains

! ************************************************************************** !

function InversionSubsurfaceCreate(driver)
  !
  ! Allocates and initializes a new subsurface inversion object
  !
  ! Author: Glenn Hammond
  ! Date: 10/25/21
  !
  use Driver_class

  class(driver_type), pointer :: driver

  class(inversion_subsurface_type), pointer :: InversionSubsurfaceCreate

  allocate(InversionSubsurfaceCreate)
  call InversionSubsurfaceCreate%Init(driver)

end function InversionSubsurfaceCreate

! ************************************************************************** !

subroutine InversionSubsurfaceInit(this,driver)
  !
  ! Initializes a new inversion object
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
  !
  use Driver_class

  class(inversion_subsurface_type) :: this
  class(driver_type), pointer :: driver

  call InversionBaseInit(this,driver)

  this%quantity_of_interest = PETSC_NULL_VEC
  this%n_qoi_per_cell = UNINITIALIZED_INTEGER
  this%measurement_vec = PETSC_NULL_VEC
  this%dist_measurement_vec = PETSC_NULL_VEC
  this%parameter_vec = PETSC_NULL_VEC
  this%dist_parameter_vec = PETSC_NULL_VEC
  this%ref_parameter_vec = PETSC_NULL_VEC
  this%ref_quantity_of_interest = PETSC_NULL_VEC
  this%ref_qoi_dataset_name = ''
  this%forward_simulation_filename = ''
  this%print_sensitivity_jacobian = PETSC_FALSE
  this%debug_adjoint = PETSC_FALSE
  this%debug_verbosity = UNINITIALIZED_INTEGER
  this%scatter_measure_to_dist_measure = PETSC_NULL_VECSCATTER
  this%scatter_param_to_dist_param = PETSC_NULL_VECSCATTER
  this%scatter_global_to_dist_param = PETSC_NULL_VECSCATTER
  this%dist_measurement_offset = UNINITIALIZED_INTEGER
  this%parameter_offset = UNINITIALIZED_INTEGER
  this%num_parameters_local = UNINITIALIZED_INTEGER
  this%qoi_is_full_vector = PETSC_FALSE
  this%first_inversion_interation = PETSC_TRUE
  this%annotate_output = PETSC_FALSE
  this%perturbation_risk_acknowledged = PETSC_FALSE

  nullify(this%local_measurement_values)
  nullify(this%local_dobs_dunknown_values)
  nullify(this%local_dobs_dparam_values)
  nullify(this%local_measurement_map)

  nullify(this%measurements)
  nullify(this%parameters)
  nullify(this%perturbation)
  nullify(this%inversion_coupled_aux)

  nullify(this%forward_simulation)
  nullify(this%realization)
  nullify(this%inversion_aux)

  ! initialize measurement reporting verbosity
  inv_meas_reporting_verbosity = 1

end subroutine InversionSubsurfaceInit

! ************************************************************************** !

function InvSuburfPerturbationCreate()
  !
  ! Allocates and initializes a new perturbation object
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  type(perturbation_type), pointer :: InvSuburfPerturbationCreate

  allocate(InvSuburfPerturbationCreate)
  InvSuburfPerturbationCreate%base_parameter_vec = PETSC_NULL_VEC
  InvSuburfPerturbationCreate%base_measurement_vec = PETSC_NULL_VEC

  InvSuburfPerturbationCreate%ndof = 0
  InvSuburfPerturbationCreate%idof_pert = 0
  InvSuburfPerturbationCreate%pert = 0.d0
  InvSuburfPerturbationCreate%base_value = 0.d0
  InvSuburfPerturbationCreate%tolerance = 1.d-6
  nullify(InvSuburfPerturbationCreate%select_cells)

end function InvSuburfPerturbationCreate

! ************************************************************************** !

subroutine InvSuburfSkipThisOnly(this)
  !
  ! A dummy subroutine that is skipped
  !
  ! Author: Glenn Hammond
  ! Date: 03/25/22
  !
  class(inversion_subsurface_type) :: this

  if (this%maximum_iteration /= -1) then
    call this%driver%PrintErrMsg('All inversion routines must be extended, &
      &even if by implementing a "skip" routines: InvSuburfSkipThisOnly')
  endif

end subroutine InvSuburfSkipThisOnly

! ************************************************************************** !

subroutine InversionSubsurfReadBlock(this,input,option)

  use Input_Aux_module
  use Option_module
  use String_module

  class(inversion_subsurface_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  error_string = 'Test Sensitivity Jacobian'

  ! just testing the sensitivity Jacobian
  this%maximum_iteration = -1

  input%ierr = 0
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    found = PETSC_TRUE
    call InversionSubsurfReadSelectCase(this,input,keyword,found, &
                                        error_string,option)
    if (found) cycle

    select case(trim(keyword))
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select

  enddo
  call InputPopBlock(input,option)

end subroutine InversionSubsurfReadBlock

! ************************************************************************** !

subroutine InversionSubsurfReadSelectCase(this,input,keyword,found, &
                                          error_string,option)

  use Input_Aux_module
  use Option_module
  use String_module
  use Units_module
  use Utility_module

  class(inversion_subsurface_type) :: this
  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: string
  type(inversion_measurement_aux_type), pointer :: new_measurement
  type(inversion_measurement_aux_type), pointer :: first_measurement
  type(inversion_measurement_aux_type), pointer :: last_measurement
  type(inversion_parameter_type), pointer :: new_parameter
  type(inversion_parameter_type), pointer :: first_parameter
  type(inversion_parameter_type), pointer :: last_parameter
  PetscInt :: i
  PetscInt :: observed_variable
  PetscReal :: measurement_time
  character(len=4) :: measurement_time_units
  character(len=MAXWORDLENGTH) :: internal_units

  nullify(new_measurement)
  nullify(last_measurement)
  nullify(new_parameter)
  nullify(last_parameter)

  found = PETSC_TRUE
  call InversionBaseReadSelectCase(this,input,keyword,found, &
                                   error_string,option)
  if (found) return

  found = PETSC_TRUE
  select case(trim(keyword))
    case('FORWARD_SIMULATION_FILENAME')
      call InputReadFilename(input,option,this%forward_simulation_filename)
      call InputErrorMsg(input,option,keyword,error_string)
    case('OBSERVATION_FUNCTION')
      option%io_buffer = 'OBSERVATION_FUNCTION (renamed OBSERVED_VARIABLE) &
        &must be specified in the MEASUREMENT block.'
      call PrintErrMsg(option)
    case('MEASUREMENTS')
      observed_variable = UNINITIALIZED_INTEGER
      measurement_time = UNINITIALIZED_DOUBLE
      measurement_time_units = ''
      string = trim(error_string)//keyword
      input%ierr = 0
      call InputPushBlock(input,option)
      do
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadCard(input,option,keyword)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(keyword)
        select case(trim(keyword))
          case('MEASUREMENT')
            new_measurement => InversionMeasurementAuxRead(input,string,option)
            ! set the observed variable if not specified in
            ! measurement subblock
            if (Uninitialized(new_measurement%iobs_var)) then
              if (Initialized(observed_variable)) then
                new_measurement%iobs_var = observed_variable
              else
                option%io_buffer = 'The OBSERVED_VARIABLE is not specified &
                  &for a MEASUREMENT.'
                call PrintErrMsg(option)
              endif
            endif
            if (Uninitialized(new_measurement%time)) then
              if (Initialized(measurement_time)) then
                new_measurement%time = measurement_time
                new_measurement%time_units = measurement_time_units
              endif
            endif
            if (associated(last_measurement)) then
              last_measurement%next => new_measurement
              new_measurement%id = last_measurement%id + 1
            else
              first_measurement => new_measurement
              first_measurement%id = 1
            endif
            last_measurement => new_measurement
            nullify(new_measurement)
          case('TIME')
            call InputReadDouble(input,option,measurement_time)
            call InputErrorMsg(input,option,keyword,error_string)
            call InputReadWord(input,option,word,PETSC_TRUE)
            if (input%ierr /= 0) word = 'sec'
            measurement_time_units = trim(word)
            internal_units = 'sec'
            measurement_time = measurement_time* &
                               UnitsConvertToInternal(word,internal_units,option)
          case('OBSERVED_VARIABLE')
            observed_variable = &
              InvMeasAuxReadObservedVariable(input,keyword,error_string,option)
          case('REPORTING_VERBOSITY')
            call InputReadInt(input,option,inv_meas_reporting_verbosity)
            call InputErrorMsg(input,option,keyword,error_string)
          case default
            call InputKeywordUnrecognized(input,keyword,error_string,option)
        end select
      enddo
      call InputPopBlock(input,option)
      if (.not.associated(last_measurement)) then
        option%io_buffer = 'No measurement found in inversion measurement block.'
        call PrintErrMsg(option)
      else if (associated(this%measurements)) then
        option%io_buffer = 'Measurements may only be defined in a single block.'
        call PrintErrMsg(option)
      else
        allocate(this%measurements(last_measurement%id))
        do i = 1, last_measurement%id
          call InversionMeasurementAuxInit(this%measurements(i))
        enddo
        last_measurement => first_measurement
        do
          if (.not.associated(last_measurement)) exit
          call InversionMeasurementAuxCopy(last_measurement, &
                                         this%measurements(last_measurement%id))
          last_measurement => last_measurement%next
        enddo
        call InversionMeasureAuxListDestroy(first_measurement)
      endif
    case('PARAMETERS')
      string = trim(error_string)//keyword
      input%ierr = 0
      call InputPushBlock(input,option)
      do
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadCard(input,option,keyword)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(keyword)
        select case(trim(keyword))
          case('PARAMETER')
            new_parameter => InversionParameterRead(input,string,option)
          case default
            call InputKeywordUnrecognized(input,keyword,error_string,option)
        end select
        if (associated(last_parameter)) then
          last_parameter%next => new_parameter
          new_parameter%id = last_parameter%id + 1
        else
          first_parameter => new_parameter
          first_parameter%id = 1
        endif
        last_parameter => new_parameter
        nullify(new_parameter)
      enddo
      call InputPopBlock(input,option)
      if (.not.associated(last_parameter)) then
        option%io_buffer = 'No parameter found in inversion parameter block.'
        call PrintErrMsg(option)
      else if (associated(this%parameters)) then
        option%io_buffer = 'Parameters may only be defined in a single block.'
        call PrintErrMsg(option)
      else
        allocate(this%parameters(last_parameter%id))
        do i = 1, last_parameter%id
          call InversionParameterInit(this%parameters(i))
        enddo
        last_parameter => first_parameter
        do
          if (.not.associated(last_parameter)) exit
          call InversionParameterCopy(last_parameter, &
                                           this%parameters(last_parameter%id))
          last_parameter => last_parameter%next
        enddo
        call InversionParameterDestroy(first_parameter)
      endif
    case('PRINT_SENSITIVITY_JACOBIAN')
      this%print_sensitivity_jacobian = PETSC_TRUE
    case('DEBUG_ADJOINT')
      this%debug_adjoint = PETSC_TRUE
      call InputReadInt(input,option,i)
      if (input%ierr == 0) then
        this%debug_verbosity = i
      endif
    case('COUPLED_FLOW_AND_ERT')
      this%inversion_option%coupled_flow_ert = PETSC_TRUE
    case('PERTURBATION')
      string = trim(error_string)//keyword
      input%ierr = 0
      this%inversion_option%use_perturbation = PETSC_TRUE
      this%perturbation => InvSuburfPerturbationCreate()
      call InputPushBlock(input,option)
      do
        call InputReadPflotranString(input,option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadCard(input,option,keyword)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(keyword)
        select case(trim(keyword))
          case('PERTURBATION_TOLERANCE')
              call InputReadDouble(input,option,this%perturbation%tolerance)
              call InputErrorMsg(input,option,keyword,error_string)
          case('SELECT_CELLS')
            call UtilityReadArray(this%perturbation%select_cells, &
                                  ZERO_INTEGER,error_string,input,option)
          case('ACKNOWLEDGE_RISK_OF_PERTURBATION')
            this%perturbation_risk_acknowledged = PETSC_TRUE
          case('ANNOTATE_PERTURBATION_OUTPUT')
            this%annotate_output = PETSC_TRUE
          case default
            call InputKeywordUnrecognized(input,keyword,error_string,option)
        end select
      enddo
      call InputPopBlock(input,option)
    case default
      found = PETSC_FALSE
  end select

end subroutine InversionSubsurfReadSelectCase

! ************************************************************************** !

subroutine InvSubsurfInitForwardRun(this,option)
  !
  ! Initializes the forward simulation
  !
  ! Author: Glenn Hammond
  ! Date: 03/21/22
  !
  use Factory_Forward_module
  use Option_module
  use String_module

  class(inversion_subsurface_type) :: this
  type(option_type), pointer :: option

  option => OptionCreate()
  write(option%group_prefix,'(i6)') this%iteration
  option%group_prefix = 'Run' // trim(adjustl(option%group_prefix))
  if (associated(this%perturbation)) then
    if (this%annotate_output) then
      if (this%perturbation%idof_pert > 0) then
        option%group_prefix = trim(option%group_prefix) // 'P' // &
          StringWrite(this%perturbation%idof_pert)
      else if (this%perturbation%idof_pert == 0) then
        option%group_prefix = trim(option%group_prefix) // 'Base'
      else
        option%group_prefix = trim(option%group_prefix) // 'Final'
      endif
    endif
  endif
  call OptionSetDriver(option,this%driver)
  call OptionSetInversionOption(option,this%inversion_option)
  call FactoryForwardInitialize(this%forward_simulation, &
                                this%forward_simulation_filename,option)
  this%realization => this%forward_simulation%realization

end subroutine InvSubsurfInitForwardRun

! ************************************************************************** !

subroutine InvSubsurfSetupForwardRunLinkage(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  use Connection_module
  use Coupler_module
  use Discretization_module
  use Grid_module
  use Material_module
  use Option_module
  use Patch_module
  use PM_Base_class
  use PM_ERT_class
  use PM_Subsurface_Flow_class
  use String_module
  use Variables_module, only : PERMEABILITY
  use Waypoint_module
  use ZFlow_Aux_module

  class(inversion_subsurface_type) :: this

  type(patch_type), pointer :: patch
  type(inversion_forward_aux_type), pointer :: inversion_forward_aux
  type(material_property_type), pointer :: material_property
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: iflag
  PetscInt :: i
  PetscInt :: local_id
  PetscInt :: temp_int
  PetscInt :: icount
  PetscInt :: num_measurements, num_measurements_local
  PetscInt :: num_parameters
  PetscInt, allocatable :: int_array(:), int_array2(:)
  PetscReal :: tempreal
  PetscReal, pointer :: vec_ptr(:)
  Vec :: v, v2
  IS :: is_petsc
  IS :: is_measure
  PetscMPIInt :: mpi_int
  PetscInt :: max_int(2)
  IS :: is_parameter
  PetscErrorCode :: ierr

  nullify(vec_ptr)

  patch => this%realization%patch
  if (.not.associated(this%inversion_aux)) then
    ! perturbation can be problematic with certain flow/transport configurations
    ! check for these situations here
    if (associated(this%perturbation) .and. &
        .not.this%perturbation_risk_acknowledged) then
      select type(pm => &
                  this%forward_simulation%flow_process_model_coupler%pm_list)
        class is(pm_subsurface_flow_type)
          if (Initialized(pm%cfl_governor) .or. &
              Initialized(pm%pressure_change_governor) .or. &
              Initialized(pm%temperature_change_governor) .or. &
              Initialized(pm%saturation_change_governor) .or. &
              Initialized(pm%xmol_change_governor)) then
            call this%driver%PrintErrMsg('The use of timestep size governors &
              &combined with a sensitivity Jacobian calculated using &
              &perturbation can produce incorrect derivatives. Please add &
              &the card ACKNOWLEDGE_RISK_OF_PERTURBATION to the PERTURBATION &
              &block in the input file.')
          endif
          ! add other situations here
        ! add other classes here
      end select
    endif
    this%n_qoi_per_cell = 1 ! 1 perm per cell

    this%inversion_aux => InversionAuxCreate()
    num_measurements = 0
    if (associated(this%measurements)) then
      num_measurements = size(this%measurements)
    else
      call this%driver%PrintErrMsg('No inversion measurements defined.')
    endif
    if (associated(this%parameters)) then
      num_parameters = 0
      do i = 1, size(this%parameters)
        call InversionParameterMapNameToInt(this%parameters(i),this%driver)
        if (len_trim(this%parameters(i)%material_name) > 0) then
          material_property => &
            MaterialPropGetPtrFromArray(this%parameters(i)%material_name, &
                              this%realization%patch%material_property_array)
          if (.not.associated(material_property)) then
            call this%driver%PrintErrMsg('Inversion MATERIAL "' // &
                trim(this%parameters(i)%material_name) // '" not found among &
            &MATERIAL_PROPERTIES.')
          endif
          this%parameters(i)%imat = abs(material_property%internal_id)
          num_parameters = num_parameters + 1
        else
          num_parameters = num_parameters + patch%grid%nmax
        endif
      enddo
      if (.not.associated(this%perturbation)) then
        do i = 1, size(this%parameters)
          if (i == 1) then
            temp_int = this%parameters(i)%iparameter
          else
            if (temp_int /= this%parameters(i)%iparameter) then
              call this%driver%PrintErrMsg('Inversion by multiple &
                &parameters of differing type (e.g. permeability, &
                &porosity) only supported for perturbation.')
            endif
          endif
        enddo
      endif
      if (num_parameters == patch%grid%nmax) then
        this%qoi_is_full_vector = PETSC_TRUE
        this%num_parameters_local = patch%grid%nlmax*this%n_qoi_per_cell
      else
        this%num_parameters_local = PETSC_DECIDE
      endif
    else
      call this%driver%PrintErrMsg('No inversion parameters defined.')
    endif
    if (size(this%parameters) > 1 .and. this%qoi_is_full_vector) then
      call this%driver%PrintErrMsg('More than one parameter not currently &
                                   &supported for full vector inversion.')
    endif

    if (this%inversion_option%coupled_flow_ert) then
      do i = 1, size(this%parameters)
        if (InversionParameterGetIDFromName(this%parameters(i)% &
                                              parameter_name,this%driver) /= &
            PERMEABILITY) then
          call this%driver%PrintErrMsg('COUPLED_ZFLOW_ERT currently only &
                                        &supported for permeability.')
        endif
      enddo
      if (.not.(associated(this%forward_simulation% &
                            flow_process_model_coupler) .and. &
                associated(this%forward_simulation% &
                            geop_process_model_coupler))) then
        call this%driver%PrintErrMsg('Coupled ZFLOW and ERT inversion &
          &requires that both ZFLOW and ERT process models are employed.')
      endif
      this%inversion_coupled_aux => InversionCoupledAuxCreate()
      this%inversion_coupled_aux%parameters => this%parameters
    endif

    ! JsensitivityT is the transpose of the sensitivity Jacobian
    ! with num measurement columns and num parameter rows
    call MatCreateDense(this%driver%comm%mycomm, &
                        this%num_parameters_local,PETSC_DECIDE, &
                        num_parameters,num_measurements, &
                        PETSC_NULL_SCALAR,this%inversion_aux%JsensitivityT, &
                        ierr);CHKERRQ(ierr)
    call MatZeroEntries(this%inversion_aux%JsensitivityT,ierr);CHKERRQ(ierr)
    ! cannot pass in this%measurement_vec as it is initialized to
    ! PETSC_NULL_VEC and MatCreateVecs keys off that input
    call MatCreateVecs(this%inversion_aux%JsensitivityT,v,v2, &
                       ierr);CHKERRQ(ierr)
    this%dist_measurement_vec = v
    this%dist_parameter_vec = v2
    call MatGetLocalSize(this%inversion_aux%JsensitivityT,temp_int, &
                         num_measurements_local,ierr);CHKERRQ(ierr)
    if (this%num_parameters_local == PETSC_DECIDE) then
      this%num_parameters_local = temp_int
    endif
    if (temp_int /= this%num_parameters_local) then
      call this%driver%PrintErrMsg('Misalignment in MatGetLocalSize ('//&
                 trim(StringWrite(temp_int))//','//&
                 trim(StringWrite(this%num_parameters_local))//') in &
                 &InvSubsurfSetupForwardRunLinkage.')
    endif
    allocate(int_array(2),int_array2(2))
    int_array(1) = this%num_parameters_local
    int_array(2) = num_measurements_local
    int_array2(:) = 0
    call MPI_Exscan(int_array,int_array2,TWO_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                    this%driver%comm%mycomm,ierr);CHKERRQ(ierr)
    this%parameter_offset = int_array2(1)
    this%dist_measurement_offset = int_array2(2)
    deallocate(int_array,int_array2)
    call VecCreateSeq(PETSC_COMM_SELF,num_measurements,this%measurement_vec, &
                      ierr);CHKERRQ(ierr)

    if (this%qoi_is_full_vector) then
      ! is_parameter should mirror the natural vec
      allocate(int_array(this%num_parameters_local))
      do i = 1, this%num_parameters_local
        int_array(i) = i
      enddo
      int_array = this%parameter_offset + int_array - 1
      call DiscretAOApplicationToPetsc(this%realization%discretization, &
                                       int_array)
      call ISCreateGeneral(this%driver%comm%mycomm,size(int_array),int_array, &
                           PETSC_COPY_VALUES,is_petsc,ierr);CHKERRQ(ierr)
      deallocate(int_array)
      call VecScatterCreate(this%realization%field%work,is_petsc, &
                            this%dist_parameter_vec,PETSC_NULL_IS, &
                            this%scatter_global_to_dist_param, &
                            ierr);CHKERRQ(ierr)
      call ISDestroy(is_petsc,ierr);CHKERRQ(ierr)
    else
      call VecCreateSeq(PETSC_COMM_SELF,num_parameters,this%parameter_vec, &
                        ierr);CHKERRQ(ierr)
      call ISCreateStride(this%driver%comm%mycomm,num_parameters,ZERO_INTEGER, &
                          ONE_INTEGER,is_parameter,ierr);CHKERRQ(ierr)
      call VecScatterCreate(this%parameter_vec,is_parameter, &
                            this%dist_parameter_vec,PETSC_NULL_IS, &
                            this%scatter_param_to_dist_param, &
                            ierr);CHKERRQ(ierr)
      call ISDestroy(is_parameter,ierr)
    endif

    do i = 1, num_measurements
      iflag = PETSC_FALSE
      ! ensure that all observed variables are being simulated
      select case(this%measurements(i)%iobs_var)
        case(OBS_LIQUID_PRESSURE)
          if (Uninitialized(zflow_liq_flow_eq)) then
            string = 'Liquid pressure'
            iflag = PETSC_TRUE
          endif
        case(OBS_LIQUID_SATURATION)
          if (Uninitialized(zflow_liq_flow_eq)) then
            string = 'Liquid saturation'
            iflag = PETSC_TRUE
          endif
        case(OBS_SOLUTE_CONCENTRATION)
          if (Uninitialized(zflow_sol_tran_eq)) then
            string = 'Solute concentration'
            iflag = PETSC_TRUE
          endif
        case(OBS_ERT_MEASUREMENT)
          if (this%realization%option%ngeopdof == 0) then
            string = 'ERT measurement'
            iflag = PETSC_TRUE
          endif
        case default
          call this%driver%PrintErrMsg('Unknown observation type in &
            &InvSubsurfSetupForwardRunLinkage: ' // &
            trim(StringWrite(this%measurements(i)%iobs_var)))
      end select
      if (iflag) then
        call this%driver%PrintErrMsg(trim(string) // ' is specified as a &
          &measurement for inversion, but it is not being simulated.')
      endif
    enddo

    max_int = UNINITIALIZED_INTEGER
    allocate(int_array(num_measurements))
    int_array = -1
    do i = 1, num_measurements
      if (this%measurements(i)%iobs_var == OBS_ERT_MEASUREMENT) then
        max_int(1) = max(max_int(1),this%measurements(i)%cell_id)
      else
        if (Initialized(this%measurements(i)%coordinate%x)) then
          call GridGetLocalIDFromCoordinate(patch%grid, &
                                            this%measurements(i)%coordinate, &
                                            this%realization%option,local_id)
          if (Initialized(local_id)) then
            this%measurements(i)%cell_id = &
              patch%grid%nG2A(patch%grid%nL2G(local_id))
            this%measurements(i)%local_id = local_id
          endif
        endif
        max_int(2) = max(max_int(2),this%measurements(i)%cell_id)
        int_array(i) = this%measurements(i)%cell_id
      endif
    enddo
    ! ensure that cell and ert measurement ids are within bounds
    call MPI_Allreduce(MPI_IN_PLACE,max_int,TWO_INTEGER_MPI,MPIU_INTEGER, &
                       MPI_MAX,this%driver%comm%mycomm,ierr);CHKERRQ(ierr)
    if (associated(this%realization%survey)) then
      if (max_int(1) > size(this%realization%survey%dsim)) then
        call this%driver%PrintErrMsg('The ERT_MEASUREMENT_ID assigned to a &
          &measurement is greater than the number of measurements in an &
          &survey (survey size = ' // &
          trim(StringWrite(num_measurements)) // &
          ' vs. maximum ERT_MEASUREMENT_ID = ' // &
          trim(StringWrite(size(this%realization%survey%dsim))) // ').')
      endif
    endif
    if (max_int(1) > patch%grid%nmax) then
      call this%driver%PrintErrMsg('A measurement cell ID is &
        &beyond the maximum cell ID of ' // StringWrite(i))
    endif
    ! ensure that all cell ids have been found
    mpi_int = num_measurements
    call MPI_Allreduce(MPI_IN_PLACE,int_array,mpi_int,MPIU_INTEGER,MPI_MAX, &
                      this%driver%comm%mycomm,ierr);CHKERRQ(ierr)
    i = maxval(int_array)
    do i = 1, num_measurements
      if (int_array(i) > 0) then
        this%measurements(i)%cell_id = int_array(i)
      endif
      if (Uninitialized(this%measurements(i)%cell_id)) then
        string = 'Measurement ' // trim(StringWrite(i)) // &
          ' at coordinate (' // &
          trim(StringWrite(this%measurements(i)%coordinate%x)) // ',' // &
          trim(StringWrite(this%measurements(i)%coordinate%y)) // ',' // &
          trim(StringWrite(this%measurements(i)%coordinate%z)) // &
          ') not mapped properly.'
        call this%driver%PrintErrMsg(string)
      endif
    enddo

    int_array = -1
    iflag = PETSC_FALSE ! no ERT measurements
    do i = 1, num_measurements
      if (this%measurements(i)%iobs_var /= OBS_ERT_MEASUREMENT) then
        int_array(i) = this%measurements(i)%cell_id
      else
        iflag = PETSC_TRUE
      endif
    enddo
    ! throw error if ert is included as a measurement when using adjoints
    if (iflag .and. .not.associated(this%perturbation)) then
      call this%driver%PrintErrMsg('ERT measurements are currently not &
        &supported adjoint-based inversion.')
    endif
    int_array = int_array - 1
    call DiscretAOApplicationToPetsc(this%realization%discretization, &
                                     int_array)

    ! create mapping between local measurements and redundant measurement_vec
    ! find local measurements
    temp_int = 0
    call MPI_Exscan(patch%grid%nlmax,temp_int,ONE_INTEGER_MPI,MPIU_INTEGER, &
                    MPI_SUM,this%driver%comm%mycomm,ierr);CHKERRQ(ierr)
    ! Exscan does not include the temp_int from this rank
    ! local ids are between temp_int and temp_int + grid%nlmax
    ! (0 < id <= grid%nlmax on process 0; it is zero-based)
    icount = 0
    do i = 1, num_measurements
      if (this%measurements(i)%iobs_var /= OBS_ERT_MEASUREMENT) then
        if (int_array(i) >= temp_int .and. &
            int_array(i) < temp_int+patch%grid%nlmax) then
          ! it is local
          icount = icount + 1
          local_id = int_array(i) - temp_int + 1
          if (Initialized(this%measurements(i)%local_id) .and. &
              this%measurements(i)%local_id /= local_id) then
            this%realization%option%io_buffer = 'Error mapping local id &
              &for measurement ' // trim(StringWrite(i))
            call PrintErrMsgByRank(this%realization%option)
          endif
          this%measurements(i)%local_id = local_id
        endif
      else if (OptionIsIORank(this%realization%option)) then
        icount = icount + 1
        this%measurements(i)%local_id = this%measurements(i)%cell_id
      endif
    enddo
    allocate(this%local_measurement_values(icount))
    allocate(this%local_measurement_map(icount))
    if (.not.associated(this%perturbation)) then
      ! if adjoint, need to allocate array for potential partial derivatives
      allocate(this%local_dobs_dunknown_values(icount))
      allocate(this%local_dobs_dparam_values(icount))
    endif

    this%local_measurement_map = UNINITIALIZED_INTEGER
    icount = 0
    do i = 1, num_measurements
      if (this%measurements(i)%iobs_var /= OBS_ERT_MEASUREMENT) then
        if (int_array(i) >= temp_int .and. &
            int_array(i) < temp_int+patch%grid%nlmax) then
          ! it is local
          icount = icount + 1
          this%local_measurement_map(icount) = i
        endif
      else if (OptionIsIORank(this%realization%option)) then
        icount = icount + 1
        this%local_measurement_map(icount) = i
      endif
    enddo

    ! map measurement vec to distributed measurement vec
    call ISCreateStride(this%driver%comm%mycomm,num_measurements,ZERO_INTEGER, &
                        ONE_INTEGER,is_measure,ierr);CHKERRQ(ierr)
    call VecScatterCreate(this%measurement_vec,is_measure, &
                          this%dist_measurement_vec,PETSC_NULL_IS, &
                          this%scatter_measure_to_dist_measure, &
                          ierr);CHKERRQ(ierr)
    call ISDestroy(is_measure,ierr)

    this%inversion_aux%measurements => this%measurements
    this%inversion_aux%measurement_vec = this%measurement_vec

    if (associated(this%perturbation)) then
      if (this%qoi_is_full_vector) then
        if (associated(this%perturbation%select_cells)) then
          this%perturbation%ndof = size(this%perturbation%select_cells)
          if (this%perturbation%ndof > this%realization%patch%grid%nmax) then
            call this%driver%PrintErrMsg('Number of SELECT_CELLS is larger &
                                        &than the problem size: '// &
                      trim(StringWrite(this%perturbation%ndof))//' '// &
                      trim(StringWrite(this%realization%patch%grid%nmax)))
          endif
        else
          this%perturbation%ndof = this%realization%patch%grid%nmax
        endif
      else
        this%perturbation%ndof = size(this%parameters)
      endif
      if (.not.this%inversion_option%coupled_flow_ert) then
        call VecDuplicate(this%measurement_vec, &
                          this%perturbation%base_measurement_vec, &
                          ierr);CHKERRQ(ierr)
      endif
    endif

    inversion_forward_aux => InversionForwardAuxCreate()
    inversion_forward_aux%JsensitivityT_ptr = this%inversion_aux%JsensitivityT
    inversion_forward_aux%measurements => this%measurements
    inversion_forward_aux%measurement_vec = this%measurement_vec
    inversion_forward_aux%inversion_coupled_aux => this%inversion_coupled_aux
    inversion_forward_aux%local_measurement_values_ptr => &
      this%local_measurement_values
    if (.not.associated(this%perturbation)) then
      inversion_forward_aux%iparameter = this%parameters(1)%iparameter
      inversion_forward_aux%local_dobs_dunknown_values_ptr => &
        this%local_dobs_dunknown_values
      inversion_forward_aux%local_dobs_dparam_values_ptr => &
        this%local_dobs_dparam_values
    endif
    this%inversion_aux%inversion_forward_aux => inversion_forward_aux

    ! if permeability is the parameter of interest, ensure that it is
    ! isotropic or specified with a vertical anisotropy ratio
    iflag = PETSC_FALSE
    if (this%qoi_is_full_vector) then
      if (MaterialAnisotropyExists(patch%material_properties)) then
        iflag = PETSC_TRUE
      endif
    else
      do i = 1, size(this%parameters)
        if (this%parameters(i)%iparameter == PERMEABILITY) then
          material_property => &
            patch%material_property_array(this%parameters(i)%imat)%ptr
          ! if PERM_HORIZONTAL and VERTICAL_ANISOTROPY_RATIO are used,
          ! isotropic permeabilithy will be false, but tempreal will be 1.
          ! this only works with perturbation
          tempreal = material_property%permeability(3,3) / &
                     material_property%permeability(1,1) / &
                     material_property%vertical_anisotropy_ratio
          if (.not.material_property%isotropic_permeability .and. &
              (.not.associated(this%perturbation) .or. &
               .not.(tempreal > 0.999d0 .and. tempreal < 1.001d0))) then
            iflag = PETSC_TRUE
            exit
          endif
        endif
      enddo
    endif
    if (iflag) then
      call this%driver%PrintErrMsg('Anisotropic permeability is a parameter &
                      &of interest in the forward simulation and is not &
                      &supported for inversion.')
    endif

  endif

  inversion_forward_aux => this%inversion_aux%inversion_forward_aux
  this%local_measurement_values = UNINITIALIZED_DOUBLE

  if (.not.associated(this%perturbation)) then
    ! set up pointer to M matrix
    inversion_forward_aux%M_ptr = &
        this%forward_simulation%flow_process_model_coupler%timestepper%solver%M
    ! create inversion_ts_aux for first time step
    nullify(inversion_forward_aux%first) ! must pass in null object
    inversion_forward_aux%first => &
      InversionTSAuxCreate(inversion_forward_aux%first, &
                           inversion_forward_aux%M_ptr)
    inversion_forward_aux%last => inversion_forward_aux%first
    call InvTSAuxAllocate(inversion_forward_aux%first, &
                          inversion_forward_aux%M_ptr, &
                          this%realization%option%nflowdof, &
                          patch%grid%nlmax)
    this%local_dobs_dunknown_values = UNINITIALIZED_DOUBLE
    this%local_dobs_dparam_values = UNINITIALIZED_DOUBLE

  else ! this%perturbation is associated

    ! destroy all flow adjoint data structures
    call InvForwardAuxDestroyList(this%inversion_aux%inversion_forward_aux, &
                                  PETSC_FALSE)
    this%inversion_aux%inversion_forward_aux%store_adjoint = PETSC_FALSE

    if (this%inversion_option%coupled_flow_ert) then
      ! ndof is always the number of parameters
      this%perturbation%ndof = size(this%parameters)
      this%realization%option%inversion%calculate_ert = &
        this%perturbation%idof_pert <= 0
      this%realization%option%inversion%record_measurements = &
        this%perturbation%idof_pert <= 0
      this%realization%option%inversion%calculate_ert_jacobian = &
        this%perturbation%idof_pert < 0
    else ! normal perturbation
      if (this%qoi_is_full_vector .and. &
          this%perturbation%idof_pert > &
            this%realization%patch%grid%nmax) then
        call this%driver%PrintErrMsg('SELECT_CELLS ID is larger than &
                                    &the problem size: '// &
                trim(StringWrite(this%perturbation%idof_pert))//' '// &
                trim(StringWrite(this%realization%patch%grid%nmax)))
      endif
    endif

    if (Uninitialized(this%parameters(1)%iparameter) .and. &
        this%qoi_is_full_vector) then
      call this%driver%PrintErrMsg('Quantity of interest not specified in &
        &InvSubsurfSetupForwardRunLinkage.')
    endif
  endif

end subroutine InvSubsurfSetupForwardRunLinkage

! ************************************************************************** !

subroutine InvSubsurfConnectToForwardRun(this)
  !
  ! Sets up the interface between a single forward run and the outer
  ! inversion wrapper
  !
  ! Author: Glenn Hammond
  ! Date: 10/18/21
  !
  use Discretization_module
  use Factory_Subsurface_module
  use Init_Subsurface_module
  use Material_module
  use String_module
  use Utility_module
  use Waypoint_module
  use ZFlow_Aux_module

  class(inversion_subsurface_type) :: this

  PetscReal :: rmin, rmax
  PetscReal :: final_time
  PetscInt :: iqoi(2)
  PetscInt :: i, sync_count
  type(waypoint_type), pointer :: waypoint
  PetscReal, pointer :: real_array(:)
  PetscBool :: iflag
  PetscErrorCode :: ierr

  ! insert measurement times into waypoint list. this must come after the
  ! simulation is initialized to obtain the final time.
  if (.not.associated(this%inversion_aux% &
                        inversion_forward_aux%sync_times)) then
    allocate(real_array(size(this%measurements)))
    final_time = &
      WaypointListGetFinalTime(this%forward_simulation%waypoint_list_subsurface)
    iflag = PETSC_FALSE
    do i = 1, size(this%measurements)
      if (Uninitialized(this%measurements(i)%time)) then
        this%measurements(i)%time = final_time
      endif
      if (this%measurements(i)%time > final_time) then
        iflag = PETSC_TRUE
      endif
      real_array(i) = this%measurements(i)%time
    enddo
    if (iflag) then
      call this%driver%PrintErrMsg( &
            'Inversion measurement times are specified beyond the end of &
            &the final simulation time.')
    endif
    call UtilitySortArray(real_array)
    sync_count = 0
    do i = 1, size(real_array)
      iflag = PETSC_FALSE
      if (i == 1) then
        iflag = PETSC_TRUE
      else if (real_array(i) > real_array(i-1)) then
        iflag = PETSC_TRUE
      endif
      if (iflag) then
        sync_count = sync_count + 1
        real_array(sync_count) = real_array(i)
      endif
    enddo
    allocate(this%inversion_aux%inversion_forward_aux%sync_times(sync_count))
    this%inversion_aux%inversion_forward_aux%sync_times(:) = &
      real_array(1:sync_count)
    deallocate(real_array)
    nullify(real_array)

    if (associated(this%inversion_coupled_aux)) then
      allocate(this%inversion_coupled_aux%solutions(sync_count))
      do i = 1, sync_count
        call InversionCoupledSolutionInit(this%inversion_coupled_aux% &
                                            solutions(i))
        this%inversion_coupled_aux%solutions(i)%time = &
          this%inversion_aux%inversion_forward_aux%sync_times(i)
      enddo

      ! allocate any full vector measurements
      call InvCoupledAllocateSolnVecs(this%inversion_coupled_aux, &
                                      this%realization%field%work)
    endif

  endif

  real_array => this%inversion_aux%inversion_forward_aux%sync_times
  do i = 1, size(real_array)
    waypoint => WaypointCreate()
    waypoint%time = real_array(i)
    waypoint%sync = PETSC_TRUE
    call WaypointInsertInList(waypoint, &
                              this%forward_simulation%waypoint_list_outer)
!                              this%forward_simulation%waypoint_list_subsurface)
  enddo

#if 0
  ! resets the waypoint pointers to the first entry, which may have changed
  ! due to the insertions above
  call WaypointListFillIn(this%forward_simulation%waypoint_list_subsurface, &
                          this%forward_simulation%option)
  call WaypointListRemoveExtraWaypnts(this%forward_simulation% &
                                        waypoint_list_subsurface, &
                                      this%forward_simulation%option)
  call WaypointListFindDuplicateTimes(this%forward_simulation% &
                                        waypoint_list_subsurface, &
                                      this%forward_simulation%option)
  call FactorySubsurfSetPMCWaypointPtrs(this%forward_simulation)
#endif

  ! must come after insertion of waypoints and setting of pointers; otherwise,
  ! the pmc timestepper cur_waypoint pointers are pointing at the time 0
  ! waypoint, instead of the one just after time 0, which is incorrect.
  call this%forward_simulation%InitializeRun()

  this%realization%patch%aux%inversion_forward_aux => &
    this%inversion_aux%inversion_forward_aux
  call InvForwardAuxResetMeasurements(this%inversion_aux%inversion_forward_aux)

  if (this%qoi_is_full_vector) then
    iqoi = InversionParameterIntToQOIArray(this%parameters(1))
    if (this%first_inversion_interation) then
      ! on first iteration of inversion, store the values
      call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                  this%realization%field%work_loc, &
                                  iqoi(1),iqoi(2))
      call DiscretizationLocalToGlobal(this%realization%discretization, &
                                      this%realization%field%work_loc, &
                                      this%realization%field%work,ONEDOF)
      call InvSubsurfScatGlobalToDistParam(this, &
                                           this%realization%field%work, &
                                           this%dist_parameter_vec, &
                                           INVSUBSCATFORWARD)
    else
      ! on subsequent iterations, overwrite what was read from input file
      ! with latest inverted values
      call InvSubsurfScatGlobalToDistParam(this, &
                                           this%realization%field%work, &
                                           this%dist_parameter_vec, &
                                           INVSUBSCATREVERSE)
      call DiscretizationGlobalToLocal(this%realization%discretization, &
                                       this%realization%field%work, &
                                       this%realization%field%work_loc,ONEDOF)
      call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                   this%realization%field%work_loc, &
                                   iqoi(1),iqoi(2))
    endif
  else
    if (this%first_inversion_interation) then
      ! load the original parameter values
      do i = 1, size(this%parameters)
        call InvSubsurfCopyParameterValue(this,i,GET_MATERIAL_VALUE)
      enddo
      call InvSubsurfCopyParameterToFromVec(this,COPY_TO_VEC)
    else
      call InvSubsurfCopyParameterToFromVec(this,COPY_FROM_VEC)
      do i = 1, size(this%parameters)
        call InvSubsurfCopyParameterValue(this,i,OVERWRITE_MATERIAL_VALUE)
      enddo
      ! update material auxvars
      call InitSubsurfAssignMatIDsToRegns(this%realization)
      call InitSubsurfAssignMatProperties(this%realization)
    endif
  endif

  if (associated(this%perturbation)) then
    ! on first pass, store and set thereafter
    if (this%perturbation%idof_pert <= 0) then
      if (this%perturbation%base_parameter_vec == PETSC_NULL_VEC) then
        if (this%qoi_is_full_vector) then
          call VecDuplicate(this%dist_parameter_vec, &
                            this%perturbation%base_parameter_vec, &
                            ierr);CHKERRQ(ierr)
        else
          call VecDuplicate(this%parameter_vec, &
                            this%perturbation%base_parameter_vec, &
                            ierr);CHKERRQ(ierr)
        endif
      endif
      if (this%qoi_is_full_vector) then
        call VecCopy(this%dist_parameter_vec, &
                     this%perturbation%base_parameter_vec,ierr);CHKERRQ(ierr)
      else
        if (this%perturbation%idof_pert == 0) then
          call VecCopy(this%parameter_vec,this%perturbation%base_parameter_vec, &
                      ierr);CHKERRQ(ierr)
        else
        !geh: try with and without
          call VecCopy(this%perturbation%base_parameter_vec,this%parameter_vec, &
                      ierr);CHKERRQ(ierr)
        endif
      endif
    else
      ! on subsequent passes, we have to overwrite the entire
      if (this%qoi_is_full_vector) then
        call VecZeroEntries(this%dist_parameter_vec,ierr);CHKERRQ(ierr)
        if (this%driver%comm%myrank == 0) then
          call VecSetValue(this%dist_parameter_vec, &
                           this%perturbation%idof_pert-1, &
                           this%perturbation%tolerance,INSERT_VALUES, &
                           ierr);CHKERRQ(ierr)
        endif
        call VecAssemblyBegin(this%dist_parameter_vec,ierr);CHKERRQ(ierr)
        call VecAssemblyEnd(this%dist_parameter_vec,ierr);CHKERRQ(ierr)
        call VecPointwiseMult(this%dist_parameter_vec,this%dist_parameter_vec, &
                              this%perturbation%base_parameter_vec, &
                              ierr);CHKERRQ(ierr)
        call VecMax(this%dist_parameter_vec,PETSC_NULL_INTEGER,rmax, &
                    ierr);CHKERRQ(ierr)
        call VecMin(this%dist_parameter_vec,PETSC_NULL_INTEGER,rmin, &
                    ierr);CHKERRQ(ierr)
        if (rmax > 0.d0) then
          this%perturbation%pert = rmax
        else
          this%perturbation%pert = rmin
        endif
        call VecAXPY(this%dist_parameter_vec,1.d0, &
                     this%perturbation%base_parameter_vec,ierr);CHKERRQ(ierr)
        call InvSubsurfScatGlobalToDistParam(this, &
                                             this%realization%field%work, &
                                             this%dist_parameter_vec, &
                                             INVSUBSCATREVERSE)

        call DiscretizationGlobalToLocal(this%realization%discretization, &
                                        this%realization%field%work, &
                                        this%realization%field%work_loc,ONEDOF)
        iqoi = InversionParameterIntToQOIArray(this%parameters(1))
        call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                    this%realization%field%work_loc, &
                                    iqoi(1),iqoi(2))
      else
        this%perturbation%base_value = &
          this%parameters(this%perturbation%idof_pert)%value
        this%perturbation%pert = this%perturbation%base_value* &
                                 this%perturbation%tolerance
        this%parameters(this%perturbation%idof_pert)%value = &
          this%perturbation%base_value + this%perturbation%pert
        ! overwrite material property value
        call InvSubsurfCopyParameterValue(this,this%perturbation%idof_pert, &
                                          OVERWRITE_MATERIAL_VALUE)
        call InitSubsurfAssignMatIDsToRegns(this%realization)
        call InitSubsurfAssignMatProperties(this%realization)
      endif
    endif
  else ! set adjoint variable
    ! pass in first parameter as an earlier check prevents adjoint-based
    ! inversion for more than one parameter type
    call InvSubsurfSetAdjointVariable(this,this%parameters(1)%iparameter)
  endif

  iflag = PETSC_TRUE
  if (associated(this%perturbation)) then
    iflag = this%perturbation%idof_pert == 0
  endif
  if (iflag) then
    call InvSubsurfPrintCurParamValues(this)
  endif

  this%first_inversion_interation = PETSC_FALSE

end subroutine InvSubsurfConnectToForwardRun

! ************************************************************************** !

subroutine InvSubsurfSetAdjointVariable(this,iparameter)
  !
  ! Sets the adjoint variable for a process model
  !
  ! Author: Glenn Hammond
  ! Date: 03/30/22

  use String_module
  use Variables_module, only : PERMEABILITY, POROSITY, VG_ALPHA, VG_M, VG_SR
  use ZFlow_Aux_module

  class(inversion_subsurface_type) :: this
  PetscInt :: iparameter

  character(len=MAXSTRINGLENGTH) :: string

  select case(this%realization%option%iflowmode)
    case(ZFLOW_MODE)
    case default
      string = 'Flow mode "' // trim(this%realization%option%flowmode) // &
        '" not supported for inversion (InvSubsurfSetAdjointVariable).'
    call this%driver%PrintErrMsg(string)
  end select

  select case(iparameter)
    case(PERMEABILITY)
      zflow_adjoint_parameter = ZFLOW_ADJOINT_PERMEABILITY
    case(POROSITY)
      zflow_adjoint_parameter = ZFLOW_ADJOINT_POROSITY
    case(VG_ALPHA,VG_M,VG_SR)
      string = 'van Genuchten parameters are unsupported for adjoint-based &
        &inversion.'
      call this%driver%PrintErrMsg(string)
    case default
      string = 'Unrecognized variable in InvSubsurfSetAdjointVariable: ' // &
               trim(StringWrite(iparameter))
      call this%driver%PrintErrMsg(string)
  end select

end subroutine InvSubsurfSetAdjointVariable

! ************************************************************************** !

subroutine InvSubsurfCopyParameterValue(this,iparam,iflag)
  !
  ! Copies parameter values back and forth
  !
  ! Author: Glenn Hammond
  ! Date: 03/30/22

  use Material_module

  class(inversion_subsurface_type) :: this
  PetscInt :: iparam
  PetscInt :: iflag

  PetscReal :: tempreal

  if (iflag /= GET_MATERIAL_VALUE) then
    ! everything else is implicit OVERWRITE_MATERIAL_VALUE
    tempreal = this%parameters(iparam)%value
  endif

  call InvSubsurfGetSetParamValueByMat(this,tempreal, &
                                       this%parameters(iparam)%iparameter, &
                                       this%parameters(iparam)%imat,iflag)

  if (iflag == GET_MATERIAL_VALUE) then
    this%parameters(iparam)%value = tempreal
  endif

end subroutine InvSubsurfCopyParameterValue

! ************************************************************************** !

subroutine InvSubsurfGetSetParamValueByMat(this,value,iparameter,imat,iflag)
  !
  ! Copies parameter values back and forth
  !
  ! Author: Glenn Hammond
  ! Date: 03/30/22

  use Characteristic_Curves_module
  use Material_module
  use String_module
  use Utility_module
  use Variables_module, only : ELECTRICAL_CONDUCTIVITY, PERMEABILITY, &
                               POROSITY, VG_ALPHA, VG_SR, VG_M

  class(inversion_subsurface_type) :: this
  PetscReal :: value
  PetscInt :: iparameter
  PetscInt :: imat
  PetscInt :: iflag

  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: tempreal
  type(material_property_type), pointer :: material_property
  type(characteristic_curves_type), pointer :: cc

  material_property => &
    this%realization%patch%material_property_array(imat)%ptr
  select case(iparameter)
    case(ELECTRICAL_CONDUCTIVITY)
      if (iflag == GET_MATERIAL_VALUE) then
        value = material_property%electrical_conductivity
      else
        material_property%electrical_conductivity = value
      endif
    case(PERMEABILITY)
      if (iflag == GET_MATERIAL_VALUE) then
        value = material_property%permeability(1,1)
      else
        material_property%permeability(1,1) = value
        material_property%permeability(2,2) = value
        if (Initialized(material_property%vertical_anisotropy_ratio)) then
          value = value * material_property%vertical_anisotropy_ratio
        endif
        material_property%permeability(3,3) = value
      endif
    case(POROSITY)
      if (iflag == GET_MATERIAL_VALUE) then
        value = material_property%porosity
      else
        material_property%porosity = value
      endif
    case(VG_ALPHA,VG_SR,VG_M)
      cc => this%realization%patch%characteristic_curves_array( &
              material_property%saturation_function_id)%ptr
      select case(iparameter)
        case(VG_ALPHA)
          if (iflag == GET_MATERIAL_VALUE) then
            value = cc%saturation_function%GetAlpha_()
          else
            call cc%saturation_function%SetAlpha_(value)
          endif
        case(VG_M)
          if (iflag == GET_MATERIAL_VALUE) then
            value = cc%saturation_function%GetM_()
            tempreal = cc%liq_rel_perm_function%GetM_()
            if (.not.Equal(value,tempreal)) then
              string = 'For inversion, saturation and relative permeability &
                &function van Genuchten "m" values must match in &
                &characteristic curve "' // trim(cc%name)
              call this%driver%PrintErrMsg(string)
            endif
          else
            call cc%saturation_function%SetM_(value)
            call cc%liq_rel_perm_function%SetM_(value)
          endif
        case(VG_SR)
          if (iflag == GET_MATERIAL_VALUE) then
            value = cc%saturation_function%GetResidualSaturation()
            tempreal = cc%liq_rel_perm_function%GetResidualSaturation()
            if (.not.Equal(value,tempreal)) then
              string = 'For inversion, saturation and relative permeability &
                &function  saturations must match in characteristic &
                &curve "' // trim(cc%name)
              call this%driver%PrintErrMsg(string)
            endif
          else
            call cc%saturation_function%SetResidualSaturation(value)
            call cc%liq_rel_perm_function%SetResidualSaturation(value)
          endif
      end select
    case default
      string = 'Unrecognized variable in &
        &InvSubsurfCopyParamValueByMat: ' // &
        trim(StringWrite(iparameter))
      call this%driver%PrintErrMsg(string)
  end select

end subroutine InvSubsurfGetSetParamValueByMat

! ************************************************************************** !

subroutine InvSubsurfCopyParameterToFromVec(this,iflag)
  !
  ! Copies parameter values back and forth
  !
  ! Author: Glenn Hammond
  ! Date: 03/30/22

  class(inversion_subsurface_type) :: this
  PetscInt :: iflag

  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: i
  PetscErrorCode :: ierr

  call VecGetArrayF90(this%parameter_vec,vec_ptr,ierr);CHKERRQ(ierr)
  if (iflag == COPY_FROM_VEC) then
    do i = 1, size(this%parameters)
      this%parameters(i)%value = vec_ptr(i)
    enddo
  elseif (iflag == COPY_TO_VEC) then
    do i = 1, size(this%parameters)
      vec_ptr(i) = this%parameters(i)%value
    enddo
  else
    call this%driver%PrintErrMsg('Unrecogized flag in &
                                 &InvSubsurfCopyParameterToFromVec')
  endif
  call VecRestoreArrayF90(this%parameter_vec,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine InvSubsurfCopyParameterToFromVec

! ************************************************************************** !

subroutine InvSubsurfGetParamValueByCell(this,value,iparameter,imat, &
                                         material_auxvar)
  !
  ! Returns the parameter value at the cell
  !
  ! Author: Glenn Hammond
  ! Date: 11/11/22

  use Characteristic_Curves_module
  use Material_module
  use Material_Aux_module, only : material_auxvar_type, &
                                  MaterialAuxVarGetValue
  use String_module
  use Variables_module, only : ELECTRICAL_CONDUCTIVITY, &
                               PERMEABILITY, PERMEABILITY_X, &
                               POROSITY, BASE_POROSITY, &
                               VG_ALPHA, VG_SR, VG_M

  class(inversion_subsurface_type) :: this
  PetscReal :: value
  PetscInt :: iparameter
  PetscInt :: imat
  type(material_auxvar_type) :: material_auxvar

  type(material_property_type), pointer :: material_property
  type(characteristic_curves_type), pointer :: cc

  select case(iparameter)
    case(ELECTRICAL_CONDUCTIVITY)
      value = MaterialAuxVarGetValue(material_auxvar,ELECTRICAL_CONDUCTIVITY)
    case(PERMEABILITY)
      value = MaterialAuxVarGetValue(material_auxvar,PERMEABILITY_X)
    case(POROSITY)
      value = MaterialAuxVarGetValue(material_auxvar,BASE_POROSITY)
    case(VG_ALPHA,VG_SR,VG_M)
      material_property => &
        this%realization%patch%material_property_array(imat)%ptr
      cc => this%realization%patch%characteristic_curves_array( &
              material_property%saturation_function_id)%ptr
      select case(iparameter)
        case(VG_ALPHA)
          value = cc%saturation_function%GetAlpha_()
        case(VG_M)
          value = cc%saturation_function%GetM_()
        case(VG_SR)
          value = cc%saturation_function%GetResidualSaturation()
      end select
    case default
      call this%driver%PrintErrMsg('Unrecognized variable in &
                                   &InvSubsurfGetParamValueByCell: ' // &
                                   trim(StringWrite(iparameter)))
  end select

end subroutine InvSubsurfGetParamValueByCell

! ************************************************************************** !

subroutine InvSubsurfExecuteForwardRun(this)
  !
  ! Executes a forward simulation
  !
  ! Author: Glenn Hammond
  ! Date: 03/21/22

  class(inversion_subsurface_type) :: this

  if (this%realization%option%status == PROCEED) then
    call this%forward_simulation%ExecuteRun()
    call InvSubsurfPostProcMeasurements(this)
  endif

end subroutine InvSubsurfExecuteForwardRun

! ************************************************************************** !

subroutine InvSubsurfCalculateSensitivity(this)
  !
  ! Calculates sensitivity matrix Jsensitivity
  !
  ! Author: Glenn Hammond
  ! Date: 03/25/22
  !
  use Option_module
  use String_module
  use Units_module
  use Utility_module

  class(inversion_subsurface_type) :: this

  if (associated(this%perturbation)) then
    call InvSubsurfPertCalcSensitivity(this)
  else
    call InvSubsurfAdjointCalcSensitivity(this)
  endif

  call InvSubsurfOutputSensitivity(this,'')

  if (.not.this%qoi_is_full_vector) then
    call InvSubsurfScatParamToDistParam(this, &
                                        this%parameter_vec, &
                                        this%dist_parameter_vec, &
                                        INVSUBSCATFORWARD)
  endif

end subroutine InvSubsurfCalculateSensitivity

! ************************************************************************** !

subroutine InvSubsurfPostProcMeasurements(this)
  !
  ! Stores measurements in appropriate arrays and vectors
  !
  ! Author: Glenn Hammond
  ! Date: 06/22/22
  !
  use Option_module
  use String_module
  use Units_module
  use Utility_module

  class(inversion_subsurface_type) :: this

  type(option_type), pointer :: option
  PetscInt :: imeasurement
  character(len=MAXWORDLENGTH) :: word
  Vec :: dobs_dunknown_vec
  Vec :: dist_dobs_dunknown_vec
  Vec :: dobs_dparam_vec
  Vec :: dist_dobs_dparam_vec
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: vec_ptr2(:)
  PetscInt :: icount
  PetscBool :: iflag
  PetscErrorCode :: ierr

  option => this%realization%option

  ! distribute measurements to measurement objects
  call VecSet(this%measurement_vec,UNINITIALIZED_DOUBLE,ierr);CHKERRQ(ierr)
  call VecSet(this%dist_measurement_vec,-888.d0,ierr);CHKERRQ(ierr)
  icount = 0
  do imeasurement = 1, size(this%measurements)
    if (Initialized(this%measurements(imeasurement)%local_id)) then
      icount = icount + 1
      call VecSetValue(this%dist_measurement_vec, &
                       this%local_measurement_map(icount)-1, &
                       this%local_measurement_values(icount),&
                       INSERT_VALUES,ierr);CHKERRQ(ierr)
    endif
  enddo
  call VecAssemblyBegin(this%dist_measurement_vec,ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(this%dist_measurement_vec,ierr);CHKERRQ(ierr)
  call InvSubsurfScatMeasToDistMeas(this,this%measurement_vec, &
                                    this%dist_measurement_vec, &
                                    INVSUBSCATREVERSE)
  call VecGetArrayF90(this%measurement_vec,vec_ptr,ierr);CHKERRQ(ierr)
  do imeasurement = 1, size(this%measurements)
    this%measurements(imeasurement)%simulated_value = vec_ptr(imeasurement)
    this%measurements(imeasurement)%measured = PETSC_TRUE
  enddo
  call VecRestoreArrayF90(this%measurement_vec,vec_ptr,ierr);CHKERRQ(ierr)
  if (associated(this%local_dobs_dunknown_values)) then
    ! distribute derivatives to measurement objects
    ! temporary vecs for derivatives (if they exist)
    call VecDuplicate(this%measurement_vec,dobs_dunknown_vec,ierr);CHKERRQ(ierr)
    call VecDuplicate(this%dist_measurement_vec,dist_dobs_dunknown_vec,&
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(this%measurement_vec,dobs_dparam_vec,ierr);CHKERRQ(ierr)
    call VecDuplicate(this%dist_measurement_vec,dist_dobs_dparam_vec,&
                      ierr);CHKERRQ(ierr)
    call VecSet(dobs_dunknown_vec,UNINITIALIZED_DOUBLE,ierr);CHKERRQ(ierr)
    call VecSet(dist_dobs_dunknown_vec,-888.d0,ierr);CHKERRQ(ierr)
    call VecSet(dobs_dparam_vec,UNINITIALIZED_DOUBLE,ierr);CHKERRQ(ierr)
    ! dist_dobs_dparam_vec has to be UNINITIALIZED_DOUBLE, otherwise -888 will
    ! be set for all uninitialized values
    call VecSet(dist_dobs_dparam_vec,UNINITIALIZED_DOUBLE,ierr);CHKERRQ(ierr)
    icount = 0
    do imeasurement = 1, size(this%measurements)
      if (Initialized(this%measurements(imeasurement)%local_id)) then
        icount = icount + 1
        ! set the partial derivative
        select case(this%measurements(imeasurement)%iobs_var)
          case(OBS_LIQUID_SATURATION)
            call VecSetValue(dist_dobs_dunknown_vec, &
                             this%local_measurement_map(icount)-1, &
                             this%local_dobs_dunknown_values(icount),&
                             INSERT_VALUES,ierr);CHKERRQ(ierr)
          case(OBS_LIQUID_PRESSURE)
            call VecSetValue(dist_dobs_dparam_vec, &
                             this%local_measurement_map(icount)-1, &
                             this%local_dobs_dparam_values(icount),&
                             INSERT_VALUES,ierr);CHKERRQ(ierr)
        end select
      endif
    enddo
    call VecAssemblyBegin(dist_dobs_dunknown_vec,ierr);CHKERRQ(ierr)
    call VecAssemblyEnd(dist_dobs_dunknown_vec,ierr);CHKERRQ(ierr)
    call VecAssemblyBegin(dist_dobs_dparam_vec,ierr);CHKERRQ(ierr)
    call VecAssemblyEnd(dist_dobs_dparam_vec,ierr);CHKERRQ(ierr)
    call InvSubsurfScatMeasToDistMeas(this,dobs_dunknown_vec, &
                                      dist_dobs_dunknown_vec, &
                                      INVSUBSCATREVERSE)
    call InvSubsurfScatMeasToDistMeas(this,dobs_dparam_vec, &
                                      dist_dobs_dparam_vec, &
                                      INVSUBSCATREVERSE)
    call VecGetArrayF90(dobs_dunknown_vec,vec_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(dobs_dparam_vec,vec_ptr2,ierr);CHKERRQ(ierr)
    do imeasurement = 1, size(this%measurements)
      this%measurements(imeasurement)%dobs_dunknown = &
        vec_ptr(imeasurement)
      this%measurements(imeasurement)%dobs_dparam = &
        vec_ptr2(imeasurement)
    enddo
    call VecRestoreArrayF90(dobs_dunknown_vec,vec_ptr,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(dobs_dparam_vec,vec_ptr,ierr);CHKERRQ(ierr)
    call VecDestroy(dist_dobs_dunknown_vec,ierr);CHKERRQ(ierr)
    call VecDestroy(dobs_dunknown_vec,ierr);CHKERRQ(ierr)
    call VecDestroy(dist_dobs_dparam_vec,ierr);CHKERRQ(ierr)
    call VecDestroy(dobs_dparam_vec,ierr);CHKERRQ(ierr)
  endif

    ! ensure that all measurement have been recorded
  do imeasurement = 1, size(this%measurements)
    call InversionMeasurementPrintConcise(this%measurements(imeasurement), &
                                          'PostProcess',option)
    if (.not.this%measurements(imeasurement)%measured) then
      option%io_buffer = 'Measurement at cell ' // &
        trim(StringWrite(this%measurements(imeasurement)%cell_id))
      if (Initialized(this%measurements(imeasurement)%time)) then
        word = 'sec'
        option%io_buffer = trim(option%io_buffer) // &
          ' at ' // trim(StringWrite(this%measurements(imeasurement)%time/ &
          UnitsConvertToInternal(this%measurements(imeasurement)%time_units, &
                                 word,option,ierr))) // &
          ' ' // trim(this%measurements(imeasurement)%time_units)
      endif
      option%io_buffer = trim(option%io_buffer) // &
        ' with a measured value from the measurement file of ' // &
        trim(StringWrite(this%measurements(imeasurement)%value)) // &
        ' was not measured during the simulation.'
      call PrintErrMsg(option)
    endif
  enddo

  iflag = PETSC_TRUE
  if (associated(this%perturbation)) then
    iflag = this%perturbation%idof_pert == 0
  endif
  if (iflag) then
    call InvSubsurfPrintCurMeasValues(this)
  endif

end subroutine InvSubsurfPostProcMeasurements

! ************************************************************************** !

subroutine InvSubsurfAdjointCalcSensitivity(this)
  !
  ! Calculates sensitivity matrix Jsensitivity
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/21
  !
  use Option_module
  use String_module
  use Timer_class
  use Utility_module

  class(inversion_subsurface_type) :: this

  type(inversion_aux_type), pointer :: inversion_aux
  type(inversion_forward_ts_aux_type), pointer :: cur_inversion_ts_aux
  type(inversion_forward_ts_aux_type), pointer :: prev_inversion_ts_aux
  type(option_type), pointer :: option
  class(timer_type), pointer :: timer
  PetscInt :: imeasurement

  option => this%realization%option
  inversion_aux => this%inversion_aux

  timer => TimerCreate()

  call timer%Start()

  call PrintHeader('SENSITIVITY JACOBIAN',option)

  ! initialize first_lambda flag
  do imeasurement = 1, size(this%measurements)
    this%measurements(imeasurement)%first_lambda = PETSC_FALSE
  enddo

  ! go to end of list
  cur_inversion_ts_aux => inversion_aux%inversion_forward_aux%last
  if (.not.associated(cur_inversion_ts_aux)) then
    option%io_buffer = 'Inversion timestep auxiliary list is NULL.'
    call PrintErrMsg(option)
  endif

  ! the last link should be allocated, but not populated. this is by design
  if (cur_inversion_ts_aux%dResdu /= PETSC_NULL_MAT) then
    option%io_buffer = 'Last link in Inversion timestep auxiliary list &
      &is not NULL.'
    call PrintErrMsg(option)
  else
    ! remove the last link
    prev_inversion_ts_aux => cur_inversion_ts_aux%prev
    if (.not.associated(prev_inversion_ts_aux)) then
      option%io_buffer = 'Next to last link in Inversion timestep &
        &auxiliary list is NULL.'
      call PrintErrMsg(option)
    endif
    nullify(prev_inversion_ts_aux%next)
    call InversionTSAuxDestroy(cur_inversion_ts_aux)
    ! point cur_inversion_ts_aux to the end of the list
    inversion_aux%inversion_forward_aux%last => prev_inversion_ts_aux
    this%inversion_aux%max_ts = prev_inversion_ts_aux%timestep
  endif

  call InvSubsurfAdjAddSensitivities(this)

  call InvForwardAuxDestroyList(inversion_aux%inversion_forward_aux, &
                                PETSC_TRUE)

  call timer%Stop()
  option%io_buffer = '    ' // &
    trim(StringWrite('(f20.1)',timer%GetCumulativeTime())) // &
    ' seconds to build all sensitivities.'
  call PrintMsg(option)
  call TimerDestroy(timer)

end subroutine InvSubsurfAdjointCalcSensitivity

! ************************************************************************** !

subroutine InvSubsurfAdjAddSensitivities(this)
  !
  ! Calculates lambda and resulting sensitivities and adds them the
  ! sensitivity matrix
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/21
  !
  use Connection_module
  use Debug_module
  use Discretization_module
  use Grid_module
  use Option_module
  use Patch_module
  use Realization_Base_class
  use Solver_module
  use String_module
  use Timer_class
  use ZFlow_Aux_module
  use Units_module
  use Utility_module
  use Variables_module

  class(inversion_subsurface_type) :: this

  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(inversion_forward_ts_aux_type), pointer :: inversion_forward_ts_aux
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(solver_type), pointer :: solver
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: vec_ptr2(:)
  PetscReal :: tempreal
  PetscInt :: imeasurement
  PetscInt :: icell_measurement
  PetscInt :: tempint
  PetscInt :: ndof
  PetscInt :: i, j, offset, local_id
  Vec :: work  ! a 1 dof vec
  Vec :: onedof_vec
  Vec :: ndof_vec ! ndof_vec
  Vec :: p_       ! ndof_vec
  Vec :: rhs      ! ndof_vec
  ! derivative of residual at k+1 time level wrt unknown at k time level
  ! times lambda at k time level
  Vec :: dReskp1_duk_lambdak
  Vec :: natural_vec

  ! for forward loop sensitivity calculation
  Vec :: dResdKLambda
  Vec :: ndof_vec1
  Vec :: ndof_vec2
  PetscInt :: iparameter
  PetscInt :: natural_id_
  Vec, pointer :: lambda(:)

  class(timer_type), pointer :: timer
  PetscReal, parameter :: tol = 1.d-6
  PetscErrorCode :: ierr

  nullify(vec_ptr)
  nullify(vec_ptr2)
  nullify(lambda)

  solver => this%forward_simulation%flow_process_model_coupler% &
              timestepper%solver
  option => this%realization%option
  discretization => this%realization%discretization
  patch => this%realization%patch
  grid => patch%grid
  ndof = option%nflowdof

  timer => TimerCreate()

  call timer%Start()

  work = this%realization%field%work ! DO NOT DESTROY!
  ndof_vec = this%realization%field%flow_xx ! DO NOT DESTROY!

  dResdKLambda = PETSC_NULL_VEC
  ndof_vec1 = PETSC_NULL_VEC
  ndof_vec2 = PETSC_NULL_VEC

  call VecDuplicate(work,onedof_vec,ierr);CHKERRQ(ierr)
  call VecDuplicate(ndof_vec,p_,ierr);CHKERRQ(ierr)
  call VecDuplicate(ndof_vec,rhs,ierr);CHKERRQ(ierr)
  call VecDuplicate(ndof_vec,dReskp1_duk_lambdak,ierr);CHKERRQ(ierr)
  if (this%qoi_is_full_vector) then
    call VecDuplicate(ndof_vec,dResdKLambda,ierr);CHKERRQ(ierr)
  else
    call VecDuplicate(ndof_vec,ndof_vec1,ierr);CHKERRQ(ierr)
    call VecDuplicate(ndof_vec1,ndof_vec2,ierr);CHKERRQ(ierr)
  endif
  ! must use natural vec, not parameter vec since offset is based on measurement
  call DiscretizationCreateVector(discretization,ONEDOF, &
                                  natural_vec,NATURAL,option)

  call VecDuplicateVecsF90(ndof_vec,size(this%measurements),lambda, &
                          ierr);CHKERRQ(ierr)
  do i = 1, size(lambda)
    call VecZeroEntries(lambda(i),ierr);CHKERRQ(ierr)
  enddo

  call MatZeroEntries(this%inversion_aux%JsensitivityT,ierr);CHKERRQ(ierr)

  inversion_forward_ts_aux => this%inversion_aux%inversion_forward_aux%last
  do
    if (.not.associated(inversion_forward_ts_aux)) exit

    call KSPSetOperators(solver%ksp,inversion_forward_ts_aux%dResdu, &
                        inversion_forward_ts_aux%dResdu,ierr);CHKERRQ(ierr)
    do imeasurement = 1, size(this%measurements)

      ! skip lambdas after the measurement time
      if (Initialized(this%measurements(imeasurement)%time) .and. &
          inversion_forward_ts_aux%time > &
            this%measurements(imeasurement)%time + tol) then
        cycle
      endif

      ! backward loop contribution (calculating lambda for the ts)
      if (.not.this%measurements(imeasurement)%first_lambda) then
        this%measurements(imeasurement)%first_lambda = PETSC_TRUE

        ! if we have a dobs_dparam, add that value into Jsens
        if (Initialized(this%measurements(imeasurement)%dobs_dparam)) then
          if (this%qoi_is_full_vector .or. size(this%parameters) > 0) then
            option%io_buffer = 'Need to refactor dobs_dparam for more than &
              &one parameter.'
            call PrintErrMsg(option)
          endif
          if (option%comm%myrank == option%driver%io_rank) then
            tempreal = -this%measurements(imeasurement)%dobs_dparam
            iparameter = 1
            call MatSetValue(this%inversion_aux%JsensitivityT,iparameter-1, &
                            imeasurement-1,tempreal,ADD_VALUES, &
                            ierr);CHKERRQ(ierr)
          endif
        endif

        ! begin lambda calculations
        call VecZeroEntries(natural_vec,ierr);CHKERRQ(ierr)
        if (option%myrank == 0) then
          icell_measurement = this%measurements(imeasurement)%cell_id
          tempreal = -1.d0
          ! scale by partial derivative if measurement is not the
          ! primary dependent variable
          select case(this%measurements(imeasurement)%iobs_var)
            case(OBS_LIQUID_SATURATION)
              tempreal = tempreal * &
                        this%measurements(imeasurement)%dobs_dunknown
          end select
          call VecSetValue(natural_vec,icell_measurement-1,tempreal, &
                          INSERT_VALUES,ierr);CHKERRQ(ierr)
        endif
        call InversionMeasurementPrintConcise( &
               this%measurements(imeasurement),'',option)
        call VecAssemblyBegin(natural_vec,ierr);CHKERRQ(ierr)
        call VecAssemblyEnd(natural_vec,ierr);CHKERRQ(ierr)
        call DiscretizationNaturalToGlobal(discretization,natural_vec, &
                                          onedof_vec,ONEDOF)
        select case(this%measurements(imeasurement)%iobs_var)
          case(OBS_LIQUID_PRESSURE,OBS_LIQUID_SATURATION)
            tempint = zflow_liq_flow_eq
          case(OBS_SOLUTE_CONCENTRATION)
            tempint = zflow_sol_tran_eq
          case(OBS_ERT_MEASUREMENT)
            option%io_buffer = 'ERT measurement are currently not supported &
              &adjoint-based inversion.'
            call PrintErrMsg(option)
          case default
            option%io_buffer = 'Unknown observation type in &
              &InvSubsurfCalcLambda'
            call PrintErrMsg(option)
        end select
        if (Uninitialized(tempint)) then
          option%io_buffer = 'The observed state variable is not being &
            &modeled.'
          call PrintErrMsg(option)
        endif
        call VecZeroEntries(p_,ierr);CHKERRQ(ierr)
        call VecStrideScatter(onedof_vec,tempint-1,p_,INSERT_VALUES, &
                              ierr);CHKERRQ(ierr)
        call VecZeroEntries(dReskp1_duk_lambdak,ierr);CHKERRQ(ierr)
      else
        call VecZeroEntries(p_,ierr);CHKERRQ(ierr)
        call VecZeroEntries(dReskp1_duk_lambdak,ierr);CHKERRQ(ierr)
        call VecGetArrayF90(dReskp1_duk_lambdak,vec_ptr,ierr);CHKERRQ(ierr)
        call VecGetArrayF90(lambda(imeasurement),vec_ptr2,ierr);CHKERRQ(ierr)
        do local_id = 1, grid%nlmax
          offset = (local_id-1)*ndof
          do i = 1, ndof
            tempreal = 0.d0
            do j = 1, ndof
              tempreal = tempreal + &
                ! this is a transpose block matmult: i,j -> j,i
                inversion_forward_ts_aux%next%dRes_du_k(j,i,local_id) * &
                vec_ptr2(offset+j)
            enddo
            vec_ptr(offset+i) = tempreal
          enddo
        enddo
        call VecRestoreArrayF90(dReskp1_duk_lambdak,vec_ptr,ierr);CHKERRQ(ierr)
        call VecRestoreArrayF90(lambda(imeasurement), &
                                vec_ptr2,ierr);CHKERRQ(ierr)
      endif
      call VecWAXPY(rhs,-1.d0,dReskp1_duk_lambdak,p_,ierr);CHKERRQ(ierr)
      call KSPSolveTranspose(solver%ksp,rhs,lambda(imeasurement), &
                            ierr);CHKERRQ(ierr)

      ! forward loop (adding sensitivity contributions)
      if (this%qoi_is_full_vector) then
        call MatMultTranspose(inversion_forward_ts_aux%dResdparam, &
                              lambda(imeasurement), &
                              dResdKLambda,ierr);CHKERRQ(ierr)
        call VecGetArrayF90(dResdKLambda,vec_ptr,ierr);CHKERRQ(ierr)
        do iparameter = 1, grid%nlmax
          natural_id_ = grid%nG2A(grid%nL2G(iparameter))
          offset = (iparameter-1)*option%nflowdof
          call MatSetValue(this%inversion_aux%JsensitivityT,natural_id_-1, &
                           imeasurement-1,vec_ptr(offset+1),ADD_VALUES, &
                           ierr);CHKERRQ(ierr)
        enddo
        call VecRestoreArrayF90(dResdKLambda,vec_ptr,ierr);CHKERRQ(ierr)
      else
        do iparameter = 1, size(this%parameters)
          call VecZeroEntries(ndof_vec1,ierr);CHKERRQ(ierr)
          call VecGetArrayF90(ndof_vec1,vec_ptr,ierr);CHKERRQ(ierr)
          do local_id = 1, grid%nlmax
            if (patch%imat(grid%nL2G(local_id)) == &
                this%parameters(iparameter)%imat) then
              offset = (local_id-1)*option%nflowdof
              vec_ptr(offset+1:offset+option%nflowdof) = 1.d0
            endif
          enddo
          call VecRestoreArrayF90(ndof_vec1,vec_ptr,ierr);CHKERRQ(ierr)
          call MatMult(inversion_forward_ts_aux%dResdparam,ndof_vec1, &
                       ndof_vec2,ierr);CHKERRQ(ierr)
          call VecDot(ndof_vec2,lambda(imeasurement), &
                      tempreal,ierr);CHKERRQ(ierr)
          if (option%comm%myrank == option%driver%io_rank) then
            call MatSetValue(this%inversion_aux%JsensitivityT,iparameter-1, &
                             imeasurement-1,tempreal,ADD_VALUES, &
                             ierr);CHKERRQ(ierr)
          endif
        enddo
      endif
    enddo
    inversion_forward_ts_aux => inversion_forward_ts_aux%prev
  enddo

  call MatAssemblyBegin(this%inversion_aux%JsensitivityT,MAT_FINAL_ASSEMBLY, &
                        ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(this%inversion_aux%JsensitivityT,MAT_FINAL_ASSEMBLY, &
                      ierr);CHKERRQ(ierr)

  call VecDestroyVecs(size(lambda),lambda,ierr);CHKERRQ(ierr)
  call VecDestroy(onedof_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(p_,ierr);CHKERRQ(ierr)
  call VecDestroy(rhs,ierr);CHKERRQ(ierr)
  call VecDestroy(dReskp1_duk_lambdak,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)

  if (dResdKLambda /= PETSC_NULL_VEC) then
    call VecDestroy(dResdKLambda,ierr);CHKERRQ(ierr)
  endif
  if (ndof_vec1 /= PETSC_NULL_VEC) then
    call VecDestroy(ndof_vec1,ierr);CHKERRQ(ierr)
  endif
  if (ndof_vec2 /= PETSC_NULL_VEC) then
    call VecDestroy(ndof_vec2,ierr);CHKERRQ(ierr)
  endif

  call timer%Stop()
  option%io_buffer = '    ' // &
    trim(StringWrite('(f20.1)',timer%GetCumulativeTime()))
  option%io_buffer = trim(option%io_buffer) // &
    ' seconds to calculate lambdas and sensitivities.'
  call PrintMsg(option)
  call TimerDestroy(timer)

end subroutine InvSubsurfAdjAddSensitivities

! ************************************************************************** !

subroutine InvSubsurfPertCalcSensitivity(this)
  !
  ! Calculates sensitivity matrix Jsensitivity using perturbation
  !
  ! Author: Glenn Hammond
  ! Date: 03/21/22
  !
  use Discretization_module
  use Material_module
  use Option_module

  class(inversion_subsurface_type) :: this

  type(option_type), pointer :: option
  PetscInt :: iteration
  PetscInt :: iqoi(2)
  PetscInt :: i
  PetscErrorCode :: ierr

  ! destroy non-perturbed forward run
  iteration = 0
  ! InvSubsurfPerturbationFillRow performs setup on iteration 0
  if (this%inversion_option%coupled_flow_ert) then
    call InvSubsurfFVCalcPartialJs(this,iteration)
  else
    call InvSubsurfPerturbationFillRow(this,iteration)
  endif
  call this%DestroyForwardRun()
  iteration = 1
  do
    if (associated(this%perturbation%select_cells)) then
      this%perturbation%idof_pert = this%perturbation%select_cells(iteration)
    else
      this%perturbation%idof_pert = iteration
    endif
    call this%InitializeForwardRun(option)
    call InvSubsurfSetupForwardRunLinkage(this) ! do not call mapped version
    call this%ConnectToForwardRun()
    call this%ExecuteForwardRun()
    if (this%inversion_option%coupled_flow_ert) then
      call InvSubsurfFVCalcPartialJs(this,iteration)
    else
      call InvSubsurfPerturbationFillRow(this,iteration)
    endif
    iteration = iteration + 1
    if (iteration > this%perturbation%ndof) exit
    ! the last forward run will be destroyed after any output of
    ! sensitivity matrices
    call this%DestroyForwardRun()
  enddo

  ! Build coupled J once partial Js are calculated
  if (this%inversion_option%coupled_flow_ert) then
    ! destroy last from loop above (we are not calculating Jsense)
    call this%DestroyForwardRun()
    ! -1 is a non-perturbed forward run after perturbation is complete
    this%perturbation%idof_pert = -1
    call this%InitializeForwardRun(option)
    call InvSubsurfSetupForwardRunLinkage(this) ! do not call mapped version
    call this%ConnectToForwardRun()
    call this%ExecuteForwardRun()
    ! the last forward run will be destroyed after any output of
    ! sensitivity matrices
    ! Finalize assembly of coupled sensitivity matrix
    call MatAssemblyBegin(this%inversion_aux%JsensitivityT, &
                          MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(this%inversion_aux%JsensitivityT, &
                        MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  endif

  ! must reset dof back to zero
  this%perturbation%idof_pert = 0

  ! reset measurement vectors to the base model
  if (this%perturbation%base_measurement_vec /= PETSC_NULL_VEC) then
    call VecCopy(this%perturbation%base_measurement_vec, &
                 this%measurement_vec,ierr);CHKERRQ(ierr)
  endif
  call InvSubsurfScatMeasToDistMeas(this, &
                                    this%measurement_vec, &
                                    this%dist_measurement_vec, &
                                    INVSUBSCATFORWARD)

  ! reset parameters to base copy
  if (this%qoi_is_full_vector) then
    iqoi = InversionParameterIntToQOIArray(this%parameters(1))
    call VecCopy(this%perturbation%base_parameter_vec,this%dist_parameter_vec, &
                 ierr);CHKERRQ(ierr)
    call InvSubsurfScatGlobalToDistParam(this, &
                                         this%realization%field%work, &
                                         this%dist_parameter_vec, &
                                         INVSUBSCATREVERSE)
    call DiscretizationGlobalToLocal(this%realization%discretization, &
                                     this%realization%field%work, &
                                     this%realization%field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 iqoi(1),iqoi(2))
  else
    call InvSubsurfCopyParameterToFromVec(this,COPY_FROM_VEC)
    do i = 1, size(this%parameters)
      call InvSubsurfCopyParameterValue(this,i,OVERWRITE_MATERIAL_VALUE)
    enddo
  endif

end subroutine InvSubsurfPertCalcSensitivity

! ************************************************************************** !

subroutine InvSubsurfPerturbationFillRow(this,iteration)
  !
  ! Fills a row (actually column since we store the transpose) of the
  ! Jacobian created through perurbation
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  use Debug_module
  use Realization_Base_class
  use String_module

  class(inversion_subsurface_type) :: this
  PetscInt :: iteration

  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: i
  PetscErrorCode :: ierr

  if (this%perturbation%idof_pert < 0) then
    call this%driver%PrintErrMsg('InvSubsurfPerturbationFillRow() called  &
      &with an idof_pert < 0')
  endif

  call VecGetArrayF90(this%measurement_vec,vec_ptr,ierr);CHKERRQ(ierr)
  do i = 1, size(this%measurements)
    vec_ptr(i) = this%measurements(i)%simulated_value
  enddo
  call VecRestoreArrayF90(this%measurement_vec,vec_ptr,ierr);CHKERRQ(ierr)

  if (this%perturbation%idof_pert == 0) then
    call VecCopy(this%measurement_vec,this%perturbation%base_measurement_vec, &
                 ierr);CHKERRQ(ierr)
    call MatZeroEntries(this%inversion_aux%JsensitivityT,ierr);CHKERRQ(ierr)
  else
    call VecAXPY(this%measurement_vec,-1.d0, &
                 this%perturbation%base_measurement_vec,ierr);CHKERRQ(ierr)
    call VecScale(this%measurement_vec,1.d0/this%perturbation%pert, &
                  ierr);CHKERRQ(ierr)
  endif

  if (this%perturbation%idof_pert == 0) return

  ! don't need to use the distributed vec, but why not
  call InvSubsurfScatMeasToDistMeas(this, &
                                    this%measurement_vec, &
                                    this%dist_measurement_vec, &
                                    INVSUBSCATFORWARD)
  call VecGetArrayF90(this%dist_measurement_vec,vec_ptr,ierr);CHKERRQ(ierr)
  do i = 1, size(vec_ptr)
    call MatSetValue(this%inversion_aux%JsensitivityT, &
                     this%perturbation%idof_pert-1, &
                     this%dist_measurement_offset+i-1,vec_ptr(i), &
                     INSERT_VALUES,ierr);CHKERRQ(ierr)
  enddo
  call VecRestoreArrayF90(this%dist_measurement_vec,vec_ptr, &
                          ierr);CHKERRQ(ierr)

  if (iteration == this%perturbation%ndof) then
    call MatAssemblyBegin(this%inversion_aux%JsensitivityT,MAT_FINAL_ASSEMBLY, &
                          ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(this%inversion_aux%JsensitivityT,MAT_FINAL_ASSEMBLY, &
                        ierr);CHKERRQ(ierr)
  endif

  if (.not.this%qoi_is_full_vector) then
    ! revert back to base value
    this%parameters(this%perturbation%idof_pert)%value = &
      this%perturbation%base_value
  endif

end subroutine InvSubsurfPerturbationFillRow

! ************************************************************************** !

subroutine InvSubsurfFVCalcPartialJs(this,iteration)
  !
  ! Fills entries in partial Jacobian vectors for full-vector inversion
  !
  ! Author: Glenn Hammond
  ! Date: 09/23/22
  !
  use Realization_Base_class

  class(inversion_subsurface_type) :: this
  PetscInt :: iteration

  PetscInt :: i

  if (this%perturbation%idof_pert < 0) then
    call this%driver%PrintErrMsg('InvSubsurfFVCalcPartialJs() called  &
      &with an idof_pert < 0')
  endif

  do i = 1, size(this%inversion_coupled_aux%solutions)
    call InvCoupledUpdateSolnVecs(this%perturbation%idof_pert, &
      this%inversion_coupled_aux%solutions(i)%perturbed_saturation_solution, &
      this%inversion_coupled_aux%solutions(i)%original_saturation_solution, &
      this%inversion_coupled_aux%solutions(i)%dsaturation_dparameter, &
      this%perturbation%pert)
    if (this%inversion_coupled_aux%solutions(i)%perturbed_solute_solution /= &
        PETSC_NULL_VEC) then
      call InvCoupledUpdateSolnVecs(this%perturbation%idof_pert, &
        this%inversion_coupled_aux%solutions(i)%perturbed_solute_solution, &
        this%inversion_coupled_aux%solutions(i)%original_solute_solution, &
        this%inversion_coupled_aux%solutions(i)%dsolute_dparameter, &
        this%perturbation%pert)
    endif
  enddo
  call InversionCoupledAuxReset(this%inversion_coupled_aux)

  !if (iteration == this%perturbation%ndof) then
  !  do i = 1, size(this%inversion_coupled_aux%solutions)
  !    do iparam = 1, size(this%parameters)
  !      print *, i, iparam, 'saturation '
  !      call VecView(this%inversion_coupled_aux%solutions(i)% &
  !                     dsaturation_dparameter(iparam), &
  !                   PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
  !      print *, i, iparam, 'solute '
  !      call VecView(this%inversion_coupled_aux%solutions(i)% &
  !                     dsolute_dparameter(iparam), &
  !                   PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
  !    enddo
  !  enddo
  !endif

  if (.not.this%qoi_is_full_vector .and. &
      this%perturbation%idof_pert > 0) then
    ! revert back to base value
    this%parameters(this%perturbation%idof_pert)%value = &
      this%perturbation%base_value
  endif

end subroutine InvSubsurfFVCalcPartialJs

! ************************************************************************** !

subroutine InvSubsurfScaleSensitivity(this)
  !
  ! Scales sensitivity Jacobian for ln(K)
  !
  ! Author: Glenn Hammond
  ! Date: 10/11/21
  !
  class(inversion_subsurface_type) :: this

  call this%driver%PrintErrMsg('InvSubsurfScaleSensitivity &
    &must be extended from the Subsurface implementation.')

end subroutine InvSubsurfScaleSensitivity

! ************************************************************************** !

subroutine InvSubsurfOutputSensitivity(this,suffix)
  !
  ! Writes sensitivity Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 10/11/21
  !
  use String_module
  class(inversion_subsurface_type) :: this
  character(len=*) :: suffix

  character(len=MAXSTRINGLENGTH) :: filename_prefix

  if (.not.this%print_sensitivity_jacobian) return

  filename_prefix = trim(this%driver%global_prefix) // '_Jsense'
  if (associated(this%perturbation)) then
    filename_prefix = trim(filename_prefix) // '_num'
  else
    filename_prefix = trim(filename_prefix) // '_anal'
  endif
  if (this%annotate_output) then
    filename_prefix = trim(filename_prefix) // '_i' // &
      StringWrite(this%iteration)
  endif
  if (len_trim(suffix) > 0) filename_prefix = trim(filename_prefix) // &
                            '_' // suffix
  call InvSubsurfOutputSensitivityASCII(this,this%inversion_aux% &
                                                    JsensitivityT, &
                                        filename_prefix)
  if (this%qoi_is_full_vector) then
    call InvSubsurfOutputSensitivityHDF5(this, &
                                         this%inversion_aux%JsensitivityT, &
                                         filename_prefix)
  endif

end subroutine InvSubsurfOutputSensitivity

! ************************************************************************** !

subroutine InvSubsurfOutputSensitivityASCII(this,JsensitivityT,filename_prefix)
  !
  ! Writes sensitivity Jacobian to an ASCII output file
  !
  ! Author: Glenn Hammond
  ! Date: 10/11/21
  !
  use Realization_Subsurface_class

  class(inversion_subsurface_type) :: this
  Mat :: JsensitivityT
  character(len=*) :: filename_prefix

  character(len=MAXSTRINGLENGTH) :: string
  PetscViewer :: viewer
  PetscErrorCode :: ierr

  if (.not.associated(this%realization)) then
    call this%driver%PrintErrMsg('InvSubsurfOutputSensitivityASCII must be &
           &called before the forward simulation is destroyed.')
  endif

  string = trim(filename_prefix) // '.txt'
  call PetscViewerASCIIOpen(this%realization%option%mycomm,string,viewer, &
                            ierr);CHKERRQ(ierr)
  call MatView(JsensitivityT,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

end subroutine InvSubsurfOutputSensitivityASCII

! ************************************************************************** !

subroutine InvSubsurfOutputSensitivityHDF5(this,JsensitivityT,filename_prefix)
  !
  ! Writes sensitivity Jacobian to an HDF5 output file
  !
  ! Author: Glenn Hammond
  ! Date: 10/18/21
  !
  use hdf5
  use HDF5_module
  use Output_HDF5_module
  use String_module

  class(inversion_subsurface_type) :: this
  Mat :: JsensitivityT
  character(len=*) :: filename_prefix

  Vec :: row_vec
  PetscMPIInt, parameter :: ON=1, OFF=0
  PetscInt :: imeasurement
  PetscInt :: num_measurement
  PetscErrorCode :: ierr

  integer(HID_T) :: file_id
  integer(HID_T) :: grp_id
  integer(HID_T) :: prop_id
  character(len=MAXSTRINGLENGTH) :: string
  integer :: hdf5_err

  if (.not.associated(this%realization)) then
    call this%driver%PrintErrMsg('InvSubsurfOutputSensitivityHDF5 must be &
           &called before the forward simulation is destroyed.')
  endif

  !HDF5 formatted output
  call h5pcreate_f(H5P_FILE_ACCESS_F,prop_id,hdf5_err)
#ifndef SERIAL_HDF5
  call h5pset_fapl_mpio_f(prop_id,this%realization%option%mycomm, &
                          MPI_INFO_NULL,hdf5_err)
#endif
  string = trim(filename_prefix) // '.h5'
  call h5fcreate_f(string,H5F_ACC_TRUNC_F,file_id,hdf5_err, &
                   H5P_DEFAULT_F,prop_id)
  call h5pclose_f(prop_id,hdf5_err)

  call OutputHDF5WriteStructCoordGroup(file_id, &
                                       this%realization%discretization, &
                                       this%realization%patch%grid, &
                                       this%realization%option)
  ! create a group for the data set
  this%realization%option%time = 0.d0
  write(string,'(''Time:'',es13.5,x,a1)') &
        this%realization%option%time/this%realization%output_option%tconv, &
        this%realization%output_option%tunit
  call h5eset_auto_f(OFF,hdf5_err)
  call h5gopen_f(file_id,string,grp_id,hdf5_err)
  if (hdf5_err /= 0) then
    call h5gcreate_f(file_id,string,grp_id,hdf5_err,OBJECT_NAMELEN_DEFAULT_F)
  endif
  call h5eset_auto_f(ON,hdf5_err)

  num_measurement = size(this%measurements)
  call VecDuplicate(this%dist_measurement_vec,row_vec,ierr);CHKERRQ(ierr)
  do imeasurement = 1, num_measurement
    call VecZeroEntries(row_vec,ierr);CHKERRQ(ierr)
    if (this%realization%option%myrank == 0) then
      call VecSetValue(row_vec,imeasurement-1,1.d0,INSERT_VALUES, &
                       ierr);CHKERRQ(ierr)
    endif
    call VecAssemblyBegin(row_vec,ierr);CHKERRQ(ierr)
    call VecAssemblyEnd(row_vec,ierr);CHKERRQ(ierr)
    call MatMult(JsensitivityT,row_vec,this%realization%field%work, &
                 ierr);CHKERRQ(ierr)
    string = 'Measurement ' // &
      trim(StringWrite(this%measurements(imeasurement)%cell_id))
    call HDF5WriteStructDataSetFromVec(string,this%realization, &
                                       this%realization%field%work,grp_id, &
                                       H5T_NATIVE_DOUBLE)
  enddo
  call VecDestroy(row_vec,ierr);CHKERRQ(ierr)
  call h5gclose_f(grp_id,hdf5_err)
  call OutputHDF5CloseFile(this%realization%option,file_id)

end subroutine InvSubsurfOutputSensitivityHDF5

! ************************************************************************** !

subroutine InvSubsurfScatGlobalToDistParam(this,global_,dist_parameter_vec, &
                                           direction)
  !
  ! Scatters from work to dist_parameter_vec
  !
  ! Author: Glenn Hammond
  ! Date: 04/01/22
  !
  class(inversion_subsurface_type) :: this
  Vec :: global_
  Vec :: dist_parameter_vec
  PetscInt :: direction

  PetscErrorCode :: ierr

  if (direction == INVSUBSCATFORWARD) then
    call VecScatterBegin(this%scatter_global_to_dist_param,global_, &
                         dist_parameter_vec,INSERT_VALUES,SCATTER_FORWARD, &
                         ierr);CHKERRQ(ierr)
    call VecScatterEnd(this%scatter_global_to_dist_param,global_, &
                       dist_parameter_vec,INSERT_VALUES,SCATTER_FORWARD, &
                       ierr);CHKERRQ(ierr)
  else ! INVSUBSCATREVERSE
    call VecScatterBegin(this%scatter_global_to_dist_param,dist_parameter_vec, &
                         global_,INSERT_VALUES,SCATTER_REVERSE, &
                         ierr);CHKERRQ(ierr)
    call VecScatterEnd(this%scatter_global_to_dist_param,dist_parameter_vec, &
                       global_,INSERT_VALUES,SCATTER_REVERSE, &
                       ierr);CHKERRQ(ierr)
  endif

end subroutine InvSubsurfScatGlobalToDistParam

! ************************************************************************** !

subroutine InvSubsurfScatParamToDistParam(this,parameter_vec, &
                                          dist_parameter_vec,direction)
  !
  ! Scatters from parameter_vec to dist_parameter_vec
  !
  ! Author: Glenn Hammond
  ! Date: 04/11/22
  !
  class(inversion_subsurface_type) :: this
  Vec :: parameter_vec
  Vec :: dist_parameter_vec
  PetscInt :: direction

  PetscErrorCode :: ierr

  if (direction == INVSUBSCATFORWARD) then
    ! the parameter_vec is full on each process
    call VecScatterBegin(this%scatter_param_to_dist_param,parameter_vec, &
                         dist_parameter_vec,INSERT_VALUES, &
                         SCATTER_FORWARD_LOCAL,ierr);CHKERRQ(ierr)
    call VecScatterEnd(this%scatter_param_to_dist_param,parameter_vec, &
                       dist_parameter_vec,INSERT_VALUES,SCATTER_FORWARD_LOCAL, &
                       ierr);CHKERRQ(ierr)
  else ! INVSUBSCATREVERSE
    call VecScatterBegin(this%scatter_param_to_dist_param,dist_parameter_vec, &
                         parameter_vec,INSERT_VALUES,SCATTER_REVERSE, &
                         ierr);CHKERRQ(ierr)
    call VecScatterEnd(this%scatter_param_to_dist_param,dist_parameter_vec, &
                       parameter_vec,INSERT_VALUES,SCATTER_REVERSE, &
                       ierr);CHKERRQ(ierr)
  endif

end subroutine InvSubsurfScatParamToDistParam

! ************************************************************************** !

subroutine InvSubsurfScatMeasToDistMeas(this,measurement_vec, &
                                        dist_measurement_vec,direction)
  !
  ! Scatters from measurement_vec to dist_measurement_vec
  !
  ! Author: Glenn Hammond
  ! Date: 04/01/22
  !
  class(inversion_subsurface_type) :: this
  Vec :: measurement_vec
  Vec :: dist_measurement_vec
  PetscInt :: direction

  PetscErrorCode :: ierr

  if (direction == INVSUBSCATFORWARD) then
    ! the measurement_vec is full on each process
    call VecScatterBegin(this%scatter_measure_to_dist_measure,measurement_vec, &
                         dist_measurement_vec,INSERT_VALUES, &
                         SCATTER_FORWARD_LOCAL,ierr);CHKERRQ(ierr)
    call VecScatterEnd(this%scatter_measure_to_dist_measure,measurement_vec, &
                       dist_measurement_vec,INSERT_VALUES, &
                       SCATTER_FORWARD_LOCAL,ierr);CHKERRQ(ierr)
  else ! INVSUBSCATREVERSE
    call VecScatterBegin(this%scatter_measure_to_dist_measure, &
                         dist_measurement_vec,measurement_vec,INSERT_VALUES, &
                         SCATTER_REVERSE,ierr);CHKERRQ(ierr)
    call VecScatterEnd(this%scatter_measure_to_dist_measure, &
                       dist_measurement_vec,measurement_vec,INSERT_VALUES, &
                       SCATTER_REVERSE,ierr);CHKERRQ(ierr)
  endif

end subroutine InvSubsurfScatMeasToDistMeas

! ************************************************************************** !

subroutine InvSubsurfPrintCurMeasValues(this)
  !
  ! Prints the current values of parameters in the list of parameters
  ! Doesn't work for full vector parameter sets
  !
  ! Author: Glenn Hammond
  ! Date: 11/17/22
  !
  use Driver_class

  class(inversion_subsurface_type) :: this

  PetscInt :: i
  PetscInt :: num_measurements

  num_measurements = size(this%measurements)
  if (this%driver%PrintToScreen()) then
    do i = 1, num_measurements
      call InvMeasurePrintComparison(STDOUT_UNIT, &
                                     this%measurements(i),i==1, &
                                     i==num_measurements, &
                                     this%realization%option)
    enddo
  endif
  if (this%driver%PrintToFile()) then
    do i = 1, num_measurements
      call InvMeasurePrintComparison(this%driver%fid_out, &
                                     this%measurements(i),i==1, &
                                     i==num_measurements, &
                                     this%realization%option)
    enddo
  endif

end subroutine InvSubsurfPrintCurMeasValues

! ************************************************************************** !

subroutine InvSubsurfPrintCurParamValues(this)
  !
  ! Prints the current values of parameters in the list of parameters
  ! Doesn't work for full vector parameter sets
  !
  ! Author: Glenn Hammond
  ! Date: 11/17/22
  !
  use Driver_class

  class(inversion_subsurface_type) :: this

  PetscInt :: i
  PetscInt :: num_parameters

  if (this%qoi_is_full_vector) return

  num_parameters = size(this%parameters)
  if (this%driver%PrintToScreen()) then
    do i = 1, num_parameters
      call InversionParameterPrint(STDOUT_UNIT, &
                                   this%parameters(i),i==1, &
                                   i==num_parameters, &
                                   this%realization%option)
    enddo
  endif
  if (this%driver%PrintToFile()) then
    do i = 1, num_parameters
      call InversionParameterPrint(this%driver%fid_out, &
                                   this%parameters(i),i==1, &
                                   i==num_parameters, &
                                   this%realization%option)
    enddo
  endif

end subroutine InvSubsurfPrintCurParamValues

! ************************************************************************** !

subroutine InvSubsurfDestroyForwardRun(this)
  !
  ! Destroys the forward simulation
  !
  ! Author: Glenn Hammond
  ! Date: 03/21/22
  !
  class(inversion_subsurface_type) :: this

  nullify(this%realization)
  call this%forward_simulation%FinalizeRun()
  call this%forward_simulation%Strip()
  deallocate(this%forward_simulation)
  nullify(this%forward_simulation)

end subroutine InvSubsurfDestroyForwardRun

! ************************************************************************** !

subroutine InvSubsurfPerturbationStrip(perturbation)
  !
  ! Deallocates members of inversion perturbation
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  use Utility_module

  type(perturbation_type), pointer :: perturbation

  PetscErrorCode :: ierr

  if (.not.associated(perturbation)) return

  call DeallocateArray(perturbation%select_cells)
  if (perturbation%base_parameter_vec /= PETSC_NULL_VEC) then
    call VecDestroy(perturbation%base_parameter_vec,ierr);CHKERRQ(ierr)
  endif
  if (perturbation%base_measurement_vec /= PETSC_NULL_VEC) then
    call VecDestroy(perturbation%base_measurement_vec,ierr);CHKERRQ(ierr)
  endif
  deallocate(perturbation)
  nullify(perturbation)

end subroutine InvSubsurfPerturbationStrip

! ************************************************************************** !

subroutine InversionSubsurfaceStrip(this)
  !
  ! Deallocates members of inversion Subsurface
  !
  ! Author: Glenn Hammond
  ! Date: 09/20/21
  !
  use Utility_module

  class(inversion_subsurface_type) :: this

  PetscInt :: i
  PetscErrorCode :: ierr

  call InversionBaseStrip(this)
  call InvSubsurfPerturbationStrip(this%perturbation)

  nullify(this%realization)
  call InversionAuxDestroy(this%inversion_aux)
  if (associated(this%forward_simulation)) then
    print *, 'Why is forward simulation still associated in &
             &InversionSubSurfStrip?'
    stop
  endif
  nullify(this%forward_simulation)

  call DeallocateArray(this%local_measurement_values)
  call DeallocateArray(this%local_dobs_dunknown_values)
  call DeallocateArray(this%local_measurement_map)

  if (this%quantity_of_interest /= PETSC_NULL_VEC) then
    call VecDestroy(this%quantity_of_interest,ierr);CHKERRQ(ierr)
  endif
  if (this%ref_quantity_of_interest /= PETSC_NULL_VEC) then
    call VecDestroy(this%ref_quantity_of_interest,ierr);CHKERRQ(ierr)
  endif
  if (this%parameter_vec /= PETSC_NULL_VEC) then
    call VecDestroy(this%parameter_vec,ierr);CHKERRQ(ierr)
  endif
  if (this%dist_parameter_vec /= PETSC_NULL_VEC) then
    call VecDestroy(this%dist_parameter_vec,ierr);CHKERRQ(ierr)
  endif
  if (this%ref_parameter_vec /= PETSC_NULL_VEC) then
    call VecDestroy(this%ref_parameter_vec,ierr);CHKERRQ(ierr)
  endif
  if (this%measurement_vec /= PETSC_NULL_VEC) then
    call VecDestroy(this%measurement_vec,ierr);CHKERRQ(ierr)
  endif
  if (this%dist_measurement_vec /= PETSC_NULL_VEC) then
    call VecDestroy(this%dist_measurement_vec,ierr);CHKERRQ(ierr)
  endif
  if (this%scatter_measure_to_dist_measure /= PETSC_NULL_VECSCATTER) then
    call VecScatterDestroy(this%scatter_measure_to_dist_measure, &
                           ierr);CHKERRQ(ierr)
  endif
  if (this%scatter_param_to_dist_param /= PETSC_NULL_VECSCATTER) then
    call VecScatterDestroy(this%scatter_measure_to_dist_measure, &
                           ierr);CHKERRQ(ierr)
  endif
  if (this%scatter_global_to_dist_param /= PETSC_NULL_VECSCATTER) then
    call VecScatterDestroy(this%scatter_measure_to_dist_measure, &
                           ierr);CHKERRQ(ierr)
  endif

  if (associated(this%measurements)) then
    do i = 1, size(this%measurements)
      call InversionMeasurementAuxStrip(this%measurements(i))
    enddo
    deallocate(this%measurements)
  endif
  nullify(this%measurements)
  if (associated(this%parameters)) then
    do i = 1, size(this%parameters)
      call InversionParameterStrip(this%parameters(i))
    enddo
    deallocate(this%parameters)
  endif
  nullify(this%parameters)

  if (associated(this%inversion_coupled_aux)) then
    call InversionCoupledAuxDestroy(this%inversion_coupled_aux)
  endif

end subroutine InversionSubsurfaceStrip

end module Inversion_Subsurface_class
