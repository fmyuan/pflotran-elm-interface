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

  type, public, extends(inversion_base_type) :: inversion_subsurface_type
    character(len=MAXSTRINGLENGTH) :: forward_simulation_filename
    character(len=MAXSTRINGLENGTH) :: checkpoint_filename
    character(len=MAXSTRINGLENGTH) :: restart_filename
    PetscInt :: restart_iteration
    class(simulation_subsurface_type), pointer :: forward_simulation
    class(realization_subsurface_type), pointer :: realization
    type(inversion_aux_type), pointer :: inversion_aux
    PetscInt :: dist_measurement_offset
    PetscInt :: parameter_offset  ! needed?
    PetscInt :: num_parameters_local
    PetscInt :: n_qoi_per_cell
    Vec :: quantity_of_interest       ! reserved for inversion_ert
    Vec :: ref_quantity_of_interest   ! reserved for inversion_ert
    character(len=MAXWORDLENGTH) :: ref_qoi_dataset_name
    ! flags reference in objects below the inversion objects should be
    ! stored in option_inversion_type
    PetscBool :: print_sensitivity_jacobian
    PetscBool :: annotate_output
    PetscBool :: perturbation_risk_acknowledged
    PetscBool :: debug_adjoint
    PetscInt :: debug_verbosity
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
    procedure, public :: WriteIterationInfo => InvSubsurfWriteIterationInfoLoc
    procedure, public :: CalculateSensitivity => InvSubsurfCalculateSensitivity
    procedure, public :: ScaleSensitivity => InvSuburfSkipThisOnly
    procedure, public :: CalculateUpdate => InvSuburfSkipThisOnly
    procedure, public :: UpdateRegularizationParameters => &
                           InvSuburfSkipThisOnly
    procedure, public :: Strip => InversionSubsurfaceStrip
  end type inversion_subsurface_type

  public :: InversionSubsurfaceCreate, &
            InversionSubsurfaceInit, &
            InversionSubsurfReadSelectCase, &
            InvSubsurfSetupForwardRunLinkage, &
            InvSubsurfConnectToForwardRun, &
            InvSubsurfOutputSensitivity, &
            InvSubsurfPrintCurParamUpdate, &
            InvSubsurfWriteIterationInfo, &
            InversionSubsurfaceStrip

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
  this%inversion_aux => InversionAuxCreate(driver)

  this%quantity_of_interest = PETSC_NULL_VEC
  this%n_qoi_per_cell = UNINITIALIZED_INTEGER
  this%ref_quantity_of_interest = PETSC_NULL_VEC
  this%ref_qoi_dataset_name = ''
  this%forward_simulation_filename = ''
  this%checkpoint_filename = ''
  this%restart_filename = ''
  this%restart_iteration = UNINITIALIZED_INTEGER
  this%print_sensitivity_jacobian = PETSC_FALSE
  this%debug_adjoint = PETSC_FALSE
  this%debug_verbosity = UNINITIALIZED_INTEGER
  this%dist_measurement_offset = UNINITIALIZED_INTEGER
  this%parameter_offset = UNINITIALIZED_INTEGER
  this%num_parameters_local = UNINITIALIZED_INTEGER
  this%annotate_output = PETSC_FALSE
  this%perturbation_risk_acknowledged = PETSC_FALSE

  nullify(this%local_measurement_values)
  nullify(this%local_dobs_dunknown_values)
  nullify(this%local_dobs_dparam_values)
  nullify(this%local_measurement_map)

  nullify(this%forward_simulation)
  nullify(this%realization)

  ! initialize measurement reporting verbosity
  inv_meas_reporting_verbosity = 1
  call InversionParamInitBounds()

end subroutine InversionSubsurfaceInit

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
  PetscReal :: lower_bound
  PetscReal :: upper_bound
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
    case('CHECKPOINT_FILENAME')
      call InputReadFilename(input,option,this%checkpoint_filename)
      call InputErrorMsg(input,option,keyword,error_string)
    case('RESTART_FILENAME')
      call InputReadFilename(input,option,this%restart_filename)
      call InputErrorMsg(input,option,keyword,error_string)
      call InputReadInt(input,option,i)
      if (input%ierr == 0) then
        this%restart_iteration = i
      endif
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
            ! cannot use InputReadAndConvertUnits here because the
            ! time and time units may be reused for subsequent measurement
            ! yet InputReadAndConvertUnits does not return the time units.
            call InputReadDouble(input,option,measurement_time)
            call InputErrorMsg(input,option,keyword,error_string)
            call InputReadWord(input,option,word,PETSC_TRUE)
            if (input%ierr /= 0) word = 'sec'
            measurement_time_units = trim(word)
            internal_units = 'sec'
            measurement_time = measurement_time * &
                           UnitsConvertToInternal(word,internal_units, &
                                                  trim(error_string)// &
                                                  ',MEASUREMENTS,MEASUREMENT,&
                                                  &TIME',option)
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
      else if (associated(this%inversion_aux%measurements)) then
        option%io_buffer = 'Measurements may only be defined in a single block.'
        call PrintErrMsg(option)
      else
        allocate(this%inversion_aux%measurements(last_measurement%id))
        do i = 1, last_measurement%id
          call InversionMeasurementAuxInit(this%inversion_aux%measurements(i))
        enddo
        last_measurement => first_measurement
        do
          if (.not.associated(last_measurement)) exit
          call InversionMeasurementAuxCopy(last_measurement, &
                         this%inversion_aux%measurements(last_measurement%id))
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
          case('BOUNDS')
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'BOUNDS,PARAMETER_NAME', &
                               error_string)
            call InputReadDouble(input,option,lower_bound)
            call InputErrorMsg(input,option,'BOUNDS,LOWER_BOUND',error_string)
            call InputReadDouble(input,option,upper_bound)
            call InputErrorMsg(input,option,'BOUNDS,UPPER_BOUND',error_string)
            i = InversionParamGetItypeFromName(word,this%driver)
            call InversionParamSetBounds(i,lower_bound,upper_bound)
            cycle ! skip appending to parameter list below
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
      else if (associated(this%inversion_aux%parameters)) then
        option%io_buffer = 'Parameters may only be defined in a single block.'
        call PrintErrMsg(option)
      else
        allocate(this%inversion_aux%parameters(last_parameter%id))
        do i = 1, last_parameter%id
          call InversionParameterInit(this%inversion_aux%parameters(i))
        enddo
        last_parameter => first_parameter
        do
          if (.not.associated(last_parameter)) exit
          call InversionParameterCopy(last_parameter, &
                             this%inversion_aux%parameters(last_parameter%id))
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
    case('NUM_PROCESS_GROUPS')
      call InputReadInt(input,option,this%inversion_option%num_process_groups)
      call InputErrorMsg(input,option,keyword,error_string)
    case('PERTURBATION')
      string = trim(error_string)//keyword
      input%ierr = 0
      this%inversion_option%use_perturbation = PETSC_TRUE
      this%inversion_aux%perturbation => InversionAuxPerturbationCreate()
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
              call InputReadDouble(input,option, &
                                   this%inversion_aux%perturbation%tolerance)
              call InputErrorMsg(input,option,keyword,error_string)
          case('SELECT_CELLS')
            call UtilityReadArray(this%inversion_aux% &
                                    perturbation%select_cells, &
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

  if (this%inversion_option%forcomm%group_id > 1) then
    this%driver%print_flags%print_to_screen = PETSC_FALSE
  endif

  option => OptionCreate()
  call OptionSetDriver(option,this%driver)
  call OptionSetComm(option,this%inversion_option%forcomm)
  call OptionSetInversionOption(option,this%inversion_option)
  this%inversion_option%perturbation_run = PETSC_FALSE
  if (associated(this%inversion_aux%perturbation)) then
    this%inversion_option%perturbation_run = &
      (this%inversion_aux%perturbation%idof_pert /= 0)
  else if (this%inversion_option%num_process_groups > 1) then
    option%io_buffer = 'The splitting of inversion communicators is only &
      &allowed for perturbation inversion runs.'
    call PrintErrMsg(option)
  endif
  if (this%inversion_option%num_process_groups > 1 .and. &
      this%inversion_aux%qoi_is_full_vector) then
    call this%driver%PrintErrMsg('Cannot use multiple comm groups for full &
                                 &vector inversion.')
  endif

  ! this %iteration may be overwritten in InvSubsurfRestartParameters
  call InvSubsurfRestartIteration(this)
  call InvSubsurfInitSetGroupPrefix(this,option)
  call FactoryForwardInitialize(this%forward_simulation, &
                                this%forward_simulation_filename,option)
  this%realization => this%forward_simulation%realization

end subroutine InvSubsurfInitForwardRun

! ************************************************************************** !

subroutine InvSubsurfInitSetGroupPrefix(this,option)
  !
  ! Initializes the forward simulation
  !
  ! Author: Glenn Hammond
  ! Date: 03/21/22
  !
  use Option_module
  use String_module

  class(inversion_subsurface_type) :: this
  type(option_type) :: option

  option%group_prefix = 'Run' // trim(StringWrite(this%iteration))
  this%inversion_option%iteration_prefix = option%group_prefix
  if (associated(this%inversion_aux%perturbation)) then
    if (this%annotate_output) then
      if (this%inversion_aux%perturbation%idof_pert > 0) then
        option%group_prefix = trim(option%group_prefix) // 'P' // &
          StringWrite(this%inversion_aux%perturbation%idof_pert)
      else if (this%inversion_aux%perturbation%idof_pert == 0) then
        option%group_prefix = trim(option%group_prefix) // 'Base'
      else
        option%group_prefix = trim(option%group_prefix) // 'Final'
      endif
    else
      if (this%inversion_aux%perturbation%idof_pert > 0) then
        if (this%inversion_option%num_process_groups > 1) then
          ! prevent collisions in the RunTmp files
          option%group_prefix = 'RunTmp-' // &
            StringWrite(this%inversion_option%forcomm%group_id)
        else
          option%group_prefix = 'RunTmp'
        endif
      endif
    endif
  endif

end subroutine InvSubsurfInitSetGroupPrefix

! ************************************************************************** !

subroutine InvSubsurfSetupForwardRunLinkage(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  use Communicator_Aux_module
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
  use Variables_module, only : PERMEABILITY, POROSITY, &
                  ARCHIE_CEMENTATION_EXPONENT, ARCHIE_SATURATION_EXPONENT, &
                  ARCHIE_TORTUOSITY_CONSTANT
  use Waypoint_module
  use ZFlow_Aux_module

  class(inversion_subsurface_type) :: this

  type(patch_type), pointer :: patch
  type(material_property_type), pointer :: material_property
  type(comm_type), pointer :: invcomm
  type(comm_type), pointer :: forcomm
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: iflag
  PetscInt :: i
  PetscInt :: local_id
  PetscInt :: temp_int
  PetscInt :: icount
  PetscInt :: num_measurements, num_measurements_local
  PetscInt :: num_parameters
  PetscInt :: param_id
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
  this%inversion_aux%material_property_array => patch%material_property_array
  this%inversion_aux%cc_array => patch%characteristic_curves_array
  invcomm => this%inversion_option%invcomm
  forcomm => this%inversion_option%forcomm

  if (this%inversion_aux%measurement_vec == PETSC_NULL_VEC) then
    ! perturbation can be problematic with certain flow/transport configurations
    ! check for these situations here
    if (associated(this%inversion_aux%perturbation) .and. &
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

    num_measurements = 0
    if (associated(this%inversion_aux%measurements)) then
      num_measurements = size(this%inversion_aux%measurements)
    else
      call this%driver%PrintErrMsg('No inversion measurements defined.')
    endif
    if (.not.associated(this%inversion_aux%parameters)) then
      call this%driver%PrintErrMsg('No inversion parameters defined.')
    endif

    num_parameters = 0
    do i = 1, size(this%inversion_aux%parameters)
      call InversionParameterMapNameToInt(this%inversion_aux%parameters(i), &
                                          this%driver)
      if (len_trim(this%inversion_aux%parameters(i)%material_name) > 0) then
        material_property => &
            MaterialPropGetPtrFromArray(this%inversion_aux% &
                                  parameters(i)%material_name, &
                                  this%inversion_aux%material_property_array)
        if (.not.associated(material_property)) then
          call this%driver%PrintErrMsg('Inversion MATERIAL "' // &
              trim(this%inversion_aux%parameters(i)%material_name) // &
              '" not found among MATERIAL_PROPERTIES.')
        endif
        this%inversion_aux%parameters(i)%imat = &
          abs(material_property%internal_id)
        num_parameters = num_parameters + 1
      else
        num_parameters = num_parameters + patch%grid%nmax
      endif
    enddo
    if (.not.associated(this%inversion_aux%perturbation)) then
      do i = 1, size(this%inversion_aux%parameters)
        if (i == 1) then
          temp_int = this%inversion_aux%parameters(i)%itype
        else
          if (temp_int /= this%inversion_aux%parameters(i)%itype) then
            call this%driver%PrintErrMsg('Inversion by multiple &
              &parameters of differing type (e.g. permeability, &
              &porosity) only supported for perturbation.')
          endif
        endif
      enddo
    endif
    if (num_parameters == patch%grid%nmax) then
      this%inversion_aux%qoi_is_full_vector = PETSC_TRUE
      this%num_parameters_local = patch%grid%nlmax*this%n_qoi_per_cell
    else
      this%num_parameters_local = PETSC_DECIDE
    endif

    if (size(this%inversion_aux%parameters) > 1 .and. &
        this%inversion_aux%qoi_is_full_vector) then
      call this%driver%PrintErrMsg('More than one parameter not currently &
                                   &supported for full vector inversion.')
    endif

    if (this%inversion_option%coupled_flow_ert) then
      do i = 1, size(this%inversion_aux%parameters)
        param_id = InversionParamGetItypeFromName(this%inversion_aux% &
                                                     parameters(i)% &
                                                     parameter_name, &
                                                   this%driver)
        select case(param_id)
          case(PERMEABILITY,POROSITY,ARCHIE_CEMENTATION_EXPONENT, &
               ARCHIE_SATURATION_EXPONENT,ARCHIE_TORTUOSITY_CONSTANT)
          case default
            string = 'COUPLED_ZFLOW_ERT does not currently support &
              &inversion for "' // &
              trim(this%inversion_aux%parameters(i)%parameter_name) // &
              '", only PERMEABILITY and POROSITY.'
            call this%driver%PrintErrMsg(string)
        end select
      enddo
      if (.not.(associated(this%forward_simulation% &
                            flow_process_model_coupler) .and. &
                associated(this%forward_simulation% &
                            geop_process_model_coupler))) then
        call this%driver%PrintErrMsg('Coupled ZFLOW and ERT inversion &
          &requires that both ZFLOW and ERT process models are employed.')
      endif
      this%inversion_aux%coupled_aux => InversionCoupledAuxCreate()
    endif

    if (associated(invcomm)) then
      ! JsensitivityT is the transpose of the sensitivity Jacobian
      ! with num measurement columns and num parameter rows
      call MatCreateDense(invcomm%communicator, &
                          this%num_parameters_local,PETSC_DECIDE, &
                          num_parameters,num_measurements, &
                          PETSC_NULL_SCALAR, &
                          this%inversion_aux%JsensitivityT, &
                          ierr);CHKERRQ(ierr)
      call MatZeroEntries(this%inversion_aux%JsensitivityT, &
                          ierr);CHKERRQ(ierr)
      ! cannot pass in this%inversion_aux%measurement_vec as it is
      ! initialized to PETSC_NULL_VEC and MatCreateVecs keys off that input
      call MatCreateVecs(this%inversion_aux%JsensitivityT,v,v2, &
                         ierr);CHKERRQ(ierr)
      this%inversion_aux%dist_measurement_vec = v
      this%inversion_aux%dist_parameter_vec = v2
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
      call MPI_Exscan(int_array,int_array2,TWO_INTEGER_MPI,MPIU_INTEGER, &
                      MPI_SUM,invcomm%communicator,ierr);CHKERRQ(ierr)
      this%parameter_offset = int_array2(1)
      this%dist_measurement_offset = int_array2(2)
      deallocate(int_array,int_array2)
    else
      ! still need the measurement vec, but all values will reside on
      ! the 0th rank in forcomm
      i = 0
      if (forcomm%rank == 0) i = num_measurements
      call VecCreateMPI(forcomm%communicator,i,PETSC_DETERMINE, &
                        this%inversion_aux%dist_measurement_vec, &
                        ierr);CHKERRQ(ierr)
    endif

    call VecCreateSeq(PETSC_COMM_SELF,num_measurements, &
                      this%inversion_aux%measurement_vec, &
                      ierr);CHKERRQ(ierr)
    if (.not.this%inversion_aux%qoi_is_full_vector) then
      call VecCreateSeq(PETSC_COMM_SELF,num_parameters, &
                        this%inversion_aux%parameter_vec, &
                        ierr);CHKERRQ(ierr)
    endif

    if (associated(invcomm)) then
      if (this%inversion_aux%qoi_is_full_vector) then
        ! is_parameter should mirror the natural vec
        allocate(int_array(this%num_parameters_local))
        do i = 1, this%num_parameters_local
          int_array(i) = i
        enddo
        int_array = this%parameter_offset + int_array - 1
        call DiscretAOApplicationToPetsc(this%realization%discretization, &
                                         int_array)
        call ISCreateGeneral(invcomm%communicator,size(int_array),int_array, &
                             PETSC_COPY_VALUES,is_petsc,ierr);CHKERRQ(ierr)
        deallocate(int_array)
        call VecScatterCreate(this%realization%field%work,is_petsc, &
                          this%inversion_aux%dist_parameter_vec, &
                          PETSC_NULL_IS, &
                          this%inversion_aux%scatter_global_to_dist_param, &
                          ierr);CHKERRQ(ierr)
        call ISDestroy(is_petsc,ierr);CHKERRQ(ierr)
      else
        ! moved outside earlier, outside of conditional
  !      call VecCreateSeq(PETSC_COMM_SELF,num_parameters, &
  !                        this%inversion_aux%parameter_vec, &
  !                        ierr);CHKERRQ(ierr)
        call VecDuplicate(this%inversion_aux%parameter_vec, &
                          this%inversion_aux%del_parameter_vec, &
                          ierr);CHKERRQ(ierr)
        call ISCreateStride(invcomm%communicator,num_parameters,ZERO_INTEGER, &
                            ONE_INTEGER,is_parameter,ierr);CHKERRQ(ierr)
        call VecScatterCreate(this%inversion_aux%parameter_vec,is_parameter, &
                            this%inversion_aux%dist_parameter_vec, &
                            PETSC_NULL_IS, &
                            this%inversion_aux%scatter_param_to_dist_param, &
                              ierr);CHKERRQ(ierr)
        call ISDestroy(is_parameter,ierr)
      endif
    endif

    do i = 1, num_measurements
      iflag = PETSC_FALSE
      ! ensure that all observed variables are being simulated
      select case(this%inversion_aux%measurements(i)%iobs_var)
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
            trim(StringWrite(this%inversion_aux%measurements(i)%iobs_var)))
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
      if (this%inversion_aux%measurements(i)%iobs_var == &
          OBS_ERT_MEASUREMENT) then
        max_int(1) = max(max_int(1),this%inversion_aux%measurements(i)%cell_id)
      else
        if (Initialized(this%inversion_aux%measurements(i)%coordinate%x)) then
          call GridGetLocalIDFromCoordinate(patch%grid, &
                                this%inversion_aux%measurements(i)%coordinate, &
                                this%realization%option,local_id)
          if (Initialized(local_id)) then
            this%inversion_aux%measurements(i)%cell_id = &
              patch%grid%nG2A(patch%grid%nL2G(local_id))
            this%inversion_aux%measurements(i)%local_id = local_id
          endif
        endif
        max_int(2) = max(max_int(2),this%inversion_aux%measurements(i)%cell_id)
        int_array(i) = this%inversion_aux%measurements(i)%cell_id
      endif
    enddo
    ! ensure that cell and ert measurement ids are within bounds
    call MPI_Allreduce(MPI_IN_PLACE,max_int,TWO_INTEGER_MPI,MPIU_INTEGER, &
                       MPI_MAX,forcomm%communicator,ierr);CHKERRQ(ierr)
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
                       forcomm%communicator,ierr);CHKERRQ(ierr)
    i = maxval(int_array)
    do i = 1, num_measurements
      if (int_array(i) > 0) then
        this%inversion_aux%measurements(i)%cell_id = int_array(i)
      endif
      if (Uninitialized(this%inversion_aux%measurements(i)%cell_id)) then
        string = 'Measurement ' // trim(StringWrite(i)) // &
          ' at coordinate (' // &
          trim(StringWrite(this%inversion_aux% &
                             measurements(i)%coordinate%x)) // ',' // &
          trim(StringWrite(this%inversion_aux% &
                             measurements(i)%coordinate%y)) // ',' // &
          trim(StringWrite(this%inversion_aux% &
                             measurements(i)%coordinate%z)) // &
          ') not mapped properly.'
        call this%driver%PrintErrMsg(string)
      endif
    enddo

    int_array = -1
    iflag = PETSC_FALSE ! no ERT measurements
    do i = 1, num_measurements
      if (this%inversion_aux%measurements(i)%iobs_var /= &
          OBS_ERT_MEASUREMENT) then
        int_array(i) = this%inversion_aux%measurements(i)%cell_id
      else
        iflag = PETSC_TRUE
      endif
    enddo
    ! throw error if ert is included as a measurement when using adjoints
    if (iflag .and. .not.associated(this%inversion_aux%perturbation)) then
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
                    MPI_SUM,forcomm%communicator,ierr);CHKERRQ(ierr)
    ! Exscan does not include the temp_int from this rank
    ! local ids are between temp_int and temp_int + grid%nlmax
    ! (0 < id <= grid%nlmax on process 0; it is zero-based)
    icount = 0
    do i = 1, num_measurements
      if (this%inversion_aux%measurements(i)%iobs_var /= &
          OBS_ERT_MEASUREMENT) then
        if (int_array(i) >= temp_int .and. &
            int_array(i) < temp_int+patch%grid%nlmax) then
          ! it is local
          icount = icount + 1
          local_id = int_array(i) - temp_int + 1
          if (Initialized(this%inversion_aux%measurements(i)%local_id) .and. &
              this%inversion_aux%measurements(i)%local_id /= local_id) then
            this%realization%option%io_buffer = 'Error mapping local id &
              &for measurement ' // trim(StringWrite(i))
            call PrintErrMsgByRank(this%realization%option)
          endif
          this%inversion_aux%measurements(i)%local_id = local_id
        endif
      else if (OptionIsIORank(this%realization%option)) then
        icount = icount + 1
        this%inversion_aux%measurements(i)%local_id = &
          this%inversion_aux%measurements(i)%cell_id
      endif
    enddo
    allocate(this%local_measurement_values(icount))
    allocate(this%local_measurement_map(icount))
    if (.not.associated(this%inversion_aux%perturbation)) then
      ! if adjoint, need to allocate array for potential partial derivatives
      allocate(this%local_dobs_dunknown_values(icount))
      allocate(this%local_dobs_dparam_values(icount))
    endif

    this%local_measurement_map = UNINITIALIZED_INTEGER
    icount = 0
    do i = 1, num_measurements
      if (this%inversion_aux%measurements(i)%iobs_var /= &
          OBS_ERT_MEASUREMENT) then
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

    if (associated(invcomm)) then
      ! map measurement vec to distributed measurement vec
      call ISCreateStride(invcomm%communicator,num_measurements,ZERO_INTEGER, &
                          ONE_INTEGER,is_measure,ierr);CHKERRQ(ierr)
      call VecScatterCreate(this%inversion_aux%measurement_vec,is_measure, &
                        this%inversion_aux%dist_measurement_vec, &
                        PETSC_NULL_IS, &
                        this%inversion_aux%scatter_measure_to_dist_measure, &
                        ierr);CHKERRQ(ierr)
      call ISDestroy(is_measure,ierr)
    endif

    this%inversion_aux%measurements => this%inversion_aux%measurements
    this%inversion_aux%measurement_vec = this%inversion_aux%measurement_vec

    if (associated(this%inversion_aux%perturbation)) then
      if (this%inversion_aux%qoi_is_full_vector) then
        if (associated(this%inversion_aux%perturbation%select_cells)) then
          this%inversion_aux%perturbation%ndof = &
            size(this%inversion_aux%perturbation%select_cells)
          if (this%inversion_aux%perturbation%ndof > &
              this%realization%patch%grid%nmax) then
            call this%driver%PrintErrMsg('Number of SELECT_CELLS is larger &
                                        &than the problem size: '// &
              trim(StringWrite(this%inversion_aux%perturbation%ndof))//' '// &
              trim(StringWrite(this%realization%patch%grid%nmax)))
          endif
        else
          this%inversion_aux%perturbation%ndof = &
            this%realization%patch%grid%nmax
        endif
      else
        this%inversion_aux%perturbation%ndof = &
          size(this%inversion_aux%parameters)
      endif
      if (.not.this%inversion_option%coupled_flow_ert) then
        call VecDuplicate(this%inversion_aux%measurement_vec, &
                          this%inversion_aux%perturbation% &
                            base_measurement_vec, &
                          ierr);CHKERRQ(ierr)
      endif
    endif

    this%inversion_aux%local_measurement_values_ptr => &
      this%local_measurement_values
    if (.not.associated(this%inversion_aux%perturbation)) then
      this%inversion_aux%local_dobs_dunknown_values_ptr => &
        this%local_dobs_dunknown_values
      this%inversion_aux%local_dobs_dparam_values_ptr => &
        this%local_dobs_dparam_values
    endif

    ! if permeability is the parameter of interest, ensure that it is
    ! isotropic or specified with a vertical anisotropy ratio
    if (this%inversion_aux%qoi_is_full_vector) then
      if (MaterialAnisotropyExists(patch%material_properties)) then
        call this%driver%PrintErrMsg('Anisotropic permeability is not &
                      &supported for full vector inversion.')
      endif
    else
      do i = 1, size(this%inversion_aux%parameters)
        if (this%inversion_aux%parameters(i)%itype == PERMEABILITY) then
          material_property => &
            this%inversion_aux%material_property_array(this%inversion_aux% &
                                                        parameters(i)%imat)%ptr
          ! if PERM_HORIZONTAL and VERTICAL_ANISOTROPY_RATIO are used,
          ! isotropic permeabilithy will be false, but tempreal will be 1.
          ! this only works with perturbation
          tempreal = material_property%permeability(3,3) / &
                     material_property%permeability(1,1) / &
                     material_property%vertical_anisotropy_ratio
          if (.not.material_property%isotropic_permeability) then
            if (associated(this%inversion_aux%perturbation)) then
              if (.not.(tempreal > 0.999d0 .and. tempreal < 1.001d0)) then
                call this%driver%PrintErrMsg('Anisotropic permeability is &
                  &only allowed for perturbation-based inversion when &
                  &specified with a PERM_HORIZONTAL and &
                  &VERTICAL_ANISOTROPY_RATIO.')
              endif
            else
              call this%driver%PrintErrMsg('Anisotropic permeability is not &
                &supported for adjoint-based inversion.')
            endif
          endif
        endif
      enddo
    endif

  endif

  this%local_measurement_values = UNINITIALIZED_DOUBLE
  if (.not.associated(this%inversion_aux%perturbation)) then
    this%local_dobs_dunknown_values = UNINITIALIZED_DOUBLE
    this%local_dobs_dparam_values = UNINITIALIZED_DOUBLE
  endif

  if (.not.associated(this%inversion_aux%perturbation)) then
    ! set up pointer to M matrix
    this%inversion_aux%M_ptr = &
        this%forward_simulation%flow_process_model_coupler%timestepper%solver%M
    ! create inversion_ts_aux for first time step
    this%inversion_aux%first_forward_ts_aux => &
      InversionForwardTSAuxCreate(this%inversion_aux%M_ptr, &
                                  this%realization%option%nflowdof, &
                                  patch%grid%nlmax)
    this%inversion_aux%last_forward_ts_aux => &
      this%inversion_aux%first_forward_ts_aux
    this%local_dobs_dunknown_values = UNINITIALIZED_DOUBLE
    this%local_dobs_dparam_values = UNINITIALIZED_DOUBLE

  else ! this%inversion_aux%perturbation is associated

    this%inversion_aux%store_adjoint = PETSC_FALSE

    if (this%inversion_option%coupled_flow_ert) then
      ! ndof is always the number of parameters
      this%inversion_aux%perturbation%ndof = &
        size(this%inversion_aux%parameters)
      this%realization%option%inversion%calculate_ert = &
        this%inversion_aux%perturbation%idof_pert <= 0
      this%realization%option%inversion%record_measurements = &
        this%inversion_aux%perturbation%idof_pert <= 0
      this%realization%option%inversion%calculate_ert_jacobian = &
        this%inversion_aux%perturbation%idof_pert < 0
    else ! normal perturbation
      if (this%inversion_aux%qoi_is_full_vector .and. &
          this%inversion_aux%perturbation%idof_pert > &
            this%realization%patch%grid%nmax) then
        call this%driver%PrintErrMsg('SELECT_CELLS ID is larger than &
                                    &the problem size: '// &
          trim(StringWrite(this%inversion_aux%perturbation%idof_pert))//' '// &
          trim(StringWrite(this%realization%patch%grid%nmax)))
      endif
    endif

    if (Uninitialized(this%inversion_aux%parameters(1)%itype) .and. &
        this%inversion_aux%qoi_is_full_vector) then
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
  use Factory_Forward_module
  use Factory_Subsurface_module
  use Init_Subsurface_module
  use Material_module
  use String_module
  use Timestepper_Base_class, only : TS_STOP_END_SIMULATION
  use Utility_module
  use Waypoint_module
  use ZFlow_Aux_module

  class(inversion_subsurface_type) :: this

  PetscReal :: final_time
  PetscInt :: i, sync_count
  type(waypoint_type), pointer :: waypoint
  type(inversion_perturbation_type), pointer :: perturbation
  PetscReal, pointer :: real_array(:)
  PetscReal, pointer :: vec_ptr(:)
  PetscMPIInt :: mpi_int
  PetscBool :: iflag
  PetscInt :: iqoi(2)
  PetscErrorCode :: ierr

  perturbation => this%inversion_aux%perturbation

  ! the allocation of sync_times come before the call to
  ! InvCoupledAllocateSolnVecs()
  if (.not.associated(this%inversion_aux%sync_times)) then
    ! insert measurement times into waypoint list. this must come after the
    ! simulation is initialized to obtain the final time.
    allocate(real_array(size(this%inversion_aux%measurements)))
    final_time = &
      WaypointListGetFinalTime(this%forward_simulation%waypoint_list_subsurface)
    iflag = PETSC_FALSE
    do i = 1, size(this%inversion_aux%measurements)
      if (Uninitialized(this%inversion_aux%measurements(i)%time)) then
        this%inversion_aux%measurements(i)%time = final_time
      endif
      if (this%inversion_aux%measurements(i)%time > final_time) then
        iflag = PETSC_TRUE
      endif
      real_array(i) = this%inversion_aux%measurements(i)%time
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
    allocate(this%inversion_aux%sync_times(sync_count))
    this%inversion_aux%sync_times(:) = real_array(1:sync_count)
    deallocate(real_array)
    nullify(real_array)

    if (associated(this%inversion_aux%coupled_aux)) then
      allocate(this%inversion_aux%coupled_aux%solutions(sync_count))
      do i = 1, sync_count
        call InversionCoupledSolutionInit(this%inversion_aux%coupled_aux% &
                                            solutions(i))
        this%inversion_aux%coupled_aux%solutions(i)%time = &
          this%inversion_aux%sync_times(i)
      enddo
      ! allocate any full vector measurements
      call InvCoupledAllocateSolnVecs(this%inversion_aux%coupled_aux, &
                                      this%realization%field%work, &
                                      size(this%inversion_aux%parameters))
    endif
  endif
  real_array => this%inversion_aux%sync_times
  do i = 1, size(real_array)
    waypoint => WaypointCreate()
    waypoint%time = real_array(i)
    waypoint%sync = PETSC_TRUE
    call WaypointInsertInList(waypoint, &
                              this%forward_simulation%waypoint_list_outer)
  enddo


  ! store the initial set of parameters
  if (associated(this%inversion_option%invcomm) .and. &
      this%inversion_aux%startup_phase) then
    if (this%inversion_aux%qoi_is_full_vector) then
      iqoi = InversionParameterIntToQOIArray( &
                                  this%inversion_aux%parameters(1))
      ! on first iteration of inversion, store the values
      call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                  this%realization%field%work_loc, &
                                  iqoi(1),iqoi(2))
      call DiscretizationLocalToGlobal(this%realization%discretization, &
                                      this%realization%field%work_loc, &
                                      this%realization%field%work,ONEDOF)
      call InvAuxScatGlobalToDistParam(this%inversion_aux, &
                                    this%realization%field%work, &
                                    this%inversion_aux%dist_parameter_vec, &
                                    INVAUX_SCATFORWARD)
    else
      ! load the original parameter values
      call InvAuxInitializeParameterValues(this%inversion_aux)
    endif
    ! must come after the copying of parameters above
    call this%RestartReadData()
  endif
  ! this flag can be set to false after restart has been read
  this%inversion_aux%startup_phase = PETSC_FALSE

  if (associated(perturbation)) then
    ! update parameters for parallel perturbation
    if (this%inversion_option%num_process_groups > 1 .and. &
        this%inversion_aux%perturbation%idof_pert == 0) then
      ! at this point, only the parameter_vec has the most up-to-date
      ! values. do not call InvAuxMaterialToParamVec() as it will
      ! overwrite parameter_vec
      call VecGetArrayF90(this%inversion_aux%parameter_vec,vec_ptr, &
                          ierr);CHKERRQ(ierr)
      mpi_int = size(this%inversion_aux%parameters)
      ! could bcast from driver%comm%rank == 0, but this is more
      ! organized
      if (this%inversion_option%forcomm_i%group_id == 1) then
!print *, this%driver%comm%rank, ' Bcast1 forcomm_i%group_id == 1'
!call MPI_Barrier(this%inversion_option%forcomm_i%communicator,ierr);CHKERRQ(ierr)
        call MPI_Bcast(vec_ptr,mpi_int,MPI_DOUBLE_PRECISION,ZERO_INTEGER_MPI, &
                      this%inversion_option%forcomm_i%communicator, &
                      ierr);CHKERRQ(ierr)
      endif
!call MPI_Barrier(this%inversion_option%forcomm%communicator,ierr);CHKERRQ(ierr)
!print *, this%driver%comm%rank, ' Bcast2 all'
      call MPI_Bcast(vec_ptr,mpi_int,MPI_DOUBLE_PRECISION,ZERO_INTEGER_MPI, &
                    this%inversion_option%forcomm%communicator, &
                    ierr);CHKERRQ(ierr)
      call VecRestoreArrayF90(this%inversion_aux%parameter_vec,vec_ptr, &
                              ierr);CHKERRQ(ierr)
      call InvAuxParamVecToMaterial(this%inversion_aux)
    endif
    ! do not connect non-perturbation runs for processes not in invcomm
    if (.not. associated(this%inversion_option%invcomm) .and. &
        this%inversion_aux%perturbation%idof_pert <= 0) then
      ! skip the execution of the forward run
      this%realization%option%status = SKIP
      this%forward_simulation%stop_flag = TS_STOP_END_SIMULATION
      return
    endif
  endif

  this%realization%patch%aux%inversion_aux => this%inversion_aux
  call InversionAuxResetMeasurements(this%inversion_aux)

  call FactoryForwardPrerequisite(this%forward_simulation)

  ! must come after insertion of waypoints and setting of pointers; otherwise,
  ! the pmc timestepper cur_waypoint pointers are pointing at the time 0
  ! waypoint, instead of the one just after time 0, which is incorrect.
  !
  ! forward simulation parameters are overwritten within InitializeRun()
  call this%forward_simulation%InitializeRun()

  if (.not.associated(perturbation)) then
    ! pass in first parameter as an earlier check prevents adjoint-based
    ! inversion for more than one parameter type
    call InvSubsurfSetAdjointVariable(this,this%inversion_aux% &
                                             parameters(1)%itype)
  endif

end subroutine InvSubsurfConnectToForwardRun

! ************************************************************************** !

subroutine InvSubsurfSetAdjointVariable(this,iparameter_type)
  !
  ! Sets the adjoint variable for a process model
  !
  ! Author: Glenn Hammond
  ! Date: 03/30/22

  use String_module
  use Variables_module, only : PERMEABILITY, POROSITY, VG_ALPHA, VG_M, &
                               VG_SR, ARCHIE_CEMENTATION_EXPONENT, &
                               ARCHIE_SATURATION_EXPONENT, &
                               ARCHIE_TORTUOSITY_CONSTANT
  use ZFlow_Aux_module

  class(inversion_subsurface_type) :: this
  PetscInt :: iparameter_type

  character(len=MAXSTRINGLENGTH) :: string

  select case(this%realization%option%iflowmode)
    case(ZFLOW_MODE)
    case default
      string = 'Flow mode "' // trim(this%realization%option%flowmode) // &
        '" not supported for inversion (InvSubsurfSetAdjointVariable).'
    call this%driver%PrintErrMsg(string)
  end select

  select case(iparameter_type)
    case(PERMEABILITY)
      zflow_adjoint_parameter = ZFLOW_ADJOINT_PERMEABILITY
    case(POROSITY)
      zflow_adjoint_parameter = ZFLOW_ADJOINT_POROSITY
    case(VG_ALPHA,VG_M,VG_SR)
      string = 'van Genuchten parameters are unsupported for adjoint-based &
        &inversion.'
      call this%driver%PrintErrMsg(string)
    case(ARCHIE_CEMENTATION_EXPONENT,ARCHIE_SATURATION_EXPONENT, &
         ARCHIE_TORTUOSITY_CONSTANT)
      string = "Archie's parameters are unsupported for adjoint-based &
        &inversion."
      call this%driver%PrintErrMsg(string)
    case default
      string = 'Unrecognized variable in InvSubsurfSetAdjointVariable: ' // &
               trim(StringWrite(iparameter_type))
      call this%driver%PrintErrMsg(string)
  end select

end subroutine InvSubsurfSetAdjointVariable

! ************************************************************************** !

subroutine InvSubsurfExecuteForwardRun(this)
  !
  ! Executes a forward simulation
  !
  ! Author: Glenn Hammond
  ! Date: 03/21/22

  use Option_module
  use Timestepper_Base_class, only : TS_STOP_END_SIMULATION

  class(inversion_subsurface_type) :: this

  type(option_type), pointer :: option

  option => this%realization%option
  if (option%status == PROCEED) then
    call this%forward_simulation%ExecuteRun()
    if (this%forward_simulation%stop_flag == TS_STOP_END_SIMULATION) then
      call InvSubsurfPostProcMeasurements(this)
    else
      option%io_buffer = 'Inversion forward simulation "' // &
        trim(option%group_prefix) // '" failed.'
      call PrintErrMsg(option)
    endif
  endif

end subroutine InvSubsurfExecuteForwardRun

! ************************************************************************** !

subroutine InvSubsurfWriteIterationInfoLoc(this)
  !
  ! Writes inversion run info
  !
  ! Author: Glenn Hammond
  ! Date: 12/16/22
  !

  use String_module

  implicit none

  class(inversion_subsurface_type) :: this

  character(len=:), allocatable :: string
  character(len=:), allocatable :: nl
  character(len=80) :: divider

  nl = new_line('a')
  write(divider,'(40("=+"))')
  string = nl // trim(divider) // nl
  call this%driver%PrintMsg(string)
  call InvSubsurfWriteIterationInfo(this)
  string = nl // divider // nl
  call this%driver%PrintMsg(string)

end subroutine InvSubsurfWriteIterationInfoLoc

! ************************************************************************** !

subroutine InvSubsurfWriteIterationInfo(this)
  !
  ! Prints information for the current inversion iteration to the screen
  ! and/or output file.
  !
  ! Author: Glenn Hammond
  ! Date: 12/16/22

  class(inversion_subsurface_type) :: this

  call InvSubsurfPrintCurMeasValues(this)
  call InvSubsurfPrintCurParamValues(this)

end subroutine InvSubsurfWriteIterationInfo

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

  if (associated(this%inversion_aux%perturbation)) then
    call InvSubsurfPertCalcSensitivity(this)
  else
    call InvSubsurfAdjointCalcSensitivity(this)
  endif

  if (associated(this%inversion_option%invcomm)) then
    call InvSubsurfOutputSensitivity(this,'')
    if (.not.this%inversion_aux%qoi_is_full_vector) then
      call InvAuxScatParamToDistParam(this%inversion_aux, &
                                      this%inversion_aux%parameter_vec, &
                                      this%inversion_aux%dist_parameter_vec, &
                                      INVAUX_SCATFORWARD)
    endif
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
  PetscInt :: num_measurements
  PetscMPIInt :: mpi_int
  PetscErrorCode :: ierr

  option => this%realization%option

  num_measurements = size(this%inversion_aux%measurements)
  ! distribute measurements to measurement objects
  call VecSet(this%inversion_aux%measurement_vec,UNINITIALIZED_DOUBLE, &
              ierr);CHKERRQ(ierr)
  call VecSet(this%inversion_aux%dist_measurement_vec,-888.d0, &
              ierr);CHKERRQ(ierr)
  icount = 0
  do imeasurement = 1, num_measurements
    if (Initialized(this%inversion_aux%measurements(imeasurement)% &
                      local_id)) then
      icount = icount + 1
      call VecSetValue(this%inversion_aux%dist_measurement_vec, &
                       this%local_measurement_map(icount)-1, &
                       this%local_measurement_values(icount),&
                       INSERT_VALUES,ierr);CHKERRQ(ierr)
!      print *, 'setm: ', this%inversion_option%forcomm%rank, &
!               this%inversion_option%forcomm%group_id, &
!               this%local_measurement_map(icount), &
!               this%local_measurement_values(icount)
    endif
  enddo
  call VecAssemblyBegin(this%inversion_aux%dist_measurement_vec, &
                        ierr);CHKERRQ(ierr)
  call VecAssemblyEnd(this%inversion_aux%dist_measurement_vec, &
                      ierr);CHKERRQ(ierr)

!  print *, 'local_measurement_map: ', &
!    this%inversion_option%forcomm%rank, this%local_measurement_map

  if (associated(this%inversion_option%invcomm)) then
    call InvAuxScatMeasToDistMeas(this%inversion_aux, &
                                  this%inversion_aux%measurement_vec, &
                                  this%inversion_aux%dist_measurement_vec, &
                                  INVAUX_SCATREVERSE)
  else
    ! forward runs without invcomm still need to copy the values into
    ! measurement_vec
    call VecGetArrayF90(this%inversion_aux%measurement_vec,vec_ptr2, &
                        ierr);CHKERRQ(ierr)
    call VecGetArrayF90(this%inversion_aux%dist_measurement_vec,vec_ptr, &
                        ierr);CHKERRQ(ierr)
    if (this%inversion_option%forcomm%rank == 0) then
      vec_ptr2(:) = vec_ptr(:)
    endif
    call VecRestoreArrayF90(this%inversion_aux%dist_measurement_vec, &
                            vec_ptr,ierr);CHKERRQ(ierr)
    ! broadcast to perturbation ranks > 0 in forcomm
    mpi_int = 0
!print *, this%driver%comm%rank, ' Bcast3 group_id > 1'
!call MPI_Barrier(this%inversion_option%forcomm%communicator,ierr);CHKERRQ(ierr)
    call MPI_Bcast(vec_ptr2,num_measurements,MPI_DOUBLE_PRECISION,mpi_int, &
                   this%inversion_option%forcomm%communicator, &
                   ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(this%inversion_aux%measurement_vec, &
                            vec_ptr2,ierr);CHKERRQ(ierr)
  endif
  call InvAuxCopyMeasToFromMeasVec(this%inversion_aux,INVAUX_COPY_FROM_VEC)
  do imeasurement = 1, num_measurements
    this%inversion_aux%measurements(imeasurement)%measured = PETSC_TRUE
  enddo
  if (associated(this%local_dobs_dunknown_values)) then
    ! distribute derivatives to measurement objects
    ! temporary vecs for derivatives (if they exist)
    call VecDuplicate(this%inversion_aux%measurement_vec, &
                      dobs_dunknown_vec,ierr);CHKERRQ(ierr)
    call VecDuplicate(this%inversion_aux%dist_measurement_vec, &
                      dist_dobs_dunknown_vec,ierr);CHKERRQ(ierr)
    call VecDuplicate(this%inversion_aux%measurement_vec, &
                      dobs_dparam_vec,ierr);CHKERRQ(ierr)
    call VecDuplicate(this%inversion_aux%dist_measurement_vec, &
                      dist_dobs_dparam_vec,ierr);CHKERRQ(ierr)
    call VecSet(dobs_dunknown_vec,UNINITIALIZED_DOUBLE,ierr);CHKERRQ(ierr)
    call VecSet(dist_dobs_dunknown_vec,-888.d0,ierr);CHKERRQ(ierr)
    call VecSet(dobs_dparam_vec,UNINITIALIZED_DOUBLE,ierr);CHKERRQ(ierr)
    ! dist_dobs_dparam_vec has to be UNINITIALIZED_DOUBLE, otherwise -888 will
    ! be set for all uninitialized values
    call VecSet(dist_dobs_dparam_vec,UNINITIALIZED_DOUBLE,ierr);CHKERRQ(ierr)
    icount = 0
    do imeasurement = 1, num_measurements
      if (Initialized(this%inversion_aux%measurements(imeasurement)% &
                        local_id)) then
        icount = icount + 1
        ! set the partial derivative
        select case(this%inversion_aux%measurements(imeasurement)%iobs_var)
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
    call InvAuxScatMeasToDistMeas(this%inversion_aux, &
                                  dobs_dunknown_vec, &
                                  dist_dobs_dunknown_vec, &
                                  INVAUX_SCATREVERSE)
    call InvAuxScatMeasToDistMeas(this%inversion_aux, &
                                  dobs_dparam_vec, &
                                  dist_dobs_dparam_vec, &
                                  INVAUX_SCATREVERSE)
    call VecGetArrayF90(dobs_dunknown_vec,vec_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(dobs_dparam_vec,vec_ptr2,ierr);CHKERRQ(ierr)
    do imeasurement = 1, num_measurements
      this%inversion_aux%measurements(imeasurement)%dobs_dunknown = &
        vec_ptr(imeasurement)
      this%inversion_aux%measurements(imeasurement)%dobs_dparam = &
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
  do imeasurement = 1, num_measurements
    call InversionMeasurementPrintConcise(this%inversion_aux% &
                                            measurements(imeasurement), &
                                          'PostProcess',option)
    if (.not.this%inversion_aux%measurements(imeasurement)%measured) then
      option%io_buffer = 'Measurement at cell ' // &
        trim(StringWrite(this%inversion_aux%measurements(imeasurement)% &
                           cell_id))
      if (Initialized(this%inversion_aux%measurements(imeasurement)%time)) then
        word = 'sec'
        option%io_buffer = trim(option%io_buffer) // &
          ' at ' // &
          trim(StringWrite(this%inversion_aux% &
                             measurements(imeasurement)%time* &
          UnitsConvertToExternal(this%inversion_aux% &
                                   measurements(imeasurement)%time_units, &
                                 word,option))) // &
          ' ' // trim(this%inversion_aux%measurements(imeasurement)%time_units)
      endif
      option%io_buffer = trim(option%io_buffer) // &
        ' with a measured value from the measurement file of ' // &
        trim(StringWrite(this%inversion_aux% &
                           measurements(imeasurement)%value)) // &
        ' was not measured during the simulation.'
      call PrintErrMsg(option)
    endif
  enddo

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
  do imeasurement = 1, size(this%inversion_aux%measurements)
    this%inversion_aux%measurements(imeasurement)%first_lambda = PETSC_FALSE
  enddo

  ! go to end of list
  cur_inversion_ts_aux => inversion_aux%last_forward_ts_aux
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
    call InvForwardTSAuxDestroyMatrices(cur_inversion_ts_aux)
    ! point cur_inversion_ts_aux to the end of the list
    inversion_aux%last_forward_ts_aux => prev_inversion_ts_aux
  endif

  call InvSubsurfAdjAddSensitivities(this)

  call InvAuxAdjCleanupAfterForwardRun(inversion_aux)

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
  type(inversion_measurement_aux_type), pointer :: measurements(:)

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
  measurements => this%inversion_aux%measurements

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
  if (this%inversion_aux%qoi_is_full_vector) then
    call VecDuplicate(ndof_vec,dResdKLambda,ierr);CHKERRQ(ierr)
  else
    call VecDuplicate(ndof_vec,ndof_vec1,ierr);CHKERRQ(ierr)
    call VecDuplicate(ndof_vec1,ndof_vec2,ierr);CHKERRQ(ierr)
  endif
  ! must use natural vec, not parameter vec since offset is based on measurement
  call DiscretizationCreateVector(discretization,ONEDOF, &
                                  natural_vec,NATURAL,option)

  call VecDuplicateVecsF90(ndof_vec,size(measurements),lambda, &
                           ierr);CHKERRQ(ierr)
  do i = 1, size(lambda)
    call VecZeroEntries(lambda(i),ierr);CHKERRQ(ierr)
  enddo

  call MatZeroEntries(this%inversion_aux%JsensitivityT,ierr);CHKERRQ(ierr)

  inversion_forward_ts_aux => this%inversion_aux%last_forward_ts_aux
  do
    if (.not.associated(inversion_forward_ts_aux)) exit

    call KSPSetOperators(solver%ksp,inversion_forward_ts_aux%dResdu, &
                        inversion_forward_ts_aux%dResdu,ierr);CHKERRQ(ierr)
    do imeasurement = 1, size(measurements)

      ! skip lambdas after the measurement time
      if (Initialized(measurements(imeasurement)%time) .and. &
          inversion_forward_ts_aux%time > &
            measurements(imeasurement)%time + tol) then
        cycle
      endif

      ! backward loop contribution (calculating lambda for the ts)
      if (.not.measurements(imeasurement)%first_lambda) then
        measurements(imeasurement)%first_lambda = PETSC_TRUE

        ! if we have a dobs_dparam, add that value into Jsens
        if (Initialized(measurements(imeasurement)%dobs_dparam)) then
          if (this%inversion_aux%qoi_is_full_vector .or. &
              size(this%inversion_aux%parameters) > 0) then
            option%io_buffer = 'Need to refactor dobs_dparam for more than &
              &one parameter.'
            call PrintErrMsg(option)
          endif
          if (OptionIsIORank(option)) then
            tempreal = -measurements(imeasurement)%dobs_dparam
            iparameter = 1
            call MatSetValue(this%inversion_aux%JsensitivityT,iparameter-1, &
                            imeasurement-1,tempreal,ADD_VALUES, &
                            ierr);CHKERRQ(ierr)
          endif
        endif

        ! begin lambda calculations
        call VecZeroEntries(natural_vec,ierr);CHKERRQ(ierr)
        if (option%myrank == 0) then
          icell_measurement = measurements(imeasurement)%cell_id
          tempreal = -1.d0
          ! scale by partial derivative if measurement is not the
          ! primary dependent variable
          select case(measurements(imeasurement)%iobs_var)
            case(OBS_LIQUID_SATURATION)
              tempreal = tempreal * &
                        measurements(imeasurement)%dobs_dunknown
          end select
          call VecSetValue(natural_vec,icell_measurement-1,tempreal, &
                          INSERT_VALUES,ierr);CHKERRQ(ierr)
        endif
        call InversionMeasurementPrintConcise( &
               measurements(imeasurement),'',option)
        call VecAssemblyBegin(natural_vec,ierr);CHKERRQ(ierr)
        call VecAssemblyEnd(natural_vec,ierr);CHKERRQ(ierr)
        call DiscretizationNaturalToGlobal(discretization,natural_vec, &
                                          onedof_vec,ONEDOF)
        select case(measurements(imeasurement)%iobs_var)
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
      if (this%inversion_aux%qoi_is_full_vector) then
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
        do iparameter = 1, size(this%inversion_aux%parameters)
          call VecZeroEntries(ndof_vec1,ierr);CHKERRQ(ierr)
          call VecGetArrayF90(ndof_vec1,vec_ptr,ierr);CHKERRQ(ierr)
          do local_id = 1, grid%nlmax
            if (patch%imat(grid%nL2G(local_id)) == &
                this%inversion_aux%parameters(iparameter)%imat) then
              offset = (local_id-1)*option%nflowdof
              vec_ptr(offset+1:offset+option%nflowdof) = 1.d0
            endif
          enddo
          call VecRestoreArrayF90(ndof_vec1,vec_ptr,ierr);CHKERRQ(ierr)
          call MatMult(inversion_forward_ts_aux%dResdparam,ndof_vec1, &
                       ndof_vec2,ierr);CHKERRQ(ierr)
          call VecDot(ndof_vec2,lambda(imeasurement), &
                      tempreal,ierr);CHKERRQ(ierr)
          if (OptionIsIORank(option)) then
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
  PetscInt :: idof_pert
  PetscInt :: iqoi(2)
  PetscErrorCode :: ierr

  ! destroy non-perturbed forward run
  idof_pert = 0
  ! InvSubsurfPerturbationFillRow performs setup on iteration 0
  if (this%inversion_option%coupled_flow_ert) then
    call InvSubsurfFVSetOrigSoln(this)
  else
    call InvSubsurfPerturbationFillRow(this,idof_pert)
  endif
  call this%DestroyForwardRun()
  ! begin by assigning perturbations based on the group id
  idof_pert = this%inversion_option%forcomm%group_id
  do
    if (associated(this%inversion_aux%perturbation%select_cells)) then
      this%inversion_aux%perturbation%idof_pert = &
        this%inversion_aux%perturbation%select_cells(idof_pert)
    else
      this%inversion_aux%perturbation%idof_pert = idof_pert
      ! the following condition allows the perturbation forward
      ! runs to be destroyed below.
!      if (.not. associated(this%inversion_option%invcomm) .and. &
!          idof_pert > this%inversion_aux%perturbation%ndof) exit
    endif
    call this%InitializeForwardRun(option)
    call InvSubsurfSetupForwardRunLinkage(this) ! do not call mapped version
    call this%ConnectToForwardRun()
    call this%ExecuteForwardRun()
    if (this%inversion_option%coupled_flow_ert) then
      call InvSubsurfFVCalcPartialJs(this,idof_pert)
    else
      call InvSubsurfPerturbationFillRow(this,idof_pert)
    endif
    idof_pert = idof_pert + this%inversion_option%num_process_groups
!    if (associated(this%inversion_option%invcomm) .and. &
    if (idof_pert > this%inversion_aux%perturbation%ndof) exit
    ! the last forward run will be destroyed after any output of
    ! sensitivity matrices (invcomm only)
    call this%DestroyForwardRun()
  enddo

  ! Build coupled J once partial Js are calculated
  if (associated(this%inversion_option%invcomm) .and. &
      this%inversion_option%coupled_flow_ert) then
    ! destroy last from loop above (we are not calculating Jsense)
    call this%DestroyForwardRun()
    ! -1 is a non-perturbed forward run after perturbation is complete
    this%inversion_aux%perturbation%idof_pert = -1
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
  this%inversion_aux%perturbation%idof_pert = 0

  ! reset measurement vectors to the base model
  if (this%inversion_aux%perturbation%base_measurement_vec /= &
      PETSC_NULL_VEC) then
    call VecCopy(this%inversion_aux%perturbation%base_measurement_vec, &
                 this%inversion_aux%measurement_vec,ierr);CHKERRQ(ierr)
  endif
  if (associated(this%inversion_option%invcomm)) then
    call InvAuxScatMeasToDistMeas(this%inversion_aux, &
                                  this%inversion_aux%measurement_vec, &
                                  this%inversion_aux%dist_measurement_vec, &
                                  INVAUX_SCATFORWARD)
    call InvAuxCopyMeasToFromMeasVec(this%inversion_aux,INVAUX_COPY_FROM_VEC)
  endif

  ! reset parameters to base copy
  if (this%inversion_aux%qoi_is_full_vector) then
    iqoi = InversionParameterIntToQOIArray(this%inversion_aux%parameters(1))
    call VecCopy(this%inversion_aux%perturbation%base_parameter_vec, &
                 this%inversion_aux%dist_parameter_vec, &
                 ierr);CHKERRQ(ierr)
    call InvAuxScatGlobalToDistParam(this%inversion_aux, &
                                     this%realization%field%work, &
                                     this%inversion_aux%dist_parameter_vec, &
                                     INVAUX_SCATREVERSE)
    call DiscretizationGlobalToLocal(this%realization%discretization, &
                                     this%realization%field%work, &
                                     this%realization%field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 iqoi(1),iqoi(2))
  else
    call InvAuxParamVecToMaterial(this%inversion_aux)
  endif

end subroutine InvSubsurfPertCalcSensitivity

! ************************************************************************** !

subroutine InvSubsurfPerturbationFillRow(this,my_dof)
  !
  ! Fills a row (actually column since we store the transpose) of the
  ! Jacobian created through perurbation
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/21
  !
  use Debug_module
  use Option_Inversion_module
  use Inversion_Measurement_Aux_module
  use Realization_Base_class
  use String_module

  class(inversion_subsurface_type) :: this
  PetscInt :: my_dof

  type(inversion_option_type), pointer :: inversion_option
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: i
  PetscInt :: igroup
  PetscInt :: iparameter
  PetscInt :: num_measurements
  PetscMPIInt :: mpi_int
  PetscMPIInt :: status(MPI_STATUS_SIZE)
  PetscErrorCode :: ierr

  if (my_dof < 0) then
    call this%driver%PrintErrMsg('InvSubsurfPerturbationFillRow() called  &
      &with an idof_pert < 0')
  endif

  inversion_option => this%inversion_option

  num_measurements = size(this%inversion_aux%measurements)

  call InvAuxCopyMeasToFromMeasVec(this%inversion_aux,INVAUX_COPY_TO_VEC)

  if (my_dof == 0) then
    ! if parallel perturbation runs, need to broadcast the base measurement
    if (associated(this%inversion_option%forcomm_i)) then
      call VecGetArrayF90(this%inversion_aux%measurement_vec, &
                          vec_ptr,ierr);CHKERRQ(ierr)
      if (this%inversion_option%forcomm%rank == 0) then
        mpi_int = 0
        ! broadcast to perturbation ranks 0 (forcomm_0)
!print *, this%driver%comm%rank, ' Bcast4 forcomm%rank == 0'
!call MPI_Barrier(this%inversion_option%forcomm_i%communicator,ierr);CHKERRQ(ierr)
        call MPI_Bcast(vec_ptr,num_measurements,MPI_DOUBLE_PRECISION,mpi_int, &
                      this%inversion_option%forcomm_i%communicator, &
                      ierr);CHKERRQ(ierr)
      endif
      if (this%inversion_option%forcomm%group_id > 1) then
        ! broadcast among perturbation ranks > 0
        mpi_int = 0
!call MPI_Barrier(this%inversion_option%forcomm%communicator,ierr);CHKERRQ(ierr)
!print *, this%driver%comm%rank, ' Bcast5 group_id > 1'
        call MPI_Bcast(vec_ptr,num_measurements,MPI_DOUBLE_PRECISION,mpi_int, &
                       this%inversion_option%forcomm%communicator, &
                       ierr);CHKERRQ(ierr)
      endif
      call VecRestoreArrayF90(this%inversion_aux%measurement_vec, &
                              vec_ptr,ierr);CHKERRQ(ierr)
    endif
    call VecCopy(this%inversion_aux%measurement_vec, &
                 this%inversion_aux%perturbation%base_measurement_vec, &
                 ierr);CHKERRQ(ierr)
    call VecGetArrayF90(this%inversion_aux%perturbation%base_measurement_vec, &
                        vec_ptr,ierr);CHKERRQ(ierr)
!    print *, 'gehbvec: ', inversion_option%forcomm%group_id, vec_ptr
    call VecRestoreArrayF90(this%inversion_aux%perturbation%base_measurement_vec, &
                            vec_ptr,ierr);CHKERRQ(ierr)
    if (associated(this%inversion_option%invcomm)) then
      call MatZeroEntries(this%inversion_aux%JsensitivityT, &
                          ierr);CHKERRQ(ierr)
    endif
  else
    call VecGetArrayF90(this%inversion_aux%measurement_vec, &
                        vec_ptr,ierr);CHKERRQ(ierr)
!    print *, 'gehmvec: ', inversion_option%forcomm%rank, &
!                           inversion_option%forcomm%group_id, vec_ptr
    call VecRestoreArrayF90(this%inversion_aux%measurement_vec, &
                            vec_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(this%inversion_aux%perturbation%base_measurement_vec, &
                        vec_ptr,ierr);CHKERRQ(ierr)
!    print *, 'gehbvec: ', inversion_option%forcomm%rank, &
!                           inversion_option%forcomm%group_id, vec_ptr
    call VecRestoreArrayF90(this%inversion_aux%perturbation%base_measurement_vec, &
                            vec_ptr,ierr);CHKERRQ(ierr)
!    print *, 'gehpert: ', inversion_option%forcomm%rank, &
!                           inversion_option%forcomm%group_id, this%inversion_aux%perturbation%pert
    call VecAXPY(this%inversion_aux%measurement_vec,-1.d0, &
                 this%inversion_aux%perturbation%base_measurement_vec, &
                 ierr);CHKERRQ(ierr)
    call VecScale(this%inversion_aux%measurement_vec, &
                  1.d0/this%inversion_aux%perturbation%pert, &
                  ierr);CHKERRQ(ierr)
  endif

  if (my_dof == 0) return

  if (associated(this%inversion_option%invcomm)) then
    do igroup = 1, this%inversion_option%num_process_groups
      if (my_dof + igroup - 1 > &
          this%inversion_aux%perturbation%ndof) exit
      if (igroup > 1) then
!        if (this%inversion_option%forcomm_i%rank == 0) then
          ! receive derivatives from each groups process 0
          call VecGetArrayF90(this%inversion_aux%measurement_vec, &
                              vec_ptr,ierr);CHKERRQ(ierr)
          call MPI_Probe(MPI_ANY_SOURCE,MPI_ANY_TAG, &
                         this%inversion_option%forcomm_i%communicator, &
                         status,ierr);CHKERRQ(ierr)
          iparameter = status(MPI_TAG)
          call MPI_Recv(vec_ptr,num_measurements,MPI_DOUBLE_PRECISION, &
                        status(MPI_SOURCE),MPI_ANY_TAG, &
                        this%inversion_option%forcomm_i%communicator, &
                        MPI_STATUS_IGNORE,ierr);CHKERRQ(ierr)
!          print *, 'recv: ', inversion_option%forcomm%rank, &
!                             inversion_option%forcomm%group_id,vec_ptr
          call VecRestoreArrayF90(this%inversion_aux%measurement_vec, &
                                  vec_ptr,ierr);CHKERRQ(ierr)
!        endif
      else
        iparameter = my_dof
      endif

      ! don't need to use the distributed vec, but why not
      call InvAuxScatMeasToDistMeas(this%inversion_aux, &
                                    this%inversion_aux%measurement_vec, &
                                    this%inversion_aux%dist_measurement_vec, &
                                    INVAUX_SCATFORWARD)
      call VecGetArrayF90(this%inversion_aux%dist_measurement_vec,vec_ptr, &
                          ierr);CHKERRQ(ierr)
      do i = 1, size(vec_ptr)
!        print *, 'gehmat: ', inversion_option%forcomm%group_id,iparameter, &
!                this%dist_measurement_offset+i, vec_ptr(i)
        call MatSetValue(this%inversion_aux%JsensitivityT, &
                        iparameter-1, &
                        this%dist_measurement_offset+i-1,vec_ptr(i), &
                        INSERT_VALUES,ierr);CHKERRQ(ierr)
      enddo
      call VecRestoreArrayF90(this%inversion_aux%dist_measurement_vec,vec_ptr, &
                              ierr);CHKERRQ(ierr)
    enddo

  else if (associated(this%inversion_option%forcomm_i)) then
!    if (this%inversion_option%forcomm_i%group_id == 1) then
      ! send derivatives to process 0 on group 1
      if (my_dof <= this%inversion_aux%perturbation%ndof) then
        call VecGetArrayF90(this%inversion_aux%measurement_vec, &
                            vec_ptr,ierr);CHKERRQ(ierr)
        mpi_int = 0
        call MPI_Send(vec_ptr,num_measurements,MPI_DOUBLE_PRECISION, &
                      mpi_int,my_dof, &
                      this%inversion_option%forcomm_i%communicator, &
                      ierr);CHKERRQ(ierr)
!        print *, 'send: ', inversion_option%forcomm%rank, &
!                            inversion_option%forcomm%group_id,vec_ptr
        call VecRestoreArrayF90(this%inversion_aux%measurement_vec, &
                                vec_ptr,ierr);CHKERRQ(ierr)
      endif
 !   endif
  endif

  if (associated(this%inversion_option%invcomm) .and. &
      iparameter == this%inversion_aux%perturbation%ndof) then
    call MatAssemblyBegin(this%inversion_aux%JsensitivityT, &
                          MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(this%inversion_aux%JsensitivityT, &
                        MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  endif

  if (.not.this%inversion_aux%qoi_is_full_vector) then
    ! revert back to base value
    this%inversion_aux%parameters(my_dof)%value = &
      this%inversion_aux%perturbation%base_value
  endif

end subroutine InvSubsurfPerturbationFillRow

! ************************************************************************** !

subroutine InvSubsurfFVSetOrigSoln(this)
  !
  ! Fills entries in partial Jacobian vectors for full-vector inversion
  !
  ! Author: Glenn Hammond
  ! Date: 09/23/22
  !
  class(inversion_subsurface_type) :: this

  type(inversion_coupled_soln_type), pointer :: solutions(:)
  PetscInt :: i
  PetscErrorCode :: ierr

  solutions => this%inversion_aux%coupled_aux%solutions
  if (associated(this%inversion_option%invcomm)) then
    do i = 1, size(solutions)
      call VecCopy(solutions(i)%perturbed_saturation_solution, &
                   solutions(i)%original_saturation_solution, &
                   ierr);CHKERRQ(ierr)
      if (solutions(i)%perturbed_solute_solution /= &
          PETSC_NULL_VEC) then
        call VecCopy(solutions(i)%perturbed_solute_solution, &
                     solutions(i)%original_solute_solution, &
                     ierr);CHKERRQ(ierr)
      endif
    enddo
  endif
  if (associated(this%inversion_option%forcomm_i)) then
!print *, this%driver%comm%rank, ' Bcast11 all'
!      call MPI_Barrier(this%driver%comm%communicator, &
!                       ierr);CHKERRQ(ierr)
    do i = 1, size(solutions)
!print *, this%driver%comm%rank, ' Bcast12 forcomm_i'
!      call MPI_Barrier(this%inversion_option%forcomm_i%communicator, &
!                       ierr);CHKERRQ(ierr)
      call InvAuxBCastVecForCommI(this%inversion_option%forcomm_i, &
                                  solutions(i)%original_saturation_solution, &
                                  this%driver)
      if (solutions(i)%perturbed_solute_solution /= &
          PETSC_NULL_VEC) then
!print *, this%driver%comm%rank, ' Bcast13 forcomm_i'
!        call MPI_Barrier(this%inversion_option%forcomm_i%communicator, &
!                         ierr);CHKERRQ(ierr)
        call InvAuxBCastVecForCommI(this%inversion_option%forcomm_i, &
                                    solutions(i)%perturbed_solute_solution, &
                                    this%driver)
      endif
    enddo
  endif
  call InversionCoupledAuxReset(this%inversion_aux%coupled_aux)


end subroutine InvSubsurfFVSetOrigSoln

! ************************************************************************** !

subroutine InvSubsurfFVCalcPartialJs(this,idof_pert)
  !
  ! Fills entries in partial Jacobian vectors for full-vector inversion
  !
  ! Author: Glenn Hammond
  ! Date: 09/23/22
  !
  use Communicator_Aux_module

  implicit none

  class(inversion_subsurface_type) :: this
  PetscInt, intent(in) :: idof_pert ! preserve intent to catch errors

  type(inversion_coupled_soln_type), pointer :: solutions(:)
  type(comm_type), pointer :: comm
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: icomm
  PetscInt :: i
  PetscMPIInt :: recv_rank
  PetscInt :: vec_size
  PetscInt :: idof_pert_recv
  PetscErrorCode :: ierr

  if (this%inversion_aux%perturbation%idof_pert < 0) then
    call this%driver%PrintErrMsg('InvSubsurfFVCalcPartialJs() called  &
      &with an idof_pert < 0')
  endif

  solutions => this%inversion_aux%coupled_aux%solutions
  do i = 1, size(solutions)
    call InvCoupledUpdateSolnVecs(idof_pert, &
      solutions(i)%perturbed_saturation_solution, &
      solutions(i)%original_saturation_solution, &
      solutions(i)%dsaturation_dparameter, &
      this%inversion_aux%perturbation%pert)
    if (solutions(i)%perturbed_solute_solution /= &
        PETSC_NULL_VEC) then
      call InvCoupledUpdateSolnVecs(idof_pert, &
        solutions(i)%perturbed_solute_solution, &
        solutions(i)%original_solute_solution, &
        solutions(i)%dsolute_dparameter, &
        this%inversion_aux%perturbation%pert)
    endif
  enddo

  comm => this%inversion_option%forcomm_i
  if (associated(comm)) then
    call VecGetLocalSize(solutions(1)%dsaturation_dparameter(idof_pert), &
                         vec_size,ierr);CHKERRQ(ierr)
    if (comm%rank == 0) then
      do icomm = 1, min(comm%size-1, &
                        this%inversion_aux%perturbation%ndof-idof_pert)
        recv_rank = icomm
!        call MPI_Probe(icomm,mpi_int,comm%communicator, &
!                       MPI_STATUS_IGNORE,ierr);CHKERRQ(ierr)
        do i = 1, size(solutions)
          idof_pert_recv = idof_pert + icomm
          call VecGetArrayF90(solutions(i)% &
                                dsaturation_dparameter(idof_pert_recv), &
                              vec_ptr,ierr);CHKERRQ(ierr)
          call MPI_Recv(vec_ptr,vec_size,MPI_DOUBLE_PRECISION, &
                        recv_rank,MPI_ANY_TAG,comm%communicator, &
                        MPI_STATUS_IGNORE,ierr);CHKERRQ(ierr)
          call VecRestoreArrayF90(solutions(i)% &
                                    dsaturation_dparameter(idof_pert_recv), &
                                  vec_ptr,ierr);CHKERRQ(ierr)
          if (solutions(i)%perturbed_solute_solution /= &
              PETSC_NULL_VEC) then
            call VecGetArrayF90(solutions(i)% &
                                  dsolute_dparameter(idof_pert_recv), &
                                vec_ptr,ierr);CHKERRQ(ierr)
            call MPI_Recv(vec_ptr,vec_size,MPI_DOUBLE_PRECISION, &
                          recv_rank,MPI_ANY_TAG,comm%communicator, &
                          MPI_STATUS_IGNORE,ierr);CHKERRQ(ierr)
            call VecRestoreArrayF90(solutions(i)% &
                                      dsolute_dparameter(idof_pert_recv), &
                                    vec_ptr,ierr);CHKERRQ(ierr)
          endif
        enddo
      enddo
    else
      do i = 1, size(solutions)
        call VecGetArrayF90(solutions(i)%dsaturation_dparameter(idof_pert), &
                            vec_ptr,ierr);CHKERRQ(ierr)
        call MPI_Send(vec_ptr,vec_size,MPI_DOUBLE_PRECISION,ZERO_INTEGER_MPI, &
                      idof_pert,comm%communicator,ierr);CHKERRQ(ierr)
        call VecRestoreArrayF90(solutions(i)% &
                                  dsaturation_dparameter(idof_pert), &
                                vec_ptr,ierr);CHKERRQ(ierr)
        if (solutions(i)%perturbed_solute_solution /= &
            PETSC_NULL_VEC) then
          call VecGetArrayF90(solutions(i)%dsolute_dparameter(idof_pert), &
                              vec_ptr,ierr);CHKERRQ(ierr)
          call MPI_Send(vec_ptr,vec_size,MPI_DOUBLE_PRECISION, &
                        ZERO_INTEGER_MPI,idof_pert,comm%communicator, &
                        ierr);CHKERRQ(ierr)
          call VecGetArrayF90(solutions(i)%dsolute_dparameter(idof_pert), &
                              vec_ptr,ierr);CHKERRQ(ierr)
        endif
      enddo
    endif
  endif
  call InversionCoupledAuxReset(this%inversion_aux%coupled_aux)

  if (.not.this%inversion_aux%qoi_is_full_vector .and. idof_pert > 0) then
    ! revert back to base value
    this%inversion_aux% &
        parameters(this%inversion_aux%perturbation%idof_pert)%value = &
      this%inversion_aux%perturbation%base_value
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
  if (associated(this%inversion_aux%perturbation)) then
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
  if (this%inversion_aux%qoi_is_full_vector) then
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
  use HDF5_Aux_module
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
  character(len=MAXSTRINGLENGTH) :: string

  if (.not.associated(this%realization)) then
    call this%driver%PrintErrMsg('InvSubsurfOutputSensitivityHDF5 must be &
           &called before the forward simulation is destroyed.')
  endif

  !HDF5 formatted output
  string = trim(filename_prefix) // '.h5'
  call HDF5FileOpen(string,file_id,PETSC_TRUE, &
                    this%realization%option)

  call OutputHDF5WriteStructCoordGroup(file_id, &
                                       this%realization%discretization, &
                                       this%realization%patch%grid, &
                                       this%realization%option)
  ! create a group for the data set
  this%realization%option%time = 0.d0
  write(string,'(''Time:'',es13.5,x,a1)') &
        this%realization%option%time/this%realization%output_option%tconv, &
        this%realization%output_option%tunit
  call HDF5GroupOpenOrCreate(file_id,string,grp_id,this%realization%option)

  num_measurement = size(this%inversion_aux%measurements)
  call VecDuplicate(this%inversion_aux%dist_measurement_vec, &
                    row_vec,ierr);CHKERRQ(ierr)
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
      trim(StringWrite(this%inversion_aux%measurements(imeasurement)%cell_id))
    call HDF5WriteStructDataSetFromVec(string,this%realization, &
                                       this%realization%field%work,grp_id, &
                                       H5T_NATIVE_DOUBLE)
  enddo
  call VecDestroy(row_vec,ierr);CHKERRQ(ierr)
  call HDF5GroupClose(grp_id,this%driver)
  call OutputHDF5CloseFile(this%realization%option,file_id)

end subroutine InvSubsurfOutputSensitivityHDF5

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

  num_measurements = size(this%inversion_aux%measurements)
  if (this%driver%PrintToScreen()) then
    do i = 1, num_measurements
      call InvMeasurePrintComparison(STDOUT_UNIT, &
                                     this%inversion_aux%measurements(i),i==1, &
                                     i==num_measurements, &
                                     this%realization%option)
    enddo
  endif
  if (this%driver%PrintToFile()) then
    do i = 1, num_measurements
      call InvMeasurePrintComparison(this%driver%fid_out, &
                                     this%inversion_aux%measurements(i),i==1, &
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

  if (this%inversion_aux%qoi_is_full_vector) return

  num_parameters = size(this%inversion_aux%parameters)
  if (this%driver%PrintToScreen()) then
    do i = 1, num_parameters
      call InversionParameterPrint(STDOUT_UNIT, &
                                   this%inversion_aux%parameters(i),i==1, &
                                   i==num_parameters, &
                                   this%realization%option)
    enddo
  endif
  if (this%driver%PrintToFile()) then
    do i = 1, num_parameters
      call InversionParameterPrint(this%driver%fid_out, &
                                   this%inversion_aux%parameters(i),i==1, &
                                   i==num_parameters, &
                                   this%realization%option)
    enddo
  endif

end subroutine InvSubsurfPrintCurParamValues

! ************************************************************************** !

subroutine InvSubsurfPrintCurParamUpdate(this)
  !
  ! Prints the current update of parameters in the list of parameters
  ! Doesn't work for full vector parameter sets
  !
  ! Author: Glenn Hammond
  ! Date: 11/17/22
  !
  use Driver_class

  class(inversion_subsurface_type) :: this

  PetscInt :: i
  PetscInt :: num_parameters

  if (this%inversion_aux%qoi_is_full_vector) return

  call InvAuxCopyParamToFromParamVec(this%inversion_aux, &
                                     INVAUX_PARAMETER_UPDATE, &
                                     INVAUX_COPY_FROM_VEC)

  num_parameters = size(this%inversion_aux%parameters)
  if (this%driver%PrintToScreen()) then
    do i = 1, num_parameters
      call InversionParameterPrintUpdate(STDOUT_UNIT, &
                                   this%inversion_aux%parameters(i),i==1, &
                                   i==num_parameters)
    enddo
  endif
  if (this%driver%PrintToFile()) then
    do i = 1, num_parameters
      call InversionParameterPrintUpdate(this%driver%fid_out, &
                                   this%inversion_aux%parameters(i),i==1, &
                                   i==num_parameters)
    enddo
  endif

end subroutine InvSubsurfPrintCurParamUpdate

! ************************************************************************** !

subroutine InvSubsurfRestartIteration(this)
  !
  ! Reads the restart iteration from an inversion checkpoint file
  !
  ! Author: Glenn Hammond
  ! Date: 12/09/22
  !
  use hdf5
  use Driver_class
  use HDF5_Aux_module
  use String_module

  class(inversion_subsurface_type) :: this

  integer(HID_T) :: file_id
  PetscInt :: restart_iteration
  PetscInt :: last_iteration

  if (.not.this%inversion_aux%startup_phase .or. &
      len_trim(this%restart_filename) == 0) return

  call this%driver%PrintMsg('Reading restart iteration from inversion &
                            &checkpoint file "' // &
                            trim(this%restart_filename) // '".')
  call HDF5FileOpen(this%restart_filename,file_id,PETSC_FALSE,this%driver)
  call HDF5AttributeRead(file_id,H5T_NATIVE_INTEGER,'Last Iteration', &
                         last_iteration,this%driver)
  if (Initialized(this%restart_iteration)) then
    if (this%restart_iteration > last_iteration) then
      call this%driver%PrintErrMsg('Restart iteration ' // &
                              trim(StringWrite(this%restart_iteration)) // &
                              ' is beyond the last checkpoint iteration (' // &
                              trim(StringWrite(last_iteration)) // ').')
    endif
    restart_iteration = this%restart_iteration
    call this%driver%PrintMsg('Restarting at iteration (' // &
                              trim(StringWrite(restart_iteration)) // ').')
  else
    restart_iteration = last_iteration
    call this%driver%PrintMsg('Restarting at last iteration (' // &
                              trim(StringWrite(restart_iteration)) // ').')
  endif
  call HDF5FileClose(file_id,this%driver)
  this%restart_iteration = restart_iteration
  this%iteration = this%restart_iteration

end subroutine InvSubsurfRestartIteration

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

subroutine InversionSubsurfaceStrip(this)
  !
  ! Deallocates members of inversion Subsurface
  !
  ! Author: Glenn Hammond
  ! Date: 09/20/21
  !
  use Utility_module

  class(inversion_subsurface_type) :: this

  PetscErrorCode :: ierr

  call InversionBaseStrip(this)

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

end subroutine InversionSubsurfaceStrip

end module Inversion_Subsurface_class
