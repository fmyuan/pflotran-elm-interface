module Inversion_Subsurface_class

#include "petsc/finclude/petscksp.h"
  use petscksp

  use PFLOTRAN_Constants_module
  use Inversion_Aux_module
  use Inversion_TS_Aux_module
  use Inversion_Measurement_Aux_module
  use Inversion_Parameter_module
  use Inversion_Base_class
  use Realization_Subsurface_class
  use Simulation_Subsurface_class

  implicit none

  private

  type, public, extends(inversion_base_type) :: inversion_subsurface_type
    character(len=MAXSTRINGLENGTH) :: forward_simulation_filename
    class(simulation_subsurface_type), pointer :: forward_simulation
    class(realization_subsurface_type), pointer :: realization
    type(inversion_aux_type), pointer :: inversion_aux
    type(inversion_measurement_aux_type), pointer :: measurements(:)
    type(inversion_parameter_type), pointer :: parameters(:)
    PetscInt :: dist_measurement_offset
    PetscInt :: parameter_offset
    PetscInt :: iqoi(2)
    PetscInt :: iobsfunc
    PetscInt :: n_qoi_per_cell
    Vec :: quantity_of_interest
    Vec :: ref_quantity_of_interest
    Vec :: measurement_vec
    Vec :: dist_measurement_vec
    Vec :: parameter_vec
    Vec :: ref_parameter_vec
    character(len=MAXWORDLENGTH) :: ref_qoi_dataset_name
    PetscBool :: print_sensitivity_jacobian
    PetscBool :: debug_adjoint
    PetscInt :: debug_verbosity
    PetscBool :: local_adjoint
    VecScatter :: scatter_global_to_measurement
    VecScatter :: scatter_measure_to_dist_measure
    PetscBool :: qoi_is_full_vector
    PetscBool :: first_iteration
  contains
    procedure, public :: Init => InversionSubsurfaceInit
    procedure, public :: ReadBlock => InversionSubsurfReadBlock
    procedure, public :: Initialize => InversionSubsurfInitialize
    procedure, public :: Step => InversionSubsurfaceStep
    procedure, public :: ConnectToForwardRun => InvSubsurfConnectToForwardRun
    procedure, public :: CalculateSensitivity => InvSubsurfCalculateSensitivity
    procedure, public :: ScaleSensitivity => InvSubsurfScaleSensitivity
    procedure, public :: OutputSensitivity => InvSubsurfOutputSensitivity
    procedure, public :: Invert => InversionSubsurfaceInvert
    procedure, public :: Strip => InversionSubsurfaceStrip
  end type inversion_subsurface_type

  public :: InversionSubsurfaceCreate, &
            InversionSubsurfaceInit, &
            InversionSubsurfReadSelectCase, &
            InversionSubsurfInitialize, &
            InvSubsurfConnectToForwardRun, &
            InvSubsurfOutputSensitivity, &
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
  use Driver_module

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
  use Driver_module
  use ZFlow_Aux_module, only : zflow_calc_adjoint

  class(inversion_subsurface_type) :: this
  class(driver_type), pointer :: driver

  call InversionBaseInit(this,driver)

  this%quantity_of_interest = PETSC_NULL_VEC
  this%iqoi = UNINITIALIZED_INTEGER
  this%iobsfunc = UNINITIALIZED_INTEGER
  this%n_qoi_per_cell = UNINITIALIZED_INTEGER
  this%measurement_vec = PETSC_NULL_VEC
  this%dist_measurement_vec = PETSC_NULL_VEC
  this%parameter_vec = PETSC_NULL_VEC
  this%ref_parameter_vec = PETSC_NULL_VEC
  this%ref_quantity_of_interest = PETSC_NULL_VEC
  this%ref_qoi_dataset_name = ''
  this%forward_simulation_filename = ''
  this%print_sensitivity_jacobian = PETSC_FALSE
  this%debug_adjoint = PETSC_FALSE
  this%local_adjoint = PETSC_FALSE
  this%debug_verbosity = UNINITIALIZED_INTEGER
  this%scatter_global_to_measurement = PETSC_NULL_VECSCATTER
  this%dist_measurement_offset = UNINITIALIZED_INTEGER
  this%qoi_is_full_vector = PETSC_FALSE
  this%first_iteration = PETSC_TRUE

  nullify(this%measurements)
  nullify(this%parameters)

  nullify(this%forward_simulation)
  nullify(this%realization)
  nullify(this%inversion_aux)

  zflow_calc_adjoint = PETSC_TRUE

end subroutine InversionSubsurfaceInit

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

  error_string = 'Test Inversion'

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
  use Variables_module, only : ELECTRICAL_CONDUCTIVITY, &
                               PERMEABILITY, POROSITY, &
                               LIQUID_PRESSURE, LIQUID_SATURATION, &
                               SOLUTE_CONCENTRATION
  use Material_Aux_module, only : POROSITY_BASE
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
    case('QUANTITY_OF_INTEREST')
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,keyword,error_string)
      call StringToUpper(word)
      select case(word)
        case('ELECTRICAL_CONDUCTIVITY')
          this%iqoi(1) = ELECTRICAL_CONDUCTIVITY
        case('PERMEABILITY')
          this%iqoi(1) = PERMEABILITY
        case('POROSITY')
          this%iqoi(1) = POROSITY
          this%iqoi(2) = POROSITY_BASE
        case default
          call InputKeywordUnrecognized(input,word,trim(error_string)// &
                                        & ','//trim(keyword),option)
      end select
    case('REFERENCE_QUANTITY_OF_INTEREST')
      call InputReadNChars(input,option,this%ref_qoi_dataset_name, &
                           MAXWORDLENGTH,PETSC_TRUE)
      call InputErrorMsg(input,option,'DATASET NAME', &
                         keyword)
    case('OBSERVATION_FUNCTION')
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,keyword,error_string)
      call StringToUpper(word)
      select case(word)
        case('LIQUID_PRESSURE')
          this%iobsfunc = LIQUID_PRESSURE
        case('LIQUID_SATURATION')
          this%iobsfunc = LIQUID_SATURATION
        case('SOLUTE_CONCENTRATION')
          this%iobsfunc = SOLUTE_CONCENTRATION
        case default
          call InputKeywordUnrecognized(input,word,trim(error_string)// &
                                        & ','//trim(keyword),option)
      end select
    case('MEASUREMENTS')
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
          case default
            call InputKeywordUnrecognized(input,keyword,error_string,option)
        end select
        if (associated(last_measurement)) then
          last_measurement%next => new_measurement
          new_measurement%id = last_measurement%id + 1
        else
          first_measurement => new_measurement
          first_measurement%id = 1
        endif
        last_measurement => new_measurement
        nullify(new_measurement)
      enddo
      call InputPopBlock(input,option)
      if (.not.associated(last_measurement)) then
        option%io_buffer = 'No measurement found in inversion measurement block.'
        call PrintErrMsg(option)
      else if (associated(this%measurements)) then
        option%io_buffer = 'Parameters may only be defined in a single block.'
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
        option%io_buffer = 'Measurements may only be defined in a single block.'
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
    case('LOCAL_ADJOINT')
      this%local_adjoint = PETSC_TRUE
    case default
      found = PETSC_FALSE
  end select

end subroutine InversionSubsurfReadSelectCase

! ************************************************************************** !

subroutine InversionSubsurfInitialize(this)
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
  use Option_module
  use Patch_module
  use String_module
  use Variables_module, only : ELECTRICAL_CONDUCTIVITY, &
                               PERMEABILITY, POROSITY, &
                               LIQUID_PRESSURE, LIQUID_SATURATION, &
                               SOLUTE_CONCENTRATION
  use Material_Aux_module, only : POROSITY_BASE

  class(inversion_subsurface_type) :: this

  type(patch_type), pointer :: patch
  type(inversion_forward_aux_type), pointer :: inversion_forward_aux
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i
  PetscInt :: local_id
  PetscInt :: num_measurements, num_measurements_local
  PetscInt :: num_parameters, num_parameters_local
  PetscInt, allocatable :: int_array(:), int_array2(:)
  PetscReal, pointer :: vec_ptr(:)
  Vec :: v, v2
  IS :: is_petsc
  IS :: is_measure
  PetscMPIInt :: mpi_int
  PetscErrorCode :: ierr

  nullify(vec_ptr)

  patch => this%realization%patch
  if (.not.associated(this%inversion_aux)) then
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
        if (len_trim(this%parameters(i)%strata_name) > 0) then
          num_parameters = num_parameters + 1
        else
          num_parameters = num_parameters + patch%grid%nmax
        endif
      enddo
      if (num_parameters == patch%grid%nmax) then
        this%qoi_is_full_vector = PETSC_TRUE
        select case(this%parameters(1)%parameter_name)
          case('ELECTRICAL_CONDUCTIVITY')
            this%iqoi(1) = ELECTRICAL_CONDUCTIVITY
          case('PERMEABILITY')
            this%iqoi(1) = PERMEABILITY
          case('POROSITY')
            this%iqoi(1) = POROSITY
            this%iqoi(2) = POROSITY_BASE
          case default
            call this%driver%PrintErrMsg('Unrecognized parameter in &
                                         &InversionSubsurfInitialize')
        end select
        num_parameters_local = patch%grid%nlmax*this%n_qoi_per_cell
      else
        num_parameters_local = PETSC_DECIDE
      endif
    else
      call this%driver%PrintErrMsg('No inversion parameters defined.')
    endif
    ! JsensitivityT is the transpose of the sensitivity Jacobian
    ! with num measurement columns and num parameter rows
    call MatCreateDense(this%driver%comm%mycomm, &
                        num_parameters_local,PETSC_DECIDE, &
                        num_parameters,num_measurements, &
                        PETSC_NULL_SCALAR, &
                        this%inversion_aux%JsensitivityT,ierr);CHKERRQ(ierr)
    call MatZeroEntries(this%inversion_aux%JsensitivityT,ierr);CHKERRQ(ierr)
    ! cannot pass in this%measurement_vec as it is initialized to
    ! PETSC_NULL_VEC and MatCreateVecs keys off that input
    call MatCreateVecs(this%inversion_aux%JsensitivityT,v,v2, &
                       ierr);CHKERRQ(ierr)
    this%dist_measurement_vec = v
    this%parameter_vec = v2
    call MatGetLocalSize(this%inversion_aux%JsensitivityT, &
                         num_parameters_local, &
                         num_measurements_local,ierr);CHKERRQ(ierr)
    allocate(int_array(2),int_array2(2))
    int_array(1) = num_parameters_local
    int_array(2) = num_measurements_local
    int_array2(:) = 0
    call MPI_Exscan(int_array,int_array2,TWO_INTEGER_MPI,MPIU_INTEGER, &
                    MPI_SUM,this%driver%comm%mycomm,ierr);CHKERRQ(ierr)
    this%parameter_offset = int_array2(1)
    this%dist_measurement_offset = int_array2(2)
    deallocate(int_array,int_array2)
    call VecCreateSeq(PETSC_COMM_SELF,num_measurements, &
                      this%measurement_vec,ierr);CHKERRQ(ierr)

    ! map coordinates to cell ids
    allocate(int_array(num_measurements))
    int_array = -1
    do i = 1, num_measurements
      if (Initialized(this%measurements(i)%coordinate%x)) then
        call GridGetLocalIDFromCoordinate(patch%grid, &
                                          this%measurements(i)%coordinate, &
                                          this%realization%option,local_id)
        if (Initialized(local_id)) then
          int_array(i) = patch%grid%nG2A(patch%grid%nL2G(local_id))
        endif
      endif
    enddo
    mpi_int = num_measurements
    call MPI_Allreduce(MPI_IN_PLACE,int_array,mpi_int,MPIU_INTEGER,MPI_MAX, &
                       this%driver%comm%mycomm,ierr)
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

    ! map measurement vecs to the solution vector
    do i = 1, num_measurements
      int_array(i) = this%measurements(i)%cell_id
    enddo
    i = maxval(int_array)
    if (i > patch%grid%nmax) then
      this%realization%option%io_buffer = 'A measurement cell ID (' // &
        trim(StringWrite(i)) // ') is beyond the maximum cell ID.'
      call PrintErrMsg(this%realization%option)
    endif
    int_array = int_array - 1
    call DiscretAOApplicationToPetsc(this%realization%discretization, &
                                     int_array)
    call ISCreateGeneral(this%driver%comm%mycomm,num_measurements, &
                         int_array, &
                         PETSC_COPY_VALUES,is_petsc,ierr);CHKERRQ(ierr)
    deallocate(int_array)
    call VecScatterCreate(this%realization%field%work,is_petsc, &
                          this%measurement_vec,PETSC_NULL_IS, &
                          this%scatter_global_to_measurement, &
                          ierr);CHKERRQ(ierr)
    call ISDestroy(is_petsc,ierr);
    call ISCreateStride(this%driver%comm%mycomm,num_measurements, &
                        ZERO_INTEGER,ONE_INTEGER, &
                        is_measure,ierr);CHKERRQ(ierr)
    call VecScatterCreate(this%measurement_vec,is_measure, &
                          this%dist_measurement_vec,PETSC_NULL_IS, &
                          this%scatter_measure_to_dist_measure, &
                          ierr);CHKERRQ(ierr)
    call ISDestroy(is_petsc,ierr);

    this%inversion_aux%scatter_global_to_measurement = &
      this%scatter_global_to_measurement
    this%inversion_aux%measurements => this%measurements
    this%inversion_aux%measurement_vec = this%measurement_vec

    inversion_forward_aux => InversionForwardAuxCreate()
    inversion_forward_aux%iobsfunc = this%iobsfunc
    inversion_forward_aux%scatter_global_to_measurement = &
      this%scatter_global_to_measurement
    inversion_forward_aux%measurements => this%measurements
    inversion_forward_aux%measurement_vec = this%measurement_vec
    ! set up pointer to M matrix
    this%inversion_aux%inversion_forward_aux => inversion_forward_aux
  endif

  inversion_forward_aux => this%inversion_aux%inversion_forward_aux
  inversion_forward_aux%M_ptr = &
    this%forward_simulation%flow_process_model_coupler%timestepper%solver%M
! create inversion_ts_aux for first time step
  nullify(inversion_forward_aux%first) ! must pass in null object
  inversion_forward_aux%first => &
    InversionTSAuxCreate(inversion_forward_aux%first, &
                         inversion_forward_aux%M_ptr)
  inversion_forward_aux%current => inversion_forward_aux%first
  call InvTSAuxAllocate(inversion_forward_aux%first, &
                        inversion_forward_aux%M_ptr, &
                        this%realization%option%nflowdof, &
                        patch%grid%nlmax)

end subroutine InversionSubsurfInitialize

! ************************************************************************** !

subroutine InversionSubsurfaceStep(this)
  !
  ! Execute a simulation
  !
  ! Author: Glenn Hammond
  ! Date: 05/27/21

  use Option_module
  use Factory_Forward_module

  class(inversion_subsurface_type) :: this

  type(option_type), pointer :: option

  option => OptionCreate()
  write(option%group_prefix,'(i6)') this%iteration
  option%group_prefix = 'Run' // trim(adjustl(option%group_prefix))
  call OptionSetDriver(option,this%driver)
  call OptionSetInversionOption(option,this%inversion_option)
  call FactoryForwardInitialize(this%forward_simulation, &
                                this%forward_simulation_filename,option)
  this%realization => this%forward_simulation%realization
  call this%Initialize()
  call this%forward_simulation%InitializeRun()
  call this%ConnectToForwardRun()
  if (option%status == PROCEED) then
    call this%forward_simulation%ExecuteRun()
  endif
  call this%CalculateSensitivity()
  call this%OutputSensitivity('')
  nullify(this%realization)
  call this%forward_simulation%FinalizeRun()
  call this%forward_simulation%Strip()
  deallocate(this%forward_simulation)
  nullify(this%forward_simulation)
  call this%Invert()

  this%converg_flag = PETSC_FALSE
  if (this%iteration > this%maximum_iteration) this%converg_flag = PETSC_TRUE

end subroutine InversionSubsurfaceStep

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
  use Init_Subsurface_module
  use Material_module
  use Region_module
  use Strata_module
  use String_module

  class(inversion_subsurface_type) :: this

  PetscInt :: i
  type(strata_type), pointer :: strata
  type(material_property_type), pointer :: material_property
  type(region_type), pointer :: region
  PetscErrorCode :: ierr

  this%realization%patch%aux%inversion_forward_aux => &
    this%inversion_aux%inversion_forward_aux
  call InvForwardAuxResetMeasurements(this%inversion_aux%inversion_forward_aux)

  if (this%qoi_is_full_vector) then
      ! on first pass, store and set thereafter
    if (this%quantity_of_interest == PETSC_NULL_VEC) then
      ! can't move VecDuplicate earlier as it will not be null for the
      ! conditional above
      call VecDuplicate(this%realization%field%work, &
                        this%quantity_of_interest,ierr);CHKERRQ(ierr)
      call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                  this%realization%field%work_loc, &
                                  this%iqoi(1),this%iqoi(2))
      call DiscretizationLocalToGlobal(this%realization%discretization, &
                                      this%realization%field%work_loc, &
                                      this%quantity_of_interest,ONEDOF)
    endif
    call DiscretizationGlobalToLocal(this%realization%discretization, &
                                     this%quantity_of_interest, &
                                     this%realization%field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 this%iqoi(1),this%iqoi(2))
  else
    if (this%quantity_of_interest == PETSC_NULL_VEC) then
      ! can't move VecDuplicate earlier as it will not be null for the
      ! conditional above
      call VecDuplicate(this%realization%field%work, &
                        this%quantity_of_interest,ierr);CHKERRQ(ierr)
    endif
    if (this%first_iteration) then
      this%first_iteration = PETSC_FALSE
      do i = 1, size(this%parameters)
        strata => this%realization%patch%strata_list%first
        do
          if (.not.associated(strata)) exit
          if (StringCompareIgnoreCase(this%parameters(i)%strata_name, &
                                      strata%region_name)) then
            exit
          endif
          strata => strata%next
        enddo
        if (associated(strata)) then
          material_property => &
            MaterialPropGetPtrFromArray(strata%material_property_name, &
                              this%realization%patch%material_property_array)
          region => RegionGetPtrFromList(strata%region_name, &
                                        this%realization%patch%region_list)
          this%parameters(i)%imat = abs(material_property%internal_id)
          this%parameters(i)%value = material_property%permeability(1,1)
        endif
      enddo
    endif
  endif

end subroutine InvSubsurfConnectToForwardRun

! ************************************************************************** !

subroutine InvSubsurfCalculateSensitivity(this)
  !
  ! Calculates sensitivity matrix Jsensitivity
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/21
  !
  use Option_module
  use String_module
  use Timer_class
  use Units_module
  use Utility_module

  class(inversion_subsurface_type) :: this

  type(inversion_aux_type), pointer :: inversion_aux
  type(inversion_forward_ts_aux_type), pointer :: cur_inversion_ts_aux
  type(inversion_forward_ts_aux_type), pointer :: prev_inversion_ts_aux
  type(option_type), pointer :: option
  class(timer_type), pointer :: timer
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: imeasurement
  PetscErrorCode :: ierr

  option => this%realization%option
  inversion_aux => this%inversion_aux

  timer => TimerCreate()

  call timer%Start()

  call PrintHeader('SENSITIVITY JACOBIAN',option)

  ! initialize first_lambda flag
  do imeasurement = 1, size(this%measurements)
    this%measurements(imeasurement)%first_lambda = PETSC_FALSE
    if (.not.this%measurements(imeasurement)%measured) then
      option%io_buffer = 'Measurement at cell ' // &
        StringWrite(this%measurements(imeasurement)%cell_id)
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

  ! work backward through list
  call OptionPrint(' Working backward through inversion_ts_aux list &
                   &calculating lambdas.',option)
  cur_inversion_ts_aux => inversion_aux%inversion_forward_aux%last
  do
    if (.not.associated(cur_inversion_ts_aux)) exit
    if (OptionPrintToScreen(option)) then
      print *, '  call InvSubsurfCalcLambda: ', cur_inversion_ts_aux%timestep
    endif
    call InvSubsurfCalcLambda(this,cur_inversion_ts_aux)
    cur_inversion_ts_aux => cur_inversion_ts_aux%prev
  enddo

  ! work forward through list
  call OptionPrint(' Working forward through inversion_ts_aux list &
                   &calculating sensitivity coefficients.',option)
  call MatZeroEntries(inversion_aux%JsensitivityT,ierr);CHKERRQ(ierr)
  cur_inversion_ts_aux => inversion_aux%inversion_forward_aux%first
  do
    if (.not.associated(cur_inversion_ts_aux)) exit
    if (OptionPrintToScreen(option)) then
      print *, '  call InvSubsurfAddSensitivity: ', cur_inversion_ts_aux%timestep
    endif
    call InvSubsurfAddSensitivity(this,cur_inversion_ts_aux)
    cur_inversion_ts_aux => cur_inversion_ts_aux%next
  enddo
  call MatAssemblyBegin(inversion_aux%JsensitivityT, &
                        MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(inversion_aux%JsensitivityT, &
                      MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  call InvForwardAuxDestroyList(inversion_aux%inversion_forward_aux, &
                                PETSC_TRUE)

  call timer%Stop()
  option%io_buffer = '    ' // &
    trim(StringWrite('(f20.1)',timer%GetCumulativeTime())) // &
    ' seconds to build all sensitivities.'
  call PrintMsg(option)
  call TimerDestroy(timer)

end subroutine InvSubsurfCalculateSensitivity

! ************************************************************************** !

subroutine InvSubsurfCalcLambda(this,inversion_forward_ts_aux)
  !
  ! Calculates sensitivity matrix Jsensitivity
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
  use Utility_module
  use Variables_module

  use PM_Base_class
  use PM_ZFlow_class

  class(inversion_subsurface_type) :: this
  type(inversion_forward_ts_aux_type), pointer :: inversion_forward_ts_aux

  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(inversion_aux_type), pointer :: inversion_aux
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(solver_type), pointer :: solver
  type(zflow_auxvar_type), pointer :: zflow_auxvars(:,:)
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
  Vec :: p        ! ndof_vec
  Vec :: rhs      ! ndof_vec
  ! derivative of residual at k+1 time level wrt unknown at k time level
  ! times lambda at k time level
  Vec :: dReskp1_duk_lambdak
  Vec :: natural_vec
  class(timer_type), pointer :: timer
  PetscReal, parameter :: tol = 1.d-6
  PetscErrorCode :: ierr

  nullify(vec_ptr)

  solver => this%forward_simulation%flow_process_model_coupler% &
              timestepper%solver
  option => this%realization%option
  discretization => this%realization%discretization
  patch => this%realization%patch
  grid => patch%grid
  inversion_aux => this%inversion_aux
  zflow_auxvars => patch%aux%ZFlow%auxvars
  ndof = option%nflowdof

  timer => TimerCreate()

  call timer%Start()

  work = this%realization%field%work ! DO NOT DESTROY!
  call VecDuplicate(work,onedof_vec,ierr);CHKERRQ(ierr)
  ndof_vec = this%realization%field%flow_xx ! DO NOT DESTROY!
  call VecDuplicate(ndof_vec,p,ierr);CHKERRQ(ierr)
  call VecDuplicate(ndof_vec,rhs,ierr);CHKERRQ(ierr)
  call VecDuplicate(ndof_vec,dReskp1_duk_lambdak,ierr);CHKERRQ(ierr)
  call DiscretizationCreateVector(discretization,ONEDOF, &
                                  natural_vec,NATURAL,option)

  if (this%debug_verbosity > 2) then
    if (OptionPrintToScreen(option)) print *, 'M'
    call MatView(inversion_forward_ts_aux%dResdu,PETSC_VIEWER_STDOUT_WORLD, &
                 ierr);CHKERRQ(ierr)
  endif
  call VecDuplicateVecsF90(ndof_vec,size(this%measurements), &
                           inversion_forward_ts_aux%lambda,ierr);CHKERRQ(ierr)
  call KSPSetOperators(solver%ksp,inversion_forward_ts_aux%dResdu, &
                       inversion_forward_ts_aux%dResdu,ierr);CHKERRQ(ierr)
  do imeasurement = 1, size(this%measurements)
    if (Initialized(this%measurements(imeasurement)%time) .and. &
        inversion_forward_ts_aux%time > this%measurements(imeasurement)%time + tol) then
      call VecZeroEntries(inversion_forward_ts_aux%lambda(imeasurement), &
                          ierr);CHKERRQ(ierr)
      cycle
    endif
    if (.not.this%measurements(imeasurement)%first_lambda) then
!      call InversionMeasurementPrint(this%measurements(imeasurement),option)
!      print *, inversion_forward_ts_aux%time/(3600.*24.*365.)
      this%measurements(imeasurement)%first_lambda = PETSC_TRUE
      call VecZeroEntries(natural_vec,ierr);CHKERRQ(ierr)
      if (option%myrank == 0) then
        icell_measurement = this%measurements(imeasurement)%cell_id
        tempreal = -1.d0
        call VecSetValue(natural_vec,icell_measurement-1,tempreal, &
                        INSERT_VALUES,ierr);CHKERRQ(ierr)
      endif
      call VecAssemblyBegin(natural_vec,ierr);CHKERRQ(ierr)
      call VecAssemblyEnd(natural_vec,ierr);CHKERRQ(ierr)
      call DiscretizationNaturalToGlobal(discretization,natural_vec, &
                                         onedof_vec,ONEDOF)
      select case(this%iobsfunc)
        case(LIQUID_PRESSURE)
          tempint = zflow_liq_flow_eq
        case(LIQUID_SATURATION)
          tempint = zflow_liq_flow_eq
          call RealizationGetVariable(this%realization,work,DERIVATIVE, &
                                      ZFLOW_LIQ_SAT_WRT_LIQ_PRES)
          call VecPointwiseMult(onedof_vec,onedof_vec,work,ierr);CHKERRQ(ierr)
        case(SOLUTE_CONCENTRATION)
          tempint = zflow_sol_tran_eq
        case default
          option%io_buffer = 'Unknown observation type in InvSubsurfCalcLambda'
          call PrintErrMsg(option)
      end select
      if (Uninitialized(tempint)) then
        option%io_buffer = 'The observed state variable is not being modeled.'
        call PrintErrMsg(option)
      endif
      call VecStrideScatter(onedof_vec,tempint-1,p,INSERT_VALUES,ierr);CHKERRQ(ierr)
      if (this%debug_verbosity > 2) then
        if (OptionPrintToScreen(option)) print *, 'p'
        call VecView(p,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      endif
      call VecZeroEntries(dReskp1_duk_lambdak,ierr);CHKERRQ(ierr)
    else
      call VecZeroEntries(p,ierr);CHKERRQ(ierr)
      call VecZeroEntries(dReskp1_duk_lambdak,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(dReskp1_duk_lambdak,vec_ptr,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(inversion_forward_ts_aux%next%lambda(imeasurement), &
                          vec_ptr2,ierr);CHKERRQ(ierr)
      do local_id = 1, grid%nlmax
        offset = (local_id-1)*ndof
        do i = 1, ndof
          tempreal = 0.d0
          do j = 1, ndof
            tempreal = tempreal + &
              ! this is a transpose block matmult: i,j -> j,i
              inversion_forward_ts_aux%next%dRes_du_k(j,i,local_id)*vec_ptr2(offset+j)
          enddo
          vec_ptr(offset+i) = tempreal
        enddo
      enddo
      call VecRestoreArrayF90(dReskp1_duk_lambdak,vec_ptr,ierr);CHKERRQ(ierr)
      call VecRestoreArrayF90(inversion_forward_ts_aux%next%lambda(imeasurement), &
                              vec_ptr2,ierr);CHKERRQ(ierr)
    endif
    call VecWAXPY(rhs,-1.d0,dReskp1_duk_lambdak,p,ierr);CHKERRQ(ierr)
    call KSPSolveTranspose(solver%ksp,rhs, &
                           inversion_forward_ts_aux%lambda(imeasurement), &
                           ierr);CHKERRQ(ierr)
    if (this%debug_verbosity > 2) then
      if (OptionPrintToScreen(option)) print *, 'lambda'
      call VecView(inversion_forward_ts_aux%lambda(imeasurement), &
                   PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      if (option%comm%mycommsize == 1) then
        call VecGetArrayF90(inversion_forward_ts_aux%lambda(imeasurement), &
                            vec_ptr,ierr);CHKERRQ(ierr)
        print *, vec_ptr(:)
        call VecRestoreArrayF90(inversion_forward_ts_aux%lambda(imeasurement), &
                                vec_ptr,ierr);CHKERRQ(ierr)
      endif
    endif
  enddo
  call VecDestroy(onedof_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(p,ierr);CHKERRQ(ierr)
  call VecDestroy(rhs,ierr);CHKERRQ(ierr)
  call VecDestroy(dReskp1_duk_lambdak,ierr);CHKERRQ(ierr)
  call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)

  call timer%Stop()
  option%io_buffer = '    ' // &
    trim(StringWrite('(f20.1)',timer%GetCumulativeTime()))
  option%io_buffer = trim(option%io_buffer) // &
    ' seconds to calculate lambdas.'
  call PrintMsg(option)
  call TimerDestroy(timer)

end subroutine InvSubsurfCalcLambda

! ************************************************************************** !

subroutine InvSubsurfAddSensitivity(this,inversion_forward_ts_aux)
  !
  ! Calculates sensitivity matrix Jsensitivity
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
  use Solver_module
  use String_module
  use Timer_class

  use PM_Base_class
  use PM_ZFlow_class

  class(inversion_subsurface_type) :: this

  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(inversion_aux_type), pointer :: inversion_aux
  type(inversion_forward_ts_aux_type), pointer :: inversion_forward_ts_aux
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(solver_type), pointer :: solver
  PetscReal, pointer :: vec_ptr(:)
  PetscReal, pointer :: vec_ptr2(:)
  PetscInt :: iparameter, imeasurement
  PetscInt :: natural_id
  PetscInt :: local_id
  PetscInt :: offset
  PetscReal :: tempreal
  Vec :: work
  Vec :: work2
  Vec :: dResdKLambda
  PetscViewer :: viewer
  class(timer_type), pointer :: timer
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr

  nullify(vec_ptr)
  nullify(vec_ptr2)

  dResdKLambda = PETSC_NULL_VEC

  solver => this%forward_simulation%flow_process_model_coupler% &
              timestepper%solver
  option => this%realization%option
  discretization => this%realization%discretization
  patch => this%realization%patch
  grid => patch%grid
  inversion_aux => this%inversion_aux

  timer => TimerCreate()

  call timer%Start()

  work = this%realization%field%work ! DO NOT DESTROY!
  call VecDuplicate(this%realization%field%flow_xx, &
                    dResdKLambda,ierr);CHKERRQ(ierr)

  if (this%debug_adjoint) then
    string = 'dResdK_ts'//trim(StringWrite(inversion_forward_ts_aux%timestep))//'.txt'
    call PetscViewerASCIIOpen(option%mycomm,string, &
                              viewer,ierr);CHKERRQ(ierr)
    call MatView(inversion_forward_ts_aux%dResdparam,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    if (this%debug_verbosity > 2) then
      if (OptionPrintToScreen(option)) print *, 'dResdK'
      call MatView(inversion_forward_ts_aux%dResdparam, &
                   PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
    endif
  endif

  do imeasurement = 1, size(this%measurements)
    if (this%debug_adjoint) then
      string = 'lambda_ts'//trim(StringWrite(inversion_forward_ts_aux%timestep)) // &
                '_' // &
                trim(StringWrite(this%measurements(imeasurement)%cell_id)) // &
                '.txt'
      call PetscViewerASCIIOpen(option%mycomm,string, &
                                viewer,ierr);CHKERRQ(ierr)
      call VecView(inversion_forward_ts_aux%lambda(imeasurement),viewer, &
                    ierr);CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
      if (this%debug_verbosity > 2) then
        if (OptionPrintToScreen(option)) print *, 'lambda ', imeasurement
        call VecView(inversion_forward_ts_aux%lambda(imeasurement), &
                     PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      endif
    endif
    if (this%qoi_is_full_vector) then
      call MatMultTranspose(inversion_forward_ts_aux%dResdparam, &
                            inversion_forward_ts_aux%lambda(imeasurement), &
                            dResdKLambda,ierr);CHKERRQ(ierr)
      if (this%debug_verbosity > 2) then
        if (OptionPrintToScreen(option)) print *, 'dGamdp ', imeasurement
        call VecView(dResdKLambda, &
                    PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
      endif
      call VecGetArrayF90(dResdKLambda,vec_ptr,ierr);CHKERRQ(ierr)
      do iparameter = 1, grid%nlmax
        natural_id = grid%nG2A(grid%nL2G(iparameter))
        offset = (iparameter-1)*option%nflowdof
        call MatSetValue(inversion_aux%JsensitivityT,natural_id-1, &
                         imeasurement-1, &
                         vec_ptr(offset+1),ADD_VALUES,ierr);CHKERRQ(ierr)
      enddo
      call VecRestoreArrayF90(dResdKLambda,vec_ptr,ierr);CHKERRQ(ierr)
    else
      call VecDuplicate(work,work2,ierr);CHKERRQ(ierr)
      do iparameter = 1, size(this%parameters)
        call VecZeroEntries(work,ierr);CHKERRQ(ierr)
        call VecGetArrayF90(work,vec_ptr,ierr);CHKERRQ(ierr)
        do local_id = 1, grid%nlmax
          if (patch%imat(grid%nL2G(local_id)) == &
              this%parameters(iparameter)%imat) then
            vec_ptr(local_id) = 1.d0
          endif
        enddo
        call VecRestoreArrayF90(work,vec_ptr,ierr);CHKERRQ(ierr)
        call MatMult(inversion_forward_ts_aux%dResdparam,work,work2, &
                     ierr);CHKERRQ(ierr)
        call VecDot(work2,inversion_forward_ts_aux%lambda(imeasurement), &
                    tempreal,ierr);CHKERRQ(ierr)
        if (option%comm%myrank == option%driver%io_rank) then
          call MatSetValue(inversion_aux%JsensitivityT,iparameter-1, &
                           imeasurement-1,tempreal,ADD_VALUES,&
                           ierr);CHKERRQ(ierr)
        endif
      enddo
      call VecDestroy(work2,ierr);CHKERRQ(ierr)
    endif
  enddo

  if (dResdKLambda /= PETSC_NULL_VEC) then
    call VecDestroy(dResdKLambda,ierr);CHKERRQ(ierr)
  endif

  call timer%Stop()
  option%io_buffer = '    ' // &
    trim(StringWrite('(f20.1)',timer%GetCumulativeTime()))
  option%io_buffer = trim(option%io_buffer) // &
    ' seconds to add contributions to Jsensitivity.'
  call PrintMsg(option)
  call TimerDestroy(timer)

end subroutine InvSubsurfAddSensitivity

! ************************************************************************** !

subroutine InvSubsurfScaleSensitivity(this)
  !
  ! Scales sensitivity Jacobian for ln(K)
  !
  ! Author: Glenn hammond
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
  ! Author: Glenn hammond
  ! Date: 10/11/21
  !

  class(inversion_subsurface_type) :: this
  character(len=*) :: suffix

  character(len=MAXSTRINGLENGTH) :: filename_prefix

  if (.not.this%print_sensitivity_jacobian) return

  filename_prefix = trim(this%driver%global_prefix) // '_Jsense'
  if (len_trim(suffix) > 0) filename_prefix = trim(filename_prefix) // '_' // &
                            suffix
  call InvSubsurfOutputSensitivityASCII(this,this%inversion_aux%JsensitivityT, &
                                        filename_prefix)
  if (this%qoi_is_full_vector) then
    call InvSubsurfOutputSensitivityHDF5(this,this%inversion_aux%JsensitivityT, &
                                         filename_prefix)
  endif

end subroutine InvSubsurfOutputSensitivity

! ************************************************************************** !

subroutine InvSubsurfOutputSensitivityASCII(this,JsensitivityT,filename_prefix)
  !
  ! Writes sensitivity Jacobian to an ASCII output file
  !
  ! Author: Glenn hammond
  ! Date: 10/11/21
  !
  use Realization_Subsurface_class
  use Variables_module, only : PERMEABILITY

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
  ! Author: Glenn hammond
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
  PetscReal, pointer :: row_ptr(:)
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
    call MatMult(JsensitivityT,row_vec, &
                 this%realization%field%work, &
                 ierr);CHKERRQ(ierr)
    string = 'Measurement ' // &
             StringWrite(this%measurements(imeasurement)%cell_id)
    call HDF5WriteStructDataSetFromVec(string,this%realization, &
                                       this%realization%field%work,grp_id, &
                                       H5T_NATIVE_DOUBLE)
  enddo
  call VecDestroy(row_vec,ierr);CHKERRQ(ierr)
  call h5gclose_f(grp_id,hdf5_err)
  call OutputHDF5CloseFile(this%realization%option,file_id)

end subroutine InvSubsurfOutputSensitivityHDF5

! ************************************************************************** !

subroutine InversionSubsurfaceInvert(this)
  !
  ! Inverts for a parameter update
  !
  ! Author: Glenn hammond
  ! Date: 09/20/21
  !

  class(inversion_subsurface_type) :: this

end subroutine InversionSubsurfaceInvert

! ************************************************************************** !

subroutine InversionSubsurfaceStrip(this)
  !
  ! Deallocates members of inversion Subsurface
  !
  ! Author: Glenn hammond
  ! Date: 09/20/21
  !
  use Utility_module

  class(inversion_subsurface_type) :: this

  PetscInt :: i
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
  if (this%quantity_of_interest /= PETSC_NULL_VEC) then
    call VecDestroy(this%quantity_of_interest,ierr);CHKERRQ(ierr)
  endif
  if (this%ref_quantity_of_interest /= PETSC_NULL_VEC) then
    call VecDestroy(this%ref_quantity_of_interest,ierr);CHKERRQ(ierr)
  endif
  if (this%parameter_vec /= PETSC_NULL_VEC) then
    call VecDestroy(this%parameter_vec,ierr);CHKERRQ(ierr)
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
  if (this%scatter_global_to_measurement /= PETSC_NULL_VECSCATTER) then
    call VecScatterDestroy(this%scatter_global_to_measurement, &
                           ierr);CHKERRQ(ierr)
  endif
  if (this%scatter_measure_to_dist_measure /= PETSC_NULL_VECSCATTER) then
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

end subroutine InversionSubsurfaceStrip

end module Inversion_Subsurface_class
