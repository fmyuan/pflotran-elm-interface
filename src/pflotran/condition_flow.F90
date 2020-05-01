module FlowCondition_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Global_Aux_module
  use Dataset_Base_class
  use Dataset_Ascii_class
  use Time_Storage_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: flow_condition_type
    PetscInt :: id                          ! id from which condition can be referenced
    PetscBool :: is_transient
    PetscBool :: sync_time_with_update
    character(len=MAXWORDLENGTH) :: name    ! name of condition (e.g. initial, recharge)
    PetscInt :: num_sub_conditions
    PetscInt :: nfluids
    PetscInt, pointer :: itype(:)           !  'fluid' sub-condition type (num_sub_conditions)
                                            ! for fluid: 1-pressure, 2-saturation, 3-rate, 4-flux, 5-flow (e.g. H) conductance,
                                            ! for thermal.: 1-temperature, 2-ethalpy, 3-rate, 4-flux, 5-T.conductance
    character(len=MAXWORDLENGTH) :: time_units
    character(len=MAXWORDLENGTH) :: length_units
    type(time_storage_type), pointer  :: default_time_storage
    class(dataset_base_type), pointer :: datum
    ! for a flow-condition, ONLY one of variables below is required, which corresponding to 'itype' above
    type(flow_sub_condition_type), pointer :: pressure   ! nfluids
    type(flow_sub_condition_type), pointer :: saturation ! nfluids
    type(flow_sub_condition_type), pointer :: rate       ! nfluids
    type(flow_sub_condition_type), pointer :: flux       ! nfluids
    type(flow_sub_condition_type), pointer :: conductance! nfluids
    PetscReal, pointer :: molarity(:,:)                  ! (nfluids, nflowspec)
    ! for a flow-condition, ONLY one of variables below is required, which corresponding to 'itype' (T. i.e. thermal) above
    type(flow_sub_condition_type), pointer :: temperature   ! one for all fluids and pored-media
    type(flow_sub_condition_type), pointer :: enthalpy
    type(flow_sub_condition_type), pointer :: energy_rate
    type(flow_sub_condition_type), pointer :: energy_flux
    type(flow_sub_condition_type), pointer :: energy_conductance
    ! any new sub conditions above must be added to FlowConditionIsTransient

    type(sub_condition_ptr_type), pointer :: sub_condition_ptr(:)  ! dim: n(flow)dof, i.e. not-derivatived

    type(flow_condition_type), pointer :: next ! pointer to next condition_type for linked-lists
  end type flow_condition_type

  type, public :: flow_sub_condition_type
    PetscInt :: itype                     ! integer describing type of sub-condition
    PetscInt :: isubtype                  ! if any
    character(len=MAXWORDLENGTH) :: ctype      ! character string describing type of sub-condition
    character(len=MAXWORDLENGTH) :: units      ! units
    character(len=MAXWORDLENGTH) :: name
    class(dataset_base_type), pointer :: gradient  !
    class(dataset_base_type), pointer :: dataset   ! time-series of data with sub-condition
  end type flow_sub_condition_type

  type, public :: sub_condition_ptr_type
    type(flow_sub_condition_type), pointer :: ptr
  end type sub_condition_ptr_type

  type, public :: condition_ptr_type
    type(flow_condition_type), pointer :: ptr
  end type condition_ptr_type

  type, public :: condition_list_type
    PetscInt :: num_conditions
    type(flow_condition_type), pointer :: first
    type(flow_condition_type), pointer :: last
    type(flow_condition_type), pointer :: array(:)
  end type condition_list_type


  public :: FlowConditionCreate, FlowConditionDestroy, FlowConditionRead, &
            FlowConditionAddToList, FlowConditionInitList, &
            FlowConditionDestroyList, &
            FlowConditionGetPtrFromList, FlowConditionUpdate, &
            FlowConditionPrint, &
            FlowConditionIsTransient, &
            ConditionReadValues, &
            GetSubConditionName, &
            FlowConditionUnknownItype

contains

! ************************************************************************** !

function FlowConditionCreate(option)
  !
  ! Creates a flow condition
  !
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! Rewritten by Fengming Yuan @ornl-ccsi/esd, 2020-04/29
  !

  use Option_module

  implicit none

  type(option_type) :: option
  type(flow_condition_type), pointer :: FlowConditionCreate

  type(flow_condition_type), pointer :: condition

  allocate(condition)
  nullify(condition%pressure)
  nullify(condition%saturation)
  nullify(condition%rate)
  nullify(condition%flux)
  nullify(condition%conductance)
  nullify(condition%molarity)
  nullify(condition%temperature)
  nullify(condition%enthalpy)
  nullify(condition%energy_rate)
  nullify(condition%energy_flux)
  nullify(condition%energy_conductance)

  !nullify(condition%sub_condition_ptr)

  nullify(condition%itype)
  nullify(condition%next)
  nullify(condition%datum)
  nullify(condition%default_time_storage)
  condition%is_transient = PETSC_FALSE
  condition%sync_time_with_update = PETSC_FALSE
  condition%time_units = ''
  condition%length_units = ''
  condition%id = 0
  condition%nfluids = 0
  condition%num_sub_conditions = 0
  condition%name = ''

  FlowConditionCreate => condition

end function FlowConditionCreate

! ************************************************************************** !

! ************************************************************************** !

function FlowSubConditionCreate(nfluids)
  !
  ! Creates a sub_condition
  !
  ! Author: Glenn Hammond
  ! Date: 02/04/08
  !

  use Dataset_Ascii_class
  use Option_module

  implicit none

  type(flow_sub_condition_type), pointer :: FlowSubConditionCreate

  PetscInt :: nfluids

  type(flow_sub_condition_type), pointer :: sub_condition
  class(dataset_ascii_type), pointer :: dataset_ascii

  allocate(sub_condition)
  sub_condition%units = ''
  sub_condition%itype = NULL_CONDITION
  sub_condition%isubtype = 0
  sub_condition%ctype = ''
  sub_condition%name = ''
  nullify(sub_condition%gradient)
  nullify(sub_condition%dataset)

  ! by default, all dataset are of type dataset_ascii_type, unless overwritten
  dataset_ascii => DatasetAsciiCreate()
  call DatasetAsciiInit(dataset_ascii)
  dataset_ascii%array_width = nfluids
  dataset_ascii%data_type = DATASET_REAL
  sub_condition%dataset => dataset_ascii
  nullify(dataset_ascii)

  FlowSubConditionCreate => sub_condition

end function FlowSubConditionCreate

! ************************************************************************** !

function GetFlowSubCondFromArrayByName(sub_condition_ptr_list,name)
  !
  ! returns a pointer to a subcondition with
  ! matching name
  !
  ! Author: Glenn Hammond
  ! Date: 06/02/08
  !

  use Input_Aux_module
  use String_module

  implicit none

  type(flow_sub_condition_type), pointer :: GetFlowSubCondFromArrayByName
  type(sub_condition_ptr_type), pointer :: sub_condition_ptr_list(:)
  character(len=MAXWORDLENGTH) :: name

  PetscInt :: idof
  PetscInt :: length

  nullify(GetFlowSubCondFromArrayByName)
  length = len_trim(name)
  do idof = 1, size(sub_condition_ptr_list)
    if (length == len_trim(sub_condition_ptr_list(idof)%ptr%name) .and. &
        StringCompare(name,sub_condition_ptr_list(idof)%ptr%name,length)) then
      GetFlowSubCondFromArrayByName => sub_condition_ptr_list(idof)%ptr
      return
    endif
  enddo

  print *, 'GetFlowSubCondFromArrayByName() needs to be updated to include &
           &the general_condition_type.'
  stop

end function GetFlowSubCondFromArrayByName

! ************************************************************************** !

subroutine FlowSubConditionVerify(option, condition, sub_condition_name, &
                                  sub_condition, default_time_storage, &
                                  destroy_if_null)
  !
  ! Verifies the data in a subcondition
  !
  ! Author: Glenn Hammond
  ! Date: 02/04/08
  !
  use Time_Storage_module
  use Option_module
  use Dataset_module

  implicit none

  type(option_type) :: option
  type(flow_condition_type) :: condition
  character(len=MAXWORDLENGTH) :: sub_condition_name
  type(flow_sub_condition_type), pointer :: sub_condition
  type(time_storage_type), pointer :: default_time_storage
  PetscBool :: destroy_if_null

  character(len=MAXSTRINGLENGTH) :: header

  if (.not.associated(sub_condition)) return

  ! dataset is not optional
  if (.not.(associated(sub_condition%dataset%rarray) .or. &
            associated(sub_condition%dataset%rbuffer) .or. &
            len_trim(sub_condition%dataset%name) > 0 .or. &
            sub_condition%itype /= NULL_CONDITION)) then
    if (destroy_if_null) call FlowSubConditionDestroy(sub_condition)
    return
  endif

  if (len_trim(sub_condition%ctype) == NULL_CONDITION) then
    option%io_buffer = 'TYPE of condition ' // trim(condition%name) // &
      ' ' // trim(sub_condition_name) // ' dataset not defined.'
    call PrintErrMsg(option)
  endif

  header = 'SUBSURFACE/FLOW_CONDITION/' // &
           trim(condition%name) // '/' // &
           trim(sub_condition_name) // '/Value(s)'
  call DatasetVerify(sub_condition%dataset,default_time_storage, &
                     header,option)
  header = 'SUBSURFACE/FLOW_CONDITION/' // &
           trim(condition%name) // '/' // &
           trim(sub_condition_name) // '/Gradient'
  call DatasetVerify(sub_condition%gradient,default_time_storage, &
                     header,option)

end subroutine FlowSubConditionVerify

! ************************************************************************** !

subroutine FlowConditionRead(condition,input,option)
  !
  ! Reads a condition from the input file
  !
  ! Author: Glenn Hammond
  ! Date: 10/31/07
  !

  use Option_module
  use Input_Aux_module
  use String_module
  use Logging_module
  use Time_Storage_module
  use Dataset_module

  implicit none

  type(flow_condition_type) :: condition
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: rate_unit_string
  character(len=MAXWORDLENGTH) :: energy_rate_unit_string
  character(len=MAXWORDLENGTH) :: internal_units
  type(flow_sub_condition_type), pointer :: pressure, saturation, &
                                            flux, rate,           &
                                            conductance,          &
                                            temperature, enthalpy,     &
                                            energy_flux,  energy_rate, &
                                            energy_conductance,        &
                                            sub_condition_ptr
  PetscReal :: default_time
  PetscInt :: idof
  type(time_storage_type), pointer :: default_time_storage
  class(dataset_ascii_type), pointer :: dataset_ascii

  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_flow_condition_read, &
                          ierr);CHKERRQ(ierr)

  default_time = 0.d0

  default_time_storage => TimeStorageCreate()
  default_time_storage%is_cyclic = PETSC_FALSE
  default_time_storage%time_interpolation_method = INTERPOLATION_STEP

  rate_unit_string = 'not_assigned'
  energy_rate_unit_string = 'not_assigned'
  internal_units = 'not_assigned'

  pressure => FlowSubConditionCreate(option%flow%nfluid)
  pressure%name = 'pressure'
  saturation => FlowSubConditionCreate(option%flow%nfluid)
  saturation%name = 'saturation'
  rate => FlowSubConditionCreate(option%flow%nfluid)
  rate%name = 'rate'
  flux => FlowSubConditionCreate(option%flow%nfluid)
  flux%name = 'flux'
  conductance => FlowSubConditionCreate(option%flow%nfluid)
  conductance%name = 'conductance'

  temperature => FlowSubConditionCreate(ONE_INTEGER)
  temperature%name = 'temperature'
  enthalpy => FlowSubConditionCreate(ONE_INTEGER)
  enthalpy%name = 'enthalpy'
  energy_rate => FlowSubConditionCreate(ONE_INTEGER)
  energy_rate%name = 'energy_rate'
  energy_flux => FlowSubConditionCreate(ONE_INTEGER)
  energy_flux%name = 'energy_flux'
  energy_conductance => FlowSubConditionCreate(ONE_INTEGER)
  energy_conductance%name = 'energy_conductance'

  condition%time_units = 'yr'
  condition%length_units = 'm'
  pressure%units = 'Pa'
  saturation%units = ' '
  rate%units = 'kg/s'
  flux%units = 'm/s'
  conductance%units = 'm/s/Pa'
  temperature%units = 'C'
  enthalpy%units = 'kJ/mol'
  energy_rate%units = 'W'
  energy_flux%units = 'W/m^2'
  energy_conductance%units = 'W/m^2/C'

  ! read the condition
  input%ierr = 0
  call InputPushBlock(input,option)
  do
  
    internal_units = 'not_assigned'

    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,'CONDITION')

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword','CONDITION')

    select case(trim(word))

      case('UNITS') ! read default units for condition arguments
        do
          call InputReadWord(input,option,word,PETSC_TRUE)
          if (InputError(input)) exit
          select case(trim(word))
            case('s','sec','min','hr','d','day','y','yr')
              condition%time_units = trim(word)
            case('mm','cm','m','met','meter','dm','km')
              condition%length_units = trim(word)
            case('Pa','KPa')
              pressure%units = trim(word)
            case('kg/s','kg/yr')
              rate%units = trim(word)
            case('m/s','m/yr')
              flux%units = trim(word)
            case('m/s/Pa','m/s/KPa')
              conductance%units = trim(word)
            case('W','J/yr')
              energy_rate%units = trim(word)
            case('W/m^2','J/m^2/yr')
              energy_flux%units = trim(word)
            case('W/m^2/C','J/m^2/yr/C')
              energy_conductance%units = trim(word)
            case('C','K')
              temperature%units = trim(word)
            case('kJ/mol')
              enthalpy%units = trim(word)
            case default
              call InputKeywordUnrecognized(input,word,'condition,units', &
                                            option)
          end select
        enddo
      case('CYCLIC')
        ! by default, is_cyclic is set to PETSC_FALSE
        default_time_storage%is_cyclic = PETSC_TRUE
      case('SYNC_TIMESTEP_WITH_UPDATE')
        condition%sync_time_with_update = PETSC_TRUE
      case('INTERPOLATION')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'INTERPOLATION','CONDITION')
        call StringToLower(word)
        select case(word)
          case('step')
            default_time_storage%time_interpolation_method = &
              INTERPOLATION_STEP
          case('linear')
            default_time_storage%time_interpolation_method = &
              INTERPOLATION_LINEAR
          case default
            call InputKeywordUnrecognized(input,word, &
                                          'condition,interpolation', &
                                          option)
        end select
      case('TYPE') ! read condition type (dirichlet, neumann, etc) for each dof
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')

          if (InputCheckExit(input,option)) exit

          if (InputError(input)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')
          call StringToUpper(word)
          select case(trim(word))
            case('PRESSURE')
              sub_condition_ptr => pressure
              sub_condition_ptr%itype = DIRICHLET_BC ! default for 'PRESSURE/SATURATTION'
            case('SATURATION')
              sub_condition_ptr => saturation
              sub_condition_ptr%itype = DIRICHLET_BC ! default for 'PRESSURE/SATUTRATION'
            case('RATE')
              sub_condition_ptr => rate
            case('FLUX')
              sub_condition_ptr => flux
              sub_condition_ptr%itype = NEUMANN_BC   ! default for 'FLUX'
            case('CONDUCTANCE')
              sub_condition_ptr => conductance
            case('TEMPERATURE')
              sub_condition_ptr => temperature
              sub_condition_ptr%itype = DIRICHLET_BC ! default for 'TEMPERATURE/ENTHALPY'
            case('ENTHALPY')
              sub_condition_ptr => enthalpy
              sub_condition_ptr%itype = DIRICHLET_BC ! default for 'TEMPERATURE/ENTHALPY'
            case('ENERGY_RATE')
              sub_condition_ptr => energy_rate
            case('ENERGY_FLUX')
              sub_condition_ptr => energy_flux
              sub_condition_ptr%itype = NEUMANN_BC   ! default for 'FLUX'
            case('ENERGY_CONDUCTANCE')
              sub_condition_ptr => energy_conductance
            case default
              call InputKeywordUnrecognized(input,word,'condition,type',option)
          end select
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'TYPE','CONDITION')
          call StringToLower(word)
          sub_condition_ptr%ctype = word
          select case(word)
            case('dirichlet')
              sub_condition_ptr%itype = DIRICHLET_BC
            case('neumann')
              sub_condition_ptr%itype = NEUMANN_BC
            case('mass_rate')
              sub_condition_ptr%itype = MASS_RATE_SS
              rate_unit_string = 'kg/sec'
            case('energy_rate')
              sub_condition_ptr%itype = ENERGY_RATE_SS
              energy_rate_unit_string = 'MJ/sec|MW'
            case('heterogeneous_energy_rate')
              sub_condition_ptr%itype = HET_ENERGY_RATE_SS
              energy_rate_unit_string = 'MJ/sec|MW'
            case('scaled_mass_rate','scaled_volumetric_rate', &
                 'scaled_energy_rate')
              select case(word)
                case('scaled_mass_rate')
                  sub_condition_ptr%itype = SCALED_MASS_RATE_SS
                  rate_unit_string = 'kg/sec'
                case('scaled_volumetric_rate')
                  sub_condition_ptr%itype = SCALED_VOLUMETRIC_RATE_SS
                  rate_unit_string = 'm^3/sec'
                case('scaled_energy_rate')
                  sub_condition_ptr%itype = SCALED_ENERGY_RATE_SS
                  energy_rate_unit_string = 'MW|MJ/sec'
              end select

              ! store name of type for error messaging below.
              string = word
              call InputReadWord(input,option,word,PETSC_TRUE)
              if (input%ierr == 0) then
                call InputPushCard(input,word,option)
                call StringToLower(word)
                sub_condition_ptr%ctype = trim(sub_condition_ptr%ctype) // word
                select case(word)
                  case('neighbor_perm')
                    sub_condition_ptr%isubtype = SCALE_BY_NEIGHBOR_PERM
                  case('volume')
                    sub_condition_ptr%isubtype = SCALE_BY_VOLUME
                  case('perm')
                    sub_condition_ptr%isubtype = SCALE_BY_PERM
                  case default
                    string = 'flow condition "' // trim(condition%name) // &
                      '" ' // trim(string)
                    call InputKeywordUnrecognized(input,word,string,option)
                end select
              else
                option%io_buffer = 'Specify one of NEIGHBOR_PERM, &
                  &VOLUME, PERM subtypes in flow condition "' // &
                  trim(condition%name) // '" ' // trim(string)
                call PrintErrMsg(option)
                endif
            case('hydrostatic')
              sub_condition_ptr%itype = HYDROSTATIC_BC
            case('conductance')
              sub_condition_ptr%itype = HYDROSTATIC_CONDUCTANCE_BC
            case('zero_gradient')
              sub_condition_ptr%itype = ZERO_GRADIENT_BC
            case('seepage')
              sub_condition_ptr%itype = HYDROSTATIC_SEEPAGE_BC
            case('dirichlet_seepage')
              sub_condition_ptr%itype = DIRICHLET_SEEPAGE_BC
            case('dirichlet_conductance')
              sub_condition_ptr%itype = DIRICHLET_CONDUCTANCE_BC
            case('volumetric_rate')
              sub_condition_ptr%itype = VOLUMETRIC_RATE_SS
              rate_unit_string = 'm^3/sec'
            case('equilibrium')
              sub_condition_ptr%itype = EQUILIBRIUM_SS
            case('unit_gradient')
              if (.not.associated(sub_condition_ptr,pressure)) then
                option%io_buffer = 'unit_gradient flow condition type may &
                  &only be associated with a PRESSURE flow condition.'
                call PrintErrMsg(option)
              endif
              sub_condition_ptr%itype = UNIT_GRADIENT_BC
            case('heterogeneous_volumetric_rate')
              sub_condition_ptr%itype = HET_VOL_RATE_SS
              rate_unit_string = 'm^3/sec'
            case('heterogeneous_mass_rate')
              sub_condition_ptr%itype = HET_MASS_RATE_SS
              rate_unit_string = 'kg/sec'
            case('heterogeneous_dirichlet')
              sub_condition_ptr%itype = HET_DIRICHLET_BC
            case('heterogeneous_seepage')
              sub_condition_ptr%itype = HET_HYDROSTATIC_SEEPAGE_BC
            case('heterogeneous_conductance')
              sub_condition_ptr%itype = HET_HYDROSTATIC_CONDUCTANCE_BC
            case('heterogeneous_surface_seepage')
              sub_condition_ptr%itype = HET_SURF_HYDROSTATIC_SEEPAGE_BC
            case('spillover')
              sub_condition_ptr%itype = SPILLOVER_BC
            case default
              call InputKeywordUnrecognized(input,word,'condition bc type',option)
          end select
        enddo
        call InputPopBlock(input,option)
      case('DATUM')
        dataset_ascii => DatasetAsciiCreate()
        call DatasetAsciiInit(dataset_ascii)
        dataset_ascii%array_width = 3
        dataset_ascii%data_type = DATASET_REAL
        condition%datum => dataset_ascii
        nullify(dataset_ascii)
        internal_units = 'meter'
        call ConditionReadValues(input,option,word, &
                                 condition%datum,word,internal_units)
      case('GRADIENT','GRAD')
        call InputPushBlock(input,option)
        do
          internal_units = 'not_assigned'
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'CONDITION')

          if (InputCheckExit(input,option)) exit

          if (InputError(input)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword','CONDITION,TYPE')
          select case(trim(word))
            case('PRES','PRESS','PRESSURE')
              sub_condition_ptr => pressure
              internal_units = 'Pa/meter'
            case('RATE')
              sub_condition_ptr => rate
              internal_units = 'kg/sec-meter'
            case('ENERGY_RATE')
              sub_condition_ptr => energy_rate
              internal_units = 'MW/meter|MJ/sec-meter'
            case('FLUX')
              sub_condition_ptr => flux
              internal_units = 'm/sec-m|unitless/sec'
            case('SATURATION')
              sub_condition_ptr => saturation
              internal_units = 'unitless/meter'
            case('TEMP','TEMPERATURE')
              sub_condition_ptr => temperature
              internal_units = 'temperature/m'
            case('H','ENTHALPY')
              sub_condition_ptr => enthalpy
              internal_units = 'kJ/mol-meter'
            case default
              call InputKeywordUnrecognized(input,word, &
                     'FLOW CONDITION,GRADIENT,TYPE',option)
          end select
          dataset_ascii => DatasetAsciiCreate()
          call DatasetAsciiInit(dataset_ascii)
          dataset_ascii%array_width = 3
          dataset_ascii%data_type = DATASET_REAL
          sub_condition_ptr%gradient => dataset_ascii
          nullify(dataset_ascii)
          call ConditionReadValues(input,option,word, &
                                   sub_condition_ptr%gradient, &
                                   word,internal_units)
          nullify(sub_condition_ptr)
        enddo
        call InputPopBlock(input,option)
      case('TEMPERATURE','TEMP')
        internal_units = 'C'
        call ConditionReadValues(input,option,word, &
                                 temperature%dataset, &
                                 temperature%units,internal_units)
      case('ENTHALPY','H')
        internal_units = 'kJ/mol'
        call ConditionReadValues(input,option,word, &
                                 enthalpy%dataset, &
                                 enthalpy%units,internal_units)
      case('PRESSURE','PRES','PRESS')
        internal_units = 'Pa'
        call ConditionReadValues(input,option,word, &
                                 pressure%dataset, &
                                 pressure%units,internal_units)
      case('RATE')
        internal_units = rate_unit_string
        call ConditionReadValues(input,option,word, &
                                 rate%dataset, &
                                 rate%units,internal_units)
      case('ENERGY_FLUX')
        input%force_units = PETSC_TRUE
        internal_units = 'MW/m^2|MJ/m^2-sec'
        call ConditionReadValues(input,option,word, &
                                 energy_flux%dataset, &
                                 energy_flux%units,internal_units)
        input%force_units = PETSC_FALSE
      case('ENERGY_CONDUCTANCE')
        input%force_units = PETSC_TRUE
        internal_units = 'W/m^2/C'
        call ConditionReadValues(input,option,word, &
                                 energy_conductance%dataset, &
                                 energy_conductance%units,internal_units)
        input%force_units = PETSC_FALSE
      case('ENERGY_RATE')
        input%force_units = PETSC_TRUE
        internal_units = energy_rate_unit_string
        input%err_buf = word
        call ConditionReadValues(input,option,word, &
                                 energy_rate%dataset, &
                                 energy_rate%units,internal_units)
        input%force_units = PETSC_FALSE
      case('FLUX','VELOCITY','VEL')
        internal_units = 'meter/sec'
        call ConditionReadValues(input,option,word, &
                                 flux%dataset, &
                                 flux%units,internal_units)
      case('SAT','SATURATION')
        internal_units = 'unitless'
        call ConditionReadValues(input,option,word, &
                                 saturation%dataset, &
                                 saturation%units,internal_units)
      case('CONDUCTANCE')
        internal_units = 'meter/sec/Pa'
        call ConditionReadValues(input,option,word, &
                                 conductance%dataset, &
                                 conductance%units,internal_units)
      case default
        call InputKeywordUnrecognized(input,word,'flow condition',option)
    end select

  enddo
  call InputPopBlock(input,option)

  ! datum is not required
  string = trim(condition%name) // '/' // 'Datum'
  call DatasetVerify(condition%datum,default_time_storage,string,option)

  ! check to ensure that a rate condition is not of type pressure
  if (associated(rate)) then
    select case(rate%itype)
      case(DIRICHLET_BC,NEUMANN_BC,HYDROSTATIC_BC,UNIT_GRADIENT_BC, &
           HYDROSTATIC_CONDUCTANCE_BC,ZERO_GRADIENT_BC,HYDROSTATIC_SEEPAGE_BC, &
           DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC,SURFACE_DIRICHLET, &
           SURFACE_SPILLOVER,HET_DIRICHLET_BC,HET_HYDROSTATIC_SEEPAGE_BC,&
           HET_HYDROSTATIC_CONDUCTANCE_BC)
        option%io_buffer = 'RATE condition must not be of type: dirichlet, &
          &neumann, zero_gradient, dirichlet_zero_gradient, hydrostatic, &
          &seepage, or conductance".'
        call PrintErrMsg(option)
    end select
  endif
  ! check to ensure that a pressure condition is not of type rate
  if (associated(pressure)) then
    select case(pressure%itype)
      case(MASS_RATE_SS,SCALED_MASS_RATE_SS,VOLUMETRIC_RATE_SS, &
           SCALED_VOLUMETRIC_RATE_SS,EQUILIBRIUM_SS)
        option%io_buffer = 'PRESSURE or FLUX condition must not be of type: &
          &mass_rate, scaled_mass_rate, volumetric_rate, &
          &scaled_volumetric_rate, equilibrium, or production_well.'
        call PrintErrMsg(option)
    end select
  endif

  ! verify the datasets
  word = 'pressure'
  call FlowSubConditionVerify(option,condition,word,pressure, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'saturation'
  call FlowSubConditionVerify(option,condition,word,saturation, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'flux'
  call FlowSubConditionVerify(option,condition,word,flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'rate'
  call FlowSubConditionVerify(option,condition,word,rate, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'energy_flux'
  call FlowSubConditionVerify(option,condition,word,energy_flux, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'energy_rate'
  call FlowSubConditionVerify(option,condition,word,energy_rate, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'temperature'
  call FlowSubConditionVerify(option,condition,word,temperature, &
                              default_time_storage, &
                              PETSC_TRUE)
  word = 'enthalpy'
  call FlowSubConditionVerify(option,condition,word,enthalpy, &
                              default_time_storage, &
                              PETSC_TRUE)

  select case(option%iflowmode)
    case default
      option%io_buffer = 'The flow mode not supported in original &
        &FlowConditionRead.'
      call PrintMsg(option)

    case(MPFLOW_MODE)
      if (.not.associated(pressure) &
           .and. .not.associated(saturation) &
           .and. .not.associated(flux) &
           .and. .not.associated(rate)) then
        option%io_buffer = 'at least one of pressure, saturation, conductance, &
          &flux, or rate condition NOT null in &
          &condition: ' // trim(condition%name)
        call PrintErrMsg(option)
      endif

      if (associated(pressure)) then
        condition%pressure => pressure
      endif
      if (associated(saturation)) then
        condition%saturation => saturation
      endif
      if (associated(flux)) then
        condition%flux => flux
      endif
      if (associated(rate)) then
        condition%rate => rate
      endif

      if (.not.associated(temperature) .and. .not.associated(energy_rate) &
          .and. .not.associated(energy_flux)) then
        option%io_buffer = 'temperature, energy_flux, and energy_rate &
          &condition null in condition: ' // trim(condition%name)
        call PrintErrMsg(option)
      endif
      if (associated(temperature) .and. associated(energy_rate) ) then
        option%io_buffer = 'Both, temperature and energy_rate cannot be &
                           &specified in condition: ' // trim(condition%name)
        call PrintErrMsg(option)
      endif
      if (associated(temperature)) condition%temperature => temperature
      if (associated(enthalpy)) condition%enthalpy => enthalpy
      if (associated(energy_flux)) condition%energy_flux => energy_flux
      if (associated(energy_rate)) condition%energy_rate => energy_rate
      if (associated(energy_conductance)) condition%energy_conductance => energy_conductance

      condition%num_sub_conditions = option%nflowdof
      allocate(condition%sub_condition_ptr(condition%num_sub_conditions))
      do idof = 1, condition%num_sub_conditions
        nullify(condition%sub_condition_ptr(idof)%ptr)
      enddo

      ! must be in this order, which matches the dofs i problem
      if (associated(pressure)) &
        condition%sub_condition_ptr(ONE_INTEGER)%ptr => pressure
      if (associated(flux)) &
        condition%sub_condition_ptr(ONE_INTEGER)%ptr => flux
      if (associated(rate)) &
        condition%sub_condition_ptr(ONE_INTEGER)%ptr => rate
      if (associated(saturation)) &
        condition%sub_condition_ptr(ONE_INTEGER)%ptr => saturation

      if ( associated(temperature)) &
        condition%sub_condition_ptr(TWO_INTEGER)%ptr => temperature
      if (associated(energy_flux)) &
        condition%sub_condition_ptr(TWO_INTEGER)%ptr => energy_flux
      if (associated(energy_rate)) &
        condition%sub_condition_ptr(TWO_INTEGER)%ptr => energy_rate

      allocate(condition%itype(condition%num_sub_conditions))
      condition%itype = 0
      if (associated(pressure)) &
        condition%itype(ONE_INTEGER) = pressure%itype
      if (associated(flux)) &
        condition%itype(ONE_INTEGER) = flux%itype
      if (associated(rate)) &
        condition%itype(ONE_INTEGER) = rate%itype
      if (associated(saturation)) &
        condition%itype(ONE_INTEGER) = saturation%itype

      if (associated(temperature)) &
        condition%itype(TWO_INTEGER) = temperature%itype
      if (associated(energy_flux)) &
        condition%itype(TWO_INTEGER) = energy_flux%itype
      if (associated(energy_rate)) &
        condition%itype(TWO_INTEGER) = energy_rate%itype

  end select

  condition%default_time_storage => default_time_storage

  call PetscLogEventEnd(logging%event_flow_condition_read,ierr);CHKERRQ(ierr)

end subroutine FlowConditionRead

! ************************************************************************** !

subroutine ConditionReadValues(input,option,keyword,dataset_base, &
                               data_external_units,data_internal_units)
  ! 
  ! Read the value(s) of a condition variable
  !
  ! Author: Glenn Hammond
  ! Date: 10/31/07
  !

  use Input_Aux_module
  use String_module
  use Option_module
  use Logging_module
  use HDF5_Aux_module
  use Units_module
  use Dataset_module
  use Dataset_Base_class
  use Dataset_Ascii_class
#if defined(PETSC_HAVE_HDF5)
  use hdf5
#endif

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  class(dataset_base_type), pointer :: dataset_base
  character(len=*) :: data_external_units
  character(len=*) :: data_internal_units

  character(len=MAXSTRINGLENGTH), pointer :: internal_unit_strings(:)
  class(dataset_ascii_type), pointer :: dataset_ascii
  character(len=MAXSTRINGLENGTH) :: string2, filename, hdf5_path
  character(len=MAXWORDLENGTH) :: word, realization_word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscInt :: length, i
  PetscErrorCode :: ierr

  call PetscLogEventBegin(logging%event_flow_condition_read_values, &
                          ierr);CHKERRQ(ierr)

  ! dataset_base, though of type dataset_base_type, should always be created
  ! as dataset_ascii_type.
  dataset_ascii => DatasetAsciiCast(dataset_base)
  if (.not.associated(dataset_ascii)) then
    ! The dataset was not of type dataset_asci and was likely set to a different
    ! type.  There is a bug in the input file.
    option%io_buffer = 'Dataset associated with ' // trim(keyword) // &
      ' in the input file is already associated with a different dataset &
      &type.  Check for duplicate definitions of ' // trim(keyword) // '.'
    call PrintErrMsg(option)
  endif

  filename = ''
  realization_word = ''
  hdf5_path = ''

  internal_unit_strings => StringSplit(data_internal_units,',')

  input%ierr = 0
  string2 = trim(input%buf)
  call InputReadWord(input,option,word,PETSC_TRUE)
  call InputErrorMsg(input,option,'file or value','CONDITION')
  call StringToLower(word)
  length = len_trim(word)
  if (StringStartsWithAlpha(word)) then
    call InputPushCard(input,word,option)
    if (length == FOUR_INTEGER .and. &
        StringCompare(word,'file',FOUR_INTEGER)) then 
      input%err_buf2 = trim(keyword) // ', FILE'
      input%err_buf = 'keyword'
      call InputReadNChars(input,option,string2,MAXSTRINGLENGTH,PETSC_TRUE)
      if (input%ierr == 0) then
        filename = string2
      else
        option%io_buffer = 'The ability to read realization dependent &
          &datasets outside the DATASET block is no longer supported'
        call PrintErrMsg(option)
      endif

      if (len_trim(filename) < 2) then
        option%io_buffer = 'No filename listed under Flow_Condition: ' // &
                           trim(keyword)
        call PrintErrMsg(option)
      endif

      if (index(filename,'.h5') > 0) then
        write(option%io_buffer,'("Reading of HDF5 datasets for flow ", &
                                 &"conditions not currently supported.")')
        call PrintErrMsg(option)
      else
        i = index(filename,'.',PETSC_TRUE)
        if (i > 2) then
          filename = filename(1:i-1) // trim(realization_word) // filename(i:)
        else
          filename = trim(filename) // trim(realization_word)
        endif
        error_string = 'CONDITION,' // trim(keyword) // ',FILE'
        call DatasetAsciiReadFile(dataset_ascii,filename,data_external_units, &
                                  data_internal_units,error_string,option)
        dataset_ascii%filename = filename
      endif
    else if (StringCompare(word,'dataset')) then
      call InputReadWord(input,option,word,PETSC_TRUE)
      input%err_buf2 = trim(keyword) // ', DATASET'
      input%err_buf = 'dataset name'
      call InputErrorMsg(input,option)
      call DatasetDestroy(dataset_base)
      dataset_base => DatasetBaseCreate()
      dataset_base%name = word
    else if (length==FOUR_INTEGER .and. StringCompare(word,'list',length)) then 
      error_string = 'CONDITION,' // trim(keyword) // ',LIST'
      call DatasetAsciiReadList(dataset_ascii,input,data_external_units, &
                                data_internal_units,error_string,option)
    else if (StringCompare(word,'dbase_value')) then
      input%buf = trim(string2)
      error_string = 'CONDITION,' // trim(keyword) // ',SINGLE'
      call DatasetAsciiReadSingle(dataset_ascii,input,data_external_units, &
                                  data_internal_units,error_string,option)
    else
      option%io_buffer = 'Keyword "' // trim(word) // &
        '" not recognized in when reading condition values for "' // &
        trim(keyword) // '".'
      call PrintErrMsg(option)
    endif
  else
    input%buf = trim(string2)
    error_string = 'CONDITION,' // trim(keyword) // ',SINGLE'
    call DatasetAsciiReadSingle(dataset_ascii,input,data_external_units, &
                                data_internal_units,error_string,option)
  endif

  deallocate(internal_unit_strings)
  nullify(internal_unit_strings)

  call PetscLogEventEnd(logging%event_flow_condition_read_values, &
                        ierr);CHKERRQ(ierr)

end subroutine ConditionReadValues

! ************************************************************************** !

subroutine FlowConditionPrint(condition,option)
  !
  ! Prints flow condition info
  !
  ! Author: Glenn Hammond
  ! Date: 12/04/08
  !

  use Option_module
  use Dataset_module

  implicit none

  type(flow_condition_type) :: condition
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i

99 format(/,80('-'))

  write(option%fid_out,'(/,2x,''Flow Condition: '',a)') trim(condition%name)

  if (condition%sync_time_with_update) then
    string = 'yes'
  else
    string = 'no'
  endif
  write(option%fid_out,'(4x,''Synchronize time with update: '', a)') trim(string)
  write(option%fid_out,'(4x,''Time units: '', a)') trim(condition%time_units)
  write(option%fid_out,'(4x,''Length units: '', a)') trim(condition%length_units)

100 format(6x,a)
  write(option%fid_out,100) 'Datum:'
  if (associated(condition%datum)) then
    call DatasetPrint(condition%datum,option)
  endif

  do i=1, condition%num_sub_conditions
    call FlowConditionPrintSubCondition(condition%sub_condition_ptr(i)%ptr, &
                                        option)
  enddo
  write(option%fid_out,99)

end subroutine FlowConditionPrint

! ************************************************************************** !

subroutine FlowConditionPrintSubCondition(subcondition,option)
  !
  ! Prints flow subcondition info
  !
  ! Author: Glenn Hammond
  ! Date: 12/04/08
  !

  use Option_module
  use Dataset_module

  implicit none

  type(flow_sub_condition_type) :: subcondition
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  write(option%fid_out,'(/,4x,''Sub Condition: '',a)') trim(subcondition%name)
  string = GetSubConditionName(subcondition%itype)

  105 format(6x,'Type: ',a)
  write(option%fid_out,105) trim(string)

  110 format(6x,a)

  write(option%fid_out,110) 'Gradient:'
  if (associated(subcondition%gradient)) then
    call DatasetPrint(subcondition%gradient,option)
  endif

  write(option%fid_out,110) 'Data:'
  if (associated(subcondition%dataset)) then
    call DatasetPrint(subcondition%dataset,option)
  endif

end subroutine FlowConditionPrintSubCondition

! ************************************************************************** !

function GetSubConditionName(subcon_itype)
  !
  ! SubConditionName: Return name of subcondition
  !
  ! Author: Gautam Bisht
  ! Date: 10/16/13
  !

  implicit none

  PetscInt :: subcon_itype

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: GetSubConditionName

  select case(subcon_itype)
    case(DIRICHLET_BC)
      string = 'dirichlet'
    case(NEUMANN_BC)
      string = 'neumann'
    case(DIRICHLET_ZERO_GRADIENT_BC)
      string = 'dirichlet-zero gradient'
    case(MASS_RATE_SS)
      string = 'mass_rate'
    case(HYDROSTATIC_BC)
      string = 'hydrostatic'
    case(HYDROSTATIC_CONDUCTANCE_BC)
      string = 'conductance'
    case(ZERO_GRADIENT_BC)
      string = 'zero gradient'
    case(HYDROSTATIC_SEEPAGE_BC)
      string = 'seepage'
    case(DIRICHLET_SEEPAGE_BC)
      string = 'dirichlet seepage'
    case(DIRICHLET_CONDUCTANCE_BC)
      string = 'dirichlet conductance'
    case(VOLUMETRIC_RATE_SS)
      string = 'volumetric rate'
    case(EQUILIBRIUM_SS)
      string = 'equilibrium'
    case(UNIT_GRADIENT_BC)
      string = 'unit gradient'
    case(SCALED_MASS_RATE_SS)
      string = 'scaled mass rate'
    case(SCALED_VOLUMETRIC_RATE_SS)
      string = 'scaled volumetric rate'
    case(HET_VOL_RATE_SS)
      string = 'heterogeneous volumetric rate'
    case(HET_MASS_RATE_SS)
      string = 'heterogeneous mass rate'
    case(HET_DIRICHLET_BC)
      string = 'heterogeneous dirichlet'
    case(HET_HYDROSTATIC_SEEPAGE_BC)
      string = 'heterogeneous seepage'
    case(HET_HYDROSTATIC_CONDUCTANCE_BC)
      string = 'heterogeneous conductance'
    case(ENERGY_RATE_SS)
      string = 'energy rate'
    case(SCALED_ENERGY_RATE_SS)
      string = 'scaled energy rate'
    case(HET_ENERGY_RATE_SS)
      string = 'heterogeneous energy rate'
    case(HET_SURF_HYDROSTATIC_SEEPAGE_BC)
      string = 'heterogeneous surface seepage'
    case(SPILLOVER_BC)
      string = 'spillover'
    case(SURFACE_DIRICHLET)
      string = 'surface_dirichlet'
    case(SURFACE_ZERO_GRADHEIGHT)
      string = 'surface_zero_gradheight'
    case(SURFACE_SPILLOVER)
      string = 'surface_spillover'
    end select

  GetSubConditionName = trim(string)

end function GetSubConditionName

! ************************************************************************** !

subroutine FlowConditionUpdate(condition_list,option)
  ! 
  ! Updates a transient condition
  !
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  !

  use Option_module
  use Dataset_module

  implicit none

  type(condition_list_type) :: condition_list
  type(option_type) :: option

  type(flow_condition_type), pointer :: condition
  type(flow_sub_condition_type), pointer :: sub_condition
  PetscInt :: isub_condition

  condition => condition_list%first
  do
    if (.not.associated(condition)) exit
    call DatasetUpdate(condition%datum,option)
    do isub_condition = 1, condition%num_sub_conditions

      sub_condition => condition%sub_condition_ptr(isub_condition)%ptr

      if (associated(sub_condition)) then
        call DatasetUpdate(sub_condition%dataset,option)
        call DatasetUpdate(sub_condition%gradient,option)
      endif

    enddo

    condition => condition%next

  enddo

end subroutine FlowConditionUpdate

! ************************************************************************** !

! ************************************************************************** !

subroutine FlowConditionInitList(list)
  !
  ! Initializes a condition list
  !
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  !

  implicit none

  type(condition_list_type) :: list

  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_conditions = 0

end subroutine FlowConditionInitList

! ************************************************************************** !

subroutine FlowConditionAddToList(new_condition,list)
  !
  ! Adds a new condition to a condition list
  !
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  !

  implicit none

  type(flow_condition_type), pointer :: new_condition
  type(condition_list_type) :: list

  list%num_conditions = list%num_conditions + 1
  new_condition%id = list%num_conditions
  if (.not.associated(list%first)) list%first => new_condition
  if (associated(list%last)) list%last%next => new_condition
  list%last => new_condition

end subroutine FlowConditionAddToList

! ************************************************************************** !

function FlowConditionGetPtrFromList(condition_name,condition_list)
  !
  ! Returns a pointer to the condition matching &
  ! condition_name
  !
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  !

  use String_module

  implicit none

  type(flow_condition_type), pointer :: FlowConditionGetPtrFromList
  character(len=MAXWORDLENGTH) :: condition_name
  type(condition_list_type) :: condition_list

  PetscInt :: length
  type(flow_condition_type), pointer :: condition

  nullify(FlowConditionGetPtrFromList)
  condition => condition_list%first

  do
    if (.not.associated(condition)) exit
    length = len_trim(condition_name)
    if (length == len_trim(condition%name) .and. &
        StringCompare(condition%name,condition_name, &
                      length)) then
      FlowConditionGetPtrFromList => condition
      return
    endif
    condition => condition%next
  enddo

end function FlowConditionGetPtrFromList

! ************************************************************************** !

function FlowConditionIsTransient(condition)
  !
  ! Returns PETSC_TRUE
  !
  ! Author: Glenn Hammond
  ! Date: 01/12/11
  !

  use Dataset_module

  implicit none

  type(flow_condition_type) :: condition

  PetscBool :: FlowConditionIsTransient

  FlowConditionIsTransient = PETSC_FALSE

  if (DatasetIsTransient(condition%datum) .or. &
      FlowSubConditionIsTransient(condition%pressure) .or. &
      FlowSubConditionIsTransient(condition%saturation) .or. &
      FlowSubConditionIsTransient(condition%flux) .or. &
      FlowSubConditionIsTransient(condition%rate) .or. &
      FlowSubConditionIsTransient(condition%temperature) .or. &
      FlowSubConditionIsTransient(condition%enthalpy) .or. &
      FlowSubConditionIsTransient(condition%energy_rate) .or. &
      FlowSubConditionIsTransient(condition%energy_flux) .or. &
      FlowSubConditionIsTransient(condition%energy_conductance) &
      ) then
    FlowConditionIsTransient = PETSC_TRUE
  endif

end function FlowConditionIsTransient

! ************************************************************************** !

function FlowSubConditionIsTransient(sub_condition)
  !
  ! Returns PETSC_TRUE
  !
  ! Author: Glenn Hammond
  ! Date: 10/26/11
  !

  use Dataset_module

  implicit none

  type(flow_sub_condition_type), pointer :: sub_condition

  PetscBool :: FlowSubConditionIsTransient

  FlowSubConditionIsTransient = PETSC_FALSE

  if (associated(sub_condition)) then
    if (DatasetIsTransient(sub_condition%dataset) .or. &
        DatasetIsTransient(sub_condition%gradient)) then
      FlowSubConditionIsTransient = PETSC_TRUE
    endif
  endif

end function FlowSubConditionIsTransient

! ************************************************************************** !

function FlowConditionUnknownItype(condition,message,type_name)
  !
  ! Returns a string indicating which flow condition has a wrong type.
  !
  ! Author: Glenn Hammond
  ! Date: 03/21/16
  !
  implicit none

  type(flow_condition_type) :: condition
  character(len=*) :: message
  character(len=*) :: type_name

  character(len=MAXSTRINGLENGTH) :: FlowConditionUnknownItype

  FlowConditionUnknownItype = 'Unknown TYPE (' // trim(type_name) // &
    ') for ' // trim(message) // ' within FLOW_CONDITION "' // &
    trim(condition%name) // '".'

end function FlowConditionUnknownItype

! **************************************************************************** !

subroutine FlowConditionDestroyList(condition_list)
  !
  ! Deallocates a list of conditions
  !
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  !

  implicit none

  type(condition_list_type), pointer :: condition_list

  type(flow_condition_type), pointer :: condition, prev_condition

  if (.not.associated(condition_list)) return

  condition => condition_list%first
  do
    if (.not.associated(condition)) exit
    prev_condition => condition
    condition => condition%next
    call FlowConditionDestroy(prev_condition)
  enddo

  condition_list%num_conditions = 0
  nullify(condition_list%first)
  nullify(condition_list%last)
  if (associated(condition_list%array)) deallocate(condition_list%array)
  nullify(condition_list%array)

  deallocate(condition_list)
  nullify(condition_list)

end subroutine FlowConditionDestroyList

! ************************************************************************** !

subroutine FlowConditionDestroy(condition)
  !
  ! Deallocates a condition
  !
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  !

  use Dataset_module
  use Dataset_Ascii_class
  use Utility_module

  implicit none

  type(flow_condition_type), pointer :: condition

  class(dataset_ascii_type), pointer :: dataset_ascii
  PetscInt :: i

  if (.not.associated(condition)) return

  ! if dataset_ascii_type, destroy.  Otherwise, they are in another list
  dataset_ascii => DatasetAsciiCast(condition%datum)
  ! dataset_ascii will be NULL if not dataset_ascii_type
  call DatasetAsciiDestroy(dataset_ascii)

  if (associated(condition%sub_condition_ptr)) then
    ! only nullify the pointers; don't destroy as they are destroyed below
    do i=1,condition%num_sub_conditions
      nullify(condition%sub_condition_ptr(i)%ptr)
    enddo
    deallocate(condition%sub_condition_ptr)
    nullify(condition%sub_condition_ptr)
  endif

  call DeallocateArray(condition%itype)

  call FlowSubConditionDestroy(condition%pressure)
  call FlowSubConditionDestroy(condition%saturation)
  call FlowSubConditionDestroy(condition%rate)
  call FlowSubConditionDestroy(condition%flux)
  call FlowSubConditionDestroy(condition%conductance)
  call FlowSubConditionDestroy(condition%temperature)
  call FlowSubConditionDestroy(condition%enthalpy)
  call FlowSubConditionDestroy(condition%energy_rate)
  call FlowSubConditionDestroy(condition%energy_flux)
  call FlowSubConditionDestroy(condition%energy_conductance)

  call TimeStorageDestroy(condition%default_time_storage)

  nullify(condition%next)

  deallocate(condition)
  nullify(condition)

end subroutine FlowConditionDestroy


! ************************************************************************** !

subroutine FlowSubConditionDestroy(sub_condition)
  !
  ! Destroys a sub_condition
  !
  ! Author: Glenn Hammond
  ! Date: 02/04/08
  !

  use Dataset_module
  use Dataset_Ascii_class

  implicit none

  type(flow_sub_condition_type), pointer :: sub_condition

  class(dataset_ascii_type), pointer :: dataset_ascii

  if (.not.associated(sub_condition)) return

  ! if dataset_ascii_type, destroy.  Otherwise, they are in another list
  dataset_ascii => DatasetAsciiCast(sub_condition%dataset)
  ! dataset_ascii will be NULL if not dataset_ascii_type
  call DatasetAsciiDestroy(dataset_ascii)
  dataset_ascii => DatasetAsciiCast(sub_condition%gradient)
  call DatasetAsciiDestroy(dataset_ascii)

  deallocate(sub_condition)
  nullify(sub_condition)

end subroutine FlowSubConditionDestroy

! ************************************************************************** !

end module FlowCondition_module
