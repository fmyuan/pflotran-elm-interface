module Inversion_Measurement_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Geometry_module
  use Option_Inversion_module

  implicit none

  private

  PetscInt, public :: inv_meas_reporting_verbosity

  character(len=MAXWORDLENGTH), parameter :: OBS_LIQUID_PRESSURE_STRING = &
                                               'LIQUID_PRESSURE'
  character(len=MAXWORDLENGTH), parameter :: OBS_LIQUID_SATURATION_STRING = &
                                               'LIQUID_SATURATION'
  character(len=MAXWORDLENGTH), parameter :: OBS_SOLUTE_CONCENTRATION_STRING = &
                                               'SOLUTE_CONCENTRATION'
  character(len=MAXWORDLENGTH), parameter :: OBS_ERT_MEASUREMENT_STRING = &
                                               'ERT_MEASUREMENT'

  PetscInt, parameter, public :: OBS_LIQUID_PRESSURE = 1
  PetscInt, parameter, public :: OBS_LIQUID_SATURATION = 2
  PetscInt, parameter, public :: OBS_SOLUTE_CONCENTRATION = 3
  PetscInt, parameter, public :: OBS_ERT_MEASUREMENT = 4

  type, public :: inversion_measurement_aux_type
    PetscInt :: id
    PetscReal :: time
    character(len=4) :: time_units
    PetscInt :: cell_id
    PetscInt :: local_id
    PetscInt :: iobs_var
    PetscReal :: value
    PetscReal :: weight
    PetscReal :: simulated_derivative
    PetscReal :: simulated_value
    PetscBool :: first_lambda
    PetscBool :: measured
    type(point3d_type) :: coordinate
    type(inversion_measurement_aux_type), pointer :: next
  end type inversion_measurement_aux_type

  public :: InversionMeasurementAuxCreate, &
            InversionMeasurementAuxInit, &
            InversionMeasurementAuxReset, &
            InversionMeasurementAuxCopy, &
            InversionMeasurementPrint, &
            InversionMeasurementPrintConcise, &
            InvMeasAnnounceToString, &
            InvMeasAuxObsVarIDToString, &
            InversionMeasurementMeasure, &
            InversionMeasurementAuxRead, &
            InvMeasAuxReadObservedVariable, &
            InversionMeasureAuxListDestroy, &
            InversionMeasurementAuxStrip, &
            InversionMeasurementAuxDestroy

contains

! ************************************************************************** !

function InversionMeasurementAuxCreate()
  !
  ! Allocate and initialize auxiliary inversion measurement object
  !
  ! Author: Glenn Hammond
  ! Date: 02/07/22
  !

  type(inversion_measurement_aux_type), pointer :: &
    InversionMeasurementAuxCreate

  type(inversion_measurement_aux_type), pointer :: aux

  allocate(aux)
  call InversionMeasurementAuxInit(aux)

  InversionMeasurementAuxCreate => aux

end function InversionMeasurementAuxCreate

! ************************************************************************** !

subroutine InversionMeasurementAuxInit(measurement)
  !
  ! Initializes auxiliary inversion measurement object
  !
  ! Author: Glenn Hammond
  ! Date: 02/07/22
  !

  type(inversion_measurement_aux_type) :: measurement

  measurement%id = UNINITIALIZED_INTEGER
  measurement%time = UNINITIALIZED_DOUBLE
  measurement%time_units = ''
  measurement%cell_id = UNINITIALIZED_INTEGER
  measurement%local_id = UNINITIALIZED_INTEGER
  measurement%iobs_var = UNINITIALIZED_INTEGER
  measurement%value = UNINITIALIZED_DOUBLE
  measurement%weight = UNINITIALIZED_DOUBLE
  measurement%simulated_derivative = UNINITIALIZED_DOUBLE
  measurement%simulated_value = UNINITIALIZED_DOUBLE
  measurement%first_lambda = PETSC_FALSE
  measurement%measured = PETSC_FALSE

  call GeometryInitCoordinate(measurement%coordinate)

  nullify(measurement%next)

end subroutine InversionMeasurementAuxInit

! ************************************************************************** !

subroutine InversionMeasurementAuxReset(measurement)
  !
  ! Resets measurement data at beginning of inversion iteration
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22
  !
  type(inversion_measurement_aux_type) :: measurement

  measurement%simulated_value = UNINITIALIZED_DOUBLE
  measurement%first_lambda = PETSC_FALSE
  measurement%measured = PETSC_FALSE

end subroutine InversionMeasurementAuxReset

! ************************************************************************** !

subroutine InversionMeasurementAuxCopy(measurement,measurement2)
  !
  ! Copies auxiliary inversion measurement object
  !
  ! Author: Glenn Hammond
  ! Date: 02/07/22
  !

  type(inversion_measurement_aux_type) :: measurement
  type(inversion_measurement_aux_type) :: measurement2

  measurement2%id = measurement%id
  measurement2%time = measurement%time
  measurement2%time_units = measurement%time_units
  measurement2%cell_id = measurement%cell_id
  measurement2%local_id = measurement%local_id
  measurement2%iobs_var = measurement%iobs_var
  measurement2%value = measurement%value
  measurement2%weight = measurement%weight
  measurement2%simulated_derivative = measurement%simulated_derivative
  measurement2%simulated_value = measurement%simulated_value
  call GeometryCopyCoordinate(measurement%coordinate,measurement2%coordinate)

end subroutine InversionMeasurementAuxCopy

! ************************************************************************** !

function InversionMeasurementAuxRead(input,error_string,option)
  !
  ! Reads measurements and appends to the list
  !
  ! Author: Glenn Hammond
  ! Date: 02/07/22
  !
  use Input_Aux_module
  use Option_module
  use String_module
  use Units_module

  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option

  type(inversion_measurement_aux_type), pointer :: InversionMeasurementAuxRead

  character(len=MAXWORDLENGTH) :: keyword
  type(inversion_measurement_aux_type), pointer :: new_measurement
  PetscReal :: units_conversion
  PetscReal :: sd
  character(len=MAXWORDLENGTH) :: internal_units, word
  character(len=MAXSTRINGLENGTH) :: string

  new_measurement => InversionMeasurementAuxCreate()

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
      case('TIME')
        call InputReadDouble(input,option,new_measurement%time)
        call InputErrorMsg(input,option,keyword,error_string)
        call InputReadWord(input,option,word,PETSC_TRUE)
        if (input%ierr /= 0) word = 'sec'
        new_measurement%time_units = trim(word)
        internal_units = 'sec'
        units_conversion = UnitsConvertToInternal(word,internal_units,option)
        new_measurement%time = new_measurement%time*units_conversion
      case('CELL_ID','ERT_MEASUREMENT_ID')
        call InputReadInt(input,option,new_measurement%cell_id)
        call InputErrorMsg(input,option,keyword,error_string)
      case('COORDINATE')
        string = trim(error_string) // ',' // keyword
        call GeometryReadCoordinate(input,option,new_measurement%coordinate, &
                                    string)
      case('VALUE')
        call InputReadDouble(input,option,new_measurement%value)
        call InputErrorMsg(input,option,keyword,error_string)
      case('STANDARD_DEVIATION')
        call InputReadDouble(input,option,sd)
        call InputErrorMsg(input,option,keyword,error_string)
        if (sd <= 0) sd = 1.d16
        new_measurement%weight = 1 / sd
      case('OBSERVED_VARIABLE')
        new_measurement%iobs_var = &
          InvMeasAuxReadObservedVariable(input,keyword,error_string,option)
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select

  enddo
  call InputPopBlock(input,option)

  if (UnInitialized(new_measurement%cell_id) .and. &
      UnInitialized(new_measurement%coordinate%x)) then
    option%io_buffer = 'CELL_ID or COORDINATE not specified for measurement.'
    call PrintErrMsg(option)
  endif
  if (UnInitialized(new_measurement%value)) then
    option%io_buffer = 'VALUE not specified for measurement.'
    call PrintErrMsg(option)
  endif
  if (UnInitialized(new_measurement%weight)) then
    sd = 0.05 * new_measurement%value
    new_measurement%weight = 1 / sd
  endif

  InversionMeasurementAuxRead => new_measurement

end function InversionMeasurementAuxRead

! ************************************************************************** !

function InvMeasAuxReadObservedVariable(input,keyword,error_string,option)
  !
  ! Reads the observed variable and returns a corresponding integer ID
  !
  ! Author: Glenn Hammond
  ! Date: 06/20/22
  !
  use Input_Aux_module
  use Option_module
  use String_module

  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option

  PetscInt :: InvMeasAuxReadObservedVariable

  character(len=MAXWORDLENGTH) :: word

  call InputReadWord(input,option,word,PETSC_TRUE)
  call InputErrorMsg(input,option,keyword,error_string)
  call StringToUpper(word)
  InvMeasAuxReadObservedVariable = InvMeasAuxObsVarStringToInt(word)
  if (Uninitialized(InvMeasAuxReadObservedVariable)) then
    call InputKeywordUnrecognized(input,word,trim(error_string)// &
                                  & ','//trim(word),option)
  endif

end function InvMeasAuxReadObservedVariable

! ************************************************************************** !

function InvMeasAuxObsVarStringToInt(string)
  !
  ! Maps an observation variable string to its integer ID
  !
  ! Author: Glenn Hammond
  ! Date: 07/15/22
  !
  character(len=*) :: string

  PetscInt :: InvMeasAuxObsVarStringToInt

  InvMeasAuxObsVarStringToInt = UNINITIALIZED_INTEGER
  select case(string)
    case(OBS_LIQUID_PRESSURE_STRING)
      InvMeasAuxObsVarStringToInt = OBS_LIQUID_PRESSURE
    case(OBS_LIQUID_SATURATION_STRING)
      InvMeasAuxObsVarStringToInt = OBS_LIQUID_SATURATION
    case(OBS_SOLUTE_CONCENTRATION_STRING)
      InvMeasAuxObsVarStringToInt = OBS_SOLUTE_CONCENTRATION
    case(OBS_ERT_MEASUREMENT_STRING)
      InvMeasAuxObsVarStringToInt = OBS_ERT_MEASUREMENT
  end select

end function InvMeasAuxObsVarStringToInt

! ************************************************************************** !

function InvMeasAuxObsVarIDToString(id,option)
  !
  ! Maps an observation variable integer ID to string
  !
  ! Author: Glenn Hammond
  ! Date: 07/15/22
  !
  use Option_module
  use String_module

  PetscInt :: id
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: InvMeasAuxObsVarIDToString

  InvMeasAuxObsVarIDToString = ''
  select case(id)
    case(OBS_LIQUID_PRESSURE)
      InvMeasAuxObsVarIDToString = OBS_LIQUID_PRESSURE_STRING
    case(OBS_LIQUID_SATURATION)
      InvMeasAuxObsVarIDToString = OBS_LIQUID_SATURATION_STRING
    case(OBS_SOLUTE_CONCENTRATION)
      InvMeasAuxObsVarIDToString = OBS_SOLUTE_CONCENTRATION_STRING
    case(OBS_ERT_MEASUREMENT)
      InvMeasAuxObsVarIDToString = OBS_ERT_MEASUREMENT_STRING
    case default
      option%io_buffer = 'Unknown measurement variable integer ID in &
        &InvMeasAuxObsVarIDToString: ' // trim(StringWrite(id))
      call PrintErrMsg(option)
  end select

end function InvMeasAuxObsVarIDToString

! ************************************************************************** !

subroutine InversionMeasurementMeasure(time,measurement,value_,option)
  !
  ! Copies the value into measurement
  !
  ! Author: Glenn Hammond
  ! Date: 02/21/22

  use Option_module
  use String_module
  use Utility_module

  PetscReal :: time
  type(inversion_measurement_aux_type) :: measurement
  PetscReal :: value_
  type(option_type) :: option

  PetscBool :: measure

  measure = PETSC_FALSE
!  if (Uninitialized(measurement%time)) then
!    measure = PETSC_TRUE
!  else
    measure = .not.measurement%measured .and. Equal(measurement%time,time)
!  endif

  if (measure) then
    measurement%simulated_value = value_
    measurement%measured = PETSC_TRUE
    if (inv_meas_reporting_verbosity > 0) then
      option%io_buffer = '  Recording measurement #' // &
        trim(StringWrite(measurement%id))
      if (inv_meas_reporting_verbosity > 1) then
        option%io_buffer = trim(option%io_buffer) // &
          ' sim value = ' // &
          trim(StringWrite('(es13.6)',measurement%simulated_value))
      endif
      if (inv_meas_reporting_verbosity > 2) then
        option%io_buffer = trim(option%io_buffer) // &
          ', orig value = ' // &
          trim(StringWrite('(es13.6)',measurement%value))
      endif
      call PrintMsg(option)
    endif
  endif

end subroutine InversionMeasurementMeasure

! ************************************************************************** !

subroutine InversionMeasurementPrintConcise(measurement,optional_string, &
                                            option)
  !
  ! Print contents of measurement object for debugging
  !
  ! Author: Glenn Hammond
  ! Date: 02/14/22
  !
  use Option_module
  use String_module
  use Units_module

  type(inversion_measurement_aux_type) :: measurement
  character(len=*) :: optional_string
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  PetscErrorCode :: ierr

  if (OptionPrintToScreen(option)) then
    word = 'sec'
    option%io_buffer = 'Measurement #' // &
      trim(StringWrite(measurement%id)) // ', Time: ' // &
      trim(StringWrite(measurement%time / &
                       UnitsConvertToInternal(measurement%time_units,word, &
                                              option,ierr))) // ' ' // &
      trim(measurement%time_units) // ', Var: ' // &
      trim(InvMeasAuxObsVarIDToString(measurement%iobs_var,option)) // &
      ', Cell: ' // &
      trim(StringWrite(measurement%cell_id)) // ', Value: ' // &
      trim(StringWrite(measurement%simulated_value)) // ', Deriv: ' // &
      trim(StringWrite(measurement%simulated_derivative))
    if (len_trim(optional_string) > 0) then
      option%io_buffer = trim(optional_string) // ' : ' // &
        trim(option%io_buffer)
    endif
    call PrintMsg(option)
  endif

end subroutine InversionMeasurementPrintConcise

! ************************************************************************** !

subroutine InversionMeasurementPrint(measurement,option)
  !
  ! Print contents of measurement object for debugging
  !
  ! Author: Glenn Hammond
  ! Date: 02/14/22
  !
  use Option_module
  use String_module
  use Units_module

  type(inversion_measurement_aux_type) :: measurement
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  PetscErrorCode :: ierr

  if (OptionPrintToScreen(option)) then
    print *, 'Measurement: ' // trim(StringWrite(measurement%id))
    word = 'sec'
    print *, '                Time: ' // &
      trim(StringWrite(measurement%time / &
                       UnitsConvertToInternal(measurement%time_units,word, &
                                              option,ierr))) // ' ' // &
      measurement%time_units
    print *, '             Cell ID: ' // trim(StringWrite(measurement%cell_id))
    print *, '            Variable: ' // &
      trim(InvMeasAuxObsVarIDToString(measurement%iobs_var,option))
    print *, '               Value: ' // trim(StringWrite(measurement%value))
    print *, '     Simulated Value: ' // &
      trim(StringWrite(measurement%simulated_value))
    print *, 'Simulated Derivative: ' // &
      trim(StringWrite(measurement%simulated_derivative))
  endif

end subroutine InversionMeasurementPrint

! ************************************************************************** !

function InvMeasAnnounceToString(measurement,rvalue,option)
  !
  ! Announces the recording of a measurement for inversion
  !
  ! Author: Glenn Hammond
  ! Date: 07/15/22
  !
  use Option_module
  use String_module
  use Units_module

  type(inversion_measurement_aux_type) :: measurement
  PetscReal :: rvalue
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: InvMeasAnnounceToString

  character(len=MAXWORDLENGTH), parameter :: word = 'sec'
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr

  string = 'Measurement #' // trim(StringWrite(measurement%id)) // &
    ' at Cell ID ' // trim(StringWrite(measurement%cell_id))
  if (Initialized(measurement%coordinate%x)) then
    string = trim(string) // ' (' // &
      trim(adjustl(StringWriteF('(f20.2)',measurement%coordinate%x))) // ',' // &
      trim(adjustl(StringWriteF('(f20.2)',measurement%coordinate%y))) // ',' // &
      trim(adjustl(StringWriteF('(f20.2)',measurement%coordinate%z))) // ')'
  endif
  string = trim(string) // ' for variable "' // &
           trim(InvMeasAuxObsVarIDToString(measurement%iobs_var,option))
  if (Initialized(rvalue)) then
    string = trim(string) // '" recorded as ' // &
      trim(StringWrite('(es22.14)',rvalue))
  endif
  string = trim(string) // ' at ' // &
           trim(StringWrite(measurement%time / &
                      UnitsConvertToInternal(measurement%time_units,word, &
                                             option,ierr))) // ' ' // &
           measurement%time_units
  InvMeasAnnounceToString = string

end function InvMeasAnnounceToString

! ************************************************************************** !

subroutine InversionMeasureAuxListDestroy(inv_measure_aux_list)
  !
  ! Deallocates a inversion auxiliary measurement list object
  !
  ! Author: Glenn Hammond
  ! Date: 02/07/22
  !
  type(inversion_measurement_aux_type), pointer :: inv_measure_aux_list

  type(inversion_measurement_aux_type), pointer :: cur_inv_measure_aux
  type(inversion_measurement_aux_type), pointer :: next_inv_measure_aux

  cur_inv_measure_aux => inv_measure_aux_list
  do
    if (.not.associated(cur_inv_measure_aux)) exit
    next_inv_measure_aux => cur_inv_measure_aux%next
    call InversionMeasurementAuxDestroy(cur_inv_measure_aux)
    cur_inv_measure_aux => next_inv_measure_aux
  enddo

  nullify(inv_measure_aux_list)

end subroutine InversionMeasureAuxListDestroy

! ************************************************************************** !

subroutine InversionMeasurementAuxStrip(aux)
  !
  ! Deallocates members of measurement auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 02/07/22
  !

  type(inversion_measurement_aux_type) :: aux

  nullify(aux%next)

end subroutine InversionMeasurementAuxStrip

! ************************************************************************** !

subroutine InversionMeasurementAuxDestroy(aux)
  !
  ! Deallocates a inversion measurement auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 02/07/22
  !
  use Utility_module, only : DeallocateArray

  type(inversion_measurement_aux_type), pointer :: aux

  if (.not.associated(aux)) return

  call InversionMeasurementAuxStrip(aux)

  deallocate(aux)
  nullify(aux)

end subroutine InversionMeasurementAuxDestroy

end module Inversion_Measurement_Aux_module
