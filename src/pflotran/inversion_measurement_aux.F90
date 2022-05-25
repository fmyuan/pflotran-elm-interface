module Inversion_Measurement_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Geometry_module

  implicit none

  private

  type, public :: inversion_measurement_aux_type
    PetscInt :: id
    PetscReal :: time
    character(len=4) :: time_units
    PetscInt :: cell_id
    PetscReal :: value
    PetscReal :: simulated_value
    PetscBool :: first_lambda
    PetscBool :: measured
    type(point3d_type) :: coordinate
    type(inversion_measurement_aux_type), pointer :: next
  end type inversion_measurement_aux_type

  public :: InversionMeasurementAuxCreate, &
            InversionMeasurementAuxInit, &
            InversionMeasurementAuxCopy, &
            InversionMeasurementPrint, &
            InversionMeasurementMeasure, &
            InversionMeasurementAuxRead, &
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
  measurement%value = UNINITIALIZED_DOUBLE
  measurement%simulated_value = UNINITIALIZED_DOUBLE
  measurement%first_lambda = PETSC_FALSE
  measurement%measured = PETSC_FALSE

  call GeometryInitCoordinate(measurement%coordinate)

  nullify(measurement%next)

end subroutine InversionMeasurementAuxInit

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
  measurement2%value = measurement%value
  measurement2%simulated_value = measurement%simulated_value
  call GeometryCopyCoordinate(measurement%coordinate,measurement2%coordinate)

end subroutine InversionMeasurementAuxCopy

! ************************************************************************** !

subroutine InversionMeasurementMeasure(time,measurement,value_)
  !
  ! Copies the value into measurement
  !
  ! Author: Glenn Hammond
  ! Date: 02/21/22
  !
  PetscReal :: time
  type(inversion_measurement_aux_type) :: measurement
  PetscReal :: value_

  PetscReal, parameter :: tol = 1.d0
  PetscBool :: measure

  measure = PETSC_FALSE
  if (Uninitialized(measurement%time)) then
    measure = PETSC_TRUE
  else
    measure = .not.measurement%measured .and. measurement%time <= time + tol
  endif

  if (measure) then
    measurement%simulated_value = value_
    measurement%measured = PETSC_TRUE
  endif

end subroutine InversionMeasurementMeasure

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

  if (optionPrintToScreen(option)) then
    print *, 'Measurment: ' // trim(StringWrite(measurement%id))
    word = 'sec'
    print *, '             Time: ' // &
      trim(StringWrite(measurement%time / &
                       UnitsConvertToInternal(measurement%time_units,word, &
                                              option,ierr))) // ' ' // &
      measurement%time_units
    print *, '          Cell ID: ' // trim(StringWrite(measurement%cell_id))
    print *, '            Value: ' // trim(StringWrite(measurement%value))
    print *, '  Simulated Value: ' // &
      trim(StringWrite(measurement%simulated_value))
  endif

end subroutine InversionMeasurementPrint

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

  InversionMeasurementAuxRead => new_measurement

end function InversionMeasurementAuxRead

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

  PetscInt :: iaux
  PetscErrorCode :: ierr

  if (.not.associated(aux)) return

  call InversionMeasurementAuxStrip(aux)

  deallocate(aux)
  nullify(aux)

end subroutine InversionMeasurementAuxDestroy

end module Inversion_Measurement_Aux_module
