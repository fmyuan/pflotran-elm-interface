module Inversion_Measurement_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: inversion_measurement_aux_type
    PetscInt :: id
    PetscReal :: time
    PetscInt :: cell_id
    PetscReal :: value
    type(inversion_measurement_aux_type), pointer :: next
  end type inversion_measurement_aux_type

  public :: InversionMeasurementAuxCreate, &
            InversionMeasurementAuxInit, &
            InversionMeasurementAuxCopy, &
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
  measurement%cell_id = UNINITIALIZED_INTEGER
  measurement%value = UNINITIALIZED_DOUBLE

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
  measurement2%cell_id = measurement%cell_id
  measurement2%value = measurement%value

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

  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option

  type(inversion_measurement_aux_type), pointer :: InversionMeasurementAuxRead

  character(len=MAXWORDLENGTH) :: keyword
  type(inversion_measurement_aux_type), pointer :: new_measurement

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
      case('CELL_ID')
        call InputReadInt(input,option,new_measurement%cell_id)
        call InputErrorMsg(input,option,keyword,error_string)
      case('VALUE')
        call InputReadDouble(input,option,new_measurement%value)
        call InputErrorMsg(input,option,keyword,error_string)
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select

  enddo
  call InputPopBlock(input,option)

  if (UnInitialized(new_measurement%cell_id)) then
    option%io_buffer = 'CELL_ID not specified for measurement.'
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
