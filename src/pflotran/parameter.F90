module Parameter_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: parameter_type
    character(len=MAXWORDLENGTH) :: name
    PetscInt :: id
    PetscReal :: value
    character(len=MAXWORDLENGTH) :: dataset_name
    type(parameter_type), pointer :: next
  end type parameter_type

  public :: ParameterCreate, ParameterDestroy, &
            ParameterRead, ParameterAddToList, &
            ParameterSetup, ParameterGetIDFromName, &
            ParameterQualityCheck

contains

! ************************************************************************** !

function ParameterCreate()
  !
  ! Creates a parameter object
  !
  ! Author: Glenn Hammond
  ! Date: 01/04/23

  type(parameter_type), pointer :: ParameterCreate

  type(parameter_type), pointer :: parameter

  allocate(parameter)
  parameter%name = ''
  parameter%id = UNINITIALIZED_INTEGER
  parameter%value = UNINITIALIZED_DOUBLE
  parameter%dataset_name = ''
  nullify(parameter%next)
  ParameterCreate => parameter

end function ParameterCreate

! ************************************************************************** !

subroutine ParameterRead(parameter,input,option)
  !
  ! Reads in contents of a parameter card
  !
  ! Author: Glenn Hammond
  ! Date: 01/04/23

  use Option_module
  use Input_Aux_module
  use String_module
  use Units_module

  type(parameter_type) :: parameter
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_str

  error_str = 'PARAMETER'

  input%ierr = 0
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_str)
    call StringToUpper(keyword)

    select case(trim(keyword))

      case('NAME')
        call InputReadWord(input,option,parameter%name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_str)
      case('VALUE')
        call InputReadDouble(input,option,parameter%value)
        call InputErrorMsg(input,option,keyword,error_str)
      case('DATASET')
        call InputReadWord(input,option,parameter%dataset_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_str)
      case default
        call InputKeywordUnrecognized(input,keyword,error_str,option)
    end select

  enddo
  call InputPopBlock(input,option)

  if (len_trim(parameter%name) == 0) then
    option%io_buffer = 'No NAME specified for a PARAMETER in the input file.'
    call PrintErrMsg(option)
  endif

  if (Uninitialized(parameter%value) .and. &
     len_trim(parameter%dataset_name) == 0) then
    option%io_buffer = 'No VALUE or DATASET specified for a PARAMETER "' // &
      trim(parameter%name) // '" in the input file.'
    call PrintErrMsg(option)
  endif

end subroutine ParameterRead

! ************************************************************************** !

subroutine ParameterAddToList(parameter,list)
  !
  ! Adds a parameter to linked list
  !
  ! Author: Glenn Hammond
  ! Date: 01/04/23

  type(parameter_type), pointer :: parameter
  type(parameter_type), pointer :: list

  type(parameter_type), pointer :: cur_parameter

  if (associated(list)) then
    cur_parameter => list
    ! loop to end of list
    do
      if (.not.associated(cur_parameter%next)) exit
      cur_parameter => cur_parameter%next
    enddo
    cur_parameter%next => parameter
  else
    list => parameter
  endif

end subroutine ParameterAddToList

! ************************************************************************** !

subroutine ParameterSetup(parameter_list,option)
  !
  ! Initializes the parameter module
  !
  ! Author: Glenn Hammond
  ! Date: 01/0523

  use Option_module
  use Option_Parameter_module

  type(parameter_type), pointer :: parameter_list
  type(option_type) :: option

  type(parameter_option_type), pointer :: parameter_option
  type(parameter_type), pointer :: cur_parameter
  PetscInt :: iparameter

  parameter_option => option%parameter

  iparameter = 0
  cur_parameter => parameter_list
  do
    if (.not.associated(cur_parameter)) exit
    iparameter = iparameter + 1
    cur_parameter%id = iparameter
    cur_parameter => cur_parameter%next
  enddo
  parameter_option%num_parameters = iparameter
  allocate(parameter_option%parameter_names(iparameter))
  cur_parameter => parameter_list
  do
    if (.not.associated(cur_parameter)) exit
    parameter_option%parameter_names(cur_parameter%id) = cur_parameter%name
    cur_parameter => cur_parameter%next
  enddo

end subroutine ParameterSetup

! ************************************************************************** !

subroutine ParameterQualityCheck(parameter_list,option)
  !
  ! Initializes the parameter module
  !
  ! Author: Glenn Hammond
  ! Date: 01/0523

  use Option_module
  use String_module

  type(parameter_type), pointer :: parameter_list
  type(option_type) :: option

  type(parameter_type), pointer :: cur_parameter
  PetscInt :: iparameter

  iparameter = 0
  cur_parameter => parameter_list
  do
    if (.not.associated(cur_parameter)) exit
    iparameter = iparameter + 1
    cur_parameter => cur_parameter%next
  enddo
  if (iparameter /= option%parameter%num_parameters) then
    option%io_buffer = &
      StringWrite(iparameter-option%parameter%num_parameters) // &
      ' parameters were added to the parameter list after it was "setup". &
      &Please send your input deck to pflotran-dev@googlegroups.com'
    call PrintErrMsg(option)
  endif

end subroutine ParameterQualityCheck

! ************************************************************************** !

function ParameterGetIDFromName(parameter_name,option)
  !
  ! Returns the parameter id/index given a name
  !
  ! Author: Glenn Hammond
  ! Date: 01/0523

  use Option_module
  use String_module

  character(len=*) :: parameter_name
  type(option_type) :: option

  PetscInt :: ParameterGetIDFromName

  PetscInt :: i
  character(len=MAXWORDLENGTH), pointer :: parameter_names(:)

  parameter_names => option%parameter%parameter_names
  ParameterGetIDFromName = 0
  do i = 1, option%parameter%num_parameters
    if (StringCompare(parameter_name,parameter_names(i))) then
      ParameterGetIDFromName = i
      return
    endif
  enddo
  if (option%parameter%num_parameters > 0) then
    option%io_buffer = 'Parameter "' // trim(parameter_name) // '" was not &
      &found among available parameters: ' // StringsMerge(parameter_names,',')
  else
    option%io_buffer = 'Parameter "' // trim(parameter_name) // '" was &
      &requested in ParameterGetIDFromName(), but no PARAMETERs were &
      &included in the input deck.'
  endif
  call PrintErrMsg(option)

end function ParameterGetIDFromName

! ************************************************************************** !

recursive subroutine ParameterDestroy(parameter)
  !
  ! Destroys a parameter
  !
  ! Author: Glenn Hammond
  ! Date: 01/04/23

  type(parameter_type), pointer :: parameter

  if (.not.associated(parameter)) return

  call ParameterDestroy(parameter%next)

  deallocate(parameter)
  nullify(parameter)

end subroutine ParameterDestroy

end module Parameter_module
