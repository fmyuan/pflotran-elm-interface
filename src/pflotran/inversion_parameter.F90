module Inversion_Parameter_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: inversion_parameter_type
    PetscInt :: id
    PetscInt :: iparameter
    PetscInt :: imat
    PetscReal :: value
    character(len=MAXWORDLENGTH) :: parameter_name
    character(len=MAXWORDLENGTH) :: material_name
    type(inversion_parameter_type), pointer :: next
  end type inversion_parameter_type

  public :: InversionParameterCreate, &
            InversionParameterInit, &
            InversionParameterRead, &
            InversionParameterCopy, &
            InversionParameterMapNameToInt, &
            InversionParameterIntToQOIArray, &
            InversionParameterPrint, &
            InversionParameterStrip, &
            InversionParameterDestroy

contains

! ************************************************************************** !

function InversionParameterCreate()
  !
  ! Allocate and initialize inversion parameter object
  !
  ! Author: Glenn Hammond
  ! Date: 03/18/22
  !

  type(inversion_parameter_type), pointer :: InversionParameterCreate

  type(inversion_parameter_type), pointer :: inversion_parameter

  allocate(inversion_parameter)
  call InversionParameterInit(inversion_parameter)

  InversionParameterCreate => inversion_parameter

end function InversionParameterCreate

! ************************************************************************** !

subroutine InversionParameterInit(inversion_parameter)
  !
  ! Initializes auxiliary inversion parameter object
  !
  ! Author: Glenn Hammond
  ! Date: 03/18/22
  !

  type(inversion_parameter_type) :: inversion_parameter

  inversion_parameter%id = UNINITIALIZED_INTEGER
  inversion_parameter%iparameter = UNINITIALIZED_INTEGER
  inversion_parameter%imat = UNINITIALIZED_INTEGER
  inversion_parameter%value = UNINITIALIZED_DOUBLE
  inversion_parameter%parameter_name = ''
  inversion_parameter%material_name = ''

  nullify(inversion_parameter%next)

end subroutine InversionParameterInit

! ************************************************************************** !

subroutine InversionParameterCopy(inversion_parameter,inversion_parameter2)
  !
  ! Copies auxiliary inversion parameter object
  !
  ! Author: Glenn Hammond
  ! Date: 03/18/22
  !

  type(inversion_parameter_type) :: inversion_parameter
  type(inversion_parameter_type) :: inversion_parameter2

  inversion_parameter2%id = inversion_parameter%id
  inversion_parameter2%iparameter = inversion_parameter%iparameter
  inversion_parameter2%imat = inversion_parameter%imat
  inversion_parameter2%value = inversion_parameter%value
  inversion_parameter2%parameter_name = inversion_parameter%parameter_name
  inversion_parameter2%material_name = inversion_parameter%material_name

end subroutine InversionParameterCopy

! ************************************************************************** !

subroutine InversionParameterPrint(fid,inversion_parameter, &
                                   print_header,print_footer,option)
  !
  ! Print contents of inversion parameter object
  !
  ! Author: Glenn Hammond
  ! Date: 03/18/22
  !
  use Option_module
  use String_module
  use Units_module

  PetscInt :: fid
  type(inversion_parameter_type) :: inversion_parameter
  PetscBool :: print_header
  PetscBool :: print_footer
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (print_header) then
!                12345678901234567890123456789012345678901234567890
    write(fid,'(40("=+"),//, &
              &" Current values of inversion parameters:",//, &
              &"    # &
              &Parameter Name      &
              &Material Name       &
              & Value",/,&
              &"    - &
              &--------------      &
              &-------------       &
              & -----")')
  endif
  write(string,'(i4," ",2a20,es13.6)') &
    inversion_parameter%id, &
    inversion_parameter%parameter_name, &
    inversion_parameter%material_name, &
    inversion_parameter%value
  write(fid,*) trim(string)
  if (print_footer) then
    write(fid,'(/,40("=+"))')
  endif

end subroutine InversionParameterPrint

! ************************************************************************** !

function InversionParameterRead(input,error_string,option)
  !
  ! Reads measurements and appends to the list
  !
  ! Author: Glenn Hammond
  ! Date: 03/18/22
  !
  use Input_Aux_module
  use Option_module
  use String_module
  use Units_module

  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option

  type(inversion_parameter_type), pointer :: InversionParameterRead

  character(len=MAXWORDLENGTH) :: keyword
  type(inversion_parameter_type), pointer :: new_inversion_parameter
  character(len=MAXWORDLENGTH) :: word

  new_inversion_parameter => InversionParameterCreate()

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
      case('NAME')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_string)
        call StringtoUpper(word)
        new_inversion_parameter%parameter_name = word
      case('MATERIAL')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_string)
        new_inversion_parameter%material_name = word
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select

  enddo
  call InputPopBlock(input,option)

  if (len_trim(new_inversion_parameter%parameter_name) == 0) then
    option%io_buffer = 'Parameter name not specified for inversion &
      &parameter block.'
    call PrintErrMsg(option)
  endif

  InversionParameterRead => new_inversion_parameter

end function InversionParameterRead

! ************************************************************************** !

subroutine InversionParameterMapNametoInt(inversion_parameter,driver)
  !
  ! Maps an inverion parameter to subsurface model parameter id
  !
  ! Author: Glenn Hammond
  ! Date: 03/25/22
  !
  use Driver_module
  use String_module
  use Variables_module, only : ELECTRICAL_CONDUCTIVITY, &
                               PERMEABILITY, POROSITY, &
                               VG_SR, VG_ALPHA, VG_M

  type(inversion_parameter_type) :: inversion_parameter
  type(driver_type) :: driver

  PetscInt :: i

  select case(inversion_parameter%parameter_name)
    case('ELECTRICAL_CONDUCTIVITY')
      i = ELECTRICAL_CONDUCTIVITY
    case('PERMEABILITY')
      i = PERMEABILITY
    case('POROSITY')
      i = POROSITY
    case('ALPHA')
      i = VG_ALPHA
    case('RESIDUAL_SATURATION')
      i = VG_SR
    case('M')
      i = VG_M
    case default
      call driver%PrintErrMsg('Unrecognized parameter in &
                              &InversionParameterMap: ' // &
                              StringWrite(inversion_parameter%parameter_name))
  end select
  inversion_parameter%iparameter = i

end subroutine InversionParameterMapNametoInt

! ************************************************************************** !

function InversionParameterIntToQOIArray(inversion_parameter)
  !
  ! Maps an inverion parameter to subsurface model parameter id
  !
  ! Author: Glenn Hammond
  ! Date: 03/25/22
  !
  use String_module
  use Variables_module, only : POROSITY
  use Material_Aux_module, only : POROSITY_BASE

  type(inversion_parameter_type) :: inversion_parameter

  PetscInt :: InversionParameterIntToQOIArray(2)

  InversionParameterIntToQOIArray(1) = inversion_parameter%iparameter
  select case(inversion_parameter%iparameter)
    case(POROSITY)
      InversionParameterIntToQOIArray(2) = POROSITY_BASE
    case default
      InversionParameterIntToQOIArray(2) = ZERO_INTEGER
  end select

end function InversionParameterIntToQOIArray

! ************************************************************************** !

subroutine InversionParameterStrip(inversion_parameter)
  !
  ! Deallocates members of an inversion parameter object
  !
  ! Author: Glenn Hammond
  ! Date: 03/18/22
  !

  type(inversion_parameter_type) :: inversion_parameter

!  nullify(inversion_parameter%next)

end subroutine InversionParameterStrip

! ************************************************************************** !

recursive subroutine InversionParameterDestroy(inversion_parameter)
  !
  ! Deallocates a inversion parameter object
  !
  ! Author: Glenn Hammond
  ! Date: 03/18/22
  !
  type(inversion_parameter_type), pointer :: inversion_parameter

  if (.not.associated(inversion_parameter)) return

  if (associated(inversion_parameter%next)) then
    call InversionParameterDestroy(inversion_parameter%next)
  endif

  call InversionParameterStrip(inversion_parameter)

  deallocate(inversion_parameter)
  nullify(inversion_parameter)

end subroutine InversionParameterDestroy

end module Inversion_Parameter_module
