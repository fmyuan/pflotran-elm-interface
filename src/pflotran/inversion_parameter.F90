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
    PetscReal :: update
    character(len=MAXWORDLENGTH) :: parameter_name
    character(len=MAXWORDLENGTH) :: material_name
    type(inversion_parameter_type), pointer :: next
  end type inversion_parameter_type

  public :: InversionParameterCreate, &
            InversionParameterInit, &
            InversionParameterRead, &
            InversionParameterCopy, &
            InversionParameterMapNameToInt, &
            InversionParamGetItypeFromName, &
            InversionParamGetNameFromItype, &
            InversionParameterIntToQOIArray, &
            InversionParameterPrint, &
            InversionParameterPrintUpdate, &
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
  inversion_parameter%update = UNINITIALIZED_DOUBLE
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
  inversion_parameter2%update = inversion_parameter%update
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
    write(fid,'(/, &
              &" Current values of inversion parameters:",//, &
              &"      # &
              &Parameter Name                    &
              &Material Name       &
              & Value",/,&
              &"      - &
              &--------------                    &
              &-------------       &
              & -----")')
  endif
  write(string,'(i6," ",a32,2x,a20,es13.6)') &
    inversion_parameter%id, &
    inversion_parameter%parameter_name, &
    inversion_parameter%material_name, &
    inversion_parameter%value
  write(fid,*) trim(string)
  if (print_footer) then
!    write(fid,'(/,40("=+"))')
  endif

end subroutine InversionParameterPrint

! ************************************************************************** !

subroutine InversionParameterPrintUpdate(fid,inversion_parameter, &
                                         print_header,print_footer)
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

  character(len=MAXSTRINGLENGTH) :: string

  if (print_header) then
    write(fid,'(/, &
              &" Current values of inversion parameter updates:",//, &
              &"      # &
              &Parameter Name      &
              &Material Name       &
              & Update",/,&
              &"      - &
              &--------------      &
              &-------------       &
              & -----")')
  endif
  write(string,'(i6," ",2a20,es13.6)') &
    inversion_parameter%id, &
    inversion_parameter%parameter_name, &
    inversion_parameter%material_name, &
    inversion_parameter%update
  write(fid,*) trim(string)
  if (print_footer) then
!    write(fid,'(/,40("=+"))')
  endif

end subroutine InversionParameterPrintUpdate

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

subroutine InversionParameterMapNametoInt(inversion_parameter,driver, &
                                          inversion_option)
  !
  ! Maps an inversion parameter to subsurface model parameter id
  !
  ! Author: Glenn Hammond
  ! Date: 03/25/22
  !
  use Driver_class
  use Option_Inversion_module

  type(inversion_parameter_type) :: inversion_parameter
  class(driver_type) :: driver
  type(inversion_option_type) :: inversion_option

  inversion_parameter%iparameter = &
    InversionParamGetItypeFromName(inversion_parameter%parameter_name, &
                                    driver,inversion_option)

end subroutine InversionParameterMapNametoInt

! ************************************************************************** !

function InversionParamGetItypeFromName(name_,driver,inversion_option)
  !
  ! Maps an inversion parameter_name to subsurface model parameter id
  !
  ! Author: Glenn Hammond
  ! Date: 11/23/22
  !
  use Driver_class
  use Option_Inversion_module
  use Variables_module, only : ELECTRICAL_CONDUCTIVITY, &
                               PERMEABILITY, POROSITY, &
                               VG_SR, VG_ALPHA, VG_M, &
                               ARCHIE_CEMENTATION_EXPONENT, &
                               ARCHIE_SATURATION_EXPONENT, &
                               ARCHIE_TORTUOSITY_CONSTANT

  character(len=MAXWORDLENGTH) :: name_
  class(driver_type) :: driver
  type(inversion_option_type) :: inversion_option

  PetscInt :: InversionParamGetItypeFromName

  PetscInt :: i

  select case(name_)
    case('ELECTRICAL_CONDUCTIVITY')
      i = ELECTRICAL_CONDUCTIVITY
      inversion_option%invert_for_elec_cond = PETSC_TRUE
    case('PERMEABILITY')
      i = PERMEABILITY
      inversion_option%invert_for_permeability = PETSC_TRUE
    case('POROSITY')
      i = POROSITY
      inversion_option%invert_for_porosity = PETSC_TRUE
    case('ALPHA')
      i = VG_ALPHA
      inversion_option%invert_for_vg_alpha = PETSC_TRUE
    case('RESIDUAL_SATURATION')
      i = VG_SR
      inversion_option%invert_for_vg_sr = PETSC_TRUE
    case('M')
      i = VG_M
      inversion_option%invert_for_vg_m = PETSC_TRUE
    case('ARCHIE_CEMENTATION_EXPONENT')
      i = ARCHIE_CEMENTATION_EXPONENT
      inversion_option%invert_for_arch_cement_exp = PETSC_TRUE
    case('ARCHIE_SATURATION_EXPONENT')
      i = ARCHIE_SATURATION_EXPONENT
      inversion_option%invert_for_arch_sat_exp = PETSC_TRUE
    case('ARCHIE_TORTUOSITY_CONSTANT')
      i = ARCHIE_TORTUOSITY_CONSTANT
      inversion_option%invert_for_arch_tort_const = PETSC_TRUE
    case default
      call driver%PrintErrMsg('Unrecognized parameter in &
                              &InversionParamGetItypeFromName: ' // &
                              trim(name_))
  end select

  InversionParamGetItypeFromName = i

end function InversionParamGetItypeFromName

! ************************************************************************** !

function InversionParamGetNameFromItype(itype,driver,inversion_option)
  !
  ! Maps an inversion parameter id to subsurface model parameter name
  !
  ! Author: Glenn Hammond
  ! Date: 06/07/23
  !
  use Driver_class
  use Option_Inversion_module
  use String_module
  use Variables_module, only : ELECTRICAL_CONDUCTIVITY, &
                               PERMEABILITY, POROSITY, &
                               VG_SR, VG_ALPHA, VG_M, &
                               ARCHIE_CEMENTATION_EXPONENT, &
                               ARCHIE_SATURATION_EXPONENT, &
                               ARCHIE_TORTUOSITY_CONSTANT

  PetscInt :: itype
  class(driver_type) :: driver
  type(inversion_option_type) :: inversion_option

  character(len=MAXWORDLENGTH) :: InversionParamGetNameFromItype

  character(len=MAXWORDLENGTH) :: word

  select case(itype)
    case(ELECTRICAL_CONDUCTIVITY)
      word = 'ELECTRICAL_CONDUCTIVITY'
    case(PERMEABILITY)
      word = 'PERMEABILITY'
    case(POROSITY)
      word = 'POROSITY'
    case(VG_ALPHA)
      word = 'ALPHA'
    case(VG_SR)
      word = 'RESIDUAL_SATURATION'
    case(VG_M)
      word = 'M'
    case(ARCHIE_CEMENTATION_EXPONENT)
      word = 'ARCHIE_CEMENTATION_EXPONENT'
    case(ARCHIE_SATURATION_EXPONENT)
      word = 'ARCHIE_SATURATION_EXPONENT'
    case(ARCHIE_TORTUOSITY_CONSTANT)
      word = 'ARCHIE_TORTUOSITY_CONSTANT'
    case default
      call driver%PrintErrMsg('Unrecognized parameter in &
                              &InversionParamGetNameFromItype: ' // &
                              StringWrite(itype))
  end select

  InversionParamGetNameFromItype = word

end function InversionParamGetNameFromItype

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
