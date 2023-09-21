module Inversion_Parameter_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Variables_module, only : ELECTRICAL_CONDUCTIVITY, &
                               PERMEABILITY, POROSITY, &
                               VG_SR, VG_ALPHA, VG_M, &
                               ARCHIE_CEMENTATION_EXPONENT, &
                               ARCHIE_SATURATION_EXPONENT, &
                               ARCHIE_TORTUOSITY_CONSTANT, &
                               SURFACE_ELECTRICAL_CONDUCTIVITY, &
                               WAXMAN_SMITS_CLAY_CONDUCTIVITY, &
                               VERTICAL_PERM_ANISOTROPY_RATIO

  implicit none

  private

  PetscInt, parameter :: MAP_ELEC_COND = 1
  PetscInt, parameter :: MAP_PERM = 2
  PetscInt, parameter :: MAP_POR = 3
  PetscInt, parameter :: MAP_VG_SR = 4
  PetscInt, parameter :: MAP_VG_ALPHA = 5
  PetscInt, parameter :: MAP_VG_M = 6
  PetscInt, parameter :: MAP_ACE = 7
  PetscInt, parameter :: MAP_ASE = 8
  PetscInt, parameter :: MAP_ATC = 9
  PetscInt, parameter :: MAP_SURF_ELEC_COND = 10
  PetscInt, parameter :: MAP_WS_CLAY_COND = 11
  PetscInt, parameter :: MAP_VERT_PERM_ANISO_RATIO = 12
  ! MAX_MAP must be updated when new MAP_XXX parameters are added
  PetscInt, parameter :: MAX_MAP = 12

  ! second index is MAX_MAP above
  PetscReal :: parameter_bounds(2,MAX_MAP) = UNINITIALIZED_DOUBLE

  type, public :: inversion_parameter_type
    PetscInt :: id
    PetscInt :: itype
    PetscInt :: imat
    PetscReal :: value
    PetscReal :: update
    PetscReal :: bounds(2)
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
            InversionParamInitBounds, &
            InversionParamSetGlobalBounds, &
            InversionParamGetGlobalBounds, &
            InversionParameterBoundParameter, &
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
  inversion_parameter%itype = UNINITIALIZED_INTEGER
  inversion_parameter%imat = UNINITIALIZED_INTEGER
  inversion_parameter%value = UNINITIALIZED_DOUBLE
  inversion_parameter%update = UNINITIALIZED_DOUBLE
  inversion_parameter%bounds = UNINITIALIZED_DOUBLE
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
  inversion_parameter2%itype = inversion_parameter%itype
  inversion_parameter2%imat = inversion_parameter%imat
  inversion_parameter2%value = inversion_parameter%value
  inversion_parameter2%update = inversion_parameter%update
  inversion_parameter2%bounds = inversion_parameter%bounds
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
  write(string,'(i6," ",a32,2x,a20,2es13.6)') &
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
              &Parameter Name                    &
              &Material Name       &
              & Update",/,&
              &"      - &
              &--------------                    &
              &-------------       &
              & -----")')
  endif
  write(string,'(i6," ",a32,2x,a20,2es13.6)') &
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
      case('INITIAL_VALUE')
        call InputReadDouble(input,option,new_inversion_parameter%value)
        call InputErrorMsg(input,option,keyword,error_string)
      case('BOUNDS')
        call InputReadDouble(input,option,new_inversion_parameter%bounds(1))
        call InputErrorMsg(input,option,'BOUNDS,LOWER_BOUND',error_string)
        call InputReadDouble(input,option,new_inversion_parameter%bounds(2))
        call InputErrorMsg(input,option,'BOUNDS,UPPER_BOUND',error_string)
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
  ! Maps an inversion parameter to subsurface model parameter id
  !
  ! Author: Glenn Hammond
  ! Date: 03/25/22
  !
  use Driver_class

  type(inversion_parameter_type) :: inversion_parameter
  class(driver_type) :: driver

  inversion_parameter%itype = &
    InversionParamGetItypeFromName(inversion_parameter%parameter_name,driver)

end subroutine InversionParameterMapNametoInt

! ************************************************************************** !

function InversionParamGetItypeFromName(name_,driver)
  !
  ! Maps an inversion parameter_name to subsurface model parameter id
  !
  ! Author: Glenn Hammond
  ! Date: 11/23/22
  !
  use Driver_class

  character(len=MAXWORDLENGTH) :: name_
  class(driver_type) :: driver

  PetscInt :: InversionParamGetItypeFromName

  PetscInt :: i

  select case(name_)
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
    case('ARCHIE_CEMENTATION_EXPONENT')
      i = ARCHIE_CEMENTATION_EXPONENT
    case('ARCHIE_SATURATION_EXPONENT')
      i = ARCHIE_SATURATION_EXPONENT
    case('ARCHIE_TORTUOSITY_CONSTANT')
      i = ARCHIE_TORTUOSITY_CONSTANT
    case('SURFACE_ELECTRICAL_CONDUCTIVITY')
      i = SURFACE_ELECTRICAL_CONDUCTIVITY
    case('WAXMAN_SMITS_CLAY_CONDUCTIVITY')
      i = WAXMAN_SMITS_CLAY_CONDUCTIVITY
    case('VERTICAL_PERM_ANISOTROPY_RATIO')
      i = VERTICAL_PERM_ANISOTROPY_RATIO
    case default
      call driver%PrintErrMsg('Unrecognized parameter in &
                              &InversionParamGetItypeFromName: ' // &
                              trim(name_))
  end select

  InversionParamGetItypeFromName = i

end function InversionParamGetItypeFromName

! ************************************************************************** !

function InversionParamGetNameFromItype(itype,driver)
  !
  ! Maps an inversion parameter id to subsurface model parameter name
  !
  ! Author: Glenn Hammond
  ! Date: 06/07/23
  !
  use Driver_class
  use String_module

  PetscInt :: itype
  class(driver_type) :: driver

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
    case(SURFACE_ELECTRICAL_CONDUCTIVITY)
      word = 'SURFACE_ELECTRICAL_CONDUCTIVITY'
    case(WAXMAN_SMITS_CLAY_CONDUCTIVITY)
      word = 'WAXMAN_SMITS_CLAY_CONDUCTIVITY'
    case(VERTICAL_PERM_ANISOTROPY_RATIO)
      word = 'VERTICAL_PERM_ANISOTROPY_RATIO'
    case default
      call driver%PrintErrMsg('Unrecognized parameter in &
                              &InversionParamGetNameFromItype: ' // &
                              StringWrite(itype))
  end select

  InversionParamGetNameFromItype = word

end function InversionParamGetNameFromItype

! ************************************************************************** !

function InvParamItypeToItypeInternal(itype)
  !
  ! Maps an inversion parameter id to the internal id in this module
  !
  ! Author: Glenn Hammond
  ! Date: 06/14/23
  !
  PetscInt :: itype

  PetscInt :: InvParamItypeToItypeInternal

  PetscInt :: i

  select case(itype)
    case(ELECTRICAL_CONDUCTIVITY)
      i = MAP_ELEC_COND
    case(PERMEABILITY)
      i = MAP_PERM
    case(POROSITY)
      i = MAP_POR
    case(VG_ALPHA)
      i = MAP_VG_ALPHA
    case(VG_SR)
      i = MAP_VG_SR
    case(VG_M)
      i = MAP_VG_M
    case(ARCHIE_CEMENTATION_EXPONENT)
      i = MAP_ACE
    case(ARCHIE_SATURATION_EXPONENT)
      i = MAP_ASE
    case(ARCHIE_TORTUOSITY_CONSTANT)
      i = MAP_ATC
    case(SURFACE_ELECTRICAL_CONDUCTIVITY)
      i = MAP_SURF_ELEC_COND
    case(WAXMAN_SMITS_CLAY_CONDUCTIVITY)
      i = MAP_WS_CLAY_COND
    case(VERTICAL_PERM_ANISOTROPY_RATIO)
      i = MAP_VERT_PERM_ANISO_RATIO
    case default
      print *, 'Unrecognized parameter in &
               &InvParamItypeToItypeInternal: ', itype
      stop
  end select

  InvParamItypeToItypeInternal = i

end function InvParamItypeToItypeInternal

! ************************************************************************** !

subroutine InversionParamInitBounds()
  !
  ! Sets the global upper and lower bounds for each variable
  !
  ! Author: Glenn Hammond
  ! Date: 06/14/23
  !
  PetscReal, parameter :: default_lower_bound = 0.d0
  PetscReal, parameter :: default_upper_bound = 1.d20

  parameter_bounds(:,:) = UNINITIALIZED_DOUBLE
  call InversionParamSetGlobalBounds(ELECTRICAL_CONDUCTIVITY, &
                                     default_lower_bound,default_upper_bound)
  call InversionParamSetGlobalBounds(PERMEABILITY,1.d-30,1.d-7)
  call InversionParamSetGlobalBounds(POROSITY,0.d0,1.d0)
  call InversionParamSetGlobalBounds(VG_ALPHA,1.d-6,1.d-3)
  call InversionParamSetGlobalBounds(VG_SR,0.d0,0.6d0) ! not based on data
  call InversionParamSetGlobalBounds(VG_M,0.2d0,0.9d0)
  call InversionParamSetGlobalBounds(ARCHIE_CEMENTATION_EXPONENT, &
                                     default_lower_bound,default_upper_bound)
  call InversionParamSetGlobalBounds(ARCHIE_SATURATION_EXPONENT, &
                                     default_lower_bound,default_upper_bound)
  call InversionParamSetGlobalBounds(ARCHIE_TORTUOSITY_CONSTANT, &
                                     default_lower_bound,default_upper_bound)
  call InversionParamSetGlobalBounds(SURFACE_ELECTRICAL_CONDUCTIVITY, &
                                     0.d0,0.1d0)
  call InversionParamSetGlobalBounds(WAXMAN_SMITS_CLAY_CONDUCTIVITY, &
                                     0.d0,0.1d0)
  call InversionParamSetGlobalBounds(VERTICAL_PERM_ANISOTROPY_RATIO, &
                                     0.01d0,100.d0)

end subroutine InversionParamInitBounds

! ************************************************************************** !

subroutine InversionParamSetGlobalBounds(itype,lower_bound,upper_bound)
  !
  ! Sets the global upper and lower bounds for each variable
  !
  ! Author: Glenn Hammond
  ! Date: 06/14/23
  !

  PetscInt :: itype
  PetscReal :: lower_bound
  PetscReal :: upper_bound

  parameter_bounds(:,InvParamItypeToItypeInternal(itype)) = &
    [lower_bound,upper_bound]

end subroutine InversionParamSetGlobalBounds

! ************************************************************************** !

subroutine InversionParamGetGlobalBounds(itype,lower_bound,upper_bound)
  !
  ! Gets the global upper and lower bounds for each variable
  !
  ! Author: Glenn Hammond
  ! Date: 06/14/23
  !

  PetscInt :: itype
  PetscReal :: lower_bound
  PetscReal :: upper_bound

  PetscInt i

  i = InvParamItypeToItypeInternal(itype)
  lower_bound = parameter_bounds(1,i)
  upper_bound = parameter_bounds(2,i)

end subroutine InversionParamGetGlobalBounds

! ************************************************************************** !

subroutine InversionParameterBoundParameter(inversion_parameter,value)
  !
  ! Truncates parameter to within the lower and upper bound
  !
  ! Author: Glenn Hammond
  ! Date: 06/14/23
  !
  type(inversion_parameter_type) :: inversion_parameter
  PetscReal :: value

  PetscReal :: lower_bound
  PetscReal :: upper_bound

  if (Initialized(inversion_parameter%bounds(1))) then
    lower_bound = inversion_parameter%bounds(1)
    upper_bound = inversion_parameter%bounds(2)
  else
    call InversionParamGetGlobalBounds(inversion_parameter%itype, &
                                       lower_bound,upper_bound)
  endif
  value = max(min(value,upper_bound),lower_bound)

end subroutine InversionParameterBoundParameter

! ************************************************************************** !

function InversionParameterIntToQOIArray(inversion_parameter)
  !
  ! Maps an inverion parameter to subsurface model parameter id
  !
  ! Author: Glenn Hammond
  ! Date: 03/25/22
  !
  use String_module
  use Material_Aux_module, only : POROSITY_BASE

  type(inversion_parameter_type) :: inversion_parameter

  PetscInt :: InversionParameterIntToQOIArray(2)

  InversionParameterIntToQOIArray(1) = inversion_parameter%itype
  select case(inversion_parameter%itype)
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
