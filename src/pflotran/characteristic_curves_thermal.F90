module Characteristic_Curves_Thermal_module
  !
  ! thermal conductivity as a function of saturation and temperature
  !
  ! Author: Kris Kuhlman
  ! Date: 05/04/20
  !

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  !---------------------------------------------------------------------------
  type, public :: thermal_conductivity_base_type
    PetscBool :: analytical_derivative_available
  contains
    procedure, public :: Init => TCFBaseInit
    procedure, public :: Verify => TCFBaseVerify
    procedure, public :: Test => TCFBaseTest
    procedure, public :: Conductivity => TCFBaseConductivity
  end type thermal_conductivity_base_type
  !---------------------------------------------------------------------------
  type, public, extends(thermal_conductivity_base_type) :: kT_default_type
    PetscReal :: kT_wet, kT_dry
  contains
    procedure, public :: Verify => TCFDefaultVerify
    procedure, public :: Conductivity => TCFDefaultConductivity
  end type kT_default_type
  !---------------------------------------------------------------------------
  type, public, extends(thermal_conductivity_base_type) :: kT_constant_type
    PetscReal :: constant_thermal_conductivity
  contains
    procedure, public :: Verify => TCFConstantVerify
    procedure, public :: Conductivity => TCFConstantConductivity
  end type kT_constant_type
  !---------------------------------------------------------------------------
  type, public, extends(kT_default_type) :: kT_power_type
    PetscReal :: ref_temp
    PetscReal :: gamma
  contains
    procedure, public :: Init => TCFPowerInit
    procedure, public :: Verify => TCFPowerVerify
    procedure, public :: Conductivity => TCFPowerConductivity
  end type kT_power_type
  !---------------------------------------------------------------------------
  type, public, extends(kT_default_type) :: kT_cubic_polynomial_type
    PetscReal :: ref_temp
    PetscReal(3) :: beta
  contains
    procedure, public :: Init => TCFCubicPolyInit
    procedure, public :: Verify => TCFCubicPolyVerify
    procedure, public :: Conductivity => TCFCubicPolyConductivity
  end type kT_cubic_polynomial_type
  !---------------------------------------------------------------------------
  type, public :: cc_thermal_type
    character(len=MAXWORDLENGTH) :: name
    PetscBool :: print_me
    PetscBool :: test
    class(thermal_conductivity_base_type), pointer :: thermal_conductivity_function
    class(cc_thermal_type), pointer :: next
  end type cc_thermal_type
  !---------------------------------------------------------------------------
  type, public :: cc_thermal_ptr_type
    class(cc_thermal_type), pointer :: ptr
  end type cc_thermal_ptr_type

contains

! ************************************************************************** !
! ************************************************************************** !

subroutine TCFBaseVerify(this,name,option)

  use Option_module

  implicit none

  class(thermal_conductivity_func_base_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  if ((.not.this%analytical_derivative_available) .and. &
       (.not.option%flow%numerical_derivatives)) then
    option%io_buffer = 'Analytical derivatives are not available for the &
         &thermal conductivity function chosen: ' // &
         trim(name)
    call PrintErrMsg(option)
  endif

end subroutine TCFBaseVerify

! ************************************************************************** !

subroutine TCFBaseThermalConductivity(this,liquid_saturation,temperature, &
     thermal_conductivity,dk_dsatl,dk_dtemp,option)

  use Option_module

  implicit none

  class(thermal_conductivity_func_base_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dk_dsatl, dk_dtemp
  type(option_type), intent(inout) :: option

  option%io_buffer = 'TCFBaseThermalConductivity must be extended.'
  call PrintErrMsg(option)

end subroutine TCFBaseThermalConductivity

! ************************************************************************** !

subroutine TCFBaseTest(this,tcc_name,option)

  use Option_module
  use Material_Aux_class

  implicit none

  class(sat_func_base_type) :: this
  character(len=MAXWORDLENGTH) :: tcc_name
  type(option_type), intent(inout) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt, parameter :: nt = 51
  PetscInt, parameter :: ns = 11
  PetscReal :: temp, deltaTemp, sat, deltaSat
  PetscReal :: temp_vec(nt)
  PetscReal :: sat_vec(ns)
  PetscReal :: kT(nt,ns)
  PetscReal :: dkT_dsat(nt,ns)
  PetscReal :: dkT_dsat_numerical(nt,ns)
  PetscReal :: dkT_dT(nt,ns)
  PetscReal :: dkT_dT_numerical(nt,ns)
  PetscReal :: temp_pert, sat_pert, unused
  PetscReal :: temp_min, temp_max, sat_min, sat_max
  PetscInt :: i,j

  ! thermal conductivity as a function of temp. and liq. sat.
  temp_min = 0.0d0 ! Celsius
  temp_max = 250.0d0
  sat_min = 0.0d0 ! relative saturation
  sat_max = 1.0d0

  deltaTemp = (temp_max - temp_min)/(nt - 1)
  deltaSat = (sat_max - sat_min)/(ns - 1)

  temp_vec = [(temp_min + i*deltaTemp, i,0,nt-1)]
  sat_vec = [(sat_min + i*deltaSat, i,0,ns-1)]

  dTemp = 1.0D-6 ! perturbations
  dSat = 1.0D-6

  do i = 1,nt
    do j = 1,ns
      call this%Conductivity(sat_vec(j),temp_vec(i), &
           kT(i,j),dkT_dsat(i,j),dkT_dT(i,j),option)

      ! calculate numerical derivative dkT_dsatl_numerical
      call this%Conductivity(satl_vec(j),temp_vec(i)+dTemp, &
           temp_pert,unused,unused,option)

      dkT_dT_numerical(i,j) = (kT(i,j) - temp_pert)/(temp_vec(i)*dTemp)

      call this%Conductivity(sat_vec(j)+dSat,temp_vec(i), &
           sat_pert,unused,unused,option)

      dkT_dsat_numerical(i,j) = (kT(i,j) - sat_pert)/(sat_vec(j)*dSat)
    end do
  end do

  write(string,*) tcc_name
  string = trim(cc_name) // '_kT_vs_sat_and_temp.dat'
  open(unit=86,file=string)
  write(86,*) '"temperature [C]", "liquid saturation [-]", "kT [W/m*K]", &
       &"dkT/dsat", "dkT/dT", "dkT/dsat_numerical", "dkT/dT_numerical"'
  do i = 1,nt
    do j = 1,ns
      write(86,'(7(ES14.6))') temp_vec(i), sat_vec(j), &
           kT(i,j), dkT_dsat(i,j), dkT_dT(i,j), &
           dkT_dsat_numerical(i,j), dkT_dT_numerical(i,j)
    end do
  end do
  close(86)

end subroutine TCFBaseTest

! ************************************************************************** !

subroutine ThermalConductivityFunctionDestroy(tcf)

  implicit none

  class(thermal_conductivity_base_type), pointer :: tcf

  if (.not.associated(tcf)) return
  deallocate(tcf)
  nullify(tcf)

end subroutine ThermalConductivityFunctionDestroy

! ************************************************************************** !
! ************************************************************************** !

function TCF_Default_Create()

  implicit none

  class(kT_default_type), pointer :: TCF_Default_Create

  allocate(TCF_Default_Create)
  TCF_Default_Create%kT_wet = 0.d0
  TCF_Default_Create%kT_dry = 0.d0
  TCF_Default_Create%analytical_derivative_available = PETSC_TRUE

end function TCF_Default_Create

! ************************************************************************** !

subroutine TCFDefaultVerify(this,name,option)

  use Option_module

  implicit none

  class(kT_default_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'THERMAL_CONDUCTIVITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'THERMAL_CONDUCTIVITY_FUNCTION,DEFAULT'
  endif
  call TCFBaseVerify(this,string,option)
  if (Uninitialized(this%kT_wet)) then
    option%io_buffer = UninitializedMessage('THERMAL_CONDUCTIVITY_WET',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%kT_dry)) then
    option%io_buffer = UninitializedMessage('THERMAL_CONDUCTIVITY_DRY',string)
    call PrintErrMsg(option)
  endif

end subroutine TCFDefaultVerify

! ************************************************************************** !

subroutine TCFDefaultConductivity(this,liquid_saturation,temperature, &
     thermal_conductivity,dk_dsatl,dk_dtemp,option)

  use Option_module

  implicit none

  class(kT_default_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dk_dsatl, dk_dtemp
  type(option_type), intent(inout) :: option

  PetscReal :: tempreal

  ! based on Somerton et al., 1974:
  ! k_eff = k_dry + sqrt(s_l)*(k_wet-k_dry)

  dk_dtemp = 0.d0 ! only a function of saturation

  if (liquid_saturation > 0.d0) then
    tempreal = sqrt(liquid_saturation) * (this%kT_wet - this%kT_dry)
    thermal_conductivity = this%kT_dry + tempreal
    dk_dsatl = 0.5d0 * tempreal / liquid_saturation
  else
    thermal_conductivity = this%kT_dry
    dk_dsatl = 0.d0
  end if

end subroutine TCFDefaultConductivity

! ************************************************************************** !
! ************************************************************************** !

function TCF_Constant_Create()

  implicit none

  class(kT_constant_type), pointer :: TCF_Constant_Create

  allocate(TCF_Constant_Create)
  TCF_Constant_Create%constant_thermal_conductivity = 0.d0
  TCF_Constant_Create%analytical_derivative_available = PETSC_TRUE

end function TCF_Constant_Create

! ************************************************************************** !

subroutine TCFConstantVerify(this,name,option)

  use Option_module

  implicit none

  class(kT_constant_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'THERMAL_CONDUCTIVITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'THERMAL_CONDUCTIVITY_FUNCTION,CONSTANT'
  endif
  call TCFBaseVerify(this,string,option)
  if (Uninitialized(this%constant_thermal_conductivity)) then
    option%io_buffer = UninitializedMessage( &
         'THERMAL_CONDUCTIVITY_CONSTANT',string)
    call PrintErrMsg(option)
  endif

end subroutine TCFConstantVerify

! ************************************************************************** !

subroutine TCFDefaultConductivity(this,liquid_saturation,temperature, &
     thermal_conductivity,dk_dsatl,dk_dtemp,option)

  use Option_module

  implicit none

  class(kT_constant_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dk_dsatl, dk_dtemp
  type(option_type), intent(inout) :: option

  thermal_conductivity = this%constant_thermal_conductivity
  dk_dsatl = 0.d0
  dk_dtemp = 0.d0

end subroutine TCFDefaultConductivity

! ************************************************************************** !
! ************************************************************************** !

function TCF_Power_Create()

  implicit none

  class(kT_power_type), pointer :: TCF_Power_Create

  allocate(TCF_Power_Create)
  TCF_Power_Create%kT_wet = 0.d0
  TCF_Power_Create%kT_dry = 0.d0
  TCF_Power_Create%gamma = 0.0d0
  TCF_Power_Create%analytical_derivative_available = PETSC_TRUE

end function TCF_Power_Create

! ************************************************************************** !

subroutine TCFPowerVerify(this,name,option)

  use Option_module

  implicit none

  class(kT_default_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'THERMAL_CONDUCTIVITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'THERMAL_CONDUCTIVITY_FUNCTION,POWER'
  endif
  call TCFBaseVerify(this,string,option)
  if (Uninitialized(this%kT_wet)) then
    option%io_buffer = UninitializedMessage('THERMAL_CONDUCTIVITY_WET',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%kT_dry)) then
    option%io_buffer = UninitializedMessage('THERMAL_CONDUCTIVITY_DRY',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%ref_temp)) then
    option%io_buffer = UninitializedMessage('REFERENCE_TEMPERATURE',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%gamma)) then
    option%io_buffer = UninitializedMessage('EXPONENT',string)
    call PrintErrMsg(option)
  endif

end subroutine TCFPowerVerify

! ************************************************************************** !

subroutine TCFPowerConductivity(this,liquid_saturation,temperature, &
     thermal_conductivity,dk_dsatl,dk_dtemp,option)

  use Option_module

  implicit none

  class(kT_power_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dk_dsatl, dk_dtemp
  type(option_type), intent(inout) :: option

  class(kT_default_type) :: that
  PetscReal :: tempreal, kT_base, shifted_temp

  ! behavior WRT saturation from default function
  call TCFDefaultConductivity(this%kT_default_type, &
       liquid_saturation,0.d0,kT_base,dk_dsatl,dk_dtemp,option)

  shifted_temp = temperature - this%ref_temp

  tempreal = kT_base * shifted_temp ** this%gamma
  thermal_conductivity = tempreal
  dk_temp = this%gamma * tempreal / shifted_temp

end subroutine TCFPowerConductivity

! ************************************************************************** !
! ************************************************************************** !

function TCF_Cubic_Polynomial_Create()

  implicit none

  class(kT_cubic_polynomial_type), pointer :: TCF_Cubic_Polynomial_Create

  allocate(TCF_Cubic_Polynomial_Create)
  TCF_Cubic_Polynomial_Create%kT_wet = 0.d0
  TCF_Cubic_Polynomial_Create%kT_dry = 0.d0
  TCF_Cubic_Polynomial_Create%ref_temp = 0.d0
  TCF_Cubic_Polynomial_Create%beta = [ 0.d0, 0.d0, 0.d0 ]
  TCF_Cubic_Polynomial_Create%analytical_derivative_available = PETSC_TRUE

end function TCF_Cubic_Polynomial_Create

! ************************************************************************** !

subroutine TCFCubicPolyVerify(this,name,option)

  use Option_module

  implicit none

  class(kT_default_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'THERMAL_CONDUCTIVITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'THERMAL_CONDUCTIVITY_FUNCTION,CUBIC_POLYNOMIAL'
  endif
  call TCFBaseVerify(this,string,option)
  if (Uninitialized(this%kT_wet)) then
    option%io_buffer = UninitializedMessage('THERMAL_CONDUCTIVITY_WET',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%kT_dry)) then
    option%io_buffer = UninitializedMessage('THERMAL_CONDUCTIVITY_DRY',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%ref_temp)) then
    option%io_buffer = UninitializedMessage('COEFFCIENTS',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%beta)) then
    option%io_buffer = UninitializedMessage('COEFFCIENTS',string)
    call PrintErrMsg(option)
  endif

end subroutine TCFCubicPolyVerify

! ************************************************************************** !

subroutine TCFCubicPolyConductivity(this,liquid_saturation,temperature, &
     thermal_conductivity,dk_dsatl,dk_dtemp,option)

  use Option_module

  implicit none

  class(kT_default_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dk_dsatl, dk_dtemp
  type(option_type), intent(inout) :: option

  PetscReal :: kT_base, shifted_temp

  ! behavior WRT saturation from default function
  call TCFDefaultConductivity(this%kT_default_type, &
       liquid_saturation,0.d0,kT_base,dk_dsatl,dk_dtemp,option)

  shifted_temp = temperature - this%ref_temp

  ! kT = k0 * (1 + beta1*T + beta2*T^2 + beta3*T^3)
  thermal_conductivity = kT_base * &
       (1.d0 + shifted_temp*(this%beta(1) &
       + shifted_temp*(this%beta(2) &
       + shifted_temp* this%beta(3))))

  dk_temp = kT_base * (this%beta(1) &
       + shifted_temp*(2.d0*this%beta(2) &
       + shifted_temp* 3.d0*this%beta(3)))

end subroutine TCFCubicPolyConductivity

! ************************************************************************** !
! ************************************************************************** !

function CharacteristicCurvesThermalCreate()

  implicit none

  class(cc_thermal_type), pointer :: CharacteristicCurvesThermalCreate
  class(cc_thermal_type), pointer :: characteristic_thermal_curves

  allocate(characteristic_thermal_curves)
  characteristic_thermal_curves%name = ''
  characteristic_thermal_curves%print_me = PETSC_FALSE
  characteristic_thermal_curves%test = PETSC_FALSE
  nullify(characteristic_thermal_curves%thermal_conductivity_function)
  nullify(characteristic_thermal_curves%next)

  CharacteristicCurvesThermalCreate => characteristic_thermal_curves

end function CharacteristicCurvesThermalCreate

! ************************************************************************** !

subroutine CharacteristicCurvesThermalRead(this,input,option)

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  class(cc_thermal_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: error_string, verify_string
  class(thermal_conductivity_base_type), &
       pointer :: thermal_conductivity_function_ptr

  nullify(thermal_conductivity_function_ptr)

  input%ierr = 0
  error_string = 'THERMAL_CHARACTERISTIC_CURVES'
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    select case(trim(keyword))
    case('THERMAL_CONDUCTIVITY_FUNCTION')
      call InputErrorMsg(input,option, &
           'THERMAL_CONDUCTIVITY_FUNCTION',error_string)
      call StringToUpper(word)
      select case(word)
      case('CONSTANT')
        this%thermal_conductivity_function => TCF_Constant_Create()
      case('DEFAULT')
        this%thermal_conductivity_function => TCF_Default_Create()
      case('POWER')
        this%thermal_conductivity_function => RCF_Power_Create()
      case('CUBIC_POLYNOMIAL')
        this%thermal_conductivity_function => TCF_Cubic_Polynomial_Create()
      case default
        call InputKeywordUnrecognized(input,word, &
             'THERMAL_CONDUCTIVITY_FUNCTION',option)
      end select
      call ThermalConductivityFunctionRead( &
           this%thermal_conductivity_function,input,option)
    case('TEST')
      this%test = PETSC_TRUE
    case default
      call InputKeywordUnrecognized(input,keyword, &
           'THERMAL_CHARACTERISTIC_CURVES',option)
    end select
  enddo
  call InputPopBlock(input,option)

  verify_string = 'THERMAL_CHARACTERISTIC_CURVES(' // trim(this%name) // '),'

  if (associated(this%thermal_conductivity_function)) then
    call this%thermal_conductivity_function%Verify(verify_string,option)
  else
    option%io_buffer = 'A thermal conductivity function has &
         &not been set under THERMAL_CHARACTERISTIC_CURVES "' // &
         trim(this%name) // '". A THERMAL_CONDUCTIVITY_FUNCTION &
         &block must be specified.'
  endif

end subroutine CharacteristicCurvesThermalRead

!-----------------------------------------------------------------------

subroutine ThermalConductivityFunctionRead( &
     thermal_conductivity_function,input,option)
  !
  ! Reads in contents of a THERMAL_CONDUCTIVITY_FUNCTION block
  !
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  class(thermal_conductivity_base_type) :: thermal_conductivity_function
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  input%ierr = 0
  error_string = 'THERMAL_CHARACTERISTIC_CURVES,&
       &THERMAL_CONDUCTIVITY_FUNCTION,'
  select type(tcf => saturation_function)
  class is(kT_constant_type)
    error_string = trim(error_string) // 'CONSTANT'
  class is(kT_default_type)
    error_string = trim(error_string) // 'DEFAULT'
  class is(kT_power_type)
    error_string = trim(error_string) // 'POWER'
  class is(kT_cubic_polynomial_type)
    error_string = trim(error_string) // 'CUBIC_POLYNOMIAL'
  end select

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    select type(tcf => thermal_conductivity_function)
      !------------------------------------------
    class is(kT_constant_type)
      select case(keyword)
      case('CONSTANT_THERMAL_CONDUCTIVITY')
        call InputReadDouble(input,option,tcf%constant_thermal_conductivity)
        call InputErrorMsg(input,option,'constant thermal conductivity', &
             error_string)
      case default
        call InputKeywordUnrecognized(input,keyword, &
             'constant thermal conductivity',option)
      end select
      !------------------------------------------
    class is(kT_default_type)
      select case(keyword)
      case('THERMAL_CONDUCTIVITY_WET')
        call InputReadDouble(input,option,tcf%kT_wet)
        call InputErrorMsg(input,option,'thermal conductivity wet', &
             error_string)
      case('THERMAL_CONDUCTIVITY_DRY')
        call InputReadDouble(input,option,tcf%kT_dry)
        call InputErrorMsg(input,option,'thermal conductivity dry', &
             error_string)
      case default
        call InputKeywordUnrecognized(input,keyword, &
             'saturation-dependent thermal conductivity',option)
      end select
      !------------------------------------------
    class is(kT_power_type)
      select case(keyword)
      case('THERMAL_CONDUCTIVITY_WET')
        call InputReadDouble(input,option,tcf%kT_wet)
        call InputErrorMsg(input,option,'thermal conductivity wet', &
             error_string)
      case('THERMAL_CONDUCTIVITY_DRY')
        call InputReadDouble(input,option,tcf%kT_dry)
        call InputErrorMsg(input,option,'thermal conductivity dry', &
             error_string)
      case('REFERENCE_TEMPERATURE')
        call InputReadDouble(input,option,tcf%ref_temp)
        call InputErrorMsg(input,option,'reference temperature', &
             error_string)
      case('EXPONENT')
        call InputReadDouble(input,option,tcf%gamma)
        call InputErrorMsg(input,option,'thermal conductivity exponent', &
             error_string)
      case default
        call InputKeywordUnrecognized(input,keyword, &
             'temp-dependent (power) thermal conductivity',option)
      end select
      !------------------------------------------
    class is(kT_cubic_polynomial_type)
      select case(keyword)
      case('THERMAL_CONDUCTIVITY_WET')
        call InputReadDouble(input,option,tcf%kT_wet)
        call InputErrorMsg(input,option,'thermal conductivity wet', &
             error_string)
      case('THERMAL_CONDUCTIVITY_DRY')
        call InputReadDouble(input,option,tcf%kT_dry)
        call InputErrorMsg(input,option,'thermal conductivity dry', &
             error_string)
      case('REFERENCE_TEMPERATURE')
        call InputReadDouble(input,option,tcf%ref_temp)
        call InputErrorMsg(input,option,'reference temperature', &
             error_string)
      case('CUBIC_POLYNOMIAL_COEFFICIENTS')
        call InputReadNDoubles(input,option,tcf%beta,3)
        call InputErrorMsg(input,option, &
             'thermal conductivity polynomial coefficients',error_string)
      case default
        call InputKeywordUnrecognized(input,keyword, &
             'temp-dependent (cubic polynomial) thermal conductivity',option)
      end select
    class default
      option%io_buffer = 'Read routine not implemented for ' &
           // trim(error_string) // '.'
      call PrintErrMsg(option)
    end select
  end do
  call InputPopBlock(input,option)

end subroutine ThermalConductivityFunctionRead

! ************************************************************************** !

subroutine CharacteristicCurvesThermalAddToList(new_thermal_cc,list)

  implicit none

  class(cc_thermal_type), pointer :: new_thermal_cc
  class(cc_thermal_type), pointer :: list

  class(cc_thermal_type), pointer :: cur_thermal_cc

  if (associated(list)) then
    cur_thermal_cc => list
    ! loop to end of list
    do
      if (.not.associated(cur_thermal_cc%next)) exit
      cur_thermal_cc => cur_thermal_cc%next
    enddo
    cur_thermal_cc%next => new_thermal_cc
  else
    list => new_thermal_cc
  endif

end subroutine CharacteristicCurvesThermalAddToList

! ************************************************************************** !

subroutine CharCurvesThermalConvertListToArray(list,array,option)

  use String_module
  use Option_module

  implicit none

  class(cc_thermal_type), pointer :: list
  type(cc_thermal_ptr_type), pointer :: array(:)
  type(option_type) :: option

  class(cc_thermal_type), pointer :: cur_thermal_cc
  PetscInt :: count

  count = 0
  cur_thermal_cc => list
  do
    if (.not.associated(cur_thermal_cc)) exit
    count = count + 1
    cur_thermal_cc => cur_thermal_cc%next
  enddo

  if (associated(array)) deallocate(array)
  allocate(array(count))

  count = 0
  cur_thermal_cc => list
  do
    if (.not.associated(cur_thermal_cc)) exit
    count = count + 1
    array(count)%ptr => cur_thermal_cc
    if (cur_thermal_cc%test) then
      call OptionSetBlocking(option,PETSC_FALSE)
      if (option%myrank == option%io_rank) then
        if (associated(cur_thermal_cc%thermal_conductivity_function)) then
          call cur_thermal_cc%thermal_conductivity_function%Test( &
               cur_thermal_cc%name,option)
        end if
        call CharacteristicCurvesThermalTest(cur_thermal_cc,option)
      endif
      call OptionSetBlocking(option,PETSC_TRUE)
      call OptionCheckNonBlockingError(option)
    endif
    cur_thermal_cc => cur_thermal_cc%next
  enddo

end subroutine CharCurvesThermalConvertListToArray

! ************************************************************************** !

function CharacteristicCurvesThermalGetID(cc_thermal_array, &
     cc_thermal_name, &
     material_property_name, option)

  use Option_module
  use String_module

  type(cc_thermal_ptr_type), pointer :: cc_thermal_array(:)
  character(len=MAXWORDLENGTH) :: cc_thermal_name
  character(len=MAXWORDLENGTH) :: material_property_name
  type(option_type) :: option

  PetscInt :: CharacteristicCurvesGetID

  CharacteristicCurvesGetID = 0
  do CharacteristicCurvesGetID = 1, size(cc_thermal_array)
    if (StringCompare(cc_thermal_name, &
         cc_thermal_array(CharacteristicCurvesGetID)%ptr%name)) then
      return
    endif
  enddo
  option%io_buffer = 'Thermal characteristic curves "' // &
       trim(cc_thermal_name) // &
       '" in material property "' // &
       trim(material_property_name) // &
       '" not found among available thermal characteristic curves.'
  call PrintErrMsg(option)

end function CharacteristicCurvesThermalGetID

! ************************************************************************** !

subroutine CharCurvesThermalInputRecord(cc_thermal_list)

  implicit none

  class(cc_thermal_type), pointer :: cc_thermal_list

  class(cc_thermal_type), pointer :: cur_thermal_ccurve
  character(len=MAXWORDLENGTH) :: word1
  PetscInt :: id = INPUT_RECORD_UNIT

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
       &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'THERMAL CHARACTERISTIC CURVES'

  cur_thermal_ccurve => cc_thermal_list
  do
    if (.not.associated(cur_thermal_ccurve)) exit

    write(id,'(a29)',advance='no') 'thermal characteristic curve name: '
    write(id,'(a)') adjustl(trim(cur_thermal_ccurve%name))

    if (associated(cur_thermal_ccurve%thermal_conductivity_function)) then
      write(id,'(a29)',advance='no') 'thermal conductivity function: '
      select type (tcf => cur_thermal_ccurve%thermal_conductivity_function)
        !---------------------------------
      class is (kT_default_type)
        write(id,'(a)') 'default'
        write(id,'(a29)',advance='no') 'kT_wet: '
        write(word1,*) tcf%kT_wet
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'kT_dry: '
        write(word1,*) tcf%kT_dry
        write(id,'(a)') adjustl(trim(word1))
        !---------------------------------
      class is (kT_constant_type)
        write(id,'(a)') 'constant'
        write(id,'(a29)',advance='no') 'kT: '
        write(word1,*) tcf%constant_thermal_conductivity
        write(id,'(a)') adjustl(trim(word1))
        !---------------------------------
      class is (kT_power_type)
        write(id,'(a)') 'temperature-dependent (power)'
        write(id,'(a29)',advance='no') 'kT_wet: '
        write(word1,*) tcf%kT_wet
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'kT_dry: '
        write(word1,*) tcf%kT_dry
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'reference temp.: '
        write(word1,*) tcf%ref_temp
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'exponent: '
        write(word1,*) tcf%gamma
        write(id,'(a)') adjustl(trim(word1))
        !---------------------------------
      class is (kT_cubic_polynomial_type)
        write(id,'(a)') 'temperature-dependent (cubic polynomial)'
        write(id,'(a29)',advance='no') 'kT_wet: '
        write(word1,*) tcf%kT_wet
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'kT_dry: '
        write(word1,*) tcf%kT_dry
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'reference temp.: '
        write(word1,*) tcf%ref_temp
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'beta(1): '
        write(word1,*) tcf%beta(1)
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'beta(2): '
        write(word1,*) tcf%beta(2)
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'beta(3): '
        write(word1,*) tcf%beta(3)
        write(id,'(a)') adjustl(trim(word1))
      end select
    endif

    write(id,'(a29)') '---------------------------: '
    cur_thermal_ccurve => cur_thermal_ccurve%next
  enddo

end subroutine CharCurvesThermalInputRecord


! ************************************************************************** !

recursive subroutine ThermalCharacteristicCurvesDestroy(tcc)

  implicit none

  class(cc_thermal_type), pointer :: tcc

  if (.not.associated(tcc)) return

  call CharacteristicCurvesDestroy(ttcc%next)

  call ThermalConductivityFunctionDestroy(tcc%thermal_conductivity_function)

  deallocate(tcc)
  nullify(tcc)

end subroutine ThermalCharacteristicCurvesDestroy
end module Characteristic_Curves_Thermal_module
