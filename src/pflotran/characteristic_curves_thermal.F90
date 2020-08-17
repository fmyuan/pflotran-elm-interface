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
  contains
    procedure, public :: Verify => TCFBaseVerify
    procedure, public :: Test => TCFBaseTest
    procedure, public :: CalculateTCond => TCFBaseConductivity
    procedure, public :: TensorOp => ThermalConductivityTensorToScalar
  end type thermal_conductivity_base_type
  !---------------------------------------------------------------------------
  type, public, extends(thermal_conductivity_base_type) :: kT_constant_type
    PetscReal :: constant_thermal_conductivity
  contains
    procedure, public :: Verify => TCFConstantVerify
    procedure, public :: CalculateTCond => TCFConstantConductivity
  end type kT_constant_type
  !---------------------------------------------------------------------------
  type, public, extends(thermal_conductivity_base_type) :: kT_default_type
    PetscReal :: kT_wet, kT_dry
    PetscReal :: kT_X, kT_Y, kT_Z, kT_XY, kT_XZ, kT_YZ
    PetscReal, dimension(2,3,3) :: kT   ! thermal conductivity tensor
    PetscReal, dimension(3,3)   :: kTf  ! anisotropy ratio tensor    
    PetscBool :: isotropic
    PetscBool :: full_tensor
  contains
    procedure, public :: Verify => TCFDefaultVerify
    procedure, public :: CalculateTCond => TCFDefaultConductivity
  end type kT_default_type
  !---------------------------------------------------------------------------
  type, public, extends(kT_default_type) :: kT_power_type
    PetscReal :: ref_temp
    PetscReal :: gamma  ! T^gamma
  contains
    procedure, public :: Verify => TCFPowerVerify
    procedure, public :: CalculateTCond => TCFPowerConductivity
  end type kT_power_type
  !---------------------------------------------------------------------------
  type, public, extends(kT_default_type) :: kT_cubic_polynomial_type
    PetscReal :: ref_temp
    PetscReal :: beta(3) ! (T,T^2,T^3)
  contains
    procedure, public :: Verify => TCFCubicPolyVerify
    procedure, public :: CalculateTCond => TCFCubicPolyConductivity
  end type kT_cubic_polynomial_type
  !---------------------------------------------------------------------------
  type, public, extends(kT_default_type) :: kT_linear_resistivity_type
    PetscReal :: ref_temp
    PetscReal :: a(2)  ! 1/(a(1) + a(2)*T)
  contains
    procedure, public :: Verify => TCFLinearResistivityVerify
    procedure, public :: CalculateTCond => TCFLinearResistivityConductivity
  end type kT_linear_resistivity_type
  !---------------------------------------------------------------------------
  type, public :: cc_thermal_type
    character(len=MAXWORDLENGTH) :: name
    PetscBool :: print_me
    PetscBool :: test
    class(thermal_conductivity_base_type), pointer :: &
         thermal_conductivity_function
    class(cc_thermal_type), pointer :: next
  end type cc_thermal_type
  !---------------------------------------------------------------------------
  type, public :: cc_thermal_ptr_type
    class(cc_thermal_type), pointer :: ptr
  end type cc_thermal_ptr_type

  public :: CharCurvesThermalCreate, &
            CharCurvesThermalGetID, &
            CharCurvesThermalRead, &
            CharCurvesThermalAddToList, &
            CharCurvesThermalConvertListToArray, &
            CharCurvesThermalInputRecord, &
            CharCurvesThermalDestroy, &
            TCFDefaultCreate, &
            TCFConstantCreate, &
            TCFPowerCreate, &
            TCFCubicPolynomialCreate, &
            TCFLinearResistivityCreate, &
            TCFAssignDefault

contains

! ************************************************************************** !

subroutine TCFBaseVerify(this,name,option)

  use Option_module

  implicit none

  class(thermal_conductivity_base_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

end subroutine TCFBaseVerify

! ************************************************************************** !

subroutine TCFBaseConductivity(this,liquid_saturation,temperature, &
     thermal_conductivity,dkT_dsatl,dkT_dtemp,option)

  use Option_module

  implicit none

  class(thermal_conductivity_base_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dkT_dsatl, dkT_dtemp
  type(option_type), intent(inout) :: option

  thermal_conductivity = 0.d0
  dkT_dsatl = 0.d0
  dkT_dtemp = 0.d0

  option%io_buffer = 'TCFBaseThermalConductivity must be extended.'
  call PrintErrMsg(option)

end subroutine TCFBaseConductivity

! ************************************************************************** !

subroutine TCFBaseTest(this,tcc_name,option)

  use Option_module

  implicit none

  class(thermal_conductivity_base_type) :: this
  character(len=MAXWORDLENGTH) :: tcc_name
  type(option_type), intent(inout) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt, parameter :: nt = 71
  PetscInt, parameter :: ns = 31
  PetscReal, parameter :: perturbation = 1.0D-6
  PetscReal :: deltaTemp, deltaSat
  PetscReal :: temp_vec(nt)
  PetscReal :: sat_vec(ns)
  PetscReal :: kT(nt,ns)
  PetscReal :: dkT_dsat(nt,ns)
  PetscReal :: dkT_dsat_numerical(nt,ns)
  PetscReal :: dkT_dtemp(nt,ns)
  PetscReal :: dkT_dtemp_numerical(nt,ns)
  PetscReal :: perturbed_temp, perturbed_sat
  PetscReal :: kT_temp_pert, kT_sat_pert, unused1, unused2
  PetscReal :: temp_min, temp_max, sat_min, sat_max
  PetscInt :: i,j

  ! thermal conductivity as a function of temp. and liq. sat.
  temp_min = 1.0d0 ! Celsius
  temp_max = 250.0d0
  sat_min = 1.0d-3
  sat_max = 1.0d0

  deltaTemp = (temp_max - temp_min)/(nt - 1)
  deltaSat = (sat_max - sat_min)/(ns - 1)

  temp_vec = [(temp_min + i*deltaTemp, i=0,nt-1)]
  sat_vec = [(sat_min + i*deltaSat, i=0,ns-1)]

  do i = 1,nt
    do j = 1,ns
      ! base case with analytical derivatives
      call this%CalculateTCond(sat_vec(j),temp_vec(i), &
           kT(i,j),dkT_dsat(i,j),dkT_dtemp(i,j),option)

      ! calculate numerical derivatives via finite differences
      perturbed_temp = temp_vec(i) * (1.d0 + perturbation)
      call this%CalculateTCond(sat_vec(j),perturbed_temp, &
           kT_temp_pert,unused1,unused2,option)

      dkT_dtemp_numerical(i,j) = (kT_temp_pert - kT(i,j))/(temp_vec(i)*perturbation)

      perturbed_sat = sat_vec(j) * (1.d0 + perturbation)
      call this%CalculateTCond(perturbed_sat,temp_vec(i), &
           kT_sat_pert,unused1,unused2,option)

      dkT_dsat_numerical(i,j) = (kT_sat_pert - kT(i,j))/(sat_vec(j)*perturbation)
    enddo
  enddo

  write(string,*) tcc_name
  string = trim(tcc_name) // '_kT_vs_sat_and_temp.dat'
  open(unit=86,file=string)
  write(86,*) '"temperature [C]", "liquid saturation [-]", "kT [W/m*K]", &
       &"dkT/dsat", "dkT/dT", "dkT/dsat_numerical", "dkT/dT_numerical"'
  do i = 1,nt
    do j = 1,ns
      write(86,'(7(ES14.6))') temp_vec(i), sat_vec(j), &
           kT(i,j), dkT_dsat(i,j), dkT_dtemp(i,j), &
           dkT_dsat_numerical(i,j), dkT_dtemp_numerical(i,j)
    enddo
  enddo
  close(86)

end subroutine TCFBaseTest

! ************************************************************************** !

subroutine TCFDestroy(tcf)

  implicit none

  class(thermal_conductivity_base_type), pointer :: tcf

  if (.not.associated(tcf)) return
  deallocate(tcf)
  nullify(tcf)

end subroutine TCFDestroy

! ************************************************************************** !

function TCFDefaultCreate()

  implicit none

  class(kT_default_type), pointer :: TCFDefaultCreate

  allocate(TCFDefaultCreate)
  TCFDefaultCreate%isotropic  = PETSC_TRUE
  TCFDefaultCreate%full_tensor= PETSC_FALSE
  TCFDefaultCreate%kT_wet = UNINITIALIZED_DOUBLE
  TCFDefaultCreate%kT_dry = UNINITIALIZED_DOUBLE
  TCFDefaultCreate%kT_X   = UNINITIALIZED_DOUBLE
  TCFDefaultCreate%kT_Y   = UNINITIALIZED_DOUBLE
  TCFDefaultCreate%kT_Z   = UNINITIALIZED_DOUBLE
  TCFDefaultCreate%kT_XY  = UNINITIALIZED_DOUBLE
  TCFDefaultCreate%kT_XZ  = UNINITIALIZED_DOUBLE
  TCFDefaultCreate%kT_YZ  = UNINITIALIZED_DOUBLE

end function TCFDefaultCreate

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
     thermal_conductivity,dkT_dsatl,dkT_dtemp,option)

  use Option_module

  implicit none

  class(kT_default_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dkT_dsatl, dkT_dtemp
  type(option_type), intent(inout) :: option

  PetscReal :: tempreal

  ! based on Somerton et al., 1974:
  ! k_eff = k_dry + sqrt(s_l)*(k_wet-k_dry)

  dkT_dtemp = 0.d0 ! only a function of saturation
  
  if (liquid_saturation > 0.d0) then
    tempreal = sqrt(liquid_saturation) * &
                   (this%kT_wet - this%kT_dry)
    thermal_conductivity = this%kT_dry + tempreal
    dkT_dsatl = 0.5d0 * tempreal / liquid_saturation
  else
    thermal_conductivity = this%kT_dry
    dkT_dsatl = 0.d0
  endif

end subroutine TCFDefaultConductivity

! ************************************************************************** !

function TCFConstantCreate()

  implicit none

  class(kT_constant_type), pointer :: TCFConstantCreate

  allocate(TCFConstantCreate)
  TCFConstantCreate%constant_thermal_conductivity = UNINITIALIZED_DOUBLE

end function TCFConstantCreate

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

subroutine TCFConstantConductivity(this,liquid_saturation,temperature, &
     thermal_conductivity,dkT_dsatl,dkT_dtemp,option)

  use Option_module

  implicit none

  class(kT_constant_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dkT_dsatl, dkT_dtemp
  type(option_type), intent(inout) :: option

  thermal_conductivity = this%constant_thermal_conductivity
  dkT_dsatl = 0.d0
  dkT_dtemp = 0.d0

end subroutine TCFConstantConductivity

! ************************************************************************** !

function TCFPowerCreate()

  implicit none

  class(kT_power_type), pointer :: TCFPowerCreate

  allocate(TCFPowerCreate)
  TCFPowerCreate%kT_wet = UNINITIALIZED_DOUBLE
  TCFPowerCreate%kT_dry = UNINITIALIZED_DOUBLE
  TCFPowerCreate%ref_temp = -273.15d0
  TCFPowerCreate%gamma = UNINITIALIZED_DOUBLE
  TCFPowerCreate%isotropic   = PETSC_TRUE
  TCFPowerCreate%full_tensor = PETSC_FALSE
  TCFPowerCreate%kT_X   = UNINITIALIZED_DOUBLE
  TCFPowerCreate%kT_Y   = UNINITIALIZED_DOUBLE
  TCFPowerCreate%kT_Z   = UNINITIALIZED_DOUBLE
  TCFPowerCreate%kT_XY  = UNINITIALIZED_DOUBLE
  TCFPowerCreate%kT_XZ  = UNINITIALIZED_DOUBLE
  TCFPowerCreate%kT_YZ  = UNINITIALIZED_DOUBLE

end function TCFPowerCreate

! ************************************************************************** !

subroutine TCFPowerVerify(this,name,option)

  use Option_module

  implicit none

  class(kT_power_type) :: this
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
  if (Uninitialized(this%gamma)) then
    option%io_buffer = UninitializedMessage('EXPONENT',string)
    call PrintErrMsg(option)
  endif

end subroutine TCFPowerVerify

! ************************************************************************** !

subroutine TCFPowerConductivity(this,liquid_saturation,temperature, &
     thermal_conductivity,dkT_dsatl,dkT_dtemp,option)

  use Option_module

  implicit none

  class(kT_power_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dkT_dsatl, dkT_dtemp
  type(option_type), intent(inout) :: option

  PetscReal :: tempreal, kT_base, shifted_temp, unused

  ! saturation behavior from default function
  call this%kT_default_type%CalculateTCond(liquid_saturation,0.d0, &
       kT_base,dkT_dsatl,unused,option)

  shifted_temp = temperature - this%ref_temp

  ! kT = kT_0 * (T/300) ** gamma [T in Kelvin]
  tempreal = (shifted_temp / 3.0D+2) ** this%gamma
  thermal_conductivity = kT_base * tempreal
  dkT_dtemp = this%gamma * kT_base * tempreal / shifted_temp
  dkT_dsatl = dkT_dsatl * tempreal

end subroutine TCFPowerConductivity

! ************************************************************************** !

function TCFCubicPolynomialCreate()

  implicit none

  class(kT_cubic_polynomial_type), pointer :: TCFCubicPolynomialCreate

  allocate(TCFCubicPolynomialCreate)
  TCFCubicPolynomialCreate%kT_wet = UNINITIALIZED_DOUBLE
  TCFCubicPolynomialCreate%kT_dry = UNINITIALIZED_DOUBLE
  TCFCubicPolynomialCreate%ref_temp = 0.d0
  TCFCubicPolynomialCreate%beta = [ UNINITIALIZED_DOUBLE, &
                                       UNINITIALIZED_DOUBLE, &
                                       UNINITIALIZED_DOUBLE ]
  TCFCubicPolynomialCreate%isotropic  = PETSC_TRUE
  TCFCubicPolynomialCreate%full_tensor= PETSC_FALSE
  TCFCubicPolynomialCreate%kT_X   = UNINITIALIZED_DOUBLE
  TCFCubicPolynomialCreate%kT_Y   = UNINITIALIZED_DOUBLE
  TCFCubicPolynomialCreate%kT_Z   = UNINITIALIZED_DOUBLE
  TCFCubicPolynomialCreate%kT_XY  = UNINITIALIZED_DOUBLE
  TCFCubicPolynomialCreate%kT_XZ  = UNINITIALIZED_DOUBLE
  TCFCubicPolynomialCreate%kT_YZ  = UNINITIALIZED_DOUBLE

end function TCFCubicPolynomialCreate

! ************************************************************************** !

subroutine TCFCubicPolyVerify(this,name,option)

  use Option_module

  implicit none

  class(kT_cubic_polynomial_type) :: this
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
  if (Uninitialized(this%beta(1)).or. &
      Uninitialized(this%beta(2)) .or. &
      Uninitialized(this%beta(3))) then
    option%io_buffer = UninitializedMessage( &
         'CUBIC_POLYNOMIAL_COEFFCIENTS (1+b(1)*T+b(2)*T^2+b(3)*T^3)',string)
    call PrintErrMsg(option)
  endif

end subroutine TCFCubicPolyVerify

! ************************************************************************** !

subroutine TCFCubicPolyConductivity(this,liquid_saturation,temperature, &
     thermal_conductivity,dkT_dsatl,dkT_dtemp,option)

  use Option_module

  implicit none

  class(kT_cubic_polynomial_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dkT_dsatl, dkT_dtemp
  type(option_type), intent(inout) :: option

  PetscReal :: kT_base, shifted_temp, unused, tempreal

  ! saturation behavior from default function
  call this%kT_default_type%CalculateTCond(liquid_saturation,0.d0, &
       kT_base,dkT_dsatl,unused,option)

  shifted_temp = temperature - this%ref_temp

  ! kT = k0 * (1 + beta1*T + beta2*T^2 + beta3*T^3)
  tempreal = (1.d0 + shifted_temp*(this%beta(1) &
       + shifted_temp*(this%beta(2) &
       + shifted_temp* this%beta(3))))
  thermal_conductivity = kT_base * tempreal

  dkT_dtemp = kT_base * (this%beta(1) &
       + shifted_temp*(2.d0*this%beta(2) &
       + shifted_temp* 3.d0*this%beta(3)))

  dkT_dsatl = dkT_dsatl * tempreal

end subroutine TCFCubicPolyConductivity

! ************************************************************************** !

function TCFLinearResistivityCreate()

  implicit none

  class(kT_linear_resistivity_type), pointer :: TCFLinearResistivityCreate

  allocate(TCFLinearResistivityCreate)
  TCFLinearResistivityCreate%kT_wet = UNINITIALIZED_DOUBLE
  TCFLinearResistivityCreate%kT_dry = UNINITIALIZED_DOUBLE
  TCFLinearResistivityCreate%ref_temp = 0.d0
  TCFLinearResistivityCreate%a = [ UNINITIALIZED_DOUBLE, &
                                      UNINITIALIZED_DOUBLE]
  TCFLinearResistivityCreate%isotropic  = PETSC_TRUE
  TCFLinearResistivityCreate%full_tensor= PETSC_FALSE
  TCFLinearResistivityCreate%kT_X   = UNINITIALIZED_DOUBLE
  TCFLinearResistivityCreate%kT_Y   = UNINITIALIZED_DOUBLE
  TCFLinearResistivityCreate%kT_Z   = UNINITIALIZED_DOUBLE
  TCFLinearResistivityCreate%kT_XY  = UNINITIALIZED_DOUBLE
  TCFLinearResistivityCreate%kT_XZ  = UNINITIALIZED_DOUBLE
  TCFLinearResistivityCreate%kT_YZ  = UNINITIALIZED_DOUBLE

end function TCFLinearResistivityCreate

! ************************************************************************** !

subroutine TCFLinearResistivityVerify(this,name,option)

  use Option_module

  implicit none

  class(kT_linear_resistivity_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'THERMAL_CONDUCTIVITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'THERMAL_CONDUCTIVITY_FUNCTION,LINEAR_RESISTIVITY'
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
  if (Uninitialized(this%a(1)) .or. Uninitialized(this%a(2))) then
     option%io_buffer = UninitializedMessage( &
          'RESISTIVITY_COEFFCIENTS 1/(a(1) + a(2)*T)',string)
     call PrintErrMsg(option)
  endif

end subroutine TCFLinearResistivityVerify

! ************************************************************************** !

subroutine TCFLinearResistivityConductivity(this,liquid_saturation,temperature, &
     thermal_conductivity,dkT_dsatl,dkT_dtemp,option)

  use Option_module

  implicit none

  class(kT_linear_resistivity_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dkT_dsatl, dkT_dtemp
  type(option_type), intent(inout) :: option

  PetscReal :: relkT, shifted_temp, tempreal, unused

  ! saturation behavior from default function
  call this%kT_default_type%CalculateTCond(liquid_saturation,0.d0, &
       relkT,dkT_dsatl,unused,option)

  shifted_temp = temperature - this%ref_temp

  ! linear thermal resistivity idea from Birch & Clark (1940) and
  ! used by Blesch, Kulacki & Christensen (1983), ONWI-495
  ! kT = relkT/(a1 + a2*T)
  tempreal = this%a(1) + shifted_temp*this%a(2)
  thermal_conductivity = relkT / tempreal

  dkT_dtemp = -this%a(2) * thermal_conductivity / tempreal
  dkT_dsatl = dkT_dsatl / tempreal

end subroutine TCFLinearResistivityConductivity

! ************************************************************************** !

function CharCurvesThermalCreate()

  implicit none

  class(cc_thermal_type), pointer :: CharCurvesThermalCreate
  class(cc_thermal_type), pointer :: characteristic_curves_thermal

  allocate(characteristic_curves_thermal)
  characteristic_curves_thermal%name = ''
  characteristic_curves_thermal%print_me = PETSC_FALSE
  characteristic_curves_thermal%test = PETSC_FALSE
  nullify(characteristic_curves_thermal%thermal_conductivity_function)
  nullify(characteristic_curves_thermal%next)

  CharCurvesThermalCreate => characteristic_curves_thermal

end function CharCurvesThermalCreate

! ************************************************************************** !

subroutine CharCurvesThermalRead(this,input,option)

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
      call InputReadCard(input,option,word)
      call InputErrorMsg(input,option, &
           'THERMAL_CONDUTIVITY_FUNCTION',error_string)
      call StringToUpper(word)
      select case(word)
      case('CONSTANT')
        this%thermal_conductivity_function => TCFConstantCreate()
      case('DEFAULT')
        this%thermal_conductivity_function => TCFDefaultCreate()
      case('POWER')
        this%thermal_conductivity_function => TCFPowerCreate()
      case('CUBIC_POLYNOMIAL')
        this%thermal_conductivity_function => TCFCubicPolynomialCreate()
      case('LINEAR_RESISTIVITY')
        this%thermal_conductivity_function => TCFLinearResistivityCreate()
      case default
        call InputKeywordUnrecognized(input,word, &
             'THERMAL_CONDUCTIVITY_FUNCTION',option)
      end select
      call TCFRead( &
           this%thermal_conductivity_function,input,option)
      call TCFCheckAnisotropy( &
           this%thermal_conductivity_function,option)
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

end subroutine CharCurvesThermalRead

! ************************************************************************** !

subroutine TCFAssignDefault(thermal_conductivity_function,&
                                                 kwet,kdry,option)

  use Option_module

  implicit none

  class(thermal_conductivity_base_type) :: thermal_conductivity_function
  PetscReal :: kwet,kdry
  type(option_type) :: option

  select type(tcf => thermal_conductivity_function)
      !------------------------------------------
    class is(kT_default_type)
      tcf%kT_dry = kdry
      tcf%kT_wet = kwet
  end select

end subroutine TCFAssignDefault

! ************************************************************************** !

subroutine TCFRead(thermal_conductivity_function,input,option)
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

  input%ierr = 0
  error_string = 'THERMAL_CHARACTERISTIC_CURVES,&
       & THERMAL_CONDUCTIVITY_FUNCTION,'
  select type(tcf => thermal_conductivity_function)
  class is(kT_constant_type)
    error_string = trim(error_string) // 'CONSTANT'
  class is(kT_default_type)
    error_string = trim(error_string) // 'DEFAULT'
  class is(kT_power_type)
    error_string = trim(error_string) // 'POWER'
  class is(kT_cubic_polynomial_type)
    error_string = trim(error_string) // 'CUBIC_POLYNOMIAL'
  class is(kT_linear_resistivity_type)
    error_string = trim(error_string) // 'LINEAR_RESISTIVITY'
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
        call InputReadAndConvertUnits(input, &
             tcf%constant_thermal_conductivity,'W/m-C', &
             'CHARACTERISTIC_CURVES_THERMAL,constant thermal conductivity', &
             option)
      case default
        call InputKeywordUnrecognized(input,keyword, &
             'constant thermal conductivity',option)
      end select
      !------------------------------------------
    class is(kT_default_type)
      call TCFDefaultRead(tcf,input,keyword,error_string,'default',option)
      !------------------------------------------
    class is(kT_power_type)
      select case(keyword)
      case('REFERENCE_TEMPERATURE')
        call InputReadDouble(input,option,tcf%ref_temp)
        call InputErrorMsg(input,option,'reference temperature', &
             error_string)
        call InputReadAndConvertUnits(input,tcf%ref_temp, &
             'C','CHARACTERISTIC_CURVES_THERMAL,reference temperature',option)
      case('EXPONENT')
        call InputReadDouble(input,option,tcf%gamma)
        call InputErrorMsg(input,option,'thermal conductivity exponent', &
             error_string)
      case default
        call TCFDefaultRead(tcf,input,keyword,error_string,'power',option)
      end select
      !------------------------------------------
    class is(kT_cubic_polynomial_type)
      select case(keyword)
      case('REFERENCE_TEMPERATURE')
        call InputReadDouble(input,option,tcf%ref_temp)
        call InputErrorMsg(input,option,'reference temperature', &
             error_string)
        call InputReadAndConvertUnits(input,tcf%ref_temp, &
             'C','CHARACTERISTIC_CURVES_THERMAL,reference temperature',option)
      case('CUBIC_POLYNOMIAL_COEFFICIENTS')
        call InputReadNDoubles(input,option,tcf%beta,3)
        call InputErrorMsg(input,option, &
             'thermal conductivity polynomial coefficients',error_string)
      case default
        call TCFDefaultRead(tcf,input,keyword,error_string, &
          'cubic polynomial',option)
      end select
      !------------------------------------------
    class is(kT_linear_resistivity_type)
      select case(keyword)      
      case('REFERENCE_TEMPERATURE')
        call InputReadDouble(input,option,tcf%ref_temp)
        call InputErrorMsg(input,option,'reference temperature', &
             error_string)
        call InputReadAndConvertUnits(input,tcf%ref_temp, &
             'C','CHARACTERISTIC_CURVES_THERMAL,reference temperature',option)
      case('LINEAR_RESISTIVITY_COEFFICIENTS')
        call InputReadNDoubles(input,option,tcf%a,2)
        call InputErrorMsg(input,option, &
             'linear thermal resistivity coefficients',error_string)
      case default
        call TCFDefaultRead(tcf,input,keyword,error_string, &
          'linear resistivity',option)
      end select

    class default
      option%io_buffer = 'Read routine not implemented for ' &
           // trim(error_string) // '.'
      call PrintErrMsg(option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine TCFRead

! ************************************************************************** !

subroutine TCFDefaultRead(tcf,input,keyword,error_string,kind,option)
  !
  ! Reads in contents of THERMAL_CONDUCTIVITY_FUNCTION block for dervived 
  ! types of thermal conductivity default base class
  !
  use Option_module
  use Input_Aux_module
  use String_module
  
  class(kT_default_type) :: tcf
  type(input_type) :: input
  character(len=MAXWORDLENGTH)   :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=*)  :: kind
  type(option_type) :: option
  
  select case(keyword)
  case('THERMAL_CONDUCTIVITY_WET')
    call InputReadDouble(input,option,tcf%kT_wet)
    call InputErrorMsg(input,option,'thermal conductivity wet', &
         error_string)
    call InputReadAndConvertUnits(input,tcf%kT_wet,'W/m-C', &
         'CHARACTERISTIC_CURVES_THERMAL,thermal conductivity wet',option)
  case('THERMAL_CONDUCTIVITY_DRY')
    call InputReadDouble(input,option,tcf%kT_dry)
    call InputErrorMsg(input,option,'thermal conductivity dry', &
         error_string)
    call InputReadAndConvertUnits(input,tcf%kT_dry,'W/m-C', &
         'CHARACTERISTIC_CURVES_THERMAL,thermal conductivity dry',option)
  case('THERMAL_CONDUCTIVITY_X')
    call InputReadDouble(input,option,tcf%kT_X)
    call InputErrorMsg(input,option, & 
       'anisotropic thermal conductivity X component', error_string)
  case('THERMAL_CONDUCTIVITY_Y')
    call InputReadDouble(input,option,tcf%kT_Y)
    call InputErrorMsg(input,option, & 
       'anisotropic thermal conductivity Y component', error_string)
  case('THERMAL_CONDUCTIVITY_Z')
    call InputReadDouble(input,option,tcf%kT_Z)
    call InputErrorMsg(input,option, & 
       'anisotropic thermal conductivity Z component', error_string)
  case('THERMAL_CONDUCTIVITY_XY')
    call InputReadDouble(input,option,tcf%kT_XY)
    call InputErrorMsg(input,option, & 
       'anisotropic thermal conductivity XY component', error_string)
  case('THERMAL_CONDUCTIVITY_XZ')
    call InputReadDouble(input,option,tcf%kT_XZ)
    call InputErrorMsg(input,option, & 
       'anisotropic thermal conductivity XZ component', error_string)
  case('THERMAL_CONDUCTIVITY_YZ')
    call InputReadDouble(input,option,tcf%kT_YZ)
    call InputErrorMsg(input,option, & 
       'anisotropic thermal conductivity YZ component', error_string)
  case default
     call InputKeywordUnrecognized(input,keyword, &
          'temp-dependent ('//trim(kind)//') thermal conductivity',option)
  end select
  
end subroutine

! ************************************************************************** !

subroutine TCFCheckAnisotropy(thermal_conductivity_function,option)
  !
  ! Check anisotropy parameters
  !
  use Option_module

  implicit none

  class(thermal_conductivity_base_type) :: thermal_conductivity_function
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'THERMAL_CHARACTERISTIC_CURVES,&
       & THERMAL_CONDUCTIVITY_FUNCTION,'
  
  select type(tcf => thermal_conductivity_function)
  class is(kT_default_type)
    error_string = trim(error_string) // ' DEFAULT'
    
    ! check if wet and dry thermal conductivities are initialized
    if (.not. Initialized(tcf%kT_dry) .or. &
        .not. Initialized(tcf%kT_wet)) then
      ! wet and dry values must be specified per anisotropic component
      option%io_buffer = 'Must specify wet and dry thermal conductivity ' &
        //'values in order to use anisotropy ratios in ' &
        // trim(error_string) // '.'
      call PrintErrMsg(option)
    else if (Initialized(tcf%kT_X) .or. Initialized(tcf%kT_Y) &
      .or.   Initialized(tcf%kT_Z)) then
      ! inputs must be anisotropy ratios between zero and one
      
      ! check diagonal components first, as tensor must at least be diagonal
      if (Initialized(tcf%kT_X)) then
        if (tcf%kT_X < 0.0d0 .or. tcf%kT_X > 1.0d0) then
          option%io_buffer = 'Anisotropy ratio for X must lie between 0 and ' &
               //'1 in '// trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        tcf%kTf(1,1) = tcf%kT_X
        tcf%kTf(1,1) = tcf%kT_X
        tcf%kT(1,1,1) = tcf%kT_dry * tcf%kT_X
        tcf%kT(2,1,1) = tcf%kT_wet * tcf%kT_X
      else
        option%io_buffer = 'Anisotropy ratio for X uninitialized in ' &
             // trim(error_string) // '.'
        call PrintErrMsg(option)
      endif
      
      if (Initialized(tcf%kT_Y)) then
        if (tcf%kT_Y < 0.0d0 .or. tcf%kT_Y > 1.0d0) then
          option%io_buffer = 'Anisotropy ratio for Y must lie between 0 and ' &
               //'1 in '// trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        tcf%kTf(2,2) = tcf%kT_Y
        tcf%kTf(2,2) = tcf%kT_Y
        tcf%kT(1,2,2) = tcf%kT_dry * tcf%kT_Y
        tcf%kT(2,2,2) = tcf%kT_wet * tcf%kT_Y
      else
        option%io_buffer = 'Anisotropy ratio for Y uninitialized in ' &
             // trim(error_string) // '.'
        call PrintErrMsg(option)
      endif
      
      if (Initialized(tcf%kT_Z)) then
        if (tcf%kT_Z < 0.0d0 .or. tcf%kT_Z > 1.0d0) then
          option%io_buffer = 'Anisotropy ratio for Z must lie between 0 and ' &
               //'1 in '// trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        tcf%kTf(3,3) = tcf%kT_Z
        tcf%kTf(3,3) = tcf%kT_Z
        tcf%kT(1,3,3) = tcf%kT_dry * tcf%kT_Z
        tcf%kT(2,3,3) = tcf%kT_wet * tcf%kT_Z
      else
        option%io_buffer = 'Anisotropy ratio for Z uninitialized in ' &
             // trim(error_string) // '.'
        call PrintErrMsg(option)
      endif
      
      ! check off-diagonal components next; if one is given, so must the others
      if (Initialized(tcf%kT_XY)) then
        if (tcf%kT_XY < 0.0d0 .or. tcf%kT_XY > 1.0d0) then
          option%io_buffer = 'Anisotropy ratio for XY must lie between 0 and ' &
               //'1 in '// trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        if (.not. Initialized(tcf%kT_XZ) .or. & 
          .not. Initialized(tcf%kT_YZ)) then
          option%io_buffer = 'All off-diagonal components must be specified' &
               // ' if XY ratio is provided in '&
               // trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        tcf%isotropic = PETSC_FALSE
        tcf%kTf(1,2) = tcf%kT_XY
        tcf%kTf(1,2) = tcf%kT_XY
        tcf%kTf(2,1) = tcf%kT_XY
        tcf%kTf(2,1) = tcf%kT_XY
        tcf%kT(1,1,2) = tcf%kT_dry * tcf%kT_XY
        tcf%kT(2,1,2) = tcf%kT_wet * tcf%kT_XY
        tcf%kT(1,2,1) = tcf%kT_dry * tcf%kT_XY
        tcf%kT(2,2,1) = tcf%kT_wet * tcf%kT_XY
      endif
      
      if (Initialized(tcf%kT_XZ)) then
        if (tcf%kT_XZ < 0.0d0 .or. tcf%kT_XZ > 1.0d0) then
          option%io_buffer = 'Anisotropy ratio for XZ must lie between 0 and ' &
               //'1 in '// trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        if (.not. Initialized(tcf%kT_XY) .or. & 
          .not. Initialized(tcf%kT_YZ)) then
          option%io_buffer = 'All off-diagonal components must be specified' &
               // ' if XZ ratio is provided in '&
               // trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        tcf%isotropic = PETSC_FALSE
        tcf%kTf(1,3) = tcf%kT_XZ
        tcf%kTf(1,3) = tcf%kT_XZ
        tcf%kTf(3,1) = tcf%kT_XZ
        tcf%kTf(3,1) = tcf%kT_XZ
        tcf%kT(1,1,3) = tcf%kT_dry * tcf%kT_XZ
        tcf%kT(2,1,3) = tcf%kT_wet * tcf%kT_XZ
        tcf%kT(1,3,1) = tcf%kT_dry * tcf%kT_XZ
        tcf%kT(2,3,1) = tcf%kT_wet * tcf%kT_XZ
      endif
      
      if (Initialized(tcf%kT_YZ)) then
        if (tcf%kT_YZ < 0.0d0 .or. tcf%kT_YZ > 1.0d0) then
          option%io_buffer = 'Anisotropy ratio for YZ must lie between 0 and ' &
               //'1 in '// trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        if (.not. Initialized(tcf%kT_XY) .or. & 
          .not. Initialized(tcf%kT_XZ)) then
          option%io_buffer = 'All off-diagonal components must be specified' &
               // ' if YZ ratio is provided in '&
               // trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        tcf%isotropic = PETSC_FALSE
        tcf%kTf(2,3) = tcf%kT_YZ
        tcf%kTf(2,3) = tcf%kT_YZ
        tcf%kTf(3,2) = tcf%kT_YZ
        tcf%kTf(3,2) = tcf%kT_YZ
        tcf%kT(1,2,3) = tcf%kT_dry * tcf%kT_YZ
        tcf%kT(2,2,3) = tcf%kT_wet * tcf%kT_YZ
        tcf%kT(1,3,2) = tcf%kT_dry * tcf%kT_YZ
        tcf%kT(2,3,2) = tcf%kT_wet * tcf%kT_YZ
      endif
      
      ! check for isotropy and fully initialize tensor
      if (tcf%kT_X == tcf%kT_Y .and. tcf%kT_Y == tcf%kT_Z) then
        if (Initialized(tcf%kT_XY) .or. Initialized(tcf%kT_XZ) &
          .or. Initialized(tcf%kT_YZ)) then
          tcf%isotropic = PETSC_FALSE
          tcf%full_tensor = PETSC_TRUE
        else
          tcf%isotropic = PETSC_TRUE
          tcf%kT(:,:,:) = 0.0d0
          tcf%kT(1,1,1) = tcf%kT_dry
          tcf%kT(1,2,2) = tcf%kT_dry
          tcf%kT(1,3,3) = tcf%kT_dry
          tcf%kT(2,1,1) = tcf%kT_wet
          tcf%kT(2,2,2) = tcf%kT_wet
          tcf%kT(2,3,3) = tcf%kT_wet
          option%io_buffer = 'Thermal conductivity will be treated as' &
               // ' isotropic in '&
               // trim(error_string) // '.'
          call PrintMsg(option)
        endif
      else
        tcf%isotropic = PETSC_FALSE
        if (Initialized(tcf%kT_XY) &
          .or. Initialized(tcf%kT_XZ) &
          .or. Initialized(tcf%kT_YZ)) then
          ! full thermal conductivity tensor
          tcf%full_tensor = PETSC_TRUE
        else 
          ! diagonal thermal conductivity tensor
          tcf%kT(:,1,2) = 0.0d0
          tcf%kT(:,1,3) = 0.0d0
          tcf%kT(:,2,3) = 0.0d0
          tcf%kT(:,2,1) = 0.0d0
          tcf%kT(:,3,1) = 0.0d0
          tcf%kT(:,3,2) = 0.0d0
          tcf%kTf(1,2) = 0.0d0
          tcf%kTf(1,3) = 0.0d0
          tcf%kTf(2,3) = 0.0d0
          tcf%kTf(2,1) = 0.0d0
          tcf%kTf(3,1) = 0.0d0
          tcf%kTf(3,2) = 0.0d0
        endif
      endif
      
    else if (Initialized(tcf%kT_XY) &
      .or. Initialized(tcf%kT_XZ) &
      .or. Initialized(tcf%kT_YZ)) then
      if (.not. Initialized(tcf%kT_X) .or. & 
        .not. Initialized(tcf%kT_Y) .or. &
        .not. Initialized(tcf%kT_Z)) then
        option%io_buffer = 'Diagonal components of thermal conductivity ' &
             // 'must be specified if off-diagonal components are '&
             // 'provided in '&
             // trim(error_string) // '.'
        call PrintErrMsg(option)
      endif
    else
      tcf%isotropic = PETSC_TRUE
      tcf%kT(:,:,:) = 0.0d0
      tcf%kT(1,1,1) = tcf%kT_dry
      tcf%kT(1,2,2) = tcf%kT_dry
      tcf%kT(1,3,3) = tcf%kT_dry
      tcf%kT(2,1,1) = tcf%kT_wet
      tcf%kT(2,2,2) = tcf%kT_wet
      tcf%kT(2,3,3) = tcf%kT_wet
      tcf%kTf(:,:) = 0.0d0
      tcf%kTf(1,1) = 1.0d0
      tcf%kTf(2,2) = 1.0d0
      tcf%kTf(3,3) = 1.0d0
      tcf%kTf(1,1) = 1.0d0
      tcf%kTf(2,2) = 1.0d0
      tcf%kTf(3,3) = 1.0d0
    endif
    
    if (tcf%full_tensor) then
      call FullTensorCheckEigenvalues(tcf%kTf,option)
    endif
    
  end select

end subroutine TCFCheckAnisotropy

! ************************************************************************** !

subroutine ThermalConductivityTensorToScalar(this,dist,option)
  !
  ! Transform thermal conductivity tensor to a scalar via dot product
  !
  use Option_module

  implicit none

  class(thermal_conductivity_base_type) :: this
  type(option_type) :: option
  ! -1 = fraction upwind
  ! 0 = magnitude
  ! 1 = unit x-dir
  ! 2 = unit y-dir
  ! 3 = unit z-dir
  PetscReal, intent(in) :: dist(-1:3)

  PetscReal :: kxw, kyw, kzw, kxyw, kxzw, kyzw
  PetscReal :: kxd, kyd, kzd, kxyd, kxzd, kyzd
  
  select type(tcf => this)
  class is(kT_default_type)
    if (tcf%isotropic) then 
      return
    endif
    
    kxd = tcf%kT(1,1,1)
    kyd = tcf%kT(1,2,2)
    kzd = tcf%kT(1,3,3)
    kxw = tcf%kT(2,1,1)
    kyw = tcf%kT(2,2,2)
    kzw = tcf%kT(2,3,3)
    
    if (tcf%full_tensor) then
      kxyd = tcf%kT(1,1,2)
      kxzd = tcf%kT(1,1,3)
      kyzd = tcf%kT(1,2,3)
      kxyw = tcf%kT(2,1,2)
      kxzw = tcf%kT(2,1,3)
      kyzw = tcf%kT(2,2,3)
      tcf%kT_dry = &
        FullThermalConductivityTensorToScalar(kxd,kyd,kzd,kxyd,kxzd,kyzd,&
        dist,option) 
      tcf%kT_wet = &
        FullThermalConductivityTensorToScalar(kxw,kyw,kzw,kxyw,kxzw,kyzw,&
        dist,option) 
    else if (.not. tcf%isotropic) then
      tcf%kT_dry = &
        DiagThermalConductivityTensorToScalar(kxd,kyd,kzd,&
        dist,option) 
      tcf%kT_wet = &
        DiagThermalConductivityTensorToScalar(kxw,kyw,kzw,&
        dist,option) 
    endif
    
  end select
  
end subroutine ThermalConductivityTensorToScalar

! ************************************************************************** !

function FullThermalConductivityTensorToScalar(kx,ky,kz,kxy,kxz,kyz,dist,&
  option)
  !
  ! Full tensor directional thermal conductivity
  !
  use Option_module

  implicit none

  type(option_type) :: option
  ! -1 = fraction upwind
  ! 0 = magnitude
  ! 1 = unit x-dir
  ! 2 = unit y-dir
  ! 3 = unit z-dir
  PetscReal :: FullThermalConductivityTensorToScalar
  PetscReal, intent(in)  :: dist(-1:3)

  PetscReal :: kx,ky,kz,kxy,kxz,kyz
  
  FullThermalConductivityTensorToScalar = &
                                kx*dabs(dist(1))**2.0 + &
                                ky*dabs(dist(2))**2.0 + &
                                kz*dabs(dist(3))**2.0 + &
                                2*kxy*dist(1)*dist(2) + &
                                2*kxz*dist(1)*dist(3) + &
                                2*kyz*dist(2)*dist(3)
  
end function FullThermalConductivityTensorToScalar

! ************************************************************************** !

function DiagThermalConductivityTensorToScalar(kx,ky,kz,dist,&
  option)
  !
  ! Diagonal tensor directional thermal conductivity
  !
  use Option_module

  implicit none

  type(option_type) :: option
  ! -1 = fraction upwind
  ! 0 = magnitude
  ! 1 = unit x-dir
  ! 2 = unit y-dir
  ! 3 = unit z-dir
  PetscReal :: DiagThermalConductivityTensorToScalar
  PetscReal, intent(in)  :: dist(-1:3)

  PetscReal :: kx,ky,kz
  
  DiagThermalConductivityTensorToScalar = kx*dabs(dist(1))**2.0 + &
                                          ky*dabs(dist(2))**2.0 + &
                                          kz*dabs(dist(3))**2.0
  
end function DiagThermalConductivityTensorToScalar

! ************************************************************************** !

subroutine FullTensorCheckEigenvalues(kTR,option)
  !
  ! Check if full tensor is positive semi-definite
  !
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  PetscReal, dimension(3,3) :: kTR  ! anistropy ratios
  
  PetscReal :: a, b, c, d           ! characteristic function coefficients
  PetscReal :: kx,ky,kz,kxy,kxz,kyz ! symmetric tensor components
  PetscReal :: l(3)   ! eigenvalues
  PetscReal :: check  ! check eigenvalues in original characteristic function
  PetscInt  :: i
  
  ! symmetric tensor components
  kx  = kTR(1,1)
  ky  = kTR(2,2)
  kz  = kTR(3,3)
  kxy = kTR(1,2)
  kxz = kTR(1,3)
  kyz = kTR(2,3)
  
  l = 0.0d0
  
  ! cubic polynomial coefficients for characteristic function of tensor
  ! det(K - I.l) = f(l) = a l^3 + b l^2 + c l + d = 0
  ! l is an eigenvalue
  a = -1.0d0
  b = kx + ky + kz
  c = (kxy**2 + kxz**2 - kx*ky + kxy*kyz - kx*kz - ky*kz)
  d = (kxy**2)*kxz - (kxz**2)*ky - kx*kxy*kyz + kxy*kxz*kyz - (kxy**2)*kz + &
    kx*ky*kz
  
  call  CubicFormula(l,a,b,c,d,option)
  
  do i=1,3
    check = FullTensorCharacteristicFunction(l(i),kTR)
    if (l(i) < 0.0d0) then
      option%io_buffer = 'Thermal conductivity tensor is not positive'&
                       //' semi-definite.'
      call PrintErrMsg(option)
    endif
    if (l(i) > 1.0d0) then
      option%io_buffer = 'Eigenvalues of thermal conductivity tensor '&
                       //' indicate user input values may be exceeded '&
                       //' along principal axes.'
      call PrintMsg(option)
    endif
    if (check > 1.0d-1) then
      option%io_buffer = 'Could not determine positive semi-definite'&
                       //' thermal conductivity tensor.'
      call PrintMsg(option)
    endif
  enddo
  
end subroutine FullTensorCheckEigenvalues

! ************************************************************************** !

function FullTensorCharacteristicFunction(l,kT)
  !
  ! Check eigenvalues in characteristic function of symmetric tensor
  !
  implicit none
  PetscReal, dimension(3,3) :: kT
  PetscReal :: l,kx,ky,kz,kxy,kxz,kyz ! symmetric tensor components
  PetscReal :: FullTensorCharacteristicFunction
  
  ! symmetric tensor components
  kx  = kT(1,1)
  ky  = kT(2,2)
  kz  = kT(3,3)
  kxy = kT(1,2)
  kxz = kT(1,3)
  kyz = kT(2,3)
  
  ! characteristic function of tensor
  ! det(K - I.l) = f(l) = a l^3 + b l^2 + c l + d = 0
  
  FullTensorCharacteristicFunction = & 
    kxy**2*kxz - kxz**2*ky - kx*kxy*kyz + kxy*kxz*kyz - kxy**2*kz + kx*ky*kz + &
    (kxy**2 + kxz**2 - kx*ky + kxy*kyz - kx*kz - ky*kz)*l + &
    (kx + ky + kz)*l**2 - l**3
  
end function FullTensorCharacteristicFunction

! ************************************************************************** !

subroutine CubicFormula(roots,a,b,c,d,option)
  !
  ! Find real solutions of general cubic formula
  !
  use Option_module
  use PFLOTRAN_Constants_module

  implicit none

  type(option_type) :: option
  PetscReal, intent(in)    :: a, b, c, d
  PetscReal, intent(inout) :: roots(3)  
  PetscReal :: p, q, t, x, y
  PetscInt  :: k
  
  ! François Viète's trigonometric formula for three real roots 
  ! of cubic polynomial
  
  p = (3.0d0*a*c - b**2)/(3.0d0*(a**2))
  q = (2.0d0*(b**3) - 9.0d0*a*b*c + 27.0d0*d*(a**2))/(27.0d0*(a**3))
  
  ! Check if three real roots are applicable
  ! This must be true for the thermal conductivity tensor
  y = 4.0d0*(p**3) + 27.0d0*(q**2)
  
  if ( y >= 0.0d0 ) then
    option%io_buffer = 'Thermal conductivity tensor does not have three'&
                     //' real eigenvalues.'
    call PrintErrMsg(option)
  endif
  
  do k = 0,2,1
    t = 2.0d0*sqrt(-1*p/3.0d0)*cos((1.0d0/3.0d0)* &
      acos((3.0d0*q/(2.0d0*p))*sqrt(-3.0d0/p))-(2.0d0*pi*k/3.0d0))
    x = t - (b/(3.0d0*a))
    roots(k+1) = x
  enddo
  
end subroutine CubicFormula

! ************************************************************************** !

subroutine CharCurvesThermalAddToList(new_thermal_cc,list)

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

end subroutine CharCurvesThermalAddToList

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
        endif
      endif
      call OptionSetBlocking(option,PETSC_TRUE)
      call OptionCheckNonBlockingError(option)
    endif
    cur_thermal_cc => cur_thermal_cc%next
  enddo

end subroutine CharCurvesThermalConvertListToArray

! ************************************************************************** !

function CharCurvesThermalGetID(cc_thermal_array, &
     cc_thermal_name, material_property_name, option)

  use Option_module
  use String_module

  type(cc_thermal_ptr_type), pointer :: cc_thermal_array(:)
  character(len=MAXWORDLENGTH) :: cc_thermal_name
  character(len=MAXWORDLENGTH) :: test1, test2
  character(len=MAXWORDLENGTH) :: material_property_name
  type(option_type) :: option

  PetscInt :: CharCurvesThermalGetID
  PetscInt :: i, j
  
  do i = 1, size(cc_thermal_array)
      test1 = cc_thermal_array(i)%ptr%name
      do j = 1, size(cc_thermal_array)
        if (i == j) cycle
        test2 = cc_thermal_array(j)%ptr%name
        if (test1 == test2) then
          option%io_buffer = 'Duplicate thermal characteristic curve '//&
                             trim(test2)//&
                             ' has been detected.'
          call PrintErrMsg(option)
        endif
      enddo
  enddo

  CharCurvesThermalGetID = 0
  do CharCurvesThermalGetID = 1, size(cc_thermal_array)
    if (StringCompare(cc_thermal_name, &
         cc_thermal_array(CharCurvesThermalGetID)%ptr%name)) then
      return
    endif
  enddo
  option%io_buffer = 'Thermal characteristic curves "' // &
       trim(cc_thermal_name) // &
       '" in material property "' // &
       trim(material_property_name) // &
       '" not found among available thermal characteristic curves.'
  call PrintErrMsg(option)

end function CharCurvesThermalGetID

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

    write(id,'(a29)',advance='no') 'thermal char. curve name: '
    write(id,'(a)') adjustl(trim(cur_thermal_ccurve%name))

    if (associated(cur_thermal_ccurve%thermal_conductivity_function)) then
      write(id,'(a29)',advance='no') 'thermal conductivity func.: '
      select type (tcf => cur_thermal_ccurve%thermal_conductivity_function)
        !---------------------------------
      class is (kT_default_type)
        write(id,'(a)') 'only saturation-dependent (default)'
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
        write(id,'(a)') 'sat.- and temp.-dependent (power)'
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
        write(id,'(a)') 'sat.- and temp.-dependent (cubic polynomial)'
        write(id,'(a29)',advance='no') 'kT_wet: '
        write(word1,*) tcf%kT_wet
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'kT_dry: '
        write(word1,*) tcf%kT_dry
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'reference temp.: '
        write(word1,*) tcf%ref_temp
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'T coefficient: '
        write(word1,*) tcf%beta(1)
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'T^2 coefficient: '
        write(word1,*) tcf%beta(2)
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'T^3 coefficient: '
        write(word1,*) tcf%beta(3)
        write(id,'(a)') adjustl(trim(word1))
        !---------------------------------
      class is (kT_linear_resistivity_type)
        write(id,'(a)') 'sat.- and temp.-dependent (lin. resistivity)'
        write(id,'(a29)',advance='no') 'kT_wet: '
        write(word1,*) tcf%kT_wet
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'kT_dry: '
        write(word1,*) tcf%kT_dry
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'reference temp.: '
        write(word1,*) tcf%ref_temp
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'const. coefficient: '
        write(word1,*) tcf%a(1)
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'T coefficient: '
        write(word1,*) tcf%a(2)
        write(id,'(a)') adjustl(trim(word1))
      end select
    endif

    write(id,'(a29)') '---------------------------: '
    cur_thermal_ccurve => cur_thermal_ccurve%next
  enddo

end subroutine CharCurvesThermalInputRecord

! ************************************************************************** !

recursive subroutine CharCurvesThermalDestroy(tcc)

  implicit none

  class(cc_thermal_type), pointer :: tcc

  if (.not.associated(tcc)) return

  call CharCurvesThermalDestroy(tcc%next)

  call TCFDestroy(tcc%thermal_conductivity_function)

  deallocate(tcc)
  nullify(tcc)

end subroutine CharCurvesThermalDestroy

end module Characteristic_Curves_Thermal_module
