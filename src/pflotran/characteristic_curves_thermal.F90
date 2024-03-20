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
    PetscReal :: alpha ! exponent for soil Kersten number (TH mode)
  contains
    procedure, public :: Verify => TCFBaseVerify
    procedure, public :: Test => TCFBaseTest
    procedure, public :: CalculateTCond => TCFBaseConductivity
    procedure, public :: CalculateFTCond => TCFBaseConductivity2
    procedure, public :: TCondTensorToScalar
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
    PetscReal :: kT_x, kT_y, kT_z, kT_xy, kT_xz, kT_yz
    PetscReal :: kT(3,3,3)   ! thermal conductivity tensor
    PetscReal :: kTf(3,3)    ! anisotropy ratio tensor
    PetscBool :: isotropic
    PetscBool :: full_tensor
  contains
    procedure, public :: Verify => TCFDefaultVerify
    procedure, public :: CalculateTCond => TCFDefaultConductivity
  end type kT_default_type
  !---------------------------------------------------------------------------
  type, public, extends(kT_default_type) :: kT_linear_type
  contains
    procedure, public :: Verify => TCFLinearVerify
    procedure, public :: CalculateTCond => TCFLinearConductivity
  end type kT_linear_type
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
    PetscReal :: ref_temp, ref_por, por_exp
    PetscReal :: a(2)  ! 1/(a(1) + a(2)*T)
    PetscReal :: b(2)  ! b(1) + b(2)*T
    PetscBool :: porosity_effect
  contains
    procedure, public :: Verify => TCFLinearResistivityVerify
    procedure, public :: CalculateTCond => TCFLinearResistivityConductivity
  end type kT_linear_resistivity_type
  !---------------------------------------------------------------------------
  type, public, extends(kT_default_type) :: kT_frozen_type
    PetscReal :: kT_frozen  ! frozen thermal conductivity
    PetscReal :: alpha_fr   ! exponent for frozen soil Kersten number
    PetscInt  :: ice_model  ! indicator of ice model
  contains
    procedure, public :: Verify => TCFFrozenVerify
    procedure, public :: Test => TCFFrozenTest  ! test with ice saturation
    procedure, public :: CalculateTCond => TCFFrozenConductivity1   ! freezing inactive
    procedure, public :: CalculateFTCond => TCFFrozenConductivity2  ! freezing active
  end type kT_frozen_type
  !---------------------------------------------------------------------------
  type, public, extends(kT_default_type) :: kT_ASM_dry_type
    PetscReal :: dry_T_coeff ! coefficient of T-dependent term
    PetscReal :: dry_T_power ! exponent of temperature in T-dependent term
  contains
    procedure, public :: Verify => TCFASMDryVerify
    procedure, public :: CalculateTCond => TCFASMDryConductivity
  end type kT_ASM_dry_type
  !---------------------------------------------------------------------------
  type, public, extends(kT_default_type) :: kT_ASM_water_filled_type
    PetscReal :: kT_water     ! thermal conductivity of water
    PetscReal :: kT_solid     ! thermal conductivity of solid components
    PetscReal :: porosity_asm ! porosity of the assembly
  contains
    procedure, public :: Verify => TCFASMWaterFilledVerify
    procedure, public :: CalculateTCond => TCFASMWaterFilledConductivity
  end type kT_ASM_water_filled_type
  !---------------------------------------------------------------------------
  type, public, extends(kT_default_type) :: kT_ASM_radial_type
    PetscReal :: dry_T_coeff  ! coefficient of T-dependent term
    PetscReal :: dry_T_power  ! exponent of temperature in T-dependent term
    PetscReal :: kT_water     ! thermal conductivity of water
    PetscReal :: kT_solid     ! thermal conductivity of solid components
    PetscReal :: porosity_asm ! porosity of the assembly
  contains
    procedure, public :: Verify => TCFASMRadialVerify
    procedure, public :: CalculateTCond => TCFASMRadialConductivity
  end type kT_ASM_radial_type
  !---------------------------------------------------------------------------
  type, public, extends(kT_default_type) :: kT_ASM_axial_type
    PetscReal :: kT_water     ! thermal conductivity of water
    PetscReal :: kT_solid     ! thermal conductivity of solid components
    PetscReal :: porosity_asm ! porosity of the assembly
  contains
    procedure, public :: Verify => TCFASMAxialVerify
    procedure, public :: CalculateTCond => TCFASMAxialConductivity
  end type kT_ASM_axial_type
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
  !---------------------------------------------------------------------------
  type, public, extends(kT_default_type) :: kT_composite_type
    PetscReal :: dist(-1:3)
    ! -1 = fraction upwind
    ! 0 = magnitude
    ! 1 = unit x-dir
    ! 2 = unit y-dir
    ! 3 = unit z-dir
    character(len=MAXWORDLENGTH) :: lkT(3) ! thermal conductivity function names
    type(cc_thermal_ptr_type) :: ckT(3)    ! thermal conductivity functions
  contains
    procedure, public :: Test => TCFCompositeTest
    procedure, public :: Verify => TCFCompositeVerify
    procedure, public :: CalculateTCond => TCFCompositeConductivity
  end type kT_composite_type
  !---------------------------------------------------------------------------

  public :: CharCurvesThermalCreate, &
            CharCurvesThermalGetID, &
            CharCurvesThermalRead, &
            CharCurvesThermalAddToList, &
            CharCurvesThermalConvertListToArray, &
            CharCurvesThermalInputRecord, &
            CharCurvesThermalDestroy, &
            CompositeTCCList, &
            TCFDefaultCreate, &
            TCFConstantCreate, &
            TCFPowerCreate, &
            TCFCubicPolynomialCreate, &
            TCFLinearResistivityCreate, &
            TCFASMDryCreate, &
            TCFASMWaterFilledCreate, &
            TCFASMRadialCreate, &
            TCFASMAxialCreate, &
            TCFCompositeCreate, &
            TCFFrozenCreate, &
            TCFAssignDefault, &
            TCFAssignFrozen

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

subroutine TCFBaseConductivity(this,liquid_saturation,temperature,porosity, &
     thermal_conductivity,dkT_dsatl,dkT_dtemp,option)

  use Option_module

  implicit none

  class(thermal_conductivity_base_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature, porosity
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

subroutine TCFBaseConductivity2(this,liquid_saturation,ice_saturation, &
     temperature,porosity,thermal_conductivity,dkT_dsatl,dkT_dsati,dkT_dtemp,option)

  use Option_module

  implicit none

  class(thermal_conductivity_base_type) :: this
  PetscReal, intent(in) :: liquid_saturation, ice_saturation, temperature, porosity
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dkT_dsatl, dkT_dsati, dkT_dtemp
  type(option_type), intent(inout) :: option

  thermal_conductivity = 0.d0
  dkT_dsatl = 0.d0
  dkT_dsati = 0.d0
  dkT_dtemp = 0.d0

  option%io_buffer = 'Base thermal conductivity must be extended for ' &
                   //'frozen parameters.'
  call PrintErrMsg(option)

end subroutine TCFBaseConductivity2

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
  PetscReal :: deltaTemp, deltaSat, deltaPor
  PetscReal :: temp0, sat0, por0, pmult
  PetscReal :: kT, dkT_dsat, dkT_dsat_numerical
  PetscReal :: dkT_dtemp, dkT_dtemp_numerical
  PetscReal :: dkT_dpor_numerical
  PetscReal :: perturbed_temp, perturbed_sat, perturbed_por
  PetscReal :: kT_temp_pert, kT_sat_pert, kT_por_pert, unused1, unused2
  PetscReal :: temp_min, temp_max, sat_min, sat_max, por_min, por_max
  PetscInt :: i,j,k, np

  ! thermal conductivity as a function of temp., liq. sat., and porosity
  temp_min = 1.0d0 ! Celsius
  temp_max = 250.0d0
  sat_min = 1.0d-3
  sat_max = 1.0d0

  select type(this)
  class is(kT_linear_resistivity_type)
     ! this is the only tcc hwere porosity-variation is implemented
     por_max = this%ref_por
     por_min = min(por_max / 10.0, 1.0D-4)
     np = 21
     deltaPor = (por_max - por_min)/(np - 1)
  class default
     por_min = 1.0D-4
     por_max = por_min
     deltaPor = 0.d0
     np = 1
  end select

  deltaTemp = (temp_max - temp_min)/(nt - 1)
  deltaSat = (sat_max - sat_min)/(ns - 1)

  write(string,*) tcc_name
  string = trim(tcc_name) // '_kT_vs_sat_temp_and_por.dat'
  open(unit=86,file=string)
  write(86,*) '"temperature [C]", "liquid saturation [-]", "porosity [-]",' &
       //'"kT [W/m*K]", "dkT/dsat", "dkT/dT", "dkT/dsat_numerical",' &
       //'"dkT/dT_numerical", "dkT/dpor_numerical"'

  do i = 1,nt
     temp0 = temp_min + deltaTemp * (i - 1)
     do j = 1,ns
        sat0 = sat_min + deltaSat * (j - 1)
        do k = 1,np
          por0 = por_min + deltaPor * (k - 1)

          ! base case with analytical derivatives
          call this%CalculateTCond(sat0,temp0,por0, &
               kT,dkT_dsat,dkT_dtemp,option)

          ! calculate numerical derivatives via finite differences
          perturbed_temp = temp0 * (1.d0 + perturbation)
          call this%CalculateTCond(sat0,perturbed_temp,por0, &
               kT_temp_pert,unused1,unused2,option)

          dkT_dtemp_numerical = (kT_temp_pert - kT)/(temp0 * &
               perturbation)

          if (j == ns) then
             pmult = -1.d0  ! backwards diff
          else
             pmult =1.d0
          endif
          perturbed_sat = sat0 * (1.d0 + pmult * perturbation)
          call this%CalculateTCond(perturbed_sat,temp0,por0, &
               kT_sat_pert,unused1,unused2,option)

          dkT_dsat_numerical = (pmult * kT_sat_pert - pmult * kT)/(sat0 * &
               perturbation)

          if (k > 1 .and. k == np) then
             pmult = -1.d0
          else
             pmult = 1.d0
          endif
          perturbed_por = por0 * (1.d0 + pmult * perturbation)
          call this%CalculateTCond(sat0,temp0,perturbed_por, &
               kT_por_pert,unused1,unused2,option)

          dkT_dpor_numerical = (pmult * kT_por_pert - pmult * kT)/(por0 * &
               perturbation)

          write(86,'(9(ES14.6,1X))') temp0, sat0, por0, &
               kT, dkT_dsat, dkT_dtemp, dkT_dsat_numerical, &
               dkT_dtemp_numerical, dkT_dpor_numerical
       enddo
    enddo
 enddo
 close(86)

end subroutine TCFBaseTest

! ************************************************************************** !

subroutine TCFFrozenTest(this,tcc_name,option)

  use Option_module

  implicit none

  class(kT_frozen_type) :: this
  character(len=MAXWORDLENGTH) :: tcc_name
  type(option_type), intent(inout) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt, parameter :: nt = 28
  PetscInt, parameter :: ns = 12
  PetscInt, parameter :: ni = 12
  PetscReal, parameter :: perturbation = 1.0D-6
  PetscReal :: deltaTemp, deltaSat, deltaIce
  PetscReal :: temp_vec(nt)
  PetscReal :: sat_vec(ns)
  PetscReal :: ice_vec(ni)
  PetscReal :: kT(nt,ns,ni)
  PetscReal :: dkT_dsat(nt,ns,ni)
  PetscReal :: dkT_dsat_numerical(nt,ns,ni)
  PetscReal :: dkT_dice(nt,ns,ni)
  PetscReal :: dkT_dice_numerical(nt,ns,ni)
  PetscReal :: dkT_dtemp(nt,ns,ni)
  PetscReal :: dkT_dtemp_numerical(nt,ns,ni)
  PetscReal :: perturbed_temp, perturbed_sat, perturbed_ice
  PetscReal :: kT_temp_pert, kT_sat_pert, kT_ice_pert
  PetscReal :: unused1, unused2, unused3
  PetscReal :: temp_min, temp_max, sat_min, sat_max, ice_min, ice_max
  PetscInt :: i,j,k

  ! resort to regular test if frozen thermal conductivity not initialized
  if (Uninitialized(this%kT_frozen)) then
    call TCFBaseTest(this,tcc_name,option)
    return
  endif

  ! thermal conductivity as a function of temp. and liq. sat.
  temp_min = 1.0d0 ! Celsius
  temp_max = 250.0d0
  sat_min = 1.0d-3
  sat_max = 1.0d0
  ice_min = 1.0d-3
  ice_max = 1.0d0

  deltaTemp = (temp_max - temp_min)/(nt - 1)
  deltaSat = (sat_max - sat_min)/(ns - 1)
  deltaIce = (ice_max - ice_min)/(ni - 1)

  temp_vec = [(temp_min + i*deltaTemp, i=0,nt-1)]
  sat_vec = [(sat_min + i*deltaSat, i=0,ns-1)]
  ice_vec = [(ice_min + i*deltaIce, i=0,ni-1)]

  do i = 1,nt
    do j = 1,ns
      do k = 1,ni
        ! base case with analytical derivatives
        call this%CalculateFTCond(sat_vec(j),ice_vec(k),temp_vec(i), &
           0.d0,kT(i,j,k),dkT_dsat(i,j,k),dkT_dice(i,j,k),dkT_dtemp(i,j,k),option)

        ! calculate numerical derivatives via finite differences
        perturbed_temp = temp_vec(i) * (1.d0 + perturbation)
        call this%CalculateFTCond(sat_vec(j),ice_vec(k),perturbed_temp, &
             0.d0,kT_temp_pert,unused1,unused2,unused3,option)

        dkT_dtemp_numerical(i,j,k) = (kT_temp_pert - kT(i,j,k))/ &
                                   (temp_vec(i)*perturbation)

        perturbed_sat = sat_vec(j) * (1.d0 + perturbation)
        call this%CalculateFTCond(perturbed_sat,ice_vec(k),temp_vec(i), &
             0.d0,kT_sat_pert,unused1,unused2,unused3,option)

        dkT_dsat_numerical(i,j,k) = (kT_sat_pert - kT(i,j,k))/ &
                                  (sat_vec(j)*perturbation)

        perturbed_ice = ice_vec(k) * (1.d0 + perturbation)
        call this%CalculateFTCond(sat_vec(j),perturbed_ice,temp_vec(i), &
             0.d0,kT_ice_pert,unused1,unused2,unused3,option)

        dkT_dice_numerical(i,j,k) = (kT_ice_pert - kT(i,j,k))/ &
                                   (ice_vec(k)*perturbation)
      enddo
    enddo
  enddo

  write(string,*) tcc_name
  string = trim(tcc_name) // '_kT_vs_sat_and_temp.dat'
  open(unit=86,file=string)
  write(86,*) '"temperature [C]", "liquid saturation [-]",' &
               //'"ice saturation [-]", "kT [W/m*K]", "dkT/dsatl", "dkT/dsati",' &
               //'"dkT/dT", "dkT/dsatl_numerical", "dkT/dsati_numerical",' &
               //'"dkT/dT_numerical"'
  do i = 1,nt
    do j = 1,ns
      do k = 1,ni
        write(86,'(10(ES14.6))') temp_vec(i), sat_vec(j), ice_vec(k), &
             kT(i,j,k), dkT_dsat(i,j,k), dkT_dice(i,j,k), dkT_dtemp(i,j,k), &
             dkT_dsat_numerical(i,j,k), dkT_dice_numerical(i,j,k), &
             dkT_dtemp_numerical(i,j,k)
      enddo
    enddo
  enddo
  close(86)

end subroutine TCFFrozenTest

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
  TCFDefaultCreate%isotropic   = PETSC_TRUE
  TCFDefaultCreate%full_tensor = PETSC_FALSE
  TCFDefaultCreate%kT_wet = UNINITIALIZED_DOUBLE
  TCFDefaultCreate%kT_dry = UNINITIALIZED_DOUBLE
  TCFDefaultCreate%alpha  = 1.0d0
  TCFDefaultCreate%kT     = UNINITIALIZED_DOUBLE
  TCFDefaultCreate%kTf    = UNINITIALIZED_DOUBLE
  TCFDefaultCreate%kT_x   = UNINITIALIZED_DOUBLE
  TCFDefaultCreate%kT_y   = UNINITIALIZED_DOUBLE
  TCFDefaultCreate%kT_z   = UNINITIALIZED_DOUBLE
  TCFDefaultCreate%kT_xy  = UNINITIALIZED_DOUBLE
  TCFDefaultCreate%kT_xz  = UNINITIALIZED_DOUBLE
  TCFDefaultCreate%kT_yz  = UNINITIALIZED_DOUBLE

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

subroutine TCFDefaultConductivity(this,liquid_saturation,temperature,porosity, &
     thermal_conductivity,dkT_dsatl,dkT_dtemp,option)

  use Option_module

  implicit none

  class(kT_default_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature, porosity
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

function TCFLinearCreate()

  implicit none

  class(kT_linear_type), pointer :: TCFLinearCreate

  allocate(TCFLinearCreate)
  TCFLinearCreate%isotropic   = PETSC_TRUE
  TCFLinearCreate%full_tensor = PETSC_FALSE
  TCFLinearCreate%kT_wet = UNINITIALIZED_DOUBLE
  TCFLinearCreate%kT_dry = UNINITIALIZED_DOUBLE
  TCFLinearCreate%alpha  = 1.0d0
  TCFLinearCreate%kT     = UNINITIALIZED_DOUBLE
  TCFLinearCreate%kTf    = UNINITIALIZED_DOUBLE
  TCFLinearCreate%kT_x   = UNINITIALIZED_DOUBLE
  TCFLinearCreate%kT_y   = UNINITIALIZED_DOUBLE
  TCFLinearCreate%kT_z   = UNINITIALIZED_DOUBLE
  TCFLinearCreate%kT_xy  = UNINITIALIZED_DOUBLE
  TCFLinearCreate%kT_xz  = UNINITIALIZED_DOUBLE
  TCFLinearCreate%kT_yz  = UNINITIALIZED_DOUBLE

end function TCFLinearCreate

! ************************************************************************** !

subroutine TCFLinearVerify(this,name,option)

  use Option_module

  implicit none

  class(kT_linear_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'THERMAL_CONDUCTIVITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'THERMAL_CONDUCTIVITY_FUNCTION,LINEAR'
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

end subroutine TCFLinearVerify

! ************************************************************************** !

subroutine TCFLinearConductivity(this,liquid_saturation,temperature,porosity, &
                                 thermal_conductivity,dkT_dsatl,dkT_dtemp, &
                                 option)

  use Option_module

  implicit none

  class(kT_linear_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature, porosity
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dkT_dsatl, dkT_dtemp
  type(option_type), intent(inout) :: option

  PetscReal :: tempreal

  dkT_dtemp = 0.d0 ! only a function of saturation

  if (liquid_saturation > 0.d0) then
    tempreal = liquid_saturation * &
               (this%kT_wet - this%kT_dry)
    thermal_conductivity = this%kT_dry + tempreal
    dkT_dsatl = 0.5d0 * tempreal / liquid_saturation
  else
    thermal_conductivity = this%kT_dry
    dkT_dsatl = 0.d0
  endif

end subroutine TCFLinearConductivity

! ************************************************************************** !

function TCFASMDryCreate()

  implicit none

  class(kT_ASM_dry_type), pointer :: TCFASMDryCreate

  allocate(TCFASMDryCreate)
  ! User doesn't need to initialize wet thermal conditivity for this case
  TCFASMDryCreate%kT_wet         = 0.0d0
  TCFASMDryCreate%alpha          = 1.0d0
  TCFASMDryCreate%kT_dry         = UNINITIALIZED_DOUBLE
  TCFASMDryCreate%dry_T_coeff    = UNINITIALIZED_DOUBLE
  TCFASMDryCreate%dry_T_power    = UNINITIALIZED_DOUBLE
  TCFASMDryCreate%kT             = UNINITIALIZED_DOUBLE
  TCFASMDryCreate%kTf            = UNINITIALIZED_DOUBLE
  TCFASMDryCreate%kT_x           = UNINITIALIZED_DOUBLE
  TCFASMDryCreate%kT_y           = UNINITIALIZED_DOUBLE
  TCFASMDryCreate%kT_z           = UNINITIALIZED_DOUBLE
  TCFASMDryCreate%kT_xy          = UNINITIALIZED_DOUBLE
  TCFASMDryCreate%kT_xz          = UNINITIALIZED_DOUBLE
  TCFASMDryCreate%kT_yz          = UNINITIALIZED_DOUBLE
  TCFASMDryCreate%isotropic      = PETSC_TRUE
  TCFASMDryCreate%full_tensor    = PETSC_FALSE

end function TCFASMDryCreate

! ************************************************************************** !

subroutine TCFASMDryVerify(this,name,option)

  use Option_module

  implicit none

  class(kT_ASM_dry_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'THERMAL_CONDUCTIVITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'THERMAL_CONDUCTIVITY_FUNCTION,ASM_DRY'
  endif
  call TCFBaseVerify(this,string,option)
  if (Uninitialized(this%dry_T_coeff)) then
    option%io_buffer = UninitializedMessage('ASM_DRY_COEFFICIENT',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%dry_T_power)) then
    option%io_buffer = UninitializedMessage('ASM_DRY_EXPONENT',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%kT_dry)) then
    option%io_buffer = UninitializedMessage('THERMAL_CONDUCTIVITY_DRY',string)
    call PrintErrMsg(option)
  endif

end subroutine TCFASMDryVerify

! ************************************************************************** !

subroutine TCFASMDryConductivity(this,liquid_saturation,temperature,porosity, &
     thermal_conductivity,dkT_dsatl,dkT_dtemp,option)

  use Option_module

  implicit none

  class(kT_ASM_dry_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature, porosity
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dkT_dsatl, dkT_dtemp
  type(option_type), intent(inout) :: option

  PetscReal :: T, tempreal, k1, k2

  T = temperature

  ! Reference: Spent Nuclear Fuel Effective Thermal Conductivity
  !            Report DI: BBA000000-017 17-5705-00010 REV 00

  ! If user did not initialize wet thermal conductivty, remove saturation
  !   dependence
  if (this%kT_wet == 0.0d0) then
    this%kT_wet = this%kT_dry
  end if

  if (T < 0.0d0) then
    T = 0.0d0 ! avoid complex number in analysis
  endif

  k1 = this%dry_T_coeff * (T ** this%dry_T_power)

  k2 = this%kT_dry + k1

  dkT_dtemp = this%dry_T_coeff * this%dry_T_power * &
    (T ** (this%dry_T_power - 1.0d0))

  if (liquid_saturation > 0.d0) then
    tempreal = sqrt(liquid_saturation) * &
                   (this%kT_wet - k2)
    thermal_conductivity = k2 + tempreal
    dkT_dsatl = 0.5d0 * tempreal / liquid_saturation
  else
    thermal_conductivity = k2
    dkT_dsatl = 0.d0
  endif

end subroutine TCFASMDryConductivity

! ************************************************************************** !

function TCFASMWaterFilledCreate()

  implicit none

  class(kT_ASM_water_filled_type), pointer :: &
    TCFASMWaterFilledCreate

  allocate(TCFASMWaterFilledCreate)
  ! User doesn't need to initialize dry or wet thermal conditivity for this case
  TCFASMWaterFilledCreate%kT_dry         = 0.0d0
  TCFASMWaterFilledCreate%kT_wet         = 0.0d0
  TCFASMWaterFilledCreate%alpha          = 1.0d0
  TCFASMWaterFilledCreate%kT_water       = UNINITIALIZED_DOUBLE
  TCFASMWaterFilledCreate%kT_solid       = UNINITIALIZED_DOUBLE
  TCFASMWaterFilledCreate%porosity_asm   = UNINITIALIZED_DOUBLE
  TCFASMWaterFilledCreate%kT             = UNINITIALIZED_DOUBLE
  TCFASMWaterFilledCreate%kTf            = UNINITIALIZED_DOUBLE
  TCFASMWaterFilledCreate%kT_x           = UNINITIALIZED_DOUBLE
  TCFASMWaterFilledCreate%kT_y           = UNINITIALIZED_DOUBLE
  TCFASMWaterFilledCreate%kT_z           = UNINITIALIZED_DOUBLE
  TCFASMWaterFilledCreate%kT_xy          = UNINITIALIZED_DOUBLE
  TCFASMWaterFilledCreate%kT_xz          = UNINITIALIZED_DOUBLE
  TCFASMWaterFilledCreate%kT_yz          = UNINITIALIZED_DOUBLE
  TCFASMWaterFilledCreate%isotropic      = PETSC_TRUE
  TCFASMWaterFilledCreate%full_tensor    = PETSC_FALSE

end function TCFASMWaterFilledCreate

! ************************************************************************** !

subroutine TCFASMWaterFilledVerify(this,name,option)

  use Option_module

  implicit none

  class(kT_ASM_water_filled_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'THERMAL_CONDUCTIVITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // &
      'THERMAL_CONDUCTIVITY_FUNCTION,ASM_WATER_FILLED'
  endif
  call TCFBaseVerify(this,string,option)
  if (Uninitialized(this%kT_water)) then
    option%io_buffer = UninitializedMessage('THERMAL_CONDUCTIVITY_WATER',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%kT_solid)) then
    option%io_buffer = UninitializedMessage('THERMAL_CONDUCTIVITY_SOLID',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%porosity_asm)) then
    option%io_buffer = UninitializedMessage('POROSITY_ASSEMBLY',string)
    call PrintErrMsg(option)
  endif

end subroutine TCFASMWaterFilledVerify

! ************************************************************************** !

subroutine TCFASMWaterFilledConductivity(this,liquid_saturation, &
     temperature,porosity,thermal_conductivity,dkT_dsatl,dkT_dtemp,option)

  use Option_module

  implicit none

  class(kT_ASM_water_filled_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature, porosity
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dkT_dsatl, dkT_dtemp
  type(option_type), intent(inout) :: option

  PetscReal :: tempreal, lamda, v1, v2, ratio

  ! Reference: equation 15 in Cheng and Hsu, 1999

  lamda = this%kT_water / this%kT_solid

  v1 = sqrt(1 - this%porosity_asm)

  v2 = 1 + (lamda - 1) * v1

  ratio = 1 - v1 + ( v1 / v2 )

  this%kT_wet = ratio * this%kT_water

  ! Allow user to impart saturation dependence
  if (this%kT_dry == 0.0d0) then
    this%kT_dry = this%kT_wet ! remove saturation dependence if user input N/A
  end if

  if (liquid_saturation > 0.d0) then
    tempreal = sqrt(liquid_saturation) * &
                   (this%kT_wet - this%kT_dry)
    thermal_conductivity = this%kT_dry + tempreal
    dkT_dsatl = 0.5d0 * tempreal / liquid_saturation
  else
    thermal_conductivity = this%kT_dry
    dkT_dsatl = 0.d0
  endif

  dkT_dtemp = 0.0d0

end subroutine TCFASMWaterFilledConductivity

! ************************************************************************** !

function TCFASMRadialCreate()

  implicit none

  class(kT_ASM_radial_type), pointer :: TCFASMRadialCreate

  allocate(TCFASMRadialCreate)
  ! User doesn't need to initialize wet thermal conditivity for this case
  TCFASMRadialCreate%kT_wet         = 0.0d0
  TCFASMRadialCreate%alpha          = 1.0d0
  TCFASMRadialCreate%kT_dry         = UNINITIALIZED_DOUBLE
  TCFASMRadialCreate%dry_T_coeff    = UNINITIALIZED_DOUBLE
  TCFASMRadialCreate%dry_T_power    = UNINITIALIZED_DOUBLE
  TCFASMRadialCreate%kT_water       = UNINITIALIZED_DOUBLE
  TCFASMRadialCreate%kT_solid       = UNINITIALIZED_DOUBLE
  TCFASMRadialCreate%porosity_asm   = UNINITIALIZED_DOUBLE
  TCFASMRadialCreate%kT             = UNINITIALIZED_DOUBLE
  TCFASMRadialCreate%kTf            = UNINITIALIZED_DOUBLE
  TCFASMRadialCreate%kT_x           = UNINITIALIZED_DOUBLE
  TCFASMRadialCreate%kT_y           = UNINITIALIZED_DOUBLE
  TCFASMRadialCreate%kT_z           = UNINITIALIZED_DOUBLE
  TCFASMRadialCreate%kT_xy          = UNINITIALIZED_DOUBLE
  TCFASMRadialCreate%kT_xz          = UNINITIALIZED_DOUBLE
  TCFASMRadialCreate%kT_yz          = UNINITIALIZED_DOUBLE
  TCFASMRadialCreate%isotropic      = PETSC_TRUE
  TCFASMRadialCreate%full_tensor    = PETSC_FALSE

end function TCFASMRadialCreate

! ************************************************************************** !

subroutine TCFASMRadialVerify(this,name,option)

  use Option_module

  implicit none

  class(kT_ASM_radial_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'THERMAL_CONDUCTIVITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'THERMAL_CONDUCTIVITY_FUNCTION,ASM_RADIAL'
  endif
  call TCFBaseVerify(this,string,option)
  if (Uninitialized(this%dry_T_coeff)) then
    option%io_buffer = UninitializedMessage('ASM_DRY_COEFFICIENT',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%dry_T_power)) then
    option%io_buffer = UninitializedMessage('ASM_DRY_EXPONENT',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%kT_dry)) then
    option%io_buffer = UninitializedMessage('THERMAL_CONDUCTIVITY_DRY',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%kT_water)) then
    option%io_buffer = UninitializedMessage('THERMAL_CONDUCTIVITY_WATER',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%kT_solid)) then
    option%io_buffer = UninitializedMessage('THERMAL_CONDUCTIVITY_SOLID',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%porosity_asm)) then
    option%io_buffer = &
      UninitializedMessage('POROSITY_ASSEMBLY',string)
    call PrintErrMsg(option)
  end if

end subroutine TCFASMRadialVerify

! ************************************************************************** !

subroutine TCFASMRadialConductivity(this,liquid_saturation,temperature, &
     porosity,thermal_conductivity,dkT_dsatl,dkT_dtemp,option)

  use Option_module

  implicit none

  class(kT_ASM_radial_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature, porosity
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dkT_dsatl, dkT_dtemp
  type(option_type), intent(inout) :: option

  PetscReal :: T, tempreal, lamda, v1, v2, ratio, k1, k2

  T = temperature

  ! Calculation of wet thermal conductivity

  ! Reference: equation 15 in Cheng and Hsu, 1999

  lamda = this%kT_water / this%kT_solid

  v1 = sqrt(1 - this%porosity_asm)

  v2 = 1 + (lamda - 1) * v1

  ratio = 1 - v1 + ( v1 / v2 )

  this%kT_wet = ratio * this%kT_water

  ! Calculation of dry thermal conductivity

  ! Reference: Spent Nuclear Fuel Effective Thermal Conductivity
  !            Report DI: BBA000000-017 17-5705-00010 REV 00

  if (this%kT_wet <= 0.0d0) then
    option%io_buffer = 'Error in wet thermal conductivity derivation' &
                     //' for assembly conditions'
    call PrintErrMsg(option)
  end if

  if (T < 0.0d0) then
    T = 0.0d0 ! avoid complex number in analysis
  endif

  k1 = this%dry_T_coeff * (T ** this%dry_T_power)

  k2 = this%kT_dry + k1

  dkT_dtemp = this%dry_T_coeff * this%dry_T_power * &
    (T ** (this%dry_T_power - 1.0d0))

  if (liquid_saturation > 0.d0) then
    tempreal = sqrt(liquid_saturation) * &
                   (this%kT_wet - k2)
    thermal_conductivity = k2 + tempreal
    dkT_dsatl = 0.5d0 * tempreal / liquid_saturation
  else
    thermal_conductivity = k2
    dkT_dsatl = 0.d0
  endif

end subroutine TCFASMRadialConductivity

! ************************************************************************** !

function TCFASMAxialCreate()

  implicit none

  class(kT_ASM_axial_type), pointer :: TCFASMAxialCreate

  allocate(TCFASMAxialCreate)
  ! User doesn't need to initialize dry or wet thermal conditivity for this case
  TCFASMAxialCreate%kT_dry        = 0.0d0
  TCFASMAxialCreate%kT_wet        = 0.0d0
  TCFASMAxialCreate%alpha         = 1.0d0
  TCFASMAxialCreate%kT_solid      = UNINITIALIZED_DOUBLE
  TCFASMAxialCreate%kT_water      = UNINITIALIZED_DOUBLE
  TCFASMAxialCreate%porosity_asm  = UNINITIALIZED_DOUBLE
  TCFASMAxialCreate%kT            = UNINITIALIZED_DOUBLE
  TCFASMAxialCreate%kTf           = UNINITIALIZED_DOUBLE
  TCFASMAxialCreate%kT_x          = UNINITIALIZED_DOUBLE
  TCFASMAxialCreate%kT_y          = UNINITIALIZED_DOUBLE
  TCFASMAxialCreate%kT_z          = UNINITIALIZED_DOUBLE
  TCFASMAxialCreate%kT_xy         = UNINITIALIZED_DOUBLE
  TCFASMAxialCreate%kT_xz         = UNINITIALIZED_DOUBLE
  TCFASMAxialCreate%kT_yz         = UNINITIALIZED_DOUBLE
  TCFASMAxialCreate%isotropic     = PETSC_TRUE
  TCFASMAxialCreate%full_tensor   = PETSC_FALSE

end function TCFASMAxialCreate

! ************************************************************************** !

subroutine TCFASMAxialVerify(this,name,option)

  use Option_module

  implicit none

  class(kT_ASM_axial_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'THERMAL_CONDUCTIVITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'THERMAL_CONDUCTIVITY_FUNCTION,ASM_AXIAL'
  endif
  call TCFBaseVerify(this,string,option)
  if (Uninitialized(this%kT_water)) then
    option%io_buffer = UninitializedMessage('THERMAL_CONDUCTIVITY_WATER',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%kT_solid)) then
    option%io_buffer = UninitializedMessage('THERMAL_CONDUCTIVITY_SOLID',string)
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%porosity_asm)) then
    option%io_buffer = UninitializedMessage('POROSITY_ASSEMBLY',string)
    call PrintErrMsg(option)
  endif

end subroutine TCFASMAxialVerify

! ************************************************************************** !

subroutine TCFASMAxialConductivity(this,liquid_saturation,temperature, &
     porosity,thermal_conductivity,dkT_dsatl,dkT_dtemp,option)

  use Option_module

  implicit none

  class(kT_ASM_axial_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature, porosity
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dkT_dsatl, dkT_dtemp
  type(option_type), intent(inout) :: option

  PetscReal :: k1, k2

  ! Assumes parallel heat condition of DPC internals and water

  k1 = this%kT_solid * (1 - this%porosity_asm)

  k2 = this%porosity_asm * this%kT_water * liquid_saturation

  thermal_conductivity = k1 + k2

  dkT_dtemp = 0.d0

  dkT_dsatl = this%porosity_asm * this%kT_water

end subroutine TCFASMAxialConductivity

! ************************************************************************** !

function TCFConstantCreate()

  implicit none

  class(kT_constant_type), pointer :: TCFConstantCreate

  allocate(TCFConstantCreate)
  TCFConstantCreate%constant_thermal_conductivity = UNINITIALIZED_DOUBLE
  TCFConstantCreate%alpha = 1.0d0

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
     porosity,thermal_conductivity,dkT_dsatl,dkT_dtemp,option)

  use Option_module

  implicit none

  class(kT_constant_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature, porosity
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
  TCFPowerCreate%ref_temp = -T273K
  TCFPowerCreate%gamma = UNINITIALIZED_DOUBLE
  TCFPowerCreate%alpha = 1.0d0
  TCFPowerCreate%isotropic   = PETSC_TRUE
  TCFPowerCreate%full_tensor = PETSC_FALSE
  TCFPowerCreate%kT     = UNINITIALIZED_DOUBLE
  TCFPowerCreate%kTf    = UNINITIALIZED_DOUBLE
  TCFPowerCreate%kT_x   = UNINITIALIZED_DOUBLE
  TCFPowerCreate%kT_y   = UNINITIALIZED_DOUBLE
  TCFPowerCreate%kT_z   = UNINITIALIZED_DOUBLE
  TCFPowerCreate%kT_xy  = UNINITIALIZED_DOUBLE
  TCFPowerCreate%kT_xz  = UNINITIALIZED_DOUBLE
  TCFPowerCreate%kT_yz  = UNINITIALIZED_DOUBLE

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
     porosity,thermal_conductivity,dkT_dsatl,dkT_dtemp,option)

  use Option_module

  implicit none

  class(kT_power_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature, porosity
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dkT_dsatl, dkT_dtemp
  type(option_type), intent(inout) :: option

  PetscReal :: tempreal, kT_base, shifted_temp, unused

  ! saturation behavior from default function
  call this%kT_default_type%CalculateTCond(liquid_saturation,0.d0,0.d0, &
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
  TCFCubicPolynomialCreate%alpha = 1.0d0
  TCFCubicPolynomialCreate%isotropic   = PETSC_TRUE
  TCFCubicPolynomialCreate%full_tensor = PETSC_FALSE
  TCFCubicPolynomialCreate%kT     = UNINITIALIZED_DOUBLE
  TCFCubicPolynomialCreate%kTf    = UNINITIALIZED_DOUBLE
  TCFCubicPolynomialCreate%kT_x   = UNINITIALIZED_DOUBLE
  TCFCubicPolynomialCreate%kT_y   = UNINITIALIZED_DOUBLE
  TCFCubicPolynomialCreate%kT_z   = UNINITIALIZED_DOUBLE
  TCFCubicPolynomialCreate%kT_xy  = UNINITIALIZED_DOUBLE
  TCFCubicPolynomialCreate%kT_xz  = UNINITIALIZED_DOUBLE
  TCFCubicPolynomialCreate%kT_yz  = UNINITIALIZED_DOUBLE

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
     porosity,thermal_conductivity,dkT_dsatl,dkT_dtemp,option)

  use Option_module

  implicit none

  class(kT_cubic_polynomial_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature, porosity
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dkT_dsatl, dkT_dtemp
  type(option_type), intent(inout) :: option

  PetscReal :: kT_base, shifted_temp, unused, tempreal

  ! saturation behavior from default function
  call this%kT_default_type%CalculateTCond(liquid_saturation,0.d0,0.d0, &
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
  TCFLinearResistivityCreate%alpha = 1.0d0
  TCFLinearResistivityCreate%isotropic   = PETSC_TRUE
  TCFLinearResistivityCreate%full_tensor = PETSC_FALSE
  TCFLinearResistivityCreate%kT     = UNINITIALIZED_DOUBLE
  TCFLinearResistivityCreate%kTf    = UNINITIALIZED_DOUBLE
  TCFLinearResistivityCreate%kT_x   = UNINITIALIZED_DOUBLE
  TCFLinearResistivityCreate%kT_y   = UNINITIALIZED_DOUBLE
  TCFLinearResistivityCreate%kT_z   = UNINITIALIZED_DOUBLE
  TCFLinearResistivityCreate%kT_xy  = UNINITIALIZED_DOUBLE
  TCFLinearResistivityCreate%kT_xz  = UNINITIALIZED_DOUBLE
  TCFLinearResistivityCreate%kT_yz  = UNINITIALIZED_DOUBLE
  TCFLinearResistivityCreate%porosity_effect = PETSC_FALSE
  TCFLinearResistivityCreate%ref_por = UNINITIALIZED_DOUBLE
  TCFLinearResistivityCreate%por_exp = UNINITIALIZED_DOUBLE
  TCFLinearResistivityCreate%b = [ UNINITIALIZED_DOUBLE, &
                                   UNINITIALIZED_DOUBLE]

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
  ! default: porosity_effect is false, all variables uninitialized
  if ((Initialized(this%ref_por)) .and. &
      (Initialized(this%por_exp)) .and. &
      (Initialized(this%b(1))) .and. &
      (Initialized(this%b(2)))) then
     ! all required variables are set: turn on porosity effect
     this%porosity_effect = PETSC_TRUE
  else if ((Initialized(this%ref_por)) .or. &
           (Initialized(this%por_exp)) .or. &
           (Initialized(this%b(1))) .or. &
           (Initialized(this%b(2)))) then
     ! some are initialized, but not all: error message
     if (Uninitialized(this%ref_por)) then
        option%io_buffer = UninitializedMessage('REFERENCE_POROSITY',string)
        call PrintErrMsg(option)
     endif
     if (Uninitialized(this%por_exp)) then
        option%io_buffer = UninitializedMessage('POROSITY_EXPONENT',string)
        call PrintErrMsg(option)
     endif
     if (Uninitialized(this%b(1)) .or. Uninitialized(this%b(2))) then
        option%io_buffer = UninitializedMessage( &
             'INITIAL_LINEAR_COEFFICIENTS b(1) + b(2)*T',string)
        call PrintErrMsg(option)
     endif
  endif

end subroutine TCFLinearResistivityVerify

! ************************************************************************** !

subroutine TCFLinearResistivityConductivity(this,liquid_saturation, &
     temperature,porosity,thermal_conductivity,dkT_dsatl,dkT_dtemp,option)

  use Option_module

  implicit none

  class(kT_linear_resistivity_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature, porosity
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dkT_dsatl, dkT_dtemp
  type(option_type), intent(inout) :: option

  PetscReal :: dkT_from_por(2)
  PetscReal :: relkT, shifted_temp, tempreal, unused, scaled_por, solid_term

  ! saturation behavior from default function
  call this%kT_default_type%CalculateTCond(liquid_saturation,0.d0,0.d0, &
       relkT,dkT_dsatl,unused,option)

  shifted_temp = temperature - this%ref_temp

  ! linear thermal resistivity idea from Birch & Clark (1940) and
  ! used by Blesch, Kulacki & Christensen (1983), ONWI-495
  ! kT = relkT/(a1 + a2*T)
  tempreal = this%a(1) + shifted_temp*this%a(2)
  thermal_conductivity = relkT / tempreal

  if (this%porosity_effect) then
     ! from GRS-281 (2012) Tab B.4 for granular salt (Salzgrus) with
     ! air-filled porosity at Gorleben (VSG)
     ! kT = (1 - por/ref_por)^por_exp * kT_intact + ...
     !      (por/ref_por) * (b(1) + b(2)*T)
     scaled_por = min(porosity / this%ref_por, 1.d0)

     solid_term = (1.d0 - scaled_por) ** this%por_exp
     thermal_conductivity = thermal_conductivity * solid_term
     thermal_conductivity = thermal_conductivity + scaled_por * &
          (this%b(1) + this%b(2)*shifted_temp)

     dkT_from_por(1) = solid_term
     dkT_from_por(2) = scaled_por * this%b(2)
  else
     dkT_from_por(1) = 1.d0
     dkT_from_por(2) = 0.d0
  end if

  dkT_dtemp = -this%a(2) * thermal_conductivity / tempreal * &
       dkT_from_por(1) + dkT_from_por(2)
  dkT_dsatl = dkT_dsatl * dkT_from_por(1) / tempreal

end subroutine TCFLinearResistivityConductivity

! ************************************************************************** !

function TCFFrozenCreate()

  implicit none

  class(kT_frozen_type), pointer :: TCFFrozenCreate

  allocate(TCFFrozenCreate)
  TCFFrozenCreate%kT_wet    = UNINITIALIZED_DOUBLE
  TCFFrozenCreate%kT_dry    = UNINITIALIZED_DOUBLE
  TCFFrozenCreate%kT_frozen = UNINITIALIZED_DOUBLE
  TCFFrozenCreate%alpha     = UNINITIALIZED_DOUBLE
  TCFFrozenCreate%alpha_fr  = UNINITIALIZED_DOUBLE
  TCFFrozenCreate%ice_model = UNINITIALIZED_INTEGER
  TCFFrozenCreate%isotropic   = PETSC_TRUE
  TCFFrozenCreate%full_tensor = PETSC_FALSE
  TCFFrozenCreate%kT     = UNINITIALIZED_DOUBLE
  TCFFrozenCreate%kTf    = UNINITIALIZED_DOUBLE
  TCFFrozenCreate%kT_x   = UNINITIALIZED_DOUBLE
  TCFFrozenCreate%kT_y   = UNINITIALIZED_DOUBLE
  TCFFrozenCreate%kT_z   = UNINITIALIZED_DOUBLE
  TCFFrozenCreate%kT_xy  = UNINITIALIZED_DOUBLE
  TCFFrozenCreate%kT_xz  = UNINITIALIZED_DOUBLE
  TCFFrozenCreate%kT_yz  = UNINITIALIZED_DOUBLE

end function TCFFrozenCreate

! ************************************************************************** !

subroutine TCFFrozenVerify(this,name,option)

  use Option_module

  implicit none

  class(kT_frozen_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'THERMAL_CONDUCTIVITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'THERMAL_CONDUCTIVITY_FUNCTION,FROZEN'
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
  if (Uninitialized(this%alpha)) then
    option%io_buffer = UninitializedMessage('EXPONENT',string)
    call PrintErrMsg(option)
  endif
  ! Freezing is optional, but related parameters must be initialized to use it
  if (Initialized(this%kT_frozen)) then
    if (Uninitialized(this%alpha_fr)) then
      option%io_buffer = UninitializedMessage('FROZEN EXPONENT (MUST BE '&
                                            //'SPECIFIED WITH FROZEN THERMAL '&
                                            //'CONDUCTIVITY)',string)
      call PrintErrMsg(option)
    elseif (Uninitialized(this%ice_model)) then
      option%io_buffer = UninitializedMessage('ICE MODEL (MUST BE '&
                                            //'SPECIFIED WITH FROZEN THERMAL '&
                                            //'CONDUCTIVITY)',string)
      call PrintErrMsg(option)
    else
      option%flow%th_freezing = PETSC_TRUE
      ! Outside of TH mode, frozen parameters aren't actually used
      if (.not. option%iflowmode == TH_MODE .and. &
          .not. option%iflowmode == TH_TS_MODE) then
        option%io_buffer = 'FREEZING MODEL ONLY UTILIZED IN TH MODE. ONLY ' &
                         //'NON-FROZEN PARAMETERS WILL BE EMPLOYED FOR ' &
                         //'THERMAL CONDUCTIVITY CALCULATION.'
        call PrintWrnMsg(option)
      endif
    endif
  endif

end subroutine TCFFrozenVerify

! ************************************************************************** !

subroutine TCFFrozenConductivity1(this,liquid_saturation,temperature, &
     porosity,thermal_conductivity,dkT_dsatl,dkT_dtemp,option)

  use Option_module

  implicit none

  class(kT_frozen_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature, porosity
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dkT_dsatl, dkT_dtemp
  type(option_type), intent(inout) :: option

  PetscReal, parameter :: epsilon = 1.d-6
  PetscReal :: Ke

  ! Soil Kersten numbers
  Ke = (liquid_saturation + epsilon)**(this%alpha) ! unfrozen

  ! Do not use freezing
  thermal_conductivity = this%kT_dry + (this%kT_wet - this%kT_dry)*Ke
  dkT_dtemp = 0.0d0
  dkT_dsatl = (this%kT_wet - this%kT_dry) * this%alpha * &
              liquid_saturation**(this%alpha - 1)

end subroutine TCFFrozenConductivity1

! ************************************************************************** !

subroutine TCFFrozenConductivity2(this,liquid_saturation,ice_saturation,   &
     temperature,porosity,thermal_conductivity, &
     dkT_dsatl,dkT_dsati,dkT_dtemp,option)

  use Option_module

  implicit none

  class(kT_frozen_type) :: this
  PetscReal, intent(in) :: liquid_saturation, ice_saturation, temperature
  PetscReal, intent(in) :: porosity
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dkT_dsatl, dkT_dsati, dkT_dtemp
  type(option_type), intent(inout) :: option

  PetscReal, parameter :: epsilon = 1.d-6
  PetscReal :: Ke, Ke_fr

  ! Soil Kersten numbers
  Ke = (liquid_saturation + epsilon)**(this%alpha)    ! unfrozen
  Ke_fr = (ice_saturation + epsilon)**(this%alpha_fr) ! frozen

  ! Use freezing
  thermal_conductivity = this%kT_wet*Ke + this%kT_frozen*Ke_fr + &
                         (1.d0 - Ke - Ke_fr)*this%kT_dry
  dkT_dtemp = 0.0d0
  dkT_dsatl = (this%kT_wet - this%kT_dry) * this%alpha * &
              liquid_saturation**(this%alpha - 1)
  dkT_dsati = (this%kT_frozen - this%kT_dry) * this%alpha_fr * &
              ice_saturation**(this%alpha_fr - 1)

end subroutine TCFFrozenConductivity2

! ************************************************************************** !

function TCFCompositeCreate()

  implicit none

  class(kT_composite_type), pointer :: TCFCompositeCreate

  allocate(TCFCompositeCreate)
  TCFCompositeCreate%isotropic   = PETSC_TRUE
  TCFCompositeCreate%full_tensor = PETSC_FALSE
  TCFCompositeCreate%kT_wet = 0.0d0 ! not used, deferred to sub-functions
  TCFCompositeCreate%kT_dry = 0.0d0 ! not used, deferred to sub-functions
  TCFCompositeCreate%alpha  = 1.0d0 ! not used, deferred to sub-functions
  TCFCompositeCreate%kT     = UNINITIALIZED_DOUBLE
  TCFCompositeCreate%kTf    = UNINITIALIZED_DOUBLE
  TCFCompositeCreate%kT_x   = UNINITIALIZED_DOUBLE
  TCFCompositeCreate%kT_y   = UNINITIALIZED_DOUBLE
  TCFCompositeCreate%kT_z   = UNINITIALIZED_DOUBLE
  TCFCompositeCreate%kT_xy  = UNINITIALIZED_DOUBLE
  TCFCompositeCreate%kT_xz  = UNINITIALIZED_DOUBLE
  TCFCompositeCreate%kT_yz  = UNINITIALIZED_DOUBLE

  TCFCompositeCreate%dist   = UNINITIALIZED_DOUBLE
  TCFCompositeCreate%lkT = [ "UNINITIALIZED", &
                             "UNINITIALIZED", &
                             "UNINITIALIZED" ]

end function TCFCompositeCreate

! ************************************************************************** !

subroutine TCFCompositeVerify(this,name,option)

  use Option_module

  implicit none

  class(kT_composite_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'THERMAL_CONDUCTIVITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'THERMAL_CONDUCTIVITY_FUNCTION,COMPOSITE'
  endif
  call TCFBaseVerify(this,string,option)
  if (this%lkT(1) == "UNINITIALIZED") then
    option%io_buffer = UninitializedMessage( &
          'THERMAL CONDUCTIVITY SUB-FUNCTION FOR X',string)
    call PrintErrMsg(option)
  endif
  if (this%lkT(2) == "UNINITIALIZED") then
    option%io_buffer = UninitializedMessage( &
          'THERMAL CONDUCTIVITY SUB-FUNCTION FOR Y',string)
    call PrintErrMsg(option)
  endif
  if (this%lkT(3) == "UNINITIALIZED") then
    option%io_buffer = UninitializedMessage( &
          'THERMAL CONDUCTIVITY SUB-FUNCTION FOR Z',string)
    call PrintErrMsg(option)
  endif

end subroutine TCFCompositeVerify

! ************************************************************************** !

subroutine TCFCompositeTest(this,tcc_name,option)

  use Option_module

  implicit none

  class(kT_composite_type) :: this
  character(len=MAXWORDLENGTH) :: tcc_name
  type(option_type), intent(inout) :: option

  return

end subroutine TCFCompositeTest

! ************************************************************************** !

subroutine TCFCompositeConductivity(this,liquid_saturation, &
     temperature,porosity,thermal_conductivity,dkT_dsatl,dkT_dtemp,option)

  use Option_module

  implicit none

  class(kT_composite_type) :: this
  PetscReal, intent(in) :: liquid_saturation, temperature, porosity
  PetscReal, intent(out) :: thermal_conductivity
  PetscReal, intent(out) :: dkT_dsatl, dkT_dtemp
  type(option_type), intent(inout) :: option

  PetscInt  :: i
  PetscReal :: dist(-1:3)       ! unit vectors
  PetscReal :: mag              ! magnitude
  PetscReal :: scaling(3)       ! scaling factor
  PetscReal :: tmpkT(3)         ! thermal conductivity along certain axis
  PetscReal :: tmpdkT_dsatl(3)  ! differential conductivity wrt saturation
  PetscReal :: tmpdkT_dtemp(3)  ! differential conductivity wrt temperature

  dist = this%dist  ! get unit vectors

  thermal_conductivity = 0.0d0
  dkT_dtemp = 0.0d0
  dkT_dsatl = 0.0d0

  tmpkT = 0.0d0
  tmpdkT_dsatl = 0.0d0
  tmpdkT_dtemp = 0.0d0
  scaling = 1.0d0

  mag = sqrt( dist(1)**2 + dist(2)**2 + dist(3)**2 )

  do i = 1, 3 ! iterate over X, Y, and Z
    select type(tcf => this%ckT(i)%ptr%thermal_conductivity_function)
    class is(thermal_conductivity_base_type)
      scaling(i) = dabs(dist(i))**2/mag
      call tcf%CalculateTCond(liquid_saturation,temperature,0.d0,tmpkT(i), &
                              tmpdkT_dsatl(i),tmpdkT_dtemp(i),option)

      thermal_conductivity = thermal_conductivity + &
                             scaling(i)*tmpkT(i)*this%kTf(i,i)

      dkT_dtemp = dkT_dtemp + scaling(i)*tmpdkT_dsatl(i)*this%kTf(i,i)
      dkT_dsatl = dkT_dsatl + scaling(i)*tmpdkT_dtemp(i)*this%kTf(i,i)
    end select
  enddo

end subroutine TCFCompositeConductivity

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
      case('LINEAR')
        this%thermal_conductivity_function => TCFLinearCreate()
      case('POWER')
        this%thermal_conductivity_function => TCFPowerCreate()
      case('CUBIC_POLYNOMIAL')
        this%thermal_conductivity_function => TCFCubicPolynomialCreate()
      case('LINEAR_RESISTIVITY')
        this%thermal_conductivity_function => TCFLinearResistivityCreate()
      case('FROZEN')
        this%thermal_conductivity_function => TCFFrozenCreate()
      case('ASM_DRY')
        this%thermal_conductivity_function => TCFASMDryCreate()
      case('ASM_WATER_FILLED')
        this%thermal_conductivity_function => TCFASMWaterFilledCreate()
      case('ASM_RADIAL')
        this%thermal_conductivity_function => TCFASMRadialCreate()
      case('ASM_AXIAL')
        this%thermal_conductivity_function => TCFASMAxialCreate()
      case('COMPOSITE')
        this%thermal_conductivity_function => TCFCompositeCreate()
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
                            kwet,kdry,alpha,option)

  use Option_module

  implicit none

  class(thermal_conductivity_base_type) :: thermal_conductivity_function
  PetscReal :: kwet,kdry,alpha
  type(option_type) :: option

  select type(tcf => thermal_conductivity_function)
      !------------------------------------------
    class is(kT_default_type)
      tcf%kT_dry = kdry
      tcf%kT_wet = kwet
      tcf%alpha  = alpha
  end select

end subroutine TCFAssignDefault

! ************************************************************************** !

subroutine TCFAssignFrozen(thermal_conductivity_function,&
                           kwet,kdry,kfrozen,alpha,alpha_fr,icemod,option)

  use Option_module

  implicit none

  class(thermal_conductivity_base_type) :: thermal_conductivity_function
  PetscReal :: kwet,kdry,kfrozen,alpha,alpha_fr
  PetscInt :: icemod
  type(option_type) :: option

  select type(tcf => thermal_conductivity_function)
      !------------------------------------------
    class is(kT_frozen_type)
      tcf%kT_dry    = kdry
      tcf%kT_wet    = kwet
      tcf%kT_frozen = kfrozen
      tcf%alpha     = alpha
      tcf%alpha_fr  = alpha_fr
      tcf%ice_model = icemod
  end select

end subroutine TCFAssignFrozen

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
  class is(kT_linear_type)
    error_string = trim(error_string) // 'LINEAR'
  class is(kT_power_type)
    error_string = trim(error_string) // 'POWER'
  class is(kT_cubic_polynomial_type)
    error_string = trim(error_string) // 'CUBIC_POLYNOMIAL'
  class is(kT_linear_resistivity_type)
    error_string = trim(error_string) // 'LINEAR_RESISTIVITY'
  class is(kT_frozen_type)
    error_string = trim(error_string) // 'FROZEN'
  class is(kT_ASM_dry_type)
    error_string = trim(error_string) // 'ASM_DRY'
  class is(kT_ASM_water_filled_type)
    error_string = trim(error_string) // 'ASM_WATER_FILLED'
  class is(kT_ASM_radial_type)
    error_string = trim(error_string) // 'ASM_RADIAL'
  class is(kT_ASM_axial_type)
    error_string = trim(error_string) // 'ASM_AXIAL'
  class is(kT_composite_type)
    error_string = trim(error_string) // 'COMPOSITE'
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
      case('KERSTEN_EXPONENT')
        call InputReadDouble(input,option,tcf%alpha)
        call InputErrorMsg(input,option,'Kersten exponent', &
             error_string)
      case default
        call InputKeywordUnrecognized(input,keyword, &
             'constant thermal conductivity',option)
      end select
      !------------------------------------------
    class is(kT_default_type)
      call TCFDefaultRead(tcf,input,keyword,error_string,'default',option)
      !------------------------------------------
    class is(kT_linear_type)
      call TCFDefaultRead(tcf,input,keyword,error_string,'linear',option)
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
    class is(kT_ASM_dry_type)
      select case(keyword)
      case('ASM_DRY_COEFFICIENT')
        call InputReadDouble(input,option,tcf%dry_T_coeff)
        call InputErrorMsg(input,option, &
             'assembly dry conditions thermal conductivity temperature '// &
             'coefficient',error_string)
      case('ASM_DRY_EXPONENT')
        call InputReadDouble(input,option,tcf%dry_T_power)
        call InputErrorMsg(input,option, &
             'assembly dry conditions thermal conductivity temperature '// &
             'exponent',error_string)
      case default
        call TCFDefaultRead(tcf,input,keyword,error_string, &
                            'assembly dry conditions',option)
      end select
      !------------------------------------------
    class is(kT_ASM_water_filled_type)
      select case(keyword)
      case('THERMAL_CONDUCTIVITY_WATER')
        call InputReadDouble(input,option,tcf%kT_water)
        call InputErrorMsg(input,option,'thermal conductivity water', &
             error_string)
        call InputReadAndConvertUnits(input,tcf%kT_water,'W/m-C', &
             'CHARACTERISTIC_CURVES_THERMAL,thermal conductivity water',option)
      case('THERMAL_CONDUCTIVITY_SOLID')
        call InputReadDouble(input,option,tcf%kT_solid)
        call InputErrorMsg(input,option,'thermal conductivity solid', &
             error_string)
        call InputReadAndConvertUnits(input,tcf%kT_solid,'W/m-C', &
             'CHARACTERISTIC_CURVES_THERMAL,thermal conductivity solid',option)
      case('POROSITY_ASSEMBLY')
        call InputReadDouble(input,option,tcf%porosity_asm)
        call InputErrorMsg(input,option, &
             'assembly porosity', &
             error_string)
      case default
        call TCFDefaultRead(tcf,input,keyword,error_string, &
                            'assembly water-filled conditions',option)
      end select
      !------------------------------------------
    class is(kT_ASM_radial_type)
      select case(keyword)
      case('ASM_DRY_COEFFICIENT')
        call InputReadDouble(input,option,tcf%dry_T_coeff)
        call InputErrorMsg(input,option, &
             'assembly dry conditions thermal conductivity temperature '// &
             'coefficient',error_string)
      case('ASM_DRY_EXPONENT')
        call InputReadDouble(input,option,tcf%dry_T_power)
        call InputErrorMsg(input,option, &
             'assembly dry conditions thermal conductivity temperature '// &
             'exponent',error_string)
      case('THERMAL_CONDUCTIVITY_WATER')
        call InputReadDouble(input,option,tcf%kT_water)
        call InputErrorMsg(input,option,'thermal conductivity water', &
             error_string)
        call InputReadAndConvertUnits(input,tcf%kT_water,'W/m-C', &
             'CHARACTERISTIC_CURVES_THERMAL,thermal conductivity water',option)
      case('THERMAL_CONDUCTIVITY_SOLID')
        call InputReadDouble(input,option,tcf%kT_solid)
        call InputErrorMsg(input,option,'thermal conductivity solid', &
             error_string)
        call InputReadAndConvertUnits(input,tcf%kT_solid,'W/m-C', &
             'CHARACTERISTIC_CURVES_THERMAL,thermal conductivity solid',option)
      case('POROSITY_ASSEMBLY')
        call InputReadDouble(input,option,tcf%porosity_asm)
        call InputErrorMsg(input,option, &
             'assembly porosity', &
             error_string)
      case default
        call TCFDefaultRead(tcf,input,keyword,error_string,&
                            'assembly radial',option)
      end select
      !------------------------------------------
    class is(kT_ASM_axial_type)
      select case(keyword)
      case('THERMAL_CONDUCTIVITY_WATER')
        call InputReadDouble(input,option,tcf%kT_water)
        call InputErrorMsg(input,option,'thermal conductivity water', &
             error_string)
        call InputReadAndConvertUnits(input,tcf%kT_water,'W/m-C', &
             'CHARACTERISTIC_CURVES_THERMAL,thermal conductivity water',option)
      case('THERMAL_CONDUCTIVITY_SOLID')
        call InputReadDouble(input,option,tcf%kT_solid)
        call InputErrorMsg(input,option,'thermal conductivity solid', &
             error_string)
        call InputReadAndConvertUnits(input,tcf%kT_solid,'W/m-C', &
             'CHARACTERISTIC_CURVES_THERMAL,thermal conductivity solid',option)
      case('POROSITY_ASSEMBLY')
        call InputReadDouble(input,option,tcf%porosity_asm)
        call InputErrorMsg(input,option, &
             'assembly porosity', &
             error_string)
      case default
        call TCFDefaultRead(tcf,input,keyword,error_string,&
                            'assembly axial',option)
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
      case('REFERENCE_POROSITY')
        call InputReadDouble(input,option,tcf%ref_por)
        call InputErrorMsg(input,option, 'reference porosity', &
             error_string)
      case('POROSITY_EXPONENT')
        call InputReadDouble(input,option,tcf%por_exp)
        call InputErrorMsg(input,option, 'porosity exponent', &
             error_string)
      case('INITIAL_LINEAR_COEFFICIENTS')
        call InputReadNDoubles(input,option,tcf%b,2)
        call InputErrorMsg(input,option, 'initial thermal conductivity linear coefficients', &
             error_string)
      case default
        call TCFDefaultRead(tcf,input,keyword,error_string, &
          'linear resistivity',option)
      end select
      !------------------------------------------
    class is(kT_frozen_type)
      select case(keyword)
      case('THERMAL_CONDUCTIVITY_FROZEN')
        call InputReadDouble(input,option,tcf%kT_frozen)
        call InputErrorMsg(input,option,'thermal conductivity frozen', &
             error_string)
        call InputReadAndConvertUnits(input,tcf%kT_frozen,'W/m-C', &
             'CHARACTERISTIC_CURVES_THERMAL,thermal conductivity frozen', &
             option)
      case('KERSTEN_EXPONENT_FROZEN')
        call InputReadDouble(input,option,tcf%alpha_fr)
        call InputErrorMsg(input,option,'Kersten exponent - frozen', &
             error_string)
      case('ICE_MODEL')
        call InputReadCard(input,option,keyword,PETSC_FALSE)
        call StringToUpper(keyword)
        select case (trim(keyword))
        case ('PAINTER_EXPLICIT')
          tcf%ice_model = PAINTER_EXPLICIT
        case ('PAINTER_KARRA_IMPLICIT')
          tcf%ice_model = PAINTER_KARRA_IMPLICIT
        case ('PAINTER_KARRA_EXPLICIT')
          tcf%ice_model = PAINTER_KARRA_EXPLICIT
        case ('PAINTER_KARRA_EXPLICIT_NOCRYO')
          tcf%ice_model = PAINTER_KARRA_EXPLICIT_NOCRYO
        case ('DALL_AMICO')
          tcf%ice_model = DALL_AMICO
        case default
          option%io_buffer = 'Cannot identify the specificed ice model. &
           &Specify PAINTER_EXPLICIT, PAINTER_KARRA_IMPLICIT, &
           &PAINTER_KARRA_EXPLICIT, PAINTER_KARRA_EXPLICIT_NOCRYO, &
           &or DALL_AMICO.'
          call PrintErrMsg(option)
        end select
      case default
        call TCFDefaultRead(tcf,input,keyword,error_string, &
          'frozen',option)
      end select
      !------------------------------------------
    class is(kT_composite_type)
      select case(keyword)
      case('COMPOSITE_X')
        call InputReadWord(input,option,tcf%lkT(1),PETSC_TRUE)
        call InputErrorMsg(input,option,'fn_x','COMPOSITE')
      case('COMPOSITE_Y')
        call InputReadWord(input,option,tcf%lkT(2),PETSC_TRUE)
        call InputErrorMsg(input,option,'fn_y','COMPOSITE')
      case('COMPOSITE_Z')
        call InputReadWord(input,option,tcf%lkT(3),PETSC_TRUE)
        call InputErrorMsg(input,option,'fn_z','COMPOSITE')
      case default
        call TCFDefaultRead(tcf,input,keyword,error_string, &
          'composite',option)
      end select
      !------------------------------------------
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
  case('KERSTEN_EXPONENT')
     call InputReadDouble(input,option,tcf%alpha)
     call InputErrorMsg(input,option,'Kersten exponent', &
          error_string)
  case('ANISOTROPY_RATIO_X')
    call InputReadDouble(input,option,tcf%kT_x)
    call InputErrorMsg(input,option, &
       'anisotropic thermal conductivity X component', error_string)
  case('ANISOTROPY_RATIO_Y')
    call InputReadDouble(input,option,tcf%kT_y)
    call InputErrorMsg(input,option, &
       'anisotropic thermal conductivity Y component', error_string)
  case('ANISOTROPY_RATIO_Z')
    call InputReadDouble(input,option,tcf%kT_z)
    call InputErrorMsg(input,option, &
       'anisotropic thermal conductivity Z component', error_string)
  case('ANISOTROPY_RATIO_XY')
    call InputReadDouble(input,option,tcf%kT_xy)
    call InputErrorMsg(input,option, &
       'anisotropic thermal conductivity XY component', error_string)
  case('ANISOTROPY_RATIO_XZ')
    call InputReadDouble(input,option,tcf%kT_xz)
    call InputErrorMsg(input,option, &
       'anisotropic thermal conductivity XZ component', error_string)
  case('ANISOTROPY_RATIO_YZ')
    call InputReadDouble(input,option,tcf%kT_yz)
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

  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'THERMAL_CHARACTERISTIC_CURVES, THERMAL_CONDUCTIVITY_FUNCTION,'

  select type(tcf => thermal_conductivity_function)
  !------------------------------------------
  class is(kT_default_type)
    error_string = trim(error_string) // ' ANISOTROPIC DEFAULT TYPE'

    ! check if wet and dry thermal conductivities are initialized
    if (.not. Initialized(tcf%kT_dry) .or. &
        .not. Initialized(tcf%kT_wet)) then
      ! wet and dry values must be specified per anisotropic component
      option%io_buffer = 'Must specify wet and dry thermal conductivity ' &
                       //'values in order to use anisotropy ratios in ' &
                       // trim(error_string) // '.'
      call PrintErrMsg(option)
    elseif (Initialized(tcf%kT_x) .or. Initialized(tcf%kT_y) &
            .or. Initialized(tcf%kT_z)) then
      ! inputs must be anisotropy ratios between zero and one

      ! check diagonal components first, as tensor must at least be diagonal
      if (Initialized(tcf%kT_x)) then
        if (tcf%kT_x < 0.0d0 .or. tcf%kT_x > 1.0d0) then
          option%io_buffer = 'Anisotropy ratio for X must lie between 0 and ' &
                           //'1 in '// trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        tcf%kTf(1,1) = tcf%kT_x
        tcf%kT(1,1,1) = tcf%kT_dry * tcf%kT_x
        tcf%kT(1,1,2) = tcf%kT_wet * tcf%kT_x
      else
        option%io_buffer = 'Anisotropy ratio for X uninitialized in ' &
                           // trim(error_string) // '.'
        call PrintErrMsg(option)
      endif

      if (Initialized(tcf%kT_y)) then
        if (tcf%kT_y < 0.0d0 .or. tcf%kT_y > 1.0d0) then
          option%io_buffer = 'Anisotropy ratio for Y must lie between 0 and ' &
                           //'1 in '// trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        tcf%kTf(2,2) = tcf%kT_y
        tcf%kT(2,2,1) = tcf%kT_dry * tcf%kT_y
        tcf%kT(2,2,2) = tcf%kT_wet * tcf%kT_y
      else
        option%io_buffer = 'Anisotropy ratio for Y uninitialized in ' &
                           // trim(error_string) // '.'
        call PrintErrMsg(option)
      endif

      if (Initialized(tcf%kT_z)) then
        if (tcf%kT_z < 0.0d0 .or. tcf%kT_z > 1.0d0) then
          option%io_buffer = 'Anisotropy ratio for Z must lie between 0 and ' &
                           //'1 in '// trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        tcf%kTf(3,3) = tcf%kT_z
        tcf%kT(3,3,1) = tcf%kT_dry * tcf%kT_z
        tcf%kT(3,3,2) = tcf%kT_wet * tcf%kT_z
      else
        option%io_buffer = 'Anisotropy ratio for Z uninitialized in ' &
                           // trim(error_string) // '.'
        call PrintErrMsg(option)
      endif

      ! check off-diagonal components next; if one is given, so must the others
      if (Initialized(tcf%kT_xy)) then
        if (tcf%kT_xy < 0.0d0 .or. tcf%kT_xy > 1.0d0) then
          option%io_buffer = 'Anisotropy ratio for XY must lie between 0 and ' &
                           //'1 in '// trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        if (.not. Initialized(tcf%kT_xz) .or. &
            .not. Initialized(tcf%kT_yz)) then
          option%io_buffer = 'All off-diagonal components must be specified ' &
                          // 'if XY ratio is provided in ' &
                          // trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        tcf%isotropic = PETSC_FALSE
        tcf%kTf(1,2) = tcf%kT_xy
        tcf%kTf(2,1) = tcf%kT_xy
        tcf%kT(1,2,1) = tcf%kT_dry * tcf%kT_xy
        tcf%kT(2,1,2) = tcf%kT_wet * tcf%kT_xy
        tcf%kT(1,2,1) = tcf%kT_dry * tcf%kT_xy
        tcf%kT(2,1,2) = tcf%kT_wet * tcf%kT_xy
      endif

      if (Initialized(tcf%kT_xz)) then
        if (tcf%kT_xz < 0.0d0 .or. tcf%kT_xz > 1.0d0) then
          option%io_buffer = 'Anisotropy ratio for XZ must lie between 0 and ' &
                           //'1 in '// trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        if (.not. Initialized(tcf%kT_xy) .or. &
            .not. Initialized(tcf%kT_yz)) then
          option%io_buffer = 'All off-diagonal components must be specified ' &
                          // 'if XZ ratio is provided in ' &
                          // trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        tcf%isotropic = PETSC_FALSE
        tcf%kTf(1,3) = tcf%kT_xz
        tcf%kTf(3,1) = tcf%kT_xz
        tcf%kT(1,3,1) = tcf%kT_dry * tcf%kT_xz
        tcf%kT(1,3,2) = tcf%kT_wet * tcf%kT_xz
        tcf%kT(3,1,1) = tcf%kT_dry * tcf%kT_xz
        tcf%kT(3,1,2) = tcf%kT_wet * tcf%kT_xz
      endif

      if (Initialized(tcf%kT_yz)) then
        if (tcf%kT_yz < 0.0d0 .or. tcf%kT_yz > 1.0d0) then
          option%io_buffer = 'Anisotropy ratio for YZ must lie between 0 and ' &
                           //'1 in '// trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        if (.not. Initialized(tcf%kT_xy) .or. &
            .not. Initialized(tcf%kT_xz)) then
          option%io_buffer = 'All off-diagonal components must be specified ' &
                          // 'if YZ ratio is provided in ' &
                          // trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        tcf%isotropic = PETSC_FALSE
        tcf%kTf(2,3) = tcf%kT_yz
        tcf%kTf(3,2) = tcf%kT_yz
        tcf%kT(2,3,1) = tcf%kT_dry * tcf%kT_yz
        tcf%kT(2,3,2) = tcf%kT_wet * tcf%kT_yz
        tcf%kT(3,2,1) = tcf%kT_dry * tcf%kT_yz
        tcf%kT(3,2,2) = tcf%kT_wet * tcf%kT_yz
      endif

      ! check for isotropy and fully initialize tensor
      if (tcf%kT_x == tcf%kT_y .and. tcf%kT_y == tcf%kT_z) then
        if (Initialized(tcf%kT_xy) .or. Initialized(tcf%kT_xz) &
            .or. Initialized(tcf%kT_yz)) then
          tcf%isotropic = PETSC_FALSE
          tcf%full_tensor = PETSC_TRUE
        else
          tcf%isotropic = PETSC_TRUE
          tcf%kT(:,:,:) = 0.0d0
          tcf%kT(1,1,1) = tcf%kT_dry
          tcf%kT(2,2,1) = tcf%kT_dry
          tcf%kT(3,3,1) = tcf%kT_dry
          tcf%kT(1,1,2) = tcf%kT_wet
          tcf%kT(2,2,2) = tcf%kT_wet
          tcf%kT(3,3,2) = tcf%kT_wet
          option%io_buffer = 'Thermal conductivity will be treated as' &
                          // ' isotropic in ' &
                          // trim(error_string) // '.'
          call PrintMsg(option)
        endif
      else
        tcf%isotropic = PETSC_FALSE
        if (Initialized(tcf%kT_xy) &
            .or. Initialized(tcf%kT_xz) &
            .or. Initialized(tcf%kT_yz)) then
          ! full thermal conductivity tensor
          tcf%full_tensor = PETSC_TRUE
        else
          ! diagonal thermal conductivity tensor
          tcf%kT(1,2,:) = 0.0d0
          tcf%kT(1,3,:) = 0.0d0
          tcf%kT(2,3,:) = 0.0d0
          tcf%kT(2,1,:) = 0.0d0
          tcf%kT(3,1,:) = 0.0d0
          tcf%kT(3,2,:) = 0.0d0
          tcf%kTf(1,2) = 0.0d0
          tcf%kTf(1,3) = 0.0d0
          tcf%kTf(2,3) = 0.0d0
          tcf%kTf(2,1) = 0.0d0
          tcf%kTf(3,1) = 0.0d0
          tcf%kTf(3,2) = 0.0d0
        endif
      endif

    elseif (Initialized(tcf%kT_xy) &
            .or. Initialized(tcf%kT_xz) &
            .or. Initialized(tcf%kT_yz)) then
      if (.not. Initialized(tcf%kT_x) .or. &
          .not. Initialized(tcf%kT_y) .or. &
          .not. Initialized(tcf%kT_z)) then
        option%io_buffer = 'Diagonal components of thermal conductivity ' &
                        // 'must be specified if off-diagonal components are ' &
                        // 'provided in '// trim(error_string) // '.'
        call PrintErrMsg(option)
      endif
    else
      tcf%isotropic = PETSC_TRUE
      tcf%kT(:,:,:) = 0.0d0
      tcf%kT(1,1,1) = tcf%kT_dry
      tcf%kT(2,2,1) = tcf%kT_dry
      tcf%kT(3,3,1) = tcf%kT_dry
      tcf%kT(1,1,2) = tcf%kT_wet
      tcf%kT(2,2,2) = tcf%kT_wet
      tcf%kT(3,3,2) = tcf%kT_wet
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
  !------------------------------------------
  class is(kT_frozen_type)
    error_string = trim(error_string) // ' ANISOTROPIC FROZEN TYPE'

    ! check if wet and dry thermal conductivities are initialized
    if (.not. Initialized(tcf%kT_dry) .or. &
        .not. Initialized(tcf%kT_wet)) then
      ! wet and dry values must be specified per anisotropic component
      option%io_buffer = 'Must specify wet and dry thermal conductivity ' &
                       //'values in order to use anisotropy ratios in ' &
                       // trim(error_string) // '.'
      call PrintErrMsg(option)
    elseif (Initialized(tcf%kT_x) .or. Initialized(tcf%kT_y) &
            .or. Initialized(tcf%kT_z)) then
      ! inputs must be anisotropy ratios between zero and one

      ! check diagonal components first, as tensor must at least be diagonal
      if (Initialized(tcf%kT_x)) then
        if (tcf%kT_x < 0.0d0 .or. tcf%kT_x > 1.0d0) then
          option%io_buffer = 'Anisotropy ratio for X must lie between 0 and ' &
                           //'1 in '// trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        tcf%kTf(1,1) = tcf%kT_x
        tcf%kT(1,1,1) = tcf%kT_dry * tcf%kT_x
        tcf%kT(1,1,2) = tcf%kT_wet * tcf%kT_x
        if (Initialized(tcf%kT_frozen)) then
          tcf%kT(1,1,3) = tcf%kT_frozen * tcf%kT_x
        endif
      else
        option%io_buffer = 'Anisotropy ratio for X uninitialized in ' &
                           // trim(error_string) // '.'
        call PrintErrMsg(option)
      endif

      if (Initialized(tcf%kT_y)) then
        if (tcf%kT_y < 0.0d0 .or. tcf%kT_y > 1.0d0) then
          option%io_buffer = 'Anisotropy ratio for Y must lie between 0 and ' &
                           //'1 in '// trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        tcf%kTf(2,2) = tcf%kT_y
        tcf%kT(2,2,1) = tcf%kT_dry * tcf%kT_y
        tcf%kT(2,2,2) = tcf%kT_wet * tcf%kT_y
        if (Initialized(tcf%kT_frozen)) then
          tcf%kT(2,2,3) = tcf%kT_frozen * tcf%kT_y
        endif
      else
        option%io_buffer = 'Anisotropy ratio for Y uninitialized in ' &
                           // trim(error_string) // '.'
        call PrintErrMsg(option)
      endif

      if (Initialized(tcf%kT_z)) then
        if (tcf%kT_z < 0.0d0 .or. tcf%kT_z > 1.0d0) then
          option%io_buffer = 'Anisotropy ratio for Z must lie between 0 and ' &
                           //'1 in '// trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        tcf%kTf(3,3) = tcf%kT_z
        tcf%kT(3,3,1) = tcf%kT_dry * tcf%kT_z
        tcf%kT(3,3,2) = tcf%kT_wet * tcf%kT_z
        if (Initialized(tcf%kT_frozen)) then
          tcf%kT(3,3,3) = tcf%kT_frozen * tcf%kT_z
        endif
      else
        option%io_buffer = 'Anisotropy ratio for Z uninitialized in ' &
                           // trim(error_string) // '.'
        call PrintErrMsg(option)
      endif

      ! check off-diagonal components next; if one is given, so must the others
      if (Initialized(tcf%kT_xy)) then
        if (tcf%kT_xy < 0.0d0 .or. tcf%kT_xy > 1.0d0) then
          option%io_buffer = 'Anisotropy ratio for XY must lie between 0 and ' &
                           //'1 in '// trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        if (.not. Initialized(tcf%kT_xz) .or. &
            .not. Initialized(tcf%kT_yz)) then
          option%io_buffer = 'All off-diagonal components must be specified ' &
                          // 'if XY ratio is provided in ' &
                          // trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        tcf%isotropic = PETSC_FALSE
        tcf%kTf(1,2) = tcf%kT_xy
        tcf%kTf(2,1) = tcf%kT_xy
        tcf%kT(1,2,1) = tcf%kT_dry * tcf%kT_xy
        tcf%kT(2,1,2) = tcf%kT_wet * tcf%kT_xy
        tcf%kT(1,2,1) = tcf%kT_dry * tcf%kT_xy
        tcf%kT(2,1,2) = tcf%kT_wet * tcf%kT_xy
        if (Initialized(tcf%kT_frozen)) then
          tcf%kT(1,2,3) = tcf%kT_frozen * tcf%kT_xy
          tcf%kT(2,1,3) = tcf%kT_frozen * tcf%kT_xy
        endif
      endif

      if (Initialized(tcf%kT_xz)) then
        if (tcf%kT_xz < 0.0d0 .or. tcf%kT_xz > 1.0d0) then
          option%io_buffer = 'Anisotropy ratio for XZ must lie between 0 and ' &
                           //'1 in '// trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        if (.not. Initialized(tcf%kT_xy) .or. &
            .not. Initialized(tcf%kT_yz)) then
          option%io_buffer = 'All off-diagonal components must be specified ' &
                          // 'if XZ ratio is provided in ' &
                          // trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        tcf%isotropic = PETSC_FALSE
        tcf%kTf(1,3) = tcf%kT_xz
        tcf%kTf(3,1) = tcf%kT_xz
        tcf%kT(1,3,1) = tcf%kT_dry * tcf%kT_xz
        tcf%kT(1,3,2) = tcf%kT_wet * tcf%kT_xz
        tcf%kT(3,1,1) = tcf%kT_dry * tcf%kT_xz
        tcf%kT(3,1,2) = tcf%kT_wet * tcf%kT_xz
        if (Initialized(tcf%kT_frozen)) then
          tcf%kT(1,3,3) = tcf%kT_frozen * tcf%kT_xz
          tcf%kT(3,1,3) = tcf%kT_frozen * tcf%kT_xz
        endif
      endif

      if (Initialized(tcf%kT_yz)) then
        if (tcf%kT_yz < 0.0d0 .or. tcf%kT_yz > 1.0d0) then
          option%io_buffer = 'Anisotropy ratio for YZ must lie between 0 and ' &
                           //'1 in '// trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        if (.not. Initialized(tcf%kT_xy) .or. &
            .not. Initialized(tcf%kT_xz)) then
          option%io_buffer = 'All off-diagonal components must be specified ' &
                          // 'if YZ ratio is provided in ' &
                          // trim(error_string) // '.'
          call PrintErrMsg(option)
        endif
        tcf%isotropic = PETSC_FALSE
        tcf%kTf(2,3) = tcf%kT_yz
        tcf%kTf(3,2) = tcf%kT_yz
        tcf%kT(2,3,1) = tcf%kT_dry * tcf%kT_yz
        tcf%kT(2,3,2) = tcf%kT_wet * tcf%kT_yz
        tcf%kT(3,2,1) = tcf%kT_dry * tcf%kT_yz
        tcf%kT(3,2,2) = tcf%kT_wet * tcf%kT_yz
        if (Initialized(tcf%kT_frozen)) then
          tcf%kT(2,3,3) = tcf%kT_frozen * tcf%kT_yz
          tcf%kT(3,2,3) = tcf%kT_frozen * tcf%kT_yz
        endif
      endif

      ! check for isotropy and fully initialize tensor
      if (tcf%kT_x == tcf%kT_y .and. tcf%kT_y == tcf%kT_z) then
        if (Initialized(tcf%kT_xy) .or. Initialized(tcf%kT_xz) &
            .or. Initialized(tcf%kT_yz)) then
          tcf%isotropic = PETSC_FALSE
          tcf%full_tensor = PETSC_TRUE
        else
          tcf%isotropic = PETSC_TRUE
          tcf%kT(:,:,:) = 0.0d0
          tcf%kT(1,1,1) = tcf%kT_dry
          tcf%kT(2,2,1) = tcf%kT_dry
          tcf%kT(3,3,1) = tcf%kT_dry
          tcf%kT(1,1,2) = tcf%kT_wet
          tcf%kT(2,2,2) = tcf%kT_wet
          tcf%kT(3,3,2) = tcf%kT_wet
          tcf%kT(1,1,3) = tcf%kT_frozen
          tcf%kT(2,2,3) = tcf%kT_frozen
          tcf%kT(3,3,3) = tcf%kT_frozen
          option%io_buffer = 'Thermal conductivity will be treated as' &
                          // ' isotropic in ' &
                          // trim(error_string) // '.'
          call PrintMsg(option)
        endif
      else
        tcf%isotropic = PETSC_FALSE
        if (Initialized(tcf%kT_xy) &
            .or. Initialized(tcf%kT_xz) &
            .or. Initialized(tcf%kT_yz)) then
          ! full thermal conductivity tensor
          tcf%full_tensor = PETSC_TRUE
        else
          ! diagonal thermal conductivity tensor
          tcf%kT(1,2,:) = 0.0d0
          tcf%kT(1,3,:) = 0.0d0
          tcf%kT(2,3,:) = 0.0d0
          tcf%kT(2,1,:) = 0.0d0
          tcf%kT(3,1,:) = 0.0d0
          tcf%kT(3,2,:) = 0.0d0
          tcf%kTf(1,2) = 0.0d0
          tcf%kTf(1,3) = 0.0d0
          tcf%kTf(2,3) = 0.0d0
          tcf%kTf(2,1) = 0.0d0
          tcf%kTf(3,1) = 0.0d0
          tcf%kTf(3,2) = 0.0d0
        endif
      endif

    elseif (Initialized(tcf%kT_xy) &
            .or. Initialized(tcf%kT_xz) &
            .or. Initialized(tcf%kT_yz)) then
      if (.not. Initialized(tcf%kT_x) .or. &
          .not. Initialized(tcf%kT_y) .or. &
          .not. Initialized(tcf%kT_z)) then
        option%io_buffer = 'Diagonal components of thermal conductivity ' &
                        // 'must be specified if off-diagonal components are ' &
                        // 'provided in '// trim(error_string) // '.'
        call PrintErrMsg(option)
      endif
    else
      tcf%isotropic = PETSC_TRUE
      tcf%kT(:,:,:) = 0.0d0
      tcf%kT(1,1,1) = tcf%kT_dry
      tcf%kT(2,2,1) = tcf%kT_dry
      tcf%kT(3,3,1) = tcf%kT_dry
      tcf%kT(1,1,2) = tcf%kT_wet
      tcf%kT(2,2,2) = tcf%kT_wet
      tcf%kT(3,3,2) = tcf%kT_wet
      tcf%kT(1,1,3) = tcf%kT_frozen
      tcf%kT(2,2,3) = tcf%kT_frozen
      tcf%kT(3,3,3) = tcf%kT_frozen
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

subroutine TCondTensorToScalar(this,dist,option)
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

  PetscReal :: kTd(3,3) ! dry thermal conductivity tensor
  PetscReal :: kTw(3,3) ! wet thermal conductivity tensor
  PetscReal :: kTf(3,3) ! frozen thermal conductivity tensor
  PetscInt :: i

  select type(tcf => this)
  ! -------------------------
  class is(kT_default_type)
    if (tcf%isotropic) then
      return
    endif

    kTd = tcf%kT(:,:,1)
    kTw = tcf%kT(:,:,2)

    if (tcf%full_tensor) then
      tcf%kT_dry = FullTCondTensorToScalar(kTd,dist,option)
      tcf%kT_wet = FullTCondTensorToScalar(kTw,dist,option)
    elseif (.not. tcf%isotropic) then
      tcf%kT_dry = DiagTCondTensorToScalar(kTd,dist,option)
      tcf%kT_wet = DiagTCondTensorToScalar(kTw,dist,option)
    endif

    if (option%iflowmode == TH_MODE .or. option%iflowmode == TH_TS_MODE) then
      tcf%kT_dry = tcf%kT_dry * option%scale
      tcf%kT_wet = tcf%kT_wet * option%scale
    endif
  ! -------------------------
  class is(kT_frozen_type)
    if (tcf%isotropic) then
      return
    endif

    kTd = tcf%kT(:,:,1)
    kTw = tcf%kT(:,:,2)
    if (Initialized(tcf%kT_frozen)) then
      kTf = tcf%kT(:,:,3)
    endif

    if (tcf%full_tensor) then
      tcf%kT_dry = FullTCondTensorToScalar(kTd,dist,option)
      tcf%kT_wet = FullTCondTensorToScalar(kTw,dist,option)
      if (Initialized(tcf%kT_frozen)) then
        tcf%kT_frozen = FullTCondTensorToScalar(kTf,dist,option)
      endif
    elseif (.not. tcf%isotropic) then
      tcf%kT_dry = DiagTCondTensorToScalar(kTd,dist,option)
      tcf%kT_wet = DiagTCondTensorToScalar(kTw,dist,option)
      if (Initialized(tcf%kT_frozen)) then
        tcf%kT_frozen = DiagTCondTensorToScalar(kTf,dist,option)
      endif
    endif

    if (option%iflowmode == TH_MODE .or. option%iflowmode == TH_TS_MODE) then
      tcf%kT_dry = tcf%kT_dry * option%scale
      tcf%kT_wet = tcf%kT_wet * option%scale
        if (Initialized(tcf%kT_frozen)) then
          tcf%kT_frozen = tcf%kT_frozen * option%scale
        endif
    endif
  ! -------------------------
  type is(kT_composite_type)

    tcf%dist = dist

    do i = 1, size(tcf%ckT)
      select type(stcf => tcf%ckT(i)%ptr%thermal_conductivity_function)
      class is(kT_default_type)

        if (stcf%isotropic) then
          cycle
        endif

        ! Apply user-defined anisotropy to sub-function
        kTd = stcf%kT(:,:,1)
        kTw = stcf%kT(:,:,2)

        if (stcf%full_tensor) then
          stcf%kT_dry = FullTCondTensorToScalar(kTd,dist,option)
          stcf%kT_wet = FullTCondTensorToScalar(kTw,dist,option)
        elseif (.not. stcf%isotropic) then
          stcf%kT_dry = DiagTCondTensorToScalar(kTd,dist,option)
          stcf%kT_wet = DiagTCondTensorToScalar(kTw,dist,option)
        endif

      end select
    enddo

  end select

end subroutine TCondTensorToScalar

! ************************************************************************** !

function FullTCondTensorToScalar(kT,dist,option)
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
  PetscReal, intent(in)  :: dist(-1:3)
  PetscReal, intent(in)  :: kT(3,3)
  PetscReal :: FullTCondTensorToScalar

  PetscReal :: kx,ky,kz,kxy,kxz,kyz

  kx  = kT(1,1)
  ky  = kT(2,2)
  kz  = kT(3,3)
  kxy = kT(1,2)
  kxz = kT(1,3)
  kyz = kT(2,3)

  FullTCondTensorToScalar = kx*dabs(dist(1))**2 + &
                            ky*dabs(dist(2))**2 + &
                            kz*dabs(dist(3))**2 + &
                            2*kxy*dist(1)*dist(2) + &
                            2*kxz*dist(1)*dist(3) + &
                            2*kyz*dist(2)*dist(3)

end function FullTCondTensorToScalar

! ************************************************************************** !

function DiagTCondTensorToScalar(kT,dist,option)
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
  PetscReal, intent(in)  :: dist(-1:3)
  PetscReal, intent(in)  :: kT(3,3)
  PetscReal :: DiagTCondTensorToScalar

  PetscReal :: kx,ky,kz

  kx  = kT(1,1)
  ky  = kT(2,2)
  kz  = kT(3,3)

  DiagTCondTensorToScalar = kx*dabs(dist(1))**2 + &
                            ky*dabs(dist(2))**2 + &
                            kz*dabs(dist(3))**2

end function DiagTCondTensorToScalar

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

  do k = 0,2
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
      if (OptionIsIORank(option)) then
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

subroutine CompositeTCCList(list,tcf,option)

  use String_module
  use Option_module

  implicit none

  class(cc_thermal_type), pointer :: list
  class(kT_composite_type) :: tcf
  type(option_type) :: option

  class(cc_thermal_type), pointer :: cur_thermal_cc
  PetscInt :: i

  do i = 1, 3 ! iterate over X, Y, and Z
    cur_thermal_cc => list
    do
      if (.not.associated(cur_thermal_cc)) exit
      if (trim(cur_thermal_cc%name) == trim(tcf%lkT(i))) then
        tcf%ckT(i)%ptr => cur_thermal_cc
        exit
      else
        cur_thermal_cc => cur_thermal_cc%next
      endif
    enddo
    if (.not. associated(tcf%ckT(i)%ptr)) then
      option%io_buffer = 'Thermal conductivity sub-function "' &
                         // trim(tcf%lkT(i)) // &
                         '" not found in list of thermal characteristic curves.'
      call PrintErrMsg(option)
    endif
  enddo

end subroutine CompositeTCCList

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
      class is (kT_linear_type)
        write(id,'(a)') 'only saturation-dependent (linear)'
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
      class is (kT_ASM_dry_type)
        write(id,'(a)') 'temp.-dependent/optional sat. (asm. dry)'
        write(id,'(a29)',advance='no') 'kT_wet: '
        write(word1,*) tcf%kT_wet
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'kT_dry: '
        write(word1,*) tcf%kT_dry
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'temperature coefficient: '
        write(word1,*) tcf%dry_T_coeff
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'temperature exponent: '
        write(word1,*) tcf%dry_T_power
        write(id,'(a)') adjustl(trim(word1))
        !---------------------------------
      class is (kT_ASM_water_filled_type)
        write(id,'(a)') 'not sat.- or temp.-dependent (asm. water-filled)'
        write(id,'(a29)',advance='no') 'kT_dry: '
        write(word1,*) tcf%kT_dry
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'kT_water: '
        write(word1,*) tcf%kT_water
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'kT_solid: '
        write(word1,*) tcf%kT_solid
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'assembly porosity: '
        write(word1,*) tcf%porosity_asm
        write(id,'(a)') adjustl(trim(word1))
        !---------------------------------
      class is (kT_ASM_radial_type)
        write(id,'(a)') 'temp.- and sat.-dependent. (asm. radial)'
        write(id,'(a29)',advance='no') 'kT_dry: '
        write(word1,*) tcf%kT_dry
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'kT_water: '
        write(word1,*) tcf%kT_water
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'kT_solid: '
        write(word1,*) tcf%kT_solid
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'temperature coefficient: '
        write(word1,*) tcf%dry_T_coeff
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'temperature exponent: '
        write(word1,*) tcf%dry_T_power
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'assembly porosity: '
        write(word1,*) tcf%porosity_asm
        write(id,'(a)') adjustl(trim(word1))
        !---------------------------------
      class is (kT_ASM_axial_type)
        write(id,'(a)') 'sat.-dependent. (asm. axial)'
        write(id,'(a29)',advance='no') 'kT_water: '
        write(word1,*) tcf%kT_water
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'kT_solid: '
        write(word1,*) tcf%kT_solid
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'assembly porosity: '
        write(word1,*) tcf%porosity_asm
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
        write(id,'(a29)',advance='no') 'porosity effect: '
        write(word1,*) tcf%porosity_effect
        write(id,'(a)') adjustl(trim(word1))
        if (tcf%porosity_effect) then
           write(id,'(a29)',advance='no') 'reference por.: '
           write(word1,*) tcf%ref_por
           write(id,'(a)') adjustl(trim(word1))
           write(id,'(a29)',advance='no') 'por. exponent: '
           write(word1,*) tcf%por_exp
           write(id,'(a)') adjustl(trim(word1))
           write(id,'(a29)',advance='no') 'initial TC, constant term: '
           write(word1,*) tcf%b(1)
           write(id,'(a)') adjustl(trim(word1))
           write(id,'(a29)',advance='no') 'initial TC, linear term: '
           write(word1,*) tcf%b(2)
           write(id,'(a)') adjustl(trim(word1))
        end if
        !---------------------------------
      class is (kT_frozen_type)
        write(id,'(a)') 'liquid and ice sat.-dependent (frozen)'
        write(id,'(a29)',advance='no') 'kT_wet: '
        write(word1,*) tcf%kT_wet
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'kT_dry: '
        write(word1,*) tcf%kT_dry
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'kersten exponent: '
        write(word1,*) tcf%alpha
        write(id,'(a)') adjustl(trim(word1))
        if (Initialized(tcf%kT_frozen)) then
          write(id,'(a29)',advance='no') 'kT_frozen: '
          write(word1,*) tcf%kT_frozen
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'kersten exponent (frozen): '
          write(word1,*) tcf%alpha_fr
          write(id,'(a)') adjustl(trim(word1))
          write(id,'(a29)',advance='no') 'ice model index: '
          write(word1,*) tcf%ice_model
          write(id,'(a)') adjustl(trim(word1))
        endif
        !---------------------------------
      class is (kT_composite_type)
        write(id,'(a)') 'composite thermal characteristic curve'
        write(id,'(a29)',advance='no') 'sub-function X: '
        write(word1,*) trim(tcf%lkT(1))
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'sub-function Y: '
        write(word1,*) trim(tcf%lkT(2))
        write(id,'(a)') adjustl(trim(word1))
        write(id,'(a29)',advance='no') 'sub-function Z: '
        write(word1,*) trim(tcf%lkT(3))
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
