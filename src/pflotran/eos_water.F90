module EOS_Water_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Geometry_module
  use EOS_Database_module

  implicit none

  private

  ! module variables
  PetscReal :: constant_density
  PetscReal :: constant_enthalpy
  PetscReal :: constant_viscosity
  PetscReal :: constant_steam_density
  PetscReal :: constant_steam_enthalpy
  PetscReal :: constant_salinity
  PetscReal :: surface_density_kg

  ! exponential pressure
  PetscReal :: exp_p_reference_density
  PetscReal :: exp_p_reference_pressure
  PetscReal :: exp_p_water_compressibility

  ! exponential pressure temperature
  PetscReal :: exp_pt_reference_density
  PetscReal :: exp_pt_reference_pressure
  PetscReal :: exp_pt_water_compressibility
  PetscReal :: exp_pt_thermal_expansion
  PetscReal :: exp_pt_reference_temperature

  ! planes for planar eos
  type(plane_type) :: water_density_tp_plane
  type(plane_type) :: water_enthalpy_tp_plane
  type(plane_type) :: steam_density_tp_plane
  type(plane_type) :: steam_enthalpy_tp_plane


  ! quadratic
  PetscReal :: quadratic_reference_density
  PetscReal :: quadratic_reference_pressure
  PetscReal :: quadratic_wat_compressibility

  ! linear
  PetscReal :: linear_reference_density
  PetscReal :: linear_reference_pressure
  PetscReal :: linear_water_compressibility

  ! halite saturated brine
  PetscBool :: halite_saturated_brine

  ! PVT tables - eos_tables
  class(eos_table_type), pointer :: pvt_table => null()

  ! In order to support generic EOS subroutines, we need the following:
  ! 1. An interface declaration that defines the argument list (best to have
  !    "Dummy" appended.
  ! 2. A procedure pointer that is initially set to null.  This pointer is
  !    pointed to the appropriate subroutine later on (e.g. EOSWaterInit())
  ! 3. An interface for derivative/non-derivative versions

  ! procedure pointer declarations
  ! standard versions
  procedure(EOSWaterViscosityDummy), pointer :: EOSWaterViscosityPtr => null()
  procedure(EOSWaterSatPressDummy), pointer :: &
    EOSWaterSaturationPressurePtr => null()
  procedure(EOSWaterDensityDummy), pointer :: EOSWaterDensityPtr => null()
  procedure(EOSWaterEnthalpyDummy), pointer :: EOSWaterEnthalpyPtr => null()
  procedure(EOSWaterSteamDenEnthDummy), pointer :: &
    EOSWaterSteamDensityEnthalpyPtr => null()
  procedure(EOSWaterDensityIceDummy), pointer :: &
    EOSWaterDensityIcePtr => null()
  ! extended versions
  procedure(EOSWaterSatPressExtDummy), pointer :: &
    EOSWaterSaturationPressureExtPtr => null()
  procedure(EOSWaterViscosityExtDummy), pointer :: &
    EOSWaterViscosityExtPtr => null()
  procedure(EOSWaterDensityExtDummy), pointer :: &
    EOSWaterDensityExtPtr => null()
  procedure(EOSWaterEnthalpyExtDummy), pointer :: &
    EOSWaterEnthalpyExtPtr => null()
  procedure(EOSWaterInternalEnergyIceDummy), pointer :: &
    EOSWaterInternalEnergyIcePtr => null()

  ! interface blocks
  interface
  ! standard versions
    subroutine EOSWaterViscosityDummy(T, P, PS, dPS_dT, &
                                      calculate_derivatives, VW, &
                                      dVW_dT, dVW_dP, ierr,table_idxs)
      implicit none
      PetscReal, intent(in) :: T, P, PS, dPS_dT
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: VW
      PetscReal, intent(out) :: dVW_dT, dVW_dP
      PetscErrorCode, intent(inout) :: ierr
      PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
    end subroutine EOSWaterViscosityDummy
    subroutine EOSWaterSatPressDummy(T, calculate_derivatives, &
                                     PS, dPS_dT, ierr)
      implicit none
      PetscReal, intent(in) :: T
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: PS, dPS_dT
      PetscErrorCode, intent(inout) :: ierr
    end subroutine EOSWaterSatPressDummy
    subroutine EOSWaterDensityDummy(t,p,calculate_derivatives, &
                                    dw,dwmol,dwp,dwt,ierr,table_idxs)
      implicit none
      PetscReal, intent(in) :: t
      PetscReal, intent(in) :: p
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: dw,dwmol,dwp,dwt
      PetscErrorCode, intent(inout) :: ierr
      PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
    end subroutine EOSWaterDensityDummy
    subroutine EOSWaterEnthalpyDummy(t,p,calculate_derivatives, &
                                     hw,hwp,hwt,ierr)
      implicit none
      PetscReal, intent(in) :: t
      PetscReal, intent(in) :: p
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: hw,hwp,hwt
      PetscErrorCode, intent(inout) :: ierr
    end subroutine EOSWaterEnthalpyDummy
    subroutine EOSWaterSteamDenEnthDummy(t,p,calculate_derivatives, &
                                         dg,dgmol,hg,dgp,dgt,hgp,hgt,ierr)
      implicit none
      PetscReal, intent(in) :: t
      PetscReal, intent(in) :: p
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: dg,dgmol,dgp,dgt
      PetscReal, intent(out) :: hg,hgp,hgt
      PetscErrorCode, intent(inout) :: ierr
    end subroutine EOSWaterSteamDenEnthDummy
    subroutine EOSWaterDensityIceDummy(t,p,calculate_derivatives, &
                                       dw,dwp,dwt,ierr)
      implicit none
      PetscReal, intent(in) :: t
      PetscReal, intent(in) :: p
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: dw,dwp,dwt
      PetscErrorCode, intent(inout) :: ierr
    end subroutine EOSWaterDensityIceDummy
    ! Extended versions
    subroutine EOSWaterSatPressExtDummy(T, aux, calculate_derivatives,&
                                        PS, dPS_dT, ierr)
      implicit none
      PetscReal, intent(in) :: T
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(in) :: aux(*)
      PetscReal, intent(out) :: PS, dPS_dT
      PetscErrorCode, intent(inout) :: ierr
    end subroutine EOSWaterSatPressExtDummy
    subroutine EOSWaterViscosityExtDummy(T, P, PS, dPS_dT, aux, &
                                         calculate_derivatives, VW, &
                                         dVW_dT, dVW_dP, ierr)
      implicit none
      PetscReal, intent(in) :: T, P, PS, dPS_dT
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(in) :: aux(*)
      PetscReal, intent(out) :: VW
      PetscReal, intent(out) :: dVW_dT, dVW_dP
      PetscErrorCode, intent(inout) :: ierr
    end subroutine EOSWaterViscosityExtDummy
    subroutine EOSWaterDensityExtDummy(t,p,aux,calculate_derivatives, &
                                       dw,dwmol,dwp,dwt,ierr)
      implicit none
      PetscReal, intent(in) :: t
      PetscReal, intent(in) :: p
      PetscReal, intent(in) :: aux(*)
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: dw,dwmol,dwp,dwt
      PetscErrorCode, intent(inout) :: ierr
    end subroutine EOSWaterDensityExtDummy
    subroutine EOSWaterEnthalpyExtDummy(t,p,aux,calculate_derivatives, &
                                        hw,hwp,hwt,ierr)
      implicit none
      PetscReal, intent(in) :: t
      PetscReal, intent(in) :: p
      PetscReal, intent(in) :: aux(*)
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: hw,hwp,hwt
      PetscErrorCode, intent(inout) :: ierr
    end subroutine EOSWaterEnthalpyExtDummy
    subroutine EOSWaterInternalEnergyIceDummy(t,u,calculate_derivatives,&
                                              du_ice_dT,du_ice_dP,ierr)
      implicit none
      PetscReal, intent(in) :: t
      PetscReal, intent(out) :: u
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: du_ice_dT
      PetscReal, intent(out) :: du_ice_dP
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSWaterInternalEnergyIceDummy
  end interface

  ! interfaces for derivative/non-derivative versions that are visible outside
  ! the module.
  ! standard versions
  interface EOSWaterViscosity
    procedure EOSWaterViscosityNoDerive
    procedure EOSWaterViscosityDerive
  end interface
  interface EOSWaterSaturationPressure
    procedure EOSWaterSatPresNoDerive
    procedure EOSWaterSatPresDerive
  end interface
  interface EOSWaterDensity
    procedure EOSWaterDensityNoDerive
    procedure EOSWaterDensityDerive
  end interface
  interface EOSWaterEnthalpy
    procedure EOSWaterEnthalpyNoDerive
    procedure EOSWaterEnthalpyDerive
  end interface
  interface EOSWaterSteamDensityEnthalpy
    procedure EOSWaterSteamDenEnthNoDerive
    procedure EOSWaterSteamDenEnthDerive
  end interface
  interface EOSWaterDensityIce
    procedure EOSWaterDensityIceNoDerive
    procedure EOSWaterDensityIceDerive
  end interface
  ! Extended versions
  interface EOSWaterSaturationPressureExt
    procedure EOSWaterSatPresExtNoDerive
    procedure EOSWaterSatPresExtDerive
  end interface
  interface EOSWaterViscosityExt
    procedure EOSWaterViscosityExtNoDerive
    procedure EOSWaterViscosityExtDerive
  end interface
  interface EOSWaterDensityExt
    procedure EOSWaterDensityExtNoDerive
    procedure EOSWaterDensityExtDerive
!geh: very useful for debuggin
!    procedure EOSWaterDensityExtNumericalDerive
  end interface
  interface EOSWaterEnthalpyExt
    procedure EOSWaterEnthalpyExtNoDerive
    procedure EOSWaterEnthalpyExtDerive
  end interface
  interface EOSWaterInternalEnergyIce
    procedure EOSWaterInternalEnergyIceNoDerive
    procedure EOSWaterInternalEnergyIceDerive
  end interface

  ! the "public" definition that makes subroutines visible outside.
  public :: EOSWaterInit, &
            EOSWaterVerify, &
            EOSWaterViscosity, &
            EOSWaterSaturationPressure, &
            EOSWaterDensity, &
            EOSWaterEnthalpy, &
            EOSWaterSteamDensityEnthalpy, &
            EOSWaterDuanMixture, &
            EOSWaterViscosityNaCl, &
            EOSWaterInternalEnergyIce, &
            EOSWaterDensityIcePainter, &
            EOSWaterSaturationTemperature, &
            EOSWaterDensityIce, &
            EOSWaterDensityTGDPB01, &
            EOSWaterDensityBRAGFLO, &
            EOSWaterSaturationPressureExt, &
            EOSWaterViscosityExt, &
            EOSWaterDensityExt, &
            EOSWaterEnthalpyExt, &
            EOSWaterInputRecord, &
            EOSWaterSetSalinity, &
            EOSWaterSolubility, &
            EOSWaterComputeSalinity

  public :: EOSWaterSetDensity, &
            EOSWaterSetEnthalpy, &
            EOSWaterSetViscosity, &
            EOSWaterSetSaturationPressure, &
            EOSWaterSetSteamDensity, &
            EOSWaterSetSteamEnthalpy, &
            EOSWaterSetIceInternalEnergy, &
            EOSWaterSetWaterTab, &
            EOSWaterSetSurfaceDensity, &
            EOSWaterGetSurfaceDensity, &
            EOSWaterTableProcess, &
            EOSWaterSurfaceTension, &
            EOSWaterKelvin

  public :: TestEOSWaterBatzleAndWang, &
            EOSWaterTest, &
            EOSWaterSteamTest

  contains

! ************************************************************************** !

subroutine EOSWaterInit()

  implicit none

  constant_density = UNINITIALIZED_DOUBLE
  constant_viscosity = UNINITIALIZED_DOUBLE
  constant_enthalpy = UNINITIALIZED_DOUBLE
  constant_steam_density = UNINITIALIZED_DOUBLE
  constant_steam_enthalpy = UNINITIALIZED_DOUBLE
  constant_salinity = UNINITIALIZED_DOUBLE
  surface_density_kg = UNINITIALIZED_DOUBLE
  exp_p_reference_density = UNINITIALIZED_DOUBLE
  exp_p_reference_pressure = UNINITIALIZED_DOUBLE
  exp_p_water_compressibility = UNINITIALIZED_DOUBLE
  exp_pt_reference_density = UNINITIALIZED_DOUBLE
  exp_pt_reference_pressure = UNINITIALIZED_DOUBLE
  exp_pt_reference_temperature = UNINITIALIZED_DOUBLE
  exp_pt_water_compressibility = UNINITIALIZED_DOUBLE
  exp_pt_thermal_expansion = UNINITIALIZED_DOUBLE
  quadratic_reference_density = UNINITIALIZED_DOUBLE
  quadratic_reference_pressure = UNINITIALIZED_DOUBLE
  quadratic_wat_compressibility = UNINITIALIZED_DOUBLE
  halite_saturated_brine = PETSC_FALSE

  ! standard versions
  EOSWaterDensityPtr => EOSWaterDensityIFC67
  EOSWaterEnthalpyPtr => EOSWaterEnthalpyIFC67
  EOSWaterViscosityPtr => EOSWaterViscosity1
  EOSWaterSaturationPressurePtr => EOSWaterSaturationPressureIFC67
  EOSWaterSteamDensityEnthalpyPtr => EOSWaterSteamDensityEnthalpyIFC67
  EOSWaterDensityIcePtr => EOSWaterDensityIcePainter
  EOSWaterInternalEnergyIcePtr => EOSWaterInternalEnergyIceDefault

  ! extended versions
  EOSWaterSaturationPressureExtPtr => EOSWaterSaturationPressureHaasExt
  EOSWaterViscosityExtPtr => EOSWaterViscosityKestinExt
  EOSWaterDensityExtPtr => EOSWaterDensityBatzleAndWangExt
  EOSWaterEnthalpyExtPtr => EOSWaterEnthalpyDriesnerExt

end subroutine EOSWaterInit

! ************************************************************************** !

subroutine EOSWaterVerify(ierr,error_string)

  implicit none

  PetscErrorCode, intent(inout) :: ierr
  character(len=MAXSTRINGLENGTH), intent(out) :: error_string

  ierr = 0
  error_string = ''
  if ((associated(EOSWaterDensityPtr,EOSWaterDensityIFC67) .and. &
        Initialized(constant_density)) .or. &
      (associated(EOSWaterEnthalpyPtr,EOSWaterEnthalpyIFC67) .and. &
        Initialized(constant_enthalpy)) &
     ) then
    ierr = 1
  endif

  if (associated(EOSWaterDensityPtr,EOSWaterDensityConstant) .and. &
      Uninitialized(constant_density)) then
    error_string = trim(error_string) // &
      ' CONSTANT density not set.'
    ierr = 1
  endif

  if (associated(EOSWaterEnthalpyPtr,EOSWaterEnthalpyConstant) .and. &
      Uninitialized(constant_enthalpy)) then
    error_string = trim(error_string) // &
      ' CONSTANT enthalpy not set.'
    ierr = 1
  endif

  if (associated(EOSWaterDensityPtr,EOSWaterDensityExpPressure) .and. &
      (Uninitialized(exp_p_reference_density) .or. &
       Uninitialized(exp_p_reference_pressure) .or. &
       Uninitialized(exp_p_water_compressibility))) then
    error_string = trim(error_string) // &
      ' Exponential parameters incorrect.'
    ierr = 1
  endif

  if (associated(EOSWaterDensityPtr,EOSWaterDensityExpPressureTemp) .and. &
      (Uninitialized(exp_pt_reference_density) .or. &
       Uninitialized(exp_pt_reference_pressure) .or. &
       Uninitialized(exp_pt_reference_temperature) .or. &
       Uninitialized(exp_pt_water_compressibility) .or. &
       Uninitialized(exp_pt_thermal_expansion))) then
    error_string = trim(error_string) // &
      ' Exponential Temperature parameters incorrect.'
    ierr = 1
  endif

  if (associated(EOSWaterDensityPtr,EOSWaterDensityLinear) .and. &
      (Uninitialized(linear_reference_density) .or. &
       Uninitialized(linear_reference_pressure) .or. &
       Uninitialized(linear_water_compressibility))) then
    error_string = trim(error_string) // &
      ' Linear parameters incorrect.'
     ierr = 1
  endif

  if (associated(EOSWaterDensityPtr,EOSWaterDensityBRAGFLO) .and. &
      (Uninitialized(exp_p_reference_density) .or. &
       Uninitialized(exp_p_reference_pressure) .or. &
       Uninitialized(exp_p_water_compressibility))) then
    error_string = trim(error_string) // &
      ' BRAGFLO parameters incorrect.'
    ierr = 1
  endif

  if ((associated(EOSWaterViscosityPtr, &
                  EOSWaterViscosityConstant) .and. &
       Uninitialized(constant_viscosity)) .or. &
      (associated(EOSWaterViscosityPtr, &
                  EOSWaterViscosity1) .and. &
       Initialized(constant_viscosity))) then
    ierr = 1
  endif

  if (associated(EOSWaterSteamDensityEnthalpyPtr, &
                 EOSWaterSteamDenEnthConstant) .and. &
      (Uninitialized(constant_steam_density) .or. &
       Uninitialized(constant_steam_enthalpy))) then
    if (Uninitialized(constant_steam_density)) then
      error_string = trim(error_string) // &
        ' CONSTANT steam density not set.'
    endif
    if (Uninitialized(constant_steam_enthalpy)) then
      error_string = trim(error_string) // &
        ' CONSTANT steam enthalpy not set.'
    endif
    ierr = 1
  endif

end subroutine EOSWaterVerify

! ************************************************************************** !

subroutine EOSWaterSetDensity(keyword,aux)

  implicit none

  character(len=*) :: keyword
  PetscReal, optional :: aux(*)

  select case(keyword)
    case('CONSTANT')
      constant_density = aux(1)
      EOSWaterDensityPtr => EOSWaterDensityConstant
    case('DEFAULT','IFC67')
      EOSWaterDensityPtr => EOSWaterDensityIFC67
      EOSWaterDensityExtPtr => EOSWaterDensityBatzleAndWangExt
    case('IF97')
      EOSWaterDensityPtr => EOSWaterDensityIF97
    case('EXPONENTIAL','EXPONENTIAL_PRESSURE')
      exp_p_reference_density = aux(1)
      exp_p_reference_pressure = aux(2)
      exp_p_water_compressibility = aux(3)
      EOSWaterDensityPtr => EOSWaterDensityExpPressure
    case('EXPONENTIAL_PRESSURE_TEMPERATURE')
      exp_pt_reference_density = aux(1)
      exp_pt_reference_pressure = aux(2)
      exp_pt_reference_temperature = aux(3)
      exp_pt_water_compressibility = aux(4)
      exp_pt_thermal_expansion = aux(5)
      EOSWaterDensityPtr => EOSWaterDensityExpPressureTemp
    case('LINEAR')
      linear_reference_density = aux(1)
      linear_reference_pressure = aux(2)
      linear_water_compressibility = aux(3)
      EOSWaterDensityPtr => EOSWaterDensityLinear
    case('BRAGFLO')
      exp_p_reference_density = aux(1)
      exp_p_reference_pressure = aux(2)
      exp_p_water_compressibility = aux(3)
      EOSWaterDensityPtr => EOSWaterDensityBRAGFLO
    case('QUADRATIC')
      if (Initialized(aux(1))) then
        quadratic_reference_density = 999.014d0 !kg/m3
      else
        quadratic_reference_density = aux(1)
      end if
      if (Initialized(aux(2))) then
        quadratic_reference_pressure = 1.0d5 !Pa
      else
        quadratic_reference_pressure = aux(2)
      end if
      if (Initialized(aux(3))) then
        quadratic_wat_compressibility = 3.94769306686405d-10 !1/Pa
      else
        quadratic_wat_compressibility = aux(3)
      end if
      EOSWaterDensityPtr => EOSWaterDensityQuadratic
    case('TRANGENSTEIN')
      EOSWaterDensityPtr => EOSWaterDensityTrangenstein
    case('PLANAR')
      EOSWaterDensityPtr => EOSWaterDensityTPPlanar
      call EOSWaterDensityTPPlanarSetup()
    case('TGDPB01')
      EOSWaterDensityPtr => EOSWaterDensityTGDPB01
    case('PAINTER')
      EOSWaterDensityPtr => EOSWaterDensityPainter
    case('BATZLE_AND_WANG')
      EOSWaterDensityPtr => EOSWaterDensityBatzleAndWang
      EOSWaterDensityExtPtr => EOSWaterDensityBatzleAndWangExt
    case('SPARROW')
      EOSWaterDensityExtPtr => EOSWaterDensitySparrowExt
    case('DRIESNER')
      EOSWaterDensityExtPtr => EOSWaterDensityDriesnerExt
    case default
      print *, 'Unknown pointer type "' // trim(keyword) // &
        '" in EOSWaterSetDensity().'
      stop
  end select

end subroutine EOSWaterSetDensity

! ************************************************************************** !

subroutine EOSWaterSetSalinity(input,word,option)

  use Input_Aux_module
  use Option_module

  implicit none

  type(input_type), pointer :: input
  character(len=*) :: word
  type(option_type) :: option

  select case (word)
    case('CONSTANT')
      call InputReadDouble(input,option,constant_salinity)
      call InputErrorMsg(input,option,'CONSTANT','EOS WATER, SALINITY')
      call InputReadAndConvertUnits(input,constant_salinity,'g/g',&
                         'EOS,WATER,HALITE_SATURATED_BRINE,SALINITY,CONSTANT',&
                          option)
    case('TEMPERATURE_CONTROLLED')
      halite_saturated_brine = PETSC_TRUE
      constant_salinity = 0.d0 ! Initialize it
    case default
      call InputKeywordUnrecognized(input,word,'EOS,WATER,SALINITY',option)
  end select
  option%flow%density_depends_on_salinity = PETSC_TRUE

end subroutine EOSWaterSetSalinity

! ************************************************************************** !

subroutine EOSWaterSetEnthalpy(keyword,aux)

  implicit none

  character(len=*) :: keyword
  PetscReal, optional :: aux(*)

  select case(keyword)
    case('CONSTANT')
      constant_enthalpy = aux(1)
      EOSWaterEnthalpyPtr => EOSWaterEnthalpyConstant
    case('DEFAULT','IFC67')
      EOSWaterEnthalpyPtr => EOSWaterEnthalpyIFC67
    case('IF97')
      EOSWaterEnthalpyPtr => EOSWaterEnthalpyIF97
    case('PLANAR')
      EOSWaterEnthalpyPtr => EOSWaterEnthalpyTPPlanar
      call EOSWaterEnthalpyTPPlanarSetup()
    case('PAINTER')
      EOSWaterEnthalpyPtr => EOSWaterEnthalpyPainter
    case('SPARROW')
      EOSWaterEnthalpyExtPtr => EOSWaterEnthalpySparrowExt
    case('DRIESNER')
      EOSWaterEnthalpyExtPtr => EOSWaterEnthalpyDriesnerExt
    case default
      print *, 'Unknown pointer type "' // trim(keyword) // &
        '" in EOSWaterSetEnthalpy().'
      stop
  end select

end subroutine EOSWaterSetEnthalpy

! ************************************************************************** !

subroutine EOSWaterSetViscosity(keyword,aux)

  implicit none

  character(len=*), intent(in) :: keyword
  PetscReal, intent(in), optional :: aux(*)

  select case(keyword)
    case('CONSTANT')
      constant_viscosity = aux(1)
      EOSWaterViscosityPtr => EOSWaterViscosityConstant
      EOSWaterViscosityExtPtr => EOSWaterViscosityConstantExt
    case('DEFAULT')
      EOSWaterViscosityPtr => EOSWaterViscosity1
      EOSWaterViscosityExtPtr => EOSWaterViscosityKestinExt
    case('BATZLE_AND_WANG')
      EOSWaterViscosityPtr => EOSWaterViscosityBatzleAndWang
      EOSWaterViscosityExtPtr => EOSWaterViscosityBatzleAndWangExt
    case('KESTIN')
      EOSWaterViscosityExtPtr => EOSWaterViscosityKestinExt
    case('GRABOWSKI')
      EOSWaterViscosityPtr => EOSWaterViscosityGrabowski
    case default
      print *, 'Unknown pointer type "' // trim(keyword) // &
        '" in EOSWaterSetViscosity().'
      stop
  end select

end subroutine EOSWaterSetViscosity

! ************************************************************************** !

subroutine EOSWaterSetSaturationPressure(keyword,aux)

  implicit none

  character(len=*), intent(in) :: keyword
  PetscReal, intent(in), optional :: aux(*)

  select case(keyword)
    case('IFC67')
      EOSWaterSaturationPressurePtr => EOSWaterSaturationPressureIFC67
    case('IF97')
      EOSWaterSaturationPressurePtr => EOSWaterSaturationPressureIF97
    case('WAGNER_AND_PRUSS')
      EOSWaterSaturationPressurePtr => EOSWaterSatPresWagnerPruss
    case('HAAS')
      EOSWaterSaturationPressureExtPtr => EOSWaterSaturationPressureHaasExt
    case('SPARROW')
      EOSWaterSaturationPressureExtPtr => EOSWaterSatPressSparrowExt
    case('HUANG-ICE','ICE')
      EOSWaterSaturationPressurePtr => EOSWaterSaturationPressureIce
    case default
      print *, 'Unknown pointer type "' // trim(keyword) // &
        '" in EOSWaterSetSaturationPressure().'
      stop
  end select

end subroutine EOSWaterSetSaturationPressure

! ************************************************************************** !

subroutine EOSWaterSetSteamDensity(keyword,aux)

  implicit none

  character(len=*), intent(in) :: keyword
  PetscReal, intent(in), optional :: aux(*)

  select case(keyword)
    case('CONSTANT')
      constant_steam_density = aux(1)
      EOSWaterSteamDensityEnthalpyPtr => EOSWaterSteamDenEnthConstant
    case('PLANAR')
      EOSWaterSteamDensityEnthalpyPtr => EOSWaterSteamDenEnthTPPlanar
      call EOSWaterSteamDenEnthTPPlanarSetup()
    case('IFC67')
      EOSWaterSteamDensityEnthalpyPtr => EOSWaterSteamDensityEnthalpyIFC67
    case('IF97')
      EOSWaterSteamDensityEnthalpyPtr => EOSWaterSteamDensityEnthalpyIF97
    case default
      print *, 'Unknown pointer type "' // trim(keyword) // &
        '" in EOSWaterSetSteamDensity().'
      stop
  end select

end subroutine EOSWaterSetSteamDensity

! ************************************************************************** !

subroutine EOSWaterSetSteamEnthalpy(keyword,aux)

  implicit none

  character(len=*), intent(in) :: keyword
  PetscReal, intent(in), optional :: aux(*)

  select case(keyword)
    case('CONSTANT')
      constant_steam_enthalpy = aux(1)
      EOSWaterSteamDensityEnthalpyPtr => EOSWaterSteamDenEnthConstant
    case('PLANAR')
      EOSWaterSteamDensityEnthalpyPtr => EOSWaterSteamDenEnthTPPlanar
      call EOSWaterSteamDenEnthTPPlanarSetup()
    case('DEFAULT','IFC67')
      EOSWaterSteamDensityEnthalpyPtr => EOSWaterSteamDensityEnthalpyIFC67
    case('IF97')
      EOSWaterSteamDensityEnthalpyPtr => EOSWaterSteamDensityEnthalpyIF97
    case default
      print *, 'Unknown pointer type "' // trim(keyword) // &
        '" in EOSWaterSetSteamEnthalpy().'
      stop
  end select

end subroutine EOSWaterSetSteamEnthalpy

! ************************************************************************** !

subroutine EOSWaterSetIceInternalEnergy(keyword,aux)

  implicit none
  
character(len=*), intent(in) :: keyword
  PetscReal, intent(in), optional :: aux(*)

  select case(keyword)
    case('DEFAULT')
      EOSWaterInternalEnergyIcePtr => EOSWaterInternalEnergyIceDefault
    case('FUKUSAKO')
      EOSWaterInternalEnergyIcePtr => EOSWaterInternalEnergyIceFukusako
    case default
      print *, 'Unknown pointer type "' // trim(keyword) // &
        '" in EOSWaterSetIceInternalEnergy().'
      stop
   end select

end subroutine EOSWaterSetIceInternalEnergy
   
! ************************************************************************** !

subroutine EOSWaterSetWaterTab(input,option)
  !
  ! Author: Paolo Orsini
  ! Date: 03/20/19
  !
  ! Set up a Water Table for density and viscosity

  use Option_module
  use Input_Aux_module
  use Lookup_Table_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option

  type(lookup_table_var_type), pointer :: db_var => null()
  character(len=MAXWORDLENGTH) :: internal_units, user_units
  PetscInt :: data_idx

  pvt_table => EOSTableCreate('WATERTAB',option)

  pvt_table%num_prop = 2

  ! units initially assing default values - overwritten by units specified
  ! in the table input
  internal_units = '' !assign default value by SetDefaultInternalUnits
  user_units = ''     !assign default value by SetMetricUnits

  !adding FVF
  data_idx = 1 !position of FVF in the table (after pressure)
  db_var => CreateLookupTableVar(EOS_FVF,internal_units,user_units,data_idx)
  call pvt_table%AddEOSProp(db_var,option)
  nullify(db_var)

  !adding VISCOSITY
  data_idx = 2 !position of viscosity in the table (after FVF)
  db_var => CreateLookupTableVar(EOS_VISCOSITY,internal_units,user_units, &
                                 data_idx)
  call pvt_table%AddEOSProp(db_var,option)
  nullify(db_var)

  !set Default internal must be called before Set Metric
  call pvt_table%SetDefaultInternalUnits(option)
  call pvt_table%SetMetricUnits(option)

  call pvt_table%Read(input,option)

  call EOSTableAddToList(pvt_table,eos_table_list)

  EOSWaterViscosityPtr => EOSWaterViscosityTable
  EOSWaterDensityPtr => EOSWaterDensityTable

end subroutine EOSWaterSetWaterTab

! ************************************************************************** !

subroutine EOSWaterSetSurfaceDensity(input_ref_density)

  implicit none

  PetscReal, intent(in) :: input_ref_density

  surface_density_kg = input_ref_density

end subroutine EOSWaterSetSurfaceDensity


! ************************************************************************** !
function EOSWaterGetSurfaceDensity()

  implicit none

  PetscReal :: EOSWaterGetSurfaceDensity

  EOSWaterGetSurfaceDensity= surface_density_kg

end function EOSWaterGetSurfaceDensity

! ************************************************************************** !

subroutine EOSWaterDensityTable(T,P,calculate_derivatives,dw,dwmol,dwp,dwt, &
                                ierr,table_idxs)
  !
  ! Author: Paolo Orsini
  ! Date: 03/20/19
  !
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscBool, intent(in) :: calculate_derivatives ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: dw     ! water density [kg/m^3]
  PetscReal, intent(out) :: dwmol     ! water density [kmol/m^3]
  PetscReal, intent(out) :: dwp ! derivative wrt table pressure [kmol/Pa]
  PetscReal, intent(out) :: dwt ! derivative wrt table temperature [kmol/C]
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  !ierr initialised in EOSPropGrad
  !Rho from pvt table is already in kmol/m3
  call pvt_table%EOSPropGrad(T,P,EOS_DENSITY,dwmol,dwt,dwp,ierr,table_idxs)
  !
  !dw = dwmol * fmw_wat !might consider to have a general formula weight
  !
  dw = dwmol * FMWH2O !kg/m^3

end subroutine EOSWaterDensityTable

! ************************************************************************** !

subroutine EOSWaterViscosityTable(T,P,PS,dPS_dT,&
                                  calculate_derivatives,VW, &
                                  dVW_dT, dVW_dP, ierr,table_idxs)
  !
  ! Author: Paolo Orsini
  ! Date: 03/20/19
  !
  implicit none

  PetscReal, intent(in) :: T       ! C
  PetscReal, intent(in) :: P       ! Pa
  PetscReal, intent(in) :: PS      ! Pa
  PetscReal, intent(in) :: dPS_dT  ! Pa/C
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: VW     ! Pa-s
  PetscReal, intent(out) :: dVW_dT !derivative viscosity wrt table temperature [Pa-s/C]
  PetscReal, intent(out) :: dVW_dP !derivative viscosity wrt table Pressure [Pa-s/Pa]
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  call pvt_table%EOSPropGrad(T,P,EOS_VISCOSITY,VW,dVW_dT,dVW_dP, &
                             ierr,table_idxs)

end subroutine EOSWaterViscosityTable

! ************************************************************************** !

subroutine EOSWaterTableProcess(option)
  !
  ! Author: Paolo Orsini
  ! Date: 03/21/19
  !
  ! Processes water pvt table - once the entire input deck has been read

  use Option_module

  implicit none

  type(option_type) :: option

  if (.not.associated(pvt_table)) return

  select case(pvt_table%name)
    case("WATERTAB")
      call pvt_table%ConvertFVFtoMolarDensity(FMWH2O,surface_density_kg)
  end select

end subroutine EOSWaterTableProcess

! ************************************************************************** !

subroutine EOSWaterViscosityNoDerive(T, P, PS, VW, ierr, table_idxs)

  implicit none

  PetscReal, intent(in) :: T, P, PS ! temperature, pressure, saturation_press
  PetscReal, intent(out) :: VW ! water viscosity
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: dPS_dT ! derivative of PS with respect to temp
  PetscReal :: dum1, dum2

  ierr = 0
  dPS_dT = 0.d0
  call EOSWaterViscosityPtr(T, P, PS, dPS_dT, PETSC_FALSE, VW, &
                            dum1, dum2, ierr, table_idxs)
end subroutine EOSWaterViscosityNoDerive

! ************************************************************************** !

subroutine EOSWaterViscosityDerive(T, P, PS, dPS_dT, VW, dVW_dT, &
                                   dVW_dP, ierr, table_idxs)

  implicit none

  PetscReal, intent(in) :: T, P, PS ! temperature, pressure, saturation_press
  PetscReal, intent(in) :: dPS_dT ! derivative of PS with respect to temp
  PetscReal, intent(out) :: VW ! water viscosity
  PetscReal, intent(out) :: dVW_dT, dVW_dP ! derivatives
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  ierr = 0
  call EOSWaterViscosityPtr(T, P, PS, dPS_dT, PETSC_TRUE, VW, &
                            dVW_dT, dVW_dP, ierr, table_idxs)

end subroutine EOSWaterViscosityDerive

! ************************************************************************** !

subroutine EOSWaterSatPresNoDerive(T, PS, ierr)

  implicit none

  PetscReal, intent(in) :: T ! temperature
  PetscReal, intent(out) :: PS ! Saturation pres. and derivative
  PetscErrorCode, intent(inout) :: ierr

  PetscReal :: dummy

  ierr = 0
  call EOSWaterSaturationPressurePtr(T, PETSC_FALSE, PS, dummy, ierr)

end subroutine EOSWaterSatPresNoDerive

! ************************************************************************** !

subroutine EOSWaterSatPresDerive(T, PS, dPS_dT, ierr)

  implicit none

  PetscReal, intent(in) :: T ! temperature
  PetscReal, intent(out) :: PS, dPS_dT ! Saturation pres. and derivative
  PetscErrorCode, intent(inout) :: ierr

  ierr = 0
  call EOSWaterSaturationPressurePtr(T, PETSC_TRUE, PS, dPS_dT, ierr)

end subroutine EOSWaterSatPresDerive

! ************************************************************************** !

subroutine EOSWaterSatPresExtNoDerive(T, aux, PS, ierr)

  implicit none

  PetscReal, intent(in) :: T ! temperature
  PetscReal, intent(in) :: aux(*)
  PetscReal, intent(out) :: PS ! Saturation pres. and derivative
  PetscErrorCode, intent(inout) :: ierr

  PetscReal :: dummy

  ierr = 0
  call EOSWaterSaturationPressureExtPtr(T, aux, PETSC_FALSE, PS, dummy, ierr)

end subroutine EOSWaterSatPresExtNoDerive

! ************************************************************************** !

subroutine EOSWaterSatPresExtDerive(T, aux, PS, dPS_dT, ierr)

  implicit none

  PetscReal, intent(in) :: T ! temperature
  PetscReal, intent(in) :: aux(*)
  PetscReal, intent(out) :: PS, dPS_dT ! Saturation pres. and derivative
  PetscErrorCode, intent(inout) :: ierr

  ierr = 0
  call EOSWaterSaturationPressureExtPtr(T, aux, PETSC_TRUE, PS, dPS_dT, ierr)

end subroutine EOSWaterSatPresExtDerive

! ************************************************************************** !

subroutine EOSWaterDensityNoDerive(t,p,dw,dwmol,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscReal, intent(out) :: dw,dwmol
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: dum1, dum2

  ierr = 0
  call EOSWaterDensityPtr(t,p,PETSC_FALSE,dw,dwmol,dum1,dum2,ierr,table_idxs)

end subroutine EOSWaterDensityNoDerive

! ************************************************************************** !

subroutine EOSWaterDensityDerive(t,p,dw,dwmol,dwp,dwt,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  ierr = 0
  call EOSWaterDensityPtr(t,p,PETSC_TRUE,dw,dwmol,dwp,dwt,ierr,table_idxs)

end subroutine EOSWaterDensityDerive

! ************************************************************************** !

subroutine EOSWaterEnthalpyNoDerive(t,p,hw,ierr)

  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscReal, intent(out) :: hw
  PetscErrorCode, intent(inout) :: ierr

  PetscReal :: dum1, dum2

  ierr = 0
  call EOSWaterEnthalpyPtr(t,p,PETSC_FALSE,hw,dum1,dum2,ierr)

end subroutine EOSWaterEnthalpyNoDerive

! ************************************************************************** !

subroutine EOSWaterEnthalpyDerive(t,p,hw,hwp,hwt,ierr)
  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscReal, intent(out) :: hw,hwp,hwt
  PetscErrorCode, intent(inout) :: ierr

  ierr = 0
  call EOSWaterEnthalpyPtr(t,p,PETSC_TRUE,hw,hwp,hwt,ierr)

end subroutine EOSWaterEnthalpyDerive

! ************************************************************************** !

subroutine EOSWaterInternalEnergyIceNoDerive(t,u,ierr)
  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(out) :: u
  PetscErrorCode, intent(out) :: ierr

  PetscReal :: dum1, dum2

  ierr = 0
  call EOSWaterInternalEnergyIcePtr(t,u,PETSC_FALSE,dum1,dum2,ierr)

end subroutine EOSWaterInternalEnergyIceNoDerive

! ************************************************************************** !

subroutine EOSWaterInternalEnergyIceDerive(t,u,du_ice_dT,du_ice_dP,ierr)
  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(out) :: u
  PetscReal, intent(out) :: du_ice_dT, du_ice_dP
  PetscErrorCode, intent(out) :: ierr

  ierr = 0
  call EOSWaterInternalEnergyIcePtr(t,u,PETSC_TRUE,du_ice_dT,du_ice_dP,ierr)

end subroutine EOSWaterInternalEnergyIceDerive

! ************************************************************************** !

subroutine EOSWaterViscosityExtNoDerive(T, P, PS, aux, VW, ierr)

  implicit none

  PetscReal, intent(in) :: T, P, PS ! temperature, pressure, saturation_press
  PetscReal, intent(in) :: aux(*)
  PetscReal, intent(out) :: VW ! water viscosity
  PetscErrorCode, intent(inout) :: ierr

  PetscReal :: dPS_dT ! derivative of PS with respect to temp
  PetscReal :: dum1, dum2

  ierr = 0
  dPS_dT = 0.d0
  call EOSWaterViscosityExtPtr(T, P, PS, dPS_dT, aux, PETSC_FALSE, VW, &
                               dum1, dum2, ierr)

end subroutine EOSWaterViscosityExtNoDerive

! ************************************************************************** !

subroutine EOSWaterViscosityExtDerive(T, P, PS, dPS_dT, aux, VW, dVW_dT, &
                                      dVW_dP, ierr)

  implicit none

  PetscReal, intent(in) :: T, P, PS ! temperature, pressure, saturation_press
  PetscReal, intent(in) :: dPS_dT ! derivative of PS with respect to temp
  PetscReal, intent(in) :: aux(*)
  PetscReal, intent(out) :: VW ! water viscosity
  PetscReal, intent(out) :: dVW_dT, dVW_dP ! derivatives
  PetscErrorCode, intent(inout) :: ierr

  ierr = 0
  call EOSWaterViscosityExtPtr(T, P, PS, dPS_dT, aux, PETSC_TRUE, VW, &
                               dVW_dT, dVW_dP, ierr)

end subroutine EOSWaterViscosityExtDerive

! ************************************************************************** !

subroutine EOSWaterDensityExtNoDerive(t,p,aux,dw,dwmol,ierr)

  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscReal, intent(in) :: aux(*)
  PetscReal, intent(out) :: dw,dwmol
  PetscErrorCode, intent(inout) :: ierr

  PetscReal :: dum1, dum2

  ierr = 0
  call EOSWaterDensityExtPtr(t,p,aux,PETSC_FALSE,dw,dwmol,dum1,dum2,ierr)

end subroutine EOSWaterDensityExtNoDerive

! ************************************************************************** !

subroutine EOSWaterDensityExtDerive(t,p,aux,dw,dwmol,dwp,dwt,ierr)
  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscReal, intent(in) :: aux(*)
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(inout) :: ierr

  ierr = 0
  call EOSWaterDensityExtPtr(t,p,aux,PETSC_TRUE,dw,dwmol,dwp,dwt,ierr)

end subroutine EOSWaterDensityExtDerive

! ************************************************************************** !

subroutine EOSWaterEnthalpyExtNoDerive(t,p,aux,hw,ierr)

  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscReal, intent(in) :: aux(*)
  PetscReal, intent(out) :: hw
  PetscErrorCode, intent(inout) :: ierr

  PetscReal :: dum1, dum2

  ierr = 0
  call EOSWaterEnthalpyExtPtr(t,p,aux,PETSC_FALSE,hw,dum1,dum2,ierr)

end subroutine EOSWaterEnthalpyExtNoDerive

! ************************************************************************** !

subroutine EOSWaterEnthalpyExtDerive(t,p,aux,hw,hwp,hwt,ierr)
  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscReal, intent(in) :: aux(*)
  PetscReal, intent(out) :: hw,hwp,hwt
  PetscErrorCode, intent(inout) :: ierr

  ierr = 0
  call EOSWaterEnthalpyExtPtr(t,p,aux,PETSC_TRUE,hw,hwp,hwt,ierr)

end subroutine EOSWaterEnthalpyExtDerive

! ************************************************************************** !

subroutine EOSWaterViscosity1(T, P, PS, dPS_dT, calculate_derivatives, &
                              VW, dVW_dT, dVW_dP, ierr,table_idxs)

! Calculates the viscosity of water and derivatives as a function of
! temperature, pressure, and saturation pressure.

  implicit none

  PetscReal, intent(in) :: T, P, PS ! temperature, pressure, saturation_press
  PetscReal, intent(in) :: dPS_dT ! derivative of PS with respect to temp
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: VW ! water viscosity
  PetscReal, intent(out) :: dVW_dT, dVW_dP ! derivatives
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: EX, PHI, AM, pwr, aln10

  EX  = 247.8d0/(T+133.15d0)
  PHI = 1.0467d0*(T-31.85d0)
  !geh: added max(P,PS)-PS in place of P-PS
  !geh: here P should be the maximum of Pl and Pg
  AM  = 1.d0+PHI*(max(P,PS)-PS)*1.d-11
  pwr = 10.d0**EX
  VW = 1.d-7*AM*241.4d0*pwr

  if (calculate_derivatives) then
    aln10 = log(10.d0)
    dVW_dT = VW/AM*1.d-11* &
            ! dAM_PHI_dT       dAM_PS_dT
            (1.0467d0*(P-PS) - PHI*dPS_dT) - &
            ! dpwr_EX_dT
            VW*aln10*247.8d0/(T+133.15d0)**2
    dVW_dP = VW/AM*PHI*1.d-11
  else
    dVW_dT = UNINITIALIZED_DOUBLE
    dVW_dP = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterViscosity1

! ************************************************************************** !

subroutine EOSWaterViscosityConstant(T, P, PS, dPS_dT, &
                                     calculate_derivatives, &
                                      VW, dVW_dT, dVW_dP, ierr,table_idxs)

! Calculates the viscosity of water and derivatives as a function of
! temperature, pressure, and saturation pressure.

  implicit none

  PetscReal, intent(in) :: T, P, PS ! temperature, pressure, saturation_press
  PetscReal, intent(in) :: dPS_dT ! derivative of PS with respect to temp
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: VW ! water viscosity
  PetscReal, intent(out) :: dVW_dT, dVW_dP ! derivatives
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  VW = constant_viscosity

  dVW_dT = 0.d0
  dVW_dP = 0.d0

end subroutine EOSWaterViscosityConstant

! ************************************************************************** !

subroutine EOSWaterViscosityConstantExt(T, P, PS, dPS_dT, aux,&
                                        calculate_derivatives, &
                                        VW, dVW_dT, dVW_dP, ierr)

! Calculates the viscosity of water and derivatives as a function of
! temperature, pressure, and saturation pressure.

  implicit none

  PetscReal, intent(in) :: T, P, PS ! temperature, pressure, saturation_press
  PetscReal, intent(in) :: dPS_dT ! derivative of PS with respect to temp
  PetscReal, intent(in) :: aux(*)
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: VW ! water viscosity
  PetscReal, intent(out) :: dVW_dT, dVW_dP ! derivatives
  PetscErrorCode, intent(inout) :: ierr

  VW = constant_viscosity

  dVW_dT = 0.d0
  dVW_dP = 0.d0

end subroutine EOSWaterViscosityConstantExt

! ************************************************************************** !

subroutine EOSWaterSaturationPressureIFC67(T, calculate_derivatives, &
                                           PS, dPS_dT, ierr)

  implicit none

  PetscReal, intent(in) :: T ! temperature
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: PS, dPS_dT ! Saturation pres. and derivative
  PetscErrorCode, intent(inout) :: ierr

  PetscReal, save, dimension(9) :: A(9)
  PetscReal :: TC, SC, PCAP, E1, E2
  PetscReal :: one_m_tc, one_m_tc_sq, E2_bottom
  PetscReal :: dTC_dT, dSC_dTC, dE1_dTC, dE2_dTC, dPC_dSC, dPC_dTC

  DATA A/ &
    -7.691234564d0,-2.608023696d1,-1.681706546d2,6.423285504d1, &
    -1.189646225d2,4.167117320d0,2.097506760E1,1.d9,6.d0/

  if (T .GT. 500.d0) then
    ierr = 1
    return
  end if
  TC = (T+273.15d0)/H2O_CRITICAL_TEMPERATURE
  one_m_tc = 1.d0-TC
  one_m_tc_sq = one_m_tc*one_m_tc
  SC = A(1)*one_m_tc+A(2)*one_m_tc_sq+A(3)*one_m_tc**3.d0+ &
        A(4)*one_m_tc**4.d0+A(5)*one_m_tc**5.d0
  E1 = TC*(1.d0+A(6)*one_m_tc+A(7)*one_m_tc_sq)
  E2_bottom = A(8)*one_m_tc_sq+A(9)
  E2 = one_m_tc/E2_bottom
  PCAP = EXP(SC/E1-E2)

  PS = PCAP*H2O_CRITICAL_PRESSURE

  if (calculate_derivatives) then
    dTC_dT = 1.d0/H2O_CRITICAL_TEMPERATURE
    dSC_dTC = -A(1)-2.d0*A(2)*one_m_tc-3.d0*A(3)*one_m_tc_sq- &
              4.d0*A(4)*one_m_tc**3.-5.d0*A(5)*one_m_tc**4.
    dE1_dTC = (1.d0+A(6)*one_m_tc+A(7)*one_m_tc_sq)+ &
              TC*(-A(6)-2.d0*A(7)*one_m_tc)
    dE2_dTC = -1.d0/E2_bottom+one_m_tc/(E2_bottom*E2_bottom)*2.d0*one_m_tc
    dPC_dTC = (-SC/(E1*E1)*dE1_dTC-dE2_dTC)*PCAP
    dPC_dSC = 1.d0/E1*PCAP
    dPS_dT = (dPC_dSC*dSC_dTC+dPC_dTC)*dTC_dT*H2O_CRITICAL_PRESSURE
  else
    dPS_dT = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterSaturationPressureIFC67

! ************************************************************************** !

subroutine EOSWaterSaturationPressureIF97(T, calculate_derivatives, &
                                           PS, dPS_dT, ierr)

  !Author: Michael Nole
  !Date: 01/20/19
  !Water saturation Pressure as f(T) from IF97 standard, valid between 273.15K and
  !623.15K (Region 4). At temperatures above 623.15K, a quadratic function covers
  !the region where superheated steam properties are accurate.


  implicit none

  PetscReal, intent(in) :: T ! temperature
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: PS, dPS_dT ! Saturation pres. and derivative
  PetscErrorCode, intent(inout) :: ierr

  PetscReal, parameter :: n(10) = [0.11670521452767d4,-0.72421316703206d6, &
               -0.17073846940092d2,0.12020824702470d5,-0.32325550322333d7, &
               0.14915108613530d2,-0.48232657361591d4,0.40511340542057d6, &
               -0.23855557567849d0,0.65017534844798d3]
  PetscReal :: theta, T_temp, A, B, C
  PetscReal :: B2_4AC, dtheta_dt, da_dt, db_dt, dc_dt

  T_temp = T + 273.15d0

  if (T_temp > 623.15d0) then
    PS = 0.34805185628969d3 - 0.11671859879975d1*T_temp + &
            0.10192970039326d-2 * T_temp*T_temp
    if (calculate_derivatives) then
      dPS_DT = -0.11671859879975d1 + 2.d0 * 0.10192970039326d-2 *T_temp
    else
      dPS_DT = UNINITIALIZED_DOUBLE
    endif
  else

    if (T_temp < 273.15d0) T_temp = 273.15d0

    theta = T_temp + n(9)/(T_temp-n(10))
    A = theta*theta + n(1)*theta +n(2)
    B = n(3)*theta*theta + n(4)*theta + n(5)
    C = n(6)*theta*theta + n(7)*theta + n(8)
    B2_4AC = (B*B-4.d0*A*C)**(0.5d0)
    PS = 1.d6*(2.d0*C/((B2_4AC)-B))**4
    if (calculate_derivatives) then
      dtheta_dt = 1.d0 - n(9)/((T_temp - n(10))*(T_temp - n(10)))
      da_dt = 2.d0*theta*dtheta_dt +n(1)*dtheta_dt
      db_dt = 2.d0*n(3)*theta*dtheta_dt + n(4)*dtheta_dt
      dc_dt = 2.d0*n(6)*theta*dtheta_dt + n(7)*dtheta_dt

      dPS_DT = 1.d6*4.d0*(2.d0*C/(B2_4AC-B))**3 * (2.d0*dc_dt/(B2_4AC-B) + &
               2.d0*C/(-(B2_4AC-B)*(B2_4AC-B))*(-db_dt + 0.5d0/B2_4AC * &
               (2.d0*B*db_dt - 4.d0*(A*dc_dt + C*da_dt))))
    else
      dPS_DT = UNINITIALIZED_DOUBLE
    endif
  endif

end subroutine EOSWaterSaturationPressureIF97


! ************************************************************************** !
subroutine EOSWaterSaturationPressureHaasExt(T, aux, calculate_derivatives, PS, &
                                          dPS_dT, ierr)

  ! Water saturation pressure with dissolved halite
  ! Haas, 1976
  ! Physical properties of the coexisting phases and thermochemical properties of
  ! the H2O component in boiling NaCl solutions
  !
  ! Author: David Fukuyama
  ! Date: 12/20/21

  implicit none

  PetscReal, intent(in) :: T
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(in) :: aux(*)
  PetscReal, intent(out) :: PS, dPS_dT
  PetscErrorCode, intent(inout) :: ierr

  PetscReal :: ln_p !bars
  PetscReal :: T0, Tx, x
  PetscReal :: dT0_dT,da_dT,db_dT,dm_dT,dln_T0_dT,dz_dT,dy_dT,dw_dT
  PetscReal :: dln_p_dT
  PetscReal :: ln_T0
  PetscReal :: a, b, m
  PetscReal :: w, y, z
  PetscReal, parameter :: a1 = 5.93582d-6
  PetscReal, parameter :: a2 = -5.19386d-5
  PetscReal, parameter :: a3 = 1.23156d-5
  PetscReal, parameter :: b1 = 1.15420d-6
  PetscReal, parameter :: b2 = 1.141254d-7
  PetscReal, parameter :: b3 = -1.92476d-8
  PetscReal, parameter :: b4 = -1.70717d-9
  PetscReal, parameter :: b5 = 1.05390d-10
  PetscReal, parameter :: e0 = 12.50849d0
  PetscReal, parameter :: e1 = -4.616913d3
  PetscReal, parameter :: e2 = 3.193455d-4
  PetscReal, parameter :: e3 = 1.1965d-11
  PetscReal, parameter :: e4 = -1.013137d-2
  PetscReal, parameter :: e5 = -5.7148d-3
  PetscReal, parameter :: e6 = 2.9370d5

  Tx = T+273.15d0

  x = aux(1) ! mass fraction
  x = 1.d3*(x*1.d2)/(58.442d0*(1.d2-(x*1.d2))) !mol/kg
  ! ln(T0) = m*ln(Tx)
  ! Tx: brine temperature
  ! T0: temperature of H2O at 0 molal, at the same pressure as Tx

  a = 1.d0 + x*(a1 + x*(a2 + a3*x))
  b = x*(b1 + x*(b2 + x*(b3 + x*(b4 + b5*x))))
  m = 1.d0/(a+b*Tx)
  ln_T0 = m*log(Tx)
  T0 = exp(ln_T0)

  z = T0 + 1.d-2
  y = 647.27d0 - T0
  w = z**2 - e6

  ! Note- there is a typo in Eq. 6 of Haas (1976) that adds the third and
  ! fourth term instead of multiplying.
  ln_p = e0+e1/z+e2*w/z*(1.d1**(e3*w**2)-1.d0)+e4*1.d1**(e5*y**1.25d0)
  PS = exp(ln_p)*1.d5

  if (calculate_derivatives) then
    da_dT = 0.d0
    db_dT = 0.d0
    dm_dT = -(b)/(a+b*Tx)**2
    dln_T0_dT = m/Tx+dm_dT*log(Tx)
    dT0_dT = exp(ln_T0)*dln_T0_dT
    dz_dT = dT0_dT
    dy_dT = -dT0_dT
    dw_dT = 2*z*dZ_dT
    dln_p_dT = -e1*z**(-2.d0)*dz_dT+ &
               (e2*1.d1**(e3*w**2-1.d0)* &
                 (z*dw_dT*e3*log(100.d0)*w**2+1.d0)-w*dz_dT)/(z**2)+&
               e2*w/z*(e3*5.d0**(3.d0*w**2-1.d0)*8.d0**(w**2.d0)* &
                         dw_dT*log(10.d0))+ &
               2.87823d0*1.d1**(e5*y**1.25d0)*e4*e5*y**0.25d0*dy_dT
    dPS_dT = 1.d5*exp(ln_p)*dln_p_dT
  endif

end subroutine EOSWaterSaturationPressureHaasExt
! ************************************************************************** !

subroutine EOSWaterSatPresWagnerPruss(T, calculate_derivatives, &
                                      PS, dPS_dT, ierr)
  !
  ! Calculates the saturation pressure of water as a function of temperature
  ! based on Wagner W. and A. Pruss (1993) "International Equations for the
  ! Saturation Properties of Ordinary Water Substance. Revised According to
  ! the International Temperature Scale of 1990. Addendum to J. Phys. Chem.
  ! Ref. Data 16, 893 (1987)".
  !
  ! Author: Glenn Hammond
  ! Date: 11/07/16
  !
  implicit none

  PetscReal, intent(in) :: T ! temperature
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: PS, dPS_dT ! Saturation pres. and derivative
  PetscErrorCode, intent(inout) :: ierr

  PetscReal, parameter :: a(6) = [-7.85951783d0,1.84408259d0,-11.7866497d0, &
                                  22.6807411d0,-15.9618719d0,1.80122502d0]
  PetscReal, parameter :: Tc = 647.096d0 ! [K]
  PetscReal, parameter :: Pc = 22.064d6  ! [Pa]
  PetscReal :: tau
  PetscReal :: T_over_Tc
  PetscReal :: polynomial
  PetscReal :: dpolynomial_dtau, dtau_dT

  T_over_Tc = (T+273.15d0)/Tc
  tau = 1.d0 - T_over_Tc
  polynomial = a(1)*tau+a(2)*tau**1.5d0+a(3)*tau**3.d0+ &
               a(4)*tau**3.5d0+a(5)*tau**4.d0+a(6)*tau**7.5d0
  PS = Pc*exp(polynomial/T_over_Tc)
  if (calculate_derivatives) then
    dtau_dT = -1.d0/Tc
    dpolynomial_dtau = a(1)+1.5d0*a(2)*sqrt(tau)+3.d0*a(3)*tau*tau+ &
                       3.5d0*a(4)*tau**2.5d0+4.d0*a(5)*tau**3.d0+ &
                       7.5d0*a(6)*tau**6.5d0
    dPS_dT = PS*(dpolynomial_dtau/T_over_Tc*dtau_dT+polynomial*Tc)
  else
    dPS_dT = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterSatPresWagnerPruss

subroutine EOSWaterSaturationPressureIce(T, calculate_derivatives, &
                                      PS, dPS_dT, ierr)
  !
  ! Calculates the saturation pressure of water as a function of temperature
  ! above and below the freezing point of water based on Huang, J. (2018). 
  ! A simple accurate formula for calculating saturation vapor pressure of 
  ! water and ice. Journal of Applied Meteorology and Climatology, 57(6), 
  ! 1265-1272.
  !
  ! Author: Michael Nole
  ! Date: 12/07/2022
  !
  implicit none

  PetscReal, intent(in) :: T ! temperature
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: PS, dPS_dT ! Saturation pres. and derivative
  PetscErrorCode, intent(inout) :: ierr

  if (T > 0.d0) then
    PS = exp(34.494d0 - 4924.99d0/(T+237.1d0))/(T+105.d0)**1.57d0
    if (calculate_derivatives) then
      dPS_dT = -1.57d0*(T+105.d0)**(-2.57d0)*exp(34.494d0-4924.99d0/(T+237.1d0))
      dPS_dT = dPS_dT +(T+105.d0)**(-1.57d0)*exp(34.494d0-4924.99d0/ &
               (T+237.1d0)) * 4924.99d0*(T+237.1d0)**(-2.d0)
    endif
  else
    PS = exp(43.494d0 - 6545.8d0/(T+278.d0))/(T+868.d0)**2
    if (calculate_derivatives) then
      dPS_dT = -2.d0*(T+868.d0)**(-3.d0)*exp(43.494d0 - 6545.8d0/(T+278.d0))
      dPS_dT = dPS_dT + (T+868.d0)**(-2.d0)*exp(43.494d0 - 6545.8d0/ &
               (T+278.d0))*6545.8d0*(T+278.d0)**(-2.d0)
    endif
  endif
!print *, T
!print *, PS
end subroutine EOSWaterSaturationPressureIce

! ************************************************************************** !

subroutine EOSWaterDensityIFC67(t,p,calculate_derivatives,dw,dwmol, &
                                dwp,dwt,ierr,table_idxs)

!  This subroutine calculates water and steam-gas mixture properties.
!  The water and steam properties are valid in the range of:
!
!            0 < p < 165.4 * 10^5 pascals (165.4 bars)
!            0 < t < 350 centigrade (623.15 Kelvin)
!
!  The properties cover densities, enthalpies, internal energies,
!  and partial derivatives of these quanties with respect to
!  pressure and temperature.
!
!  For saturated fluid, it will also calculate water saturation
!  temperature for the specified pressure using Newton-Raphson and
!  the derivative dts/dp (=tsp) or Ps for a given temperature.
!
!  Ref.: International Formulation Committee of the Sixth International
!       Conference on Properties of Steam (1967).

  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal, save :: aa(0:22)
  PetscReal, save :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12

  PetscReal :: beta,beta2x,beta4,theta,theta2x,theta18,theta20
  PetscReal :: xx,yy,zz
  PetscReal :: u0,u1,u2,u3,u4,u5,u6,u7,u8,u9
!  PetscReal :: v0_1, v1_1, v2_1, v3_1, v4_1
!  PetscReal :: v1_2, v2_2, v3_2, v4_2, v20_2, v40_2
!  PetscReal :: v1_3, v2_3, v3_3, v4_3
!  PetscReal :: v1_4, v2_4, v3_4
!  PetscReal :: v1_5, v2_5
!  PetscReal :: v1_6
!  PetscReal :: term1,term2,term2t,term3,term3t,term3p,term4,term4t,term4p, &
!               term5,term5t,term5p,term6,term6t,term6p,term7,term7t,term7p
!  PetscReal :: dv2t,dv2p,dv3t
  PetscReal :: vr,ypt,zpt,zpp,vrpt,vrpp,cnv
  PetscReal :: tc1,pc1,vc1,utc1,upc1,vc1mol
  PetscReal :: d2z_dp2    ! 2nd derivative of z w.r.t. pressure
  PetscReal :: d2vr_dp2   ! 2nd derivative of vr w.r.t. pressure
  PetscReal :: d2wmol_dp2 ! 2nd derivative of density w.r.t. pressure
  PetscReal, parameter :: zero = 0.d0
  PetscReal, parameter :: one = 1.d0
  PetscReal, parameter :: two = 2.d0
  PetscReal, parameter :: three = 3.d0
  PetscReal, parameter :: four = 4.d0
  PetscReal, parameter :: five = 5.d0
  PetscReal, parameter :: six = 6.d0
  PetscReal, parameter :: ten = 10.d0

  data aa/ &
!-----data aa0,aa1,aa2,aa3/
        6.824687741d03,-5.422063673d02,-2.096666205d04, 3.941286787d04, &
!-----data aa4,aa5,aa6,aa7/
        -6.733277739d04, 9.902381028d04,-1.093911774d05, 8.590841667d04, &
!-----data aa8,aa9,aa10,aa11/
        -4.511168742d04, 1.418138926d04,-2.017271113d03, 7.982692717d00, &
!-----data aa12,aa13,aa14,aa15/
        -2.616571843d-2, 1.522411790d-3, 2.284279054d-2, 2.421647003d02, &
!-----data aa16,aa17,aa18,aa19/
        1.269716088d-10,2.074838328d-7, 2.174020350d-8, 1.105710498d-9, &
!-----data aa20,aa21,aa22/
        1.293441934d01, 1.308119072d-5, 6.047626338d-14/

  data a1,a2,a3,a4/ &
  8.438375405d-1, 5.362162162d-4, 1.720000000d00, 7.342278489d-2/
  data a5,a6,a7,a8/ &
  4.975858870d-2, 6.537154300d-1, 1.150000000d-6, 1.510800000d-5/
  data a9,a10,a11,a12/ &
  1.418800000d-1, 7.002753165d00, 2.995284926d-4, 2.040000000d-1/

  tc1 = H2O_CRITICAL_TEMPERATURE    ! K
  pc1 = H2O_CRITICAL_PRESSURE       ! Pa
  vc1 = 0.00317d0  ! m^3/kg
  utc1 = one/tc1   ! 1/C
  upc1 = one/pc1   ! 1/Pa
  vc1mol = vc1*FMWH2O ! m^3/kmol

  theta = (t+273.15d0)*utc1
  theta2x = theta*theta
  theta18 = theta**18.d0
  theta20 = theta18*theta2x

  beta = p*upc1
  beta2x = beta*beta
  beta4  = beta2x*beta2x

  yy = one-a1*theta2x-a2*theta**(-6.d0)
  xx = a3*yy*yy-two*(a4*theta-a5*beta)

!   Note: xx may become negative near the critical point-pcl.
  if (xx.gt.zero) then
    xx = sqrt(xx)
  else
    write(*,*) 'Warning: negative term in density (eos_water.F90:&
      &EOSWaterDensityIFC67):'
    write(*,*) 't= ',t,' p= ',p,' xx= ',xx
    ierr = 1
    xx = 1.d-6               !set arbitrarily
  end if
  zz = yy + xx
  u0 = -five/17.d0
  u1 = aa(11)*a5*zz**u0
  u2 = one/(a8+theta**11.d0)
  u3 = aa(17)+(two*aa(18)+three*aa(19)*beta)*beta
  u4 = one/(a7+theta18*theta)
  u5 = (a10+beta)**(-4.d0)
  u6 = a11-three*u5
  u7 = aa(20)*theta18*(a9+theta2x)
  u8 = aa(15)*(a6-theta)**9.d0

  vr = u1+aa(12)+theta*(aa(13)+aa(14)*theta)+u8*(a6-theta) &
        +aa(16)*u4-u2*u3-u6*u7+(three*aa(21)*(a12-theta) &
        +four*aa(22)*beta/theta20)*beta2x

  dwmol = one/(vr*vc1mol) ! kmol/m^3
  dw = one/(vr*vc1) ! kg/m^3

  ! ypt used for enthalpy even if derivative not calculated
  ypt = six*a2*theta**(-7.d0)-two*a1*theta

  !---calculate derivatives for water density
  if (calculate_derivatives) then
    zpt = ypt+(a3*yy*ypt-a4)/xx
    zpp = a5/xx
    u9 = u0*u1/zz
    vrpt = u9*zpt+aa(13)+two*aa(14)*theta-ten*u8 &
        -19.d0*aa(16)*u4*u4*theta18+11.d0*u2*u2*u3*theta**10.d0 &
        -aa(20)*u6*(18.d0*a9*theta18+20.d0*theta20)/theta &
        -(three*aa(21)+80.d0*aa(22)*beta/(theta20*theta))*beta2x

    vrpp = u9*zpp-u2*(two*aa(18)+six*aa(19)*beta)-12.d0*u7*u5/ &
        (a10+beta)+(six*aa(21)*(a12-theta)+12.d0*aa(22)*beta/ &
        theta20)*beta

    cnv = -one/(vc1mol*vr*vr)
    dwt = cnv*vrpt*utc1 ! kmol/m^3/C
    dwp = cnv*vrpp*upc1 ! kmol/m^3/Pa

    ! 2nd derivative w.r.t pressure
    d2z_dp2 = -a5*a5/xx/xx/xx

    d2vr_dp2 = u9*(u0-one)/zz*zpp*zpp + &
               u0*u1/zz*d2z_dp2 + &
               six*u2*aa(19) + &
               60.d0*u7*u5/(a10+beta)/(a10+beta) + &
               six*(a12 - theta) + &
               24.d0*aa(22)*beta/theta20

    d2wmol_dp2 = -cnv*upc1*upc1*(2.d0/vr*vrpp*vrpp -  d2vr_dp2) ! kmol/m^3/Pa^2

  else
    dwt = UNINITIALIZED_DOUBLE
    dwp = UNINITIALIZED_DOUBLE
    d2wmol_dp2 = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterDensityIFC67

subroutine EOSWaterDensityIF97(T,P,calculate_derivatives,dw,dwmol, &
                                dwp,dwt,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T   ! Temperature in centigrade
  PetscReal, intent(in) :: P   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal, parameter :: Tf = 273.15d0 !K
  PetscReal, parameter :: p_ref = 16.53d6 !Pa
  PetscReal, parameter :: T_ref = 1386.d0  ! K
  PetscReal, parameter :: R = 0.461526d0 ! kJ/kg-K
  PetscReal, parameter :: n_i(34) = [ 1.4632971213167d-1, &
          -8.4548187169114d-1, -0.37563603672040d1, 0.33855169168385d1, &
          -9.5791963387872d-1, 1.5772038513228d-1, -0.16616417199501d-1, &
          0.81214629983568d-3, 0.28319080123804d-3, -0.60706301565874d-3, &
          -0.18990068218419d-1, -0.32529748770505d-1, -0.21841717175414d-1, &
          -0.52838357969930d-4, -0.47184321073267d-3, -0.30001780793026d-3, &
          0.47661393906987d-4, -0.44141845330846d-5, -0.72694996297594d-15, &
          -0.31679644845054d-4, -0.28270797985312d-5, -0.85205128120103d-9, &
          -0.22425281908000d-5, -0.65171222895601d-6, -0.14341729937924d-12, &
          -0.40516996860117d-6, -0.12734301741641d-8, -0.17424871230634d-9, &
          -0.68762131295531d-18, 0.14478307828521d-19, 0.26335781662795d-22, &
          -0.11947622640071d-22, 0.18228094581404d-23, -0.93537087292458d-25 ]
  PetscInt, parameter :: I_i(34) = [0,0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3, &
                                    3,3,4,4,4,5,8,8,21,23,29,30,31,32]
  PetscInt, parameter :: J_i(34) = [-2,-1,0,1,2,3,4,5,-9,-7,-1,0,1,3,-3,0,1, &
                                     3,17,-4,0,6,-5,-2,10,-8,-11,-6,-29,-31, &
                                     -38,-39,-40,-41]
  PetscReal :: pi, tao, g_pi, T_temp
  PetscReal :: dg_pi_dT, dv_dt, dg_pi_dp, dv_dp


  if (Tf <= 623.15d0) then
    ! Region 1: Valid from 273.15K to 623.15 K, Ps(T) to 100MPa

    T_temp = T+Tf
    pi = P/p_ref
    tao = T_ref/T_temp

    g_pi = sum((-n_i*I_i*(7.1d0-pi)**(I_i-1))*(tao-1.222d0)**(J_i))
    dw = g_pi * pi*R*T_temp/P * 1.d3

    dw = 1.d0/dw
    dwmol = dw/FMWH2O

    if (calculate_derivatives) then
      dg_pi_dT = T_ref/(T_temp*T_temp) * sum(n_i*I_i*(7.1d0-pi)**(I_i-1)* &
                                             J_i*(tao-1.222d0)**(J_i-1))
      dv_dt = R*pi/P * (g_pi+T_temp*dg_pi_dT)
      dwt = -1.d3/FMWH2O*dw*dw * dv_dt

      dg_pi_dp = sum(n_i*I_i*(I_i-1)*(7.1d0-pi)**(I_i-2)* &
                     (tao-1.222d0)**(J_i)) / p_ref
      dv_dp = R*T_temp*(-g_pi*pi/(P*P) + (g_pi/p_ref + dg_pi_dp*pi)/P)
      dwp = -1.d3/FMWH2O*dw*dw *dv_dp
    else
      dwp = UNINITIALIZED_DOUBLE
      dwt = UNINITIALIZED_DOUBLE
    endif
  else
    ! Region 3: Valid in "wedge" >623.15K, >Ps(T), and 100MPa
    ! put into separate routine, to minimize weighing down this
    ! much more commonly used routine (Region 1)
    call EOSWaterDensityIF97Region3(T,P,calculate_derivatives,dw,dwmol, &
                                    dwp,dwt)
  end if

end subroutine EOSWaterDensityIF97

subroutine EOSWaterDensityIF97Region3(T,P,calculate_derivatives,dw,dwmol, &
                                     dwp,dwt)
  ! only ever called by EOSWaterDensityIF97()
  ! implements "backward equations" and "auxiliary equations" for
  ! specific volume as a funciton of p,T near the critical point
  ! in report IAPWS SR5-05(2016).

  implicit none
  PetscReal, intent(in) :: T   ! Temperature [C]
  PetscReal, intent(in) :: P   ! Pressure [Pa]
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode :: ierr

  PetscReal :: dumb, psat_64315, psat_62315, T_temp, nu
  PetscInt, parameter :: T3ab=1, T3cd=2, T3gh=3
  PetscInt, parameter :: T3ij=4, T3jk=5, T3mn=6
  PetscInt, parameter :: T3op=7, T3qu=8, T3rx=9
  PetscInt, parameter :: T3ef=12, T3uv=10, T3wx=11

  PetscReal, parameter :: Tf = 273.15d0 !K
  PetscReal, parameter :: R = 0.461526d0 ! kJ/kg-K

  ! Region 3:  Region 3: Valid in "wedge" >623.15K, >Ps(T), and 100MPa
  ! first, determine which of 24 sub-regions we fall into:

  T_temp = T+Tf

  ! 21.0434 MPa
  call EOSWaterSaturationPressureIF97(370.d0, .false., &
                                      psat_64315, dumb, ierr)
  ! 16.5292 MPa
  call EOSWaterSaturationPressureIF97(350.d0, .false., &
                                      psat_62315, dumb, ierr)

  ! ranges given in Table 2 and Table 10 (near critical point)
  nu = UNINITIALIZED_DOUBLE
  if (P <= 100.0D+6 .and. P > 40.0D+6) then
    if (T_temp > T3bdry(P,T3ab)) then
      nu = IF97_subregion_3b(P,T_temp)
    else
      nu = IF97_subregion_3a(P,T_temp)
    end if
  else if (P > 25.0D+6) then
    if (T_temp  > T3bdry(P,T3ef)) then
      nu = IF97_subregion_3f(P,T_temp)
    else if (T_temp > T3bdry(P,T3ab)) then
      nu = IF97_subregion_3e(P,T_temp)
    else if (T_temp > T3bdry(P,T3cd)) then
      nu = IF97_subregion_3d(P,T_temp)
    else
      nu = IF97_subregion_3c(P,T_temp)
    end if
  else if (P > 23.5D+6) then
    if (T_temp > T3bdry(P,T3jk)) then
      nu = IF97_subregion_3k(P,T_temp)
    else if (T_temp > T3bdry(P,T3ij)) then
      nu = IF97_subregion_3j(P,T_temp)
    else if (T_temp > T3bdry(P,T3ef)) then
      nu = IF97_subregion_3i(P,T_temp)
    else if (T_temp > T3bdry(P,T3gh)) then
      nu = IF97_subregion_3h(P,T_temp)
    else if (T_temp > T3bdry(P,T3cd)) then
      nu = IF97_subregion_3g(P,T_temp)
    else
      nu = IF97_subregion_3c(P,T_temp)
    end if
  else if (P > 23.0D+6) then
    if (T_temp > T3bdry(P,T3jk)) then
      nu = IF97_subregion_3k(P,T_temp)
    else if (T_temp > T3bdry(P,T3ij)) then
      nu = IF97_subregion_3j(P,T_temp)
    else if (T_temp > T3bdry(P,T3ef)) then
      nu = IF97_subregion_3i(P,T_temp)
    else if (T_temp > T3bdry(P,T3gh)) then
      nu = IF97_subregion_3h(P,T_temp)
    else if (T_temp > T3bdry(P,T3cd)) then
      nu = IF97_subregion_3l(P,T_temp)
    else
      nu = IF97_subregion_3c(P,T_temp)
    end if
  else if (P > 22.5D+6) then
    if (T_temp > T3bdry(P,T3jk)) then
      nu = IF97_subregion_3k(P,T_temp)
    else if (T_temp > T3bdry(P,T3ij)) then
      nu = IF97_subregion_3j(P,T_temp)
    else if (T_temp > T3bdry(P,T3op)) then
      nu = IF97_subregion_3p(P,T_temp)
    else if (T_temp > T3bdry(P,T3ef)) then
      nu = IF97_subregion_3o(P,T_temp)
    else if (T_temp > T3bdry(P,T3mn)) then
      nu = IF97_subregion_3n(P,T_temp)
    else if (T_temp > T3bdry(P,T3gh)) then
      nu = IF97_subregion_3m(P,T_temp)
    else if (T_temp > T3bdry(P,T3cd)) then
      nu = IF97_subregion_3l(P,T_temp)
    else
      nu = IF97_subregion_3c(P,T_temp)
    end if
  else if (P > 21.11D+6) then
    if (T_temp > T3bdry(P,T3jk)) then
      nu = IF97_subregion_3k(P,T_temp)
    else if (T_temp > T3bdry(P,T3rx)) then
      nu = IF97_subregion_3r(P,T_temp)
    else if (T_temp > T3bdry(P,T3wx)) then
      nu = IF97_subregion_3x(P,T_Temp)
    else if (T_temp > T3bdry(P,T3ef)) then
      nu = IF97_subregion_3w(P,T_Temp)
    else if (T_temp > T3bdry(P,T3uv)) then
      nu = IF97_subregion_3v(P,T_Temp)
    else if (T_temp > T3bdry(P,T3qu)) then
      nu = IF97_subregion_3u(P,T_Temp)
    else if (T_temp > T3bdry(P,T3cd)) then
      nu = IF97_subregion_3q(P,T_Temp)
    else
      nu = IF97_subregion_3c(P,T_temp)
    end if
  else if (P > 22.064D+6) then ! P_crit
    if (T_temp > T3bdry(P,T3jk)) then
      nu = IF97_subregion_3k(P,T_temp)
    else if (T_temp > T3bdry(P,T3rx)) then
      nu = IF97_subregion_3r(P,T_temp)
    else if (T_temp > T3bdry(P,T3wx)) then
      nu = IF97_subregion_3x(P,T_Temp)
    else if (T_temp > T3bdry(P,T3ef)) then
      nu = IF97_subregion_3z(P,T_Temp)
    else if (T_temp > T3bdry(P,T3uv)) then
      nu = IF97_subregion_3y(P,T_Temp)
    else if (T_temp > T3bdry(P,T3qu)) then
      nu = IF97_subregion_3u(P,T_Temp)
    else if (T_temp > T3bdry(P,T3cd)) then
      nu = IF97_subregion_3q(P,T_Temp)
    else
      nu = IF97_subregion_3c(P,T_temp)
    end if
  else if (P > 21.93161551d+6) then
    if (T_temp > T3bdry(P,T3jk)) then
      nu = IF97_subregion_3k(P,T_temp)
    else if (T_temp > T3bdry(P,T3rx)) then
      nu = IF97_subregion_3r(P,T_temp)
    else if (T_temp > T3bdry(P,T3wx)) then
      nu = IF97_subregion_3x(P,T_Temp)
    else if (T_temp > IF97SaturationTemperature(P)) then
      nu = IF97_subregion_3z(P,T_Temp)
    else if (T_temp > T3bdry(P,T3uv)) then
      nu = IF97_subregion_3y(P,T_Temp)
    else if (T_temp > T3bdry(P,T3qu)) then
      nu = IF97_subregion_3u(P,T_Temp)
    else if (T_temp > T3bdry(P,T3cd)) then
      nu = IF97_subregion_3q(P,T_Temp)
    else
      nu = IF97_subregion_3c(P,T_temp)
    end if
  else if (P > 21.90096265d+6) then
    if (T_temp > T3bdry(P,T3jk)) then
      nu = IF97_subregion_3k(P,T_temp)
    else if (T_temp > T3bdry(P,T3rx)) then
      nu = IF97_subregion_3r(P,T_temp)
    else if (T_temp > T3bdry(P,T3wx)) then
      nu = IF97_subregion_3x(P,T_Temp)
    else if (T_temp > IF97SaturationTemperature(P)) then
      nu = IF97_subregion_3z(P,T_Temp)
    else if (T_temp > T3bdry(P,T3qu)) then
      nu = IF97_subregion_3u(P,T_Temp)
    else if (T_temp > T3bdry(P,T3cd)) then
      nu = IF97_subregion_3q(P,T_Temp)
    else
      nu = IF97_subregion_3c(P,T_temp)
    end if
  else if (P > psat_64315) then ! 21.0434 MPa
    if (T_temp > T3bdry(P,T3jk)) then
      nu = IF97_subregion_3k(P,T_temp)
    else if (T_temp > T3bdry(P,T3rx)) then
      nu = IF97_subregion_3r(P,T_temp)
    else if (T_temp > IF97SaturationTemperature(P)) then
      nu = IF97_subregion_3x(P,T_temp)
    else if (T_temp > T3bdry(P,T3qu)) then
      nu = IF97_subregion_3u(P,T_temp)
    else if (T_temp > T3bdry(P,T3cd)) then
      nu = IF97_subregion_3q(P,T_temp)
    else
      nu = IF97_subregion_3c(P,T_temp)
    end if
  else if (P > 20.5D+6) then
    if (T_temp > T3bdry(P,T3jk)) then
      nu = IF97_subregion_3k(P,T_temp)
    else if (T_temp > IF97SaturationTemperature(p)) then
      nu = IF97_subregion_3r(P,T_temp)
    else if (T_temp > T3bdry(P,T3cd)) then
      nu = IF97_subregion_3s(P,T_temp)
    else
      nu = IF97_subregion_3c(P,T_temp)
    end if
  else if (P > 1.900881189173929D+1) then
    if (T_temp > IF97SaturationTemperature(P)) then
      nu = IF97_subregion_3t(P,T_temp)
    else if (T_temp > T3bdry(P,T3cd)) then
      nu = IF97_subregion_3s(P,T_temp)
    else
      nu = IF97_subregion_3c(P,T_temp)
    end if
  else if (P > psat_62315) then ! 16.5292 MPa
    if (T_temp >= IF97SaturationTemperature(P)) then
      nu = IF97_subregion_3t(P,T_temp)
    else
      nu = IF97_subregion_3c(P,T_temp)
    end if
  else
    stop 'IF97 Region 3 ERROR: either > 100 MPa or < 16.5292 MPa'
  end if

  dw = 1.d0 / (nu * R * T_temp / P * 1.0D+3)

  if (calculate_derivatives) then
    stop 'IF97 region 3 water density derivative not implemented yet'
  else
    dwp = UNINITIALIZED_DOUBLE
    dwt = UNINITIALIZED_DOUBLE
    dwmol = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterDensityIF97Region3

function T3bdry(P,idx)
  implicit none
  PetscReal, parameter :: n_i(11,5) = reshape([ &
        0.154793642129415D4, -0.187661219490113D3, 0.213144632222113D2, & !T3ab
       -0.191887498864292D4, 0.918419702359447D3, &
        0.585276966696349D3, 0.278233532206915D1, -0.127283549295878D-1, & !T3cd
        0.159090746562729D-3, 0.0d0, &
       -0.249284240900418D5, 0.428143584791546D4, -0.269029173140130D3, & !T3gh
        0.751608051114157D1, -0.787105249910383D-1, &
        0.584814781649163D3, -0.616179320924617D0, 0.260763050899562D0, & !T3ij
       -0.587071076864459D-2, 0.515308185433082D-4, &
        0.617229772068439D3, -0.770600270141675D1, 0.697072596851896D0, & !T3jk
       -0.157391839848015D-1, 0.137897492684194D-3, &
        0.535339483742384D3, 0.761978122720128D1, -0.158365725441648D0, & !T3mn
        0.192871054508108D-2, 0.0d0, &
        0.969461372400213D3, -0.332500170441278D3, 0.642859598466067D2, & !T3op
        0.773845935768222D3, -0.152313732937084D4, &
        0.565603648239126D3, 0.529062258221222D1, -0.102020639611016D0, & !T3qu
        0.122240301070145D-2, 0.0d0, &
        0.584561202520006D3, -0.102961025163669D1, 0.243293362700452D0, & !T3rx
       -0.294905044740799D-2, 0.0d0, &
        0.528199646263062D3, 0.890579602135307D1, -0.222814134903755D0, & !T3uv
        0.286791682263697D-2, 0.0d0, &
        0.728052609145380D1, 0.973505869861952D2, 0.147370491183191D2, & !T3wx
        0.329196213998375D3, 0.873371668682417D3],shape=[11,5],order=[2,1])
  PetscInt, parameter :: I_i(11,5) = reshape([0,1,2,-1,-2, 0,1,2,3,0, 0,1,2,3,4, &
       0,1,2,3,4, 0,1,2,3,4, 0,1,2,3,0, 0,1,2,-1,-2, 0,1,2,3,0, 0,1,2,3,0, 0,1,2,3,0, &
       0,1,2,-1,-2],shape=[11,5],order=[2,1])
  PetscReal, intent(in) :: P
  PetscInt, intent(in) :: idx
  PetscReal :: T3bdry, pi

  pi = P * 1.0D-6
  select case(idx)
  case(1,7,11)
    T3bdry = sum(n_i(idx,:) * log(pi) ** I_i(idx,:)) ! eqn 2
  case(12)
    T3bdry = 3.727888004d0 * (pi - 22.064d0) + 647.096d0 ! eqn 3
  case(2:6,8:10)
    T3bdry = sum(n_i(idx,:) * pi ** I_i(idx,:)) ! eqn 1
  case default
    stop "ERROR: IF97 T3bdry(). Invalid select case"
  end select
end function T3bdry

function IF97_subregion_3a(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(30) = [0.110879558823853D-2, &
          0.572616740810616D3, -0.767051948380852D5, -0.253321069529674D-1, &
          0.628008049345689D4, 0.234105654131876D6, 0.216867826045856D0, &
          -0.156237904341963D3, -0.269893956176613D5, -0.180407100085505D-3, &
          0.116732227668261D-2, 0.266987040856040D2, 0.282776617243286D5, &
          -0.242431520029523D4, 0.435217323022733D-3, -0.122494831387441D-1, &
          0.179357604019989D1, 0.442729521058314D2, -0.593223489018342D-2, &
          0.453186261685774D0, 0.135825703129140D1, 0.408748415856745D-1, &
          0.474686397863312D0, 0.118646814997915D1, 0.546987265727549D0, &
          0.195266770452643D0, -0.502268790869663D-1, -0.369645308193377D0, &
          0.633828037528420D-2, 0.797441793901017D-1]
  PetscInt, parameter :: I_i(30) = [-12,-12,-12,-10,-10,-10,-8,-8,-8,-6,-5,-5,-5, &
                                    -4,-3,-3,-3,-3,-2,-2,-2,-1,-1,-1,0,0,1,1,2, &
                                    2]
  PetscInt, parameter :: J_i(30) = [5,10,12,5,10,12,5,8,10,1,1,5,10,8,0,1,3,6,0, &
                                    2,3,0,1,2,0,1,0,2,0,2]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0024d0, p_star = 100.d0, T_star = 760.d0
  PetscReal, parameter :: a = 0.085d0, b = 0.817d0
  PetscReal, parameter :: c = 1.0d0, d = 1.00d0, e = 1.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3a

function IF97_subregion_3b(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(32) = [-0.827670470003621D-1, &
          0.416887126010565D2, 0.483651982197059D-1, -0.291032084950276D5, &
          -0.111422582236948D3, -0.202300083904014D-1, 0.294002509338515D3, &
          0.140244997609658D3, -0.344384158811459D3, 0.361182452612149D3, &
          -0.140699677420738D4, -0.202023902676481D-2, 0.171346792457471D3, &
          -0.425597804058632D1, 0.691346085000334D-5, 0.151140509678925D-2, &
          -0.416375290166236D-1, -0.413754957011042D2, -0.506673295721637D2, &
          -0.572212965569023D-3, 0.608817368401785D1, 0.239600660256161D2, &
          0.122261479925384D-1, 0.216356057692938D1, 0.398198903368642D0, &
          -0.116892827834085D0, -0.102845919373532D0, -0.492676637589284D0, &
          0.655540456406790D-1, -0.240462535078530D0, -0.269798180310075D-1, &
          0.128369435967012D0]
  PetscInt, parameter :: I_i(32) = [-12,-12,-10,-10,-8,-6,-6,-6,-5,-5,-5,-4,-4, &
                                    -4,-3,-3,-3,-3,-3,-2,-2,-2,-1,-1,0,0,1,1,2, &
                                    3,4,4]
  PetscInt, parameter :: J_i(32) = [10,12,8,14,8,5,6,8,5,8,10,2,4,5,0,1,2,3,5,0, &
                                    2,5,0,2,0,1,0,2,0,2,0,1]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0041d0, p_star = 100.d0, T_star = 860.d0
  PetscReal, parameter :: a = 0.280d0, b = 0.779d0
  PetscReal, parameter :: c = 1.0d0, d = 1.00d0, e = 1.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3b

function IF97_subregion_3c(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(35) = [0.311967788763030D1, &
          0.276713458847564D5, 0.322583103403269D8, -0.342416065095363D3, &
          -0.899732529907377D6, -0.793892049821251D8, 0.953193003217388D2, &
          0.229784742345072D4, 0.175336675322499D6, 0.791214365222792D7, &
          0.319933345844209D-4, -0.659508863555767D2, -0.833426563212851D6, &
          0.645734680583292D-1, -0.382031020570813D7, 0.406398848470079D-4, &
          0.310327498492008D2, -0.892996718483724D-3, 0.234604891591616D3, &
          0.377515668966951D4, 0.158646812591361D-1, 0.707906336241843D0, &
          0.126016225146570D2, 0.736143655772152D0, 0.676544268999101D0, &
          -0.178100588189137D2, -0.156531975531713D0, 0.117707430048158D2, &
          0.840143653860447D-1, -0.186442467471949D0, -0.440170203949645D2, &
          0.123290423502494D7, -0.240650039730845D-1, -0.107077716660869D7, &
          0.438319858566475D-1]
  PetscInt, parameter :: I_i(35) = [-12,-12,-12,-10,-10,-10,-8,-8,-8,-6,-5,-5,-5, &
                                    -4,-4,-3,-3,-2,-2,-2,-1,-1,-1,0,0,0,1,1,2,2, &
                                    2,2,3,3,8]
  PetscInt, parameter :: J_i(35) = [6,8,10,6,8,10,5,6,7,8,1,4,7,2,8,0,3,0,4,5,0, &
                                    1,2,0,1,2,0,2,0,1,3,7,0,7,1]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0022d0, p_star = 40.d0, T_star = 690.d0
  PetscReal, parameter :: a = 0.259d0, b = 0.903d0
  PetscReal, parameter :: c = 1.0d0, d = 1.00d0, e = 1.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3c

function IF97_subregion_3d(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(38) = [-0.452484847171645D-9, &
          0.315210389538801D-4, -0.214991352047545D-2, 0.508058874808345D3, &
          -0.127123036845932D8, 0.115371133120497D13, -0.197805728776273D-15, &
          0.241554806033972D-10, -0.156481703640525D-5, 0.277211346836625D-2, &
          -0.203578994462286D2, 0.144369489909053D7, -0.411254217946539D11, &
          0.623449786243773D-5, -0.221774281146038D2, -0.689315087933158D5, &
          -0.195419525060713D8, 0.316373510564015D4, 0.224040754426988D7, &
          -0.436701347922356D-5, -0.404213852833996D-3, -0.348153203414663D3, &
          -0.385294213555289D6, 0.135203700099403D-6, 0.134648383271089D-3, &
          0.125031835351736D6, 0.968123678455841D-1, 0.225660517512438D3, &
          -0.190102435341872D-3, -0.299628410819229D-1, 0.500833915372121D-2, &
          0.387842482998411D0, -0.138535367777182D4, 0.870745245971773D0, &
          0.171946252068742D1, -0.326650121426383D-1, 0.498044171727877D4, &
          0.551478022765087D-2]
  PetscInt, parameter :: I_i(38) = [-12,-12,-12,-12,-12,-12,-10,-10,-10,-10,-10, &
                                    -10,-10,-8,-8,-8,-8,-6,-6,-5,-5,-5,-5,-4,-4, &
                                    -4,-3,-3,-2,-2,-1,-1,-1,0,0,1,1,3]
  PetscInt, parameter :: J_i(38) = [4,6,7,10,12,16,0,2,4,6,8,10,14,3,7,8,10,6,8, &
                                    1,2,5,7,0,1,7,2,4,0,1,0,1,5,0,2,0,6,0]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0029d0, p_star = 40.d0, T_star = 690.d0
  PetscReal, parameter :: a = 0.559d0, b = 0.939d0
  PetscReal, parameter :: c = 1.0d0, d = 1.00d0, e = 4.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3d

function IF97_subregion_3e(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(29) = [0.715815808404721D9, &
          -0.114328360753449D12, 0.376531002015720D-11, -0.903983668691157D-4, &
          0.665695908836252D6, 0.535364174960127D10, 0.794977402335603D11, &
          0.922230563421437D2, -0.142586073991215D6, -0.111796381424162D7, &
          0.896121629640760D4, -0.669989239070491D4, 0.451242538486834D-2, &
          -0.339731325977713D2, -0.120523111552278D1, 0.475992667717124D5, &
          -0.266627750390341D6, -0.153314954386524D-3, 0.305638404828265D0, &
          0.123654999499486D3, -0.104390794213011D4, -0.157496516174308D-1, &
          0.685331118940253D0, 0.178373462873903D1, -0.544674124878910D0, &
          0.204529931318843D4, -0.228342359328752D5, 0.413197481515899D0, &
          -0.341931835910405D2]
  PetscInt, parameter :: I_i(29) = [-12,-12,-10,-10,-10,-10,-10,-8,-8,-8,-6,-5, &
                                    -4,-4,-3,-3,-3,-2,-2,-2,-2,-1,0,0,1,1,1,2,2 &
                                    ]
  PetscInt, parameter :: J_i(29) = [14,16,3,6,10,14,16,7,8,10,6,6,2,4,2,6,7,0,1, &
                                    3,4,0,0,1,0,4,6,0,2]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0032d0, p_star = 40.d0, T_star = 710.d0
  PetscReal, parameter :: a = 0.587d0, b = 0.918d0
  PetscReal, parameter :: c = 1.0d0, d = 1.00d0, e = 1.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3e

function IF97_subregion_3f(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(42) = [-0.251756547792325D-7, &
          0.601307193668763D-5, -0.100615977450049D-2, 0.999969140252192D0, &
          0.214107759236486D1, -0.165175571959086D2, -0.141987303638727D-2, &
          0.269251915156554D1, 0.349741815858722D2, -0.300208695771783D2, &
          -0.131546288252539D1, -0.839091277286169D1, 0.181545608337015D-9, &
          -0.591099206478909D-3, 0.152115067087106D1, 0.252956470663225D-4, &
          0.100726265203786D-14, -0.149774533860650D1, -0.793940970562969D-9, &
          -0.150290891264717D-3, 0.151205531275133D1, 0.470942606221652D-5, &
          0.195049710391712D-12, -0.911627886266077D-8, 0.604374640201265D-3, &
          -0.225132933900136D-15, 0.610916973582981D-11, -0.303063908043404D-6, &
          -0.137796070798409D-4, -0.919296736666106D-3, 0.639288223132545D-9, &
          0.753259479898699D-6, -0.400321478682929D-12, 0.756140294351614D-8, &
          -0.912082054034891D-11, -0.237612381140539D-7, 0.269586010591874D-4, &
          -0.732828135157839D-10, 0.241995578306660D-9, -0.405735532730322D-3, &
          0.189424143498011D-9, -0.486632965074563D-9]
  PetscInt, parameter :: I_i(42) = [0,0,0,0,0,0,1,1,1,1,2,2,3,3,3,4,5,5,6,7,7,10, &
                                    12,12,12,14,14,14,14,14,16,16,18,18,20,20,20, &
                                    22,24,24,28,32]
  PetscInt, parameter :: J_i(42) = [-3,-2,-1,0,1,2,-1,1,2,3,0,1,-5,-2,0,-3,-8,1, &
                                    -6,-4,1,-6,-10,-8,-4,-12,-10,-8,-6,-4,-10,-8, &
                                    -12,-10,-12,-10,-6,-12,-12,-4,-12,-12]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0064d0, p_star = 40.d0, T_star = 730.d0
  PetscReal, parameter :: a = 0.587d0, b = 0.891d0
  PetscReal, parameter :: c = 0.5d0, d = 1.00d0, e = 4.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3f

function IF97_subregion_3g(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(38) = [0.412209020652996D-4, &
          -0.114987238280587D7, 0.948180885032080D10, -0.195788865718971D18, &
          0.496250704871300D25, -0.105549884548496D29, -0.758642165988278D12, &
          -0.922172769596101D23, 0.725379072059348D30, -0.617718249205859D2, &
          0.107555033344858D5, -0.379545802336487D8, 0.228646846221831D12, &
          -0.499741093010619D7, -0.280214310054101D31, 0.104915406769586D7, &
          0.613754229168619D28, 0.802056715528378D32, -0.298617819828065D8, &
          -0.910782540134681D2, 0.135033227281565D6, -0.712949383408211D19, &
          -0.104578785289542D37, 0.304331584444093D2, 0.593250797959445D10, &
          -0.364174062110798D28, 0.921791403532461D0, -0.337693609657471D0, &
          -0.724644143758508D2, -0.110480239272601D0, 0.536516031875059D1, &
          -0.291441872156205D4, 0.616338176535305D40, -0.120889175861180D39, &
          0.818396024524612D23, 0.940781944835829D9, -0.367279669545448D5, &
          -0.837513931798655D16]
  PetscInt, parameter :: I_i(38) = [-12,-12,-12,-12,-12,-12,-10,-10,-10,-8,-8,-8, &
                                    -8,-6,-6,-5,-5,-4,-3,-2,-2,-2,-2,-1,-1,-1,0, &
                                    0,0,1,1,1,3,5,6,8,10,10]
  PetscInt, parameter :: J_i(38) = [7,12,14,18,22,24,14,20,24,7,8,10,12,8,22,7, &
                                    20,22,7,3,5,14,24,2,8,18,0,1,2,0,1,3,24,22, &
                                    12,3,0,6]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0027d0, p_star = 25.d0, T_star = 660.d0
  PetscReal, parameter :: a = 0.872d0, b = 0.971d0
  PetscReal, parameter :: c = 1.0d0, d = 1.00d0, e = 4.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3g

function IF97_subregion_3h(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(29) = [0.561379678887577D-1, &
          0.774135421587083D10, 0.111482975877938D-8, -0.143987128208183D-2, &
          0.193696558764920D4, -0.605971823585005D9, 0.171951568124337D14, &
          -0.185461154985145D17, 0.387851168078010D-16, -0.395464327846105D-13, &
          -0.170875935679023D3, -0.212010620701220D4, 0.177683337348191D8, &
          0.110177443629575D2, -0.234396091693313D6, -0.656174421999594D7, &
          0.156362212977396D-4, -0.212946257021400D1, 0.135249306374858D2, &
          0.177189164145813D0, 0.139499167345464D4, -0.703670932036388D-2, &
          -0.152011044389648D0, 0.981916922991113D-4, 0.147199658618076D-2, &
          0.202618487025578D2, 0.899345518944240D0, -0.211346402240858D0, &
          0.249971752957491D2]
  PetscInt, parameter :: I_i(29) = [-12,-12,-10,-10,-10,-10,-10,-10,-8,-8,-8,-8, &
                                    -8,-6,-6,-6,-5,-5,-5,-4,-4,-3,-3,-2,-1,-1,0, &
                                    1,1]
  PetscInt, parameter :: J_i(29) = [8,12,4,6,8,10,14,16,0,1,6,7,8,4,6,8,2,3,4,2, &
                                    4,1,2,0,0,2,0,0,2]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0032d0, p_star = 25.d0, T_star = 660.d0
  PetscReal, parameter :: a = 0.898d0, b = 0.983d0
  PetscReal, parameter :: c = 1.0d0, d = 1.00d0, e = 4.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3h

function IF97_subregion_3i(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(42) = [0.106905684359136D1, &
          -0.148620857922333D1, 0.259862256980408D15, -0.446352055678749D-11, &
          -0.566620757170032D-6, -0.235302885736849D-2, -0.269226321968839D0, &
          0.922024992944392D1, 0.357633505503772D-11, -0.173942565562222D2, &
          0.700681785556229D-5, -0.267050351075768D-3, -0.231779669675624D1, &
          -0.753533046979752D-12, 0.481337131452891D1, -0.223286270422356D22, &
          -0.118746004987383D-4, 0.646412934136496D-2, -0.410588536330937D-9, &
          0.422739537057241D20, 0.313698180473812D-12, 0.164395334345040D-23, &
          -0.339823323754373D-5, -0.135268639905021D-1, -0.723252514211625D-14, &
          0.184386437538366D-8, -0.463959533752385D-1, -0.992263100376750D14, &
          0.688169154439335D-16, -0.222620998452197D-10, -0.540843018624083D-7, &
          0.345570606200257D-2, 0.422275800304086D11, -0.126974478770487D-14, &
          0.927237985153679D-9, 0.612670812016489D-13, -0.722693924063497D-11, &
          -0.383669502636822D-3, 0.374684572410204D-3, -0.931976897511086D5, &
          -0.247690616026922D-1, 0.658110546759474D2]
  PetscInt, parameter :: I_i(42) = [0,0,0,1,1,1,1,2,3,3,4,4,4,5,5,5,7,7,8,8,10, &
                                    12,12,12,14,14,14,14,18,18,18,18,18,20,20,22, &
                                    24,24,32,32,36,36]
  PetscInt, parameter :: J_i(42) = [0,1,10,-4,-2,-1,0,0,-5,0,-3,-2,-1,-6,-1,12, &
                                    -4,-3,-6,10,-8,-12,-6,-4,-10,-8,-4,5,-12,-10, &
                                    -8,-6,2,-12,-10,-12,-12,-8,-10,-5,-10,-8]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0041d0, p_star = 25.d0, T_star = 660.d0
  PetscReal, parameter :: a = 0.910d0, b = 0.984d0
  PetscReal, parameter :: c = 0.5d0, d = 1.00d0, e = 4.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3i

function IF97_subregion_3j(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(29) = [-0.111371317395540D-3, &
          0.100342892423685D1, 0.530615581928979D1, 0.179058760078792D-5, &
          -0.728541958464774D-3, -0.187576133371704D2, 0.199060874071849D-2, &
          0.243574755377290D2, -0.177040785499444D-3, -0.259680385227130D-2, &
          -0.198704578406823D3, 0.738627790224287D-4, -0.236264692844138D-2, &
          -0.161023121314333D1, 0.622322971786473D4, -0.960754116701669D-8, &
          -0.510572269720488D-10, 0.767373781404211D-2, 0.663855469485254D-14, &
          -0.717590735526745D-9, 0.146564542926508D-4, 0.309029474277013D-11, &
          -0.464216300971708D-15, -0.390499637961161D-13, -0.236716126781431D-9, &
          0.454652854268717D-11, -0.422271787482497D-2, 0.283911742354706D-10, &
          0.270929002720228D1]
  PetscInt, parameter :: I_i(29) = [0,0,0,1,1,1,2,2,3,4,4,5,5,5,6,10,12,12,14,14, &
                                    14,16,18,20,20,24,24,28,28]
  PetscInt, parameter :: J_i(29) = [-1,0,1,-2,-1,1,-1,1,-2,-2,2,-3,-2,0,3,-6,-8, &
                                    -3,-10,-8,-5,-10,-12,-12,-10,-12,-6,-12,-5]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0054d0, p_star = 25.d0, T_star = 670.d0
  PetscReal, parameter :: a = 0.875d0, b = 0.964d0
  PetscReal, parameter :: c = 0.5d0, d = 1.00d0, e = 4.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3j

function IF97_subregion_3k(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(34) = [-0.401215699576099D9, &
          0.484501478318406D11, 0.394721471363678D-14, 0.372629967374147D5, &
          -0.369794374168666D-29, -0.380436407012452D-14, 0.475361629970233D-6, &
          -0.879148916140706D-3, 0.844317863844331D0, 0.122433162656600D2, &
          -0.104529634830279D3, 0.589702771277429D3, -0.291026851164444D14, &
          0.170343072841850D-5, -0.277617606975748D-3, -0.344709605486686D1, &
          0.221333862447095D2, -0.194646110037079D3, 0.808354639772825D-15, &
          -0.180845209145470D-10, -0.696664158132412D-5, -0.181057560300994D-2, &
          0.255830298579027D1, 0.328913873658481D4, -0.173270241249904D-18, &
          -0.661876792558034D-6, -0.395688923421250D-2, 0.604203299819132D-17, &
          -0.400879935920517D-13, 0.160751107464958D-8, 0.383719409025556D-4, &
          -0.649565446702457D-14, -0.149095328506000D-11, 0.541449377329581D-8]
  PetscInt, parameter :: I_i(34) = [-2,-2,-1,-1,0,0,0,0,0,0,0,0,0,1,1,1,1,1,2,2, &
                                    2,2,2,2,5,5,5,6,6,6,6,8,10,12]
  PetscInt, parameter :: J_i(34) = [10,12,-5,6,-12,-6,-2,-1,0,1,2,3,14,-3,-2,0, &
                                    1,2,-8,-6,-3,-2,0,4,-12,-6,-3,-12,-10,-8,-5, &
                                    -12,-12,-10]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0077d0, p_star = 25.d0, T_star = 680.d0
  PetscReal, parameter :: a = 0.802d0, b = 0.935d0
  PetscReal, parameter :: c = 1.0d0, d = 1.00d0, e = 1.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3k

function IF97_subregion_3l(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(43) = [0.260702058647537D10, &
          -0.188277213604704D15, 0.554923870289667D19, -0.758966946387758D23, &
          0.413865186848908D27, -0.815038000738060D12, -0.381458260489955D33, &
          -0.123239564600519D-1, 0.226095631437174D8, -0.495017809506720D12, &
          0.529482996422863D16, -0.444359478746295D23, 0.521635864527315D35, &
          -0.487095672740742D55, -0.714430209937547D6, 0.127868634615495D0, &
          -0.100752127917598D2, 0.777451437960990D7, -0.108105480796471D25, &
          -0.357578581169659D-5, -0.212857169423484D1, 0.270706111085238D30, &
          -0.695953622348829D33, 0.110609027472280D0, 0.721559163361354D2, &
          -0.306367307532219D15, 0.265839618885530D-4, 0.253392392889754D-1, &
          -0.214443041836579D3, 0.937846601489667D0, 0.223184043101700D1, &
          0.338401222509191D2, 0.494237237179718D21, -0.198068404154428D0, &
          -0.141415349881140D31, -0.993862421613651D2, 0.125070534142731D3, &
          -0.996473529004439D3, 0.473137909872765D5, 0.116662121219322D33, &
          -0.315874976271533D16, -0.445703369196945D33, 0.642794932373694D33]
  PetscInt, parameter :: I_i(43) = [-12,-12,-12,-12,-12,-10,-10,-8,-8,-8,-8,-8, &
                                    -8,-8,-6,-5,-5,-4,-4,-3,-3,-3,-3,-2,-2,-2,-1, &
                                    -1,-1,0,0,0,0,1,1,2,4,5,5,6,10,10,14]
  PetscInt, parameter :: J_i(43) = [14,16,18,20,22,14,24,6,10,12,14,18,24,36,8, &
                                    4,5,7,16,1,3,18,20,2,3,10,0,1,3,0,1,2,12,0, &
                                    16,1,0,0,1,14,4,12,10]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0026d0, p_star = 24.d0, T_star = 650.d0
  PetscReal, parameter :: a = 0.908d0, b = 0.989d0
  PetscReal, parameter :: c = 1.0d0, d = 1.00d0, e = 4.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3l

function IF97_subregion_3m(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(40) = [0.811384363481847D0, &
          -0.568199310990094D4, -0.178657198172556D11, 0.795537657613427D32, &
          -0.814568209346872D5, -0.659774567602874D8, -0.152861148659302D11, &
          -0.560165667510446D12, 0.458384828593949D6, -0.385754000383848D14, &
          0.453735800004273D8, 0.939454935735563D12, 0.266572856432938D28, &
          -0.547578313899097D10, 0.200725701112386D15, 0.185007245563239D13, &
          0.185135446828337D9, -0.170451090076385D12, 0.157890366037614D15, &
          -0.202530509748774D16, 0.368193926183570D60, 0.170215539458936D18, &
          0.639234909918741D42, -0.821698160721956D15, -0.795260241872306D24, &
          0.233415869478510D18, -0.600079934586803D23, 0.594584382273384D25, &
          0.189461279349492D40, -0.810093428842645D46, 0.188813911076809D22, &
          0.111052244098768D36, 0.291133958602503D46, -0.329421923951460D22, &
          -0.137570282536696D26, 0.181508996303902D28, -0.346865122768353D30, &
          -0.211961148774260D38, -0.128617899887675D49, 0.479817895699239D65]
  PetscInt, parameter :: I_i(40) = [0,3,8,20,1,3,4,5,1,6,2,4,14,2,5,3,0,1,1,1,28, &
                                    2,16,0,5,0,3,4,12,16,1,8,14,0,2,3,4,8,14,24 &
                                    ]
  PetscInt, parameter :: J_i(40) = [0,0,0,2,5,5,5,5,6,6,7,8,8,10,10,12,14,14,18, &
                                    20,20,22,22,24,24,28,28,28,28,28,32,32,32,36, &
                                    36,36,36,36,36,36]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0028d0, p_star = 23.d0, T_star = 650.d0
  PetscReal, parameter :: a = 1.000d0, b = 0.997d0
  PetscReal, parameter :: c = 1.0d0, d = 0.25d0, e = 1.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3m

function IF97_subregion_3n(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(39) = [0.280967799943151D-38, &
          0.614869006573609D-30, 0.582238667048942D-27, 0.390628369238462D-22, &
          0.821445758255119D-20, 0.402137961842776D-14, 0.651718171878301D-12, &
          -0.211773355803058D-7, 0.264953354380072D-2, -0.135031446451331D-31, &
          -0.607246643970893D-23, -0.402352115234494D-18, -0.744938506925544D-16, &
          0.189917206526237D-12, 0.364975183508473D-5, 0.177274872361946D-25, &
          -0.334952758812999D-18, -0.421537726098389D-8, -0.391048167929649D-1, &
          0.541276911564176D-13, 0.705412100773699D-11, 0.258585887897486D-8, &
          -0.493111362030162D-10, -0.158649699894543D-5, -0.525037427886100D0, &
          0.220019901729615D-2, -0.643064132636925D-2, 0.629154149015048D2, &
          0.135147318617061D3, 0.240560808321713D-6, -0.890763306701305D-3, &
          -0.440209599407714D4, -0.302807107747776D3, 0.159158748314599D4, &
          0.232534272709876D6, -0.792681207132600D6, -0.869871364662769D11, &
          0.354542769185671D12, 0.400849240129329D15]
  PetscInt, parameter :: I_i(39) = [0,3,4,6,7,10,12,14,18,0,3,5,6,8,12,0,3,7,12, &
                                    2,3,4,2,4,7,4,3,5,6,0,0,3,1,0,1,0,1,0,1]
  PetscInt, parameter :: J_i(39) = [-12,-12,-12,-12,-12,-12,-12,-12,-12,-10,-10, &
                                    -10,-10,-10,-10,-8,-8,-8,-8,-6,-6,-6,-5,-5, &
                                    -5,-4,-3,-3,-3,-2,-1,-1,0,1,1,2,4,5,6]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0031d0, p_star = 23.d0, T_star = 650.d0
  PetscReal, parameter :: a = 0.976d0, b = 0.997d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = exp(sum(n_i * (pi - a)**I_i * (theta - b)**J_i)) * v_star
end function IF97_subregion_3n

function IF97_subregion_3o(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(24) = [0.128746023979718D-34, &
          -0.735234770382342D-11, 0.289078692149150D-2, 0.244482731907223D0, &
          0.141733492030985D-23, -0.354533853059476D-28, -0.594539202901431D-17, &
          -0.585188401782779D-8, 0.201377325411803D-5, 0.138647388209306D1, &
          -0.173959365084772D-4, 0.137680878349369D-2, 0.814897605805513D-14, &
          0.425596631351839D-25, -0.387449113787755D-17, 0.139814747930240D-12, &
          -0.171849638951521D-2, 0.641890529513296D-21, 0.118960578072018D-10, &
          -0.155282762571611D-17, 0.233907907347507D-7, -0.174093247766213D-12, &
          0.377682649089149D-8, -0.516720236575302D-10]
  PetscInt, parameter :: I_i(24) = [0,0,0,2,3,4,4,4,4,4,5,5,6,7,8,8,8,10,10,14, &
                                    14,20,20,24]
  PetscInt, parameter :: J_i(24) = [-12,-4,-1,-1,-10,-12,-8,-5,-4,-1,-4,-3,-8,-12, &
                                    -10,-8,-4,-12,-8,-12,-8,-12,-10,-12]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0034d0, p_star = 23.d0, T_star = 650.d0
  PetscReal, parameter :: a = 0.974d0, b = 0.996d0
  PetscReal, parameter :: c = 0.5d0, d = 1.00d0, e = 1.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3o

function IF97_subregion_3p(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(27) = [-0.982825342010366D-4, &
          0.105145700850612D1, 0.116033094095084D3, 0.324664750281543D4, &
          -0.123592348610137D4, -0.561403450013495D-1, 0.856677401640869D-7, &
          0.236313425393924D3, 0.972503292350109D-2, -0.103001994531927D1, &
          -0.149653706199162D-8, -0.215743778861592D-4, -0.834452198291445D1, &
          0.586602660564988D0, 0.343480022104968D-25, 0.816256095947021D-5, &
          0.294985697916798D-2, 0.711730466276584D-16, 0.400954763806941D-9, &
          0.107766027032853D2, -0.409449599138182D-6, -0.729121307758902D-5, &
          0.677107970938909D-8, 0.602745973022975D-7, -0.382323011855257D-10, &
          0.179946628317437D-2, -0.345042834640005D-3]
  PetscInt, parameter :: I_i(27) = [0,0,0,0,1,2,3,3,4,6,7,7,8,10,12,12,12,14,14, &
                                    14,16,18,20,22,24,24,36]
  PetscInt, parameter :: J_i(27) = [-1,0,1,2,1,-1,-3,0,-2,-2,-5,-4,-2,-3,-12,-6, &
                                    -5,-10,-8,-3,-8,-8,-10,-10,-12,-8,-12]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0041d0, p_star = 23.d0, T_star = 650.d0
  PetscReal, parameter :: a = 0.972d0, b = 0.997d0
  PetscReal, parameter :: c = 0.5d0, d = 1.00d0, e = 1.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3p

function IF97_subregion_3q(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(24) = [-0.820433843259950D5, &
          0.473271518461586D11, -0.805950021005413D-1, 0.328600025435980D2, &
          -0.356617029982490D4, -0.172985781433335D10, 0.351769232729192D8, &
          -0.775489259985144D6, 0.710346691966018D-4, 0.993499883820274D5, &
          -0.642094171904570D0, -0.612842816820083D4, 0.232808472983776D3, &
          -0.142808220416837D-4, -0.643596060678456D-2, -0.428577227475614D1, &
          0.225689939161918D4, 0.100355651721510D-2, 0.333491455143516D0, &
          0.109697576888873D1, 0.961917379376452D0, -0.838165632204598D-1, &
          0.247795908411492D1, -0.319114969006533D4]
  PetscInt, parameter :: I_i(24) = [-12,-12,-10,-10,-10,-10,-8,-6,-5,-5,-4,-4,-3, &
                                    -2,-2,-2,-2,-1,-1,-1,0,1,1,1]
  PetscInt, parameter :: J_i(24) = [10,12,6,7,8,10,8,6,2,5,3,4,3,0,1,2,4,0,1,2, &
                                    0,0,1,3]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0022d0, p_star = 23.d0, T_star = 650.d0
  PetscReal, parameter :: a = 0.848d0, b = 0.983d0
  PetscReal, parameter :: c = 1.0d0, d = 1.00d0, e = 4.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3q

function IF97_subregion_3r(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(27) = [0.144165955660863D-2, &
          -0.701438599628258D13, -0.830946716459219D-16, 0.261975135368109D0, &
          0.393097214706245D3, -0.104334030654021D5, 0.490112654154211D9, &
          -0.147104222772069D-3, 0.103602748043408D1, 0.305308890065089D1, &
          -0.399745276971264D7, 0.569233719593750D-11, -0.464923504407778D-1, &
          -0.535400396512906D-17, 0.399988795693162D-12, -0.536479560201811D-6, &
          0.159536722411202D-1, 0.270303248860217D-14, 0.244247453858506D-7, &
          -0.983430636716454D-5, 0.663513144224454D-1, -0.993456957845006D1, &
          0.546491323528491D3, -0.143365406393758D5, 0.150764974125511D6, &
          -0.337209709340105D-9, 0.377501980025469D-8]
  PetscInt, parameter :: I_i(27) = [-8,-8,-3,-3,-3,-3,-3,0,0,0,0,3,3,8,8,8,8,10, &
                                    10,10,10,10,10,10,10,12,14]
  PetscInt, parameter :: J_i(27) = [6,14,-3,3,4,5,8,-1,0,1,5,-6,-2,-12,-10,-8,-5, &
                                    -12,-10,-8,-6,-5,-4,-3,-2,-12,-12]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0054d0, p_star = 23.d0, T_star = 650.d0
  PetscReal, parameter :: a = 0.874d0, b = 0.982d0
  PetscReal, parameter :: c = 1.0d0, d = 1.00d0, e = 1.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3r

function IF97_subregion_3s(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(29) = [-0.532466612140254D23, &
          0.100415480000824D32, -0.191540001821367D30, 0.105618377808847D17, &
          0.202281884477061D59, 0.884585472596134D8, 0.166540181638363D23, &
          -0.313563197669111D6, -0.185662327545324D54, -0.624942093918942D-1, &
          -0.504160724132590D10, 0.187514491833092D5, 0.121399979993217D-2, &
          0.188317043049455D1, -0.167073503962060D4, 0.965961650599775D0, &
          0.294885696802488D1, -0.653915627346115D5, 0.604012200163444D50, &
          -0.198339358557937D0, -0.175984090163501D58, 0.356314881403987D1, &
          -0.575991255144384D3, 0.456213415338071D5, -0.109174044987829D8, &
          0.437796099975134D34, -0.616552611135792D46, 0.193568768917797D10, &
          0.950898170425042D54]
  PetscInt, parameter :: I_i(29) = [-12,-12,-10,-8,-6,-5,-5,-4,-4,-3,-3,-2,-1,-1, &
                                    -1,0,0,0,0,1,1,3,3,3,4,4,4,5,14]
  PetscInt, parameter :: J_i(29) = [20,24,22,14,36,8,16,6,32,3,8,4,1,2,3,0,1,4, &
                                    28,0,32,0,1,2,3,18,24,4,24]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0022d0, p_star = 21.d0, T_star = 640.d0
  PetscReal, parameter :: a = 0.886d0, b = 0.990d0
  PetscReal, parameter :: c = 1.0d0, d = 1.00d0, e = 4.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3s

function IF97_subregion_3t(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(33) = [0.155287249586268D1, &
          0.664235115009031D1, -0.289366236727210D4, -0.385923202309848D13, &
          -0.291002915783761D1, -0.829088246858083D12, 0.176814899675218D1, &
          -0.534686695713469D9, 0.160464608687834D18, 0.196435366560186D6, &
          0.156637427541729D13, -0.178154560260006D1, -0.229746237623692D16, &
          0.385659001648006D8, 0.110554446790543D10, -0.677073830687349D14, &
          -0.327910592086523D31, -0.341552040860644D51, -0.527251339709047D21, &
          0.245375640937055D24, -0.168776617209269D27, 0.358958955867578D29, &
          -0.656475280339411D36, 0.355286045512301D39, 0.569021454413270D58, &
          -0.700584546433113D48, -0.705772623326374D65, 0.166861176200148D53, &
          -0.300475129680486D61, -0.668481295196808D51, 0.428432338620678D69, &
          -0.444227367758304D72, -0.281396013562745D77]
  PetscInt, parameter :: I_i(33) = [0,0,0,0,1,1,2,2,2,3,3,4,4,7,7,7,7,7,10,10,10, &
                                    10,10,18,20,22,22,24,28,32,32,32,36]
  PetscInt, parameter :: J_i(33) = [0,1,4,12,0,10,0,6,14,3,8,0,10,3,4,7,20,36,10, &
                                    12,14,16,22,18,32,22,36,24,28,22,32,36,36]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0088d0, p_star = 20.d0, T_star = 650.d0
  PetscReal, parameter :: a = 0.803d0, b = 1.020d0
  PetscReal, parameter :: c = 1.0d0, d = 1.00d0, e = 1.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3t

function IF97_subregion_3u(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(38) = [0.122088349258355D18, &
          0.104216468608488D10, -0.882666931564652D16, 0.259929510849499D20, &
          0.222612779142211D15, -0.878473585050085D18, -0.314432577551552D22, &
          -0.216934916996285D13, 0.159079648196849D21, -0.339567617303423D3, &
          0.884387651337836D13, -0.843405926846418D21, 0.114178193518022D2, &
          -0.122708229235641D-3, -0.106201671767107D3, 0.903443213959313D25, &
          -0.693996270370852D28, 0.648916718965575D-8, 0.718957567127851D4, &
          0.105581745346187D-2, -0.651903203602581D15, -0.160116813274676D25, &
          -0.510254294237837D-8, -0.152355388953402D0, 0.677143292290144D12, &
          0.276378438378930D15, 0.116862983141686D-1, -0.301426947980171D14, &
          0.169719813884840D-7, 0.104674840020929D27, -0.108016904560140D5, &
          -0.990623601934295D-12, 0.536116483602738D7, 0.226145963747881D22, &
          -0.488731565776210D-9, 0.151001548880670D-4, -0.227700464643920D5, &
          -0.781754507698846D28]
  PetscInt, parameter :: I_i(38) = [-12,-10,-10,-10,-8,-8,-8,-6,-6,-5,-5,-5,-3, &
                                    -1,-1,-1,-1,0,0,1,2,2,3,5,5,5,6,6,8,8,10,12, &
                                    12,12,14,14,14,14]
  PetscInt, parameter :: J_i(38) = [14,10,12,14,10,12,14,8,12,4,8,12,2,-1,1,12, &
                                    14,-3,1,-2,5,10,-5,-4,2,3,-5,2,-8,8,-4,-12, &
                                    -4,4,-12,-10,-6,6]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0026d0, p_star = 23.d0, T_star = 650.d0
  PetscReal, parameter :: a = 0.902d0, b = 0.988d0
  PetscReal, parameter :: c = 1.0d0, d = 1.00d0, e = 1.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3u

function IF97_subregion_3v(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(39) = [-0.415652812061591D-54, &
          0.177441742924043D-60, -0.357078668203377D-54, 0.359252213604114D-25, &
          -0.259123736380269D2, 0.594619766193460D5, -0.624184007103158D11, &
          0.313080299915944D17, 0.105006446192036D-8, -0.192824336984852D-5, &
          0.654144373749937D6, 0.513117462865044D13, -0.697595750347391D19, &
          -0.103977184454767D29, 0.119563135540666D-47, -0.436677034051655D-41, &
          0.926990036530639D-29, 0.587793105620748D21, 0.280375725094731D-17, &
          -0.192359972440634D23, 0.742705723302738D27, -0.517429682450605D2, &
          0.820612048645469D7, -0.188214882341448D-8, 0.184587261114837D-1, &
          -0.135830407782663D-5, -0.723681885626348D17, -0.223449194054124D27, &
          -0.111526741826431D-34, 0.276032601145151D-28, 0.134856491567853D15, &
          0.652440293345860D-9, 0.510655119774360D17, -0.468138358908732D32, &
          -0.760667491183279D16, -0.417247986986821D-18, 0.312545677756104D14, &
          -0.100375333864186D15, 0.247761392329058D27]
  PetscInt, parameter :: I_i(39) = [-10,-8,-6,-6,-6,-6,-6,-6,-5,-5,-5,-5,-5,-5, &
                                    -4,-4,-4,-4,-3,-3,-3,-2,-2,-1,-1,0,0,0,1,1, &
                                    3,4,4,4,5,8,10,12,14]
  PetscInt, parameter :: J_i(39) = [-8,-12,-12,-3,5,6,8,10,1,2,6,8,10,14,-12,-10, &
                                    -6,10,-3,10,12,2,4,-2,0,-2,6,10,-12,-10,3,-6, &
                                    3,10,2,-12,-2,-3,1]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0031d0, p_star = 23.d0, T_star = 650.d0
  PetscReal, parameter :: a = 0.960d0, b = 0.995d0
  PetscReal, parameter :: c = 1.0d0, d = 1.00d0, e = 1.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3v

function IF97_subregion_3w(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(35) = [-0.586219133817016D-7, &
          -0.894460355005526D11, 0.531168037519774D-30, 0.109892402329239D0, &
          -0.575368389425212D-1, 0.228276853990249D5, -0.158548609655002D19, &
          0.329865748576503D-27, -0.634987981190669D-24, 0.615762068640611D-8, &
          -0.961109240985747D8, -0.406274286652625D-44, -0.471103725498077D-12, &
          0.725937724828145D0, 0.187768525763682D-38, -0.103308436323771D4, &
          -0.662552816342168D-1, 0.579514041765710D3, 0.237416732616644D-26, &
          0.271700235739893D-14, -0.907886213483600D2, -0.171242509570207D-36, &
          0.156792067854621D3, 0.923261357901470D0, -0.597865988422577D1, &
          0.321988767636389D7, -0.399441390042203D-29, 0.493429086046981D-7, &
          0.812036983370565D-19, -0.207610284654137D-11, -0.340821291419719D-6, &
          0.542000573372233D-17, -0.856711586510214D-12, 0.266170454405981D-13, &
          0.858133791857099D-5]
  PetscInt, parameter :: I_i(35) = [-12,-12,-10,-10,-8,-8,-8,-6,-6,-6,-6,-5,-4, &
                                    -4,-3,-3,-2,-2,-1,-1,-1,0,0,1,2,2,3,3,5,5,5, &
                                    8,8,10,10]
  PetscInt, parameter :: J_i(35) = [8,14,-1,8,6,8,14,-4,-3,2,8,-10,-1,3,-10,3,1, &
                                    2,-8,-4,1,-12,1,-1,-1,2,-12,-5,-10,-8,-6,-12, &
                                    -10,-12,-8]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0039d0, p_star = 23.d0, T_star = 650.d0
  PetscReal, parameter :: a = 0.959d0, b = 0.995d0
  PetscReal, parameter :: c = 1.0d0, d = 1.00d0, e = 4.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3w

function IF97_subregion_3x(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(36) = [0.377373741298151D19, &
          -0.507100883722913D13, -0.103363225598860D16, 0.184790814320773D-5, &
          -0.924729378390945D-3, -0.425999562292738D24, -0.462307771873973D-12, &
          0.107319065855767D22, 0.648662492280682D11, 0.244200600688281D1, &
          -0.851535733484258D10, 0.169894481433592D22, 0.215780222509020D-26, &
          -0.320850551367334D0, -0.382642448458610D17, -0.275386077674421D-28, &
          -0.563199253391666D6, -0.326068646279314D21, 0.397949001553184D14, &
          0.100824008584757D-6, 0.162234569738433D5, -0.432355225319745D11, &
          -0.592874245598610D12, 0.133061647281106D1, 0.157338197797544D7, &
          0.258189614270853D14, 0.262413209706358D25, -0.920011937431142D-1, &
          0.220213765905426D-2, -0.110433759109547D2, 0.847004870612087D7, &
          -0.592910695762536D9, -0.183027173269660D-4, 0.181339603516302D0, &
          -0.119228759669889D4, 0.430867658061468D7]
  PetscInt, parameter :: I_i(36) = [-8,-6,-5,-4,-4,-4,-3,-3,-1,0,0,0,1,1,2,3,3, &
                                    3,4,5,5,5,6,8,8,8,8,10,12,12,12,12,14,14,14, &
                                    14]
  PetscInt, parameter :: J_i(36) = [14,10,10,1,2,14,-2,12,5,0,4,10,-10,-1,6,-12, &
                                    0,8,3,-6,-2,1,1,-6,-3,1,8,-8,-10,-8,-5,-4,-12, &
                                    -10,-8,-6]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0049d0, p_star = 23.d0, T_star = 650.d0
  PetscReal, parameter :: a = 0.910d0, b = 0.988d0
  PetscReal, parameter :: c = 1.0d0, d = 1.00d0, e = 1.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3x

function IF97_subregion_3y(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(20) = [-0.525597995024633D-9, &
          0.583441305228407D4, -0.134778968457925D17, 0.118973500934212D26, &
          -0.159096490904708D27, -0.315839902302021D-6, 0.496212197158239D3, &
          0.327777227273171D19, -0.527114657850696D22, 0.210017506281863D-16, &
          0.705106224399834D21, -0.266713136106469D31, -0.145370512554562D-7, &
          0.149333917053130D28, -0.149795620287641D8, -0.381881906271100D16, &
          0.724660165585797D-4, -0.937808169550193D14, 0.514411468376383D10, &
          -0.828198594040141D5]
  PetscInt, parameter :: I_i(20) = [0,0,0,0,1,2,2,2,2,3,3,3,4,4,5,5,8,8,10,12]
  PetscInt, parameter :: J_i(20) = [-3,1,5,8,8,-4,-1,4,5,-8,4,8,-6,6,-2,1,-8,-2, &
                                    -5,-8]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0031d0, p_star = 22.d0, T_star = 650.d0
  PetscReal, parameter :: a = 0.996d0, b = 0.994d0
  PetscReal, parameter :: c = 1.0d0, d = 1.00d0, e = 4.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3y

function IF97_subregion_3z(p,T) result (v)
  implicit none
  PetscReal, parameter :: n_i(23) = [0.244007892290650D-10, &
          -0.463057430331242D7, 0.728803274777712D10, 0.327776302858856D16, &
          -0.110598170118409D10, -0.323899915729957D13, 0.923814007023245D16, &
          0.842250080413712D-12, 0.663221436245506D12, -0.167170186672139D15, &
          0.253749358701391D4, -0.819731559610523D-20, 0.328380587890663D12, &
          -0.625004791171543D8, 0.803197957462023D21, -0.204397011338353D-10, &
          -0.378391047055938D4, 0.972876545938620D-2, 0.154355721681459D2, &
          -0.373962862928643D4, -0.682859011374572D11, -0.248488015614543D-3, &
          0.394536049497068D7]
  PetscInt, parameter :: I_i(23) = [-8,-6,-5,-5,-4,-4,-4,-3,-3,-3,-2,-1,0,1,2,3, &
                                    3,6,6,6,6,8,8]
  PetscInt, parameter :: J_i(23) = [3,6,6,8,5,6,8,-2,5,6,2,-6,3,1,6,-6,-2,-6,-5, &
                                    -4,-1,-8,-4]
  PetscReal :: p, T, v, pi, theta
  PetscReal, parameter :: v_star = 0.0038d0, p_star = 22.d0, T_star = 650.d0
  PetscReal, parameter :: a = 0.993d0, b = 0.994d0
  PetscReal, parameter :: c = 1.0d0, d = 1.00d0, e = 4.d0
  pi = p/(p_star * 1.0D+6)
  theta = (T + 273.15d0)/T_star
  v = sum(n_i * ((pi - a)**c)**I_i * ((theta - b)**d)**J_i)**e * v_star
end function IF97_subregion_3z

function IF97SaturationTemperature(Ps) result(Ts)
  ! section 8.2 of IAPWS R7-97(2012)
  ! this is the inverse of EOSWaterSaturationPressureIF97()
  implicit none
  PetscReal, parameter :: n(10) = [0.11670521452767D4, -0.72421316703206D6, &
       -0.17073846940092D2, 0.12020824702470D5, -0.32325550322333D7, &
       0.14915108613530D2, -0.48232657361591D4, 0.40511340542057D6, &
       -0.23855557567849D0, 0.65017534844798D3]
  PetscReal :: Ps ! saturation pressure (Pa)
  PetscReal :: Ts ! saturation temperature (K)
  PetscReal :: beta, betasq, D, E, F, G
  betasq = sqrt(Ps/1.0D+6) ! eqn 29a
  beta = sqrt(betasq)
  E = betasq + n(3)*beta + n(6)
  F = n(1)*betasq + n(4)*beta + n(7)
  G = n(2)*betasq + n(5)*beta + n(8)
  D = 2.d0 * G / (-F-sqrt(F**2 - 4.d0*E*G))
  Ts = (n(10) + D - sqrt((n(10) + D)**2 - 4.d0*(n(9) + n(10)*D)))/2.d0 ! eqn 31

end function IF97SaturationTemperature

! ************************************************************************** !

subroutine EOSWaterEnthalpyIFC67(t,p,calculate_derivatives,hw, &
                                 hwp,hwt,ierr)

!  This subroutine calculates water and steam-gas mixture properties.
!  The water and steam properties are valid in the range of:
!
!            0 < p < 165.4 * 10^5 pascals (165.4 bars)
!            0 < t < 350 centigrade (623.15 Kelvin)
!
!  The properties cover densities, enthalpies, internal energies,
!  and partial derivatives of these quanties with respect to
!  pressure and temperature.
!
!  For saturated fluid, it will also calculate water saturation
!  temperature for the specified pressure using Newton-Raphson and
!  the derivative dts/dp (=tsp) or Ps for a given temperature.
!
!  Ref.: International Formulation Committee of the Sixth International
!       Conference on Properties of Steam (1967).

  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: hw,hwp,hwt
  PetscErrorCode, intent(inout) :: ierr

  PetscInt :: i

  PetscReal, save :: aa(0:22)
  PetscReal, save :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12

  PetscReal :: beta,beta2x,beta4,theta,utheta,theta2x,theta18,theta20
  PetscReal :: xx,yy,zz
  PetscReal :: u0,u1!,u2,u3,u4,u5,u6,u7,u8,u9
  PetscReal :: tempreal
  PetscReal :: v0_1, v1_1, v2_1, v3_1, v4_1
  PetscReal :: v1_2, v2_2, v3_2, v4_2, v20_2, v40_2
  PetscReal :: v1_3, v2_3, v3_3, v4_3
  PetscReal :: v1_4, v2_4, v3_4
  PetscReal :: v1_5, v2_5
  PetscReal :: v1_6
  PetscReal :: term1,term2,term2t,term3,term3t,term3p,term4,term4t,term4p, &
               term5,term5t,term5p,term6,term6t,term6p,term7,term7t,term7p
  PetscReal :: dv2t,dv2p,dv3t
!  PetscReal :: vr,ypt,yptt,zpt,zpp,vrpt,vrpp,cnv
  PetscReal :: ypt,yptt,zpt,zpp
  PetscReal :: tc1,pc1,vc1,utc1,upc1,vc1mol
  PetscReal, parameter :: zero = 0.d0
  PetscReal, parameter :: one = 1.d0
  PetscReal, parameter :: two = 2.d0
  PetscReal, parameter :: three = 3.d0
  PetscReal, parameter :: five = 5.d0
  PetscReal, parameter :: six = 6.d0
  PetscReal, parameter :: nine = 9.d0
  PetscReal, parameter :: ten = 10.d0

  data aa/ &
!-----data aa0,aa1,aa2,aa3/
        6.824687741d03,-5.422063673d02,-2.096666205d04, 3.941286787d04, &
!-----data aa4,aa5,aa6,aa7/
        -6.733277739d04, 9.902381028d04,-1.093911774d05, 8.590841667d04, &
!-----data aa8,aa9,aa10,aa11/
        -4.511168742d04, 1.418138926d04,-2.017271113d03, 7.982692717d00, &
!-----data aa12,aa13,aa14,aa15/
        -2.616571843d-2, 1.522411790d-3, 2.284279054d-2, 2.421647003d02, &
!-----data aa16,aa17,aa18,aa19/
        1.269716088d-10,2.074838328d-7, 2.174020350d-8, 1.105710498d-9, &
!-----data aa20,aa21,aa22/
        1.293441934d01, 1.308119072d-5, 6.047626338d-14/

  data a1,a2,a3,a4/ &
  8.438375405d-1, 5.362162162d-4, 1.720000000d00, 7.342278489d-2/
  data a5,a6,a7,a8/ &
  4.975858870d-2, 6.537154300d-1, 1.150000000d-6, 1.510800000d-5/
  data a9,a10,a11,a12/ &
  1.418800000d-1, 7.002753165d00, 2.995284926d-4, 2.040000000d-1/

  tc1 = H2O_CRITICAL_TEMPERATURE    ! K
  pc1 = H2O_CRITICAL_PRESSURE     ! Pa
  vc1 = 0.00317d0  ! m^3/kg
  utc1 = one/tc1   ! 1/C
  upc1 = one/pc1   ! 1/Pa
  vc1mol = vc1*FMWH2O ! m^3/kmol

  theta = (t+273.15d0)*utc1
  theta2x = theta*theta
  theta18 = theta**18.d0
  theta20 = theta18*theta2x

  beta = p*upc1
  beta2x = beta*beta
  beta4  = beta2x*beta2x

  yy = one-a1*theta2x-a2*theta**(-6.d0)
  xx = a3*yy*yy-two*(a4*theta-a5*beta)
!   Note: xx may become negative near the critical point-pcl.
  if (xx.gt.zero) then
    xx = sqrt(xx)
  else
    write(*,*) 'Warning: negative term in density (eos_water.F90:&
      &EOSWaterEnthalpyIFC67):'
    write(*,*) 't= ',t,' p= ',p,' xx= ',xx
    ierr = 1
    xx = 1.d-6               !set arbitrarily
  end if
  zz = yy + xx
  u0 = -five/17.d0
  u1 = aa(11)*a5*zz**u0
#if 0
  u2 = one/(a8+theta**11)
  u3 = aa(17)+(two*aa(18)+three*aa(19)*beta)*beta
  u4 = one/(a7+theta18*theta)
  u5 = (a10+beta)**(-4)
  u6 = a11-three*u5
  u7 = aa(20)*theta18*(a9+theta2x)
  u8 = aa(15)*(a6-theta)**9

  vr = u1+aa(12)+theta*(aa(13)+aa(14)*theta)+u8*(a6-theta) &
        +aa(16)*u4-u2*u3-u6*u7+(three*aa(21)*(a12-theta) &
        +four*aa(22)*beta/theta20)*beta2x

  dwmol = one/(vr*vc1mol) ! kmol/m^3
  dw = one/(vr*vc1) ! kg/m^3

  !---calculate derivatives for water density
  if (calculate_derivatives) then
    u9 = u0*u1/zz
    vrpt = u9*zpt+aa(13)+two*aa(14)*theta-ten*u8 &
        -19.d0*aa(16)*u4*u4*theta18+11.d0*u2*u2*u3*theta**10.d0 &
        -aa(20)*u6*(18.d0*a9*theta18+20.d0*theta20)/theta &
        -(three*aa(21)+80.d0*aa(22)*beta/(theta20*theta))*beta2x

    vrpp = u9*zpp-u2*(two*aa(18)+six*aa(19)*beta)-12.d0*u7*u5/ &
        (a10+beta)+(six*aa(21)*(a12-theta)+12.d0*aa(22)*beta/ &
        theta20)*beta

    cnv = -one/(vc1mol*vr*vr)
    dwt = cnv*vrpt*utc1 ! kmol/m^3/C
    dwp = cnv*vrpp*upc1 ! kmol/m^3/Pa
  else
    dwt = UNINITIALIZED_DOUBLE
    dwp = UNINITIALIZED_DOUBLE
  endif
#endif

  ! ypt used for enthalpy even if derivative not calculated
  ypt = six*a2*theta**(-7.d0)-two*a1*theta


!---compute enthalpy internal energy and derivatives for water
  utheta = one/theta
  term1 = aa(0)*theta
  term2 = -aa(1)
  ! term2t is part of the derivative calc., but left here to avoid
  ! recomputing the expensive do loop
  term2t = zero
  do i = 3,10
    tempreal = dfloat(i-2)*aa(i)*theta**(i-1)
    term2t = term2t+tempreal*utheta*dfloat(i-1)
    term2 = term2+tempreal
  end do

  ! "v" section 1
  v0_1 = u1/a5
  v2_1 = 17.d0*(zz/29.d0-yy/12.d0)+five*theta*ypt/12.d0
  v3_1 = a4*theta-(a3-one)*theta*yy*ypt
  v1_1 = zz*v2_1+v3_1
  term3 = v0_1*v1_1

  ! block 1 removed from here

  ! "v" section 2
  v1_2 = nine*theta+a6
  v20_2 = (a6-theta)
  v2_2 = v20_2**9.d0
  v3_2 = a7+20.d0*theta**19.d0
  v40_2 = a7+theta**19.d0
  v4_2 = one/(v40_2*v40_2)
  ! term4p is a derivative, but left due to dependency in term4
  term4p = aa(12)-aa(14)*theta2x+aa(15)*v1_2*v2_2+aa(16)*v3_2*v4_2
  term4 = term4p*beta

  ! block 2 removed from here

  ! "v" section 3
  v1_3 = beta*(aa(17)+aa(18)*beta+aa(19)*beta2x)
  v2_3 = 12.d0*theta**11.d0+a8
  v4_3 = one/(a8+theta**11.d0)
  v3_3 = v4_3*v4_3
  term5 = v1_3*v2_3*v3_3

  ! block 3 removed from here

  ! "v" section 4
  v1_4 = (a10+beta)**(-3.d0)+a11*beta
  v3_4 = (17.d0*a9+19.d0*theta2x)
  v2_4 = aa(20)*theta18*v3_4
  term6 = v1_4*v2_4

  ! block 4 removed from here

  ! "v" section 5
  v1_5 = 21.d0*aa(22)/theta20*beta4
  v2_5 = aa(21)*a12*beta2x*beta
  term7 = v1_5+v2_5

  ! "v" section 6
  v1_6 = pc1*vc1mol
  hw = (term1-term2+term3+term4-term5+term6+term7)*v1_6

  if (calculate_derivatives) then

    zpt = ypt+(a3*yy*ypt-a4)/xx
    zpp = a5/xx

    ! block 1
    yptt = -two*a1-42.d0*a2/theta**8.d0
    dv2t = 17.d0*(zpt/29.d0-ypt/12.d0)+five/12.d0*(ypt+theta*yptt)
    dv3t = a4-(a3-one)*(theta*yy*yptt+yy*ypt+theta*ypt*ypt)
    dv2p = 17.d0*zpp/29.d0
    v4_1 = five*v1_1/(17.d0*zz)
    term3t = v0_1*(zz*dv2t+(v2_1-v4_1)*zpt+dv3t)
    term3p = v0_1*(zz*dv2p+(v2_1-v4_1)*zpp)

    ! block 2
    term4t = (-two*aa(14)*theta+nine*aa(15)*(v2_2-v1_2*v2_2/v20_2) &
             +38.d0*theta18*aa(16)*(ten*v4_2-v3_2*v4_2/v40_2))*beta

    ! block 3
    term5p = v3_3*v2_3*(aa(17)+two*aa(18)*beta+three*aa(19)*beta2x)
    term5t = v1_3*(132.d0*v3_3*theta**10.d0-22.d0*v2_3*v3_3*v4_3*theta**10.d0)

    ! block 4
    term6p = v2_4*(a11-three*(a10+beta)**(-4.d0))
    term6t = v1_4*aa(20)*theta18*(18.d0*v3_4*utheta+38.d0*theta)

    ! block 5
    term7p = beta2x*(three*aa(21)*a12+84.d0*aa(22)*beta/theta20)
    term7t = -420.d0*aa(22)*beta4/(theta20*theta)

    hwp = (term3p+term4p-term5p+term6p+term7p)*vc1mol
    hwt = (aa(0)-term2t+term3t+term4t-term5t+term6t+term7t)*v1_6*utc1
  else
    hwp = UNINITIALIZED_DOUBLE
    hwt = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterEnthalpyIFC67

subroutine EOSWaterEnthalpyIF97(T,P,calculate_derivatives,hw, &
                                 hwp,hwt,ierr)
  implicit none

  PetscReal, intent(in) :: T   ! Temperature in centigrade
  PetscReal, intent(in) :: P   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: hw,hwp,hwt
  PetscErrorCode, intent(inout) :: ierr

  PetscReal, parameter :: Tf = 273.15d0 !K
  PetscReal, parameter :: p_ref = 16.53d6 !Pa
  PetscReal, parameter :: T_ref = 1386.d0  ! K
  PetscReal, parameter :: R = 0.461526d0 ! kJ/kg-K
  PetscReal, parameter :: n_i(34) = [1.4632971213167d-1, &
          -8.4548187169114d-1, -0.37563603672040d1, 0.33855169168385d1, &
          -9.5791963387872d-1, 1.5772038513228d-1, -0.16616417199501d-1, &
          0.81214629983568d-3, 0.28319080123804d-3, -0.60706301565874d-3, &
          -0.18990068218419d-1, -0.32529748770505d-1, -0.21841717175414d-1, &
          -0.52838357969930d-4, -0.47184321073267d-3, -0.30001780793026d-3, &
          0.47661393906987d-4, -0.44141845330846d-5, -0.72694996297594d-15, &
          -0.31679644845054d-4, -0.28270797985312d-5, -0.85205128120103d-9, &
          -0.22425281908000d-5, -0.65171222895601d-6, -0.14341729937924d-12, &
          -0.40516996860117d-6, -0.12734301741641d-8, -0.17424871230634d-9, &
          -0.68762131295531d-18, 0.14478307828521d-19, 0.26335781662795d-22, &
          -0.11947622640071d-22, 0.18228094581404d-23, -0.93537087292458d-25]
  PetscInt, parameter :: I_i(34) = [0,0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3, &
                                    3,3,4,4,4,5,8,8,21,23,29,30,31,32]
  PetscInt, parameter :: J_i(34) = [-2,-1,0,1,2,3,4,5,-9,-7,-1,0,1,3,-3,0,1, &
                                     3,17,-4,0,6,-5,-2,10,-8,-11,-6,-29,-31, &
                                     -38,-39,-40,-41]
  PetscReal :: pi, tao,  g_tao, T_temp

  T_temp = T+Tf
  pi = P/p_ref
  tao = T_ref/T_temp

  if (Tf <= 623.15d0) then
    ! Region 1: Valid from 273.15K to 623.15 K, Ps(T) to 100MPa
    g_tao = sum((n_i*(7.1d0-pi)**(I_i))*J_i*(tao-1.222d0)**(J_i-1))

    hw = g_tao *T_ref*R
    hw = hw*FMWH2O * 1.d3

    if (calculate_derivatives) then
      hwp = T_ref*R/p_ref * sum(-n_i*I_i*(7.1d0-pi)**(I_i-1) * &
           J_i*(tao-1.222d0)**(J_i-1))
      hwt = -T_ref*T_ref*R/(T_temp*T_temp) * sum(n_i*(7.1d0-pi)**(I_i)* &
           J_i*(J_i-1)*(tao-1.222d0)**(J_i-2))
      hwp = hwp*FMWH2O * 1.d3
      hwt = hwt*FMWH2O * 1.d3
    else
      hwp = UNINITIALIZED_DOUBLE
      hwt = UNINITIALIZED_DOUBLE
    endif
  else
    ! Region 3: Valid in "wedge" >623.15K, >Ps(T), and 100MPa
    ! moved to subroutine to avoid weighting down this more-used routine
    call EOSWaterEnthalpyIF97Region3(T,P,.false.,hw,hwp,hwt)
  end if

end subroutine EOSWaterEnthalpyIF97

subroutine EOSWaterEnthalpyIF97Region3(T,P,calculate_derivatives,hw,hwp,hwt)

  ! Region 3 is in terms of Helmholtz free energy, rather than
  ! Gibbs free energy. Enthalpy is in terms of density and temperature.

  implicit none

  PetscReal, intent(in) :: T   ! Temperature in centigrade
  PetscReal, intent(in) :: P   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: hw,hwp,hwt

  PetscReal, parameter :: Tf = 273.15d0 !K
  PetscReal, parameter :: T_c = 647.096d0 ! K
  PetscReal, parameter :: rho_c = 322.d0 ! kg/m^3
  PetscReal, parameter :: R = 0.461526d0 ! kJ/kg-K
  PetscReal, parameter :: n_i(40) = [0.10658070028513D1, &
       -0.15732845290239D2, 0.20944395974307D2, -0.76867707878716D1, &
       0.26185947787954D1, -0.28080781148620D1, 0.12053369696517D1, &
       -0.84566812812502D-2, -0.12654315477714D1, -0.11524407806681D1, &
       0.88521043984318d0, -0.64207765181607d0, 0.38493460186671d0, &
       -0.85214708824206d0, 0.48972281541877D1, -0.30502617256965d1, &
       0.39420536879154d-1, 0.12558408424308d0, -0.27999329698710d0, &
       0.13899799569460D1, -0.20189915023570D1, -0.82147637173963D-2, &
       -0.47596035734923d0, 0.43984074473500D-1, -0.44476435428739d0, &
       0.90572070719733d0, 0.70522450087967d0, 0.10770512626332d0, &
       -0.32913623258954d0, -0.50871062041158d0, -0.22175400873096D-1, &
       0.94260751665092D-1, 0.16436278447961d0, -0.13503372241348D-1, &
       -0.14834345352472D-1, 0.57922953628084D-3, 0.32308904703711D-2, &
       0.80964802996215D-4, -0.16557679795037D-3, -0.44923899061815D-4]
  PetscInt, parameter :: I_i(40) = [0,0,0,0,0, 0,0,0,1,1, 1,1,2,2,2, &
       2,2,2,3,3, 3,3,3,4,4, 4,4,5,5,5, 6,6,6,7,8, 9,9,10,10,11]
  PetscInt, parameter :: J_i(40) = [0,0,1,2,7, 10,12,23,2,6, 15,17,0,2,6, &
       7,22,26,0,2, 4,16,26,0,2, 4,26,1,3,26, 0,2,26,2,26, 2,26,0,1,26]
  PetscReal :: tau,  phi_tau, phi_delta, T_temp, delta, dw, dumb

  T_temp = T+Tf
  ! P only used to compute density
  call EOSWaterDensityIF97Region3(T,P,.false.,dw,dumb,dumb,dumb)
  delta = dw/rho_c
  tau = T_c/T_temp

  ! Table 32
  phi_tau = sum(n_i(2:) * delta ** I_i(2:) * &
                J_i(2:) * tau ** (J_i(2:) - 1))
  phi_delta = n_i(1)/delta + sum(n_i(2:) * I_i(2:) * &
                                 delta ** (I_i(2:) - 1) * tau ** J_i(2:))

  ! Table 31
  hw = (tau * phi_tau + delta * phi_delta) * R * T_temp

  if (calculate_derivatives) then
    stop 'IF97 region 3 water enthalpy derivative not implemented yet'
  else
    hwt = UNINITIALIZED_DOUBLE
    hwp = UNINITIALIZED_DOUBLE
  end if

end subroutine EOSWaterEnthalpyIF97Region3

! ************************************************************************** !
subroutine EOSWaterEnthalpyDriesnerExt(T,P,aux,calculate_derivatives,hw,hwp,&
     hwt,ierr)
  !
  ! Calculates brine enthalpy based on the formuatlion in Driesner, 2007.
  !
  ! Author: David Fukuyama
  ! Date: 05/31/21
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: P
  PetscReal, intent(in) :: aux(*)
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: hw,hwp,hwt
  PetscErrorCode, intent(inout) :: ierr

  PetscReal, parameter :: Pa_to_bar = 1.d-5

  PetscReal :: T_h, p_bar, t_c
  PetscReal :: q1,q2, q1x1, q2x1
  PetscReal :: q10,q11,q12,q20,q21,q22,q23
  PetscReal :: dq1_dp,dq2_dp,dq1x1_dp,dq2x1_dp
  PetscReal :: dq10_dp,dq11_dp,dq12_dp,dq20_dp,dq21_dp,dq22_dp,dq23_dp,dT_h_dp
  PetscReal :: s,x, brine_molar_mass

  t_c = T
  p_bar = P*Pa_to_bar

  s = aux(1) !mass frac
  brine_molar_mass = 1.d0 / (s / FMWNACL + (1.d0 - s) / FMWH2O)
  x = s * brine_molar_mass / FMWNACL

  q11 = -32.1724d0 + 0.0621255d0 * p_bar  ! table 5
  q21 = -1.69513d0 - 4.52781d-4 * p_bar - 6.04279d-8  * p_bar**2
  q22 = 0.0612567d0 + 1.88082d-5 * p_bar

  q1x1 = 47.9048d0 - 9.36994d-3 * p_bar + 6.51059d-6 * p_bar**2 ! eq 25
  q2x1 = 0.241022d0 + 3.45087d-5 * p_bar - 4.28356d-9 * p_bar**2 ! eq 26

  q12 = -q11 - q1x1
  q10 = q1x1

  q20 = 1.d0 - q21 * sqrt(q22)
  q23 = q2x1 - q20 - q21 * sqrt(1.d0 + q22)

  q1 = q10 + q11 * (1.d0 - x) + q12 * (1.d0 - x)**2  ! eq 23
  q2 = q20 + q21 * sqrt(x + q22) + (q23 * x)       ! eq 24
  T_h = q1 + q2 * t_c                              ! eq 22

  call EOSWaterEnthalpyIF97(T_h,P,calculate_derivatives,hw,hwp,hwt,ierr)

  if (calculate_derivatives) then
    dq11_dp = 0.0621255d0
    dq21_dp = -4.52781d-4 - 2.d0 * 6.04279d-8  * p_bar
    dq22_dp = 1.88082d-5

    dq1x1_dp = -9.36994d-3 + 2.d0*6.51059d-6 * p_bar
    dq2x1_dp = 3.45087d-5 -2.d0 * 4.28356d-9 * p_bar

    dq12_dp = -dq11_dp - dq1x1_dp
    dq10_dp = dq1x1_dp

    dq20_dp = -q21 * dq22_dp/(2.d0*sqrt(q22)) - dq21_dp
    dq23_dp = dq2x1_dp - dq20_dp - dq21_dp - q21*dq22_dp/(2.d0*sqrt(q22)) - dq21_dp*sqrt(q22)

    dq1_dp = dq10_dp + dq11_dp * (1.d0 - x) + dq12_dp * (1.d0 - x)**2
    dq2_dp = dq20_dp + dq21_dp * sqrt(x + q22) + q21 * (dq22_dp)/(2.d0*sqrt(x+q22)) + dq23_dp * x
    dT_h_dp = dq1_dp + dq2_dp * t_c
    ! d(h_if97)/dT*dT/dP+d(h_if97)/dP
    hwp = hwt*dT_h_dp+hwp
    ! d(h_if97/dP*dP/dT+d(h_if97)/dT
    hwt = hwp/dT_h_dp+hwt

  endif
end subroutine EOSWaterEnthalpyDriesnerExt

! ************************************************************************** !
subroutine EOSWaterDensityDriesnerExt(T,P, aux, &
                                   calculate_derivatives, &
                                   dw, dwmol, dwp, dwt, ierr)
  !
  ! Water density calculation from Driesner (2007)
  ! doi:10.1016/j.gca.2007.05.026
  !
  ! Author: David Fukuyama
  ! Date: 05/26/21
  !
  implicit none

  PetscReal, intent(in) :: T   ! Temperature in centigrade
  PetscReal, intent(in) :: P   ! Pressure in Pascal
  PetscReal, intent(in) :: aux(*)    ! solute mole fraction
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw ! kg/m^3
  PetscReal, intent(out) :: dwmol ! kmol/m^3
  PetscReal, intent(out) :: dwp ! kmol/m^3-Pa
  PetscReal, intent(out) :: dwt ! kmol/m^3-C
  PetscErrorCode, intent(inout) :: ierr

  PetscReal, parameter :: Pa_to_bar = 1.d-5

  PetscReal :: T_v ! Scaled temperature for volumetric correlation
  PetscReal :: x,s
  PetscReal :: n1,n2,n11,n12,n20,n21,n22,n23,n1x1,n2x1
  PetscReal :: D_T, n30,n31,n300,n301,n302,n310,n311,n312

  PetscReal :: t_C ! temperature in Celcius
  PetscReal :: p_bar
  PetscReal :: brine_molar_mass

  t_C = T
  p_bar = P*Pa_to_bar

  s = aux(1) !mass frac
  brine_molar_mass = 1.d0 / (s / FMWNACL + (1.d0 - s) / FMWH2O)
  x = s * brine_molar_mass / FMWNACL
  
  n11 = -54.2958d0 - 45.7623d0 * exp(-9.44785d-4 * p_bar) ! Table 4
  n21 = -2.6142d0 - 0.000239092d0 * p_bar
  n22 = 0.0356828d0 + (4.37235d-6 + 2.0566d-9 * p_bar) * p_bar
  n1x1 = 330.47d0 + 0.942876d0 * sqrt(p_bar) + p_bar * (0.0817193d0 &
         + p_bar * (-2.47556d-8 + p_bar * 3.45052d-10)) ! eq 11
  n2x1 = -0.0370751d0 + 0.00237723d0 * sqrt(p_bar) + p_bar * (5.42049d-5 &
        + p_bar * (5.84709d-9 - p_bar * 5.99373d-13))  ! eq 12
  n12 = -n1x1 - n11
  n20 = 1.d0 - n21 * sqrt(n22)
  
  n23 = n2x1 - n20 - n21 * sqrt(1.d0 + n22)
  n1  = n1x1 + (1.d0 - x) * n11  + n12 * (1.d0 - x)**2 ! eq 9
  n2  = n20 + n21 * sqrt(x + n22) + n23 * x            ! eq 10
  n300 = 7.60664d6 / (p_bar + 472.051d0)**2
  n301 = -50.d0 - 86.1446d0*exp(-6.21128d-4 * p_bar)
  n302 = 294.318d0 * exp(-5.66735d-3 * p_bar)
  n310 = -0.0732761d0 * exp(-2.3772d-3 * p_bar)  - 5.2948d-5 * p_bar
  n311 = -47.2747d0 + 24.3653d0 * exp(-1.25533d-3 * p_bar)
  n312 = -0.278529d0 - 0.00081381d0 * p_bar
  n30 = n300 * (exp(n301*x)-1.d0) + n302*x
  n31 = n310 * (exp(n311*x)) + n312*x
  D_T = n30 * exp(n31*T)
  T_v = n1 + n2 * t_c + D_T! doesn't include D(T) correction
  call EOSWaterDensityIF97(T_v,P,calculate_derivatives,dw,dwmol,dwp,dwt,ierr)
  if (calculate_derivatives) then
    ! calculate derivatives

  else
     dwp = UNINITIALIZED_DOUBLE
     dwt = UNINITIALIZED_DOUBLE
  endif

  
  dw = dw * brine_molar_mass / FMWH2O

end subroutine EOSWaterDensityDriesnerExt

! ************************************************************************** !

subroutine EOSWaterEnthalpySparrowExt(T,P,aux,calculate_derivatives,hw,hwp,&
     hwt,ierr)
  !
  ! Calculates brine enthalpy based on Sparrow, 2003.
  ! Valid from 0<T<300C
  !
  ! Author: David Fukuyama
  ! Date: 06/01/21
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: P
  PetscReal, intent(in) :: aux(*)
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: hw,hwp,hwt
  PetscErrorCode, intent(inout) :: ierr

  PetscReal :: A, B, C, D, E
  PetscReal :: s, molal
  PetscReal, parameter :: mol_to_kmol = 1.d-3
  PetscReal, parameter :: kJ_to_J = 1.d3

  s = aux(1) ! mass fraction
  ! eq 8
  A = (0.0005d0  +s*(0.0378d0  +s*(-0.3682d0 +s*(-0.6529d0 +s*2.89d0))))*1.d3
  B = 4.145d0    +s*(-4.973d0  +s*(4.482d0   +s*(18.31d0   +s*(-46.41d0))))
  C = 0.0007d0   +s*(-0.0059d0 +s*(0.0854d0  +s*(-0.4951d0 +s*0.8255d0)))
  D = (-0.0048d0 +s*(0.0639d0  +s*(-0.714d0  +s*(3.273d0   +s*(-4.85d0)))))*1.d-3
  E = (0.0202d0  +s*(-0.2432d0 +s*(2.054d0   +s*(-8.211d0  +s*11.43d0))))*1.d-6
  hw = A+T*(B+T*(C+T*(D+T*E))) !kJ/kg

  molal = (1.d3*(s*100.d0)/(58.442d0*(1.d2 - (s*100.d0))))*1.d2 !mol/kg

  hw = hw*kJ_to_J/((molal+55.508435d0)*mol_to_kmol) !J/kmol
  if (calculate_derivatives) then
    hwt = (B+T*(C+T*(D+T*E)))*kJ_to_J/((molal+55.508435d0)*mol_to_kmol)
    hwp = 0.d0
  else
    hwt = UNINITIALIZED_DOUBLE
    hwp = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterEnthalpySparrowExt

! ************************************************************************** !

subroutine EOSWaterDensityConstant(t,p,calculate_derivatives,dw,dwmol, &
                                   dwp,dwt,ierr,table_idxs)
  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  dw = constant_density ! kg/m^3
  dwmol = dw/FMWH2O ! kmol/m^3

  dwp = 0.d0
  dwt = 0.d0

end subroutine EOSWaterDensityConstant

! ************************************************************************** !

subroutine EOSWaterEnthalpyConstant(t,p,calculate_derivatives, &
                                    hw,hwp,hwt,ierr)
  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: hw,hwp,hwt
  PetscErrorCode, intent(inout) :: ierr

  hw = constant_enthalpy ! J/kmol

  hwp = 0.d0
  hwt = 0.d0

end subroutine EOSWaterEnthalpyConstant

! ************************************************************************** !

subroutine EOSWaterDensityExpPressure(t,p,calculate_derivatives, &
                                      dw,dwmol,dwp,dwt,ierr,table_idxs)
  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  ! kg/m^3
  dw = exp_p_reference_density*exp(exp_p_water_compressibility* &
                                     (p-exp_p_reference_pressure))
  dwmol = dw/FMWH2O ! kmol/m^3

  if (calculate_derivatives) then
    dwp = dwmol*exp_p_water_compressibility !kmol/m^3/Pa
  else
    dwp = UNINITIALIZED_DOUBLE
  endif
  dwt = 0.d0

end subroutine EOSWaterDensityExpPressure

! ************************************************************************** !

subroutine EOSWaterDensityExpPressureTemp(t,p,calculate_derivatives, &
                                          dw,dwmol,dwp,dwt,ierr,table_idxs)

  ! Calculate water as an exponential function of temperature
  ! and pressure

  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  ! kg/m^3
  dw = exp_pt_reference_density*(1.d0 / (exp(exp_pt_thermal_expansion* &
       (t - exp_pt_reference_temperature)) * exp(exp_pt_water_compressibility* &
       (exp_pt_reference_pressure - p))))
  dwmol = dw/FMWH2O ! kmol/m^3

  if (calculate_derivatives) then
    dwp = dwmol*exp_pt_water_compressibility !kmol/m^3/Pa
    dwt = -dwmol*exp_pt_thermal_expansion
  else
    dwp = UNINITIALIZED_DOUBLE
    dwt = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterDensityExpPressureTemp

! ************************************************************************** !

subroutine EOSWaterDensityLinear(t,p,calculate_derivatives, &
                                      dw,dwmol,dwp,dwt,ierr,table_idxs)
  !
  ! Water density linear model
  !
  ! Author: Satish Karra
  ! Date: 06/19/17

  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  ! kg/m^3
  dw = linear_reference_density*(1.d0 + &
         linear_water_compressibility*(p-linear_reference_pressure))

  dwmol = dw/FMWH2O ! kmol/m^3

  if (calculate_derivatives) then
    dwp = linear_reference_density*linear_water_compressibility/FMWH2O
  else
    dwp = UNINITIALIZED_DOUBLE
  endif
  dwt = 0.d0

end subroutine EOSWaterDensityLinear

! ************************************************************************** !

subroutine EOSWaterDensityBRAGFLO(t,p,calculate_derivatives, &
                                  dw,dwmol,dwp,dwt,ierr,table_idxs)
  !
  ! Water density based on formulation in BRAGFLO.  The BRAGFLO user manual
  ! is incorrect as it does not include the truncation in the code (see
  ! subroutine DENO near line 6848 of Bragflo.f
  !
  ! Author: Glenn Hammond
  ! Date: 06/02/17

  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade (ignored)
  PetscReal, intent(in) :: p   ! Pressure in Pascal
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: p_adjust

  ! kg/m^3
  p_adjust = max(p,0.d0)
  dw = exp_p_reference_density* &
         exp(min(80.d0,exp_p_water_compressibility* &
                       (p_adjust-exp_p_reference_pressure)))

  dwmol = dw/FMWH2O ! kmol/m^3

  if (calculate_derivatives) then
    print *, 'Analytical derivatives in EOSWaterDensityBRAGFLO() &
      &currently not possible due to discontinuities in the formulation.'
    stop
    !dwp = dwmol*exponent_water_compressibility !kmol/m^3/Pa
  else
    dwp = UNINITIALIZED_DOUBLE
  endif
  dwt = 0.d0

end subroutine EOSWaterDensityBRAGFLO

! ************************************************************************** !

subroutine EOSWaterDensityQuadratic(t,p,calculate_derivatives, &
                                      dw,dwmol,dwp,dwt,ierr,table_idxs)
  !
  ! Water density quadratic model
  !
  ! Author: Paolo Orsini
  ! Date: 02/10/17

  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: X_pr

  ! [-] =                         1/Pa *  Pa
  X_pr = quadratic_wat_compressibility * (p - quadratic_reference_pressure )

  ! kg/m^3
  dw = quadratic_reference_density * (1.0 + X_pr + (X_pr**2)/2.0d0 )

  dwmol = dw/FMWH2O ! kmol/m^3

  if (calculate_derivatives) then
          ! kg/m^3 * 1/Pa * [-] / (kg/kmol) = kmol/m^3/Pa
    dwp = quadratic_reference_density * quadratic_wat_compressibility * &
          (1.0 + X_pr) / FMWH2O
  else
    dwp = UNINITIALIZED_DOUBLE
  endif
  dwt = 0.d0

end subroutine EOSWaterDensityQuadratic

! ************************************************************************** !

subroutine EOSWaterDensityTrangenstein(t,p,calculate_derivatives, &
                                       dw,dwmol,dwp,dwt,ierr,table_idxs)
  !
  ! Water density model taken from Trangenstein
  ! Trangenstein, J. A. “Analysis of a Model and Sequential Numerical Method
  ! for Thermal Reservoir Simulation,” The Mathematics of Oil Recovery,
  ! King, P.R. (ed.), Oxford, UK, Clarendon Press (1992), 359.
  !
  ! Author: Paolo Orsini
  ! Date: 02/22/17

  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal, parameter :: a0 = 9.99839520d+02
  PetscReal, parameter :: a1 = 1.69551760d+01
  PetscReal, parameter :: a2 = -7.98700000d-03
  PetscReal, parameter :: a3 = -4.61704610d-05
  PetscReal, parameter :: a4 = 1.05563020d-07
  PetscReal, parameter :: a5 = -2.80543530d-10
  PetscReal, parameter :: a6 = 1.68798500d-02
  PetscReal, parameter :: a7 = 10.20

  PetscReal, parameter :: cpw = 4.00d-005  !1/atm
  ! conversion 1/atm -> 1/Pa
  ! cpw_Pa = cpw * 1.01325*1.0d-1

  PetscReal :: t_numer, d_t_numer_dt, t_denom, d_t_denom_dt,  p_exponent, d_p_exponent_dp

  dw = (a0 + a1*t + a2 * t**2.0 + a3 * t**3.0 + a4 * t**4.0 + a5 * t**5.0 ) / &
       ( 1.0 + a6 * t ) * &           !Pa -> MPa
       dexp( cpw / (1.01325*1.0d-1) * (p * 1.0d-6 - a7 ) )
               !1/atm -> 1/MPa

  dwmol = dw/FMWH2O ! kmol/m^3

  if (calculate_derivatives) then

    !!! DS - FIRST ATTEMPT at these derivatives
    t_numer = (a0 + a1*t + a2 * t**2.0 + a3 * t**3.0 + a4 * t**4.0 + a5 * t**5.0 )
    d_t_numer_dt = (a1 + 2.0* a2 * t +  3.0 * a3 * t**2.0 + 4.0 * a4 * t**3.0 + 5.0 * a5 * t**4.0 )

    t_denom = ( 1.0 + a6 * t )
    d_t_denom_dt = a6

    p_exponent = cpw / (1.01325*1.0d-1) * (p * 1.0d-6 - a7 )
    d_p_exponent_dp = cpw / (1.01325*1.0d-1) * 1.0d-6

    dwt = (d_t_numer_dt / t_denom) - (t_numer * d_t_denom_dt / t_denom / t_denom  )
    dwt = dwt * dexp(p_exponent)

    dwp = dw * d_p_exponent_dp

    !! derivatives are undestood to be in mol units so
    dwt = dwt/FMWH2O
    dwp = dwp/FMWH2O

    !PO TODO add derivatives
    !print *, 'Derivatives not set up in EOSWaterDensityTrangenstein.'
    !stop
  else
    dwp = UNINITIALIZED_DOUBLE
    dwt = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterDensityTrangenstein

! ************************************************************************** !

subroutine EOSWaterViscosityGrabowski(T, P, PS, dPS_dT, &
                                      calculate_derivatives, VW, &
                                      dVW_dT, dVW_dP, ierr,table_idxs)
  !
  ! Grabowski, J. W. and Rubin, B. “A Preliminary Numerical Simulation
  ! Study of In-situ Combustion in a Cold Lake Oil Sands Reservoir”
  ! Journal of Canadian Petroleum Technology, (1981) 20, No. 2, 79-89.
  !
  ! Author: Paolo Orsini
  ! Date: 02/26/17
  !
  implicit none
  PetscReal, intent(in) :: T       ! C
  PetscReal, intent(in) :: P       ! Pa
  PetscReal, intent(in) :: PS      ! Pa
  PetscReal, intent(in) :: dPS_dT  ! Pa/C
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: VW     ! Pa-s
  PetscReal, intent(out) :: dVW_dT, dVW_dP
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  ! convert from centipoise to Pa-s (1 cP = 1.d-3 Pa-s)
  PetscReal, parameter :: centipoise_to_Pa_s = 1.d-3
  PetscReal, parameter :: a = 2.1850
  PetscReal, parameter :: b = 0.04012
  PetscReal, parameter :: c = 5.1547d-6

  PetscReal :: t_F ! temperature in F
  PetscReal :: denom, d_denom_d_t_F, d_VW_d_t_F, d_t_F_d_T

  t_F = T*(9.0/5.0)+32.0

  VW = a / (-1.0  + b * t_F + c * t_F**2.0 )
  VW = VW * centipoise_to_Pa_s

  if (calculate_derivatives) then

    !!! DS - FIRST ATTEMPT at these derivatives
    denom = (-1.0  + b * t_F + c * t_F**2.0 )
    d_denom_d_t_F = ( b + 2.0 * c * t_F )
    !! simple chain rule for this:
    d_VW_d_t_F = -1.0 * a * d_denom_d_t_F / denom / denom

    d_t_F_d_T = 9.0/5.0
    !! and again simple chain rule:
    dVW_dT = d_t_F_d_T * d_VW_d_t_F
    !! scaling needed:
    dVW_dT = dVW_dT * centipoise_to_Pa_s



    dVW_dP = 0.d0

    !PO TODO add derivatives
    !print *, 'Derivatives not set up in EOSWaterViscosityGrabowski'
    !stop
    !dVW_dT =
  else
    dVW_dP = UNINITIALIZED_DOUBLE
    dVW_dT = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterViscosityGrabowski

! ************************************************************************** !

subroutine EOSWaterDensityTPPlanarSetup()
  !
  ! Setups up plane for density of water as a function of temperature and
  ! pressure using a simple plane equation
  !
  ! Author: Glenn Hammond
  ! Date: 1/31/17
  !
  implicit none

  PetscReal, parameter :: p0 = 1.d6
  PetscReal, parameter :: t0 = 30.d0
  PetscReal, parameter :: drho_dp = 4.42d-7
  PetscReal, parameter :: drho_dT = -0.312d0
  PetscReal, parameter :: rho_reference = 996.d0

  call GeomComputePlaneWithGradients(water_density_tp_plane,p0,t0, &
                                     rho_reference,drho_dp,drho_dT)

end subroutine EOSWaterDensityTPPlanarSetup

! ************************************************************************** !

subroutine EOSWaterDensityTPPlanar(t,p,calculate_derivatives, &
                                   dw,dwmol,dwp,dwt,ierr,table_idxs)
  !
  ! Calculates the density of water as a function of temperature and pressure
  ! using a simple plane equation
  !
  ! Author: Glenn Hammond
  ! Date: 11/04/16
  !
  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  ! kg/m^3
  dw = GeometryGetPlaneZIntercept(water_density_tp_plane,p,t)
  dwmol = dw/FMWH2O ! kmol/m^3

  if (calculate_derivatives) then
    call GeomGetPlaneGradientinXandY(water_density_tp_plane,dwp,dwt)
    dwp = dwp/FMWH2O
    dwt = dwt/FMWH2O
  else
    dwp = UNINITIALIZED_DOUBLE
    dwt = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterDensityTPPlanar

! ************************************************************************** !

subroutine EOSWaterEnthalpyTPPlanarSetup()
  !
  ! Setups up plane for enthalpy of water as a function of temperature and
  ! pressure using a simple plane equation
  !
  ! Author: Glenn Hammond
  ! Date: 1/31/17
  !
  implicit none

  PetscReal, parameter :: p0 = 1.d6
  PetscReal, parameter :: t0 = 30.d0
  PetscReal, parameter :: dh_dp = 1.65d-2
  PetscReal, parameter :: dh_dT = 7.5d4
  PetscReal, parameter :: h_reference = 2.27d6 ! J/kmol

  call GeomComputePlaneWithGradients(water_enthalpy_tp_plane,p0,t0, &
                                     h_reference,dh_dp,dh_dT)

end subroutine EOSWaterEnthalpyTPPlanarSetup

! ************************************************************************** !

subroutine EOSWaterEnthalpyTPPlanar(t,p,calculate_derivatives, &
                                    hw,hwp,hwt,ierr)
  !
  ! Calculates the enthalpy of water as a function of temperature and pressure
  ! using a simple plane equation
  !
  ! Author: Glenn Hammond
  ! Date: 11/04/16
  !
  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: hw,hwp,hwt
  PetscErrorCode, intent(inout) :: ierr

  ! kg/m^3
  hw = GeometryGetPlaneZIntercept(water_enthalpy_tp_plane,p,t)

  if (calculate_derivatives) then
    call GeomGetPlaneGradientinXandY(water_enthalpy_tp_plane,hwp,hwt)
  else
    hwp = UNINITIALIZED_DOUBLE
    hwt = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterEnthalpyTPPlanar

! ************************************************************************** !

subroutine EOSWaterSteamDenEnthNoDerive(t,p,dg,dgmol,hg,ierr)

  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade.
  PetscReal, intent(in) :: p   ! Vapor Pressure in Pascals.
  PetscReal, intent(out) :: dg,dgmol
  PetscReal, intent(out) :: hg
  PetscErrorCode, intent(inout) :: ierr

  PetscReal :: dum1, dum2, dum3, dum4

  call EOSWaterSteamDensityEnthalpyPtr(t,p,PETSC_FALSE,dg,dgmol,hg, &
                                       dum1,dum2,dum3,dum4,ierr)

end subroutine EOSWaterSteamDenEnthNoDerive

! ************************************************************************** !

subroutine EOSWaterSteamDenEnthDerive(t,pv,dg,dgmol,hg, &
                                      dgp,dgt,hgp,hgt,ierr)
  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade.
  PetscReal, intent(in) :: pv  ! Vapor Pressure in Pascals.
  PetscReal, intent(out) :: dg,dgmol,dgp,dgt
  PetscReal, intent(out) :: hg,hgp,hgt
  PetscErrorCode, intent(inout) :: ierr

  call EOSWaterSteamDensityEnthalpyPtr(t,pv,PETSC_TRUE,dg,dgmol,hg,dgp, &
                                       dgt,hgp,hgt,ierr)

end subroutine EOSWaterSteamDenEnthDerive

! ************************************************************************** !

subroutine EOSWaterSteamDensityEnthalpyIFC67(t,pv,calculate_derivatives, &
                                             dg,dgmol,hg, &
                                             dgp,dgt,hgp,hgt,ierr)
! t/C  p/Pa dgmol/(mol/m^3)  h/J/kmol
  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade.
  PetscReal, intent(in) :: pv  ! Vapor Pressure in Pascals.
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dg,dgmol,dgp,dgt
  PetscReal, intent(out) :: hg,hgp,hgt
  PetscErrorCode, intent(inout) :: ierr

  PetscInt, save :: n(8),ll(8),x(8,2),z(8,3)
  PetscReal, save :: xi1,xl0,xl1,xl2
  PetscReal, save :: b(8,2),bb(0:9,0:6)
  PetscReal :: sumbx(8),sumbxt(8)

  PetscInt :: i,j
  PetscReal, save :: delp,delt
  PetscReal :: beta,betap,betal,betalp,betalp1,betal1,ubeta,theta, &
            thetap,utheta
  PetscReal :: xx,xxp
  PetscReal :: f1,f2,fim1,fim2,sum,sumt,sum1,sum1t,sum1p, &
            sum2,sum2t
  PetscReal :: u1,u1t,u1p,u2,u2p,u2t,u3,u3p,u3t,v1,v1t
  PetscReal :: term,term1,term1t,term1p,term2,term2t,term2p, &
            term3,term3t,term3p,term4,term4t,term4p,term5,term5t,term5p
  PetscReal :: hr,hrp,hrt,hrpt,hrpp
  PetscReal :: vr,vrpt,vrpp
  PetscReal :: tc1,pc1,vc1,utc1,upc1,vc1mol
  PetscReal, parameter :: zero = 0.d0
  PetscReal, parameter :: one = 1.d0
  PetscReal, parameter :: two = 2.d0
  PetscReal, parameter :: three = 3.d0
  PetscReal, parameter :: four = 4.d0
  PetscReal, parameter :: five = 5.d0
  PetscReal, parameter :: six = 6.d0
  PetscReal, parameter :: ten = 10.d0

  data delt,delp/1.d-6,1.d-6/

  data n/2,3,2,2,3,2,2,2/
  data ll/0,0,0,0,0,1,1,2/
  data x/0,0,0,0,0,14,19,54, &
          0,0,0,0,0, 0, 0,27/
  data z/13,18,18,25,32,12,24,24, &
          3, 2,10,14,28,11,18,14, &
          0, 1, 0, 0,24, 0, 0, 0/

  data b/7.6333333333d-1,0.d0,0.d0,0.d0,0.d0,4.006073948d-1, &
          8.636081627d-2,-8.532322921d-1, &
          0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,3.460208861d-1/

!---sub-region 2
  data bb/ &
!---bb(0-9,0)
          1.683599274d1,  8*0.d0        , 1.936587558d2, &
!---bb(0-9,1)
          2.856067796d01, 6.670375918d-2, 8.390104328d-2, &
          4.520918904d-1,-5.975336707d-1, 5.958051609d-1, &
          1.190610271d-1, 1.683998803d-1, 6.552390126d-3, &
          -1.388522425d+3, &
!---bb(0-9,2)
          -5.438923329d01, 1.388983801d0 , 2.614670893d-2, &
          1.069036614d-1,-8.847535804d-2,-5.159303373d-1, &
          -9.867174132d-2,-5.809438001d-2, 5.710218649d-4, &
          4.126607219d03, &
!---bb(0-9,3)
          4.330662834d-1, 0.d0          ,-3.373439453d-2, &
          0.d0          , 0.d0          , 2.075021122d-1, &
          0.d0          , 0.d0          , 0.d0          , &
          -6.508211677d03, &
!---bb(0-9,4)
          -6.547711697d-1, 8*0.d0        , 5.745984054d03, &
!---bb(0-9,5)
          8.565182058d-2, 8*0.d0        ,-2.693088365d03, &
!---bb(0-9,6)
                          9*0.d0        , 5.235718623d02/

  data xi1/4.260321148d0/
  data xl0,xl1,xl2/15.74373327d0,-34.17061978d0,19.31380707d0/

  tc1 = H2O_CRITICAL_TEMPERATURE       ! K
  pc1 = H2O_CRITICAL_PRESSURE        ! Pa
  vc1 = 0.00317d0     ! m^3/kg
  utc1 = one/tc1
  upc1 = one/pc1
  vc1mol = vc1*FMWH2O ! m^3/kmol

  theta  = (t+273.15d0)*utc1
  beta   = pv*upc1
  ubeta  = one/beta
  utheta = one/theta
  xx = exp(b(1,1)*(one-theta))

!---compute steam density and derivatives

  term1 = xi1*theta*ubeta
  term1t = xi1*ubeta
  term1p = -term1*ubeta

  do i = 1,8
    sum = zero
    sumt = zero
    do j = 1,n(i)
      u1 = bb(i,j)*xx**z(i,j)
      sumt = sumt-b(1,1)*z(i,j)*u1
      sum = sum + u1
    end do
    sumbx(i) = sum
    sumbxt(i) = sumt
  end do

  term2  = sumbx(1)+beta*(two*sumbx(2)+beta*(three*sumbx(3) &
            +beta*(four*sumbx(4)+five*beta*sumbx(5))))
  term2t = sumbxt(1)+beta*(two*sumbxt(2)+beta*(three*sumbxt(3) &
            +beta*(four*sumbxt(4)+five*beta*sumbxt(5))))
  term2p = two*sumbx(2)+beta*(six*sumbx(3)+beta*(12.d0*sumbx(4) &
            +20.d0*beta*sumbx(5)))

  term3  = zero
  term3p = zero
  term3t = zero

  do i = 6,8
    fim1 = dfloat(i-1)
    fim2 = fim1-one

    sum2 = zero
    sum2t = zero
    do j = 1,ll(i)
      u1 = b(i,j)*xx**x(i,j)
      sum2t = sum2t-b(1,1)*x(i,j)*u1
      sum2 = sum2 + u1
    end do

    f1 = fim2*beta**(1-i)
    f2 = beta**(2-i)
    u1 = one/(f2+sum2)**2
    term = u1*f1*sumbx(i)
    term3 = term3+term
    u2 = two*term/(f2+sum2)
    term3t = term3t+u1*(f1*sumbxt(i)+u2*sum2t)
    term3p = term3p-u1*(f1*fim1*sumbx(i)+u2*fim2*f2)*ubeta
  end do

  term4 = bb(9,0)
  term4t = zero
  do i = 1,6
    u1 = bb(9,i)*xx**i
    term4t = term4t-b(1,1)*dfloat(i)*u1
    term4 = term4 + u1
  end do

  betal = xl0+theta*(xl1+xl2*theta)
  betalp = xl1+ two*xl2*theta
  u1 = 11.d0*(beta/betal)**10
  term4t = u1*(term4t-10.d0*betalp*term4/betal)
  term4 = term4*u1
  term4p = 10.d0*term4*ubeta

  vr = term1-term2-term3+term4
  vrpt = term1t-term2t-term3t+term4t
  vrpp = term1p-term2p-term3p+term4p

  u1 = -one/(vc1mol*vr)
  dgmol = -u1
  dg = one/(vc1*vr)

  if (calculate_derivatives) then
    u1 = u1/vr
    dgt = u1*vrpt*utc1
    dgp = u1*vrpp*upc1
  else
    dgt = UNINITIALIZED_DOUBLE
    dgp = UNINITIALIZED_DOUBLE
  endif

!---compute steam enthalpy, internal energy, and derivatives
  thetap = theta+delt
  betap  = beta +delp
  xxp    = exp(b(1,1)*(one-thetap))

  term1  = bb(0,0)*theta
  term1t = bb(0,0)*thetap

  term2  = zero
  term2t = zero
  do i = 1,5
    u1 = bb(0,i)*(dfloat(i)-two)
    term2t = term2t+u1*thetap**(i-1)
    term2  = term2+u1*theta**(i-1)
  end do

  term3  = zero
  term3t = zero
  term3p = zero
  u1 = one
  u1p = one
  do i = 1,5
    u1 = u1*beta
    u1p = u1p*betap

    sum = zero
    sumt = zero
    do j = 1,n(i)
      sumt = sumt+bb(i,j)*(one+z(i,j)*b(1,1)*thetap)*xxp**z(i,j)
      sum  = sum+bb(i,j)*(one+z(i,j)*b(1,1)*theta)*xx**z(i,j)
    end do

    term3t = term3t+u1*sumt
    term3p = term3p+u1p*sum
    term3  = term3+u1*sum
  end do

  term4  = zero
  term4t = zero
  term4p = zero
  do i = 6,8

    sum1  = zero
    sum2  = zero
    sum1t = zero
    sum2t = zero

    do j = 1,ll(i)
      u1 = b(i,j)*xxp**x(i,j)
      sum1t = sum1t+x(i,j)*u1
      sum2t = sum2t+u1
      u1 = b(i,j)*xx**x(i,j)
      sum1 = sum1+x(i,j)*u1
      sum2 = sum2+u1
    end do

    u1 = one/(beta**(2-i)+sum2)
    u2 = one-b(1,1)*theta*sum1*u1
    u3 = b(1,1)*theta

    u1t = one/(beta**(2-i)+sum2t)
    u2t = one-b(1,1)*thetap*sum1t*u1t
    u3t = b(1,1)*thetap

    u1p = one/(betap**(2-i)+sum2)
    u2p = one-b(1,1)*theta*sum1*u1p
    u3p = u3

    sum1 = zero
    sum1t = zero
    sum1p = zero
    do j = 1,n(i)
      sum1t = sum1t + bb(i,j)*xxp**z(i,j)*(u2t+z(i,j)*u3t)
      sum1p = sum1p + bb(i,j)*xx **z(i,j)*(u2p+z(i,j)*u3p)
      sum1  = sum1  + bb(i,j)*xx **z(i,j)*(u2 +z(i,j)*u3)
    end do

    term4t = term4t+sum1t*u1t
    term4p = term4p+sum1p*u1p
    term4  = term4+sum1*u1
  end do

  u1 = ten*betalp/betal
  term5 = (one+theta*u1)*bb(9,0)

  betal1  = xl0+thetap*(xl1+xl2*thetap)
  betalp1 = xl1+two*xl2*thetap
  u1t     = ten*betalp1/betal1
  term5t  = (one+thetap*u1t)*bb(9,0)

  do i = 1,6
    v1     = one+theta*(u1+dfloat(i)*b(1,1))
    v1t    = one+thetap*(u1t+dfloat(i)*b(1,1))
    term5t = v1t*bb(9,i)*xxp**i+term5t
    term5  =  v1*bb(9,i)*xx **i+term5
  end do

  term5  = term5*beta*(beta/betal)**10
  term5t = term5t*beta*(beta/betal1)**10
  term5p = term5*(betap/beta)**11

  hr   = term1 -term2 -term3 -term4 +term5
  hrt  = term1t-term2t-term3t-term4t+term5t
  hrp  = term1 -term2 -term3p-term4p+term5p
  hrpt = (hrt-hr)/delt
  hrpp = (hrp-hr)/delp

  v1 = pc1*vc1mol  ! Pa = (nRT/V) = J/m^3 : J/m^3 * m^3/kmol = J/kmol
  hg = hr*v1       ! J/kmol

  if (calculate_derivatives) then
    hgt = hrpt*v1*utc1
    hgp = hrpp*vc1mol
  else
    hgt = UNINITIALIZED_DOUBLE
    hgp = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterSteamDensityEnthalpyIFC67

! ************************************************************************** !

subroutine EOSWaterSteamDenEnthConstant(t,pv,calculate_derivatives, &
                                        dg,dgmol,hg,dgp,dgt,hgp,hgt,ierr)
! t/C  p/Pa dgmol/(mol/m^3)  h/J/kmol
  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade.
  PetscReal, intent(in) :: pv  ! Vapor Pressure in Pascals.
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dg,dgmol,dgp,dgt
  PetscReal, intent(out) :: hg,hgp,hgt
  PetscErrorCode, intent(inout) :: ierr

  dg = constant_steam_density
  dgmol = dg/FMWH2O
  hg = constant_steam_enthalpy

  if (calculate_derivatives) then
    dgt = 0.d0
    dgp = 0.d0
    hgt = 0.d0
    hgp = 0.d0
  else
    dgt = UNINITIALIZED_DOUBLE
    dgp = UNINITIALIZED_DOUBLE
    hgt = UNINITIALIZED_DOUBLE
    hgp = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterSteamDenEnthConstant

! ************************************************************************** !

subroutine EOSWaterSteamDenEnthTPPlanarSetup()
  !
  ! Setups up plane for density adn enthalpy of steam as a function of
  ! temperature and pressure using a simple plane equation
  !
  ! Author: Glenn Hammond
  ! Date: 1/31/17
  !
  implicit none

  PetscReal, parameter :: p0 = 1.d5
  PetscReal, parameter :: t0 = 30.d0
  PetscReal, parameter :: drho_dp = 7.8d-6
  PetscReal, parameter :: drho_dT = -3.54d-3
  PetscReal, parameter :: rho_reference = 0.846d0 ! kg/m^3
  PetscReal, parameter :: dh_dp = -5.2d0
  PetscReal, parameter :: dh_dT = 4.d4
  PetscReal, parameter :: h_reference = 4.54d7 ! J/kmol

  call GeomComputePlaneWithGradients(steam_density_tp_plane,p0,t0, &
                                     rho_reference,drho_dp,drho_dT)
  call GeomComputePlaneWithGradients(steam_enthalpy_tp_plane,p0,t0, &
                                     h_reference,dh_dp,dh_dT)

end subroutine EOSWaterSteamDenEnthTPPlanarSetup

! ************************************************************************** !

subroutine EOSWaterSteamDenEnthTPPlanar(t,pv,calculate_derivatives, &
                                        dv,dvmol,hv,dvp,dvt,hvp,hvt,ierr)
! t/C  p/Pa dgmol/(mol/m^3)  h/J/kmol
                                          !
  ! Calculates the enthalpy of water as a function of temperature and pressure
  ! using a simple plane equation
  !
  ! Author: Glenn Hammond
  ! Date: 11/04/16
  !
  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade.
  PetscReal, intent(in) :: pv  ! Vapor Pressure in Pascals.
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dv,dvmol,dvp,dvt ! kmol/m^3
  PetscReal, intent(out) :: hv,hvp,hvt       ! J/kmol
  PetscErrorCode, intent(inout) :: ierr

  ! FMWAIR = 28.96d0

  ! kg/m^3
  dv = GeometryGetPlaneZIntercept(steam_density_tp_plane,pv,t)
  dvmol = dv/FMWH2O ! kmol/m^3
  hv = GeometryGetPlaneZIntercept(steam_enthalpy_tp_plane,pv,t)

  if (calculate_derivatives) then
    call GeomGetPlaneGradientinXandY(steam_density_tp_plane,dvp,dvt)
    dvp = dvp/FMWH2O
    dvt = dvt/FMWH2O
    call GeomGetPlaneGradientinXandY(steam_enthalpy_tp_plane,hvp,hvt)
  else
    dvp = UNINITIALIZED_DOUBLE
    dvt = UNINITIALIZED_DOUBLE
    hvp = UNINITIALIZED_DOUBLE
    hvt = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterSteamDenEnthTPPlanar

! ************************************************************************** !

subroutine EOSWaterSteamDensityEnthalpyIF97(T, Pv, calculate_derivatives, &
                                             dg, dgmol, hg, dgp, dgt, hgp, &
                                             hgt, ierr)

  !Author: Michael Nole
  !Date: 01/20/19
  !Superheated steam EOS from IF97
  !Region 2, valid on: {273.15K <= T <= 623.15K; 0 < P < Ps(T)}
  !                    {623.15K < T <= 863.15K; 0 < P <= P(T) 2-3 boundary fn}
  !                    {863.15K < T <= 1073.15K; 0 < P < 100MPa}
  !Region 5, valid on: {1073.15K < T < 2273.15K; 0 < P < 50MPa}
  ! IAPWS R7-97(2012)

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: Pv
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dg,dgmol,dgp,dgt
  PetscReal, intent(out) :: hg,hgp,hgt
  PetscErrorCode, intent(inout) :: ierr

  PetscReal, parameter :: p_ref = 1.d6 !16.53d6 !Pa
  PetscReal, parameter :: T_ref = 540d0 !1386  ! K
  PetscReal, parameter :: R = 0.461526d0 ! kJ/kg-K
  PetscReal, parameter :: Tf = 273.15d0
  PetscReal :: pi, tao, T_temp
  PetscReal :: gamma_0_pi, gamma_r_pi, gamma_0_tao, gamma_r_tao
  ! region 2, Tables 10 and 11
  PetscReal, parameter :: n_i0(9) = [-0.96927686500217d1, 0.10086655968018d2, &
          -0.56087911283020d-2,  0.71452738081455d-1, -0.40710498223928d0, &
          0.14240819171444d1, -0.43839511319450d1, -0.28408632460772d0, &
          0.21268463753307d-1]
  PetscReal, parameter :: n_i(43) = [-0.17731742473213d-2, &
          -0.17834862292358d-1, -0.45996013696365d-1, -0.57581259083432d-1, &
          -0.50325278727930d-1, -0.33032641670203d-4, -0.18948987516315d-3, &
          -0.39392777243355d-2, -0.43797295650573d-1, -0.26674547914087d-4, &
          0.20481737692309d-7, 0.43870667284435d-6, -0.32277677238570d-4, &
          -0.15033924542148d-2, -0.40668253562649d-1, -0.78847309559367d-9, &
          0.12790717852285d-7, 0.48225372718507d-6, 0.22922076337661d-5, &
          -0.16714766451061d-10, -0.21171472321355d-2, -0.23895741934104d2, &
          -0.59059564324270d-17, -0.12621808899101d-5, -0.38946842435739d-1, &
          0.11256211360459d-10, -0.82311340897998d1, 0.19809712802088d-7, &
          0.10406965210174d-18, -0.10234747095929d-12, -0.10018179379511d-8, &
          -0.80882908646985d-10, 0.10693031879409d0, -0.33662250574171d0, &
          0.89185845355421d-24, 0.30629316876232d-12, -0.42002467698208d-5, &
          -0.59056029685639d-25, 0.37826947613457d-5, -0.12768608934681d-14, &
          0.73087610595061d-28, 0.55414715350778d-16, -0.94369707241210d-6]
  PetscInt, parameter :: J_i0(9) = [0,1,-5,-4,-3,-2,-1,2,3]
  PetscInt, parameter :: I_i(43) = [1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,5,6, &
                                    6,6,7,7,7,8,8,9,10,10,10,16,16,18,20,20, &
                                    20,21,22,23,24,24,24]
  PetscInt, parameter :: J_i(43) = [0,1,2,3,6,1,2,4,7,36,0,1,3,6,35,1,2,3,7, &
                                    3,16,35,0,11,25,8,36,13,4,10,14,29,50,57 &
                                    ,20,35,48,21,53,39,26,40,58]
  ! region 5, Tables 37 and 38
  PetscReal, parameter :: n5_i0(6) = [-0.13179983674201D2, &
       0.68540841634434D1, -0.24805148933466D-1, 0.36901534980333d0, &
       -0.31161318213925D1, -0.32961626538917d0]
  PetscReal, parameter :: n5_i(6) = [0.15736404855259D-2, 0.90153761673944D-3, &
       -0.50270077677648D-2, 0.22440037409485D-5, -0.41163275453471D-5, &
       0.37919454822955D-7]
  PetscInt, parameter :: J5_i0(6) = [0,1,-3,-2,-1,2]
  PetscInt, parameter :: I5_i(6) = [1,1,1,2,2,3]
  PetscInt, parameter :: J5_i(6) = [1,2,3,3,9,7]

  T_temp = T+Tf
  pi = Pv/p_ref
  tao = T_ref/T_temp

  if (T_temp <= 1073.15d0) then
    ! region 2, tables 13 & 14
    gamma_0_pi = 1.d0/pi
    gamma_r_pi = sum(n_i*I_i*(pi**(I_i-1))*(tao-0.5)**J_i)
    gamma_0_tao = sum(n_i0*J_i0*tao**(J_i0-1))
    gamma_r_tao = sum(n_i*pi**(I_i)*J_i*(tao-0.5)**(J_i-1))
    dg = R*T_temp/Pv * pi * (gamma_0_pi + gamma_r_pi) *1.d3
    hg = R*T_temp * tao * (gamma_0_tao + gamma_r_tao)

    dg = 1/dg
    dgmol = dg/FMWH2O

    if (calculate_derivatives) then
      dgt = R/p_ref*((gamma_0_pi+gamma_r_pi)-T_ref/(T_temp)* &
           sum(n_i*I_i*pi**(I_i-1)*J_i*(tao-0.5d0)**(J_i-1)))
      dgp = R*T_temp/(p_ref*p_ref)*(-1.d0/(pi*pi)+ &
           sum(n_i*I_i*(I_i-1)*pi**(I_i-2)*(tao-0.5d0)**(J_i)))
      hgt = -tao*tao*R*(sum(n_i0*J_i0*(J_i0-1)*tao**(J_i0-2))+ &
           sum(n_i*pi**(I_i)*J_i*(J_i-1)*(tao-0.5d0)**(J_i-2)))
      hgp = T_ref*R/p_ref *(sum(n_i*I_i*pi**(I_i-1)*J_i*(tao-0.5d0)**(J_i-1)))

      hgt = hgt*FMWH2O * 1.d3
      hgp = hgp*FMWH2O * 1.d3
      dgt = -dgt*dg*dg*1.d3/FMWH2O
      dgp = -dgp*dg*dg*1.d3/FMWH2O
    else
      dgt = UNINITIALIZED_DOUBLE
      dgp = UNINITIALIZED_DOUBLE
      hgt = UNINITIALIZED_DOUBLE
      hgp = UNINITIALIZED_DOUBLE
    endif
    hg = hg*FMWH2O * 1.d3
  else if (T_temp <= 2273.15d0) then
    ! region 5: table 40
    gamma_0_pi = 1.d0/pi
    gamma_0_tao = sum(n5_i0*J5_i0*tao**(J5_i0-1))
    ! table 41
    gamma_r_pi = sum(n5_i*I5_i*(pi**(I5_i-1))*tao**J5_i)
    gamma_r_tao = sum(n5_i*pi**(I5_i)*J5_i*tao**(J5_i-1))
    dg = R*T_temp/Pv * pi * (gamma_0_pi + gamma_r_pi) * 1.d3
    hg = R*T_temp * tao * (gamma_0_tao + gamma_r_tao)

    dg = 1/dg
    if (calculate_derivatives) then
    stop 'IF97 region 5 steam density/enthalpy derivative not implemented yet'
    else
      dgt = UNINITIALIZED_DOUBLE
      dgp = UNINITIALIZED_DOUBLE
      hgt = UNINITIALIZED_DOUBLE
      hgp = UNINITIALIZED_DOUBLE
    end if
  else
    stop 'wow. much hotness.'
  end if

end subroutine EOSWaterSteamDensityEnthalpyIF97

subroutine EOSWaterDuanMixture(t,p,xmol,y_nacl,avgmw,dw_kg,denmix)

! Duan et al. (2008) Energy and Fuels, v 22, 1666-1674.

  implicit none

  PetscReal :: t,tk,p,xco2,xmol,x1,y_nacl,vphi_a1,vphi_a2,vphi,denmix,pw_kg,dw_kg,avgmw

  PetscReal :: fmwh2o = 18.01534d0
  PetscReal :: fmwco2 = 44.0098d0
  PetscReal :: fmwnacl = 58.44277d0
  PetscReal :: dummy
  PetscErrorCode :: ierr

  !duan mixing **************************
  tk = t + 273.15D0; xco2 = xmol;
  call EOSWaterDensity(t,p,pw_kg,dummy,ierr)
  x1 = 1.D0-xco2;
  vphi_a1 = (0.3838402D-3*tk - 0.5595385D0)*tk + 0.30429268D3 + &
            (-0.72044305D5 + 0.63003388D7/tk)/tk;
  vphi_a2 = (-0.57709332D-5*tk + 0.82764653D-2)*tk - 0.43813556D1 + &
            (0.10144907D4 - 0.86777045D5/tk)/tk;
  vphi = (1.D0 + vphi_a1 + vphi_a2*p*1.D-6)*(fmwh2o*1.D-3/pw_kg);
  vphi = x1*((1.D0 - y_nacl)*fmwh2o + y_nacl*fmwnacl)*1.D-3/dw_kg + xco2*vphi;
  denmix = (x1*((1.D0 - y_nacl)*fmwh2o + y_nacl*fmwnacl) + xco2*fmwco2)*1.D-3/vphi;
  denmix = denmix/avgmw

end subroutine EOSWaterDuanMixture

! ************************************************************************** !

subroutine EOSWaterViscosityNaCl (t,p_Pa,xnacl,visnacl)

  !viscosity: Kestin et al. (1981)

  implicit none

  PetscReal, intent(in) :: t        ! [C]
  PetscReal, intent(in) :: p_Pa     ! [Pa]
  PetscReal, intent(in) :: xnacl    ! [-]
  PetscReal, intent(out) :: visnacl ! [Pa-s]


  PetscReal, save :: a1,a2,a3,b1,b2,b3,c1,c2,c3,c4,wnacl
  PetscReal :: ak,bk,ck
  PetscReal :: beta,betap,betas,betaw
  PetscReal :: tt,mnacl,fac,mu0,ms
  PetscReal :: p_GPa

  data a1,a2,a3 / 3.324d-2, 3.624d-3, -1.879d-4 /
  data b1,b2,b3 / -3.96d-2, 1.02d-2, -7.02d-4 /
  data c1,c2,c3,c4 / 1.2378d0, -1.303d-3, 3.06d-6, 2.55d-8 /

  data wnacl / 58.44277d-3 / ! (kg/mol NaCl)

  !convert pressure to GPa
  p_GPa = p_Pa*1.d-9

  mnacl = xnacl/(1.d0-xnacl)/wnacl

  tt = 20.d0-t
  ck = (c1 + (c2 + (c3+c4*tt)*tt)*tt)*tt/(96.d0+t)
  ak = (a1 + (a2 + a3*mnacl)*mnacl)*mnacl
  bk = (b1 + (b2 + b3*mnacl)*mnacl)*mnacl

  ms = 6.044d0 + (2.8d-3 + 3.6d-5*t)*t
  fac = mnacl/ms
  betaw = -1.297d0 + (5.74d-2 + (-6.97d-4 + (4.47d-6 - 1.05d-8*t)*t)*t)*t
  betas = 0.545d0 + 2.8d-3 * t - betaw
  betap = (2.5d0 + (-2.d0 + 0.5d0*fac)*fac)*fac
  beta = betas*betap + betaw

  mu0 = 1001.74d-6 * 10.d0**(ak + ck*(bk + 1.d0))

  visnacl = mu0*(1.d0 + beta*p_GPa)

end subroutine EOSWaterViscosityNaCl

! ************************************************************************** !

subroutine EOSWaterViscosityKestinExt(T, P, PS, dPS_dT, aux, &
                                      calculate_derivatives, VW, &
                                      dVW_dT, dVW_dP, ierr)

  !viscosity: Kestin et al. (1981)

  implicit none

  PetscReal, intent(in) :: T   ! C
  PetscReal, intent(in) :: P   ! Pa
  PetscReal, intent(in) :: PS  ! Pa
  PetscReal, intent(in) :: dPS_dT
  PetscReal, intent(in) :: aux(*)
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: VW ! Pa-s
  PetscReal, intent(out) :: dVW_dT, dVW_dP
  PetscErrorCode, intent(inout) :: ierr

  PetscReal, save :: a1,a2,a3,b1,b2,b3,c1,c2,c3,c4,wnacl
  PetscReal :: xnacl    ! [-]
  PetscReal :: ak,bk,ck
  PetscReal :: beta,betap,betas,betaw
  PetscReal :: tt,mnacl,fac,mu0,ms
  PetscReal :: p_GPa

  data a1,a2,a3 / 3.324d-2, 3.624d-3, -1.879d-4 /
  data b1,b2,b3 / -3.96d-2, 1.02d-2, -7.02d-4 /
  data c1,c2,c3,c4 / 1.2378d0, -1.303d-3, 3.06d-6, 2.55d-8 /

  data wnacl / 58.44277d-3 / ! (kg/mol NaCl)

  if (calculate_derivatives) then
    print *, 'Derivatives not set up in EOSWaterViscosityKestinExt().'
    stop
  endif

  !convert pressure to GPa
  p_GPa = P*1.d-9

  xnacl = aux(1)
  mnacl = xnacl/(1.d0-xnacl)/wnacl
  tt = 20.d0-t
  ck = (c1 + (c2 + (c3+c4*tt)*tt)*tt)*tt/(96.d0+t)
  ak = (a1 + (a2 + a3*mnacl)*mnacl)*mnacl
  bk = (b1 + (b2 + b3*mnacl)*mnacl)*mnacl

  ms = 6.044d0 + (2.8d-3 + 3.6d-5*t)*t
  fac = mnacl/ms
  betaw = -1.297d0 + (5.74d-2 + (-6.97d-4 + (4.47d-6 - 1.05d-8*t)*t)*t)*t
  betas = 0.545d0 + 2.8d-3 * t - betaw
  betap = (2.5d0 + (-2.d0 + 0.5d0*fac)*fac)*fac
  beta = betas*betap + betaw

  mu0 = 1001.74d-6 * 10.d0**(ak + ck*(bk + 1.d0))

  VW = mu0*(1.d0 + beta*p_GPa)

end subroutine EOSWaterViscosityKestinExt

! ************************************************************************** !

subroutine EOSWaterSaturationTemperature(ps,ts_guess,ts,t_ps,ierr)

  !  This function calculates saturation temperature for a given Ps c
  !  Ref.: International Formulation Committee of the Sixth International
  !       Conference on Properties of Steam (1967).

  !    ps  = saturation pressure (pascals)
  !    ts  = saturation temperature (deg. C)
  !    tsp = estimated ts on entry and dT/dps on return

  implicit none

  PetscReal, intent(in) :: ps
  PetscReal, intent(in) :: ts_guess
  PetscReal, intent(out) :: ts
  PetscReal, intent(out) :: t_ps
  PetscErrorCode :: ierr


  PetscReal, parameter :: epsilon = 1.d-10
  PetscReal, parameter :: tc1 = H2O_CRITICAL_TEMPERATURE
  PetscReal, parameter :: pc1 = H2O_CRITICAL_PRESSURE

  PetscReal :: theta, beta, u1, err
  PetscReal :: t1num, t1nump
  PetscReal :: t1, t1den, t1denp, term1, term1p, t2, term2, term2p
  PetscReal :: f, fp
  PetscReal :: kn(9)
  PetscInt :: itr

  data kn / -7.691234564d0,-2.608023696d1,-1.681706546d2, &
            6.423285504d1,-1.189646225d2, 4.167117320d0, &
            2.097506760d1, 1.d9         , 6.d0/

!geh  if (ipvtab.eq.0 .or. tsp.gt.369.d0) then

!-------newton-raphson iteration for calculating ts by analytical funcs

!-------on entry, ts_guess = estimated ts
  theta = (ts_guess+273.15d0)/tc1
  beta  = ps/pc1

  u1  = 1.d0-theta
  itr = 0

  do
    itr = itr + 1

    t1num  = u1*(kn(1)+u1*(kn(2)+u1*(kn(3)+u1*(kn(4)+u1*kn(5)))))
    t1nump = u1*(2.d0*kn(2)+u1*(3.d0*kn(3)+u1*(4.d0*kn(4)+5.d0* &
             u1*kn(5))))+kn(1)
    t1     = 1.d0+u1*(kn(6)+kn(7)*u1)
    t1den  = 1.d0/(theta*t1)
    t1denp = theta*(kn(6)+2.d0*kn(7)*u1)-t1
    term1  = t1num*t1den
    term1p = (t1nump-term1*t1denp)*t1den

    t2     = 1.d0/(kn(8)*u1*u1+kn(9))
    term2  = u1*t2
    term2p = t2*(1.d0-2.d0*kn(8)*u1*u1*t2)
    f      = exp(term1-term2)-beta
    fp     = (f+beta)*(term1p-term2p)
    err    = f/fp
    u1     = u1-err
    theta  = 1.d0-u1

    if (dabs(err) <= epsilon .or. itr >= 20) exit

  enddo

  ts = theta*tc1-273.15d0

!-------Note-(dbeta/dtheta) = -fp ; tsp = dT/dps
  t_ps = -tc1/(pc1*fp)

end subroutine EOSWaterSaturationTemperature

! ************************************************************************** !

subroutine EOSWaterDensityIcePainter(T, P, calculate_derivatives, &
                                     den_ice, dden_ice_dT, dden_ice_dP, ierr)
  ! Subroutine to calculate the density of ice at given temperature
  ! and pressure
  ! T is in deg C, P is in Pa, density is in kmol/m3
  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: P
  PetscBool, intent(in) :: calculate_derivatives

  PetscReal, intent(out) :: den_ice
  PetscReal, intent(out) :: dden_ice_dT
  PetscReal, intent(out) :: dden_ice_dP
  PetscErrorCode, intent(inout) :: ierr

  PetscReal, parameter :: P_ref = 1.d5
  PetscReal, parameter :: alpha = 3.3d-10
  PetscReal, parameter :: beta = 1.53d-4

  den_ice = 5.09424d1*(1.d0 + alpha*(P - P_ref) - beta*(T)) !in Kmol/m3
  if (calculate_derivatives) then
    dden_ice_dT = 5.09424d1*(-beta)
    dden_ice_dP = 5.09424d1*alpha
  else
    dden_ice_dT = UNINITIALIZED_DOUBLE
    dden_ice_dP = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterDensityIcePainter

! ************************************************************************** !

subroutine EOSWaterInternalEnergyIceDefault(T, u_ice, calculate_derivatives, &
                                            du_ice_dT, du_ice_dP, ierr)
  ! Subroutine to calculate the internal energy of ice at given temperature and
  ! pressure
  ! T is in deg C, internal energy is in J/mol
  implicit none

  PetscReal, intent(in) :: T
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: u_ice
  PetscReal, intent(out) :: du_ice_dT
  PetscReal, intent(out) :: du_ice_dP
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal, parameter :: a = -10.6644d0
  PetscReal, parameter :: b = 0.1698d0
  PetscReal, parameter :: c = 198148.d0
  PetscReal, parameter :: T_ref = 273.15d0

  ! from Maier-Kelly type fit (integrated tref to t)
  ! in J/mol

  u_ice = a*(T) + b/2.d0*((T + T_ref)**(2.d0) - T_ref**(2.d0)) + &
          c*(1.d0/T_ref - 1.d0/(T + T_ref))
  u_ice = u_ice - HEAT_OF_FUSION*FMWH2O*1.d-3   ! kJ/kmol

  if (calculate_derivatives) then
    du_ice_dT = a + b*(T + T_ref) + c/((T + T_ref)**(2.d0)) !kJ/kmol/K
    du_ice_dP = 0.d0
  endif

end subroutine EOSWaterInternalEnergyIceDefault

! ************************************************************************** !

subroutine EOSWaterInternalEnergyIceFukusako(T, u_ice, calculate_derivatives, &
                                             du_ice_dT, du_ice_dP,ierr)
  
  !Internal energy of ice as f(Temperature) (Fukusako and Yamada, 1993)
  !
  ! T is in deg C, internal energy is in J/mol
  !Author: Michael Nole
  !Date: 04/04/19
  !
  implicit none

  PetscReal, intent(in) :: T
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: u_ice
  PetscReal, intent(out) :: du_ice_dT, du_ice_dP
  PetscErrorCode, intent(out) :: ierr  

  PetscReal, parameter :: Lw = -3.34110d5 ! Latent heat of fusion, J/kg 
  PetscReal :: T_temp

  T_temp = T + 273.15d0

  if (T_temp >= 90.d0) then
    u_ice = Lw + 185.d0 * (T_temp-273.15d0) + 3.445 * &
                 (T_temp**2 - 273.15d0**2) 
  else
    u_ice = Lw + 4.475 * (T_temp**2 - 273.15d0**2)
  endif

  ! J/kg to J/mol
  u_ice = u_ice * FMWH2O * 1.d-3

  if (calculate_derivatives) then
    if (T_temp >= 90.d0) then
      du_ice_dT = 185.d0 + 2 * 3.445 * T_temp
    else
      du_ice_dT = 4.475 * 2 * T_temp
    endif
    du_ice_dP = 0.d0
  endif

end subroutine EOSWaterInternalEnergyIceFukusako

! ************************************************************************** !

subroutine EOSWaterDensityPainter(t,p,calculate_derivatives,dw,dwmol, &
                                  dwp,dwt,ierr,table_idxs)

! wateos_simple: Simple water equation of state from Scott Painter
! Author: Satish Karra, LANL
! Date: 02/1/12
! T in C, P in Pa
  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal, parameter :: a = 999.915d0
  PetscReal, parameter :: b = 0.0416516d0
  PetscReal, parameter :: c = -0.0100836d0
  PetscReal, parameter :: d = 0.000206355
  PetscReal, parameter :: alpha = 5.0d-10     ! in Pa^(-1)
  PetscReal, parameter :: T_ref = 273.15d0    ! in K
  PetscReal, parameter :: P_ref = 1.0d5       ! in Pa

  PetscReal :: den_w_one_bar, T_K
  PetscReal :: u_J_kg, h_J_kg

  ! Density of water
  T_K = T + T_ref    ! convert to Kelvin
  den_w_one_bar = a + b*(T_K - T_ref) + c*(T_K - T_ref)**(2.d0) + &
                  d*(T_K - T_ref)**(3.d0)
  dw = den_w_one_bar*(1 + alpha*(P - P_ref))
  dwmol = dw/FMWH2O     ! in mol

  ! Internal energy
  u_J_kg = 4.217*1.0d3*(T_K - T_ref)    ! in J/kg
  h_J_kg = u_J_kg + P/dw    ! in J/kg

  if (calculate_derivatives) then
    ! Derivatives of density
    dwp = 1/FMWH2O*den_w_one_bar*alpha    ! in Kmol/Pa
    dwt = 1/FMWH2O*(1 + alpha*(P - P_ref))*(b + 2.d0*c*(T_K - T_ref) + &
                              3.d0*d*(T_K - T_ref)**(2.d0))      ! in Kmol/K
  else
    dwp = UNINITIALIZED_DOUBLE
    dwt = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterDensityPainter

! ************************************************************************** !

subroutine EOSWaterEnthalpyPainter(T, P, calculate_derivatives, &
                                   h_J_kmol, dh_dp, dh_dt, ierr)

! wateos_simple: Simple water equation of state from Scott Painter
! Author: Satish Karra, LANL
! Date: 02/1/12
! T in C, P in Pa
  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: P
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: h_J_kmol
  PetscReal :: den_water_kg, den_water_kmol
  PetscReal :: dden_water_dp, dden_water_dt
  PetscReal, intent(out) :: dh_dp, dh_dt

  PetscErrorCode, intent(inout) :: ierr

  PetscReal, parameter :: a = 999.915d0
  PetscReal, parameter :: b = 0.0416516d0
  PetscReal, parameter :: c = -0.0100836d0
  PetscReal, parameter :: d = 0.000206355
  PetscReal, parameter :: alpha = 5.0d-10     ! in Pa^(-1)
  PetscReal, parameter :: T_ref = 273.15d0    ! in K
  PetscReal, parameter :: P_ref = 1.0d5       ! in Pa

  PetscReal :: den_w_one_bar, T_K
  PetscReal :: u_J_kg, h_J_kg
  PetscReal :: du_dt

  ! Density of water
  T_K = T + T_ref    ! convert to Kelvin
  den_w_one_bar = a + b*(T_K - T_ref) + c*(T_K - T_ref)**(2.d0) + &
                  d*(T_K - T_ref)**(3.d0)
  den_water_kg = den_w_one_bar*(1 + alpha*(P - P_ref))
  den_water_kmol = den_water_kg/FMWH2O     ! in mol

  ! Internal energy
  u_J_kg = 4.217*1.0d3*(T_K - T_ref)    ! in J/kg
  h_J_kg = u_J_kg + P/den_water_kg    ! in J/kg
  h_J_kmol = h_J_kg*FMWH2O     ! in J/kmol

  if (calculate_derivatives) then
    ! Derivatives of density
    dden_water_dp = 1/FMWH2O*den_w_one_bar*alpha    ! in Kmol/Pa
    dden_water_dt = 1/FMWH2O*(1 + alpha*(P - P_ref))*(b + 2.d0*c*(T_K - T_ref) + &
                              3.d0*d*(T_K - T_ref)**(2.d0))      ! in Kmol/K

    ! Derivatives of enthalpy
    dh_dp = FMWH2O/den_water_kg   ! in J/kmol/Pa
    du_dt = 4.217*1.d3                  ! in J/kg/K
    dh_dt = FMWH2O*(du_dt + P*(-1.d0/den_water_kg**(2.d0))* &
                    dden_water_dt*FMWH2O)    ! in MJ/kmol/K
  else
    dden_water_dp = UNINITIALIZED_DOUBLE
    dden_water_dp = UNINITIALIZED_DOUBLE
    dh_dp = UNINITIALIZED_DOUBLE
    du_dt = UNINITIALIZED_DOUBLE
    dh_dt = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterEnthalpyPainter

! ************************************************************************** !

subroutine EOSWaterDensityIceNoDerive(t,p,dw,ierr)

  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscReal, intent(out) :: dw
  PetscErrorCode, intent(inout) :: ierr

  PetscReal :: dum1, dum2

  call EOSWaterDensityIcePtr(t,p,PETSC_FALSE,dw,dum1,dum2,ierr)

end subroutine EOSWaterDensityIceNoDerive

! ************************************************************************** !

subroutine EOSWaterDensityIceDerive(t,p,dw,dwp,dwt,ierr)
  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscReal, intent(out) :: dw,dwp,dwt
  PetscErrorCode, intent(inout) :: ierr

  call EOSWaterDensityIcePtr(t,p,PETSC_TRUE,dw,dwp,dwt,ierr)

end subroutine EOSWaterDensityIceDerive

! ************************************************************************** !

subroutine EOSWaterDensityTGDPB01(t, p, calculate_derivatives, &
                                  dw, dwmol, dwp, dwt, ierr,table_idxs)

  !
  ! Tanaka M. , G. Girard, R. Davis, A. Peuto, and N. Bignell. 2001.
  ! Recommended table for the density of water between 0 °C
  ! and 40 °C based on recent experimental reports. Metrologia,
  ! 38:301-309 [doi:10.1088/0026-1394/38/4/3].
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 24/06/15
  !
  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal,parameter :: a1 = -3.983035d0     ! [degC]
  PetscReal,parameter :: a2 = 301.797d0       ! [degC]
  PetscReal,parameter :: a3 = 522528.9d0      ! [degC^{2}]
  PetscReal,parameter :: a4 = 69.34881d0      ! [degC]
  PetscReal,parameter :: a5 = 999.974950d0    ! [kg m^{-3}]
  PetscReal,parameter :: k0 = 50.74d-11       ! [Pa^{-1}]
  PetscReal,parameter :: k1 = -0.326d-11      ! [Pa^{-1} degC^{-1}]
  PetscReal,parameter :: k2 = 0.00416d-11     ! [Pa^{-1} degC^{-2}]
  PetscReal,parameter :: p0 = 101325.d0       ! [Pa]
  PetscReal :: dent
  PetscReal :: kappa
  PetscReal :: ddent_dt
  PetscReal :: ddent_dt_1
  PetscReal :: ddent_dt_2
  PetscReal :: ddent_dt_3
  PetscReal :: ddent_dp
  PetscReal :: dkappa_dp
  PetscReal :: dkappa_dt

  ! Density of water as function of temperature
  dent = a5*(1.d0 - ((t + a1)**2.d0)*(t + a2)/a3/(t + a4))

  ! Compressibility of water
  kappa = (1.d0 + (k0 + k1*t + k2*t**2.d0)*(p - p0))

  ! Density of water
  dw    = dent*kappa ! [kg m^{-3}]
  dwmol = dw/FMWH2O  ! [kmol m^{-3}]

  if (calculate_derivatives) then
    ! Derivative
    ddent_dp = 0.d0
    ddent_dt_1 = -((t + a1)**2.d0)/a3/(t + a4)
    ddent_dt_2 = -2.d0*(t + a1)*(t + a2)/a3/(t + a4)
    ddent_dt_3 =  ((t + a1)**2.d0)*(t + a2)/a3/((t + a4)**2.d0)
    ddent_dt   = a5*(ddent_dt_1 + ddent_dt_2 + ddent_dt_3)

    dkappa_dp = (k0 + k1*t + k2*t**2.d0)
    dkappa_dt = (k1 + 2.d0*k2*t)*(p - p0)

    dwt = (ddent_dt*kappa + dent*dkappa_dt)/FMWH2O ! [kmol m^{-3} degC^{-1}]
    dwp = (ddent_dp*kappa + dent*dkappa_dp)/FMWH2O ! [kmol m^{-3} Pa^{-1}]
  else
    dwt = UNINITIALIZED_DOUBLE
    dwp = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterDensityTGDPB01

! ************************************************************************** !

subroutine EOSWaterDensityBatzleAndWang(tin, pin, calculate_derivatives, &
                                        dw, dwmol, dwp, dwt, ierr,table_idxs)

  !
  ! From Batlze M. and Z. Wang (1992) Seismic properties of fluids, Geophysics,
  ! Vol. 57, No. 11, Pg. 1396-1408.
  !
  ! Equation 27a
  !
  ! Author: Glenn Hammond
  ! Date: 02/08/16
  !
  implicit none

  PetscReal, intent(in) :: tin   ! Temperature in centigrade
  PetscReal, intent(in) :: pin   ! Pressure in Pascal
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw ! kg/m^3
  PetscReal, intent(out) :: dwmol ! kmol/m^3
  PetscReal, intent(out) :: dwp ! kmol/m^3-Pa
  PetscReal, intent(out) :: dwt ! kmol/m^3-C
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal, parameter :: g_cm3_to_kg_m3 = 1.d3
  PetscReal, parameter :: Pa_to_MPa = 1.d-6
  PetscReal :: t ! temperature in Celcius
  PetscReal :: t_sq, t_cub
  PetscReal :: p_MPa, p_MPa_sq

  t = tin
  t_sq = t*t
  t_cub = t*t_sq
  p_MPa = pin*Pa_to_MPa
  p_MPa_sq = p_MPa*p_MPa

  ! temperature is in C and pressure in MPa
  ! g/cm^3
  ! Eq. 27a
  dw = 1.d0 + 1.d-6*(-80.d0*t - 3.3d0*t_sq + 1.75d-3*t_cub + &
                     489.d0*p_MPa - 2.d0*t*p_MPa + 1.6d-2*t_sq*p_MPa - &
                     1.3d-5*t_cub*p_MPa - 3.33d-1*p_MPa_sq - 2.d-3*t*p_MPa_sq)
  ! convert from g/cm^3 to kg/m^3
  dw = dw * g_cm3_to_kg_m3
  dwmol = dw/FMWH2O ! kmol/m^3

  if (calculate_derivatives) then
    dwp = 1.d-6*(489.d0 - 2.d0*t + 1.6d-2*t_sq - 1.3d-5*t_cub - &
                 6.66d-1*p_MPa - 4.d-3*t*p_MPa) *  &
          g_cm3_to_kg_m3/FMWH2O
    ! convert from kmol/m^3-MPa to kmol/m^3-Pa
    dwp = dwp*Pa_to_MPa
    dwt = 1.d-6*(-80.d0 - 6.6d0*t + 5.25d-3*t_sq - 2.d0*p_MPa + &
                 3.2d-2*t*p_MPa - 3.9d-5*t_sq*p_MPa - 2.d-3*p_MPa_sq) * &
          g_cm3_to_kg_m3/FMWH2O
  else
    dwp = UNINITIALIZED_DOUBLE
    dwt = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterDensityBatzleAndWang

! ************************************************************************** !

subroutine EOSWaterDensityBatzleAndWangExt(tin, pin, aux, &
                                           calculate_derivatives, &
                                           dw, dwmol, dwp, dwt, ierr)

  !
  ! From Batlze M. and Z. Wang (1992) Seismic properties of fluids, Geophysics,
  ! Vol. 57, No. 11, Pg. 1396-1408.
  !
  ! Equation 27b
  !
  ! Author: Glenn Hammond
  ! Date: 02/08/16
  !
  implicit none

  PetscReal, intent(in) :: tin   ! Temperature in centigrade
  PetscReal, intent(in) :: pin   ! Pressure in Pascal
  PetscReal, intent(in) :: aux(*)
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw ! kg/m^3
  PetscReal, intent(out) :: dwmol ! kmol/m^3
  PetscReal, intent(out) :: dwp ! kmol/m^3-Pa
  PetscReal, intent(out) :: dwt ! kmol/m^3-C
  PetscErrorCode, intent(inout) :: ierr

  PetscReal, parameter :: g_cm3_to_kg_m3 = 1.d3
  PetscReal, parameter :: Pa_to_MPa = 1.d-6
  PetscReal :: t_C ! temperature in Celcius
  PetscReal :: p_MPa
  PetscReal :: s

  t_C = tin
  p_MPa = pin*Pa_to_MPa

  s = aux(1)
  call EOSWaterDensityPtr(tin, pin, calculate_derivatives, &
                          dw, dwmol, dwp, dwt, ierr)

  ! temperature is in C and pressure in MPa
  ! kg/m^3
  ! Eq. 27b
  dw = dw + &
       s*(0.668d0 + 0.44d0*s + &
          1.d-6*(300.d0*p_MPa - 2400.d0*p_MPa*s + &
                 t_C*(80.d0 + 3.d0*t_C - 3300.d0*s - 13.d0*p_MPa + &
                      47.d0*p_Mpa*s))) * &
       g_cm3_to_kg_m3

  ! molar density H2O = solution density (kg/m3) * (1-mass frac salt) / (kg/kmol water)
  !dwmol = dw * (1-s) / (FMWH2O)
  dwmol = dw / FMWH2O

  if (calculate_derivatives) then
        ! v - this dwp is in the correct units of kmol/m^3-Pa
    dwp = dwp + &
          s*(1.d-6*(300.d0 - 2400.d0*s + t_C*(-13.d0 + 47.d0*s))) * &
                                 ! v - convert from kmol/m^3-MPa to kmol/m^3-Pa
          g_cm3_to_kg_m3/FMWH2O*Pa_to_MPa
    dwt = dwt + &
          s*(1.d-6*(80.d0 + 6.d0*t_C - 3300.d0*s - 13.d0*p_MPa + &
                    47.d0*p_Mpa*s)) * &
          g_cm3_to_kg_m3/FMWH2O
  else
    dwp = UNINITIALIZED_DOUBLE
    dwt = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterDensityBatzleAndWangExt

! ************************************************************************** !

subroutine EOSWaterViscosityBatzleAndWang(T, P, PS, dPS_dT, &
                                          calculate_derivatives, VW, &
                                          dVW_dT, dVW_dP, ierr,table_idxs)
  !
  ! From Batlze M. and Z. Wang (1992) Seismic properties of fluids, Geophysics,
  ! Vol. 57, No. 11, Pg. 1396-1408.
  !
  ! Equation 32
  !
  ! Author: Glenn Hammond
  ! Date: 02/08/16
  !
  implicit none
  PetscReal, intent(in) :: T       ! C
  PetscReal, intent(in) :: P       ! Pa
  PetscReal, intent(in) :: PS      ! Pa
  PetscReal, intent(in) :: dPS_dT  ! Pa/C
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: VW     ! Pa-s
  PetscReal, intent(out) :: dVW_dT, dVW_dP
  PetscErrorCode, intent(inout) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  ! convert from centipoise to Pa-s (1 cP = 1.d-3 Pa-s)
  PetscReal, parameter :: centipoise_to_Pa_s = 1.d-3
  PetscReal :: t_C
  PetscReal :: exponential_term
  PetscReal :: temperature_term

  t_C = T

  ! this is Eq. 32 without all the salt terms.
  ! -0.057138d0 = -1.d0*0.42d0*(-0.17d0)**2.d0+0.045d0
  exponential_term = -0.057138d0*t_C**0.8d0
  temperature_term = 1.65d0*exp(exponential_term)
  VW = 0.1d0 + temperature_term
  VW = VW * centipoise_to_Pa_s

  if (calculate_derivatives) then
    dVW_dP = 0.d0
    dVW_dT = 0.8d0*temperature_term*exponential_term/t_C*centipoise_to_Pa_s
  else
    dVW_dP = UNINITIALIZED_DOUBLE
    dVW_dT = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterViscosityBatzleAndWang

! ************************************************************************** !

subroutine EOSWaterViscosityBatzleAndWangExt(T, P, PS, dPS_dT, aux, &
                                             calculate_derivatives, VW, &
                                             dVW_dT, dVW_dP, ierr)
  !
  ! From Batlze M. and Z. Wang (1992) Seismic properties of fluids, Geophysics,
  ! Vol. 57, No. 11, Pg. 1396-1408.
  !
  ! Equation 32
  !
  ! Author: Glenn Hammond
  ! Date: 02/08/16
  !
  implicit none

  PetscReal, intent(in) :: T, P, PS, dPS_dT
  PetscReal, intent(in) :: aux(*)
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: VW
  PetscReal, intent(out) :: dVW_dT, dVW_dP
  PetscErrorCode, intent(inout) :: ierr

  ! convert from centipoise to Pa-s (1 cP = 1.d-3 Pa-s)
  PetscReal, parameter :: centipoise_to_Pa_s = 1.d-3
  PetscReal :: t_C
  PetscReal :: s
  PetscReal :: exponential_term
  PetscReal :: temperature_term

  s = aux(1)
  t_C = T

  exponential_term = -1.d0*(0.42d0*(s**0.8d0-0.17d0)**2.d0 + 0.045d0)* &
                     t_C**0.8d0
  temperature_term = (1.65d0 + 91.9d0*s**3.d0)*exp(exponential_term)
  VW = 0.1d0 + 0.333d0*s + temperature_term
  VW = VW * centipoise_to_Pa_s

  if (calculate_derivatives) then
    dVW_dP = 0.d0
    dVW_dT = 0.8d0*temperature_term*exponential_term/t_C*centipoise_to_Pa_s
  else
    dVW_dP = UNINITIALIZED_DOUBLE
    dVW_dT = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterViscosityBatzleAndWangExt

! ************************************************************************** !

subroutine EOSWaterSatPressSparrowExt(T,aux,calculate_derivatives, &
                                   PS, dPS_dT, ierr)
  !
  ! Calculates water saturation pressure based on Sparrow, 2003.
  ! Valid from 0<T<300 C
  !
  ! Author: David Fukuyama
  ! 06/01/21
  !
  implicit none

  PetscReal, intent(in) :: T  ! Temperature
  PetscReal, intent(in) :: aux(*) ! NaCl mole fraction
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: PS, dPS_dT ! Saturation pres. and derivative
  PetscErrorCode, intent(inout) :: ierr

  PetscReal, parameter :: MPa_to_Pa = 1.d6
  PetscReal :: A,B,C,D,E,s

  s = aux(1)
  if (T <= 150.d0) then ! eq 6
    A = (0.9083d0  +s*(-0.569d0 +s*(0.1945d0  +s*(-3.736d0 +s*2.82d0))))*1.d-3
    B = (-0.0669d0 +s*(0.0582d0 +s*(-0.1668d0 +s*(0.6761d0 +s*(-2.091d0)))))*1.d-3
    C = (7.541d0   +s*(-5.143d0 +s*(6.482d0   +s*(-52.62d0 +s*115.7d0))))*1.d-6
    D = (-0.0922d0 +s*(0.0649d0 +s*(-0.1313d0 +s*(0.8024d0 +s*(-1.986d0)))))*1.d-6
    E = (1.237d0   +s*(-0.753d0 +s*(0.1448d0  +s*(-6.964d0 +s*14.61d0))))*1.d-9
  elseif (T > 150.d0) then
    A = -3.248d0   +s*(7.081d0   +s*(-49.93d0 +s*(219.6d0  +s*(-308.5d0))))
    B = 0.0610d0   +s*(-0.1185d0 +s*(0.7916d0 +s*(-3.474d0 +s*4.882d0)))
    C = (-0.4109d0 +s*(0.6789d0  +s*(-4.155d0 +s*(18.34d0  +s*(-25.89d0)))))*1.d-3
    D = (1.13d0    +s*(-1.432d0  +s*(7.169d0  +s*(-33.17d0 +s*47.45d0))))*1.d-6
    E = 0.d0
  endif

  PS = A+T*(B+T*(C+T*(D+T*E)))
  PS = PS*MPa_to_Pa

  if (calculate_derivatives) then
    dPS_dT = B+T*(C+T*(D+T*E)) * MPa_to_Pa
  endif

end subroutine EOSWaterSatPressSparrowExt

! ************************************************************************** !

subroutine EOSWaterDensitySparrowExt(T,P, aux, &
                                  calculate_derivatives, &
                                  dw, dwmol, dwp, dwt, ierr)
  !
  ! Brine density calculation from Sparrow (2003)
  ! Values valid from 0<T<300 C
  !
  ! Author: David Fukuyama
  ! Date: 05/31/21
  !

  implicit none

  PetscReal, intent(in) :: T ! Celsius
  PetscReal, intent(in) :: P ! Pa
  PetscReal, intent(in) :: aux(*)
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(inout) :: ierr

  PetscReal :: A, B, C, D, E, s

  s = aux(1)
  ! eq 7 (0 C <= T <= 300 C)
  A = (1.001d0   +s*(0.7666d0 +s*(-0.0149d0 +s*(0.2663d0 +s*0.8845d0))))*1.d3
  B = -0.0214d0  +s*(-3.496d0 +s*(10.02d0   +s*(-6.56d0  +s*(-31.37d0))))
  C = (-5.263d0  +s*(39.87d0  +s*(-176.2d0  +s*(363.5d0  +s*(-7.784d0)))))*1.d-3
  D = (15.42d0   +s*(-167.d0  +s*(980.7d0   +s*(-2573.d0 +s*876.6d0))))*1.d-6
  E = (-0.0276d0 +s*(0.2978d0 +s*(-2.017d0  +s*(6.345d0  +s*(-3.914d0)))))*1.d-6

  dw = A+T*(B+T*(C+T*(D+E*T))) !kg/m^3
  !dwmol = dw*(1-s)/FMWH2O ! kmol/m^3
  dwmol = dw/FMWH2O

end subroutine EOSWaterDensitySparrowExt

! ************************************************************************** !

subroutine EOSWaterSolubility(T, solubility)
  !
  ! Determines solubility of NaCl in water based on Sparrow, 2003
  ! https://doi.org/10.1016/S0011-9164(03)90068-3
  ! returns solubility (mass fraction [g NaCl/g H2O])
  !
  ! Author: David Fukuyama
  ! Date: 12/17/21
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(out) :: solubility

  ! 0 C <= T <= 450 C
  solubility = 0.2628d0 + 62.75d-6 * T  + 1.084d-6 * T**2 ! eq 5

end subroutine EOSWaterSolubility

! ************************************************************************** !

subroutine EOSWaterComputeSalinity(T, x)
  !
  ! Computes salinity of NaCl in water based on EOS block inputs
  !
  ! Author: David Fukuyama
  ! Date: 4/6/22
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(out) :: x

  if (.not. halite_saturated_brine) then
     ! salinity [mass %]
    if (Initialized(constant_salinity)) then
      x = constant_salinity
    else
      x = 0.d0
    endif
  else !compute solubility
    call EOSWaterSolubility(T,x)
  endif

end subroutine EOSWaterComputeSalinity

! ************************************************************************** !
  subroutine TestEOSWaterBatzleAndWang()
  !
  ! From Batlze M. and Z. Wang (1992) Seismic properties of fluids, Geophysics,
  ! Vol. 57, No. 11, Pg. 1396-1408.
  !
  ! Code for testing
  !
  ! Author: Glenn Hammond
  ! Date: 02/15/16
  !
  PetscReal :: p, t, dw, dwmol, dwp, dwt, ps, dps_dt, vw, dvw_dt, dvw_dp
  PetscReal :: t2, dw2, p2, vw2
  PetscReal :: aux(1)
  PetscErrorCode :: ierr

  p = 0.101325d6
  t = 20.d0
  ps = -999.d0
  dps_dt = -999.d0
  aux(1) = 0.14d0

  call EOSWaterSetDensity('BATZLE_AND_WANG')
  call EOSWaterSetViscosity('BATZLE_AND_WANG')

  call EOSWaterDensityBatzleAndWang(t,p, PETSC_TRUE, &
                                    dw, dwmol, dwp, dwt, ierr)
  print *, 'dw:    ', dw
  print *, 'dwmol: ', dwmol
  print *, 'dwp:   ', dwp
  print *, 'dwt:   ', dwt

  t2 = t + 1.d-6*t
  call EOSWaterDensityBatzleAndWang(t2,p, PETSC_TRUE, &
                                    dw2, dwmol, dwp, dwt, ierr)
  print *
  print *, 'dw2:   ', dw2
  print *, 'dw(t): ', (dw2-dw)/(t2-t)

  p2 = p + 1.d-6*p
  call EOSWaterDensityBatzleAndWang(t,p2, PETSC_TRUE, &
                                    dw2, dwmol, dwp, dwt, ierr)
  print *
  print *, 'dw2:   ', dw2
  print *, 'dw(p): ', (dw2-dw)/(p2-p)


  call EOSWaterDensityBatzleAndWangExt(t,p,aux, PETSC_TRUE, &
                                    dw, dwmol, dwp, dwt, ierr)
  print *, 'Density-Ext'
  print *, 'dw:    ', dw
  print *, 'dwmol: ', dwmol
  print *, 'dwp:   ', dwp
  print *, 'dwt:   ', dwt

  call EOSWaterViscosityBatzleAndWang(t, p, PS, dPS_dT, &
                                      PETSC_TRUE, vw, &
                                      dvw_dt, dvw_dp, ierr)
  print *
  print *, 'vw:      ', vw
  print *, 'dvw_dp:  ', dvw_dp
  print *, 'dvw_dt:  ', dvw_dt


  call EOSWaterViscosityBatzleAndWangExt(t, p, PS, dPS_dT, aux, &
                                         PETSC_TRUE, vw, &
                                         dvw_dt, dvw_dp, ierr)
  print *, 'Ext-'
  print *, 'vw:      ', vw
  print *, 'dvw_dp:  ', dvw_dp
  print *, 'dvw(t)t:  ', dvw_dt

  call EOSWaterViscosityBatzleAndWangExt(t2, p, PS, dPS_dT, aux, &
                                         PETSC_TRUE, vw2, &
                                         dvw_dt, dvw_dp, ierr)

  print *, 'Ext-numerical'
  print *, 'vw:      ', vw2
  print *, 'dvw(t)t:  ', (vw2-vw)/(t2-t)

  call EOSWaterViscosityBatzleAndWangExt(t, p, PS, dPS_dT, aux, &
                                         PETSC_TRUE, vw, &
                                         dvw_dt, dvw_dp, ierr)
  print *, 'Ext-S'
  print *, 'vw:      ', vw
  print *, 'dvw_dp:  ', dvw_dp
  print *, 'dvw(t)t:  ', dvw_dt

end subroutine TestEOSWaterBatzleAndWang

! ************************************************************************** !

subroutine EOSWaterDensityExtNumericalDerive(t,p,aux,dw,dwmol,dwp,dwt,ierr)

  implicit none

  PetscReal, intent(in) :: t     ! Temperature in centigrade
  PetscReal, intent(in) :: p     ! Pressure in Pascal
  PetscReal, intent(in) :: aux(*)
  PetscReal, intent(out) :: dw ! kg/m^3
  PetscReal, intent(out) :: dwmol ! kmol/m^3
  PetscReal, intent(out) :: dwp ! kmol/m^3-Pa
  PetscReal, intent(out) :: dwt ! kmol/m^3-C
  PetscErrorCode, intent(inout) :: ierr

  PetscReal :: dwp_analytical, dwt_analytical
  PetscReal :: dwp_numerical, dwt_numerical
  PetscReal :: dum1, dum2, dum3
  PetscReal :: dwmol_tpert, dwmol_ppert
  PetscReal :: tpert, ppert, t_plus_tpert, p_plus_ppert
  PetscReal :: salinity(1)
  PetscReal :: pert_tol = 1.d-6

  tpert = t*pert_tol
  ppert = p*pert_tol
  t_plus_tpert = t + tpert
  p_plus_ppert = p + ppert
  salinity(1) = aux(1)

#if 0
  ! test against non-extended version
  call EOSWaterDensityPtr(t,p,PETSC_TRUE,dw,dwmol,dwp_analytical, &
                          dwt_analytical,ierr)
  dwp = dwp_analytical
  dwt = dwt_analytical
#else
  call EOSWaterDensityExtPtr(t,p,salinity,PETSC_TRUE, &
                             dw,dwmol,dwp_analytical,dwt_analytical,ierr)
  dwp = dwp_analytical
  dwt = dwt_analytical
  call EOSWaterDensityExtPtr(t_plus_tpert,p,salinity,PETSC_FALSE, &
                             dum1,dwmol_tpert,dum2,dum3,ierr)
  call EOSWaterDensityExtPtr(t,p_plus_ppert,salinity,PETSC_FALSE, &
                             dum1,dwmol_ppert,dum2,dum3,ierr)

  dwp_numerical = (dwmol_ppert-dwmol)/ppert
  dwt_numerical = (dwmol_tpert-dwmol)/tpert

  if (.not.PETSC_TRUE) then
    dwp = dwp_numerical
    dwt = dwt_numerical
  else
    dwp = dwp_analytical
    dwt = dwt_analytical
  endif

  if (dabs((dwp_numerical-dwp_analytical)/dwp_numerical) > 1.d-4) then
    print *, p, t, salinity(1), dw, dwmol, dwp, dwp_analytical, dwp_numerical
  endif
#endif

end subroutine EOSWaterDensityExtNumericalDerive

! **************************************************************************** !

subroutine EOSWaterInputRecord()
  !
  ! Prints ingested equation of state information to the input record file.
  !
  ! Author: Jenn Frederick
  ! Date: 05/04/2016
  !

  implicit none

  character(len=MAXWORDLENGTH) :: word1
  PetscInt :: id = INPUT_RECORD_UNIT

  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'WATER'

  ! water density [kg/m^3]
  if (associated(EOSWaterDensityPtr,EOSWaterDensityConstant)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(word1,*) constant_density
    write(id,'(a)') 'constant, ' // adjustl(trim(word1)) // ' kg/m^3'
  endif
  if (associated(EOSWaterDensityPtr,EOSWaterDensityExpPressure)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(id,'(a)') 'exponential'
    write(id,'(a29)',advance='no') 'exp. ref. density: '
    write(word1,*) exp_p_reference_density
    write(id,'(a)') adjustl(trim(word1)) // ' kg/m^3'
    write(id,'(a29)',advance='no') 'exp. ref. pressure: '
    write(word1,*) exp_p_reference_pressure
    write(id,'(a)') adjustl(trim(word1)) // ' Pa'
    write(id,'(a29)',advance='no') 'exp. water compressibility: '
    write(word1,*) exp_p_water_compressibility
    write(id,'(a)') adjustl(trim(word1)) // ' 1/Pa'
  endif
  if (associated(EOSWaterDensityPtr,EOSWaterDensityExpPressureTemp)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(id,'(a)') 'exponential temperature'
    write(id,'(a29)',advance='no') 'temp. ref. density: '
    write(word1,*) exp_pt_reference_density
    write(id,'(a)') adjustl(trim(word1)) // ' kg/m^3'
    write(id,'(a29)',advance='no') 'temp. ref. pressure: '
    write(word1,*) exp_pt_reference_pressure
    write(id,'(a)') adjustl(trim(word1)) // ' Pa'
    write(id,'(a29)',advance='no') 'temp. ref. pressure: '
    write(word1,*) exp_pt_reference_temperature
    write(id,'(a)') adjustl(trim(word1)) // ' C'
    write(id,'(a29)',advance='no') 'temp. water compressibility: '
    write(word1,*) exp_pt_water_compressibility
    write(id,'(a)') adjustl(trim(word1)) // ' 1/Pa'
    write(id,'(a29)',advance='no') 'temp. thermal expansion: '
    write(word1,*) exp_pt_thermal_expansion
    write(id,'(a)') adjustl(trim(word1)) // ' 1/K'
  endif
  if (associated(EOSWaterDensityPtr,EOSWaterDensityLinear)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(id,'(a)') 'linear'
    write(id,'(a29)',advance='no') 'linear ref. density: '
    write(word1,*) linear_reference_density
    write(id,'(a)') adjustl(trim(word1)) // ' kg/m^3'
    write(id,'(a29)',advance='no') 'linear ref. pressure: '
    write(word1,*) linear_reference_pressure
    write(id,'(a)') adjustl(trim(word1)) // ' Pa'
    write(id,'(a29)',advance='no') 'linear water compressibility: '
    write(word1,*) linear_water_compressibility
    write(id,'(a)') adjustl(trim(word1)) // ' 1/Pa'
  endif
  if (associated(EOSWaterDensityPtr,EOSWaterDensityBRAGFLO)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(id,'(a)') 'BRAGFLO'
    write(id,'(a29)',advance='no') 'exp. ref. density: '
    write(word1,*) exp_p_reference_density
    write(id,'(a)') adjustl(trim(word1)) // ' kg/m^3'
    write(id,'(a29)',advance='no') 'exp. ref. pressure: '
    write(word1,*) exp_p_reference_pressure
    write(id,'(a)') adjustl(trim(word1)) // ' Pa'
    write(id,'(a29)',advance='no') 'exp. water compressibility: '
    write(word1,*) exp_p_water_compressibility
    write(id,'(a)') adjustl(trim(word1)) // ' 1/Pa'
  endif
  if (associated(EOSWaterDensityPtr,EOSWaterDensityIFC67)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(id,'(a)') 'default, IFC67'
  endif
  if (associated(EOSWaterDensityPtr,EOSWaterDensityTGDPB01)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(id,'(a)') 'TGDPB01'
  endif
  if (associated(EOSWaterDensityPtr,EOSWaterDensityPainter)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(id,'(a)') 'PAINTER'
  endif
  if (associated(EOSWaterDensityPtr,EOSWaterDensityBatzleAndWang) .and. &
      associated(EOSWaterDensityExtPtr,EOSWaterDensityBatzleAndWangExt)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(id,'(a)') 'Batzle and Wang'
  endif
  if (associated(EOSWaterDensityPtr,EOSWaterDensityQuadratic)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(id,'(a)') 'Quadratic'
    write(id,'(a29)',advance='no') 'quad. ref. density: '
    write(word1,*) quadratic_reference_density
    write(id,'(a)') adjustl(trim(word1)) // ' kg/m^3'
    write(id,'(a29)',advance='no') 'quad. ref. pressure: '
    write(word1,*) quadratic_reference_pressure
    write(id,'(a)') adjustl(trim(word1)) // ' Pa'
    write(id,'(a29)',advance='no') 'quad. water compressibility: '
    write(word1,*) quadratic_wat_compressibility
    write(id,'(a)') adjustl(trim(word1)) // ' 1/Pa'
  end if
  if (associated(EOSWaterDensityPtr,EOSWaterDensityTrangenstein)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(id,'(a)') 'Trangenstein'
  end if

  ! water viscosity [Pa-s]
  if (associated(EOSWaterViscosityPtr,EOSWaterViscosityConstant)) then
    write(id,'(a29)',advance='no') 'water viscosity: '
    write(word1,*) constant_viscosity
    write(id,'(a)') 'constant, ' // trim(word1) // ' Pa-sec'
  endif
  if (associated(EOSWaterViscosityPtr,EOSWaterViscosity1)) then
    write(id,'(a29)',advance='no') 'water viscosity: '
    write(id,'(a)') 'default'
  endif
  if (associated(EOSWaterViscosityPtr,EOSWaterViscosityBatzleAndWang)) then
    write(id,'(a29)',advance='no') 'water viscosity: '
    write(id,'(a)') 'Batzle and Wang'
  endif
  if (associated(EOSWaterViscosityPtr,EOSWaterViscosityGrabowski)) then
    write(id,'(a29)',advance='no') 'water viscosity: '
    write(id,'(a)') 'Grabowski Batzle'
  endif

  ! water enthalpy [J/kmol]
  if (associated(EOSWaterViscosityPtr,EOSWaterEnthalpyConstant)) then
    write(id,'(a29)',advance='no') 'water enthalpy: '
    write(word1,*) constant_enthalpy
    write(id,'(a)') 'constant, ' // trim(word1) // ' J/kmol'
  endif
  if (associated(EOSWaterViscosityPtr,EOSWaterEnthalpyIFC67)) then
    write(id,'(a29)',advance='no') 'water enthalpy: '
    write(id,'(a)') 'default, IFC67'
  endif
  if (associated(EOSWaterViscosityPtr,EOSWaterEnthalpyPainter)) then
    write(id,'(a29)',advance='no') 'water enthalpy: '
    write(id,'(a)') 'PAINTER'
  endif

  ! steam density
  if (associated(EOSWaterSteamDensityEnthalpyPtr, &
                 EOSWaterSteamDenEnthConstant)) then
    write(id,'(a29)',advance='no') 'steam density: '
    write(word1,*) constant_steam_density
    write(id,'(a)') 'constant, ' // trim(word1) // ' kg/m^3'
  endif
  if (associated(EOSWaterSteamDensityEnthalpyPtr, &
                 EOSWaterSteamDensityEnthalpyIFC67)) then
    write(id,'(a29)',advance='no') 'steam density: '
    write(id,'(a)') 'default, IFC67'
  endif

  ! steam enthalpy [J/kmol]
  if (associated(EOSWaterSteamDensityEnthalpyPtr, &
                 EOSWaterSteamDenEnthConstant)) then
    write(id,'(a29)',advance='no') 'steam enthalpy: '
    write(word1,*) constant_steam_enthalpy
    write(id,'(a)') 'constant, ' // trim(word1) // ' J/kmol'
  endif
  if (associated(EOSWaterSteamDensityEnthalpyPtr, &
                 EOSWaterSteamDensityEnthalpyIFC67)) then
    write(id,'(a29)',advance='no') 'steam enthalpy: '
    write(id,'(a)') 'default, IFC67'
  endif

  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'

end subroutine EOSWaterInputRecord

! ************************************************************************** !

subroutine EOSWaterTest(temp_low,temp_high,pres_low,pres_high, &
                        ntemp,npres,uniform_temp,uniform_pres,filename)

  use Utility_module, only : InitToNan

  implicit none

  PetscReal :: temp_low
  PetscReal :: temp_high
  PetscReal :: pres_low
  PetscReal :: pres_high
  PetscInt :: npres
  PetscInt :: ntemp
  PetscBool :: uniform_temp
  PetscBool :: uniform_pres
  character(len=MAXWORDLENGTH) :: filename

  PetscReal, allocatable :: temp(:)
  PetscReal, allocatable :: pres(:)
  PetscReal, allocatable :: density_kg(:,:)
  PetscReal, allocatable :: enthalpy(:,:)
  PetscReal, allocatable :: viscosity(:,:)
  PetscReal, allocatable :: saturation_pressure_array(:)
  PetscReal :: dum1, dum2, dum3
  PetscInt :: itemp, ipres
  PetscReal :: ln_low, ln_high
  PetscReal :: saturation_pressure
  PetscReal :: NaN
  character(len=MAXWORDLENGTH) :: eos_density_name
  character(len=MAXWORDLENGTH) :: eos_enthalpy_name
  character(len=MAXWORDLENGTH) :: eos_viscosity_name
  character(len=MAXWORDLENGTH) :: eos_saturation_pressure_name
  character(len=MAXSTRINGLENGTH) :: header, string

  PetscErrorCode :: ierr

  NaN = InitToNan()

  allocate(temp(ntemp))
  temp = UNINITIALIZED_DOUBLE
  allocate(pres(ntemp))
  pres = UNINITIALIZED_DOUBLE
  allocate(density_kg(npres,ntemp))
  density_kg = UNINITIALIZED_DOUBLE
  allocate(viscosity(npres,ntemp))
  viscosity = UNINITIALIZED_DOUBLE
  allocate(enthalpy(npres,ntemp))
  enthalpy = UNINITIALIZED_DOUBLE
  allocate(saturation_pressure_array(ntemp))
  saturation_pressure_array = UNINITIALIZED_DOUBLE

  if (uniform_pres) then
    do ipres = 1, npres
      pres(ipres) = &
        (pres_high-pres_low)/max(dble(npres-1),1.d0) * (ipres-1) + pres_low
    enddo
  else
    ln_high = log(pres_high)
    ln_low = log(pres_low)
    do ipres = 1, npres
      pres(ipres) = &
        exp((ln_high-ln_low)/max(dble(npres-1),1.d0) * (ipres-1) + ln_low)
    enddo
  endif

  if (uniform_temp) then
    do itemp = 1, ntemp
      temp(itemp) = &
        (temp_high-temp_low)/max(dble(ntemp-1),1.d0) * (itemp-1) + temp_low
    enddo
  else
    ln_high = log(temp_high)
    ln_low = log(temp_low)
    do itemp = 1, ntemp
      temp(itemp) = &
        exp((ln_high-ln_low)/max(dble(ntemp-1),1.d0) * (itemp-1) + ln_low)
    enddo
  endif

  ! density
  if (associated(EOSWaterDensityPtr,EOSWaterDensityConstant)) then
    eos_density_name = 'Constant'
  else if (associated(EOSWaterDensityPtr,EOSWaterDensityExpPressure)) then
    eos_density_name = 'Exponential Pressure'
  else if (associated(EOSWaterDensityPtr,EOSWaterDensityExpPressureTemp)) then
    eos_density_name = 'Exponential Pressure Temperature'
  else if (associated(EOSWaterDensityPtr,EOSWaterDensityLinear)) then
    eos_density_name = 'Linear'
  else if (associated(EOSWaterDensityPtr,EOSWaterDensityBRAGFLO)) then
    eos_density_name = 'BRAGFLO'
  else if (associated(EOSWaterDensityPtr,EOSWaterDensityIFC67)) then
    eos_density_name = 'IFC67'
  else if (associated(EOSWaterDensityPtr,EOSWaterDensityTGDPB01)) then
    eos_density_name = 'TGDPB01'
  else if (associated(EOSWaterDensityPtr,EOSWaterDensityPainter)) then
    eos_density_name = 'Painter'
  else if (associated(EOSWaterDensityPtr,EOSWaterDensityBatzleAndWang)) then
    eos_density_name = 'Batzle and Wang'
  else if (associated(EOSWaterDensityPtr,EOSWaterDensityQuadratic)) then
    eos_density_name = 'Quadratic'
  else if (associated(EOSWaterDensityPtr,EOSWaterDensityTrangenstein)) then
    eos_density_name = 'Trangenstein'
  else
    eos_density_name = 'Unknown'
  endif

  ! enthalpy
  if (associated(EOSWaterEnthalpyPtr,EOSWaterEnthalpyConstant)) then
    eos_enthalpy_name = 'Constant'
  else if (associated(EOSWaterEnthalpyPtr,EOSWaterEnthalpyIFC67)) then
    eos_enthalpy_name = 'IFC67'
  else if (associated(EOSWaterEnthalpyPtr,EOSWaterEnthalpyPainter)) then
    eos_enthalpy_name = 'Painter'
  else
    eos_enthalpy_name = 'Unknown'
  endif

  ! viscosity
  if (associated(EOSWaterViscosityPtr,EOSWaterViscosityConstant)) then
    eos_viscosity_name = 'Constant'
  else if (associated(EOSWaterViscosityPtr,EOSWaterViscosity1)) then
    eos_viscosity_name = 'Default'
  else if (associated(EOSWaterViscosityPtr,EOSWaterViscosityBatzleAndWang)) then
    eos_viscosity_name = 'Batzle and Wang'
  else if (associated(EOSWaterViscosityPtr,EOSWaterViscosityGrabowski)) then
    eos_viscosity_name = 'Grabowski'
  else
    eos_viscosity_name = 'Unknown'
  endif

  ! saturation pressure
  if (associated(EOSWaterSaturationPressurePtr, &
                 EOSWaterSaturationPressureIFC67)) then
    eos_saturation_pressure_name = 'IFC67'
  elseif (associated(EOSWaterSaturationPressurePtr, &
                     EOSWaterSatPresWagnerPruss)) then
    eos_saturation_pressure_name = 'WagnerAndPruss'
  elseif (associated(EOSWaterSaturationPressurePtr, &
                     EOSWaterSaturationPressureIce)) then
    eos_saturation_pressure_name = 'Huang-Ice'
  else
    eos_saturation_pressure_name = 'Unknown'
  endif

  do itemp = 1, ntemp
    do ipres = 1, npres
      ! EOSWaterSaturationPressurePtr() must come before call to
      ! EOSWaterViscosityPtr()
      call EOSWaterSaturationPressurePtr(temp(itemp),PETSC_FALSE, &
                                         saturation_pressure,dum1,ierr)
      if (ipres == 1) &
        saturation_pressure_array(itemp) = saturation_pressure
      call EOSWaterDensityPtr(temp(itemp),pres(ipres),PETSC_FALSE, &
                              density_kg(ipres,itemp), &
                              dum1,dum2,dum3,ierr)
      call EOSWaterEnthalpyPtr(temp(itemp),pres(ipres),PETSC_FALSE, &
                               enthalpy(ipres,itemp),dum1,dum2,ierr)
      call EOSWaterViscosityPtr(temp(itemp),pres(ipres),saturation_pressure, &
                                dum1,PETSC_FALSE,viscosity(ipres,itemp), &
                                dum2,dum3,ierr)
    enddo
  enddo

100 format(100es16.8)
  if (len_trim(filename) == 0) then
    string = 'eos_water_test.txt'
  else
    string = filename
  endif
  open(unit=IUNIT_TEMP,file=string)
  header = 'T[C], P[Pa], &
    &Density (' // trim(eos_density_name) // ') [kg/m^3], &
    &Enthalpy (' // trim(eos_enthalpy_name) // ') [J/kmol], &
    &Viscosity (' // trim(eos_viscosity_name) // ') [Pa-s], &
    &Saturation Pressure (' // trim(eos_saturation_pressure_name) // ') [Pa]'
  write(IUNIT_TEMP,'(a)') trim(header)
  write(IUNIT_TEMP,'(100i9)') ntemp, npres
  do itemp = 1, ntemp
    do ipres = 1, npres
      write(IUNIT_TEMP,100) temp(itemp), pres(ipres), &
            density_kg(ipres,itemp), enthalpy(ipres,itemp), &
            viscosity(ipres,itemp), saturation_pressure_array(itemp)
    enddo
  enddo
  close(IUNIT_TEMP)

  deallocate(temp)
  deallocate(pres)
  deallocate(density_kg)
  deallocate(enthalpy)
  deallocate(viscosity)
  deallocate(saturation_pressure_array)

end subroutine EOSWaterTest

! ************************************************************************** !

subroutine EOSWaterSteamTest(temp_low,temp_high,pres_low,pres_high, &
                             ntemp,npres,uniform_temp,uniform_pres,filename)

  use Utility_module, only : InitToNan

  implicit none

  PetscReal :: temp_low
  PetscReal :: temp_high
  PetscReal :: pres_low
  PetscReal :: pres_high
  PetscInt :: npres
  PetscInt :: ntemp
  PetscBool :: uniform_temp
  PetscBool :: uniform_pres
  character(len=MAXWORDLENGTH) :: filename

  PetscReal, allocatable :: temp(:)
  PetscReal, allocatable :: pres(:)
  PetscReal, allocatable :: density_kg(:,:)
  PetscReal, allocatable :: density_kg_psat(:)
  PetscReal, allocatable :: enthalpy(:,:)
  PetscReal, allocatable :: enthalpy_psat(:)
  PetscReal, allocatable :: viscosity(:,:)
  PetscReal, allocatable :: saturation_pressure_array(:)
  PetscReal :: dum1, dum2, dum3, dum4, dum5
  PetscInt :: itemp, ipres
  PetscReal :: ln_low, ln_high
  PetscReal :: saturation_pressure
  PetscReal :: NaN
  character(len=MAXWORDLENGTH) :: eos_density_name
  character(len=MAXWORDLENGTH) :: eos_enthalpy_name
  character(len=MAXWORDLENGTH) :: eos_saturation_pressure_name
  character(len=MAXSTRINGLENGTH) :: header, string

  PetscErrorCode :: ierr

  NaN = InitToNan()

  allocate(temp(ntemp))
  temp = UNINITIALIZED_DOUBLE
  allocate(pres(ntemp))
  pres = UNINITIALIZED_DOUBLE
  allocate(density_kg(npres,ntemp))
  density_kg = UNINITIALIZED_DOUBLE
  allocate(density_kg_psat(ntemp))
  density_kg_psat = UNINITIALIZED_DOUBLE
  allocate(viscosity(npres,ntemp))
  viscosity = UNINITIALIZED_DOUBLE
  allocate(enthalpy(npres,ntemp))
  enthalpy = UNINITIALIZED_DOUBLE
  allocate(enthalpy_psat(ntemp))
  enthalpy_psat = UNINITIALIZED_DOUBLE
  allocate(saturation_pressure_array(ntemp))
  saturation_pressure_array = UNINITIALIZED_DOUBLE

  if (uniform_pres) then
    do ipres = 1, npres
      pres(ipres) = (pres_high-pres_low)/dble(npres-1) * (ipres-1) + pres_low
    enddo
  else
    ln_high = log(pres_high)
    ln_low = log(pres_low)
    do ipres = 1, npres
      pres(ipres) = exp((ln_high-ln_low)/dble(npres-1) * (ipres-1) + ln_low)
    enddo
  endif

  if (uniform_temp) then
    do itemp = 1, ntemp
      temp(itemp) = (temp_high-temp_low)/dble(ntemp-1) * (itemp-1) + temp_low
    enddo
  else
    ln_high = log(temp_high)
    ln_low = log(temp_low)
    do itemp = 1, ntemp
      temp(itemp) = exp((ln_high-ln_low)/dble(ntemp-1) * (itemp-1) + ln_low)
    enddo
  endif

  ! density
  if (associated(EOSWaterSteamDensityEnthalpyPtr,EOSWaterSteamDenEnthConstant)) then
    eos_density_name = 'Constant'
  else if (associated(EOSWaterSteamDensityEnthalpyPtr, &
                      EOSWaterSteamDensityEnthalpyIFC67)) then
    eos_density_name = 'IFC67'
  else
    eos_density_name = 'Unknown'
  endif

  ! enthalpy
  if (associated(EOSWaterSteamDensityEnthalpyPtr,EOSWaterSteamDenEnthConstant)) then
    eos_enthalpy_name = 'Constant'
  else if (associated(EOSWaterSteamDensityEnthalpyPtr, &
                      EOSWaterSteamDensityEnthalpyIFC67)) then
    eos_enthalpy_name = 'IFC67'
  else
    eos_enthalpy_name = 'Unknown'
  endif

  ! saturation pressure
  if (associated(EOSWaterSaturationPressurePtr, &
                 EOSWaterSaturationPressureIFC67)) then
    eos_saturation_pressure_name = 'IFC67'
  elseif (associated(EOSWaterSaturationPressurePtr, &
                     EOSWaterSatPresWagnerPruss)) then
    eos_saturation_pressure_name = 'WagnerAndPruss'
  else
    eos_saturation_pressure_name = 'Unknown'
  endif

  do itemp = 1, ntemp
    do ipres = 1, npres
      ! EOSWaterSaturationPressurePtr() must come before call to
      ! EOSWaterViscosityPtr()
      call EOSWaterSaturationPressurePtr(temp(itemp),PETSC_FALSE, &
                                         saturation_pressure,dum1,ierr)
      if (ipres == 1) then
        saturation_pressure_array(itemp) = saturation_pressure
        call EOSWaterSteamDensityEnthalpyPtr(temp(itemp),saturation_pressure, &
                                             PETSC_FALSE, &
                                             density_kg_psat(itemp),dum1, &
                                             enthalpy_psat(itemp), &
                                             dum2,dum3,dum4,dum5,ierr)
      endif
      call EOSWaterSteamDensityEnthalpyPtr(temp(itemp),pres(ipres),PETSC_FALSE, &
                                           density_kg(ipres,itemp),dum1, &
                                           enthalpy(ipres,itemp), &
                                           dum2,dum3,dum4,dum5,ierr)
    enddo
  enddo

100 format(100es16.8)
  if (len_trim(filename) == 0) then
    string = 'eos_water_steam_test.txt'
  else
    string = filename
  endif
  open(unit=IUNIT_TEMP,file=string)
  header = 'T[C], P[Pa], &
    &Steam Density (' // trim(eos_density_name) // ') [kg/m^3], &
    &Steam Enthalpy (' // trim(eos_enthalpy_name) // ') [J/kmol]'
  write(IUNIT_TEMP,'(a)') trim(header)
  write(IUNIT_TEMP,'(100i9)') ntemp, npres
  do itemp = 1, ntemp
    do ipres = 1, npres
      write(IUNIT_TEMP,100) temp(itemp), pres(ipres), &
            density_kg(ipres,itemp), enthalpy(ipres,itemp)
    enddo
  enddo
  close(IUNIT_TEMP)

  if (len_trim(filename) == 0) then
    string = 'eos_water_steam_psat_test.txt'
  else
    itemp = index(filename,'.')
    string = filename
    string(itemp:itemp+4) = '_psat'
    string(itemp+5:) = filename(itemp:)
  endif
  open(unit=IUNIT_TEMP,file=string)
  header = 'T[C], &
    &Water Saturation Pressure (' // trim(eos_saturation_pressure_name) // ') [Pa], &
    &Steam Density @ Saturation Pressure (' // trim(eos_density_name) // ') [kg/m^3], &
    &Steam Enthalpy @ Saturation Pressure (' // trim(eos_enthalpy_name) // ') [J/kmol]'
  write(IUNIT_TEMP,'(a)') trim(header)
  write(IUNIT_TEMP,'(100i9)') ntemp
  do itemp = 1, ntemp
    write(IUNIT_TEMP,100) temp(itemp), &
          saturation_pressure_array(itemp), &
          density_kg_psat(itemp), enthalpy_psat(itemp)
  enddo
  close(IUNIT_TEMP)

  deallocate(temp)
  deallocate(pres)
  deallocate(density_kg)
  deallocate(enthalpy)
  deallocate(density_kg_psat)
  deallocate(enthalpy_psat)
  deallocate(saturation_pressure_array)

end subroutine EOSWaterSteamTest

! ************************************************************************** !

subroutine EOSWaterSurfaceTension(T,sigma)
  !
  ! Surface tension of water equation from Revised Release on Surface
  ! Tension of Ordinary Water Substance, June 2014. Valid from -25C to
  ! 373 C
  !
  ! Author: Michael Nole
  !

  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(out) :: sigma

  PetscReal, parameter :: Tc = 647.096d0
  PetscReal, parameter :: B = 235.8d0
  PetscReal, parameter :: b_2 = -0.625d0
  PetscReal, parameter :: mu = 1.256d0
  PetscReal, parameter :: sigma_base = 0.073d0
  PetscReal :: Temp
  PetscReal :: tao

  Temp=T+273.15d0

  if (T <= 373.d0) then
    tao = 1.d0-Temp/Tc
    sigma = B*(tao**mu)*(1+b_2*tao)
    sigma = sigma * 1.d-3
  else
    sigma = 0.d0
  endif
  sigma= sigma/sigma_base

  !TOUGH3 way (not pressure-dependent)
  !if (Temp >= 101) sigma = 0

end subroutine EOSWaterSurfaceTension



! ************************************************************************** !

subroutine EOSWaterKelvin(Pc,rhow,T,Psat,Pv)
  !
  ! Kelvin equation for vapor pressure lowering
  ! vis-a-vis the Young-Laplace equation where
  ! Kelvin Eqn: ln(P/P_sat) = 2*sigma*Vm/(rRT)
  ! Y-L Eqn: Pc = 2*sigma*cos(theta)/r
  ! cos(theta) -> 1
  !
  ! Author: Michael Nole
  ! Date: 02/04/22
  !

  use PFLOTRAN_Constants_module

  implicit none

  PetscReal, intent(in) :: Pc, T, Psat, rhow
  PetscReal, intent(out) :: Pv

  PetscReal :: vp_factor, T_temp

  T_temp = T + 273.15d0

  ! For water:
  vp_factor = Pc / (rhow * 1000.d0 * IDEAL_GAS_CONSTANT * T_temp)

  Pv = exp(vp_factor) * Psat

end subroutine EOSWaterKelvin

! ************************************************************************** !

end module EOS_Water_module
