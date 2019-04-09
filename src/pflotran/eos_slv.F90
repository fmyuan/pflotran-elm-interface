module EOS_Slv_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use EOSData_module

  implicit none

  private

  ! module variables
  PetscReal :: fmw_slv           !kg/Kmol
  PetscReal :: constant_density
  PetscReal :: constant_viscosity
  PetscReal :: reference_density_kg

  PetscReal :: exponent_reference_density
  PetscReal :: exponent_reference_pressure
  PetscReal :: exponent_slv_compressibility

  ! This is the offset added to temperature [C] used to calculate the energy
  ! equation of state.  Switching between 0. and 273.15 greatly changes results.
#if defined(MATCH_TOUGH2)
  PetscReal, parameter :: T_energy_offset = 0.d0
#else
  PetscReal, parameter :: T_energy_offset = 273.15d0
#endif

  ! EOS databases
  class(eos_database_type), pointer :: eos_dbase => null()

! PVT tables - eos_tables
  class(eos_table_type), pointer :: pvt_table => null()

! Set null initial values for proceedure pointers of three types

  procedure(EOSSlvViscosityDummy),pointer::EOSSlvViscosityPtr=>null()
  procedure(EOSSlvDensityDummy)  ,pointer::EOSSlvDensityPtr  =>null()
  procedure(EOSSlvEnergyDummy)   ,pointer::EOSSlvEnergyPtr   =>null()

! interface blocks
  interface
    subroutine EOSSlvViscosityDummy(T, P, V_mix, &
                                    calculate_derivative, dV_dT, &
                                    dV_dP, ierr, table_idxs)
      implicit none
      PetscReal, intent(in) :: T            ! temperature [C]
      PetscReal, intent(in) :: P            ! pressure    [Pa]
      PetscReal, intent(out) :: V_mix       ! mixture viscosity
      PetscBool, intent(in) :: calculate_derivative ! Request derivatives
      PetscReal, intent(out) :: dV_dT       ! derivative wrt temperature
      PetscReal, intent(out) :: dV_dP       ! derivative wrt slv pressure
      PetscErrorCode, intent(out) :: ierr
      PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
    end subroutine EOSSlvViscosityDummy

    subroutine EOSSlvDensityDummy(T,P,Rho_slv,dRho_dT,dRho_dP,ierr,table_idxs)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: P        ! pressure [Pa]
      PetscReal, intent(out) :: Rho_slv ! slv density [kmol/m^3]
      PetscReal, intent(out) :: dRho_dT ! derivative slv density wrt temperature
      PetscReal, intent(out) :: dRho_dP ! derivative slv density wrt pressuret
      PetscErrorCode, intent(out) :: ierr
      PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
    end subroutine EOSSlvDensityDummy

    subroutine EOSSlvEnergyDummy(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
      implicit none
      PetscReal, intent(in)  :: T       ! temperature [C]
      PetscReal, intent(in)  :: P       ! pressure [Pa]
      PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
      PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
      PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
      PetscReal, intent(out) :: U       ! internal energy [J/kmol]
      PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
      PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSSlvEnergyDummy

  end interface

! Overloads for derivative/non-deriv versions called from outside module.

  interface EOSSlvViscosity
    procedure EOSSlvViscosityNoDerive
    procedure EOSSlvViscosityDerive
  end interface
  interface EOSSlvDensity
    procedure EOSSlvDensityNoDerive
    procedure EOSSlvDensityDerive
  end interface
  interface EOSSlvEnergy
    procedure EOSSlvEnergyNoDerive
    procedure EOSSlvEnergyDerive
  end interface
  interface EOSSlvDensityEnergy
    procedure EOSSlvDenEnthNoDerive
    procedure EOSSlvDenEnthDerive
  end interface

! The "public" definition that makes subroutines visible outside.

  public :: EOSSlvInit, &
            EOSSlvVerify, &
            EOSSlvViscosity, &
            EOSSlvDensity, &
            EOSSlvEnergy, &
            EOSSlvDensityEnergy, &
            EOSSlvInputRecord, &
            EOSSlvTest

  public :: EOSSlvSetDensityIdeal, &
            EOSSlvSetEnergyIdeal, &
            EOSSlvSetFMW, &
            EOSSlvGetFMW, &
            EOSSlvSetViscosityConstant, &
            EOSSlvSetReferenceDensity, &
            EOSSlvGetReferenceDensity, &
            EOSSlvSetPVDS, &
            EOSSlvSetEOSDBase, &
            EOSSlvTableProcess, &
            EOSSlvDBaseDestroy

  contains

! ************************************************************************** !

subroutine EOSSlvInit()

  implicit none

  fmw_slv = UNINITIALIZED_DOUBLE
  constant_density = UNINITIALIZED_DOUBLE
  constant_viscosity = UNINITIALIZED_DOUBLE

  reference_density_kg = UNINITIALIZED_DOUBLE

  exponent_reference_density = UNINITIALIZED_DOUBLE
  exponent_reference_pressure = UNINITIALIZED_DOUBLE
  exponent_slv_compressibility = UNINITIALIZED_DOUBLE

  EOSSlvDensityPtr  =>EOSSlvDensityTable
  EOSSlvEnergyPtr   =>EOSSlvEnergyIdeal
  EOSSlvViscosityPtr=>EOSSlvViscosityTable

end subroutine EOSSlvInit

! ************************************************************************** !

subroutine EOSSlvVerify(ierr,error_string)

  implicit none

  PetscErrorCode, intent(out) :: ierr
  character(len=MAXSTRINGLENGTH), intent(out) :: error_string

  ierr = 0

! Ideal density but constant density set
  error_string = ''
  if ( associated(EOSSlvDensityPtr,EOSSlvDensityIdeal) .and. &
       Initialized(constant_density)) then
    ierr = 1
  endif

! Constant density but constant density unset
  if (associated(EOSSlvDensityPtr,EOSSlvDensityConstant) .and. &
      Uninitialized(constant_density)) then
    error_string = trim(error_string)//' CONSTANT density not set.'
    ierr = 1
  endif

! Constant viscosity but constant viscosity not set
  if ( associated(EOSSlvViscosityPtr,EOSSlvViscosityConstant) .and. &
       Uninitialized(constant_viscosity) ) then
    ierr = 1
  endif

! Table density but table not set
  if(       associated(EOSSlvDensityPtr,EOSSlvDensityTable).and. &
     (.not.associated(pvt_table) ) ) then
    ierr=1
  endif

  if (Uninitialized(fmw_slv)) then
    fmw_slv = FMWCO2
  end if

end subroutine EOSSlvVerify

! ************************************************************************** !

subroutine EOSSlvSetDensityIdeal()

  implicit none

  EOSSlvDensityPtr=> EOSSlvDensityIdeal

end subroutine EOSSlvSetDensityIdeal

! ************************************************************************** !

subroutine EOSSlvSetEOSDBase(filename,option)

!------------------------------------------------------------------------------
! Request use of database for solvent characterisation
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : May 2018
!------------------------------------------------------------------------------

  use Option_module

  implicit none

  character(len=MAXWORDLENGTH) :: filename
  type(option_type) :: option

  eos_dbase => EOSDatabaseCreate(filename,'slv_database')
  call eos_dbase%Read(option)

! Set up property function pointers

  EOSSlvDensityPtr => EOSSlvDensityEOSDBase
  EOSSlvEnergyPtr => EOSSlvEnergyEOSDBase
  EOSSlvViscosityPtr => EOSSlvViscosityEOSDBase

end subroutine EOSSlvSetEOSDBase

! *****************************************************************************

subroutine EOSSlvEnergyEOSDBase(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)

!------------------------------------------------------------------------------
! Energy and Enthalpy look-up from table
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : May 2018
!------------------------------------------------------------------------------

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
  PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
  PetscErrorCode, intent(out) :: ierr

  ierr=0

  call eos_dbase%EOSPropGrad(T,P,EOS_ENTHALPY,H,dH_dT,dH_dP,ierr)
  call eos_dbase%EOSPropGrad(T,P,EOS_INTERNAL_ENERGY,U,dU_dT,dU_dP,ierr)

! J/kg * kg/Kmol = J/Kmol
  H = H  * fmw_slv
  dH_dT = dH_dT * fmw_slv
  dH_dP = dH_dP * fmw_slv

! J/kg * kg/Kmol = J/Kmol
  U = U  * fmw_slv
  dU_dT = dU_dT * fmw_slv
  dU_dP = dU_dP * fmw_slv

end subroutine EOSSlvEnergyEOSDBase

! *****************************************************************************

subroutine EOSSlvViscosityEOSDBase(T, P, V_mix, &
                                   calculate_derivative, dV_dT,dV_dP, &
                                   ierr, table_idxs)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! solvent pressure [Pa]
  PetscReal, intent(out) :: V_mix   ! mixture viscosity
  PetscBool, intent(in) :: calculate_derivative ! Request derivatives
  PetscReal, intent(out) :: dV_dT   ! derivative wrt temperature
  PetscReal, intent(out) :: dV_dP   ! derivative wrt gas pressure 
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  ierr=0

  call eos_dbase%EOSPropGrad(T,P,EOS_VISCOSITY,V_mix,dV_dT,dV_dP,ierr)

end subroutine EOSSlvViscosityEOSDBase

! ************************************************************************** !

subroutine EOSSlvSetEnergyIdeal()

  implicit none

  EOSSlvEnergyPtr=> EOSSlvEnergyIdeal

end subroutine EOSSlvSetEnergyIdeal

! ************************************************************************** !

subroutine EOSSlvSetFMW(input_fmw_slv)

  implicit none

  PetscReal,intent(in) :: input_fmw_slv

  fmw_slv = input_fmw_slv

end subroutine EOSSlvSetFMW

! ************************************************************************** !

function EOSSlvGetFMW()

  implicit none

  PetscReal :: EOSSlvGetFMW

  EOSSlvGetFMW = fmw_slv

end function EOSSlvGetFMW

! ************************************************************************** !

subroutine EOSSlvSetReferenceDensity(input_ref_density_kg)

  implicit none

  PetscReal,intent(in) :: input_ref_density_kg

  reference_density_kg = input_ref_density_kg

end subroutine EOSSlvSetReferenceDensity

! ************************************************************************** !

function EOSSlvGetReferenceDensity()

  implicit none

  PetscReal :: EOSSlvGetReferenceDensity

  EOSSlvGetReferenceDensity = reference_density_kg

end function EOSSlvGetReferenceDensity

! ************************************************************************** !

subroutine EOSSlvSetViscosityConstant(viscosity)

  implicit none

  PetscReal,intent(in) :: viscosity

  constant_viscosity = viscosity
  EOSSlvViscosityPtr => EOSSlvViscosityConstant

end subroutine EOSSlvSetViscosityConstant

! ************************************************************************** !

subroutine EOSSlvViscosityNoDerive(T,P,V_mix,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: V_mix   ! mixture viscosity
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: dum1, dum2

  call EOSSlvViscosityPtr(T,P,V_mix,PETSC_FALSE,dum1,dum2,ierr,table_idxs)

end subroutine EOSSlvViscosityNoDerive

! ************************************************************************** !
subroutine EOSSlvViscosityDerive(T,P,V_mix,dV_dT,dV_dP,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: V_mix   ! mixture viscosity
  PetscReal, intent(out) :: dV_dT,dv_dP   ! derivatives
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  call EOSSlvViscosityPtr(T,P,V_mix,PETSC_TRUE,dV_dT,dV_dP,ierr,table_idxs)

end subroutine EOSSlvViscosityDerive

! ************************************************************************** !

subroutine EOSSlvViscosityConstant(T, P, V, &
                                   calculate_derivative, dV_dT, dV_dP, &
                                   ierr, table_idxs)
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: V       ! mixture viscosity
  PetscBool, intent(in) :: calculate_derivative
  PetscReal, intent(out) :: dV_dT   ! derivative wrt temperature
  PetscReal, intent(out) :: dV_dP   ! derivative wrt pressure

  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  ierr=0

  V = constant_viscosity

  dV_dT=0.d0
  dV_dP=0.d0
  !if (calculate_derivative) call throwDerivativeError()

end subroutine EOSSlvViscosityConstant

! ************************************************************************** !

subroutine EOSSlvViscosityTable(T,P,V, &
                                   calculate_derivative,dV_dT,dV_dP,ierr,table_idxs)
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: V       ! mixture viscosity
  PetscBool, intent(in) :: calculate_derivative ! Request derivatives
  PetscReal, intent(out) :: dV_dT       ! derivative wrt temperature
  PetscReal, intent(out) :: dV_dP       ! derivative wrt pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  ierr=0

  call pvt_table%EOSPropGrad(T,P,EOS_VISCOSITY,V,dV_dT,dV_dP, &
                             ierr,table_idxs)

end subroutine EOSSlvViscosityTable

! ************************************************************************** !

subroutine EOSSlvDensityNoDerive(T,P,Rho_slv,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_slv ! slv density [kmol/m^3]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: dum1, dum2

  call EOSSlvDensityPtr(T,P,Rho_slv,dum1,dum2,ierr,table_idxs)

end subroutine EOSSlvDensityNoDerive

! ************************************************************************** !

subroutine EOSSlvDensityDerive(T,P,Rho_slv,dRho_dT,dRho_dP,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_slv ! slv density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative slv density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative slv density wrt pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  call EOSSlvDensityPtr(T,P,Rho_slv,dRho_dT,dRho_dP,ierr,table_idxs)

end subroutine EOSSlvDensityDerive

! ************************************************************************** !

subroutine EOSSlvEnergyNoDerive(T,P,H,U,ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscErrorCode, intent(out) :: ierr

  PetscReal :: dum1, dum2, dum3, dum4

  call EOSSlvEnergyPtr(T,P,H,dum1,dum2,U,dum3,dum4,ierr)

end subroutine EOSSlvEnergyNoDerive

! ************************************************************************** !

subroutine EOSSlvEnergyDerive(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
  PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
  PetscErrorCode, intent(out) :: ierr

  call EOSSlvEnergyPtr(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)

end subroutine EOSSlvEnergyDerive

! ************************************************************************** !

subroutine EOSSlvDenEnthNoDerive(T,P,Rho_slv,H,U,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_slv ! slv density [kmol/m^3]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: dum1, dum2, dum3, dum4, dum5, dum6

  call EOSSlvDensityEnergy1(T,P,Rho_slv,dum1,dum2, &
                            H,dum3,dum4,U,dum5,dum6,ierr,table_idxs)

end subroutine EOSSlvDenEnthNoDerive

! ************************************************************************** !

subroutine EOSSlvDenEnthDerive(T,P,Rho_slv,dRho_dT,dRho_dP, &
                               H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_slv ! slv density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative slv density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative slv density wrt pressuret
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
  PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  call EOSSlvDensityEnergy1(T,P,Rho_slv,dRho_dT,dRho_dP, &
                              H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr,table_idxs)

end subroutine EOSSlvDenEnthDerive

! ************************************************************************** !

subroutine EOSSlvDensityEnergy1(T,P,Rho_slv,dRho_dT,dRho_dP, &
                                H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr,&
                                table_idxs)
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_slv ! slv density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative slv density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative slv density wrt pressuret
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
  PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  call EOSSlvDensityPtr(T,P,Rho_slv,dRho_dT,dRho_dP,ierr,table_idxs)
  call EOSSlvEnergyPtr(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)

end subroutine EOSSlvDensityEnergy1

! ************************************************************************** !

subroutine EOSSlvDensityIdeal(T,P,Rho_slv,dRho_dT,dRho_dP,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_slv ! slv density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative slv density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative slv density wrt pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal  T_kelvin

  T_kelvin = T + 273.15d0
  Rho_slv = P / T_kelvin / IDEAL_GAS_CONSTANT * 1.d-3 ! mol/m^3 -> kmol/m^3

  dRho_dP =  Rho_slv / P
  dRho_dT = -Rho_slv / T_kelvin

  ierr=0

end subroutine EOSSlvDensityIdeal

! ************************************************************************** !

subroutine EOSSlvDensityEOSDBase(T,P,Rhs_slv,dRhs_dT,dRhs_dP,ierr,table_idxs)

!------------------------------------------------------------------------------
! Do density lookup for solvent
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : May 2018
!------------------------------------------------------------------------------

  implicit none

  PetscReal, intent(in)  :: T       ! temperature [C ]
  PetscReal, intent(in)  :: P       ! pressure    [Pa]
  PetscReal, intent(out) :: Rhs_slv ! solvent density [kmol/m^3]
  PetscReal, intent(out) :: dRhs_dT ! derivative solvent density wrt temperature
  PetscReal, intent(out) :: dRhs_dP ! derivative solvent density wrt pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  ierr=0
  call eos_dbase%EOSPropGrad(T,P,EOS_DENSITY,Rhs_slv,dRhs_dT,dRhs_dP,ierr)

  Rhs_slv= Rhs_slv / fmw_slv ! kmol/m^3

  dRhs_dT = dRhs_dT / fmw_slv
  dRhs_dP = dRhs_dP / fmw_slv

end subroutine EOSSlvDensityEOSDBase

subroutine EOSSlvDensityTable(T, P, Rho_slv, dRho_dT, dRho_dP, ierr, &
                                 table_idxs)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_slv ! oil density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative oil density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative oil density wrt pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  ierr = 0
  call pvt_table%EOSPropGrad(T,P,EOS_DENSITY,Rho_slv,dRho_dT,dRho_dP, &
                             ierr,table_idxs)

end subroutine EOSSlvDensityTable

! ************************************************************************** !

subroutine EOSSlvEnergyIdeal(T,P,H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
  PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
  PetscErrorCode, intent(out) :: ierr

! Cv_co2 units: J/mol-K, from Wiki at (15C,1 atm)
  PetscReal, parameter :: Cv_co2 = 28.075
  PetscReal :: T_energy
  PetscReal :: T_k

! T_energy is either T or T + 273.15
! do not change below
  T_energy = T + T_energy_offset
  T_k = T + 273.15d0
  U = Cv_co2 * T_energy * 1.d3  ! J/mol -> J/kmol
  H = U + IDEAL_GAS_CONSTANT * T_k * 1.d3 ! J/mol -> J/kmol

  dU_dP = 0.d0
  dU_dT = Cv_co2 * 1.d3
  dH_dP = 0.d0
  dH_dT = dU_dT + IDEAL_GAS_CONSTANT * 1.d3

  ierr=0

end subroutine EOSSlvEnergyIdeal

! ************************************************************************** !

subroutine EOSSlvDensityConstant(T,P,Rho_slv,dRho_dT,dRho_dP,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho_slv ! slv density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative slv density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative slv density wrt pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  ierr = 0
  Rho_slv = constant_density ! kmol/m^3

  dRho_dT = 0.d0
  dRho_dP = 0.d0

end subroutine EOSSlvDensityConstant

! **************************************************************************** !

subroutine EOSSlvSetPVDS(input,option)
  !
  ! Author: Dave Ponting
  ! Date  : May 2018
  !
  ! Set up a PVDS table

  use Option_module
  use Input_Aux_module
  use Lookup_Table_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option

  type(lookup_table_var_type), pointer :: db_var => null()
  character(len=MAXWORDLENGTH) :: internal_units, user_units
  PetscInt :: data_idx

  pvt_table => EOSTableCreate('PVDS',option)

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

  EOSSlvDensityPtr => EOSSlvDensityTable
  EOSSlvViscosityPtr => EOSSlvViscosityTable

end subroutine EOSSlvSetPVDS

! ************************************************************************** !

subroutine EOSSlvTableProcess(option)

  ! Processes slv pvt table - once the entire input deck has been read

  use Option_module

  implicit none

  type(option_type) :: option

  if (.not.associated(pvt_table)) return

  select case(pvt_table%name)
  case("PVDS")
      call pvt_table%ConvertFVFtoMolarDensity(fmw_slv,reference_density_kg)
  end select

end subroutine EOSSlvTableProcess

! ************************************************************************** !

subroutine EOSSlvInputRecord()
  ! 
  ! Prints ingested equation of state information to the input record file.
  ! 

  implicit none

  character(len=MAXWORDLENGTH) :: word1, word2
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'SOLVENT'

  ! slv density [kg/m^3]

  if (associated(EOSSlvDensityPtr,EOSSlvDensityIdeal)) then
    write(id,'(a29)',advance='no') 'slv density: '
    write(id,'(a)') 'ideal co2'
  endif
  if (associated(EOSSlvDensityPtr,EOSSlvDensityConstant)) then
    write(id,'(a29)',advance='no') 'slv density: '
    write(id,'(a)') 'constant'
  endif
  if (associated(EOSSlvDensityPtr,EOSSlvDensityConstant)) then
    write(id,'(a29)',advance='no') 'slv density: '
    write(id,'(a)') 'table'
  endif

  ! slv viscosity [Pa-s]
  if (associated(EOSSlvViscosityPtr,EOSSlvViscosityConstant)) then
    write(id,'(a29)',advance='no') 'slv viscosity: '
    write(word1,*) constant_viscosity
    write(id,'(a)') 'constant, ' // trim(word1) // ' Pa-sec'
  endif
  if (associated(EOSSlvViscosityPtr,EOSSlvViscosityTable)) then
    write(id,'(a29)',advance='no') 'slv viscosity: '
    write(id,'(a)') 'table'
  endif

  ! slv enthalpy [J/kmol]
  if ( associated(EOSSlvEnergyPtr,EOSSlvEnergyIdeal)) then
    write(id,'(a29)',advance='no') 'slv enthalpy: '
    write(id,'(a)') 'default, ideal co2'
  endif

  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'

end subroutine EOSSlvInputRecord

! ************************************************************************** !

subroutine EOSSlvTest(temp_low,temp_high,pres_low,pres_high, &
                      ntemp,npres,uniform_temp,uniform_pres,filename)

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
  PetscReal, allocatable :: internal_energy(:,:)
  PetscReal, allocatable :: viscosity(:,:)

  PetscInt :: itemp, ipres
  PetscReal :: ln_low, ln_high

  character(len=MAXWORDLENGTH) :: eos_density_name
  character(len=MAXWORDLENGTH) :: eos_energy_name
  character(len=MAXWORDLENGTH) :: eos_viscosity_name
  character(len=MAXSTRINGLENGTH) :: header, string

  allocate(temp(ntemp))
  temp = UNINITIALIZED_DOUBLE
  allocate(pres(npres))
  pres = UNINITIALIZED_DOUBLE
  allocate(density_kg(npres,ntemp))
  density_kg = UNINITIALIZED_DOUBLE
  allocate(viscosity(npres,ntemp))
  viscosity = UNINITIALIZED_DOUBLE
  allocate(enthalpy(npres,ntemp))
  enthalpy = UNINITIALIZED_DOUBLE
  allocate(internal_energy(npres,ntemp))
  internal_energy = UNINITIALIZED_DOUBLE

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
  if (associated(EOSSlvDensityPtr,EOSSlvDensityConstant)) then
    eos_density_name = 'Constant'
  else if (associated(EOSSlvDensityPtr,EOSSlvDensityIdeal)) then
    eos_density_name = 'Ideal Gas Law for CO2'
  else if (associated(EOSSlvDensityPtr,EOSSlvDensityTable)) then
    eos_density_name = 'PVDS table'
  else 
    eos_density_name = 'Unknown'
  endif

  ! energy
  if (associated(EOSSlvEnergyPtr,EOSSlvEnergyIdeal)) then
    eos_energy_name = 'Ideal Gas Law for CO2'
  else
    eos_energy_name = 'Unknown'
  endif

  ! viscosity
  if (associated(EOSSlvViscosityPtr,EOSSlvViscosityConstant)) then
    eos_viscosity_name = 'Constant'
  else if (associated(EOSSlvViscosityPtr,EOSSlvViscosityTable)) then
    eos_viscosity_name = 'PVDS table'
  else
    eos_viscosity_name = 'Unknown'
  endif

100 format(100es16.8)
  if (len_trim(filename) == 0) then
    string = 'eos_slv_test.txt'
  else
    string = filename
  endif
  open(unit=IUNIT_TEMP,file=string)
  header = 'T[C], P[Pa], &
    &Density (' // trim(eos_density_name) // ') [kmol/m^3], &
    &Enthalpy (' // trim(eos_energy_name) // ') [J/kmol], &
    &Internal Energy (' // trim(eos_energy_name) // ') [J/kmol], &
    &Viscosity (' // trim(eos_viscosity_name) // ') [Pa-s]'
  write(IUNIT_TEMP,'(a)') trim(header)
  write(IUNIT_TEMP,'(100i9)') ntemp, npres
  do itemp = 1, ntemp
    do ipres = 1, npres
      write(IUNIT_TEMP,100) temp(itemp), pres(ipres), &
            density_kg(ipres,itemp), enthalpy(ipres,itemp), &
            internal_energy(ipres,itemp), viscosity(ipres,itemp)
    enddo
  enddo
  close(IUNIT_TEMP)

  deallocate(temp)
  deallocate(pres)
  deallocate(density_kg)
  deallocate(enthalpy)
  deallocate(internal_energy)
  deallocate(viscosity)

end subroutine EOSSlvTest

! ************************************************************************** !

subroutine throwDerivativeError()

!------------------------------------------------------------------------------
! Issue error if derivatives requestsed (not available at present)
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : May 2018
!------------------------------------------------------------------------------

  use Option_module

  implicit none

  type(option_type) :: option

  option%io_buffer = 'Derivatives are not available'
  call printErrMsg(option)

end subroutine throwDerivativeError

! ************************************************************************** !

subroutine EOSSlvDBaseDestroy()

  implicit none

  !destroy databases
  call EOSDatabaseDestroy(eos_dbase)

  nullify(pvt_table)

end subroutine EOSSlvDBaseDestroy

! ************************************************************************** !

end module EOS_Slv_module
