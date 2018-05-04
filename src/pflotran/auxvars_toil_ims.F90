module AuxVars_TOilIms_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  use AuxVars_Base_module
  use AuxVars_Flow_module
  use AuxVars_FlowEnergy_module

  implicit none

  private

  PetscBool, public :: toil_analytical_derivatives = PETSC_FALSE
  PetscBool, public :: toil_analytical_derivatives_compare = PETSC_FALSE
  PetscReal, public :: toil_dcomp_tol = 1.d0
  PetscBool, public :: toil_GP = PETSC_FALSE

  type, public :: toil_derivative_auxvar_type
    PetscReal, pointer :: dp_dsat(:) !! pressure w.r.t. saturation, will take into account 
                                     !! dependence of liquid pressure on saturation
    PetscReal, pointer :: dsat_dp(:,:) !! saturation w.r.t. pressure !! USED?
    PetscReal, pointer :: dden_dp(:,:) !! density w.r.t. pressure
    PetscReal, pointer :: dden_dt(:)   !! density w.r.t. temperature
    PetscReal, pointer :: dsat_dt(:)   !! saturattion w.r.t. temperature
    PetscReal, pointer :: dU_dp(:)     !! internal energy w.r.t. pressure !! USED?
    PetscReal, pointer :: dU_dt(:)
    PetscReal, pointer :: dH_dp(:)
    PetscReal, pointer :: dH_dt(:)

    PetscReal, pointer :: dmobility(:,:) !! this seems to be the only one
                                         !! that qualifies for the full 
                                         !! square matrix treatment
                                         !! (2x3 here)

    PetscReal, pointer :: dpor_dp

  contains
    procedure, public :: Init => AuxVarDerivsTOilImsInit
    procedure, public :: Strip => AuxVarDerivsTOilImsStrip
  end type toil_derivative_auxvar_type

  type, public, extends(auxvar_flow_energy_type) :: auxvar_toil_ims_type
    ! no data at the moment
    PetscBool :: hasDerivatives
    type(toil_derivative_auxvar_type), pointer :: d
  contains
    procedure, public :: Init => AuxVarTOilImsInit
    procedure, public :: Strip => AuxVarTOilImsStrip
  end type auxvar_toil_ims_type

  public :: AuxVarTOilImsStrip

contains

! ************************************************************************** !

subroutine AuxVarTOilImsInit(this,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: PAolo Orsini
  ! Date: 5/27/16
  !

  use Option_module

  implicit none

  class(auxvar_toil_ims_type) :: this
  type(option_type) :: option

  this%effective_porosity = 0.d0
  this%pert = 0.d0

  if (toil_analytical_derivatives) then
    this%hasDerivatives = PETSC_TRUE

#if 0
    allocate(this%df) 
    allocate(this%de) 
    call this%df%Init(option)
    call this%de%Init(option)
#endif
    allocate(this%d)
    call this%d%Init(option)
  else
    this%hasDerivatives = PETSC_FALSE
  endif

  !PO to do:
  !allign TOIL_IMS AuxVarCompute to new pc array (see AuxVarFlowInit),
  !then remove the following iniitlisation block, and replace with a call to
  !AuxVarFlowInit
  !the indended part below is the block to replace with AuxVarFlowInit
    !two phases (water,oil) and capillary pressure
    ! allocate(this%pres(option%nphase+ONE_INTEGER))
    ! this%pres = 0.d0
    ! allocate(this%sat(option%nphase))
    ! this%sat = 0.d0
    ! allocate(this%den(option%nphase))
    ! this%den = 0.d0
    ! allocate(this%den_kg(option%nphase))
    ! this%den_kg = 0.d0
    ! allocate(this%mobility(option%nphase))
    ! this%mobility = 0.d0
    !if visc needed only for toil_ims, use AuxVarFlowInit to allocate sat,
    !den, den_kg, mobility then allocate viscosity in here
    !more likley to allocate viscosity for all flow module though
    ! allocate(this%viscosity(option%nphase))
    ! this%viscosity = 0.d0
  ! end block to replace with AuxVarFlowInit
  call AuxVarFlowInit(this,option)

  !PO to do:
  !replace the indended block below with AuxVarFlowEnergyInit
    ! this%temp = 0.d0
    ! allocate(this%H(option%nphase))
    ! this%H = 0.d0
    ! allocate(this%U(option%nphase))
    ! this%U = 0.d0
  !end block to replace with AuxVarFlowEnergyInit
  call AuxVarFlowEnergyInit(this,option)


end subroutine AuxVarTOilImsInit

! ************************************************************************** !

subroutine AuxVarDerivsTOilImsInit(this,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: PAolo Orsini
  ! Date: 5/27/16
  ! 

  use Option_module

  implicit none
  
  class(toil_derivative_auxvar_type) :: this
  type(option_type) :: option

  allocate(this%dp_dsat(option%nphase))
  this%dp_dsat= 0.d0

  allocate(this%dsat_dp(option%nphase, 1))
  this%dsat_dp = 0.d0
  allocate(this%dden_dp(option%nphase, 1))
  this%dden_dp = 0.d0
  allocate(this%dsat_dt(option%nphase))
  this%dsat_dt = 0.d0
  allocate(this%dden_dt(option%nphase))
  this%dden_dt = 0.d0

  allocate(this%dU_dt(option%nphase))
  this%dU_dt = 0.d0
  allocate(this%dU_dp(option%nphase))
  this%dU_dP = 0.d0

  allocate(this%dH_dt(option%nphase))
  this%dH_dt = 0.d0
  allocate(this%dH_dp(option%nphase))
  this%dH_dP = 0.d0

  allocate(this%dmobility(option%nphase, 3))
  this%dmobility = 0.d0

  allocate(this%dpor_dp)
  this%dpor_dp = 0.d0

end subroutine AuxVarDerivsTOilImsInit

! ************************************************************************** !

subroutine AuxVarTOilImsStrip(this)
  !
  ! TOilImsAuxVarDestroy: Deallocates a toil_ims auxiliary object
  !
  ! Author: Paolo Orsini
  ! Date: 10/30/16
  !
  use Utility_module, only : DeallocateArray

  implicit none

  class(auxvar_toil_ims_type) :: this

  call AuxVarFlowStrip(this)

  call AuxVarFlowEnergyStrip(this)

end subroutine AuxVarTOilImsStrip

! ************************************************************************** !

subroutine AuxVarDerivsTOilImsStrip(this)
 ! 
  use Utility_module, only : DeallocateArray

  implicit none

  class(toil_derivative_auxvar_type) :: this

  call DeallocateArray(this%dp_dsat)

  call DeallocateArray(this%dsat_dp)
  call DeallocateArray(this%dden_dp)
  call DeallocateArray(this%dsat_dt)
  call DeallocateArray(this%dden_dt)

  call DeallocateArray(this%dU_dt)
  call DeallocateArray(this%dU_dP)
  call DeallocateArray(this%dH_dt)
  call DeallocateArray(this%dH_dP)

  call deallocatearray(this%dmobility)

end subroutine AuxVarDerivsTOilImsStrip

! ************************************************************************** !

subroutine AuxVarDerivsTOilIMSCopy(auxvar,auxvar2,option)
  ! 
  ! Copies an auxiliary derivatives variable
  ! 

  use Option_module

  implicit none
  
  class(toil_derivative_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%dp_dsat= auxvar%dp_dsat
  auxvar2%dsat_dp = auxvar%dsat_dp
  auxvar2%dden_dp = auxvar%dden_dp
  auxvar2%dden_dt = auxvar%dden_dt
  auxvar2%dsat_dt = auxvar%dsat_dt
  auxvar2%dU_dp = auxvar%dU_dp
  auxvar2%dU_dt = auxvar%dU_dt
  auxvar2%dH_dp = auxvar%dH_dp
  auxvar2%dH_dt = auxvar%dH_dt
  auxvar2%dmobility = auxvar%dmobility
  auxvar2%dpor_dp = auxvar%dpor_dp


end subroutine AuxVarDerivsTOilIMSCopy

! ************************************************************************** !

end module AuxVars_TOilIms_module
