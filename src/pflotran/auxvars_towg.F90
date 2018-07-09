module AuxVars_TOWG_module

  use PFLOTRAN_Constants_module

  use AuxVars_Base_module
  use AuxVars_Flow_module
  use AuxVars_FlowEnergy_module
  use petscsys

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  PetscBool, public :: towg_analytical_derivatives = PETSC_FALSE
  PetscBool, public :: towg_analytical_derivatives_compare = PETSC_FALSE
  PetscReal, public :: towg_dcomp_tol = 1.d0

#if 0
  type, public :: towg_derivative_auxvar_type
    ! _dp is always with respect to oil pressure, no other pressures are used as solution variables.
    ! _dt is with respect to temperature
    ! _dsat_[x] is with respect to saturtion of phase [x] 

    PetscReal, pointer :: dp_dsat_oil(:) !! pressure w.r.t. saturation, will take into account 
                                         !! dependence of liquid pressure on saturation
    PetscReal, pointer :: dp_dsat_gas(:) 
    PetscReal, pointer :: dsat_dp(:) !! saturation w.r.t. pressure
    PetscReal, pointer :: dden_dp(:) !! density w.r.t. pressure
    !! dden_dsat_oil
    !! dden_dsat_gas       etc
    PetscReal, pointer :: dden_dt(:)   !! density w.r.t. temperature
    PetscReal, pointer :: dsat_dt(:)   !! saturattion w.r.t. temperature
    PetscReal, pointer :: dU_dp(:)     !! internal energy w.r.t. pressure
    PetscReal, pointer :: dU_dt(:)
    PetscReal, pointer :: dH_dp(:)
    PetscReal, pointer :: dH_dt(:)

    PetscReal, pointer :: dmobility(:,:) !! this seems to be the only one
                                         !! that qualifies for the full 
                                         !! square matrix treatment
                                         !! (2x3 here)

    PetscReal, pointer :: dpor_dp_oil

  contains
    procedure, public :: Init => AuxVarDerivsTowgInit
    procedure, public :: Strip => AuxVarDerivsTowgStrip
  end type towg_derivative_auxvar_type
#endif


  type, public, extends(auxvar_flow_energy_type) :: auxvar_towg_type
    PetscInt :: istate_store(2) ! 1 = previous timestep; 2 = previous iteration
    type(tl_auxvar_type), pointer :: tl=>null()
    type(bo_auxvar_type), pointer :: bo=>null()
    !PetscBool :: hasDerivatives
    !type(towg_derivative_auxvar_type), pointer :: d
  contains
    procedure, public :: Init => AuxVarTOWGInit
    procedure, public :: Strip => AuxVarTOWGStrip
    procedure, public :: InitTL
    procedure, public :: StripTL
    procedure, public :: InitBO
    procedure, public :: StripBO
  end type auxvar_towg_type

  type, public ::  tl_auxvar_type
    PetscReal :: den_oil_eff_kg
    PetscReal :: den_gas_eff_kg
  end type tl_auxvar_type

  type, public ::  bo_auxvar_type
    PetscReal :: bubble_point
    PetscReal :: xo
    PetscReal :: xg
    PetscBool :: hasDerivatives
    type(bo_derivative_auxvar_type), pointer :: d
  end type bo_auxvar_type

  type, public ::  bo_derivative_auxvar_type
    !!! derivatives w.r.t. bubble point
    !!! ...
    PetscReal :: ddeno_dpb !! oil density w.r.t. bubble pt
    PetscReal :: dHo_dpb,dUo_dpb

    !!! deivaties of xo, xg with respect to various
    !!! solution variables
    !!! ...
    PetscReal :: dxo_dpb !! xo w.r.t. bubble pt
    PetscReal :: dxg_dpb !! xg w.r.t. bubble pt
    PetscReal :: dxo_dt  !! xo w.r.t. temperature
    PetscReal :: dxg_dt  !! xg w.r.t. temperature
  contains
    procedure, public :: Init => AuxVarDerivsBOTOWGInit
    procedure, public :: Strip => AuxVarDerivsBOTOWGStrip
  end type bo_derivative_auxvar_type

  public :: AuxVarTOWGStrip

contains

! ************************************************************************** !

subroutine AuxVarTOWGInit(this,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: PAolo Orsini
  ! Date: 11/07/16
  ! 

  use Option_module

  implicit none
  
  class(auxvar_towg_type) :: this
  type(option_type) :: option

  this%effective_porosity = 0.d0
  this%pert = 0.d0

#if 0
  if (towg_analytical_derivatives) then
    this%hasDerivatives = PETSC_TRUE
    allocate(this%d)
    call this%d%Init(option)
  else
    this%hasDerivatives = PETSC_FALSE
  endif
#endif
  
  call AuxVarFlowInit(this,option)

  call AuxVarFlowEnergyInit(this,option)

  this%istate_store = 0

end subroutine AuxVarTOWGInit

! ************************************************************************** !

#if 0
subroutine AuxVarDerivsTowgInit(this,option)
  ! 
  ! 

  use Option_module

  implicit none
  
  class(towg_derivative_auxvar_type) :: this
  type(option_type) :: option

  allocate(this%dp_dsat_oil(option%nphase))
  this%dp_dsat_oil = 0.d0
  allocate(this%dp_dsat_gas(option%nphase))
  this%dp_dsat_gas = 0.d0

  allocate(this%dsat_dp(option%nphase))
  this%dsat_dp = 0.d0
  allocate(this%dden_dp(option%nphase))
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

  allocate(this%dpor_dp_oil)
  this%dpor_dp_oil = 0.d0

end subroutine AuxVarDerivsTowgInit
#endif

! ************************************************************************** !

subroutine InitTL(this,option)
  ! 
  ! Initialize auxiliary object or Todd Longstaff model
  ! 
  ! Author: Paolo Orsini
  ! Date: 07/06/17
  ! 

  use Option_module

  implicit none
  
  class(auxvar_towg_type) :: this
  type(option_type) :: option

  allocate(this%tl)

  this%tl%den_oil_eff_kg = 0.0
  this%tl%den_gas_eff_kg = 0.0

end subroutine InitTL

!--Routine to initialise the black oil substructure----------------------------

subroutine InitBO(this,option)

!------------------------------------------------------------------------------
! Used in TOWG_BLACK_OIL and TOWG_SOLVENT_TL, initialises the bo sub-structure
! of auxvars. Contains the bubble point and the oil mole fractions.
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Sep 2017
!------------------------------------------------------------------------------

  use Option_module

  implicit none

  class(auxvar_towg_type) :: this
  type(option_type) :: option

  allocate(this%bo)

  this%bo%bubble_point=0.0
  this%bo%xg          =0.0
  this%bo%xo          =0.0

  if (towg_analytical_derivatives) then
    !this%hasDerivatives = PETSC_TRUE
    this%bo%hasDerivatives = PETSC_TRUE
    allocate(this%bo%d)
    call this%bo%d%Init(option)
  else
    !this%hasDerivatives = PETSC_FALSE
    this%bo%hasDerivatives = PETSC_FALSE
  endif

end subroutine InitBO

!--Routine to initialise the TOWG substructure-------------------------------

subroutine AuxVarTOWGStrip(this)
  ! 
  ! AuxVarTOWGStrip: Deallocates a towg auxiliary object
  ! 
  ! Author: Paolo Orsini
  ! Date: 10/30/16
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  class(auxvar_towg_type) :: this

  call AuxVarFlowStrip(this)

  call AuxVarFlowEnergyStrip(this)

  if (associated(this%tl)) call this%StripTL()
  if (associated(this%bo)) call this%StripBO()

end subroutine AuxVarTOWGStrip


! ************************************************************************** !

#if 0
subroutine AuxVarDerivsTowgStrip(this)
 ! 
  use Utility_module, only : DeallocateArray

  implicit none

  class(towg_derivative_auxvar_type) :: this

  call DeallocateArray(this%dp_dsat_oil)
  call DeallocateArray(this%dp_dsat_gas)

  call DeallocateArray(this%dsat_dp)
  call DeallocateArray(this%dden_dp)
  call DeallocateArray(this%dsat_dt)
  call DeallocateArray(this%dden_dt)

  call DeallocateArray(this%dU_dt)
  call DeallocateArray(this%dU_dP)
  call DeallocateArray(this%dH_dt)
  call DeallocateArray(this%dH_dP)

  call deallocatearray(this%dmobility)

end subroutine AuxVarDerivsTowgStrip
#endif

! ************************************************************************** !

subroutine StripTL(this)
  ! 
  ! StripTL: Deallocates a the Todd Longstaff component of 
  !          the towg auxiliary object
  ! 
  ! Author: Paolo Orsini
  ! Date: 07/06/17
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  class(auxvar_towg_type) :: this

  deallocate(this%tl)

end subroutine StripTL

!--Routine to strip the black oil substructure---------------------------------

subroutine StripBO(this)

!------------------------------------------------------------------------------
! Used in TOWG_BLACK_OIL and TOWG_SOLVENT_TL.
! Deallocate the bo sub-structure of auxvars
!------------------------------------------------------------------------------
! Author: Dave Ponting
! Date  : Sep 2017
!------------------------------------------------------------------------------

  use Utility_module, only : DeallocateArray

  implicit none

  class(auxvar_towg_type) :: this

  deallocate(this%bo)

end subroutine StripBO

! ************************************************************************** !

subroutine AuxVarDerivsBOTowgInit(this,option)
  ! 
  ! 

  use Option_module

  implicit none
  
  class(bo_derivative_auxvar_type) :: this
  type(option_type) :: option

   this%ddeno_dpb = 0.d0 !! oil density w.r.t. bubble pt
   this%dHo_dpb = 0.d0
   this%dUo_dpb = 0.d0

   this%dxo_dpb = 0.d0!! xo w.r.t. bubble pt
   this%dxg_dpb = 0.d0!! xg w.r.t. bubble pt
   this%dxo_dt = 0.d0 !! xo w.r.t. temperature
   this%dxg_dt = 0.d0 !! xg w.r.t. temperature



end subroutine AuxVarDerivsBOTowgInit

! ************************************************************************** !

subroutine AuxVarDerivsBOTOWGStrip(this,option)
  use Option_module

  implicit none
  
  class(bo_derivative_auxvar_type) :: this
  type(option_type) :: option
end subroutine AuxVarDerivsBOTOWGStrip


end module AuxVars_TOWG_module

