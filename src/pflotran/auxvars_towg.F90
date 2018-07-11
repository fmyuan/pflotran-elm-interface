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


  type, public, extends(auxvar_flow_energy_type) :: auxvar_towg_type
    PetscInt :: istate_store(2) ! 1 = previous timestep; 2 = previous iteration
    type(tl_auxvar_type), pointer :: tl=>null()
    type(bo_auxvar_type), pointer :: bo=>null()
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

    PetscBool :: has_derivs
    ! new generalized versions of derivatives:
    PetscReal, pointer :: D_xo(:)   ! (idof)
    PetscReal, pointer :: D_xg(:)   ! (idof)

    !type(bo_derivative_auxvar_type), pointer :: d
  end type bo_auxvar_type

  type, public ::  bo_derivative_auxvar_type

    ! deivaties of xo, xg with respect to various
    ! solution variables
    PetscReal :: dxo_dpb ! xo w.r.t. bubble pt
    PetscReal :: dxg_dpb ! xg w.r.t. bubble pt
    PetscReal :: dxo_dt  ! xo w.r.t. temperature
    PetscReal :: dxg_dt  ! xg w.r.t. temperature

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

  call AuxVarFlowInit(this,option)

  call AuxVarFlowEnergyInit(this,option)

  this%istate_store = 0

end subroutine AuxVarTOWGInit

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

  !if (towg_analytical_derivatives) then
  if (.NOT. option%flow%numerical_derivatives) then
    !this%hasDerivatives = PETSC_TRUE
    this%bo%has_derivs= PETSC_TRUE
    !allocate(this%bo%d)
    !call this%bo%d%Init(option)

    !!! new general derivs:
    allocate(this%bo%D_xo(option%nflowdof))
    this%bo%D_xo = 0.d0
    allocate(this%bo%D_xg(option%nflowdof))
    this%bo%D_xg = 0.d0

  else
    !this%hasDerivatives = PETSC_FALSE
    this%bo%has_derivs = PETSC_FALSE
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

  if (this%has_derivs) then
    call DeallocateArray(this%bo%D_xo)
    call DeallocateArray(this%bo%D_xg)
  endif

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

