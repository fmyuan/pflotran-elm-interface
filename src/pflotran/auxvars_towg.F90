module AuxVars_TOWG_module

  use PFLOTRAN_Constants_module

  use AuxVars_Base_module
  use AuxVars_Flow_module
  use AuxVars_FlowEnergy_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

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
  end type bo_auxvar_type

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

  deallocate(this%bo)

end subroutine StripBO

end module AuxVars_TOWG_module

