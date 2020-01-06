module AuxVars_Flow_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  use AuxVars_Base_module

  implicit none

  PetscReal,public :: flow_aux_debug_tol = 1.d-1
  PetscReal,public :: flow_aux_debug_reltol = 1.d-1
  PetscBool,public :: flow_aux_use_GP = PETSC_FALSE

  private

  type, public, extends(auxvar_base_type) :: auxvar_flow_type
    PetscReal, pointer :: sat(:)    ! (nfluids+1) unitless: the '+1' is for fluid in solid phase (e.g. ice)
    PetscReal, pointer :: den(:)    ! (nfluids+1) kmol/m^3:
    PetscReal, pointer :: den_kg(:) ! (nfluids+1) kg/m^3  :
    PetscReal, pointer :: pres(:)      ! (nfluids): exlcuding solid phase
    PetscReal, pointer :: mobility(:)  ! relative permissivity/dynamic viscosity
    PetscReal, pointer :: viscosity(:) ! dynamic viscosity
    PetscReal, pointer :: xmol(:,:)    ! (nflowspec,nfluids) transportants in fluids
    PetscReal          :: por          ! porosity of porous media for flow
    PetscReal          :: pc           ! liq. fluid capillary pressure only (specific variables useful)

    ! derivatives
    PetscBool :: has_derivs
    PetscReal, pointer :: D_sat(:,:)    ! (nfluids+1, nflowdof)
    PetscReal, pointer :: D_den(:,:)    ! (nfluids+1, nflowdof) kmol/m^3 phase
    PetscReal, pointer :: D_den_kg(:,:) ! (nfluids+1, nflowdof) kg/m^3 phase
    PetscReal, pointer :: D_pres(:,:)   ! (nfluids, nflowdof)
    PetscReal, pointer :: D_mobility(:,:)  ! (nfluids, nflowdof) relative perm./visc.
    PetscReal, pointer :: D_viscosity(:,:) ! (nfluids, nflowdof) dynamic viscosity
    PetscReal, pointer :: D_por(:)         ! (nflowdof)
    PetscReal, pointer :: D_pc(:)          ! liq. capillary pressure (nflowdof)

  contains
    !
  end type auxvar_flow_type

  public :: AuxVarFlowInit, AuxVarFlowCopy, AuxVarFlowStrip

contains

! ************************************************************************** !
subroutine AuxVarFlowInit(this,option)
  !
  ! initialize flow auxiliary variables
  !
  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019
  !

  use Option_module

  implicit none

  class(auxvar_flow_type) :: this
  type(option_type) :: option

  nullify(this%pres)
  nullify(this%sat)
  nullify(this%den)
  nullify(this%den_kg)
  nullify(this%mobility)
  nullify(this%viscosity)

  nullify(this%D_pres)
  nullify(this%D_sat)
  nullify(this%D_pc)
  nullify(this%D_den)
  nullify(this%D_den_kg)
  nullify(this%D_mobility)
  nullify(this%D_por)

  !
  allocate(this%pres(option%nfluids))
  this%pres = 0.d0
  this%pc = 0.0d0
  allocate(this%sat(option%nfluids))
  this%sat = 0.d0
  allocate(this%den(option%nfluids))
  this%den = 0.d0
  allocate(this%den_kg(option%nfluids))
  this%den_kg = 0.d0
  allocate(this%mobility(option%nfluids))
  this%mobility = 0.d0
  allocate(this%viscosity(option%nfluids))
  this%viscosity = 0.d0
  this%por = 0.d0

  this%has_derivs = option%flow%numerical_derivatives
  if (.not.option%flow%numerical_derivatives) then

    this%has_derivs = PETSC_TRUE

    allocate(this%D_pres(option%nfluids,option%nflowdof))
    this%D_pres = 0.d0
    allocate(this%D_sat(option%nfluids,option%nflowdof))
    this%D_sat = 0.d0
    allocate(this%D_pc(option%nflowdof))
    this%D_pc = 0.d0
    allocate(this%D_den(option%nfluids,option%nflowdof))
    this%D_den = 0.d0
    allocate(this%D_den_kg(option%nfluids,option%nflowdof))
    this%D_den_kg = 0.d0
    allocate(this%D_mobility(option%nfluids,option%nflowdof))
    this%D_mobility = 0.d0
    allocate(this%D_viscosity(option%nfluids,option%nflowdof))
    this%D_viscosity = 0.d0
    allocate(this%D_por(option%nflowdof))
    this%D_por = 0.d0
  endif

end subroutine AuxVarFlowInit

! ************************************************************************** !
subroutine AuxVarFlowCopy(auxvar, auxvar2)
  !
  ! dupliacate flow auxiliary variables
  !
  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019
  !


  implicit none

  class (auxvar_flow_type) :: auxvar
  class (auxvar_flow_type) :: auxvar2

  auxvar2%pres = auxvar%pres
  auxvar2%pc = auxvar%pc
  auxvar2%sat = auxvar%sat
  auxvar2%den = auxvar%den
  auxvar2%den_kg = auxvar%den_kg
  auxvar2%mobility = auxvar%mobility
  auxvar2%viscosity = auxvar%viscosity
  if (auxvar%has_derivs) then
    auxvar2%D_pres = auxvar%D_pres
    auxvar2%D_sat  = auxvar%D_sat
    auxvar2%D_pc   = auxvar%D_pc
    auxvar2%D_den  = auxvar%D_den
    auxvar2%D_den_kg   = auxvar%D_den_kg
    auxvar2%D_mobility = auxvar%D_mobility
    auxvar2%D_viscosity = auxvar%D_viscosity
    auxvar2%D_por      = auxvar%D_por
  endif

end subroutine AuxVarFlowCopy

! ************************************************************************** !

subroutine AuxVarFlowStrip(this)
  !
  ! AuxVarFlowDestroy: Deallocates a flow auxiliary object
  !
  ! Author: Paolo Orsini
  ! Date: 8/5/16
  !
  use Utility_module, only : DeallocateArray

  implicit none

  class(auxvar_flow_type) :: this

  call DeallocateArray(this%pres)
  call DeallocateArray(this%sat)
  call DeallocateArray(this%den)
  call DeallocateArray(this%den_kg)
  call DeallocateArray(this%mobility)
  call DeallocateArray(this%viscosity)

  if (this%has_derivs) then 
    call DeallocateArray(this%D_pres)
    call DeallocateArray(this%D_sat)
    call DeallocateArray(this%D_pc)
    call DeallocateArray(this%D_den)
    call DeallocateArray(this%D_den_kg)
    call DeallocateArray(this%D_mobility)
    call DeallocateArray(this%D_viscosity)
    call DeallocateArray(this%D_por)
  endif

end subroutine AuxVarFlowStrip
! ************************************************************************** !

end module AuxVars_Flow_module
