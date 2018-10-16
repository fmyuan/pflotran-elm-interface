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
    PetscReal, pointer :: pres(:)   ! (iphase)
    PetscReal, pointer :: sat(:)    ! (iphase)
    PetscReal, pointer :: pc(:)     ! capillary pressure (iphase-1)
    PetscReal, pointer :: den(:)    ! (iphase) kmol/m^3 phase
    PetscReal, pointer :: den_kg(:) ! (iphase) kg/m^3 phase
    PetscReal, pointer :: mobility(:) ! relative perm / dynamic viscosity
    PetscReal, pointer :: viscosity(:) ! dynamic viscosity
    PetscInt, pointer :: table_idx(:)

    ! derivatives
    PetscBool :: has_derivs
    PetscReal, pointer :: D_pres(:,:)   ! (iphase)
    PetscReal, pointer :: D_sat(:,:)    ! (iphase)
    PetscReal, pointer :: D_pc(:,:)     ! capillary pressure (iphase-1)
    PetscReal, pointer :: D_den(:,:)    ! (iphase) kmol/m^3 phase
    PetscReal, pointer :: D_den_kg(:,:) ! (iphase) kg/m^3 phase
    PetscReal, pointer :: D_mobility(:,:) ! relative perm / dynamic viscosity
    PetscReal, pointer :: D_por(:)

  contains
    !procedure, public :: Init => InitAuxVarFlow
  end type auxvar_flow_type

  public :: AuxVarFlowInit, AuxVarFlowStrip

contains

! ************************************************************************** !
subroutine AuxVarFlowInit(this,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: PAolo Orsini
  ! Date: 5/27/16
  !

  use Option_module

  implicit none

  class(auxvar_flow_type) :: this
  type(option_type) :: option

  allocate(this%pres(option%nphase))
  this%pres = 0.d0
  allocate(this%pc(option%nphase - ONE_INTEGER))
  this%pc = 0.0d0
  allocate(this%sat(option%nphase))
  this%sat = 0.d0
  allocate(this%den(option%nphase))
  this%den = 0.d0
  allocate(this%den_kg(option%nphase))
  this%den_kg = 0.d0
  allocate(this%mobility(option%nphase))
  this%mobility = 0.d0
  allocate(this%viscosity(option%nphase))
  this%viscosity = 0.d0
  if (option%num_table_indices > 0) then
    allocate(this%table_idx(option%num_table_indices))
    this%table_idx = 1
  else
    nullify(this%table_idx)
  endif

  this%has_derivs = PETSC_FALSE
  if (.not.option%flow%numerical_derivatives) then

    this%has_derivs = PETSC_TRUE

    allocate(this%D_pres(option%nphase,option%nflowdof))
    this%D_pres = 0.d0
    allocate(this%D_sat(option%nphase,option%nflowdof))
    this%D_sat = 0.d0
    allocate(this%D_pc(option%nphase - ONE_INTEGER,option%nflowdof))
    this%D_pc = 0.d0
    allocate(this%D_den(option%nphase,option%nflowdof))
    this%D_den = 0.d0
    allocate(this%D_den_kg(option%nphase,option%nflowdof))
    this%D_den_kg = 0.d0
    allocate(this%D_mobility(option%nphase,option%nflowdof))
    this%D_mobility = 0.d0
    allocate(this%D_por(option%nflowdof))
    this%D_por = 0.d0
  endif


end subroutine AuxVarFlowInit

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
  call DeallocateArray(this%pc)
  call DeallocateArray(this%sat)
  call DeallocateArray(this%den)
  call DeallocateArray(this%den_kg)
  call DeallocateArray(this%mobility)
  call DeallocateArray(this%viscosity)
  call DeallocateArray(this%table_idx)

  if (this%has_derivs) then 
    call DeallocateArray(this%D_pres)
    call DeallocateArray(this%D_sat)
    call DeallocateArray(this%D_pc)
    call DeallocateArray(this%D_den)
    call DeallocateArray(this%D_den_kg)
    call DeallocateArray(this%D_mobility)
    call DeallocateArray(this%D_por)
  endif

end subroutine AuxVarFlowStrip
! ************************************************************************** !

end module AuxVars_Flow_module
