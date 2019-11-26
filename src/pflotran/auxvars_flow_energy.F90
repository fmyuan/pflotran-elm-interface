module AuxVars_Flow_Energy_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  use AuxVars_Base_module
  use AuxVars_Flow_module

  implicit none
  
  private 

  type, public, extends(auxvar_flow_type) :: auxvar_flow_energy_type
    PetscReal :: TK            ! temperature for all phases in Kelvin
    PetscReal, pointer :: H(:) ! (nfluids+1) MJ/kmol
    PetscReal, pointer :: U(:) ! (nfluids+1) MJ/kmol

    ! derivatives
    PetscReal, pointer :: D_H(:,:) ! MJ/kmol
    PetscReal, pointer :: D_U(:,:) ! MJ/kmol
  contains
   !..............
  end type auxvar_flow_energy_type

  public :: AuxVarFlowEnergyInit, &
            AuxVarFlowEnergyCopy, &
            AuxVarFlowEnergyStrip

contains

! ************************************************************************** !
subroutine AuxVarFlowEnergyInit(this,option)
  ! 
  ! Initialize energy auxiliary variables, in addition to flow's
  !
  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019
 !

  use Option_module

  implicit none

  class (auxvar_flow_energy_type) :: this
  type(option_type) :: option

  this%TK = 0.d0
  allocate(this%H(option%nfluids+ONE_INTEGER))  ! 'ONE_INTEGER' is for fluid's solid phase
  this%H = 0.d0
  allocate(this%U(option%nfluids+ONE_INTEGER))
  this%U = 0.d0

  nullify(this%D_H)
  nullify(this%D_U)

  this%has_derivs = option%flow%numerical_derivatives
  if (.not.option%flow%numerical_derivatives) then
    this%has_derivs = PETSC_TRUE
    allocate(this%D_H(option%nfluids+ONE_INTEGER,option%nflowdof))
    this%D_H = 0.d0
    allocate(this%D_U(option%nfluids+ONE_INTEGER,option%nflowdof))
    this%D_U = 0.d0
  endif

  call AuxVarFlowInit(this, option)

end subroutine AuxVarFlowEnergyInit

! ************************************************************************** !
subroutine AuxVarFlowEnergyCopy(auxvar, auxvar2)
  !
  ! dupliacate energy auxiliary variables, in addition to flow's
  !
  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019
  !

  use Option_module

  implicit none

  class(auxvar_flow_energy_type) :: auxvar
  class(auxvar_flow_energy_type) :: auxvar2

  auxvar2%TK= auxvar%TK
  auxvar2%H = auxvar%H
  auxvar2%U = auxvar%U
  auxvar2%has_derivs = auxvar%has_derivs
  if(auxvar%has_derivs) then
    auxvar2%D_H = auxvar%D_H
    auxvar2%D_U = auxvar%D_U
  endif

  call AuxVarFlowCopy(auxvar, auxvar2)

end subroutine AuxVarFlowEnergyCopy

! ************************************************************************** !

subroutine AuxVarFlowEnergyStrip(this)
  ! 
  ! AuxVarFlowDestroy: Deallocates a toil_ims auxiliary object
  ! 
  ! Author: Paolo Orsini
  ! Date: 8/5/16
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  class(auxvar_flow_energy_type) :: this

  call DeallocateArray(this%H)  
  call DeallocateArray(this%U)  

  call DeallocateArray(this%D_H)  
  call DeallocateArray(this%D_U)

  call AuxVarFlowStrip(this)

end subroutine AuxVarFlowEnergyStrip

! ************************************************************************** !

end module AuxVars_Flow_Energy_module

