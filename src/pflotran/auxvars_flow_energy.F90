module AuxVars_Flow_Energy_module

  ! variables of flow model Extended for flow with thermal processes

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  use AuxVars_Flow_module

  implicit none
  
  private 

  type, public, extends(auxvar_flow_type) :: auxvar_flow_energy_type
    ! note: temperature in oC (tc) derived from 'global_auxvar_type'
    PetscReal :: TK            ! in Kevin, assuming ONE value for all 'fluids'
    PetscReal, pointer :: H(:) ! (nfluid+1) enthalpy in MJ/kmol , with 1 more dimension for when (liq_ or air_)fluid phase-changes to the solid or immobile
    PetscReal, pointer :: U(:) ! (nfluid+1) internal energy in MJ/kmol

    PetscReal, pointer :: D_H(:,:) !  (nfluid+1, nflowdof) MJ/kmol
    PetscReal, pointer :: D_U(:,:) !  (nfluid+1, nflowdof) MJ/kmol
  contains
   !..............
  end type auxvar_flow_energy_type

  public :: AuxVarFlowEnergyInit, &
            AuxVarFlowEnergyCopy, &
            AuxVarFlowEnergyStrip

contains

! ************************************************************************** !
subroutine AuxVarFlowEnergyInit(auxvar,option)
  ! 
  ! Initialise energy auxiliary variables
  ! 
  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019
  !
  ! 

  use Option_module

  implicit none

  class(auxvar_flow_energy_type) :: auxvar
  type(option_type) :: option
  PetscInt:: nfluid, nflowdof

  call AuxVarFlowInit(auxvar, option)

  auxvar%TK = auxvar%tc+TC2TK

  nullify(auxvar%H)
  nullify(auxvar%U)

  nfluid = option%flow%nfluid
  if(option%iflowmode > 0 .and. nfluid > 0) then
    allocate(auxvar%H(nfluid+1))
    auxvar%H = 0.d0
    allocate(auxvar%U(nfluid+1))
    auxvar%U = 0.d0
  endif

  nullify(auxvar%D_H)
  nullify(auxvar%D_U)
  nflowdof = option%nflowdof
  if (.not.option%flow%numerical_derivatives) then
    auxvar%has_derivs = PETSC_TRUE    ! not needed, but just in case

    allocate(auxvar%D_H(nfluid+1,nflowdof))
    auxvar%D_H = 0.d0
    allocate(auxvar%D_U(nfluid+1,nflowdof))
    auxvar%D_U = 0.d0
  endif


end subroutine AuxVarFlowEnergyInit

! ************************************************************************** !
subroutine AuxVarFlowEnergyCopy(auxvar, auxvar2)
  !
  ! dupliacate energy auxiliary variables
  !
  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019
  !  !

  implicit none

  class(auxvar_flow_energy_type) :: auxvar
  class(auxvar_flow_energy_type) :: auxvar2

  ! flow variables
  call AuxVarFlowCopy(auxvar, auxvar2)

  !
  auxvar2%TK = auxvar%TK
  if (associated(auxvar%H)) &
    auxvar2%H  = auxvar%H
  if (associated(auxvar%U)) &
    auxvar2%U  = auxvar%U

  auxvar2%has_derivs = auxvar%has_derivs

  if (associated(auxvar%D_H)) &
    auxvar2%D_H = auxvar%D_H
  if (associated(auxvar%D_U)) &
    auxvar2%D_U = auxvar%D_U


end subroutine AuxVarFlowEnergyCopy

! ************************************************************************** !

subroutine AuxVarFlowEnergyStrip(auxvar)
  ! 
  ! AuxVarFlowDestroy: Deallocates a toil_ims auxiliary object
  ! 
  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019
  !
  ! 
  use Utility_module, only : DeallocateArray

  implicit none

  class(auxvar_flow_energy_type) :: auxvar

  call DeallocateArray(auxvar%H)
  call DeallocateArray(auxvar%U)

  if (auxvar%has_derivs) then
    call DeallocateArray(auxvar%D_H)
    call DeallocateArray(auxvar%D_U)
  endif

  call AuxVarFlowStrip(auxvar)

end subroutine AuxVarFlowEnergyStrip

! ************************************************************************** !

end module AuxVars_Flow_Energy_module

