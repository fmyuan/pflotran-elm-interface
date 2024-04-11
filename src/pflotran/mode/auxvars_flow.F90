module AuxVars_Flow_module
  ! variables globally shared (accessible) Extended for flow process model (PM)

  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019
  !

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  use Option_Flow_module
  use Global_Aux_module

  implicit none

  PetscReal,public :: flow_aux_debug_tol = 1.d-1
  PetscReal,public :: flow_aux_debug_reltol = 1.d-1
  PetscBool,public :: flow_aux_use_GP = PETSC_FALSE

  private

  type, public, extends(global_auxvar_type) :: auxvar_flow_type

    PetscReal          :: pc                   ! liq. fluid capillary pressure only in (negative) Pa (specific variable very useful)

    ! note: here 'sat','den','pres' are derived from 'global_auxvar_type'
    PetscReal, pointer :: mobility(:)          ! (nfluid): relative permissivity/dynamic viscosity
    PetscReal, pointer :: viscosity(:)         ! (nfluid): dynamic viscosity

    PetscReal, pointer :: den_kg(:)            ! (nfluid+1) kg/m^3, with 1 more dimension for when (liq_ or air_)fluid phase-changes to the solid or immobile

    ! dim: option_flow%nspecgas, for all dissolved-gas solutes in a liq. fuild solution (when nfluid>=2, disgassing-disolving process may be needed)
    PetscReal, pointer :: fugacoeff(:)         ! fugacity coefficent in Henry's Law
    PetscReal, pointer :: xmass(:)             ! for convertion from molarity (moles/L solution) to molality(moles/kg solvent)
    PetscReal, pointer :: liq_molarity(:)      ! (nflowgas) aq. gas species in liq. fluid in moles/L solution
    PetscReal, pointer :: gas_molefraction(:)  ! (nflowgas) gas species in gas fluid

    ! dim: (nfluid+1, option_flow%nspecflow). NOTE that this is for fluid species in mixture or itself, NOT transportants in it.
    PetscReal, pointer :: mass_balance(:,:)       ! kmol
    PetscReal, pointer :: mass_balance_delta(:,:) ! kmol/timestep

    ! derivatives
    PetscBool :: has_derivs

    PetscReal, pointer :: D_por(:)             ! (nflowdof)
    PetscReal, pointer :: D_pc(:)              ! (nflowdof)

    PetscReal, pointer :: D_sat(:,:)           ! (nfluid, nflowdof)
    PetscReal, pointer :: D_den(:,:)           ! (nfluid, nflowdof) kmol/m^3 phase
    PetscReal, pointer :: D_pres(:,:)          ! (nfluid, nflowdof)
    PetscReal, pointer :: D_mobility(:,:)      ! (nfluid, nflowdof) relative perm./visc.
    PetscReal, pointer :: D_viscosity(:,:)     ! (nfluid, nflowdof) dynamic viscosity

    PetscReal, pointer :: D_den_kg(:,:)        ! (nfluid+1, nflowdof) kg/m^3 phase

    PetscReal, pointer :: D_liqmolarity(:,:)      ! (nspecgas, nflowdof) aq. gas in liq. fluid
    PetscReal, pointer :: D_gasmolefraction(:,:)  ! (nspecgas, nflowdof) gas in gas fluid

  contains
    !
  end type auxvar_flow_type

  public :: AuxVarFlowInit, AuxVarFlowCopy, AuxVarFlowStrip

contains

! ************************************************************************** !
subroutine AuxVarFlowInit(auxvar,option_flow)
  !
  ! initialize flow auxiliary variables
  !
  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019
  !

  use Option_module

  implicit none

  class(auxvar_flow_type) :: auxvar
  type(option_type) :: option
  type(option_flow_type) :: option_flow

  PetscInt:: nflowgas, nfluid, nflowspec, nflowdof

  ! base auxiliary variables
  call GlobalAuxVarInit(auxvar, option)

  ! specific flow relevant variables
  auxvar%pc = 0.0d0

  nullify(auxvar%mobility)
  nullify(auxvar%viscosity)
  nullify(auxvar%den_kg)

  nfluid = option%nfluid
  if(option%iflowmode > 0 .and. nfluid > 0) then
    allocate(auxvar%mobility(nfluid))
    auxvar%mobility = 0.d0
    allocate(auxvar%viscosity(nfluid))
    auxvar%viscosity = 0.d0
    allocate(auxvar%den_kg(nfluid+1))
    auxvar%den_kg = 0.d0
  endif

  nullify(auxvar%fugacoeff)
  nullify(auxvar%xmass)
  nullify(auxvar%liq_molarity)
  nullify(auxvar%gas_molefraction)

  nflowgas   = option_flow%nspecgas
  if (option%iflowmode > 0 .and. nflowgas > 0) then
      allocate(auxvar%fugacoeff(nflowgas))
      auxvar%fugacoeff = 1.d0
      allocate(auxvar%xmass(nflowgas))
      auxvar%xmass = 1.d0
      allocate(auxvar%liq_molarity(nflowgas))
      auxvar%liq_molarity = 0.d0
      allocate(auxvar%gas_molefraction(nflowgas))
      auxvar%gas_molefraction = 0.d0
  endif

  nullify(auxvar%mass_balance)
  nullify(auxvar%mass_balance_delta)
  ! for fluid(s) only, including its phase-changing mass (NOT transportants in it)
  nflowspec  = option_flow%nspecflow

  if (option%iflag/=0 .and. option%compute_mass_balance_new .and. nflowspec>0) then
    allocate(auxvar%mass_balance(nfluid+1, nflowspec))
    auxvar%mass_balance = 0.d0
    allocate(auxvar%mass_balance_delta(nfluid+1, nflowspec))
    auxvar%mass_balance_delta = 0.d0
  endif

  !
  nullify(auxvar%D_pc)
  nullify(auxvar%D_por)
  nullify(auxvar%D_pres)
  nullify(auxvar%D_sat)
  nullify(auxvar%D_den)
  nullify(auxvar%D_mobility)
  nullify(auxvar%D_viscosity)
  nullify(auxvar%D_den_kg)
  nullify(auxvar%D_liqmolarity)
  nullify(auxvar%D_gasmolefraction)

  nflowdof = option%nflowdof
  auxvar%has_derivs = PETSC_FALSE
  if (.not.option%numerical_derivatives) then
    auxvar%has_derivs = PETSC_TRUE

    allocate(auxvar%D_pc(nflowdof))
    auxvar%D_pc = 0.d0
    allocate(auxvar%D_por(nflowdof))
    auxvar%D_por = 0.d0

    allocate(auxvar%D_pres(nfluid,nflowdof))
    auxvar%D_pres = 0.d0
    allocate(auxvar%D_sat(nfluid,nflowdof))
    auxvar%D_sat = 0.d0
    allocate(auxvar%D_den(nfluid,nflowdof))
    auxvar%D_den = 0.d0
    allocate(auxvar%D_mobility(nfluid,nflowdof))
    auxvar%D_mobility = 0.d0
    allocate(auxvar%D_viscosity(nfluid,nflowdof))
    auxvar%D_viscosity = 0.d0

    allocate(auxvar%D_den_kg(nfluid+1,nflowdof))
    auxvar%D_den_kg = 0.d0

    allocate(auxvar%D_liqmolarity(nflowgas,nflowdof))
    auxvar%D_liqmolarity = 0.d0
    allocate(auxvar%D_gasmolefraction(nflowgas,nflowdof))
    auxvar%D_gasmolefraction = 0.d0

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

  ! base variables
  call GlobalAuxVarCopy(auxvar,auxvar2)

  ! specific variables for flow
  auxvar2%pc  = auxvar%pc
  if (associated(auxvar%sat)) &
    auxvar2%sat = auxvar%sat
  if (associated(auxvar%den)) &
    auxvar2%den = auxvar%den
  if (associated(auxvar%den_kg)) &
    auxvar2%den_kg = auxvar%den_kg
  if (associated(auxvar%mobility)) &
    auxvar2%mobility = auxvar%mobility
  if (associated(auxvar%viscosity)) &
    auxvar2%viscosity = auxvar%viscosity

  if (associated(auxvar%fugacoeff)) &
    auxvar2%fugacoeff = auxvar%fugacoeff
  if (associated(auxvar%xmass)) &
    auxvar2%xmass = auxvar%xmass
  if (associated(auxvar%liq_molarity)) &
    auxvar2%liq_molarity = auxvar%liq_molarity
  if (associated(auxvar%gas_molefraction)) &
    auxvar2%gas_molefraction = auxvar%gas_molefraction

  if (associated(auxvar%mass_balance) .and. &
      associated(auxvar2%mass_balance)) then
    auxvar2%mass_balance = auxvar%mass_balance
   endif
  if (associated(auxvar%mass_balance_delta) .and. &
      associated(auxvar2%mass_balance_delta)) then
    auxvar2%mass_balance_delta = auxvar%mass_balance_delta
  endif

  if (associated(auxvar%D_por)) &
    auxvar2%D_por  = auxvar%D_por
  if (associated(auxvar%D_pc)) &
    auxvar2%D_pc   = auxvar%D_pc

  if (associated(auxvar%D_pres)) &
    auxvar2%D_pres = auxvar%D_pres
  if (associated(auxvar%D_sat)) &
    auxvar2%D_sat  = auxvar%D_sat
  if (associated(auxvar%D_den)) &
    auxvar2%D_den  = auxvar%D_den
  if (associated(auxvar%D_den_kg)) &
    auxvar2%D_den_kg   = auxvar%D_den_kg
  if (associated(auxvar%D_mobility)) &
    auxvar2%D_mobility = auxvar%D_mobility
  if (associated(auxvar%D_viscosity)) &
    auxvar2%D_viscosity= auxvar%D_viscosity
  if (associated(auxvar%D_liqmolarity)) &
    auxvar2%D_viscosity= auxvar%D_liqmolarity
  if (associated(auxvar%D_gasmolefraction)) &
    auxvar2%D_viscosity= auxvar%D_gasmolefraction

end subroutine AuxVarFlowCopy

! ************************************************************************** !

subroutine AuxVarFlowStrip(auxvar)
  !
  ! AuxVarFlowDestroy: Deallocates a flow auxiliary object
  !
  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019
  !
  use Utility_module, only : DeallocateArray

  implicit none

  class(auxvar_flow_type) :: auxvar

  call DeallocateArray(auxvar%pres)
  call DeallocateArray(auxvar%sat)
  call DeallocateArray(auxvar%den)
  call DeallocateArray(auxvar%den_kg)
  call DeallocateArray(auxvar%mobility)
  call DeallocateArray(auxvar%viscosity)

  if (associated(auxvar%fugacoeff)) &
    call DeallocateArray(auxvar%fugacoeff)
  if (associated(auxvar%xmass)) &
    call DeallocateArray(auxvar%xmass)
  if (associated(auxvar%liq_molarity)) &
    call DeallocateArray(auxvar%liq_molarity)
  if (associated(auxvar%gas_molefraction)) &
   call DeallocateArray(auxvar%gas_molefraction)

  if (associated(auxvar%mass_balance)) &
   call DeallocateArray(auxvar%mass_balance)
  if (associated(auxvar%mass_balance_delta)) &
   call DeallocateArray(auxvar%mass_balance_delta)

  if (auxvar%has_derivs) then
    call DeallocateArray(auxvar%D_pres)
    call DeallocateArray(auxvar%D_sat)
    call DeallocateArray(auxvar%D_pc)
    call DeallocateArray(auxvar%D_den)
    call DeallocateArray(auxvar%D_den_kg)
    call DeallocateArray(auxvar%D_mobility)
    call DeallocateArray(auxvar%D_viscosity)
    call DeallocateArray(auxvar%D_por)
  endif

end subroutine AuxVarFlowStrip
! ************************************************************************** !

end module AuxVars_Flow_module
