module Global_Aux_module
  ! variables globally shared (accessible) for all process models (PMs)

  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019
  !

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none
  
  private 

  ! NOTES on a few concenpts -
  !  - 'fluid(s)':
  !       flow matter or mixture, like liq. water or solution ('liq_fluid'), air('air_fluid'), oil('oil_fluid') etc.
  !       that can transport itself and matter (species) in it.
  !  - 'phase(s)':
  !       a matter in physical forms or states, usually as liq., gas, and solid. Maybe for one fluid or species in it.
  !  - 'flow species':
  !       the individual gas/liq/oil species, including primary fluid itself, in a fluid (either 'liq_fluid', 'air_fluid', 'oil_fluid').
  !  - 'dof':
  !       'Degree of Freedom', usually here means independent variables or states.

  type, public :: global_auxvar_type
    PetscReal          :: por               ! porosity of porous media for fluid(s) to flow
    PetscReal          :: tc                ! oC
    ! dim: constants: option%flow%nfluid
    PetscReal, pointer :: pres(:)           ! Pa
    PetscReal, pointer :: sat(:)            ! fraction
    PetscReal, pointer :: den(:)            ! kmol/m^3
    ! dim2: timesteps (currently 2) to be stored
    PetscReal, pointer :: por_store(:)      ! (timesteps)
    PetscReal, pointer :: tc_store(:)       !
    PetscReal, pointer :: pres_store(:,:)   ! (nfluid, timesteps)
    PetscReal, pointer :: sat_store(:,:)
    PetscReal, pointer :: den_store(:,:)

  end type global_auxvar_type
  
  ! the following data-type is for domain-wide and all processe-modes (including flows/RTs/etc.) uses
  type, public :: global_type
    PetscReal :: time_t, time_tpdt
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(global_auxvar_type), pointer :: auxvars(:)
    type(global_auxvar_type), pointer :: auxvars_bc(:)
    type(global_auxvar_type), pointer :: auxvars_ss(:)
  end type global_type
  
  interface GlobalAuxVarDestroy
    module procedure GlobalAuxVarSingleDestroy
    module procedure GlobalAuxVarArrayDestroy
  end interface GlobalAuxVarDestroy
  
  public :: GlobalAuxCreate, GlobalAuxDestroy, &
            GlobalAuxVarInit, GlobalAuxVarCopy, &
            GlobalAuxVarDestroy, GlobalAuxVarStrip

contains

! ************************************************************************** !

function GlobalAuxCreate()
  ! 
  ! Allocate and initialize auxiliary object
  ! 
  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019
  !
  ! 

  implicit none
  
  type(global_type), pointer :: GlobalAuxCreate
  
  type(global_type), pointer :: aux

  allocate(aux) 
  aux%time_t = 0.d0
  aux%time_tpdt = 0.d0
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)

  GlobalAuxCreate => aux
  
end function GlobalAuxCreate

! ************************************************************************** !

subroutine GlobalAuxVarInit(auxvar,option)
  ! 
  ! Initialize auxiliary object
  ! 
  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019
  !

  use Option_module

  implicit none
  
  class(global_auxvar_type) :: auxvar
  type(option_type) :: option

  PetscInt :: nfluid
  
  !
  auxvar%por = option%reference_porosity
  auxvar%tc  = 0.d0
  nullify(auxvar%pres)
  nullify(auxvar%sat)
  nullify(auxvar%den)

  nullify(auxvar%por_store)
  nullify(auxvar%tc_store)
  nullify(auxvar%pres_store)
  nullify(auxvar%sat_store)
  nullify(auxvar%den_store)

  nfluid = max(1,option%flow%nfluid)  ! at least liq. fluid exists, for purpose of transport with known velocity (data-driven)
  if(nfluid>0) then
    allocate(auxvar%pres(nfluid))
    auxvar%pres = 0.d0
    allocate(auxvar%sat(nfluid))
    auxvar%sat = 0.d0
    allocate(auxvar%den(nfluid))
    auxvar%den = 0.d0
  endif

  ! stored variables only for 2 timesteps, if flow mode on
  if (option%iflowmode > 0) then
    allocate(auxvar%por_store(TWO_INTEGER))
    auxvar%por_store = 0.d0

    allocate(auxvar%tc_store(TWO_INTEGER))
    auxvar%tc_store = 0.d0

    allocate(auxvar%pres_store(nfluid,TWO_INTEGER))
    auxvar%pres_store = option%reference_pressure

    allocate(auxvar%sat_store(nfluid,TWO_INTEGER))
    auxvar%sat_store = 0.d0

    allocate(auxvar%den_store(nfluid,TWO_INTEGER))
    auxvar%den_store = 0.d0
  endif
 
end subroutine GlobalAuxVarInit

! ************************************************************************** !

subroutine GlobalAuxVarCopy(auxvar,auxvar2)
  ! 
  ! Copies an auxiliary variable
  ! 
  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019
  !

  implicit none
  
  class(global_auxvar_type) :: auxvar, auxvar2

  auxvar2%por  = auxvar%por
  auxvar2%tc   = auxvar%tc
  auxvar2%pres = auxvar%pres
  auxvar2%sat  = auxvar%sat
  auxvar2%den  = auxvar%den

  if (associated(auxvar2%por_store)) then
    auxvar2%tc_store = auxvar%por_store
  endif
  if (associated(auxvar2%tc_store)) then
    auxvar2%tc_store = auxvar%tc_store
  endif
  if (associated(auxvar2%pres_store)) then
    auxvar2%pres_store = auxvar%pres_store  
  endif
  if (associated(auxvar2%den_store)) then
    auxvar2%den_store = auxvar%den_store  
  endif
  if (associated(auxvar2%sat_store)) then
    auxvar2%sat_store = auxvar%sat_store  
  endif

end subroutine GlobalAuxVarCopy

! ************************************************************************** !

subroutine GlobalAuxVarSingleDestroy(auxvar)
  ! 
  ! Deallocates a mode auxiliary object

  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019
  !


  implicit none

  class(global_auxvar_type), pointer :: auxvar
  
  if (associated(auxvar)) then
    call GlobalAuxVarStrip(auxvar)
    deallocate(auxvar)
  endif
  nullify(auxvar)

end subroutine GlobalAuxVarSingleDestroy

! ************************************************************************** !

subroutine GlobalAuxVarArrayDestroy(auxvars)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019
  !

  implicit none

  type(global_auxvar_type), pointer :: auxvars(:)
  
  PetscInt :: iaux
  
  if (associated(auxvars)) then
    do iaux = 1, size(auxvars)
      call GlobalAuxVarStrip(auxvars(iaux))
    enddo  
    deallocate(auxvars)
  endif
  nullify(auxvars)

end subroutine GlobalAuxVarArrayDestroy

! ************************************************************************** !

subroutine GlobalAuxVarStrip(auxvar)
  ! 
  ! Deallocates all members of single auxiliary object
  ! 
  ! Author: Fengming Yuan @CCSI/ORNL
  ! Date: 11/07/2019
  !

  use Utility_module, only: DeallocateArray

  implicit none

  type(global_auxvar_type) :: auxvar
  
  call DeallocateArray(auxvar%por_store)
  call DeallocateArray(auxvar%tc_store)
  call DeallocateArray(auxvar%pres)
  call DeallocateArray(auxvar%sat)
  call DeallocateArray(auxvar%den)
  call DeallocateArray(auxvar%pres_store)
  call DeallocateArray(auxvar%sat_store)
  call DeallocateArray(auxvar%den_store)

end subroutine GlobalAuxVarStrip

! ************************************************************************** !

subroutine GlobalAuxDestroy(aux)
  ! 
  ! Deallocates a mode auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  ! 

  implicit none

  type(global_type), pointer :: aux
  
  if (.not.associated(aux)) return
  
  call GlobalAuxVarDestroy(aux%auxvars)
  call GlobalAuxVarDestroy(aux%auxvars_bc)
  call GlobalAuxVarDestroy(aux%auxvars_ss)
  
  deallocate(aux)
  nullify(aux)
  
end subroutine GlobalAuxDestroy

end module Global_Aux_module
