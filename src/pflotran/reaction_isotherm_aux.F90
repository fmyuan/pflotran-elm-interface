module Reaction_Isotherm_Aux_module
  
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  !kd units
  PetscInt, parameter, public :: KD_UNIT_KG_M3_BULK = 0
  PetscInt, parameter, public :: KD_UNIT_MLW_GSOIL = 1

  type, public :: isotherm_linklist_type
    PetscInt :: id
    PetscInt :: itype
    character(len=MAXWORDLENGTH) :: species_name
    character(len=MAXWORDLENGTH) :: kd_mineral_name
    PetscReal :: Kd
    PetscReal :: Langmuir_B
    PetscReal :: Freundlich_n
    type(isotherm_linklist_type), pointer :: next
  end type isotherm_linklist_type

  type, public :: isotherm_rxn_type
    PetscReal, pointer :: eqisothermcoefficient(:)
    PetscReal, pointer :: eqisothermlangmuirb(:)
    PetscReal, pointer :: eqisothermfreundlichn(:)
  end type isotherm_rxn_type

  type, public :: isotherm_type
    type(isotherm_linklist_type), pointer :: isotherm_list
    type(isotherm_linklist_type), pointer :: multicontinuum_isotherm_list
    type(isotherm_rxn_type), pointer :: isotherm_rxn
    type(isotherm_rxn_type), pointer :: multicontinuum_isotherm_rxn
    PetscInt :: kd_unit
    PetscInt, pointer :: eqisothermtype(:)
    PetscInt, pointer :: eqkdspecid(:)
    PetscInt, pointer :: eqkdmineral(:)
    PetscInt :: neqkdrxn
  end type isotherm_type

  public :: IsothermCreate, &
            IsothermRxnCreate, &
            IsothermDestroy

contains
  
! ************************************************************************** !
 
function IsothermCreate()
  ! 
  ! Allocate and initialize isotherm reaction object
  ! 

  implicit none

  type(isotherm_type), pointer :: IsothermCreate

  type(isotherm_type), pointer :: isotherm

  allocate(isotherm)

  nullify(isotherm%isotherm_list)
  nullify(isotherm%multicontinuum_isotherm_list)
  nullify(isotherm%eqisothermtype)
  nullify(isotherm%eqkdspecid)
  nullify(isotherm%eqkdmineral)
  nullify(isotherm%isotherm_rxn)
  nullify(isotherm%multicontinuum_isotherm_rxn)

  isotherm%kd_unit = -999
  isotherm%neqkdrxn = 0

  IsothermCreate => isotherm

end function IsothermCreate

! ************************************************************************** !

function IsothermRxnCreate()

  ! Allocate and initialize an isotherm sorption reaction
  ! 
  ! 
  implicit none
  
  type(isotherm_linklist_type), pointer :: IsothermRxnCreate
  
  type(isotherm_linklist_type), pointer :: rxn

  allocate(rxn)
  rxn%id = 0
  rxn%itype = 0
  rxn%species_name = ''
  rxn%kd_mineral_name = ''
  rxn%Kd = 0.d0
  rxn%Langmuir_B = 0.d0
  rxn%Freundlich_n = 0.d0
  nullify(rxn%next)

  IsothermRxnCreate => rxn

end function IsothermRxnCreate

! ************************************************************************** !

subroutine IsothermRxnDestroy(rxn)

  ! 
  ! Deallocates an isotherm reaction
  ! 


  implicit none

  type(isotherm_linklist_type), pointer :: rxn

  if (.not.associated(rxn)) return

  deallocate(rxn)  
  nullify(rxn)

end subroutine IsothermRxnDestroy

! ************************************************************************** !

subroutine IsothermDestroy(isotherm,option)

  ! 
  ! Deallocates an isotherm object
  ! 

  use Utility_module, only: DeallocateArray
  use Option_module
  
  implicit none

  type(isotherm_type) :: isotherm
  type(option_type) :: option
  type(isotherm_linklist_type), pointer :: kd_rxn, prev_kd_rxn

  kd_rxn => isotherm%isotherm_list

  do
    if (.not.associated(kd_rxn)) exit
    prev_kd_rxn => kd_rxn
    kd_rxn => kd_rxn%next
    call IsothermRxnDestroy(prev_kd_rxn)
  enddo
  nullify(isotherm%isotherm_list)

  ! secondary continuum
  if (option%use_mc) then
    kd_rxn => isotherm%multicontinuum_isotherm_list
    do
      if (.not.associated(kd_rxn)) exit
      prev_kd_rxn => kd_rxn
      kd_rxn => kd_rxn%next
      call IsothermRxnDestroy(prev_kd_rxn)
    enddo
    nullify(isotherm%multicontinuum_isotherm_list)
  endif
 
  if (associated(isotherm%isotherm_rxn)) then
    call DeallocateArray(isotherm%isotherm_rxn%eqisothermcoefficient)
    call DeallocateArray(isotherm%isotherm_rxn%eqisothermlangmuirb)
    call DeallocateArray(isotherm%isotherm_rxn%eqisothermfreundlichn)
  endif

  nullify(isotherm%isotherm_rxn)

  if (associated(isotherm%multicontinuum_isotherm_rxn)) then
    call DeallocateArray(isotherm%multicontinuum_isotherm_rxn%eqisothermcoefficient)
    call DeallocateArray(isotherm%multicontinuum_isotherm_rxn%eqisothermlangmuirb)
    call DeallocateArray(isotherm%multicontinuum_isotherm_rxn%eqisothermfreundlichn)
  endif

  nullify(isotherm%multicontinuum_isotherm_rxn)

  call DeallocateArray(isotherm%eqisothermtype)
  call DeallocateArray(isotherm%eqkdspecid)
  call DeallocateArray(isotherm%eqkdmineral)
  
end subroutine IsothermDestroy

end module Reaction_Isotherm_Aux_module
