module Reaction_Isotherm_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  !kd units
  PetscInt, parameter, public :: KD_UNIT_KG_M3_BULK = 0
  PetscInt, parameter, public :: KD_UNIT_MLW_GSOIL = 1

  type, public :: isotherm_link_type
    PetscInt :: id
    PetscInt :: itype
    character(len=MAXWORDLENGTH) :: species_name
    character(len=MAXWORDLENGTH) :: kd_mineral_name
    PetscReal :: Kd
    PetscReal :: Langmuir_B
    PetscReal :: Freundlich_n
    type(isotherm_link_type), pointer :: next
  end type isotherm_link_type

  type, public :: isotherm_rxn_type
    PetscReal, pointer :: eqisothermcoeff(:)
    PetscReal, pointer :: eqisothermlangmuirb(:)
    PetscReal, pointer :: eqisothermfreundlichn(:)
  end type isotherm_rxn_type

  type, public :: isotherm_type
    type(isotherm_link_type), pointer :: isotherm_list
    type(isotherm_link_type), pointer :: multicontinuum_isotherm_list
    type(isotherm_rxn_type), pointer :: isotherm_rxn
    type(isotherm_rxn_type), pointer :: multicontinuum_isotherm_rxn
    PetscInt :: ikd_units
    PetscInt, pointer :: eqisothermtype(:)
    PetscInt, pointer :: eqkdspecid(:)
    PetscInt, pointer :: eqkdmineral(:)
    PetscInt :: neqkdrxn
  end type isotherm_type

  public :: ReactionIsothermCreateObject, &
            ReactionIsothermCreateLink, &
            ReactionIsothermCreateRxn, &
            ReactionIsothermDestroyObject

contains

! ************************************************************************** !

function ReactionIsothermCreateObject()
  !
  ! Allocate and initialize isotherm reaction object
  !

  implicit none

  type(isotherm_type), pointer :: ReactionIsothermCreateObject

  type(isotherm_type), pointer :: isotherm

  allocate(isotherm)

  nullify(isotherm%isotherm_list)
  nullify(isotherm%multicontinuum_isotherm_list)
  nullify(isotherm%eqisothermtype)
  nullify(isotherm%eqkdspecid)
  nullify(isotherm%eqkdmineral)
  nullify(isotherm%isotherm_rxn)
  nullify(isotherm%multicontinuum_isotherm_rxn)

  isotherm%ikd_units = UNINITIALIZED_INTEGER
  isotherm%neqkdrxn = 0

  ReactionIsothermCreateObject => isotherm

end function ReactionIsothermCreateObject

! ************************************************************************** !

function ReactionIsothermCreateLink()

  ! Allocate and initialize an isotherm sorption reaction
  !
  !
  implicit none

  type(isotherm_link_type), pointer :: ReactionIsothermCreateLink

  type(isotherm_link_type), pointer :: rxn

  allocate(rxn)
  rxn%id = 0
  rxn%itype = 0
  rxn%species_name = ''
  rxn%kd_mineral_name = ''
  rxn%Kd = 0.d0
  rxn%Langmuir_B = 0.d0
  rxn%Freundlich_n = 0.d0
  nullify(rxn%next)

  ReactionIsothermCreateLink => rxn

end function ReactionIsothermCreateLink

! ************************************************************************** !

subroutine ReactionIsothermCreateRxn(isotherm_rxn, isotherm)

  implicit none

  type(isotherm_type), pointer :: isotherm

  type(isotherm_rxn_type), pointer :: isotherm_rxn

  allocate(isotherm_rxn)
  ! allocate arrays
  allocate(isotherm_rxn%eqisothermcoeff(isotherm%neqkdrxn))
  isotherm_rxn%eqisothermcoeff = 0.d0
  allocate(isotherm_rxn%eqisothermlangmuirb(isotherm%neqkdrxn))
  isotherm_rxn%eqisothermlangmuirb = 0.d0
  allocate(isotherm_rxn%eqisothermfreundlichn(isotherm%neqkdrxn))
  isotherm_rxn%eqisothermfreundlichn = 0.d0

end subroutine ReactionIsothermCreateRxn

! ************************************************************************** !

subroutine IsothermRxnLinkDestroy(link)
  !
  ! Deallocates an isotherm reaction
  !
  implicit none

  type(isotherm_link_type), pointer :: link

  if (.not.associated(link)) return

  deallocate(link)
  nullify(link)

end subroutine IsothermRxnLinkDestroy

! ************************************************************************** !

subroutine IsothermRxnDestroy(isotherm_rxn)
  !
  ! Deallocates an isotherm reaction
  !
  use Utility_module, only: DeallocateArray

  implicit none

  type(isotherm_rxn_type), pointer :: isotherm_rxn

  if (.not.associated(isotherm_rxn)) return

  call DeallocateArray(isotherm_rxn%eqisothermcoeff)
  call DeallocateArray(isotherm_rxn%eqisothermlangmuirb)
  call DeallocateArray(isotherm_rxn%eqisothermfreundlichn)

  deallocate(isotherm_rxn)
  nullify(isotherm_rxn)

end subroutine IsothermRxnDestroy

! ************************************************************************** !

subroutine ReactionIsothermDestroyObject(isotherm,option)
  !
  ! Deallocates an isotherm object
  !

  use Utility_module, only: DeallocateArray
  use Option_module

  implicit none

  type(isotherm_type), pointer :: isotherm
  type(option_type) :: option
  type(isotherm_link_type), pointer :: isotherm_rxn_link, &
                                       prev_isotherm_rxn_link

  isotherm_rxn_link => isotherm%isotherm_list

  do
    if (.not.associated(isotherm_rxn_link)) exit
    prev_isotherm_rxn_link => isotherm_rxn_link
    isotherm_rxn_link => isotherm_rxn_link%next
    call IsothermRxnLinkDestroy(prev_isotherm_rxn_link)
  enddo
  nullify(isotherm%isotherm_list)

  ! secondary continuum
  if (option%use_sc) then
    isotherm_rxn_link => isotherm%multicontinuum_isotherm_list
    do
      if (.not.associated(isotherm_rxn_link)) exit
      prev_isotherm_rxn_link => isotherm_rxn_link
      isotherm_rxn_link => isotherm_rxn_link%next
      call IsothermRxnLinkDestroy(prev_isotherm_rxn_link)
    enddo
    nullify(isotherm%multicontinuum_isotherm_list)
  endif

  call IsothermRxnDestroy(isotherm%isotherm_rxn)
  call IsothermRxnDestroy(isotherm%multicontinuum_isotherm_rxn)

  call DeallocateArray(isotherm%eqisothermtype)
  call DeallocateArray(isotherm%eqkdspecid)
  call DeallocateArray(isotherm%eqkdmineral)

  deallocate(isotherm)
  nullify(isotherm)

end subroutine ReactionIsothermDestroyObject

end module Reaction_Isotherm_Aux_module
