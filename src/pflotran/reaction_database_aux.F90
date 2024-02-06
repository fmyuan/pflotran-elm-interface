module Reaction_Database_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Reaction_Equation_module

  implicit none

  private

  type, public :: database_rxn_type
    type(reaction_equation_type), pointer :: reaction_equation
    PetscReal, pointer :: logK(:)
    PetscReal, pointer :: logKCoeff_hpt(:)
  end type database_rxn_type

  interface ReactionDBCreateRxn
    module procedure ReactionDBCreateRxn1
    module procedure ReactionDBCreateRxn2
  end interface ReactionDBCreateRxn

  public :: ReactionDBAlignSpeciesInRxn, &
            ReactionDBSubSpecInRxn, &
            ReactionDBCreateRxn, &
            ReactionDBCheckLegitLogKs, &
            ReactionDBDestroyRxn

contains

! ************************************************************************** !

function ReactionDBCreateRxn1()
  !
  ! Allocate and initialize an equilibrium reaction
  !
  ! Author: Glenn Hammond
  ! Date: 09/01/08
  !

  implicit none

  type(database_rxn_type), pointer :: ReactionDBCreateRxn1

  type(database_rxn_type), pointer :: dbaserxn

  allocate(dbaserxn)
  nullify(dbaserxn%reaction_equation)
  nullify(dbaserxn%logK)
  nullify(dbaserxn%logKCoeff_hpt)

  ReactionDBCreateRxn1 => dbaserxn

end function ReactionDBCreateRxn1

! ************************************************************************** !

function ReactionDBCreateRxn2(num_species,num_logKs)
  !
  ! Allocate and initialize an equilibrium reaction
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !

  implicit none

  PetscInt :: num_species
  PetscInt :: num_logKs

  type(database_rxn_type), pointer :: ReactionDBCreateRxn2

  type(database_rxn_type), pointer :: dbaserxn

  dbaserxn => ReactionDBCreateRxn()
  dbaserxn%reaction_equation => ReactionEquationCreate(num_species)
  allocate(dbaserxn%logK(num_logKs))
  dbaserxn%logK = UNINITIALIZED_DOUBLE

  ReactionDBCreateRxn2 => dbaserxn

end function ReactionDBCreateRxn2

! ************************************************************************** !

subroutine ReactionDBAlignSpeciesInRxn(num_basis_species,basis_names, &
                                       reaction_equation,species_name,option)
  !
  ! Aligns the ordering of species in reaction with
  ! the current basis
  !
  ! Author: Glenn Hammond
  ! Date: 10/07/08
  !
  use Option_module
  use String_module
  use Utility_module, only : DeallocateArray

  implicit none

  PetscInt :: num_basis_species
  character(len=MAXWORDLENGTH) :: basis_names(num_basis_species), species_name
  type(reaction_equation_type) :: reaction_equation
  type(option_type) :: option

  PetscInt :: i_rxn_species
  PetscInt :: i_basis_species
  PetscInt :: num_species
  PetscReal :: stoich_new(num_basis_species)
  PetscBool :: found

  ! reallocate specid to proper size
  call DeallocateArray(reaction_equation%specid)
  num_species = size(reaction_equation%spec_name)
  reaction_equation%nspec = num_species
  allocate(reaction_equation%specid(num_species))
  reaction_equation%specid(:) = UNINITIALIZED_INTEGER

  stoich_new = 0.d0
  do i_rxn_species = 1, reaction_equation%nspec
    found = PETSC_FALSE
    do i_basis_species = 1, num_basis_species
      if (StringCompare(reaction_equation%spec_name(i_rxn_species), &
                        basis_names(i_basis_species), &
                        MAXWORDLENGTH)) then
        stoich_new(i_basis_species) = reaction_equation%stoich(i_rxn_species)
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      option%io_buffer = &
        trim(reaction_equation%spec_name(i_rxn_species)) // &
        ' not found in basis (ReactionDBAlignSpeciesInRxn) for species ' // &
        trim(species_name)
      call PrintErrMsg(option)
    endif
  enddo

  ! zero everthing out
  reaction_equation%spec_name = ''
  reaction_equation%stoich = 0.d0
  reaction_equation%specid(1:num_species) = 0

  ! fill in
  i_rxn_species = 0
  do i_basis_species = 1, num_basis_species
    if (dabs(stoich_new(i_basis_species)) > 1.d-40) then
      i_rxn_species = i_rxn_species + 1
      reaction_equation%spec_name(i_rxn_species) = &
        basis_names(i_basis_species)
      reaction_equation%stoich(i_rxn_species) = stoich_new(i_basis_species)
      reaction_equation%specid(i_rxn_species) = i_basis_species
    endif
  enddo

  if (i_rxn_species /= reaction_equation%nspec) then
    write(option%io_buffer,*) &
                   'Number of reaction species does not match original:', &
                    i_rxn_species, reaction_equation%nspec
    call PrintErrMsg(option)
  endif

end subroutine ReactionDBAlignSpeciesInRxn

! ************************************************************************** !

subroutine ReactionDBSubSpecInRxn(name1,dbaserxn1,dbaserxn2)
  !
  ! Swaps out a chemical species in a chemical reaction, replacing it with
  ! the species in a secondary reaction (swaps 1 into 2)
  !
  ! Author: Glenn Hammond
  ! Date: 10/06/08
  !
  use String_module

  implicit none

  character(len=MAXWORDLENGTH) :: name1
  type(database_rxn_type) :: dbaserxn1
  type(database_rxn_type) :: dbaserxn2

  PetscReal :: scale

  call ReactionEquationSubSpecInRxn(name1,dbaserxn1%reaction_equation, &
                                    dbaserxn2%reaction_equation,scale)
  dbaserxn2%logK = dbaserxn2%logK + scale*dbaserxn1%logK

end subroutine ReactionDBSubSpecInRxn

! ************************************************************************** !

function ReactionDBCheckLegitLogKs(dbaserxn,species_name,temperatures, &
                                   option)
  !
  ! Checks whether legitimate log Ks exist for
  ! all database temperatures if running
  ! non-isothermal
  !
  ! Author: Glenn Hammond
  ! Date: 01/07/13
  !
  use Option_module
  use Utility_module, only : Equal

  implicit none

  type(database_rxn_type), pointer :: dbaserxn
  character(len=MAXWORDLENGTH) :: species_name
  PetscReal :: temperatures(:)
  type(option_type) :: option

  PetscBool :: ReactionDBCheckLegitLogKs

  PetscInt :: itemp
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word

  ReactionDBCheckLegitLogKs = PETSC_TRUE

  if (.not.associated(dbaserxn)) return
  if (option%use_isothermal .and. &
      Equal(option%flow%reference_temperature,25.d0)) return

  string = ''
  do itemp = 1, size(dbaserxn%logK)
    if (Equal(dabs(dbaserxn%logK(itemp)),500.d0)) then
      write(word,'(f5.1)') temperatures(itemp)
      string = trim(string) // ' ' // word
      ReactionDBCheckLegitLogKs = PETSC_FALSE
    endif
  enddo

  if (.not.ReactionDBCheckLegitLogKs) then
    option%io_buffer = ' ERROR: Undefined log Ks for temperatures (' // &
                       trim(adjustl(string)) // ') for species "' // &
                       trim(species_name) // '" in database.'
    call PrintMsg(option)
  endif

end function ReactionDBCheckLegitLogKs

! ************************************************************************** !

subroutine ReactionDBDestroyRxn(dbaserxn)
  !
  ! Deallocates a database reaction
  !
  ! Author: Glenn Hammond
  ! Date: 05/29/08
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(database_rxn_type), pointer :: dbaserxn

  if (.not.associated(dbaserxn)) return

  call ReactionEquationDestroy(dbaserxn%reaction_equation)
  call DeallocateArray(dbaserxn%logK)
  call DeallocateArray(dbaserxn%logKCoeff_hpt)

  deallocate(dbaserxn)
  nullify(dbaserxn)

end subroutine ReactionDBDestroyRxn

end module Reaction_Database_Aux_module
