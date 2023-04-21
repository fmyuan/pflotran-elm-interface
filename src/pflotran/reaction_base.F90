module Reaction_Base_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  private

  type, public :: reaction_base_type
    ! flag for solving for the change in the log of the concentration
    PetscBool :: use_log_formulation
    PetscBool :: equilibrate_at_each_cell
    PetscInt, pointer :: print_cells(:)
  end type reaction_base_type

  public :: ReactionBaseInit, &
            ReactionBaseStrip

contains

! ************************************************************************** !

subroutine ReactionBaseInit(reaction_base)
  !
  ! Initialize base reaction object
  !
  ! Author: Glenn Hammond
  ! Date: 10/21/19
  !
  implicit none

  class(reaction_base_type) :: reaction_base

  reaction_base%use_log_formulation = PETSC_FALSE
  reaction_base%equilibrate_at_each_cell = PETSC_FALSE
  nullify(reaction_base%print_cells)

end subroutine ReactionBaseInit

! ************************************************************************** !

subroutine ReactionBaseStrip(reaction_base)
  !
  ! Deallocates members of a base reaction object
  !
  ! Author: Glenn Hammond
  ! Date: 10/21/19
  !
  use Utility_module

  implicit none

  class(reaction_base_type) :: reaction_base

  call DeallocateArray(reaction_base%print_cells)

end subroutine ReactionBaseStrip

end module Reaction_Base_module
