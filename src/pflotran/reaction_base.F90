module Reaction_Base_module
  
#include "petsc/finclude/petscsys.h"
  use petscsys

  private 

  type, public :: reaction_base_type
    ! flag for solving for the change in the log of the concentration
    PetscBool :: use_log_formulation 
    PetscBool :: equilibrate_at_each_cell
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

end subroutine ReactionBaseInit

! ************************************************************************** !

subroutine ReactionBaseStrip(reaction_base)
  ! 
  ! Deallocates members of a base reaction object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/21/19
  ! 
  implicit none

  class(reaction_base_type) :: reaction_base

end subroutine ReactionBaseStrip

end module Reaction_Base_module
