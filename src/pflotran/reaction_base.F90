module Reaction_Base_module
  
#include "petsc/finclude/petscsys.h"
  use petscsys

  private 

  type, public :: reaction_base_type
  end type reaction_base_type

  public :: ReactionBaseInit, &
            ReactionBaseStrip
             
contains

! ************************************************************************** !

subroutine ReactionBaseInit(reaction)
  ! 
  ! Initialize base reaction object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/21/19
  ! 
  implicit none
  
  class(reaction_base_type) :: reaction

end subroutine ReactionBaseInit

! ************************************************************************** !

subroutine ReactionBaseStrip(reaction)
  ! 
  ! Deallocates members of a base reaction object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/21/19
  ! 
  implicit none

  class(reaction_base_type) :: reaction

end subroutine ReactionBaseStrip

end module Reaction_Base_module
