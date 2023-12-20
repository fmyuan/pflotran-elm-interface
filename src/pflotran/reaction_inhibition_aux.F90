module Reaction_Inhibition_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module

  implicit none

  private

  ! inhibition parameters
  PetscInt, parameter, public :: INHIBITION_MONOD = 1
  PetscInt, parameter, public :: INHIBITION_THRESHOLD = 2
  PetscInt, parameter, public :: INHIBITION_SMOOTHSTEP = 3

  type, public :: inhibition_type
    PetscInt :: id
    PetscInt :: itype
    character(len=MAXWORDLENGTH) :: species_name
    PetscReal :: inhibition_constant
    PetscReal :: inhibition_constant2
    type(inhibition_type), pointer :: next
  end type inhibition_type

  interface ReactionInhibitionThreshold
    module procedure ReactionInhibitionThreshold1
    module procedure ReactionInhibitionThreshold2
  end interface

  interface ReactionInhibitionSmoothstep
    module procedure ReactionInhibitionSmoothstep1
    module procedure ReactionInhibitionSmoothstep2
  end interface

  public :: ReactionInhibitionCreate, &
            ReactionInhibitionMonod, &
            ReactionInhibitionThreshold, &
            ReactionInhibitionSmoothstep, &
            ReactionInhibitionDestroy

contains

! ************************************************************************** !

function ReactionInhibitionCreate()
  !
  ! Allocate and initialize a inhibition object
  !
  ! Author: Glenn Hammond
  ! Date: 10/30/12, 11/21/23
  !

  implicit none

  type(inhibition_type), pointer :: ReactionInhibitionCreate

  type(inhibition_type), pointer :: inhibition

  allocate(inhibition)
  inhibition%id = 0
  inhibition%itype = 0
  inhibition%species_name = ''
  inhibition%inhibition_constant = UNINITIALIZED_DOUBLE
  inhibition%inhibition_constant2 = UNINITIALIZED_DOUBLE
  nullify(inhibition%next)

  ReactionInhibitionCreate => inhibition

end function ReactionInhibitionCreate

! ************************************************************************** !

subroutine ReactionInhibitionMonod(concentration,threshold_concentration, &
                                   inhibition_factor,derivative)
  !
  ! Calculates inhibition through the Monod term
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/23
  !
  implicit none

  PetscReal :: concentration
  PetscReal :: threshold_concentration
  PetscReal :: inhibition_factor
  PetscReal :: derivative

  PetscReal :: local_threshold_concentration
  PetscReal :: denominator

  local_threshold_concentration = dabs(threshold_concentration)
  denominator = local_threshold_concentration + concentration
  if (threshold_concentration < 0.d0) then ! inhibit above
    inhibition_factor = local_threshold_concentration / denominator
    derivative = -1.d0 * local_threshold_concentration / &
                         (denominator*denominator)
  else ! inhibit below
    inhibition_factor = concentration / denominator
    derivative = (denominator - concentration) / (denominator*denominator)
  endif

end subroutine ReactionInhibitionMonod

! ************************************************************************** !

subroutine ReactionInhibitionThreshold1(concentration, &
                                        threshold_concentration, &
                                        inhibition_factor,derivative)
  !
  ! Calculates threshold inhibition using the arc tangent function
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/23
  !
  implicit none

  PetscReal :: concentration
  PetscReal :: threshold_concentration
  PetscReal :: inhibition_factor
  PetscReal :: derivative

  PetscReal :: threshold_f

  threshold_f = 1.d5/dabs(threshold_concentration)
  call ReactionInhibitionThreshold2(concentration,threshold_concentration, &
                                    threshold_f,inhibition_factor,derivative)

end subroutine ReactionInhibitionThreshold1

 ! ************************************************************************** !

subroutine ReactionInhibitionThreshold2(concentration, &
                                        threshold_concentration, &
                                        threshold_constant, &
                                        inhibition_factor,derivative)
  !
  ! Calculates threshold inhibition using the arc tangent function
  !
  ! Author: Glenn Hammond
  ! Date: 05/17/23
  !

  implicit none

  PetscReal :: concentration
  PetscReal :: threshold_concentration
  PetscReal :: threshold_constant
  PetscReal :: inhibition_factor
  PetscReal :: derivative

  PetscReal :: tempreal

  tempreal = (concentration-dabs(threshold_concentration))*threshold_constant
  ! derivative of atan(X) = 1 / (1 + X^2) dX
  derivative = threshold_constant / (1.d0+tempreal*tempreal) / PI
  if (threshold_concentration < 0.d0) then ! INHIBIT_ABOVE_THRESHOLD
    inhibition_factor = 0.5d0 - atan(tempreal)/PI
    derivative = -1.d0 * derivative
  else ! INHIBIT_BELOW_THRESHOLD
    inhibition_factor = 0.5d0 + atan(tempreal)/PI
  endif

end subroutine ReactionInhibitionThreshold2

! ************************************************************************** !

subroutine ReactionInhibitionSmoothstep1(concentration, &
                                         threshold_concentration, &
                                         inhibition_factor,derivative)
  !
  ! Calculates threshold inhibition using sigmoid function
  !
  ! Author: Glenn Hammond
  ! Date: 11/27/23
  !
  implicit none

  PetscReal :: concentration
  PetscReal :: threshold_concentration
  PetscReal :: inhibition_factor
  PetscReal :: derivative

  PetscReal, parameter :: log10_interval = 3.d0

  call ReactionInhibitionSmoothstep2(concentration,threshold_concentration, &
                                     log10_interval,inhibition_factor, &
                                     derivative)

end subroutine ReactionInhibitionSmoothstep1

! ************************************************************************** !

subroutine ReactionInhibitionSmoothstep2(concentration, &
                                         threshold_concentration, &
                                         log10_interval,inhibition_factor, &
                                         derivative)
  !
  ! Calculates threshold inhibition using sigmoid function
  !
  ! Author: Glenn Hammond - based on Peishi Jiang implementation
  ! Date: 11/27/23
  !
  use Utility_module

  implicit none

  PetscReal :: concentration
  PetscReal :: threshold_concentration
  PetscReal :: log10_interval
  PetscReal :: inhibition_factor
  PetscReal :: derivative

  PetscReal :: log_inhibition
  PetscReal :: log_concentration
  PetscReal :: lower_bound

  log_inhibition = log10(dabs(threshold_concentration))
  log_concentration = log10(concentration)
  lower_bound = log_inhibition - 0.5d0 * log10_interval

  call Smoothstep(log_concentration,lower_bound, &
                  lower_bound+log10_interval,inhibition_factor, &
                  derivative)

  ! inhibition
  if (threshold_concentration < 0.d0) then ! INHIBIT_ABOVE_THRESHOLD
    inhibition_factor = 1.d0 - inhibition_factor
    derivative = -1.d0 * derivative / (concentration*LOG_TO_LN)
  else ! INHIBIT_BELOW_THRESHOLD
      derivative = derivative / (concentration*LOG_TO_LN)
  endif

end subroutine ReactionInhibitionSmoothstep2

! ************************************************************************** !

recursive subroutine ReactionInhibitionDestroy(inhibition)
  !
  ! Deallocates a inhibition object
  !
  ! Author: Glenn Hammond
  ! Date: 10/30/12
  !

  implicit none

  type(inhibition_type), pointer :: inhibition

  if (.not. associated(inhibition)) return

  call ReactionInhibitionDestroy(inhibition%next)

  deallocate(inhibition)
  nullify(inhibition)

end subroutine ReactionInhibitionDestroy

end module Reaction_Inhibition_Aux_module
