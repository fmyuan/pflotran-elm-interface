module Upwind_Direction_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  implicit none

  private

  PetscBool, public :: fix_upwind_direction = PETSC_TRUE
  PetscBool, public :: update_upwind_direction = PETSC_FALSE
  PetscBool, public :: count_upwind_direction_flip = PETSC_FALSE
  PetscInt, public :: upwind_dir_update_freq = 9999

  ! variables that track the number of times the upwind direction changes
  ! during the residual and Jacobian calculations.
  PetscInt, public :: liq_upwind_flip_count_by_res
  PetscInt, public :: gas_upwind_flip_count_by_res
  PetscInt, public :: liq_bc_upwind_flip_count_by_res
  PetscInt, public :: gas_bc_upwind_flip_count_by_res
  PetscInt, public :: liq_upwind_flip_count_by_jac
  PetscInt, public :: gas_upwind_flip_count_by_jac
  PetscInt, public :: liq_bc_upwind_flip_count_by_jac
  PetscInt, public :: gas_bc_upwind_flip_count_by_jac

  public :: UpwindDirectionInit, &
            UpwindDirection, &
            UpwindDirectionPrintStats

contains

! ************************************************************************** !

subroutine UpwindDirectionInit()
  !
  ! Allocates and initializes variables for upwind direction
  !
  ! Author: Glenn Hammond
  ! Date: 10/19/18
  !
  implicit none

  update_upwind_direction = PETSC_FALSE

  ! these counters are used for performance analysis only.  they will not
  ! affect simulation results.
  liq_upwind_flip_count_by_res = 0
  gas_upwind_flip_count_by_res = 0
  liq_bc_upwind_flip_count_by_res = 0
  gas_bc_upwind_flip_count_by_res = 0
  liq_upwind_flip_count_by_jac = 0
  gas_upwind_flip_count_by_jac = 0
  liq_bc_upwind_flip_count_by_jac = 0
  gas_bc_upwind_flip_count_by_jac = 0

end subroutine UpwindDirectionInit

! ************************************************************************** !

function UpwindDirection(upwind_direction,delta_pressure, &
                         derivative_call, &
                         count_upwind_direction_flip_, &
                         upwind_flip_count_by_res, &
                         upwind_flip_count_by_jac)
  !
  ! Calculates the upwind direction
  !
  ! Author: Glenn Hammond
  ! Date: 10/19/18
  !
  implicit none

  PetscInt :: upwind_direction
  PetscReal :: delta_pressure
  PetscBool :: upwind
  PetscBool :: derivative_call
  PetscBool :: count_upwind_direction_flip_
  PetscInt :: upwind_flip_count_by_res
  PetscInt :: upwind_flip_count_by_jac

  PetscBool :: UpwindDirection

  PetscInt :: prev_upwind_direction
  PetscInt :: new_upwind_direction

  if (fix_upwind_direction) then
    if (update_upwind_direction .or. count_upwind_direction_flip_) then
      prev_upwind_direction = upwind_direction
      if (delta_pressure >= 0.d0) then
        ! positive means upstream
        new_upwind_direction = iabs(prev_upwind_direction)
      else
        ! negative means downstream
        new_upwind_direction = -iabs(prev_upwind_direction)
      endif
      if (count_upwind_direction_flip_) then
        if (new_upwind_direction /= prev_upwind_direction) then
          if (derivative_call) then
            upwind_flip_count_by_jac = upwind_flip_count_by_jac + 1
          else
            upwind_flip_count_by_res = upwind_flip_count_by_res + 1
          endif
        endif
      endif
      if (update_upwind_direction) then
        upwind_direction = new_upwind_direction
      endif
    endif
    UpwindDirection = (upwind_direction > 0)
  else
    UpwindDirection = (delta_pressure >= 0.d0)
  endif

end function UpwindDirection

! ************************************************************************** !
subroutine UpwindDirectionPrintStats(option)
  !
  ! Prints statistic recorded for upwind direction
  !
  ! Author: Glenn Hammond
  ! Date: 10/22/18
  !
  use Option_module

  implicit none

  type(option_type) :: option

  if (fix_upwind_direction .and. &
      count_upwind_direction_flip .and. &
      OptionPrintToScreen(option)) then
    print *, 'Res upwind dir. flip (liq): ', liq_upwind_flip_count_by_res
    print *, 'Res upwind dir. flip (gas): ', gas_upwind_flip_count_by_res
    print *, 'Res upwind dir. flip (liq bc): ', liq_bc_upwind_flip_count_by_res
    print *, 'Res upwind dir. flip (gas bc): ', gas_bc_upwind_flip_count_by_res
    print *, 'Jac upwind dir. flip (liq): ', liq_upwind_flip_count_by_jac
    print *, 'Jac upwind dir. flip (gas): ', gas_upwind_flip_count_by_jac
    print *, 'Jac upwind dir. flip (liq bc): ', liq_bc_upwind_flip_count_by_jac
    print *, 'Jac upwind dir. flip (gas bc): ', gas_bc_upwind_flip_count_by_jac
  endif

end subroutine UpwindDirectionPrintStats

! ************************************************************************** !

end module Upwind_Direction_module
