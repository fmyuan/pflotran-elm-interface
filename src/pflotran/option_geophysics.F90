module Option_Geophysics_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: geophysics_option_type

     PetscInt :: num_electrodes

  end type geophysics_option_type

  public :: OptionGeophysicsCreate, &
            OptionGeophysicsInitAll, &
            OptionGeophysicsInitRealization, &
            OptionGeophysicsDestroy

contains

! ************************************************************************** !

function OptionGeophysicsCreate()
  !
  ! Allocates and initializes a new Option object
  !
  ! Author: Glenn Hammond
  ! Date: 02/12/21
  !

  implicit none

  type(geophysics_option_type), pointer :: OptionGeophysicsCreate

  type(geophysics_option_type), pointer :: option

  allocate(option)

  ! DO NOT initialize members of the option type here.  One must decide
  ! whether the member needs initialization once for all stochastic
  ! simulations or initialization for every realization (e.g. within multiple
  ! stochastic simulations).  This is done in OptionInitAll() and
  ! OptionInitRealization()
  call OptionGeophysicsInitAll(option)
  OptionGeophysicsCreate => option

end function OptionGeophysicsCreate

! ************************************************************************** !

subroutine OptionGeophysicsInitAll(option)
  !
  ! Initializes all option variables
  !
  ! Author: Glenn Hammond
  ! Date: 02/12/21
  !

  implicit none

  type(geophysics_option_type) :: option

  ! These variables should only be initialized once at the beginning of a
  ! PFLOTRAN run (regardless of whether stochastic)

  call OptionGeophysicsInitRealization(option)

end subroutine OptionGeophysicsInitAll

! ************************************************************************** !

subroutine OptionGeophysicsInitRealization(option)
  !
  ! Initializes option variables specific to a single
  ! realization
  !
  ! Author: Glenn Hammond
  ! Date: 02/12/21
  !

  implicit none

  type(geophysics_option_type) :: option

  ! These variables should be initialized once at the beginning of every
  ! PFLOTRAN realization or simulation of a single realization

  option%num_electrodes = UNINITIALIZED_INTEGER

end subroutine OptionGeophysicsInitRealization

! ************************************************************************** !

subroutine OptionGeophysicsDestroy(option)
  !
  ! Deallocates an option
  !
  ! Author: Glenn Hammond
  ! Date: 02/12/21
  !

  implicit none

  type(geophysics_option_type), pointer :: option

  if (.not.associated(option)) return

  deallocate(option)
  nullify(option)

end subroutine OptionGeophysicsDestroy

end module Option_Geophysics_module
