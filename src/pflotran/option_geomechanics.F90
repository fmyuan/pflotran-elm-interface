module Option_Geomechanics_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: geomechanics_option_type

    PetscBool :: geomech_initial ! initial_solve
    PetscReal :: geomech_time
    PetscInt :: subsurf_coupling
    PetscReal :: gravity(3)
    PetscInt :: split_scheme

  end type geomechanics_option_type

  public :: OptionGeomechanicsCreate, &
            OptionGeomechanicsInitAll, &
            OptionGeomechanicsInitRealization, &
            OptionGeomechanicsDestroy

contains

! ************************************************************************** !

function OptionGeomechanicsCreate()
  !
  ! Allocates and initializes a new Option object
  !
  ! Author: Jumanah Al Kubaisy
  ! Date: 2/3/25
  !

  implicit none

  type(geomechanics_option_type), pointer :: OptionGeomechanicsCreate

  type(geomechanics_option_type), pointer :: option

  allocate(option)

  ! DO NOT initialize members of the option type here.  One must decide
  ! whether the member needs initialization once for all stochastic
  ! simulations or initialization for every realization (e.g. within multiple
  ! stochastic simulations).  This is done in OptionInitAll() and
  ! OptionInitRealization()
  call OptionGeomechanicsInitAll(option)
  OptionGeomechanicsCreate => option

end function OptionGeomechanicsCreate

! ************************************************************************** !

subroutine OptionGeomechanicsInitAll(option)
  !
  ! Initializes all option variables
  !
  ! Author: Jumanah Al Kubaisy
  ! Date: 2/3/25
  !

  implicit none

  type(geomechanics_option_type) :: option

  ! These variables should only be initialized once at the beginning of a
  ! PFLOTRAN run (regardless of whether stochastic)

  call OptionGeomechanicsInitRealization(option)

end subroutine OptionGeomechanicsInitAll

! ************************************************************************** !

subroutine OptionGeomechanicsInitRealization(option)
  !
  ! Initializes option variables specific to a single
  ! realization
  !
  ! Author: Jumanah Al Kubaisy
  ! Date: 2/3/25
  !

  implicit none

  type(geomechanics_option_type) :: option

  ! These variables should be initialized once at the beginning of every
  ! PFLOTRAN realization or simulation of a single realization

  option%geomech_initial = PETSC_FALSE
  option%geomech_time = 0.d0
  option%subsurf_coupling = 0
  option%split_scheme = 0
  option%gravity(:) = 0.d0
  option%gravity(3) = -1.d0*EARTH_GRAVITY    ! m/s^2

end subroutine OptionGeomechanicsInitRealization

! ************************************************************************** !

subroutine OptionGeomechanicsDestroy(option)
  !
  ! Deallocates an option
  !
  ! Author: Jumanah Al Kubaisy
  ! Date: 2/3/25
  !

  implicit none

  type(geomechanics_option_type), pointer :: option

  if (.not.associated(option)) return

  deallocate(option)
  nullify(option)

end subroutine OptionGeomechanicsDestroy

end module Option_Geomechanics_module
