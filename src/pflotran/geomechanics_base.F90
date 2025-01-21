module Geomechanics_base_module

#include "petsc/finclude/petscvec.h"
  use Geomechanics_Regression_module
  use PMC_Geomechanics_class
  use Geomechanics_Realization_class

  implicit none

  private

  type, public :: geomechanics_base_type
    class(pmc_geomechanics_type), pointer :: process_model_coupler
    class(realization_geomech_type), pointer :: realization
    type(geomechanics_regression_type), pointer :: regression
  end type geomechanics_base_type

  public :: GeomechCreate, &
            GeomechDestroy

contains

! ************************************************************************** !

function GeomechCreate()

  ! Create a geomech object

  implicit none

  type (geomechanics_base_type),pointer :: GeomechCreate

  type (geomechanics_base_type),pointer :: geomech

  allocate(geomech)
  nullify(geomech%process_model_coupler)
  nullify(geomech%realization)
  nullify(geomech%regression)

  GeomechCreate => geomech

end function GeomechCreate

! ************************************************************************** !

subroutine GeomechDestroy(geomech)

  ! Destroys geomech object

  implicit none

  type (geomechanics_base_type),pointer :: geomech

  if (.not.associated(geomech)) return

  deallocate(geomech)
  nullify(geomech)

end subroutine GeomechDestroy

end module Geomechanics_base_module
