module Geomechanics_Attr_module

#include "petsc/finclude/petscvec.h"
  use Geomechanics_Regression_module
  use PMC_Geomechanics_class
  use Geomechanics_Realization_class

  implicit none

  private

  type, public :: geomechanics_attr_type
    class(pmc_geomechanics_type), pointer :: process_model_coupler
    class(realization_geomech_type), pointer :: realization
    type(geomechanics_regression_type), pointer :: regression
  end type geomechanics_attr_type

  public :: GeomechAttrCreate, &
            GeomechAttrDestroy

contains

! ************************************************************************** !

function GeomechAttrCreate()

  ! Create a geomech object

  implicit none

  type (geomechanics_attr_type),pointer :: GeomechAttrCreate

  type (geomechanics_attr_type),pointer :: geomech

  allocate(geomech)
  nullify(geomech%process_model_coupler)
  nullify(geomech%realization)
  nullify(geomech%regression)

  GeomechAttrCreate => geomech

end function GeomechAttrCreate

! ************************************************************************** !

subroutine GeomechAttrDestroy(geomech)

  ! Destroys geomech object

  implicit none

  type (geomechanics_attr_type),pointer :: geomech

  if (.not.associated(geomech)) return

  deallocate(geomech)
  nullify(geomech)

end subroutine GeomechAttrDestroy

! ************************************************************************** !

end module Geomechanics_Attr_module
