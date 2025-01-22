module Geomechanics_linkage_module

#include "petsc/finclude/petscvec.h"
  use Geomechanics_Regression_module
  use PMC_Geomechanics_class
  use Geomechanics_Realization_class

  implicit none

  private

  type, public :: geomechanics_linkage_type
    class(pmc_geomechanics_type), pointer :: process_model_coupler
    class(realization_geomech_type), pointer :: realization
    type(geomechanics_regression_type), pointer :: regression
  end type geomechanics_linkage_type

  public :: GeomechLinkageCreate, &
            GeomechLinkageDestroy

contains

! ************************************************************************** !

function GeomechLinkageCreate()

  ! Create a geomech object

  implicit none

  type (geomechanics_linkage_type),pointer :: GeomechLinkageCreate

  type (geomechanics_linkage_type),pointer :: geomech

  allocate(geomech)
  nullify(geomech%process_model_coupler)
  nullify(geomech%realization)
  nullify(geomech%regression)

  GeomechLinkageCreate => geomech

end function GeomechLinkageCreate

! ************************************************************************** !

subroutine GeomechlinkageDestroy(geomech)

  ! Destroys geomech object

  implicit none

  type (geomechanics_linkage_type),pointer :: geomech

  if (.not.associated(geomech)) return

  deallocate(geomech)
  nullify(geomech)

end subroutine GeomechLinkageDestroy

! ************************************************************************** !

end module Geomechanics_linkage_module
