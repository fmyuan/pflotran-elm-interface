module Mapping_module

  use clm_pflotran_interface_type

  implicit none

#include "definitions.h"
#include "finclude/petsclog.h"

  private

  type, public :: inside_each_pflotran_cell

     PetscInt           :: num_clm_cells
     PetscInt,  pointer :: id_clm_cells(:)
     PetscReal, pointer :: perc_vol_overlap(:)
     PetscReal          :: total_vol_overlap
  end type inside_each_pflotran_cell

  type, public :: inside_each_clm_cell

     PetscInt           :: num_pflotran_cells
     PetscInt,  pointer :: id_pflotran_cells(:)
     PetscReal, pointer :: perc_vol_overlap(:)
     PetscReal          :: total_vol_overlap 
  end type inside_each_clm_cell

  type, public  :: mapping_type

     type(inside_each_pflotran_cell), pointer :: pf2clm(:)
     type(inside_each_clm_cell),      pointer :: clm2pf(:)

  end type mapping_type

  public :: MappingCreate

contains

  function MappingCreate()

    implicit none

    type(mapping_type), pointer :: MappingCreate

    type(mapping_type), pointer :: mapping

    allocate(mapping)
    nullify(mapping%pf2clm)
    nullify(mapping%clm2pf)

    MappingCreate => mapping

  end function MappingCreate

end module Mapping_module
