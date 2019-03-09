module NW_Transport_Aux_module

  ! this module cannot depend on any other modules besides Option_module
  ! and Matrix_Block_Aux_module
  use Matrix_Block_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"

  PetscReal, public :: nwt_itol_scaled_res = UNINITIALIZED_DOUBLE
  PetscReal, public :: nwt_itol_rel_update = UNINITIALIZED_DOUBLE
  PetscReal, public :: nwt_min_saturation = 0.d0

end module NW_Transport_Aux_module