
module pflotran_clm_interface_type

  implicit none

#include "definitions.h"

  private

  type, public :: pflotran_clm_type

     PetscReal, pointer :: zwt(:)
     PetscReal, pointer :: alpha(:)
     PetscReal, pointer :: lambda(:)
     PetscReal, pointer :: qsrc_flx(:)
     PetscReal, pointer :: sat_new(:)

  end type pflotran_clm_type



  type(pflotran_clm_type) , public, target , save :: pf_clm_data

end module pflotran_clm_interface_type
