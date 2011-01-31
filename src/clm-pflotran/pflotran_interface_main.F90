program pflotran_interface_main

  use pflotran_model_module

  implicit none

#include "finclude/petscsysdef.h"

  PetscInt :: time

  type(pflotran_model_type),pointer :: pflotran_m

  allocate(pflotran_m)

  pflotran_m => pflotranModelCreate()

  call pflotranModelInitMapping2(pflotran_m)

  !call pflotranModelStepperRunInit(pflotran_m)


  !do time = 1,48

   !  call pflotranModelStepperRunTillPauseTime(pflotran_m,time * 1800.0d0)

  !enddo

  !call pflotranModelStepperRunTillPauseTime(pflotran_m,4.0d0)
  !call pflotranModelStepperRunTillPauseTime(pflotran_m,7.0d0)
  
  !call pflotranModelStepperRunFinalize(pflotran_m)

  call pflotranModelDestroy(pflotran_m)

end program pflotran_interface_main
