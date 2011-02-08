program pflotran_interface_main

  use pflotran_model_module
  use clm_pflotran_interface_data
  implicit none

#include "finclude/petscsysdef.h"

  PetscInt :: time

  type(pflotran_model_type),pointer :: pflotran_m
  type(clm_pflotran_data),pointer :: clmpf_data

  allocate(pflotran_m)
  allocate(clmpf_data)

  pflotran_m => pflotranModelCreate()
  clmpf_data => clm_pf_data_create()

  call pflotranModelInitMapping2(pflotran_m, clmpf_data)

  !call pflotranModelStepperRunInit(pflotran_m)


  !do time = 1,48

   !  call pflotranModelStepperRunTillPauseTime(pflotran_m,time * 1800.0d0)

  !enddo

  !call pflotranModelStepperRunTillPauseTime(pflotran_m,4.0d0)
  !call pflotranModelStepperRunTillPauseTime(pflotran_m,7.0d0)
  
  !call pflotranModelStepperRunFinalize(pflotran_m)

  call pflotranModelDestroy(pflotran_m)

end program pflotran_interface_main
