program pflotran_interface_main

  use pflotran_model_module
  use clm_pflotran_interface_data
  implicit none

#include "finclude/petscsysdef.h"

  PetscInt :: time

  type(pflotran_model_type),pointer :: pflotran_m
  type(clm_pflotran_data),pointer :: clmpf_data
  PetscInt, pointer                  :: clm_cell_ids(:)
  PetscInt                           :: clm_npts, ii

  allocate(pflotran_m)
  allocate(clmpf_data)

  pflotran_m => pflotranModelCreate()
  clmpf_data => clm_pf_data_create()

  !clm_npts = 1
  !allocate(clm_cell_ids(clm_npts))
  !clm_cell_ids(1) = pflotran_m%option%myrank

  clm_npts = 40/pflotran_m%option%mycommsize
  !clm_npts = 2
  allocate(clm_cell_ids(clm_npts))
  do ii = 1,clm_npts
    clm_cell_ids(ii) = clm_npts*pflotran_m%option%myrank + ii -1
  	if(pflotran_m%option%myrank.eq.0) write(*,*),clm_cell_ids(ii)
  enddo

  call pflotranModelInitMapping2(pflotran_m,clmpf_data,clm_cell_ids,clm_npts)

  !call pflotranModelStepperRunInit(pflotran_m)


  !do time = 1,48

   !  call pflotranModelStepperRunTillPauseTime(pflotran_m,time * 1800.0d0)

  !enddo

  !call pflotranModelStepperRunTillPauseTime(pflotran_m,4.0d0)
  !call pflotranModelStepperRunTillPauseTime(pflotran_m,7.0d0)
  
  !call pflotranModelStepperRunFinalize(pflotran_m)

  call pflotranModelDestroy(pflotran_m)

end program pflotran_interface_main
