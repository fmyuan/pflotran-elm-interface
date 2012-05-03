program pflotran_interface_main

  use pflotran_model_module
  !use clm_pflotran_interface_data
  use Mapping_module
  use Input_module
  use HDF5_aux_module
  use hdf5
  use Grid_module
  use Field_module
  use Richards_Aux_module
  use Richards_module
  use Discretization_module
  
  implicit none

#include "finclude/petscsysdef.h"

  type(pflotran_model_type),pointer  :: pflotran_m
  !type(clm_pflotran_interface_data_type),pointer    :: clm_pf_idata

  
  PetscInt :: time

  PetscInt, pointer                  :: clm_cell_ids(:)
  PetscInt                           :: clm_npts, ii,fileid,num_u_a,jj
  PetscInt                           :: npts
  PetscScalar, pointer :: hksat_x_loc(:) ! hydraulic conductivity in x-dir at saturation (mm H2O /s)

  ! To read HDF5 soil properties  
  character(len=MAXSTRINGLENGTH)     :: filename
  character(len=MAXSTRINGLENGTH)     :: group_name
  character(len=MAXSTRINGLENGTH)     :: dataset_name
  character(len=MAXWORDLENGTH)       :: card

  PetscInt :: PRINT_RANK    
  PRINT_RANK = 0

  ! A
  allocate(pflotran_m)
  !allocate(clm_pf_idata)
  
  ! Create the model and CLM data
  pflotran_m => pflotranModelCreate()

  call pflotranModelStepperRunInit(pflotran_m)
  do time = 1,10
     call pflotranModelStepperRunTillPauseTime(pflotran_m,time * 3600.0d0)
  enddo
  call pflotranModelStepperRunFinalize(pflotran_m)

  call pflotranModelDestroy(pflotran_m)


end program pflotran_interface_main
