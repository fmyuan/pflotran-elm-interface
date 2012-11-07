program pflotran_interface_main

  use pflotran_model_module
  use clm_pflotran_interface_data
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
  pflotran_m => pflotranModelCreate(1)

#if 1
  if (pflotran_m%option%mycommsize == 1) then
    clm_npts = 30*30*10
    allocate (clm_cell_ids(clm_npts))
    do ii = 1,clm_npts
      clm_cell_ids(ii) = ii-1
      !clm_cell_ids(ii) = ii
    enddo
  else
    clm_npts = 30*30*10/2
    allocate (clm_cell_ids(clm_npts))
    do ii = 1,clm_npts
      !clm_cell_ids(ii) = ii-1 + 30*30*10/2*pflotran_m%option%myrank
      clm_cell_ids(ii) = ii + 30*30*10/2*pflotran_m%option%myrank
    enddo
  endif
#endif

#if 0
  if (pflotran_m%option%mycommsize == 1) then
    clm_npts = 8
    allocate(clm_cell_ids(clm_npts))
    do ii = 1,clm_npts
      clm_cell_ids(ii) = ii-1
    enddo
  else
    if(pflotran_m%option%myrank == 0) then
      clm_npts = 5
      allocate(clm_cell_ids(clm_npts))
      clm_cell_ids(1) = 4
      clm_cell_ids(2) = 5
      clm_cell_ids(3) = 6
      clm_cell_ids(4) = 7
      clm_cell_ids(5) = 8
    else
      clm_npts = 3
      allocate(clm_cell_ids(clm_npts))
      clm_cell_ids(1) = 3
      clm_cell_ids(2) = 2
      clm_cell_ids(3) = 1
    endif
  endif
#endif
  
    clm_pf_idata%nlclm = clm_npts
    clm_pf_idata%ngclm = clm_npts
    clm_pf_idata%nlpf  = pflotran_m%realization%patch%grid%nlmax
    clm_pf_idata%ngpf  = pflotran_m%realization%patch%grid%ngmax
  call clm_pflotran_interface_data_allocate_memory(1)
    call pflotranModelInitMapping3(pflotran_m, clm_cell_ids,clm_npts, CLM2PF_FLUX_MAP_ID)
    call pflotranModelInitMapping3(pflotran_m, clm_cell_ids,clm_npts, CLM2PF_SOIL_MAP_ID)
    call pflotranModelInitMapping3(pflotran_m, clm_cell_ids,clm_npts, PF2CLM_FLUX_MAP_ID)
    call pflotranModelSetSoilProp3(pflotran_m)
  !call pflotranModelInitMapping3(pflotran_m, clm_cell_ids, clm_npts,CLM2PF_FLUX_MAP_ID)


  !write(*,*), 'Done pflotranModelCreate()'

  call pflotranModelStepperRunInit(pflotran_m)
  do time = 1,1
     call pflotranModelStepperRunTillPauseTime(pflotran_m,time * 3600.0d0)
  enddo
  call pflotranModelStepperRunFinalize(pflotran_m)

  call pflotranModelDestroy(pflotran_m)


end program pflotran_interface_main
