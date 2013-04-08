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
  use Option_module
  
  implicit none

#include "finclude/petscsysdef.h"

  type(pflotran_model_type),pointer  :: pflotran_m

  
  PetscErrorCode :: ierr
  PetscInt :: time

  PetscInt, pointer                  :: clm_cell_ids(:)
  PetscInt                           :: clm_npts, ii,fileid,num_u_a,jj
  PetscInt                           :: npts
  PetscInt           :: ntimes

  ! To read HDF5 soil properties  
  character(len=MAXSTRINGLENGTH)     :: filename
  character(len=MAXSTRINGLENGTH)     :: group_name
  character(len=MAXSTRINGLENGTH)     :: dataset_name
  character(len=MAXWORDLENGTH)       :: card

  PetscInt :: PRINT_RANK    
  PRINT_RANK = 0

  call MPI_Init(ierr)

  ! A
  allocate(pflotran_m)
  
  ! Create the model
  pflotran_m => pflotranModelCreate(MPI_COMM_WORLD)

  ! Set up CLM cell ids
  if (pflotran_m%option%mycommsize == 1) then
    clm_npts = 200*10
    allocate (clm_cell_ids(clm_npts))
    do ii = 1,clm_npts
      clm_cell_ids(ii) = ii-1
    enddo
  else
    if (pflotran_m%option%mycommsize == 2) then
      clm_npts = 200*10/2
      allocate (clm_cell_ids(clm_npts))
      do ii = 1,clm_npts
        clm_cell_ids(ii) = ii-1 + 200*10/2*pflotran_m%option%myrank
      enddo
    else
      pflotran_m%option%io_buffer = 'The example can only run with max 2 procs.'
      call printErrMsg(pflotran_m%option)
    endif
  endif

  clm_pf_idata%nlclm = clm_npts
  clm_pf_idata%ngclm = clm_npts
  clm_pf_idata%nlpf  = pflotran_m%realization%patch%grid%nlmax
  clm_pf_idata%ngpf  = pflotran_m%realization%patch%grid%ngmax

  ! Allocate memory for CLM-PFLOTRAN data transfer
  call CLMPFLOTRANIdataCreate(MPI_COMM_WORLD)
  
  ! Set mapping between CLM and PFLOTRAN
  call pflotranModelInitMapping(pflotran_m, clm_cell_ids,clm_npts, CLM2PF_FLUX_MAP_ID)
  call pflotranModelInitMapping(pflotran_m, clm_cell_ids,clm_npts, CLM2PF_SOIL_MAP_ID)
  call pflotranModelInitMapping(pflotran_m, clm_cell_ids,clm_npts, PF2CLM_FLUX_MAP_ID)
! call pflotranModelSetSoilProp(pflotran_m)

  ! Initialize PFLOTRAN Stepper
  call pflotranModelStepperRunInit(pflotran_m)
  
  ! Run PFLOTRAN 'ntimes'. For each time run PFLOTRAN for 3600s and PFLOTRAN
  ! can take multiple smaller steps to reach the 3600s interval.
  ntimes = 10
  do time = 1,ntimes
     
     ! When coupled with CLM:
     ! GetSourceSinkFromCLM()
  
     ! Run PFLOTRAN
     call pflotranModelStepperRunTillPauseTime(pflotran_m,time * 3600.0d0)
     
     ! When coupled with CLM
     ! PassSaturationValuesToCLM()
     
  enddo
  
  ! Finalize PFLOTRAN Stepper
  call pflotranModelStepperRunFinalize(pflotran_m)

  call pflotranModelDestroy(pflotran_m)

  call MPI_Finalize(ierr)

end program pflotran_interface_main
