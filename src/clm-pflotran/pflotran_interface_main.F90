program pflotran_interface_main

  use pflotran_model_module, only : pflotran_model_type, pflotranModelCreate, &
       pflotranModelInitMapping, pflotranModelStepperRunInit, &
       pflotranModelStepperRunTillPauseTime, pflotranModelDestroy, &
       CLM_SRF_TO_PF_SRF, PF_SRF_TO_CLM_SRF
  use clm_pflotran_interface_data
  use Mapping_module
  use Input_Aux_module
  use Option_module
  
  use Simulation_Base_class, only : simulation_base_type
  use Subsurface_Simulation_class, only : subsurface_simulation_type
  use Surface_Simulation_class, only : surface_simulation_type
  use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
  use Realization_Base_class, only : realization_base_type
  use Surface_Realization_class, only : surface_realization_type

  use PFLOTRAN_Constants_module

  implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"

  type(pflotran_model_type), pointer  :: pflotran_m
  class(realization_base_type), pointer :: realization
  class(surface_realization_type), pointer :: surf_realization

  
  PetscErrorCode :: ierr
  PetscInt :: time

  PetscInt, pointer                  :: clm_cell_ids(:), clm_surf_cell_ids(:)
  PetscInt                           :: clm_npts, clm_surf_npts, ii, fileid, num_u_a, jj
  PetscInt                           :: npts
  PetscInt           :: ntimes
  Vec :: myvec

  ! To read HDF5 soil properties  
  character(len=MAXSTRINGLENGTH)     :: filename
  character(len=MAXSTRINGLENGTH)     :: group_name
  character(len=MAXSTRINGLENGTH)     :: dataset_name
  character(len=MAXWORDLENGTH)       :: card

  PetscInt :: PRINT_RANK    
  PRINT_RANK = 0

  call MPI_Init(ierr)

  ! Create the model
  filename = 'pflotran'
  pflotran_m => pflotranModelCreate(MPI_COMM_WORLD, filename)

  select type (simulation => pflotran_m%simulation)
    class is (subsurface_simulation_type)
       realization => simulation%realization
       nullify(surf_realization)
    class is (surfsubsurface_simulation_type)
       realization => simulation%realization
       surf_realization => simulation%surf_realization
    class is (surface_simulation_type)
       nullify(realization)
       surf_realization => simulation%surf_realization
    class default
       nullify(realization)
       nullify(surf_realization)
       pflotran_m%option%io_buffer = "ERROR: pflotran model only works on combinations of subsurface and surface simulations."
       call printErrMsg(pflotran_m%option)
   end select

  ! Set up CLM cell ids
  if (pflotran_m%option%mycommsize == 1) then
    clm_npts = 5000*10
    clm_surf_npts = 5000
    allocate (clm_cell_ids(clm_npts))
    allocate (clm_surf_cell_ids(clm_surf_npts))
    do ii = 1, clm_npts
      clm_cell_ids(ii) = ii-1
    enddo
    do ii = 1, clm_surf_npts
      clm_surf_cell_ids(ii) = (ii-1)*10
    enddo
  else
    if (pflotran_m%option%mycommsize == 2) then
      clm_surf_npts = 5000/2
      clm_npts       = clm_surf_npts*10
      allocate (clm_cell_ids(clm_npts))
      allocate (clm_surf_cell_ids(clm_surf_npts))
      do ii = 1, clm_npts
        clm_cell_ids(ii) = ii-1 + clm_npts*pflotran_m%option%myrank
      enddo
      do ii = 1, clm_surf_npts
        clm_surf_cell_ids(ii) = (ii-1)*10 + clm_npts*pflotran_m%option%myrank
      enddo
    else
      pflotran_m%option%io_buffer = 'The example can only run with max 2 procs.'
      call printErrMsg(pflotran_m%option)
    endif
  endif

  call CLMPFLOTRANIDataInit()

  clm_pf_idata%nlclm_3d = clm_npts
  clm_pf_idata%ngclm_3d = clm_npts
  clm_pf_idata%nlpf_3d  = realization%patch%grid%nlmax
  clm_pf_idata%ngpf_3d  = realization%patch%grid%ngmax

  ! Allocate memory for CLM-PFLOTRAN data transfer
  call CLMPFLOTRANIDataCreateVec(MPI_COMM_WORLD)

  ! Set mapping between CLM and PFLOTRAN
  call pflotranModelInitMapping(pflotran_m, clm_cell_ids, clm_npts, 1)
  !call pflotranModelInitMapping(pflotran_m, clm_cell_ids, clm_npts, CLM_SUB_TO_PF_EXTENDED_SUB)
  !call pflotranModelInitMapping(pflotran_m, clm_cell_ids, clm_npts, PF_SUB_TO_CLM_SUB)
  !call pflotranModelInitMapping(pflotran_m, clm_surf_cell_ids, clm_surf_npts, CLM_SRF_TO_PF_2DSUB)
#ifdef SURFACE_FLOW
  clm_pf_idata%nlclm_2d = clm_surf_npts
  clm_pf_idata%ngclm_2d = clm_surf_npts
  if(pflotran_m%option%nsurfflowdof>0) then
    clm_pf_idata%nlpf_2d  = surf_realization%patch%grid%nlmax
    clm_pf_idata%ngpf_2d  = surf_realization%patch%grid%ngmax
    call pflotranModelInitMapping(pflotran_m, clm_surf_cell_ids, clm_surf_npts, PF_SRF_TO_CLM_SRF)
    call pflotranModelInitMapping(pflotran_m, clm_surf_cell_ids, clm_surf_npts, CLM_SRF_TO_PF_SRF)
  endif
#endif

! call pflotranModelSetSoilProp(pflotran_m)

  ! Initialize PFLOTRAN Stepper
  call pflotranModelStepperRunInit(pflotran_m)
  
  ! Run PFLOTRAN 'ntimes'. For each time run PFLOTRAN for 3600s and PFLOTRAN
  ! can take multiple smaller steps to reach the 3600s interval.
  ntimes = 10
  do time = 1, ntimes
     
     ! When coupled with CLM:
     ! GetSourceSinkFromCLM()
  
     ! Run PFLOTRAN
     call pflotranModelStepperRunTillPauseTime(pflotran_m, time * 3600.0d0)
     
     ! When coupled with CLM
     ! PassSaturationValuesToCLM()
     
  enddo
  
  ! Finalize PFLOTRAN Stepper
  !call pflotranModelStepperRunFinalize(pflotran_m)

  call CLMPFLOTRANIDataDestroy()

  call pflotranModelDestroy(pflotran_m)

  call MPI_Finalize(ierr)

end program pflotran_interface_main
