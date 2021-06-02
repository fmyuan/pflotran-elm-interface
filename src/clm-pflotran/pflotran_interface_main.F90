program pflotran_interface_main

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"

  use petscsys
  use petscvec
  use pflotran_model_module         , only : pflotran_model_type, pflotranModelCreate, &
       pflotranModelInitMapping, pflotranModelStepperRunInit, &
       pflotranModelStepperRunTillPauseTime, pflotranModelDestroy
  use clm_pflotran_interface_data
  use Mapping_module
  use Input_Aux_module
  use Option_module
  use Input_Aux_module
  use String_module
  
  use Simulation_Base_class         , only : simulation_base_type
  use Simulation_Subsurface_class   , only : simulation_subsurface_type, &
                                             SimSubsurfCast
  use Realization_Base_class        , only : realization_base_type
  use Timestepper_Base_class        , only : TS_STOP_END_SIMULATION

  use PFLOTRAN_Constants_module

  implicit none


  type(pflotran_model_type)       , pointer :: pflotran_m
  class(realization_base_type)    , pointer :: realization
  class(simulation_subsurface_type), pointer :: simulation
  
  PetscErrorCode                            :: ierr
  PetscInt                                  :: time

  PetscInt                        , pointer :: clm_cell_ids(:), clm_surf_cell_ids(:)
  PetscInt                                  :: clm_npts, clm_surf_npts, ii
  PetscInt                                  :: ntimes

  ! To read HDF5 soil properties  
  character(len=MAXSTRINGLENGTH)            :: filename
  character(len=MAXSTRINGLENGTH)            :: string
  PetscBool                                 :: pflotranin_option_found
  PetscBool                                 :: input_prefix_option_found
  character(len=MAXSTRINGLENGTH)  , pointer :: strings(:)

  PetscInt                                  :: PRINT_RANK    
  PRINT_RANK = 0

  call MPI_Init(ierr)
 
  ! Create the model
  pflotran_m => pflotranModelCreate(MPI_COMM_WORLD, filename)
  simulation => SimSubsurfCast(pflotran_m%simulation)

  select type (simulation)
    class is (simulation_subsurface_type)
       realization => simulation%realization
    class default
       nullify(realization)
       pflotran_m%option%io_buffer = "ERROR: pflotran model only works on combinations of subsurface and surface simulations."
       call PrintErrMsg(pflotran_m%option)
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
      call PrintErrMsg(pflotran_m%option)
    endif
  endif

  call CLMPFLOTRANIDataInit()

  clm_pf_idata%nlclm_sub = clm_npts
  clm_pf_idata%ngclm_sub = clm_npts
  clm_pf_idata%nlpf_sub  = realization%patch%grid%nlmax
  clm_pf_idata%ngpf_sub  = realization%patch%grid%ngmax

  ! Allocate memory for CLM-PFLOTRAN data transfer
  call CLMPFLOTRANIDataCreateVec(MPI_COMM_WORLD)

  ! Set mapping between CLM and PFLOTRAN
  call pflotranModelInitMapping(pflotran_m, clm_cell_ids, clm_npts, 1)
  !call pflotranModelInitMapping(pflotran_m, clm_cell_ids, clm_npts, CLM_SUB_TO_PF_EXTENDED_SUB)
  !call pflotranModelInitMapping(pflotran_m, clm_cell_ids, clm_npts, PF_SUB_TO_CLM_SUB)
  !call pflotranModelInitMapping(pflotran_m, clm_surf_cell_ids, clm_surf_npts, CLM_SRF_TO_PF_2DSUB)
#ifdef SURFACE_FLOW
  clm_pf_idata%nlclm_srf = clm_surf_npts
  clm_pf_idata%ngclm_srf = clm_surf_npts
  if(pflotran_m%option%nsurfflowdof>0) then
    clm_pf_idata%nlpf_srf  = surf_realization%patch%grid%nlmax
    clm_pf_idata%ngpf_srf  = surf_realization%patch%grid%ngmax
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

  ! flag ensures shutdown due to successful run.
  simulation%stop_flag = TS_STOP_END_SIMULATION
  
  ! Finalize PFLOTRAN Stepper
  !call pflotranModelStepperRunFinalize(pflotran_m)

  call CLMPFLOTRANIDataDestroy()

  call pflotranModelDestroy(pflotran_m)

  call MPI_Finalize(ierr)

end program pflotran_interface_main
