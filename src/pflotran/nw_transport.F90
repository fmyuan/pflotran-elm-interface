module NW_Transport_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Global_Aux_module
  use Material_Aux_class
  use PM_NWT_class
  use NW_Transport_Aux_module
  
  use PFLOTRAN_Constants_module
  
  implicit none
  
  private
  
  PetscReal, parameter :: perturbation_tolerance = 1.d-5
  
  public :: NWTTimeCut, &
            NWTSetup
            
contains

! ************************************************************************** !

subroutine NWTTimeCut(realization)
  ! 
  ! Resets arrays for a time step cut.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/12/2019
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Global_module
  !use Secondary_Continuum_module, only : SecondaryRTTimeCut
 
  implicit none
  
  type(realization_subsurface_type) :: realization
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  
  PetscErrorCode :: ierr

  field => realization%field
  option => realization%option
 
  ! copy previous solution back to current solution
  call VecCopy(field%tran_yy,field%tran_xx,ierr);CHKERRQ(ierr)
  
  ! set densities and saturations to t+dt
  if (realization%option%nflowdof > 0) then
    call GlobalWeightAuxVars(realization, &
                             realization%option%transport%tran_weight_t1)
  endif

  !if (option%use_mc) then
  !  call SecondaryRTTimeCut(realization)
  !endif
 
end subroutine NWTTimeCut

! ************************************************************************** !

subroutine NWTSetup(realization)
  ! 
  ! Sets up the nuclear waste transport realization.
  ! Author: Jenn Frederick
  ! Date: 03/12/2019
  ! 
  
  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Grid_module
  use Material_module
  use Material_Aux_class
  use Coupler_module
  use Condition_module
  use Connection_module
  use Fluid_module
  use Transport_Constraint_module
  use Output_Aux_module
  
  implicit none

  type(realization_subsurface_type) :: realization
  
  type(patch_type), pointer :: patch
  type(nw_trans_realization_type), pointer :: nw_trans
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(material_property_type), pointer :: cur_material_property
  type(fluid_property_type), pointer :: cur_fluid_property
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(output_variable_list_type), pointer :: list
  
  PetscInt :: ghosted_id, local_id
  PetscInt :: iconn, sum_connection 
  PetscInt :: iphase
  PetscInt :: flag(3)
  PetscInt :: ncomp, nphase
  
  patch => realization%patch
  grid => patch%grid
  nw_trans => realization%nw_trans
  option => realization%option
  ncomp = nw_trans%params%ncomp
  nphase = option%transport%nphase
  
  patch%aux%NWT => NWTAuxCreate(ncomp,nphase)
  
  cur_material_property => realization%material_properties
  do                                      
    if (.not.associated(cur_material_property)) exit
    if (maxval(cur_material_property%dispersivity(2:3)) > 0.d0) then
      nw_trans%params%calculate_transverse_dispersion = PETSC_TRUE
      exit
    endif
    cur_material_property => cur_material_property%next
  enddo
  
  material_auxvars => patch%aux%Material%auxvars
  flag = 0
  !TODO(geh): change to looping over ghosted ids once the legacy code is 
  !           history and the communicator can be passed down.
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)

    ! Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    
    if (material_auxvars(ghosted_id)%volume < 0.d0 .and. flag(1) == 0) then
      flag(1) = 1
      option%io_buffer = 'Non-initialized cell volume.'
      call printMsg(option)
    endif
    if (material_auxvars(ghosted_id)%porosity < 0.d0 .and. flag(2) == 0) then
      flag(2) = 1
      option%io_buffer = 'Non-initialized porosity.'
      call printMsg(option)
    endif
    if (material_auxvars(ghosted_id)%tortuosity < 0.d0 .and. flag(3) == 0) then
      flag(3) = 1
      option%io_buffer = 'Non-initialized tortuosity.'
      call printMsg(option)
    endif 
  
  enddo 
  
  if (maxval(flag) > 0) then
    option%io_buffer = &
      'Material property errors found in NWTSetup (Nuclear Waste Transport).'
    call printErrMsg(option)
  endif
  
  ! jenn:todo Should we make this compatible with secondary continuum?
  ! Look at reactive_transport.F90 lines 254-271.
  
  ! allocate auxvar data structures for all grid cells
#ifdef COMPUTE_INTERNAL_MASS_FLUX
  option%iflag = 1 ! allocate mass_balance array
#else  
  option%iflag = 0 ! be sure not to allocate mass_balance array
#endif
  allocate(patch%aux%NWT%auxvars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call NWTAuxVarInit(patch%aux%NWT%auxvars(ghosted_id),nw_trans,option)
  enddo
  patch%aux%NWT%num_aux = grid%ngmax
  
  ! count the number of boundary connections and allocate
  ! auxvar data structures for them
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_condition_list)
  if (sum_connection > 0) then
    option%iflag = 1 ! enable allocation of mass_balance array 
    allocate(patch%aux%NWT%auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call NWTAuxVarInit(patch%aux%NWT%auxvars_bc(iconn),nw_trans,option)
    enddo
  endif
  patch%aux%NWT%num_aux_bc = sum_connection
  
  ! count the number of source/sink connections and allocate
  ! auxvar data structures for them
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (sum_connection > 0) then
    option%iflag = 1 ! enable allocation of mass_balance array 
    allocate(patch%aux%NWT%auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call NWTAuxVarInit(patch%aux%NWT%auxvars_ss(iconn),nw_trans,option)
    enddo
  endif
  patch%aux%NWT%num_aux_ss = sum_connection
  option%iflag = 0
  
  ! initialize parameters
  allocate(nw_trans%diffusion_coefficient(ncomp,nphase))
  nw_trans%diffusion_coefficient = 1.d-9
  allocate(nw_trans%diffusion_activation_energy(ncomp,nphase))
  nw_trans%diffusion_activation_energy = 0.d0
  allocate(nw_trans%species_print(ncomp))
  nw_trans%species_print = PETSC_FALSE
  
  cur_fluid_property => realization%fluid_properties
  do 
    if (.not.associated(cur_fluid_property)) exit
    iphase = cur_fluid_property%phase_id
    ! setting of phase diffusion coefficients must come before individual
    ! species below
    ! jenn:todo Problem: nw_trans%diffusion_coefficient is zero! Is that right?
    if (iphase <= nphase) then
      nw_trans%diffusion_coefficient(:,iphase) = &
        cur_fluid_property%diffusion_coefficient
      nw_trans%diffusion_activation_energy(:,iphase) = &
        cur_fluid_property%diffusion_activation_energy
    endif
    cur_fluid_property => cur_fluid_property%next
  enddo
  
  ! jenn:todo Will we support species-dependent diffusion coefficients?
  ! If so, look at reactive_transport.F90, beginning line 338 in RTSetup().
  
  ! setup output
  list => realization%output_option%output_snap_variable_list
  call PMNWTSetPlotVariables(list,nw_trans,option, &
                             realization%output_option%tunit)
  if (.not.associated(realization%output_option%output_snap_variable_list, &
                      realization%output_option%output_obs_variable_list)) then
    list => realization%output_option%output_obs_variable_list
    call PMNWTSetPlotVariables(list,nw_trans,option, &
                            realization%output_option%tunit)
  endif
  
end subroutine NWTSetup



end module NW_Transport_module
