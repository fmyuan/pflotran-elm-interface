module Hydrate_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Hydrate_Aux_module
  use Hydrate_Common_module
  use Global_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

  public :: HydrateSetup, &
            HydrateInitializeTimestep, &
            HydrateUpdateSolution, &
            HydrateTimeCut,&
            HydrateUpdateAuxVars, &
            HydrateUpdateFixedAccum, &
            HydrateComputeMassBalance, &
            HydrateResidual, &
            HydrateJacobian, &
            HydrateGetTecplotHeader, &
            HydrateSetPlotVariables, &
            HydrateMapBCAuxVarsToGlobal, &
            HydrateDestroy

contains

! ************************************************************************** !

subroutine HydrateSetup(realization)
  ! 
  ! Creates arrays for auxiliary variables
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
  use Fluid_module
  use Material_Aux_class
  use Output_Aux_module
  use Matrix_Zeroing_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(output_variable_list_type), pointer :: list
  type(coupler_type), pointer :: boundary_condition
  type(material_parameter_type), pointer :: material_parameter

  PetscInt :: ghosted_id, iconn, sum_connection, local_id
  PetscInt :: i, idof, count, ndof
  PetscBool :: error_found
  PetscInt :: flag(10)
  PetscErrorCode :: ierr
                                                ! extra index for derivatives
  type(hydrate_auxvar_type), pointer :: hyd_auxvars(:,:)
  type(hydrate_auxvar_type), pointer :: hyd_auxvars_bc(:)
  type(hydrate_auxvar_type), pointer :: hyd_auxvars_ss(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(fluid_property_type), pointer :: cur_fluid_property

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  
  patch%aux%Hydrate => HydrateAuxCreate(option)

  hydrate_analytical_derivatives = .not.option%flow%numerical_derivatives

  ! ensure that material properties specific to this module are properly
  ! initialized
  material_parameter => patch%aux%Material%material_parameter
  error_found = PETSC_FALSE
  
  if (minval(material_parameter%soil_heat_capacity(:)) < 0.d0) then
    option%io_buffer = 'ERROR: Non-initialized soil heat capacity.'
    call PrintMsgByRank(option)
    error_found = PETSC_TRUE
  endif
  if (minval(material_parameter%soil_thermal_conductivity(:,:)) < 0.d0) then
    option%io_buffer = 'ERROR: Non-initialized soil thermal conductivity.'
    call PrintMsgByRank(option)
    error_found = PETSC_TRUE
  endif
  
  material_auxvars => patch%aux%Material%auxvars
  flag = 0
  !TODO(geh): change to looping over ghosted ids once the legacy code is 
  !           history and the communicator can be passed down.
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    if (material_auxvars(ghosted_id)%volume < 0.d0 .and. flag(1) == 0) then
      flag(1) = 1
      option%io_buffer = 'ERROR: Non-initialized cell volume.'
      call PrintMsgByRank(option)
    endif
    if (material_auxvars(ghosted_id)%porosity < 0.d0 .and. flag(2) == 0) then
      flag(2) = 1
      option%io_buffer = 'ERROR: Non-initialized porosity.'
      call PrintMsgByRank(option)
    endif
    if (material_auxvars(ghosted_id)%tortuosity < 0.d0 .and. flag(3) == 0) then
      flag(3) = 1
      option%io_buffer = 'ERROR: Non-initialized tortuosity.'
      call PrintMsgByRank(option)
    endif
    if (material_auxvars(ghosted_id)%soil_particle_density < 0.d0 .and. &
        flag(4) == 0) then
      flag(4) = 1
      option%io_buffer = 'ERROR: Non-initialized soil particle density.'
      call PrintMsgByRank(option)
    endif
    if (minval(material_auxvars(ghosted_id)%permeability) < 0.d0 .and. &
        flag(5) == 0) then
      option%io_buffer = 'ERROR: Non-initialized permeability.'
      call PrintMsgByRank(option)
      flag(5) = 1
    endif
  enddo
  
  error_found = error_found .or. (maxval(flag) > 0)
  call MPI_Allreduce(MPI_IN_PLACE,error_found,ONE_INTEGER_MPI,MPI_LOGICAL, &
                     MPI_LOR,option%mycomm,ierr)
  if (error_found) then
    option%io_buffer = 'Material property errors found in HydrateSetup.'
    call PrintErrMsg(option)
  endif
  
  ! allocate auxvar data structures for all grid cells  
  ndof = option%nflowdof
  allocate(hyd_auxvars(0:2*ndof,grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    do idof = 0, 2 * ndof
      call HydrateAuxVarInit(hyd_auxvars(idof,ghosted_id), &
                         (hydrate_analytical_derivatives .and. idof==0), &
                          option)
    enddo
  enddo
  patch%aux%Hydrate%auxvars => hyd_auxvars
  patch%aux%Hydrate%num_aux = grid%ngmax

  ! count the number of boundary connections and allocate
  ! auxvar data structures for them 
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_condition_list)
  if (sum_connection > 0) then
    allocate(hyd_auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call HydrateAuxVarInit(hyd_auxvars_bc(iconn),PETSC_FALSE,option)
    enddo
    patch%aux%Hydrate%auxvars_bc => hyd_auxvars_bc
  endif
  patch%aux%Hydrate%num_aux_bc = sum_connection

  ! count the number of source/sink connections and allocate
  ! auxvar data structures for them  
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (sum_connection > 0) then
    allocate(hyd_auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call HydrateAuxVarInit(hyd_auxvars_ss(iconn),PETSC_FALSE,option)
    enddo
    patch%aux%Hydrate%auxvars_ss => hyd_auxvars_ss
  endif
  patch%aux%Hydrate%num_aux_ss = sum_connection

  ! initialize parameters
  cur_fluid_property => realization%fluid_properties
  do 
    if (.not.associated(cur_fluid_property)) exit
    patch%aux%Hydrate%hydrate_parameter% &
      diffusion_coefficient(cur_fluid_property%phase_id) = &
        cur_fluid_property%diffusion_coefficient
    cur_fluid_property => cur_fluid_property%next
  enddo  
  ! check whether diffusion coefficients are initialized.
  if (Uninitialized(patch%aux%Hydrate%hydrate_parameter% &
      diffusion_coefficient(LIQUID_PHASE))) then
    option%io_buffer = &
      UninitializedMessage('Liquid phase diffusion coefficient','')
    call PrintErrMsg(option)
  endif
  if (Uninitialized(patch%aux%Hydrate%hydrate_parameter% &
      diffusion_coefficient(GAS_PHASE))) then
    option%io_buffer = &
      UninitializedMessage('Gas phase diffusion coefficient','')
    call PrintErrMsg(option)
  endif

  list => realization%output_option%output_snap_variable_list
  call HydrateSetPlotVariables(realization,list)
  list => realization%output_option%output_obs_variable_list
  call HydrateSetPlotVariables(realization,list)
  
  hydrate_ts_count = 0
  hydrate_ts_cut_count = 0
  hydrate_ni_count = 0


  call PatchSetupUpwindDirection(patch,option)

end subroutine HydrateSetup

! ************************************************************************** !

subroutine HydrateInitializeTimestep(realization)
  ! 
  ! Update data in module prior to time step
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  use Realization_Subsurface_class
  use Upwind_Direction_module

  
  implicit none
  
  type(realization_subsurface_type) :: realization
  
  if (hydrate_restrict_state_chng) then
    realization%patch%aux%Hydrate%auxvars(:,:)%istatechng = PETSC_FALSE
  endif
  
  hydrate_newton_iteration_number = 0
  update_upwind_direction = PETSC_TRUE
  call HydrateUpdateFixedAccum(realization)
  
  hydrate_ni_count = 0

end subroutine HydrateInitializeTimestep

! ************************************************************************** !
subroutine HydrateUpdateSolution(realization)
  ! 
  ! Updates data in module after a successful time
  ! step
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module
  use Grid_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(hydrate_auxvar_type), pointer :: hyd_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  PetscInt :: local_id, ghosted_id
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  hyd_auxvars => patch%aux%Hydrate%auxvars  
  global_auxvars => patch%aux%Global%auxvars
  
  if (realization%option%compute_mass_balance_new) then
    call HydrateUpdateMassBalance(realization)
  endif
  
  ! update stored state
  do ghosted_id = 1, grid%ngmax
    hyd_auxvars(ZERO_INTEGER,ghosted_id)%istate_store(PREV_TS) = &
      global_auxvars(ghosted_id)%istate
  enddo
  hydrate_ts_count = hydrate_ts_count + 1
  hydrate_ts_cut_count = 0
  hydrate_ni_count = 0
 
end subroutine HydrateUpdateSolution

! ************************************************************************** !

subroutine HydrateTimeCut(realization)
  ! 
  ! Resets arrays for time step cut
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Patch_module
  use Discretization_module
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  type(hydrate_auxvar_type), pointer :: hyd_auxvars(:,:)
  
  PetscInt :: local_id, ghosted_id
  PetscErrorCode :: ierr

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  global_auxvars => patch%aux%Global%auxvars
  hyd_auxvars => patch%aux%Hydrate%auxvars

  ! restore stored state
  do ghosted_id = 1, grid%ngmax
    global_auxvars(ghosted_id)%istate = &
      hyd_auxvars(ZERO_INTEGER,ghosted_id)%istate_store(PREV_TS)
  enddo
  
  hydrate_ts_cut_count = hydrate_ts_cut_count + 1

  call HydrateInitializeTimestep(realization)  

end subroutine HydrateTimeCut

! ************************************************************************** !

subroutine HydrateNumericalJacobianTest(xx,realization,B)
  ! 
  ! Computes the a test numerical jacobian
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Grid_module
  use Field_module

  implicit none

  Vec :: xx
  type(realization_subsurface_type) :: realization
  Mat :: B

  Vec :: xx_pert
  Vec :: res
  Vec :: res_pert
  Mat :: A
  PetscViewer :: viewer
  PetscErrorCode :: ierr

  PetscReal, pointer :: vec_p(:), vec2_p(:)

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  PetscReal :: derivative, perturbation
  PetscReal :: perturbation_tolerance = 1.d-6
  PetscInt, save :: icall = 0
  character(len=MAXWORDLENGTH) :: word, word2

  PetscInt :: idof, idof2, icell

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  icall = icall + 1
  call VecDuplicate(xx,xx_pert,ierr);CHKERRQ(ierr)
  call VecDuplicate(xx,res,ierr);CHKERRQ(ierr)
  call VecDuplicate(xx,res_pert,ierr);CHKERRQ(ierr)

  call MatCreate(option%mycomm,A,ierr);CHKERRQ(ierr)
  call MatSetType(A,MATAIJ,ierr);CHKERRQ(ierr)
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,grid%nlmax*option%nflowdof, &
                   grid%nlmax*option%nflowdof, &
                   ierr);CHKERRQ(ierr)
  call MatSeqAIJSetPreallocation(A,27,PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
  call MatSetFromOptions(A,ierr);CHKERRQ(ierr)
  call MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE, &
                    ierr);CHKERRQ(ierr)

  call VecZeroEntries(res,ierr);CHKERRQ(ierr)
  call HydrateResidual(PETSC_NULL_SNES,xx,res,realization,ierr)
#if 0
  word  = 'num_0.dat'
  call PetscViewerASCIIOpen(option%mycomm,word,viewer,ierr);CHKERRQ(ierr)
  call VecView(res,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
  call VecGetArrayF90(res,vec2_p,ierr);CHKERRQ(ierr)
  do icell = 1,grid%nlmax
    if (patch%imat(grid%nL2G(icell)) <= 0) cycle
    do idof = (icell-1)*option%nflowdof+1,icell*option%nflowdof 
      call VecCopy(xx,xx_pert,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(xx_pert,vec_p,ierr);CHKERRQ(ierr)
      perturbation = vec_p(idof)*perturbation_tolerance
      vec_p(idof) = vec_p(idof)+perturbation
      call VecRestoreArrayF90(xx_pert,vec_p,ierr);CHKERRQ(ierr)
      call VecZeroEntries(res_pert,ierr);CHKERRQ(ierr)
      call HydrateResidual(PETSC_NULL_SNES,xx_pert,res_pert,realization,ierr)
#if 0
      write(word,*) idof
      word  = 'num_' // trim(adjustl(word)) // '.dat'
      call PetscViewerASCIIOpen(option%mycomm,word,viewer,ierr);CHKERRQ(ierr)
      call VecView(res_pert,viewer,ierr);CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
      call VecGetArrayF90(res_pert,vec_p,ierr);CHKERRQ(ierr)
      do idof2 = 1, grid%nlmax*option%nflowdof
        derivative = (vec_p(idof2)-vec2_p(idof2))/perturbation
        if (dabs(derivative) > 1.d-30) then
          call MatSetValue(A,idof2-1,idof-1,derivative,INSERT_VALUES, &
                           ierr);CHKERRQ(ierr)
        endif
      enddo
      call VecRestoreArrayF90(res_pert,vec_p,ierr);CHKERRQ(ierr)
    enddo
  enddo
  call VecRestoreArrayF90(res,vec2_p,ierr);CHKERRQ(ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

#if 1
  write(word,*) icall
  word = 'numerical_jacobian-' // trim(adjustl(word)) // '.out'
  call PetscViewerASCIIOpen(option%mycomm,word,viewer,ierr);CHKERRQ(ierr)
  call MatView(A,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

!geh: uncomment to overwrite numerical Jacobian
!  call MatCopy(A,B,DIFFERENT_NONZERO_PATTERN,ierr)
  call MatDestroy(A,ierr);CHKERRQ(ierr)

  call VecDestroy(xx_pert,ierr);CHKERRQ(ierr)
  call VecDestroy(res,ierr);CHKERRQ(ierr)
  call VecDestroy(res_pert,ierr);CHKERRQ(ierr)

end subroutine HydrateNumericalJacobianTest

! ************************************************************************** !

subroutine HydrateComputeMassBalance(realization,mass_balance)
  ! 
  ! Initializes mass balance
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
  use Material_Aux_class
 
  implicit none
  
  type(realization_subsurface_type) :: realization
  PetscReal :: mass_balance(realization%option%nflowspec, &
                            realization%option%nphase)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(hydrate_auxvar_type), pointer :: hydrate_auxvars(:,:)
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscErrorCode :: ierr
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iphase, icomp
  PetscReal :: vol_phase

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  hydrate_auxvars => patch%aux%Hydrate%auxvars
  material_auxvars => patch%aux%Material%auxvars

  mass_balance = 0.d0

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    do iphase = 1, option%nphase
      ! volume_phase = saturation*porosity*volume
      vol_phase = &
        hydrate_auxvars(ZERO_INTEGER,ghosted_id)%sat(iphase)* &
        hydrate_auxvars(ZERO_INTEGER,ghosted_id)%effective_porosity* &
        material_auxvars(ghosted_id)%volume
      ! mass = volume_phase*density
      do icomp = 1, option%nflowspec
        mass_balance(icomp,iphase) = mass_balance(icomp,iphase) + &
          hydrate_auxvars(ZERO_INTEGER,ghosted_id)%den(iphase)* &
          hydrate_auxvars(ZERO_INTEGER,ghosted_id)%xmol(icomp,iphase) * &
          fmw_comp(icomp)*vol_phase
      enddo
    enddo
  enddo

end subroutine HydrateComputeMassBalance

! ************************************************************************** !

subroutine HydrateZeroMassBalanceDelta(realization)
  ! 
  ! Zeros mass balance delta array
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

  do iconn = 1, patch%aux%Hydrate%num_aux_bc
    global_auxvars_bc(iconn)%mass_balance_delta = 0.d0
  enddo
  do iconn = 1, patch%aux%Hydrate%num_aux_ss
    global_auxvars_ss(iconn)%mass_balance_delta = 0.d0
  enddo

end subroutine HydrateZeroMassBalanceDelta

! ************************************************************************** !

subroutine HydrateUpdateMassBalance(realization)
  ! 
  ! Updates mass balance
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  
  PetscInt :: iconn
  PetscInt :: icomp

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

  do iconn = 1, patch%aux%Hydrate%num_aux_bc
    do icomp = 1, option%nflowspec
      global_auxvars_bc(iconn)%mass_balance(icomp,:) = &
        global_auxvars_bc(iconn)%mass_balance(icomp,:) + &
        global_auxvars_bc(iconn)%mass_balance_delta(icomp,:)* &
        fmw_comp(icomp)*option%flow_dt
    enddo
  enddo
  do iconn = 1, patch%aux%Hydrate%num_aux_ss
    do icomp = 1, option%nflowspec
      global_auxvars_ss(iconn)%mass_balance(icomp,:) = &
        global_auxvars_ss(iconn)%mass_balance(icomp,:) + &
        global_auxvars_ss(iconn)%mass_balance_delta(icomp,:)* &
        fmw_comp(icomp)*option%flow_dt
    enddo
  enddo

end subroutine HydrateUpdateMassBalance

! ************************************************************************** !

subroutine HydrateUpdateAuxVars(realization,update_state)
  ! 
  ! Updates the auxiliary variables associated with the Hydrate problem
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  use Material_Aux_class
  use EOS_Water_module
  use Saturation_Function_module
  use Hydrate_Aux_module
  
  implicit none

  type(realization_subsurface_type) :: realization
  PetscBool :: update_state
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(hydrate_auxvar_type), pointer :: hyd_auxvars(:,:), hyd_auxvars_bc(:), &
                                        hyd_auxvars_ss(:)
  type(hydrate_auxvar_type) :: hyd_auxvar, hyd_auxvar_ss
  type(global_auxvar_type) :: global_auxvar_ss, global_auxvar
  type(global_auxvar_type), pointer :: global_auxvars(:), &
                                       global_auxvars_bc(:), global_auxvars_ss(:)
  
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn, natural_id
  PetscInt :: ghosted_start, ghosted_end, i
  PetscInt :: iphasebc, iphase
  PetscInt :: offset
  PetscInt :: istate
  PetscInt :: wat_comp_id, air_comp_id
  PetscReal :: gas_pressure, capillary_pressure, liquid_saturation
  PetscReal :: saturation_pressure, temperature
  PetscReal :: qsrc(3)
  PetscInt :: real_index, variable, flow_src_sink_type
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof), & 
               xxss(realization%option%nflowdof)
  PetscReal :: cell_pressure,qsrc_vol(2),scale
  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  hyd_auxvars => patch%aux%Hydrate%auxvars
  hyd_auxvars_bc => patch%aux%Hydrate%auxvars_bc
  hyd_auxvars_ss => patch%aux%Hydrate%auxvars_ss
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars
    
  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
     
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    ghosted_end = ghosted_id*option%nflowdof
    ghosted_start = ghosted_end - option%nflowdof + 1
    ! HYDRATE_UPDATE_FOR_ACCUM indicates call from non-perturbation
    option%iflag = HYDRATE_UPDATE_FOR_ACCUM
    natural_id = grid%nG2A(ghosted_id)
    if (grid%nG2L(ghosted_id) == 0) natural_id = -natural_id
    call HydrateAuxVarCompute(xx_loc_p(ghosted_start:ghosted_end), &
                       hyd_auxvars(ZERO_INTEGER,ghosted_id), &
                       global_auxvars(ghosted_id), &
                       material_auxvars(ghosted_id), &
                       patch%characteristic_curves_array( &
                         patch%cc_id(ghosted_id))%ptr, &
                       natural_id, &
                       option)
    if (update_state) then
      call HydrateAuxVarUpdateState(xx_loc_p(ghosted_start:ghosted_end), &
                                    hyd_auxvars(ZERO_INTEGER,ghosted_id), &
                                    global_auxvars(ghosted_id), &
                                    material_auxvars(ghosted_id), &
                                    patch%characteristic_curves_array( &
                                      patch%cc_id(ghosted_id))%ptr, &
                                    natural_id, option)
    endif
  enddo

  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      !geh: negate to indicate boundary connection, not actual cell
      natural_id = -grid%nG2A(ghosted_id) 
      offset = (ghosted_id-1)*option%nflowdof
      if (patch%imat(ghosted_id) <= 0) cycle

      xxbc(:) = xx_loc_p(offset+1:offset+option%nflowdof)
      istate = boundary_condition%flow_aux_int_var(HYDRATE_STATE_INDEX,iconn)
      if (istate == HYD_ANY_STATE) then
        istate = global_auxvars(ghosted_id)%istate
        select case(istate)
          case(L_STATE,G_STATE)
            do idof = 1, option%nflowdof
              select case(boundary_condition%flow_bc_type(idof))
                case(DIRICHLET_BC,HYDROSTATIC_BC)
                  real_index = boundary_condition%flow_aux_mapping(dof_to_primary_variable(idof,istate))
                  xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
              end select   
            enddo
          case(HA_STATE) !MAN: Testing HA_STATE
            do idof = 1, option%nflowdof
              select case(boundary_condition%flow_bc_type(idof))
                case(DIRICHLET_BC,HYDROSTATIC_BC)
                  real_index = boundary_condition%flow_aux_mapping( &
                          dof_to_primary_variable(idof,GA_STATE))
                  xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
              end select
            enddo
          case(GA_STATE)
            do idof = 1, option%nflowdof
              select case(boundary_condition%flow_bc_type(idof))
                case(HYDROSTATIC_BC)
                  real_index = boundary_condition%flow_aux_mapping(dof_to_primary_variable(idof,istate))
                  xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                case(DIRICHLET_BC)
                  variable = dof_to_primary_variable(idof,istate)
                  select case(variable)
                    ! for gas pressure dof
                    case(HYDRATE_GAS_PRESSURE_INDEX)
                      real_index = boundary_condition%flow_aux_mapping(variable)
                      if (real_index /= 0) then
                        xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                      else
                        option%io_buffer = 'Mixed FLOW_CONDITION "' // &
                          trim(boundary_condition%flow_condition%name) // &
                          '" needs gas pressure defined.'
                        call PrintErrMsg(option)
                      endif
                    ! for air pressure dof
                    case(HYDRATE_AIR_PRESSURE_INDEX)
                      real_index = boundary_condition%flow_aux_mapping(variable)
                      if (real_index == 0) then ! air pressure not found
                        ! if air pressure is not available, let's try temperature 
                        real_index = boundary_condition%flow_aux_mapping(HYDRATE_TEMPERATURE_INDEX)
                        if (real_index /= 0) then
                          temperature = boundary_condition%flow_aux_real_var(real_index,iconn)
                          call EOSWaterSaturationPressure(temperature,saturation_pressure,ierr)
                          ! now verify whether gas pressure is provided through BC
                          if (boundary_condition%flow_bc_type(ONE_INTEGER) == NEUMANN_BC) then
                            gas_pressure = xxbc(ONE_INTEGER)
                          else
                            real_index = boundary_condition%flow_aux_mapping(HYDRATE_GAS_PRESSURE_INDEX)
                            if (real_index /= 0) then
                              gas_pressure = boundary_condition%flow_aux_real_var(real_index,iconn)
                            else
                              option%io_buffer = 'Mixed FLOW_CONDITION "' // &
                                trim(boundary_condition%flow_condition%name) // &
                                '" needs gas pressure defined to calculate air ' // &
                                'pressure from temperature.'
                              call PrintErrMsg(option)
                            endif
                          endif
                          xxbc(idof) = gas_pressure - saturation_pressure
                        else
                          option%io_buffer = 'Cannot find boundary constraint for air pressure.'
                          call PrintErrMsg(option)
                        endif
                      else
                        xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                      endif
                    ! for gas saturation dof
                    case(HYDRATE_GAS_SATURATION_INDEX)
                      real_index = boundary_condition%flow_aux_mapping(variable)
                      if (real_index /= 0) then
                        xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                      else
!geh: should be able to use the saturation within the cell
!                        option%io_buffer = 'Mixed FLOW_CONDITION "' // &
!                          trim(boundary_condition%flow_condition%name) // &
!                          '" needs saturation defined.'
!                        call PrintErrMsg(option)
                      endif
                    case(HYDRATE_TEMPERATURE_INDEX)
                      real_index = boundary_condition%flow_aux_mapping(variable)
                      if (real_index /= 0) then
                        xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
                      else
                        option%io_buffer = 'Mixed FLOW_CONDITION "' // &
                          trim(boundary_condition%flow_condition%name) // &
                          '" needs temperature defined.'
                        call PrintErrMsg(option)
                      endif
                  end select
                case(NEUMANN_BC)
                case default
                  option%io_buffer = 'Unknown BC type in HydrateUpdateAuxVars().'
                  call PrintErrMsg(option)
              end select
            enddo  
        end select
      else
        ! we do this for all BCs; Neumann bcs will be set later
        do idof = 1, option%nflowdof
          if (istate > 3) then
            real_index = boundary_condition%flow_aux_mapping(&
                    dof_to_primary_variable(idof,GA_STATE))
          else
            real_index = boundary_condition%flow_aux_mapping(&
                    dof_to_primary_variable(idof,istate))
          endif
          if (real_index > 0) then
            xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
          else
            option%io_buffer = 'Error setting up boundary condition in HydrateUpdateAuxVars'
            call PrintErrMsg(option)
          endif
        enddo
      endif
          
      ! set this based on data given
      !MAN fix this
      global_auxvars_bc(sum_connection)%istate = istate
      ! HYDRATE_UPDATE_FOR_BOUNDARY indicates call from non-perturbation
      option%iflag = HYDRATE_UPDATE_FOR_BOUNDARY

      call HydrateAuxVarCompute(xxbc,hyd_auxvars_bc(sum_connection), &
                                global_auxvars_bc(sum_connection), &
                                material_auxvars(ghosted_id), &
                                patch%characteristic_curves_array( &
                                  patch%cc_id(ghosted_id))%ptr, &
                                natural_id, &
                                option)
      ! update state and update aux var; this could result in two update to
      ! the aux var as update state updates if the state changes
      call HydrateAuxVarUpdateState(xxbc,hyd_auxvars_bc(sum_connection), &
                                    global_auxvars_bc(sum_connection), &
                                    material_auxvars(ghosted_id), &
                                    patch%characteristic_curves_array( &
                                      patch%cc_id(ghosted_id))%ptr, &
                                    natural_id,option)
    enddo
    boundary_condition => boundary_condition%next
  enddo
  
  wat_comp_id = option%water_id
  air_comp_id = option%air_id
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do
  
    if (.not.associated(source_sink)) exit
      
    qsrc = source_sink%flow_condition%hydrate%rate%dataset%rarray(:)
    cur_connection_set => source_sink%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle

      flow_src_sink_type = source_sink%flow_condition%hydrate%rate%itype

      global_auxvar = global_auxvars(ghosted_id)
      hyd_auxvar = hyd_auxvars(ZERO_INTEGER, ghosted_id)
      hyd_auxvar_ss = hyd_auxvars_ss(sum_connection)
      global_auxvar_ss = global_auxvars_ss(sum_connection)
    
      flow_src_sink_type = source_sink%flow_condition%hydrate%rate%itype
    
      if (associated(source_sink%flow_condition%hydrate%temperature)) then
        hyd_auxvar_ss%temp = source_sink%flow_condition%hydrate% &
                           temperature%dataset%rarray(1)
      else
        hyd_auxvar_ss%temp = hyd_auxvar%temp
      endif
    
      ! Check if liquid pressure is set
      if (associated(source_sink%flow_condition%hydrate%liquid_pressure)) then
        hyd_auxvar_ss%pres(wat_comp_id) = source_sink%flow_condition% &
                                    hydrate%liquid_pressure%dataset%rarray(1)
      else
        hyd_auxvar_ss%pres(wat_comp_id) = hyd_auxvar%pres(option%liquid_phase)
      endif
    
      ! Check if gas pressure is set
      if (associated(source_sink%flow_condition%hydrate%gas_pressure)) then
        hyd_auxvar_ss%pres(air_comp_id) = source_sink%flow_condition% &
                               hydrate%gas_pressure%dataset%rarray(1)
      else
        hyd_auxvar_ss%pres(air_comp_id) = hyd_auxvar%pres(option%gas_phase)
      endif
    
      select case(flow_src_sink_type)
      case(MASS_RATE_SS)
        qsrc_vol(air_comp_id) = qsrc(air_comp_id)/(fmw_comp(air_comp_id)* &
                              hyd_auxvar%den(air_comp_id))
        qsrc_vol(wat_comp_id) = qsrc(wat_comp_id)/(fmw_comp(wat_comp_id)* &
                              hyd_auxvar%den(wat_comp_id))
      case(SCALED_MASS_RATE_SS)                       ! kg/sec -> kmol/sec
        qsrc_vol(air_comp_id) = qsrc(air_comp_id)/(fmw_comp(air_comp_id)* &
                              hyd_auxvar%den(air_comp_id))*scale 
        qsrc_vol(wat_comp_id) = qsrc(wat_comp_id)/(fmw_comp(wat_comp_id)* &
                              hyd_auxvar%den(wat_comp_id))*scale 
      case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
      ! qsrc1 = m^3/sec             ! den = kmol/m^3
        qsrc_vol(air_comp_id) = qsrc(air_comp_id)*scale
        qsrc_vol(wat_comp_id) = qsrc(wat_comp_id)*scale
      end select
    
      xxss(1) = maxval(hyd_auxvar_ss%pres(option% &
                     liquid_phase:option%gas_phase))
      if (dabs(qsrc_vol(wat_comp_id)) < 1.d-40 .and. &
          dabs(qsrc_vol(air_comp_id)) < 1.d-40) then
        xxss(2) = 0.d0
      else
        xxss(2) = qsrc_vol(air_comp_id)/(qsrc_vol(wat_comp_id) &
                + qsrc_vol(air_comp_id))
      endif
      xxss(3) = hyd_auxvar_ss%temp
    
      cell_pressure = maxval(hyd_auxvar%pres(option% &
                           liquid_phase:option%gas_phase))    
    
      if (cell_pressure>xxss(1) .or. qsrc(wat_comp_id)<0 .or. &
         qsrc(air_comp_id)<0.d0) then
        xxss(1) = cell_pressure
        xxss(2) = hyd_auxvar%sat(air_comp_id)
        xxss(3) = hyd_auxvar%temp
      endif
    
      if (dabs(qsrc(wat_comp_id)) > 0.d0 .and. &
          dabs(qsrc(air_comp_id)) > 0.d0) then
        global_auxvar_ss%istate = GA_STATE
      elseif (dabs(qsrc(wat_comp_id)) > 0.d0) then
        global_auxvar_ss%istate = L_STATE
      elseif (dabs(qsrc(air_comp_id)) > 0.d0) then
        global_auxvar_ss%istate = G_STATE
      else
        global_auxvar_ss%istate = GA_STATE
      endif
    
      if (global_auxvar_ss%istate /= global_auxvar%istate) then
        global_auxvar_ss%istate = GA_STATE
      endif
    
      allocate(global_auxvar_ss%m_nacl(1))
      global_auxvar_ss%m_nacl(1) = 0.d0
      option%iflag = HYDRATE_UPDATE_FOR_SS
    
      ! Compute state variables 
      call HydrateAuxVarCompute(xxss,hyd_auxvar_ss, &
                                global_auxvar_ss, &
                                material_auxvars(ghosted_id), &
                                patch%characteristic_curves_array( &
                                patch%cc_id(source_sink%region% &
                                cell_ids(1)))%ptr, &
                                source_sink%region%cell_ids(1), &
                                option)
    enddo
    source_sink => source_sink%next
  enddo
  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  patch%aux%Hydrate%auxvars_up_to_date = PETSC_TRUE

end subroutine HydrateUpdateAuxVars

! ************************************************************************** !

subroutine HydrateUpdateFixedAccum(realization)
  ! 
  ! Updates the fixed portion of the
  ! accumulation term
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Material_Aux_class
  use Hydrate_Aux_module
  use Hydrate_Common_module

  implicit none
  
  type(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(hydrate_auxvar_type), pointer :: hyd_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(material_parameter_type), pointer :: material_parameter

  type(hydrate_parameter_type), pointer :: hydrate_parameter
  PetscInt :: ghosted_id, local_id, local_start, local_end, natural_id
  PetscInt :: imat
  PetscReal, pointer :: xx_p(:)
  PetscReal, pointer :: accum_p(:), accum_p2(:)
  PetscReal :: Jac_dummy(realization%option%nflowdof, &
                         realization%option%nflowdof)
                          
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  hyd_auxvars => patch%aux%Hydrate%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  material_parameter => patch%aux%Material%material_parameter
  hydrate_parameter => patch%aux%Hydrate%hydrate_parameter

  call VecGetArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    natural_id = grid%nG2A(ghosted_id)
    local_end = local_id*option%nflowdof
    local_start = local_end - option%nflowdof + 1
    ! HYDRATE_UPDATE_FOR_FIXED_ACCUM indicates call from non-perturbation
    option%iflag = HYDRATE_UPDATE_FOR_FIXED_ACCUM

     call HydrateAuxVarCompute(xx_p(local_start:local_end), &
                              hyd_auxvars(ZERO_INTEGER,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              patch%characteristic_curves_array( &
                                patch%cc_id(ghosted_id))%ptr, &
                              natural_id, &
                              option)
    call HydrateAccumulation(hyd_auxvars(ZERO_INTEGER,ghosted_id), &
                             global_auxvars(ghosted_id), &
                             material_auxvars(ghosted_id),patch%grid% &
                             z(ghosted_id),grid%z_max_global, &
                             hydrate_parameter,&
                             material_parameter%soil_heat_capacity(imat), &
                             option,accum_p(local_start:local_end), &
                             Jac_dummy,PETSC_FALSE, &
                             local_id == hydrate_debug_cell_id)
  enddo
  
  
  call VecRestoreArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  
end subroutine HydrateUpdateFixedAccum

! ************************************************************************** !

subroutine HydrateResidual(snes,xx,r,realization,ierr)
  ! 
  ! Computes the residual equation
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module

  use Connection_module
  use Grid_module
  use Coupler_module  
  use Debug_module
  use Material_Aux_class
  use Upwind_Direction_module
  use Hydrate_Common_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_subsurface_type) :: realization
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  Mat, parameter :: null_mat = tMat(0)
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(material_parameter_type), pointer :: material_parameter
  type(hydrate_parameter_type), pointer :: hydrate_parameter
  type(hydrate_auxvar_type), pointer :: hyd_auxvars(:,:), hyd_auxvars_bc(:), &
                                        hyd_auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: iconn
  PetscInt :: iphase
  PetscReal :: scale
  PetscReal :: ss_flow_vol_flux(realization%option%nphase)
  PetscInt :: sum_connection
  PetscInt :: local_start, local_end
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
  PetscInt :: i, imat, imat_up, imat_dn
  PetscInt, save :: iplot = 0
  PetscInt :: flow_src_sink_type

  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum_p(:), accum_p2(:)
  PetscReal, pointer :: vec_p(:)
  
  PetscReal :: qsrc(3)
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word

  PetscInt :: icc_up, icc_dn
  PetscReal :: Res(realization%option%nflowdof)
  PetscReal :: Jac_dummy(realization%option%nflowdof, &
                         realization%option%nflowdof)
  PetscReal :: v_darcy(realization%option%nphase)
  

  discretization => realization%discretization
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  hyd_auxvars => patch%aux%Hydrate%auxvars
  hyd_auxvars_bc => patch%aux%Hydrate%auxvars_bc
  hyd_auxvars_ss => patch%aux%Hydrate%auxvars_ss
  hydrate_parameter => patch%aux%Hydrate%hydrate_parameter
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars
  
  
  hydrate_newton_iteration_number = hydrate_newton_iteration_number + 1
  ! bragflo uses the following logic, update when
  !   it == 1, before entering iteration loop
  !   it > 1 and mod(it-1,frequency) == 0
  ! the first is set in HydrateInitializeTimestep, the second is set here
  if (hydrate_newton_iteration_number > 1 .and. &
      mod(hydrate_newton_iteration_number-1, &
          upwind_dir_update_freq) == 0) then
    update_upwind_direction = PETSC_TRUE
  endif

  ! Communication -----------------------------------------
  ! These 3 must be called before HydrateUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
 
  ! do update state
  hydrate_high_temp_ts_cut = PETSC_FALSE
  ! MAN: add Newton-TR compatibility
  !hydrate_allow_state_change = PETSC_TRUE
  !hydrate_state_changed = PETSC_FALSE
  !if (hydrate_sub_newton_iter_num > 1 .and. hydrate_using_newtontr) then
  !  ! when newtonTR is active and has inner iterations to re-evaluate the 
  !  ! residual,primary variables must not change. -hdp
  !  hydrate_allow_state_change = PETSC_FALSE
  !endif
                                             ! do update state
  call HydrateUpdateAuxVars(realization,hydrate_allow_state_change)

  ! override flags since they will soon be out of date
  patch%aux%Hydrate%auxvars_up_to_date = PETSC_FALSE 

  ! always assume variables have been swapped; therefore, must copy back
  call VecLockPop(xx,ierr); CHKERRQ(ierr)
  call DiscretizationLocalToGlobal(discretization,field%flow_xx_loc,xx, &
                                   NFLOWDOF)
  call VecLockPush(xx,ierr); CHKERRQ(ierr)

  if (option%compute_mass_balance_new) then
    call HydrateZeroMassBalanceDelta(realization)
  endif

  option%iflag = 1
  ! now assign access pointer to local variables
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)

  ! Accumulation terms ------------------------------------
  ! accumulation at t(k) (doesn't change during Newton iteration)
  call VecGetArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  r_p = -accum_p
  call VecRestoreArrayReadF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

  ! accumulation at t(k+1)
  call VecGetArrayF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    local_end = local_id * option%nflowdof
    local_start = local_end - option%nflowdof + 1
    call HydrateAccumulation(hyd_auxvars(ZERO_INTEGER,ghosted_id), &
                             global_auxvars(ghosted_id), &
                             material_auxvars(ghosted_id),patch%grid% &
                             z(ghosted_id),0.d0,&
                             hydrate_parameter,&
                             material_parameter%soil_heat_capacity(imat), &
                             option,Res,Jac_dummy,&
                             hydrate_analytical_derivatives, &
                             local_id == hydrate_debug_cell_id)
    r_p(local_start:local_end) =  r_p(local_start:local_end) + Res(:)
    accum_p2(local_start:local_end) = Res(:)
  enddo
  call VecRestoreArrayF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)

  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      imat_up = patch%imat(ghosted_id_up) 
      imat_dn = patch%imat(ghosted_id_dn) 
      if (imat_up <= 0 .or. imat_dn <= 0) cycle

      icc_up = patch%cc_id(ghosted_id_up)
      icc_dn = patch%cc_id(ghosted_id_dn)

      call HydrateFlux(hyd_auxvars(ZERO_INTEGER,ghosted_id_up), &
                       global_auxvars(ghosted_id_up), &
                       material_auxvars(ghosted_id_up), &
                       material_parameter%soil_thermal_conductivity(:,imat_up), &
                       hyd_auxvars(ZERO_INTEGER,ghosted_id_dn), &
                       global_auxvars(ghosted_id_dn), &
                       material_auxvars(ghosted_id_dn), &
                       material_parameter%soil_thermal_conductivity(:,imat_dn), &
                       cur_connection_set%area(iconn), &
                       cur_connection_set%dist(:,iconn), &
                       patch%flow_upwind_direction(:,iconn), &
                       hydrate_parameter,option,v_darcy,Res, &
                       Jac_dummy,Jac_dummy, &
                       hydrate_analytical_derivatives, &
                       update_upwind_direction, &
                       count_upwind_direction_flip, &
                       (local_id_up == hydrate_debug_cell_id .or. &
                        local_id_dn == hydrate_debug_cell_id))

      patch%internal_velocities(:,sum_connection) = v_darcy
      if (associated(patch%internal_flow_fluxes)) then
        patch%internal_flow_fluxes(:,sum_connection) = Res(:)
      endif
      
      if (local_id_up > 0) then
        local_end = local_id_up * option%nflowdof
        local_start = local_end - option%nflowdof + 1
        r_p(local_start:local_end) = r_p(local_start:local_end) + Res(:)
      endif
         
      if (local_id_dn > 0) then
        local_end = local_id_dn * option%nflowdof
        local_start = local_end - option%nflowdof + 1
        r_p(local_start:local_end) = r_p(local_start:local_end) - Res(:)
      endif
    enddo

    cur_connection_set => cur_connection_set%next
  enddo    

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      imat_dn = patch%imat(ghosted_id)
      if (imat_dn <= 0) cycle

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      icc_dn = patch%cc_id(ghosted_id)

      call HydrateBCFlux(boundary_condition%flow_bc_type, &
                     boundary_condition%flow_aux_mapping, &
                     boundary_condition%flow_aux_real_var(:,iconn), &
                     hyd_auxvars_bc(sum_connection), &
                     global_auxvars_bc(sum_connection), &
                     hyd_auxvars(ZERO_INTEGER,ghosted_id), &
                     global_auxvars(ghosted_id), &
                     material_auxvars(ghosted_id), &
                     material_parameter%soil_thermal_conductivity(:,imat_dn), &
                     cur_connection_set%area(iconn), &
                     cur_connection_set%dist(:,iconn), &
                     patch%flow_upwind_direction_bc(:,iconn), &
                     hydrate_parameter,option, &
                     v_darcy,Res,Jac_dummy, &
                     hydrate_analytical_derivatives, &
                     update_upwind_direction, &
                     count_upwind_direction_flip, &
                     local_id == hydrate_debug_cell_id)
      patch%boundary_velocities(:,sum_connection) = v_darcy
      if (associated(patch%boundary_flow_fluxes)) then
        patch%boundary_flow_fluxes(:,sum_connection) = Res(:)
      endif
      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_bc(sum_connection)%mass_balance_delta(1:2,1) = &
          global_auxvars_bc(sum_connection)%mass_balance_delta(1:2,1) - &
          Res(1:2)
      endif

      local_end = local_id * option%nflowdof
      local_start = local_end - option%nflowdof + 1
      r_p(local_start:local_end)= r_p(local_start:local_end) - Res(:)

    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first 
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      local_end = local_id * option%nflowdof
      local_start = local_end - option%nflowdof + 1

      if (associated(source_sink%flow_aux_real_var)) then
        scale = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
      else
        scale = 1.d0
      endif

      qsrc=source_sink%flow_condition%hydrate%rate%dataset%rarray(:)
      flow_src_sink_type=source_sink%flow_condition%hydrate%rate%itype
      
      call HydrateSrcSink(option,qsrc,flow_src_sink_type, &
                          hyd_auxvars_ss(sum_connection), &
                          hyd_auxvars(ZERO_INTEGER,ghosted_id), &
                          global_auxvars(ghosted_id), &
                          ss_flow_vol_flux, &
                          scale,Res,Jac_dummy, &
                          hydrate_analytical_derivatives, &
                          local_id == hydrate_debug_cell_id)

      r_p(local_start:local_end) =  r_p(local_start:local_end) - Res(:)

      if (associated(patch%ss_flow_vol_fluxes)) then
        patch%ss_flow_vol_fluxes(:,sum_connection) = ss_flow_vol_flux
      endif      
      if (associated(patch%ss_flow_fluxes)) then
        patch%ss_flow_fluxes(:,sum_connection) = Res(:)
      endif      
      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_ss(sum_connection)%mass_balance_delta(1:2,1) = &
          global_auxvars_ss(sum_connection)%mass_balance_delta(1:2,1) - &
          Res(1:2)
      endif

    enddo
    source_sink => source_sink%next
  enddo
  
  if (patch%aux%Hydrate%inactive_cells_exist) then
    do i=1,patch%aux%Hydrate%matrix_zeroing%n_inactive_rows
      r_p(patch%aux%Hydrate%matrix_zeroing%inactive_rows_local(i)) = 0.d0
    enddo
  endif

  if (hydrate_high_temp_ts_cut) then
    r_p(:) = 1.d20
  endif
  
  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  
  call HydrateSSSandbox(r,null_mat,PETSC_FALSE,grid,material_auxvars, &
                        hyd_auxvars,option)

  ! Mass Transfer
  if (field%flow_mass_transfer /= PETSC_NULL_VEC) then
    ! scale by -1.d0 for contribution to residual.  A negative contribution
    ! indicates mass being added to system.
    !call VecGetArrayF90(field%flow_mass_transfer,vec_p,ierr);CHKERRQ(ierr)
    !call VecRestoreArrayF90(field%flow_mass_transfer,vec_p,ierr);CHKERRQ(ierr)
    call VecAXPY(r,-1.d0,field%flow_mass_transfer,ierr);CHKERRQ(ierr)
  endif                      
                        
  if (Initialized(hydrate_debug_cell_id)) then
    call VecGetArrayReadF90(r, r_p, ierr);CHKERRQ(ierr)
    do local_id = hydrate_debug_cell_id-1, hydrate_debug_cell_id+1
      write(*,'(''  residual   : '',i2,10es12.4)') local_id, &
        r_p((local_id-1)*option%nflowdof+1:(local_id-1)*option%nflowdof+2), &
        r_p(local_id*option%nflowdof)*1.d6
    enddo
    call VecRestoreArrayReadF90(r, r_p, ierr);CHKERRQ(ierr)
  endif
  
  if (realization%debug%vecview_residual) then
    call DebugWriteFilename(realization%debug,string,'Gresidual','', &
                            hydrate_ts_count,hydrate_ts_cut_count, &
                            hydrate_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif
  if (realization%debug%vecview_solution) then
    call DebugWriteFilename(realization%debug,string,'Gxx','', &
                            hydrate_ts_count,hydrate_ts_cut_count, &
                            hydrate_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif

  update_upwind_direction = PETSC_FALSE

end subroutine HydrateResidual

! ************************************************************************** !

subroutine HydrateJacobian(snes,xx,A,B,realization,ierr)
  ! 
  ! Computes the Jacobian
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Connection_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Material_Aux_class
  use Upwind_Direction_module
  use Hydrate_Aux_module

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr

  Mat :: J
  MatType :: mat_type
  PetscReal :: norm
  PetscViewer :: viewer

  PetscInt :: icc_up,icc_dn
  PetscReal :: qsrc, scale
  PetscInt :: imat, imat_up, imat_dn
  PetscInt :: local_id, ghosted_id, natural_id
  PetscInt :: irow
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  Vec, parameter :: null_vec = tVec(0)
  
  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof), &
               Jdn(realization%option%nflowdof,realization%option%nflowdof)
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection 
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity 
  PetscInt, pointer :: zeros(:)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(material_parameter_type), pointer :: material_parameter
  type(hydrate_parameter_type), pointer :: hydrate_parameter
  type(hydrate_auxvar_type), pointer :: hyd_auxvars(:,:), &
                                        hyd_auxvars_bc(:), &
                                        hyd_auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:) 
  class(material_auxvar_type), pointer :: material_auxvars(:)
  
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  hyd_auxvars => patch%aux%Hydrate%auxvars
  hyd_auxvars_bc => patch%aux%Hydrate%auxvars_bc
  hyd_auxvars_ss => patch%aux%Hydrate%auxvars_ss
  hydrate_parameter => patch%aux%Hydrate%hydrate_parameter
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars

  hydrate_force_iteration = PETSC_FALSE

  call MatGetType(A,mat_type,ierr);CHKERRQ(ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  else
    J = A
  endif

  call MatZeroEntries(J,ierr);CHKERRQ(ierr)

  if (.not.hydrate_analytical_derivatives) then
    ! Perturb aux vars
    do ghosted_id = 1, grid%ngmax  ! For each local node do...
      if (patch%imat(ghosted_id) <= 0) cycle
      natural_id = grid%nG2A(ghosted_id)
      call HydrateAuxVarPerturb(hyd_auxvars(:,ghosted_id), &
                                global_auxvars(ghosted_id), &
                                material_auxvars(ghosted_id), &
                                patch%characteristic_curves_array( &
                                  patch%cc_id(ghosted_id))%ptr, &
                                natural_id,option)
    enddo
  endif
  
  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    call HydrateAccumDerivative(hyd_auxvars(:,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              grid%z(ghosted_id),grid%z_max_global, &
                              hydrate_parameter, &
                              material_parameter%soil_heat_capacity(imat), &
                              option,Jup) 
    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                  ADD_VALUES,ierr);CHKERRQ(ierr)
  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call DebugWriteFilename(realization%debug,string,'Gjacobian_accum','', &
                            hydrate_ts_count,hydrate_ts_cut_count, &
                            hydrate_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif


  ! Interior Flux Terms -----------------------------------  
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0    
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      imat_up = patch%imat(ghosted_id_up)
      imat_dn = patch%imat(ghosted_id_dn)
      if (imat_up <= 0 .or. imat_dn <= 0) cycle

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   
   
      icc_up = patch%cc_id(ghosted_id_up)
      icc_dn = patch%cc_id(ghosted_id_dn)
                              
      call HydrateFluxDerivative(hyd_auxvars(:,ghosted_id_up), &
                     global_auxvars(ghosted_id_up), &
                     material_auxvars(ghosted_id_up), &
                     material_parameter%soil_thermal_conductivity(:,imat_up), &
                     hyd_auxvars(:,ghosted_id_dn), &
                     global_auxvars(ghosted_id_dn), &
                     material_auxvars(ghosted_id_dn), &
                     material_parameter%soil_thermal_conductivity(:,imat_dn), &
                     cur_connection_set%area(iconn), &
                     cur_connection_set%dist(:,iconn), &
                     patch%flow_upwind_direction(:,iconn), &
                     hydrate_parameter,option,&
                     Jup,Jdn)
      if (local_id_up > 0) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                      Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                      Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
      if (local_id_dn > 0) then
        Jup = -Jup
        Jdn = -Jdn
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                                      Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                                      Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call DebugWriteFilename(realization%debug,string,'Gjacobian_flux','', &
                            hydrate_ts_count,hydrate_ts_cut_count, &
                            hydrate_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      imat_dn = patch%imat(ghosted_id)
      if (imat_dn <= 0) cycle

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      icc_dn = patch%cc_id(ghosted_id)

      call HydrateBCFluxDerivative(boundary_condition%flow_bc_type, &
                      boundary_condition%flow_aux_mapping, &
                      boundary_condition%flow_aux_real_var(:,iconn), &
                      hyd_auxvars_bc(sum_connection), &
                      global_auxvars_bc(sum_connection), &
                      hyd_auxvars(:,ghosted_id), &
                      global_auxvars(ghosted_id), &
                      material_auxvars(ghosted_id), &
                      material_parameter%soil_thermal_conductivity(:,imat_dn), &
                      cur_connection_set%area(iconn), &
                      cur_connection_set%dist(:,iconn), &
                      patch%flow_upwind_direction_bc(:,iconn), &
                      hydrate_parameter,option, &
                      Jdn)

      Jdn = -Jdn
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call DebugWriteFilename(realization%debug,string,'Gjacobian_bcflux','', &
                            hydrate_ts_count,hydrate_ts_cut_count, &
                            hydrate_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif

  ! Source/sinks
  source_sink => patch%source_sink_list%first 
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      
      if (patch%imat(ghosted_id) <= 0) cycle

      if (associated(source_sink%flow_aux_real_var)) then
        scale = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
      else
        scale = 1.d0
      endif
      
      Jup = 0.d0
      call HydrateSrcSinkDerivative(option,source_sink,hyd_auxvars_ss( &
                        sum_connection), &
                        hyd_auxvars(:,ghosted_id), &
                        global_auxvars(ghosted_id), &
                        scale,Jup)

      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)

    enddo
    source_sink => source_sink%next
  enddo
  
!  call HydrateSSSandbox(null_vec,A,PETSC_TRUE,grid,material_auxvars, &
!                        hyd_auxvars,option)

  if (realization%debug%matview_Jacobian_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call DebugWriteFilename(realization%debug,string,'Gjacobian_srcsink','', &
                            hydrate_ts_count,hydrate_ts_cut_count, &
                            hydrate_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif
  
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  ! zero out isothermal and inactive cells
  if (patch%aux%Hydrate%inactive_cells_exist) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    call MatZeroRowsLocal(A,patch%aux%Hydrate%matrix_zeroing%n_inactive_rows, &
                          patch%aux%Hydrate%matrix_zeroing% &
                            inactive_rows_local_ghosted, &
                          qsrc,PETSC_NULL_VEC,PETSC_NULL_VEC, &
                          ierr);CHKERRQ(ierr)
  endif
  
  if (realization%debug%matview_Jacobian) then
    call DebugWriteFilename(realization%debug,string,'Gjacobian','', &
                            hydrate_ts_count,hydrate_ts_cut_count, &
                            hydrate_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(J,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif
  if (realization%debug%norm_Jacobian) then
    option => realization%option
    call MatNorm(J,NORM_1,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("1 norm: ",es11.4)') norm
    call PrintMsg(option)
    call MatNorm(J,NORM_FROBENIUS,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("2 norm: ",es11.4)') norm
    call PrintMsg(option)
    call MatNorm(J,NORM_INFINITY,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("inf norm: ",es11.4)') norm
    call PrintMsg(option)
  endif

!  call MatView(J,PETSC_VIEWER_STDOUT_WORLD,ierr)

#if 0
  imat = 1
  if (imat == 1) then
    call HydrateNumericalJacobianTest(xx,realization,J) 
  endif
#endif

  ! update after evaluations to ensure zero-based index to match screen output
  hydrate_ni_count = hydrate_ni_count + 1

end subroutine HydrateJacobian

! ************************************************************************** !

function HydrateGetTecplotHeader(realization,icolumn)
  ! 
  ! Returns Hydrate Lite contribution to
  ! Tecplot file header
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 
  
  use Realization_Subsurface_class
  use Option_module
  use Field_module
    
  implicit none
  
  character(len=MAXSTRINGLENGTH) :: HydrateGetTecplotHeader
  type(realization_subsurface_type) :: realization
  PetscInt :: icolumn
  
  character(len=MAXSTRINGLENGTH) :: string, string2
  type(option_type), pointer :: option
  type(field_type), pointer :: field  
  PetscInt :: i

  option => realization%option
  field => realization%field
  
  string = ''
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-T [C]"'')') icolumn
  else
    write(string2,'('',"T [C]"'')')
  endif
  string = trim(string) // trim(string2)
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-P [Pa]"'')') icolumn
  else
    write(string2,'('',"P [Pa]"'')')
  endif
  string = trim(string) // trim(string2)
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-State"'')') icolumn
  else
    write(string2,'('',"State"'')')
  endif
  string = trim(string) // trim(string2)
  
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Sat(l)"'')') icolumn
  else
    write(string2,'('',"Sat(l)"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Sat(g)"'')') icolumn
  else
    write(string2,'('',"Sat(g)"'')')
  endif
  string = trim(string) // trim(string2)
    
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Rho(l)"'')') icolumn
  else
    write(string2,'('',"Rho(l)"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Rho(g)"'')') icolumn
  else
    write(string2,'('',"Rho(g)"'')')
  endif
  string = trim(string) // trim(string2)
    
  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-U(l)"'')') icolumn
  else
    write(string2,'('',"U(l)"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-U(g)"'')') icolumn
  else
    write(string2,'('',"U(g)"'')')
  endif
  string = trim(string) // trim(string2)
  
  do i=1,option%nflowspec
    if (icolumn > -1) then
      icolumn = icolumn + 1
      write(string2,'('',"'',i2,''-Xl('',i2,'')"'')') icolumn, i
    else
      write(string2,'('',"Xl('',i2,'')"'')') i
    endif
    string = trim(string) // trim(string2)
  enddo

  do i=1,option%nflowspec
    if (icolumn > -1) then
      icolumn = icolumn + 1
      write(string2,'('',"'',i2,''-Xg('',i2,'')"'')') icolumn, i
    else
      write(string2,'('',"Xg('',i2,'')"'')') i
    endif
    string = trim(string) // trim(string2)
  enddo
 
  HydrateGetTecplotHeader = string

end function HydrateGetTecplotHeader

! ************************************************************************** !

subroutine HydrateSetPlotVariables(realization,list)
  ! 
  ! Adds variables to be printed to list
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 
  
  use Realization_Subsurface_class
  use Output_Aux_module
  use Variables_module
    
  implicit none
  
  type(realization_subsurface_type) :: realization
  type(output_variable_list_type), pointer :: list

  character(len=MAXWORDLENGTH) :: name, units
  type(output_variable_type), pointer :: output_variable

  if (associated(list%first)) then
    return
  endif
  
  if (list%flow_vars) then
  
    name = 'Liquid Pressure'
    units = 'Pa'
    call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                                LIQUID_PRESSURE)

    name = 'Gas Pressure'
    units = 'Pa'
    call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                                GAS_PRESSURE)

    name = 'Liquid Saturation'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                                LIQUID_SATURATION)
    
    name = 'Gas Saturation'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                                GAS_SATURATION)
    
    name = 'Liquid Density'
    units = 'kg/m^3'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                LIQUID_DENSITY)
    
    name = 'Gas Density'
    units = 'kg/m^3'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                GAS_DENSITY)
    
    name = 'X_g^l'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                LIQUID_MOLE_FRACTION, &
                                realization%option%air_id)
    
    name = 'X_l^l'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                LIQUID_MOLE_FRACTION, &
                                realization%option%water_id)
    
    name = 'X_g^g'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                GAS_MOLE_FRACTION, &
                                realization%option%air_id)
    
    name = 'X_l^g'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                GAS_MOLE_FRACTION, &
                                realization%option%water_id)
  
  endif
  
  if (list%energy_vars) then
  
    name = 'Temperature'
    units = 'C'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                TEMPERATURE)
    
    name = 'Liquid Energy'
    units = 'MJ/kmol'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                LIQUID_ENERGY)
    
    name = 'Gas Energy'
    units = 'MJ/kmol'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                GAS_ENERGY)
    
    name = 'Thermodynamic State'
    units = ''
    output_variable => OutputVariableCreate(name,OUTPUT_DISCRETE,units,STATE)
    output_variable%plot_only = PETSC_TRUE ! toggle output off for observation
    output_variable%iformat = 1 ! integer
    call OutputVariableAddToList(list,output_variable)   
  
  endif
  
end subroutine HydrateSetPlotVariables

! ************************************************************************** !

function HydrateAverageDensity(iphase,istate_up,istate_dn, &
                               density_up,density_dn,dden_up,dden_dn)
  ! 
  ! Averages density, using opposite cell density if phase non-existent
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  implicit none

  PetscInt :: iphase
  PetscInt :: istate_up, istate_dn
  PetscReal :: density_up(:), density_dn(:)
  PetscReal :: dden_up, dden_dn

  PetscReal :: HydrateAverageDensity

  dden_up = 0.d0
  dden_dn = 0.d0
  if (iphase == LIQUID_PHASE) then
    if (istate_up == G_STATE) then
      HydrateAverageDensity = density_dn(iphase)
      dden_dn = 1.d0
    else if (istate_dn == G_STATE) then
      HydrateAverageDensity = density_up(iphase)
      dden_up = 1.d0
    else
      HydrateAverageDensity = 0.5d0*(density_up(iphase)+density_dn(iphase))
      dden_up = 0.5d0
      dden_dn = 0.5d0
    endif
  else if (iphase == GAS_PHASE) then
    if (istate_up == L_STATE) then
      HydrateAverageDensity = density_dn(iphase)
      dden_dn = 1.d0      
    else if (istate_dn == L_STATE) then
      HydrateAverageDensity = density_up(iphase)
      dden_up = 1.d0      
    else
      HydrateAverageDensity = 0.5d0*(density_up(iphase)+density_dn(iphase))
      dden_up = 0.5d0
      dden_dn = 0.5d0      
    endif
  endif

end function HydrateAverageDensity

! ************************************************************************** !

subroutine HydrateSSSandbox(residual,Jacobian,compute_derivative, &
                            grid,material_auxvars,hydrate_auxvars,option)
  ! 
  ! Evaluates source/sink term storing residual and/or Jacobian
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 
#include "petsc/finclude/petscmat.h"
  use petscmat
  use Option_module
  use Grid_module
  use Material_Aux_class, only: material_auxvar_type
  use SrcSink_Sandbox_module
  use SrcSink_Sandbox_Base_class
  
  implicit none

  PetscBool :: compute_derivative
  Vec :: residual
  Mat :: Jacobian
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(hydrate_auxvar_type), pointer :: hydrate_auxvars(:,:)
  
  type(grid_type) :: grid
  type(option_type) :: option
  
  PetscReal, pointer :: r_p(:)
  PetscReal :: res(option%nflowdof)
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  class(srcsink_sandbox_base_type), pointer :: cur_srcsink
  PetscInt :: local_id, ghosted_id, istart, iend, irow, idof
  PetscReal :: res_pert(option%nflowdof)
  PetscReal :: aux_real(10)
  PetscErrorCode :: ierr
  
  if (.not.compute_derivative) then
    call VecGetArrayF90(residual,r_p,ierr);CHKERRQ(ierr)
  endif
  
  cur_srcsink => ss_sandbox_list
  do
    if (.not.associated(cur_srcsink)) exit
    aux_real = 0.d0
    local_id = cur_srcsink%local_cell_id
    ghosted_id = grid%nL2G(local_id)
    res = 0.d0
    Jac = 0.d0
    call HydrateSSSandboxLoadAuxReal(cur_srcsink,aux_real, &
                      hydrate_auxvars(ZERO_INTEGER,ghosted_id),option)
    call cur_srcsink%Evaluate(res,Jac,PETSC_FALSE, &
                              material_auxvars(ghosted_id), &
                              aux_real,option)
    if (compute_derivative) then
      do idof = 1, option%nflowdof
        res_pert = 0.d0
        call HydrateSSSandboxLoadAuxReal(cur_srcsink,aux_real, &
                                    hydrate_auxvars(idof,ghosted_id),option)
        call cur_srcsink%Evaluate(res_pert,Jac,PETSC_FALSE, &
                                  material_auxvars(ghosted_id), &
                                  aux_real,option)
        do irow = 1, option%nflowdof
          Jac(irow,idof) = (res_pert(irow)-res(irow)) / &
                            hydrate_auxvars(idof,ghosted_id)%pert
        enddo
      enddo
      call MatSetValuesBlockedLocal(Jacobian,1,ghosted_id-1,1, &
                                    ghosted_id-1,Jac,ADD_VALUES, &
                                    ierr);CHKERRQ(ierr)
    else
      iend = local_id*option%nflowdof
      istart = iend - option%nflowdof + 1
      r_p(istart:iend) = r_p(istart:iend) - res
    endif
    cur_srcsink => cur_srcsink%next
  enddo
  
  if (.not.compute_derivative) then
    call VecRestoreArrayF90(residual,r_p,ierr);CHKERRQ(ierr)
  endif

end subroutine HydrateSSSandbox

! ************************************************************************** !

subroutine HydrateSSSandboxLoadAuxReal(srcsink,aux_real,hyd_auxvar,option)

  use Option_module
  use SrcSink_Sandbox_Base_class
  use SrcSink_Sandbox_WIPP_Gas_class
  use SrcSink_Sandbox_WIPP_Well_class

  implicit none

  class(srcsink_sandbox_base_type) :: srcsink
  PetscReal :: aux_real(:)
  type(hydrate_auxvar_type) hyd_auxvar
  type(option_type) :: option
  
  aux_real = 0.d0
  select type(srcsink)
    class is(srcsink_sandbox_wipp_gas_type)
      aux_real(WIPP_GAS_WATER_SATURATION_INDEX) = &
        hyd_auxvar%sat(option%liquid_phase)
      aux_real(WIPP_GAS_TEMPERATURE_INDEX) = &
        hyd_auxvar%temp
    class is(srcsink_sandbox_wipp_well_type)
      aux_real(WIPP_WELL_LIQUID_MOBILITY) = &
        hyd_auxvar%mobility(option%liquid_phase)
      aux_real(WIPP_WELL_GAS_MOBILITY) = &
        hyd_auxvar%mobility(option%gas_phase)
      aux_real(WIPP_WELL_LIQUID_PRESSURE) = &
        hyd_auxvar%pres(option%liquid_phase)
      aux_real(WIPP_WELL_GAS_PRESSURE) = &
        hyd_auxvar%pres(option%gas_phase)
      aux_real(WIPP_WELL_LIQUID_ENTHALPY) = &
        hyd_auxvar%H(option%liquid_phase)
      aux_real(WIPP_WELL_GAS_ENTHALPY) = &
        hyd_auxvar%H(option%gas_phase)
      aux_real(WIPP_WELL_XMOL_AIR_IN_LIQUID) = &
        hyd_auxvar%xmol(option%air_id,option%liquid_phase)
      aux_real(WIPP_WELL_XMOL_WATER_IN_GAS) = &
        hyd_auxvar%xmol(option%water_id,option%gas_phase)
      aux_real(WIPP_WELL_LIQUID_DENSITY) = &
        hyd_auxvar%den(option%liquid_phase)
      aux_real(WIPP_WELL_GAS_DENSITY) = &
        hyd_auxvar%den(option%gas_phase)
  end select
  
end subroutine HydrateSSSandboxLoadAuxReal

! ************************************************************************** !

subroutine HydrateMapBCAuxVarsToGlobal(realization)
  ! 
  ! Maps variables in hydrate auxvar to global equivalent.
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Coupler_module
  use Connection_module

  implicit none

  type(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(hydrate_auxvar_type), pointer :: hyd_auxvars_bc(:)  
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)  

  PetscInt :: sum_connection, iconn
  
  option => realization%option
  patch => realization%patch

  if (option%ntrandof == 0) return ! no need to update
  
  hyd_auxvars_bc => patch%aux%Hydrate%auxvars_bc
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      global_auxvars_bc(sum_connection)%sat = &
        hyd_auxvars_bc(sum_connection)%sat
      global_auxvars_bc(sum_connection)%den_kg = &
        hyd_auxvars_bc(sum_connection)%den_kg
      global_auxvars_bc(sum_connection)%temp = &
        hyd_auxvars_bc(sum_connection)%temp
    enddo
    boundary_condition => boundary_condition%next
  enddo
  
end subroutine HydrateMapBCAuxVarsToGlobal

! ************************************************************************** !

subroutine HydrateDestroy(realization)
  ! 
  ! Deallocates variables associated with Hydrate
  ! 
  ! Author: Michael Nole
  ! Date: 07/23/19
  ! 

  use Realization_Subsurface_class

  implicit none

  type(realization_subsurface_type) :: realization
  
  ! place anything that needs to be freed here.
  ! auxvars are deallocated in auxiliary.F90.

end subroutine HydrateDestroy

end module Hydrate_module
