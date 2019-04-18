module NW_Transport_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Global_Aux_module
  use Material_Aux_class
  use NW_Transport_Aux_module
  
  use PFLOTRAN_Constants_module
  
  implicit none
  
  private
  
  PetscReal, parameter :: perturbation_tolerance = 1.d-5
  
  public :: NWTTimeCut, &
            NWTSetup, &
            NWTUpdateAuxVars, &
            NWTAuxVarCompute, &
            NWTInitializeTimestep, &
            NWTUpdateTransportCoefs, &
            NWTUpdateFixedAccumulation
            
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
  ! jenn:todo Calling these routine requires "use PM_NWT" which causes a
  ! circular dependency. Need to figure out how to set up plotting stuff
  ! outside of the PM, or have the PM do it independently in it's Setup().
  !call PMNWTSetPlotVariables(list,nw_trans,option, &
  !                           realization%output_option%tunit)
  if (.not.associated(realization%output_option%output_snap_variable_list, &
                      realization%output_option%output_obs_variable_list)) then
    list => realization%output_option%output_obs_variable_list
  !  call PMNWTSetPlotVariables(list,nw_trans,option, &
  !                             realization%output_option%tunit)
  endif
  
end subroutine NWTSetup

! ************************************************************************** !

subroutine NWTUpdateAuxVars(realization,update_cells,update_bcs)
  ! 
  ! Updates the auxiliary variables associated with
  ! nuclear waste transport mode.
  ! 
  ! Author: Jenn Frederick
  ! Date: 04/02/2019
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Option_module
  use Field_module
  use Logging_module
  use Global_Aux_module
  
  implicit none

  type(realization_subsurface_type) :: realization
  PetscBool :: update_bcs
  PetscBool :: update_cells
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(nw_trans_realization_type), pointer :: nw_trans
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: ghosted_id, local_id, sum_connection, iconn
  PetscInt :: istart, iend, istart_loc, iend_loc
  PetscReal, pointer :: xx_loc_p(:)
  PetscInt, parameter :: iphase = 1
  PetscInt :: offset
  PetscErrorCode :: ierr
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars(:)
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars_bc(:)
  PetscInt, save :: icall
  
  data icall/0/

  option => realization%option
  patch => realization%patch  
  grid => patch%grid
  field => realization%field
  nw_trans => realization%nw_trans
  nwt_auxvars => patch%aux%NWT%auxvars
  nwt_auxvars_bc => patch%aux%NWT%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars

  
  call VecGetArrayReadF90(field%tran_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

  if (update_cells) then

    ! jenn:todo make a logging%event_nwt_auxvars?
    call PetscLogEventBegin(logging%event_rt_auxvars,ierr);CHKERRQ(ierr)
  
    do ghosted_id = 1, grid%ngmax
      if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle

      offset = (ghosted_id-1)*nw_trans%params%ncomp
      istart = offset + 1
      iend = offset + nw_trans%params%ncomp
      
      nwt_auxvars(ghosted_id)%molality = xx_loc_p(istart:iend)

      call NWTAuxVarCompute(nwt_auxvars(ghosted_id), &
                            global_auxvars(ghosted_id), &
                            nw_trans,option)

    enddo

    call PetscLogEventEnd(logging%event_rt_auxvars,ierr);CHKERRQ(ierr)
  endif

  if (update_bcs) then

    call PetscLogEventBegin(logging%event_rt_auxvars_bc,ierr);CHKERRQ(ierr)

    boundary_condition => patch%boundary_condition_list%first
    sum_connection = 0    
    do 
      if (.not.associated(boundary_condition)) exit
      cur_connection_set => boundary_condition%connection_set

      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        
        if (patch%imat(ghosted_id) <= 0) cycle

        offset = (ghosted_id-1)*nw_trans%params%ncomp
        istart_loc =  1
        iend_loc = nw_trans%params%ncomp
        istart = offset + istart_loc
        iend = offset + iend_loc
        
        select case(boundary_condition%tran_condition%itype)
          case(CONCENTRATION_SS,DIRICHLET_BC,NEUMANN_BC)
            ! don't need to do anything as the constraint below provides all
            ! the concentrations, etc.
              
            ! jenn:todo Is this kludge still needed?
            !geh: terrible kludge, but should work for now.
            !geh: the problem is that ...%pri_molal() on first call is 
            !      zero and PETSC_TRUE is passed into 
            !      ReactionEquilibrateConstraint() below for 
            !      use_prev_soln_as_guess.  If the previous solution is 
            !      zero, the code will crash.
            if (nwt_auxvars_bc(sum_connection)%molality(1) < 1.d-200) then
              nwt_auxvars_bc(sum_connection)%molality = xx_loc_p(istart:iend)
            endif
          case(DIRICHLET_ZERO_GRADIENT_BC)
            if (patch%boundary_velocities(iphase,sum_connection) >= 0.d0) then
                  ! don't need to do anything as the constraint below 
                  ! provides all the concentrations, etc.
              if (nwt_auxvars_bc(sum_connection)%molality(1) < 1.d-200) then
                nwt_auxvars_bc(sum_connection)%molality = xx_loc_p(istart:iend)
              endif
            else
              ! same as zero_gradient below
              nwt_auxvars_bc(sum_connection)%molality = xx_loc_p(istart:iend)
            endif
          case(ZERO_GRADIENT_BC)
            nwt_auxvars_bc(sum_connection)%molality = xx_loc_p(istart:iend)               
        end select
          
      enddo ! iconn
      boundary_condition => boundary_condition%next
    enddo

    call PetscLogEventEnd(logging%event_rt_auxvars_bc,ierr);CHKERRQ(ierr)

  endif 

  call VecRestoreArrayReadF90(field%tran_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)
  icall = icall+ 1
  
end subroutine NWTUpdateAuxVars

! ************************************************************************** !

subroutine NWTAuxVarCompute(nwt_auxvar,global_auxvar,nw_trans,option)
  ! 
  ! Do I actually need this? I'm just going to keep it as a placeholder.
  ! 
  ! Author: Jenn Frederick
  ! Date: 04/02/2019
  ! 

  use Option_module
  
  implicit none
  
  type(nw_transport_auxvar_type) :: nwt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(nw_trans_realization_type) :: nw_trans
  type(option_type) :: option

  PetscReal :: ln_conc(nw_trans%params%ncomp)

  ln_conc = log(nwt_auxvar%molality)

end subroutine NWTAuxVarCompute

! ************************************************************************** !

subroutine NWTInitializeTimestep(realization)
  ! 
  ! Author: Jenn Frederick
  ! Date: 04/18/2019
  ! 

  use Realization_Subsurface_class

  type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr
  
  ! copying of solution to tran_yy for temporary storage in case of 
  ! time step cut must be performed here as tran_xx change outside of
  ! transport (e.g. pm_ufd_decay)
  call VecCopy(realization%field%tran_xx,realization%field%tran_yy, &
               ierr);CHKERRQ(ierr)
  call NWTUpdateFixedAccumulation(realization)

end subroutine NWTInitializeTimestep

! ************************************************************************** !

subroutine NWTUpdateTransportCoefs(realization)
  ! 
  ! Calculates coefficients for transport matrix
  ! 
  ! Author: Jenn Frederick
  ! Date: 04/18/2019
  ! 

  use Realization_Subsurface_class
  use Discretization_module
  use Patch_module
  use Connection_module
  use Coupler_module
  use Option_module
  use Field_module  
  use Grid_module  
  use Transport_module

  implicit none
  
  type(realization_subsurface_type) :: realization
  
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(nw_trans_realization_type), pointer :: nw_trans
  PetscInt :: local_id, ghosted_id, ghosted_face_id, id
  
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set  
  PetscInt :: sum_connection, iconn, num_connections
  PetscInt :: ghosted_id_up, ghosted_id_dn, local_id_up, local_id_dn
  PetscReal, allocatable :: cell_centered_Darcy_velocities(:,:)
  PetscReal, allocatable :: cell_centered_Darcy_velocities_ghosted(:,:,:)
  PetscReal ::  local_Darcy_velocities_up(3,2)
  PetscReal ::  local_Darcy_velocities_dn(3,2)
  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: i
  PetscInt :: iphase
  PetscInt :: nphase
  PetscErrorCode :: ierr
    
  option => realization%option
  field => realization%field
  patch => realization%patch
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars
  grid => patch%grid
  nw_trans => realization%nw_trans
  nphase = nw_trans%params%nphase
  
  local_Darcy_velocities_up = UNINITIALIZED_DOUBLE
  local_Darcy_velocities_dn = UNINITIALIZED_DOUBLE

  if (nw_trans%params%calculate_transverse_dispersion) then
    allocate(cell_centered_Darcy_velocities_ghosted(3,nphase, &
                                                    patch%grid%ngmax))
    cell_centered_Darcy_velocities_ghosted = 0.d0
    allocate(cell_centered_Darcy_velocities(3,patch%grid%nlmax))
    do iphase = 1, nphase
      call PatchGetCellCenteredVelocities(patch,iphase, &
                                          cell_centered_Darcy_velocities)
      ! at this point, velocities are at local cell centers; we need 
      ! ghosted too.
      do i=1,3
        call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
        vec_ptr(:) = cell_centered_Darcy_velocities(i,:)
        call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
        call DiscretizationGlobalToLocal(realization%discretization, &
                                         field%work, &
                                         field%work_loc,ONEDOF)
        call VecGetArrayF90(field%work_loc,vec_ptr,ierr);CHKERRQ(ierr)
        cell_centered_Darcy_velocities_ghosted(i,iphase,:) = vec_ptr(:)
        call VecRestoreArrayF90(field%work_loc,vec_ptr,ierr);CHKERRQ(ierr)
      enddo
    enddo
    deallocate(cell_centered_Darcy_velocities)
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

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping   

      if (patch%imat(ghosted_id_up) <= 0 .or.  &
          patch%imat(ghosted_id_dn) <= 0) cycle
          
      ! have to use temporary array since unallocated arrays cannot be
      ! indexed in call to subroutine.
      if (allocated(cell_centered_Darcy_velocities_ghosted)) then
        local_Darcy_velocities_up(:,1:nphase) = &
          cell_centered_Darcy_velocities_ghosted(:,1:nphase,ghosted_id_up)
        local_Darcy_velocities_dn(:,1:nphase) = &
          cell_centered_Darcy_velocities_ghosted(:,1:nphase,ghosted_id_dn)
      endif

      call TDispersion(global_auxvars(ghosted_id_up), &
                      material_auxvars(ghosted_id_up), &
                      local_Darcy_velocities_up, &
                      patch%material_property_array(patch%imat(ghosted_id_up))% &
                        ptr%dispersivity, &
                      global_auxvars(ghosted_id_dn), &
                      material_auxvars(ghosted_id_dn), &
                      local_Darcy_velocities_dn, &
                      patch%material_property_array(patch%imat(ghosted_id_dn))% &
                        ptr%dispersivity, &
                      cur_connection_set%dist(:,iconn), &
                      nw_trans%params%ncomp,nw_trans%params%nphase, &
                      realization,option, &
                      patch%internal_velocities(:,sum_connection), &
                      patch%internal_tran_coefs(:,:,sum_connection))
                                           
    enddo
    cur_connection_set => cur_connection_set%next
  enddo    
  
! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
 
    cur_connection_set => boundary_condition%connection_set
    num_connections = cur_connection_set%num_connections
    do iconn = 1, num_connections
      sum_connection = sum_connection + 1
  
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      if (allocated(cell_centered_Darcy_velocities_ghosted)) then
        local_Darcy_velocities_up(:,1:nphase) = &
          cell_centered_Darcy_velocities_ghosted(:,1:nphase,ghosted_id)
      endif
      
      call TDispersionBC(boundary_condition%tran_condition%itype, &
                        global_auxvars_bc(sum_connection), &
                        global_auxvars(ghosted_id), &
                        material_auxvars(ghosted_id), &
                        local_Darcy_velocities_up, &
                        patch%material_property_array(patch%imat(ghosted_id))% &
                          ptr%dispersivity, &
                        cur_connection_set%dist(:,iconn), &
                        nw_trans%params%ncomp,nw_trans%params%nphase, &
                        realization,option, &
                        patch%boundary_velocities(:,sum_connection), &
                        patch%boundary_tran_coefs(:,:,sum_connection))
    enddo
    boundary_condition => boundary_condition%next
  enddo
  
  if (allocated(cell_centered_Darcy_velocities_ghosted)) &
    deallocate(cell_centered_Darcy_velocities_ghosted)

end subroutine NWTUpdateTransportCoefs

! ************************************************************************** !

subroutine NWTUpdateFixedAccumulation(realization)
  ! 
  ! Computes the derivative??? of the accumulation term in 
  ! the residual function
  ! 
  ! Author: Jenn Frederick
  ! Date: 04/18/2019
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use NW_Transport_Aux_module
  use Option_module
  use Field_module  
  use Grid_module

  implicit none
  
  type(realization_subsurface_type) :: realization
  
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(nw_trans_realization_type), pointer :: nw_trans
  PetscReal, pointer :: xx_p(:), accum_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: dof_offset, istart, iend
  PetscErrorCode :: ierr
  PetscReal :: vol_frac_prim
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  nwt_auxvars => patch%aux%NWT%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  grid => patch%grid
  nw_trans => realization%nw_trans

  ! cannot use tran_xx_loc vector here as it has not yet been updated.
  call VecGetArrayReadF90(field%tran_xx,xx_p, ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%tran_accum, accum_p, ierr);CHKERRQ(ierr)
  
  vol_frac_prim = 1.d0 ! what is this?
  
! Do not use NWTUpdateAuxVars() as it loops over ghosted ids

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    
    ! compute offset in solution vector for first dof in grid cell
    dof_offset = (local_id-1)*nw_trans%params%ncomp
    
    ! calculate range of species
    istart = dof_offset + 1
    iend = dof_offset + nw_trans%params%ncomp

    ! copy primary aqueous species
    nwt_auxvars(ghosted_id)%molality = xx_p(istart:iend)
    
    ! DO NOT RECOMPUTE THE ACTIVITY COEFFICIENTS BEFORE COMPUTING THE
    ! FIXED PORTION OF THE ACCUMULATION TERM - geh
    call NWTAuxVarCompute(nwt_auxvars(ghosted_id), &
                          global_auxvars(ghosted_id), &
                          nw_trans,option)
    call NWTAccumulation(nwt_auxvars(ghosted_id), &
                         global_auxvars(ghosted_id), &
                         material_auxvars(ghosted_id), &
                         nw_trans,option, &
                         accum_p(istart:iend)) 
  enddo

  call VecRestoreArrayReadF90(field%tran_xx,xx_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%tran_accum, accum_p, ierr);CHKERRQ(ierr)

end subroutine NWTUpdateFixedAccumulation

! ************************************************************************** !

subroutine NWTAccumulation(nwt_auxvar,global_auxvar,material_auxvar, &
                           nw_trans,option,Res)
  ! 
  ! Computes the accumulation term in the residual function
  ! 
  ! Author: Jenn Frederick
  ! Date: 04/18/2019
  ! 

  use Option_module

  implicit none
  
  type(nw_transport_auxvar_type) :: nwt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  type(nw_trans_realization_type), pointer :: nw_trans
  PetscReal :: Res(nw_trans%params%ncomp)
  
  PetscInt :: iphase
  PetscInt :: istart, iend
  PetscReal :: psv_t
  
  iphase = 1
  Res = 0.d0
  
  ! units = (mol solute/L water)*(m^3 por/m^3 bulk)*(m^3 water/m^3 por)*
  !         (m^3 bulk)*(1000L water/m^3 water)/(sec) = mol/sec
  ! 1000.d0 converts vol from m^3 -> L
  ! all residual entries should be in mol/sec
  psv_t = material_auxvar%porosity*global_auxvar%sat(iphase)*1000.d0* &
          material_auxvar%volume
  istart = 1
  iend = nw_trans%params%ncomp
  ! Your accumulation term goes here:
  Res(istart:iend) = 0.d0
  !Res(istart:iend) = psv_t*nwt_auxvar%total(:,iphase) 

end subroutine NWTAccumulation

end module NW_Transport_module
