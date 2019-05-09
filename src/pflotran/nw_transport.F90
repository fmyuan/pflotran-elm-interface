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
  
  public :: NWTMaxChange, &
            NWTTimeCut, &
            NWTSetup, &
            NWTUpdateAuxVars, &
            NWTAuxVarCompute, &
            NWTInitializeTimestep, &
            NWTUpdateTransportCoefs, &
            NWTUpdateFixedAccumulation, &
            NWTResidualAccum, &
            NWTResidual
            
contains

! ************************************************************************** !

subroutine NWTMaxChange(realization,dcmax)
  ! 
  ! Computes the maximum change in the solution vector
  ! 
  ! Author: Jenn Frederick
  ! Date: 04/18/2019
  ! 

  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Patch_module
  use Grid_module
  
  implicit none
  
  type(realization_subsurface_type) :: realization
  PetscReal :: dcmax(:)
  
  type(field_type), pointer :: field 
  PetscErrorCode :: ierr
  
  field => realization%field

  dcmax = 0.d0
  
  call VecWAXPY(field%tran_dxx,-1.d0,field%tran_xx,field%tran_yy, &
                ierr);CHKERRQ(ierr)
  
  call VecStrideNormAll(field%tran_dxx,NORM_INFINITY,dcmax,ierr);CHKERRQ(ierr)
      
end subroutine NWTMaxChange

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
  PetscInt :: nspecies, nphase
  
  patch => realization%patch
  grid => patch%grid
  nw_trans => realization%nw_trans
  option => realization%option
  nspecies = nw_trans%params%nspecies
  nphase = option%transport%nphase
  
  patch%aux%NWT => NWTAuxCreate()
  
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
  allocate(nw_trans%diffusion_coefficient(nspecies,nphase))
  nw_trans%diffusion_coefficient = 1.d-9
  allocate(nw_trans%diffusion_activation_energy(nspecies,nphase))
  nw_trans%diffusion_activation_energy = 0.d0
  allocate(nw_trans%species_print(nspecies))
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

    call PetscLogEventBegin(logging%event_nwt_auxvars,ierr);CHKERRQ(ierr)
  
    do ghosted_id = 1, grid%ngmax
      if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
      !geh - Ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle

      offset = (ghosted_id-1)*nw_trans%params%nspecies
      istart = offset + 1
      iend = offset + nw_trans%params%nspecies
      
      nwt_auxvars(ghosted_id)%total_bulk_conc = xx_loc_p(istart:iend)

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

        offset = (ghosted_id-1)*nw_trans%params%nspecies
        istart_loc =  1
        iend_loc = nw_trans%params%nspecies
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
            if (nwt_auxvars_bc(sum_connection)%total_bulk_conc(1) < 1.d-200) then
              nwt_auxvars_bc(sum_connection)%total_bulk_conc = &
                                                          xx_loc_p(istart:iend)
            endif
          case(DIRICHLET_ZERO_GRADIENT_BC)
            if (patch%boundary_velocities(iphase,sum_connection) >= 0.d0) then
                  ! don't need to do anything as the constraint below 
                  ! provides all the concentrations, etc.
              if (nwt_auxvars_bc(sum_connection)%total_bulk_conc(1) < 1.d-200) then
                nwt_auxvars_bc(sum_connection)%total_bulk_conc = &
                                                          xx_loc_p(istart:iend)
              endif
            else
              ! same as zero_gradient below
              nwt_auxvars_bc(sum_connection)%total_bulk_conc = &
                                                          xx_loc_p(istart:iend)
            endif
          case(ZERO_GRADIENT_BC)
            nwt_auxvars_bc(sum_connection)%total_bulk_conc = &
                                                          xx_loc_p(istart:iend)               
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

  PetscReal :: ln_conc(nw_trans%params%nspecies)

  ln_conc = log(nwt_auxvar%total_bulk_conc)

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

      ! jenn:todo Make your own version of TDispersion and return the 
      ! subroutine as it was since I don't need it anymore.
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
                      nw_trans%params%nspecies,nw_trans%params%nphase, &
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
                        nw_trans%params%nspecies,nw_trans%params%nphase, &
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

subroutine NWTResidual(snes,xx,r,realization,ierr)
  ! 
  ! Computes the residual equation
  ! 
  ! Author: Jennifer Frederick
  ! Date: 04/18/2019
  ! 

  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Coupler_module
  use Connection_module
  use Discretization_module
  use Option_module
  use Grid_module
  use Logging_module
  use Debug_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_subsurface_type) :: realization
  PetscReal, pointer :: xx_p(:), log_xx_p(:)
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: r_p(:), fixed_accum_p(:), vec_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: nphase, iphase, nspecies
  PetscInt :: istartall, iendall, offset
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(nw_trans_realization_type), pointer :: nw_trans
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars(:)
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn, sum_connection
  PetscReal :: Res(realization%nw_trans%params%nspecies)
  PetscViewer :: viewer  
  
  character(len=MAXSTRINGLENGTH) :: string

  call PetscLogEventBegin(logging%event_nwt_residual,ierr);CHKERRQ(ierr)

  patch => realization%patch
  field => realization%field
  discretization => realization%discretization
  option => realization%option
  nw_trans => patch%nw_trans
  nwt_auxvars => patch%aux%NWT%auxvars
  nwt_auxvars_ss => patch%aux%NWT%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars
  nphase = nw_trans%params%nphase
  nspecies = nw_trans%params%nspecies

  ! Communication -----------------------------------------
  if (realization%nw_trans%use_log_formulation) then
    ! have to convert the log concentration to non-log form
    call VecGetArrayF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(xx,log_xx_p,ierr);CHKERRQ(ierr)
    xx_p(:) = exp(log_xx_p(:))
    call VecRestoreArrayF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(xx,log_xx_p,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(discretization,field%tran_xx, &
                                     field%tran_xx_loc,NTRANDOF)
  else
    call DiscretizationGlobalToLocal(discretization,xx,field%tran_xx_loc, &
                                     NTRANDOF)
  endif
  
  ! Get pointer to residual Vector data
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%tran_accum, fixed_accum_p, ierr);CHKERRQ(ierr)
  
  !== Accumulation Terms ======================================
  if (.not.option%steady_state) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      ! ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
      
      offset = (local_id-1)*nspecies
      istartall = offset + 1
      iendall = offset + nspecies
      
      call NWTResidualAccum(nwt_auxvars(ghosted_id), &
                            global_auxvars(ghosted_id), &
                            material_auxvars(ghosted_id), &
                            nw_trans,Res)
      r_p(istartall:iendall) = r_p(istartall:iendall) + &
        (Res(1:nspecies) - fixed_accum_p(istartall:iendall))/option%tran_dt
      
    enddo
  endif
  
  !== Source/Sink Terms =======================================
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    do iconn = 1, cur_connection_set%num_connections 
      sum_connection = sum_connection + 1     
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      ! ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
      
      offset = (local_id-1)*nspecies
      istartall = offset + 1
      iendall = offset + nspecies
      
      call NWTResidualSrcSink(nwt_auxvars(ghosted_id), &
                              source_sink,patch,sum_connection, &
                              nw_trans,Res)
      r_p(istartall:iendall) = r_p(istartall:iendall) - Res(1:nspecies)
      
      if (associated(patch%ss_tran_fluxes)) then
        patch%ss_tran_fluxes(:,sum_connection) = - Res(:)
      endif
      if (option%compute_mass_balance_new) then
        ! contribution to boundary 
        iphase = LIQUID_PHASE
        nwt_auxvars_ss(sum_connection)%mass_balance_delta(:,iphase) = &
          nwt_auxvars_ss(sum_connection)%mass_balance_delta(:,iphase) - Res
      endif
    enddo
    source_sink => source_sink%next
  enddo

  !== Decay and Ingrowth ======================================
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    ! ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    
    offset = (local_id-1)*nspecies
    istartall = offset + 1
    iendall = offset + nspecies
    
    call NWTResidualRx(nwt_auxvars(ghosted_id), &
                       material_auxvars(ghosted_id), &
                       nw_trans,Res)
    r_p(istartall:iendall) = r_p(istartall:iendall) - Res(1:nspecies)
    
  enddo
  
  
  ! pass #2 for internal and boundary flux terms
  !call NWTResidualFlux(snes,xx,r,realization,ierr)

  
  ! Restore residual Vector data
  call VecRestoreArrayF90(field%tran_accum, fixed_accum_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)

  if (realization%debug%vecview_residual) then
    string = 'NWTresidual'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%vecview_solution) then
    string = 'NWTxx'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(field%tran_xx,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  
  call PetscLogEventEnd(logging%event_nwt_residual,ierr);CHKERRQ(ierr)

end subroutine NWTResidual

! ************************************************************************** !

subroutine NWTUpdateFixedAccumulation(realization)
  ! 
  ! Computes the fixed portion of the accumulation term in 
  ! the residual function (the accumulation at t=t).
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
  PetscReal, pointer :: xx_p(:), fixed_accum_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: dof_offset, istart, iend
  PetscErrorCode :: ierr
  
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

  call VecGetArrayF90(field%tran_accum, fixed_accum_p, ierr);CHKERRQ(ierr)
    
! Do not use NWTUpdateAuxVars() as it loops over ghosted ids

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    
    ! compute offset in solution vector for first dof in grid cell
    dof_offset = (local_id-1)*nw_trans%params%nspecies
    
    ! calculate range of species
    istart = dof_offset + 1
    iend = dof_offset + nw_trans%params%nspecies

    ! copy primary dependent variable
    nwt_auxvars(ghosted_id)%total_bulk_conc = xx_p(istart:iend)
    
    call NWTAuxVarCompute(nwt_auxvars(ghosted_id), &
                          global_auxvars(ghosted_id), &
                          nw_trans,option)
    call NWTResidualAccum(nwt_auxvars(ghosted_id), &
                          global_auxvars(ghosted_id), &
                          material_auxvars(ghosted_id), &
                          nw_trans,fixed_accum_p(istart:iend)) 
  enddo

  call VecRestoreArrayReadF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%tran_accum,fixed_accum_p,ierr);CHKERRQ(ierr)

end subroutine NWTUpdateFixedAccumulation

! ************************************************************************** !

subroutine NWTResidualAccum(nwt_auxvar,global_auxvar,material_auxvar, &
                            nw_trans,Res)
  ! 
  ! Computes the accumulation term in the residual function.
  ! All residual entries should be in [mol/sec].
  ! 
  ! Author: Jenn Frederick
  ! Date: 04/18/2019
  ! 

  use Option_module

  implicit none
  
  type(nw_transport_auxvar_type) :: nwt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(nw_trans_realization_type), pointer :: nw_trans
  PetscReal :: Res(nw_trans%params%nspecies)
  
  PetscInt :: iphase
  PetscInt :: istart, iend
  PetscReal :: pv
  
  iphase = LIQUID_PHASE
  Res = 0.d0
  
  ! All residual entries for accumulation should be in [mol-species].
  ! Dividing by dt will occur later in NWTResidual.
  
  ! porosity in [m^3-void/m^3-bulk]
  ! volume in [m^3-bulk]
  ! pv in [m^3-void]
  pv = material_auxvar%porosity*material_auxvar%volume
      
  istart = 1
  iend = nw_trans%params%nspecies
  
  ! -- Aqueous-Component ----------------------------------------
  ! saturation in [m^3-liq/m^3-void]
  ! aqueous conc in [mol-species/m^3-liq]
  Res(istart:iend) = pv*global_auxvar%sat(iphase)* &
                     nwt_auxvar%aqueous_eq_conc(:)
                     
  ! -- Precipitated-Component -----------------------------------
  ! mineral volume fraction in [m^3-mnrl/m^3-void]
  ! precipitated conc in [mol-species/m^3-mnrl]
  Res(istart:iend) = Res(istart:iend) + &
                     pv*nwt_auxvar%mnrl_vol_frac(:)* &
                     nwt_auxvar%mnrl_eq_conc(:)
                  
  ! -- Sorbed-Component -----------------------------------------
  ! jenn:todo Add sorbed component into NWTResidualAccum().

end subroutine NWTResidualAccum

! ************************************************************************** !

subroutine NWTResidualSrcSink(nwt_auxvar,source_sink,patch, &
                              sum_connection,nw_trans,Res)
  ! 
  ! Computes the source/sink terms in the residual function.
  ! All residual entries should be in [mol-species/sec].
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/08/2019
  ! 

  use Patch_module
  use Coupler_module

  implicit none
  
  type(nw_transport_auxvar_type) :: nwt_auxvar
  type(coupler_type), pointer :: source_sink
  type(patch_type), pointer :: patch
  PetscInt :: sum_connection
  type(nw_trans_realization_type), pointer :: nw_trans
  PetscReal :: Res(nw_trans%params%nspecies)
  
  PetscInt :: istart, iend, iphase
  PetscReal :: qsrc
  PetscReal :: coef_in, coef_out
  
  iphase = LIQUID_PHASE
  Res = 0.d0
  
  if (associated(patch%ss_flow_vol_fluxes)) then
    ! qsrc = [m^3-liq/sec] 
    qsrc = patch%ss_flow_vol_fluxes(LIQUID_PHASE,sum_connection)
  endif
      
  istart = 1
  iend = nw_trans%params%nspecies
  
  ! -- Aqueous-Component ----------------------------------------
  select case(source_sink%tran_condition%itype)
    case(EQUILIBRIUM_SS)
      ! jenn:todo What is EQUILIBRIUM_SS option?
    case(MASS_RATE_SS)
      ! jenn:todo What is MASS_RATE_SS option?
    case default
      if (qsrc > 0.d0) then ! source of fluid flux
        ! represents inside of the domain
        coef_in = 0.d0
        ! represents outside of the domain
        coef_out = qsrc
      else                  ! sink of fluid flux
        ! represents inside of the domain
        coef_in = qsrc
        ! represents outside of the domain
        coef_out = 0.d0
      endif
  end select
  ! units of coef = [m^3-liq/sec]
  ! units of aqueous_eq_conc = [mol-species/m^3-liq]
  ! units of residual entries = [mol-species/sec]
  Res(istart:iend) = (coef_in*nwt_auxvar%aqueous_eq_conc(:)) + &
                     (coef_out*source_sink%tran_condition% &
                               cur_constraint_coupler%nwt_auxvar% &
                               aqueous_eq_conc(:))
                               
  ! -- Precipitated-Component -----------------------------------
  ! There is no contribution from precipitated components.
  
  ! -- Sorbed-Component -----------------------------------------
  ! There is no contribution from sorbed components.

end subroutine NWTResidualSrcSink

! ************************************************************************** !

subroutine NWTResidualRx(nwt_auxvar,material_auxvar,nw_trans,Res)
  ! 
  ! Computes the decay/ingrowth term in the residual function.
  ! All residual entries should be in [mol/sec].
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/09/2019
  ! 

  use Option_module

  implicit none
  
  type(nw_transport_auxvar_type) :: nwt_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(nw_trans_realization_type), pointer :: nw_trans
  PetscReal :: Res(nw_trans%params%nspecies)
  
  PetscInt :: iphase, parent_id
  PetscReal :: vol
  type(species_type), pointer :: species
  type(radioactive_decay_rxn_type), pointer :: rad_rxn
  
  iphase = LIQUID_PHASE
  Res = 0.d0
  
  ! All residual entries for decay/ingrowth should be in [mol-species].
  
  ! volume in [m^3-bulk]
  vol = material_auxvar%volume
  
  species => nw_trans%species_list
  do 
    if (.not.associated(species)) exit
    
    if (species%radioactive) then
      ! Find the reaction object associated with this species
      rad_rxn => nw_trans%rad_decay_rxn_list
      do
        if (.not.associated(rad_rxn)) exit
        if (rad_rxn%species_id == species%id) exit
        rad_rxn => rad_rxn%next
      enddo
      ! Add in species decay
      Res(species%id) = -(rad_rxn%rate_constant * &
                          nwt_auxvar%total_bulk_conc(species%id))
      
      ! Add in contribution from parent (if exists)
      if (rad_rxn%parent_id > 0.d0) then
        parent_id = rad_rxn%parent_id
        ! Find the reaction object associated with the parent species
        rad_rxn => nw_trans%rad_decay_rxn_list
        do
          if (.not.associated(rad_rxn)) exit
          if (rad_rxn%species_id == parent_id) exit
          rad_rxn => rad_rxn%next
        enddo
        Res(species%id) = Res(species%id) + (rad_rxn%rate_constant * &
                          nwt_auxvar%total_bulk_conc(parent_id))
      endif
    endif
    
    species => species%next
  enddo
                    
end subroutine NWTResidualRx

! ************************************************************************** !

subroutine NWTResidualFlux(snes,xx,r,realization,ierr)
  ! 
  ! Computes the flux terms in the residual function for
  ! nuclear waste transport
  ! 
  ! Author: Jenn Frederick
  ! Date: 02/14/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Transport_module
  use Option_module
  use Field_module
  use Grid_module
  use Connection_module
  use Coupler_module  
  use Debug_module
  
  implicit none

  type :: flux_ptrs
    PetscReal, dimension(:), pointer :: flux_p 
  end type

  type (flux_ptrs), dimension(0:2) :: fluxes
  SNES, intent(in) :: snes
  Vec, intent(inout) :: xx
  Vec, intent(inout) :: r
  type(realization_subsurface_type) :: realization  
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: r_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt, parameter :: iphase = 1
  PetscInt :: i, istart, iend                        
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(nw_trans_realization_type), pointer :: nw_trans
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars(:), nwt_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:) 
  
  PetscReal, pointer :: face_fluxes_p(:)

  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn
  PetscInt :: ghosted_id_up, ghosted_id_dn, local_id_up, local_id_dn
  PetscReal :: fraction_upwind, distance, dist_up, dist_dn
  PetscInt :: axis, side, nlx, nly, nlz, ngx, ngxy, pstart, pend, flux_id
  PetscInt :: direction, max_x_conn, max_y_conn
    
  PetscReal :: coef_up(realization%nw_trans%params%nspecies, &
                       realization%option%transport%nphase)
  PetscReal :: coef_dn(realization%nw_trans%params%nspecies, &
                       realization%option%transport%nphase)
  PetscReal :: Res(realization%nw_trans%params%nspecies)

  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  nwt_auxvars => patch%aux%NWT%auxvars
  nwt_auxvars_bc => patch%aux%NWT%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  
  if (option%compute_mass_balance_new) then
    ! jenn:todo Create NWTZeroMassBalanceDelta(realization)
    !call RTZeroMassBalanceDelta(realization)
  endif
  
  ! Get pointer to Vector data
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
 
  r_p = 0.d0

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

      ! TFluxCoef will eventually be moved to another routine where it should be
      ! called only once per flux interface at the beginning of a transport
      ! time step.
      
      call NWTFluxCoef(nw_trans,&
                       global_auxvars(ghosted_id_up), &
                       global_auxvars(ghosted_id_dn), &
                       option,cur_connection_set%area(iconn), &
                       patch%internal_velocities(:,sum_connection), &
                       cur_connection_set%dist(-1,iconn), &
                       coef_up,coef_dn)
                      
      !call TFlux(rt_parameter, &
      !            rt_auxvars(ghosted_id_up), &
      !            global_auxvars(ghosted_id_up), &
      !            rt_auxvars(ghosted_id_dn), &
      !            global_auxvars(ghosted_id_dn), &
      !            coef_up,coef_dn,option,Res)

#ifdef COMPUTE_INTERNAL_MASS_FLUX
      nwt_auxvars(local_id_up)%mass_balance_delta(:,iphase) = &
        nwt_auxvars(local_id_up)%mass_balance_delta(:,iphase) - Res        
#endif
      if (local_id_up>0) then
        iend = local_id_up*nw_trans%params%nspecies
        istart = iend-nw_trans%params%nspecies+1
        r_p(istart:iend) = r_p(istart:iend) + Res(1:nw_trans%params%nspecies)
      endif
      
      if (local_id_dn>0) then
        iend = local_id_dn*nw_trans%params%nspecies
        istart = iend-nw_trans%params%nspecies+1
        r_p(istart:iend) = r_p(istart:iend) - Res(1:nw_trans%params%nspecies)
      endif

      if (associated(patch%internal_tran_fluxes)) then
        patch%internal_tran_fluxes(1:nw_trans%params%nspecies,iconn) = &
            Res(1:nw_trans%params%nspecies)
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

      if (patch%imat(ghosted_id) <= 0) cycle
      
      !call TFluxCoef(rt_parameter, &
      !            global_auxvars_bc(sum_connection), &
      !            global_auxvars(ghosted_id), &
      !            option,cur_connection_set%area(iconn), &
      !            patch%boundary_velocities(:,sum_connection), &
      !            patch%boundary_tran_coefs(:,:,sum_connection), &
      !            0.5d0, &
      !            coef_up,coef_dn)
                  
      !call TFlux(rt_parameter, &
      !           rt_auxvars_bc(sum_connection), &
      !           global_auxvars_bc(sum_connection), &
      !           rt_auxvars(ghosted_id), &
      !           global_auxvars(ghosted_id), &
      !           coef_up,coef_dn,option,Res)
                  
      iend = local_id*nw_trans%params%nspecies
      istart = iend-nw_trans%params%nspecies+1
      r_p(istart:iend)= r_p(istart:iend) - Res(1:nw_trans%params%nspecies)

      if (option%compute_mass_balance_new) then
      ! contribution to boundary 
        nwt_auxvars_bc(sum_connection)%mass_balance_delta(:,iphase) = &
          nwt_auxvars_bc(sum_connection)%mass_balance_delta(:,iphase) - Res
!     ! contribution to internal 
!       nwt_auxvars(ghosted_id)%mass_balance_delta(:,iphase) = &
!         nwt_auxvars(ghosted_id)%mass_balance_delta(:,iphase) + Res
      endif  
                 
      if (associated(patch%boundary_tran_fluxes)) then
        patch%boundary_tran_fluxes(1:nw_trans%params%nspecies,sum_connection) = &
            Res(1:nw_trans%params%nspecies)
      endif
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Restore vectors
  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
 
end subroutine NWTResidualFlux

! ************************************************************************** !

subroutine NWTFluxCoef(nw_trans, &
                       global_auxvar_up,global_auxvar_dn, &
                       option,area,velocity, &
                       fraction_upwind,T_up,T_dn)
  ! 
  ! Computes flux coefficients for residual terms.
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/08/2019
  ! 

  use Option_module

  implicit none
  
  type(nw_trans_realization_type), pointer :: nw_trans 
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn 
  type(option_type) :: option
  PetscReal :: area
  PetscReal :: velocity(*)
  PetscReal :: fraction_upwind
  PetscReal :: T_up(nw_trans%params%nspecies,nw_trans%params%nphase)
  PetscReal :: T_dn(nw_trans%params%nspecies,nw_trans%params%nphase)

  PetscInt :: icomp, ncomp, iphase
  PetscInt :: unit_n_up, unit_n_dn
  PetscReal :: coef_up(nw_trans%params%nspecies)
  PetscReal :: coef_dn(nw_trans%params%nspecies)
  PetscReal :: q
  
  ncomp = nw_trans%params%ncomp

  T_up(:,:) = 0.d0
  T_dn(:,:) = 0.d0

  do icomp = 1, ncomp  ! loop thru (1) aq, (2) ppt, and (3) sorbed components
  
    if (icomp == 1) then
      q = velocity(LIQUID_PHASE)  ! liquid is the only mobile phase
    else
      q = 0.d0
    endif

    ! upstream weighting
    if (q > 0.d0) then
      unit_n_up = -1
      unit_n_dn = +1
      coef_up(:) = (q*unit_n_up*area)
      coef_dn(:) = (q*unit_n_dn*area)
    else
      unit_n_up = -1
      unit_n_dn = +1
      coef_up(:) = (q*unit_n_up*area)
      coef_dn(:) = (q*unit_n_dn*area)
    endif
 
    T_up(:,iphase) = coef_up
    T_dn(:,iphase) = coef_dn
    
  enddo
    
end subroutine NWTFluxCoef

! ************************************************************************** !

end module NW_Transport_module
