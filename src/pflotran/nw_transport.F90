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
            NWTSetup, &
            NWTUpdateAuxVars, &
            NWTAuxVarCompute, &
            NWTInitializeTimestep, &
            NWTUpdateFixedAccumulation, &
            NWTResidual, &
            NWTJacobian, &
            NWTComputeMassBalance, &
            NWTUpdateMassBalance, &
            NWTDestroy
            
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
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars(:)
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars_bc(:)
  PetscInt, save :: icall
  
  data icall/0/

  option => realization%option
  patch => realization%patch  
  grid => patch%grid
  field => realization%field
  nw_trans => realization%nw_trans
  material_auxvars => patch%aux%Material%auxvars
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
                            material_auxvars(ghosted_id), &
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

subroutine NWTAuxVarCompute(nwt_auxvar,global_auxvar,material_auxvar, &
                            nw_trans,option)
  ! 
  ! Computes the secondary variables from the primary dependent variable.
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/29/2019
  ! 

  use Option_module
  use Material_Aux_class
  
  implicit none
  
  type(nw_transport_auxvar_type) :: nwt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(nw_trans_realization_type) :: nw_trans
  type(option_type) :: option

  PetscReal :: ln_conc(nw_trans%params%nspecies)

  ln_conc = log(nwt_auxvar%total_bulk_conc)
  
  ! aqueous concentration (equilibrium)
  nwt_auxvar%aqueous_eq_conc(:) = nwt_auxvar%total_bulk_conc(:) / &
                     (global_auxvar%sat(LIQUID_PHASE)*material_auxvar%porosity)

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
  PetscInt :: ghosted_id, ghosted_id_up, ghosted_id_dn
  PetscInt :: local_id, local_id_up, local_id_dn
  PetscInt :: nphase, iphase, nspecies
  PetscInt :: istart, iend, offset
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(nw_trans_realization_type), pointer :: nw_trans
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars(:)
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars_ss(:)
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars_bc(:)
  type(coupler_type), pointer :: source_sink
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(connection_set_list_type), pointer :: connection_set_list
  PetscInt :: iconn, sum_connection
  PetscReal :: Res(realization%nw_trans%params%nspecies)
  PetscReal :: Res_up(realization%nw_trans%params%nspecies)
  PetscReal :: Res_dn(realization%nw_trans%params%nspecies)
  PetscViewer :: viewer  
  
  character(len=MAXSTRINGLENGTH) :: string

  call PetscLogEventBegin(logging%event_nwt_residual,ierr);CHKERRQ(ierr)

  patch => realization%patch
  field => realization%field
  discretization => realization%discretization
  option => realization%option
  grid => patch%grid
  nw_trans => patch%nw_trans
  nwt_auxvars => patch%aux%NWT%auxvars
  nwt_auxvars_ss => patch%aux%NWT%auxvars_ss
  nwt_auxvars_bc => patch%aux%NWT%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars
  ! note: there is no patch%aux%Material%auxvars_bc
  material_auxvars_bc => patch%aux%Material%auxvars
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
  
#if 1
  !== Accumulation Terms ======================================
  if (.not.option%steady_state) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      ! ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
            
      call NWTResidualAccum(nwt_auxvars(ghosted_id), &
                            global_auxvars(ghosted_id), &
                            material_auxvars(ghosted_id), &
                            nw_trans,Res)
                            
      offset = (local_id-1)*nspecies
      istart = offset + 1
      iend = offset + nspecies
      r_p(istart:iend) = r_p(istart:iend) + &
        (Res(1:nspecies) - fixed_accum_p(istart:iend))/option%tran_dt
      
    enddo
  endif
#endif
  
#if 1
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
            
      call NWTResidualSrcSink(nwt_auxvars(ghosted_id), &
                              source_sink,patch,sum_connection, &
                              nw_trans,Res)
      
      offset = (local_id-1)*nspecies
      istart = offset + 1
      iend = offset + nspecies
      r_p(istart:iend) = r_p(istart:iend) - Res(1:nspecies)
      
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
#endif

#if 1
  !== Decay and Ingrowth ======================================
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    ! ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
        
    call NWTResidualRx(nwt_auxvars(ghosted_id), &
                       material_auxvars(ghosted_id), &
                       nw_trans,Res)
    
    offset = (local_id-1)*nspecies
    istart = offset + 1
    iend = offset + nspecies
    r_p(istart:iend) = r_p(istart:iend) - Res(1:nspecies)
    
  enddo
#endif

#if 0
  !== Fluxes ==================================================
  if (option%compute_mass_balance_new) then
    ! jenn:todo Create NWTZeroMassBalanceDelta(realization)
    !call RTZeroMassBalanceDelta(realization)
  endif
  
  ! Interior Flux Terms ---------------------------------------
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
      local_id_dn = grid%nG2L(ghosted_id_dn) ! ghost to local mapping
      
      ! ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id_up) <= 0 .or.  &
          patch%imat(ghosted_id_dn) <= 0) cycle
          
      call NWTResidualFlux(nwt_auxvars(ghosted_id_up), &
                           nwt_auxvars(ghosted_id_dn), &
                           global_auxvars(ghosted_id_up), &
                           global_auxvars(ghosted_id_dn), &
                           material_auxvars(ghosted_id_up), &
                           material_auxvars(ghosted_id_dn), &
                           cur_connection_set%area(iconn), &
                           cur_connection_set%dist(:,iconn), &
                           patch%internal_velocities(:,sum_connection), &
                           nw_trans,option,Res_up,Res_dn)
                            
      if (local_id_up>0) then
        offset = (local_id_up-1)*nspecies
        istart = offset + 1
        iend = offset + nspecies
        r_p(istart:iend) = r_p(istart:iend) + Res_up(1:nspecies)
      endif
      
      if (local_id_dn>0) then
        offset = (local_id_dn-1)*nspecies
        istart = offset + 1
        iend = offset + nspecies
        r_p(istart:iend) = r_p(istart:iend) + Res_dn(1:nspecies)
      endif
      
      if (associated(patch%internal_tran_fluxes)) then
      ! jenn:todo Not sure how to handle internal_tran_fluxes = Res, because
      ! I have a Res_up and Res_dn, not just Res.
      !  patch%internal_tran_fluxes(1:nw_trans%params%nspecies,iconn) = &
      !      Res(1:nw_trans%params%nspecies)
      endif
      
    enddo
    cur_connection_set => cur_connection_set%next
  enddo
  
  ! Boundary Flux Terms ---------------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do
    if (.not.associated(boundary_condition)) exit
    sum_connection = sum_connection + 1
    
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      ! ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
      
      call NWTResidualFlux(nwt_auxvars_bc(sum_connection), &
                           nwt_auxvars(ghosted_id), &
                           global_auxvars_bc(sum_connection), &
                           global_auxvars(ghosted_id), &
                           material_auxvars_bc(sum_connection), &
                           material_auxvars(ghosted_id), &
                           cur_connection_set%area(iconn), &
                           cur_connection_set%dist(:,iconn), &
                           patch%internal_velocities(:,sum_connection), &
                           nw_trans,option,Res_up,Res_dn)
                            
      offset = (local_id-1)*nspecies
      istart = offset + 1
      iend = offset + nspecies
      r_p(istart:iend) = r_p(istart:iend) + Res_dn(1:nspecies)
      ! note: Don't need to worry about Res_up because that is outside of
      ! the domain, and doesn't have a place in r_p.
      
      if (option%compute_mass_balance_new) then
      ! contribution to boundary
      ! jenn:todo Not sure how to handle mass_balance_delta = Res, because
      ! I have a Res_up and Res_dn, not just Res.
      !  nwt_auxvars_bc(sum_connection)%mass_balance_delta(:,iphase) = &
      !    nwt_auxvars_bc(sum_connection)%mass_balance_delta(:,iphase) - Res
      endif  
                 
      if (associated(patch%boundary_tran_fluxes)) then
      ! jenn:todo Not sure how to handle boundary_tran_fluxes = Res, because
      ! I have a Res_up and Res_dn, not just Res.
      !  patch%boundary_tran_fluxes(1:nw_trans%params%nspecies,sum_connection) = &
      !      Res(1:nw_trans%params%nspecies)
      endif
  
    boundary_condition => boundary_condition%next
  enddo
#endif
  
  ! multiply residual by (-1.0) because Newton's Method is (J)(dC) = (-R)
  r_p = -1.0d0*r_p
  
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
                          material_auxvars(ghosted_id), &
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
  
  PetscInt :: istart, iend
  PetscReal :: pv
  
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
  Res(istart:iend) = pv*global_auxvar%sat(LIQUID_PHASE)* &
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

subroutine NWTResidualFlux(nwt_auxvar_up,nwt_auxvar_dn, &
                           global_auxvar_up,global_auxvar_dn, &
                           material_auxvar_up,material_auxvar_dn, &
                           area,dist,velocity,nw_trans,option,Res_up,Res_dn)
  ! 
  ! Computes the flux terms in the residual function.
  ! All residual entries should be in [mol/sec].
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/09/2019
  ! 

  use Option_module
  use Connection_module

  implicit none
  
  type(nw_transport_auxvar_type) :: nwt_auxvar_up, nwt_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  PetscReal :: area, dist(-1:3)
  PetscReal :: velocity(*)
  type(nw_trans_realization_type), pointer :: nw_trans
  type(option_type) :: option
  PetscReal :: Res_up(nw_trans%params%nspecies)
  PetscReal :: Res_dn(nw_trans%params%nspecies)
  
  PetscInt :: unit_n_up, unit_n_dn
  PetscInt :: nspecies
  PetscReal :: q
  PetscReal :: sat_up, sat_dn
  PetscReal :: dist_up, dist_dn
  PetscReal, pointer :: diffusivity_up(:)
  PetscReal, pointer :: diffusivity_dn(:)
  PetscReal, pointer :: molecular_diffusion_up(:)
  PetscReal, pointer :: molecular_diffusion_dn(:)
  PetscReal, pointer :: diffusion_coefficient(:,:)
  PetscReal :: distance_gravity, upwind_weight ! both are dummy variables
  PetscReal :: harmonic_D_over_dist(nw_trans%params%nspecies)
  PetscReal :: diffusive_flux(nw_trans%params%nspecies)
  
  Res_up = 0.d0
  Res_dn = 0.d0
  nspecies = nw_trans%params%nspecies
  
  allocate(molecular_diffusion_up(nspecies))
  allocate(molecular_diffusion_dn(nspecies))
  allocate(diffusivity_up(nspecies))
  allocate(diffusivity_dn(nspecies))
  harmonic_D_over_dist(:) = 0.d0
  diffusive_flux(:) = 0.d0
  
  sat_up = global_auxvar_up%sat(LIQUID_PHASE)
  sat_dn = global_auxvar_dn%sat(LIQUID_PHASE)
  
  diffusion_coefficient => nw_trans%diffusion_coefficient
  molecular_diffusion_up(:) = diffusion_coefficient(:,LIQUID_PHASE)
  molecular_diffusion_dn(:) = diffusion_coefficient(:,LIQUID_PHASE)
  
  ! get dist_up and dist_dn from dist and dummy variables
  call ConnectionCalculateDistances(dist,option%gravity,dist_up, &
                                    dist_dn,distance_gravity, &
                                    upwind_weight)
  
  diffusivity_up(:) = max(sat_up * material_auxvar_up%porosity * &
                          material_auxvar_up%tortuosity * &
                          molecular_diffusion_up(:), 1.d-40)
  diffusivity_dn(:) = max(sat_dn * material_auxvar_dn%porosity * &
                          material_auxvar_dn%tortuosity * &
                          molecular_diffusion_dn(:), 1.d-40)
                          
  ! weighted harmonic average of diffusivity divided by distance
   harmonic_D_over_dist(:) = (diffusivity_up(:)*diffusivity_dn(:))/ &
                      (diffusivity_up(:)*dist_up + diffusivity_dn(:)*dist_dn)
  
  ! All residual entries for flux terms should be in [mol-species].
  
  ! Diffusive fluxes:
  diffusive_flux(:) = harmonic_D_over_dist(:) * &
                      (nwt_auxvar_up%aqueous_eq_conc(:) &
                       - nwt_auxvar_dn%aqueous_eq_conc(:)) 
                       
  ! Note: For dispersion, do a git pull - Glenn updated transport.F90 
  ! When adding dispersion, look at TDispersion() and the routine that
  ! calls it, UpdateTransportCoefs(), because you need to do something 
  ! with the cell centered velocities. Also, the boundary cells may need
  ! their own calculation for dispersion (There is a TDispersionBC).
                
  ! units of q = [m/s]
  q = velocity(LIQUID_PHASE)  ! liquid is the only mobile phase
  ! units of unit_n = [-] unitless
  unit_n_up = -1 
  unit_n_dn = +1
  
  ! upstream weighting
  if (q > 0.d0) then ! q flows from _up to _dn (think: upstream to downstream)
    Res_up(:) = (unit_n_up*area) * &
                 (q*nwt_auxvar_up%aqueous_eq_conc(:) - diffusive_flux(:))
    Res_dn(:) = (unit_n_dn*area) * &
                 (q*nwt_auxvar_up%aqueous_eq_conc(:) - diffusive_flux(:))
  else               ! q flows from _dn to _up (think: downstream to upstream)
    Res_up(:) = (unit_n_up*area) * &
                 (q*nwt_auxvar_dn%aqueous_eq_conc(:) - diffusive_flux(:))
    Res_dn(:) = (unit_n_dn*area) * &
                 (q*nwt_auxvar_dn%aqueous_eq_conc(:) - diffusive_flux(:))
  endif

  deallocate(diffusivity_dn)
  deallocate(diffusivity_up)
  deallocate(molecular_diffusion_dn)
  deallocate(molecular_diffusion_up)
                    
end subroutine NWTResidualFlux

! ************************************************************************** !

subroutine NWTJacobian(snes,xx,A,B,realization,ierr)
  ! 
  ! Computes the Jacobian matrix for Nuclear Waste Transport.
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/14/2019
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module
  use Logging_module
  use Debug_module
  use Connection_module
  use Coupler_module  

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr

  Mat :: J
  MatType :: mat_type
  PetscViewer :: viewer  
  type(option_type), pointer :: option
  type(grid_type),  pointer :: grid
  type(patch_type), pointer :: patch
  type(nw_trans_realization_type), pointer :: nw_trans
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars(:)
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars_bc(:)
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(connection_set_list_type), pointer :: connection_set_list
  type(coupler_type), pointer :: boundary_condition
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: iconn, sum_connection
  PetscInt :: iphase
  PetscReal :: Jac(realization%nw_trans%params%nspecies, &
                   realization%nw_trans%params%nspecies)
  PetscReal :: JacUp(realization%nw_trans%params%nspecies, &
                     realization%nw_trans%params%nspecies)
  PetscReal :: JacDn(realization%nw_trans%params%nspecies, &
                     realization%nw_trans%params%nspecies)
  PetscReal :: rdum
    
  option => realization%option
  patch => realization%patch  
  grid => patch%grid
  nw_trans => realization%nw_trans
  
  nwt_auxvars => patch%aux%NWT%auxvars
  nwt_auxvars_bc => patch%aux%NWT%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars
  ! note: there is no patch%aux%Material%auxvars_bc
  material_auxvars_bc => patch%aux%Material%auxvars
  
  iphase = LIQUID_PHASE

  call PetscLogEventBegin(logging%event_nwt_jacobian,ierr);CHKERRQ(ierr)

  call MatGetType(A,mat_type,ierr);CHKERRQ(ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  else
    J = A
  endif
    
  call MatZeroEntries(J,ierr);CHKERRQ(ierr)
  
#if 1
  !== Accumulation Terms ======================================
  if (.not.option%steady_state) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      ! ignore inactive cells with inactive materials
      if (patch%imat(ghosted_id) <= 0) cycle
      
      call NWTJacobianAccum(material_auxvars(ghosted_id), &
                            nw_trans,option,Jac) 
              
      ! PETSc uses 0-based indexing so the position must be (ghosted_id-1)              
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jac, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)
  
    enddo
  endif
#endif

#if 1
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
      
      call NWTJacobianSrcSink(material_auxvars(ghosted_id), &
                              global_auxvars(ghosted_id),source_sink, &
                              patch,sum_connection,nw_trans,Jac) 
                                
      ! PETSc uses 0-based indexing so the position must be (ghosted_id-1)
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jac, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)
    
    enddo
        
    source_sink => source_sink%next
  enddo
#endif

#if 1  
  !== Decay and Ingrowth ======================================
  do local_id = 1, grid%nlmax  
    ghosted_id = grid%nL2G(local_id)
    ! ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    
    call NWTJacobianRx(material_auxvars(ghosted_id), &
                       nw_trans,Jac)
    
    ! PETSc uses 0-based indexing so the position must be (ghosted_id-1)              
    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jac, &
                                  ADD_VALUES,ierr);CHKERRQ(ierr)
    
  enddo
#endif
  
#if 0
  !== Fluxes ==================================================
    
  ! Interior Flux Terms ---------------------------------------
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
      local_id_dn = grid%nG2L(ghosted_id_dn) ! ghost to local mapping
      
      if (patch%imat(ghosted_id_up) <= 0 .or.  &
          patch%imat(ghosted_id_dn) <= 0) cycle
          
      call NWTJacobianFlux(nwt_auxvars(ghosted_id_up), &
                           nwt_auxvars(ghosted_id_dn), &
                           global_auxvars(ghosted_id_up), &
                           global_auxvars(ghosted_id_dn), &
                           material_auxvars(ghosted_id_up), &
                           material_auxvars(ghosted_id_dn), &
                           cur_connection_set%area(iconn), &
                           cur_connection_set%dist(:,iconn), &
                           patch%internal_velocities(:,sum_connection), &
                           nw_trans,option,JacUp,JacDn)
          
      if (local_id_up>0) then
        ! PETSc uses 0-based indexing so the position must be (ghosted_id-1)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                      JacUp,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
   
      if (local_id_dn>0) then
        ! PETSc uses 0-based indexing so the position must be (ghosted_id-1)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                                      JacDn,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
      
    enddo
    
  cur_connection_set => cur_connection_set%next
  enddo

  ! Boundary Flux Terms ---------------------------------------
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
          
      call NWTJacobianFlux(nwt_auxvars_bc(sum_connection), &
                           nwt_auxvars(ghosted_id), &
                           global_auxvars_bc(sum_connection), &
                           global_auxvars(ghosted_id), &
                           material_auxvars_bc(sum_connection), &
                           material_auxvars(ghosted_id), &
                           cur_connection_set%area(iconn), &
                           cur_connection_set%dist(:,iconn), &
                           patch%internal_velocities(:,sum_connection), &
                           nw_trans,option,JacUp,JacDn)
      
      ! PETSc uses 0-based indexing so the position must be (ghosted_id-1)              
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,JacDn, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)
      ! note: Don't need to worry about JacUp because that is outside of
      ! the domain, and doesn't have a place in A.
      
    enddo
    
  boundary_condition => boundary_condition%next
  enddo
#endif
    
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)  
    
  if (realization%debug%matview_Jacobian) then
    string = 'NWTjacobian'
    call DebugCreateViewer(realization%debug,string,realization%option,viewer)
    call MatView(J,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  if (realization%nw_trans%use_log_formulation) then
    call MatDiagonalScaleLocal(J,realization%field%tran_work_loc, &
                               ierr);CHKERRQ(ierr)

    if (realization%debug%matview_Jacobian) then
      string = 'NWTjacobianLog'
      call DebugCreateViewer(realization%debug,string,realization%option,viewer)
      call MatView(J,viewer,ierr);CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    endif
    
  endif

  call PetscLogEventEnd(logging%event_nwt_jacobian,ierr);CHKERRQ(ierr)
  
end subroutine NWTJacobian

! ************************************************************************** !

subroutine NWTJacobianAccum(material_auxvar,nw_trans,option,Jac)
  ! 
  ! Computes the accumulation terms in the Jacobian matrix.
  ! All Jacobian entries should be in [m^3-bulk/sec].
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/14/2019
  !                             
  
  use Option_module
  
  implicit none
  
  class(material_auxvar_type) :: material_auxvar
  type(nw_trans_realization_type) :: nw_trans
  type(option_type) :: option
  PetscReal :: Jac(nw_trans%params%nspecies,nw_trans%params%nspecies)
  
  PetscReal :: vol_dt
  PetscInt :: istart, iend, ispecies
  
  Jac = 0.d0
  
  ! units of volume = [m^3-bulk]
  ! units of tran_dt = [sec]
  vol_dt = material_auxvar%volume/option%tran_dt
  
  istart = 1
  iend = nw_trans%params%nspecies
  do ispecies=istart,iend
    ! units of Jac = [m^3-bulk/sec]
    Jac(ispecies,ispecies) = vol_dt
  enddo

end subroutine NWTJacobianAccum

! ************************************************************************** !

subroutine NWTJacobianSrcSink(material_auxvar,global_auxvar,source_sink, &
                              patch,sum_connection,nw_trans,Jac)
  ! 
  ! Computes the source/sink terms in the Jacobian matrix.
  ! All Jacobian entries should be in [m^3-bulk/sec].
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/14/2019
  ! 

  use Patch_module
  use Coupler_module

  implicit none
  
  class(material_auxvar_type) :: material_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(coupler_type), pointer :: source_sink
  type(patch_type), pointer :: patch
  PetscInt :: sum_connection
  type(nw_trans_realization_type), pointer :: nw_trans
  PetscReal :: Jac(nw_trans%params%nspecies,nw_trans%params%nspecies)
  
  PetscInt :: istart, iend, ispecies
  PetscReal :: qsrc, u, vol
  PetscReal :: coef_in
  
  Jac = 0.d0
  
  if (associated(patch%ss_flow_vol_fluxes)) then
    ! qsrc = [m^3-liq/sec] 
    qsrc = patch%ss_flow_vol_fluxes(LIQUID_PHASE,sum_connection)
  endif
      
  istart = 1
  iend = nw_trans%params%nspecies
  
  ! transform qsrc into pore velocity volumetric flow
  ! units of porosity = [m^3-void/m^3-bulk]
  ! units of saturation = [m^3-liq/m^3-void]
  ! units of u = [m^3-bulk/sec]
  u = qsrc / (material_auxvar%porosity * global_auxvar%sat(LIQUID_PHASE))
  
  vol = material_auxvar%volume
  
  ! -- Aqueous-Component ----------------------------------------
  select case(source_sink%tran_condition%itype)
    case(EQUILIBRIUM_SS)
      ! jenn:todo What is EQUILIBRIUM_SS option?
    case(MASS_RATE_SS)
      ! jenn:todo What is MASS_RATE_SS option?
    case default
      ! Note: We only care about coef_in here, because the Jac is a derivative
      ! w.r.t. total_bulk_conc, which only exists in the inside of the domain.
      ! On the outside of the domain, we have a specified conc, which is not a
      ! fn(total_bulk_conc), thus the derivative is zero.
      if (u > 0.d0) then ! source of fluid flux
        ! represents inside of the domain
        coef_in = 0.d0
      else               ! sink of fluid flux
        ! represents inside of the domain
        coef_in = u
      endif
  end select
  
  ! units of coef = [m^3-bulk/sec]
  ! units of volume = [m^3-bulk]
  istart = 1
  iend = nw_trans%params%nspecies
  do ispecies=istart,iend
    ! units of Jac = [m^3-bulk/sec]
    Jac(ispecies,ispecies) = -1.d0 * vol * (coef_in/vol)
    ! Note: I multiply and then divide by volume to be consistent with the 
    ! details provided in the theory guide for this transport mode.
    ! Note: I multiply by -1 because src/sinks are subtracted in the residual.
  enddo
                               
  ! -- Precipitated-Component -----------------------------------
  ! There is no contribution from precipitated components.
  
  ! -- Sorbed-Component -----------------------------------------
  ! There is no contribution from sorbed components.

end subroutine NWTJacobianSrcSink

! ************************************************************************** !

subroutine NWTJacobianRx(material_auxvar,nw_trans,Jac)
  ! 
  ! Computes the radioactive decay/ingrowth terms in the Jacobian matrix.
  ! All Jacobian entries should be in [m^3-bulk/sec].
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/14/2019
  !                             
  
  use Option_module
  
  implicit none
  
  class(material_auxvar_type) :: material_auxvar
  type(nw_trans_realization_type) :: nw_trans
  PetscReal :: Jac(nw_trans%params%nspecies,nw_trans%params%nspecies)
  
  type(radioactive_decay_rxn_type), pointer :: rad_rxn
  PetscReal :: vol
  PetscReal :: decay_rate,parent_decay_rate
  PetscInt :: parent_id, ispecies
  PetscInt :: istart, iend
  PetscBool :: has_parent
  
  Jac = 0.d0
  
  ! units of volume = [m^3-bulk]
  vol = material_auxvar%volume
  
  istart = 1
  iend = nw_trans%params%nspecies
  do ispecies=istart,iend
    decay_rate = 0.d0
    parent_decay_rate = 0.d0
    has_parent = PETSC_FALSE
    ! find the decay rate for species_id = ispecies (if its radioactive)
    rad_rxn => nw_trans%rad_decay_rxn_list
    do
      if (.not.associated(rad_rxn)) exit
      if (rad_rxn%species_id == ispecies) then
        decay_rate = rad_rxn%rate_constant
        parent_id = rad_rxn%parent_id
        ! check if the species has a parent
        if (parent_id > 0) then
          has_parent = PETSC_TRUE
          parent_decay_rate = rad_rxn%rate_constant_parent          
        endif
        exit
      endif
      rad_rxn => rad_rxn%next
    enddo
    
    ! fill in the diagonal of the Jacobian first
    ! units of Jac = [m^3-bulk/sec]
    Jac(ispecies,ispecies) = -1.d0 * (-1.d0*vol*decay_rate)
    ! Note: I multiply by -1 because rx's are subtracted in the residual.
    
    ! fill in the off-diagonal associated with ingrowth from a parent
    if (has_parent) then
      Jac(ispecies,parent_id) = -1.d0 * (vol*parent_decay_rate)
      ! Note: I multiply by -1 because rx's are subtracted in the residual.
    endif
    
  enddo

end subroutine NWTJacobianRx

! ************************************************************************** !

subroutine NWTJacobianFlux(nwt_auxvar_up,nwt_auxvar_dn, &
                           global_auxvar_up,global_auxvar_dn, &
                           material_auxvar_up,material_auxvar_dn, &
                           area,dist,velocity,nw_trans,option,Jac_up,Jac_dn)
  ! 
  ! Computes the flux terms in the Jacobian matrix.
  ! All Jacobian entries should be in [m^3-bulk/sec].
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/14/2019
  ! 

  use Option_module
  use Connection_module

  implicit none
  
  type(nw_transport_auxvar_type) :: nwt_auxvar_up, nwt_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  PetscReal :: area, dist(-1:3)
  PetscReal :: velocity(*) ! at connection, not at cell center
  type(nw_trans_realization_type), pointer :: nw_trans
  type(option_type) :: option
  PetscReal :: Jac_up(nw_trans%params%nspecies,nw_trans%params%nspecies)
  PetscReal :: Jac_dn(nw_trans%params%nspecies,nw_trans%params%nspecies)
  
  PetscInt :: unit_n_up, unit_n_dn
  PetscInt :: nspecies, ispecies
  PetscReal :: q, u
  PetscReal :: harmonic_porosity, harmonic_saturation
  PetscReal :: sat_up, sat_dn
  PetscReal :: por_up, por_dn
  PetscReal :: dist_up, dist_dn
  PetscReal, pointer :: diffusivity_up(:)
  PetscReal, pointer :: diffusivity_dn(:)
  PetscReal, pointer :: molecular_diffusion_up(:)
  PetscReal, pointer :: molecular_diffusion_dn(:)
  PetscReal, pointer :: diffusion_coefficient(:,:)
  PetscReal :: distance_gravity, upwind_weight ! both are dummy variables
  PetscReal :: harmonic_D_over_dist(nw_trans%params%nspecies)
  
  Jac_up = 0.d0
  Jac_dn = 0.d0
  nspecies = nw_trans%params%nspecies
  
  allocate(molecular_diffusion_up(nspecies))
  allocate(molecular_diffusion_dn(nspecies))
  allocate(diffusivity_up(nspecies))
  allocate(diffusivity_dn(nspecies))
  harmonic_D_over_dist(:) = 0.d0
  
  sat_up = global_auxvar_up%sat(LIQUID_PHASE)
  sat_dn = global_auxvar_dn%sat(LIQUID_PHASE)
  por_up = material_auxvar_up%porosity
  por_dn = material_auxvar_dn%porosity
  
  diffusion_coefficient => nw_trans%diffusion_coefficient
  molecular_diffusion_up(:) = diffusion_coefficient(:,LIQUID_PHASE)
  molecular_diffusion_dn(:) = diffusion_coefficient(:,LIQUID_PHASE)
  
  ! get dist_up and dist_dn from dist and dummy variables
  call ConnectionCalculateDistances(dist,option%gravity,dist_up, &
                                    dist_dn,distance_gravity, &
                                    upwind_weight)
                                    
  diffusivity_up(:) = max(material_auxvar_up%tortuosity * &
                          molecular_diffusion_up(:), 1.d-40)
  diffusivity_dn(:) = max(material_auxvar_dn%tortuosity * &
                          molecular_diffusion_dn(:), 1.d-40)
                          
  ! weighted harmonic average of diffusivity divided by distance
   harmonic_D_over_dist(:) = (diffusivity_up(:)*diffusivity_dn(:))/ &
                      (diffusivity_up(:)*dist_up + diffusivity_dn(:)*dist_dn)
                       
  ! Note: For dispersion, do a git pull - Glenn updated transport.F90 
  ! When adding dispersion, look at TDispersion() and the routine that
  ! calls it, UpdateTransportCoefs(), because you need to do something 
  ! with the cell centered velocities. Also, the boundary cells may need
  ! their own calculation for dispersion (There is a TDispersionBC).
                
  ! units of q = [m-liq/s]
  ! units of u = [m-bulk/s]
  q = velocity(LIQUID_PHASE)  ! liquid is the only mobile phase
  ! weighted harmonic average
  harmonic_porosity = (por_up*por_dn)/(por_up*dist_up + por_dn*dist_dn)* &
                      (dist_up+dist_dn)
  harmonic_saturation = (sat_up*sat_dn)/(sat_up*dist_up + sat_dn*dist_dn)* &
                        (dist_up+dist_dn)
  u = q / (harmonic_porosity*harmonic_saturation)
  
  ! units of unit_n = [-] unitless
  unit_n_up = -1 
  unit_n_dn = +1
  
  do ispecies=1,nspecies
    Jac_up(ispecies,ispecies) = (unit_n_up*area) * &
                                (u - harmonic_D_over_dist(ispecies))
    Jac_dn(ispecies,ispecies) = (unit_n_dn*area) * &
                                (u - harmonic_D_over_dist(ispecies))
  enddo

  deallocate(diffusivity_dn)
  deallocate(diffusivity_up)
  deallocate(molecular_diffusion_dn)
  deallocate(molecular_diffusion_up)
                    
end subroutine NWTJacobianFlux

! ************************************************************************** !

subroutine NWTComputeMassBalance(realization,max_size,sum_mol)
  !
  ! Sums up the amount of moles in each component.
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/27/2019
  ! 

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module


  type(realization_subsurface_type) :: realization
  PetscInt :: max_size
  PetscReal :: sum_mol(max_size,4)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(nw_trans_realization_type), pointer :: nw_trans

  PetscReal :: sum_mol_tot(max_size)
  PetscReal :: sum_mol_aq(max_size)
  PetscReal :: sum_mol_sb(max_size)
  PetscReal :: sum_mol_mnrl(max_size)

  PetscErrorCode :: ierr
  PetscInt :: local_id, ghosted_id
  PetscInt :: ncomp
  PetscReal :: liquid_saturation, porosity, volume

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  nw_trans => realization%nw_trans

  nwt_auxvars => patch%aux%NWT%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars

  sum_mol = 0.d0
  sum_mol_tot = 0.d0
  sum_mol_aq = 0.d0
  sum_mol_sb = 0.d0
  sum_mol_mnrl = 0.d0

  ncomp = nw_trans%params%ncomp

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    ! ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    
    liquid_saturation = global_auxvars(ghosted_id)%sat(LIQUID_PHASE)
    porosity = material_auxvars(ghosted_id)%porosity
    volume = material_auxvars(ghosted_id)%volume ! [m^3]
    
    ! aqueous (sum_mol_aq) [mol]
    sum_mol_aq(1:ncomp) = sum_mol_aq(1:ncomp) + &
           nwt_auxvars(ghosted_id)%aqueous_eq_conc(:) * &  ! [mol/m^3-liq]
           liquid_saturation*porosity*volume               ! [m^3-liq]

    ! equilibrium sorption (sum_mol_sb) [mol]
    sum_mol_sb(1:ncomp) = sum_mol_sb(1:ncomp) + &
            nwt_auxvars(ghosted_id)%sorb_eq_conc(:) * &    ! [mol/m^3-sorb]
            volume                                         ! [m^3-sorb]
    ! jenn:tod Is this the correct calc for sum_mol_sb? Is volume right?

    ! mineral volume fractions (sum_mol_mnrl) [mol]
    sum_mol_mnrl(1:ncomp) = sum_mol_mnrl(1:ncomp) + &
            nwt_auxvars(ghosted_id)%mnrl_eq_conc(:) * &    ! [mol/m^3-mnrl]
            nwt_auxvars(ghosted_id)%mnrl_vol_frac(:) * &   ! [m^3-mnrl/m^3-void]
            porosity*volume                     ! [m^3-void/m^3-bulk * m^3-bulk]
  enddo

  sum_mol_tot = sum_mol_aq + sum_mol_sb + sum_mol_mnrl     ! [mol]

  sum_mol(:,1) = sum_mol_tot
  sum_mol(:,2) = sum_mol_aq
  sum_mol(:,3) = sum_mol_sb
  sum_mol(:,4) = sum_mol_mnrl

end subroutine NWTComputeMassBalance

! ************************************************************************** !

subroutine NWTUpdateMassBalance(realization)
  ! 
  ! Updates the mass balance.
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/27/2019
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars_bc(:)
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  nwt_auxvars_bc => patch%aux%NWT%auxvars_bc
  nwt_auxvars_ss => patch%aux%NWT%auxvars_ss

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%NWT%num_aux
    patch%aux%NWT%auxvars(iconn)%mass_balance = &
      patch%aux%NWT%auxvars(iconn)%mass_balance + &
      patch%aux%NWT%auxvars(iconn)%mass_balance_delta*option%tran_dt
  enddo
#endif

  do iconn = 1, patch%aux%NWT%num_aux_bc
    nwt_auxvars_bc(iconn)%mass_balance = &
      nwt_auxvars_bc(iconn)%mass_balance + &
      nwt_auxvars_bc(iconn)%mass_balance_delta*option%tran_dt
  enddo

  do iconn = 1, patch%aux%NWT%num_aux_ss
    nwt_auxvars_ss(iconn)%mass_balance = &
      nwt_auxvars_ss(iconn)%mass_balance + &
      nwt_auxvars_ss(iconn)%mass_balance_delta*option%tran_dt
  enddo

end subroutine NWTUpdateMassBalance

! ************************************************************************** !

subroutine NWTDestroy(realization)
  !
  ! Destroys objects in the NW Transport module.
  !
  ! Author: Jenn Frederick
  ! Date: 05/27/2019
  !
  
  use Realization_Subsurface_class
  
  implicit none
  
  type(realization_subsurface_type) :: realization
  
  ! placeholder, does nothing at the moment
  ! note: the aux objects are destroyed from auxiliary.F90 
  !       and the realization nw_trans object is destroyed in 
  !       realization_subsurface.F90, RealizationStrip().

end subroutine NWTDestroy

! ************************************************************************** !

end module NW_Transport_module
