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
            NWTZeroMassBalanceDelta, &
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
  use Option_module
  use Grid_module
  use Material_module
  use Material_Aux_class
  use Coupler_module
  use Condition_module
  use Connection_module
  use Fluid_module
  use Output_Aux_module
  
  implicit none

  type(realization_subsurface_type) :: realization
  
  class(reaction_nw_type), pointer :: reaction_nw
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
  
  grid => realization%patch%grid
  reaction_nw => realization%reaction_nw
  option => realization%option
  nspecies = reaction_nw%params%nspecies
  nphase = option%transport%nphase
  
  realization%patch%aux%NWT => NWTAuxCreate()
  
  cur_material_property => realization%material_properties
  do                                      
    if (.not.associated(cur_material_property)) exit
    if (maxval(cur_material_property%dispersivity(2:3)) > 0.d0) then
      reaction_nw%params%calculate_transverse_dispersion = PETSC_TRUE
      exit
    endif
    cur_material_property => cur_material_property%next
  enddo
  
  material_auxvars => realization%patch%aux%Material%auxvars
  flag = 0
  !TODO(geh): change to looping over ghosted ids once the legacy code is 
  !           history and the communicator can be passed down.
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)

    ! Ignore inactive cells with inactive materials
    if (realization%patch%imat(ghosted_id) <= 0) cycle
    
    if (material_auxvars(ghosted_id)%volume < 0.d0 .and. flag(1) == 0) then
      flag(1) = 1
      option%io_buffer = 'Non-initialized cell volume.'
      call PrintMsg(option)
    endif
    if (material_auxvars(ghosted_id)%porosity < 0.d0 .and. flag(2) == 0) then
      flag(2) = 1
      option%io_buffer = 'Non-initialized porosity.'
      call PrintMsg(option)
    endif
    if (material_auxvars(ghosted_id)%tortuosity < 0.d0 .and. flag(3) == 0) then
      flag(3) = 1
      option%io_buffer = 'Non-initialized tortuosity.'
      call PrintMsg(option)
    endif 
  
  enddo 
  
  if (maxval(flag) > 0) then
    option%io_buffer = &
      'Material property errors found in NWTSetup (Nuclear Waste Transport).'
    call PrintErrMsg(option)
  endif
    
  ! allocate auxvar data structures for all grid cells
#ifdef COMPUTE_INTERNAL_MASS_FLUX
  option%iflag = 1 ! allocate mass_balance array
#else  
  option%iflag = 0 ! be sure not to allocate mass_balance array
#endif
  allocate(realization%patch%aux%NWT%auxvars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call NWTAuxVarInit(realization%patch%aux%NWT%auxvars(ghosted_id),reaction_nw, &
                      option)
  enddo
  realization%patch%aux%NWT%num_aux = grid%ngmax
  
  ! count the number of boundary connections and allocate
  ! auxvar data structures for them
  sum_connection = CouplerGetNumConnectionsInList( &
                                     realization%patch%boundary_condition_list)
  if (sum_connection > 0) then
    option%iflag = 1 ! enable allocation of mass_balance array 
    allocate(realization%patch%aux%NWT%auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call NWTAuxVarInit(realization%patch%aux%NWT%auxvars_bc(iconn),reaction_nw, &
                         option)
    enddo
  endif
  realization%patch%aux%NWT%num_aux_bc = sum_connection
  
  ! count the number of source/sink connections and allocate
  ! auxvar data structures for them
  sum_connection = CouplerGetNumConnectionsInList( &
                                            realization%patch%source_sink_list)
  if (sum_connection > 0) then
    option%iflag = 1 ! enable allocation of mass_balance array 
    allocate(realization%patch%aux%NWT%auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call NWTAuxVarInit(realization%patch%aux%NWT%auxvars_ss(iconn),reaction_nw, &
                         option)
    enddo
  endif
  realization%patch%aux%NWT%num_aux_ss = sum_connection
  option%iflag = 0
  
  ! initialize parameters
  allocate(reaction_nw%diffusion_coefficient(nspecies,nphase))
  reaction_nw%diffusion_coefficient = 1.d-9
  allocate(reaction_nw%diffusion_activation_energy(nspecies,nphase))
  reaction_nw%diffusion_activation_energy = 0.d0
  allocate(reaction_nw%species_print(nspecies))
  reaction_nw%species_print = PETSC_FALSE
  
  cur_fluid_property => realization%fluid_properties
  do 
    if (.not.associated(cur_fluid_property)) exit
    iphase = cur_fluid_property%phase_id
    ! setting of phase diffusion coefficients must come before individual
    ! species below
    if (iphase <= nphase) then
      reaction_nw%diffusion_coefficient(:,iphase) = &
        cur_fluid_property%diffusion_coefficient
      reaction_nw%diffusion_activation_energy(:,iphase) = &
        cur_fluid_property%diffusion_activation_energy
    endif
    cur_fluid_property => cur_fluid_property%next
  enddo
  
  !TODO(jenn) Will we support species-dependent diffusion coefficients?
  ! If so, look at reactive_transport.F90, beginning line 338 in RTSetup().
  
  ! setup output
  list => realization%output_option%output_snap_variable_list
  call NWTSetPlotVariables(list,reaction_nw,option, &
                           realization%output_option%tunit)
  if (.not.associated(realization%output_option%output_snap_variable_list, &
                      realization%output_option%output_obs_variable_list)) then
    list => realization%output_option%output_obs_variable_list
    call NWTSetPlotVariables(list,reaction_nw,option, &
                             realization%output_option%tunit)
  endif
  
end subroutine NWTSetup

! ************************************************************************** !

subroutine NWTProcessConstraint(reaction_nw,constraint_name, &
                                nwt_species_constraint,option)
  ! 
  ! Ensures ordering of species is consistant between the reaction_nw object
  ! and the constraint object. 
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/22/2019
  ! 
  use Option_module
  use String_module
  use Utility_module
  
  implicit none
  
  class(reaction_nw_type), pointer :: reaction_nw
  character(len=MAXWORDLENGTH) :: constraint_name
  type(nwt_species_constraint_type), pointer :: nwt_species_constraint
  type(option_type) :: option
  
  PetscBool :: found
  PetscInt :: ispecies, jspecies
  PetscReal :: constraint_conc(reaction_nw%params%nspecies)
  PetscInt :: constraint_type(reaction_nw%params%nspecies)
  character(len=MAXWORDLENGTH) :: constraint_species_names( &
                                                     reaction_nw%params%nspecies)
  
  constraint_conc = 0.d0
  constraint_type = 0
  constraint_species_names = ''
  
  do ispecies = 1, reaction_nw%params%nspecies
    found = PETSC_FALSE
    do jspecies = 1, reaction_nw%params%nspecies
      if (StringCompare(nwt_species_constraint%names(ispecies), &
                        reaction_nw%species_names(jspecies),MAXWORDLENGTH)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      option%io_buffer = &
               'Species ' // trim(nwt_species_constraint%names(ispecies)) // &
               ' from CONSTRAINT ' // trim(constraint_name) // &
               ' not found among species.'
      call PrintErrMsg(option)
    else
      constraint_conc(jspecies) = &
                               nwt_species_constraint%constraint_conc(ispecies)
      constraint_type(jspecies) = &
                               nwt_species_constraint%constraint_type(ispecies)
      constraint_species_names(jspecies) = &
                                         nwt_species_constraint%names(ispecies)
    endif
  enddo
  
  ! place ordered constraint parameters back in original arrays
  nwt_species_constraint%constraint_conc = constraint_conc
  nwt_species_constraint%constraint_type = constraint_type
  nwt_species_constraint%names = constraint_species_names
    
end subroutine NWTProcessConstraint

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
  use Grid_module
  use Coupler_module
  use Connection_module
  use Option_module
  use Field_module
  use Logging_module
  use Global_Aux_module
  use Transport_Constraint_NWT_module
  
  implicit none

  type(realization_subsurface_type) :: realization
  PetscBool :: update_bcs
  PetscBool :: update_cells
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  class(reaction_nw_type), pointer :: reaction_nw
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: ghosted_id, local_id, sum_connection, iconn
  PetscInt :: istart, iend
  PetscReal, pointer :: xx_loc_p(:)
  PetscInt, parameter :: iphase = 1
  PetscInt :: offset
  PetscErrorCode :: ierr
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)  
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars(:)
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars_bc(:)
  type(nw_transport_auxvar_type), pointer :: nwt_auxvar
  class(tran_constraint_coupler_nwt_type), pointer :: constraint_coupler
  PetscInt, save :: icall
  
  data icall/0/

  option => realization%option
  grid => realization%patch%grid
  field => realization%field
  reaction_nw => realization%reaction_nw
  material_auxvars => realization%patch%aux%Material%auxvars
  nwt_auxvars => realization%patch%aux%NWT%auxvars
  nwt_auxvars_bc => realization%patch%aux%NWT%auxvars_bc
  global_auxvars => realization%patch%aux%Global%auxvars
  global_auxvars_bc => realization%patch%aux%Global%auxvars_bc

  
  call VecGetArrayReadF90(field%tran_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

  if (update_cells) then

    call PetscLogEventBegin(logging%event_nwt_auxvars,ierr);CHKERRQ(ierr)
  
    do ghosted_id = 1, grid%ngmax
      if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
      !geh - Ignore inactive cells with inactive materials
      if (realization%patch%imat(ghosted_id) <= 0) cycle

      offset = (ghosted_id-1)*reaction_nw%params%nspecies
      istart = offset + 1
      iend = offset + reaction_nw%params%nspecies
      
      nwt_auxvars(ghosted_id)%total_bulk_conc = xx_loc_p(istart:iend)

      call NWTAuxVarCompute(nwt_auxvars(ghosted_id), &
                            global_auxvars(ghosted_id), &
                            material_auxvars(ghosted_id), &
                            reaction_nw,option)

    enddo

    call PetscLogEventEnd(logging%event_rt_auxvars,ierr);CHKERRQ(ierr)
  endif

  if (update_bcs) then

    call PetscLogEventBegin(logging%event_rt_auxvars_bc,ierr);CHKERRQ(ierr)

    boundary_condition => realization%patch%boundary_condition_list%first
    sum_connection = 0    
    do 
      if (.not.associated(boundary_condition)) exit
      cur_connection_set => boundary_condition%connection_set
      nwt_auxvar => &
        TranConstraintNWTGetAuxVar(boundary_condition%tran_condition% &
                                   cur_constraint_coupler)
      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1 
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        
!geh: Since a minimum precipitate concentration of 1.d-20 is always present
!     (see NWTEqDissPrecipSorb), we must NWTAuxVarCompute() based on the
!     total bulk concentration at each boundary connection just like we do
!     for each internal grid cell above. Otherwise, the precipitate
!     concentation of 1.d-20 is not factored into the boundary concentration
!     for pure aqueous boundaries, and this generates error. To prove this
!     change #if 0 -> #if 1 below and run a transport simulation with a
!     single aqueous constraint (not concentration gradient). You will see
!     that the Newton solve struggles because slightly different aqueous
!     concentrations are assigned to the boundary faces than the internal
!     cell centers.

!     ***I propose that we do away with the minimum sorbed and precipitate
!        concentrations altogether as they add artificial mass to the system.***

#if 1
        nwt_auxvars_bc(sum_connection)%total_bulk_conc(:) = &
                       nwt_auxvar%total_bulk_conc(:)
        call NWTAuxVarCompute(nwt_auxvars_bc(sum_connection), &
                              global_auxvars_bc(sum_connection), &
                              material_auxvars(ghosted_id), &
                              reaction_nw,option)
#else
        nwt_auxvars_bc(sum_connection)%total_bulk_conc(:) = &
                       nwt_auxvar%total_bulk_conc(:)
        nwt_auxvars_bc(sum_connection)%aqueous_eq_conc(:) = &
                       nwt_auxvar%aqueous_eq_conc(:)
        nwt_auxvars_bc(sum_connection)%sorb_eq_conc(:) = &
                       nwt_auxvar%sorb_eq_conc(:)
        nwt_auxvars_bc(sum_connection)%mnrl_eq_conc(:) = &
                       nwt_auxvar%mnrl_eq_conc(:)
        nwt_auxvars_bc(sum_connection)%mnrl_vol_frac(:) = &
                       nwt_auxvar%mnrl_vol_frac(:)
#endif
          
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
                            reaction_nw,option)
  ! 
  ! Computes the secondary variables from the primary dependent variable.
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/29/2019
  ! 

  use Option_module
  use Material_Aux_class
  use NWT_Equilibrium_module
  
  implicit none
  
  type(nw_transport_auxvar_type) :: nwt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  class(reaction_nw_type) :: reaction_nw
  type(option_type) :: option
  
  type(species_type), pointer :: cur_species
  PetscReal :: solubility(reaction_nw%params%nspecies)  ! [mol/m^3-liq]
  PetscReal :: aq_mass(reaction_nw%params%nspecies)     ! [mol/m^3-liq]
  PetscReal :: ppt_mass(reaction_nw%params%nspecies)    ! [mol/m^3-bulk]
  PetscReal :: sorb_mass(reaction_nw%params%nspecies)   ! [mol/m^3-bulk]
  PetscReal :: mnrl_molar_density(reaction_nw%params%nspecies)  ! [mol/m^3-mnrl]
  PetscReal :: ele_kd(reaction_nw%params%nspecies)      ! [m^3-water/m^3-bulk]
  PetscBool :: dry_out
  PetscInt :: ispecies
  PetscReal :: sat, por

  sat = max(MIN_LIQ_SAT,global_auxvar%sat(LIQUID_PHASE))
  por = material_auxvar%porosity
  
  if (sat > 0.d0) then
    dry_out = PETSC_FALSE
  else
    dry_out = PETSC_TRUE
  endif
  
  cur_species => reaction_nw%species_list
  do 
    if (.not.associated(cur_species)) exit
    solubility(cur_species%id) = cur_species%solubility_limit
    mnrl_molar_density(cur_species%id) = cur_species%mnrl_molar_density
    ele_kd(cur_species%id) = cur_species%ele_kd
    cur_species => cur_species%next
  enddo
  
  ! calculate the equilibrium state and partition the total_bulk_conc
  do ispecies = 1,reaction_nw%params%nspecies
    call NWTEqDissPrecipSorb(solubility(ispecies),material_auxvar, &
                             global_auxvar,dry_out,ele_kd(ispecies), &
                             nwt_auxvar%total_bulk_conc(ispecies), &
                             aq_mass(ispecies),ppt_mass(ispecies), &
                             sorb_mass(ispecies))
  enddo

  !-------aqueous concentration (equilibrium)
  nwt_auxvar%aqueous_eq_conc(:) = aq_mass(:)
                     
  !-------sorbed concentration (equilibrium)
  nwt_auxvar%sorb_eq_conc(:) = sorb_mass(:)
  
  !-------precipitated concentration (equilibrium)
  nwt_auxvar%mnrl_eq_conc(:) = ppt_mass(:)
  nwt_auxvar%mnrl_vol_frac(:) = nwt_auxvar%mnrl_eq_conc(:)/ &
                              (material_auxvar%porosity*mnrl_molar_density(:))

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
  use Coupler_module
  use Connection_module
  use Discretization_module
  use Option_module
  use Grid_module
  use Logging_module
  use Debug_module
  use WIPP_Flow_Aux_module

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
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  class(reaction_nw_type), pointer :: reaction_nw
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
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
  PetscBool :: bc
  PetscReal :: Res(realization%reaction_nw%params%nspecies)
  PetscReal :: Res_up(realization%reaction_nw%params%nspecies)
  PetscReal :: Res_dn(realization%reaction_nw%params%nspecies)
  PetscReal :: area
  PetscViewer :: viewer  
  
  character(len=MAXSTRINGLENGTH) :: string

  call PetscLogEventBegin(logging%event_nwt_residual,ierr);CHKERRQ(ierr)

  field => realization%field
  discretization => realization%discretization
  option => realization%option
  grid => realization%patch%grid
  reaction_nw => realization%patch%reaction_nw
  nwt_auxvars => realization%patch%aux%NWT%auxvars
  nwt_auxvars_ss => realization%patch%aux%NWT%auxvars_ss
  nwt_auxvars_bc => realization%patch%aux%NWT%auxvars_bc
  global_auxvars => realization%patch%aux%Global%auxvars
  global_auxvars_bc => realization%patch%aux%Global%auxvars_bc
  material_auxvars => realization%patch%aux%Material%auxvars
  ! note: there is no realization%patch%aux%Material%auxvars_bc
  material_auxvars_bc => realization%patch%aux%Material%auxvars
  nphase = reaction_nw%params%nphase
  nspecies = reaction_nw%params%nspecies
  
  ! Strictly used to check if alpha needs to be considered in the area
  nullify(wippflo_auxvars)
  if (associated(realization%patch%aux%WIPPFlo)) then
    wippflo_auxvars => realization%patch%aux%WIPPFlo%auxvars
  endif

  ! Communication -----------------------------------------
  if (realization%reaction_nw%use_log_formulation) then
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
  
  ! Zero out the residual pointer
  r_p = 0.d0
  
  ! Update the auxiliary variables
  call NWTUpdateAuxVars(realization,PETSC_TRUE,PETSC_TRUE)
    
#if 1
  !== Accumulation Terms ======================================
  if (.not.option%steady_state) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      ! ignore inactive cells with inactive materials
      if (realization%patch%imat(ghosted_id) <= 0) cycle

      !WRITE(*,*)  '  auxvars_aqc = ', nwt_auxvars(ghosted_id)%aqueous_eq_conc(:) 
       
      call NWTResidualAccum(nwt_auxvars(ghosted_id), &
                            global_auxvars(ghosted_id), &
                            material_auxvars(ghosted_id), &
                            reaction_nw,Res)
                            
      offset = (local_id-1)*nspecies
      istart = offset + 1
      iend = offset + nspecies
      r_p(istart:iend) = r_p(istart:iend) - &
        (Res(1:nspecies) - fixed_accum_p(istart:iend))/option%tran_dt
      
    enddo
  endif
#endif

  !WRITE(*,*)  '     ResAccum = ', Res(:)
  !WRITE(*,*)  '       r_p(2) = ', r_p(723), '     loc = 723'
  !WRITE(*,*)  '       r_p(2) = ', r_p(791), '     loc = 791'
  
#if 1
  !== Source/Sink Terms =======================================
  source_sink => realization%patch%source_sink_list%first
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    do iconn = 1, cur_connection_set%num_connections 
      sum_connection = sum_connection + 1     
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      ! ignore inactive cells with inactive materials
      if (realization%patch%imat(ghosted_id) <= 0) cycle
            
      call NWTResidualSrcSink(nwt_auxvars(ghosted_id), &
                              global_auxvars(ghosted_id), &
                              source_sink,realization%patch%ss_flow_vol_fluxes, &
                              sum_connection,reaction_nw,Res)
      !WRITE(*,*)  '  SrcSinkName = ', source_sink%name                       
      !WRITE(*,*)  '   ResSrcSink = ', Res(:)
      
      offset = (local_id-1)*nspecies
      istart = offset + 1
      iend = offset + nspecies
      r_p(istart:iend) = r_p(istart:iend) + Res(1:nspecies)
      
      if (associated(realization%patch%ss_tran_fluxes)) then
        realization%patch%ss_tran_fluxes(:,sum_connection) = - Res(:)
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

  !WRITE(*,*)  '       r_p(3) = ', r_p(:)

#if 1
  !== Decay and Ingrowth ======================================
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    ! ignore inactive cells with inactive materials
    if (realization%patch%imat(ghosted_id) <= 0) cycle
        
    call NWTResidualRx(nwt_auxvars(ghosted_id), &
                       material_auxvars(ghosted_id), &
                       reaction_nw,Res)
    
    offset = (local_id-1)*nspecies
    istart = offset + 1
    iend = offset + nspecies
    r_p(istart:iend) = r_p(istart:iend) + Res(1:nspecies)  
    
  enddo
#endif

  !WRITE(*,*)  '        ResRx = ', Res(:)
  !WRITE(*,*)  '       r_p(4) = ', r_p(:)

#if 1
  !== Fluxes ==================================================
  if (option%compute_mass_balance_new) then
    call NWTZeroMassBalanceDelta(realization)
  endif
  
  ! Interior Flux Terms ---------------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0 
  bc = PETSC_FALSE
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      
      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! ghost to local mapping
      
      area = cur_connection_set%area(iconn)
      if (associated(wippflo_auxvars)) then
        area = area * 0.5d0 * &
               (wippflo_auxvars(ZERO_INTEGER,ghosted_id_up)%alpha + &
                wippflo_auxvars(ZERO_INTEGER,ghosted_id_dn)%alpha)
      endif
            
      ! ignore inactive cells with inactive materials
      if (realization%patch%imat(ghosted_id_up) <= 0 .or.  &
          realization%patch%imat(ghosted_id_dn) <= 0) cycle
          
      call NWTResidualFlux(nwt_auxvars(ghosted_id_up), &
                      nwt_auxvars(ghosted_id_dn), &
                      global_auxvars(ghosted_id_up), &
                      global_auxvars(ghosted_id_dn), &
                      material_auxvars(ghosted_id_up), &
                      material_auxvars(ghosted_id_dn), &
                      area, &
                      cur_connection_set%dist(:,iconn), &
                      realization%patch%internal_velocities(:,sum_connection), &
                      reaction_nw,option,bc,Res_up,Res_dn)
                            
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
      
      if (associated(realization%patch%internal_tran_fluxes)) then
      ! NOTE: (jenn)
      ! Res_dn and Res_up are equal and opposite. I have chosen to consistently
      ! assign Res_dn for the internal transport flux value for the integral
      ! flux calculation.
        realization%patch%internal_tran_fluxes(1:nspecies,iconn) = &
            Res_dn(1:nspecies)
      endif
      
    enddo
    cur_connection_set => cur_connection_set%next
  enddo
   
  !WRITE(*,*)  '       r_p(5) = ', r_p(723), '     loc = 723'
  !WRITE(*,*)  '       r_p(5) = ', r_p(791), '     loc = 791'
  !WRITE(*,*)  '     int_flux = ', r_p(723)*option%tran_dt, '     loc = 723'
  !WRITE(*,*)  '     int_flux = ', r_p(791)*option%tran_dt, '     loc = 791'

  !WRITE(*,*)  '   g_sat(792) = ', (1.0 - global_auxvars(792)%sat(LIQUID_PHASE))
  !WRITE(*,*)  ' aq_conc(792) = ', nwt_auxvars(792)%aqueous_eq_conc(:)
  
  ! Boundary Flux Terms ---------------------------------------
  boundary_condition => realization%patch%boundary_condition_list%first
  sum_connection = 0    
  bc = PETSC_TRUE
  do
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
    
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      
      area = cur_connection_set%area(iconn)

      ! ignore inactive cells with inactive materials
      if (realization%patch%imat(ghosted_id) <= 0) cycle
      
      call NWTResidualFlux(nwt_auxvars_bc(sum_connection), &
                      nwt_auxvars(ghosted_id), &
                      global_auxvars_bc(sum_connection), &
                      global_auxvars(ghosted_id), &
                      material_auxvars_bc(sum_connection), &
                      material_auxvars(ghosted_id), &
                      area, &
                      cur_connection_set%dist(:,iconn), &
                      realization%patch%boundary_velocities(:,sum_connection), &
                      reaction_nw,option,bc,Res_up,Res_dn)
                            
      offset = (local_id-1)*nspecies
      istart = offset + 1
      iend = offset + nspecies
      r_p(istart:iend) = r_p(istart:iend) + Res_dn(1:nspecies) 
      ! note: Don't need to worry about Res_up because that is outside of
      ! the domain, and doesn't have a place in r_p.
      
      if (option%compute_mass_balance_new) then
      ! contribution to boundary
        iphase = LIQUID_PHASE
        nwt_auxvars_bc(sum_connection)%mass_balance_delta(:,iphase) = &
          nwt_auxvars_bc(sum_connection)%mass_balance_delta(:,iphase) +  &
          Res_dn(1:nspecies)
      endif  
                 
      if (associated(realization%patch%boundary_tran_fluxes)) then
      ! NOTE: (jenn)
      ! Res_dn and Res_up are equal and opposite. I have chosen to consistently
      ! assign Res_dn for the internal transport flux value for the integral
      ! flux calculation.
        realization%patch%boundary_tran_fluxes(1:nspecies,sum_connection) = &
            Res_dn(1:nspecies)
      endif
    enddo
    boundary_condition => boundary_condition%next
  enddo
#endif


  !WRITE(*,*)  '       r_p(6) = ', r_p(:)
  
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
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  class(reaction_nw_type), pointer :: reaction_nw
  PetscReal, pointer :: xx_p(:), fixed_accum_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: dof_offset, istart, iend
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  nwt_auxvars => realization%patch%aux%NWT%auxvars
  global_auxvars => realization%patch%aux%Global%auxvars
  material_auxvars => realization%patch%aux%Material%auxvars
  grid => realization%patch%grid
  reaction_nw => realization%reaction_nw

  ! cannot use tran_xx_loc vector here as it has not yet been updated.
  call VecGetArrayReadF90(field%tran_xx,xx_p, ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%tran_accum, fixed_accum_p, ierr);CHKERRQ(ierr)
    
! Do not use NWTUpdateAuxVars() as it loops over ghosted ids

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (realization%patch%imat(ghosted_id) <= 0) cycle
    
    ! compute offset in solution vector for first dof in grid cell
    dof_offset = (local_id-1)*reaction_nw%params%nspecies
    
    ! calculate range of species
    istart = dof_offset + 1
    iend = dof_offset + reaction_nw%params%nspecies

    ! copy primary dependent variable
    nwt_auxvars(ghosted_id)%total_bulk_conc = xx_p(istart:iend)
    
    call NWTAuxVarCompute(nwt_auxvars(ghosted_id), &
                          global_auxvars(ghosted_id), &
                          material_auxvars(ghosted_id), &
                          reaction_nw,option)
    call NWTResidualAccum(nwt_auxvars(ghosted_id), &
                          global_auxvars(ghosted_id), &
                          material_auxvars(ghosted_id), &
                          reaction_nw,fixed_accum_p(istart:iend)) 
  enddo

  call VecRestoreArrayReadF90(field%tran_xx,xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%tran_accum,fixed_accum_p,ierr);CHKERRQ(ierr)

end subroutine NWTUpdateFixedAccumulation

! ************************************************************************** !

subroutine NWTResidualAccum(nwt_auxvar,global_auxvar,material_auxvar, &
                            reaction_nw,Res)
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
  class(reaction_nw_type), pointer :: reaction_nw
  PetscReal :: Res(reaction_nw%params%nspecies)
  
  PetscInt :: istart, iend
  PetscReal :: por, sat, volume
  
  Res = 0.d0
  
  ! All residual entries for accumulation should be in [mol-species].
  ! Dividing by dt will occur later in NWTResidual, so that the units match
  ! the rest of the residual units of [mol/sec].
  
  ! porosity in [m^3-void/m^3-bulk]
  ! saturation in [m^3-liq/m^3-void]
  ! volume in [m^3-bulk]
  por = material_auxvar%porosity
  sat = max(MIN_LIQ_SAT,global_auxvar%sat(LIQUID_PHASE))
  volume = material_auxvar%volume
      
  istart = 1
  iend = reaction_nw%params%nspecies
  
  ! -- Aqueous-Component ----------------------------------------
  
  ! aqueous conc in [mol-species/m^3-liq]
  Res(istart:iend) = volume*por*sat*nwt_auxvar%aqueous_eq_conc(:)
                     
  ! -- Precipitated-Component -----------------------------------
  ! precipitated conc in [mol-species/m^3-bulk]
  Res(istart:iend) = Res(istart:iend) + &
                     volume*nwt_auxvar%mnrl_eq_conc(:)
                  
  ! -- Sorbed-Component -----------------------------------------
  ! sorbed conc in [mol-species/m^3-bulk]
  Res(istart:iend) = Res(istart:iend) + &
                     volume*nwt_auxvar%sorb_eq_conc(:)

end subroutine NWTResidualAccum

! ************************************************************************** !

subroutine NWTResidualSrcSink(nwt_auxvar,global_auxvar,source_sink,&
                              ss_flow_vol_fluxes,sum_connection,reaction_nw,Res)
  ! 
  ! Computes the source/sink terms in the residual function.
  ! All residual entries should be in [mol-species/sec].
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/08/2019
  ! 
  use Coupler_module
  use Transport_Constraint_NWT_module

  implicit none
  
  type(nw_transport_auxvar_type) :: nwt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(coupler_type), pointer :: source_sink
  PetscReal, pointer :: ss_flow_vol_fluxes(:,:)
  PetscInt :: sum_connection
  class(reaction_nw_type), pointer :: reaction_nw
  type(nw_transport_auxvar_type), pointer :: nwt_auxvar_out
  PetscReal :: Res(reaction_nw%params%nspecies)
  
  PetscInt :: istart, iend
  PetscReal :: qsrc, sat
  PetscReal :: coef_in, coef_out
  
  sat = max(MIN_LIQ_SAT,global_auxvar%sat(LIQUID_PHASE))
  Res = 0.d0

  nwt_auxvar_out => &
    TranConstraintNWTGetAuxVar(source_sink%tran_condition% &
                                   cur_constraint_coupler)
  
  if (associated(ss_flow_vol_fluxes)) then
    ! qsrc = [m^3-liq/sec] 
    qsrc = ss_flow_vol_fluxes(LIQUID_PHASE,sum_connection)
  endif
      
  istart = 1
  iend = reaction_nw%params%nspecies
  
  ! -- Aqueous-Component ----------------------------------------
  select case(source_sink%tran_condition%itype)
    case(EQUILIBRIUM_SS)
      !TODO(jenn) What is EQUILIBRIUM_SS option?
    case(MASS_RATE_SS)
      !TODO(jenn) What is MASS_RATE_SS option?
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
  Res(istart:iend) = (coef_in*sat*nwt_auxvar%aqueous_eq_conc(:)) + &
                     (coef_out*sat*nwt_auxvar_out%aqueous_eq_conc(:))
  !WRITE(*,*)  '      coef_in = ', coef_in
  !WRITE(*,*)  '   in aq_conc = ', nwt_auxvar%aqueous_eq_conc(:)
  !WRITE(*,*)  '     coef_out = ', coef_out
  !WRITE(*,*)  '  out aq_conc = ', source_sink%tran_condition% &
  !                                cur_constraint_coupler%nwt_auxvar% &
  !                                aqueous_eq_conc(:)
                               
  ! -- Precipitated-Component -----------------------------------
  ! There is no contribution from precipitated components.
  
  ! -- Sorbed-Component -----------------------------------------
  ! There is no contribution from sorbed components.

end subroutine NWTResidualSrcSink

! ************************************************************************** !

subroutine NWTResidualRx(nwt_auxvar,material_auxvar,reaction_nw,Res)
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
  class(reaction_nw_type), pointer :: reaction_nw
  PetscReal :: Res(reaction_nw%params%nspecies)
  
  PetscInt :: iphase, parent_id
  PetscReal :: vol
  type(species_type), pointer :: species
  type(radioactive_decay_rxn_type), pointer :: rad_rxn
  
  iphase = LIQUID_PHASE
  Res = 0.d0
  
  ! All residual entries for decay/ingrowth should be in [mol-species/sec].
  
  ! volume in [m^3-bulk]
  vol = material_auxvar%volume
  
  !TODO(jenn): convert to compressed arrays instead of tranversing a linked
  !            list, which is much less efficient for memory access. see
  !            'do icplx' loop in RTotalAqueous for an example.
  species => reaction_nw%species_list
  do 
    if (.not.associated(species)) exit
    
    if (species%radioactive) then
      ! Find the reaction object associated with this species
      rad_rxn => reaction_nw%rad_decay_rxn_list
      do
        if (.not.associated(rad_rxn)) exit
        if (rad_rxn%species_id == species%id) exit
        rad_rxn => rad_rxn%next
      enddo
      ! Add in species decay
      Res(species%id) = -(rad_rxn%rate_constant * &
                          nwt_auxvar%total_bulk_conc(species%id))
      ! units are [mol-species/m^3-bulk/sec] right now
      
      ! Add in contribution from parent (if exists)
      if (rad_rxn%parent_id > 0.d0) then
        parent_id = rad_rxn%parent_id
        ! Find the reaction object associated with the parent species
        rad_rxn => reaction_nw%rad_decay_rxn_list
        do
          if (.not.associated(rad_rxn)) exit
          if (rad_rxn%species_id == parent_id) exit
          rad_rxn => rad_rxn%next
        enddo
        Res(species%id) = Res(species%id) + (rad_rxn%rate_constant * &
                          nwt_auxvar%total_bulk_conc(parent_id))
        ! units are [mol-species/m^3-bulk/sec] right now
      endif
    endif
    
    species => species%next
  enddo
  
  Res = Res*vol
  ! units are now [mol-species/sec]
                    
end subroutine NWTResidualRx

! ************************************************************************** !

subroutine NWTResidualFlux(nwt_auxvar_up,nwt_auxvar_dn, &
                           global_auxvar_up,global_auxvar_dn, &
                           material_auxvar_up,material_auxvar_dn, &
                           area,dist,velocity,reaction_nw,option,bc,Res_up,Res_dn)
  ! 
  ! Computes the flux terms in the residual function.
  ! All residual entries should be in [mol-species/sec].
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/09/2019
  ! 

  use Option_module
  use Connection_module

  implicit none
  
  ! if bc, _up is actually _bc:
  type(nw_transport_auxvar_type) :: nwt_auxvar_up, nwt_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  PetscReal :: area, dist(-1:3)
  PetscReal :: velocity(*)
  class(reaction_nw_type), pointer :: reaction_nw
  type(option_type) :: option
  PetscBool :: bc
  PetscReal :: Res_up(reaction_nw%params%nspecies)
  PetscReal :: Res_dn(reaction_nw%params%nspecies)
  
  PetscInt :: unit_n_up, unit_n_dn
  PetscInt :: nspecies
  PetscReal :: q, u
  PetscReal :: harmonic_ps
  PetscReal :: sat_up, sat_dn
  PetscReal :: por_up, por_dn
  PetscReal :: ps_up, ps_dn
  PetscReal :: dist_up, dist_dn
  PetscReal, pointer :: diffusivity_up(:)
  PetscReal, pointer :: diffusivity_dn(:)
  PetscReal, pointer :: molecular_diffusion_up(:)
  PetscReal, pointer :: molecular_diffusion_dn(:)
  PetscReal, pointer :: diffusion_coefficient(:,:)
  PetscReal :: distance_gravity, upwind_weight ! both are dummy variables
  PetscReal :: harmonic_D_over_dist(reaction_nw%params%nspecies)
  PetscReal :: diffusive_flux(reaction_nw%params%nspecies)
  
  Res_up = 0.d0
  Res_dn = 0.d0
  nspecies = reaction_nw%params%nspecies
  
  allocate(molecular_diffusion_up(nspecies))
  allocate(molecular_diffusion_dn(nspecies))
  allocate(diffusivity_up(nspecies))
  allocate(diffusivity_dn(nspecies))
  harmonic_D_over_dist(:) = 0.d0
  diffusive_flux(:) = 0.d0
  
  sat_up = max(MIN_LIQ_SAT,global_auxvar_up%sat(LIQUID_PHASE))
  sat_dn = max(MIN_LIQ_SAT,global_auxvar_dn%sat(LIQUID_PHASE))
  por_up = material_auxvar_up%porosity
  por_dn = material_auxvar_dn%porosity
  ps_up = por_up*sat_up
  ps_dn = por_dn*sat_dn
  
  diffusion_coefficient => reaction_nw%diffusion_coefficient
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
                      (nwt_auxvar_dn%aqueous_eq_conc(:) &
                       - nwt_auxvar_up%aqueous_eq_conc(:))
                       
  ! Note: For dispersion, do a git pull - Glenn updated transport.F90 
  ! When adding dispersion, look at TDispersion() and the routine that
  ! calls it, UpdateTransportCoefs(), because you need to do something 
  ! with the cell centered velocities. Also, the boundary cells may need
  ! their own calculation for dispersion (There is a TDispersionBC).
                
  ! units of q = [m/s]
  q = velocity(LIQUID_PHASE)  ! liquid is the only mobile phase

  ! units of unit_n = [-] unitless
  unit_n_up = -1  ! original
  unit_n_dn = +1  ! original

  ! upstream weighting
  if (.not.bc) then
    if (q > 0.d0) then ! q flows from _up to _dn (think: upstream to downstream)
      Res_up(:) = (unit_n_up*area) * &
                   (q*sat_up*nwt_auxvar_up%aqueous_eq_conc(:) &
                    - diffusive_flux(:))
      Res_dn(:) = (unit_n_dn*area) * &
                   (q*sat_up*nwt_auxvar_up%aqueous_eq_conc(:) &
                    - diffusive_flux(:))
    else               ! q flows from _dn to _up (think: downstream to upstream)
      Res_up(:) = (unit_n_up*area) * &
                   (q*sat_dn*nwt_auxvar_dn%aqueous_eq_conc(:) &
                    - diffusive_flux(:))
      Res_dn(:) = (unit_n_dn*area) * &
                   (q*sat_dn*nwt_auxvar_dn%aqueous_eq_conc(:) &
                    - diffusive_flux(:))
    endif
  else ! boundary calculation and there is only Res_dn(:)
    if (q > 0.d0) then ! q flows into domain
      Res_dn(:) = (unit_n_dn*area) * &
                   (q*sat_up*nwt_auxvar_up%aqueous_eq_conc(:) &
                    - diffusive_flux(:))
    else               ! q flows out of domain
      Res_dn(:) = (unit_n_dn*area) * &
                   (q*sat_dn*nwt_auxvar_dn%aqueous_eq_conc(:) &
                    - diffusive_flux(:))
    endif
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

#include "petsc/finclude/petscmat.h"
  use petscmat
  use Realization_Subsurface_class
  use Grid_module
  use Option_module
  use Field_module
  use Logging_module
  use Debug_module
  use Connection_module
  use Coupler_module 
  use WIPP_Flow_Aux_module 

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr

  Mat :: J, J_num, J_ana, J_flux
  MatType :: mat_type
  PetscViewer :: viewer, viewer2 
  Vec :: xx_pert, res, res_pert, fixed_accum
  PetscReal, pointer :: res_pert_p(:), res_p(:), xx_pert_p(:)
  PetscReal, pointer :: fixed_accum_p(:)
  PetscReal, pointer :: derivative(:), derivative_accum(:), derivative_flux(:) 
  PetscReal, pointer :: derivative_rxn(:), derivative_srcsink(:) 
  PetscReal, pointer :: perturbation(:)
  PetscReal, pointer :: residual(:)
  PetscReal, pointer :: residual_pert(:)
  PetscReal, pointer :: dRes(:), dRes_accum(:), dRes_rxn(:), dRes_flux(:)
  PetscReal, pointer :: work_loc_p(:) 
  PetscReal :: res_at_cell(realization%reaction_nw%params%nspecies)
  PetscReal :: Res_up(realization%reaction_nw%params%nspecies)
  PetscReal :: Res_dn(realization%reaction_nw%params%nspecies)
  type(option_type), pointer :: option
  type(grid_type),  pointer :: grid
  type(field_type), pointer :: field
  class(reaction_nw_type), pointer :: reaction_nw
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
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
  PetscInt :: istart, iend, offset
  PetscInt :: i, k, idof, ispecies, nspecies
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: iconn, sum_connection
  PetscInt :: cell_number
  PetscInt :: iphase
  PetscReal :: Jac_accum(realization%reaction_nw%params%nspecies, &
                   realization%reaction_nw%params%nspecies)
  PetscReal :: Jac_srcsink(realization%reaction_nw%params%nspecies, &
                   realization%reaction_nw%params%nspecies)
  PetscReal :: Jac_rxn(realization%reaction_nw%params%nspecies, &
                   realization%reaction_nw%params%nspecies)
  PetscReal :: JacUp(realization%reaction_nw%params%nspecies, &
                     realization%reaction_nw%params%nspecies)
  PetscReal :: JacDn(realization%reaction_nw%params%nspecies, &
                     realization%reaction_nw%params%nspecies)
  PetscReal :: area, value, unA_up, unA_dn
  PetscReal :: compare_accum, compare_rxn, compare_flux
  PetscBool :: numerical
    
  option => realization%option
  grid => realization%patch%grid
  field => realization%field
  reaction_nw => realization%reaction_nw
  
  nwt_auxvars => realization%patch%aux%NWT%auxvars
  nwt_auxvars_bc => realization%patch%aux%NWT%auxvars_bc
  global_auxvars => realization%patch%aux%Global%auxvars
  global_auxvars_bc => realization%patch%aux%Global%auxvars_bc
  material_auxvars => realization%patch%aux%Material%auxvars
  ! note: there is no realization%patch%aux%Material%auxvars_bc
  material_auxvars_bc => realization%patch%aux%Material%auxvars
  
  ! Strictly used to check if alpha needs to be considered in the area
  nullify(wippflo_auxvars)
  if (associated(realization%patch%aux%WIPPFlo)) then
    wippflo_auxvars => realization%patch%aux%WIPPFlo%auxvars
  endif
  
  iphase = LIQUID_PHASE
  nspecies = realization%reaction_nw%params%nspecies

  ! For debugging only. Set the type of jacobian here:
  numerical = PETSC_FALSE

  call PetscLogEventBegin(logging%event_nwt_jacobian,ierr);CHKERRQ(ierr)

  call MatGetType(A,mat_type,ierr);CHKERRQ(ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  else
    J = A
  endif
    
  ! Zero out the Jacobian matrix
  call MatZeroEntries(J,ierr);CHKERRQ(ierr)

  ! For debugging viewer:
  call MatCreate(option%mycomm,J_ana,ierr);CHKERRQ(ierr)
  call MatSetSizes(J_ana,PETSC_DECIDE,PETSC_DECIDE, &
                   grid%nlmax*option%ntrandof, &
                   grid%nlmax*option%ntrandof,ierr);CHKERRQ(ierr)
  call MatSetType(J_ana,MATAIJ,ierr);CHKERRQ(ierr)
  call MatSetFromOptions(J_ana,ierr);CHKERRQ(ierr)
  call MatSetUp(J_ana,ierr);CHKERRQ(ierr)
  call MatZeroEntries(J_ana,ierr);CHKERRQ(ierr)
  ! For debugging jacobian entries:
  call MatCreate(option%mycomm,J_flux,ierr);CHKERRQ(ierr)
  call MatSetSizes(J_flux,PETSC_DECIDE,PETSC_DECIDE, &
                   grid%nlmax*option%ntrandof, &
                   grid%nlmax*option%ntrandof,ierr);CHKERRQ(ierr)
  call MatSetType(J_flux,MATAIJ,ierr);CHKERRQ(ierr)
  call MatSetFromOptions(J_flux,ierr);CHKERRQ(ierr)
  call MatSetUp(J_flux,ierr);CHKERRQ(ierr)
  call MatZeroEntries(J_flux,ierr);CHKERRQ(ierr)
  
#if 1
  !== Accumulation Terms ======================================
  if (.not.option%steady_state) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      ! ignore inactive cells with inactive materials
      if (realization%patch%imat(ghosted_id) <= 0) cycle
      
      call NWTJacobianAccum(material_auxvars(ghosted_id), &
                            reaction_nw,option,Jac_accum) 
              
      ! PETSc uses 0-based indexing so the position must be (ghosted_id-1)              
      call MatSetValuesBlockedLocal(J,1,ghosted_id-1,1,ghosted_id-1,Jac_accum, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)
      call MatSetValuesBlocked(J_ana,1,ghosted_id-1,1,ghosted_id-1,Jac_accum, &
                               ADD_VALUES,ierr);CHKERRQ(ierr)
  
    enddo
  endif
#endif

#if 1
  !== Source/Sink Terms =======================================
  source_sink => realization%patch%source_sink_list%first 
  sum_connection = 0
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      ! ignore inactive cells with inactive materials
      if (realization%patch%imat(ghosted_id) <= 0) cycle
      
      call NWTJacobianSrcSink(material_auxvars(ghosted_id),source_sink, &
                              realization%patch%ss_flow_vol_fluxes, &
                              sum_connection,reaction_nw,Jac_srcsink) 
                              
      !WRITE(*,*)  '   JacSrcSink = ', Jac(:,:)
                                
      ! PETSc uses 0-based indexing so the position must be (ghosted_id-1)
      call MatSetValuesBlockedLocal(J,1,ghosted_id-1,1,ghosted_id-1,Jac_srcsink, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)
      call MatSetValuesBlocked(J_ana,1,ghosted_id-1,1,ghosted_id-1,Jac_srcsink, &
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
    if (realization%patch%imat(ghosted_id) <= 0) cycle
    
    call NWTJacobianRx(material_auxvars(ghosted_id), &
                       reaction_nw,Jac_rxn)
    
    ! PETSc uses 0-based indexing so the position must be (ghosted_id-1)              
    call MatSetValuesBlockedLocal(J,1,ghosted_id-1,1,ghosted_id-1,Jac_rxn, &
                                  ADD_VALUES,ierr);CHKERRQ(ierr)
    call MatSetValuesBlocked(J_ana,1,ghosted_id-1,1,ghosted_id-1,Jac_rxn, &
                             ADD_VALUES,ierr);CHKERRQ(ierr)
    
  enddo
#endif
  
#if 1
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
      
      area = cur_connection_set%area(iconn)
      if (associated(wippflo_auxvars)) then
        area = area * 0.5d0 * &
               (wippflo_auxvars(ZERO_INTEGER,ghosted_id_up)%alpha + &
                wippflo_auxvars(ZERO_INTEGER,ghosted_id_dn)%alpha)
      endif
      
      if (realization%patch%imat(ghosted_id_up) <= 0 .or.  &
          realization%patch%imat(ghosted_id_dn) <= 0) cycle
          
      call NWTJacobianFlux(nwt_auxvars(ghosted_id_up), &
                      nwt_auxvars(ghosted_id_dn), &
                      global_auxvars(ghosted_id_up), &
                      global_auxvars(ghosted_id_dn), &
                      material_auxvars(ghosted_id_up), &
                      material_auxvars(ghosted_id_dn), &
                      area, &
                      cur_connection_set%dist(:,iconn), &
                      realization%patch%internal_velocities(:,sum_connection), &
                      reaction_nw,option,JacUp,JacDn)
          
      if (local_id_up>0) then
        ! PETSc uses 0-based indexing so the position must be (ghosted_id-1)
        call MatSetValuesBlockedLocal(J,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                      JacUp,ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesBlocked(J_ana,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                 JacUp,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
   
      if (local_id_dn>0) then
        ! PETSc uses 0-based indexing so the position must be (ghosted_id-1)
        call MatSetValuesBlockedLocal(J,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                                      JacDn,ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesBlocked(J_ana,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                                 JacDn,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
      
    enddo
    
  cur_connection_set => cur_connection_set%next
  enddo

  ! Boundary Flux Terms ---------------------------------------
  boundary_condition => realization%patch%boundary_condition_list%first
  sum_connection = 0 
  do 
    if (.not.associated(boundary_condition)) exit
    
    cur_connection_set => boundary_condition%connection_set
    
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      
      area = cur_connection_set%area(iconn)
      
      if (realization%patch%imat(ghosted_id) <= 0) cycle
      
      call NWTJacobianFlux(nwt_auxvars_bc(sum_connection), &
                      nwt_auxvars(ghosted_id), &
                      global_auxvars_bc(sum_connection), &
                      global_auxvars(ghosted_id), &
                      material_auxvars_bc(sum_connection), &
                      material_auxvars(ghosted_id), &
                      area, &
                      cur_connection_set%dist(:,iconn), &
                      realization%patch%boundary_velocities(:,sum_connection), &
                      reaction_nw,option,JacUp,JacDn)
      
      ! PETSc uses 0-based indexing so the position must be (ghosted_id-1)              
      call MatSetValuesBlockedLocal(J,1,ghosted_id-1,1,ghosted_id-1,JacDn, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)
      call MatSetValuesBlocked(J_ana,1,ghosted_id-1,1,ghosted_id-1,JacDn, &
                               ADD_VALUES,ierr);CHKERRQ(ierr)
      ! note: Don't need to worry about JacUp because that is outside of
      ! the domain, and doesn't have a place in A.
      
    enddo
    
  boundary_condition => boundary_condition%next
  enddo
#endif

  call MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr) 

  call MatAssemblyBegin(J_ana,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(J_ana,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  !== ?????? ==================================================
  
  if (reaction_nw%use_log_formulation) then
    call VecGetArrayF90(realization%field%tran_work_loc,work_loc_p, &
                        ierr);CHKERRQ(ierr)
    do ghosted_id = 1, grid%ngmax  
      offset = (ghosted_id-1)*reaction_nw%params%ncomp
      istart = offset + 1
      iend = offset + reaction_nw%params%ncomp
      
      if (realization%patch%imat(ghosted_id) <= 0) then
        work_loc_p(istart:iend) = 1.d0
      else
        work_loc_p(istart:iend) = nwt_auxvars(ghosted_id)%total_bulk_conc(:)
      endif
    enddo
    call VecRestoreArrayF90(realization%field%tran_work_loc, work_loc_p,  &
                            ierr);CHKERRQ(ierr)
  endif 
    
  if (realization%debug%matview_Jacobian) then
    string = 'NWTjacobianAN'
    call DebugCreateViewer(realization%debug,string,realization%option,viewer)
    call MatView(J_ana,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  if (realization%reaction_nw%use_log_formulation) then
    call MatDiagonalScaleLocal(J,realization%field%tran_work_loc, &
                               ierr);CHKERRQ(ierr)
    if (realization%debug%matview_Jacobian) then
      string = 'NWTjacobianAN'
      call DebugCreateViewer(realization%debug,string,realization%option,viewer)
      call MatView(J,viewer,ierr);CHKERRQ(ierr)
      call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
    endif
  endif

  ! numerical jacobian:
  ! for debugging ONLY!! assumes no decay, and single processor 
#if 1
  call MatCreate(option%mycomm,J_num,ierr);CHKERRQ(ierr)
  call MatSetSizes(J_num,PETSC_DECIDE,PETSC_DECIDE, &
                   grid%nlmax*option%ntrandof, &
                   grid%nlmax*option%ntrandof,ierr);CHKERRQ(ierr)
  call MatSetType(J_num,MATAIJ,ierr);CHKERRQ(ierr)
  call MatSetFromOptions(J_num,ierr);CHKERRQ(ierr)
  call MatSetUp(J_num,ierr);CHKERRQ(ierr)
  call MatZeroEntries(J_num,ierr);CHKERRQ(ierr)

  ! Prepare Vecs for numerical Jacobian
  call VecDuplicate(field%tran_xx,xx_pert,ierr);CHKERRQ(ierr)
  call VecDuplicate(field%tran_xx,res,ierr);CHKERRQ(ierr)
  call VecDuplicate(field%tran_xx,res_pert,ierr);CHKERRQ(ierr)
  call VecDuplicate(field%tran_accum,fixed_accum,ierr);CHKERRQ(ierr)

  allocate(derivative(grid%nlmax*option%ntrandof))
  allocate(derivative_accum(grid%nlmax*option%ntrandof))
  allocate(derivative_flux(grid%nlmax*option%ntrandof))
  allocate(derivative_rxn(grid%nlmax*option%ntrandof))
  allocate(derivative_srcsink(grid%nlmax*option%ntrandof))
  allocate(perturbation(grid%nlmax*option%ntrandof))
  allocate(residual(grid%nlmax*option%ntrandof))
  allocate(residual_pert(grid%nlmax*option%ntrandof))
  allocate(dRes(grid%nlmax*option%ntrandof))
  allocate(dRes_accum(grid%nlmax*option%ntrandof))
  allocate(dRes_rxn(grid%nlmax*option%ntrandof))
  allocate(dRes_flux(grid%nlmax*option%ntrandof))

  cell_number = 3

  ! get the dM value:
  call VecCopy(field%tran_xx,xx_pert,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(xx_pert,xx_pert_p,ierr);CHKERRQ(ierr)
  do i=1,size(xx_pert_p)
    perturbation(i) = max(xx_pert_p(i)*perturbation_tolerance,1.0d-30)
    xx_pert_p(i) = xx_pert_p(i)+perturbation(i)
  enddo

  call VecCopy(field%tran_accum,fixed_accum,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(fixed_accum,fixed_accum_p,ierr);CHKERRQ(ierr)
  residual(:) = 0.d0
  residual_pert(:) = 0.d0
  do i=1,grid%nlmax
    offset = (i-1)*realization%reaction_nw%params%nspecies
    istart = offset + 1
    iend = offset + realization%reaction_nw%params%nspecies
    ! note: auxvars should be current 
    ! get the unperturbed residual:
    call NWTResidualAccum(nwt_auxvars(i),global_auxvars(i), &
                          material_auxvars(i),reaction_nw,res_at_cell)
    residual(istart:iend) = &
      (res_at_cell(1:realization%reaction_nw%params%nspecies) &
        - fixed_accum_p(istart:iend))/option%tran_dt 
    ! use dM to get a perturbed residual:
    nwt_auxvars(i)%total_bulk_conc = xx_pert_p(istart:iend)
    call NWTAuxVarCompute(nwt_auxvars(i),global_auxvars(i), &
                          material_auxvars(i),reaction_nw,option)
    call NWTResidualAccum(nwt_auxvars(i),global_auxvars(i), &
                          material_auxvars(i),reaction_nw,res_at_cell)
    residual_pert(istart:iend) = &
      (res_at_cell(1:realization%reaction_nw%params%nspecies) &
        - fixed_accum_p(istart:iend))/option%tran_dt 
    ! calculate the numerical Jacobian entry:
    dRes_accum(istart:iend) = residual_pert(istart:iend) &
                              - residual(istart:iend)
    derivative_accum(istart:iend) = dRes_accum(istart:iend) / &
                                    perturbation(istart:iend)
  enddo

  ! call Residual routine with original unperturbed M to set auxvars back
  call NWTResidual(PETSC_NULL_SNES,field%tran_xx,res,realization,ierr)
  residual(:) = 0.d0
  residual_pert(:) = 0.d0
  do i=1,grid%nlmax
    offset = (i-1)*realization%reaction_nw%params%nspecies
    istart = offset + 1
    iend = offset + realization%reaction_nw%params%nspecies
    ! note: auxvars should be current 
    ! get the unperturbed residual:
    call NWTResidualRx(nwt_auxvars(i),material_auxvars(i), &
                       reaction_nw,res_at_cell)
    residual(istart:iend) = &
                      res_at_cell(1:realization%reaction_nw%params%nspecies)
    ! use dM to get a perturbed residual:
    nwt_auxvars(i)%total_bulk_conc = xx_pert_p(istart:iend)
    call NWTAuxVarCompute(nwt_auxvars(i),global_auxvars(i), &
                          material_auxvars(i),reaction_nw,option)
    call NWTResidualRx(nwt_auxvars(i),material_auxvars(i),reaction_nw, &
                       res_at_cell)
    residual_pert(istart:iend) = &
                      res_at_cell(1:realization%reaction_nw%params%nspecies) 
    ! calculate the numerical Jacobian entry:
    dRes_rxn(istart:iend) = residual_pert(istart:iend) &
                            - residual(istart:iend)
    derivative_rxn(istart:iend) = dRes_rxn(istart:iend) / &
                                  perturbation(istart:iend)     
  enddo

  ! call Residual routine with original unperturbed M to set auxvars back
  call NWTResidual(PETSC_NULL_SNES,field%tran_xx,res,realization,ierr)
  residual(:) = 0.d0
  residual_pert(:) = 0.d0
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
      
      area = cur_connection_set%area(iconn)
      if (associated(wippflo_auxvars)) then
        area = area * 0.5d0 * &
               (wippflo_auxvars(ZERO_INTEGER,ghosted_id_up)%alpha + &
                wippflo_auxvars(ZERO_INTEGER,ghosted_id_dn)%alpha)
      endif
          
      call NWTResidualFlux(nwt_auxvars(ghosted_id_up), &
                      nwt_auxvars(ghosted_id_dn), &
                      global_auxvars(ghosted_id_up), &
                      global_auxvars(ghosted_id_dn), &
                      material_auxvars(ghosted_id_up), &
                      material_auxvars(ghosted_id_dn), &
                      area, &
                      cur_connection_set%dist(:,iconn), &
                      realization%patch%internal_velocities(:,sum_connection), &
                      reaction_nw,option,PETSC_FALSE,Res_up,Res_dn)
                            
      if (local_id_up>0) then
        offset = (local_id_up-1)*realization%reaction_nw%params%nspecies
        istart = offset + 1
        iend = offset + realization%reaction_nw%params%nspecies
        residual(istart:iend) = residual(istart:iend) &
                              + Res_up(1:realization%reaction_nw%params%nspecies)
      endif
      
      if (local_id_dn>0) then
        offset = (local_id_dn-1)*realization%reaction_nw%params%nspecies
        istart = offset + 1
        iend = offset + realization%reaction_nw%params%nspecies
        residual(istart:iend) = residual(istart:iend) &
                              + Res_dn(1:realization%reaction_nw%params%nspecies)
      endif

    enddo
    cur_connection_set => cur_connection_set%next
  enddo
  ! use dM to get a perturbed residual:
  do i=1,grid%nlmax
    offset = (i-1)*realization%reaction_nw%params%nspecies
    istart = offset + 1
    iend = offset + realization%reaction_nw%params%nspecies
    nwt_auxvars(i)%total_bulk_conc = xx_pert_p(istart:iend)
    call NWTAuxVarCompute(nwt_auxvars(i),global_auxvars(i), &
                          material_auxvars(i),reaction_nw,option)
  enddo
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
      
      area = cur_connection_set%area(iconn)
      if (associated(wippflo_auxvars)) then
        area = area * 0.5d0 * &
               (wippflo_auxvars(ZERO_INTEGER,ghosted_id_up)%alpha + &
                wippflo_auxvars(ZERO_INTEGER,ghosted_id_dn)%alpha)
      endif
          
      call NWTResidualFlux(nwt_auxvars(ghosted_id_up), &
                      nwt_auxvars(ghosted_id_dn), &
                      global_auxvars(ghosted_id_up), &
                      global_auxvars(ghosted_id_dn), &
                      material_auxvars(ghosted_id_up), &
                      material_auxvars(ghosted_id_dn), &
                      area, &
                      cur_connection_set%dist(:,iconn), &
                      realization%patch%internal_velocities(:,sum_connection), &
                      reaction_nw,option,PETSC_FALSE,Res_up,Res_dn)
                            
      if (local_id_up>0) then
        offset = (local_id_up-1)*realization%reaction_nw%params%nspecies
        istart = offset + 1
        iend = offset + realization%reaction_nw%params%nspecies
        residual_pert(istart:iend) = residual_pert(istart:iend) &
                                   + Res_up(1:realization%reaction_nw%params%nspecies)
      endif
      
      if (local_id_dn>0) then
        offset = (local_id_dn-1)*realization%reaction_nw%params%nspecies
        istart = offset + 1
        iend = offset + realization%reaction_nw%params%nspecies
        residual_pert(istart:iend) = residual_pert(istart:iend) &
                                   + Res_dn(1:realization%reaction_nw%params%nspecies)
      endif

    enddo
    cur_connection_set => cur_connection_set%next
  enddo
  ! calculate the numerical Jacobian entry:
  do i=1,grid%nlmax
    offset = (i-1)*nspecies
    istart = offset + 1
    iend = offset + nspecies
    dRes_flux(istart:iend) = residual_pert(istart:iend) &
                             - residual(istart:iend)
    derivative_flux(istart:iend) = dRes_flux(istart:iend) / &
                                   perturbation(istart:iend) 
    if (i == cell_number) WRITE(*,*) 'flux dRes = ', dRes_flux(istart)
    if (i == cell_number) WRITE(*,*) 'flux dM = ', perturbation(istart)   
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  residual(:) = 0.d0
  residual_pert(:) = 0.d0
  ! get the unperturbed residual:
  call NWTResidual(PETSC_NULL_SNES,field%tran_xx,res,realization,ierr)
  call VecGetArrayF90(res,res_p,ierr);CHKERRQ(ierr)
  ! get the dM value:
  do i=1,grid%nlmax
    nwt_auxvars(i)%total_bulk_conc = xx_pert_p(istart:iend)
  enddo
  ! use dM to get a perturbed residual:
  call NWTResidual(PETSC_NULL_SNES,xx_pert,res_pert,realization,ierr)
  call VecGetArrayF90(res_pert,res_pert_p,ierr);CHKERRQ(ierr)
  ! calculate the numerical Jacobian entry:
  do i=1,size(res_pert_p)
    dRes(i) = res_pert_p(i)-res_p(i)
    derivative(i) = dRes(i)/perturbation(i)
  enddo
  ! call Residual routine with original unperturbed M to set auxvars back
  call NWTResidual(PETSC_NULL_SNES,field%tran_xx,res,realization,ierr)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
      
      area = cur_connection_set%area(iconn)
      if (associated(wippflo_auxvars)) then
        area = area * 0.5d0 * &
               (wippflo_auxvars(ZERO_INTEGER,ghosted_id_up)%alpha + &
                wippflo_auxvars(ZERO_INTEGER,ghosted_id_dn)%alpha)
      endif
                
      call NWTJacobianFlux(nwt_auxvars(ghosted_id_up), &
                      nwt_auxvars(ghosted_id_dn), &
                      global_auxvars(ghosted_id_up), &
                      global_auxvars(ghosted_id_dn), &
                      material_auxvars(ghosted_id_up), &
                      material_auxvars(ghosted_id_dn), &
                      area, &
                      cur_connection_set%dist(:,iconn), &
                      realization%patch%internal_velocities(:,sum_connection), &
                      reaction_nw,option,JacUp,JacDn)
      
      if (ghosted_id_up == cell_number) then
        if (realization%patch%internal_velocities(LIQUID_PHASE,sum_connection) > 0.d0) then
          unA_up = realization%patch%internal_velocities(LIQUID_PHASE,sum_connection)/ &
                material_auxvars(ghosted_id_up)%porosity * area * -1.d0
        else 
          unA_up = realization%patch%internal_velocities(LIQUID_PHASE,sum_connection)/ &
                material_auxvars(ghosted_id_dn)%porosity * area * -1.d0
        endif
        WRITE(*,*) 'u*A.n = ', unA_up
        !WRITE(*,*) 'JacUp = ', JacUp(1,1)
      endif
      if (ghosted_id_dn == cell_number) then
        if (realization%patch%internal_velocities(LIQUID_PHASE,sum_connection) > 0.d0) then
          unA_dn = realization%patch%internal_velocities(LIQUID_PHASE,sum_connection)/ &
                material_auxvars(ghosted_id_up)%porosity * area * 1.d0
        else 
          unA_dn = realization%patch%internal_velocities(LIQUID_PHASE,sum_connection)/ &
                material_auxvars(ghosted_id_dn)%porosity * area * 1.d0
        endif
        WRITE(*,*) 'u*A.n = ', unA_dn
        !WRITE(*,*) 'JacDn = ', JacDn(1,1)
      endif

      k = 0
      if (local_id_up>0) then
        offset = (local_id_up-1)*nspecies
        istart = offset + 1
        iend = offset + nspecies
        do i=istart,iend
          k = k + 1
          !WRITE(*,*) 'Adding value', JacUp(k,k), 'in diagonal position', i-1
          ! PETSc uses 0-based indexing so the position must be (ghosted_id-1)
          call MatSetValues(J_flux,1,i-1,1,i-1,JacUp(k,k),ADD_VALUES,ierr);CHKERRQ(ierr)
        enddo
      endif
      k = 0
      if (local_id_dn>0) then
        offset = (local_id_dn-1)*nspecies
        istart = offset + 1
        iend = offset + nspecies
        do i=istart,iend
          k = k + 1
          !WRITE(*,*) 'Adding value', JacUp(k,k), 'in diagonal position', i-1
          ! PETSc uses 0-based indexing so the position must be (ghosted_id-1)
          call MatSetValues(J_flux,1,i-1,1,i-1,JacDn(k,k),ADD_VALUES,ierr);CHKERRQ(ierr)
        enddo
      endif
      
    enddo

  cur_connection_set => cur_connection_set%next
  enddo
  call MatAssemblyBegin(J_flux,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(J_flux,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  WRITE(*,*) 'net sum = ', unA_dn + unA_up

  do i=1,grid%nlmax
    offset = (i-1)*nspecies
    istart = offset + 1
    iend = offset + nspecies
    call MatGetValues(J_ana,1,istart-1,1,istart-1,value,ierr);CHKERRQ(ierr)
    if (i == cell_number) WRITE(*,*) '[num]jac_whole =', derivative(istart), &
     '[an]compare = ', value, '% diff =', (derivative(istart)-value)/value*100.d0
    call NWTJacobianAccum(material_auxvars(i),reaction_nw,option,Jac_accum)
    compare_accum = Jac_accum(1,1)
    if (i == cell_number) WRITE(*,*) '[num]jac_accum =', derivative_accum(istart), &
     '[an]compare = ', compare_accum, '% diff =', &
     (derivative_accum(cell_number)-compare_accum)/compare_accum*100.d0
    call NWTJacobianRx(material_auxvars(i),reaction_nw,Jac_rxn)
    compare_rxn = Jac_rxn(1,1)
    if (i == cell_number) WRITE(*,*) '[num]jac_rxn =', derivative_rxn(istart), &
     '[an]compare = ', compare_rxn, '% diff =', &
     (derivative_rxn(istart)-compare_rxn)/compare_rxn*100.d0
    call MatGetValues(J_flux,1,istart-1,1,istart-1,value,ierr);CHKERRQ(ierr)
    compare_flux = value
    if (i == cell_number) WRITE(*,*) '[num]jac_flux =', derivative_flux(istart), &
     '[an]compare = ', compare_flux, '% diff =', &
     (derivative_flux(istart)-compare_flux)/compare_flux*100.d0  
  enddo

  call VecRestoreArrayF90(xx_pert,xx_pert_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(res_pert,res_pert_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(res,res_p,ierr);CHKERRQ(ierr)

  if (numerical) call MatZeroEntries(J,ierr);CHKERRQ(ierr)

  ! note: off-diagonal terms from decay are not included!
  !       need to figure out how to do this
  do local_id = 1,grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    ! ignore inactive cells with inactive materials
    if (realization%patch%imat(ghosted_id) <= 0) cycle
    do ispecies = 1,option%ntrandof
      idof = (local_id-1)*ispecies + ispecies
      
      if (dabs(derivative(idof)) > 1.d-50) then
        call MatSetValue(J_num,idof-1,idof-1,derivative(idof),insert_values, &
                         ierr);CHKERRQ(ierr)
      endif
       if (numerical) then
         if (dabs(derivative(idof)) > 1.d-50) then
           call MatSetValue(J,idof-1,idof-1,derivative(idof),insert_values, &
                            ierr);CHKERRQ(ierr)
         else
           call MatSetValue(J,idof-1,idof-1,0.d0,insert_values, &
                            ierr);CHKERRQ(ierr)
         endif
       endif

    enddo
  enddo

  if (numerical) then
    call MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  endif
  call MatAssemblyBegin(J_num,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(J_num,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  call VecDestroy(xx_pert,ierr);CHKERRQ(ierr)
  call VecDestroy(res,ierr);CHKERRQ(ierr)
  call VecDestroy(res_pert,ierr);CHKERRQ(ierr)

  if (realization%debug%matview_Jacobian) then
    string = 'NWTjacobianNUM'
    call DebugCreateViewer(realization%debug,string,realization%option,viewer2)
    call MatView(J_num,viewer2,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer2,ierr);CHKERRQ(ierr)
  endif

  if (realization%debug%matview_Jacobian) then
    string = 'NWTjacobianNUMFLUX'
    call DebugCreateViewer(realization%debug,string,realization%option,viewer2)
    call MatView(J_flux,viewer2,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer2,ierr);CHKERRQ(ierr)
  endif

  call MatDestroy(J_num,ierr);CHKERRQ(ierr)
  call MatDestroy(J_flux,ierr);CHKERRQ(ierr)

#endif

  call MatDestroy(J_ana,ierr);CHKERRQ(ierr)
  call PetscLogEventEnd(logging%event_nwt_jacobian,ierr);CHKERRQ(ierr)
  
end subroutine NWTJacobian

! ************************************************************************** !

subroutine NWTJacobianAccum(material_auxvar,reaction_nw,option,Jac)
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
  class(reaction_nw_type) :: reaction_nw
  type(option_type) :: option
  PetscReal :: Jac(reaction_nw%params%nspecies,reaction_nw%params%nspecies)
  
  PetscReal :: vol_dt
  PetscInt :: istart, iend, ispecies
  
  Jac = 0.d0
  
  ! units of volume = [m^3-bulk]
  ! units of tran_dt = [sec]
  vol_dt = material_auxvar%volume/option%tran_dt
  
  istart = 1
  iend = reaction_nw%params%nspecies
  do ispecies=istart,iend
    ! units of Jac = [m^3-bulk/sec]
    Jac(ispecies,ispecies) = vol_dt
  enddo

end subroutine NWTJacobianAccum

! ************************************************************************** !

subroutine NWTJacobianSrcSink(material_auxvar,source_sink,ss_flow_vol_fluxes, &
                              sum_connection,reaction_nw,Jac)
  ! 
  ! Computes the source/sink terms in the Jacobian matrix.
  ! All Jacobian entries should be in [m^3-bulk/sec].
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/14/2019
  ! 

  use Coupler_module

  implicit none
  
  class(material_auxvar_type) :: material_auxvar
  type(coupler_type), pointer :: source_sink
  PetscReal, pointer :: ss_flow_vol_fluxes(:,:)
  PetscInt :: sum_connection
  class(reaction_nw_type), pointer :: reaction_nw
  PetscReal :: Jac(reaction_nw%params%nspecies,reaction_nw%params%nspecies)
  
  PetscInt :: istart, iend, ispecies
  PetscReal :: qsrc, u, vol
  PetscReal :: coef_in
  PetscBool :: dry_out
  
  Jac = 0.d0
  
  if (associated(ss_flow_vol_fluxes)) then
    ! qsrc = [m^3-liq/sec] 
    qsrc = ss_flow_vol_fluxes(LIQUID_PHASE,sum_connection)
  endif
  
  ! transform qsrc into pore velocity volumetric flow
  ! units of porosity = [m^3-void/m^3-bulk]
  ! units of u = [m^3-bulk/sec]
  u = qsrc/material_auxvar%porosity
  
  ! units of vol = [m^3-bulk]
  vol = material_auxvar%volume
  
  ! -- Aqueous-Component ----------------------------------------
  select case(source_sink%tran_condition%itype)
    case(EQUILIBRIUM_SS)
      !TODO(jenn) What is EQUILIBRIUM_SS option?
    case(MASS_RATE_SS)
      !TODO(jenn) What is MASS_RATE_SS option?
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
  ! units of vol = [m^3-bulk]
  istart = 1
  iend = reaction_nw%params%nspecies
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

subroutine NWTJacobianRx(material_auxvar,reaction_nw,Jac)
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
  class(reaction_nw_type) :: reaction_nw
  PetscReal :: Jac(reaction_nw%params%nspecies,reaction_nw%params%nspecies)
  
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
  iend = reaction_nw%params%nspecies
  !TODO(jenn): convert from linked list to array format
  do ispecies=istart,iend
    decay_rate = 0.d0
    parent_decay_rate = 0.d0
    has_parent = PETSC_FALSE
    ! find the decay rate for species_id = ispecies (if its radioactive)
    rad_rxn => reaction_nw%rad_decay_rxn_list
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
    Jac(ispecies,ispecies) = 1.d0*(-1.d0*vol*decay_rate)
    
    ! fill in the off-diagonal associated with ingrowth from a parent
    if (has_parent) then
      Jac(ispecies,parent_id) = 1.d0*(vol*parent_decay_rate)
    endif
    
  enddo

end subroutine NWTJacobianRx

! ************************************************************************** !

subroutine NWTJacobianFlux(nwt_auxvar_up,nwt_auxvar_dn, &
                           global_auxvar_up,global_auxvar_dn, &
                           material_auxvar_up,material_auxvar_dn, &
                           area,dist,velocity,reaction_nw,option,Jac_up,Jac_dn)
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
  class(reaction_nw_type), pointer :: reaction_nw
  type(option_type) :: option
  PetscReal :: Jac_up(reaction_nw%params%nspecies,reaction_nw%params%nspecies)
  PetscReal :: Jac_dn(reaction_nw%params%nspecies,reaction_nw%params%nspecies)
  
  PetscInt :: unit_n_up, unit_n_dn
  PetscInt :: nspecies, ispecies
  PetscBool :: dry_out_up, dry_out_dn
  PetscReal :: q, u
  PetscReal :: harmonic_ps
  PetscReal :: sat_up, sat_dn
  PetscReal :: por_up, por_dn
  PetscReal :: ps_up, ps_dn
  PetscReal :: dist_up, dist_dn
  PetscReal, pointer :: diffusivity_up(:)
  PetscReal, pointer :: diffusivity_dn(:)
  PetscReal, pointer :: molecular_diffusion_up(:)
  PetscReal, pointer :: molecular_diffusion_dn(:)
  PetscReal, pointer :: diffusion_coefficient(:,:)
  PetscReal :: distance_gravity, upwind_weight ! both are dummy variables
  PetscReal :: harmonic_D_over_dist(reaction_nw%params%nspecies)
  
  Jac_up = 0.d0
  Jac_dn = 0.d0
  nspecies = reaction_nw%params%nspecies
  
  allocate(molecular_diffusion_up(nspecies))
  allocate(molecular_diffusion_dn(nspecies))
  allocate(diffusivity_up(nspecies))
  allocate(diffusivity_dn(nspecies))
  harmonic_D_over_dist(:) = 0.d0
  
  sat_up = max(MIN_LIQ_SAT,global_auxvar_up%sat(LIQUID_PHASE))
  sat_dn = max(MIN_LIQ_SAT,global_auxvar_dn%sat(LIQUID_PHASE))
  por_up = material_auxvar_up%porosity
  por_dn = material_auxvar_dn%porosity
  ps_up = por_up*sat_up
  ps_dn = por_dn*sat_dn
  
  diffusion_coefficient => reaction_nw%diffusion_coefficient
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
   
  if (sat_up > 0.d0) then
    dry_out_up = PETSC_FALSE
  else
    dry_out_up = PETSC_TRUE
  endif
 if (sat_dn > 0.d0) then
    dry_out_dn = PETSC_FALSE
  else
    dry_out_dn = PETSC_TRUE
  endif
  
  ! units of unit_n = [-] unitless
  unit_n_up = -1  ! original
  unit_n_dn = +1  ! original 
                
  ! units of q = [m-liq/s]
  ! units of u = [m-bulk/s]
  q = velocity(LIQUID_PHASE)  ! liquid is the only mobile phase
  
  if (q == 0.d0) then
    u = 0.d0
  else
    ! upwinding for porosity
    if (q > 0.d0) then ! q flows from _up to _dn
      u = q/por_up
    else               ! q flows from _dn to _up
      u = q/por_dn
    endif
    !if (dry_out_up .and. dry_out_dn) then
    !  ! this situation means both cells at the connection face are dry
    !  ! this situation will probably never get entered, because q should be 
    !  ! zero, and it would have got caught above
    !  u = 0.d0
    !elseif ((.not.dry_out_up) .and. (.not.dry_out_dn)) then
    !  ! this situation is the typical one, where both cells have some liquid
    !  harmonic_ps = (ps_up*ps_dn)/(ps_up*dist_up + ps_dn*dist_dn)* &
    !                (dist_up+dist_dn)   ! weighted harmonic average
    !  u = q/harmonic_ps
    !else
    !  ! one of the cells on either side of the face is dried out
    !  ! calculate u with upwinding, not harmonic average
    !  if (q > 0.d0) then ! q flows from _up to _dn
    !    if (dry_out_up) then
    !      u = 0.d0
    !    else
    !      u = q/(por_up*sat_up)
    !    endif
    !  else               ! q flows from _dn to _up
    !    if (dry_out_dn) then
    !      u = 0.d0
    !    else
    !      u = q/(por_dn*sat_dn)
    !    endif
    !  endif
    !endif
  endif
  
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

subroutine NWTZeroMassBalanceDelta(realization)
  ! 
  ! Zero's out the mass balance delta array.
  ! 
  ! Author: Jennifer Frederick
  ! Date: 06/30/2020
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
    patch%aux%NWT%auxvars(iconn)%mass_balance_delta = 0.d0
  enddo
#endif

  do iconn = 1, patch%aux%NWT%num_aux_bc
    nwt_auxvars_bc(iconn)%mass_balance_delta = 0.d0
  enddo

  do iconn = 1, patch%aux%NWT%num_aux_ss
    nwt_auxvars_ss(iconn)%mass_balance_delta = 0.d0
  enddo

end subroutine NWTZeroMassBalanceDelta

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
  use Field_module
  use Grid_module


  type(realization_subsurface_type) :: realization
  PetscInt :: max_size
  PetscReal :: sum_mol(max_size,4)
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  class(reaction_nw_type), pointer :: reaction_nw

  PetscReal :: sum_mol_tot(max_size)
  PetscReal :: sum_mol_aq(max_size)
  PetscReal :: sum_mol_sb(max_size)
  PetscReal :: sum_mol_mnrl(max_size)

  PetscErrorCode :: ierr
  PetscInt :: local_id, ghosted_id
  PetscInt :: nspecies
  PetscReal :: liquid_saturation, porosity, volume

  option => realization%option
  grid => realization%patch%grid
  field => realization%field

  reaction_nw => realization%reaction_nw

  nwt_auxvars => realization%patch%aux%NWT%auxvars
  global_auxvars => realization%patch%aux%Global%auxvars
  material_auxvars => realization%patch%aux%Material%auxvars

  sum_mol = 0.d0
  sum_mol_tot = 0.d0
  sum_mol_aq = 0.d0
  sum_mol_sb = 0.d0
  sum_mol_mnrl = 0.d0

  nspecies = reaction_nw%params%nspecies

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    ! ignore inactive cells with inactive materials
    if (realization%patch%imat(ghosted_id) <= 0) cycle
    
    liquid_saturation = global_auxvars(ghosted_id)%sat(LIQUID_PHASE)
    porosity = material_auxvars(ghosted_id)%porosity
    volume = material_auxvars(ghosted_id)%volume ! [m^3]
    
    ! aqueous (sum_mol_aq) [mol]
    sum_mol_aq(1:nspecies) = sum_mol_aq(1:nspecies) + &
           nwt_auxvars(ghosted_id)%aqueous_eq_conc(:) * &  ! [mol/m^3-liq]
           liquid_saturation*porosity*volume               ! [m^3-liq]

    ! equilibrium sorption (sum_mol_sb) [mol]
    sum_mol_sb(1:nspecies) = sum_mol_sb(1:nspecies) + &
            nwt_auxvars(ghosted_id)%sorb_eq_conc(:) * &    ! [mol/m^3-bulk]
            volume                                         ! [m^3-bulk]

    ! mineral volume fractions (sum_mol_mnrl) [mol]
    sum_mol_mnrl(1:nspecies) = sum_mol_mnrl(1:nspecies) + &
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
  use Grid_module
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars_bc(:)
  type(nw_transport_auxvar_type), pointer :: nwt_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option

  nwt_auxvars_bc => realization%patch%aux%NWT%auxvars_bc
  nwt_auxvars_ss => realization%patch%aux%NWT%auxvars_ss

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, realization%patch%aux%NWT%num_aux
    realization%patch%aux%NWT%auxvars(iconn)%mass_balance = &
      realization%patch%aux%NWT%auxvars(iconn)%mass_balance + &
      realization%patch%aux%NWT%auxvars(iconn)%mass_balance_delta*option%tran_dt
  enddo
#endif

  do iconn = 1, realization%patch%aux%NWT%num_aux_bc
    nwt_auxvars_bc(iconn)%mass_balance = &
      nwt_auxvars_bc(iconn)%mass_balance + &
      nwt_auxvars_bc(iconn)%mass_balance_delta*option%tran_dt
  enddo

  do iconn = 1, realization%patch%aux%NWT%num_aux_ss
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
  !       and the realization reaction_nw object is destroyed in 
  !       realization_subsurface.F90, RealizationStrip().

end subroutine NWTDestroy

! ************************************************************************** !

end module NW_Transport_module
