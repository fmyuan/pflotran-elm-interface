module WIPP_Flow_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use WIPP_Flow_Aux_module
  use WIPP_Flow_Common_module
  use Global_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

  public :: WIPPFloSetup, &
            WIPPFloInitializeTimestep, &
            WIPPFloUpdateSolution, &
            WIPPFloTimeCut,&
            WIPPFloUpdateAuxVars, &
            WIPPFloUpdateFixedAccum, &
            WIPPFloComputeMassBalance, &
            WIPPFloResidual, &
            WIPPFloJacobian, &
            WIPPFloSetPlotVariables, &
            WIPPFloMapBCAuxVarsToGlobal, &
            WIPPFloDestroy

contains

! ************************************************************************** !

subroutine WIPPFloSetup(realization)
  ! 
  ! Creates arrays for auxiliary variables
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
  use Material_Aux_class
  use Output_Aux_module
 
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
  PetscReal :: dist(3)
  PetscBool :: error_found
  PetscInt :: flag(10)
                                                ! extra index for derivatives
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars_bc(:)
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars_ss(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  
  patch%aux%WIPPFlo => WIPPFloAuxCreate(option)

  wippflo_analytical_derivatives = .not.option%flow%numerical_derivatives

  ! ensure that material properties specific to this module are properly
  ! initialized
  material_parameter => patch%aux%Material%material_parameter
  error_found = PETSC_FALSE
  if (minval(material_parameter%soil_residual_saturation(:,:)) < 0.d0) then
    option%io_buffer = 'ERROR: Non-initialized soil residual saturation.'
    call printMsg(option)
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
      call printMsg(option)
    endif
    if (material_auxvars(ghosted_id)%porosity < 0.d0 .and. flag(2) == 0) then
      flag(2) = 1
      option%io_buffer = 'ERROR: Non-initialized porosity.'
      call printMsg(option)
    endif
    if (minval(material_auxvars(ghosted_id)%permeability) < 0.d0 .and. &
        flag(5) == 0) then
      option%io_buffer = 'ERROR: Non-initialized permeability.'
      call printMsg(option)
      flag(5) = 1
    endif
  enddo
  
  if (error_found .or. maxval(flag) > 0) then
    option%io_buffer = 'Material property errors found in WIPPFloSetup.'
    call printErrMsg(option)
  endif
  
  ! allocate auxvar data structures for all grid cells  
  if (wippflo_analytical_derivatives) then
    ndof = 0
  else
    ndof = option%nflowdof
  endif
  allocate(wippflo_auxvars(0:ndof,grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    do idof = 0, ndof
      call WIPPFloAuxVarInit(wippflo_auxvars(idof,ghosted_id), &
                         (wippflo_analytical_derivatives .and. idof==0), &
                          option)
    enddo
  enddo
  patch%aux%WIPPFlo%auxvars => wippflo_auxvars
  patch%aux%WIPPFlo%num_aux = grid%ngmax

  ! count the number of boundary connections and allocate
  ! auxvar data structures for them 
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_condition_list)
  if (sum_connection > 0) then
    allocate(wippflo_auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call WIPPFloAuxVarInit(wippflo_auxvars_bc(iconn),PETSC_FALSE,option)
    enddo
    patch%aux%WIPPFlo%auxvars_bc => wippflo_auxvars_bc
  endif
  patch%aux%WIPPFlo%num_aux_bc = sum_connection

  ! count the number of source/sink connections and allocate
  ! auxvar data structures for them  
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (sum_connection > 0) then
    allocate(wippflo_auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call WIPPFloAuxVarInit(wippflo_auxvars_ss(iconn),PETSC_FALSE,option)
    enddo
    patch%aux%WIPPFlo%auxvars_ss => wippflo_auxvars_ss
  endif
  patch%aux%WIPPFlo%num_aux_ss = sum_connection

  ! allocate arrays for upwind directions
  ! internal connections
  connection_set_list => grid%internal_connection_set_list
  sum_connection = ConnectionGetNumberInList(connection_set_list)
  allocate(patch%aux%WIPPFlo%upwind_direction(option%nphase,sum_connection))
  cur_connection_set => connection_set_list%first
  sum_connection = 0  
  do 
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      dist(1:3) = dabs(cur_connection_set%dist(1:3,iconn))
      if (dist(1) > dist(2) .and. dist(1) > dist(3)) then
        i = X_DIRECTION
      else if (dist(2) > dist(1) .and. dist(2) > dist(3)) then
        i = Y_DIRECTION
      else if (dist(3) > dist(1) .and. dist(3) > dist(2)) then
        i = Z_DIRECTION
      endif
      patch%aux%WIPPFlo%upwind_direction(:,sum_connection) = i
    enddo
    cur_connection_set => cur_connection_set%next
  enddo
  ! boundary connections
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_condition_list)
  boundary_condition => patch%boundary_condition_list%first
  allocate(patch%aux%WIPPFlo%upwind_direction_bc(option%nphase,sum_connection))
  sum_connection = 0  
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      dist(1:3) = dabs(cur_connection_set%dist(1:3,iconn))
      if (dist(1) > dist(2) .and. dist(1) > dist(3)) then
        i = X_DIRECTION
      else if (dist(2) > dist(1) .and. dist(2) > dist(3)) then
        i = Y_DIRECTION
      else if (dist(3) > dist(1) .and. dist(3) > dist(2)) then
        i = Z_DIRECTION
      endif
      patch%aux%WIPPFlo%upwind_direction_bc(:,sum_connection) = i
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! create array for zeroing Jacobian entries if inactive
  allocate(patch%aux%WIPPFlo%row_zeroing_array(grid%nlmax))
  patch%aux%WIPPFlo%row_zeroing_array = 0
  
  list => realization%output_option%output_snap_variable_list
  call WIPPFloSetPlotVariables(realization,list)
  list => realization%output_option%output_obs_variable_list
  call WIPPFloSetPlotVariables(realization,list)
  
  ! these counters are used for performance analysis only.  they will not 
  ! affect simulation results.
  liq_upwind_flip_count_by_res = 0
  gas_upwind_flip_count_by_res = 0
  liq_bc_upwind_flip_count_by_res = 0
  gas_bc_upwind_flip_count_by_res = 0
  liq_upwind_flip_count_by_jac = 0
  gas_upwind_flip_count_by_jac = 0
  liq_bc_upwind_flip_count_by_jac = 0
  gas_bc_upwind_flip_count_by_jac = 0

end subroutine WIPPFloSetup

! ************************************************************************** !

subroutine WIPPFloInitializeTimestep(realization)
  ! 
  ! Update data in module prior to time step
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Realization_Subsurface_class
  
  implicit none
  
  type(realization_subsurface_type) :: realization

  wippflo_update_upwind_direction = PETSC_TRUE
  call WIPPFloUpdateFixedAccum(realization)
  
end subroutine WIPPFloInitializeTimestep

! ************************************************************************** !

subroutine WIPPFloUpdateSolution(realization)
  ! 
  ! Updates data in module after a successful time
  ! step
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Realization_Subsurface_class
  
  implicit none
  
  type(realization_subsurface_type) :: realization

  if (realization%option%compute_mass_balance_new) then
    call WIPPFloUpdateMassBalance(realization)
  endif
  
end subroutine WIPPFloUpdateSolution

! ************************************************************************** !

subroutine WIPPFloTimeCut(realization)
  ! 
  ! Resets arrays for time step cut
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
  use Realization_Subsurface_class
 
  implicit none
  
  type(realization_subsurface_type) :: realization

  call WIPPFloInitializeTimestep(realization)  

end subroutine WIPPFloTimeCut

! ************************************************************************** !

subroutine WIPPFloComputeMassBalance(realization,mass_balance)
  ! 
  ! Initializes mass balance
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
  use Material_Aux_class
 
  implicit none
  
  type(realization_subsurface_type) :: realization
  PetscReal :: mass_balance(realization%option%nphase)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
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

  wippflo_auxvars => patch%aux%WIPPFlo%auxvars
  material_auxvars => patch%aux%Material%auxvars

  mass_balance = 0.d0

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    do iphase = 1, option%nphase
      ! volume_phase = saturation*porosity*volume
      vol_phase = &
        wippflo_auxvars(ZERO_INTEGER,ghosted_id)%sat(iphase)* &
        wippflo_auxvars(ZERO_INTEGER,ghosted_id)%effective_porosity* &
        material_auxvars(ghosted_id)%volume
      ! mass = volume_phase*density
      mass_balance(iphase) = mass_balance(iphase) + &
        wippflo_auxvars(ZERO_INTEGER,ghosted_id)%den(iphase)* &
        fmw_comp(iphase)*vol_phase
    enddo
  enddo

end subroutine WIPPFloComputeMassBalance

! ************************************************************************** !

subroutine WIPPFloZeroMassBalanceDelta(realization)
  ! 
  ! Zeros mass balance delta array
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
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

  do iconn = 1, patch%aux%WIPPFlo%num_aux_bc
    global_auxvars_bc(iconn)%mass_balance_delta = 0.d0
  enddo
  do iconn = 1, patch%aux%WIPPFlo%num_aux_ss
    global_auxvars_ss(iconn)%mass_balance_delta = 0.d0
  enddo

end subroutine WIPPFloZeroMassBalanceDelta

! ************************************************************************** !

subroutine WIPPFloUpdateMassBalance(realization)
  ! 
  ! Updates mass balance
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
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

  do iconn = 1, patch%aux%WIPPFlo%num_aux_bc
    do icomp = 1, option%nflowspec
      global_auxvars_bc(iconn)%mass_balance(icomp,:) = &
        global_auxvars_bc(iconn)%mass_balance(icomp,:) + &
        global_auxvars_bc(iconn)%mass_balance_delta(icomp,:)* &
        fmw_comp(icomp)*option%flow_dt
    enddo
  enddo
  do iconn = 1, patch%aux%WIPPFlo%num_aux_ss
    do icomp = 1, option%nflowspec
      global_auxvars_ss(iconn)%mass_balance(icomp,:) = &
        global_auxvars_ss(iconn)%mass_balance(icomp,:) + &
        global_auxvars_ss(iconn)%mass_balance_delta(icomp,:)* &
        fmw_comp(icomp)*option%flow_dt
    enddo
  enddo

end subroutine WIPPFloUpdateMassBalance

! ************************************************************************** !

subroutine WIPPFloUpdateAuxVars(realization)
  ! 
  ! Updates the auxiliary variables associated with the WIPPFlo problem
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
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
  use General_Aux_module, only : ANY_STATE, TWO_PHASE_STATE
  
  implicit none

  type(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars_bc(:)  
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)  
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn, natural_id
  PetscInt :: ghosted_start, ghosted_end
  PetscInt :: offset
  PetscInt :: istate
  PetscInt :: real_index, variable
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscErrorCode :: ierr
  
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  wippflo_auxvars => patch%aux%WIPPFlo%auxvars
  wippflo_auxvars_bc => patch%aux%WIPPFlo%auxvars_bc
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
    ! WIPPFLO_UPDATE_FOR_ACCUM indicates call from non-perturbation
    option%iflag = WIPPFLO_UPDATE_FOR_ACCUM
    natural_id = grid%nG2A(ghosted_id)
    if (grid%nG2L(ghosted_id) == 0) natural_id = -natural_id
    call WIPPFloAuxVarCompute(xx_loc_p(ghosted_start:ghosted_end), &
                       wippflo_auxvars(ZERO_INTEGER,ghosted_id), &
                       global_auxvars(ghosted_id), &
                       material_auxvars(ghosted_id), &
                       patch%characteristic_curves_array( &
                         patch%sat_func_id(ghosted_id))%ptr, &
                       natural_id, &
                       option)
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
      istate = boundary_condition%flow_aux_int_var(WIPPFLO_STATE_INDEX,iconn)
      if (istate == ANY_STATE) then
        do idof = 1, option%nflowdof
          select case(boundary_condition%flow_bc_type(idof))
            case(HYDROSTATIC_BC)
              real_index = boundary_condition% &
                             flow_aux_mapping(dof_to_primary_variable(idof))
              xxbc(idof) = boundary_condition% &
                             flow_aux_real_var(real_index,iconn)
            case(DIRICHLET_BC)
              variable = dof_to_primary_variable(idof)
              select case(variable)
                ! for liquid pressure dof
                case(WIPPFLO_LIQUID_PRESSURE_INDEX)
                  real_index = boundary_condition%flow_aux_mapping(variable)
                  if (real_index /= 0) then
                    xxbc(idof) = boundary_condition% &
                                   flow_aux_real_var(real_index,iconn)
                  else
                    option%io_buffer = 'Mixed FLOW_CONDITION "' // &
                      trim(boundary_condition%flow_condition%name) // &
                      '" needs gas pressure defined.'
                    call printErrMsg(option)
                  endif
                ! for gas saturation dof
                case(WIPPFLO_GAS_SATURATION_INDEX)
                  real_index = boundary_condition%flow_aux_mapping(variable)
                  if (real_index /= 0) then
                    xxbc(idof) = boundary_condition% &
                                   flow_aux_real_var(real_index,iconn)
                  endif
              end select
            case(NEUMANN_BC)
            case default
              option%io_buffer = 'Unknown BC type in WIPPFloUpdateAuxVars().'
              call printErrMsg(option)
          end select
        enddo  
      else
        ! we do this for all BCs; Neumann bcs will be set later
        do idof = 1, option%nflowdof
          real_index = boundary_condition% &
                         flow_aux_mapping(dof_to_primary_variable(idof))
          if (real_index > 0) then
            xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
          else
            option%io_buffer = 'Error setting up boundary condition in &
                               &WIPPFloUpdateAuxVars'
            call printErrMsg(option)
          endif
        enddo
      endif
          
      ! set this based on data given 
      global_auxvars_bc(sum_connection)%istate = TWO_PHASE_STATE
      ! WIPPFLO_UPDATE_FOR_BOUNDARY indicates call from non-perturbation
      option%iflag = WIPPFLO_UPDATE_FOR_BOUNDARY
      call WIPPFloAuxVarCompute(xxbc,wippflo_auxvars_bc(sum_connection), &
                                global_auxvars_bc(sum_connection), &
                                material_auxvars(ghosted_id), &
                                patch%characteristic_curves_array( &
                                  patch%sat_func_id(ghosted_id))%ptr, &
                                natural_id, &
                                option)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  patch%aux%WIPPFlo%auxvars_up_to_date = PETSC_TRUE

end subroutine WIPPFloUpdateAuxVars

! ************************************************************************** !

subroutine WIPPFloUpdateFixedAccum(realization)
  ! 
  ! Updates the fixed portion of the
  ! accumulation term
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Material_Aux_class

  implicit none
  
  type(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(material_parameter_type), pointer :: material_parameter

  PetscInt :: ghosted_id, local_id, local_start, local_end, natural_id
  PetscInt :: imat
  PetscReal, pointer :: xx_p(:), iphase_loc_p(:)
  PetscReal, pointer :: accum_p(:)
  PetscReal :: Jac_dummy(realization%option%nflowdof, &
                         realization%option%nflowdof)
                          
  PetscErrorCode :: ierr
  
  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  wippflo_auxvars => patch%aux%WIPPFlo%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  material_parameter => patch%aux%Material%material_parameter
    
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
    ! WIPPFLO_UPDATE_FOR_FIXED_ACCUM indicates call from non-perturbation
    option%iflag = WIPPFLO_UPDATE_FOR_FIXED_ACCUM
    call WIPPFloAuxVarCompute(xx_p(local_start:local_end), &
                              wippflo_auxvars(ZERO_INTEGER,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              patch%characteristic_curves_array( &
                                patch%sat_func_id(ghosted_id))%ptr, &
                              natural_id, &
                              option)
    call WIPPFloAccumulation(wippflo_auxvars(ZERO_INTEGER,ghosted_id), &
                             global_auxvars(ghosted_id), &
                             material_auxvars(ghosted_id), &
                             material_parameter%soil_heat_capacity(imat), &
                             option,accum_p(local_start:local_end), &
                             Jac_dummy,PETSC_FALSE, &
                             PETSC_FALSE)
  enddo
  
  call VecRestoreArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)
  
end subroutine WIPPFloUpdateFixedAccum

! ************************************************************************** !

subroutine WIPPFloResidual(snes,xx,r,realization,ierr)
  ! 
  ! Computes the residual equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
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

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_subsurface_type) :: realization
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  Mat, parameter :: null_mat = tMat(-1)
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(material_parameter_type), pointer :: material_parameter
  type(wippflo_parameter_type), pointer :: wippflo_parameter
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: iconn
  PetscReal :: scale
  PetscReal :: ss_flow_vol_flux(realization%option%nphase)
  PetscInt :: sum_connection
  PetscInt :: local_start, local_end
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
  PetscInt :: i, imat, imat_up, imat_dn
  PetscInt, pointer :: upwind_direction(:,:), upwind_direction_bc(:,:)

  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum_p(:), accum_p2(:)
  PetscReal, pointer :: vec_p(:)
  
  character(len=MAXSTRINGLENGTH) :: string

  PetscInt :: icap_up, icap_dn
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
  wippflo_auxvars => patch%aux%WIPPFlo%auxvars
  wippflo_auxvars_bc => patch%aux%WIPPFlo%auxvars_bc
  wippflo_parameter => patch%aux%WIPPFlo%wippflo_parameter
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars
  upwind_direction => patch%aux%WIPPFlo%upwind_direction
  upwind_direction_bc => patch%aux%WIPPFlo%upwind_direction_bc
  
  ! Communication -----------------------------------------
  ! must be called before WIPPFloUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  call WIPPFloUpdateAuxVars(realization)

  ! override flags since they will soon be out of date
  patch%aux%WIPPFlo%auxvars_up_to_date = PETSC_FALSE 

  if (option%compute_mass_balance_new) then
    call WIPPFloZeroMassBalanceDelta(realization)
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
  call VecGetArrayReadF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    local_end = local_id * option%nflowdof
    local_start = local_end - option%nflowdof + 1
    call WIPPFloAccumulation(wippflo_auxvars(ZERO_INTEGER,ghosted_id), &
                             global_auxvars(ghosted_id), &
                             material_auxvars(ghosted_id), &
                             material_parameter%soil_heat_capacity(imat), &
                             option,Res,Jac_dummy, &
                             wippflo_analytical_derivatives, &
                             PETSC_FALSE)
    r_p(local_start:local_end) =  r_p(local_start:local_end) + Res(:)
    accum_p2(local_start:local_end) = Res(:)
  enddo
  call VecRestoreArrayReadF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)

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

      icap_up = patch%sat_func_id(ghosted_id_up)
      icap_dn = patch%sat_func_id(ghosted_id_dn)

      call WIPPFloFlux(wippflo_auxvars(ZERO_INTEGER,ghosted_id_up), &
                       global_auxvars(ghosted_id_up), &
                       material_auxvars(ghosted_id_up), &
                       wippflo_auxvars(ZERO_INTEGER,ghosted_id_dn), &
                       global_auxvars(ghosted_id_dn), &
                       material_auxvars(ghosted_id_dn), &
                       cur_connection_set%area(iconn), &
                       cur_connection_set%dist(:,iconn), &
                       upwind_direction(:,sum_connection), &
                       wippflo_parameter,option,v_darcy,Res, &
                       Jac_dummy,Jac_dummy, &
                       PETSC_FALSE, & ! derivative call
                       wippflo_fix_upwind_direction, &
                       wippflo_update_upwind_direction, &
                       wippflo_count_upwind_dir_flip, & 
                       wippflo_analytical_derivatives, &
                       PETSC_FALSE)

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

      icap_dn = patch%sat_func_id(ghosted_id)

      call WIPPFloBCFlux(boundary_condition%flow_bc_type, &
                     boundary_condition%flow_aux_mapping, &
                     boundary_condition%flow_aux_real_var(:,iconn), &
                     wippflo_auxvars_bc(sum_connection), &
                     global_auxvars_bc(sum_connection), &
                     wippflo_auxvars(ZERO_INTEGER,ghosted_id), &
                     global_auxvars(ghosted_id), &
                     material_auxvars(ghosted_id), &
                     cur_connection_set%area(iconn), &
                     cur_connection_set%dist(:,iconn), &
                     upwind_direction_bc(:,sum_connection), &
                     wippflo_parameter,option, &
                     v_darcy,Res,Jac_dummy, &
                     PETSC_FALSE, & ! derivative call
                     wippflo_fix_upwind_direction, &
                     wippflo_update_upwind_direction, &
                     wippflo_count_upwind_dir_flip, & 
                     wippflo_analytical_derivatives, &
                     PETSC_FALSE)
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
      
      call WIPPFloSrcSink(option,source_sink%flow_condition%general%rate% &
                                  dataset%rarray(:), &
                          source_sink%flow_condition%general%rate%itype, &
                          wippflo_auxvars(ZERO_INTEGER,ghosted_id), &
                          global_auxvars(ghosted_id), &
                          ss_flow_vol_flux, &
                          scale,Res,Jac_dummy, &
                          wippflo_analytical_derivatives, &
                          PETSC_FALSE)

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

  if (patch%aux%WIPPFlo%inactive_cells_exist) then
    do i=1,patch%aux%WIPPFlo%n_inactive_rows
      r_p(patch%aux%WIPPFlo%inactive_rows_local(i)) = 0.d0
    enddo
  endif
  
  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
  
  call WIPPFloSSSandbox(r,null_mat,PETSC_FALSE,grid,material_auxvars, &
                        wippflo_auxvars,option)

  ! Mass Transfer
  if (field%flow_mass_transfer /= PETSC_NULL_VEC) then
    ! scale by -1.d0 for contribution to residual.  A negative contribution
    ! indicates mass being added to system.
    !call VecGetArrayF90(field%flow_mass_transfer,vec_p,ierr);CHKERRQ(ierr)
    !call VecRestoreArrayF90(field%flow_mass_transfer,vec_p,ierr);CHKERRQ(ierr)
    call VecAXPY(r,-1.d0,field%flow_mass_transfer,ierr);CHKERRQ(ierr)
  endif                      
                        
  if (realization%debug%vecview_residual) then
    string = 'WFresidual'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%vecview_solution) then
    string = 'WFxx'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  wippflo_update_upwind_direction = PETSC_FALSE

end subroutine WIPPFloResidual

! ************************************************************************** !

subroutine WIPPFloJacobian(snes,xx,A,B,realization,ierr)
  ! 
  ! Computes the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
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

  PetscInt :: icap_up,icap_dn
  PetscReal :: qsrc, scale
  PetscInt :: imat, imat_up, imat_dn
  PetscInt :: local_id, ghosted_id, natural_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  Vec, parameter :: null_vec = tVec(-1)
  
  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof), &
               Jdn(realization%option%nflowdof,realization%option%nflowdof)
  
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection  
  PetscInt, pointer :: zeros(:)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option 
  type(field_type), pointer :: field 
  type(material_parameter_type), pointer :: material_parameter
  type(wippflo_parameter_type), pointer :: wippflo_parameter
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:) 
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt, pointer :: upwind_direction(:,:), upwind_direction_bc(:,:)
  
  character(len=MAXSTRINGLENGTH) :: string
  
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  wippflo_auxvars => patch%aux%WIPPFlo%auxvars
  wippflo_auxvars_bc => patch%aux%WIPPFlo%auxvars_bc
  wippflo_parameter => patch%aux%WIPPFlo%wippflo_parameter
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars
  upwind_direction => patch%aux%WIPPFlo%upwind_direction
  upwind_direction_bc => patch%aux%WIPPFlo%upwind_direction_bc

  call MatGetType(A,mat_type,ierr);CHKERRQ(ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  else
    J = A
  endif

  call MatZeroEntries(J,ierr);CHKERRQ(ierr)

  if (.not.wippflo_analytical_derivatives) then
    ! Perturb aux vars
    do ghosted_id = 1, grid%ngmax  ! For each local node do...
      if (patch%imat(ghosted_id) <= 0) cycle
      natural_id = grid%nG2A(ghosted_id)
      call WIPPFloAuxVarPerturb(wippflo_auxvars(:,ghosted_id), &
                                global_auxvars(ghosted_id), &
                                material_auxvars(ghosted_id), &
                                patch%characteristic_curves_array( &
                                  patch%sat_func_id(ghosted_id))%ptr, &
                                natural_id,option)
    enddo
  endif
  
  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    call WIPPFloAccumDerivative(wippflo_auxvars(:,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              material_parameter%soil_heat_capacity(imat), &
                              option, &
                              Jup) 
    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                  ADD_VALUES,ierr);CHKERRQ(ierr)
  enddo

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
   
      icap_up = patch%sat_func_id(ghosted_id_up)
      icap_dn = patch%sat_func_id(ghosted_id_dn)
                              
      call WIPPFloFluxDerivative(wippflo_auxvars(:,ghosted_id_up), &
                     global_auxvars(ghosted_id_up), &
                     material_auxvars(ghosted_id_up), &
                     wippflo_auxvars(:,ghosted_id_dn), &
                     global_auxvars(ghosted_id_dn), &
                     material_auxvars(ghosted_id_dn), &
                     cur_connection_set%area(iconn), &
                     cur_connection_set%dist(:,iconn), &
                     upwind_direction(:,sum_connection), &
                     wippflo_parameter,option, &
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

      icap_dn = patch%sat_func_id(ghosted_id)

      call WIPPFloBCFluxDerivative(boundary_condition%flow_bc_type, &
                      boundary_condition%flow_aux_mapping, &
                      boundary_condition%flow_aux_real_var(:,iconn), &
                      wippflo_auxvars_bc(sum_connection), &
                      global_auxvars_bc(sum_connection), &
                      wippflo_auxvars(:,ghosted_id), &
                      global_auxvars(ghosted_id), &
                      material_auxvars(ghosted_id), &
                      cur_connection_set%area(iconn), &
                      cur_connection_set%dist(:,iconn), &
                      upwind_direction_bc(:,sum_connection), &
                      wippflo_parameter,option, &
                      Jdn)

      Jdn = -Jdn
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Source/sinks
  source_sink => patch%source_sink_list%first 
  do 
    if (.not.associated(source_sink)) exit
    
    cur_connection_set => source_sink%connection_set
    
    do iconn = 1, cur_connection_set%num_connections      
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      if (associated(source_sink%flow_aux_real_var)) then
        scale = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
      else
        scale = 1.d0
      endif
      
      Jup = 0.d0
      call WIPPFloSrcSinkDerivative(option, &
                        source_sink%flow_condition%general%rate% &
                                  dataset%rarray(:), &
                        source_sink%flow_condition%general%rate%itype, &
                        wippflo_auxvars(:,ghosted_id), &
                        global_auxvars(ghosted_id), &
                        scale,Jup)

      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)

    enddo
    source_sink => source_sink%next
  enddo
  
  call WIPPFloSSSandbox(null_vec,A,PETSC_TRUE,grid,material_auxvars, &
                        wippflo_auxvars,option)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  ! zero out inactive cells
  if (patch%aux%WIPPFlo%inactive_cells_exist) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    call MatZeroRowsLocal(A,patch%aux%WIPPFlo%n_inactive_rows, &
                          patch%aux%WIPPFlo%inactive_rows_local_ghosted, &
                          qsrc,PETSC_NULL_VEC,PETSC_NULL_VEC, &
                          ierr);CHKERRQ(ierr)
  endif

  if (realization%debug%matview_Jacobian) then
    string = 'WFjacobian'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(J,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (realization%debug%norm_Jacobian) then
    option => realization%option
    call MatNorm(J,NORM_1,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("1 norm: ",es11.4)') norm
    call printMsg(option) 
    call MatNorm(J,NORM_FROBENIUS,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("2 norm: ",es11.4)') norm
    call printMsg(option) 
    call MatNorm(J,NORM_INFINITY,norm,ierr);CHKERRQ(ierr)
    write(option%io_buffer,'("inf norm: ",es11.4)') norm
    call printMsg(option) 
  endif

end subroutine WIPPFloJacobian

! ************************************************************************** !

subroutine WIPPFloSetPlotVariables(realization,list)
  ! 
  ! Adds variables to be printed to list
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
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
  
end subroutine WIPPFloSetPlotVariables

! ************************************************************************** !

subroutine WIPPFloSSSandbox(residual,Jacobian,compute_derivative, &
                            grid,material_auxvars,wippflo_auxvars,option)
  ! 
  ! Evaluates source/sink term storing residual and/or Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
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
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  
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
    call WIPPFloSSSandboxLoadAuxReal(cur_srcsink,aux_real, &
                      wippflo_auxvars(ZERO_INTEGER,ghosted_id),option)
    call cur_srcsink%Evaluate(res,Jac,PETSC_FALSE, &
                              material_auxvars(ghosted_id), &
                              aux_real,option)
    if (compute_derivative) then
      do idof = 1, option%nflowdof
        res_pert = 0.d0
        call WIPPFloSSSandboxLoadAuxReal(cur_srcsink,aux_real, &
                                    wippflo_auxvars(idof,ghosted_id),option)
        call cur_srcsink%Evaluate(res_pert,Jac,PETSC_FALSE, &
                                  material_auxvars(ghosted_id), &
                                  aux_real,option)
        do irow = 1, option%nflowdof
          Jac(irow,idof) = (res_pert(irow)-res(irow)) / &
                            wippflo_auxvars(idof,ghosted_id)%pert
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

end subroutine WIPPFloSSSandbox

! ************************************************************************** !

subroutine WIPPFloSSSandboxLoadAuxReal(srcsink,aux_real,wippflo_auxvar,option)

  use Option_module
  use SrcSink_Sandbox_Base_class
  use SrcSink_Sandbox_WIPP_Gas_class
  use SrcSink_Sandbox_WIPP_Well_class

  implicit none

  class(srcsink_sandbox_base_type) :: srcsink
  PetscReal :: aux_real(:)
  type(wippflo_auxvar_type) wippflo_auxvar
  type(option_type) :: option
  
  aux_real = 0.d0
  select type(srcsink)
    class is(srcsink_sandbox_wipp_gas_type)
      aux_real(WIPP_GAS_WATER_SATURATION_INDEX) = &
        wippflo_auxvar%sat(option%liquid_phase)
    class is(srcsink_sandbox_wipp_well_type)
      aux_real(WIPP_WELL_LIQUID_MOBILITY) = &
        wippflo_auxvar%mobility(option%liquid_phase)
      aux_real(WIPP_WELL_GAS_MOBILITY) = &
        wippflo_auxvar%mobility(option%gas_phase)
      aux_real(WIPP_WELL_LIQUID_PRESSURE) = &
        wippflo_auxvar%pres(option%liquid_phase)
      aux_real(WIPP_WELL_GAS_PRESSURE) = &
        wippflo_auxvar%pres(option%gas_phase)
      aux_real(WIPP_WELL_LIQUID_DENSITY) = &
        wippflo_auxvar%den(option%liquid_phase)
      aux_real(WIPP_WELL_GAS_DENSITY) = &
        wippflo_auxvar%den(option%gas_phase)
  end select
  
end subroutine WIPPFloSSSandboxLoadAuxReal

! ************************************************************************** !

subroutine WIPPFloMapBCAuxVarsToGlobal(realization)
  ! 
  ! Maps variables in general auxvar to global equivalent.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
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
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars_bc(:)  
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)  

  PetscInt :: sum_connection, iconn
  
  option => realization%option
  patch => realization%patch

  if (option%ntrandof == 0) return ! no need to update
  
  wippflo_auxvars_bc => patch%aux%WIPPFlo%auxvars_bc
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      global_auxvars_bc(sum_connection)%sat = &
        wippflo_auxvars_bc(sum_connection)%sat
      global_auxvars_bc(sum_connection)%den_kg = &
        wippflo_auxvars_bc(sum_connection)%den_kg
    enddo
    boundary_condition => boundary_condition%next
  enddo
  
end subroutine WIPPFloMapBCAuxVarsToGlobal

! ************************************************************************** !

subroutine WIPPFloDestroy(realization)
  ! 
  ! Deallocates variables associated with Richard
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Realization_Subsurface_class
  use Option_module

  implicit none

  type(realization_subsurface_type) :: realization
  
  ! place anything that needs to be freed here.
  ! auxvars are deallocated in auxiliary.F90.

  if (wippflo_fix_upwind_direction .and. &
      wippflo_count_upwind_dir_flip .and. &
      OptionPrintToScreen(realization%option)) then
    print *, 'Res upwind dir. flip (liq): ', liq_upwind_flip_count_by_res
    print *, 'Res upwind dir. flip (gas): ', gas_upwind_flip_count_by_res
    print *, 'Res upwind dir. flip (liq bc): ', liq_bc_upwind_flip_count_by_res
    print *, 'Res upwind dir. flip (gas bc): ', gas_bc_upwind_flip_count_by_res
    print *, 'Jac upwind dir. flip (liq): ', liq_upwind_flip_count_by_jac
    print *, 'Jac upwind dir. flip (gas): ', gas_upwind_flip_count_by_jac
    print *, 'Jac upwind dir. flip (liq bc): ', liq_bc_upwind_flip_count_by_jac
    print *, 'Jac upwind dir. flip (gas bc): ', gas_bc_upwind_flip_count_by_jac
  endif

end subroutine WIPPFloDestroy

end module WIPP_Flow_module
