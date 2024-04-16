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
            WIPPFloZeroMassBalanceDelta, &
            WIPPFloResidual, &
            WIPPFloJacobian, &
            WIPPFloSetPlotVariables, &
            WIPPFloMapBCAuxVarsToGlobal, &
            WIPPFloCreepShutDown, &
            WIPPFloSSSandbox, &
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
  use Material_Aux_module
  use Output_Aux_module
  use Characteristic_Curves_module
  use WIPP_Characteristic_Curve_module
  use Matrix_Zeroing_module
!  use PM_Base_class
!  use PM_WIPP_Flow_class

  implicit none

  class(realization_subsurface_type) :: realization
!  class(pm_base_type), pointer :: pm

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(output_variable_list_type), pointer :: list
  type(material_parameter_type), pointer :: material_parameter

  PetscInt :: ghosted_id, iconn, sum_connection, local_id
  PetscInt :: i, idof, ndof
  PetscBool :: error_found
  PetscInt :: flag(10)
  PetscErrorCode :: ierr
                                                ! extra index for derivatives
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars_bc(:)
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars_ss(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  class(characteristic_curves_type), pointer :: cc

  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%aux%WIPPFlo => WIPPFloAuxCreate(option)

  ! ensure that material properties specific to this module are properly
  ! initialized
  material_parameter => patch%aux%Material%material_parameter
  error_found = PETSC_FALSE

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
    if (material_auxvars(ghosted_id)%porosity_base < 0.d0 .and. &
        flag(2) == 0) then
      flag(2) = 1
      option%io_buffer = 'ERROR: Non-initialized porosity.'
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
                     MPI_LOR,option%mycomm,ierr);CHKERRQ(ierr)
  if (error_found) then
    option%io_buffer = 'Material property errors found in WIPPFloSetup.'
    call PrintErrMsg(option)
  endif

  ndof = option%nflowdof
  allocate(wippflo_auxvars(0:ndof,grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    do idof = 0, ndof
      call WIPPFloAuxVarInit(wippflo_auxvars(idof,ghosted_id),option)
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
      call WIPPFloAuxVarInit(wippflo_auxvars_bc(iconn),option)
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
      call WIPPFloAuxVarInit(wippflo_auxvars_ss(iconn),option)
    enddo
    patch%aux%WIPPFlo%auxvars_ss => wippflo_auxvars_ss
  endif
  patch%aux%WIPPFlo%num_aux_ss = sum_connection

  list => realization%output_option%output_snap_variable_list
  call WIPPFloSetPlotVariables(realization,list)
  list => realization%output_option%output_obs_variable_list
  call WIPPFloSetPlotVariables(realization,list)

  if (wippflo_use_lumped_harm_flux) then
    XXFlux => WIPPFloFluxLumpedHarmonic
    XXBCFlux => WIPPFloBCFluxLumpedHarmonic
  else
    XXFlux => WIPPFloFluxHarmonicPermOnly
    XXBCFlux => WIPPFloBCFluxHarmonicPermOnly
  endif

  if (wippflo_use_bragflo_cc) then
    do i = 1, size(realization%patch%characteristic_curves_array)
      cc => realization%patch%characteristic_curves_array(i)%ptr
      call WIPPCCVerify(cc%saturation_function, &
                        cc%liq_rel_perm_function, &
                        cc%gas_rel_perm_function, &
                        option)
    enddo
  endif

  wippflo_ts_count = 0
  wippflo_ts_cut_count = 0
  wippflo_ni_count = 0

  call PatchSetupUpwindDirection(patch,option)

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
  use Upwind_Direction_module

  implicit none

  class(realization_subsurface_type) :: realization

  wippflo_newton_iteration_number = 0
  update_upwind_direction = PETSC_TRUE
  call WIPPFloUpdateFixedAccum(realization)

  wippflo_ni_count = 0
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

  class(realization_subsurface_type) :: realization

  if (realization%option%compute_mass_balance_new) then
    call WIPPFloUpdateMassBalance(realization)
  endif

  call WIPPFloCreepShutDown(realization)

  wippflo_ts_count = wippflo_ts_count + 1
  wippflo_ts_cut_count = 0
  wippflo_ni_count = 0

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

  class(realization_subsurface_type) :: realization

  wippflo_ts_cut_count = wippflo_ts_cut_count + 1

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
  use Material_Aux_module

  implicit none

  class(realization_subsurface_type) :: realization
  PetscReal :: mass_balance(realization%option%nphase)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iphase
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

  class(realization_subsurface_type) :: realization

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

  class(realization_subsurface_type) :: realization

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
  use Material_Aux_module
  use General_Aux_module, only : ANY_STATE, TWO_PHASE_STATE

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

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

  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

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
                         patch%cc_id(ghosted_id))%ptr, &
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
                    call PrintErrMsg(option)
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
              call PrintErrMsg(option)
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
            call PrintErrMsg(option)
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
                                  patch%cc_id(ghosted_id))%ptr, &
                                natural_id, &
                                option)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

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
  use Material_Aux_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(material_parameter_type), pointer :: material_parameter

  PetscInt :: ghosted_id, local_id, local_start, local_end, natural_id
  PetscInt :: imat
  PetscReal, pointer :: xx_p(:)
  PetscReal, pointer :: accum_p(:)

  PetscErrorCode :: ierr

  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  wippflo_auxvars => patch%aux%WIPPFlo%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  material_parameter => patch%aux%Material%material_parameter

  call VecGetArrayReadF90(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)

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
                                patch%cc_id(ghosted_id))%ptr, &
                              natural_id, &
                              option)
    call WIPPFloAccumulation(wippflo_auxvars(ZERO_INTEGER,ghosted_id), &
                             global_auxvars(ghosted_id), &
                             material_auxvars(ghosted_id), &
                             option,accum_p(local_start:local_end), &
                             PETSC_FALSE)
    call WIPPFloConvertUnitsToBRAGFlo(accum_p(local_start:local_end), &
                                      material_auxvars(ghosted_id),option)
  enddo

  call VecRestoreArrayReadF90(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)

end subroutine WIPPFloUpdateFixedAccum

! ************************************************************************** !

subroutine WIPPFloNumericalJacobianTest(xx,B,realization,pmwss_ptr,pmwell_ptr)
  !
  ! Computes the a test numerical jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 03/03/15
  !

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Grid_module
  use Field_module
  use PM_WIPP_SrcSink_class
  use PM_Well_class

  implicit none

  Vec :: xx
  class(realization_subsurface_type) :: realization
  Mat :: B
  class(pm_wipp_srcsink_type), pointer :: pmwss_ptr
  class(pm_well_type), pointer :: pmwell_ptr

  Vec :: xx_pert
  Vec :: res
  Vec :: res_pert
  Mat :: A
  PetscErrorCode :: ierr

  PetscReal, pointer :: vec_p(:), vec2_p(:)

  PetscViewer :: viewer
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  PetscReal :: derivative, perturbation
  PetscInt, save :: icall = 0
  character(len=MAXWORDLENGTH) :: word

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
                   grid%nlmax*option%nflowdof,ierr);CHKERRQ(ierr)
  call MatSeqAIJSetPreallocation(A,27,PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
  call MatSetFromOptions(A,ierr);CHKERRQ(ierr)
  call MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE, &
                    ierr);CHKERRQ(ierr)

  call VecZeroEntries(res,ierr);CHKERRQ(ierr)
  print *, 'Unperturbed Residual'
  call WIPPFloResidual(PETSC_NULL_SNES,xx,res,realization,pmwss_ptr, &
                       pmwell_ptr,ierr)
#if 0
  word  = 'num_0.dat'
  call PetscViewerASCIIOpen(option%mycomm,word,viewer,ierr);CHKERRQ(ierr)
  call VecView(res,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
  call VecGetArrayF90(res,vec2_p,ierr);CHKERRQ(ierr)
  print *, 'Perturbed Residual'
  do icell = 1,grid%nlmax
    if (patch%imat(grid%nL2G(icell)) <= 0) cycle
    do idof = (icell-1)*option%nflowdof+1,icell*option%nflowdof

      call VecCopy(xx,xx_pert,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(xx_pert,vec_p,ierr);CHKERRQ(ierr)
      if (mod(idof,2) == 1) then ! liquid pressure
        perturbation = max(wippflo_pres_rel_pert * vec_p(idof), &
                           wippflo_pres_min_pert)
      else ! gas saturation
        perturbation = max(wippflo_sat_rel_pert * vec_p(idof), &
                           wippflo_sat_min_pert)
        if (1.d0 - vec_p(idof) - &
            patch%characteristic_curves_array( &
              patch%cc_id(icell))%ptr% &
              saturation_function%Sr < 0.d0) then
          if (vec_p(idof) + perturbation > 1.d0) then
            perturbation = -1.d0 * perturbation
          endif
        else
          perturbation = -1.d0 * perturbation
          if (vec_p(idof) + perturbation < 0.d0) then
            perturbation = -1.d0 * perturbation
          endif
        endif
      endif
      vec_p(idof) = vec_p(idof) + perturbation
      call VecRestoreArrayF90(xx_pert,vec_p,ierr);CHKERRQ(ierr)

      call VecZeroEntries(res_pert,ierr);CHKERRQ(ierr)
      write(*,'("icell: ",i4," idof: ",i2," overall dof: ",i4,&
                &" pert: ",1pe12.5)') &
        icell,2-mod(idof,option%nflowdof),idof,perturbation
      call WIPPFloResidual(PETSC_NULL_SNES,xx_pert,res_pert,realization, &
                           pmwss_ptr,pmwell_ptr,ierr)
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
        if (dabs(derivative) > 1.d-40) then
          if (idof2 == wippflo_jacobian_test_rdof) then
            print *, 'dof: ', idof2, idof, '| cell: ', (idof2+1)/2, (idof+1)/2
            print *, vec_p(idof2), vec2_p(idof2), perturbation
          endif
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

end subroutine WIPPFloNumericalJacobianTest

! ************************************************************************** !

subroutine WIPPFloResidual(snes,xx,r,realization,pmwss_ptr,pmwell_ptr,ierr)
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
  use Debug_module

  use Connection_module
  use Grid_module
  use Coupler_module
  use Material_Aux_module
  use PM_WIPP_SrcSink_class
  use PM_Well_class
  use Upwind_Direction_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  class(realization_subsurface_type) :: realization
  class(pm_wipp_srcsink_type), pointer :: pmwss_ptr
  class(pm_well_type), pointer :: pmwell_ptr
  PetscErrorCode :: ierr
  PetscViewer :: viewer

  Mat, parameter :: null_mat = tMat(0)
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
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt, pointer :: upwind_direction(:,:)
  PetscInt, pointer :: upwind_direction_bc(:,:)

  PetscInt :: iconn
  PetscReal :: scale
  PetscReal :: ss_flow_vol_flux(realization%option%nphase)
  PetscInt :: sum_connection
  PetscInt :: local_start, local_end
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
  PetscInt :: i, imat, imat_up, imat_dn

  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum_p(:), accum_p2(:)
  PetscReal, pointer :: xx_p(:), scaled_xx_p(:)
  PetscReal, pointer :: vec_p(:)
  PetscBool :: debug_connection

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: irow

  PetscInt :: icc_up, icc_dn
  PetscReal :: Res(realization%option%nflowdof)
  PetscReal :: temp_Res(realization%option%nflowdof)
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
  upwind_direction => patch%flow_upwind_direction
  upwind_direction_bc => patch%flow_upwind_direction_bc

  if (.not.wippflo_jacobian_test_active) then
    wippflo_newton_iteration_number = wippflo_newton_iteration_number + 1
    ! bragflo uses the following logic, update when
    !   it == 1, before entering iteration loop
    !   it > 1 and mod(it-1,frequency) == 0
    ! the first is set in WIPPFloInitializeTimestep, the second is set here
    if (wippflo_newton_iteration_number > 1 .and. &
        mod(wippflo_newton_iteration_number-1, &
            upwind_dir_update_freq) == 0) then
      update_upwind_direction = PETSC_TRUE
    endif
  endif

  ! Communication -----------------------------------------
  ! must be called before WIPPFloUpdateAuxVars()
  if (option%flow%scale_all_pressure) then
    ! have to convert the log concentration to non-log form
    call VecGetArrayF90(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(xx,scaled_xx_p,ierr);CHKERRQ(ierr)
    do irow = 1, size(xx_p), 2
      xx_p(irow) = scaled_xx_p(irow) * option%flow%pressure_scaling_factor
      xx_p(irow+1) = scaled_xx_p(irow+1)
    enddo
    call VecRestoreArrayF90(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(xx,scaled_xx_p,ierr);CHKERRQ(ierr)
    call DiscretizationGlobalToLocal(discretization,field%flow_xx, &
                                     field%flow_xx_loc,NFLOWDOF)
  else
    call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,&
                                     NFLOWDOF)
  endif

  call WIPPFloUpdateAuxVars(realization)

  ! override flags since they will soon be out of date
  patch%aux%WIPPFlo%auxvars_up_to_date = PETSC_FALSE

  if (option%compute_mass_balance_new) then
    call WIPPFloZeroMassBalanceDelta(realization)
  endif

  option%iflag = 1
  ! now assign access pointer to local variables
  call VecGetArrayF90(r,r_p,ierr);CHKERRQ(ierr)

  ! Accumulation terms ------------------------------------
  ! accumulation at t(k) (doesn't change during Newton iteration)
  if (wippflo_calc_accum) then
  call VecGetArrayReadF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)
  r_p = -accum_p
  if (wippflo_residual_test .and. &
      wippflo_residual_test_cell  > 0) then
    local_end = wippflo_residual_test_cell * option%nflowdof
    local_start = local_end - option%nflowdof + 1
    write(*,'(" Aold: ",2es12.4,i5)') &
      -1.d0*r_p(local_start:local_end)/option%flow_dt, &
      wippflo_residual_test_cell
  endif
  call VecRestoreArrayReadF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)

  ! accumulation at t(k+1)
  call VecGetArrayF90(field%flow_accum2,accum_p2,ierr);CHKERRQ(ierr)
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
                             option,Res,PETSC_FALSE)
    call WIPPFloConvertUnitsToBRAGFlo(Res,material_auxvars(ghosted_id),option)
    r_p(local_start:local_end) =  r_p(local_start:local_end) + Res(:)
    accum_p2(local_start:local_end) = Res(:)
    if (wippflo_residual_test .and. &
        wippflo_residual_test_cell == local_id) then
      write(*,'(" DT[y]: ",es12.4)') option%flow_dt/3600.d0/24.d0/365.d0
      write(*,'(" A(calc): ",i5,8es12.4)') &
        wippflo_residual_test_cell, &
        wippflo_auxvars(ZERO_INTEGER,ghosted_id)%effective_porosity, &
        wippflo_auxvars(ZERO_INTEGER,ghosted_id)%den_kg(1), &
        wippflo_auxvars(ZERO_INTEGER,ghosted_id)%sat(1), &
        wippflo_auxvars(ZERO_INTEGER,ghosted_id)%den_kg(2), &
        wippflo_auxvars(ZERO_INTEGER,ghosted_id)%sat(2)
      write(*,'(" A: ",2es12.4,i5)') &
        Res(:)/option%flow_dt, wippflo_residual_test_cell
    endif
  enddo
  call VecRestoreArrayF90(field%flow_accum2,accum_p2,ierr);CHKERRQ(ierr)
  else
    r_p = 0.d0
  endif

  if (wippflo_calc_flux) then
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

      if (wippflo_jacobian_test) then
        !if (wippflo_jacobian_test_rdof > 0 .and. &
        !    .not.(ghosted_id_up == &
        !          int((wippflo_jacobian_test_rdof+1)/2) .and. &
        !          ghosted_id_dn == &
        !          int((wippflo_jacobian_test_rdof+1)/2))) cycle
        if (wippflo_jacobian_test_rdof > 0 .and. &
            .not.(ghosted_id_up == &
                  int((wippflo_jacobian_test_rdof+1)/2) .or. &
                  ghosted_id_dn == &
                  int((wippflo_jacobian_test_rdof+1)/2))) cycle
      endif
      debug_connection = PETSC_FALSE
      if (wippflo_residual_test .and. &
          (wippflo_residual_test_cell == local_id_up .or. &
           wippflo_residual_test_cell == local_id_dn)) then
        debug_connection = PETSC_TRUE
      endif
      if (wippflo_print_oscillatory_behavior) then
        if (((wippflo_prev_liq_res_cell(1) == local_id_up .and. &
              wippflo_prev_liq_res_cell(2) == local_id_dn) .or. &
             (wippflo_prev_liq_res_cell(2) == local_id_up .and. &
              wippflo_prev_liq_res_cell(1) == local_id_dn)) .or. &
            (wippflo_prev_liq_res_cell(1) == &
             wippflo_prev_liq_res_cell(2) .and. &
             (wippflo_prev_liq_res_cell(1) == local_id_up .or. &
              wippflo_prev_liq_res_cell(1) == local_id_dn))) then
          write(*,'("debug flux: ",2i5)') local_id_up, local_id_dn
          debug_connection = PETSC_TRUE
        endif
      endif
      call XXFlux(wippflo_auxvars(ZERO_INTEGER,ghosted_id_up), &
                       global_auxvars(ghosted_id_up), &
                       material_auxvars(ghosted_id_up), &
                       wippflo_auxvars(ZERO_INTEGER,ghosted_id_dn), &
                       global_auxvars(ghosted_id_dn), &
                       material_auxvars(ghosted_id_dn), &
                       cur_connection_set%area(iconn), &
                       cur_connection_set%dist(:,iconn), &
                       upwind_direction(:,sum_connection), &
                       wippflo_parameter,option,v_darcy,Res, &
                       PETSC_FALSE, & ! derivative call
                       update_upwind_direction, &
                       count_upwind_direction_flip, &
                       debug_connection)

      patch%internal_velocities(:,sum_connection) = v_darcy
      if (associated(patch%internal_flow_fluxes)) then
        patch%internal_flow_fluxes(:,sum_connection) = Res(:)
      endif

      if (local_id_up > 0) then
        local_end = local_id_up * option%nflowdof
        local_start = local_end - option%nflowdof + 1
        temp_Res = Res
        call WIPPFloConvertUnitsToBRAGFlo(temp_Res, &
                                         material_auxvars(ghosted_id_up), &
                                         option)
        r_p(local_start:local_end) = r_p(local_start:local_end) + temp_Res(:)
        if ((wippflo_jacobian_test .and. wippflo_jacobian_test_rdof > 0) .or. &
            (wippflo_residual_test .and. &
             wippflo_residual_test_cell  == local_id_up)) then
          write(*,'(" Fup: ",2es12.4,2i5)') -1.d0*temp_Res/option%flow_dt, &
            local_id_up, local_id_dn
!          write(*,'("      ",8es12.4)') cur_connection_set%area(iconn), &
!            2.d0*(cur_connection_set%dist(0,iconn)* &
!                  cur_connection_set%dist(-1,iconn)), &
!            wippflo_auxvars(ZERO_INTEGER,ghosted_id_up)%alpha, &
!            material_auxvars(ghosted_id_up)%volume
        endif
      endif

      if (local_id_dn > 0) then
        local_end = local_id_dn * option%nflowdof
        local_start = local_end - option%nflowdof + 1
        temp_Res = Res
        call WIPPFloConvertUnitsToBRAGFlo(temp_Res, &
                                         material_auxvars(ghosted_id_dn), &
                                         option)
        r_p(local_start:local_end) = r_p(local_start:local_end) - temp_Res(:)
        if ((wippflo_jacobian_test .and. wippflo_jacobian_test_rdof > 0) .or. &
            (wippflo_residual_test .and. &
             wippflo_residual_test_cell  == local_id_dn)) then
          write(*,'(" Fdn: ",2es12.4,2i5)') temp_Res/option%flow_dt, &
            local_id_up, local_id_dn
!          write(*,'("      ",8es12.4)') cur_connection_set%area(iconn), &
!            2.d0*(cur_connection_set%dist(0,iconn)* &
!                  (1.d0-cur_connection_set%dist(-1,iconn))), &
!            wippflo_auxvars(ZERO_INTEGER,ghosted_id_dn)%alpha, &
!            material_auxvars(ghosted_id_dn)%volume
        endif
      endif
    enddo

    cur_connection_set => cur_connection_set%next
  enddo
  endif

  if (wippflo_calc_bcflux) then
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

      icc_dn = patch%cc_id(ghosted_id)

      call XXBCFlux(boundary_condition%flow_bc_type, &
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
                     v_darcy,Res, &
                     PETSC_FALSE, & ! derivative call
                     update_upwind_direction, &
                     count_upwind_direction_flip, &
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
      call WIPPFloConvertUnitsToBRAGFlo(Res, &
                                       material_auxvars(ghosted_id), &
                                       option)
      r_p(local_start:local_end)= r_p(local_start:local_end) - Res(:)
      if ((wippflo_jacobian_test .and. wippflo_jacobian_test_rdof > 0) .or. &
          (wippflo_residual_test .and. &
           wippflo_residual_test_cell  == local_id)) then
        write(*,'(" BCF: ",2es12.4,i5 )') Res/option%flow_dt, local_id
      endif

    enddo
    boundary_condition => boundary_condition%next
  enddo
  endif

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

      if (Initialized(wippflo_auxvars(ZERO_INTEGER,ghosted_id)%well%pl)) then
        if (dabs(wippflo_auxvars(ZERO_INTEGER,ghosted_id)%well%dpl) < &
            1.d-15) then
          scale = 0.d0
        else
          ! jmfrede 09/14/2022 Getting rid of this scale factor because it
          ! changes wildly within Newton iterations, which seems to make
          ! WIPP_FLOW have a harder time converging. When it does converge,
          ! the scale seems to be ~ 1 anyways. Uncomment the WRITE statement
          ! to quickly see the value of scale printed to screen.
          scale = dabs(wippflo_auxvars(ZERO_INTEGER,ghosted_id)% &
                  pres(ONE_INTEGER)-wippflo_auxvars(ZERO_INTEGER,ghosted_id)% &
                  well%pl)/dabs(wippflo_auxvars(ZERO_INTEGER,ghosted_id)% &
                  well%dpl)
          scale = 1.d0
        endif
        !WRITE(*,*) 'SCALE = ', scale
      endif

      call WIPPFloSrcSink(option,source_sink%flow_condition%general%rate% &
                                  dataset%rarray(:), &
                          source_sink%flow_condition%general%rate%itype, &
                          wippflo_auxvars(ZERO_INTEGER,ghosted_id), &
                          global_auxvars(ghosted_id), &
                          material_auxvars(ghosted_id), &
                          ss_flow_vol_flux, &
                          scale,Res, &
                          PETSC_FALSE)

      if (associated(patch%ss_flow_vol_fluxes)) then
        patch%ss_flow_vol_fluxes(:,sum_connection) = ss_flow_vol_flux
      endif
      if (associated(patch%ss_flow_fluxes)) then
        patch%ss_flow_fluxes(:,sum_connection) = Res(:)
      endif
      if (option%compute_mass_balance_new) then
        if (associated(global_auxvars_ss)) then
        ! contribution to boundary
          global_auxvars_ss(sum_connection)%mass_balance_delta(1:2,1) = &
            global_auxvars_ss(sum_connection)%mass_balance_delta(1:2,1) - &
            Res(1:2)
        endif
      endif
      call WIPPFloConvertUnitsToBRAGFlo(Res, &
                                       material_auxvars(ghosted_id), &
                                       option)
      r_p(local_start:local_end) =  r_p(local_start:local_end) - Res(:)

    enddo
    source_sink => source_sink%next
  enddo

  if (wippflo_calc_chem) then
  ! WIPP gas/brine generation process model source/sinks
  if (associated(pmwss_ptr)) then
    if (pmwss_ptr%rate_update_frequency == LAG_NEWTON_ITERATION .or. &
        pmwss_ptr%rate_update_frequency == NO_LAG) then
      call PMWSSUpdateRates(pmwss_ptr,PETSC_FALSE,ierr)
    endif
    call PMWSSCalcResidualValues(pmwss_ptr,r_p,ss_flow_vol_flux)
  endif
  endif

  ! Compute WIPP well model source/sinks for the quasi-implicitly coupled well
  ! model approach
  if (wippflo_well_quasi_imp_coupled) then
    if (associated(pmwell_ptr)) then
      if (any(pmwell_ptr%well_grid%h_rank_id == option%myrank)) then
        call pmwell_ptr%UpdateFlowRates(ZERO_INTEGER,ZERO_INTEGER,-999,ierr)
        if (pmwell_ptr%well_force_ts_cut == ZERO_INTEGER) then
          call pmwell_ptr%ModifyFlowResidual(r_p,ss_flow_vol_flux)
        endif
      endif
    endif
  endif

  if (patch%aux%WIPPFlo%inactive_cells_exist) then
    do i=1,patch%aux%WIPPFlo%matrix_zeroing%n_inactive_rows
      r_p(patch%aux%WIPPFlo%matrix_zeroing%inactive_rows_local(i)) = 0.d0
    enddo
  endif

  call VecRestoreArrayF90(r,r_p,ierr);CHKERRQ(ierr)

  call WIPPFloSSSandbox(r,null_mat,PETSC_FALSE,grid,material_auxvars, &
                        wippflo_auxvars,option)

  ! Mass Transfer
  if (field%flow_mass_transfer /= PETSC_NULL_VEC) then
    ! scale by -1.d0 for contribution to residual.  A negative contribution
    ! indicates mass being added to system.
    call VecGetArrayF90(r,r_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%flow_mass_transfer,vec_p,ierr);CHKERRQ(ierr)
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      imat = patch%imat(ghosted_id)
      if (imat <= 0) cycle
      local_end = local_id * option%nflowdof
      local_start = local_end - option%nflowdof + 1
      Res = vec_p(local_start:local_end)
      call WIPPFloConvertUnitsToBRAGFlo(Res,material_auxvars(ghosted_id),option)
      r_p(local_start:local_end) = r_p(local_start:local_end) - Res
    enddo
    call VecRestoreArrayF90(r,r_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%flow_mass_transfer,vec_p, &
                            ierr);CHKERRQ(ierr)
!geh: due to the potential for units conversion, cannot VecAXPY
!    call VecAXPY(r,-1.d0,field%flow_mass_transfer,ierr);CHKERRQ(ierr)
  endif

  update_upwind_direction = PETSC_FALSE

  if (realization%debug%vecview_residual) then
    call DebugWriteFilename(realization%debug,string,'WFresidual','', &
                            wippflo_ts_count,wippflo_ts_cut_count, &
                            wippflo_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif
  if (realization%debug%vecview_solution) then
    call DebugWriteFilename(realization%debug,string,'WFxx','', &
                            wippflo_ts_count,wippflo_ts_cut_count, &
                            wippflo_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif

end subroutine WIPPFloResidual

! ************************************************************************** !

subroutine WIPPFloJacobian(snes,xx,A,B,realization,pmwss_ptr,pmwell_ptr,ierr)
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
  use Material_Aux_module
  use PM_WIPP_SrcSink_class
  use PM_Well_class
  use Upwind_Direction_module
  use Debug_module

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  class(realization_subsurface_type) :: realization
  class(pm_wipp_srcsink_type), pointer :: pmwss_ptr
  class(pm_well_type), pointer :: pmwell_ptr
  PetscErrorCode :: ierr
  PetscViewer :: viewer

  Mat :: J
  MatType :: mat_type
  PetscReal :: norm

  PetscInt :: icc_up,icc_dn
  PetscReal :: qsrc, scale
  PetscInt :: imat, imat_up, imat_dn
  PetscInt :: local_id, ghosted_id, natural_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  Vec, parameter :: null_vec = tVec(0)

  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof), &
               Jdn(realization%option%nflowdof,realization%option%nflowdof)
  PetscReal :: Jtmp(realization%option%nflowdof,realization%option%nflowdof)

  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(material_parameter_type), pointer :: material_parameter
  type(wippflo_parameter_type), pointer :: wippflo_parameter
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  PetscInt, pointer :: upwind_direction(:,:)
  PetscInt, pointer :: upwind_direction_bc(:,:)

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
  upwind_direction => patch%flow_upwind_direction
  upwind_direction_bc => patch%flow_upwind_direction_bc

  call MatGetType(A,mat_type,ierr);CHKERRQ(ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  else
    J = A
  endif

  call MatZeroEntries(J,ierr);CHKERRQ(ierr)

  ! Perturb aux vars
  do ghosted_id = 1, grid%ngmax  ! For each local node do...
    if (patch%imat(ghosted_id) <= 0) cycle
    natural_id = grid%nG2A(ghosted_id)
    call WIPPFloAuxVarPerturb(wippflo_auxvars(:,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              patch%characteristic_curves_array( &
                                patch%cc_id(ghosted_id))%ptr, &
                              natural_id,option)
  enddo

  if (wippflo_calc_accum) then
  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    call WIPPFloAccumDerivative(wippflo_auxvars(:,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              option, &
                              Jup)
    call WIPPFloConvertUnitsToBRAGFlo(Jup,material_auxvars(ghosted_id), &
                                      option)
    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                  ADD_VALUES,ierr);CHKERRQ(ierr)
  enddo

  if (realization%debug%matview_Matrix_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call DebugWriteFilename(realization%debug,string,'WFjacobian_accum','', &
                            wippflo_ts_count,wippflo_ts_cut_count, &
                            wippflo_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif

  endif

  if (wippflo_calc_flux) then
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

      if (wippflo_jacobian_test) then
        !if (wippflo_jacobian_test_rdof > 0 .and. &
        !   .not.(ghosted_id_up == &
        !          int((wippflo_jacobian_test_rdof+1)/2) .and. &
        !          ghosted_id_dn == &
        !          int((wippflo_jacobian_test_rdof+1)/2))) cycle
        if (wippflo_jacobian_test_rdof > 0 .and. &
            .not.(ghosted_id_up == &
                  int((wippflo_jacobian_test_rdof+1)/2) .or. &
                  ghosted_id_dn == &
                  int((wippflo_jacobian_test_rdof+1)/2))) cycle
      endif
      call XXFluxDerivative(wippflo_auxvars(:,ghosted_id_up), &
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
        Jtmp = Jup
        call WIPPFloConvertUnitsToBRAGFlo(Jtmp, &
                                          material_auxvars(ghosted_id_up), &
                                          option)
        if (wippflo_jacobian_test .and. wippflo_jacobian_test_rdof > 0) then
          print *, 'up-up: ',Jtmp, local_id_up
        endif
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                      Jtmp,ADD_VALUES,ierr);CHKERRQ(ierr)
        Jtmp = Jdn
        call WIPPFloConvertUnitsToBRAGFlo(Jtmp, &
                                          material_auxvars(ghosted_id_up), &
                                          option)
        if (wippflo_jacobian_test .and. wippflo_jacobian_test_rdof > 0) then
          print *, 'up-dn: ',Jtmp, local_id_up
        endif
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                      Jtmp,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
      if (local_id_dn > 0) then
        Jup = -Jup
        Jdn = -Jdn
        Jtmp = Jdn
        call WIPPFloConvertUnitsToBRAGFlo(Jtmp, &
                                          material_auxvars(ghosted_id_dn), &
                                          option)
        if (wippflo_jacobian_test) then
          print *, 'dn-dn: ',Jtmp, local_id_dn
        endif
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                                      Jtmp,ADD_VALUES,ierr);CHKERRQ(ierr)
        Jtmp = Jup
        call WIPPFloConvertUnitsToBRAGFlo(Jtmp, &
                                          material_auxvars(ghosted_id_dn), &
                                          option)
        if (wippflo_jacobian_test) then
          print *, 'dn-up: ',Jtmp, local_id_dn
        endif
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                                      Jtmp,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  if (realization%debug%matview_Matrix_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call DebugWriteFilename(realization%debug,string,'WFjacobian_flux','', &
                            wippflo_ts_count,wippflo_ts_cut_count, &
                            wippflo_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif

  endif

  if (wippflo_calc_bcflux) then
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

      icc_dn = patch%cc_id(ghosted_id)

      call XXBCFluxDerivative(boundary_condition%flow_bc_type, &
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
      call WIPPFloConvertUnitsToBRAGFlo(Jdn,material_auxvars(ghosted_id), &
                                        option)
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  if (realization%debug%matview_Matrix_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call DebugWriteFilename(realization%debug,string,'WFjacobian_bcflux','', &
                            wippflo_ts_count,wippflo_ts_cut_count, &
                            wippflo_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif

  endif

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
                        material_auxvars(ghosted_id), &
                        scale,Jup)

      call WIPPFloConvertUnitsToBRAGFlo(Jup,material_auxvars(ghosted_id), &
                                        option)
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)

    enddo
    source_sink => source_sink%next
  enddo

  if (realization%debug%matview_Matrix_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call DebugWriteFilename(realization%debug,string,'WFjacobian_srcsink','', &
                            wippflo_ts_count,wippflo_ts_cut_count, &
                            wippflo_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif



  if (wippflo_calc_chem) then
  ! WIPP gas/brine generation process model source/sinks
  if (associated(pmwss_ptr)) then
    if (pmwss_ptr%rate_update_frequency == NO_LAG) then
      call PMWSSUpdateRates(pmwss_ptr,PETSC_TRUE,ierr)
      call PMWSSCalcJacobianValues(pmwss_ptr,A,ierr)
    endif
  endif
  endif

  ! Compute WIPP well model source/sinks for the quasi-implicitly coupled well
  ! model approach
  if (wippflo_well_quasi_imp_coupled) then
  if (associated(pmwell_ptr)) then
    if (any(pmwell_ptr%well_grid%h_rank_id == option%myrank)) then
      call pmwell_ptr%UpdateFlowRates(ONE_INTEGER,ONE_INTEGER,-999,ierr)
      if (pmwell_ptr%well_force_ts_cut == ZERO_INTEGER) then
        call pmwell_ptr%UpdateFlowRates(TWO_INTEGER,TWO_INTEGER,-999,ierr)
        if (pmwell_ptr%well_force_ts_cut == ZERO_INTEGER) then
          call pmwell_ptr%ModifyFlowJacobian(A,ierr)
        endif
      endif
    endif
  endif
  endif

  call WIPPFloSSSandbox(null_vec,A,PETSC_TRUE,grid,material_auxvars, &
                        wippflo_auxvars,option)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  ! zero out inactive cells
  if (patch%aux%WIPPFlo%inactive_cells_exist) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    call MatZeroRowsLocal(A,patch%aux%WIPPFlo%matrix_zeroing%n_inactive_rows, &
                          patch%aux%WIPPFlo%matrix_zeroing% &
                            inactive_rows_local_ghosted, &
                          qsrc,PETSC_NULL_VEC,PETSC_NULL_VEC, &
                          ierr);CHKERRQ(ierr)
  endif

  if (realization%debug%matview_Matrix) then
    call DebugWriteFilename(realization%debug,string,'WFjacobian','', &
                            wippflo_ts_count,wippflo_ts_cut_count, &
                            wippflo_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(J,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif
  if (realization%debug%norm_Matrix) then
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


  if (wippflo_jacobian_test) then
    wippflo_jacobian_test_active = PETSC_TRUE
    call WIPPFloNumericalJacobianTest(xx,A,realization,pmwss_ptr,pmwell_ptr)
    wippflo_jacobian_test_active = PETSC_FALSE
  endif

  wippflo_ni_count = wippflo_ni_count + 1

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

  class(realization_subsurface_type) :: realization
  type(output_variable_list_type), pointer :: list

  character(len=MAXWORDLENGTH) :: name, units

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

  endif

end subroutine WIPPFloSetPlotVariables

! ************************************************************************** !

subroutine WIPPFloCreepShutDown(realization)
  !
  ! Author: Jennifer M. Frederick
  ! Date: 04/18/2018
  !
  use Realization_Subsurface_class
  use Grid_module
  use Option_module
  use WIPP_module
  use Creep_Closure_module
  use Material_Aux_module, only: material_auxvar_type, MaterialAuxVarSetValue
  use Variables_module, only : SOIL_REFERENCE_PRESSURE

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  class(creep_closure_type), pointer :: creep_closure
  PetscInt :: ghosted_id
  PetscInt :: creep_closure_id
  PetscReal :: cell_pressure

  option => realization%option
  grid => realization%patch%grid
  material_auxvars => realization%patch%aux%Material%auxvars
  wippflo_auxvars => realization%patch%aux%WIPPFlo%auxvars

  do ghosted_id = 1, grid%ngmax
    creep_closure_id = material_auxvars(ghosted_id)%creep_closure_id
    creep_closure => wipp%creep_closure_tables_array(creep_closure_id)%ptr
    if (associated(creep_closure)) then
      cell_pressure = wippflo_auxvars(ZERO_INTEGER,ghosted_id)%&
                                                   pres(option%liquid_phase)
      if (option%time > creep_closure%time_datamax .or. &
          option%time > creep_closure%time_closeoff .or. &
          cell_pressure > creep_closure%shutdown_pressure) then
        material_auxvars(ghosted_id)%porosity_base = &
          wippflo_auxvars(ZERO_INTEGER,ghosted_id)%effective_porosity
        call MaterialAuxVarSetValue(material_auxvars(ghosted_id), &
                                    SOIL_REFERENCE_PRESSURE,cell_pressure)
        ! index 1 of wipp%creep_closure_tables_array is a null pointer
        ! which will shut down creep closure permanently since the pointer
        ! creep_closure => wipp%creep_closure_tables_array(creep_closure_id)%ptr
        ! will no longer be associated in future conditionals
        material_auxvars(ghosted_id)%creep_closure_id = 1
      endif
    endif
  enddo

end subroutine WIPPFloCreepShutDown

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
  use Material_Aux_module, only: material_auxvar_type
  use SrcSink_Sandbox_module
  use SrcSink_Sandbox_Base_class

  implicit none

  PetscBool :: compute_derivative
  Vec :: residual
  Mat :: Jacobian
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)

  type(grid_type) :: grid
  type(option_type) :: option

  PetscReal, pointer :: r_p(:)
  PetscReal :: res(option%nflowdof)
  PetscReal :: Jac(option%nflowdof,option%nflowdof)
  class(srcsink_sandbox_base_type), pointer :: cur_srcsink
  PetscInt :: local_id, ghosted_id, istart, iend, irow, idof, icell
  PetscReal :: res_pert(option%nflowdof)
  PetscReal :: aux_real(10)
  PetscErrorCode :: ierr

  if (.not.compute_derivative) then
    call VecGetArrayF90(residual,r_p,ierr);CHKERRQ(ierr)
  endif

  cur_srcsink => ss_sandbox_list
  do
    if (.not.associated(cur_srcsink)) exit
    do icell = 1, size(cur_srcsink%local_cell_ids)
      local_id = cur_srcsink%local_cell_ids(icell)
      ghosted_id = grid%nL2G(local_id)
      aux_real = 0.d0
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
        call MatSetValuesBlockedLocal(Jacobian,1,ghosted_id-1,1,ghosted_id-1, &
                                      Jac,ADD_VALUES,ierr);CHKERRQ(ierr)
      else
        iend = local_id*option%nflowdof
        istart = iend - option%nflowdof + 1
        r_p(istart:iend) = r_p(istart:iend) - res
      endif
    enddo
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

  class(realization_subsurface_type) :: realization

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

  class(realization_subsurface_type) :: realization

  ! place anything that needs to be freed here.
  ! auxvars are deallocated in auxiliary.F90.

end subroutine WIPPFloDestroy

end module WIPP_Flow_module
