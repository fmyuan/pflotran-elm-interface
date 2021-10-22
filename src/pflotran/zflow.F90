module ZFlow_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use ZFlow_Aux_module
  use ZFlow_Common_module
  use Global_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: ZFlowSetup, &
            ZFlowInitializeTimestep, &
            ZFlowUpdateSolution, &
            ZFlowTimeCut,&
            ZFlowUpdateAuxVars, &
            ZFlowUpdateFixedAccum, &
            ZFlowComputeMassBalance, &
            ZFlowZeroMassBalanceDelta, &
            ZFlowResidual, &
            ZFlowSetPlotVariables, &
            ZFlowMapBCAuxVarsToGlobal, &
            ZFlowDestroy

contains

! ************************************************************************** !

subroutine ZFlowSetup(realization)
  !
  ! Creates arrays for auxiliary variables
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
  use Material_Aux_class
  use Output_Aux_module
  use Characteristic_Curves_module
  use Matrix_Zeroing_module
  use EOS_Water_module

  implicit none

  type(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(output_variable_list_type), pointer :: list
  type(material_parameter_type), pointer :: material_parameter

  PetscInt :: ghosted_id, iconn, sum_connection, local_id
  PetscBool :: error_found
  PetscInt :: flag(10)
  PetscInt :: temp_int, idof
  PetscErrorCode :: ierr
                                                ! extra index for derivatives
  type(zflow_auxvar_type), pointer :: zflow_auxvars(:,:)
  type(zflow_auxvar_type), pointer :: zflow_auxvars_bc(:)
  type(zflow_auxvar_type), pointer :: zflow_auxvars_ss(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)

  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%aux%ZFlow => ZFlowAuxCreate(option)

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
                     MPI_LOR,option%mycomm,ierr)
  if (error_found) then
    option%io_buffer = 'Material property errors found in ZFlowSetup.'
    call PrintErrMsg(option)
  endif

  zflow_numerical_derivatives = option%flow%numerical_derivatives
  temp_int = 0
  if (zflow_numerical_derivatives) temp_int = 1
  allocate(zflow_auxvars(0:temp_int,grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    do idof = 0, temp_int
      call ZFlowAuxVarInit(zflow_auxvars(idof,ghosted_id),option)
    enddo
  enddo
  patch%aux%ZFlow%auxvars => zflow_auxvars
  patch%aux%ZFlow%num_aux = grid%ngmax

  ! count the number of boundary connections and allocate
  ! auxvar data structures for them
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_condition_list)
  if (sum_connection > 0) then
    allocate(zflow_auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call ZFlowAuxVarInit(zflow_auxvars_bc(iconn),option)
    enddo
    patch%aux%ZFlow%auxvars_bc => zflow_auxvars_bc
  endif
  patch%aux%ZFlow%num_aux_bc = sum_connection

  ! count the number of source/sink connections and allocate
  ! auxvar data structures for them
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (sum_connection > 0) then
    allocate(zflow_auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call ZFlowAuxVarInit(zflow_auxvars_ss(iconn),option)
    enddo
    patch%aux%ZFlow%auxvars_ss => zflow_auxvars_ss
  endif
  patch%aux%ZFlow%num_aux_ss = sum_connection

  list => realization%output_option%output_snap_variable_list
  call ZFlowSetPlotVariables(realization,list)
  list => realization%output_option%output_obs_variable_list
  call ZFlowSetPlotVariables(realization,list)

  XXFlux => ZFlowFluxHarmonicPermOnly
  XXBCFlux => ZFlowBCFluxHarmonicPermOnly

  zflow_ts_count = 0
  zflow_ts_cut_count = 0
  zflow_ni_count = 0

end subroutine ZFlowSetup

! ************************************************************************** !

subroutine ZFlowInitializeTimestep(realization)
  !
  ! Update data in module prior to time step
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !
  use Realization_Subsurface_class
  use Upwind_Direction_module

  implicit none

  type(realization_subsurface_type) :: realization

  call ZFlowUpdateFixedAccum(realization)
  zflow_ni_count = 0

end subroutine ZFlowInitializeTimestep

! ************************************************************************** !

subroutine ZFlowUpdateSolution(realization)
  !
  ! Updates data in module after a successful time
  ! step
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use Realization_Subsurface_class

  implicit none

  type(realization_subsurface_type) :: realization

  if (realization%option%compute_mass_balance_new) then
    call ZFlowUpdateMassBalance(realization)
  endif

  zflow_ts_count = zflow_ts_count + 1
  zflow_ts_cut_count = 0
  zflow_ni_count = 0

end subroutine ZFlowUpdateSolution

! ************************************************************************** !

subroutine ZFlowTimeCut(realization)
  !
  ! Resets arrays for time step cut
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !
  use Realization_Subsurface_class

  implicit none

  type(realization_subsurface_type) :: realization

  zflow_ts_cut_count = zflow_ts_cut_count + 1

  call ZFlowInitializeTimestep(realization)

end subroutine ZFlowTimeCut

! ************************************************************************** !

subroutine ZFlowComputeMassBalance(realization,mass_balance)
  !
  ! Initializes mass balance
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
  use Material_Aux_class

  implicit none

  type(realization_subsurface_type) :: realization
  PetscReal :: mass_balance(1)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(zflow_auxvar_type), pointer :: zflow_auxvars(:,:)
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt, parameter :: iphase = 1

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  zflow_auxvars => patch%aux%ZFlow%auxvars
  material_auxvars => patch%aux%Material%auxvars

  mass_balance = 0.d0

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    ! volume_phase = saturation*porosity*volume
    mass_balance(iphase) = mass_balance(iphase) + &
        zflow_auxvars(ZERO_INTEGER,ghosted_id)%sat* &
        zflow_auxvars(ZERO_INTEGER,ghosted_id)%effective_porosity* &
        material_auxvars(ghosted_id)%volume
  enddo

end subroutine ZFlowComputeMassBalance

! ************************************************************************** !

subroutine ZFlowZeroMassBalanceDelta(realization)
  !
  ! Zeros mass balance delta array
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
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

  do iconn = 1, patch%aux%ZFlow%num_aux_bc
    global_auxvars_bc(iconn)%mass_balance_delta = 0.d0
  enddo
  do iconn = 1, patch%aux%ZFlow%num_aux_ss
    global_auxvars_ss(iconn)%mass_balance_delta = 0.d0
  enddo

end subroutine ZFlowZeroMassBalanceDelta

! ************************************************************************** !

subroutine ZFlowUpdateMassBalance(realization)
  !
  ! Updates mass balance
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
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

  do iconn = 1, patch%aux%ZFlow%num_aux_bc
    global_auxvars_bc(iconn)%mass_balance(1,1) = &
      global_auxvars_bc(iconn)%mass_balance(1,1) + &
      global_auxvars_bc(iconn)%mass_balance_delta(1,1)*option%flow_dt
  enddo
  do iconn = 1, patch%aux%ZFlow%num_aux_ss
    global_auxvars_ss(iconn)%mass_balance(1,1) = &
      global_auxvars_ss(iconn)%mass_balance(1,1) + &
      global_auxvars_ss(iconn)%mass_balance_delta(1,1)*option%flow_dt
  enddo

end subroutine ZFlowUpdateMassBalance

! ************************************************************************** !

subroutine ZFlowUpdateAuxVars(realization)
  !
  ! Updates the auxiliary variables associated with the ZFlow problem
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
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
  type(zflow_auxvar_type), pointer :: zflow_auxvars(:,:)
  type(zflow_auxvar_type), pointer :: zflow_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: ghosted_id, local_id, sum_connection, iconn, natural_id
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscErrorCode :: ierr

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  zflow_auxvars => patch%aux%ZFlow%auxvars
  zflow_auxvars_bc => patch%aux%ZFlow%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars

  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells

    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    ! ZFLOW_UPDATE_FOR_ACCUM indicates call from non-perturbation
    option%iflag = ZFLOW_UPDATE_FOR_ACCUM
    natural_id = grid%nG2A(ghosted_id)
    if (grid%nG2L(ghosted_id) == 0) natural_id = -natural_id
    call ZFlowAuxVarCompute(xx_loc_p(ghosted_id:ghosted_id), &
                            zflow_auxvars(ZERO_INTEGER,ghosted_id), &
                            global_auxvars(ghosted_id), &
                            material_auxvars(ghosted_id), &
                            patch%characteristic_curves_array( &
                              patch%cc_id(ghosted_id))%ptr, &
                            natural_id, &
                            PETSC_TRUE,option)
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
      if (patch%imat(ghosted_id) <= 0) cycle

      select case(boundary_condition%flow_condition% &
                    itype(ZFLOW_PRESSURE_DOF))
        case(DIRICHLET_BC, DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
             HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC)
          xxbc(1) = boundary_condition% &
                      flow_aux_real_var(ZFLOW_PRESSURE_DOF,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC,UNIT_GRADIENT_BC, &
             SURFACE_ZERO_GRADHEIGHT)
          xxbc(1) = xx_loc_p(local_id)
        case default
          option%io_buffer = 'boundary itype not set up in ZFlowUpdateAuxVars'
          call PrintErrMsg(option)
      end select

      ! ZFLOW_UPDATE_FOR_BOUNDARY indicates call from non-perturbation
      option%iflag = ZFLOW_UPDATE_FOR_BOUNDARY
      call ZFlowAuxVarCompute(xxbc,zflow_auxvars_bc(sum_connection), &
                              global_auxvars_bc(sum_connection), &
                              material_auxvars(ghosted_id), &
                              patch%characteristic_curves_array( &
                                patch%cc_id(ghosted_id))%ptr, &
                              natural_id, &
                              PETSC_FALSE,option)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  patch%aux%ZFlow%auxvars_up_to_date = PETSC_TRUE

end subroutine ZFlowUpdateAuxVars

! ************************************************************************** !

subroutine ZFlowUpdateFixedAccum(realization)
  !
  ! Updates the fixed portion of the
  ! accumulation term
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
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
  type(zflow_auxvar_type), pointer :: zflow_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(material_parameter_type), pointer :: material_parameter

  PetscInt :: ghosted_id, local_id, local_start, local_end, natural_id
  PetscInt :: imat
  PetscReal, pointer :: xx_p(:)
  PetscReal, pointer :: accum_p(:)
  PetscReal :: Jdum(1,1)

  PetscErrorCode :: ierr

  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  zflow_auxvars => patch%aux%ZFlow%auxvars
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
    local_end = local_id * option%nflowdof
    local_start = local_end - option%nflowdof + 1
    natural_id = grid%nG2A(ghosted_id)
    ! ZFLOW_UPDATE_FOR_FIXED_ACCUM indicates call from non-perturbation
    option%iflag = ZFLOW_UPDATE_FOR_FIXED_ACCUM
    call ZFlowAuxVarCompute(xx_p(local_start:local_end), &
                            zflow_auxvars(ZERO_INTEGER,ghosted_id), &
                            global_auxvars(ghosted_id), &
                            material_auxvars(ghosted_id), &
                            patch%characteristic_curves_array( &
                              patch%cc_id(ghosted_id))%ptr, &
                            natural_id, &
                            PETSC_TRUE,option)
    call ZFlowAccumulation(zflow_auxvars(ZERO_INTEGER,ghosted_id), &
                           global_auxvars(ghosted_id), &
                           material_auxvars(ghosted_id), &
                           option,accum_p(local_id:local_id), &
                           Jdum,PETSC_FALSE)
  enddo

  call VecRestoreArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

end subroutine ZFlowUpdateFixedAccum

! ************************************************************************** !

subroutine ZFlowResidual(snes,xx,r,A,realization,ierr)
  !
  ! Computes the residual equation
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
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
  use Material_Aux_class
  use Upwind_Direction_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  Mat :: A
  type(realization_subsurface_type) :: realization
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
  type(zflow_parameter_type), pointer :: zflow_parameter
  type(zflow_auxvar_type), pointer :: zflow_auxvars(:,:)
  type(zflow_auxvar_type), pointer :: zflow_auxvars_bc(:)
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
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
  PetscInt :: i, imat, imat_up, imat_dn

  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum_p(:), accum_p2(:)
  PetscReal, pointer :: vec_p(:)

  character(len=MAXSTRINGLENGTH) :: string

  PetscInt :: icc_up, icc_dn
  PetscReal :: Res(1)
  PetscReal :: Jup(1,1),Jdn(1,1)
  PetscReal :: dJupdKup(1),dJupdKdn(1),dJdndKup(1),dJdndKdn(1)
  PetscReal :: drhsdKup(1),drhsdKdn(1)
  Mat, pointer :: dMdK(:)
  Vec, pointer :: dbdK(:)
  PetscReal :: v_darcy(1)
  PetscReal :: dM(2)

  dJupdKup = 0.d0
  dJupdKdn = 0.d0
  dJdndKdn = 0.d0
  dJdndKup = 0.d0
  drhsdKup = 0.d0
  drhsdKdn = 0.d0

  discretization => realization%discretization
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  zflow_auxvars => patch%aux%ZFlow%auxvars
  zflow_auxvars_bc => patch%aux%ZFlow%auxvars_bc
  zflow_parameter => patch%aux%ZFlow%zflow_parameter
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars

  if (associated(patch%aux%inversion_aux)) then
    dMdK => patch%aux%inversion_aux%dMdK
    dbdK => patch%aux%inversion_aux%dbdK
  else
    nullify(dMdK)
    nullify(dbdK)
  endif

  if (zflow_simult_function_evals) then
    call MatZeroEntries(A,ierr);CHKERRQ(ierr)
    if (associated(dMdK)) then
      do i = 1, size(dMdK)
        call MatZeroEntries(dMdK(i),ierr);CHKERRQ(ierr)
      enddo
      do i = 1, size(dbdK)
        call VecZeroEntries(dbdK(i),ierr);CHKERRQ(ierr)
      enddo
    endif
  endif

  ! Communication -----------------------------------------
  ! must be called before ZFlowUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  call ZFlowUpdateAuxVars(realization)

  ! override flags since they will soon be out of date
  patch%aux%ZFlow%auxvars_up_to_date = PETSC_FALSE

  if (zflow_numerical_derivatives) then
        ! Perturb aux vars
    do ghosted_id = 1, grid%ngmax  ! For each local node do...
      if (patch%imat(ghosted_id) <= 0) cycle
      call ZFlowAuxVarPerturb(zflow_auxvars(:,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              patch%characteristic_curves_array( &
                                patch%cc_id(ghosted_id))%ptr, &
                              grid%nG2A(ghosted_id),option)
    enddo
  endif

  if (option%compute_mass_balance_new) then
    call ZFlowZeroMassBalanceDelta(realization)
  endif

  option%iflag = 1
  ! now assign access pointer to local variables
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)

  ! Accumulation terms ------------------------------------
  ! accumulation at t(k) (doesn't change during Newton iteration)
  if (zflow_calc_accum) then
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
      call ZFlowAccumulation(zflow_auxvars(ZERO_INTEGER,ghosted_id), &
                             global_auxvars(ghosted_id), &
                             material_auxvars(ghosted_id), &
                             option,Res,Jup, &
                             zflow_simult_function_evals)
      if (zflow_numerical_derivatives) then
        call ZFlowAccumDerivative(zflow_auxvars(:,ghosted_id), &
                                  global_auxvars(ghosted_id), &
                                  material_auxvars(ghosted_id), &
                                  option,Jup)
      endif
      r_p(local_id) =  r_p(local_id) + Res(1)
      accum_p2(local_id) = Res(1)
      if (zflow_simult_function_evals) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                      ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
    enddo
    call VecRestoreArrayF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  else
    r_p = 0.d0
  endif

  if (zflow_calc_flux) then
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

        call XXFlux(zflow_auxvars(ZERO_INTEGER,ghosted_id_up), &
                    global_auxvars(ghosted_id_up), &
                    material_auxvars(ghosted_id_up), &
                    zflow_auxvars(ZERO_INTEGER,ghosted_id_dn), &
                    global_auxvars(ghosted_id_dn), &
                    material_auxvars(ghosted_id_dn), &
                    cur_connection_set%area(iconn), &
                    cur_connection_set%dist(:,iconn), &
                    zflow_parameter,option,v_darcy, &
                    Res,Jup,Jdn, &
                    dJupdKup,dJupdKdn,dJdndKup,dJdndKdn, &
                    drhsdKup,drhsdKdn, &
                    zflow_simult_function_evals)
        if (zflow_numerical_derivatives) then
          call XXFluxDerivative(zflow_auxvars(:,ghosted_id_up), &
                                global_auxvars(ghosted_id_up), &
                                material_auxvars(ghosted_id_up), &
                                zflow_auxvars(:,ghosted_id_dn), &
                                global_auxvars(ghosted_id_dn), &
                                material_auxvars(ghosted_id_dn), &
                                cur_connection_set%area(iconn), &
                                cur_connection_set%dist(:,iconn), &
                                zflow_parameter,option, &
                                Jup,Jdn)
        endif
        patch%internal_velocities(:,sum_connection) = v_darcy
        if (associated(patch%internal_flow_fluxes)) then
          patch%internal_flow_fluxes(:,sum_connection) = Res(:)
        endif

        if (zflow_simult_function_evals .and. &
            associated(patch%aux%inversion_aux)) then
          patch%aux%inversion_aux%dFluxdIntConn(:,sum_connection) = &
            [dJupdKup(1),dJupdKdn(1),dJdndKup(1),dJdndKdn(1), &
             drhsdKup(1),drhsdKdn(1)]
        endif

        if (local_id_up > 0) then
          r_p(local_id_up) = r_p(local_id_up) + Res(1)
          if (zflow_simult_function_evals) then
            call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1, &
                                          1,ghosted_id_up-1, &
                                          Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
            call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1, &
                                          1,ghosted_id_dn-1, &
                                          Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
            if (associated(dMdK)) then
              ! flip sign for rhs of equation
              drhsdKup = -drhsdKup
              drhsdKdn = -drhsdKdn
              call MatSetValuesBlockedLocal(dMdK(local_id_up),1,ghosted_id_up-1, &
                                            1,ghosted_id_up-1, &
                                            dJupdKup,ADD_VALUES,ierr);CHKERRQ(ierr)
              call MatSetValuesBlockedLocal(dMdK(local_id_up),1,ghosted_id_up-1, &
                                            1,ghosted_id_dn-1, &
                                            dJdndKup,ADD_VALUES,ierr);CHKERRQ(ierr)
              call MatSetValuesBlockedLocal(dMdK(local_id_dn),1,ghosted_id_up-1, &
                                            1,ghosted_id_up-1, &
                                            dJupdKdn,ADD_VALUES,ierr);CHKERRQ(ierr)
              call MatSetValuesBlockedLocal(dMdK(local_id_dn),1,ghosted_id_up-1, &
                                            1,ghosted_id_dn-1, &
                                            dJdndKdn,ADD_VALUES,ierr);CHKERRQ(ierr)
              call VecSetValue(dbdK(local_id_up),ghosted_id_up-1, &
                               drhsdKup(1),ADD_VALUES,ierr);CHKERRQ(ierr)
              call VecSetValue(dbdK(local_id_dn),ghosted_id_up-1, &
                               drhsdKdn(1),ADD_VALUES,ierr);CHKERRQ(ierr)
            endif
          endif
        endif

        if (local_id_dn > 0) then
          r_p(local_id_dn) = r_p(local_id_dn) - Res(1)
          if (zflow_simult_function_evals) then
            Jup = -Jup
            Jdn = -Jdn
            call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1, &
                                          1,ghosted_id_dn-1, &
                                          Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
            call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1, &
                                          1,ghosted_id_up-1, &
                                          Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
            if (associated(dMdK)) then
              dJupdKup = -dJupdKup
              dJdndKup = -dJdndKup
              dJupdKdn = -dJupdKdn
              dJdndKdn = -dJdndKdn
              drhsdKup = -drhsdKup
              drhsdKdn = -drhsdKdn
              call MatSetValuesBlockedLocal(dMdK(local_id_up),1,ghosted_id_dn-1, &
                                            1,ghosted_id_dn-1, &
                                            dJdndKup,ADD_VALUES,ierr);CHKERRQ(ierr)
              call MatSetValuesBlockedLocal(dMdK(local_id_up),1,ghosted_id_dn-1, &
                                            1,ghosted_id_up-1, &
                                            dJupdKup,ADD_VALUES,ierr);CHKERRQ(ierr)
              call MatSetValuesBlockedLocal(dMdK(local_id_dn),1,ghosted_id_dn-1, &
                                            1,ghosted_id_dn-1, &
                                            dJdndKdn,ADD_VALUES,ierr);CHKERRQ(ierr)
              call MatSetValuesBlockedLocal(dMdK(local_id_dn),1,ghosted_id_dn-1, &
                                            1,ghosted_id_up-1, &
                                            dJupdKdn,ADD_VALUES,ierr);CHKERRQ(ierr)
              call VecSetValue(dbdK(local_id_up),ghosted_id_dn-1, &
                               drhsdKup(1),ADD_VALUES,ierr);CHKERRQ(ierr)
              call VecSetValue(dbdK(local_id_dn),ghosted_id_dn-1, &
                               drhsdKdn(1),ADD_VALUES,ierr);CHKERRQ(ierr)
           endif
          endif
        endif
      enddo

      cur_connection_set => cur_connection_set%next
    enddo
  endif

  if (zflow_calc_bcflux) then
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
                      zflow_auxvars_bc(sum_connection), &
                      global_auxvars_bc(sum_connection), &
                      zflow_auxvars(ZERO_INTEGER,ghosted_id), &
                      global_auxvars(ghosted_id), &
                      material_auxvars(ghosted_id), &
                      cur_connection_set%area(iconn), &
                      cur_connection_set%dist(:,iconn), &
                      zflow_parameter,option, &
                      v_darcy,Res,Jdn, &
                      dJdndKdn,drhsdKdn, &
                      zflow_simult_function_evals)
        if (zflow_numerical_derivatives) then
          call XXBCFluxDerivative(boundary_condition%flow_bc_type, &
                                  boundary_condition%flow_aux_mapping, &
                                  boundary_condition% &
                                    flow_aux_real_var(:,iconn), &
                                  zflow_auxvars_bc(sum_connection), &
                                  global_auxvars_bc(sum_connection), &
                                  zflow_auxvars(:,ghosted_id), &
                                  global_auxvars(ghosted_id), &
                                  material_auxvars(ghosted_id), &
                                  cur_connection_set%area(iconn), &
                                  cur_connection_set%dist(:,iconn), &
                                  zflow_parameter,option, &
                                  Jdn)
        endif
        patch%boundary_velocities(:,sum_connection) = v_darcy
        if (associated(patch%boundary_flow_fluxes)) then
          patch%boundary_flow_fluxes(:,sum_connection) = Res(:)
        endif
        if (option%compute_mass_balance_new) then
          ! contribution to boundary
          global_auxvars_bc(sum_connection)%mass_balance_delta(:,1) = &
            global_auxvars_bc(sum_connection)%mass_balance_delta(:,1) - &
            Res(:)
        endif
        r_p(local_id)= r_p(local_id) - Res(1)
        if (zflow_simult_function_evals) then
          Jdn = -Jdn
          call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn, &
                                        ADD_VALUES,ierr);CHKERRQ(ierr)
          if (associated(patch%aux%inversion_aux)) then
            dJdndKdn = -dJdndKdn
            ! no need to flip sign on rhs since downwind
            patch%aux%inversion_aux%dFluxdBCConn(:,sum_connection) = &
              [dJdndKdn(1),drhsdKdn(1)]
          endif
          if (associated(dMdK)) then
            call MatSetValuesBlockedLocal(dMdK(local_id),1,ghosted_id-1, &
                                          1,ghosted_id-1, &
                                          dJdndKdn,ADD_VALUES,ierr);CHKERRQ(ierr)
            call VecSetValue(dbdK(local_id),ghosted_id-1, &
                             drhsdKdn(1),ADD_VALUES,ierr);CHKERRQ(ierr)
          endif
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

      if (associated(source_sink%flow_aux_real_var)) then
        scale = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
      else
        scale = 1.d0
      endif

      call ZFlowSrcSink(option,source_sink%flow_condition%rate% &
                                  dataset%rarray(:), &
                        source_sink%flow_condition%rate%itype, &
                        zflow_auxvars(ZERO_INTEGER,ghosted_id), &
                        global_auxvars(ghosted_id), &
                        material_auxvars(ghosted_id), &
                        ss_flow_vol_flux, &
                        scale,Res,Jdn, &
                        zflow_simult_function_evals)
      if (zflow_numerical_derivatives) then
        call ZFlowSrcSinkDerivative(option, &
                                    source_sink%flow_condition%rate% &
                                      dataset%rarray(:), &
                                    source_sink%flow_condition%rate%itype, &
                                    zflow_auxvars(:,ghosted_id), &
                                    global_auxvars(ghosted_id), &
                                    material_auxvars(ghosted_id), &
                                    scale,Jdn)
      endif
      if (associated(patch%ss_flow_vol_fluxes)) then
        patch%ss_flow_vol_fluxes(:,sum_connection) = ss_flow_vol_flux
      endif
      if (associated(patch%ss_flow_fluxes)) then
        patch%ss_flow_fluxes(:,sum_connection) = Res(:)
      endif
      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_ss(sum_connection)%mass_balance_delta(:,1) = &
          global_auxvars_ss(sum_connection)%mass_balance_delta(:,1) - &
          Res(:)
      endif
      r_p(local_id) =  r_p(local_id) - Res(1)
      if (zflow_simult_function_evals) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn, &
                                      ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
    enddo
    source_sink => source_sink%next
  enddo

  if (patch%aux%ZFlow%inactive_cells_exist) then
    do i=1,patch%aux%ZFlow%matrix_zeroing%n_inactive_rows
      r_p(patch%aux%ZFlow%matrix_zeroing%inactive_rows_local(i)) = 0.d0
    enddo
  endif

  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)

  if (zflow_simult_function_evals) then

    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      ! zero out inactive cells

    if (associated(dMdK)) then
      do i = 1, size(dMdK)
        call MatAssemblyBegin(dMdK(i),MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
        call MatAssemblyEnd(dMdK(i),MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      enddo
      do i = 1, size(dbdK)
        call VecAssemblyBegin(dbdK(i),ierr);CHKERRQ(ierr)
        call VecAssemblyEnd(dbdK(i),ierr);CHKERRQ(ierr)
      enddo
    endif

    if (patch%aux%ZFlow%inactive_cells_exist) then
      scale = 1.d0 ! solely a temporary variable in this conditional
      call MatZeroRowsLocal(A,patch%aux%ZFlow%matrix_zeroing%n_inactive_rows, &
                            patch%aux%ZFlow%matrix_zeroing% &
                              inactive_rows_local_ghosted, &
                            scale,PETSC_NULL_VEC,PETSC_NULL_VEC, &
                            ierr);CHKERRQ(ierr)
    endif

    if (realization%debug%matview_Matrix) then
      call DebugWriteFilename(realization%debug,string,'ZFjacobian','', &
                              zflow_ts_count,zflow_ts_cut_count, &
                              zflow_ni_count)
      call DebugCreateViewer(realization%debug,string,option,viewer)
      call MatView(A,viewer,ierr);CHKERRQ(ierr)
      call DebugViewerDestroy(realization%debug,viewer)
    endif
  endif

  ! Mass Transfer
  if (field%flow_mass_transfer /= PETSC_NULL_VEC) then
    ! scale by -1.d0 for contribution to residual.  A negative contribution
    ! indicates mass being added to system.
    call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%flow_mass_transfer,vec_p,ierr);CHKERRQ(ierr)
    ! geh: leave in expanded do loop form instead of VecAXPY for flexibility
    !      in the future
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      imat = patch%imat(ghosted_id)
      if (imat <= 0) cycle
      r_p(local_id) = r_p(local_id) - vec_p(local_id)
    enddo
    call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%flow_mass_transfer,vec_p,ierr);CHKERRQ(ierr)
!geh: due to the potential for units conversion, cannot VecAXPY
!    call VecAXPY(r,-1.d0,field%flow_mass_transfer,ierr);CHKERRQ(ierr)
  endif

  if (realization%debug%vecview_residual) then
    call DebugWriteFilename(realization%debug,string,'ZFresidual','', &
                            zflow_ts_count,zflow_ts_cut_count, &
                            zflow_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif
  if (realization%debug%vecview_solution) then
    call DebugWriteFilename(realization%debug,string,'ZFxx','', &
                            zflow_ts_count,zflow_ts_cut_count, &
                            zflow_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif

end subroutine ZFlowResidual

! ************************************************************************** !

subroutine ZFlowSetPlotVariables(realization,list)
  !
  ! Adds variables to be printed to list
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use Realization_Subsurface_class
  use Output_Aux_module
  use Variables_module

  implicit none

  type(realization_subsurface_type) :: realization
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

    name = 'Liquid Saturation'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                                 LIQUID_SATURATION)
  endif

end subroutine ZFlowSetPlotVariables

! ************************************************************************** !

subroutine ZFlowMapBCAuxVarsToGlobal(realization)
  !
  ! Maps variables in zflow auxvar to global equivalent.
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
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
  type(zflow_auxvar_type), pointer :: zflow_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)

  PetscInt :: sum_connection, iconn

  option => realization%option
  patch => realization%patch

  if (option%ntrandof == 0) return ! no need to update

  zflow_auxvars_bc => patch%aux%ZFlow%auxvars_bc
  global_auxvars_bc => patch%aux%Global%auxvars_bc

  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      global_auxvars_bc(sum_connection)%sat = &
        zflow_auxvars_bc(sum_connection)%sat
    enddo
    boundary_condition => boundary_condition%next
  enddo

end subroutine ZFlowMapBCAuxVarsToGlobal

! ************************************************************************** !

subroutine ZFlowDestroy(realization)
  !
  ! Deallocates variables associated with Richard
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !

  use Realization_Subsurface_class
  use Option_module

  implicit none

  type(realization_subsurface_type) :: realization

  ! place anything that needs to be freed here.
  ! auxvars are deallocated in auxiliary.F90.

end subroutine ZFlowDestroy

end module ZFlow_module
