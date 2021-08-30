module PNF_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PNF_Aux_module
  use Global_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: PNFSetup, &
            PNFInitializeTimestep, &
            PNFUpdateSolution, &
            PNFTimeCut,&
            PNFUpdateAuxVars, &
            PNFUpdateFixedAccum, &
            PNFComputeMassBalance, &
            PNFZeroMassBalanceDelta, &
            PNFResidual, &
            PNFSetPlotVariables, &
            PNFDestroy

contains

! ************************************************************************** !

subroutine PNFSetup(realization)
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
  PetscErrorCode :: ierr
                                                ! extra index for derivatives
  type(pnf_auxvar_type), pointer :: pnf_auxvars(:)
  type(pnf_auxvar_type), pointer :: pnf_auxvars_bc(:)
  type(pnf_auxvar_type), pointer :: pnf_auxvars_ss(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)

  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%aux%PNF => PNFAuxCreate(option)

  allocate(pnf_auxvars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call PNFAuxVarInit(pnf_auxvars(ghosted_id),option)
  enddo
  patch%aux%PNF%auxvars => pnf_auxvars
  patch%aux%PNF%num_aux = grid%ngmax

  ! count the number of boundary connections and allocate
  ! auxvar data structures for them
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_condition_list)
  if (sum_connection > 0) then
    allocate(pnf_auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call PNFAuxVarInit(pnf_auxvars_bc(iconn),option)
    enddo
    patch%aux%PNF%auxvars_bc => pnf_auxvars_bc
  endif
  patch%aux%PNF%num_aux_bc = sum_connection

  ! count the number of source/sink connections and allocate
  ! auxvar data structures for them
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (sum_connection > 0) then
    allocate(pnf_auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call PNFAuxVarInit(pnf_auxvars_ss(iconn),option)
    enddo
    patch%aux%PNF%auxvars_ss => pnf_auxvars_ss
  endif
  patch%aux%PNF%num_aux_ss = sum_connection

  list => realization%output_option%output_snap_variable_list
  call PNFSetPlotVariables(realization,list)
  list => realization%output_option%output_obs_variable_list
  call PNFSetPlotVariables(realization,list)

  PNF_ts_count = 0
  PNF_ts_cut_count = 0
  PNF_ni_count = 0

end subroutine PNFSetup

! ************************************************************************** !

subroutine PNFInitializeTimestep(realization)
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

  call PNFUpdateFixedAccum(realization)
  PNF_ni_count = 0

end subroutine PNFInitializeTimestep

! ************************************************************************** !

subroutine PNFUpdateSolution(realization)
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
    call PNFUpdateMassBalance(realization)
  endif

  PNF_ts_count = PNF_ts_count + 1
  PNF_ts_cut_count = 0
  PNF_ni_count = 0

end subroutine PNFUpdateSolution

! ************************************************************************** !

subroutine PNFTimeCut(realization)
  !
  ! Resets arrays for time step cut
  !
  ! Author: Glenn Hammond
  ! Date: 08/13/21
  !
  use Realization_Subsurface_class

  implicit none

  type(realization_subsurface_type) :: realization

  PNF_ts_cut_count = PNF_ts_cut_count + 1

  call PNFInitializeTimestep(realization)

end subroutine PNFTimeCut

! ************************************************************************** !

subroutine PNFComputeMassBalance(realization,mass_balance)
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

end subroutine PNFComputeMassBalance

! ************************************************************************** !

subroutine PNFZeroMassBalanceDelta(realization)
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

  do iconn = 1, patch%aux%PNF%num_aux_bc
    global_auxvars_bc(iconn)%mass_balance_delta = 0.d0
  enddo
  do iconn = 1, patch%aux%PNF%num_aux_ss
    global_auxvars_ss(iconn)%mass_balance_delta = 0.d0
  enddo

end subroutine PNFZeroMassBalanceDelta

! ************************************************************************** !

subroutine PNFUpdateMassBalance(realization)
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

  do iconn = 1, patch%aux%PNF%num_aux_bc
    global_auxvars_bc(iconn)%mass_balance(1,1) = &
      global_auxvars_bc(iconn)%mass_balance(1,1) + &
      global_auxvars_bc(iconn)%mass_balance_delta(1,1)*option%flow_dt
  enddo
  do iconn = 1, patch%aux%PNF%num_aux_ss
    global_auxvars_ss(iconn)%mass_balance(1,1) = &
      global_auxvars_ss(iconn)%mass_balance(1,1) + &
      global_auxvars_ss(iconn)%mass_balance_delta(1,1)*option%flow_dt
  enddo

end subroutine PNFUpdateMassBalance

! ************************************************************************** !

subroutine PNFUpdateAuxVars(realization)
  !
  ! Updates the auxiliary variables associated with the PNF problem
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
  type(pnf_auxvar_type), pointer :: pnf_auxvars(:)
  type(pnf_auxvar_type), pointer :: pnf_auxvars_bc(:)
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

  pnf_auxvars => patch%aux%PNF%auxvars
  pnf_auxvars_bc => patch%aux%PNF%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars

  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells

    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    ! PNF_UPDATE_FOR_ACCUM indicates call from non-perturbation
    option%iflag = PNF_UPDATE_FOR_ACCUM
    natural_id = grid%nG2A(ghosted_id)
    if (grid%nG2L(ghosted_id) == 0) natural_id = -natural_id
    call PNFAuxVarCompute(xx_loc_p(ghosted_id:ghosted_id), &
                            pnf_auxvars(ghosted_id), &
                            global_auxvars(ghosted_id), &
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
      if (patch%imat(ghosted_id) <= 0) cycle

      select case(boundary_condition%flow_condition% &
                    itype(RICHARDS_PRESSURE_DOF))
        case(DIRICHLET_BC, DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
             HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC)
          xxbc(1) = boundary_condition% &
                      flow_aux_real_var(RICHARDS_PRESSURE_DOF,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC,UNIT_GRADIENT_BC, &
             SURFACE_ZERO_GRADHEIGHT)
          xxbc(1) = xx_loc_p(local_id)
        case default
          option%io_buffer = 'boundary itype not set up in PNFUpdateAuxVars'
          call PrintErrMsg(option)
      end select

      ! PNF_UPDATE_FOR_BOUNDARY indicates call from non-perturbation
      option%iflag = PNF_UPDATE_FOR_BOUNDARY
      call PNFAuxVarCompute(xxbc,pnf_auxvars_bc(sum_connection), &
                              global_auxvars_bc(sum_connection), &
                              natural_id, &
                              option)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p, ierr);CHKERRQ(ierr)

  patch%aux%PNF%auxvars_up_to_date = PETSC_TRUE

end subroutine PNFUpdateAuxVars

! ************************************************************************** !

subroutine PNFUpdateFixedAccum(realization)
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
  type(pnf_auxvar_type), pointer :: pnf_auxvars(:)
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

  pnf_auxvars => patch%aux%PNF%auxvars
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
    ! PNF_UPDATE_FOR_FIXED_ACCUM indicates call from non-perturbation
    option%iflag = PNF_UPDATE_FOR_FIXED_ACCUM
    call PNFAuxVarCompute(xx_p(local_start:local_end), &
                            pnf_auxvars(ghosted_id), &
                            global_auxvars(ghosted_id), &
                            natural_id, &
                            option)
    ! call PNFAccumulation(pnf_auxvars(ghosted_id), &
    !                        global_auxvars(ghosted_id), &
    !                        material_auxvars(ghosted_id), &
    !                        option,accum_p(local_id:local_id), &
    !                        Jdum,PETSC_FALSE)
  enddo

  call VecRestoreArrayReadF90(field%flow_xx,xx_p, ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum, accum_p, ierr);CHKERRQ(ierr)

end subroutine PNFUpdateFixedAccum

! ************************************************************************** !

subroutine PNFResidual(snes,xx,r,A,realization,ierr)
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
  type(PNF_parameter_type), pointer :: PNF_parameter
  type(pnf_auxvar_type), pointer :: pnf_auxvars(:)
  type(pnf_auxvar_type), pointer :: pnf_auxvars_bc(:)
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
  PetscReal :: v_darcy(1)

  discretization => realization%discretization
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  pnf_auxvars => patch%aux%PNF%auxvars
  pnf_auxvars_bc => patch%aux%PNF%auxvars_bc
  PNF_parameter => patch%aux%PNF%PNF_parameter
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars

  if (PNF_simult_function_evals) then
    call MatZeroEntries(A,ierr);CHKERRQ(ierr)
  endif

  ! Communication -----------------------------------------
  ! must be called before PNFUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  call PNFUpdateAuxVars(realization)

  ! override flags since they will soon be out of date
  patch%aux%PNF%auxvars_up_to_date = PETSC_FALSE

  if (PNF_numerical_derivatives) then
    !     ! Perturb aux vars
    ! do ghosted_id = 1, grid%ngmax  ! For each local node do...
    !   if (patch%imat(ghosted_id) <= 0) cycle
    !   call PNFAuxVarPerturb(pnf_auxvars(:,ghosted_id), &
    !                           global_auxvars(ghosted_id), &
    !                           material_auxvars(ghosted_id), &
    !                           patch%characteristic_curves_array( &
    !                             patch%cc_id(ghosted_id))%ptr, &
    !                           grid%nG2A(ghosted_id),option)
    ! enddo
  endif

  if (option%compute_mass_balance_new) then
    call PNFZeroMassBalanceDelta(realization)
  endif

  option%iflag = 1
  ! now assign access pointer to local variables
  call VecGetArrayF90(r, r_p, ierr);CHKERRQ(ierr)

  ! Accumulation terms ------------------------------------
  ! accumulation at t(k) (doesn't change during Newton iteration)
  if (PNF_calc_accum) then
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
      ! call PNFAccumulation(pnf_auxvars(ghosted_id), &
      !                        global_auxvars(ghosted_id), &
      !                        material_auxvars(ghosted_id), &
      !                        option,Res,Jup, &
      !                        PNF_simult_function_evals)
      ! if (PNF_numerical_derivatives) then
      !   call PNFAccumDerivative(pnf_auxvars(:,ghosted_id), &
      !                             global_auxvars(ghosted_id), &
      !                             material_auxvars(ghosted_id), &
      !                             option,Jup)
      ! endif
      r_p(local_id) =  r_p(local_id) + Res(1)
      accum_p2(local_id) = Res(1)
      if (PNF_simult_function_evals) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                      ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
    enddo
    call VecRestoreArrayF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)
  else
    r_p = 0.d0
  endif

  if (PNF_calc_flux) then
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

        ! call XXFlux(pnf_auxvars(ghosted_id_up), &
        !             global_auxvars(ghosted_id_up), &
        !             material_auxvars(ghosted_id_up), &
        !             pnf_auxvars(ghosted_id_dn), &
        !             global_auxvars(ghosted_id_dn), &
        !             material_auxvars(ghosted_id_dn), &
        !             cur_connection_set%area(iconn), &
        !             cur_connection_set%dist(:,iconn), &
        !             PNF_parameter,option,v_darcy, &
        !             Res,Jup,Jdn, &
        !             PNF_simult_function_evals)
        ! if (PNF_numerical_derivatives) then
        !   call XXFluxDerivative(pnf_auxvars(:,ghosted_id_up), &
        !                         global_auxvars(ghosted_id_up), &
        !                         material_auxvars(ghosted_id_up), &
        !                         pnf_auxvars(:,ghosted_id_dn), &
        !                         global_auxvars(ghosted_id_dn), &
        !                         material_auxvars(ghosted_id_dn), &
        !                         cur_connection_set%area(iconn), &
        !                         cur_connection_set%dist(:,iconn), &
        !                         PNF_parameter,option, &
        !                         Jup,Jdn)
        ! endif
        patch%internal_velocities(:,sum_connection) = v_darcy
        if (associated(patch%internal_flow_fluxes)) then
          patch%internal_flow_fluxes(:,sum_connection) = Res(:)
        endif

        if (local_id_up > 0) then
          r_p(local_id_up) = r_p(local_id_up) + Res(1)
          if (PNF_simult_function_evals) then
            call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                          Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
            call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                          Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
          endif
        endif

        if (local_id_dn > 0) then
          r_p(local_id_dn) = r_p(local_id_dn) - Res(1)
          if (PNF_simult_function_evals) then
            Jup = -Jup
            Jdn = -Jdn
            call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                                          Jdn,ADD_VALUES,ierr);CHKERRQ(ierr)
            call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                                          Jup,ADD_VALUES,ierr);CHKERRQ(ierr)
          endif
        endif
      enddo

      cur_connection_set => cur_connection_set%next
    enddo
  endif

  if (PNF_calc_bcflux) then
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

        ! call XXBCFlux(boundary_condition%flow_bc_type, &
        !               boundary_condition%flow_aux_mapping, &
        !               boundary_condition%flow_aux_real_var(:,iconn), &
        !               pnf_auxvars_bc(sum_connection), &
        !               global_auxvars_bc(sum_connection), &
        !               pnf_auxvars(ghosted_id), &
        !               global_auxvars(ghosted_id), &
        !               material_auxvars(ghosted_id), &
        !               cur_connection_set%area(iconn), &
        !               cur_connection_set%dist(:,iconn), &
        !               PNF_parameter,option, &
        !               v_darcy,Res,Jdn, &
        !               PNF_simult_function_evals)
        ! if (PNF_numerical_derivatives) then
        !   call XXBCFluxDerivative(boundary_condition%flow_bc_type, &
        !                           boundary_condition%flow_aux_mapping, &
        !                           boundary_condition% &
        !                             flow_aux_real_var(:,iconn), &
        !                           pnf_auxvars_bc(sum_connection), &
        !                           global_auxvars_bc(sum_connection), &
        !                           pnf_auxvars(:,ghosted_id), &
        !                           global_auxvars(ghosted_id), &
        !                           material_auxvars(ghosted_id), &
        !                           cur_connection_set%area(iconn), &
        !                           cur_connection_set%dist(:,iconn), &
        !                           PNF_parameter,option, &
        !                           Jdn)
        ! endif
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
        if (PNF_simult_function_evals) then
          Jdn = -Jdn
          call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn, &
                                        ADD_VALUES,ierr);CHKERRQ(ierr)
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

      ! call PNFSrcSink(option,source_sink%flow_condition%rate% &
      !                             dataset%rarray(:), &
      !                   source_sink%flow_condition%rate%itype, &
      !                   pnf_auxvars(ghosted_id), &
      !                   global_auxvars(ghosted_id), &
      !                   material_auxvars(ghosted_id), &
      !                   ss_flow_vol_flux, &
      !                   scale,Res,Jdn, &
      !                   PNF_simult_function_evals)
      ! if (PNF_numerical_derivatives) then
      !   call PNFSrcSinkDerivative(option, &
      !                               source_sink%flow_condition%rate% &
      !                                 dataset%rarray(:), &
      !                               source_sink%flow_condition%rate%itype, &
      !                               pnf_auxvars(:,ghosted_id), &
      !                               global_auxvars(ghosted_id), &
      !                               material_auxvars(ghosted_id), &
      !                               scale,Jdn)
      ! endif
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
      if (PNF_simult_function_evals) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn, &
                                      ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
    enddo
    source_sink => source_sink%next
  enddo

  if (patch%aux%PNF%inactive_cells_exist) then
    do i=1,patch%aux%PNF%matrix_zeroing%n_inactive_rows
      r_p(patch%aux%PNF%matrix_zeroing%inactive_rows_local(i)) = 0.d0
    enddo
  endif

  call VecRestoreArrayF90(r, r_p, ierr);CHKERRQ(ierr)

  if (PNF_simult_function_evals) then

    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      ! zero out inactive cells

    if (patch%aux%PNF%inactive_cells_exist) then
      scale = 1.d0 ! solely a temporary variable in this conditional
      call MatZeroRowsLocal(A,patch%aux%PNF%matrix_zeroing%n_inactive_rows, &
                            patch%aux%PNF%matrix_zeroing% &
                              inactive_rows_local_ghosted, &
                            scale,PETSC_NULL_VEC,PETSC_NULL_VEC, &
                            ierr);CHKERRQ(ierr)
    endif

    if (realization%debug%matview_Matrix) then
      call DebugWriteFilename(realization%debug,string,'PNFjacobian','', &
                              PNF_ts_count,PNF_ts_cut_count, &
                              PNF_ni_count)
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
    call DebugWriteFilename(realization%debug,string,'PNFresidual','', &
                            PNF_ts_count,PNF_ts_cut_count, &
                            PNF_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif
  if (realization%debug%vecview_solution) then
    call DebugWriteFilename(realization%debug,string,'PNFxx','', &
                            PNF_ts_count,PNF_ts_cut_count, &
                            PNF_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif

end subroutine PNFResidual

! ************************************************************************** !

subroutine PNFSetPlotVariables(realization,list)
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

end subroutine PNFSetPlotVariables

! ************************************************************************** !

subroutine PNFDestroy(realization)
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

end subroutine PNFDestroy

end module PNF_module
