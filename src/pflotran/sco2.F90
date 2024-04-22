module SCO2_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use SCO2_Aux_module
  use SCO2_Common_module
  use Global_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  public :: SCO2Setup, &
            SCO2InitializeTimestep, &
            SCO2UpdateSolution, &
            SCO2TimeCut,&
            SCO2UpdateAuxVars, &
            SCO2UpdateFixedAccum, &
            SCO2ComputeMassBalance, &
            SCO2ComputeComponentMassBalance, &
            SCO2Residual, &
            SCO2Jacobian, &
            SCO2SetPlotVariables, &
            SCO2MapBCAuxVarsToGlobal, &
            SCO2Destroy
contains

! ************************************************************************** !

subroutine SCO2Setup(realization)
  !
  ! Creates arrays for auxiliary variables
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
  !

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Grid_module
  use Fluid_module
  use Material_Aux_module
  use Output_Aux_module
  use Matrix_Zeroing_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(output_variable_list_type), pointer :: list
  type(material_parameter_type), pointer :: material_parameter

  PetscInt :: ghosted_id, iconn, sum_connection, local_id
  PetscInt :: lid, gid, sid, wid, co2_id
  PetscInt :: i, idof, ndof
  PetscBool :: error_found
  PetscInt :: flag(10)
  PetscBool, allocatable :: dof_is_active(:)
  PetscErrorCode :: ierr

  type(sco2_auxvar_type), pointer :: sco2_auxvars(:,:)
  type(sco2_auxvar_type), pointer :: sco2_auxvars_bc(:)
  type(sco2_auxvar_type), pointer :: sco2_auxvars_ss(:,:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(fluid_property_type), pointer :: cur_fluid_property

  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%aux%SCO2 => SCO2AuxCreate(option)

  material_parameter => patch%aux%Material%material_parameter
  error_found = PETSC_FALSE

  lid = option%liquid_phase
  gid = option%gas_phase
  sid = option%salt_id
  wid = option%water_id
  co2_id = option%co2_id

  if (minval(material_parameter%soil_heat_capacity(:)) < 0.d0) then
    option%io_buffer = 'ERROR: Non-initialized soil heat capacity.'
    call PrintMsgByRank(option)
    error_found = PETSC_TRUE
  endif
  if ( .not. associated(realization%characteristic_curves_thermal)) then
    if (minval(material_parameter%soil_thermal_conductivity(:,:)) < 0.d0)then
      option%io_buffer = 'ERROR: Non-initialized soil thermal conductivity.'
      call PrintMsg(option)
      error_found = PETSC_TRUE
    endif
  endif

  material_auxvars => patch%aux%Material%auxvars
  flag = 0
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
                     MPI_LOR,option%mycomm,ierr);CHKERRQ(ierr)
  if (error_found) then
    option%io_buffer = 'Material property errors found in SCO2Setup.'
    call PrintErrMsg(option)
  endif

  ! allocate auxvar data structures for all grid cells
  ndof = option%nflowdof
  if (sco2_central_diff_jacobian) then
    allocate(sco2_auxvars(0:2*ndof,grid%ngmax))
  else
    allocate(sco2_auxvars(0:ndof,grid%ngmax))
  endif
  do ghosted_id = 1, grid%ngmax
    if (sco2_central_diff_jacobian) then
      do idof = 0, 2*ndof
        call SCO2AuxVarInit(sco2_auxvars(idof,ghosted_id),option)
      enddo
    else
      do idof = 0, ndof
        call SCO2AuxVarInit(sco2_auxvars(idof,ghosted_id),option)
      enddo
    endif
  enddo
  patch%aux%SCO2%auxvars => sco2_auxvars
  patch%aux%SCO2%num_aux = grid%ngmax

  ! count the number of boundary connections and allocate
  ! auxvar data structures for them
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_condition_list)
  if (sum_connection > 0) then
    allocate(sco2_auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call SCO2AuxVarInit(sco2_auxvars_bc(iconn),option)
    enddo
    patch%aux%SCO2%auxvars_bc => sco2_auxvars_bc
  endif
  patch%aux%SCO2%num_aux_bc = sum_connection

  ! count the number of source/sink connections and allocate
  ! auxvar data structures for them
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (sum_connection > 0) then
    allocate(sco2_auxvars_ss(0:ONE_INTEGER,sum_connection))
    do iconn = 1, sum_connection
      do i = 0,ONE_INTEGER
        ! Index 0 contains user-provided information.
        ! Index 1 contains state variables used for perturbed values
        ! in the source/sink term
        call SCO2AuxVarInit(sco2_auxvars_ss(i,iconn),option)
      enddo
    enddo
    patch%aux%SCO2%auxvars_ss => sco2_auxvars_ss
  endif
  patch%aux%SCO2%num_aux_ss = sum_connection


  ! initialize parameters
  ! MAN: Need to check this holds diffusion coefficient for salt.
  !      Should index by component, not phase?
  cur_fluid_property => realization%fluid_properties
  do
    if (.not.associated(cur_fluid_property)) exit
    if (cur_fluid_property%phase_id == LIQUID_PHASE) then
      patch%aux%SCO2%sco2_parameter% &
      diffusion_coefficient(co2_id,lid) = &
        cur_fluid_property%diffusion_coefficient
      patch%aux%SCO2%sco2_parameter% &
      diffusion_coefficient(sid,lid) = &
        cur_fluid_property%salt_diffusion_coefficient
      ! Not used
      patch%aux%SCO2%sco2_parameter% &
      diffusion_coefficient(wid,lid) = 0.d0
    elseif (cur_fluid_property%phase_id == GAS_PHASE) then
      patch%aux%SCO2%sco2_parameter% &
      diffusion_coefficient(:,gid) = &
        cur_fluid_property%gas_diffusion_coefficient
      patch%aux%SCO2%sco2_parameter% &
      diffusion_coefficient(sid,gid) = 0.d0
    endif
    cur_fluid_property => cur_fluid_property%next
  enddo
  ! check whether diffusion coefficients are initialized.
  if (Uninitialized(patch%aux%SCO2%sco2_parameter% &
      diffusion_coefficient(co2_id,lid))) then
    option%io_buffer = &
      UninitializedMessage('Liquid phase diffusion coefficient','')
    call PrintErrMsg(option)
  endif
  if (Uninitialized(patch%aux%SCO2%sco2_parameter% &
      diffusion_coefficient(co2_id,gid))) then
    option%io_buffer = &
      UninitializedMessage('Gas phase diffusion coefficient','')
    call PrintErrMsg(option)
  endif

  list => realization%output_option%output_snap_variable_list
  call SCO2SetPlotVariables(realization,list)
  list => realization%output_option%output_obs_variable_list
  call SCO2SetPlotVariables(realization,list)

  allocate(dof_is_active(option%nflowdof))
  dof_is_active = PETSC_TRUE
  if (option%coupled_well) then
    dof_is_active(option%nflowdof) = PETSC_FALSE
  endif
  call PatchCreateZeroArray(patch,dof_is_active, &
                            patch%aux%SCO2%matrix_zeroing, &
                            patch%aux%SCO2%inactive_cells_exist,option)
  deallocate(dof_is_active)

  call PatchSetupUpwindDirection(patch,option)

  sco2_ts_count = 0
  sco2_ts_cut_count = 0
  sco2_ni_count = 0

end subroutine SCO2Setup

! ************************************************************************** !

subroutine SCO2InitializeTimestep(realization)
  !
  ! Update data in module prior to time step
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
  !

  use Realization_Subsurface_class
  use Upwind_Direction_module

  implicit none

  class(realization_subsurface_type) :: realization

  if (sco2_restrict_state_chng) then
    realization%patch%aux%SCO2%auxvars(:,:)%istatechng = PETSC_FALSE
  endif

  sco2_newton_iteration_number = -999
  sco2_sub_newton_iter_num = 0
  update_upwind_direction = PETSC_TRUE
  call SCO2UpdateFixedAccum(realization)

  sco2_ni_count = 0

end subroutine SCO2InitializeTimestep

! ************************************************************************** !

subroutine SCO2UpdateSolution(realization)
  !
  ! Updates data in module after a successful time
  ! step
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
  !

  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module
  use Grid_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(sco2_auxvar_type), pointer :: sco2_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  PetscInt :: ghosted_id, lid, tgid

  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  sco2_auxvars => patch%aux%SCO2%auxvars
  global_auxvars => patch%aux%Global%auxvars

  lid = option%liquid_phase
  tgid = option%trapped_gas_phase

  if (realization%option%compute_mass_balance_new) then
    call SCO2UpdateMassBalance(realization)
  endif

  ! update stored state
  do ghosted_id = 1, grid%ngmax
    sco2_auxvars(ZERO_INTEGER,ghosted_id)%istate_store(PREV_TS) = &
      global_auxvars(ghosted_id)%istate
    if (sco2_auxvars(ZERO_INTEGER,ghosted_id)%sat(lid) < &
        sco2_auxvars(ZERO_INTEGER,ghosted_id)%sl_min) then

      sco2_auxvars(ZERO_INTEGER,ghosted_id)%sl_min = &
      sco2_auxvars(ZERO_INTEGER,ghosted_id)%sat(lid)

      if (global_auxvars(ghosted_id)%istate == SCO2_TRAPPED_GAS_STATE) then
        sco2_auxvars(ZERO_INTEGER,ghosted_id)%sg_trapped = &
                      sco2_auxvars(ZERO_INTEGER,ghosted_id)%sat(tgid)
      else
        sco2_auxvars(ZERO_INTEGER,ghosted_id)%sg_trapped = 0.d0
      endif

    endif

    global_auxvars(ghosted_id)%sat(1:option%nphase) = &
                                sco2_auxvars(ZERO_INTEGER,ghosted_id)% &
                                sat(1:option%nphase)
    global_auxvars(ghosted_id)%den(1:option%nphase) = &
                                sco2_auxvars(ZERO_INTEGER,ghosted_id)% &
                                den(1:option%nphase)
    global_auxvars(ghosted_id)%den_kg(1:option%nphase) = &
                                sco2_auxvars(ZERO_INTEGER,ghosted_id)% &
                                den_kg(1:option%nphase)

  enddo

  sco2_ts_count = sco2_ts_count + 1
  sco2_ts_cut_count = 0
  sco2_ni_count = 0

end subroutine SCO2UpdateSolution

! ************************************************************************** !

subroutine SCO2TimeCut(realization)
  !
  ! Resets arrays for time step cut
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
  !
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Patch_module
  use Discretization_module
  use Grid_module

  implicit none

  class(realization_subsurface_type) :: realization
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(sco2_auxvar_type), pointer :: sco2_auxvars(:,:)

  PetscInt :: ghosted_id

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  global_auxvars => patch%aux%Global%auxvars
  sco2_auxvars => patch%aux%SCO2%auxvars

  ! restore stored state
  do ghosted_id = 1, grid%ngmax
    global_auxvars(ghosted_id)%istate = &
      sco2_auxvars(ZERO_INTEGER,ghosted_id)%istate_store(PREV_TS)
  enddo

  sco2_ts_cut_count = sco2_ts_cut_count + 1

  call SCO2InitializeTimestep(realization)

end subroutine SCO2TimeCut

! ************************************************************************** !

subroutine SCO2ComputeMassBalance(realization,mass_balance,mass_trapped)
  !
  ! Initializes mass balance
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
  !

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
  use Material_Aux_module

  implicit none

  class(realization_subsurface_type) :: realization
  PetscReal :: mass_balance(realization%option%nflowspec, &
                            realization%option%nphase)
  PetscReal :: mass_trapped(realization%option%nphase)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(sco2_auxvar_type), pointer :: sco2_auxvars(:,:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iphase, icomp
  PetscReal :: vol_phase
  PetscInt :: tgid, gid, co2_id

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  sco2_auxvars => patch%aux%SCO2%auxvars
  material_auxvars => patch%aux%Material%auxvars

  tgid = option%trapped_gas_phase
  gid = option%gas_phase
  co2_id = option%co2_id

  mass_balance = 0.d0
  mass_trapped = 0.d0

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    do iphase = 1, option%nphase
      ! volume_phase = saturation*porosity*volume
      vol_phase = &
        sco2_auxvars(ZERO_INTEGER,ghosted_id)%sat(iphase)* &
        sco2_auxvars(ZERO_INTEGER,ghosted_id)%effective_porosity* &
        material_auxvars(ghosted_id)%volume
      ! mass = volume_phase*density
      do icomp = 1, option%nflowspec
        mass_balance(icomp,iphase) = mass_balance(icomp,iphase) + &
          sco2_auxvars(ZERO_INTEGER,ghosted_id)%den_kg(iphase)* &
          sco2_auxvars(ZERO_INTEGER,ghosted_id)%xmass(icomp,iphase) * &
          vol_phase
      enddo
    enddo
    ! Trapped gas mass
    vol_phase = &
        sco2_auxvars(ZERO_INTEGER,ghosted_id)%sat(tgid)* &
        sco2_auxvars(ZERO_INTEGER,ghosted_id)%effective_porosity* &
        material_auxvars(ghosted_id)%volume
    mass_trapped(co2_id) = mass_trapped(co2_id) + &
        sco2_auxvars(ZERO_INTEGER,ghosted_id)%den_kg(gid)* &
        sco2_auxvars(ZERO_INTEGER,ghosted_id)%xmass(co2_id,gid) * &
        vol_phase
  enddo

end subroutine SCO2ComputeMassBalance

! ************************************************************************** !

subroutine SCO2ComputeComponentMassBalance(realization,num_cells,num_comp, &
                                           num_phase,sum_kg,cell_ids)
  !
  ! Author: Michael Nole
  ! Date: 02/29/2024
  !

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
  use Material_Aux_module


  class(realization_subsurface_type) :: realization
  PetscInt :: num_cells
  PetscInt :: num_comp
  PetscInt :: num_phase
  PetscReal :: sum_kg(num_comp,num_phase)
  PetscInt, pointer :: cell_ids(:)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(sco2_auxvar_type), pointer :: sco2_auxvars(:,:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: icomp, iphase, k
  PetscReal :: porosity, volume

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  sco2_auxvars => patch%aux%SCO2%auxvars
  material_auxvars => patch%aux%Material%auxvars

  sum_kg = 0.d0

  do k=1, num_cells
    local_id = cell_ids(k)
    ghosted_id = grid%nL2G(local_id)

    if (patch%imat(ghosted_id) <= 0) cycle
    volume = material_auxvars(ghosted_id)%volume


    do iphase = 1,num_phase
      do icomp = 1,num_comp
        if (iphase /=3) then
          porosity = sco2_auxvars(ZERO_INTEGER,ghosted_id)%effective_porosity
        else
          porosity = material_auxvars(ghosted_id)%porosity
        endif
        sum_kg(icomp,iphase) = sum_kg(icomp,iphase) + &
                    sco2_auxvars(ZERO_INTEGER,ghosted_id)%xmass(icomp,iphase) * &
                    sco2_auxvars(ZERO_INTEGER,ghosted_id)%den_kg(iphase) * &
                    sco2_auxvars(ZERO_INTEGER,ghosted_id)%sat(iphase) * &
                    porosity * volume
      enddo
    enddo
  enddo

end subroutine SCO2ComputeComponentMassBalance
! ************************************************************************** !

subroutine SCO2ZeroMassBalanceDelta(realization)
  !
  ! Zeros mass balance delta array
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
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

  do iconn = 1, patch%aux%SCO2%num_aux_bc
    global_auxvars_bc(iconn)%mass_balance_delta = 0.d0
  enddo
  do iconn = 1, patch%aux%SCO2%num_aux_ss
    global_auxvars_ss(iconn)%mass_balance_delta = 0.d0
  enddo

end subroutine SCO2ZeroMassBalanceDelta

! ************************************************************************** !

subroutine SCO2UpdateMassBalance(realization)
  !
  ! Updates mass balance
  !
  ! Author: Glenn Hammond
  ! Date: 01/26/24
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

  do iconn = 1, patch%aux%SCO2%num_aux_bc
    do icomp = 1, option%nflowspec
      global_auxvars_bc(iconn)%mass_balance(icomp,:) = &
        global_auxvars_bc(iconn)%mass_balance(icomp,:) + &
        global_auxvars_bc(iconn)%mass_balance_delta(icomp,:)* &
        fmw_comp(icomp)*option%flow_dt
    enddo
  enddo
  do iconn = 1, patch%aux%SCO2%num_aux_ss
    do icomp = 1, option%nflowspec
      global_auxvars_ss(iconn)%mass_balance(icomp,:) = &
        global_auxvars_ss(iconn)%mass_balance(icomp,:) + &
        global_auxvars_ss(iconn)%mass_balance_delta(icomp,:)* &
        fmw_comp(icomp)*option%flow_dt
    enddo
  enddo

end subroutine SCO2UpdateMassBalance

! ************************************************************************** !

subroutine SCO2UpdateAuxVars(realization,update_state,update_state_bc)
  !
  ! Updates the SCO2 mode auxiliary variables.
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
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
  use EOS_Water_module
  use Saturation_Function_module

  implicit none

  class(realization_subsurface_type) :: realization
  PetscBool :: update_state
  PetscBool :: update_state_bc

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition, source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(sco2_auxvar_type), pointer :: sco2_auxvars(:,:), sco2_auxvars_bc(:), &
                                     sco2_auxvars_ss(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:), &
                                       global_auxvars_bc(:), global_auxvars_ss(:)
  type(sco2_parameter_type), pointer :: sco2_parameter
  type(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: ghosted_id, local_id, sum_connection, idof, iconn, natural_id
  PetscInt :: ghosted_start, ghosted_end
  PetscInt :: offset
  PetscInt :: istate

  PetscInt :: wid, co2_id, sid, lid, gid, co2_pressure_id
  PetscReal :: qsrc(realization%option%nflowdof)
  PetscInt :: real_index, flow_src_sink_type
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)

  PetscReal :: cell_pressure, scale

  PetscReal :: Res_dummy(realization%option%nflowdof)
  PetscErrorCode :: ierr

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  sco2_auxvars => patch%aux%SCO2%auxvars
  sco2_auxvars_bc => patch%aux%SCO2%auxvars_bc
  sco2_auxvars_ss => patch%aux%SCO2%auxvars_ss
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars
  sco2_parameter => patch%aux%SCO2%sco2_parameter

  wid = option%water_id
  co2_id = option%co2_id
  sid = option%salt_id
  lid = option%liquid_phase
  gid = option%gas_phase
  co2_pressure_id = option%co2_pressure_id

  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    if (patch%imat(ghosted_id) <= 0) cycle
    ghosted_end = ghosted_id*option%nflowdof
    ghosted_start = ghosted_end - option%nflowdof + 1
    ! SCO2_UPDATE_FOR_ACCUM indicates call from non-perturbation
    option%iflag = SCO2_UPDATE_FOR_ACCUM
    natural_id = grid%nG2A(ghosted_id)
    if (sco2_well_coupling == SCO2_FULLY_IMPLICIT_WELL) then
      ! MAN: this is a hack to get well to initialize properly
      if (xx_loc_p(ghosted_end) /= &
          sco2_auxvars(ZERO_INTEGER,ghosted_id)%pres(gid)) then
        sco2_auxvars(ZERO_INTEGER,ghosted_id)%well%bh_p = xx_loc_p(ghosted_end)
      endif
    endif
    if (grid%nG2L(ghosted_id) == 0) natural_id = -natural_id
    call SCO2AuxVarCompute(xx_loc_p(ghosted_start:ghosted_end), &
                           sco2_auxvars(ZERO_INTEGER,ghosted_id), &
                           global_auxvars(ghosted_id), &
                           material_auxvars(ghosted_id), &
                           patch%characteristic_curves_array( &
                           patch%cc_id(ghosted_id))%ptr, &
                           sco2_parameter, natural_id, &
                           option)

    if (update_state) then
      call SCO2AuxVarUpdateState(xx_loc_p(ghosted_start:ghosted_end), &
                                 sco2_auxvars(ZERO_INTEGER,ghosted_id), &
                                 global_auxvars(ghosted_id), &
                                 material_auxvars(ghosted_id), &
                                 patch%characteristic_curves_array( &
                                 patch%cc_id(ghosted_id))%ptr, &
                                 sco2_parameter,natural_id,option)

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
      !Negate to indicate boundary connection, not actual cell
      natural_id = -grid%nG2A(ghosted_id)
      offset = (ghosted_id-1)*option%nflowdof
      if (patch%imat(ghosted_id) <= 0) cycle

      xxbc(:) = xx_loc_p(offset+1:offset+option%nflowdof)
      istate = boundary_condition%flow_aux_int_var(ONE_INTEGER,iconn)
      if (istate == SCO2_ANY_STATE) then
        istate = global_auxvars(ghosted_id)%istate
        select case(istate)
          case(SCO2_LIQUID_STATE,SCO2_GAS_STATE,SCO2_TRAPPED_GAS_STATE)
            do idof = 1, option%nflowdof
              select case(boundary_condition%flow_bc_type(idof))
                case(DIRICHLET_BC,HYDROSTATIC_BC)
                  real_index = boundary_condition% &
                     flow_aux_mapping(dof_to_primary_variable(idof,istate))
                  xxbc(idof) = boundary_condition% &
                     flow_aux_real_var(real_index,iconn)
              end select
            enddo
          case(SCO2_LIQUID_GAS_STATE)
            do idof = 1, option%nflowdof
              select case(boundary_condition%flow_bc_type(idof))
                case(HYDROSTATIC_BC)
                  real_index = boundary_condition% &
                        flow_aux_mapping(dof_to_primary_variable(idof,istate))
                  xxbc(idof) = boundary_condition% &
                        flow_aux_real_var(real_index,iconn)
                case(DIRICHLET_BC)
                    option%io_buffer = 'Mixed FLOW_CONDITION "' // &
                          trim(boundary_condition%flow_condition%name) // &
                          '" not fully supported yet for SCO2 Mode.'
                    call PrintErrMsg(option)
                case(NEUMANN_BC)
                case default
                  option%io_buffer = 'Unknown BC type in SCO2UpdateAuxVars().'
                  call PrintErrMsg(option)
              end select
            enddo
        end select
      else
        ! we do this for all BCs; Neumann bcs will be set later
        do idof = 1, option%nflowdof
          real_index = boundary_condition%flow_aux_mapping(&
                       dof_to_primary_variable(idof,istate))
          if (real_index > 0) then
            xxbc(idof) = boundary_condition%flow_aux_real_var(real_index,iconn)
          else
            option%io_buffer = &
                     'Error setting up boundary condition in SCO2UpdateAuxVars'
            call PrintErrMsg(option)
          endif
        enddo
      endif

      ! set this based on data given
      global_auxvars_bc(sum_connection)%istate = istate

      ! SCO2_UPDATE_FOR_BOUNDARY indicates call from non-perturbation
      option%iflag = SCO2_UPDATE_FOR_BOUNDARY
      call SCO2AuxVarCompute(xxbc,sco2_auxvars_bc(sum_connection), &
                             global_auxvars_bc(sum_connection), &
                             material_auxvars(ghosted_id), &
                             patch%characteristic_curves_array( &
                             patch%cc_id(ghosted_id))%ptr, &
                             sco2_parameter,natural_id, option)
      if (update_state_bc) then
        ! update state and update aux var; this could result in two updates to
        ! the aux var as update state updates if the state changes
        call SCO2AuxVarUpdateState(xxbc,sco2_auxvars_bc(sum_connection), &
                                   global_auxvars_bc(sum_connection), &
                                   material_auxvars(ghosted_id), &
                                   patch%characteristic_curves_array( &
                                   patch%cc_id(ghosted_id))%ptr, &
                                   sco2_parameter,natural_id,option)
        if (.not. any(boundary_condition%flow_bc_type == NEUMANN_BC)) then
            boundary_condition%flow_aux_int_var(ONE_INTEGER,iconn) = &
                          global_auxvars_bc(sum_connection)%istate
            select case(boundary_condition%flow_aux_int_var(ONE_INTEGER,iconn))
              case(SCO2_LIQUID_STATE)
                boundary_condition%flow_aux_real_var( &
                               SCO2_LIQUID_PRESSURE_DOF,iconn) = &
                               sco2_auxvars_bc(sum_connection)%pres(lid)
                boundary_condition%flow_aux_real_var( &
                               SCO2_CO2_MASS_FRAC_DOF,iconn) = &
                               sco2_auxvars_bc(sum_connection)%xmass(co2_id,lid)
                boundary_condition%flow_aux_real_var( &
                               SCO2_SALT_MASS_FRAC_DOF,iconn) = &
                               sco2_auxvars_bc(sum_connection)%m_salt(1)
              case(SCO2_GAS_STATE)
                boundary_condition%flow_aux_real_var( &
                               SCO2_GAS_PRESSURE_DOF,iconn) = &
                               sco2_auxvars_bc(sum_connection)%pres(gid)
                boundary_condition%flow_aux_real_var( &
                               SCO2_CO2_PRESSURE_DOF,iconn) = &
                               sco2_auxvars_bc(sum_connection)% &
                               pres(co2_pressure_id)
                boundary_condition%flow_aux_real_var( &
                               SCO2_SALT_MASS_FRAC_DOF,iconn) = &
                               sco2_auxvars_bc(sum_connection)%m_salt(2)
              case(SCO2_TRAPPED_GAS_STATE)
                boundary_condition%flow_aux_real_var( &
                               SCO2_LIQUID_PRESSURE_DOF,iconn) = &
                               sco2_auxvars_bc(sum_connection)%pres(lid)
                boundary_condition%flow_aux_real_var( &
                               SCO2_GAS_SATURATION_DOF,iconn) = &
                               sco2_auxvars_bc(sum_connection)%sat(gid)
                boundary_condition%flow_aux_real_var( &
                               SCO2_SALT_MASS_FRAC_DOF,iconn) = &
                               sco2_auxvars_bc(sum_connection)%m_salt(1)
              case(SCO2_LIQUID_GAS_STATE)
                boundary_condition%flow_aux_real_var( &
                               SCO2_LIQUID_PRESSURE_DOF,iconn) = &
                               sco2_auxvars_bc(sum_connection)%pres(lid)
                boundary_condition%flow_aux_real_var( &
                               SCO2_GAS_PRESSURE_DOF,iconn) = &
                               sco2_auxvars_bc(sum_connection)%pres(gid)
                boundary_condition%flow_aux_real_var( &
                               SCO2_SALT_MASS_FRAC_DOF,iconn) = &
                               sco2_auxvars_bc(sum_connection)%m_salt(1)
            end select
        endif
      endif
    enddo
    boundary_condition => boundary_condition%next
  enddo

  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit
    if (associated(source_sink%flow_condition%well)) then
      source_sink => source_sink%next
      cycle
    endif

    qsrc = source_sink%flow_condition%sco2%rate%dataset%rarray(:)
    cur_connection_set => source_sink%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0) cycle

      flow_src_sink_type = source_sink%flow_condition%sco2%rate%itype

      cell_pressure = maxval(sco2_auxvars(ZERO_INTEGER,ghosted_id)% &
                             pres(option%liquid_phase:option%gas_phase))

      ! Need to make sure all possible primary variables are covered before
      ! sending to AuxVarCompute

      if (associated(source_sink%flow_condition%sco2%temperature)) then
        sco2_auxvars_ss(ZERO_INTEGER,sum_connection)%temp = &
          source_sink%flow_condition%sco2%temperature%dataset%rarray(1)
      else
        sco2_auxvars_ss(ZERO_INTEGER,sum_connection)%temp = &
          sco2_auxvars(ZERO_INTEGER,ghosted_id)%temp
      endif

      ! Check if liquid pressure is set
      if (associated(source_sink%flow_condition%sco2%liquid_pressure)) then
        sco2_auxvars_ss(ZERO_INTEGER,sum_connection)%pres(lid) = &
          source_sink%flow_condition%sco2%liquid_pressure%dataset%rarray(1)
      else
        sco2_auxvars_ss(ZERO_INTEGER,sum_connection)%pres(lid) = &
          sco2_auxvars(ZERO_INTEGER,ghosted_id)%pres(option%liquid_phase)
      endif

      ! Check if gas pressure is set
      if (associated(source_sink%flow_condition%sco2%gas_pressure)) then
        sco2_auxvars_ss(ZERO_INTEGER,sum_connection)%pres(gid) = &
          source_sink%flow_condition%sco2%gas_pressure%dataset%rarray(1)
      else
        sco2_auxvars_ss(ZERO_INTEGER,sum_connection)%pres(gid) = &
                                                            cell_pressure
      endif

      ! Salt mass
      if (associated(source_sink%flow_condition%sco2%salt_mass)) then
        sco2_auxvars_ss(ZERO_INTEGER,sum_connection)%xmass(sid,lid) = &
          source_sink%flow_condition%sco2%salt_mass%dataset%rarray(1)
      else
        sco2_auxvars_ss(ZERO_INTEGER,sum_connection)%xmass(sid,lid) = 0.d0
      endif

      sco2_auxvars_ss(ZERO_INTEGER,sum_connection)%m_salt(:) = &
        sco2_auxvars_ss(ZERO_INTEGER,sum_connection)%xmass(sid,lid)

      ! CO2 mass fraction
      if (associated(source_sink%flow_condition%sco2%co2_mass_fraction)) then
        sco2_auxvars_ss(ZERO_INTEGER,sum_connection)%xmass(co2_id,lid) = &
          source_sink%flow_condition%sco2%co2_mass_fraction%dataset%rarray(1)
      else
        sco2_auxvars_ss(ZERO_INTEGER,sum_connection)%xmass(co2_id,lid) = 0.d0
      endif

      ! CO2 partial pressure
      if (associated(source_sink%flow_condition%sco2%co2_pressure)) then
        sco2_auxvars_ss(ZERO_INTEGER,sum_connection)%pres(co2_pressure_id) = &
          source_sink%flow_condition%sco2%co2_pressure%dataset%rarray(1)
      else
        sco2_auxvars_ss(ZERO_INTEGER,sum_connection)%pres(co2_pressure_id) = &
          cell_pressure
      endif

      ! MAN: can we infer the state instead, or does this matter?
      !      This may get weird with trapped gas state, need to
      !      check the min liq sat applied to flow conditions...
      if (dabs(qsrc(wid)) > 0.d0 .and. &
          dabs(qsrc(co2_id)) > 0.d0) then
        global_auxvars_ss(sum_connection)%istate = SCO2_LIQUID_GAS_STATE
      elseif (dabs(qsrc(wid)) > 0.d0) then
        global_auxvars_ss(sum_connection)%istate = SCO2_LIQUID_STATE
      elseif (dabs(qsrc(co2_id)) > 0.d0) then
        global_auxvars_ss(sum_connection)%istate = SCO2_GAS_STATE
      else
        global_auxvars_ss(sum_connection)%istate = SCO2_LIQUID_GAS_STATE
      endif

      !if (global_auxvars_ss(sum_connection)%istate /= &
      !    global_auxvars(ghosted_id)%istate) then
      !  global_auxvars_ss(sum_connection)%istate = SCO2_LIQUID_GAS_STATE
      !endif

      option%iflag = SCO2_UPDATE_FOR_SS

      ! Compute state variables
      call SCO2AuxVarComputeAndSrcSink(option,qsrc,flow_src_sink_type, &
                          sco2_auxvars_ss(ZERO_INTEGER,sum_connection), &
                          sco2_auxvars(ZERO_INTEGER,ghosted_id), &
                          global_auxvars(ghosted_id), &
                          global_auxvars_ss(sum_connection), &
                          material_auxvars(ghosted_id), &
                          patch%characteristic_curves_array( &
                          patch%cc_id(ghosted_id))%ptr, &
                          sco2_parameter, grid%nG2A(ghosted_id), &
                          scale, Res_dummy, PETSC_TRUE) ! aux_var_compute_only

    enddo
    source_sink => source_sink%next
  enddo
  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

  patch%aux%SCO2%auxvars_up_to_date = PETSC_TRUE

end subroutine SCO2UpdateAuxVars

! ************************************************************************** !

subroutine SCO2UpdateFixedAccum(realization)
  !
  ! Updates the fixed portion of the accumulation term
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
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
  type(sco2_auxvar_type), pointer :: sco2_auxvars(:,:)
  type(sco2_parameter_type), pointer :: sco2_parameter
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

  sco2_auxvars => patch%aux%SCO2%auxvars
  sco2_parameter => patch%aux%SCO2%sco2_parameter
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  material_parameter => patch%aux%Material%material_parameter

  call VecGetArrayReadF90(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    natural_id = grid%nG2A(ghosted_id)
    local_end = local_id*option%nflowdof
    local_start = local_end - option%nflowdof + 1
    ! SCO2_UPDATE_FOR_FIXED_ACCUM indicates call from non-perturbation
    option%iflag = SCO2_UPDATE_FOR_FIXED_ACCUM


    call SCO2AuxVarCompute(xx_p(local_start:local_end), &
                           sco2_auxvars(ZERO_INTEGER,ghosted_id), &
                           global_auxvars(ghosted_id), &
                           material_auxvars(ghosted_id), &
                           patch%characteristic_curves_array( &
                           patch%cc_id(ghosted_id))%ptr, &
                           sco2_parameter, natural_id, option)

    call SCO2Accumulation(sco2_auxvars(ZERO_INTEGER,ghosted_id), &
                             global_auxvars(ghosted_id), &
                             material_auxvars(ghosted_id), &
                             material_parameter%soil_heat_capacity(imat), &
                             option,accum_p(local_start:local_end))
  enddo


  call VecRestoreArrayReadF90(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)

end subroutine SCO2UpdateFixedAccum

! ************************************************************************** !

subroutine SCO2Residual(snes,xx,r,realization,pm_well,ierr)
  !
  ! Computes the residual
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
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
  use Material_Aux_module
  use Upwind_Direction_module
  use PM_Well_class

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  class(realization_subsurface_type) :: realization
  class(pm_well_type), pointer :: pm_well
  PetscViewer :: viewer
  PetscErrorCode :: ierr

  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(material_parameter_type), pointer :: material_parameter
  type(sco2_parameter_type), pointer :: sco2_parameter
  type(sco2_auxvar_type), pointer :: sco2_auxvars(:,:), sco2_auxvars_bc(:), &
                                        sco2_auxvars_ss(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
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
  PetscInt :: flow_src_sink_type
  PetscInt :: co2_id, sid, wid

  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum_p(:), accum_p2(:)

  PetscReal :: qsrc(realization%option%nflowdof)

  character(len=MAXSTRINGLENGTH) :: string

  PetscInt :: icct_up, icct_dn
  PetscReal :: Res(realization%option%nflowdof)
  PetscReal :: v_darcy(realization%option%nphase)

  discretization => realization%discretization
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  sco2_auxvars => patch%aux%SCO2%auxvars
  sco2_auxvars_bc => patch%aux%SCO2%auxvars_bc
  sco2_auxvars_ss => patch%aux%SCO2%auxvars_ss
  sco2_parameter => patch%aux%SCO2%sco2_parameter
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars

  wid = option%water_id
  co2_id = option%co2_id
  sid = option%salt_id

  !MAN: check if we need this
  if (sco2_newton_iteration_number > 1 .and. &
      mod(sco2_newton_iteration_number-1, &
          upwind_dir_update_freq) == 0) then
    update_upwind_direction = PETSC_TRUE
  endif

  ! Communication -----------------------------------------
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)

  ! update state
  sco2_allow_state_change = PETSC_TRUE
  sco2_state_changed = PETSC_FALSE

  if (sco2_sub_newton_iter_num > 0 .and. option%flow%using_newtontrdc .and. &
      sco2_newtontrdc_hold_inner) then
    ! when newtonTR is active and has inner iterations to re-evaluate the residual,
    ! primary variables must not change.
    sco2_allow_state_change = PETSC_FALSE
  endif
                                            ! do update state
  call SCO2UpdateAuxVars(realization,sco2_allow_state_change, &
                         sco2_allow_state_change)

  ! override flags since they will soon be out of date
  patch%aux%SCO2%auxvars_up_to_date = PETSC_FALSE

  ! always assume variables have been swapped; therefore, must copy back
  call VecLockPop(xx,ierr);CHKERRQ(ierr)
  call DiscretizationLocalToGlobal(discretization,field%flow_xx_loc,xx, &
                                   NFLOWDOF)
  call VecLockPush(xx,ierr);CHKERRQ(ierr)

  if (option%compute_mass_balance_new) then
    call SCO2ZeroMassBalanceDelta(realization)
  endif

  option%iflag = SCO2_UPDATE_FOR_ACCUM
  ! now assign access pointer to local variables
  call VecGetArrayF90(r,r_p,ierr);CHKERRQ(ierr)

  ! Accumulation terms ------------------------------------
  ! accumulation at t(k) (doesn't change during Newton iteration)
  call VecGetArrayReadF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)
  r_p = -accum_p
  call VecRestoreArrayReadF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)

  ! accumulation at t(k+1)
  call VecGetArrayF90(field%flow_accum2,accum_p2,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    local_end = local_id * option%nflowdof
    local_start = local_end - option%nflowdof + 1
    call SCO2Accumulation(sco2_auxvars(ZERO_INTEGER,ghosted_id), &
                             global_auxvars(ghosted_id), &
                             material_auxvars(ghosted_id), &
                             material_parameter%soil_heat_capacity(imat), &
                             option,Res)
    r_p(local_start:local_end) =  r_p(local_start:local_end) + Res(:)
    accum_p2(local_start:local_end) = Res(:)
  enddo
  if (sco2_well_coupling == SCO2_FULLY_IMPLICIT_WELL) then
    if (associated(pm_well)) then
      if (pm_well%well_grid%h_rank_id(1) == option%myrank) then
        if (dabs(pm_well%well%th_qg) > 0.d0) then
          accum_p2(local_end) = pm_well%well%th_qg
        elseif (dabs(pm_well%well%th_ql) > 0.d0) then
          accum_p2(local_end) = pm_well%well%th_ql
        else
          accum_p2(local_end) = 0.d0
        endif
        r_p(local_end) = 0.d0
      endif
    endif
  endif
  call VecRestoreArrayF90(field%flow_accum2, accum_p2, ierr);CHKERRQ(ierr)

  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0
  do
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      if (.not. Initialized(cur_connection_set%face_id(iconn))) cycle
      if (cur_connection_set%area(iconn) <= 0.d0) cycle
      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping

      imat_up = patch%imat(ghosted_id_up)
      imat_dn = patch%imat(ghosted_id_dn)
      if (imat_up <= 0 .or. imat_dn <= 0) cycle

      icct_up = patch%cct_id(ghosted_id_up)
      icct_dn = patch%cct_id(ghosted_id_dn)

      call SCO2Flux(sco2_auxvars(ZERO_INTEGER,ghosted_id_up), &
                       global_auxvars(ghosted_id_up), &
                       material_auxvars(ghosted_id_up), &
                       patch%char_curves_thermal_array(icct_up)%ptr, &
                       sco2_auxvars(ZERO_INTEGER,ghosted_id_dn), &
                       global_auxvars(ghosted_id_dn), &
                       material_auxvars(ghosted_id_dn), &
                       patch%char_curves_thermal_array(icct_dn)%ptr, &
                       cur_connection_set%area(iconn), &
                       cur_connection_set%dist(:,iconn), &
                       patch%flow_upwind_direction(:,iconn), &
                       option,v_darcy,Res, &
                       update_upwind_direction, &
                       count_upwind_direction_flip)

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

      icct_dn = patch%cct_id(ghosted_id)

      call SCO2BCFlux(boundary_condition%flow_bc_type, &
                     boundary_condition%flow_aux_mapping, &
                     boundary_condition%flow_aux_real_var(:,iconn), &
                     sco2_auxvars_bc(sum_connection), &
                     global_auxvars_bc(sum_connection), &
                     sco2_auxvars(ZERO_INTEGER,ghosted_id), &
                     global_auxvars(ghosted_id), &
                     material_auxvars(ghosted_id), &
                     patch%char_curves_thermal_array(icct_dn)%ptr, &
                     cur_connection_set%area(iconn), &
                     cur_connection_set%dist(:,iconn), &
                     patch%flow_upwind_direction_bc(:,iconn), &
                     option,v_darcy,Res, &
                     update_upwind_direction, &
                     count_upwind_direction_flip)
      patch%boundary_velocities(:,sum_connection) = v_darcy

      if (associated(patch%boundary_flow_fluxes)) then
        patch%boundary_flow_fluxes(:,sum_connection) = Res(:)
      endif

      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_bc(sum_connection)%mass_balance_delta(wid,1) = &
          global_auxvars_bc(sum_connection)%mass_balance_delta(wid,1) - &
          Res(SCO2_WATER_EQUATION_INDEX)
        global_auxvars_bc(sum_connection)%mass_balance_delta(co2_id,1) = &
          global_auxvars_bc(sum_connection)%mass_balance_delta(co2_id,1) - &
          Res(SCO2_CO2_EQUATION_INDEX)
        global_auxvars_bc(sum_connection)%mass_balance_delta(sid,1) = &
          global_auxvars_bc(sum_connection)%mass_balance_delta(sid,1) - &
          Res(SCO2_SALT_EQUATION_INDEX)
      endif

      local_end = local_id * option%nflowdof
      local_start = local_end - option%nflowdof + 1
      r_p(local_start:local_end)= r_p(local_start:local_end) - Res(:)

    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Update well source/sink terms
  if (sco2_well_coupling == SCO2_FULLY_IMPLICIT_WELL) then
    if (associated(pm_well)) then
      if (any(pm_well%well_grid%h_rank_id == option%myrank)) then
        call pm_well%UpdateFlowRates(ZERO_INTEGER,ZERO_INTEGER,-999,ierr)
        call pm_well%ModifyFlowResidual(r_p,ss_flow_vol_flux)
      endif
    endif
  endif

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit
    if (associated(source_sink%flow_condition%well)) then
      source_sink => source_sink%next
      cycle
    endif

    cur_connection_set => source_sink%connection_set

    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle


      local_end = local_id * option%nflowdof
      local_start = local_end - option%nflowdof + 1

      if (sco2_well_coupling == SCO2_FULLY_IMPLICIT_WELL) then
        ! nflowdof is padded with an extra dof for potential wells
        local_end = local_end - 1
      endif

      if (associated(source_sink%flow_aux_real_var)) then
        scale = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
      else
        scale = 1.d0
      endif

      qsrc=source_sink%flow_condition%sco2%rate%dataset%rarray(:)
      flow_src_sink_type=source_sink%flow_condition%sco2%rate%itype

      ! Index 0 contains user-specified conditions
      ! Index 1 contains auxvars to be used in src/sink calculations
      call SCO2AuxVarComputeAndSrcSink(option,qsrc,flow_src_sink_type, &
                            sco2_auxvars_ss(ZERO_INTEGER,sum_connection), &
                            sco2_auxvars(ZERO_INTEGER,ghosted_id), &
                            global_auxvars(ghosted_id), &
                            global_auxvars_ss(sum_connection), &
                            material_auxvars(ghosted_id), &
                            patch%characteristic_curves_array( &
                            patch%cc_id(ghosted_id))%ptr, &
                            sco2_parameter, &
                            grid%nG2A(ghosted_id), &
                            scale,Res,PETSC_FALSE)

      r_p(local_start:local_end) =  r_p(local_start:local_end) - Res(:)

      if (associated(patch%ss_flow_vol_fluxes)) then
        patch%ss_flow_vol_fluxes(:,sum_connection) = ss_flow_vol_flux
      endif
      if (associated(patch%ss_flow_fluxes)) then
        patch%ss_flow_fluxes(:,sum_connection) = Res(:)
      endif
      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_ss(sum_connection)%mass_balance_delta(wid,1) = &
          global_auxvars_ss(sum_connection)%mass_balance_delta(wid,1) - &
          Res(SCO2_WATER_EQUATION_INDEX)
        global_auxvars_ss(sum_connection)%mass_balance_delta(co2_id,1) = &
          global_auxvars_ss(sum_connection)%mass_balance_delta(co2_id,1) - &
          Res(SCO2_CO2_EQUATION_INDEX)
        global_auxvars_ss(sum_connection)%mass_balance_delta(sid,1) = &
          global_auxvars_ss(sum_connection)%mass_balance_delta(sid,1) - &
          Res(SCO2_SALT_EQUATION_INDEX)
      endif

    enddo
    source_sink => source_sink%next
  enddo

  if (patch%aux%SCO2%inactive_cells_exist) then
    do i=1,patch%aux%SCO2%matrix_zeroing%n_inactive_rows
      r_p(patch%aux%SCO2%matrix_zeroing%inactive_rows_local(i)) = 0.d0
    enddo
  endif

  call VecRestoreArrayF90(r,r_p,ierr);CHKERRQ(ierr)

  ! Mass Transfer
  if (field%flow_mass_transfer /= PETSC_NULL_VEC) then
    ! scale by -1.d0 for contribution to residual.  A negative contribution
    ! indicates mass being added to system.
    call VecAXPY(r,-1.d0,field%flow_mass_transfer,ierr);CHKERRQ(ierr)
  endif

  if (realization%debug%vecview_residual) then
    call DebugWriteFilename(realization%debug,string,'Sresidual','', &
                            sco2_ts_count,sco2_ts_cut_count, &
                            sco2_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif
  if (realization%debug%vecview_solution) then
    call DebugWriteFilename(realization%debug,string,'Sxx','', &
                            sco2_ts_count,sco2_ts_cut_count, &
                            sco2_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif

  update_upwind_direction = PETSC_FALSE

end subroutine SCO2Residual

! ************************************************************************** !

subroutine SCO2Jacobian(snes,xx,A,B,realization,pm_well,ierr)
  !
  ! Computes the Jacobian
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
  !

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Connection_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Material_Aux_module
  use Upwind_Direction_module
  use PM_Well_class

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  class(realization_subsurface_type) :: realization
  class(pm_well_type), pointer :: pm_well
  PetscErrorCode :: ierr

  Mat :: J
  MatType :: mat_type
  PetscReal :: norm
  PetscViewer :: viewer

  PetscReal :: qsrc, scale
  PetscInt :: imat, imat_up, imat_dn
  PetscInt :: local_id, ghosted_id, natural_id
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
  ! PetscInt, pointer :: zeros(:)
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(material_parameter_type), pointer :: material_parameter
  type(sco2_parameter_type), pointer :: sco2_parameter
  type(sco2_auxvar_type), pointer :: sco2_auxvars(:,:), &
                                        sco2_auxvars_bc(:), &
                                        sco2_auxvars_ss(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:), &
                                       global_auxvars_ss(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: well_dof
  PetscInt :: deactivate_row

  well_dof = ZERO_INTEGER
  deactivate_row = UNINITIALIZED_INTEGER
  if (associated(pm_well)) well_dof = ONE_INTEGER

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  sco2_auxvars => patch%aux%SCO2%auxvars
  sco2_auxvars_bc => patch%aux%SCO2%auxvars_bc
  sco2_auxvars_ss => patch%aux%SCO2%auxvars_ss
  sco2_parameter => patch%aux%SCO2%sco2_parameter
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars

  call SNESGetIterationNumber(snes,sco2_newton_iteration_number, &
                              ierr);CHKERRQ(ierr)
  sco2_newton_iteration_number = sco2_newton_iteration_number + 1

  sco2_sub_newton_iter_num = 0
  sco2_force_iteration = PETSC_FALSE

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

    call SCO2AuxVarPerturb(sco2_auxvars(:,ghosted_id), &
                           global_auxvars(ghosted_id), &
                           material_auxvars(ghosted_id), &
                           patch%characteristic_curves_array( &
                           patch%cc_id(ghosted_id))%ptr, &
                           sco2_parameter,natural_id,option)
  enddo

  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    call SCO2AccumDerivative(sco2_auxvars(:,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              material_parameter%soil_heat_capacity(imat), &
                              option,well_dof,Jup)
    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                  ADD_VALUES,ierr);CHKERRQ(ierr)
  enddo

  if (realization%debug%matview_Matrix_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call DebugWriteFilename(realization%debug,string,'Sjacobian_accum','', &
                            sco2_ts_count,sco2_ts_cut_count, &
                            sco2_ni_count)
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
      if (.not. Initialized(cur_connection_set%face_id(iconn))) cycle
      if (cur_connection_set%area(iconn) <= 0.d0) cycle

      sum_connection = sum_connection + 1

      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)

      imat_up = patch%imat(ghosted_id_up)
      imat_dn = patch%imat(ghosted_id_dn)
      if (imat_up <= 0 .or. imat_dn <= 0) cycle

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping

      call SCO2FluxDerivative(sco2_auxvars(:,ghosted_id_up), &
                     global_auxvars(ghosted_id_up), &
                     material_auxvars(ghosted_id_up), &
                     patch%char_curves_thermal_array( &
                       patch%cct_id(ghosted_id_up))%ptr, &
                     sco2_auxvars(:,ghosted_id_dn), &
                     global_auxvars(ghosted_id_dn), &
                     material_auxvars(ghosted_id_dn), &
                     patch%char_curves_thermal_array( &
                       patch%cct_id(ghosted_id_dn))%ptr, &
                     cur_connection_set%area(iconn), &
                     cur_connection_set%dist(:,iconn), &
                     patch%flow_upwind_direction(:,iconn), &
                     option,well_dof,Jup,Jdn)
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

  if (realization%debug%matview_Matrix_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call DebugWriteFilename(realization%debug,string,'Sjacobian_flux','', &
                            sco2_ts_count,sco2_ts_cut_count, &
                            sco2_ni_count)
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

      call SCO2BCFluxDerivative(boundary_condition%flow_bc_type, &
                      boundary_condition%flow_aux_mapping, &
                      boundary_condition%flow_aux_real_var(:,iconn), &
                      sco2_auxvars_bc(sum_connection), &
                      global_auxvars_bc(sum_connection), &
                      sco2_auxvars(:,ghosted_id), &
                      global_auxvars(ghosted_id), &
                      material_auxvars(ghosted_id), &
                      patch%char_curves_thermal_array( &
                      patch%cct_id(ghosted_id))%ptr, &
                      cur_connection_set%area(iconn), &
                      cur_connection_set%dist(:,iconn), &
                      patch%flow_upwind_direction_bc(:,iconn), &
                      option,well_dof,Jdn)

      Jdn = -Jdn
      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  if (realization%debug%matview_Matrix_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call DebugWriteFilename(realization%debug,string,'Sjacobian_bcflux','', &
                            sco2_ts_count,sco2_ts_cut_count, &
                            sco2_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif

  ! Source/sinks
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit
    if (associated(source_sink%flow_condition%well)) then
      source_sink => source_sink%next
      cycle
    endif

    cur_connection_set => source_sink%connection_set

    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1

      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)

      if (patch%imat(ghosted_id) <= 0 .or. &
          associated(source_sink%flow_condition%well)) cycle

      if (associated(source_sink%flow_aux_real_var)) then
        scale = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
      else
        scale = 1.d0
      endif

      Jup = 0.d0
      call SCO2SrcSinkDerivative(option,source_sink, &
                        sco2_auxvars_ss(:,sum_connection), &
                        sco2_auxvars(:,ghosted_id), &
                        global_auxvars(ghosted_id), &
                        global_auxvars_ss(sum_connection), &
                        patch%characteristic_curves_array( &
                          patch%cc_id(ghosted_id))%ptr, &
                        sco2_parameter, grid%nG2A(ghosted_id),&
                        material_auxvars(ghosted_id), &
                        scale,well_dof,Jup)

      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)

    enddo
    source_sink => source_sink%next
  enddo

  ! Well Terms
  ! Need to update all well source/sink terms wrt
  ! perturbation in bottom pressure
  if (sco2_well_coupling == FULLY_IMPLICIT_WELL) then
    if (associated(pm_well)) then
      if (any(pm_well%well_grid%h_rank_id == option%myrank)) then
        ! Perturb the well and well's reservoir variables.
        call pm_well%Perturb()

        ! Go through and update the well contributions to the Jacobian:
        ! dRi/d(P_well), dRwell/d(P_well), dRi/dxi, and dRwell,dxi
        call pm_well%ModifyFlowJacobian(A,ierr)

        ! Deactivate unused rows
        if (all(pm_well%well%liq%Q == 0.d0) .and. &
            all(pm_well%well%gas%Q == 0.d0)) then
          if (pm_well%well_grid%h_rank_id(1) == option%myrank) then
            deactivate_row = pm_well%well_grid%h_ghosted_id(1) * &
                             option%nflowdof
          endif
        endif

      endif
    endif
  endif

  if (realization%debug%matview_Matrix_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call DebugWriteFilename(realization%debug,string,'Sjacobian_srcsink','', &
                            sco2_ts_count,sco2_ts_cut_count, &
                            sco2_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  ! zero out isothermal and inactive cells
  if (patch%aux%SCO2%inactive_cells_exist) then
    qsrc = 1.d0 ! solely a temporary variable in this conditional
    call MatZeroRowsLocal(A,patch%aux%SCO2%matrix_zeroing%n_inactive_rows, &
                          patch%aux%SCO2%matrix_zeroing% &
                            inactive_rows_local_ghosted, &
                          qsrc,PETSC_NULL_VEC,PETSC_NULL_VEC, &
                          ierr);CHKERRQ(ierr)
    if (Initialized(deactivate_row)) then
      deactivate_row = deactivate_row - 1
      call MatZeroRowsLocal(A,ONE_INTEGER, deactivate_row, &
                          qsrc,PETSC_NULL_VEC,PETSC_NULL_VEC, &
                          ierr);CHKERRQ(ierr)
    endif
  endif

  if (realization%debug%matview_Matrix) then
    call DebugWriteFilename(realization%debug,string,'Sjacobian','', &
                            sco2_ts_count,sco2_ts_cut_count, &
                            sco2_ni_count)
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

  ! update after evaluations to ensure zero-based index to match screen output
  sco2_ni_count = sco2_ni_count + 1

end subroutine SCO2Jacobian

! ************************************************************************** !

subroutine SCO2SetPlotVariables(realization,list)
  !
  ! Adds variables to be printed to list
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
  !

  use Realization_Subsurface_class
  use Output_Aux_module
  use Variables_module

  implicit none

  class(realization_subsurface_type) :: realization
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

    name = 'CO2 Partial Pressure'
    units = 'Pa'
    call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                                CO2_PRESSURE)

    name = 'Vapor Pressure'
    units = 'Pa'
    call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                                VAPOR_PRESSURE)

    name = 'Liquid Saturation'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                                LIQUID_SATURATION)

    name = 'Gas Saturation'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                                GAS_SATURATION)

    name = 'Salt Precipitate Saturation'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                                PRECIPITATE_SATURATION)

    name = 'Trapped Gas Saturation'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                                TRAPPED_GAS_SATURATION)

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
    name = 'X_s^l'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                 LIQUID_MOLE_FRACTION, &
                                 realization%option%salt_id)

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

end subroutine SCO2SetPlotVariables

! ************************************************************************** !

subroutine SCO2MapBCAuxVarsToGlobal(realization)
  !
  ! Maps variables in SCO2 auxvar to global equivalent.
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
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
  type(sco2_auxvar_type), pointer :: sco2_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)

  PetscInt :: sum_connection, iconn

  option => realization%option
  patch => realization%patch

  if (option%ntrandof == 0) return ! no need to update

  sco2_auxvars_bc => patch%aux%SCO2%auxvars_bc
  global_auxvars_bc => patch%aux%Global%auxvars_bc

  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      global_auxvars_bc(sum_connection)%sat = &
        sco2_auxvars_bc(sum_connection)%sat
      global_auxvars_bc(sum_connection)%den_kg = &
        sco2_auxvars_bc(sum_connection)%den_kg
      global_auxvars_bc(sum_connection)%temp = &
        sco2_auxvars_bc(sum_connection)%temp
    enddo
    boundary_condition => boundary_condition%next
  enddo

end subroutine SCO2MapBCAuxVarsToGlobal

! ************************************************************************** !

subroutine SCO2Destroy(realization)
  !
  ! Deallocates variables associated with SCO2
  !
  ! Author: Michael Nole
  ! Date: 01/26/24
  !

  use Realization_Subsurface_class

  implicit none

  class(realization_subsurface_type) :: realization

end subroutine SCO2Destroy

end module SCO2_module
