module TH_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use TH_Aux_module
  use Global_Aux_module
  use Material_Aux_module
  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal

  implicit none

  private

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.D-8
  PetscReal, parameter :: floweps   = 1.D-24
  PetscReal, parameter :: perturbation_tolerance = 1.d-8
  PetscReal, parameter :: unit_z(3) = [0.d0,0.d0,1.d0]

  public THResidual, &
         THJacobian, &
         THUpdateFixedAccumulation, &
         THTimeCut, &
         THSetup, &
         THNumericalJacobianTest, &
         THMaxChange, &
         THUpdateSolution, &
         THGetTecplotHeader, &
         THInitializeTimestep, &
         THComputeMassBalance, &
         THResidualToMass, &
         THSecondaryHeat, &
         THSecondaryHeatJacobian, &
         THUpdateAuxVars, &
         THDestroy, &
         THAccumulation, &
         THResidualInternalConn, &
         THResidualBoundaryConn, &
         THResidualSourceSink, &
         THJacobianInternalConn, &
         THJacobianBoundaryConn, &
         THJacobianSourceSink, &
         THUpdateLocalVecs

  PetscInt, parameter :: jh2o = 1

contains

! ************************************************************************** !

subroutine THTimeCut(realization)
  !
  ! Resets arrays for time step cut
  !
  ! Author: ???
  ! Date: 12/13/07
  !

  use Realization_Subsurface_class

  implicit none

  class(realization_subsurface_type) :: realization

  TH_ts_cut_count = TH_ts_cut_count + 1
  call THInitializeTimestep(realization)

end subroutine THTimeCut

! ************************************************************************** !

subroutine THSetup(realization)
  !
  ! Author: ???
  ! Date: 02/22/08
  !

  use Realization_Subsurface_class
  use Patch_module
  use Output_Aux_module

  class(realization_subsurface_type) :: realization

  type(patch_type), pointer :: cur_patch
  type(output_variable_list_type), pointer :: list

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call THSetupPatch(realization)
    cur_patch => cur_patch%next
  enddo

  list => realization%output_option%output_snap_variable_list
  call THSetPlotVariables(realization,list)
  list => realization%output_option%output_obs_variable_list
  call THSetPlotVariables(realization,list)

  TH_ts_count = 0
  TH_ts_cut_count = 0
  TH_ni_count = 0

end subroutine THSetup

! ************************************************************************** !

subroutine THSetupPatch(realization)
  !
  ! Creates arrays for auxiliary variables
  !
  ! Author: ???
  ! Date: 02/22/08
  !

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Grid_module
  use Region_module
  use Coupler_module
  use Connection_module
  use Fluid_module
  use Secondary_Continuum_Aux_module
  use Secondary_Continuum_module
  use Characteristic_Curves_Thermal_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(TH_auxvar_type), pointer :: TH_auxvars(:), TH_auxvars_bc(:)
  type(TH_auxvar_type), pointer :: TH_auxvars_ss(:)
  type(fluid_property_type), pointer :: cur_fluid_property
  type(sec_heat_type), pointer :: TH_sec_heat_vars(:)
  type(coupler_type), pointer :: initial_condition
  class(cc_thermal_type), pointer :: thermal_cc
  character(len=MAXWORDLENGTH) :: word
  PetscReal :: area_per_vol

  PetscInt :: ghosted_id, iconn, sum_connection
  PetscInt :: i, iphase, local_id, material_id, icct
  PetscBool :: error_found
  PetscErrorCode :: ierr

  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%aux%TH => THAuxCreate(option)
  patch%aux%SC_heat => SecondaryAuxHeatCreate(option)

! option%io_buffer = 'Before TH can be run, the th_parameter object ' // &
!                    'must be initialized with the proper variables ' // &
!                    'THAuxCreate() is called anywhere.'

! call printErrMsg(option)

  if (option%flow%th_freezing) then
    allocate(patch%aux%TH%th_parameter%sir(option%nphase, &
              size(patch%saturation_function_array)))
  else
    allocate(patch%aux%TH%th_parameter%sir(option%nphase, &
              size(patch%characteristic_curves_array)))

  endif


  !Jitu, 08/04/2010: Check these allocations. Currently assumes only
  !single value in the array <modified pcl 1-13-11>
  allocate(patch%aux%TH%th_parameter%dencpr( &
             size(patch%char_curves_thermal_array)))
  allocate(patch%aux%TH%th_parameter% &
           ckwet(size(patch%char_curves_thermal_array)))
  allocate(patch%aux%TH%th_parameter% &
           ckdry(size(patch%char_curves_thermal_array)))
  allocate(patch%aux%TH%th_parameter% &
           alpha(size(patch%char_curves_thermal_array)))
  if (option%flow%th_freezing) then
   allocate(patch%aux%TH%th_parameter%ckfrozen( &
              size(patch%char_curves_thermal_array)))
   allocate(patch%aux%TH%th_parameter%alpha_fr( &
              size(patch%char_curves_thermal_array)))
  endif

  !Copy the values in the th_parameter from the global realization
  error_found = PETSC_FALSE
  do i = 1, size(patch%material_property_array)
    word = patch%material_property_array(i)%ptr%name
    if (Uninitialized(patch%material_property_array(i)%ptr%specific_heat)) then
      option%io_buffer = 'ERROR: Non-initialized HEAT_CAPACITY in material ' &
                         // trim(word)
      call PrintMsgByRank(option)
      error_found = PETSC_TRUE
    endif

    material_id = abs(patch%material_property_array(i)%ptr%internal_id)

    icct = patch%material_property_array(i)%ptr% &
                thermal_conductivity_function_id

    thermal_cc => patch%char_curves_thermal_array(icct)%ptr

    ! kg rock/m^3 rock * J/kg rock-K * 1.e-6 MJ/J = MJ/m^3-K
    patch%aux%TH%th_parameter%dencpr(icct) = &
      patch%material_property_array(i)%ptr%rock_density*option%scale* &
        patch%material_property_array(i)%ptr%specific_heat

    select type(tcf => thermal_cc%thermal_conductivity_function)
      !------------------------------------------
      type is(kT_frozen_type)
        patch%aux%TH%th_parameter%ckdry(icct) = tcf%kT_dry*option%scale
        patch%aux%TH%th_parameter%ckwet(icct) = tcf%kT_wet*option%scale
        tcf%kT_dry = tcf%kT_dry*option%scale ! apply scale to original value
        tcf%kT_wet = tcf%kT_wet*option%scale ! apply scale to original value
        patch%aux%TH%th_parameter%alpha(icct) = tcf%alpha
        if (option%flow%th_freezing) then
          patch%aux%TH%th_parameter%alpha_fr(icct) = tcf%alpha_fr
          patch%aux%TH%th_parameter%ckfrozen(icct) = &
          tcf%kT_frozen*option%scale
          tcf%kT_frozen = tcf%kT_frozen*option%scale ! apply scale to original value
        endif
      !------------------------------------------
      type is(kT_constant_type)
        patch%aux%TH%th_parameter%ckdry(icct) = &
          tcf%constant_thermal_conductivity*option%scale
        patch%aux%TH%th_parameter%ckwet(icct) = &
          tcf%constant_thermal_conductivity*option%scale
        tcf%constant_thermal_conductivity = &
          tcf%constant_thermal_conductivity*option%scale ! apply scale to original value
        patch%aux%TH%th_parameter%alpha(icct) = tcf%alpha
      !------------------------------------------
      type is(kT_default_type)
        patch%aux%TH%th_parameter%ckdry(icct) = tcf%kT_dry*option%scale
        patch%aux%TH%th_parameter%ckwet(icct) = tcf%kT_wet*option%scale
        tcf%kT_dry = tcf%kT_dry*option%scale ! apply scale to original value
        tcf%kT_wet = tcf%kT_wet*option%scale ! apply scale to original value
        patch%aux%TH%th_parameter%alpha(icct) = tcf%alpha
      !------------------------------------------
      type is(kT_linear_type)
        patch%aux%TH%th_parameter%ckdry(icct) = tcf%kT_dry*option%scale
        patch%aux%TH%th_parameter%ckwet(icct) = tcf%kT_wet*option%scale
        tcf%kT_dry = tcf%kT_dry*option%scale ! apply scale to original value
        tcf%kT_wet = tcf%kT_wet*option%scale ! apply scale to original value
        patch%aux%TH%th_parameter%alpha(icct) = tcf%alpha
      !------------------------------------------
      type is(kT_linear_resistivity_type)
        patch%aux%TH%th_parameter%ckdry(icct) = tcf%kT_dry*option%scale
        patch%aux%TH%th_parameter%ckwet(icct) = tcf%kT_wet*option%scale
        tcf%kT_dry = tcf%kT_dry*option%scale ! apply scale to original value
        tcf%kT_wet = tcf%kT_wet*option%scale ! apply scale to original value
        patch%aux%TH%th_parameter%alpha(icct) = tcf%alpha
      !------------------------------------------
      type is(kT_cubic_polynomial_type)
        patch%aux%TH%th_parameter%ckdry(icct) = tcf%kT_dry*option%scale
        patch%aux%TH%th_parameter%ckwet(icct) = tcf%kT_wet*option%scale
        tcf%kT_dry = tcf%kT_dry*option%scale ! apply scale to original value
        tcf%kT_wet = tcf%kT_wet*option%scale ! apply scale to original value
        patch%aux%TH%th_parameter%alpha(icct) = tcf%alpha
      !------------------------------------------
      type is(kT_power_type)
        patch%aux%TH%th_parameter%ckdry(icct) = tcf%kT_dry*option%scale
        patch%aux%TH%th_parameter%ckwet(icct) = tcf%kT_wet*option%scale
        tcf%kT_dry = tcf%kT_dry*option%scale ! apply scale to original value
        tcf%kT_wet = tcf%kT_wet*option%scale ! apply scale to original value
        patch%aux%TH%th_parameter%alpha(icct) = tcf%alpha
      class default
        option%io_buffer = 'ERROR: TH mode does not support thermal '&
                         //'characteristic curve "' // trim(thermal_cc%name) &
                         //'" in material: ' // trim(word)
        call PrintErrMsg(option)
    end select

    if (.not. associated(thermal_cc%thermal_conductivity_function)) then
      option%io_buffer = 'ERROR: Non-initialized thermal conductivity function &
                         &in material ' // trim(word)
      call PrintMsgByRank(option)
      error_found = PETSC_TRUE
    endif

    if (Uninitialized(patch%aux%TH%th_parameter%alpha(icct))) then
      option%io_buffer = 'ERROR: Non-initialized KERSTEN EXPONENT '&
                       //'in material ' // trim(word)
      call PrintMsgByRank(option)
      error_found = PETSC_TRUE
    endif

    if (option%flow%th_freezing) then
      if (patch%aux%TH%th_parameter%ckfrozen(icct) < 0.0d0 ) then
        option%io_buffer = 'ERROR: Non-initialized FROZEN THERMAL '&
        //'CONDUCTIVITY when freezing activated in '&
        //'material ' // trim(word)
        call PrintMsgByRank(option)
        error_found = PETSC_TRUE
      endif
      if (Uninitialized(patch%aux%TH%th_parameter%alpha_fr(icct))) then
        option%io_buffer = 'ERROR: Non-initialized FROZEN KERSTEN EXPONENT '&
                         //'when freezing activated in material ' // trim(word)
        call PrintMsgByRank(option)
        error_found = PETSC_TRUE
      endif
    endif

    if (patch%aux%TH%th_parameter%ckwet(icct) < 1.d-40 .and. &
        patch%aux%TH%th_parameter%ckdry(icct) < 1.d-40) then
      option%io_buffer = 'ERROR: Either the wet or dry thermal conductivity &
        &must be non-zero in material: ' // trim(word)
      call PrintMsgByRank(option)
      error_found = PETSC_TRUE
    endif

  enddo

  call MPI_Allreduce(MPI_IN_PLACE,error_found,ONE_INTEGER_MPI,MPI_LOGICAL, &
                     MPI_LOR,option%mycomm,ierr);CHKERRQ(ierr)
  if (error_found) then
    option%io_buffer = 'Material property errors found in THSetup.'
    call PrintErrMsg(option)
  endif

  if (option%flow%th_freezing) then
    do i = 1, size(patch%saturation_function_array)
      patch%aux%TH%th_parameter% &
        sir(:,patch%saturation_function_array(i)%ptr%id) = &
        patch%saturation_function_array(i)%ptr%Sr(:)
    enddo
  endif

  ! allocate auxvar data structures for all grid cells
  allocate(TH_auxvars(grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    call THAuxVarInit(TH_auxvars(ghosted_id),option)
  enddo


  if (option%use_sc) then
    initial_condition => patch%initial_condition_list%first
    allocate(TH_sec_heat_vars(grid%nlmax))

    do local_id = 1, grid%nlmax

      ghosted_id = grid%nL2G(local_id)
    ! Assuming the same secondary continuum for all regions (need to
    ! make it an array)
    ! S. Karra 07/18/12
      call SecondaryContinuumSetProperties( &
        TH_sec_heat_vars(local_id)%sec_continuum, &
        patch%material_property_array(1)%ptr%multicontinuum%name, &
        patch%aux%Material%auxvars(ghosted_id)%soil_properties(half_matrix_width_index), &
        patch%material_property_array(1)%ptr% &
          multicontinuum%matrix_block_size, &
        patch%material_property_array(1)%ptr% &
          multicontinuum%fracture_spacing, &
        patch%material_property_array(1)%ptr%multicontinuum%radius, &
        patch%material_property_array(1)%ptr%multicontinuum%porosity, &
        option)

      TH_sec_heat_vars(local_id)%ncells = &
        patch%material_property_array(1)%ptr%multicontinuum%ncells
      TH_sec_heat_vars(local_id)%half_aperture = &
        patch%material_property_array(1)%ptr%multicontinuum%half_aperture
      TH_sec_heat_vars(local_id)%epsilon = &
        patch%aux%Material%auxvars(ghosted_id)%soil_properties(epsilon_index)
      TH_sec_heat_vars(local_id)%log_spacing = &
        patch%material_property_array(1)%ptr%multicontinuum%log_spacing
      TH_sec_heat_vars(local_id)%outer_spacing = &
        patch%material_property_array(1)%ptr%multicontinuum%outer_spacing

      allocate(TH_sec_heat_vars(local_id)%area( &
                               TH_sec_heat_vars(local_id)%ncells))
      allocate(TH_sec_heat_vars(local_id)%vol( &
                               TH_sec_heat_vars(local_id)%ncells))
      allocate(TH_sec_heat_vars(local_id)%dm_minus( &
                               TH_sec_heat_vars(local_id)%ncells))
      allocate(TH_sec_heat_vars(local_id)%dm_plus( &
                               TH_sec_heat_vars(local_id)%ncells))
      allocate(TH_sec_heat_vars(local_id)%sec_continuum% &
             distance(TH_sec_heat_vars(local_id)%ncells))

      call SecondaryContinuumType(TH_sec_heat_vars(local_id)%sec_continuum, &
                                  TH_sec_heat_vars(local_id)%ncells, &
                                  TH_sec_heat_vars(local_id)%area, &
                                  TH_sec_heat_vars(local_id)%vol, &
                                  TH_sec_heat_vars(local_id)%dm_minus, &
                                  TH_sec_heat_vars(local_id)%dm_plus, &
                                  TH_sec_heat_vars(local_id)%half_aperture, &
                                  TH_sec_heat_vars(local_id)%epsilon, &
                                  TH_sec_heat_vars(local_id)%log_spacing, &
                                  TH_sec_heat_vars(local_id)%outer_spacing, &
                                  area_per_vol,option)

      TH_sec_heat_vars(local_id)%interfacial_area = area_per_vol* &
        (1.d0 - TH_sec_heat_vars(local_id)%epsilon)* &
        patch%material_property_array(1)%ptr% &
        multicontinuum%area_scaling

    ! Setting the initial values of all secondary node temperatures same
    ! as primary node temperatures (with initial dirichlet BC only)
    ! -- sk 06/26/12
      allocate(TH_sec_heat_vars(local_id)%sec_temp( &
                                 TH_sec_heat_vars(local_id)%ncells))

      if (option%flow%set_secondary_init_temp) then
        TH_sec_heat_vars(local_id)%sec_temp = &
          patch%material_property_array(1)%ptr%multicontinuum%init_temp
      else
        TH_sec_heat_vars(local_id)%sec_temp = &
        initial_condition%flow_condition%temperature%dataset%rarray(1)
      endif

    enddo

    patch%aux%SC_heat%sec_heat_vars => TH_sec_heat_vars

  endif


  patch%aux%TH%auxvars => TH_auxvars
  patch%aux%TH%num_aux = grid%ngmax

  ! count the number of boundary connections and allocate
  ! auxvar data structures for them
  boundary_condition => patch%boundary_condition_list%first

  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    sum_connection = sum_connection + &
                     boundary_condition%connection_set%num_connections
    boundary_condition => boundary_condition%next
  enddo

  if (sum_connection > 0) then
    allocate(TH_auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call THAuxVarInit(TH_auxvars_bc(iconn),option)
    enddo
    patch%aux%TH%auxvars_bc => TH_auxvars_bc
  endif
  patch%aux%TH%num_aux_bc = sum_connection

  ! Create aux vars for source/sink
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (sum_connection > 0) then
    allocate(TH_auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call THAuxVarInit(TH_auxvars_ss(iconn),option)
    enddo
    patch%aux%TH%auxvars_ss => TH_auxvars_ss
  endif
  patch%aux%TH%num_aux_ss = sum_connection

  ! initialize parameters
  cur_fluid_property => realization%fluid_properties
  do
    if (.not.associated(cur_fluid_property)) exit
    iphase = cur_fluid_property%phase_id
    cur_fluid_property => cur_fluid_property%next
  enddo

end subroutine THSetupPatch

! ************************************************************************** !

subroutine THComputeMassBalance(realization, mass_balance)
  !
  ! THomputeMassBalance:
  ! Adapted from RichardsComputeMassBalance: need to be checked
  !
  ! Author: Jitendra Kumar
  ! Date: 07/21/2010
  !

  use Realization_Subsurface_class
  use Patch_module

  class(realization_subsurface_type) :: realization
  PetscReal :: mass_balance(realization%option%nphase)

  type(patch_type), pointer :: cur_patch

  mass_balance = 0.d0

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call THComputeMassBalancePatch(realization, mass_balance)
    cur_patch => cur_patch%next
  enddo

end subroutine THComputeMassBalance

! ************************************************************************** !

subroutine THComputeMassBalancePatch(realization,mass_balance)
  !
  ! THomputeMassBalancePatch:
  ! Adapted from RichardsComputeMassBalancePatch: need to be checked
  !
  ! Author: Jitendra Kumar
  ! Date: 07/21/2010
  !

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
  use Material_Aux_module, only : material_auxvar_type

  implicit none

  class(realization_subsurface_type) :: realization
  PetscReal :: mass_balance(realization%option%nphase)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(TH_auxvar_type),pointer :: TH_auxvars(:)

  PetscInt :: local_id
  PetscInt :: ghosted_id

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  TH_auxvars => patch%aux%TH%auxvars

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    ! mass = volume*saturation*density

    mass_balance = mass_balance + &
      global_auxvars(ghosted_id)%den_kg* &
      global_auxvars(ghosted_id)%sat* &
      material_auxvars(ghosted_id)%porosity* &
      material_auxvars(ghosted_id)%volume

    if (option%flow%th_freezing) then
      ! mass = volume*saturation_ice*density_ice
      mass_balance = mass_balance + &
        TH_auxvars(ghosted_id)%ice%den_ice*FMWH2O* &
        TH_auxvars(ghosted_id)%ice%sat_ice* &
        material_auxvars(ghosted_id)%porosity* &
        material_auxvars(ghosted_id)%volume
    endif

  enddo

end subroutine THComputeMassBalancePatch

! ************************************************************************** !

subroutine THZeroMassBalDeltaPatch(realization)
  !
  ! Zeros mass balance delta array
  !
  ! Author: Satish Karra, LANL
  ! Date: 12/13/11
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

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%TH%num_aux
    patch%aux%Global%auxvars(iconn)%mass_balance_delta = 0.d0
  enddo
#endif

  ! Intel 10.1 on Chinook reports a SEGV if this conditional is not
  ! placed around the internal do loop - geh
  if (patch%aux%TH%num_aux_bc > 0) then
    do iconn = 1, patch%aux%TH%num_aux_bc
      global_auxvars_bc(iconn)%mass_balance_delta = 0.d0
    enddo
  endif
  if (patch%aux%TH%num_aux_ss > 0) then
    do iconn = 1, patch%aux%TH%num_aux_ss
      global_auxvars_ss(iconn)%mass_balance_delta = 0.d0
    enddo
  endif

end subroutine THZeroMassBalDeltaPatch

! ************************************************************************** !

subroutine THUpdateMassBalancePatch(realization)
  !
  ! Updates mass balance
  !
  ! Author: ???
  ! Date: 12/13/11
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

#ifdef COMPUTE_INTERNAL_MASS_FLUX
  do iconn = 1, patch%aux%TH%num_aux
    patch%aux%Global%auxvars(iconn)%mass_balance = &
      patch%aux%Global%auxvars(iconn)%mass_balance + &
      patch%aux%Global%auxvars(iconn)%mass_balance_delta*FMWH2O* &
      option%flow_dt
  enddo
#endif

  if (patch%aux%TH%num_aux_bc > 0) then
    do iconn = 1, patch%aux%TH%num_aux_bc
      global_auxvars_bc(iconn)%mass_balance = &
        global_auxvars_bc(iconn)%mass_balance + &
        global_auxvars_bc(iconn)%mass_balance_delta*FMWH2O*option%flow_dt
    enddo
  endif
  if (patch%aux%TH%num_aux_ss > 0) then
    do iconn = 1, patch%aux%TH%num_aux_ss
      global_auxvars_ss(iconn)%mass_balance = &
        global_auxvars_ss(iconn)%mass_balance + &
        global_auxvars_ss(iconn)%mass_balance_delta*FMWH2O*option%flow_dt
    enddo
  endif


end subroutine THUpdateMassBalancePatch

! ************************************************************************** !

subroutine THUpdateAuxVars(realization)
  !
  ! Updates the auxiliary variables associated with
  ! the TH problem
  !
  ! Author: ???
  ! Date: 12/10/07
  !

  use Realization_Subsurface_class
  use Patch_module

  class(realization_subsurface_type) :: realization

  type(patch_type), pointer :: cur_patch

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call THUpdateAuxVarsPatch(realization)
    cur_patch => cur_patch%next
  enddo

end subroutine THUpdateAuxVars

! ************************************************************************** !

subroutine THUpdateAuxVarsPatch(realization)
  !
  ! Updates the auxiliary variables associated with
  ! the TH problem
  !
  ! Author: ???
  ! Date: 12/10/07
  !

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(TH_auxvar_type), pointer :: TH_auxvars(:)
  type(TH_auxvar_type), pointer :: TH_auxvars_bc(:)
  type(TH_auxvar_type), pointer :: TH_auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(th_parameter_type), pointer :: th_parameter

  PetscInt :: ghosted_id, local_id, istart, iend, sum_connection, idof, iconn
  PetscInt :: iphasebc, iphase
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscReal, pointer :: xx(:)
  PetscReal :: tsrc1
  PetscErrorCode :: ierr
  PetscInt :: icct

  !!
!  PetscReal, allocatable :: gradient(:,:)
  !!

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  !!
!  allocate(gradient(grid%ngmax,3))
!  gradient = 0.d0
  !!

  TH_auxvars => patch%aux%TH%auxvars
  TH_auxvars_bc => patch%aux%TH%auxvars_bc
  TH_auxvars_ss => patch%aux%TH%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars
  th_parameter => patch%aux%TH%th_parameter

  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells

    if (patch%imat(ghosted_id) <= 0) cycle
    iend = ghosted_id*option%nflowdof
    istart = iend-option%nflowdof+1
    iphase = global_auxvars(ghosted_id)%istate
    icct = patch%cct_id(ghosted_id)

    if (option%flow%th_freezing) then
       call THAuxVarComputeFreezing(xx_loc_p(istart:iend), &
            TH_auxvars(ghosted_id),global_auxvars(ghosted_id), &
            material_auxvars(ghosted_id), &
            iphase, &
            patch%saturation_function_array(patch%cc_id(ghosted_id))%ptr, &
            patch%char_curves_thermal_array(patch%cct_id(ghosted_id))%ptr, &
            th_parameter, icct, &
            grid%nG2A(ghosted_id),PETSC_TRUE,option)
    else
       call THAuxVarComputeNoFreezing(xx_loc_p(istart:iend), &
            TH_auxvars(ghosted_id),global_auxvars(ghosted_id), &
            material_auxvars(ghosted_id), &
            iphase, &
            patch%characteristic_curves_array(patch%cc_id(ghosted_id))%ptr, &
            patch%char_curves_thermal_array(patch%cct_id(ghosted_id))%ptr, &
            th_parameter, icct, &
            grid%nG2A(ghosted_id),PETSC_TRUE,option)
    endif

    global_auxvars(ghosted_id)%istate = iphase
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
      if (patch%imat(ghosted_id) <= 0) cycle
      icct = patch%cct_id(ghosted_id)

      do idof=1,option%nflowdof
        select case(boundary_condition%flow_condition%itype(idof))
          case(DIRICHLET_BC,DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
               HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC, &
               HYDROSTATIC_CONDUCTANCE_BC, &
               HET_SURF_HYDROSTATIC_SEEPAGE_BC, &
               HET_DIRICHLET_BC,HET_HYDROSTATIC_SEEPAGE_BC, &
               HET_HYDROSTATIC_CONDUCTANCE_BC)
            xxbc(idof) = boundary_condition%flow_aux_real_var(idof,iconn)
          case(NEUMANN_BC,ZERO_GRADIENT_BC)
            xxbc(idof) = xx_loc_p((ghosted_id-1)*option%nflowdof+idof)
        end select
      enddo

      select case(boundary_condition%flow_condition%itype(TH_PRESSURE_DOF))
        case(DIRICHLET_BC,DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
             HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC, &
             HET_DIRICHLET_BC,HET_HYDROSTATIC_SEEPAGE_BC, &
             HET_HYDROSTATIC_CONDUCTANCE_BC)
          iphasebc = boundary_condition%flow_aux_int_var(1,iconn)
        case(NEUMANN_BC,ZERO_GRADIENT_BC)
          iphasebc = global_auxvars(ghosted_id)%istate
      end select

      if (option%flow%th_freezing) then
         call THAuxVarComputeFreezing(xxbc,TH_auxvars_bc(sum_connection), &
                                      global_auxvars_bc(sum_connection), &
              material_auxvars(ghosted_id), &
              iphasebc, &
              patch%saturation_function_array(patch%cc_id(ghosted_id))%ptr, &
              patch%char_curves_thermal_array(patch%cct_id(ghosted_id))%ptr, &
              th_parameter, icct, &
              -grid%nG2A(ghosted_id),PETSC_FALSE,option)
      else
         call THAuxVarComputeNoFreezing(xxbc,TH_auxvars_bc(sum_connection), &
              global_auxvars_bc(sum_connection), &
              material_auxvars(ghosted_id), &
              iphasebc, &
              patch%characteristic_curves_array(patch%cc_id(ghosted_id))%ptr, &
              patch%char_curves_thermal_array(patch%cct_id(ghosted_id))%ptr, &
              th_parameter, icct, &
              -grid%nG2A(ghosted_id),PETSC_FALSE,option)
      endif
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! source/sinks
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  allocate(xx(option%nflowdof))
  do
    if (.not.associated(source_sink)) exit
    cur_connection_set => source_sink%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      icct = patch%cct_id(ghosted_id)

      iend = ghosted_id*option%nflowdof
      istart = iend-option%nflowdof+1
      iphase = global_auxvars(ghosted_id)%istate

      select case(source_sink%flow_condition%itype(TH_TEMPERATURE_DOF))
        case (HET_DIRICHLET_BC)
          tsrc1 = source_sink%flow_aux_real_var(TWO_INTEGER,iconn)
        case (DIRICHLET_BC)
          tsrc1 = source_sink%flow_condition%temperature%dataset%rarray(1)
        case (ENERGY_RATE_SS,SCALED_ENERGY_RATE_SS,HET_ENERGY_RATE_SS, &
              ZERO_GRADIENT_BC)
          tsrc1 = xx_loc_p((ghosted_id-1)*option%nflowdof+2)
        case default
          option%io_buffer='Unsupported temperature flow condtion for ' // &
            'a source-sink in TH mode: ' // trim(source_sink%name)
          call PrintErrMsg(option)
      end select

      xx(1) = xx_loc_p(istart)
      xx(2) = tsrc1

      if (option%flow%th_freezing) then
         call THAuxVarComputeFreezing(xx, &
                                      TH_auxvars_ss(sum_connection), &
                                      global_auxvars_ss(sum_connection), &
                                      material_auxvars(ghosted_id), &
                                      iphase, &
                                      patch%saturation_function_array( &
                                        patch%cc_id(ghosted_id))%ptr, &
                                      patch%char_curves_thermal_array( &
                                        patch%cct_id(ghosted_id))%ptr, &
                                      th_parameter, icct, &
                                      -grid%nG2A(ghosted_id),PETSC_FALSE,option)
      else
         call THAuxVarComputeNoFreezing(xx, &
                                      TH_auxvars_ss(sum_connection), &
                                      global_auxvars_ss(sum_connection), &
                                      material_auxvars(ghosted_id), &
                                      iphase, &
                                      patch%characteristic_curves_array( &
                                        patch%cc_id(ghosted_id))%ptr, &
                                      patch%char_curves_thermal_array( &
                                        patch%cct_id(ghosted_id))%ptr, &
                                      th_parameter, icct, &
                                      -grid%nG2A(ghosted_id),PETSC_FALSE,option)
      endif
    enddo
    source_sink => source_sink%next
  enddo
  deallocate(xx)

  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

  patch%aux%TH%auxvars_up_to_date = PETSC_TRUE

end subroutine THUpdateAuxVarsPatch

! ************************************************************************** !

subroutine THInitializeTimestep(realization)
  !
  ! Update data in module prior to time step
  !
  ! Author: ???
  ! Date: 02/20/08
  !

  use Realization_Subsurface_class

  implicit none

  class(realization_subsurface_type) :: realization

  call THUpdateFixedAccumulation(realization)

end subroutine THInitializeTimestep

! ************************************************************************** !

subroutine THUpdateSolution(realization)
  !
  ! Updates data in module after a successful time step
  !
  ! Author: ???
  ! Date: 02/13/08
  !

  use Realization_Subsurface_class
  use Field_module
  use Patch_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(field_type), pointer :: field
  type(patch_type), pointer :: cur_patch

  field => realization%field

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call THUpdateSolutionPatch(realization)
    cur_patch => cur_patch%next
  enddo

end subroutine THUpdateSolution

! ************************************************************************** !

subroutine THUpdateSolutionPatch(realization)
  !
  ! Updates data in module after a successful time
  ! step
  !
  ! Author: Satish Karra, LANL
  ! Date: 12/13/11, 02/28/14
  !


  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Field_module
  use Secondary_Continuum_Aux_module
  use Secondary_Continuum_module

  implicit none

  class(realization_subsurface_type) :: realization
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(th_parameter_type), pointer :: th_parameter
  type(TH_auxvar_type), pointer :: auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(sec_heat_type), pointer :: TH_sec_heat_vars(:)

  PetscInt :: istart, iend
  PetscInt :: local_id, ghosted_id
  ! secondary continuum variables
  PetscReal :: sec_dencpr
  PetscInt :: icct

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  th_parameter => patch%aux%TH%th_parameter
  auxvars => patch%aux%TH%auxvars
  global_auxvars => patch%aux%Global%auxvars

  if (option%use_sc) then
    TH_sec_heat_vars => patch%aux%SC_heat%sec_heat_vars
  endif

  if (realization%option%compute_mass_balance_new) then
    call THUpdateMassBalancePatch(realization)
  endif

  if (option%use_sc) then
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1

      ! secondary rho*c_p same as primary for now
      icct = patch%cct_id(ghosted_id)
      sec_dencpr = th_parameter%dencpr(icct)

      call THSecHeatAuxVarCompute(TH_sec_heat_vars(local_id), &
                            global_auxvars(ghosted_id), &
                            th_parameter%ckwet(icct), &
                            sec_dencpr, &
                            option)

    enddo
  endif


  TH_ts_count = TH_ts_count + 1
  TH_ts_cut_count = 0
  TH_ni_count = 0

end subroutine THUpdateSolutionPatch

! ************************************************************************** !

subroutine THUpdateFixedAccumulation(realization)
  !
  ! Updates the fixed portion of the
  ! accumulation term
  !
  ! Author: ???
  ! Date: 12/10/07
  !

  use Realization_Subsurface_class
  use Patch_module

  class(realization_subsurface_type) :: realization

  type(patch_type), pointer :: cur_patch

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    realization%patch => cur_patch
    call THUpdateFixedAccumPatch(realization)
    cur_patch => cur_patch%next
  enddo

end subroutine THUpdateFixedAccumulation

! ************************************************************************** !

subroutine THUpdateFixedAccumPatch(realization)
  !
  ! Updates the fixed portion of the
  ! accumulation term
  !
  ! Author: ???
  ! Date: 12/10/07
  !

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Secondary_Continuum_Aux_module


  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(TH_auxvar_type), pointer :: TH_auxvars(:)
  type(th_parameter_type), pointer :: th_parameter
  type(sec_heat_type), pointer :: TH_sec_heat_vars(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: ghosted_id, local_id, istart, iend, iphase
  PetscReal, pointer :: xx_p(:)
  PetscReal, pointer :: accum_p(:)
  PetscReal :: vol_frac_prim
  PetscInt :: icct

  PetscErrorCode :: ierr

  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  th_parameter => patch%aux%TH%th_parameter
  TH_auxvars => patch%aux%TH%auxvars
  global_auxvars => patch%aux%Global%auxvars
  TH_sec_heat_vars => patch%aux%SC_heat%sec_heat_vars
  material_auxvars => patch%aux%Material%auxvars

  call VecGetArrayReadF90(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)


  vol_frac_prim = 1.d0

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle

    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    iphase = global_auxvars(ghosted_id)%istate
    icct = patch%cct_id(ghosted_id)

    if (option%flow%th_freezing) then
       call THAuxVarComputeFreezing(xx_p(istart:iend), &
            TH_auxvars(ghosted_id),global_auxvars(ghosted_id), &
            material_auxvars(ghosted_id), &
            iphase, &
            patch%saturation_function_array(patch%cc_id(ghosted_id))%ptr, &
            patch%char_curves_thermal_array(patch%cct_id(ghosted_id))%ptr, &
            th_parameter, icct, &
            grid%nG2A(ghosted_id),PETSC_TRUE,option)
    else
       call THAuxVarComputeNoFreezing(xx_p(istart:iend), &
            TH_auxvars(ghosted_id),global_auxvars(ghosted_id), &
            material_auxvars(ghosted_id), &
            iphase, &
            patch%characteristic_curves_array(patch%cc_id(ghosted_id))%ptr, &
            patch%char_curves_thermal_array(patch%cct_id(ghosted_id))%ptr, &
            th_parameter, icct, &
            grid%nG2A(ghosted_id),PETSC_TRUE,option)
    endif


    if (option%use_sc) then
      vol_frac_prim = TH_sec_heat_vars(local_id)%epsilon
    endif

    global_auxvars(ghosted_id)%istate = iphase
    call THAccumulation(TH_auxvars(ghosted_id),global_auxvars(ghosted_id), &
                        material_auxvars(ghosted_id), &
                        th_parameter%dencpr(patch%cct_id(ghosted_id)), &
                        option,vol_frac_prim,accum_p(istart:iend))
  enddo

  call VecRestoreArrayReadF90(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)

  call VecRestoreArrayF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)

#if 0
   call THNumericalJacobianTest(field%flow_xx,realization)
#endif

end subroutine THUpdateFixedAccumPatch

! ************************************************************************** !

subroutine THNumericalJacobianTest(xx,realization)
  !
  ! Computes the a test numerical jacobian
  !
  ! Author: ???
  ! Date: 12/13/07
  !

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Grid_module
  use Field_module

  implicit none

  Vec :: xx
  class(realization_subsurface_type) :: realization

  Vec :: xx_pert
  Vec :: res
  Vec :: res_pert
  Mat :: A
  PetscViewer :: viewer
  PetscErrorCode :: ierr

  PetscReal :: derivative, perturbation

  PetscReal, pointer :: vec_p(:), vec2_p(:)

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field

  PetscInt :: idof, idof2, icell

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  call VecDuplicate(xx,xx_pert,ierr);CHKERRQ(ierr)
  call VecDuplicate(xx,res,ierr);CHKERRQ(ierr)
  call VecDuplicate(xx,res_pert,ierr);CHKERRQ(ierr)

  call MatCreate(option%mycomm,A,ierr);CHKERRQ(ierr)
  call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,grid%nlmax*option%nflowdof, &
                   grid%nlmax*option%nflowdof,ierr);CHKERRQ(ierr)
  call MatSetType(A,MATAIJ,ierr);CHKERRQ(ierr)
  call MatSetFromOptions(A,ierr);CHKERRQ(ierr)

  call THResidual(PETSC_NULL_SNES,xx,res,realization,ierr)
  call VecGetArrayF90(res,vec2_p,ierr);CHKERRQ(ierr)
  do icell = 1,grid%nlmax
    if (patch%imat(icell) <= 0) cycle
    do idof = (icell-1)*option%nflowdof+1,icell*option%nflowdof
      call VecCopy(xx,xx_pert,ierr);CHKERRQ(ierr)
      call VecGetArrayF90(xx_pert,vec_p,ierr);CHKERRQ(ierr)
      perturbation = vec_p(idof)*perturbation_tolerance
      vec_p(idof) = vec_p(idof)+perturbation
      call VecRestoreArrayF90(xx_pert,vec_p,ierr);CHKERRQ(ierr)
      call THResidual(PETSC_NULL_SNES,xx_pert,res_pert,realization,ierr)
      call VecGetArrayF90(res_pert,vec_p,ierr);CHKERRQ(ierr)
      do idof2 = 1, grid%nlmax*option%nflowdof
        derivative = (vec_p(idof2)-vec2_p(idof2))/perturbation
        if (dabs(derivative) > 1.d-30) then
          call matsetvalue(a,idof2-1,idof-1,derivative,insert_values, &
                           ierr);CHKERRQ(ierr)
        endif
      enddo
      call VecRestoreArrayF90(res_pert,vec_p,ierr);CHKERRQ(ierr)
    enddo
  enddo
  call VecRestoreArrayF90(res,vec2_p,ierr);CHKERRQ(ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call PetscViewerASCIIOpen(option%mycomm,'numerical_jacobian.out',viewer, &
                            ierr);CHKERRQ(ierr)
  call MatView(A,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)

  call MatDestroy(A,ierr);CHKERRQ(ierr)

  call VecDestroy(xx_pert,ierr);CHKERRQ(ierr)
  call VecDestroy(res,ierr);CHKERRQ(ierr)
  call VecDestroy(res_pert,ierr);CHKERRQ(ierr)

end subroutine THNumericalJacobianTest

! ************************************************************************** !

subroutine THAccumDerivative(TH_auxvar,global_auxvar, &
                             material_auxvar, &
                             rock_dencpr, &
                             th_parameter, &
                             icct, &
                             option,sat_func, &
                             characteristic_curves, &
                             thermal_cc, &
                             vol_frac_prim,J)
  !
  ! Computes derivatives of the accumulation
  ! term for the Jacobian
  !
  ! Author: ???
  ! Date: 12/13/07
  !

  use Option_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Thermal_module
  use Saturation_Function_module
  use Material_Aux_module, only : material_auxvar_type
  use EOS_Water_module

  implicit none

  type(TH_auxvar_type) :: TH_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: vol,por,rock_dencpr
  type(th_parameter_type) :: th_parameter
  PetscInt :: icct
  class(characteristic_curves_type) :: characteristic_curves
  class(cc_thermal_type) :: thermal_cc
  type(saturation_function_type) :: sat_func
  PetscReal :: J(option%nflowdof,option%nflowdof)

  PetscReal :: porXvol

  PetscInt :: iphase, ideriv
  type(TH_auxvar_type) :: TH_auxvar_pert
  type(global_auxvar_type) :: global_auxvar_pert
  ! leave as type
  type(material_auxvar_type) :: material_auxvar_pert
  PetscReal :: x(option%nflowdof), x_pert(option%nflowdof), pert
  PetscReal :: res(option%nflowdof), res_pert(option%nflowdof)
  PetscReal :: J_pert(option%nflowdof,option%nflowdof)
  PetscReal :: vol_frac_prim
  PetscReal :: dcompressed_porosity_dp

  PetscReal :: pres, temp
  PetscReal :: sat, dsat_dp, dsat_dT
  PetscReal :: den, dden_dp, dden_dT
  PetscReal :: u, du_dp, du_dT

  ! ice variables
  PetscReal :: sat_g, den_g, mol_g, u_g
  PetscReal :: ddeng_dT, dmolg_dT, dsatg_dp, dsatg_dT, dug_dT
  PetscReal :: sat_i, den_i, u_i
  PetscReal :: dsati_dp, dsati_dT
  PetscReal :: ddeni_dp, ddeni_dT
  PetscReal :: dui_dT


  ! X = {p, T}; R = {R_p, R_T}

  vol = material_auxvar%volume
  pres = global_auxvar%pres(1)
  temp = global_auxvar%temp
  sat = global_auxvar%sat(1)
  den = global_auxvar%den(1)
  dden_dp = TH_auxvar%dden_dp
  dden_dT = TH_auxvar%dden_dT
  dsat_dp = TH_auxvar%dsat_dp
  dsat_dT = TH_auxvar%dsat_dT
  u = TH_auxvar%u
  du_dT = TH_auxvar%du_dT
  du_dp = TH_auxvar%du_dp

  por = material_auxvar%porosity
  dcompressed_porosity_dp = material_auxvar%dporosity_dp

  porXvol = por*vol

  ! d(por*sat*den)/dP * vol
  J(TH_PRESSURE_DOF,TH_PRESSURE_DOF) = (sat*dden_dp + dsat_dp*den)*porXvol + &
    dcompressed_porosity_dp*sat*den*vol

  J(TH_PRESSURE_DOF,TH_TEMPERATURE_DOF) = sat*dden_dT*porXvol
  J(TH_TEMPERATURE_DOF,TH_PRESSURE_DOF) = (dsat_dp*den*u + &
                                           sat*dden_dp*u + &
                                           sat*den*du_dp)*porXvol + &
                                           (den*sat*u - rock_dencpr*temp)* &
                                           vol*dcompressed_porosity_dp
  J(TH_TEMPERATURE_DOF,TH_TEMPERATURE_DOF) = &
            sat*(dden_dT*u + den*du_dT)*porXvol + (1.d0 - por)*vol*rock_dencpr


  if (option%flow%th_freezing) then
     ! SK, 11/17/11
     sat_g    = TH_auxvar%ice%sat_gas
     sat_i    = TH_auxvar%ice%sat_ice
     dsati_dT = TH_auxvar%ice%dsat_ice_dT
     dsati_dp = TH_auxvar%ice%dsat_ice_dp
     dsatg_dp = TH_auxvar%ice%dsat_gas_dp
     dsatg_dT = TH_auxvar%ice%dsat_gas_dT

     u_g      = TH_auxvar%ice%u_gas
     u_i      = TH_auxvar%ice%u_ice
     dug_dT   = TH_auxvar%ice%du_gas_dT
     dui_dT   = TH_auxvar%ice%du_ice_dT

     den_i    = TH_auxvar%ice%den_ice
     den_g    = TH_auxvar%ice%den_gas
     mol_g    = TH_auxvar%ice%mol_gas
     ddeni_dT = TH_auxvar%ice%dden_ice_dT
     ddeni_dp = TH_auxvar%ice%dden_ice_dp
     ddeng_dT = TH_auxvar%ice%dden_gas_dT
     dmolg_dT = TH_auxvar%ice%dmol_gas_dT

     J(TH_PRESSURE_DOF,TH_PRESSURE_DOF) = J(TH_PRESSURE_DOF,TH_PRESSURE_DOF) + &
                                          (dsatg_dp*den_g*mol_g + &
                                           dsati_dp*den_i       + &
                                           sat_i   *ddeni_dp     )*porXvol + &
                                          (sat_g   *den_g*mol_g + &
                                           sat_i   *den_i        )* &
                                          dcompressed_porosity_dp*vol

     J(TH_PRESSURE_DOF,TH_TEMPERATURE_DOF) = &
       J(TH_PRESSURE_DOF,TH_TEMPERATURE_DOF) + &
                            (TH_auxvar%dsat_dT*global_auxvar%den(1) + &
                             dsatg_dT * den_g    * mol_g            + &
                             sat_g    * ddeng_dT * mol_g            + &
                             sat_g    * den_g    * dmolg_dT         + &
                             dsati_dT * den_i                       + &
                             sat_i    * ddeni_dT                    )*porXvol

     J(TH_TEMPERATURE_DOF,TH_PRESSURE_DOF) = &
       J(TH_TEMPERATURE_DOF,TH_PRESSURE_DOF) + &
                     (dsatg_dp * den_g    * u_g + &
                      dsati_dp * den_i    * u_i + &
                      sat_i    * ddeni_dp * u_i )*porXvol + &
                     (sat_g    * den_g    * u_g + &
                      sat_i    * den_i    * u_i )*dcompressed_porosity_dp*vol

     J(TH_TEMPERATURE_DOF,TH_TEMPERATURE_DOF) = &
       J(TH_TEMPERATURE_DOF,TH_TEMPERATURE_DOF) + &
                (TH_auxvar%dsat_dT*global_auxvar%den(1)*TH_auxvar%u + &
                  dsatg_dT * den_g    * u_g                         + &
                  sat_g    * ddeng_dT * u_g                         + &
                  sat_g    * den_g    * dug_dT                      + &
                  dsati_dT * den_i    * u_i                         + &
                  sat_i    * ddeni_dT * u_i                         + &
                  sat_i    * den_i    * dui_dT                      )*porXvol
  endif

  J = J/option%flow_dt
  J(option%nflowdof,:) = vol_frac_prim*J(option%nflowdof,:)

  ! If only solving the energy equation,
  !  - Set jacobian term corresponding to mass-equation to zero, and
  !  - Set off-diagonal jacobian terms to zero.
  if (option%flow%only_energy_eq) then
    J(1,1) = 1.d0
    J(1,2) = 0.d0
    J(2,1) = 0.d0
  endif

  if (option%flow%numerical_derivatives) then
    call GlobalAuxVarInit(global_auxvar_pert,option)
    call MaterialAuxVarInit(material_auxvar_pert,option)

    call THAuxVarCopy(TH_auxvar,TH_auxvar_pert,option)
    call GlobalAuxVarCopy(global_auxvar,global_auxvar_pert,option)
    call MaterialAuxVarCopy(material_auxvar,material_auxvar_pert,option)

    x(1) = global_auxvar%pres(1)
    x(2) = global_auxvar%temp

    call THAccumulation(TH_auxvar,global_auxvar,material_auxvar, &
                         rock_dencpr,option, &
                         vol_frac_prim,res)

    do ideriv = 1,option%nflowdof
      pert = x(ideriv)*perturbation_tolerance
      x_pert = x
      if (option%flow%th_freezing) then
         if (ideriv == 1) then
            if (x_pert(ideriv) < option%flow%reference_pressure) then
               pert = - pert
            endif
            x_pert(ideriv) = x_pert(ideriv) + pert
         endif

         if (ideriv == 2) then
            if (x_pert(ideriv) < 0.d0) then
               pert = - 1.d-8
            else
               pert =  1.d-8
            endif
            x_pert(ideriv) = x_pert(ideriv) + pert
         endif
      else
         x_pert(ideriv) = x_pert(ideriv) + pert
      endif

      if (option%flow%th_freezing) then
         call THAuxVarComputeFreezing(x_pert,TH_auxvar_pert, &
                                 global_auxvar_pert,material_auxvar_pert, &
                                 iphase,sat_func,thermal_cc, &
                                 th_parameter, icct, &
                                 -999,PETSC_TRUE,option)
      else
         call THAuxVarComputeNoFreezing(x_pert,TH_auxvar_pert,&
                              global_auxvar_pert,material_auxvar_pert,&
                              iphase,characteristic_curves,thermal_cc, &
                              th_parameter,icct, &
                              -999,PETSC_TRUE,option)
      endif

      call THAccumulation(TH_auxvar_pert,global_auxvar_pert, &
                          material_auxvar_pert, &
                          rock_dencpr,option,vol_frac_prim, &
                          res_pert)
      J_pert(:,ideriv) = (res_pert(:)-res(:))/pert
    enddo

    J = J_pert
    call GlobalAuxVarStrip(global_auxvar_pert)
  endif

end subroutine THAccumDerivative

! ************************************************************************** !

subroutine THAccumulation(auxvar,global_auxvar, &
                          material_auxvar, &
                          rock_dencpr,option,vol_frac_prim,Res)
  !
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual
  !
  ! Author: ???
  ! Date: 12/13/07
  !

  use Option_module
  use Material_Aux_module, only : material_auxvar_type
  use EOS_Water_module

  implicit none

  type(TH_auxvar_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: Res(1:option%nflowdof)
  PetscReal ::rock_dencpr

  PetscReal :: vol,por
  PetscReal :: porXvol, mol(option%nflowspec), eng
  PetscReal :: vol_frac_prim

  ! ice variables
  PetscReal :: sat_g, den_g, mol_g, u_g
  PetscReal :: sat_i, den_i, u_i

  vol = material_auxvar%volume

  por = material_auxvar%porosity

  ! TechNotes, TH Mode: First term of Equation 8
  porXvol = por*vol

  mol(1) = global_auxvar%sat(1)*global_auxvar%den(1)*porXvol

! TechNotes, TH Mode: First term of Equation 9
  ! rock_dencpr [MJ/m^3 rock-K]
  eng = global_auxvar%sat(1) * &
        global_auxvar%den(1) * &
        auxvar%u * porXvol + &
        (1.d0 - por) * vol * rock_dencpr * global_auxvar%temp

  if (option%flow%th_freezing) then
     ! SK, 11/17/11
     sat_g = auxvar%ice%sat_gas
     sat_i = auxvar%ice%sat_ice
     u_i   = auxvar%ice%u_ice
     den_i = auxvar%ice%den_ice
     den_g = auxvar%ice%den_gas
     mol_g = auxvar%ice%mol_gas
     u_g = auxvar%ice%u_gas
     mol(1) = mol(1) + (sat_g*den_g*mol_g + sat_i*den_i)*porXvol
     eng = eng + (sat_g*den_g*u_g + sat_i*den_i*u_i)*porXvol
  endif

  Res(1:option%nflowdof-1) = mol(:)/option%flow_dt
  Res(option%nflowdof) = vol_frac_prim*eng/option%flow_dt

end subroutine THAccumulation

! ************************************************************************** !

subroutine THFluxDerivative(auxvar_up,global_auxvar_up, &
                            material_auxvar_up, &
                            Dk_up, &
                            icct_up, &
                            auxvar_dn,global_auxvar_dn, &
                            material_auxvar_dn, &
                            Dk_dn, &
                            icct_dn, &
                            area, &
                            dist, upweight, &
                            sir_up,sir_dn, &
                            option,sf_up,sf_dn, &
                            cc_up,cc_dn, &
                            tcc_up,tcc_dn, &
                            Dk_dry_up,Dk_dry_dn, &
                            Dk_ice_up,Dk_ice_dn, &
                            alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                            th_parameter, &
                            Jup,Jdn)

  !
  ! Computes the derivatives of the internal flux terms
  ! for the Jacobian
  !
  ! Author: ???
  ! Date: 12/13/07
  !

  use Option_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Thermal_module
  use Saturation_Function_module
  use Connection_module
  use EOS_Water_module
  use Utility_module

  implicit none

  type(TH_auxvar_type) :: auxvar_up, auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: perm_up, perm_dn
  PetscReal :: Dk_up, Dk_dn
  PetscReal :: Dk_dry_up, Dk_dry_dn
  PetscReal :: Dk_ice_up, Dk_ice_dn
  PetscReal :: alpha_up, alpha_dn
  PetscReal :: alpha_fr_up, alpha_fr_dn
  PetscInt :: icct_up, icct_dn
  PetscReal :: v_darcy, area
  PetscReal :: dist(-1:3)
  type(saturation_function_type), pointer :: sf_up, sf_dn
  class(characteristic_curves_type), pointer :: cc_up, cc_dn
  class(cc_thermal_type), pointer :: tcc_up, tcc_dn
  type(th_parameter_type) :: th_parameter
  PetscReal :: Jup(option%nflowdof,option%nflowdof)
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)

  PetscReal :: por_up, por_dn
  PetscReal :: tor_up, tor_dn
  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: fluxm,fluxe,q
  PetscReal :: uh,ukvr,DK,Dq
  PetscReal :: upweight,density_ave,gravity,dphi

  PetscReal :: dden_ave_dp_up, dden_ave_dp_dn, dden_ave_dT_up, dden_ave_dT_dn
  PetscReal :: dgravity_dden_up, dgravity_dden_dn
  PetscReal :: dphi_dp_up, dphi_dp_dn, dphi_dT_up, dphi_dT_dn
  PetscReal :: dukvr_dp_up, dukvr_dp_dn, dukvr_dT_up, dukvr_dT_dn
  PetscReal :: duh_dp_up, duh_dp_dn, duh_dT_up, duh_dT_dn
  PetscReal :: dq_dp_up, dq_dp_dn, dq_dT_up, dq_dT_dn

  PetscReal :: Dk_eff_up, Dk_eff_dn
  PetscReal :: dk_ds_up, dK_di_up, dk_dT_up
  PetscReal :: dk_ds_dn, dK_di_dn, dk_dT_dn
  PetscReal :: dKe_dT_up, dKe_dp_up
  PetscReal :: dKe_dT_dn, dKe_dp_dn
  PetscReal :: dDk_dT_up, dDk_dT_dn
  PetscReal :: dDk_dp_up, dDk_dp_dn

  PetscInt :: iphase, ideriv
  type(TH_auxvar_type) :: auxvar_pert_up, auxvar_pert_dn
  type(global_auxvar_type) :: global_auxvar_pert_up, global_auxvar_pert_dn
  PetscReal :: x_up(option%nflowdof), x_dn(option%nflowdof)
  PetscReal :: x_pert_up(option%nflowdof), x_pert_dn(option%nflowdof)
  PetscReal :: pert_up, pert_dn
  PetscReal :: res(option%nflowdof)
  PetscReal :: res_pert_up(option%nflowdof)
  PetscReal :: res_pert_dn(option%nflowdof)
  PetscReal :: J_pert_up(option%nflowdof,option%nflowdof)
  PetscReal :: J_pert_dn(option%nflowdof,option%nflowdof)
  type(material_auxvar_type), allocatable :: material_auxvar_pert_dn, &
                                              material_auxvar_pert_up

  ! ice variables
  PetscReal :: Ddiffgas_avg, Ddiffgas_up, Ddiffgas_dn
  PetscReal :: p_g
  PetscReal :: deng_up, deng_dn
  PetscReal :: psat_up, psat_dn
  PetscReal :: molg_up, molg_dn
  PetscReal :: satg_up, satg_dn
  PetscReal :: Diffg_up, Diffg_dn
  PetscReal :: ddeng_dT_up, ddeng_dT_dn
  PetscReal :: dpsat_dT_up, dpsat_dT_dn
  PetscReal :: dmolg_dT_up, dmolg_dT_dn
  PetscReal :: dDiffg_dT_up, dDiffg_dT_dn
  PetscReal :: dDiffg_dp_up, dDiffg_dp_dn
  PetscReal :: dsatg_dp_up, dsatg_dp_dn
  PetscReal :: Diffg_ref, p_ref, T_ref
  PetscErrorCode :: ierr
  PetscReal :: Ke_fr_up,Ke_fr_dn   ! frozen soil Kersten numbers
  PetscReal :: dKe_fr_dT_up, dKe_fr_dT_dn
  PetscReal :: dKe_fr_dp_up, dKe_fr_dp_dn
  PetscReal :: fv_up, fv_dn
  PetscReal :: dfv_dT_up, dfv_dT_dn
  PetscReal :: dfv_dp_up, dfv_dp_dn
  PetscReal :: dmolg_dp_up, dmolg_dp_dn
  PetscBool :: is_flowing

  call ConnectionCalculateDistances(dist,option%gravity,dd_up,dd_dn, &
                                    dist_gravity,upweight)
  call PermeabilityTensorToScalar(material_auxvar_up,dist,perm_up)
  call PermeabilityTensorToScalar(material_auxvar_dn,dist,perm_dn)

  por_up = material_auxvar_up%porosity_base
  por_dn = material_auxvar_dn%porosity_base

  tor_up = material_auxvar_up%tortuosity
  tor_dn = material_auxvar_dn%tortuosity

  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)

  fluxm = 0.D0
  fluxe = 0.D0
  v_darcy = 0.D0

  Jup = 0.d0
  Jdn = 0.d0

  dden_ave_dp_up = 0.d0
  dden_ave_dT_up = 0.d0
  dden_ave_dp_dn = 0.d0
  dden_ave_dT_dn = 0.d0
  dgravity_dden_up = 0.d0
  dgravity_dden_dn = 0.d0
  dphi_dp_up = 0.d0
  dphi_dp_dn = 0.d0
  dphi_dT_up = 0.d0
  dphi_dT_dn = 0.d0
  dukvr_dp_up = 0.d0
  dukvr_dp_dn = 0.d0
  dukvr_dT_up = 0.d0
  dukvr_dT_dn = 0.d0
  duh_dp_up = 0.d0
  duh_dp_dn = 0.d0
  duh_dT_up = 0.d0
  duh_dT_dn = 0.d0
  dq_dp_up = 0.d0
  dq_dp_dn = 0.d0
  dq_dT_up = 0.d0
  dq_dT_dn = 0.d0
  dDk_dT_up = 0.d0
  dDk_dT_dn = 0.d0
  dDk_dp_up = 0.d0
  dDk_dp_dn = 0.d0

  if (option%flow%th_freezing) then
    dfv_dT_up = 0.d0
    dfv_dT_dn = 0.d0
    dfv_dp_up = 0.d0
    dfv_dp_dn = 0.d0
    dmolg_dp_up = 0.d0
    dmolg_dp_dn = 0.d0
    dmolg_dT_up = 0.d0
    dmolg_dT_dn = 0.d0
  endif

  ! Flow term
  is_flowing = PETSC_FALSE

  if (option%flow%th_freezing) then
    if (global_auxvar_up%sat(1) > sir_up .or. &
        global_auxvar_dn%sat(1) > sir_dn) then
      is_flowing = PETSC_TRUE
    endif
  else
    if (auxvar_up%kvr > eps .or. &
        auxvar_dn%kvr > eps) then
      is_flowing = PETSC_TRUE
    endif
  endif

  if (is_flowing) then
    if (global_auxvar_up%sat(1) <eps) then
      upweight=0.d0
    else if (global_auxvar_dn%sat(1) <eps) then
      upweight=1.d0
    endif
    density_ave = upweight*global_auxvar_up%den(1)+ &
                  (1.D0-upweight)*global_auxvar_dn%den(1)
    dden_ave_dp_up = upweight*auxvar_up%dden_dp
    dden_ave_dp_dn = (1.D0-upweight)*auxvar_dn%dden_dp
    dden_ave_dT_up = upweight*auxvar_up%dden_dT
    dden_ave_dT_dn = (1.D0-upweight)*auxvar_dn%dden_dT

    gravity = (upweight*global_auxvar_up%den(1)*auxvar_up%avgmw + &
              (1.D0-upweight)*global_auxvar_dn%den(1)*auxvar_dn%avgmw) &
              * dist_gravity
    dgravity_dden_up = upweight*auxvar_up%avgmw*dist_gravity
    dgravity_dden_dn = (1.d0-upweight)*auxvar_dn%avgmw*dist_gravity

    if (th_ice_model /= DALL_AMICO) then
      dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity
      dphi_dp_up = 1.d0 + dgravity_dden_up*auxvar_up%dden_dp
      dphi_dp_dn = -1.d0 + dgravity_dden_dn*auxvar_dn%dden_dp
      dphi_dT_up = dgravity_dden_up*auxvar_up%dden_dT
      dphi_dT_dn = dgravity_dden_dn*auxvar_dn%dden_dT
    else
      dphi = auxvar_up%ice%pres_fh2o - auxvar_dn%ice%pres_fh2o + gravity
      dphi_dp_up =  auxvar_up%ice%dpres_fh2o_dp + &
                    dgravity_dden_up*auxvar_up%dden_dp
      dphi_dp_dn = -auxvar_dn%ice%dpres_fh2o_dp + &
                    dgravity_dden_dn*auxvar_dn%dden_dp
      dphi_dT_up =  auxvar_up%ice%dpres_fh2o_dT + &
                    dgravity_dden_up*auxvar_up%dden_dT
      dphi_dT_dn = -auxvar_dn%ice%dpres_fh2o_dT + &
                    dgravity_dden_dn*auxvar_dn%dden_dT
    endif

    if (dphi>=0.D0) then
      ukvr = auxvar_up%kvr
      dukvr_dp_up = auxvar_up%dkvr_dp
      dukvr_dT_up = auxvar_up%dkvr_dT

      uh = auxvar_up%h
      duh_dp_up = auxvar_up%dh_dp
      duh_dT_up = auxvar_up%dh_dT
    else
      ukvr = auxvar_dn%kvr
      dukvr_dp_dn = auxvar_dn%dkvr_dp
      dukvr_dT_dn = auxvar_dn%dkvr_dT

      uh = auxvar_dn%h
      duh_dp_dn = auxvar_dn%dh_dp
      duh_dT_dn = auxvar_dn%dh_dT
    endif

    call InterfaceApprox(auxvar_up%kvr, auxvar_dn%kvr, &
                         auxvar_up%dkvr_dp, auxvar_dn%dkvr_dp, &
                         dphi, &
                         option%flow%rel_perm_aveg, &
                         ukvr, dukvr_dp_up, dukvr_dp_dn)

    call InterfaceApprox(auxvar_up%kvr, auxvar_dn%kvr, &
                         auxvar_up%dkvr_dT, auxvar_dn%dkvr_dT, &
                         dphi, &
                         option%flow%rel_perm_aveg, &
                         ukvr, dukvr_dT_up, dukvr_dT_dn)

    if (ukvr>floweps) then
      v_darcy= Dq * ukvr * dphi

      q = v_darcy * area


      dq_dT_up = Dq*(dukvr_dT_up*dphi+ukvr*dphi_dT_up)*area
      dq_dT_dn = Dq*(dukvr_dT_dn*dphi+ukvr*dphi_dT_dn)*area

      dq_dp_up = Dq*(dukvr_dp_up*dphi+ukvr*dphi_dp_up)*area
      dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area



      ! If only solving the energy equation, ensure Jup(2,2) & Jdn(2,2)
      ! have no contribution from the mass equation
      if (option%flow%only_energy_eq) then
         v_darcy = 0.d0
         q = 0.d0
         dq_dT_up = 0.d0
         dq_dT_dn = 0.d0
      endif

      Jup(1,1) = (dq_dp_up*density_ave+q*dden_ave_dp_up)
      Jup(1,2) = (dq_dT_up*density_ave+q*dden_ave_dT_up)

      Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)
      Jdn(1,2) = (dq_dT_dn*density_ave+q*dden_ave_dT_dn)


      ! based on flux = q*density_ave*uh
      Jup(option%nflowdof,1) = (dq_dp_up*density_ave+q*dden_ave_dp_up)*uh+ &
                               q*density_ave*duh_dp_up
      Jup(option%nflowdof,2) = (dq_dT_up*density_ave+q*dden_ave_dT_up)*uh+ &
                               q*density_ave*duh_dT_up
      Jdn(option%nflowdof,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)*uh+ &
                               q*density_ave*duh_dp_dn
      Jdn(option%nflowdof,2) = (dq_dT_dn*density_ave+q*dden_ave_dT_dn)*uh+ &
                               q*density_ave*duh_dT_dn

    endif
  endif

  if (option%flow%th_freezing) then
    ! Added by Satish Karra, updated 11/11/11
    satg_up = auxvar_up%ice%sat_gas
    satg_dn = auxvar_dn%ice%sat_gas
    if ((satg_up > eps) .and. (satg_dn > eps)) then
      p_g = option%flow%reference_pressure  ! set to reference pressure
      deng_up = p_g/(IDEAL_GAS_CONSTANT* &
                (global_auxvar_up%temp + 273.15d0))*1.d-3
      deng_dn = p_g/(IDEAL_GAS_CONSTANT* &
                (global_auxvar_dn%temp + 273.15d0))*1.d-3

      Diffg_ref = 2.13D-5 ! Reference diffusivity, need to read from input file
      p_ref = 1.01325d5   ! in Pa
      T_ref = 25.d0       ! in deg C

      Diffg_up = Diffg_ref*(p_ref/p_g)*((global_auxvar_up%temp + 273.15d0)/ &
           (T_ref + 273.15d0))**(1.8)
      Diffg_dn = Diffg_ref*(p_ref/p_g)*((global_auxvar_dn%temp + 273.15d0)/ &
           (T_ref + 273.15d0))**(1.8)

      Ddiffgas_up = por_up*tor_up*satg_up*deng_up*Diffg_up
      Ddiffgas_dn = por_dn*tor_dn*satg_dn*deng_dn*Diffg_dn
      call EOSWaterSaturationPressure(global_auxvar_up%temp, psat_up, &
                                      dpsat_dT_up, ierr)
      call EOSWaterSaturationPressure(global_auxvar_dn%temp, psat_dn, &
                                      dpsat_dT_dn, ierr)

      ! vapor pressure lowering due to capillary pressure
      fv_up = exp(-auxvar_up%pc/(global_auxvar_up%den(1)* &
           IDEAL_GAS_CONSTANT*(global_auxvar_up%temp + 273.15d0)))
      fv_dn = exp(-auxvar_dn%pc/(global_auxvar_dn%den(1)* &
           IDEAL_GAS_CONSTANT*(global_auxvar_dn%temp + 273.15d0)))

      molg_up = psat_up*fv_up/p_g
      molg_dn = psat_dn*fv_dn/p_g

      dfv_dT_up = fv_up*(auxvar_up%pc/IDEAL_GAS_CONSTANT/ &
           (global_auxvar_up%den(1)* &
           (global_auxvar_up%temp + 273.15d0))**2)* &
           (auxvar_up%dden_dT*(global_auxvar_up%temp + 273.15d0) &
           + global_auxvar_up%den(1))
      dfv_dT_dn = fv_dn*(auxvar_dn%pc/IDEAL_GAS_CONSTANT/ &
           (global_auxvar_dn%den(1)* &
           (global_auxvar_dn%temp + 273.15d0))**2)* &
           (auxvar_dn%dden_dT*(global_auxvar_dn%temp + 273.15d0) &
           + global_auxvar_dn%den(1))

      dfv_dp_up = fv_up*(auxvar_up%pc/IDEAL_GAS_CONSTANT/ &
           (global_auxvar_up%den(1))**2/ &
           (global_auxvar_up%temp + 273.15d0)*auxvar_up%dden_dp &
           + 1.d0/IDEAL_GAS_CONSTANT/global_auxvar_up%den(1)/ &
           (global_auxvar_up%temp + 273.15d0))
      dfv_dp_dn = fv_dn*(auxvar_dn%pc/IDEAL_GAS_CONSTANT/ &
           (global_auxvar_dn%den(1))**2/ &
           (global_auxvar_dn%temp + 273.15d0)*auxvar_dn%dden_dp &
           + 1.d0/IDEAL_GAS_CONSTANT/global_auxvar_dn%den(1)/ &
           (global_auxvar_dn%temp + 273.15d0))

      dmolg_dT_up = (1/p_g)*dpsat_dT_up*fv_up + psat_up/p_g*dfv_dT_up
      dmolg_dT_dn = (1/p_g)*dpsat_dT_dn*fv_dn + psat_dn/p_g*dfv_dT_dn

      dmolg_dp_up = psat_up/p_g*dfv_dp_up
      dmolg_dp_dn = psat_dn/p_g*dfv_dp_dn

      ddeng_dT_up = - p_g/(IDEAL_GAS_CONSTANT*(global_auxvar_up%temp + &
           273.15d0)**2)*1.d-3
      ddeng_dT_dn = - p_g/(IDEAL_GAS_CONSTANT*(global_auxvar_dn%temp + &
           273.15d0)**2)*1.d-3

      dDiffg_dT_up = 1.8*Diffg_up/(global_auxvar_up%temp + 273.15d0)
      dDiffg_dT_dn = 1.8*Diffg_dn/(global_auxvar_dn%temp + 273.15d0)

      dDiffg_dp_up = 0.d0
      dDiffg_dp_dn = 0.d0

      dsatg_dp_up = auxvar_up%ice%dsat_gas_dp
      dsatg_dp_dn = auxvar_dn%ice%dsat_gas_dp

      if (molg_up > molg_dn) then
         upweight = 0.d0
      else
         upweight = 1.d0
      endif

      Ddiffgas_avg = upweight*Ddiffgas_up + (1.D0 - upweight)*Ddiffgas_dn

#ifndef NO_VAPOR_DIFFUION
      Jup(1,1) = Jup(1,1) + &
           (upweight*por_up*tor_up*deng_up*(Diffg_up*dsatg_dp_up &
           + satg_up*dDiffg_dp_up)* &
           (molg_up - molg_dn) + Ddiffgas_up*dmolg_dp_up)/ &
           (dd_up + dd_dn)*area

      Jup(1,2) = Jup(1,2) + (upweight*por_up*tor_up*satg_up*(Diffg_up* &
           ddeng_dT_up + deng_up*dDiffg_dT_up)*(molg_up - molg_dn) &
           + Ddiffgas_avg*dmolg_dT_up)/(dd_up + dd_dn)*area

      Jdn(1,1) = Jdn(1,1) + ((1.D0 - upweight)*por_dn*tor_dn*deng_dn* &
           (Diffg_dn*dsatg_dp_dn + satg_dn*dDiffg_dp_dn)* &
           (molg_up - molg_dn) + Ddiffgas_avg*(-dmolg_dp_dn))/ &
           (dd_up + dd_dn)*area

      Jdn(1,2) = Jdn(1,2) + &
           ((1.D0 - upweight)*por_dn*tor_dn*satg_dn*(Diffg_dn* &
           ddeng_dT_dn + deng_dn*dDiffg_dp_dn)*(molg_up - molg_dn) &
           + Ddiffgas_avg*(-dmolg_dT_dn))/(dd_up + dd_dn)*area
#endif
   endif

  endif


    dKe_dp_up = auxvar_up%dKe_dp
    dKe_dp_dn = auxvar_dn%dKe_dp

    dKe_dT_up = auxvar_up%dKe_dT
    dKe_dT_dn = auxvar_dn%dKe_dT

  if (option%flow%th_freezing) then

    call tcc_up%thermal_conductivity_function% &
         TCondTensorToScalar(dist,option)
    call tcc_up%thermal_conductivity_function%CalculateFTCond( &
         global_auxvar_up%sat(1),auxvar_up%ice%sat_ice,global_auxvar_up%temp, &
         material_auxvar_up%porosity,Dk_eff_up,dk_ds_up,dK_di_up,dk_dT_up,option)

    call tcc_dn%thermal_conductivity_function% &
         TCondTensorToScalar(dist,option)
    call tcc_dn%thermal_conductivity_function%CalculateFTCond( &
         global_auxvar_dn%sat(1),auxvar_dn%ice%sat_ice,global_auxvar_dn%temp, &
         material_auxvar_dn%porosity,Dk_eff_dn,dk_ds_dn,dK_di_dn,dk_dT_dn,option)

    Ke_fr_up = auxvar_up%ice%Ke_fr
    Ke_fr_dn = auxvar_dn%ice%Ke_fr

    dKe_fr_dT_up = auxvar_up%ice%dKe_fr_dT
    dKe_fr_dT_dn = auxvar_dn%ice%dKe_fr_dT

    dKe_fr_dp_up = auxvar_up%ice%dKe_fr_dp
    dKe_fr_dp_dn = auxvar_dn%ice%dKe_fr_dp

  else

    call tcc_up%thermal_conductivity_function% &
         TCondTensorToScalar(dist,option)
    call tcc_up%thermal_conductivity_function%CalculateTCond( &
         global_auxvar_up%sat(1),global_auxvar_up%temp, &
         material_auxvar_up%porosity,Dk_eff_up,dk_ds_up,dk_dT_up,option)

    call tcc_dn%thermal_conductivity_function% &
         TCondTensorToScalar(dist,option)
    call tcc_dn%thermal_conductivity_function%CalculateTCond( &
         global_auxvar_dn%sat(1),global_auxvar_dn%temp, &
         material_auxvar_dn%porosity,Dk_eff_dn,dk_ds_dn,dk_dT_dn,option)

  endif

  Dk = (Dk_eff_up * Dk_eff_dn) / (dd_dn*Dk_eff_up + dd_up*Dk_eff_dn)

  if (option%flow%th_freezing) then

    dDk_dT_up = Dk**2/Dk_eff_up**2*dd_up*(Dk_up*dKe_dT_up + &
        Dk_ice_up*dKe_fr_dT_up + (- dKe_dT_up - dKe_fr_dT_up)* &
        Dk_dry_up)
    dDk_dT_dn = Dk**2/Dk_eff_dn**2*dd_dn*(Dk_dn*dKe_dT_dn + &
        Dk_ice_dn*dKe_fr_dT_dn + (- dKe_dT_dn - dKe_fr_dT_dn)* &
        Dk_dry_dn)

    dDk_dp_up = Dk**2/Dk_eff_up**2*dd_up*(Dk_up*dKe_dp_up + &
        Dk_ice_up*dKe_fr_dp_up + (- dKe_dp_up - dKe_fr_dp_up)* &
        Dk_dry_up)

    dDk_dp_dn = Dk**2/Dk_eff_dn**2*dd_dn*(Dk_dn*dKe_dp_dn + &
        Dk_ice_dn*dKe_fr_dp_dn + (- dKe_dp_dn - dKe_fr_dp_dn)* &
        Dk_dry_dn)

  else

    dDk_dT_up = Dk**2/Dk_eff_up**2*dd_up*(Dk_up - Dk_dry_up)*dKe_dT_up
    dDk_dT_dn = Dk**2/Dk_eff_dn**2*dd_dn*(Dk_dn - Dk_dry_dn)*dKe_dT_dn

    dDk_dp_up = Dk**2/Dk_eff_up**2*dd_up*(Dk_up - Dk_dry_up)*dKe_dp_up
    dDk_dp_dn = Dk**2/Dk_eff_dn**2*dd_dn*(Dk_dn - Dk_dry_dn)*dKe_dp_dn

  endif

  !  cond = Dk*area*(global_auxvar_up%temp-global_auxvar_dn%temp)
  Jup(option%nflowdof,1) = Jup(option%nflowdof,1) + &
                           area*(global_auxvar_up%temp - &
                           global_auxvar_dn%temp)*dDk_dp_up
  Jdn(option%nflowdof,1) = Jdn(option%nflowdof,1) + &
                           area*(global_auxvar_up%temp - &
                           global_auxvar_dn%temp)*dDk_dp_dn

  Jup(option%nflowdof,2) = Jup(option%nflowdof,2) + Dk*area + &
                           area*(global_auxvar_up%temp - &
                           global_auxvar_dn%temp)*dDk_dT_up
  Jdn(option%nflowdof,2) = Jdn(option%nflowdof,2) + Dk*area*(-1.d0) + &
                           area*(global_auxvar_up%temp - &
                           global_auxvar_dn%temp)*dDk_dT_dn

  ! If only solving the energy equation,
  !  - Set jacobian term corresponding to mass-equation to zero, and
  !  - Set off-diagonal jacobian terms to zero.
  if (option%flow%only_energy_eq) then
    Jup(1,1) = 0.d0
    Jup(1,2) = 0.d0
    Jup(option%nflowdof,1) = 0.d0

    Jdn(1,1) = 0.d0
    Jdn(1,2) = 0.d0
    Jdn(option%nflowdof,1) = 0.d0

  endif

  ! note: Res is the flux contribution, for node up J = J + Jup
  !                                              dn J = J - Jdn

  if (option%flow%numerical_derivatives) then
    call THAuxVarCopy(auxvar_up,auxvar_pert_up,option)
    call THAuxVarCopy(auxvar_dn,auxvar_pert_dn,option)

    call GlobalAuxVarInit(global_auxvar_pert_up,option)
    call GlobalAuxVarInit(global_auxvar_pert_dn,option)
    call GlobalAuxVarCopy(global_auxvar_up,global_auxvar_pert_up,option)
    call GlobalAuxVarCopy(global_auxvar_dn,global_auxvar_pert_dn,option)

    allocate(material_auxvar_pert_up,material_auxvar_pert_dn)
    call MaterialAuxVarInit(material_auxvar_pert_up,option)
    call MaterialAuxVarInit(material_auxvar_pert_dn,option)
    call MaterialAuxVarCopy(material_auxvar_up,material_auxvar_pert_up,option)
    call MaterialAuxVarCopy(material_auxvar_dn,material_auxvar_pert_dn,option)

    x_up(1) = global_auxvar_up%pres(1)
    x_up(2) = global_auxvar_up%temp
    x_dn(1) = global_auxvar_dn%pres(1)
    x_dn(2) = global_auxvar_dn%temp

    call THFlux( &
      auxvar_up,global_auxvar_up, &
      material_auxvar_up, &
      Dk_up, tcc_up, &
      auxvar_dn,global_auxvar_dn, &
      material_auxvar_dn, &
      Dk_dn, tcc_dn, &
      area, &
      dist, upweight, &
      sir_up, sir_dn, &
      option,v_darcy,Dk_dry_up,Dk_dry_dn, &
      Dk_ice_up,Dk_ice_dn, &
      alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
      res)

    do ideriv = 1,option%nflowdof
      pert_up = x_up(ideriv)*perturbation_tolerance
      pert_dn = x_dn(ideriv)*perturbation_tolerance
      x_pert_up = x_up
      x_pert_dn = x_dn

      if (option%flow%th_freezing) then
        if (ideriv == 1) then
          if (x_pert_up(ideriv) < option%flow%reference_pressure) then
            pert_up = - pert_up
          endif
          x_pert_up(ideriv) = x_pert_up(ideriv) + pert_up

          if (x_pert_dn(ideriv) < option%flow%reference_pressure) then
            pert_dn = - pert_dn
          endif
          x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
        endif

        if (ideriv == 2) then
          if (x_pert_up(ideriv) < 0.d0) then
            pert_up = - 1.d-5
          else
            pert_up = 1.d-5
          endif
          x_pert_up(ideriv) = x_pert_up(ideriv) + pert_up

          if (x_pert_dn(ideriv) < 0.d0) then
            pert_dn = - 1.d-5
          else
            pert_dn = 1.d-5
          endif
          x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
        endif

      else
         x_pert_up(ideriv) = x_pert_up(ideriv) + pert_up
         x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn

      endif

      if (option%flow%th_freezing) then
        call THAuxVarComputeFreezing(x_pert_up,auxvar_pert_up, &
             global_auxvar_pert_up, material_auxvar_pert_up, &
             iphase,sf_up, tcc_up, &
             th_parameter,icct_up, &
             -999,PETSC_TRUE,option)
        call THAuxVarComputeFreezing(x_pert_dn,auxvar_pert_dn, &
             global_auxvar_pert_dn, material_auxvar_pert_up, &
             iphase,sf_dn, tcc_dn, &
             th_parameter,icct_up, &
             -999,PETSC_TRUE,option)
      else
        call THAuxVarComputeNoFreezing(x_pert_up,auxvar_pert_up, &
             global_auxvar_pert_up, material_auxvar_pert_up, &
             iphase,cc_up,tcc_up, &
             th_parameter,icct_up, &
             -999,PETSC_TRUE,option)
        call THAuxVarComputeNoFreezing(x_pert_dn,auxvar_pert_dn, &
             global_auxvar_pert_dn,material_auxvar_pert_dn, &
             iphase,cc_dn,tcc_dn, &
             th_parameter,icct_dn, &
             -999,PETSC_TRUE,option)
      endif

      call THFlux(auxvar_pert_up,global_auxvar_pert_up, &
                   material_auxvar_pert_up, &
                   Dk_up, tcc_up, &
                   auxvar_dn,global_auxvar_dn, &
                   material_auxvar_pert_dn, &
                   Dk_dn, tcc_dn, &
                   area, &
                   dist, upweight, &
                   sir_up,sir_dn, &
                   option,v_darcy,Dk_dry_up, &
                   Dk_dry_dn,Dk_ice_up,Dk_ice_dn, &
                   alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                   res_pert_up)
      call THFlux(auxvar_up,global_auxvar_up, &
                   material_auxvar_pert_up, &
                   Dk_up, tcc_up, &
                   auxvar_pert_dn,global_auxvar_pert_dn, &
                   material_auxvar_pert_dn, &
                   Dk_dn, tcc_dn, &
                   area, &
                   dist, upweight, &
                   sir_up,sir_dn, &
                   option,v_darcy,Dk_dry_up, &
                   Dk_dry_dn,Dk_ice_up,Dk_ice_dn, &
                   alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                   res_pert_dn)

      J_pert_up(:,ideriv) = (res_pert_up(:)-res(:))/pert_up
      J_pert_dn(:,ideriv) = (res_pert_dn(:)-res(:))/pert_dn
    enddo

    Jup = J_pert_up
    Jdn = J_pert_dn
    call GlobalAuxVarStrip(global_auxvar_pert_up)
    call GlobalAuxVarStrip(global_auxvar_pert_dn)
    call MaterialAuxVarStrip(material_auxvar_pert_up)
    call MaterialAuxVarStrip(material_auxvar_pert_dn)
  endif

end subroutine THFluxDerivative

! ************************************************************************** !
subroutine THFlux(auxvar_up,global_auxvar_up, &
                  material_auxvar_up, &
                  Dk_up, tcc_up, &
                  auxvar_dn,global_auxvar_dn, &
                  material_auxvar_dn, &
                  Dk_dn, tcc_dn, &
                  area, &
                  dist, &
                  upweight, &
                  sir_up,sir_dn, &
                  option,v_darcy,Dk_dry_up, &
                  Dk_dry_dn,Dk_ice_up,Dk_ice_dn, &
                  alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                  Res)
  !
  ! Computes the internal flux terms for the residual
  !
  ! Author: ???
  ! Date: 12/13/07
  !

  use Option_module
  use Connection_module
  use EOS_Water_module
  use Utility_module
  use Characteristic_Curves_Thermal_module

  implicit none

  type(TH_auxvar_type) :: auxvar_up, auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  class(cc_thermal_type), pointer :: tcc_up, tcc_dn
  type(option_type) :: option
  PetscReal :: sir_up, sir_dn
  PetscReal :: dd_up, dd_dn
  PetscReal :: Dk_up, Dk_dn
  PetscReal :: Dk_dry_up, Dk_dry_dn
  PetscReal :: Dk_ice_up, Dk_ice_dn
  PetscReal :: alpha_up, alpha_dn
  PetscReal :: alpha_fr_up, alpha_fr_dn
  PetscReal :: Dk_eff_up, Dk_eff_dn
  PetscReal :: dk_ds_up, dK_di_up, dk_dT_up
  PetscReal :: dk_ds_dn, dK_di_dn, dk_dT_dn
  PetscReal :: v_darcy,area
  PetscReal :: Res(1:option%nflowdof)
  PetscReal :: dist(-1:3)
  PetscReal :: fluxm,fluxe,q
  PetscReal :: uh,ukvr,DK,Dq
  PetscReal :: upweight,density_ave,cond,gravity,dphi

  PetscReal :: por_up, por_dn
  PetscReal :: tor_up, tor_dn
  PetscReal :: perm_up, perm_dn

  ! ice variables
  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: Ddiffgas_avg, Ddiffgas_up, Ddiffgas_dn
  PetscReal :: p_g
  PetscReal :: deng_up, deng_dn
  PetscReal :: psat_up, psat_dn
  PetscReal :: molg_up, molg_dn
  PetscReal :: satg_up, satg_dn
  PetscReal :: Diffg_up, Diffg_dn
  PetscReal :: Diffg_ref, p_ref, T_ref
  PetscErrorCode :: ierr
  PetscReal :: Ke_fr_up,Ke_fr_dn   ! frozen soil Kersten numbers
  PetscReal :: fv_up, fv_dn
  PetscBool :: is_flowing

  call ConnectionCalculateDistances(dist,option%gravity,dd_up,dd_dn, &
                                    dist_gravity,upweight)
  call PermeabilityTensorToScalar(material_auxvar_up,dist,perm_up)
  call PermeabilityTensorToScalar(material_auxvar_dn,dist,perm_dn)

  por_up = material_auxvar_up%porosity_base
  por_dn = material_auxvar_dn%porosity_base

  tor_up = material_auxvar_up%tortuosity
  tor_dn = material_auxvar_dn%tortuosity

  Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)

  fluxm = 0.D0
  fluxe = 0.D0
  v_darcy = 0.D0

  ! Flow term
  is_flowing = PETSC_FALSE

  if (option%flow%th_freezing) then
    if (global_auxvar_up%sat(1) > sir_up .or.  &
        global_auxvar_dn%sat(1) > sir_dn) then
      is_flowing = PETSC_TRUE
    endif
  else
    if (auxvar_up%kvr > eps .or. &
        auxvar_dn%kvr > eps) then
      is_flowing = PETSC_TRUE
    endif
  endif

  if (is_flowing) then
    if (global_auxvar_up%sat(1) < eps) then
      upweight=0.d0
    else if (global_auxvar_dn%sat(1) < eps) then
      upweight=1.d0
    endif
    density_ave = upweight*global_auxvar_up%den(1)+(1.D0-upweight)* &
                  global_auxvar_dn%den(1)

    gravity = (upweight*global_auxvar_up%den(1)*auxvar_up%avgmw + &
              (1.D0-upweight)*global_auxvar_dn%den(1)*auxvar_dn%avgmw) &
              * dist_gravity

    if (th_ice_model /= DALL_AMICO) then
      dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity
    else
      dphi = auxvar_up%ice%pres_fh2o - auxvar_dn%ice%pres_fh2o + gravity
    endif

    if (dphi >= 0.D0) then
      ukvr = auxvar_up%kvr
      uh = auxvar_up%h
    else
      ukvr = auxvar_dn%kvr
      uh = auxvar_dn%h
    endif

    call InterfaceApprox(auxvar_up%kvr, auxvar_dn%kvr, dphi, &
                         option%flow%rel_perm_aveg, ukvr)

    if (ukvr > floweps) then
      v_darcy = Dq * ukvr * dphi

      ! If only solving the energy equation, ensure Res(2) has no
      ! contribution from mass equation by setting darcy velocity
      ! to be zero
      if (option%flow%only_energy_eq) v_darcy = 0.d0

      q = v_darcy * area

      fluxm = fluxm + q*density_ave
      fluxe = fluxe + q*density_ave*uh
    endif
  endif


  if (option%flow%th_freezing) then
    ! Added by Satish Karra, 10/24/11
    satg_up = auxvar_up%ice%sat_gas
    satg_dn = auxvar_dn%ice%sat_gas
    if ((satg_up > eps) .and. (satg_dn > eps)) then
      p_g = option%flow%reference_pressure ! set to reference pressure
      deng_up = p_g/(IDEAL_GAS_CONSTANT* &
                (global_auxvar_up%temp + 273.15d0))*1.d-3
      deng_dn = p_g/(IDEAL_GAS_CONSTANT* &
                (global_auxvar_dn%temp + 273.15d0))*1.d-3

      Diffg_ref = 2.13D-5 ! Reference diffusivity, need to read from input file
      p_ref = 1.01325d5 ! in Pa
      T_ref = 25.d0 ! in deg C

      Diffg_up = Diffg_ref*(p_ref/p_g)*((global_auxvar_up%temp + 273.15d0)/ &
           (T_ref + 273.15d0))**(1.8)
      Diffg_dn = Diffg_ref*(p_ref/p_g)*((global_auxvar_dn%temp + 273.15d0)/ &
           (T_ref + 273.15d0))**(1.8)

      Ddiffgas_up = por_up*tor_up*satg_up*deng_up*Diffg_up
      Ddiffgas_dn = por_dn*tor_dn*satg_dn*deng_dn*Diffg_dn
      call EOSWaterSaturationPressure(global_auxvar_up%temp, psat_up, ierr)
      call EOSWaterSaturationPressure(global_auxvar_dn%temp, psat_dn, ierr)

      ! vapor pressure lowering due to capillary pressure
      fv_up = exp(-auxvar_up%pc/(global_auxvar_up%den(1)* &
           IDEAL_GAS_CONSTANT*(global_auxvar_up%temp + 273.15d0)))
      fv_dn = exp(-auxvar_dn%pc/(global_auxvar_dn%den(1)* &
           IDEAL_GAS_CONSTANT*(global_auxvar_dn%temp + 273.15d0)))

      molg_up = psat_up*fv_up/p_g
      molg_dn = psat_dn*fv_dn/p_g

      if (molg_up > molg_dn) then
        upweight = 0.d0
      else
        upweight = 1.d0
      endif

      Ddiffgas_avg = upweight*Ddiffgas_up + (1.D0 - upweight)*Ddiffgas_dn
#ifndef NO_VAPOR_DIFFUSION
      fluxm = fluxm + Ddiffgas_avg*area*(molg_up - molg_dn)/ &
           (dd_up + dd_dn)
#endif

    endif

  endif

  if (option%flow%th_freezing) then

    Ke_fr_up = auxvar_up%ice%Ke_fr
    Ke_fr_dn = auxvar_dn%ice%Ke_fr

    call tcc_up%thermal_conductivity_function% &
         TCondTensorToScalar(dist,option)
    call tcc_up%thermal_conductivity_function%CalculateFTCond( &
         global_auxvar_up%sat(1),auxvar_up%ice%sat_ice,global_auxvar_up%temp, &
         material_auxvar_up%porosity,Dk_eff_up,dk_ds_up,dK_di_up,dk_dT_up,option)

    call tcc_dn%thermal_conductivity_function% &
         TCondTensorToScalar(dist,option)
    call tcc_dn%thermal_conductivity_function%CalculateFTCond( &
         global_auxvar_dn%sat(1),auxvar_dn%ice%sat_ice,global_auxvar_dn%temp, &
         material_auxvar_dn%porosity,Dk_eff_dn,dk_ds_dn,dK_di_dn,dk_dT_dn,option)

  else

    call tcc_up%thermal_conductivity_function% &
         TCondTensorToScalar(dist,option)
    call tcc_up%thermal_conductivity_function%CalculateTCond( &
         global_auxvar_up%sat(1),global_auxvar_up%temp,material_auxvar_up%porosity, &
         Dk_eff_up,dk_ds_up,dk_dT_up,option)

    call tcc_dn%thermal_conductivity_function% &
         TCondTensorToScalar(dist,option)
    call tcc_dn%thermal_conductivity_function%CalculateTCond( &
         global_auxvar_dn%sat(1),global_auxvar_dn%temp,material_auxvar_dn%porosity, &
         Dk_eff_dn,dk_ds_dn,dk_dT_dn,option)

  endif

  Dk = (Dk_eff_up * Dk_eff_dn) / (dd_dn*Dk_eff_up + dd_up*Dk_eff_dn)
  cond = Dk*area*(global_auxvar_up%temp - global_auxvar_dn%temp)

  fluxe = fluxe + cond

  ! If only solving the energy equation, ensure Res(1) is zero
  if (option%flow%only_energy_eq) fluxm = 0.d0

  Res(1:option%nflowdof-1) = fluxm
  Res(option%nflowdof) = fluxe

 ! note: Res is the flux contribution, for node 1 R = R + Res_FL
 !                                              2 R = R - Res_FL


end subroutine THFlux

! ************************************************************************** !

subroutine THBCFluxDerivative(ibndtype,auxvars, &
                              auxvar_up,global_auxvar_up, &
                              auxvar_dn,global_auxvar_dn, &
                              material_auxvar_dn, &
                              Dk_dn, &
                              area, &
                              dist, &
                              sir_dn, &
                              option, &
                              sf_dn, &
                              cc_dn, &
                              tcc_dn, &
                              Dk_dry_dn, &
                              Dk_ice_dn, &
                              Jdn)
  !
  ! Computes the derivatives of the boundary flux
  ! terms for the Jacobian
  !
  ! Author: ???
  ! Date: 12/13/07
  !
  use Option_module
  use Saturation_Function_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Thermal_module
  use Connection_module
  use EOS_Water_module
  use Utility_module

  implicit none

  PetscInt :: ibndtype(:)
  type(TH_auxvar_type) :: auxvar_up, auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_dn
  type(option_type) :: option

  PetscReal :: sir_dn
  PetscBool :: is_flowing
  PetscReal :: auxvars(:) ! from aux_real_var array in boundary condition
  PetscReal :: por_dn,perm_dn,Dk_dn,tor_dn
  PetscReal :: area
  type(saturation_function_type) :: sf_dn
  class(characteristic_curves_type) :: cc_dn
  class(cc_thermal_type), pointer :: tcc_dn

  PetscReal :: Dk_dry_dn
  PetscReal :: Dk_ice_dn
  PetscReal :: Jdn(option%nflowdof,option%nflowdof)
  PetscReal :: dist(-1:3)

  PetscReal :: dist_gravity  ! distance along gravity vector

  PetscReal :: dd_dn
  PetscReal :: v_darcy
  PetscReal :: fluxm,fluxe,q,density_ave
  PetscReal :: uh,ukvr,diffdp,DK,Dq
  PetscReal :: upweight,gravity,dphi

  PetscReal :: ddiff_dp_dn, ddiff_dT_dn
  PetscReal :: dden_ave_dp_dn, dden_ave_dT_dn
  PetscReal :: dgravity_dden_dn
  PetscReal :: dphi_dp_dn, dphi_dT_dn
  PetscReal :: dukvr_dp_dn, dukvr_dT_dn
  PetscReal :: duh_dp_dn, duh_dT_dn
  PetscReal :: dq_dp_dn, dq_dT_dn
  PetscReal :: Dk_eff_dn
  PetscReal :: dk_ds_dn, dK_di_dn, dk_dT_dn
  PetscReal :: dDk_dT_dn, dDk_dp_dn
  PetscReal :: dKe_dT_dn, dKe_dp_dn
  PetscReal :: dKe_fr_dT_dn, dKe_fr_dp_dn

#if 0
  PetscInt :: iphase, ideriv
  type(TH_auxvar_type) :: auxvar_pert_dn, auxvar_pert_up
  type(global_auxvar_type) :: global_auxvar_pert_dn, global_auxvar_pert_up
  type(material_auxvar_type), allocatable :: material_auxvar_pert_dn, &
                                              material_auxvar_pert_up

  PetscReal :: perturbation
  PetscReal :: x_up(option%nflowdof), x_dn(option%nflowdof)
  PetscReal :: x_pert_up(option%nflowdof), x_pert_dn(option%nflowdof)
  PetscReal :: pert_up, pert_dn
  PetscReal :: res(option%nflowdof)
  PetscReal :: res_pert_up(option%nflowdof)
  PetscReal :: res_pert_dn(option%nflowdof)
  PetscReal :: J_pert_dn(option%nflowdof,option%nflowdof)
#endif

  ! ice variables
  PetscReal :: Ddiffgas_avg, Ddiffgas_up, Ddiffgas_dn
  PetscReal :: p_g
  PetscReal :: deng_up, deng_dn
  PetscReal :: psat_up, psat_dn
  PetscReal :: molg_up, molg_dn
  PetscReal :: satg_up, satg_dn
  PetscReal :: Diffg_up, Diffg_dn
  PetscReal :: ddeng_dT_dn
  PetscReal :: dpsat_dT_dn
  PetscReal :: dmolg_dT_dn
  PetscReal :: dDiffg_dT_dn
  PetscReal :: dDiffg_dp_dn
  PetscReal :: dsatg_dp_dn
  PetscReal :: Diffg_ref, p_ref, T_ref
  PetscErrorCode :: ierr
  PetscReal :: T_th
  PetscBool :: skip_thermal_conduction

  skip_thermal_conduction = PETSC_FALSE
  T_th  = 0.5d0

  fluxm = 0.d0
  fluxe = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0

  Jdn = 0.d0

  dden_ave_dp_dn = 0.d0
  dden_ave_dT_dn = 0.d0
  ddiff_dp_dn = 0.d0
  ddiff_dT_dn = 0.d0
  dgravity_dden_dn = 0.d0
  dphi_dp_dn = 0.d0
  dphi_dT_dn = 0.d0
  dukvr_dp_dn = 0.d0
  dukvr_dT_dn = 0.d0
  duh_dp_dn = 0.d0
  duh_dT_dn = 0.d0
  dq_dp_dn = 0.d0
  dq_dT_dn = 0.d0

  dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))
  dd_dn = dist(0)

  call PermeabilityTensorToScalar(material_auxvar_dn,dist,perm_dn)
  por_dn = material_auxvar_dn%porosity_base
  tor_dn = material_auxvar_dn%tortuosity

  ! Flow
  diffdp = por_dn*tor_dn/dd_dn*area
  select case(ibndtype(TH_PRESSURE_DOF))
    ! figure out the direction of flow
    case(DIRICHLET_BC,DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
         HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC, &
         HET_DIRICHLET_BC,HET_HYDROSTATIC_SEEPAGE_BC, &
         HET_HYDROSTATIC_CONDUCTANCE_BC)
      if (ibndtype(TH_PRESSURE_DOF) == DIRICHLET_CONDUCTANCE_BC .or. &
          ibndtype(TH_PRESSURE_DOF) == HYDROSTATIC_CONDUCTANCE_BC .or. &
          ibndtype(TH_PRESSURE_DOF) == HET_HYDROSTATIC_CONDUCTANCE_BC) then
        Dq = auxvars(TH_CONDUCTANCE_DOF)
      else
        Dq = perm_dn / dd_dn
      endif
      ! Flow term
      is_flowing = PETSC_FALSE

      if (option%flow%th_freezing) then
        if (global_auxvar_up%sat(1) > sir_dn .or.  &
            global_auxvar_dn%sat(1) > sir_dn) then
          is_flowing = PETSC_TRUE
        endif
      else
        if (auxvar_up%kvr > eps .or. &
            auxvar_dn%kvr > eps) then
          is_flowing = PETSC_TRUE
        endif
      endif

      if (is_flowing) then
        upweight=1.D0
        if (global_auxvar_up%sat(1) < eps) then
          upweight=0.d0
        else if (global_auxvar_dn%sat(1) < eps) then
          upweight=1.d0
        endif

        density_ave = upweight*global_auxvar_up%den(1)+ &
                      (1.D0-upweight)*global_auxvar_dn%den(1)
        dden_ave_dp_dn = (1.D0-upweight)*auxvar_dn%dden_dp
        dden_ave_dT_dn = (1.D0-upweight)*auxvar_dn%dden_dT

        if (ibndtype(TH_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
          dden_ave_dT_dn = dden_ave_dT_dn + upweight*auxvar_up%dden_dT
        endif

        gravity = (upweight*global_auxvar_up%den(1)*auxvar_up%avgmw + &
                  (1.D0-upweight)*global_auxvar_dn%den(1)*auxvar_dn%avgmw) &
                  * dist_gravity
        dgravity_dden_dn = (1.d0-upweight)*auxvar_dn%avgmw*dist_gravity

        if (th_ice_model /= DALL_AMICO) then
          dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity
          dphi_dp_dn = -1.d0 + dgravity_dden_dn*auxvar_dn%dden_dp
          dphi_dT_dn = dgravity_dden_dn*auxvar_dn%dden_dT
        else
          dphi = auxvar_up%ice%pres_fh2o - auxvar_dn%ice%pres_fh2o + gravity
          dphi_dp_dn = -auxvar_dn%ice%dpres_fh2o_dp + &
                        dgravity_dden_dn*auxvar_dn%dden_dp
          dphi_dT_dn = -auxvar_dn%ice%dpres_fh2o_dT + &
                        dgravity_dden_dn*auxvar_dn%dden_dT
        endif

        select case(ibndtype(TH_PRESSURE_DOF))
          case(HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC, &
               DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
               HET_HYDROSTATIC_SEEPAGE_BC,HET_HYDROSTATIC_CONDUCTANCE_BC)
            ! boundary cell is <= pref
            if (global_auxvar_up%pres(1)- &
                option%flow%reference_pressure < eps) then
              ! skip thermal conduction whenever water table is lower than cell
              skip_thermal_conduction = PETSC_TRUE
              ! flow inward
              if (dphi > 0.d0) then
                dphi = 0.d0
                dphi_dp_dn = 0.d0
                dphi_dT_dn = 0.d0
              endif
            endif
        end select

        if (ibndtype(TH_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
                     !( dgravity_dden_up                   ) (dden_dt_up)
          dphi_dT_dn = dphi_dT_dn + upweight*auxvar_up%avgmw* &
                                    dist_gravity*auxvar_up%dden_dT
        endif


        if (dphi>=0.D0) then
          ukvr = auxvar_up%kvr
          if (ibndtype(TH_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
            dukvr_dT_dn = auxvar_up%dkvr_dT
          endif
        else
          ukvr = auxvar_dn%kvr
          dukvr_dp_dn = auxvar_dn%dkvr_dp
          dukvr_dT_dn = auxvar_dn%dkvr_dT
        endif

        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
          q = v_darcy * area
          dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
          dq_dT_dn = Dq*(dukvr_dT_dn*dphi+ukvr*dphi_dT_dn)*area
        endif
      endif

    case(HET_SURF_HYDROSTATIC_SEEPAGE_BC)
      Dq = perm_dn / dd_dn

      ! Flow term
      is_flowing = PETSC_FALSE

      if (option%flow%th_freezing) then
        if (global_auxvar_up%sat(1) > sir_dn .or.  &
            global_auxvar_dn%sat(1) > sir_dn) then
          is_flowing = PETSC_TRUE
        endif
      else
        if (auxvar_up%kvr > eps .or. &
            auxvar_dn%kvr > eps) then
          is_flowing = PETSC_TRUE
        endif
      endif

      if (is_flowing) then
        upweight=1.D0
        if (global_auxvar_up%sat(1) < eps) then
          upweight=0.d0
        else if (global_auxvar_dn%sat(1) < eps) then
          upweight=1.d0
        endif

        density_ave = upweight*global_auxvar_up%den(1)+ &
                      (1.D0-upweight)*global_auxvar_dn%den(1)
        dden_ave_dp_dn = (1.D0-upweight)*auxvar_dn%dden_dp
        dden_ave_dT_dn = (1.D0-upweight)*auxvar_dn%dden_dT

        if (ibndtype(TH_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
          dden_ave_dT_dn = dden_ave_dT_dn + upweight*auxvar_up%dden_dT
        endif

        gravity = (upweight*global_auxvar_up%den(1)*auxvar_up%avgmw + &
                  (1.D0-upweight)*global_auxvar_dn%den(1)*auxvar_dn%avgmw) &
                  * dist_gravity
        dgravity_dden_dn = (1.d0-upweight)*auxvar_dn%avgmw*dist_gravity

        if (th_ice_model /= DALL_AMICO) then
          dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity
          dphi_dp_dn = -1.d0 + dgravity_dden_dn*auxvar_dn%dden_dp
          dphi_dT_dn = dgravity_dden_dn*auxvar_dn%dden_dT
        else
          dphi = auxvar_up%ice%pres_fh2o - auxvar_dn%ice%pres_fh2o + gravity
          dphi_dp_dn = -auxvar_dn%ice%dpres_fh2o_dp + &
                       dgravity_dden_dn*auxvar_dn%dden_dp
          dphi_dT_dn = -auxvar_dn%ice%dpres_fh2o_dT + &
                       dgravity_dden_dn*auxvar_dn%dden_dT
        endif

        ! flow in         ! boundary cell is <= pref
        if (dphi > 0.d0 .and. &
            global_auxvar_up%pres(1)-option%flow%reference_pressure < eps) then
          dphi = 0.d0
          dphi_dp_dn = 0.d0
          dphi_dT_dn = 0.d0
        endif

        if (ibndtype(TH_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
                        !( dgravity_dden_up                   ) (dden_dT_up)
          dphi_dT_dn = dphi_dT_dn + upweight*auxvar_up%avgmw* &
                                    dist_gravity*auxvar_up%dden_dT
        endif

        if (dphi>=0.D0) then
          ukvr = auxvar_up%kvr
          if (ibndtype(TH_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
            dukvr_dT_dn = auxvar_up%dkvr_dT
          endif
        else
          ukvr = auxvar_dn%kvr
          dukvr_dp_dn = auxvar_dn%dkvr_dp
          dukvr_dT_dn = auxvar_dn%dkvr_dT
        endif

        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
          q = v_darcy * area
          dq_dp_dn = Dq*(dukvr_dp_dn*dphi+ukvr*dphi_dp_dn)*area
          dq_dT_dn = Dq*(dukvr_dT_dn*dphi+ukvr*dphi_dT_dn)*area
        endif
      endif

    case(NEUMANN_BC)
      if (dabs(auxvars(TH_PRESSURE_DOF)) > floweps) then
        v_darcy = auxvars(TH_PRESSURE_DOF)
        if (v_darcy > 0.d0) then
          density_ave = global_auxvar_up%den(1)
          if (ibndtype(TH_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
            dden_ave_dT_dn = auxvar_up%dden_dT
          endif
        else
          density_ave = global_auxvar_dn%den(1)
          dden_ave_dp_dn = auxvar_dn%dden_dp
          dden_ave_dT_dn = auxvar_dn%dden_dT
        endif
        q = v_darcy * area
      endif

    case(ZERO_GRADIENT_BC)
      ! do nothing

  end select

  if (v_darcy >= 0.D0) then
    uh = auxvar_up%h
    if (ibndtype(TH_PRESSURE_DOF) == ZERO_GRADIENT_BC) then
      duh_dp_dn = auxvar_up%dh_dp
    endif
    if (ibndtype(TH_TEMPERATURE_DOF) == ZERO_GRADIENT_BC) then
      duh_dT_dn = auxvar_up%dh_dT
    endif
  else
    uh = auxvar_dn%h
    duh_dp_dn = auxvar_dn%dh_dp
    duh_dT_dn = auxvar_dn%dh_dT
  endif

  Jdn(1,1) = (dq_dp_dn*density_ave+q*dden_ave_dp_dn)
  Jdn(1,2) = (dq_dT_dn*density_ave+q*dden_ave_dT_dn)

  ! If only solving the energy equation, ensure Jdn(2,2) has no
  ! contribution from mass equation
  if (option%flow%only_energy_eq) then
    q = 0.d0
    dq_dT_dn = 0.d0
  endif

  ! based on flux = q*density_ave*uh
  Jdn(option%nflowdof,1) =  &
     ((dq_dp_dn*density_ave+q*dden_ave_dp_dn)*uh+q*density_ave*duh_dp_dn)
  Jdn(option%nflowdof,2) =  &
     ((dq_dT_dn*density_ave+q*dden_ave_dT_dn)*uh+q*density_ave*duh_dT_dn)

  ! Conduction term
  select case(ibndtype(TH_TEMPERATURE_DOF))
    case(DIRICHLET_BC,HET_DIRICHLET_BC)
      ! Dk =  auxvar_dn%Dk_eff / dd_dn
      !cond = Dk*area*(global_auxvar_up%temp-global_auxvar_dn%temp)

      if (skip_thermal_conduction) then
        ! skip thermal conducton when the boundary pressure is below the
        ! reference pressure (e.g. river stage is below cell center).
        dDk_dT_dn = 0.d0
        dDk_dp_dn = 0.d0
        Dk = 0.d0
      else
        if (option%flow%th_freezing) then

          call tcc_dn%thermal_conductivity_function% &
               TCondTensorToScalar(dist,option)
          call tcc_dn%thermal_conductivity_function%CalculateFTCond( &
               global_auxvar_dn%sat(1),auxvar_dn%ice%sat_ice, &
               global_auxvar_dn%temp,material_auxvar_dn%porosity, &
               Dk_eff_dn,dk_ds_dn,dK_di_dn,dk_dT_dn,option)

          dKe_dp_dn    = auxvar_dn%dKe_dp
          dKe_dT_dn    = auxvar_dn%dKe_dT
          dKe_fr_dT_dn = auxvar_dn%ice%dKe_fr_dT
          dKe_fr_dp_dn = auxvar_dn%ice%dKe_fr_dp
          Dk           = Dk_eff_dn/dd_dn

          dDk_dT_dn = Dk**2/Dk_eff_dn**2*dd_dn*(Dk_dn*dKe_dT_dn + &
              Dk_ice_dn*dKe_fr_dT_dn + (- dKe_dT_dn - dKe_fr_dT_dn)* &
              Dk_dry_dn)
          dDk_dp_dn = Dk**2/Dk_eff_dn**2*dd_dn*(Dk_dn*dKe_dp_dn + &
              Dk_ice_dn*dKe_fr_dp_dn + (- dKe_dp_dn - dKe_fr_dp_dn)* &
              Dk_dry_dn)

        else

          ! Dk_eff_dn = auxvar_dn%Dk_eff

          call tcc_dn%thermal_conductivity_function% &
               TCondTensorToScalar(dist,option)
          call tcc_dn%thermal_conductivity_function%CalculateTCond( &
               global_auxvar_dn%sat(1),global_auxvar_dn%temp, &
               material_auxvar_dn%porosity,Dk_eff_dn,dk_ds_dn,dk_dT_dn,option)

          dKe_dp_dn = auxvar_dn%dKe_dp
          dKe_dT_dn = auxvar_dn%dKe_dT
          Dk        = Dk_eff_dn/dd_dn

          dDk_dT_dn = Dk**2/Dk_eff_dn**2*dd_dn*(Dk_dn - Dk_dry_dn)*dKe_dT_dn
          dDk_dp_dn = Dk**2/Dk_eff_dn**2*dd_dn*(Dk_dn - Dk_dry_dn)*dKe_dp_dn

        endif
      endif

      Jdn(option%nflowdof,1) = Jdn(option%nflowdof,1) + &
              area*(global_auxvar_up%temp - global_auxvar_dn%temp)*dDk_dp_dn

      Jdn(option%nflowdof,2) = Jdn(option%nflowdof,2) + Dk*area*(-1.d0) + &
              area*(global_auxvar_up%temp - global_auxvar_dn%temp)*dDk_dT_dn

      if (option%flow%th_freezing) then
         ! Added by Satish Karra, 11/21/11
         satg_up = auxvar_up%ice%sat_gas
         satg_dn = auxvar_dn%ice%sat_gas
         if ((satg_up > eps) .and. (satg_dn > eps)) then
            p_g = option%flow%reference_pressure  ! set to reference pressure
            deng_up = p_g/(IDEAL_GAS_CONSTANT*(global_auxvar_up%temp + &
                 273.15d0))*1.d-3
            deng_dn = p_g/(IDEAL_GAS_CONSTANT*(global_auxvar_dn%temp + &
                 273.15d0))*1.d-3

            ! Reference diffusivity, need to read from input file
            Diffg_ref = 2.13D-5
            p_ref = 1.01325d5 ! in Pa
            T_ref = 25.d0 ! in deg C

            Diffg_up = Diffg_ref*(p_ref/p_g)*((global_auxvar_up%temp + &
                 273.15d0)/(T_ref + 273.15d0))**(1.8)
            Diffg_dn = Diffg_ref*(p_ref/p_g)*((global_auxvar_dn%temp + &
                 273.15d0)/(T_ref + 273.15d0))**(1.8)
            Ddiffgas_up = satg_up*deng_up*Diffg_up
            Ddiffgas_dn = satg_dn*deng_dn*Diffg_dn
            call EOSWaterSaturationPressure(global_auxvar_up%temp, &
                                            psat_up, ierr)
            call EOSWaterSaturationPressure(global_auxvar_dn%temp, &
                                            psat_dn, dpsat_dT_dn, ierr)
            molg_up = psat_up/p_g
            molg_dn = psat_dn/p_g
            ddeng_dT_dn = - p_g/(IDEAL_GAS_CONSTANT*(global_auxvar_dn%temp + &
                 273.15d0)**2)*1.d-3
            dmolg_dT_dn = (1/p_g)*dpsat_dT_dn
            dDiffg_dT_dn = 1.8*Diffg_dn/(global_auxvar_dn%temp + 273.15d0)
            dDiffg_dp_dn = 0.d0
            dsatg_dp_dn = auxvar_dn%ice%dsat_gas_dp

            if (molg_up > molg_dn) then
               upweight = 0.d0
            else
               upweight = 1.d0
            endif

            Ddiffgas_avg = upweight*Ddiffgas_up+(1.D0 - upweight)*Ddiffgas_dn

            Jdn(1,1) = Jdn(1,1) + por_dn*tor_dn*(1.D0 - upweight)* &
                 Ddiffgas_dn/satg_dn*dsatg_dp_dn*(molg_up - molg_dn)/dd_dn* &
                 area
            Jdn(1,2) = Jdn(1,2) + por_dn*tor_dn*(1.D0 - upweight)* &
                 (Ddiffgas_avg/deng_dn*ddeng_dT_dn + Ddiffgas_avg/Diffg_dn* &
                 dDiffg_dT_dn)*(molg_up - molg_dn)/dd_dn*area + por_dn* &
                 tor_dn*Ddiffgas_avg*(-dmolg_dT_dn)/dd_dn*area
         endif
      endif

  end select

  ! If only solving the energy equation,
  !  - Set jacobian term corresponding to mass-equation to zero, and
  !  - Set off-diagonal jacobian terms to zero.
  if (option%flow%only_energy_eq) then
    Jdn(1,1) = 0.d0
    Jdn(1,2) = 0.d0
    Jdn(option%nflowdof,1) = 0.d0
  endif

#if 0
  if (option%flow%numerical_derivatives) then
    allocate(material_auxvar_pert_up,material_auxvar_pert_dn)

    call MaterialAuxVarInit(material_auxvar_pert_up,option)
    call MaterialAuxVarInit(material_auxvar_pert_dn,option)

    call GlobalAuxVarInit(global_auxvar_pert_up,option)
    call GlobalAuxVarInit(global_auxvar_pert_dn,option)
    call THAuxVarCopy(auxvar_up,auxvar_pert_up,option)
    call THAuxVarCopy(auxvar_dn,auxvar_pert_dn,option)
    call GlobalAuxVarCopy(global_auxvar_up,global_auxvar_pert_up,option)
    call GlobalAuxVarCopy(global_auxvar_dn,global_auxvar_pert_dn,option)

    call MaterialAuxVarCopy(material_auxvar_dn,material_auxvar_pert_up, &
                            option)
    call MaterialAuxVarCopy(material_auxvar_dn,material_auxvar_pert_dn, &
                            option)

    x_up(1) = global_auxvar_up%pres(1)
    x_up(2) = global_auxvar_up%temp
    x_dn(1) = global_auxvar_dn%pres(1)
    x_dn(2) = global_auxvar_dn%temp
    do ideriv = 1,3
      if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
        x_up(ideriv) = x_dn(ideriv)
      endif
    enddo
    if (option%flow%th_freezing) then
       call THAuxVarComputeFreezing(x_dn,auxvar_dn, &
            global_auxvar_dn, &
            material_auxvar_dn, &
            iphase,sf_dn, &
            th_parameter,icct_up, &
            -999,PETSC_TRUE,option) ! do not perturb boundary porosity
       call THAuxVarComputeFreezing(x_up,auxvar_up, &
            global_auxvar_up, &
            material_auxvar_up, &
            iphase,sf_dn, &
            th_parameter,icct_up, &
            -999,PETSC_FALSE,option)
    else
       call THAuxVarComputeNoFreezing(x_dn,auxvar_dn, &
            global_auxvar_dn, &
            material_auxvar_dn, &
            iphase,cc_dn, &
            -999,PETSC_TRUE,option)
       call THAuxVarComputeNoFreezing(x_up,auxvar_up, &
            global_auxvar_up, &
            material_auxvar_up, &
            iphase,cc_dn, &
            -999,PETSC_FALSE,option) ! do not perturb boundary porosity
    endif

    call THBCFlux(ibndtype,auxvars,auxvar_up,global_auxvar_up, &
                  material_auxvar_up, &
                  auxvar_dn,global_auxvar_dn, &
                  material_auxvar_dn, &
                  Dk_dn, tcc_dn, &
                  area,dist_gravity,sir_dn,option,v_darcy, &
                  fluxe_bulk, fluxe_cond, &
                  res)
    if (ibndtype(TH_PRESSURE_DOF) == ZERO_GRADIENT_BC .or. &
        ibndtype(TH_TEMPERATURE_DOF) == ZERO_GRADIENT_BC ) then
      x_pert_up = x_up
    endif

    do ideriv = 1,option%nflowdof
      pert_dn = x_dn(ideriv)*perturbation_tolerance
      x_pert_dn = x_dn

      if (option%flow%th_freezing) then

         if (ideriv == 1) then
            if (x_pert_dn(ideriv) < option%flow%reference_pressure) then
               pert_dn = - pert_dn
            endif
            x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
         endif

         if (ideriv == 2) then
            if (x_pert_dn(ideriv) < 0.d0) then
               pert_dn = - 1.d-5
            else
               pert_dn = 1.d-5
            endif
            x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
         endif
      else
         x_pert_dn(ideriv) = x_pert_dn(ideriv) + pert_dn
      endif

      x_pert_up = x_up
      if (ibndtype(ideriv) == ZERO_GRADIENT_BC) then
        x_pert_up(ideriv) = x_pert_dn(ideriv)
      endif

      if (option%flow%th_freezing) then
         call THAuxVarComputeFreezing(x_pert_dn,auxvar_pert_dn, &
              global_auxvar_pert_dn, &
              material_auxvar_pert_dn, &
              iphase,sf_dn, &
              -999,PETSC_TRUE,option)
         call THAuxVarComputeFreezing(x_pert_up,auxvar_pert_up, &
              global_auxvar_pert_up, &
              material_auxvar_pert_up, &
              iphase,sf_dn, &
              -999,PETSC_FALSE,option) ! do not perturb boundary porosity
      else
         call THAuxVarComputeNoFreezing(x_pert_dn,auxvar_pert_dn, &
              global_auxvar_pert_dn, &
              material_auxvar_pert_dn, &
              iphase,cc_dn, &
              -999,PETSC_TRUE,option)
         call THAuxVarComputeNoFreezing(x_pert_up,auxvar_pert_up, &
              global_auxvar_pert_up, &
              material_auxvar_pert_up, &
              iphase,cc_dn, &
              -999,PETSC_FALSE,option) ! do not perturb boundary porosity
      endif

      call THBCFlux(ibndtype,auxvars,auxvar_pert_up,global_auxvar_pert_up, &
                    material_auxvar_pert_up, &
                    auxvar_pert_dn,global_auxvar_pert_dn, &
                    material_auxvar_pert_dn, &
                    Dk_dn, tcc_dn, &
                    area,dist_gravity,sir_dn,option,v_darcy, &
                    fluxe_bulk, fluxe_cond, &
                    res_pert_dn)
      J_pert_dn(:,ideriv) = (res_pert_dn(:)-res(:))/pert_dn
    enddo
    Jdn = J_pert_dn
    call GlobalAuxVarStrip(global_auxvar_pert_up)
    call GlobalAuxVarStrip(global_auxvar_pert_dn)
  endif
#endif

end subroutine THBCFluxDerivative

! ************************************************************************** !

subroutine THBCFlux(ibndtype,auxvars,auxvar_up,global_auxvar_up, &
                    auxvar_dn,global_auxvar_dn, &
                    material_auxvar_dn, &
                    Dk_dn, tcc_dn, &
                    area, &
                    dist,sir_dn, &
                    option,v_darcy, &
                    fluxe_bulk, fluxe_cond, &
                    Res)
  !
  ! Computes the  boundary flux terms for the residual
  !
  ! Author: ???
  ! Date: 12/13/07
  !
  use Option_module
  use Connection_module
  use EOS_Water_module
  use Condition_module
  use Utility_module
  use Characteristic_Curves_Thermal_module

  implicit none

  PetscInt :: ibndtype(:)
  type(TH_auxvar_type) :: auxvar_up, auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_dn
  class(cc_thermal_type), pointer :: tcc_dn
  type(option_type) :: option
  PetscReal :: sir_dn
  PetscBool :: is_flowing

  PetscReal :: auxvars(:) ! from aux_real_var array
  PetscReal :: Dk_dn
  PetscReal :: v_darcy, area
  PetscReal :: Res(1:option%nflowdof)
  PetscReal :: dist(-1:3)
  PetscReal, intent(out) :: fluxe_bulk, fluxe_cond

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dd_dn

  PetscReal :: por_dn,perm_dn,tor_dn
  PetscReal :: fluxm,fluxe,q,density_ave
  PetscReal :: uh,ukvr,diffdp,DK,Dq
  PetscReal :: upweight,cond,gravity,dphi
  PetscReal :: dk_ds_dn, dK_di_dn, dk_dT_dn

  ! ice variables
  PetscReal :: Ddiffgas_avg, Ddiffgas_dn, Ddiffgas_up
  PetscReal :: p_g
  PetscReal :: deng_dn, deng_up
  PetscReal :: psat_dn, psat_up
  PetscReal :: molg_dn, molg_up
  PetscReal :: satg_dn, satg_up
  PetscReal :: Diffg_dn, Diffg_up
  PetscReal :: Diffg_ref, p_ref, T_ref
  PetscErrorCode :: ierr
  PetscReal :: fv_up, fv_dn
  PetscReal :: T_th,fctT
  PetscBool :: skip_thermal_conduction

  skip_thermal_conduction = PETSC_FALSE
  T_th  = 0.5d0

  fluxm = 0.d0
  fluxe = 0.d0
  v_darcy = 0.d0
  density_ave = 0.d0
  q = 0.d0
  fctT = 0.d0
  fluxe_bulk = 0.d0
  fluxe_cond = 0.d0

  dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))
  dd_dn = dist(0)

  call PermeabilityTensorToScalar(material_auxvar_dn,dist,perm_dn)
  por_dn = material_auxvar_dn%porosity_base
  tor_dn = material_auxvar_dn%tortuosity

  ! Flow
  diffdp = por_dn*tor_dn/dd_dn*area
  select case(ibndtype(TH_PRESSURE_DOF))
    case(DIRICHLET_BC,DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
         HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC, &
         HET_DIRICHLET_BC,HET_HYDROSTATIC_SEEPAGE_BC, &
         HET_HYDROSTATIC_CONDUCTANCE_BC)
      if (ibndtype(TH_PRESSURE_DOF) == DIRICHLET_CONDUCTANCE_BC .or. &
          ibndtype(TH_PRESSURE_DOF) == HYDROSTATIC_CONDUCTANCE_BC .or. &
          ibndtype(TH_PRESSURE_DOF) == HET_HYDROSTATIC_CONDUCTANCE_BC) then
        Dq = auxvars(TH_CONDUCTANCE_DOF)
      else
        Dq = perm_dn / dd_dn
      endif

      ! Flow term
      is_flowing = PETSC_FALSE

      if (option%flow%th_freezing) then
        if (global_auxvar_up%sat(1) > sir_dn .or.  &
            global_auxvar_dn%sat(1) > sir_dn) then
          is_flowing = PETSC_TRUE
        endif
      else
        if (auxvar_up%kvr > eps .or. &
            auxvar_dn%kvr > eps) then
          is_flowing = PETSC_TRUE
        endif
      endif

      if (is_flowing) then
        upweight=1.D0
        if (global_auxvar_up%sat(1) < eps) then
          upweight=0.d0
        else if (global_auxvar_dn%sat(1) < eps) then
          upweight=1.d0
        endif
        density_ave = upweight*global_auxvar_up%den(1)+ &
                      (1.D0-upweight)*global_auxvar_dn%den(1)

        gravity = (upweight*global_auxvar_up%den(1)*auxvar_up%avgmw + &
                  (1.D0-upweight)*global_auxvar_dn%den(1)*auxvar_dn%avgmw) &
                  * dist_gravity

        if (th_ice_model /= DALL_AMICO) then
          dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity
        else
          dphi = auxvar_up%ice%pres_fh2o - auxvar_dn%ice%pres_fh2o + gravity
        endif

        select case(ibndtype(TH_PRESSURE_DOF))
          case(HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC, &
               DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
               HET_HYDROSTATIC_SEEPAGE_BC,HET_HYDROSTATIC_CONDUCTANCE_BC)
            ! boundary cell is <= pref
            if (global_auxvar_up%pres(1)- &
                option%flow%reference_pressure < eps) then
              ! skip thermal conduction whenever water table is lower than cell
              skip_thermal_conduction = PETSC_TRUE
              ! flow inward
              if (dphi > 0.d0) then
                dphi = 0.d0
              endif
            endif
        end select

        if (dphi>=0.D0) then
          ukvr = auxvar_up%kvr
        else
          ukvr = auxvar_dn%kvr
        endif

        call InterfaceApprox(auxvar_up%kvr, auxvar_dn%kvr, &
                             dphi, &
                             option%flow%rel_perm_aveg, &
                             ukvr)

        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
        endif
      endif

    case(HET_SURF_HYDROSTATIC_SEEPAGE_BC)
      Dq = perm_dn / dd_dn
      ! Flow term
      is_flowing = PETSC_FALSE

      if (option%flow%th_freezing) then
        if (global_auxvar_up%sat(1) > sir_dn .or.  &
            global_auxvar_dn%sat(1) > sir_dn) then
          is_flowing = PETSC_TRUE
        endif
      else
        if (auxvar_up%kvr > eps .or. &
            auxvar_dn%kvr > eps) then
          is_flowing = PETSC_TRUE
        endif
      endif

      if (is_flowing) then
        upweight=1.D0
        if (global_auxvar_up%sat(1) < eps) then
          upweight=0.d0
        else if (global_auxvar_dn%sat(1) < eps) then
          upweight=1.d0
        endif
        density_ave = upweight*global_auxvar_up%den(1)+ &
                      (1.D0-upweight)*global_auxvar_dn%den(1)

        gravity = (upweight*global_auxvar_up%den(1)*auxvar_up%avgmw + &
             (1.D0-upweight)*global_auxvar_dn%den(1)*auxvar_dn%avgmw) &
             * dist_gravity

        if (th_ice_model /= DALL_AMICO) then
          dphi = global_auxvar_up%pres(1) - global_auxvar_dn%pres(1) + gravity
        else
          dphi = auxvar_up%ice%pres_fh2o - auxvar_dn%ice%pres_fh2o + gravity
        endif

        if (dphi > 0.d0 .and. &
            global_auxvar_up%pres(1)-option%flow%reference_pressure < eps) then
          dphi = 0.d0
        endif

        if (dphi>=0.D0) then
          ukvr = auxvar_up%kvr
        else
          ukvr = auxvar_dn%kvr
        endif

        call InterfaceApprox(auxvar_up%kvr, auxvar_dn%kvr, &
             dphi, &
             option%flow%rel_perm_aveg, &
             ukvr)

        if (ukvr*Dq>floweps) then
          v_darcy = Dq * ukvr * dphi
        endif
      endif

    case(NEUMANN_BC)
      if (dabs(auxvars(TH_PRESSURE_DOF)) > floweps) then
        v_darcy = auxvars(TH_PRESSURE_DOF)
        if (v_darcy > 0.d0) then
          density_ave = global_auxvar_up%den(1)
        else
          density_ave = global_auxvar_dn%den(1)
        endif
      endif

    case(ZERO_GRADIENT_BC)
      ! do nothing needed to bypass default case

    case default
      option%io_buffer = 'BC type "' // &
        trim(GetSubConditionType(ibndtype(TH_PRESSURE_DOF))) // &
        '" not implemented in TH mode.'
      call PrintErrMsg(option)

  end select

  ! If only solving the energy equation, ensure Res(2) has no
  ! contribution from mass equation by setting darcy velocity
  ! to be zero
  if (option%flow%only_energy_eq) q = 0.d0

  q = v_darcy * area

  if (v_darcy >= 0.D0) then
    uh = auxvar_up%h
  else
    uh = auxvar_dn%h
  endif

  fluxm = fluxm + q*density_ave
  fluxe = fluxe + q*density_ave*uh
  fluxe_bulk = q*density_ave*uh

  ! Conduction term
  select case(ibndtype(TH_TEMPERATURE_DOF))
    case(DIRICHLET_BC,HET_DIRICHLET_BC)
      call tcc_dn%thermal_conductivity_function% &
      TCondTensorToScalar(dist,option)
      if (option%flow%th_freezing) then
        call tcc_dn%thermal_conductivity_function%CalculateFTCond( &
        global_auxvar_dn%sat(1),auxvar_dn%ice%sat_ice, &
        global_auxvar_dn%temp,material_auxvar_dn%porosity,auxvar_dn%Dk_eff, &
        dk_ds_dn,dK_di_dn,dk_dT_dn,option)
      else
        call tcc_dn%thermal_conductivity_function%CalculateTCond( &
         global_auxvar_dn%sat(1),global_auxvar_dn%temp,material_auxvar_dn%porosity, &
         auxvar_dn%Dk_eff,dk_ds_dn,dk_dT_dn,option)
      endif

      Dk =  auxvar_dn%Dk_eff / dd_dn

      cond = Dk*area*(global_auxvar_up%temp-global_auxvar_dn%temp)

      if (skip_thermal_conduction) then
        ! skip thermal conducton when the boundary pressure is below the
        ! reference pressure (e.g. river stage is below cell center).
        cond = 0.d0
      endif

      fluxe = fluxe + cond
      fluxe_cond = cond

      if (option%flow%th_freezing) then
         ! Added by Satish Karra,
         satg_up = auxvar_up%ice%sat_gas
         satg_dn = auxvar_dn%ice%sat_gas
         if ((satg_up > eps) .and. (satg_dn > eps)) then
            p_g = option%flow%reference_pressure ! set to reference pressure
            deng_up = p_g/(IDEAL_GAS_CONSTANT* &
                      (global_auxvar_up%temp + 273.15d0))*1.d-3
            deng_dn = p_g/(IDEAL_GAS_CONSTANT* &
                      (global_auxvar_dn%temp + 273.15d0))*1.d-3

            ! Reference diffusivity, need to read from input file
            Diffg_ref = 2.13D-5
            p_ref = 1.01325d5 ! in Pa
            T_ref = 25.d0 ! in deg C

            Diffg_up = Diffg_ref*(p_ref/p_g)*((global_auxvar_up%temp + &
                 273.15d0)/(T_ref + 273.15d0))**(1.8)
            Diffg_dn = Diffg_ref*(p_ref/p_g)*((global_auxvar_dn%temp + &
                 273.15d0)/(T_ref + 273.15d0))**(1.8)
            Ddiffgas_up = satg_up*deng_up*Diffg_up
            Ddiffgas_dn = satg_dn*deng_dn*Diffg_dn
            call EOSWaterSaturationPressure(global_auxvar_up%temp,psat_up,ierr)
            call EOSWaterSaturationPressure(global_auxvar_dn%temp,psat_dn,ierr)

            ! vapor pressure lowering due to capillary pressure
            fv_up = exp(-auxvar_up%pc/(global_auxvar_up%den(1)* &
                 IDEAL_GAS_CONSTANT*(global_auxvar_up%temp + 273.15d0)))
            fv_dn = exp(-auxvar_dn%pc/(global_auxvar_dn%den(1)* &
                 IDEAL_GAS_CONSTANT*(global_auxvar_dn%temp + 273.15d0)))

            molg_up = psat_up*fv_up/p_g
            molg_dn = psat_dn*fv_dn/p_g

            if (molg_up > molg_dn) then
               upweight = 0.d0
            else
               upweight = 1.d0
            endif

            Ddiffgas_avg = upweight*Ddiffgas_up + (1.D0 - upweight)*Ddiffgas_dn
            fluxm = fluxm + por_dn*tor_dn*Ddiffgas_avg*(molg_up - molg_dn)/ &
                 dd_dn*area
         endif
      endif

    case(NEUMANN_BC)
      !geh: default internal energy units are MJ
      !       (option%scale = 1.d-6 is for J->MJ)
      ! added by SK 10/18/11
      fluxe = fluxe + auxvars(TH_TEMPERATURE_DOF)*area*(1.d6*option%scale)
      fluxe_cond = auxvars(TH_TEMPERATURE_DOF)*area*(1.d6*option%scale)
    case(ZERO_GRADIENT_BC)
      ! No change in fluxe
    case default
      option%io_buffer = 'BC type "' // &
        trim(GetSubConditionType(ibndtype(TH_TEMPERATURE_DOF))) // &
        '" not implemented in TH mode.'
      call PrintErrMsg(option)
  end select

  ! If only solving the energy equation, set Res(1) is 0.d0
  if (option%flow%only_energy_eq) fluxm = 0.d0

  Res(1:option%nflowspec) = fluxm
  Res(option%nflowdof) = fluxe

end subroutine THBCFlux

! ************************************************************************** !

subroutine THResidual(snes,xx,r,realization,ierr)
  !
  ! Computes the residual equation
  !
  ! Author: ???
  ! Date: 12/10/07
  !

  use Realization_Subsurface_class
  use Patch_module
  use Discretization_module
  use Field_module
  use Option_module
  use Variables_module
  use Material_module
  use Debug_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  class(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr

  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string
  PetscViewer :: viewer

  field => realization%field
  discretization => realization%discretization
  option => realization%option

 ! check initial guess -----------------------------------------------
  ierr = THInitGuessCheck(xx,option)
  if (ierr<0) then
    call SNESSetFunctionDomainError(snes,ierr);CHKERRQ(ierr)
    return
  endif

  call THResidualPreliminaries(xx,r,realization,ierr)

  call THResidualInternalConn(r,realization,ierr)
  call THResidualBoundaryConn(r,realization,ierr)
  call THResidualAccumulation(r,realization,ierr)
  call THResidualSourceSink(r,realization,ierr)

  if (realization%debug%vecview_residual) then
    call DebugWriteFilename(realization%debug,string,'THresidual','', &
                            TH_ts_count,TH_ts_cut_count, &
                            TH_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(r,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif
  if (realization%debug%vecview_solution) then
    call DebugWriteFilename(realization%debug,string,'THxx','', &
                            TH_ts_count,TH_ts_cut_count, &
                            TH_ni_count)
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call VecView(xx,viewer,ierr);CHKERRQ(ierr)
    call DebugViewerDestroy(realization%debug,viewer)
  endif

end subroutine THResidual

! ************************************************************************** !

subroutine THResidualPreliminaries(xx,r,realization,ierr)
  !
  ! Perform preliminary work prior to residual computation
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/06/2019
  !

  use Connection_module
  use Realization_Subsurface_class
  use Patch_module
  use Option_module

  implicit none

  Vec, intent(inout) :: xx
  Vec, intent(inout) :: r
  class(realization_subsurface_type) :: realization

  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  PetscErrorCode :: ierr

  patch => realization%patch
  option => realization%option

  call VecZeroEntries(r,ierr);CHKERRQ(ierr)

  call THUpdateLocalVecs(xx,realization,ierr)

  call THUpdateAuxVarsPatch(realization)
  ! override flags since they will soon be out of date
  patch%aux%TH%auxvars_up_to_date = PETSC_FALSE

  if (option%compute_mass_balance_new) then
    call THZeroMassBalDeltaPatch(realization)
  endif

end subroutine THResidualPreliminaries

! ************************************************************************** !

subroutine THUpdateLocalVecs(xx,realization,ierr)
  !
  ! Updates local vectors needed for residual computation
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/06/2019
  !

  use Realization_Subsurface_class
  use Field_module
  use Discretization_module
  use Option_module
  use Logging_module
  use Material_module
  use Material_Aux_module
  use Variables_module
  use Debug_module

  implicit none

  Vec :: xx
  class(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr

  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(option_type), pointer :: option

  field => realization%field
  discretization => realization%discretization
  option => realization%option

   ! Communication -----------------------------------------
  ! These 3 must be called before THUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)

  call MaterialGetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_X,ZERO_INTEGER)
  call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_X,ZERO_INTEGER)

  call MaterialGetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Y,ZERO_INTEGER)
  call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Y,ZERO_INTEGER)

  call MaterialGetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Z,ZERO_INTEGER)
  call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(realization%patch%aux%Material,field%work_loc, &
                               PERMEABILITY_Z,ZERO_INTEGER)

end subroutine THUpdateLocalVecs

! ************************************************************************** !

subroutine THResidualInternalConn(r,realization,ierr)
  !
  ! Perform preliminary work prior to residual computation
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/06/2019
  !

  use Connection_module
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Secondary_Continuum_Aux_module
  use Secondary_Continuum_module
  use Characteristic_Curves_Thermal_module

  implicit none

  Vec :: r
  class(realization_subsurface_type) :: realization

  PetscErrorCode :: ierr
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn

  PetscReal, pointer :: r_p(:)

  PetscInt :: icc_up, icc_dn, icct_up, icct_dn
  PetscReal :: dd_up, dd_dn
           ! thermal conductivity wet constants at upstream, downstream faces.
  PetscReal :: D_up, D_dn
  PetscReal :: Dk_dry_up, Dk_dry_dn ! dry thermal conductivities
  PetscReal :: Dk_ice_up, Dk_ice_dn ! frozen soil thermal conductivities
  PetscReal :: alpha_up, alpha_dn
  PetscReal :: alpha_fr_up, alpha_fr_dn

  PetscReal :: upweight
  PetscReal :: Res(realization%option%nflowdof)

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(th_parameter_type), pointer :: th_parameter
  type(TH_auxvar_type), pointer :: auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(sec_heat_type), pointer :: TH_sec_heat_vars(:)
  class(cc_thermal_type), pointer :: tcc_dn, tcc_up
  PetscReal :: v_darcy

  PetscInt :: iconn, istart, iend
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  th_parameter => patch%aux%TH%th_parameter
  auxvars => patch%aux%TH%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  TH_sec_heat_vars => patch%aux%SC_heat%sec_heat_vars

! now assign access pointer to local variables
  call VecGetArrayF90(r,r_p,ierr);CHKERRQ(ierr)
  !print *,' Finished scattering non deriv'

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

      if (option%flow%only_vertical_flow) then
        !geh: place second conditional within first to avoid excessive
        !     dot products when .not. option%flow%only_vertical_flow
        if (dot_product(cur_connection_set%dist(1:3,iconn),unit_z) < &
            1.d-10) cycle
      endif

      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
      ! distance = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = distance * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))
      dd_up = distance*fraction_upwind
      dd_dn = distance-dd_up ! should avoid truncation error
      ! upweight could be calculated as 1.d0-fraction_upwind
      ! however, this introduces ever so slight error causing pflow-overhaul not
      ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
      upweight = dd_dn/(dd_up+dd_dn)

      icct_up = patch%cct_id(ghosted_id_up)
      icct_dn = patch%cct_id(ghosted_id_dn)
      icc_up = patch%cc_id(ghosted_id_up)
      icc_dn = patch%cc_id(ghosted_id_dn)

      D_up = th_parameter%ckwet(icct_up)
      D_dn = th_parameter%ckwet(icct_dn)

      Dk_dry_up = th_parameter%ckdry(icct_up)
      Dk_dry_dn = th_parameter%ckdry(icct_dn)

      alpha_up = th_parameter%alpha(icct_up)
      alpha_dn = th_parameter%alpha(icct_dn)

      tcc_up => patch%char_curves_thermal_array(icct_up)%ptr
      tcc_dn => patch%char_curves_thermal_array(icct_dn)%ptr

      if (option%flow%th_freezing) then
         Dk_ice_up = th_parameter%ckfrozen(icct_up)
         DK_ice_dn = th_parameter%ckfrozen(icct_dn)

         alpha_fr_up = th_parameter%alpha_fr(icct_up)
         alpha_fr_dn = th_parameter%alpha_fr(icct_dn)
      else
         Dk_ice_up = Dk_dry_up
         Dk_ice_dn = Dk_dry_dn

         alpha_fr_up = alpha_up
         alpha_fr_dn = alpha_dn
      endif

      call THFlux(auxvars(ghosted_id_up),global_auxvars(ghosted_id_up), &
                  material_auxvars(ghosted_id_up), &
                  D_up, tcc_up, &
                  auxvars(ghosted_id_dn),global_auxvars(ghosted_id_dn), &
                  material_auxvars(ghosted_id_dn), &
                  D_dn, tcc_dn, &
                  cur_connection_set%area(iconn), &
                  cur_connection_set%dist(:,iconn), &
                  upweight,th_parameter%sir(1,icc_up), &
                  th_parameter%sir(1,icc_dn), &
                  option,v_darcy,Dk_dry_up, &
                  Dk_dry_dn,Dk_ice_up,Dk_ice_dn, &
                  alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                  Res)

      patch%internal_velocities(1,sum_connection) = v_darcy
      patch%internal_flow_fluxes(:,sum_connection) = Res(:)

      if (local_id_up>0) then
        iend = local_id_up*option%nflowdof
        istart = iend-option%nflowdof+1
        r_p(istart:iend) = r_p(istart:iend) + Res(1:option%nflowdof)
      endif

      if (local_id_dn>0) then
        iend = local_id_dn*option%nflowdof
        istart = iend-option%nflowdof+1
        r_p(istart:iend) = r_p(istart:iend) - Res(1:option%nflowdof)
      endif

    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  call VecRestoreArrayF90(r,r_p,ierr);CHKERRQ(ierr)


end subroutine THResidualInternalConn

! ************************************************************************** !

subroutine THResidualBoundaryConn(r,realization,ierr)
  !
  ! Perform preliminary work prior to residual computation
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/06/2019
  !

  use Connection_module
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Secondary_Continuum_Aux_module
  use Secondary_Continuum_module
  use Characteristic_Curves_Thermal_module

  implicit none

  Vec :: r
  class(realization_subsurface_type) :: realization

  PetscErrorCode :: ierr
  PetscInt :: local_id, ghosted_id

  PetscReal, pointer :: r_p(:)

  PetscInt :: icc_dn, icct_dn
           ! thermal conductivity wet constants at upstream, downstream faces.
  PetscReal :: D_dn

  PetscReal :: Res(realization%option%nflowdof)

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(th_parameter_type), pointer :: th_parameter
  type(TH_auxvar_type), pointer :: auxvars(:), auxvars_bc(:)
  type(TH_auxvar_type), pointer :: auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(sec_heat_type), pointer :: TH_sec_heat_vars(:)
  class(cc_thermal_type), pointer :: tcc_dn
  PetscReal :: v_darcy

  PetscInt :: iconn, istart, iend
  PetscInt :: sum_connection
  PetscReal :: distance_gravity
  PetscReal :: fluxe_bulk, fluxe_cond

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  th_parameter => patch%aux%TH%th_parameter
  auxvars => patch%aux%TH%auxvars
  auxvars_bc => patch%aux%TH%auxvars_bc
  auxvars_ss => patch%aux%TH%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars
  TH_sec_heat_vars => patch%aux%SC_heat%sec_heat_vars

! now assign access pointer to local variables
  call VecGetArrayF90(r,r_p,ierr);CHKERRQ(ierr)
  !print *,' Finished scattering non deriv'


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

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      icct_dn = patch%cct_id(ghosted_id)
      D_dn = th_parameter%ckwet(icct_dn)

      distance_gravity = cur_connection_set%dist(0,iconn) * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))

      icc_dn = patch%cc_id(ghosted_id)

      tcc_dn => patch%char_curves_thermal_array(icct_dn)%ptr

      call THBCFlux(boundary_condition%flow_condition%itype, &
                    boundary_condition%flow_aux_real_var(:,iconn), &
                    auxvars_bc(sum_connection), &
                    global_auxvars_bc(sum_connection), &
                    auxvars(ghosted_id), &
                    global_auxvars(ghosted_id), &
                    material_auxvars(ghosted_id), &
                    D_dn, tcc_dn, &
                    cur_connection_set%area(iconn), &
                    cur_connection_set%dist(-1:3,iconn), &
                    th_parameter%sir(1,icc_dn), &
                    option, &
                    v_darcy, &
                    fluxe_bulk, fluxe_cond, &
                    Res)

      patch%boundary_velocities(1,sum_connection) = v_darcy
      patch%boundary_flow_fluxes(:,sum_connection) = Res(:)

      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_bc(sum_connection)%mass_balance_delta(1,1) = &
          global_auxvars_bc(sum_connection)%mass_balance_delta(1,1) - Res(1)
      endif

      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      r_p(istart:iend)= r_p(istart:iend) - Res(1:option%nflowdof)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayF90(r,r_p,ierr);CHKERRQ(ierr)

end subroutine THResidualBoundaryConn

! ************************************************************************** !

subroutine THResidualAccumulation(r,realization,ierr)
  !
  ! Perform preliminary work prior to residual computation
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/06/2019
  !

  use Connection_module
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Secondary_Continuum_Aux_module
  use Secondary_Continuum_module

  implicit none

  Vec :: r
  class(realization_subsurface_type) :: realization

  PetscErrorCode :: ierr
  PetscInt :: local_id, ghosted_id

  PetscReal, pointer :: accum_p(:)
  PetscReal, pointer :: accum2_p(:)
  PetscReal, pointer :: r_p(:)

  PetscReal :: Res(realization%option%nflowdof)

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(th_parameter_type), pointer :: th_parameter
  type(TH_auxvar_type), pointer :: auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(sec_heat_type), pointer :: TH_sec_heat_vars(:)

  PetscInt :: istart, iend
  PetscReal :: vol_frac_prim

  ! secondary continuum variables
  PetscReal :: sec_dencpr
  PetscReal :: res_sec_heat
  PetscInt :: icct

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  th_parameter => patch%aux%TH%th_parameter
  auxvars => patch%aux%TH%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  TH_sec_heat_vars => patch%aux%SC_heat%sec_heat_vars

! now assign access pointer to local variables
  call VecGetArrayF90(r,r_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_accum2,accum2_p,ierr);CHKERRQ(ierr)
  !print *,' Finished scattering non deriv'

  ! Calculating volume fractions for primary and secondary continua

  vol_frac_prim = 1.d0

  ! Accumulation terms ------------------------------------
  r_p = r_p - accum_p

  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1

    if (option%use_sc) then
      vol_frac_prim = TH_sec_heat_vars(local_id)%epsilon
    endif

    call THAccumulation(auxvars(ghosted_id),global_auxvars(ghosted_id), &
                        material_auxvars(ghosted_id), &
                        th_parameter%dencpr(patch%cct_id(ghosted_id)), &
                        option,vol_frac_prim,Res)
    r_p(istart:iend) = r_p(istart:iend) + Res
    accum2_p(istart:iend) = Res
  enddo

  ! ================== Secondary continuum heat source terms ==================
  if (option%use_sc) then
  ! Secondary continuum contribution (Added by SK 06/02/2012)
  ! only one secondary continuum for now for each primary continuum node
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1

      ! secondary rho*c_p same as primary for now
      icct = patch%cct_id(ghosted_id)
      sec_dencpr = th_parameter%dencpr(icct)

      call THSecondaryHeat(TH_sec_heat_vars(local_id), &
                          global_auxvars(ghosted_id), &
                          th_parameter%ckwet(icct), &
                          sec_dencpr, &
                          option,res_sec_heat)

      r_p(iend) = r_p(iend) - res_sec_heat*material_auxvars(ghosted_id)%volume
    enddo
  endif

  call VecRestoreArrayF90(r,r_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum2,accum2_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)


end subroutine THResidualAccumulation

! ************************************************************************** !

subroutine THResidualSourceSink(r,realization,ierr)
  !
  ! Perform preliminary work prior to residual computation
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/06/2019
  !

  use Connection_module
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Secondary_Continuum_Aux_module
  use Secondary_Continuum_module

  implicit none

  Vec :: r
  class(realization_subsurface_type) :: realization

  PetscErrorCode :: ierr
  PetscInt :: i
  PetscInt :: local_id, ghosted_id

  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: xx_loc_p(:), yy_p(:)

  PetscReal :: qsrc1, esrc1

  PetscReal :: Res_src(realization%option%nflowdof)

  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(th_parameter_type), pointer :: th_parameter
  type(TH_auxvar_type), pointer :: auxvars(:), auxvars_bc(:)
  type(TH_auxvar_type), pointer :: auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(sec_heat_type), pointer :: TH_sec_heat_vars(:)
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal, pointer :: mmsrc(:)
  PetscReal :: well_status
  PetscReal :: well_factor
  PetscReal :: pressure_bh
  PetscReal :: pressure_max
  PetscReal :: pressure_min
  PetscReal :: Dq, dphi, v_darcy, ukvr

  PetscInt :: iconn, istart, iend
  PetscInt :: sum_connection

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  th_parameter => patch%aux%TH%th_parameter
  auxvars => patch%aux%TH%auxvars
  auxvars_bc => patch%aux%TH%auxvars_bc
  auxvars_ss => patch%aux%TH%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars
  TH_sec_heat_vars => patch%aux%SC_heat%sec_heat_vars

! now assign access pointer to local variables
  call VecGetArrayF90(r,r_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%flow_yy,yy_p,ierr);CHKERRQ(ierr)
  !print *,' Finished scattering non deriv'


  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit

    cur_connection_set => source_sink%connection_set

    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      iend = local_id * option%nflowdof
      istart = iend - option%nflowdof + 1
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      if (source_sink%flow_condition%rate%itype /= HET_MASS_RATE_SS .and. &
        source_sink%flow_condition%itype(1) /= WELL_SS) &
        qsrc1 = source_sink%flow_condition%rate%dataset%rarray(1)

      Res_src = 0.d0
      select case (source_sink%flow_condition%rate%itype)
        case(MASS_RATE_SS)
          qsrc1 = qsrc1 / FMWH2O ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
        case(SCALED_MASS_RATE_SS)
          qsrc1 = qsrc1 / FMWH2O * &
                      ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
            source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case(VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc1 = m^3/sec
          qsrc1 = qsrc1*global_auxvars(ghosted_id)%den(1) ! den = kmol/m^3
        case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc1 = m^3/sec
          qsrc1 = qsrc1*global_auxvars(ghosted_id)%den(1)* & ! den = kmol/m^3
            source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case(HET_MASS_RATE_SS)
          qsrc1 = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)/FMWH2O
        case(WELL_SS) ! production well, Karra 11/10/2015
          ! if node pessure is lower than the given extraction pressure,
          ! shut it down
          !  well parameter explanation
          !   1. well status. 1 injection; -1 production; 0 shut in!
          !   2. well factor [m^3],  the effective permeability [m^2/s]
          !   3. bottomhole pressure:  [Pa]
          !   4. max pressure: [Pa]
          !   5. min pressure: [Pa]
          mmsrc => source_sink%flow_condition%well%dataset%rarray

          well_status = mmsrc(1)
          well_factor = mmsrc(2)
          pressure_bh = mmsrc(3)
          pressure_max = mmsrc(4)
          pressure_min = mmsrc(5)

          ! production well (well status = -1)
          if (dabs(well_status + 1.d0) < 1.d-1) then
            if (global_auxvars(ghosted_id)%pres(1) > pressure_min) then
              Dq = well_factor
              dphi = global_auxvars(ghosted_id)%pres(1) - pressure_bh
              if (dphi >= 0.d0) then ! outflow only
                ukvr = auxvars(ghosted_id)%kvr
                if (ukvr < 1.d-20) ukvr = 0.d0
                v_darcy = 0.d0
                if (ukvr*Dq > floweps) then
                  v_darcy = Dq * ukvr * dphi
                  ! store volumetric rate for ss_fluid_fluxes()
                  qsrc1 = -1.d0*v_darcy*global_auxvars(ghosted_id)%den(1)
                endif
              endif
            endif
          endif

        case default
          write(string,*) source_sink%flow_condition%rate%itype
          option%io_buffer = &
            'TH mode source_sink%flow_condition%rate%itype = ' // &
            trim(adjustl(string)) // ', not implemented.'
      end select

      Res_src(TH_PRESSURE_DOF) = qsrc1

      esrc1 = 0.d0
      select case(source_sink%flow_condition%itype(TH_TEMPERATURE_DOF))
        case (ENERGY_RATE_SS)
          esrc1 = source_sink%flow_condition%energy_rate%dataset%rarray(1)
        case (SCALED_ENERGY_RATE_SS)
          esrc1 = source_sink%flow_condition%energy_rate%dataset%rarray(1) * &
                  source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case (HET_ENERGY_RATE_SS)
          esrc1 = source_sink%flow_aux_real_var(TWO_INTEGER,iconn)
      end select
      ! convert J/s --> MJ/s
      !geh: default internal energy units are MJ (option%scale = 1.d-6 is for
      !     J->MJ)
      Res_src(TH_TEMPERATURE_DOF) = esrc1*1.d6*option%scale

      ! Update residual term associated with T
      if (qsrc1 > 0.d0) then ! injection
        Res_src(TH_TEMPERATURE_DOF) = Res_src(TH_TEMPERATURE_DOF) + &
          qsrc1*auxvars_ss(sum_connection)%h
      else
        ! extraction
        Res_src(TH_TEMPERATURE_DOF) = Res_src(TH_TEMPERATURE_DOF) + &
          qsrc1*auxvars(ghosted_id)%h
      endif

      r_p(istart:iend) = r_p(istart:iend) - Res_src

      if (option%compute_mass_balance_new) then
        global_auxvars_ss(sum_connection)%mass_balance_delta(1,1) = &
          global_auxvars_ss(sum_connection)%mass_balance_delta(1,1) - Res_src(1)
      endif
      if (associated(patch%ss_flow_vol_fluxes)) then
        ! fluid flux [m^3/sec] = qsrc_mol [kmol/sec] / den [kmol/m^3]
        patch%ss_flow_vol_fluxes(1,sum_connection) = qsrc1 / &
                                           global_auxvars(ghosted_id)%den(1)
      endif
      if (associated(patch%ss_flow_fluxes)) then
        patch%ss_flow_fluxes(1,sum_connection) = qsrc1

      endif


    enddo
    source_sink => source_sink%next
  enddo

  ! scale the residual by the volume
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    r_p (istart:iend)= r_p(istart:iend)/material_auxvars(ghosted_id)%volume
  enddo

  if (option%use_isothermal) then
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      istart = TWO_INTEGER + (local_id-1)*option%nflowdof
      r_p(istart)=xx_loc_p(2 + (ghosted_id-1)*option%nflowdof)-yy_p(istart-1)
    enddo
  endif

  if (patch%aux%TH%inactive_cells_exist) then
    do i=1,patch%aux%TH%matrix_zeroing%n_inactive_rows
      r_p(patch%aux%TH%matrix_zeroing%inactive_rows_local(i)) = 0.d0
    enddo
  endif

  call VecRestoreArrayF90(r,r_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(field%flow_yy,yy_p,ierr);CHKERRQ(ierr)

end subroutine THResidualSourceSink

! ************************************************************************** !

subroutine THJacobian(snes,xx,A,B,realization,ierr)
  !
  ! Computes the Jacobian
  !
  ! Author: ???
  ! Date: 12/10/07
  ! Refactored by Satish Karra, LANL 06/12/2019

  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Option_module
  use Debug_module

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  class(realization_subsurface_type) :: realization
  PetscErrorCode :: ierr

  Mat :: J
  MatType :: mat_type
  PetscViewer :: viewer
  type(option_type),  pointer :: option
  PetscReal :: norm

  character(len=MAXSTRINGLENGTH) :: string

  call MatGetType(A,mat_type,ierr);CHKERRQ(ierr)
  if (mat_type == MATMFFD) then
    J = B
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  else
    J = A
  endif

  call MatZeroEntries(J,ierr);CHKERRQ(ierr)

#if 0
   call THNumericalJacobianTest(xx,realization)
#endif

  call THJacobianInternalConn(J,realization,ierr)
  call THJacobianBoundaryConn(J,realization,ierr)
  call THJacobianAccumulation(J,realization,ierr)
  call THJacobianSourceSink(J,realization,ierr)

  if (realization%debug%matview_Matrix) then
    call DebugWriteFilename(realization%debug,string,'THjacobian','', &
                            TH_ts_count,TH_ts_cut_count, &
                            TH_ni_count)
    call DebugCreateViewer(realization%debug,string,realization%option,viewer)
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

  TH_ni_count = TH_ni_count + 1

end subroutine THJacobian

! ************************************************************************** !

subroutine THJacobianInternalConn(A,realization,ierr)
  !
  ! Computes the jacobian contribution from internal flux
  !
  ! Author: ??
  ! Date: 12/13/07
  ! Refactored by Satish Karra, LANL 06/12/2019
  !


  use Connection_module
  use Option_module
  use Grid_module
  use Realization_Subsurface_class
  use Patch_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Secondary_Continuum_Aux_module
  use Saturation_Function_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Thermal_module

  Mat :: A
  class(realization_subsurface_type) :: realization

  PetscErrorCode :: ierr
  PetscInt :: icct_up, icct_dn

  PetscReal, pointer :: xx_loc_p(:)
  PetscInt :: icc_up, icc_dn
  PetscReal :: dd_up, dd_dn
  ! thermal conductivity wet constants at upstream, downstream faces.
  PetscReal :: D_up, D_dn
  PetscReal :: Dk_dry_up, Dk_dry_dn ! dry thermal conductivities
  PetscReal :: Dk_ice_up, Dk_ice_dn ! frozen soil thermal conductivities
  PetscReal :: alpha_up, alpha_dn
  PetscReal :: alpha_fr_up, alpha_fr_dn
  PetscReal :: upweight
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn

  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof), &
               Jdn(realization%option%nflowdof,realization%option%nflowdof)

  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: distance_gravity
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(th_parameter_type), pointer :: th_parameter
  type(TH_auxvar_type), pointer :: auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(saturation_function_type), pointer :: sf_up
  type(saturation_function_type), pointer :: sf_dn
  class(characteristic_curves_type), pointer :: cc_up, cc_dn
  class(cc_thermal_type), pointer :: tcc_up, tcc_dn

  character(len=MAXSTRINGLENGTH) :: string

  PetscViewer :: viewer

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  th_parameter => patch%aux%TH%th_parameter
  auxvars => patch%aux%TH%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars

  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

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

      if (patch%imat(ghosted_id_up) <= 0 .or. &
          patch%imat(ghosted_id_dn) <= 0) cycle

      if (option%flow%only_vertical_flow) then
        !geh: place second conditional within first to avoid excessive
        !     dot products when .not. option%flow%only_vertical_flow
        if (dot_product(cur_connection_set%dist(1:3,iconn),unit_z) < &
            1.d-10) cycle
      endif

      local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
      local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping

      fraction_upwind = cur_connection_set%dist(-1,iconn)
      distance = cur_connection_set%dist(0,iconn)
      ! distance = scalar - magnitude of distance
      ! gravity = vector(3)
      ! dist(1:3,iconn) = vector(3) - unit vector
      distance_gravity = distance * &
                         dot_product(option%gravity, &
                                     cur_connection_set%dist(1:3,iconn))
      dd_up = distance*fraction_upwind
      dd_dn = distance-dd_up ! should avoid truncation error
      ! upweight could be calculated as 1.d0-fraction_upwind
      ! however, this introduces ever so slight error causing pflow-overhaul not
      ! to match pflow-orig.  This can be changed to 1.d0-fraction_upwind
      upweight = dd_dn/(dd_up+dd_dn)

      icct_up = patch%cct_id(ghosted_id_up)
      icct_dn = patch%cct_id(ghosted_id_dn)

      D_up = th_parameter%ckwet(icct_up)
      D_dn = th_parameter%ckwet(icct_dn)

      Dk_dry_up = th_parameter%ckdry(icct_up)
      Dk_dry_dn = th_parameter%ckdry(icct_dn)

      alpha_up = th_parameter%alpha(icct_up)
      alpha_dn = th_parameter%alpha(icct_dn)

      icc_up = patch%cc_id(ghosted_id_up)
      icc_dn = patch%cc_id(ghosted_id_dn)

      if (option%flow%th_freezing) then
         Dk_ice_up = th_parameter%ckfrozen(icct_up)
         DK_ice_dn = th_parameter%ckfrozen(icct_dn)

         alpha_fr_up = th_parameter%alpha_fr(icct_up)
         alpha_fr_dn = th_parameter%alpha_fr(icct_dn)

         sf_up => patch%saturation_function_array(icc_up)%ptr
         sf_dn => patch%saturation_function_array(icc_dn)%ptr
      else
         Dk_ice_up = Dk_dry_up
         Dk_ice_dn = Dk_dry_dn

         alpha_fr_up = alpha_up
         alpha_fr_dn = alpha_dn

         cc_up => patch%characteristic_curves_array(icc_up)%ptr
         cc_dn => patch%characteristic_curves_array(icc_dn)%ptr
      endif

      tcc_up => patch%char_curves_thermal_array(icct_up)%ptr
      tcc_dn => patch%char_curves_thermal_array(icct_dn)%ptr

      call THFluxDerivative(auxvars(ghosted_id_up), &
                            global_auxvars(ghosted_id_up), &
                            material_auxvars(ghosted_id_up), &
                            D_up, &
                            icct_up, &
                            auxvars(ghosted_id_dn), &
                            global_auxvars(ghosted_id_dn), &
                            material_auxvars(ghosted_id_dn), &
                            D_dn, &
                            icct_dn, &
                            cur_connection_set%area(iconn), &
                            cur_connection_set%dist(-1:3,iconn), &
                            upweight, &
                            th_parameter%sir(1,icc_up), &
                            th_parameter%sir(1,icc_dn), &
                            option, &
                            sf_up,sf_dn,cc_up,cc_dn,tcc_up,tcc_dn, &
                            Dk_dry_up,Dk_dry_dn, &
                            Dk_ice_up,Dk_ice_dn, &
                            alpha_up,alpha_dn,alpha_fr_up,alpha_fr_dn, &
                            th_parameter, &
                            Jup,Jdn)

!  scale by the volume of the cell

      if (local_id_up > 0) then
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1, &
                                      Jup/material_auxvars(ghosted_id_up)% &
                                        volume, &
                                      ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                                      Jdn/material_auxvars(ghosted_id_up)% &
                                        volume, &
                                      ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
      if (local_id_dn > 0) then
        Jup = -Jup
        Jdn = -Jdn

        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1, &
                                      Jdn/material_auxvars(ghosted_id_dn)% &
                                        volume, &
                                      ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesBlockedLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                                      Jup/material_auxvars(ghosted_id_dn)% &
                                        volume, &
                                      ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  if (realization%debug%matview_Matrix_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_flux'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif


end subroutine THJacobianInternalConn

! ************************************************************************** !

subroutine THJacobianBoundaryConn(A,realization,ierr)
  !
  ! Computes the jacobian contribution from boundary flux
  !
  ! Author: ??
  ! Date: 12/13/07
  ! Refactored by Satish Karra, LANL 06/12/2019
  !


  use Connection_module
  use Option_module
  use Grid_module
  use Realization_Subsurface_class
  use Patch_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Secondary_Continuum_Aux_module
  use Saturation_Function_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Thermal_module

  Mat :: A
  class(realization_subsurface_type) :: realization

  PetscErrorCode :: ierr
  PetscInt :: icct_dn

  PetscReal, pointer :: xx_loc_p(:)
  PetscInt :: icc_dn
  ! thermal conductivity wet constants at upstream, downstream faces.
  PetscReal :: D_dn
  PetscReal :: Dk_dry_dn ! dry thermal conductivities
  PetscReal :: Dk_ice_dn ! frozen soil thermal conductivities
  PetscReal :: alpha_dn
  PetscReal :: alpha_fr_dn
  PetscInt :: local_id, ghosted_id

  PetscReal :: Jdn(realization%option%nflowdof,realization%option%nflowdof)

  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(th_parameter_type), pointer :: th_parameter
  type(TH_auxvar_type), pointer :: auxvars_bc(:), auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:), global_auxvars_bc(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(saturation_function_type), pointer :: sf_dn
  class(characteristic_curves_type), pointer :: cc_dn
  class(cc_thermal_type), pointer :: tcc_dn

  character(len=MAXSTRINGLENGTH) :: string

  PetscViewer :: viewer

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  th_parameter => patch%aux%TH%th_parameter
  auxvars => patch%aux%TH%auxvars
  auxvars_bc => patch%aux%TH%auxvars_bc
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  material_auxvars => patch%aux%Material%auxvars


  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

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

      if (ghosted_id<=0) then
        print *, "Wrong boundary node index... STOP!!!"
        stop
      endif

      icct_dn   = patch%cct_id(ghosted_id)
      D_dn      = th_parameter%ckwet(icct_dn)
      Dk_dry_dn = th_parameter%ckdry(icct_dn)
      alpha_dn  = th_parameter%alpha(icct_dn)

      icc_dn = patch%cc_id(ghosted_id)

      if (option%flow%th_freezing) then
         DK_ice_dn = th_parameter%ckfrozen(icct_dn)
         alpha_fr_dn = th_parameter%alpha_fr(icct_dn)

         sf_dn => patch%saturation_function_array(icc_dn)%ptr
      else
         Dk_ice_dn = Dk_dry_dn
         alpha_fr_dn = alpha_dn

         cc_dn => patch%characteristic_curves_array(icc_dn)%ptr
      endif

      tcc_dn => patch%char_curves_thermal_array(icct_dn)%ptr

      call THBCFluxDerivative(boundary_condition%flow_condition%itype, &
                              boundary_condition%flow_aux_real_var(:,iconn), &
                              auxvars_bc(sum_connection), &
                              global_auxvars_bc(sum_connection), &
                              auxvars(ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              D_dn, &
                              cur_connection_set%area(iconn), &
                              cur_connection_set%dist(-1:3,iconn), &
                              th_parameter%sir(1,icc_dn), &
                              option, &
                              sf_dn,cc_dn, tcc_dn, &
                              Dk_dry_dn,Dk_ice_dn, &
                              Jdn)
      Jdn = -Jdn

      !  scale by the volume of the cell
      Jdn = Jdn/material_auxvars(ghosted_id)%volume

      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jdn, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)

    enddo
    boundary_condition => boundary_condition%next
  enddo

  if (realization%debug%matview_Matrix_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_bcflux'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

end subroutine THJacobianBoundaryConn

! ************************************************************************** !

subroutine THJacobianAccumulation(A,realization,ierr)
  !
  ! Computes the jacobian contribution from accumulation term
  !
  ! Author: ??
  ! Date: 12/13/07
  ! Refactored by Satish Karra, LANL 06/12/2019
  !


  use Connection_module
  use Option_module
  use Grid_module
  use Realization_Subsurface_class
  use Patch_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Secondary_Continuum_Aux_module
  use Saturation_Function_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Thermal_module

  Mat :: A
  class(realization_subsurface_type) :: realization

  PetscErrorCode :: ierr

  PetscReal, pointer :: xx_loc_p(:)
  PetscInt :: local_id, ghosted_id

  PetscReal :: Jup(realization%option%nflowdof,realization%option%nflowdof)

  PetscInt :: istart, iend
  PetscInt :: icc
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(th_parameter_type), pointer :: th_parameter
  type(TH_auxvar_type), pointer :: auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  type(saturation_function_type), pointer :: sat_func
  class(characteristic_curves_type), pointer :: characteristic_curves
  class(cc_thermal_type), pointer :: thermal_cc



  type(sec_heat_type), pointer :: sec_heat_vars(:)
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: icct

  PetscViewer :: viewer
  PetscReal :: vol_frac_prim

  ! secondary continuum variables
  PetscReal :: jac_sec_heat

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  th_parameter => patch%aux%TH%th_parameter
  auxvars => patch%aux%TH%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars

  sec_heat_vars => patch%aux%SC_heat%sec_heat_vars

  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

  vol_frac_prim = 1.d0

  ! Accumulation terms ------------------------------------
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    ! Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    iend = local_id*option%nflowdof
    istart = iend-option%nflowdof+1
    icc = patch%cc_id(ghosted_id)

    if (option%use_sc) then
      vol_frac_prim = sec_heat_vars(local_id)%epsilon
    endif

    icct = patch%cct_id(ghosted_id)

    if (option%flow%th_freezing) then
      sat_func => patch%saturation_function_array(icc)%ptr
    else
      characteristic_curves => patch%characteristic_curves_array(icc)%ptr
    endif

    thermal_cc => patch%char_curves_thermal_array(icct)%ptr

    call THAccumDerivative(auxvars(ghosted_id),global_auxvars(ghosted_id), &
                            material_auxvars(ghosted_id), &
                            th_parameter%dencpr(icct), &
                            th_parameter, icct, option, &
                            sat_func, characteristic_curves, &
                            thermal_cc, &
                            vol_frac_prim,Jup)

    if (option%use_sc) then
      call THSecondaryHeatJacobian(sec_heat_vars(local_id), &
                        th_parameter%ckwet(icct), &
                        th_parameter%dencpr(icct), &
                        option,jac_sec_heat)

      Jup(option%nflowdof,2) = Jup(option%nflowdof,2) - &
                               jac_sec_heat*material_auxvars(ghosted_id)%volume
    endif

    ! scale by the volume of the cell
    Jup = Jup/material_auxvars(ghosted_id)%volume

    call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jup, &
                                  ADD_VALUES,ierr);CHKERRQ(ierr)
  enddo


  if (realization%debug%matview_Matrix_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_accum'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif


  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

end subroutine THJacobianAccumulation

! ************************************************************************** !

subroutine THJacobianSourceSink(A,realization,ierr)
  !
  ! Computes the jacobian contribution from source sink
  !
  ! Author: ??
  ! Date: 12/13/07
  ! Refactored by Satish Karra, LANL 06/12/2019
  !


  use Connection_module
  use Option_module
  use Grid_module
  use Realization_Subsurface_class
  use Patch_module
  use Coupler_module
  use Field_module
  use Debug_module
  use Secondary_Continuum_Aux_module

  Mat :: A
  class(realization_subsurface_type) :: realization

  PetscErrorCode :: ierr

  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: qsrc1
  PetscInt :: local_id, ghosted_id
  PetscReal :: f_up
  PetscReal :: Jsrc(realization%option%nflowdof,realization%option%nflowdof)

  PetscInt :: istart

  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(th_parameter_type), pointer :: th_parameter
  type(TH_auxvar_type), pointer :: auxvars(:), auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  character(len=MAXSTRINGLENGTH) :: string

  PetscViewer :: viewer

  patch => realization%patch
  grid => patch%grid
  option => realization%option
  field => realization%field

  th_parameter => patch%aux%TH%th_parameter
  auxvars => patch%aux%TH%auxvars
  auxvars_ss => patch%aux%TH%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars

  call VecGetArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

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

      if (source_sink%flow_condition%rate%itype /= HET_MASS_RATE_SS .and. &
        source_sink%flow_condition%itype(1) /= WELL_SS) &
        qsrc1 = source_sink%flow_condition%rate%dataset%rarray(1)

      select case (source_sink%flow_condition%rate%itype)
        case(MASS_RATE_SS)
          qsrc1 = qsrc1 / FMWH2O ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
        case(SCALED_MASS_RATE_SS)
          qsrc1 = qsrc1 / FMWH2O * &
               ! [kg/s -> kmol/s; fmw -> g/mol = kg/kmol]
            source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case(VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc1 = m^3/sec
          qsrc1 = qsrc1*global_auxvars(ghosted_id)%den(1) ! den = kmol/m^3
        case(SCALED_VOLUMETRIC_RATE_SS)  ! assume local density for now
          ! qsrc1 = m^3/sec
          qsrc1 = qsrc1*global_auxvars(ghosted_id)%den(1)* & ! den = kmol/m^3
                   source_sink%flow_aux_real_var(ONE_INTEGER,iconn)
        case(HET_MASS_RATE_SS)
          qsrc1 = source_sink%flow_aux_real_var(ONE_INTEGER,iconn)/FMWH2O
        case default
          write(string,*) source_sink%flow_condition%rate%itype
          option%io_buffer = &
            'TH mode source_sink%flow_condition%rate%itype = ' // &
            trim(adjustl(string)) // ', not implemented.'
      end select

      Jsrc = 0.d0

      if (qsrc1 > 0.d0) then ! injection
        Jsrc(TH_TEMPERATURE_DOF,TH_PRESSURE_DOF) = &
          -qsrc1*auxvars_ss(sum_connection)%dh_dp
        ! since tsrc1 is prescribed, there is no derivative
        ! dresT_dT = -qsrc1*hw_dT
        istart = ghosted_id*option%nflowdof
      else
        ! extraction
        Jsrc(TH_TEMPERATURE_DOF,TH_PRESSURE_DOF) = &
          -qsrc1*auxvars(ghosted_id)%dh_dp
        Jsrc(TH_TEMPERATURE_DOF,TH_TEMPERATURE_DOF) = &
          -qsrc1*auxvars(ghosted_id)%dh_dT
        istart = ghosted_id*option%nflowdof
      endif

      ! scale by the volume of the cell
      Jsrc = Jsrc/material_auxvars(ghosted_id)%volume

      call MatSetValuesBlockedLocal(A,1,ghosted_id-1,1,ghosted_id-1,Jsrc, &
                                    ADD_VALUES,ierr);CHKERRQ(ierr)


    enddo
    source_sink => source_sink%next
  enddo

  if (realization%debug%matview_Matrix_detailed) then
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    string = 'jacobian_srcsink'
    call DebugCreateViewer(realization%debug,string,option,viewer)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif


  call VecRestoreArrayF90(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

! zero out isothermal and inactive cells
#ifdef ISOTHERMAL_MODE_DOES_NOT_WORK
  zero = 0.d0
  call MatZeroRowsLocal(A,patch%aux%%matrix_zeroing%n_inactive_rows, &
                        patch%aux%%matrix_zeroing% &
                          inactive_rows_local_ghosted, &
                        zero,PETSC_NULL_VEC,PETSC_NULL_VEC, &
                        ierr);CHKERRQ(ierr)
  do i=1, patch%aux%TH%matrix_zeroing%n_inactive_rows
    ii = mod(patch%aux%TH%matrix_zeroing%inactive_rows_local(i),option%nflowdof)
    ip1 = patch%aux%TH%matrix_zeroing%inactive_rows_local_ghosted(i)
    if (ii == 0) then
      ip2 = ip1-1
    else if (ii == option%nflowdof-1) then
      ip2 = ip1+1
    else
      ip2 = ip1
    endif
    call MatSetValuesLocal(A,1,ip1,1,ip2,1.d0,INSERT_VALUES, &
                           ierr);CHKERRQ(ierr)
  enddo

  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
#else
  if (patch%aux%TH%inactive_cells_exist) then
    f_up = 1.d0
    call MatZeroRowsLocal(A,patch%aux%TH%matrix_zeroing%n_inactive_rows, &
                          patch%aux%TH%matrix_zeroing% &
                            inactive_rows_local_ghosted, &
                          f_up,PETSC_NULL_VEC,PETSC_NULL_VEC, &
                          ierr);CHKERRQ(ierr)
  endif
#endif
end subroutine THJacobianSourceSink

! ************************************************************************** !

subroutine THMaxChange(realization,dpmax,dtmpmax)
  !
  ! Computes the maximum change in the solution vector
  !
  ! Author: ???
  ! Date: 01/15/08
  !

  use Realization_Subsurface_class
  use Option_module
  use Field_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(field_type), pointer :: field

  PetscReal :: dpmax, dtmpmax
  PetscErrorCode :: ierr

  option => realization%option
  field => realization%field

  dpmax = 0.d0
  dtmpmax = 0.d0

  call VecWAXPY(field%flow_dxx,-1.d0,field%flow_xx,field%flow_yy, &
                ierr);CHKERRQ(ierr)
  call VecStrideNorm(field%flow_dxx,ZERO_INTEGER,NORM_INFINITY,dpmax, &
                     ierr);CHKERRQ(ierr)
  call VecStrideNorm(field%flow_dxx,ONE_INTEGER,NORM_INFINITY,dtmpmax, &
                     ierr);CHKERRQ(ierr)

end subroutine THMaxChange

! ************************************************************************** !

subroutine THResidualToMass(realization)
  !
  ! Computes mass balance from residual equation
  !
  ! Author: ???
  ! Date: 12/10/07
  !

  use Realization_Subsurface_class
  use Patch_module
  use Discretization_module
  use Field_module
  use Option_module
  use Grid_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(field_type), pointer :: field
  type(patch_type), pointer :: cur_patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option

  PetscReal, pointer :: mass_balance_p(:)
  type(TH_auxvar_type), pointer :: auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  PetscErrorCode :: ierr
  PetscInt :: local_id, ghosted_id
  PetscInt :: istart

  option => realization%option
  field => realization%field

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit

    grid => cur_patch%grid
    auxvars => cur_patch%aux%TH%auxvars

    call VecGetArrayF90(field%flow_ts_mass_balance,mass_balance_p, &
                        ierr);CHKERRQ(ierr)

    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (cur_patch%imat(ghosted_id) <= 0) cycle

      istart = (ghosted_id-1)*option%nflowdof+1
      mass_balance_p(istart) = mass_balance_p(istart)/ &
                                global_auxvars(ghosted_id)%den(1)* &
                                global_auxvars(ghosted_id)%den_kg(1)
    enddo

    call VecRestoreArrayF90(field%flow_ts_mass_balance,mass_balance_p, &
                            ierr);CHKERRQ(ierr)

    cur_patch => cur_patch%next
  enddo

end subroutine THResidualToMass

! ************************************************************************** !

function THGetTecplotHeader(realization,icolumn)
  !
  ! THLiteGetTecplotHeader: Returns TH contribution to
  ! Tecplot file header
  !
  ! Author: ???
  ! Date: 02/13/08
  !

  use Realization_Subsurface_class
  use Option_module
  use Field_module

  implicit none

  character(len=MAXSTRINGLENGTH) :: THGetTecplotHeader
  class(realization_subsurface_type) :: realization
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
    write(string2,'('',"'',i2,''-Sl"'')') icolumn
  else
    write(string2,'('',"Sl"'')')
  endif
  string = trim(string) // trim(string2)

  if (option%flow%th_freezing) then
     if (icolumn > -1) then
        icolumn = icolumn + 1
        write(string2,'('',"'',i2,''-Sg"'')') icolumn
     else
        write(string2,'('',"Sg"'')')
     endif
     string = trim(string) // trim(string2)

     if (icolumn > -1) then
        icolumn = icolumn + 1
        write(string2,'('',"'',i2,''-Si"'')') icolumn
     else
        write(string2,'('',"Si"'')')
     endif
     string = trim(string) // trim(string2)

     if (icolumn > -1) then
        icolumn = icolumn + 1
        write(string2,'('',"'',i2,''-deni"'')') icolumn
     else
        write(string2,'('',"deni"'')')
     endif
     string = trim(string) // trim(string2)
  endif

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-denl"'')') icolumn
  else
    write(string2,'('',"denl"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-Ul"'')') icolumn
  else
    write(string2,'('',"Ul"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-visl"'')') icolumn
  else
    write(string2,'('',"visl"'')')
  endif
  string = trim(string) // trim(string2)

  if (icolumn > -1) then
    icolumn = icolumn + 1
    write(string2,'('',"'',i2,''-mobilityl"'')') icolumn
  else
    write(string2,'('',"mobilityl"'')')
  endif
  string = trim(string) // trim(string2)

  do i=1,option%nflowspec
    if (icolumn > -1) then
      icolumn = icolumn + 1
      write(string2,'('',"'',i2,''-Xl('',i2,'')"'')') icolumn,i
    else
      write(string2,'('',"Xl('',i2,'')"'')') i
    endif
    string = trim(string) // trim(string2)
  enddo

  THGetTecplotHeader = string

end function THGetTecplotHeader

! ************************************************************************** !

subroutine THSetPlotVariables(realization,list)
  !
  ! Adds variables to be printed to list
  !
  ! Author: Glenn Hammond
  ! Date: 10/15/12
  !

  use Realization_Subsurface_class
  use Output_Aux_module
  use Variables_module
  use Material_Aux_module
  use Option_module

  implicit none

  class(realization_subsurface_type) :: realization
  type(output_variable_list_type), pointer :: list

  character(len=MAXWORDLENGTH) :: name, units
  type(option_type), pointer :: option

  option => realization%option

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

    name = 'Liquid Density'
    units = 'kg/m^3'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                LIQUID_DENSITY)

    name = 'Liquid Energy'
    units = 'kJ/mol'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                LIQUID_ENERGY)

    name = 'Liquid Viscosity'
    units = 'Pa.s'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                LIQUID_VISCOSITY)

    name = 'Liquid Mobility'
    units = '1/Pa.s'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                LIQUID_MOBILITY)

  endif

  if (list%energy_vars) then

    name = 'Temperature'
    units = 'C'
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                TEMPERATURE)

  endif

  if (option%flow%th_freezing) then

    if (th_ice_model /= DALL_AMICO) then
      name = 'Gas Saturation'
      units = ''
      call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
          GAS_SATURATION)
    endif

    name = 'Ice Saturation'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
        ICE_SATURATION)

    name = 'Ice Density'
    units = 'kg/m^3'
    call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
        ICE_DENSITY)

  endif

  if (soil_compressibility_index > 0) then

    name = 'Porosity'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
                                 POROSITY)

  endif
! name = 'Phase'
! units = ''
! output_variable%iformat = 1 ! integer
! call OutputVariableAddToList(list,name,OUTPUT_GENERIC,units, &
!                              PHASE)

end subroutine THSetPlotVariables

! ************************************************************************** !

subroutine THSecondaryHeat(sec_heat_vars,global_auxvar, &
                            therm_conductivity,dencpr, &
                            option,res_heat)
  !
  ! Calculates the source term contribution due to secondary
  ! continuum in the primary continuum residual
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/2/12
  !

  use Option_module
  use Global_Aux_module
  use Secondary_Continuum_Aux_module

  implicit none

  type(sec_heat_type) :: sec_heat_vars
  type(global_auxvar_type) :: global_auxvar
  type(option_type) :: option
  PetscReal :: coeff_left(sec_heat_vars%ncells)
  PetscReal :: coeff_diag(sec_heat_vars%ncells)
  PetscReal :: coeff_right(sec_heat_vars%ncells)
  PetscReal :: rhs(sec_heat_vars%ncells)
  PetscReal :: area(sec_heat_vars%ncells)
  PetscReal :: vol(sec_heat_vars%ncells)
  PetscReal :: dm_plus(sec_heat_vars%ncells)
  PetscReal :: dm_minus(sec_heat_vars%ncells)
  PetscInt :: i, ngcells
  PetscReal :: area_fm
  PetscReal :: alpha, therm_conductivity, dencpr
  PetscReal :: temp_primary_node
  PetscReal :: m
  PetscReal :: temp_current_N
  PetscReal :: res_heat

  ngcells = sec_heat_vars%ncells
  area = sec_heat_vars%area
  vol = sec_heat_vars%vol
  dm_plus = sec_heat_vars%dm_plus
  dm_minus = sec_heat_vars%dm_minus
  area_fm = sec_heat_vars%interfacial_area
  temp_primary_node = global_auxvar%temp

  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  rhs = 0.d0

  alpha = option%flow_dt*therm_conductivity/dencpr


! Setting the coefficients
  do i = 2, ngcells-1
    coeff_left(i) = -alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i))
    coeff_diag(i) = alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i)) + &
                    alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i)) + 1.d0
    coeff_right(i) = -alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i))
  enddo

  coeff_diag(1) = alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1)) + 1.d0
  coeff_right(1) = -alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1))

  coeff_left(ngcells) = -alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells))
  coeff_diag(ngcells) = alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells)) &
                       + alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells)) &
                       + 1.d0

  ! secondary continuum values from previous time step
  rhs = sec_heat_vars%sec_temp
  rhs(ngcells) = rhs(ngcells) + &
                 alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells))* &
                 temp_primary_node

  ! Thomas algorithm for tridiagonal system
  ! Forward elimination
  do i = 2, ngcells
    m = coeff_left(i)/coeff_diag(i-1)
    coeff_diag(i) = coeff_diag(i) - m*coeff_right(i-1)
    rhs(i) = rhs(i) - m*rhs(i-1)
  enddo

  ! Back substitution
  ! We only need the temperature at the outer-most node (closest to
  ! primary node)
  temp_current_N = rhs(ngcells)/coeff_diag(ngcells)

  ! Calculate the coupling term
  res_heat = area_fm*therm_conductivity*(temp_current_N - temp_primary_node)/ &
             dm_plus(ngcells)

end subroutine THSecondaryHeat

! ************************************************************************** !

subroutine THSecondaryHeatJacobian(sec_heat_vars, &
                                    therm_conductivity, &
                                    dencpr, &
                                    option,jac_heat)
  !
  ! Calculates the source term jacobian contribution
  ! due to secondary continuum in the primary continuum residual
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/6/12
  !

  use Option_module
  use Global_Aux_module
  use Secondary_Continuum_Aux_module

  implicit none

  type(sec_heat_type) :: sec_heat_vars
  type(option_type) :: option
  PetscReal :: coeff_left(sec_heat_vars%ncells)
  PetscReal :: coeff_diag(sec_heat_vars%ncells)
  PetscReal :: coeff_right(sec_heat_vars%ncells)
  PetscReal :: rhs(sec_heat_vars%ncells)
  PetscReal :: area(sec_heat_vars%ncells)
  PetscReal :: vol(sec_heat_vars%ncells)
  PetscReal :: dm_plus(sec_heat_vars%ncells)
  PetscReal :: dm_minus(sec_heat_vars%ncells)
  PetscInt :: i, ngcells
  PetscReal :: area_fm
  PetscReal :: alpha, therm_conductivity, dencpr
  PetscReal :: m
  PetscReal :: Dtemp_N_Dtemp_prim
  PetscReal :: jac_heat

  ngcells = sec_heat_vars%ncells
  area = sec_heat_vars%area
  vol = sec_heat_vars%vol
  dm_plus = sec_heat_vars%dm_plus
  area_fm = sec_heat_vars%interfacial_area
  dm_minus = sec_heat_vars%dm_minus

  coeff_left = 0.d0
  coeff_diag = 0.d0
  coeff_right = 0.d0
  rhs = 0.d0

  alpha = option%flow_dt*therm_conductivity/dencpr

! Setting the coefficients
  do i = 2, ngcells-1
    coeff_left(i) = -alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i))
    coeff_diag(i) = alpha*area(i-1)/((dm_minus(i) + dm_plus(i-1))*vol(i)) + &
                    alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i)) + 1.d0
    coeff_right(i) = -alpha*area(i)/((dm_minus(i+1) + dm_plus(i))*vol(i))
  enddo

  coeff_diag(1) = alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1)) + 1.d0
  coeff_right(1) = -alpha*area(1)/((dm_minus(2) + dm_plus(1))*vol(1))

  coeff_left(ngcells) = -alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells))
  coeff_diag(ngcells) = alpha*area(ngcells-1)/ &
                       ((dm_minus(ngcells) + dm_plus(ngcells-1))*vol(ngcells)) &
                       + alpha*area(ngcells)/(dm_plus(ngcells)*vol(ngcells)) &
                       + 1.d0

  ! Thomas algorithm for tridiagonal system
  ! Forward elimination
  do i = 2, ngcells
    m = coeff_left(i)/coeff_diag(i-1)
    coeff_diag(i) = coeff_diag(i) - m*coeff_right(i-1)
    ! We do not have to calculate rhs terms
  enddo

  ! We need the temperature derivative at the outer-most node (closest
  ! to primary node)
  Dtemp_N_Dtemp_prim = 1.d0/coeff_diag(ngcells)*alpha*area(ngcells)/ &
                       (dm_plus(ngcells)*vol(ngcells))

  ! Calculate the jacobian term
  jac_heat = area_fm*therm_conductivity*(Dtemp_N_Dtemp_prim - 1.d0)/ &
             dm_plus(ngcells)


end subroutine THSecondaryHeatJacobian

! ************************************************************************** !

function THInitGuessCheck(xx, option)
  !
  ! Checks if the initial guess is valid.
  ! Note: Only implemented for DALL_AMICO formulation.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/04/2014
  !
  use Option_module

  Vec :: xx
  type(option_type), pointer :: option

  PetscInt :: THInitGuessCheck
  PetscInt :: idx
  PetscReal :: pres_min, pres_max
  PetscReal :: temp_min, temp_max
  PetscInt :: ipass, ipass0
  PetscErrorCode :: ierr

  ipass = 1

  if (th_ice_model /= DALL_AMICO) then
    THInitGuessCheck = ipass
    return
  endif

  call VecStrideMin(xx,ZERO_INTEGER,idx,pres_min,ierr);CHKERRQ(ierr)
  call VecStrideMin(xx,ONE_INTEGER,idx,temp_min,ierr);CHKERRQ(ierr)
  call VecStrideMax(xx,ZERO_INTEGER,idx,pres_max,ierr);CHKERRQ(ierr)
  call VecStrideMax(xx,ONE_INTEGER,idx,temp_max,ierr);CHKERRQ(ierr)

  if (pres_min < -1.d10 .or. pres_min > 1.d10 .or. &
      temp_min < -100.d0 .or. temp_max > 100.d0) then
      ipass = -1
  endif

   call MPI_Barrier(option%mycomm,ierr);CHKERRQ(ierr)
   if (option%comm%mycommsize>1)then
      call MPI_Allreduce(ipass,ipass0,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                         option%mycomm,ierr);CHKERRQ(ierr)
      if (ipass0 < option%comm%mycommsize) ipass=-1
   endif
   THInitGuessCheck = ipass

end function THInitGuessCheck

! ************************************************************************** !

subroutine EnergyToTemperatureBisection(T,TL,TR,h,energy,Cwi,Pr,option)
  !
  ! Solves the following nonlinear equation using the bisection method
  !
  ! R(T) = rho(T) Cwi hw T - energy = 0
  !
  ! Author: Nathan Collier, ORNL
  ! Date: 11/2014
  !
  use EOS_Water_module
  use Option_module

  implicit none

  PetscReal :: T,TL,TR,h,energy,Cwi,Pr
  type(option_type), pointer :: option

  PetscReal :: Tp,rho,rho_t,f,fR,fL,rtol
  PetscInt :: iter,niter
  PetscBool :: found
  PetscErrorCode :: ierr

  call EOSWaterDensity(TR,Pr,rho,rho_T,ierr)
  fR = rho*Cwi*h*(TR+273.15d0) - energy
  call EOSWaterDensity(TL,Pr,rho,rho_T,ierr)
  fL = rho*Cwi*h*(TL+273.15d0) - energy

  if (fL*fR > 0.d0) then
     print *,"[TL,TR] = ",TL,TR
     print *,"[fL,fR] = ",fL,fR
     write(option%io_buffer,'("th.F90: EnergyToTemperatureBisection -->&
                              & root is not bracketed")')
     call PrintErrMsg(option)
  endif

  T = 0.5d0*(TL+TR)
  call EOSWaterDensity(T,Pr,rho,rho_T,ierr)
  f = rho*Cwi*h*(T+273.15d0) - energy

  found = PETSC_FALSE
  niter = 200
  rtol  = 1.d-6
  do iter = 1,niter
     Tp = T
     if (fL*f < 0.d0) then
        TR = T
     else
        TL = T
     endif

     T = 0.5d0*(TL+TR)

     call EOSWaterDensity(T,Pr,rho,rho_T,ierr)
     f = rho*Cwi*h*(T+273.15d0) - energy

     if (abs((T-Tp)/(T+273.15d0)) < rtol) then
        found = PETSC_TRUE
        exit
     endif
  enddo

  if (found .eqv. PETSC_FALSE) then
     print *,"[TL,T,TR] = ",TL,T,TR
     write(option%io_buffer,'("th.F90: EnergyToTemperatureBisection -->&
                              & root not found!")')
     call PrintErrMsg(option)
  endif

end subroutine EnergyToTemperatureBisection

! ************************************************************************** !

subroutine THDestroy(patch)
  !
  ! Deallocates variables associated with Richard
  !
  ! Author: ???
  ! Date: 02/14/08
  !

  use Patch_module

  implicit none

  type(patch_type) :: patch

  ! need to free array in aux vars
  call THAuxDestroy(patch%aux%TH)

end subroutine THDestroy

end module TH_module
