module Patch_module

#include "petsc/finclude/petscvec.h"
  use petscvec

  use PFLOTRAN_Constants_module
  use Option_module
  use Grid_module
  use Coupler_module
  use Observation_module
  use Integral_Flux_module
  use Strata_module
  use Region_module
  use Reaction_Base_module
  use Reaction_Aux_module
  use NW_Transport_Aux_module
  use Dataset_Base_class
  use Material_module
  use Field_module
  use Saturation_Function_module
  use Characteristic_Curves_Thermal_module
  use Characteristic_Curves_module
  use Material_Transform_module
  use Auxiliary_module

  use General_Aux_module
  use Hydrate_Aux_module
  use TH_Aux_module

  implicit none

  private

  type, public :: patch_type

    PetscInt :: id

    ! These arrays will be used by all modes, mode-specific arrays should
    ! go in the auxiliary data stucture for that mode
    PetscInt, pointer :: imat(:)
    PetscInt, pointer :: imat_internal_to_external(:)
    PetscInt, pointer :: cc_id(:)    ! characteristic curves id
    PetscInt, pointer :: cct_id(:)   ! thermal characteristic curves id
    PetscInt, pointer :: mtf_id(:)   ! material transform id

    PetscReal, pointer :: internal_velocities(:,:)
    PetscReal, pointer :: boundary_velocities(:,:)
    PetscReal, pointer :: internal_tran_coefs(:,:,:)
    PetscReal, pointer :: boundary_tran_coefs(:,:,:)
    PetscReal, pointer :: internal_flow_fluxes(:,:)
    PetscReal, pointer :: boundary_flow_fluxes(:,:)
    ! fluid fluxes in moles/sec
    PetscReal, pointer :: ss_flow_fluxes(:,:)
    ! volumetric flux (m^3/sec) for liquid phase needed for transport
    PetscReal, pointer :: ss_flow_vol_fluxes(:,:)
    PetscReal, pointer :: internal_tran_fluxes(:,:)
    PetscReal, pointer :: boundary_tran_fluxes(:,:)
    PetscReal, pointer :: ss_tran_fluxes(:,:)

    ! for upwind direction in multiphase flow
    PetscInt, pointer :: flow_upwind_direction(:,:)
    PetscInt, pointer :: flow_upwind_direction_bc(:,:)

    type(grid_type), pointer :: grid

    type(region_list_type), pointer :: region_list

    type(coupler_list_type), pointer :: boundary_condition_list
    type(coupler_list_type), pointer :: initial_condition_list
    type(coupler_list_type), pointer :: source_sink_list

    type(material_property_type), pointer :: material_properties
    type(material_property_ptr_type), pointer :: material_property_array(:)
    type(saturation_function_type), pointer :: saturation_functions
    type(saturation_function_ptr_type), pointer :: saturation_function_array(:)
    class(characteristic_curves_type), pointer :: characteristic_curves
    type(characteristic_curves_ptr_type), pointer :: characteristic_curves_array(:)
    class(cc_thermal_type), pointer :: characteristic_curves_thermal
    type(cc_thermal_ptr_type), pointer :: char_curves_thermal_array(:)
    class(material_transform_type), pointer :: material_transform
    type(material_transform_ptr_type), pointer :: material_transform_array(:)

    type(strata_list_type), pointer :: strata_list
    type(observation_list_type), pointer :: observation_list
    type(integral_flux_list_type), pointer :: integral_flux_list

    ! Pointers to objects in mother realization object
    type(field_type), pointer :: field
    class(dataset_base_type), pointer :: datasets
    class(reaction_rt_type), pointer :: reaction
    class(reaction_nw_type), pointer :: reaction_nw
    class(reaction_base_type), pointer :: reaction_base

    type(auxiliary_type) :: aux

    type(patch_type), pointer :: next

  end type patch_type

  ! pointer data structure required for making an array of patch pointers in F90
  type, public :: patch_ptr_type
    type(patch_type), pointer :: ptr           ! pointer to the patch_type
  end type patch_ptr_type

  type, public :: patch_list_type
    PetscInt :: num_patch_objects
    type(patch_type), pointer :: first
    type(patch_type), pointer :: last
    type(patch_ptr_type), pointer :: array(:)
  end type patch_list_type

  PetscInt, parameter, public :: INT_VAR = 0
  PetscInt, parameter, public :: REAL_VAR = 1

  interface PatchGetVariable
    module procedure PatchGetVariable1
    module procedure PatchGetVariable2
  end interface

  interface PatchUnsupportedVariable
    module procedure PatchUnsupportedVariable1
    module procedure PatchUnsupportedVariable2
    module procedure PatchUnsupportedVariable3
    module procedure PatchUnsupportedVariable4
  end interface

  public :: PatchCreate, PatchDestroy, PatchCreateList, PatchDestroyList, &
            PatchAddToList, PatchConvertListToArray, PatchProcessCouplers, &
            PatchUpdateAllCouplerAuxVars, PatchInitAllCouplerAuxVars, &
            PatchLocalizeRegions, PatchUpdateUniformVelocity, &
            PatchGetVariable, PatchGetVariableValueAtCell, &
            PatchSetVariable, PatchCouplerInputRecord, &
            PatchInitConstraints, &
            PatchCountCells, PatchGetIvarsFromKeyword, &
            PatchGetVarNameFromKeyword, &
            PatchCalculateCFL1Timestep, &
            PatchGetCellCenteredVelocities, &
            PatchGetWaterMassInRegion, &
            PatchGetCompMassInRegionAssign, &
            PatchUpdateCouplerSaturation, &
            PatchSetupUpwindDirection

contains

! ************************************************************************** !

function PatchCreate()
  !
  ! Allocates and initializes a new Patch object
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  implicit none

  type(patch_type), pointer :: PatchCreate

  type(patch_type), pointer :: patch

  allocate(patch)

  patch%id = 0
  nullify(patch%imat)
  nullify(patch%imat_internal_to_external)
  nullify(patch%cc_id)
  nullify(patch%cct_id)
  nullify(patch%mtf_id)
  nullify(patch%internal_velocities)
  nullify(patch%boundary_velocities)
  nullify(patch%internal_tran_coefs)
  nullify(patch%boundary_tran_coefs)
  nullify(patch%internal_flow_fluxes)
  nullify(patch%boundary_flow_fluxes)
  nullify(patch%internal_tran_fluxes)
  nullify(patch%boundary_tran_fluxes)
  nullify(patch%ss_flow_fluxes)
  nullify(patch%ss_tran_fluxes)
  nullify(patch%ss_flow_vol_fluxes)

  nullify(patch%flow_upwind_direction)
  nullify(patch%flow_upwind_direction_bc)

  nullify(patch%grid)

  allocate(patch%region_list)
  call RegionInitList(patch%region_list)

  allocate(patch%boundary_condition_list)
  call CouplerInitList(patch%boundary_condition_list)
  allocate(patch%initial_condition_list)
  call CouplerInitList(patch%initial_condition_list)
  allocate(patch%source_sink_list)
  call CouplerInitList(patch%source_sink_list)

  nullify(patch%material_properties)
  nullify(patch%material_property_array)
  nullify(patch%saturation_functions)
  nullify(patch%saturation_function_array)
  nullify(patch%characteristic_curves)
  nullify(patch%characteristic_curves_array)
  nullify(patch%characteristic_curves_thermal)
  nullify(patch%char_curves_thermal_array)
  nullify(patch%material_transform)
  nullify(patch%material_transform_array)

  allocate(patch%observation_list)
  call ObservationInitList(patch%observation_list)
  allocate(patch%integral_flux_list)
  call IntegralFluxInitList(patch%integral_flux_list)
  allocate(patch%strata_list)
  call StrataInitList(patch%strata_list)

  call AuxInit(patch%aux)

  nullify(patch%field)
  nullify(patch%datasets)
  nullify(patch%reaction_base)
  nullify(patch%reaction)
  nullify(patch%reaction_nw)

  nullify(patch%next)

  PatchCreate => patch

end function PatchCreate

! ************************************************************************** !

function PatchCreateList()
  !
  ! PatchListCreate: Creates a patch list
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  implicit none

  type(patch_list_type), pointer :: PatchCreateList

  type(patch_list_type), pointer :: patch_list

  allocate(patch_list)
  nullify(patch_list%first)
  nullify(patch_list%last)
  nullify(patch_list%array)
  patch_list%num_patch_objects = 0

  PatchCreateList => patch_list

end function PatchCreateList

! ************************************************************************** !

subroutine PatchAddToList(new_patch,patch_list)
  !
  ! Adds a new patch to list
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  implicit none

  type(patch_type), pointer :: new_patch
  type(patch_list_type) :: patch_list

  if (associated(new_patch)) then
     patch_list%num_patch_objects = patch_list%num_patch_objects + 1
     new_patch%id = patch_list%num_patch_objects
     if (.not.associated(patch_list%first)) patch_list%first => new_patch
     if (associated(patch_list%last)) patch_list%last%next => new_patch
     patch_list%last => new_patch
  end if
end subroutine PatchAddToList

! ************************************************************************** !

subroutine PatchConvertListToArray(patch_list)
  !
  ! Creates an array of pointers to the
  ! patchs in the patch list
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  implicit none

  type(patch_list_type) :: patch_list

  PetscInt :: count
  type(patch_type), pointer :: cur_patch


  allocate(patch_list%array(patch_list%num_patch_objects))

  cur_patch => patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    patch_list%array(cur_patch%id)%ptr => cur_patch
    cur_patch => cur_patch%next
  enddo

end subroutine PatchConvertListToArray

! ************************************************************************** !

subroutine PatchLocalizeRegions(patch,regions,option)
  !
  ! Localizes regions within each patch
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  use Option_module
  use Output_Aux_module
  use Region_module

  implicit none

  type(patch_type) :: patch
  type(region_list_type) :: regions
  type(option_type) :: option

  type(region_type), pointer :: cur_region
  type(region_type), pointer :: patch_region

  cur_region => regions%first
  do
    if (.not.associated(cur_region)) exit
    patch_region => RegionCreate(cur_region)
    call RegionAddToList(patch_region,patch%region_list)
    cur_region => cur_region%next
  enddo

  !geh: All grids must be localized through GridLocalizeRegions.  Patch
  !     should not differentiate between structured/unstructured, etc.
  call GridLocalizeRegions(patch%grid,patch%region_list,option)

end subroutine PatchLocalizeRegions

! ************************************************************************** !

subroutine PatchProcessCouplers(patch,flow_conditions,transport_conditions, &
                                geophysics_conditions, option)

  !
  ! Assigns conditions and regions to couplers
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  use Option_module
  use Material_module
  use Condition_module
  use Transport_Constraint_module
  use Connection_module

  implicit none

  type(patch_type) :: patch
  type(condition_list_type) :: flow_conditions
  type(tran_condition_list_type) :: transport_conditions
  type(geop_condition_list_type) :: geophysics_conditions
  type(option_type) :: option

  type(coupler_type), pointer :: coupler
  type(coupler_list_type), pointer :: coupler_list
  type(strata_type), pointer :: strata
  type(observation_type), pointer :: observation, next_observation
  type(integral_flux_type), pointer :: integral_flux

  PetscInt :: temp_int, isub
  PetscInt :: nphase
  PetscErrorCode :: ierr

  ! boundary conditions
  coupler => patch%boundary_condition_list%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           patch%region_list)
    if (.not.associated(coupler%region)) then
      option%io_buffer = 'Region "' // trim(coupler%region_name) // &
                 '" in boundary condition "' // &
                 trim(coupler%name) // &
                 '" not found in region list'
      call PrintErrMsg(option)
    endif
    if (associated(patch%grid%structured_grid)) then
      if (coupler%region%num_cells > 0 .and. &
          (coupler%region%iface == 0 .and. &
           .not.associated(coupler%region%faces))) then
        option%io_buffer = 'Region "' // trim(coupler%region_name) // &
                 '", which is tied to a boundary condition, has not &
                 &been assigned a face in the structured grid. '
        call PrintErrMsg(option)
      endif
    endif
    ! pointer to flow condition
    if (option%nflowdof > 0) then
      if (len_trim(coupler%flow_condition_name) > 0) then
        coupler%flow_condition => &
          FlowConditionGetPtrFromList(coupler%flow_condition_name, &
                                      flow_conditions)
        if (.not.associated(coupler%flow_condition)) then
          option%io_buffer = 'Flow condition "' // &
                   trim(coupler%flow_condition_name) // &
                   '" in boundary condition "' // &
                   trim(coupler%name) // &
                   '" not found in flow condition list'
          call PrintErrMsg(option)
        endif
      else
        option%io_buffer = 'A FLOW_CONDITION must be specified in &
                           &BOUNDARY_CONDITION: ' // trim(coupler%name) // '.'
        call PrintErrMsg(option)
      endif
    endif
    ! pointer to transport condition
    if (option%ntrandof > 0) then
      if (len_trim(coupler%tran_condition_name) > 0) then
        coupler%tran_condition => &
          TranConditionGetPtrFromList(coupler%tran_condition_name, &
                                      transport_conditions)
        if (.not.associated(coupler%tran_condition)) then
           option%io_buffer = 'Transport condition "' // &
                   trim(coupler%tran_condition_name) // &
                   '" in boundary condition "' // &
                   trim(coupler%name) // &
                   '" not found in transport condition list'
          call PrintErrMsg(option)
        endif
      else
        option%io_buffer = 'A TRANSPORT_CONDITION must be specified in &
                           &BOUNDARY_CONDITION: ' // trim(coupler%name) // '.'
        call PrintErrMsg(option)
      endif
    endif
    ! pointer to geophysics condition
    if (option%ngeopdof > 0) then
      if (len_trim(coupler%geop_condition_name) > 0) then
        coupler%geop_condition => &
          GeopConditionGetPtrFromList(coupler%geop_condition_name, &
                                      geophysics_conditions)
        if (.not.associated(coupler%geop_condition)) then
           option%io_buffer = 'Geophysics condition "' // &
                   trim(coupler%geop_condition_name) // &
                   '" in boundary condition "' // &
                   trim(coupler%name) // &
                   '" not found in geophysics condition list'
          call PrintErrMsg(option)
        endif
      else
        option%io_buffer = 'A GEOPHYSICS_CONDITION must be specified in &
                           &BOUNDARY_CONDITION: ' // trim(coupler%name) // '.'
        call PrintErrMsg(option)
      endif
    endif
    coupler => coupler%next
  enddo


  ! initial conditions
  coupler => patch%initial_condition_list%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           patch%region_list)
    if (.not.associated(coupler%region)) then
      option%io_buffer = 'Region "' // trim(coupler%region_name) // &
                 '" in initial condition "' // &
                 trim(coupler%name) // &
                 '" not found in region list'
      call PrintErrMsg(option)
    endif
    ! pointer to flow condition
    if (option%nflowdof > 0) then
      if (len_trim(coupler%flow_condition_name) > 0) then
        coupler%flow_condition => &
          FlowConditionGetPtrFromList(coupler%flow_condition_name, &
                                      flow_conditions)
        if (.not.associated(coupler%flow_condition)) then
          option%io_buffer = 'Flow condition "' // &
                   trim(coupler%flow_condition_name) // &
                   '" in initial condition "' // &
                   trim(coupler%name) // &
                   '" not found in flow condition list'
          call PrintErrMsg(option)
        endif
      else
        option%io_buffer = 'A FLOW_CONDITION must be specified in ' // &
                           'INITIAL_CONDITION: ' // trim(coupler%name) // '.'
        call PrintErrMsg(option)
      endif
    endif
    ! pointer to transport condition
    if (option%ntrandof > 0) then
      if (len_trim(coupler%tran_condition_name) > 0) then
        coupler%tran_condition => &
          TranConditionGetPtrFromList(coupler%tran_condition_name, &
                                      transport_conditions)
        if (.not.associated(coupler%tran_condition)) then
          option%io_buffer = 'Transport condition "' // &
                   trim(coupler%tran_condition_name) // &
                   '" in initial condition "' // &
                   trim(coupler%name) // &
                   '" not found in transport condition list'
          call PrintErrMsg(option)
        endif
      else
        option%io_buffer = 'A TRANSPORT_CONDITION must be specified in &
                           &INITIAL_CONDITION: ' // trim(coupler%name) // '.'
        call PrintErrMsg(option)
      endif
    endif
    coupler => coupler%next
  enddo

  ! source/sinks
  coupler => patch%source_sink_list%first
  do
    if (.not.associated(coupler)) exit
    ! pointer to region
    coupler%region => RegionGetPtrFromList(coupler%region_name, &
                                           patch%region_list)
    if (.not.associated(coupler%region)) then
      option%io_buffer = 'Region "' // trim(coupler%region_name) // &
                 '" in source/sink "' // &
                 trim(coupler%name) // &
                 '" not found in region list'
      call PrintErrMsg(option)
    endif

    ! pointer to flow condition
    if (option%nflowdof > 0) then
      if (len_trim(coupler%flow_condition_name) > 0) then
        coupler%flow_condition => &
          FlowConditionGetPtrFromList(coupler%flow_condition_name, &
                                      flow_conditions)
        if (.not.associated(coupler%flow_condition)) then
          option%io_buffer = 'Flow condition "' // &
                   trim(coupler%flow_condition_name) // &
                   '" in source/sink "' // &
                   trim(coupler%name) // &
                   '" not found in flow condition list'
          call PrintErrMsg(option)
        endif
        ! check to ensure that a rate subcondition exists
        if (.not.associated(coupler%flow_condition%rate) .and. &
              .not.associated(coupler%flow_condition%well)) then
          temp_int = 0
          if (associated(coupler%flow_condition%general)) then
            if (associated(coupler%flow_condition%general%rate)) then
              temp_int = 1
            endif
          endif
          if (associated(coupler%flow_condition%hydrate)) then
            if (associated(coupler%flow_condition%hydrate%rate)) then
              temp_int = 1
            endif
          endif
          if (temp_int == 0) then
            option%io_buffer = 'FLOW_CONDITIONs associated with &
              &SOURCE_SINKs must have a RATE or WELL expression within them.'
            call PrintErrMsg(option)
          endif
        endif
      else
        option%io_buffer = 'A FLOW_CONDITION must be specified in &
                           &SOURCE_SINK: ' // trim(coupler%name) // '.'
        call PrintErrMsg(option)
      endif
    endif
    ! pointer to transport condition
    if (option%ntrandof > 0) then
      if (len_trim(coupler%tran_condition_name) > 0) then
        coupler%tran_condition => &
          TranConditionGetPtrFromList(coupler%tran_condition_name, &
                                      transport_conditions)
        if (.not.associated(coupler%tran_condition)) then
          option%io_buffer = 'Transport condition "' // &
                   trim(coupler%flow_condition_name) // &
                   '" in source/sink "' // &
                   trim(coupler%name) // &
                   '" not found in transport condition list'
          call PrintErrMsg(option)
        endif
      else
        option%io_buffer = 'A TRANSPORT_CONDITION must be specified in &
                           &SOURCE_SINK: ' // trim(coupler%name) // '.'
        call PrintErrMsg(option)
      endif
    endif
    coupler => coupler%next
  enddo

!----------------------------
! AUX

  ! strata
  ! connect pointers from strata to regions
  strata => patch%strata_list%first
  do
    if (.not.associated(strata)) exit
    ! pointer to region
    if (len_trim(strata%region_name) > 0) then
      strata%region => RegionGetPtrFromList(strata%region_name, &
                                                  patch%region_list)
      if (.not.associated(strata%region)) then
        option%io_buffer = 'Region "' // trim(strata%region_name) // &
                 '" in strata not found in region list'
        call PrintErrMsg(option)
      endif
      if (strata%active) then
        ! pointer to material
        strata%material_property => &
          MaterialPropGetPtrFromArray(strata%material_property_name, &
                                      patch%material_property_array)
        if (.not.associated(strata%material_property)) then
          option%io_buffer = 'Material "' // &
                            trim(strata%material_property_name) // &
                            '" not found in material list'
          call PrintErrMsg(option)
        endif
      endif
    else
      nullify(strata%region)
      nullify(strata%material_property)
    endif
    strata => strata%next
  enddo

  ! connectivity between initial conditions, boundary conditions,
  ! srcs/sinks, etc and grid
  call CouplerListComputeConnections(patch%grid,option, &
                                     patch%initial_condition_list)
  call CouplerListComputeConnections(patch%grid,option, &
                                     patch%boundary_condition_list)
  call CouplerListComputeConnections(patch%grid,option, &
                                     patch%source_sink_list)

  ! linkage of observation to regions and couplers must take place after
  ! connection list have been created.
  ! observation
  observation => patch%observation_list%first
  do
    if (.not.associated(observation)) exit
    next_observation => observation%next
    select case(observation%itype)
      case(OBSERVATION_SCALAR,OBSERVATION_AGGREGATE)
        ! pointer to region
        observation%region => RegionGetPtrFromList(observation%linkage_name, &
                                                    patch%region_list)
        if (.not.associated(observation%region)) then
          option%io_buffer = 'Region "' // &
                   trim(observation%linkage_name) // &
                 '" in observation point "' // &
                 trim(observation%name) // &
                 '" not found in region list'
          call PrintErrMsg(option)
        endif
        call MPI_Allreduce(observation%region%num_cells,temp_int, &
                           ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM, &
                           option%mycomm,ierr)
        if (temp_int == 0) then
          option%io_buffer = 'Region "' // trim(observation%region%name) // &
            '" is used in an observation point but lies outside the &
            &model domain.'
          call PrintErrMsg(option)
        endif
        if (observation%region%num_cells == 0 .and. observation%itype == &
            OBSERVATION_SCALAR) then
          ! remove the observation object
          call ObservationRemoveFromList(observation,patch%observation_list)
        endif
      case(OBSERVATION_FLUX)
        coupler => CouplerGetPtrFromList(observation%linkage_name, &
                                         patch%boundary_condition_list,option)
        if (associated(coupler)) then
          observation%connection_set => coupler%connection_set
        else
          option%io_buffer = 'Boundary Condition "' // &
                   trim(observation%linkage_name) // &
                   '" not found in Boundary Condition list'
          call PrintErrMsg(option)
        endif
        if (observation%connection_set%num_connections == 0) then
          ! cannot remove from list, since there must be a global reduction
          ! across all procs
          ! therefore, just nullify connection set
          nullify(observation%connection_set)
        endif
    end select
    observation => next_observation
  enddo

  ! linkage of observation to regions and couplers must take place after
  ! connection list have been created.
  ! observation
  integral_flux => patch%integral_flux_list%first
  do
    if (.not.associated(integral_flux)) exit
    call PatchGetIntegralFluxConnections(patch,integral_flux,option)
    call IntegralFluxSizeStorage(integral_flux,option)
    integral_flux => integral_flux%next
    option%flow%store_fluxes = PETSC_TRUE
    option%transport%store_fluxes = PETSC_TRUE
  enddo

  temp_int = ConnectionGetNumberInList(patch%grid%internal_connection_set_list)
  temp_int = max(temp_int,1)
  nphase = max(option%nphase,option%transport%nphase)

  ! all simulations
  allocate(patch%internal_velocities(nphase,temp_int))
  patch%internal_velocities = 0.d0

  ! flow
  if (option%nflowdof > 0) then
    if (option%flow%store_fluxes) then
      allocate(patch%internal_flow_fluxes(option%nflowdof,temp_int))
      patch%internal_flow_fluxes = 0.d0
    endif
  endif

  ! transport
  if (option%ntrandof > 0) then
    allocate(patch%internal_tran_coefs(option%ntrandof,nphase,temp_int))
    patch%internal_tran_coefs = 0.d0
    if (option%transport%store_fluxes) then
      allocate(patch%internal_tran_fluxes(option%ntrandof,temp_int))
      patch%internal_tran_fluxes = 0.d0
    endif
  endif

  temp_int = CouplerGetNumConnectionsInList(patch%boundary_condition_list)

  if (temp_int > 0) then
    ! all simulations
    allocate(patch%boundary_velocities(nphase,temp_int))
    patch%boundary_velocities = 0.d0
    ! flow
    if (option%nflowdof > 0) then
      if (option%flow%store_fluxes) then
        allocate(patch%boundary_flow_fluxes(option%nflowdof,temp_int))
        patch%boundary_flow_fluxes = 0.d0
      endif
    endif
    ! transport
    if (option%ntrandof > 0) then
      allocate(patch%boundary_tran_coefs(option%ntrandof,nphase, &
                                         temp_int))
      patch%boundary_tran_coefs = 0.d0
      if (option%transport%store_fluxes) then
        allocate(patch%boundary_tran_fluxes(option%ntrandof,temp_int))
        patch%boundary_tran_fluxes = 0.d0
      endif
    endif
  endif

  temp_int = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (temp_int > 0) then
    ! flow
    if (option%nflowdof > 0) then
      allocate(patch%ss_flow_fluxes(option%nflowdof,temp_int))
      patch%ss_flow_fluxes = 0.d0
    endif
    ! transport
    if (option%ntrandof > 0) then
      allocate(patch%ss_tran_fluxes(option%ntrandof,temp_int))
      patch%ss_tran_fluxes = 0.d0
      ! only needed by transport
      allocate(patch%ss_flow_vol_fluxes(nphase,temp_int))
      patch%ss_flow_vol_fluxes = 0.d0
    endif
  endif

end subroutine PatchProcessCouplers

! ************************************************************************** !

subroutine PatchInitAllCouplerAuxVars(patch,option)
  !
  ! Initializes coupler auxillary variables
  ! within list
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  use Option_module
  use Reaction_Aux_module

  implicit none

  type(patch_type), pointer :: patch
  type(option_type) :: option

  PetscBool :: force_update_flag = PETSC_TRUE

  call PatchInitCouplerAuxVars(patch%initial_condition_list,patch, &
                               option)
  call PatchInitCouplerAuxVars(patch%boundary_condition_list,patch, &
                               option)
  call PatchInitCouplerAuxVars(patch%source_sink_list,patch, &
                               option)

  !geh: This should not be included in PatchUpdateAllCouplerAuxVars
  ! as it will result in excessive updates to initial conditions
  ! that are not necessary after the simulation has started time stepping.
  call PatchUpdateCouplerAuxVars(patch,patch%initial_condition_list, &
                                 force_update_flag,option)
  call PatchUpdateAllCouplerAuxVars(patch,force_update_flag,option)

end subroutine PatchInitAllCouplerAuxVars

! ************************************************************************** !

subroutine PatchInitCouplerAuxVars(coupler_list,patch,option)
  !
  ! Initializes coupler auxillary variables
  ! within list
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  use Option_module
  use Connection_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use NW_Transport_Aux_module
  use Global_Aux_module
  use ZFlow_Aux_module
  use Condition_module
  use Transport_Constraint_Base_module
  use Transport_Constraint_NWT_module
  use Transport_Constraint_RT_module
  use Transport_Constraint_module
  use General_Aux_module
  use WIPP_Flow_Aux_module

  implicit none

  type(coupler_list_type), pointer :: coupler_list
  type(patch_type), pointer :: patch
  type(option_type) :: option

  PetscInt :: num_connections
  PetscBool :: force_update_flag
  PetscBool :: iflag

  type(coupler_type), pointer :: coupler
  class(tran_constraint_coupler_base_type), pointer :: cur_constraint_coupler
  PetscInt :: idof
  PetscInt :: ndof
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: temp_int

  if (.not.associated(coupler_list)) return

  coupler => coupler_list%first
  do
    if (.not.associated(coupler)) exit

    if (associated(coupler%connection_set)) then
      num_connections = coupler%connection_set%num_connections

      ! FLOW
      if (associated(coupler%flow_condition)) then
        ! determine whether flow_condition is transient
        coupler%flow_condition%is_transient = &
          FlowConditionIsTransient(coupler%flow_condition)
        if (coupler%itype == INITIAL_COUPLER_TYPE .or. &
            coupler%itype == BOUNDARY_COUPLER_TYPE) then

          if (associated(coupler%flow_condition%pressure) .or. &
              associated(coupler%flow_condition%concentration) .or. &
              associated(coupler%flow_condition%saturation) .or. &
              associated(coupler%flow_condition%temperature) .or. &
              associated(coupler%flow_condition%general) .or. &
              associated(coupler%flow_condition%hydrate)) then

            ! allocate arrays that match the number of connections
            select case(option%iflowmode)

              case(RICHARDS_MODE,RICHARDS_TS_MODE)
                temp_int = 1
                if (associated(coupler%flow_condition%pressure)) then
                  select case(coupler%flow_condition%pressure%itype)
                    case(HYDROSTATIC_CONDUCTANCE_BC, &
                         DIRICHLET_CONDUCTANCE_BC, &
                         HET_HYDROSTATIC_CONDUCTANCE_BC)
                      temp_int = temp_int + 1
                  end select
                endif
                allocate(coupler%flow_aux_real_var(temp_int,num_connections))
                allocate(coupler%flow_aux_int_var(1,num_connections))
                coupler%flow_aux_real_var = 0.d0
                coupler%flow_aux_int_var = 0

              case(ZFLOW_MODE)
                ndof = option%nflowdof
                temp_int = 0
                iflag = PETSC_FALSE
                select case(coupler%flow_condition%pressure%itype)
                  case(HYDROSTATIC_CONDUCTANCE_BC, &
                       DIRICHLET_CONDUCTANCE_BC, &
                       HET_HYDROSTATIC_CONDUCTANCE_BC)
                    temp_int = temp_int + 1
                    iflag = PETSC_TRUE
                end select
                allocate(coupler%flow_bc_type(ndof))
                allocate(coupler%flow_aux_real_var(ndof+temp_int, &
                                                   num_connections))
                !geh: don't need this
                !allocate(coupler%flow_aux_int_var(1,num_connections))
                coupler%flow_bc_type = 0
                coupler%flow_aux_real_var = 0.d0
                !coupler%flow_aux_int_var = 0
                coupler%flow_aux_mapping => ZFlowAuxMapConditionIndices(iflag)

              case(PNF_MODE)
                allocate(coupler%flow_bc_type(1))
                allocate(coupler%flow_aux_real_var(1,num_connections))
                allocate(coupler%flow_aux_int_var(1,num_connections))
                coupler%flow_bc_type = 0
                coupler%flow_aux_real_var = 0.d0
                coupler%flow_aux_int_var = 0

              case(TH_MODE,TH_TS_MODE)
                temp_int = 2
                select case(coupler%flow_condition%pressure%itype)
                  case(HYDROSTATIC_CONDUCTANCE_BC, &
                       DIRICHLET_CONDUCTANCE_BC, &
                       HET_HYDROSTATIC_CONDUCTANCE_BC)
                    temp_int = temp_int + 1
                end select
                allocate(coupler%flow_aux_real_var(temp_int,num_connections))
                allocate(coupler%flow_aux_int_var(1,num_connections))
                coupler%flow_aux_real_var = 0.d0
                coupler%flow_aux_int_var = 0

              case(MPH_MODE)
                allocate(coupler%flow_aux_real_var(option%nflowdof, &
                                                   num_connections))
                allocate(coupler%flow_aux_int_var(1,num_connections))
                coupler%flow_aux_real_var = 0.d0
                coupler%flow_aux_int_var = 0

              case(G_MODE)
                allocate(coupler%flow_aux_mapping(GENERAL_MAX_INDEX))
                allocate(coupler%flow_bc_type(THREE_INTEGER))
                allocate(coupler%flow_aux_real_var(FIVE_INTEGER, &
                                                   num_connections))
                allocate(coupler%flow_aux_int_var(ONE_INTEGER,num_connections))
                coupler%flow_aux_mapping = 0
                coupler%flow_bc_type = 0
                coupler%flow_aux_real_var = 0.d0
                coupler%flow_aux_int_var = 0

              case(H_MODE)
                allocate(coupler%flow_aux_mapping(HYDRATE_MAX_INDEX))
                !MAN: Need to fix these
                allocate(coupler%flow_bc_type(THREE_INTEGER))
                allocate(coupler%flow_aux_real_var(FIVE_INTEGER, &
                                                   num_connections))
                allocate(coupler%flow_aux_int_var(ONE_INTEGER,num_connections))
                coupler%flow_aux_mapping = 0
                coupler%flow_bc_type = 0
                coupler%flow_aux_real_var = 0.d0
                coupler%flow_aux_int_var = 0

              case(WF_MODE)
                allocate(coupler%flow_aux_mapping(WIPPFLO_MAX_INDEX))
                allocate(coupler%flow_bc_type(THREE_INTEGER))
                allocate(coupler%flow_aux_real_var(TWO_INTEGER, &
                                                   num_connections))
                allocate(coupler%flow_aux_int_var(ONE_INTEGER,num_connections))
                coupler%flow_aux_mapping = 0
                coupler%flow_bc_type = 0
                coupler%flow_aux_real_var = 0.d0
                coupler%flow_aux_int_var = 0

              case default
                option%io_buffer = 'Failed allocation for flow condition "' // &
                  trim(coupler%flow_condition%name)
                call PrintErrMsg(option)
            end select

          else if (associated(coupler%flow_condition%rate)) then
            option%io_buffer = 'Flow condition "' // &
              trim(coupler%flow_condition%name) // '" can only be used in a &
              &SOURCE_SINK since a rate is prescribed.'
            call PrintErrMsg(option)
          endif ! associated(coupler%flow_condition%pressure)

        else if (coupler%itype == SRC_SINK_COUPLER_TYPE) then

          if (associated(coupler%flow_condition%rate)) then

            if (option%iflowmode == ZFLOW_MODE) then
              ndof = option%nflowdof
              iflag = PETSC_FALSE
              temp_int = 0
              select case(coupler%flow_condition%rate%itype)
                case(SCALED_VOLUMETRIC_RATE_SS)
                  temp_int = 1
                  iflag = PETSC_TRUE
                case(SCALED_MASS_RATE_SS,MASS_RATE_SS, &
                     HET_VOL_RATE_SS,HET_MASS_RATE_SS)
                  option%io_buffer = 'Mass rate source/sinks not &
                    &supported in ZFLOW mode.'
                  call PrintErrMsg(option)
              end select
              allocate(coupler%flow_bc_type(ndof))
              allocate(coupler%flow_aux_real_var(ndof+temp_int,num_connections))
              coupler%flow_bc_type = 0
              coupler%flow_aux_real_var = 0.d0
              coupler%flow_aux_mapping => &
                ZFlowAuxMapConditionIndices(iflag)
            else
              select case(coupler%flow_condition%rate%itype)
                case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS, &
                     VOLUMETRIC_RATE_SS,MASS_RATE_SS, &
                     HET_VOL_RATE_SS,HET_MASS_RATE_SS)
                  select case(option%iflowmode)
                    case(RICHARDS_MODE,RICHARDS_TS_MODE,PNF_MODE)
                      allocate(coupler%flow_aux_real_var(1,num_connections))
                      coupler%flow_aux_real_var = 0.d0
                    case(TH_MODE,TH_TS_MODE)
                      allocate(coupler%flow_aux_real_var(option%nflowdof, &
                                                         num_connections))
                      coupler%flow_aux_real_var = 0.d0
                    case(MPH_MODE)
                      ! do nothing
                    case default
                      string = GetSubConditionName(coupler%flow_condition%&
                                                   rate%itype)
                      option%io_buffer='Source/Sink of rate%itype = "' // &
                        trim(adjustl(string)) // &
                        '", not implemented in this mode.'
                      call PrintErrMsg(option)
                  end select
                case default
                  string = GetSubConditionName(coupler%flow_condition% &
                                               rate%itype)
                  option%io_buffer = &
                    FlowConditionUnknownItype(coupler%flow_condition,'rate', &
                                              string)
                  call PrintErrMsg(option)
              end select
            endif
          ! handles source/sinks in general mode
          else if (associated(coupler%flow_condition%general)) then
            if (associated(coupler%flow_condition%general%rate)) then
              select case(coupler%flow_condition%general%rate%itype)
                case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS)
                  allocate(coupler%flow_aux_real_var(1,num_connections))
                  coupler%flow_aux_real_var = 0.d0
              end select
            endif
          endif
        endif ! coupler%itype == SRC_SINK_COUPLER_TYPE
      endif ! associated(coupler%flow_condition)
    endif ! associated(coupler%connection_set)

    ! TRANSPORT
    if (associated(coupler%tran_condition)) then
      cur_constraint_coupler => &
                          coupler%tran_condition%constraint_coupler_list
      do
        if (.not.associated(cur_constraint_coupler)) exit
        ! Setting option%iflag = 0 ensures that the "mass_balance" array
        ! is not allocated.
        option%iflag = 0
        ! Only allocate the XXX_auxvar objects if they have not been allocated.
        ! Since coupler%tran_condition is a pointer to a separate list of
        ! tran conditions, the XXX_auxvar object may already be allocated.
        if (.not.associated(cur_constraint_coupler%global_auxvar)) then
          allocate(cur_constraint_coupler%global_auxvar)
          call GlobalAuxVarInit(cur_constraint_coupler%global_auxvar,option)
        endif
        select type(cur_constraint_coupler)
          class is (tran_constraint_coupler_rt_type)
            if (.not.associated(cur_constraint_coupler%rt_auxvar)) then
              allocate(cur_constraint_coupler%rt_auxvar)
              call RTAuxVarInit(cur_constraint_coupler%rt_auxvar, &
                                patch%reaction,option)
            endif
          class is (tran_constraint_coupler_nwt_type)
            if (.not.associated(cur_constraint_coupler%nwt_auxvar)) then
              allocate(cur_constraint_coupler%nwt_auxvar)
              call NWTAuxVarInit(cur_constraint_coupler%nwt_auxvar, &
                                 NWTReactionCast(patch%reaction_base),option)
            endif
        end select
        cur_constraint_coupler => cur_constraint_coupler%next
      enddo
    endif
    coupler => coupler%next
  enddo

end subroutine PatchInitCouplerAuxVars

! ************************************************************************** !

subroutine PatchUpdateAllCouplerAuxVars(patch,force_update_flag,option)
  !
  ! Updates auxiliary variables associated
  ! with couplers in list
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  use Option_module

  implicit none

  type(patch_type) :: patch
  PetscBool :: force_update_flag
  type(option_type) :: option

  PetscInt :: iconn

  !geh: no need to update initial conditions as they only need updating
  !     once as performed in PatchInitCouplerAuxVars()
  call PatchUpdateCouplerAuxVars(patch,patch%boundary_condition_list, &
                                 force_update_flag,option)
  call PatchUpdateCouplerAuxVars(patch,patch%source_sink_list, &
                                 force_update_flag,option)

end subroutine PatchUpdateAllCouplerAuxVars

! ************************************************************************** !

subroutine PatchUpdateCouplerAuxVars(patch,coupler_list,force_update_flag, &
                                     option)
  !
  ! Updates auxiliary variables associated
  ! with couplers in list
  !
  ! Author: Glenn Hammond
  ! Date: 11/26/07
  !
  use Option_module
  use Condition_module
  use Hydrostatic_module


  use Grid_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_HDF5_class

  implicit none

  type(patch_type) :: patch
  type(coupler_list_type), pointer :: coupler_list
  PetscBool :: force_update_flag
  type(option_type) :: option

  type(coupler_type), pointer :: coupler
  type(flow_condition_type), pointer :: flow_condition

  if (.not.associated(coupler_list)) return

  coupler => coupler_list%first

  do
    if (.not.associated(coupler)) exit

    ! FLOW
    if (associated(coupler%flow_aux_real_var)) then

      flow_condition => coupler%flow_condition
      if (force_update_flag .or. flow_condition%is_transient) then
        select case(option%iflowmode)
          case(G_MODE)
            call PatchUpdateCouplerAuxVarsG(patch,coupler,option)
          case(H_MODE)
            call PatchUpdateCouplerAuxVarsH(patch,coupler,option)
          case(WF_MODE)
            call PatchUpdateCouplerAuxVarsWF(patch,coupler,option)
          case(MPH_MODE)
            call PatchUpdateCouplerAuxVarsMPH(patch,coupler,option)
          case(TH_MODE,TH_TS_MODE)
            call PatchUpdateCouplerAuxVarsTH(patch,coupler,option)
          case(RICHARDS_MODE, RICHARDS_TS_MODE)
            call PatchUpdateCouplerAuxVarsRich(patch,coupler,option)
          case(ZFLOW_MODE)
            call PatchUpdateCouplerAuxVarsZFlow(patch,coupler,option)
          case(PNF_MODE)
            call PatchUpdateCouplerAuxVarsPNF(patch,coupler,option)
        end select
      endif
    endif

    ! TRANSPORT
    ! nothing for transport at this point in time
    coupler => coupler%next
  enddo

end subroutine PatchUpdateCouplerAuxVars

! ************************************************************************** !

subroutine PatchUpdateCouplerAuxVarsWF(patch,coupler,option)
  !
  ! Updates flow auxiliary variables associated
  ! with a coupler for WF_MODE
  !
  ! Author: Glenn Hammond
  ! Date: 11/26/13
  !

  use Option_module
  use Condition_module
  use Hydrostatic_module
  use Utility_module, only : DeallocateArray

  use WIPP_Flow_Aux_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_HDF5_class
  use Dataset_Ascii_class
  use Dataset_module
  use General_Aux_module, only : LIQUID_STATE, GAS_STATE, TWO_PHASE_STATE, &
                                 ANY_STATE

  implicit none

  type(patch_type) :: patch
  type(coupler_type), pointer :: coupler
  type(option_type) :: option

  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(flow_general_condition_type), pointer :: general
  PetscBool :: update
  PetscBool :: dof1, dof2
  PetscReal :: temperature, p_sat, p_air, p_gas, p_cap, s_liq, xmol
  PetscReal :: relative_humidity
  PetscReal :: dummy_real
  PetscReal :: x(option%nflowdof)
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr

  PetscInt :: idof, num_connections,sum_connection
  PetscInt :: iconn, local_id, ghosted_id
  ! use to map flow_aux_map to the flow_aux_real_var array
  PetscInt :: real_count

  num_connections = coupler%connection_set%num_connections

  flow_condition => coupler%flow_condition

  general => flow_condition%general
  dof1 = PETSC_FALSE
  dof2 = PETSC_FALSE
  real_count = 0
  select case(flow_condition%iphase)
    case(TWO_PHASE_STATE)
      coupler%flow_aux_int_var(WIPPFLO_STATE_INDEX,1:num_connections) = &
        TWO_PHASE_STATE
      select case(general%liquid_pressure%itype)
        case(DIRICHLET_BC)
          real_count = real_count + 1
          coupler%flow_aux_mapping(WIPPFLO_LIQUID_PRESSURE_INDEX) = real_count
          select type(dataset => general%liquid_pressure%dataset)
            class is(dataset_ascii_type)
              coupler%flow_aux_real_var(real_count,1:num_connections) = &
                dataset%rarray(1)
              coupler%flow_bc_type(WIPPFLO_LIQUID_EQUATION_INDEX) = &
                DIRICHLET_BC
            class is(dataset_gridded_hdf5_type)
              call PatchUpdateCouplerGridDataset(coupler,option, &
                                                 patch%grid,dataset, &
                                                 real_count)
              coupler%flow_bc_type(WIPPFLO_LIQUID_EQUATION_INDEX) = &
                DIRICHLET_BC
            class is(dataset_common_hdf5_type)
              ! skip cell indexed datasets used in initial conditions
            class default
              call PrintMsg(option,'general%liquid_pressure%itype,DIRICHLET_BC')
              call DatasetUnknownClass(dataset,option, &
                                       'PatchUpdateCouplerAuxVarsWF')
          end select
          dof1 = PETSC_TRUE
        case(HYDROSTATIC_BC)
          ! have to increment so that saturation is correct.
          real_count = real_count + 1
          call HydrostaticUpdateCoupler(coupler,option,patch%grid)
          coupler%flow_bc_type(WIPPFLO_LIQUID_EQUATION_INDEX) = HYDROSTATIC_BC
          dof1 = PETSC_TRUE
        case default
          string = &
            GetSubConditionName(general%liquid_pressure%itype)
          option%io_buffer = &
            FlowConditionUnknownItype(coupler%flow_condition, &
              'wipp flow liquid pressure',string)
          call PrintErrMsg(option)
      end select
      ! in two-phase flow, gas saturation is second dof
      select case(general%gas_saturation%itype)
        case(DIRICHLET_BC)
          real_count = real_count + 1
          coupler%flow_aux_mapping(WIPPFLO_GAS_SATURATION_INDEX) = real_count
          select type(dataset => general%gas_saturation%dataset)
            class is(dataset_ascii_type)
              coupler%flow_aux_real_var(real_count,1:num_connections) = &
                dataset%rarray(1)
              coupler%flow_bc_type(WIPPFLO_GAS_EQUATION_INDEX) = DIRICHLET_BC
            class is(dataset_gridded_hdf5_type)
              call PatchUpdateCouplerGridDataset(coupler,option, &
                                                 patch%grid,dataset, &
                                                 real_count)
              coupler%flow_bc_type(WIPPFLO_GAS_EQUATION_INDEX) = DIRICHLET_BC
            class is(dataset_common_hdf5_type)
              ! skip cell indexed datasets used in initial conditions
            class default
              call PrintMsg(option,'general%gas_saturation%itype,DIRICHLET_BC')
              call DatasetUnknownClass(dataset,option, &
                                       'PatchUpdateCouplerAuxVarsWF')
          end select
          dof2 = PETSC_TRUE
        case default
          string = &
            GetSubConditionName(general%gas_saturation%itype)
          option%io_buffer = &
            FlowConditionUnknownItype(coupler%flow_condition, &
              'wipp flow gas saturation',string)
          call PrintErrMsg(option)
      end select
    case(LIQUID_STATE)
      option%io_buffer = 'LIQUID State not support for WIPP Flow mode.'
      call PrintErrMsg(option)
    case(GAS_STATE)
      option%io_buffer = 'GAS State not support for WIPP Flow mode.'
      call PrintErrMsg(option)
    case(ANY_STATE)
      if (associated(coupler%flow_aux_int_var)) then ! not used with rate
        coupler%flow_aux_int_var(WIPPFLO_STATE_INDEX,1:num_connections) = &
          ANY_STATE
      endif
  end select

  if (associated(general%liquid_flux)) then
    coupler%flow_bc_type(WIPPFLO_LIQUID_EQUATION_INDEX) = NEUMANN_BC
    real_count = real_count + 1
    coupler%flow_aux_mapping(WIPPFLO_LIQUID_FLUX_INDEX) = real_count
    select type(selector => general%liquid_flux%dataset)
      class is(dataset_ascii_type)
        coupler%flow_aux_real_var(real_count,1:num_connections) = &
                                           general%liquid_flux%dataset%rarray(1)
        dof1 = PETSC_TRUE
      class is(dataset_gridded_hdf5_type)
        call PatchVerifyDatasetGriddedForFlux(selector,coupler,option)
        call PatchUpdateCouplerGridDataset(coupler,option,patch%grid,selector, &
                                           real_count)
        dof1 = PETSC_TRUE
      class default
        call PrintMsg(option,'general%liquid_flux%dataset')
        call DatasetUnknownClass(selector,option, &
                                 'PatchUpdateCouplerAuxVarsWF')
    end select
  endif
  if (associated(general%gas_flux)) then
    coupler%flow_bc_type(WIPPFLO_GAS_EQUATION_INDEX) = NEUMANN_BC
    real_count = real_count + 1
    coupler%flow_aux_mapping(WIPPFLO_GAS_FLUX_INDEX) = real_count
    select type(selector => general%gas_flux%dataset)
      class is(dataset_ascii_type)
        coupler%flow_aux_real_var(real_count,1:num_connections) = &
                                              general%gas_flux%dataset%rarray(1)
        dof2 = PETSC_TRUE
      class is(dataset_gridded_hdf5_type)
        call PatchVerifyDatasetGriddedForFlux(selector,coupler,option)
        call PatchUpdateCouplerGridDataset(coupler,option,patch%grid,selector, &
                                           real_count)
        dof2 = PETSC_TRUE
      class default
        call PrintMsg(option,'general%gas_flux%dataset')
        call DatasetUnknownClass(selector,option, &
                                 'PatchUpdateCouplerAuxVarsWF')
    end select
  endif
  if (associated(general%energy_flux)) then
          option%io_buffer = 'Temperature not supported for two-phase'
          call PrintErrMsg(option)
!geh: removed for immiscible
    !coupler%flow_bc_type(WIPPFLO_ENERGY_EQUATION_INDEX) = NEUMANN_BC
    !real_count = real_count + 1
    !coupler%flow_aux_mapping(WIPPFLO_ENERGY_FLUX_INDEX) = real_count
    !select type(selector => general%energy_flux%dataset)
    !  class is(dataset_ascii_type)
    !    coupler%flow_aux_real_var(real_count,1:num_connections) = &
    !      general%energy_flux%dataset%rarray(1)
    !    dof3 = PETSC_TRUE
    !  class is(dataset_gridded_hdf5_type)
    !    call PatchVerifyDatasetGriddedForFlux(selector,coupler,option)
    !    call PatchUpdateCouplerGridDataset(coupler,option,patch%grid,selector, &
    !                                       real_count)
    !    dof3 = PETSC_TRUE
    !  class default
    !    option%io_buffer = 'Unknown dataset class for general%energy_flux.'
    !    call PrintErrMsg(option)
    !end select
  endif

  if (real_count > 2) then
    option%io_buffer = &
      'More than two dofs assigned in PatchUpdateCouplerAuxVarsWF.'
    call PrintErrMsg(option)
  endif

  if (associated(general%rate)) then
    select case(general%rate%itype)
      case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS)
        call PatchScaleSourceSink(patch,coupler,general%rate%isubtype,option)
        dof1 = PETSC_TRUE
        dof2 = PETSC_TRUE
    end select
  endif

  if (real_count == 0) then ! no need for the auxiliary arrays
    call DeallocateArray(coupler%flow_aux_mapping)
    call DeallocateArray(coupler%flow_bc_type)
    call DeallocateArray(coupler%flow_aux_real_var)
    call DeallocateArray(coupler%flow_aux_int_var)
  endif

  !geh: is this really correct, or should it be .or.
  if (.not.dof1 .or. .not.dof2) then
    option%io_buffer = 'Error with general phase boundary condition'
    call PrintErrMsg(option)
  endif

end subroutine PatchUpdateCouplerAuxVarsWF

! ************************************************************************** !

subroutine PatchUpdateCouplerAuxVarsG(patch,coupler,option)
  !
  ! Updates flow auxiliary variables associated
  ! with a coupler for G_MODE
  !
  ! Author: Glenn Hammond
  ! Date: 11/26/13
  !

  use Option_module
  use Condition_module
  use Hydrostatic_module
  use EOS_Water_module
  use Utility_module

  use General_Aux_module
  use Grid_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_HDF5_class
  use Dataset_Ascii_class
  use Dataset_module
  use String_module

  implicit none

  type(patch_type) :: patch
  type(coupler_type), pointer :: coupler
  type(option_type) :: option

  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(flow_general_condition_type), pointer :: general
  PetscBool :: update
  PetscBool :: dof1, dof2, dof3
  PetscReal :: temperature, p_sat, p_cap, s_liq, xmol
  PetscReal :: relative_humidity
  PetscReal :: gas_sat, hyd_sat, air_pressure, gas_pressure, liq_pressure
  PetscReal :: dummy_real
  PetscReal :: x(option%nflowdof)
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr

  PetscInt :: idof, num_connections,sum_connection
  PetscInt :: iconn, local_id, ghosted_id
  PetscInt :: real_count
  PetscInt :: dof_count_local(3)
  PetscInt :: dof_count_global(3)
  PetscReal, parameter :: min_two_phase_gas_pressure = 3.d3

  num_connections = coupler%connection_set%num_connections

  flow_condition => coupler%flow_condition

  general => flow_condition%general
  dof1 = PETSC_FALSE
  dof2 = PETSC_FALSE
  dof3 = PETSC_FALSE
  real_count = 0

  ! mapping of flow_aux_mapping to the flow_aux_real_var array:
  if (associated(coupler%flow_aux_mapping)) then
    ! liquid and gas pressure are set to 1st dof index
    ! liquid flux is set to 1st dof index
    coupler%flow_aux_mapping(GENERAL_GAS_PRESSURE_INDEX) = 1
    coupler%flow_aux_mapping(GENERAL_LIQUID_PRESSURE_INDEX) = 1
    coupler%flow_aux_mapping(GENERAL_LIQUID_FLUX_INDEX) = 1
    ! temperature is set to 2nd dof index
    ! energy flux is set to 2nd dof index
    coupler%flow_aux_mapping(GENERAL_TEMPERATURE_INDEX) = 2
    coupler%flow_aux_mapping(GENERAL_ENERGY_FLUX_INDEX) = 2
    ! air mole fraction, gas sat., and air pressure are set to 3rd dof index
    ! gas flux is set to 3rd dof index
    coupler%flow_aux_mapping(GENERAL_MOLE_FRACTION_INDEX) = 3
    coupler%flow_aux_mapping(GENERAL_GAS_SATURATION_INDEX) = 3
    coupler%flow_aux_mapping(GENERAL_AIR_PRESSURE_INDEX) = 3
    coupler%flow_aux_mapping(GENERAL_GAS_FLUX_INDEX) = 3
    coupler%flow_aux_mapping(GENERAL_GAS_WATER_MOL_FRAC_INDEX) = 3
  endif

  select case(flow_condition%iphase)
    case(MULTI_STATE)
      select type(dataset => general%gas_saturation%dataset)
        class is(dataset_ascii_type)
          gas_sat = general%gas_saturation%dataset%rarray(1)
          if (gas_sat > 0.d0 .and. gas_sat < 1.d0) then
            coupler%flow_aux_int_var(GENERAL_STATE_INDEX,1:num_connections) = &
              TWO_PHASE_STATE
          ! Cannot user gas_sat == 0.d0 or Equal(gas_sat,0.d0) as optimization
          ! in the Intel compiler changes the answer.
          else if (gas_sat < 0.5d0) then
            coupler%flow_aux_int_var(GENERAL_STATE_INDEX,1:num_connections) = &
              LIQUID_STATE
          else
            coupler%flow_aux_int_var(GENERAL_STATE_INDEX,1:num_connections) = &
              GAS_STATE
          endif
        class is(dataset_gridded_hdf5_type)
          ! If the gas pressure dataset is defined, we must ensure that its
          ! minimum pressure is greater than min_two_phase_gas_pressure.
          ! Otherwise, a zero gas pressure may be incorrectly used within
          ! an interpolation of gas pressure for a two phase cell mixed
          ! with a single phase liquid cell
          if (associated(general%gas_pressure)) then
            dummy_real = &
              DatasetGetMinRValue(general%gas_pressure%dataset,option)
            if (dummy_real < min_two_phase_gas_pressure) then
              option%io_buffer = 'Minimum gas pressure exceeded for &
                    &FLOW_CONDITION "' // trim(flow_condition%name) // &
                    '": ' // trim(StringFormatDouble(dummy_real)) // '.'
              call PrintErrMsg(option)
            endif
          endif
          do iconn = 1, num_connections
            call PatchGetCouplerValueFromDataset(coupler,option,patch%grid, &
                                  general%gas_saturation%dataset,iconn,gas_sat)
            if (gas_sat > 0.d0 .and. gas_sat < 1.d0) then
              coupler%flow_aux_int_var(GENERAL_STATE_INDEX,iconn) = &
                TWO_PHASE_STATE
            else if (gas_sat < 0.5d0) then
              coupler%flow_aux_int_var(GENERAL_STATE_INDEX,iconn) = LIQUID_STATE
            else
              coupler%flow_aux_int_var(GENERAL_STATE_INDEX,iconn) = GAS_STATE
            endif
          enddo
        class default
          call PrintMsg(option,'general%gas_saturation%dataset,MULTI_STATE')
          call DatasetUnknownClass(dataset,option, &
                                   'PatchUpdateCouplerAuxVarsG')
      end select
    case(TWO_PHASE_STATE)
      coupler%flow_aux_int_var(GENERAL_STATE_INDEX,1:num_connections) = &
        TWO_PHASE_STATE
        ! no need to loop in the next do loop if its all the same state, which
        ! you know from flow_condition%iphase
    case(LIQUID_STATE)
      coupler%flow_aux_int_var(GENERAL_STATE_INDEX,1:num_connections) = &
        LIQUID_STATE
      if (general%liquid_pressure%itype == HYDROSTATIC_BC) then
        if (general%mole_fraction%itype /= DIRICHLET_BC) then
          option%io_buffer = 'Hydrostatic liquid state pressure BC for &
            &flow condition "' // trim(flow_condition%name) // &
            '" requires a mole fraction BC of type DIRICHLET.'
          call PrintErrMsg(option)
        endif
        if (general%temperature%itype /= DIRICHLET_BC) then
          option%io_buffer = 'Hydrostatic liquid state pressure BC for &
            &flow condition "' // trim(flow_condition%name) // &
            '" requires a temperature BC of type DIRICHLET.'
          call PrintErrMsg(option)
        endif
        call HydrostaticUpdateCoupler(coupler,option,patch%grid)
        do iconn = 1, num_connections
          if (coupler%flow_aux_int_var(GENERAL_STATE_INDEX,iconn) /= &
              LIQUID_STATE) then
            select case(coupler%flow_aux_int_var(GENERAL_STATE_INDEX,iconn))
              case(GAS_STATE)
                string = 'gas state'
              case(TWO_PHASE_STATE)
                string = 'two phase state'
              case(ANY_STATE)
                string = 'any phase state'
            end select
            option%io_buffer = 'A ' // trim(string) // ' cell was found &
              &within a HYDROSTATIC_BC boundary condition for GENERAL mode. &
              &A hydrostatic boundary condition may not be used to set &
              &state variables in the vadose zone for GENERAL mode.'
            call PrintErrMsg(option)
          endif
        enddo
        dof1 = PETSC_TRUE; dof2 = PETSC_TRUE; dof3 = PETSC_TRUE;
      endif
    case(GAS_STATE)
      coupler%flow_aux_int_var(GENERAL_STATE_INDEX,1:num_connections) = &
        GAS_STATE
    case(ANY_STATE)
      if (associated(coupler%flow_aux_int_var)) then ! not used with rate
        coupler%flow_aux_int_var(GENERAL_STATE_INDEX,1:num_connections) = &
          ANY_STATE
      endif
  end select

  ! loop over each connection in the coupler and check its state
  ! set the flow_aux_mapping, flow_aux_real_var, etc on a connection
  ! basis rather than in coupler chunks
  ! this might be slower since we need to loop over all the connections
  ! but it makes the algorithm more general
  if (associated(coupler%flow_aux_int_var)) then
    do iconn = 1, num_connections
      select case(coupler%flow_aux_int_var(GENERAL_STATE_INDEX,iconn))
      ! ---------------------------------------------------------------------- !
        case(TWO_PHASE_STATE)
          ! gas pressure; 1st dof ------------------------ !
          select case(general%gas_pressure%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                     patch%grid,general%gas_pressure%dataset,iconn,gas_pressure)
              coupler%flow_aux_real_var(ONE_INTEGER,iconn) = gas_pressure
              dof1 = PETSC_TRUE
              coupler%flow_bc_type(GENERAL_LIQUID_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(general%gas_pressure%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'GENERAL_MODE two phase state gas pressure ',string)
              call PrintErrMsg(option)
          end select
          ! temperature; 2nd dof ------------------------- !
          select case(general%temperature%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                       patch%grid,general%temperature%dataset,iconn,temperature)
              if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
                coupler%flow_aux_real_var(TWO_INTEGER,iconn) = temperature
              else
                call EOSWaterSaturationPressure(temperature,p_sat,ierr)
                call PatchGetCouplerValueFromDataset(coupler,option, &
                     patch%grid,general%gas_pressure%dataset,iconn,gas_pressure)
                ! should it still be index = 2 here below?
                coupler%flow_aux_real_var(TWO_INTEGER,iconn) = &
                                                              gas_pressure-p_sat
              endif
              dof2 = PETSC_TRUE
              coupler%flow_bc_type(GENERAL_ENERGY_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(general%temperature%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'GENERAL_MODE two phase state temperature ',string)
              call PrintErrMsg(option)
          end select
          ! gas saturation; 3rd dof ---------------------- !
          select case(general%gas_saturation%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                        patch%grid,general%gas_saturation%dataset,iconn,gas_sat)
              coupler%flow_aux_real_var(THREE_INTEGER,iconn) = gas_sat
              dof3 = PETSC_TRUE
              coupler%flow_bc_type(GENERAL_GAS_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(general%gas_saturation%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'GENERAL_MODE two phase state gas saturation ',string)
              call PrintErrMsg(option)
          end select
      ! ---------------------------------------------------------------------- !
        case(LIQUID_STATE)
          if (general%liquid_pressure%itype == HYDROSTATIC_BC) then
  !         option%io_buffer = 'Hydrostatic BC for general phase cannot possibly ' // &
  !           'be set up correctly. - GEH'
  !         call PrintErrMsg(option)
            if (general%mole_fraction%itype /= DIRICHLET_BC) then
              option%io_buffer = 'Hydrostatic liquid state pressure BC for &
                &flow condition "' // trim(flow_condition%name) // &
                '" requires a mole fraction BC of type DIRICHLET.'
              call PrintErrMsg(option)
            endif
            if (general%temperature%itype /= DIRICHLET_BC) then
              option%io_buffer = 'Hydrostatic liquid state pressure BC for &
                &flow condition "' // trim(flow_condition%name) // &
                '" requires a temperature BC of type DIRICHLET.'
              call PrintErrMsg(option)
            endif
            ! ---> see code that just prints error
            coupler%flow_bc_type(1) = HYDROSTATIC_BC
            coupler%flow_bc_type(2:3) = DIRICHLET_BC
          else
          ! liquid pressure; 1st dof --------------------- !
            select case(general%liquid_pressure%itype)
              case(DIRICHLET_BC)
                call PatchGetCouplerValueFromDataset(coupler,option, &
                  patch%grid,general%liquid_pressure%dataset,iconn,liq_pressure)
                coupler%flow_aux_real_var(ONE_INTEGER,iconn) = liq_pressure
                dof1 = PETSC_TRUE
                coupler%flow_bc_type(GENERAL_LIQUID_EQUATION_INDEX) = &
                                                                    DIRICHLET_BC
              case default
                string = GetSubConditionName(general%liquid_pressure%itype)
                option%io_buffer = &
                  FlowConditionUnknownItype(coupler%flow_condition, &
                  'GENERAL_MODE liquid state liquid pressure ',string)
                call PrintErrMsg(option)
            end select
          ! temperature; 2nd dof ------------------------- !
            select case(general%temperature%itype)
              case(DIRICHLET_BC)
                call PatchGetCouplerValueFromDataset(coupler,option, &
                       patch%grid,general%temperature%dataset,iconn,temperature)
                coupler%flow_aux_real_var(TWO_INTEGER,iconn) = temperature
                dof2 = PETSC_TRUE
                coupler%flow_bc_type(GENERAL_ENERGY_EQUATION_INDEX) = &
                                                                    DIRICHLET_BC
              case default
                string = GetSubConditionName(general%temperature%itype)
                option%io_buffer = &
                  FlowConditionUnknownItype(coupler%flow_condition, &
                  'GENERAL_MODE liquid state temperature ',string)
                call PrintErrMsg(option)
            end select
          ! mole fraction; 3rd dof ----------------------- !
            select case(general%mole_fraction%itype)
              case(DIRICHLET_BC)
                call PatchGetCouplerValueFromDataset(coupler,option, &
                            patch%grid,general%mole_fraction%dataset,iconn,xmol)
                if (general_immiscible) then
                  xmol = GENERAL_IMMISCIBLE_VALUE
                endif
                coupler%flow_aux_real_var(THREE_INTEGER,iconn) = xmol
                dof3 = PETSC_TRUE
                coupler%flow_bc_type(GENERAL_GAS_EQUATION_INDEX) = DIRICHLET_BC
              case default
                string = GetSubConditionName(general%mole_fraction%itype)
                option%io_buffer = &
                  FlowConditionUnknownItype(coupler%flow_condition, &
                  'GENERAL_MODE liquid state mole fraction ',string)
                call PrintErrMsg(option)
            end select
          endif
      ! ---------------------------------------------------------------------- !
        case(GAS_STATE)
          gas_pressure = UNINITIALIZED_DOUBLE
          temperature = UNINITIALIZED_DOUBLE
          ! gas pressure; 1st dof ------------------------ !
          select case(general%gas_pressure%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                  patch%grid,general%gas_pressure%dataset,iconn,gas_pressure)
              coupler%flow_aux_real_var(ONE_INTEGER,iconn) = gas_pressure
              dof1 = PETSC_TRUE
              coupler%flow_bc_type(GENERAL_GAS_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(general%gas_pressure%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                'GENERAL_MODE gas state gas pressure',string)
              call PrintErrMsg(option)
          end select
          ! temperature; 2nd dof ------------------------- !
          select case(general%temperature%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                       patch%grid,general%temperature%dataset,iconn,temperature)
              coupler%flow_aux_real_var(TWO_INTEGER,iconn) = &
                temperature
              dof2 = PETSC_TRUE
              coupler%flow_bc_type(GENERAL_ENERGY_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(general%temperature%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                'GENERAL_MODE gas state temperature',string)
              call PrintErrMsg(option)
          end select
          ! air mole fraction; 3rd dof ------------------- !
          if (associated(general%mole_fraction)) then
            select case(general%mole_fraction%itype)
              case(DIRICHLET_BC)
                if (Uninitialized(gas_pressure) .or. &
                    Uninitialized(temperature)) then
                  option%io_buffer = 'GAS_PRESSURE or TEMPERATURE not set &
                    &correctly in flow condition "' // &
                    trim(flow_condition%name) // '".'
                  call PrintErrMsg(option)
                endif
                call PatchGetCouplerValueFromDataset(coupler,option, &
                            patch%grid,general%mole_fraction%dataset,iconn,xmol)
                air_pressure = xmol * gas_pressure
                if (general_immiscible) then
                  air_pressure = gas_pressure - GENERAL_IMMISCIBLE_VALUE
                endif
                call EOSWaterSaturationPressure(temperature,p_sat,ierr)
                if (gas_pressure - air_pressure >= p_sat) then
                  option%io_buffer = 'MOLE_FRACTION set in flow &
                    &condition "' // trim(flow_condition%name) // &
                    '" results in a vapor pressure exceeding the water &
                    &saturation pressure, which indicates that a two-phase &
                    &state with GAS_PRESSURE and GAS_SATURATION should be used.'
                  call PrintErrMsg(option)
                endif
                if (general_gas_air_mass_dof == GENERAL_AIR_PRESSURE_INDEX) then
                  coupler%flow_aux_real_var(THREE_INTEGER,iconn) = air_pressure
                  dof3 = PETSC_TRUE
                  coupler%flow_bc_type(GENERAL_LIQUID_EQUATION_INDEX) = &
                                                                     DIRICHLET_BC
                else
                  coupler%flow_aux_real_var(THREE_INTEGER,iconn) = 1.d0 - xmol
                  dof3 = PETSC_TRUE
                  coupler%flow_bc_type(GENERAL_LIQUID_EQUATION_INDEX) = &
                                                                    DIRICHLET_BC
                endif
              case default
                string = GetSubConditionName(general%mole_fraction%itype)
                option%io_buffer = &
                  FlowConditionUnknownItype(coupler%flow_condition, &
                  'GENERAL_MODE air mole fraction',string)
                call PrintErrMsg(option)
            end select
        ! relative humidity; 3rd dof ------------------- !
          else
            select case(general%relative_humidity%itype)
              case(DIRICHLET_BC)
                if (Uninitialized(gas_pressure) .or. &
                    Uninitialized(temperature)) then
                  option%io_buffer = 'GAS_PRESSURE or TEMPERATURE not set &
                    &correctly in flow condition "' // &
                    trim(flow_condition%name) // '".'
                  call PrintErrMsg(option)
                endif
                call PatchGetCouplerValueFromDataset(coupler,option, &
                  patch%grid,general%relative_humidity%dataset, &
                  iconn,relative_humidity)  ! relative_humidity is in percent
                if (relative_humidity < 0.d0 .or. &
                    relative_humidity > 100.d0) then
                  option%io_buffer = 'RELATIVE_HUMIDITY in flow &
                    &condition "' // trim(flow_condition%name) // '" outside &
                    &bounds of 0-100%.'
                  call PrintErrMsg(option)
                endif
                call EOSWaterSaturationPressure(temperature,p_sat,ierr)
                                  ! convert from % to fraction
                air_pressure = gas_pressure - relative_humidity*1.d-2*p_sat
                if (general_immiscible) then
                  air_pressure = gas_pressure - GENERAL_IMMISCIBLE_VALUE
                endif
                coupler%flow_aux_real_var(THREE_INTEGER,iconn) = air_pressure
                dof3 = PETSC_TRUE
                coupler%flow_bc_type(GENERAL_LIQUID_EQUATION_INDEX) = &
                                                                    DIRICHLET_BC
              case default
                string = GetSubConditionName(general%relative_humidity%itype)
                option%io_buffer = &
                  FlowConditionUnknownItype(coupler%flow_condition, &
                  'GENERAL_MODE relative humidity',string)
                call PrintErrMsg(option)
            end select
          endif
      ! ---------------------------------------------------------------------- !
        case(ANY_STATE)
          ! temperature; 2nd dof ------------------------- !
          if (associated(general%temperature)) then
            select case(general%temperature%itype)
              case(DIRICHLET_BC)
                call PatchGetCouplerValueFromDataset(coupler,option, &
                       patch%grid,general%temperature%dataset,iconn,temperature)
                coupler%flow_aux_real_var(TWO_INTEGER,iconn) = temperature
                dof2 = PETSC_TRUE
                coupler%flow_bc_type(GENERAL_ENERGY_EQUATION_INDEX) = &
                                                                    DIRICHLET_BC
              case default
                string = GetSubConditionName(general%temperature%itype)
                option%io_buffer = &
                  FlowConditionUnknownItype(coupler%flow_condition, &
                  'GENERAL_MODE gas state temperature ',string)
                call PrintErrMsg(option)
            end select
          endif
      ! ---------------------------------------------------------------------- !
      end select
    enddo
  endif

  select case(flow_condition%iphase)
    case(MULTI_STATE)
    case(TWO_PHASE_STATE)
    case(LIQUID_STATE)
    ! ---> this code just prints an error, I think:
      if (general%liquid_pressure%itype == HYDROSTATIC_BC) then
        do iconn=1,coupler%connection_set%num_connections
          if (coupler%flow_aux_int_var(ONE_INTEGER,iconn) == TWO_PHASE_STATE) then
            !geh: This cannot possibly be working.  real_count needs to be incremented
            !     but what variable is mapped?  Need to figure out how real_count
            !     factors into the hydrostatic condition
            option%io_buffer = 'Need to fix PatchUpdateCouplerAuxVarsG() ' // &
              'for a variable saturated hydrostatic condition.'
            call PrintErrMsgByRank(option)

            ! we have to remap the capillary pressure to saturation and
            ! temperature to air pressure
            local_id = coupler%connection_set%id_dn(iconn)
            ghosted_id = patch%grid%nL2G(local_id)
            ! we have to convert capillary pressure (stored in air
            ! pressure index) to a saturation
            ! index     variable
            !  1        coupler%flow_aux_mapping(GENERAL_GAS_PRESSURE_INDEX) = 1
            !        air pressure in this case hijacked for capillary pressure
            !  2        coupler%flow_aux_mapping(GENERAL_AIR_PRESSURE_INDEX) = 2
            !  3        coupler%flow_aux_mapping(GENERAL_TEMPERATURE_INDEX) = 3
            gas_pressure = coupler%flow_aux_real_var( &
                      coupler%flow_aux_mapping( &
                        GENERAL_GAS_PRESSURE_INDEX),iconn)
            p_cap = coupler%flow_aux_real_var( &
                      coupler%flow_aux_mapping( &
                        GENERAL_AIR_PRESSURE_INDEX),iconn)
            temperature = coupler%flow_aux_real_var( &
                            coupler%flow_aux_mapping( &
                              GENERAL_TEMPERATURE_INDEX),iconn)
            coupler%flow_aux_mapping(general_2ph_energy_dof) = real_count
            if (general_2ph_energy_dof == GENERAL_TEMPERATURE_INDEX) then
              coupler%flow_aux_real_var(real_count,1:num_connections) = &
                temperature
            else
              call EOSWaterSaturationPressure(temperature,p_sat,ierr)
              coupler%flow_aux_real_var( &
                coupler%flow_aux_mapping( &
                  GENERAL_AIR_PRESSURE_INDEX),iconn) = &
                    gas_pressure - p_sat ! air pressure
            endif
            call patch%characteristic_curves_array(patch%cc_id(ghosted_id))% &
                   ptr%saturation_function%Saturation(p_cap,s_liq, &
                   dummy_real,option)
            ! %flow_aux_mapping(GENERAL_GAS_SATURATION_INDEX) set to 3 in hydrostatic
            coupler%flow_aux_real_var( &
              coupler%flow_aux_mapping( &
                GENERAL_GAS_SATURATION_INDEX),iconn) = &
              1.d0 - s_liq
          endif
        enddo
        coupler%flow_bc_type(1) = HYDROSTATIC_BC
        coupler%flow_bc_type(2:3) = DIRICHLET_BC
      else
      endif
    case(GAS_STATE)
    case(ANY_STATE)
  end select

  if (associated(general%liquid_flux)) then
    coupler%flow_bc_type(GENERAL_LIQUID_EQUATION_INDEX) = NEUMANN_BC
    select type(selector => general%liquid_flux%dataset)
      class is(dataset_ascii_type)
        coupler%flow_aux_real_var(ONE_INTEGER,1:num_connections) = &
                                           general%liquid_flux%dataset%rarray(1)
        dof1 = PETSC_TRUE
      class is(dataset_gridded_hdf5_type)
        call PatchVerifyDatasetGriddedForFlux(selector,coupler,option)
        call PatchUpdateCouplerGridDataset(coupler,option,patch%grid,selector, &
                                           ONE_INTEGER)
        dof1 = PETSC_TRUE
      class default
        call PrintMsg(option,'general%liquid_flux%dataset')
        call DatasetUnknownClass(selector,option, &
                                 'PatchUpdateCouplerAuxVarsG')
    end select
  endif
  if (associated(general%energy_flux)) then
    coupler%flow_bc_type(GENERAL_ENERGY_EQUATION_INDEX) = NEUMANN_BC
    select type(selector => general%energy_flux%dataset)
      class is(dataset_ascii_type)
        coupler%flow_aux_real_var(TWO_INTEGER,1:num_connections) = &
          general%energy_flux%dataset%rarray(1)
        dof2 = PETSC_TRUE
      class is(dataset_gridded_hdf5_type)
        call PatchVerifyDatasetGriddedForFlux(selector,coupler,option)
        call PatchUpdateCouplerGridDataset(coupler,option,patch%grid,selector, &
                                           TWO_INTEGER)
        dof2 = PETSC_TRUE
      class default
        call PrintMsg(option,'general%energy_flux%dataset')
        call DatasetUnknownClass(selector,option, &
                                 'PatchUpdateCouplerAuxVarsG')
    end select
  endif
  if (associated(general%gas_flux)) then
    coupler%flow_bc_type(GENERAL_GAS_EQUATION_INDEX) = NEUMANN_BC
    select type(selector => general%gas_flux%dataset)
      class is(dataset_ascii_type)
        coupler%flow_aux_real_var(THREE_INTEGER,1:num_connections) = &
                                              general%gas_flux%dataset%rarray(1)
        dof3 = PETSC_TRUE
      class is(dataset_gridded_hdf5_type)
        call PatchVerifyDatasetGriddedForFlux(selector,coupler,option)
        call PatchUpdateCouplerGridDataset(coupler,option,patch%grid,selector, &
                                           THREE_INTEGER)
        dof3 = PETSC_TRUE
      class default
        call PrintMsg(option,'general%gas_flux%dataset')
        call DatasetUnknownClass(selector,option, &
                                 'PatchUpdateCouplerAuxVarsG')
    end select
  endif

  if (associated(general%rate)) then
    select case(general%rate%itype)
      case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS)
        call PatchScaleSourceSink(patch,coupler,general%rate%isubtype,option)
    end select
  endif

  dof_count_global = 0
  dof_count_local = 0
  if (dof1) dof_count_local(1) = 1
  if (dof2) dof_count_local(2) = 1
  if (dof3) dof_count_local(3) = 1
  call MPI_Allreduce(dof_count_local,dof_count_global,THREE_INTEGER_MPI, &
                     MPI_INTEGER,MPI_SUM,option%mycomm,ierr)
  if (dof_count_global(1) > 0) dof1 = PETSC_TRUE
  if (dof_count_global(2) > 0) dof2 = PETSC_TRUE
  if (dof_count_global(3) > 0) dof3 = PETSC_TRUE
  ! need to check if these dofs are true on any process, because the
  ! boundary condition might be split up on 2 or more processes
  if (.not.dof1 .or. .not.dof2 .or. .not.dof3) then
    if (coupler%itype .ne. SRC_SINK_COUPLER_TYPE) then
      option%io_buffer = 'Error with GENERAL_MODE phase boundary condition: &
                          &Missing dof.'
      call PrintErrMsg(option)
    endif
  endif

end subroutine PatchUpdateCouplerAuxVarsG

! ************************************************************************** !

subroutine PatchUpdateCouplerAuxVarsH(patch,coupler,option)
  !
  ! Updates flow auxiliary variables associated
  ! with a coupler for H_MODE
  !
  ! Author: Michael Nole
  ! Date: 07/22/19
  !

  use Option_module
  use Condition_module
  use Hydrostatic_module
  use EOS_Water_module
  use Utility_module

  use Grid_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_HDF5_class
  use Dataset_Ascii_class
  use Dataset_module
  use String_module
  use Hydrate_Aux_module

  implicit none

  type(patch_type) :: patch
  type(coupler_type), pointer :: coupler
  type(option_type) :: option

  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(flow_hydrate_condition_type), pointer :: hydrate
  PetscBool :: update
  PetscBool :: dof1, dof2, dof3
  PetscReal :: temperature, p_sat, p_cap, s_liq, xmol
  PetscReal :: relative_humidity
  PetscReal :: liq_sat, gas_sat, hyd_sat, ice_sat, air_pressure, gas_pressure, &
               liq_pressure
  PetscReal :: dummy_real
  PetscReal :: x(option%nflowdof)
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr

  PetscInt :: idof, num_connections,sum_connection
  PetscInt :: iconn, local_id, ghosted_id
  PetscInt :: real_count
  PetscInt :: dof_count_local(3)
  PetscInt :: dof_count_global(3)
  PetscReal, parameter :: min_two_phase_gas_pressure = 3.d3

  num_connections = coupler%connection_set%num_connections

  flow_condition => coupler%flow_condition

  hydrate => flow_condition%hydrate
  dof1 = PETSC_FALSE
  dof2 = PETSC_FALSE
  dof3 = PETSC_FALSE
  real_count = 0

  ! mapping of flow_aux_mapping to the flow_aux_real_var array:
  if (associated(coupler%flow_aux_mapping)) then
    ! liquid and gas pressure are set to 1st dof index
    ! liquid flux is set to 1st dof index
    coupler%flow_aux_mapping(HYDRATE_GAS_PRESSURE_INDEX) = 1
    coupler%flow_aux_mapping(HYDRATE_LIQUID_PRESSURE_INDEX) = 1
    coupler%flow_aux_mapping(HYDRATE_LIQUID_FLUX_INDEX) = 1
    coupler%flow_aux_mapping(HYDRATE_LIQ_SATURATION_INDEX) = 1
    ! temperature is set to 2nd dof index
    ! energy flux is set to 2nd dof index
    coupler%flow_aux_mapping(HYDRATE_TEMPERATURE_INDEX) = 2
    coupler%flow_aux_mapping(HYDRATE_ENERGY_FLUX_INDEX) = 2
    coupler%flow_aux_mapping(HYDRATE_ICE_SATURATION_INDEX) = 2
    ! air mole fraction, gas sat., and air pressure are set to 3rd dof index
    ! gas flux is set to 3rd dof index
    coupler%flow_aux_mapping(HYDRATE_LIQ_MOLE_FRACTION_INDEX) = 3
    coupler%flow_aux_mapping(HYDRATE_GAS_SATURATION_INDEX) = 3
    coupler%flow_aux_mapping(HYDRATE_AIR_PRESSURE_INDEX) = 3
    coupler%flow_aux_mapping(HYDRATE_GAS_FLUX_INDEX) = 3
    coupler%flow_aux_mapping(HYDRATE_HYD_SATURATION_INDEX) = 3

  endif

  select case(flow_condition%iphase)
    case(HYD_MULTI_STATE)
      select type(dataset => hydrate%gas_saturation%dataset)
        class is(dataset_ascii_type)
          gas_sat = hydrate%gas_saturation%dataset%rarray(1)
          if (gas_sat > 0.d0 .and. gas_sat < 1.d0) then
            coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,1:num_connections) = &
              GA_STATE
          ! Cannot user gas_sat == 0.d0 or Equal(gas_sat,0.d0) as optimization
          ! in the Intel compiler changes the answer.
          else if (gas_sat < 0.5d0) then
            coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,1:num_connections) = &
              L_STATE
          else
            coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,1:num_connections) = &
              G_STATE
          endif
        class is(dataset_gridded_hdf5_type)
          ! If the gas pressure dataset is defined, we must ensure that its
          ! minimum pressure is greater than min_two_phase_gas_pressure.
          ! Otherwise, a zero gas pressure may be incorrectly used within
          ! an interpolation of gas pressure for a two phase cell mixed
          ! with a single phase liquid cell
          if (associated(hydrate%gas_pressure)) then
            dummy_real = &
              DatasetGetMinRValue(hydrate%gas_pressure%dataset,option)
            if (dummy_real < min_two_phase_gas_pressure) then
              option%io_buffer = 'Minimum gas pressure exceeded for &
                    &FLOW_CONDITION "' // trim(flow_condition%name) // &
                    '": ' // trim(StringFormatDouble(dummy_real)) // '.'
              call PrintErrMsg(option)
            endif
          endif
          do iconn = 1, num_connections
            call PatchGetCouplerValueFromDataset(coupler,option,patch%grid, &
                                  hydrate%gas_saturation%dataset,iconn,gas_sat)
            if (gas_sat > 0.d0 .and. gas_sat < 1.d0) then
              coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,iconn) = &
                GA_STATE
            else if (gas_sat < 0.5d0) then
              coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,iconn) = L_STATE
            else
              coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,iconn) = G_STATE
            endif
          enddo
        class default
          call PrintMsg(option,'hydrate%gas_saturation%dataset,MULTI_STATE')
          call DatasetUnknownClass(dataset,option, &
                                   'PatchUpdateCouplerAuxVarsH')
      end select
    case(L_STATE)
      coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,1:num_connections) = &
        L_STATE
      if (hydrate%liquid_pressure%itype == HYDROSTATIC_BC) then
        if (hydrate%mole_fraction%itype /= DIRICHLET_BC) then
          option%io_buffer = 'Hydrostatic liquid state pressure BC for &
            &flow condition "' // trim(flow_condition%name) // &
            '" requires a mole fraction BC of type DIRICHLET.'
          call PrintErrMsg(option)
        endif
        if (hydrate%temperature%itype /= DIRICHLET_BC) then
          option%io_buffer = 'Hydrostatic liquid state pressure BC for &
            &flow condition "' // trim(flow_condition%name) // &
            '" requires a temperature BC of type DIRICHLET.'
          call PrintErrMsg(option)
        endif
        call HydrostaticUpdateCoupler(coupler,option,patch%grid)
        do iconn = 1, num_connections
          if (coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,iconn) /= &
              L_STATE) then
            select case(coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,iconn))
              case(G_STATE)
                string = 'gas state'
              case(GA_STATE)
                string = 'two phase state'
              case(HYD_ANY_STATE)
                string = 'any phase state'
            end select
            option%io_buffer = 'A ' // trim(string) // ' cell was found &
              &within a HYDROSTATIC_BC boundary condition for HYDRATE mode. &
              &A hydrostatic boundary condition may not be used to set &
              &state variables in the vadose zone for HYDRATE mode.'
            call PrintErrMsg(option)
          endif
        enddo
        dof1 = PETSC_TRUE; dof2 = PETSC_TRUE; dof3 = PETSC_TRUE;
      endif
    case(G_STATE)
      coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,1:num_connections) = &
        G_STATE
    case(H_STATE)
      coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,1:num_connections) = &
        H_STATE
    case(I_STATE)
      coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,1:num_connections) = &
        I_STATE
    case(GA_STATE)
      coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,1:num_connections) = &
        GA_STATE
        ! no need to loop in the next do loop if its all the same state, which
        ! you know from flow_condition%iphase
    case(HG_STATE)
      coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,1:num_connections) = &
        HG_STATE
    case(HA_STATE)
      coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,1:num_connections) = HA_STATE
      if (associated(hydrate%liquid_pressure) .and.  hydrate% &
                              liquid_pressure%itype == HYDROSTATIC_BC) then
        if (hydrate%temperature%itype /= DIRICHLET_BC) then
          option%io_buffer = 'Hydrostatic hydrate-aq. state pressure BC for &
            &flow condition "' // trim(flow_condition%name) // &
            '" requires a temperature BC of type DIRICHLET.'
          call PrintErrMsg(option)
        endif
        call HydrostaticUpdateCoupler(coupler,option,patch%grid)
        do iconn = 1, num_connections
          if (coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,iconn) /= &
              L_STATE .and. coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,iconn)&
              /= HA_STATE) then
            string = 'Non L or HA state'
            option%io_buffer = 'A ' // trim(string) // ' cell was found &
              &within a HYDROSTATIC_BC boundary condition for HYDRATE mode. &
              &A hydrostatic boundary condition may not be used to set &
              &state variables in the vadose zone for HYDRATE mode.'
            call PrintErrMsg(option)
          endif
        enddo
        dof1 = PETSC_TRUE; dof2 = PETSC_TRUE; dof3 = PETSC_TRUE;
      endif
    case(HI_STATE)
      coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,1:num_connections) = &
        HI_STATE
    case(GI_STATE)
      coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,1:num_connections) = &
        GI_STATE
    case(HGA_STATE)
      coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,1:num_connections) = &
        HGA_STATE
    case(HAI_STATE)
      coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,1:num_connections) = &
        HAI_STATE
    case(GAI_STATE)
      coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,1:num_connections) = &
        GAI_STATE
    case(HGI_STATE)
      coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,1:num_connections) = &
        HGI_STATE
    case(QUAD_STATE)
      coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,1:num_connections) = &
        QUAD_STATE
    case(HYD_ANY_STATE)
      if (associated(coupler%flow_aux_int_var)) then ! not used with rate
        coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,1:num_connections) = &
          HYD_ANY_STATE
      endif
  end select

  ! loop over each connection in the coupler and check its state
  ! set the flow_aux_mapping, flow_aux_real_var, etc on a connection
  ! basis rather than in coupler chunks
  ! this might be slower since we need to loop over all the connections
  ! but it makes the algorithm more general
  if (associated(coupler%flow_aux_int_var)) then
    do iconn = 1, num_connections
      select case(coupler%flow_aux_int_var(HYDRATE_STATE_INDEX,iconn))
      ! ---------------------------------------------------------------------- !
        case(L_STATE)
          if (hydrate%liquid_pressure%itype == HYDROSTATIC_BC) then
            if (hydrate%mole_fraction%itype /= DIRICHLET_BC) then
              option%io_buffer = 'Hydrostatic liquid state pressure BC for &
                &flow condition "' // trim(flow_condition%name) // &
                '" requires a mole fraction BC of type DIRICHLET.'
              call PrintErrMsg(option)
            endif
            if (hydrate%temperature%itype /= DIRICHLET_BC) then
              option%io_buffer = 'Hydrostatic liquid state pressure BC for &
                &flow condition "' // trim(flow_condition%name) // &
                '" requires a temperature BC of type DIRICHLET.'
              call PrintErrMsg(option)
            endif
            ! ---> see code that just prints error
            coupler%flow_bc_type(1) = HYDROSTATIC_BC
            coupler%flow_bc_type(2:3) = DIRICHLET_BC
          else
          ! liquid pressure; 1st dof --------------------- !
            select case(hydrate%liquid_pressure%itype)
              case(DIRICHLET_BC,DIRICHLET_SEEPAGE_BC)
                call PatchGetCouplerValueFromDataset(coupler,option, &
                  patch%grid,hydrate%liquid_pressure%dataset,iconn,liq_pressure)
                coupler%flow_aux_real_var(ONE_INTEGER,iconn) = liq_pressure
                dof1 = PETSC_TRUE
                coupler%flow_bc_type(HYDRATE_LIQUID_EQUATION_INDEX) = &
                                                                    DIRICHLET_BC
              case default
                string = GetSubConditionName(hydrate%liquid_pressure%itype)
                option%io_buffer = &
                  FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE liquid state liquid pressure ',string)
                call PrintErrMsg(option)
            end select
          ! temperature; 2nd dof ------------------------- !
            select case(hydrate%temperature%itype)
              case(DIRICHLET_BC)
                call PatchGetCouplerValueFromDataset(coupler,option, &
                       patch%grid,hydrate%temperature%dataset,iconn,temperature)
                coupler%flow_aux_real_var(TWO_INTEGER,iconn) = temperature
                dof2 = PETSC_TRUE
                coupler%flow_bc_type(HYDRATE_ENERGY_EQUATION_INDEX) = &
                                                                    DIRICHLET_BC
              case default
                string = GetSubConditionName(hydrate%temperature%itype)
                option%io_buffer = &
                  FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE liquid state temperature ',string)
                call PrintErrMsg(option)
            end select
          ! mole fraction; 3rd dof ----------------------- !
            select case(hydrate%mole_fraction%itype)
              case(DIRICHLET_BC)
                call PatchGetCouplerValueFromDataset(coupler,option, &
                            patch%grid,hydrate%mole_fraction%dataset,iconn,xmol)
                coupler%flow_aux_real_var(THREE_INTEGER,iconn) = xmol
                dof3 = PETSC_TRUE
                coupler%flow_bc_type(HYDRATE_GAS_EQUATION_INDEX) = DIRICHLET_BC
              case default
                string = GetSubConditionName(hydrate%mole_fraction%itype)
                option%io_buffer = &
                  FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE liquid state mole fraction ',string)
                call PrintErrMsg(option)
            end select
          endif
      ! ---------------------------------------------------------------------- !
        case(G_STATE)
          gas_pressure = UNINITIALIZED_DOUBLE
          temperature = UNINITIALIZED_DOUBLE
          ! gas pressure; 1st dof ------------------------ !
          select case(hydrate%gas_pressure%itype)
            case(DIRICHLET_BC,DIRICHLET_SEEPAGE_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                  patch%grid,hydrate%gas_pressure%dataset,iconn,gas_pressure)
              coupler%flow_aux_real_var(ONE_INTEGER,iconn) = gas_pressure
              dof1 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_GAS_EQUATION_INDEX) = &
                      hydrate%gas_pressure%itype !DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%gas_pressure%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                'HYDRATE MODE gas state gas pressure',string)
              call PrintErrMsg(option)
          end select
          ! temperature; 2nd dof ------------------------- !
          select case(hydrate%temperature%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                       patch%grid,hydrate%temperature%dataset,iconn,temperature)
              coupler%flow_aux_real_var(TWO_INTEGER,iconn) = &
                temperature
              dof2 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_ENERGY_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%temperature%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                'HYDRATE MODE gas state temperature',string)
              call PrintErrMsg(option)
          end select
          ! air mole fraction; 3rd dof ------------------- !
          if (associated(hydrate%mole_fraction)) then
            select case(hydrate%mole_fraction%itype)
              case(DIRICHLET_BC)
                if (Uninitialized(gas_pressure) .or. &
                    Uninitialized(temperature)) then
                  option%io_buffer = 'GAS_PRESSURE or TEMPERATURE not set &
                    &correctly in flow condition "' // &
                    trim(flow_condition%name) // '".'
                  call PrintErrMsg(option)
                endif
                call PatchGetCouplerValueFromDataset(coupler,option, &
                            patch%grid,hydrate%mole_fraction%dataset,iconn,xmol)
                air_pressure = xmol * gas_pressure
                call EOSWaterSaturationPressure(temperature,p_sat,ierr)
                if (gas_pressure - air_pressure >= p_sat) then
                  option%io_buffer = 'MOLE_FRACTION set in flow &
                    &condition "' // trim(flow_condition%name) // &
                    '" results in a vapor pressure exceeding the water &
                    &saturation pressure, which indicates that a two-phase &
                    &state with GAS_PRESSURE and GAS_SATURATION should be used.'
                  call PrintErrMsg(option)
                endif
                coupler%flow_aux_real_var(THREE_INTEGER,iconn) = air_pressure
                dof3 = PETSC_TRUE
                coupler%flow_bc_type(HYDRATE_LIQUID_EQUATION_INDEX) = &
                                                                    DIRICHLET_BC
              case default
                string = GetSubConditionName(hydrate%mole_fraction%itype)
                option%io_buffer = &
                  FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE air mole fraction',string)
                call PrintErrMsg(option)
            end select
        ! relative humidity; 3rd dof ------------------- !
          else
            select case(hydrate%relative_humidity%itype)
              case(DIRICHLET_BC)
                if (Uninitialized(gas_pressure) .or. &
                    Uninitialized(temperature)) then
                  option%io_buffer = 'GAS_PRESSURE or TEMPERATURE not set &
                    &correctly in flow condition "' // &
                    trim(flow_condition%name) // '".'
                  call PrintErrMsg(option)
                endif
                call PatchGetCouplerValueFromDataset(coupler,option, &
                  patch%grid,hydrate%relative_humidity%dataset, &
                  iconn,relative_humidity)  ! relative_humidity is in percent
                if (relative_humidity < 0.d0 .or. &
                    relative_humidity > 100.d0) then
                  option%io_buffer = 'RELATIVE_HUMIDITY in flow &
                    &condition "' // trim(flow_condition%name) // '" outside &
                    &bounds of 0-100%.'
                  call PrintErrMsg(option)
                endif
                call EOSWaterSaturationPressure(temperature,p_sat,ierr)
                                  ! convert from % to fraction
                air_pressure = gas_pressure - relative_humidity*1.d-2*p_sat
                coupler%flow_aux_real_var(THREE_INTEGER,iconn) = air_pressure
                dof3 = PETSC_TRUE
                coupler%flow_bc_type(HYDRATE_LIQUID_EQUATION_INDEX) = &
                                                                    DIRICHLET_BC
              case default
                string = GetSubConditionName(hydrate%relative_humidity%itype)
                option%io_buffer = &
                  FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE relative humidity',string)
                call PrintErrMsg(option)
            end select
          endif
    ! ---------------------------------------------------------------------- !
        case(H_STATE)
          ! gas pressure; 1st dof ------------------------ !
          select case(hydrate%gas_pressure%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                     patch%grid,hydrate%gas_pressure%dataset,iconn,gas_pressure)
              coupler%flow_aux_real_var(ONE_INTEGER,iconn) = gas_pressure
              dof1 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_LIQUID_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%gas_pressure%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE H-state gas pressure ',string)
              call PrintErrMsg(option)
          end select
          ! temperature; 2nd dof ------------------------- !
          select case(hydrate%temperature%itype)
              case(DIRICHLET_BC)
                call PatchGetCouplerValueFromDataset(coupler,option, &
                       patch%grid,hydrate%temperature%dataset,iconn,temperature)
                coupler%flow_aux_real_var(TWO_INTEGER,iconn) = temperature
                dof2 = PETSC_TRUE
                coupler%flow_bc_type(HYDRATE_ENERGY_EQUATION_INDEX) = &
                                                                    DIRICHLET_BC
              case default
                string = GetSubConditionName(hydrate%temperature%itype)
                option%io_buffer = &
                  FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE H-state state temperature ',string)
                call PrintErrMsg(option)
            end select
          ! mole fraction; 3rd dof ----------------------- !
            select case(hydrate%mole_fraction%itype)
              case(DIRICHLET_BC)
                xmol = MOL_RATIO_METH
                coupler%flow_aux_real_var(THREE_INTEGER,iconn) = xmol
                dof3 = PETSC_TRUE
                coupler%flow_bc_type(HYDRATE_GAS_EQUATION_INDEX) = DIRICHLET_BC
              case default
                string = GetSubConditionName(hydrate%mole_fraction%itype)
                option%io_buffer = &
                  FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE H-state mole fraction ',string)
                call PrintErrMsg(option)
            end select
    ! ---------------------------------------------------------------------- !
        case(I_STATE)
          ! gas pressure; 1st dof ------------------------ !
          select case(hydrate%gas_pressure%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                     patch%grid,hydrate%gas_pressure%dataset,iconn,gas_pressure)
              coupler%flow_aux_real_var(ONE_INTEGER,iconn) = gas_pressure
              dof1 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_LIQUID_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%gas_pressure%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE I-state gas pressure ',string)
              call PrintErrMsg(option)
          end select
          ! temperature; 2nd dof ------------------------- !
            select case(hydrate%temperature%itype)
              case(DIRICHLET_BC)
                call PatchGetCouplerValueFromDataset(coupler,option, &
                       patch%grid,hydrate%temperature%dataset,iconn,temperature)
                coupler%flow_aux_real_var(TWO_INTEGER,iconn) = temperature
                dof2 = PETSC_TRUE
                coupler%flow_bc_type(HYDRATE_ENERGY_EQUATION_INDEX) = &
                                                                    DIRICHLET_BC
              case default
                string = GetSubConditionName(hydrate%temperature%itype)
                option%io_buffer = &
                  FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE I-state state temperature ',string)
                call PrintErrMsg(option)
            end select
          ! mole fraction; 3rd dof ----------------------- !
            select case(hydrate%mole_fraction%itype)
              case(DIRICHLET_BC)
                xmol = 0.d0
                coupler%flow_aux_real_var(THREE_INTEGER,iconn) = xmol
                dof3 = PETSC_TRUE
                coupler%flow_bc_type(HYDRATE_GAS_EQUATION_INDEX) = DIRICHLET_BC
              case default
                string = GetSubConditionName(hydrate%mole_fraction%itype)
                option%io_buffer = &
                  FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE I-state mole fraction ',string)
                call PrintErrMsg(option)
            end select
    ! ---------------------------------------------------------------------- !
        case(GA_STATE, HG_STATE)
          ! gas pressure; 1st dof ------------------------ !
          select case(hydrate%gas_pressure%itype)
            case(DIRICHLET_BC,DIRICHLET_SEEPAGE_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                     patch%grid,hydrate%gas_pressure%dataset,iconn,gas_pressure)
              coupler%flow_aux_real_var(ONE_INTEGER,iconn) = gas_pressure
              dof1 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_LIQUID_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%gas_pressure%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE GA-State gas pressure ',string)
              call PrintErrMsg(option)
          end select
          ! temperature; 2nd dof ------------------------- !
          select case(hydrate%temperature%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                       patch%grid,hydrate%temperature%dataset,iconn,temperature)
              coupler%flow_aux_real_var(TWO_INTEGER,iconn) = temperature
              dof2 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_ENERGY_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%temperature%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE GA-State temperature ',string)
              call PrintErrMsg(option)
          end select
          ! gas saturation; 3rd dof ---------------------- !
          select case(hydrate%gas_saturation%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                        patch%grid,hydrate%gas_saturation%dataset,iconn,gas_sat)
              coupler%flow_aux_real_var(THREE_INTEGER,iconn) = gas_sat
              dof3 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_GAS_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%gas_saturation%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE GA-State gas saturation ',string)
              call PrintErrMsg(option)
          end select
    ! ---------------------------------------------------------------------- !
        case(HA_STATE, HI_STATE)
          ! gas pressure; 1st dof ------------------------ !
          if (associated(hydrate%gas_pressure)) then
            select case(hydrate%gas_pressure%itype)
              case(DIRICHLET_BC)
                call PatchGetCouplerValueFromDataset(coupler,option, &
                     patch%grid,hydrate%gas_pressure%dataset,iconn,gas_pressure)
                coupler%flow_aux_real_var(ONE_INTEGER,iconn) = gas_pressure
                dof1 = PETSC_TRUE
                coupler%flow_bc_type(HYDRATE_LIQUID_EQUATION_INDEX) = &
                                                                  DIRICHLET_BC
              case default
                string = GetSubConditionName(hydrate%gas_pressure%itype)
                option%io_buffer = &
                  FlowConditionUnknownItype(coupler%flow_condition, &
                    'HYDRATE MODE HA-state gas pressure ',string)
                call PrintErrMsg(option)
            end select
            ! temperature; 2nd dof ------------------------- !
            select case(hydrate%temperature%itype)
              case(DIRICHLET_BC)
                call PatchGetCouplerValueFromDataset(coupler,option, &
                       patch%grid,hydrate%temperature%dataset,iconn,temperature)
                coupler%flow_aux_real_var(TWO_INTEGER,iconn) = temperature
                dof3 = PETSC_TRUE
                coupler%flow_bc_type(HYDRATE_ENERGY_EQUATION_INDEX) = &
                                                                 DIRICHLET_BC
              case default
                string = GetSubConditionName(hydrate%temperature%itype)
                option%io_buffer = &
                  FlowConditionUnknownItype(coupler%flow_condition, &
                    'HYDRATE MODE HA-state temperature ',string)
                call PrintErrMsg(option)
            end select

          else
            if (hydrate%liquid_pressure%itype == HYDROSTATIC_BC) then
              if (hydrate%temperature%itype /= DIRICHLET_BC) then
                option%io_buffer = 'Hydrostatic HA state pressure BC for &
                  &flow condition "' // trim(flow_condition%name) // &
                  '" requires a temperature BC of type DIRICHLET.'
                call PrintErrMsg(option)
              endif
              ! ---> see code that just prints error
              coupler%flow_bc_type(1) = HYDROSTATIC_BC
              coupler%flow_bc_type(2:3) = DIRICHLET_BC
            else
            ! liquid pressure; 1st dof --------------------- !
              select case(hydrate%liquid_pressure%itype)
                case(DIRICHLET_BC)
                  call PatchGetCouplerValueFromDataset(coupler,option,patch% &
                      grid,hydrate%liquid_pressure%dataset,iconn,liq_pressure)
                  coupler%flow_aux_real_var(ONE_INTEGER,iconn) = liq_pressure
                  dof1 = PETSC_TRUE
                  coupler%flow_bc_type(HYDRATE_LIQUID_EQUATION_INDEX) = &
                                                                  DIRICHLET_BC
                case default
                  string = GetSubConditionName(hydrate%liquid_pressure%itype)
                  option%io_buffer = &
                    FlowConditionUnknownItype(coupler%flow_condition, &
                    'HYDRATE MODE HA-state liquid pressure ',string)
                  call PrintErrMsg(option)
              end select
            ! temperature; 2nd dof ------------------------- !
              select case(hydrate%temperature%itype)
                case(DIRICHLET_BC)
                  call PatchGetCouplerValueFromDataset(coupler,option, &
                       patch%grid,hydrate%temperature%dataset,iconn,temperature)
                  coupler%flow_aux_real_var(TWO_INTEGER,iconn) = temperature
                  dof3 = PETSC_TRUE
                  coupler%flow_bc_type(HYDRATE_ENERGY_EQUATION_INDEX) = &
                                                                    DIRICHLET_BC
                case default
                  string = GetSubConditionName(hydrate%temperature%itype)
                  option%io_buffer = &
                    FlowConditionUnknownItype(coupler%flow_condition, &
                    'HYDRATE MODE HA-state temperature ',string)
                  call PrintErrMsg(option)
              end select
            endif
          endif
          ! hydrate saturation; 3rd dof ---------------------- !
          select case(hydrate%hydrate_saturation%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                        patch%grid,hydrate%hydrate_saturation%dataset,iconn, &
                        hyd_sat)
              coupler%flow_aux_real_var(THREE_INTEGER,iconn) = hyd_sat
              dof2 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_GAS_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%hydrate_saturation%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE HA-state hydrate saturation ',string)
              call PrintErrMsg(option)
          end select
      ! ---------------------------------------------------------------------- !
        case(GI_STATE)
          ! gas pressure; 1st dof ------------------------ !
          select case(hydrate%gas_pressure%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                     patch%grid,hydrate%gas_pressure%dataset,iconn,gas_pressure)
              coupler%flow_aux_real_var(ONE_INTEGER,iconn) = gas_pressure
              dof1 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_LIQUID_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%gas_pressure%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE GI-State gas pressure ',string)
              call PrintErrMsg(option)
          end select
          ! temperature; 2nd dof ------------------------- !
          select case(hydrate%temperature%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                       patch%grid,hydrate%temperature%dataset,iconn,temperature)
              coupler%flow_aux_real_var(TWO_INTEGER,iconn) = temperature
              dof2 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_ENERGY_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%temperature%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE GI-state temperature ',string)
              call PrintErrMsg(option)
          end select
          ! ice saturation; 3rd dof ---------------------- !
          select case(hydrate%gas_saturation%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                        patch%grid,hydrate%ice_saturation%dataset,iconn,ice_sat)
              coupler%flow_aux_real_var(THREE_INTEGER,iconn) = ice_sat
              dof3 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_GAS_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%ice_saturation%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE GI-state gas saturation ',string)
              call PrintErrMsg(option)
          end select
      ! ---------------------------------------------------------------------- !
        case(AI_STATE)
          ! gas pressure; 1st dof ------------------------ !
          select case(hydrate%gas_pressure%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                     patch%grid,hydrate%gas_pressure%dataset,iconn,gas_pressure)
              coupler%flow_aux_real_var(ONE_INTEGER,iconn) = gas_pressure
              dof1 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_LIQUID_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%gas_pressure%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE AI-State gas pressure ',string)
              call PrintErrMsg(option)
          end select
          ! temperature; 2nd dof ------------------------- !
          select case(hydrate%temperature%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                       patch%grid,hydrate%temperature%dataset,iconn,temperature)
              coupler%flow_aux_real_var(TWO_INTEGER,iconn) = temperature
              dof2 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_ENERGY_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%temperature%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE AI-state temperature ',string)
              call PrintErrMsg(option)
          end select
          ! liquid saturation; 3rd dof ---------------------- !
          select case(hydrate%gas_saturation%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                        patch%grid,hydrate%liquid_saturation%dataset,iconn, &
                        liq_sat)
              coupler%flow_aux_real_var(THREE_INTEGER,iconn) = liq_sat
              dof3 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_GAS_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%liquid_saturation%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE AI-state liquid saturation ',string)
              call PrintErrMsg(option)
          end select
      ! ---------------------------------------------------------------------- !
        case(HGA_STATE)
          ! gas saturation; 1st dof ------------------------ !
          select case(hydrate%gas_saturation%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                     patch%grid,hydrate%gas_saturation%dataset,iconn,gas_sat)
              coupler%flow_aux_real_var(ONE_INTEGER,iconn) = gas_sat
              dof1 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_LIQUID_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%gas_saturation%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE HGA-State gas saturation ',string)
              call PrintErrMsg(option)
          end select
          ! temperature; 2nd dof ------------------------- !
          select case(hydrate%temperature%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                       patch%grid,hydrate%temperature%dataset,iconn,temperature)
              coupler%flow_aux_real_var(TWO_INTEGER,iconn) = temperature
              dof2 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_ENERGY_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%temperature%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE HGA-state temperature ',string)
              call PrintErrMsg(option)
          end select
          ! hydrate saturation; 3rd dof ---------------------- !
          select case(hydrate%hydrate_saturation%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                        patch%grid,hydrate%hydrate_saturation%dataset,iconn, &
                        hyd_sat)
              coupler%flow_aux_real_var(THREE_INTEGER,iconn) = hyd_sat
              dof3 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_GAS_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%hydrate_saturation%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE HGA-state hydrate saturation ',string)
              call PrintErrMsg(option)
          end select
      ! ---------------------------------------------------------------------- !
        case(HAI_STATE, GAI_STATE)
          ! gas pressure; 1st dof ------------------------ !
          select case(hydrate%gas_pressure%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                     patch%grid,hydrate%gas_pressure%dataset,iconn,gas_pressure)
              coupler%flow_aux_real_var(ONE_INTEGER,iconn) = gas_pressure
              dof1 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_LIQUID_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%gas_pressure%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE HAI-State gas pressure ',string)
              call PrintErrMsg(option)
          end select
          ! ice saturation; 2nd dof ------------------------- !
          select case(hydrate%ice_saturation%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                       patch%grid,hydrate%ice_saturation%dataset,iconn,ice_sat)
              coupler%flow_aux_real_var(TWO_INTEGER,iconn) = ice_sat
              dof2 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_ENERGY_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%ice_saturation%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE HAI-state ice saturation ',string)
              call PrintErrMsg(option)
          end select
          ! hydrate saturation; 3rd dof ---------------------- !
          select case(hydrate%hydrate_saturation%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                        patch%grid,hydrate%hydrate_saturation%dataset,iconn, &
                        hyd_sat)
              coupler%flow_aux_real_var(THREE_INTEGER,iconn) = hyd_sat
              dof3 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_GAS_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%hydrate_saturation%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE HAI-state hydrate saturation ',string)
              call PrintErrMsg(option)
          end select
      ! ---------------------------------------------------------------------- !
        case(HGI_STATE)
          ! ice saturation; 1st dof ------------------------ !
          select case(hydrate%ice_saturation%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                     patch%grid,hydrate%ice_saturation%dataset,iconn,ice_sat)
              coupler%flow_aux_real_var(ONE_INTEGER,iconn) = ice_sat
              dof1 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_LIQUID_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%ice_saturation%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE HGI-State ice saturation ',string)
              call PrintErrMsg(option)
          end select
          ! temperature; 2nd dof ------------------------- !
          select case(hydrate%temperature%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                       patch%grid,hydrate%temperature%dataset,iconn,temperature)
              coupler%flow_aux_real_var(TWO_INTEGER,iconn) = temperature
              dof2 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_ENERGY_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%temperature%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE HGI-state temperature ',string)
              call PrintErrMsg(option)
          end select
          ! hydrate saturation; 3rd dof ---------------------- !
          select case(hydrate%hydrate_saturation%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                        patch%grid,hydrate%hydrate_saturation%dataset,iconn, &
                        hyd_sat)
              coupler%flow_aux_real_var(THREE_INTEGER,iconn) = hyd_sat
              dof3 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_GAS_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%hydrate_saturation%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE HGI-state hydrate saturation ',string)
              call PrintErrMsg(option)
          end select
      ! ---------------------------------------------------------------------- !
        case(QUAD_STATE)
          ! liquid saturation; 1st dof ------------------------ !
          select case(hydrate%liquid_saturation%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                     patch%grid,hydrate%liquid_saturation%dataset,iconn,liq_sat)
              coupler%flow_aux_real_var(ONE_INTEGER,iconn) = liq_sat
              dof1 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_LIQUID_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%liquid_saturation%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE QUAD-State liquid saturation ',string)
              call PrintErrMsg(option)
          end select
          ! ice satuation; 2nd dof ------------------------- !
          select case(hydrate%ice_saturation%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                       patch%grid,hydrate%ice_saturation%dataset,iconn,ice_sat)
              coupler%flow_aux_real_var(TWO_INTEGER,iconn) = ice_sat
              dof2 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_ENERGY_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%ice_saturation%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE QUAD-state ice saturation ',string)
              call PrintErrMsg(option)
          end select
          ! gas saturation; 3rd dof ---------------------- !
          select case(hydrate%gas_saturation%itype)
            case(DIRICHLET_BC)
              call PatchGetCouplerValueFromDataset(coupler,option, &
                        patch%grid,hydrate%gas_saturation%dataset,iconn, &
                        gas_sat)
              coupler%flow_aux_real_var(THREE_INTEGER,iconn) = gas_sat
              dof3 = PETSC_TRUE
              coupler%flow_bc_type(HYDRATE_GAS_EQUATION_INDEX) = DIRICHLET_BC
            case default
              string = GetSubConditionName(hydrate%gas_saturation%itype)
              option%io_buffer = &
                FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE QUAD-state gas saturation ',string)
              call PrintErrMsg(option)
          end select
      ! ---------------------------------------------------------------------- !
        case(HYD_ANY_STATE)
          ! temperature; 2nd dof ------------------------- !
          if (associated(hydrate%temperature)) then
            select case(hydrate%temperature%itype)
              case(DIRICHLET_BC)
                call PatchGetCouplerValueFromDataset(coupler,option, &
                       patch%grid,hydrate%temperature%dataset,iconn,temperature)
                coupler%flow_aux_real_var(TWO_INTEGER,iconn) = temperature
                dof2 = PETSC_TRUE
                coupler%flow_bc_type(HYDRATE_ENERGY_EQUATION_INDEX) = &
                                                                    DIRICHLET_BC
              case default
                string = GetSubConditionName(hydrate%temperature%itype)
                option%io_buffer = &
                  FlowConditionUnknownItype(coupler%flow_condition, &
                  'HYDRATE MODE gas state temperature ',string)
                call PrintErrMsg(option)
            end select
          endif
      ! ---------------------------------------------------------------------- !
      end select
    enddo
  endif

  if (associated(hydrate%liquid_flux)) then
    coupler%flow_bc_type(HYDRATE_LIQUID_EQUATION_INDEX) = NEUMANN_BC
    select type(selector => hydrate%liquid_flux%dataset)
      class is(dataset_ascii_type)
        coupler%flow_aux_real_var(ONE_INTEGER,1:num_connections) = &
                                           hydrate%liquid_flux%dataset%rarray(1)
        dof1 = PETSC_TRUE
      class is(dataset_gridded_hdf5_type)
        call PatchVerifyDatasetGriddedForFlux(selector,coupler,option)
        call PatchUpdateCouplerGridDataset(coupler,option,patch%grid,selector, &
                                           ONE_INTEGER)
        dof1 = PETSC_TRUE
      class default
        call PrintMsg(option,'hydrate%liquid_flux%dataset')
        call DatasetUnknownClass(selector,option, &
                                 'PatchUpdateCouplerAuxVarsH')
    end select
  endif
  if (associated(hydrate%energy_flux)) then
    coupler%flow_bc_type(HYDRATE_ENERGY_EQUATION_INDEX) = NEUMANN_BC
    select type(selector => hydrate%energy_flux%dataset)
      class is(dataset_ascii_type)
        coupler%flow_aux_real_var(TWO_INTEGER,1:num_connections) = &
          hydrate%energy_flux%dataset%rarray(1)
        dof2 = PETSC_TRUE
      class is(dataset_gridded_hdf5_type)
        call PatchVerifyDatasetGriddedForFlux(selector,coupler,option)
        call PatchUpdateCouplerGridDataset(coupler,option,patch%grid,selector, &
                                           TWO_INTEGER)
        dof2 = PETSC_TRUE
      class default
        call PrintMsg(option,'hydrate%energy_flux%dataset')
        call DatasetUnknownClass(selector,option, &
                                 'PatchUpdateCouplerAuxVarsG')
    end select
  endif
  if (associated(hydrate%gas_flux)) then
    coupler%flow_bc_type(HYDRATE_GAS_EQUATION_INDEX) = NEUMANN_BC
    select type(selector => hydrate%gas_flux%dataset)
      class is(dataset_ascii_type)
        coupler%flow_aux_real_var(THREE_INTEGER,1:num_connections) = &
                                              hydrate%gas_flux%dataset%rarray(1)
        dof3 = PETSC_TRUE
      class is(dataset_gridded_hdf5_type)
        call PatchVerifyDatasetGriddedForFlux(selector,coupler,option)
        call PatchUpdateCouplerGridDataset(coupler,option,patch%grid,selector, &
                                           THREE_INTEGER)
        dof3 = PETSC_TRUE
      class default
        call PrintMsg(option,'hydrate%gas_flux%dataset')
        call DatasetUnknownClass(selector,option, &
                                 'PatchUpdateCouplerAuxVarsH')
    end select
  endif

  if (associated(hydrate%rate)) then
    select case(hydrate%rate%itype)
      case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS)
        call PatchScaleSourceSink(patch,coupler,hydrate%rate%isubtype,option)
    end select
  endif

  dof_count_global = 0
  dof_count_local = 0
  if (dof1) dof_count_local(1) = 1
  if (dof2) dof_count_local(2) = 1
  if (dof3) dof_count_local(3) = 1
  call MPI_Allreduce(dof_count_local,dof_count_global,THREE_INTEGER_MPI, &
                     MPI_INTEGER,MPI_SUM,option%mycomm,ierr)
  if (dof_count_global(1) > 0) dof1 = PETSC_TRUE
  if (dof_count_global(2) > 0) dof2 = PETSC_TRUE
  if (dof_count_global(3) > 0) dof3 = PETSC_TRUE
  ! need to check if these dofs are true on any process, because the
  ! boundary condition might be split up on 2 or more processes
  if (.not.dof1 .or. .not.dof2 .or. .not.dof3) then
    if (coupler%itype .ne. SRC_SINK_COUPLER_TYPE) then
      option%io_buffer = 'Error with HYDRATE_MODE phase boundary condition: &
                          &Missing dof.'
      call PrintErrMsg(option)
    endif
  endif

end subroutine PatchUpdateCouplerAuxVarsH

! ************************************************************************** !

subroutine PatchUpdateCouplerAuxVarsMPH(patch,coupler,option)
  !
  ! Updates flow auxiliary variables associated
  ! with a coupler for MPH_MODE
  !
  ! Author: Glenn Hammond
  ! Date: 11/26/07
  !

  use Option_module
  use Condition_module
  use Hydrostatic_module


  use Grid_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_HDF5_class

  implicit none

  type(patch_type) :: patch
  type(coupler_type), pointer :: coupler
  type(option_type) :: option

  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(flow_general_condition_type), pointer :: general
  class(dataset_common_hdf5_type), pointer :: dataset
  PetscBool :: update
  PetscBool :: dof1, dof2, dof3
  PetscReal :: temperature, p_sat
  PetscReal :: x(option%nflowdof)
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr

  PetscInt :: idof, num_connections,sum_connection
  PetscInt :: iconn, local_id, ghosted_id

  num_connections = coupler%connection_set%num_connections

  flow_condition => coupler%flow_condition

  if (associated(flow_condition%pressure)) then
    coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
                flow_condition%iphase
    select case(flow_condition%pressure%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        coupler%flow_aux_real_var(MPH_PRESSURE_DOF,1:num_connections) = &
                flow_condition%pressure%dataset%rarray(1)
      case(HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC)
        call HydrostaticUpdateCoupler(coupler,option,patch%grid)
 !  case(SATURATION_BC)
    end select
    select case(flow_condition%temperature%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        if (flow_condition%pressure%itype /= HYDROSTATIC_BC .or. &
           (flow_condition%pressure%itype == HYDROSTATIC_BC .and. &
           flow_condition%temperature%itype /= DIRICHLET_BC)) then
          coupler%flow_aux_real_var(MPH_TEMPERATURE_DOF,1:num_connections) = &
                  flow_condition%temperature%dataset%rarray(1)
        endif
    end select
    select case(flow_condition%concentration%itype)
      case(DIRICHLET_BC,ZERO_GRADIENT_BC)
        if (flow_condition%pressure%itype /= HYDROSTATIC_BC .or. &
           (flow_condition%pressure%itype == HYDROSTATIC_BC .and. &
           flow_condition%concentration%itype /= DIRICHLET_BC)) then
          coupler%flow_aux_real_var(MPH_CONCENTRATION_DOF,1:num_connections) = &
                  flow_condition%concentration%dataset%rarray(1)
        endif
    end select
  else
    select case(flow_condition%temperature%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        coupler%flow_aux_real_var(MPH_TEMPERATURE_DOF,1:num_connections) = &
                  flow_condition%temperature%dataset%rarray(1)
    end select
    select case(flow_condition%concentration%itype)
      case(DIRICHLET_BC,ZERO_GRADIENT_BC)
         coupler%flow_aux_real_var(MPH_CONCENTRATION_DOF,1:num_connections) = &
                  flow_condition%concentration%dataset%rarray(1)
    end select
  endif
  if (associated(flow_condition%rate)) then
    select case(flow_condition%rate%itype)
      case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS)
        call PatchScaleSourceSink(patch,coupler,flow_condition%rate%isubtype, &
                                  option)
    end select
  endif
  if (associated(flow_condition%saturation)) then
    call PatchUpdateCouplerSaturation(coupler,option,patch%grid, &
                                 patch%characteristic_curves_array, &
                                 patch%cc_id)
  endif

end subroutine PatchUpdateCouplerAuxVarsMPH

! ************************************************************************** !

subroutine PatchUpdateCouplerAuxVarsTH(patch,coupler,option)
  !
  ! Updates flow auxiliary variables associated
  ! with a coupler for TH_MODE
  !
  ! Author: Glenn Hammond
  ! Date: 11/26/07
  !

  use Option_module
  use Condition_module
  use Hydrostatic_module


  use Grid_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_HDF5_class
  use Dataset_Ascii_class
  use Dataset_module

  implicit none

  type(patch_type) :: patch
  type(coupler_type), pointer :: coupler
  type(option_type) :: option

  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(flow_general_condition_type), pointer :: general
  class(dataset_common_hdf5_type), pointer :: dataset
  PetscBool :: update
  PetscBool :: dof1, dof2, dof3
  PetscReal :: temperature, p_sat
  PetscReal :: x(option%nflowdof)
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr
  PetscBool :: apply_temp_cond
  PetscInt :: rate_scale_type

  PetscInt :: idof, num_connections,sum_connection
  PetscInt :: iconn, local_id, ghosted_id
  PetscInt :: iphase

  num_connections = coupler%connection_set%num_connections

  flow_condition => coupler%flow_condition

  if (associated(flow_condition%pressure)) then
    !geh: this is a fix for an Intel compiler bug. Not sure why Intel cannot
    !     access flow_condition%iphase directly....
    iphase = flow_condition%iphase
    coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = iphase
!    coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
!                                                        flow_condition%iphase
    select case(flow_condition%pressure%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC,SPILLOVER_BC, &
           DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC)
        select type(selector =>flow_condition%pressure%dataset)
          class is(dataset_ascii_type)
            coupler%flow_aux_real_var(TH_PRESSURE_DOF,1:num_connections) = &
              selector%rarray(1)
          class is(dataset_gridded_hdf5_type)
            call PatchUpdateCouplerGridDataset(coupler,option, &
                                               patch%grid,selector, &
                                               TH_PRESSURE_DOF)
          class is(dataset_common_hdf5_type)
            ! skip cell indexed datasets used in initial conditions
          class default
            call PrintMsg(option,'th%pressure%itype,DIRICHLET_BC')
            call DatasetUnknownClass(selector,option, &
                                     'PatchUpdateCouplerAuxVarsTH')
        end select
        select case(flow_condition%pressure%itype)
          case(DIRICHLET_CONDUCTANCE_BC)
            coupler%flow_aux_real_var(TH_CONDUCTANCE_DOF, &
                                      1:num_connections) = &
                                           flow_condition%pressure%aux_real(1)
        end select
      case(HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC)
        call HydrostaticUpdateCoupler(coupler,option,patch%grid)
      case(HET_DIRICHLET_BC,HET_HYDROSTATIC_SEEPAGE_BC, &
           HET_HYDROSTATIC_CONDUCTANCE_BC)
        call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                flow_condition%pressure%dataset,TH_PRESSURE_DOF,option)
        if (flow_condition%pressure%itype == &
            HET_HYDROSTATIC_CONDUCTANCE_BC) then
          coupler%flow_aux_real_var(TH_CONDUCTANCE_DOF,1:num_connections) = &
            flow_condition%pressure%aux_real(1)
        endif
      case default
        string = &
          GetSubConditionName(flow_condition%pressure%itype)
        option%io_buffer = &
          FlowConditionUnknownItype(flow_condition,'TH pressure',string)
        call PrintErrMsg(option)
    end select
    if (associated(flow_condition%temperature)) then
      select case(flow_condition%temperature%itype)
        case(DIRICHLET_BC,ZERO_GRADIENT_BC)
          select type(selector =>flow_condition%temperature%dataset)
            class is(dataset_ascii_type)
              if (flow_condition%pressure%itype /= HYDROSTATIC_BC .or. &
                 (flow_condition%pressure%itype == HYDROSTATIC_BC .and. &
                 flow_condition%temperature%itype /= DIRICHLET_BC)) then
                coupler%flow_aux_real_var(TH_TEMPERATURE_DOF, &
                                          1:num_connections) = &
                  selector%rarray(1)
              endif
            class is(dataset_gridded_hdf5_type)
              call PatchUpdateCouplerGridDataset(coupler,option, &
                                                 patch%grid,selector, &
                                                 TH_TEMPERATURE_DOF)
            class is(dataset_common_hdf5_type)
              ! skip cell indexed datasets used in initial conditions
            class default
              call PrintMsg(option,'th%temperature%itype,DIRICHLET_BC')
              call DatasetUnknownClass(selector,option, &
                                       'PatchUpdateCouplerAuxVarsTH')
          end select
        case (HET_DIRICHLET_BC)
          call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                  flow_condition%temperature%dataset, &
                  TH_TEMPERATURE_DOF,option)
        case default
          string = &
            GetSubConditionName(flow_condition%temperature%itype)
          option%io_buffer = &
            FlowConditionUnknownItype(flow_condition,'TH temperature',string)
          call PrintErrMsg(option)
      end select
    endif
    if (associated(flow_condition%energy_flux)) then
      coupler%flow_aux_real_var(TH_TEMPERATURE_DOF,1:num_connections) = &
        flow_condition%energy_flux%dataset%rarray(1)
    endif
  endif

  apply_temp_cond = PETSC_FALSE
  if (associated(flow_condition%temperature) .and. associated(flow_condition%pressure)) then
    if (flow_condition%pressure%itype /= HYDROSTATIC_BC) then
      apply_temp_cond = PETSC_TRUE
    else
      if (flow_condition%temperature%itype /= DIRICHLET_BC) then
        apply_temp_cond = PETSC_TRUE
      endif
    endif
  else
    apply_temp_cond = PETSC_TRUE
  endif

  if (associated(flow_condition%temperature) .and. apply_temp_cond) then
    select case(flow_condition%temperature%itype)
      case(DIRICHLET_BC,ZERO_GRADIENT_BC)
        select type(selector =>flow_condition%temperature%dataset)
          class is(dataset_ascii_type)
            coupler%flow_aux_real_var(TH_TEMPERATURE_DOF, &
                                      1:num_connections) = &
              selector%rarray(1)
          class is(dataset_gridded_hdf5_type)
            call PatchUpdateCouplerGridDataset(coupler,option, &
                                               patch%grid,selector, &
                                               TH_TEMPERATURE_DOF)
          class is(dataset_common_hdf5_type)
            ! skip cell indexed datasets used in initial conditions
          class default
            call PrintMsg(option,'th%pressure%itype,DIRICHLET_BC')
            call DatasetUnknownClass(selector,option, &
                                     'PatchUpdateCouplerAuxVarsTH')
        end select
      case (HET_DIRICHLET_BC)
        call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                flow_condition%temperature%dataset, &
                TH_TEMPERATURE_DOF,option)
      case default
        string = &
          GetSubConditionName(flow_condition%temperature%itype)
        option%io_buffer = &
          FlowConditionUnknownItype(flow_condition,'TH temperature',string)
        call PrintErrMsg(option)
    end select
  endif

  if (associated(flow_condition%energy_flux)) then
    select case(flow_condition%energy_flux%itype)
      case(NEUMANN_BC)
        select type(selector =>flow_condition%energy_flux%dataset)
          class is(dataset_ascii_type)
            coupler%flow_aux_real_var(TH_TEMPERATURE_DOF, &
                                      1:num_connections) = &
              selector%rarray(1)
          class is(dataset_gridded_hdf5_type)
            call PatchUpdateCouplerGridDataset(coupler,option, &
                                               patch%grid,selector, &
                                               TH_TEMPERATURE_DOF)
          class default
            call PrintMsg(option,'th%pressure%itype,NEUMANN_BC')
            call DatasetUnknownClass(selector,option, &
                                     'PatchUpdateCouplerAuxVarsTH')
        end select
      case default
        string = &
          GetSubConditionName(flow_condition%energy_flux%itype)
        option%io_buffer = &
          FlowConditionUnknownItype(flow_condition,'TH energy flux',string)
        call PrintErrMsg(option)
    end select
  endif

  !geh: we set this flag to ensure that we are not scaling mass and energy
  !     differently
  rate_scale_type = 0
  if (associated(flow_condition%rate)) then
    select case(flow_condition%rate%itype)
      case (HET_MASS_RATE_SS,HET_VOL_RATE_SS)
        call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                                            flow_condition%rate%dataset, &
                                            TH_PRESSURE_DOF,option)
      case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS)
        call PatchScaleSourceSink(patch,coupler,flow_condition%rate%isubtype, &
                                  option)
        rate_scale_type = flow_condition%rate%isubtype
      case(MASS_RATE_SS,VOLUMETRIC_RATE_SS)
      ! do nothing here
      case default
        string = &
          GetSubConditionName(flow_condition%rate%itype)
        option%io_buffer = &
          FlowConditionUnknownItype(flow_condition,'TH rate',string)
        call PrintErrMsg(option)
    end select
  endif
  if (associated(flow_condition%energy_rate)) then
    select case (flow_condition%energy_rate%itype)
      case (ENERGY_RATE_SS)
        !geh: this is pointless as %dataset%rarray(1) is reference in TH,
        !     not the flow_aux_real_var!
        coupler%flow_aux_real_var(TH_TEMPERATURE_DOF,1:num_connections) = &
                  flow_condition%energy_rate%dataset%rarray(1)
      case (SCALED_ENERGY_RATE_SS)
        if (rate_scale_type == 0) then
          call PatchScaleSourceSink(patch,coupler, &
                                    flow_condition%energy_rate%isubtype,option)
        else if (rate_scale_type == flow_condition%energy_rate%isubtype) then
          !geh: do nothing as it is taken care of later.
        else
          option%io_buffer = 'MASS and ENERGY scaling mismatch in ' // &
            'FLOW_CONDITION "' // trim(flow_condition%name) // '".'
          call PrintErrMsg(option)
        endif
        !geh: do nothing as the
      case (HET_ENERGY_RATE_SS)
        call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                flow_condition%energy_rate%dataset, &
                TH_TEMPERATURE_DOF,option)
      case default
        string = &
          GetSubConditionName(flow_condition%energy_rate%itype)
        option%io_buffer = &
          FlowConditionUnknownItype(flow_condition,'TH energy rate',string)
        call PrintErrMsg(option)
    end select
  endif
  if (associated(flow_condition%saturation)) then
    call PatchUpdateCouplerSaturation(coupler,option,patch%grid, &
                                 patch%characteristic_curves_array, &
                                 patch%cc_id)
  endif

end subroutine PatchUpdateCouplerAuxVarsTH

! ************************************************************************** !

subroutine PatchUpdateCouplerAuxVarsMIS(patch,coupler,option)
  !
  ! Updates flow auxiliary variables associated
  ! with a coupler for MIS_MODE
  !
  ! Author: Glenn Hammond
  ! Date: 11/26/07
  !

  use Option_module
  use Condition_module
  use Hydrostatic_module


  use Grid_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_HDF5_class

  implicit none

  type(patch_type) :: patch
  type(coupler_type), pointer :: coupler
  type(option_type) :: option

  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(flow_general_condition_type), pointer :: general
  class(dataset_common_hdf5_type), pointer :: dataset
  PetscBool :: update
  PetscBool :: dof1, dof2, dof3
  PetscReal :: temperature, p_sat
  PetscReal :: x(option%nflowdof)
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr

  PetscInt :: idof, num_connections,sum_connection
  PetscInt :: iconn, local_id, ghosted_id

  num_connections = coupler%connection_set%num_connections

  flow_condition => coupler%flow_condition
  if (associated(flow_condition%pressure)) then
    select case(flow_condition%pressure%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        coupler%flow_aux_real_var(MIS_PRESSURE_DOF, &
                                  1:num_connections) = &
          flow_condition%pressure%dataset%rarray(1)
      case(HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC)
        call HydrostaticUpdateCoupler(coupler,option,patch%grid)
   !  case(SATURATION_BC)
    end select
  endif
  if (associated(flow_condition%concentration)) then
    select case(flow_condition%concentration%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        if (associated(flow_condition%concentration%dataset)) then
          coupler%flow_aux_real_var(MIS_CONCENTRATION_DOF, &
                                    1:num_connections) = &
            flow_condition%concentration%dataset%rarray(1)
        endif
      case(HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC)
        call HydrostaticUpdateCoupler(coupler,option,patch%grid)
   !  case(SATURATION_BC)
    end select
  endif
  if (associated(flow_condition%rate)) then
    select case(flow_condition%rate%itype)
      case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS)
        call PatchScaleSourceSink(patch,coupler, &
                                  flow_condition%rate%isubtype,option)
    end select
  endif

end subroutine PatchUpdateCouplerAuxVarsMIS

! ************************************************************************** !

subroutine PatchUpdateCouplerAuxVarsRich(patch,coupler,option)
  !
  ! Updates flow auxiliary variables associated
  ! with a coupler for RICHARDS_MODE
  !
  ! Author: Glenn Hammond
  ! Date: 11/26/07
  !

  use Option_module
  use Condition_module
  use Hydrostatic_module


  use Grid_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_HDF5_class
  use Dataset_Ascii_class
  use Dataset_module

  implicit none

  type(patch_type) :: patch
  type(coupler_type), pointer :: coupler
  type(option_type) :: option

  type(flow_condition_type), pointer :: flow_condition
  class(dataset_common_hdf5_type), pointer :: dataset
  PetscBool :: update
  PetscBool :: dof1, dof2, dof3
  PetscReal :: temperature, p_sat
  PetscReal :: x(option%nflowdof)
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr

  PetscInt :: idof, num_connections,sum_connection
  PetscInt :: iconn, local_id, ghosted_id

  num_connections = coupler%connection_set%num_connections

  flow_condition => coupler%flow_condition
  if (associated(flow_condition%pressure)) then
    select case(flow_condition%pressure%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC,SURFACE_DIRICHLET, &
           SURFACE_SPILLOVER,DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC)
        select type(dataset => &
                    flow_condition%pressure%dataset)
          class is(dataset_ascii_type)
            coupler%flow_aux_real_var(RICHARDS_PRESSURE_DOF, &
                                      1:num_connections) = dataset%rarray(1)
          class is(dataset_gridded_hdf5_type)
            call PatchUpdateCouplerGridDataset(coupler,option, &
                                            patch%grid,dataset, &
                                            RICHARDS_PRESSURE_DOF)
          class is(dataset_common_hdf5_type)
            ! skip cell indexed datasets used in initial conditions
          class default
            call PrintMsg(option,'pressure%itype,DIRICHLET-type')
            call DatasetUnknownClass(dataset,option, &
                                     'PatchUpdateCouplerAuxVarsRich')
        end select
        select case(flow_condition%pressure%itype)
          case(DIRICHLET_CONDUCTANCE_BC)
            coupler%flow_aux_real_var(RICHARDS_CONDUCTANCE_DOF, &
                                      1:num_connections) = &
                                           flow_condition%pressure%aux_real(1)
        end select
      case(HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC)
        call HydrostaticUpdateCoupler(coupler,option,patch%grid)
   !  case(SATURATION_BC)
      case(HET_DIRICHLET_BC,HET_HYDROSTATIC_SEEPAGE_BC, &
           HET_HYDROSTATIC_CONDUCTANCE_BC)
        call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                flow_condition%pressure%dataset, &
                RICHARDS_PRESSURE_DOF,option)
        if (flow_condition%pressure%itype == &
            HET_HYDROSTATIC_CONDUCTANCE_BC) then
          coupler%flow_aux_real_var(RICHARDS_CONDUCTANCE_DOF, &
                                    1:num_connections) = &
            flow_condition%pressure%aux_real(1)
        endif
    end select
  endif
  if (associated(flow_condition%saturation)) then
    call PatchUpdateCouplerSaturation(coupler,option,patch%grid, &
                                 patch%characteristic_curves_array, &
                                 patch%cc_id)
  endif
  if (associated(flow_condition%rate)) then
    select case(flow_condition%rate%itype)
      case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS)
        call PatchScaleSourceSink(patch,coupler, &
                                  flow_condition%rate%isubtype,option)
      case (HET_VOL_RATE_SS,HET_MASS_RATE_SS)
        call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                flow_condition%rate%dataset, &
                RICHARDS_PRESSURE_DOF,option)
    end select
  endif

end subroutine PatchUpdateCouplerAuxVarsRich

! ************************************************************************** !

subroutine PatchUpdateCouplerAuxVarsZFlow(patch,coupler,option)
  !
  ! Updates flow auxiliary variables associated
  ! with a coupler for ZFLOW_MODE
  !
  ! Author: Glenn Hammond
  ! Date: 08/17/21
  !

  use Option_module
  use Condition_module
  use Hydrostatic_module
  use ZFlow_Aux_module

  use Grid_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_HDF5_class
  use Dataset_Ascii_class
  use Dataset_module

  implicit none

  type(patch_type) :: patch
  type(coupler_type), pointer :: coupler
  type(option_type) :: option

  type(flow_condition_type), pointer :: flow_condition
  class(dataset_common_hdf5_type), pointer :: dataset
  PetscInt :: water_index, water_aux_index, energy_index, solute_index
  PetscInt :: num_connections,sum_connection
  PetscInt :: iconn, local_id, ghosted_id
  PetscBool :: hydrostatic_update_called

  hydrostatic_update_called = PETSC_FALSE
  num_connections = coupler%connection_set%num_connections

  water_index = coupler%flow_aux_mapping(ZFLOW_COND_WATER_INDEX)
  water_aux_index = coupler%flow_aux_mapping(ZFLOW_COND_WATER_AUX_INDEX)
  energy_index = coupler%flow_aux_mapping(ZFLOW_COND_ENERGY_INDEX)
  solute_index = coupler%flow_aux_mapping(ZFLOW_COND_SOLUTE_INDEX)

  flow_condition => coupler%flow_condition
  if (associated(flow_condition%pressure)) then
    coupler%flow_bc_type(water_index) = flow_condition%pressure%itype
    select case(flow_condition%pressure%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC, &
           DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC)
        select type(dataset => &
                    flow_condition%pressure%dataset)
          class is(dataset_ascii_type)
            coupler%flow_aux_real_var(water_index,1:num_connections) = &
              dataset%rarray(1)
          class is(dataset_gridded_hdf5_type)
            call PatchUpdateCouplerGridDataset(coupler,option,patch%grid, &
                                               dataset,water_index)
          class is(dataset_common_hdf5_type)
            ! skip cell indexed datasets used in initial conditions
          class default
            call PrintMsg(option,'pressure%itype,DIRICHLET-type')
            call DatasetUnknownClass(dataset,option, &
                                     'PatchUpdateCouplerAuxVarsZFlow')
        end select
        select case(flow_condition%pressure%itype)
          case(DIRICHLET_CONDUCTANCE_BC)
            coupler%flow_aux_real_var(water_aux_index,1:num_connections) = &
                                           flow_condition%pressure%aux_real(1)
        end select
      case(HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC)
        hydrostatic_update_called = PETSC_TRUE
        call HydrostaticUpdateCoupler(coupler,option,patch%grid)
      case(HET_DIRICHLET_BC,HET_HYDROSTATIC_SEEPAGE_BC, &
           HET_HYDROSTATIC_CONDUCTANCE_BC)
        call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                flow_condition%pressure%dataset, &
                water_index,option)
        if (flow_condition%pressure%itype == &
            HET_HYDROSTATIC_CONDUCTANCE_BC) then
          coupler%flow_aux_real_var(water_aux_index, &
                                    1:num_connections) = &
            flow_condition%pressure%aux_real(1)
        endif
    end select
  endif

  if (associated(flow_condition%rate)) then
    coupler%flow_bc_type(water_index) = flow_condition%rate%itype
    coupler%flow_aux_real_var(water_index,1:num_connections) = &
          flow_condition%rate%dataset%rarray(1)
    select case(flow_condition%rate%itype)
      case(SCALED_VOLUMETRIC_RATE_SS)
        call PatchScaleSourceSink(patch,coupler, &
                                  flow_condition%rate%isubtype,option)
! not yet supported
!      case(HET_VOL_RATE_SS)
!        call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
!                flow_condition%rate%dataset, &
!                water_index,option)
      case default
    end select
  endif

  if (associated(flow_condition%concentration)) then
    coupler%flow_bc_type(solute_index) = flow_condition%concentration%itype
    select case(flow_condition%concentration%itype)
      case(DIRICHLET_BC,ZERO_GRADIENT_BC)
        if (.not.hydrostatic_update_called) then
          coupler%flow_aux_real_var(solute_index,1:num_connections) = &
            flow_condition%concentration%dataset%rarray(1)
        endif
      case(NEUMANN_BC)
        option%io_buffer = 'NEUMANN_BC not supported for ZFLOW concentration'
        call PrintErrMsg(option)
    end select
  endif

end subroutine PatchUpdateCouplerAuxVarsZFlow

! ************************************************************************** !

subroutine PatchUpdateCouplerAuxVarsPNF(patch,coupler,option)
  !
  ! Updates flow auxiliary variables associated
  ! with a coupler for PNF_MODE
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !

  use Option_module
  use Condition_module
  use PNF_Aux_module

  use Grid_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_HDF5_class
  use Dataset_Ascii_class
  use Dataset_module

  implicit none

  type(patch_type) :: patch
  type(coupler_type), pointer :: coupler
  type(option_type) :: option

  type(flow_condition_type), pointer :: flow_condition
  class(dataset_common_hdf5_type), pointer :: dataset
  PetscBool :: update
  PetscBool :: dof1, dof2, dof3
  PetscReal :: temperature, p_sat
  PetscReal :: x(option%nflowdof)
  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscErrorCode :: ierr

  PetscInt :: idof, num_connections,sum_connection
  PetscInt :: iconn, local_id, ghosted_id

  num_connections = coupler%connection_set%num_connections

  flow_condition => coupler%flow_condition
  if (associated(flow_condition%pressure)) then
    select case(flow_condition%pressure%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC)
        select type(dataset => &
                    flow_condition%pressure%dataset)
          class is(dataset_ascii_type)
            coupler%flow_aux_real_var(PNF_LIQUID_PRESSURE_DOF, &
                                      1:num_connections) = dataset%rarray(1)
          class is(dataset_gridded_hdf5_type)
            call PatchUpdateCouplerGridDataset(coupler,option, &
                                            patch%grid,dataset, &
                                            PNF_LIQUID_PRESSURE_DOF)
          class is(dataset_common_hdf5_type)
            ! skip cell indexed datasets used in initial conditions
          class default
            call PrintMsg(option,'pressure%itype,DIRICHLET-type')
            call DatasetUnknownClass(dataset,option, &
                                     'PatchUpdateCouplerAuxVarsPNF')
        end select
      case(HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC)
        option%io_buffer = 'HYDROSTATIC BCs not supported in PNF flow'
        call PrintErrMsg(option)
      case(HET_DIRICHLET_BC,HET_HYDROSTATIC_SEEPAGE_BC, &
           HET_HYDROSTATIC_CONDUCTANCE_BC)
        option%io_buffer = 'Heterogenenous BCs not supported in PNF flow'
        call PrintErrMsg(option)
    end select
    coupler%flow_bc_type(PNF_LIQUID_EQUATION_INDEX) = &
      flow_condition%pressure%itype
  endif
  if (associated(flow_condition%saturation)) then
    option%io_buffer = 'Saturation BCs not supported in PNF flow'
    call PrintErrMsg(option)
  endif
  if (associated(flow_condition%rate)) then
    select case(flow_condition%rate%itype)
      case(SCALED_VOLUMETRIC_RATE_SS)
        call PatchScaleSourceSink(patch,coupler, &
                                  flow_condition%rate%isubtype,option)
    end select
  endif

end subroutine PatchUpdateCouplerAuxVarsPNF

! ************************************************************************** !

subroutine PatchGetCouplerValueFromDataset(coupler,option,grid,dataset,iconn, &
                                           value)
  !
  ! Gets you an auxiliary variable from a dataset for a given connection.
  !
  ! Author: Glenn Hammond/Jennifer M. Frederick
  ! Date: 04/18/2018
  !

  use Option_module
  use Grid_module
  use Coupler_module
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_HDF5_class
  use Dataset_Ascii_class
  use Dataset_module

  implicit none

  type(coupler_type) :: coupler
  type(option_type) :: option
  type(grid_type) :: grid
  class(dataset_base_type) :: dataset
  PetscInt :: iconn
  PetscReal :: value

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscReal :: x
  PetscReal :: y
  PetscReal :: z
  PetscReal :: dist(-1:3)

  select type(dataset)
    class is(dataset_ascii_type)
      value = dataset%rarray(1)
    class is(dataset_gridded_hdf5_type)
      local_id = coupler%connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      x = grid%x(ghosted_id)
      y = grid%y(ghosted_id)
      z = grid%z(ghosted_id)
      if (associated(coupler%connection_set%dist)) then
        dist = coupler%connection_set%dist(:,iconn)
        x = x-dist(0)*dist(1)
        y = y-dist(0)*dist(2)
        z = z-dist(0)*dist(3)
      endif
      call DatasetGriddedHDF5InterpolateReal(dataset,x,y,z,value,option)
    class default
      call DatasetUnknownClass(dataset,option, &
                               'PatchGetCouplerValueFromDataset')
  end select

end subroutine PatchGetCouplerValueFromDataset

! ************************************************************************** !

subroutine PatchUpdateCouplerGridDataset(coupler,option,grid,dataset,dof)
  !
  ! Updates auxiliary variables from a gridded dataset.
  !
  ! Author: Glenn Hammond
  ! Date: 11/26/07
  !

  use Option_module
  use Grid_module
  use Coupler_module
  use Dataset_Gridded_HDF5_class

  implicit none

  type(coupler_type) :: coupler
  type(option_type) :: option
  type(grid_type) :: grid
  class(dataset_gridded_hdf5_type) :: dataset
  PetscInt :: dof

  PetscReal :: temp_real
  PetscInt :: iconn
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscReal :: x
  PetscReal :: y
  PetscReal :: z
  PetscReal :: dist(-1:3)

  do iconn = 1, coupler%connection_set%num_connections
    local_id = coupler%connection_set%id_dn(iconn)
    ghosted_id = grid%nL2G(local_id)
    x = grid%x(ghosted_id)
    y = grid%y(ghosted_id)
    z = grid%z(ghosted_id)
    if (associated(coupler%connection_set%dist)) then
      dist = coupler%connection_set%dist(:,iconn)
      x = x-dist(0)*dist(1)
      y = y-dist(0)*dist(2)
      z = z-dist(0)*dist(3)
    endif
    call DatasetGriddedHDF5InterpolateReal(dataset,x,y,z,temp_real,option)
    coupler%flow_aux_real_var(dof,iconn) = temp_real
  enddo

end subroutine PatchUpdateCouplerGridDataset

! ************************************************************************** !

subroutine PatchUpdateCouplerSaturation(coupler,option,grid,cc_array,cc_id)
  !
  ! Computes the pressures for a saturation
  ! initial/boundary condition
  !
  ! Author: Glenn Hammond
  ! Date: 06/07/11
  !
  ! Author: Heeho Park
  ! Modified Date: 10/05/20
  !

  use Option_module
  use Condition_module
  use Connection_module

  implicit none

  type(coupler_type) :: coupler
  type(option_type) :: option
  type(grid_type) :: grid
  type(characteristic_curves_ptr_type), intent(in) :: cc_array(:)
  class(characteristic_curves_type), pointer :: cc_ptr

  PetscInt :: cc_id(:)
  PetscInt :: local_id, ghosted_id, iconn
  PetscReal :: saturation
  PetscReal :: capillary_pressure
  PetscReal :: dpc_dsat
  PetscReal :: liquid_pressure

  type(flow_condition_type), pointer :: condition

  type(connection_set_type), pointer :: cur_connection_set

  condition => coupler%flow_condition

  if (option%iflowmode /= RICHARDS_MODE) then
    option%io_buffer = 'PatchUpdateCouplerSaturation is not set up for this flow mode.'
    call PrintErrMsg(option)
  endif

  ! in this case, the saturation is stored within concentration dataset
  saturation = condition%saturation%dataset%rarray(1)

  do iconn = 1, coupler%connection_set%num_connections
    local_id = coupler%connection_set%id_dn(iconn)
    ghosted_id = grid%nL2G(local_id)
    cc_ptr => cc_array(cc_id(ghosted_id))%ptr

    call cc_ptr%saturation_function%CapillaryPressure(saturation,&
                                           capillary_pressure,dpc_dsat,option)
    liquid_pressure = option%flow%reference_pressure - capillary_pressure
    coupler%flow_aux_real_var(1,iconn) = liquid_pressure
  enddo

end subroutine PatchUpdateCouplerSaturation

! ************************************************************************** !

subroutine PatchScaleSourceSink(patch,source_sink,iscale_type,option)
  !
  ! Scales select source/sinks based on perms*volume
  !
  ! Author: Glenn Hammond
  ! Date: 01/12/11
  !

#include "petsc/finclude/petscdmda.h"
  use petscdmda
  use Option_module
  use Field_module
  use Coupler_module
  use Connection_module
  use Condition_module
  use Grid_module
  use Material_Aux_module
  use ZFlow_Aux_module
  use Variables_module, only : PERMEABILITY_X, PERMEABILITY_Y, PERMEABILITY_Z

  implicit none

  type(patch_type) :: patch
  type(coupler_type) :: source_sink
  PetscInt :: iscale_type
  type(option_type) :: option

  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(connection_set_type), pointer :: cur_connection_set
  type(field_type), pointer :: field

  PetscReal, pointer :: vec_ptr(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id, neighbor_ghosted_id
  PetscInt :: iconn
  PetscReal :: scale, sum
  PetscInt :: icount, x_count, y_count, z_count
  PetscInt, parameter :: x_width = 1, y_width = 1, z_width = 0
  PetscInt :: ghosted_neighbors(27)
  PetscBool :: inactive_found
  type(material_auxvar_type), pointer :: material_auxvars(:)

  field => patch%field
  grid => patch%grid
  material_auxvars => patch%aux%Material%auxvars

  grid => patch%grid

  call VecZeroEntries(field%work,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)

  inactive_found = PETSC_FALSE
  cur_connection_set => source_sink%connection_set

  select case(iscale_type)
    case(SCALE_BY_VOLUME)
      do iconn = 1, cur_connection_set%num_connections
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        if (patch%imat(ghosted_id) <= 0) inactive_found = PETSC_TRUE
        vec_ptr(local_id) = vec_ptr(local_id) + &
          material_auxvars(ghosted_id)%volume
      enddo
    case(SCALE_BY_PERM)
      do iconn = 1, cur_connection_set%num_connections
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        if (patch%imat(ghosted_id) <= 0) inactive_found = PETSC_TRUE
        vec_ptr(local_id) = vec_ptr(local_id) + &
          material_auxvars(ghosted_id)%permeability(perm_xx_index) * &
          material_auxvars(ghosted_id)%volume
      enddo
    case(SCALE_BY_NEIGHBOR_PERM)
      do iconn = 1, cur_connection_set%num_connections
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        if (patch%imat(ghosted_id) <= 0) inactive_found = PETSC_TRUE
        !geh: kludge for 64-bit integers.
        call GridGetGhostedNeighbors(grid,ghosted_id,DMDA_STENCIL_STAR, &
                                    x_width,y_width,z_width, &
                                    x_count,y_count,z_count, &
                                    ghosted_neighbors,option)
        if (x_count + y_count + z_count == 0) then
          write(option%io_buffer,*) grid%nG2A(ghosted_id)
          option%io_buffer = 'Cell ' // trim(adjustl(option%io_buffer)) // &
            ' in FLOW_CONDITION "' // trim(source_sink%flow_condition%name) // &
            '" in SOURCE_SINK "' // trim(source_sink%name) // &
            '" has no neighbors, and therefore, NEIGHBOR_PERM cannot be used.'
          call PrintErrMsgByRank(option)
        endif
        ! ghosted neighbors is ordered first in x, then, y, then z
        icount = 0
        sum = 0.d0
        ! x-direction
        do while (icount < x_count)
          icount = icount + 1
          neighbor_ghosted_id = ghosted_neighbors(icount)
          if (patch%imat(neighbor_ghosted_id) <= 0) inactive_found = PETSC_TRUE
          sum = sum + &
                MaterialAuxVarGetValue(material_auxvars(neighbor_ghosted_id), &
                                       PERMEABILITY_X) * &
                grid%structured_grid%dy(neighbor_ghosted_id)* &
                grid%structured_grid%dz(neighbor_ghosted_id)
        enddo
        ! y-direction
        do while (icount < x_count + y_count)
          icount = icount + 1
          neighbor_ghosted_id = ghosted_neighbors(icount)
          if (patch%imat(neighbor_ghosted_id) <= 0) inactive_found = PETSC_TRUE
          sum = sum + &
                MaterialAuxVarGetValue(material_auxvars(neighbor_ghosted_id), &
                                       PERMEABILITY_Y) * &
                grid%structured_grid%dx(neighbor_ghosted_id)* &
                grid%structured_grid%dz(neighbor_ghosted_id)
        enddo
        ! z-direction
        do while (icount < x_count + y_count + z_count)
          icount = icount + 1
          neighbor_ghosted_id = ghosted_neighbors(icount)
          if (patch%imat(neighbor_ghosted_id) <= 0) inactive_found = PETSC_TRUE
          sum = sum + &
                MaterialAuxVarGetValue(material_auxvars(neighbor_ghosted_id), &
                                       PERMEABILITY_Z) * &
                grid%structured_grid%dx(neighbor_ghosted_id)* &
                grid%structured_grid%dy(neighbor_ghosted_id)
        enddo
        vec_ptr(local_id) = vec_ptr(local_id) + sum
      enddo
    case(0)
      option%io_buffer = 'Unknown scaling type in PatchScaleSourceSink ' // &
        'for FLOW_CONDITION "' // trim(source_sink%flow_condition%name) // '".'
      call PrintErrMsg(option)
  end select

  call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
  call VecNorm(field%work,NORM_1,scale,ierr);CHKERRQ(ierr)
  if (scale < 1.d-40) then
    option%io_buffer = 'Zero infinity norm in PatchScaleSourceSink for ' // &
      'FLOW_CONDITION "' // trim(source_sink%flow_condition%name) // '".'
    call PrintErrMsg(option)
  endif
  scale = 1.d0/scale
  call VecScale(field%work,scale,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(field%work,vec_ptr, ierr);CHKERRQ(ierr)
  do iconn = 1, cur_connection_set%num_connections
    local_id = cur_connection_set%id_dn(iconn)
    select case(option%iflowmode)
      !geh: This is a scaling factor that is stored that would be applied to
      !     all phases.
      case(RICHARDS_MODE,RICHARDS_TS_MODE,G_MODE,H_MODE,TH_MODE,TH_TS_MODE, &
           WF_MODE)
        source_sink%flow_aux_real_var(ONE_INTEGER,iconn) = &
          vec_ptr(local_id)
      case(ZFLOW_MODE)
        if (zflow_calc_adjoint .and. iscale_type /= SCALE_BY_VOLUME) then
          option%io_buffer = 'ZFLOW inversion must factor in derivative &
            &of scaled source/sink by permeability'
          call PrintErrMsg(option)
        endif
        source_sink%flow_aux_real_var( &
            source_sink%flow_aux_mapping(ZFLOW_COND_WATER_AUX_INDEX),&
            iconn) = vec_ptr(local_id)
      case(MPH_MODE,PNF_MODE)
        option%io_buffer = 'PatchScaleSourceSink not set up for flow mode'
        call PrintErrMsg(option)
    end select
  enddo
  call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)

  if (inactive_found) then
    option%io_buffer = 'Inactive cells found in source/sink "' // &
      trim(source_sink%name) // '" coupler region cells or neighboring cells.'
    call PrintErrMsgByRank(option)
  endif

end subroutine PatchScaleSourceSink

! ************************************************************************** !

subroutine PatchUpdateHetroCouplerAuxVars(patch,coupler,dataset_base, &
                                          isub_condition,option)
  !
  ! This subroutine updates aux vars for distributed copuler_type
  !
  ! Author: Gautam Bisht, LBL
  ! Date: 10/03/2012
  !

#include "petsc/finclude/petscdmda.h"
  use petscdmda
  use Option_module
  use Field_module
  use Coupler_module
  use Connection_module
  use Condition_module
  use Grid_module
  use Dataset_module
  use Dataset_Map_HDF5_class
  use Dataset_Base_class
  use Dataset_Ascii_class

  implicit none

  type(patch_type) :: patch
  type(coupler_type) :: coupler
  class(dataset_base_type), pointer :: dataset_base
  PetscInt :: isub_condition
  type(option_type) :: option

  type(connection_set_type), pointer :: cur_connection_set
  type(grid_type),pointer :: grid
  PetscErrorCode :: ierr
  PetscInt :: iconn
  PetscInt :: ghosted_id,local_id
  PetscInt,pointer ::cell_ids_nat(:)
  type(flow_sub_condition_type) :: flow_sub_condition

  class(dataset_map_hdf5_type), pointer :: dataset_map_hdf5
  class(dataset_ascii_type), pointer :: dataset_ascii

  grid => patch%grid

  if (isub_condition>option%nflowdof*option%nphase) then
    option%io_buffer='ERROR: PatchUpdateHetroCouplerAuxVars  '// &
      'isub_condition > option%nflowdof*option%nphase.'
    call PrintErrMsg(option)
  endif

  if (option%iflowmode/=RICHARDS_MODE .and. &
      option%iflowmode/=TH_MODE .and. &
      option%iflowmode/=TH_TS_MODE .and. &
      option%iflowmode/=RICHARDS_TS_MODE) then
    option%io_buffer='PatchUpdateHetroCouplerAuxVars only implemented '// &
      ' for RICHARDS or TH mode.'
    call PrintErrMsg(option)
  endif

  cur_connection_set => coupler%connection_set

  select type(selector=>dataset_base)
    class is(dataset_map_hdf5_type)
      dataset_map_hdf5 => selector

      ! If called for the first time, create the map
      if (dataset_map_hdf5%first_time) then
        allocate(cell_ids_nat(cur_connection_set%num_connections))
        do iconn=1,cur_connection_set%num_connections
          local_id = cur_connection_set%id_dn(iconn)
          ghosted_id = grid%nL2G(local_id)
          cell_ids_nat(iconn)=grid%nG2A(ghosted_id)
        enddo

        call PatchCreateFlowConditionDatasetMap(patch%grid,dataset_map_hdf5,&
                cell_ids_nat,cur_connection_set%num_connections,option)

        dataset_map_hdf5%first_time = PETSC_FALSE
        deallocate(cell_ids_nat)

      endif

      ! Save the data in the array
      do iconn=1,cur_connection_set%num_connections
        coupler%flow_aux_real_var(isub_condition,iconn) = &
          dataset_map_hdf5%rarray(dataset_map_hdf5%datatocell_ids(iconn))
      enddo

    class is(dataset_ascii_type)
      dataset_ascii => selector

      do iconn=1,cur_connection_set%num_connections
        coupler%flow_aux_real_var(isub_condition,iconn) = &
          dataset_ascii%rarray(1)
      enddo

    class default
      option%io_buffer = 'Incorrect dataset class for coupler "' // &
                         trim(coupler%name) // '".'
      call PrintMsg(option)
      call DatasetUnknownClass(selector,option, &
                               'PatchUpdateHetroCouplerAuxVars')
  end select

end subroutine PatchUpdateHetroCouplerAuxVars

! ************************************************************************** !

subroutine PatchCreateFlowConditionDatasetMap(grid,dataset_map_hdf5,cell_ids,ncells,option)
  !
  ! This routine creates dataset-map for flow condition
  !
  ! Author: Gautam Bisht, LBL
  ! Date: 10/26/12
  !
  use Grid_module
  use Dataset_Map_HDF5_class
  use Option_module

  implicit none

  type(grid_type) :: grid
  class(dataset_map_hdf5_type) :: dataset_map_hdf5
  type(option_type) :: option
  PetscInt,pointer :: cell_ids(:)
  PetscInt :: ncells

  PetscInt, allocatable :: int_array(:)
  PetscInt :: ghosted_id,local_id
  PetscInt :: ii,count
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr
  PetscInt :: max_id_loc, max_id_global
  PetscInt :: istart

  IS :: is_from, is_to
  Vec :: map_ids_1, map_ids_2,map_ids_3
  VecScatter ::vec_scatter
  PetscViewer :: viewer

  ! Step-1: Rearrange map dataset
  max_id_loc = maxval(dataset_map_hdf5%mapping(2,:))
  call MPI_Allreduce(max_id_loc,max_id_global,ONE_INTEGER,MPIU_INTEGER, &
                     MPI_MAX,option%mycomm,ierr)
  call VecCreateMPI(option%mycomm,dataset_map_hdf5%map_dims_local(2),&
                    PETSC_DETERMINE,map_ids_1,ierr);CHKERRQ(ierr)
  call VecCreateMPI(option%mycomm,PETSC_DECIDE,max_id_global,map_ids_2, &
                    ierr);CHKERRQ(ierr)
  call VecSet(map_ids_2,0.d0,ierr);CHKERRQ(ierr)

  istart = 0
  call MPI_Exscan(dataset_map_hdf5%map_dims_local(2), istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)

  allocate(int_array(dataset_map_hdf5%map_dims_local(2)))
  do ii=1,dataset_map_hdf5%map_dims_local(2)
    int_array(ii)=ii+istart
  enddo
  int_array=int_array-1

  call ISCreateBlock(option%mycomm,1,dataset_map_hdf5%map_dims_local(2), &
                     int_array,PETSC_COPY_VALUES,is_from,ierr);CHKERRQ(ierr)
  deallocate(int_array)

  allocate(int_array(dataset_map_hdf5%map_dims_local(2)))
  do ii=1,dataset_map_hdf5%map_dims_local(2)
    int_array(ii)=dataset_map_hdf5%mapping(2,ii)
  enddo
  int_array=int_array-1

  call ISCreateBlock(option%mycomm,1,dataset_map_hdf5%map_dims_local(2), &
                     int_array,PETSC_COPY_VALUES,is_to,ierr);CHKERRQ(ierr)
  deallocate(int_array)

  !call VecCreateSeq(PETSC_COMM_SELF,dataset_map%map_dims_global(2),map_ids_1,ierr)
  !call VecCreateSeq(PETSC_COMM_SELF,maxval(dataset_map%map(2,:)),map_ids_2,ierr)
  !call VecSet(map_ids_2,0,ierr)

  call VecScatterCreate(map_ids_1,is_from,map_ids_2,is_to,vec_scatter, &
                        ierr);CHKERRQ(ierr)
  call ISDestroy(is_from,ierr);CHKERRQ(ierr)
  call ISDestroy(is_to,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(map_ids_1,vec_ptr,ierr);CHKERRQ(ierr)
  do ii=1,dataset_map_hdf5%map_dims_local(2)
    vec_ptr(ii)=dataset_map_hdf5%mapping(1,ii)
  enddo
  call VecRestoreArrayF90(map_ids_1,vec_ptr,ierr);CHKERRQ(ierr)

  call VecScatterBegin(vec_scatter,map_ids_1,map_ids_2, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(vec_scatter,map_ids_1,map_ids_2, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(vec_scatter,ierr);CHKERRQ(ierr)

  ! Step-2: Get ids in map dataset for cells
  allocate(int_array(ncells))
  allocate(dataset_map_hdf5%cell_ids_local(ncells))
  int_array=cell_ids-1

  call ISCreateBlock(option%mycomm,1,ncells,int_array,PETSC_COPY_VALUES,is_from, &
                     ierr);CHKERRQ(ierr)

  istart = 0
  call MPI_Exscan(ncells, istart, ONE_INTEGER_MPI, &
                  MPIU_INTEGER, MPI_SUM, option%mycomm, ierr)

  do local_id=1,ncells
    int_array(local_id)=local_id+istart
  enddo
  int_array=int_array-1

  call ISCreateBlock(option%mycomm,1,ncells,int_array,PETSC_COPY_VALUES,is_to, &
                     ierr);CHKERRQ(ierr)
  deallocate(int_array)

  !call VecCreateSeq(PETSC_COMM_SELF,ncells,map_ids_3,ierr)
  call VecCreateMPI(option%mycomm,ncells,PETSC_DETERMINE,map_ids_3, &
                    ierr);CHKERRQ(ierr)

  call VecScatterCreate(map_ids_2,is_from,map_ids_3,is_to,vec_scatter, &
                        ierr);CHKERRQ(ierr)
  call ISDestroy(is_from,ierr);CHKERRQ(ierr)
  call ISDestroy(is_to,ierr);CHKERRQ(ierr)

  call VecScatterBegin(vec_scatter,map_ids_2,map_ids_3, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(vec_scatter,map_ids_2,map_ids_3, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterDestroy(vec_scatter,ierr);CHKERRQ(ierr)

  ! Step-3: Save the datatocell_ids
  allocate(dataset_map_hdf5%datatocell_ids(ncells))
  call VecGetArrayF90(map_ids_3,vec_ptr,ierr);CHKERRQ(ierr)
  do local_id=1,ncells
    dataset_map_hdf5%datatocell_ids(local_id) = int(vec_ptr(local_id))
  enddo
  call VecRestoreArrayF90(map_ids_3,vec_ptr,ierr);CHKERRQ(ierr)

  call VecDestroy(map_ids_1,ierr);CHKERRQ(ierr)
  call VecDestroy(map_ids_2,ierr);CHKERRQ(ierr)
  call VecDestroy(map_ids_3,ierr);CHKERRQ(ierr)

end subroutine PatchCreateFlowConditionDatasetMap

! ************************************************************************** !

subroutine PatchInitConstraints(patch,reaction_base,option)
  !
  ! Initializes constraint concentrations
  !
  ! Author: Glenn Hammond
  ! Date: 12/04/08
  !
  use Reaction_Base_module

  implicit none

  type(patch_type) :: patch
  class(reaction_base_type), pointer :: reaction_base
  type(option_type) :: option

  call PatchInitCouplerConstraints(patch%initial_condition_list, &
                                   reaction_base,option)

  call PatchInitCouplerConstraints(patch%boundary_condition_list, &
                                   reaction_base,option)

  call PatchInitCouplerConstraints(patch%source_sink_list, &
                                   reaction_base,option)

end subroutine PatchInitConstraints

! ************************************************************************** !

subroutine PatchInitCouplerConstraints(coupler_list,reaction_base,option)
  !
  ! Initializes constraint concentrations
  ! for a given coupler
  !
  ! Author: Glenn Hammond
  ! Date: 12/04/08
  !

  use Reaction_module
  use Reactive_Transport_Aux_module
  use NW_Transport_Aux_module
  use NWT_Equilibrium_module
  use Reaction_Aux_module
  use Global_Aux_module
  use Material_Aux_module
  use Transport_Constraint_Base_module
  use Transport_Constraint_NWT_module
  use Transport_Constraint_RT_module
  use Transport_Constraint_module
  use Dataset_Ascii_class

  use EOS_Water_module

  implicit none

  type(coupler_list_type), pointer :: coupler_list
  class(reaction_base_type), pointer :: reaction_base
  type(option_type) :: option

  type(reactive_transport_auxvar_type), pointer :: rt_auxvar
  type(nw_transport_auxvar_type), pointer :: nwt_auxvar
  type(global_auxvar_type), pointer :: global_auxvar
  type(material_auxvar_type), allocatable :: material_auxvar
  type(coupler_type), pointer :: cur_coupler
  class(tran_constraint_coupler_base_type), pointer :: cur_constraint_coupler
  class(reaction_rt_type), pointer :: reaction
  class(reaction_nw_type), pointer :: reaction_nw
  PetscReal :: dum1
  PetscErrorCode :: ierr

  nullify(global_auxvar)
  nullify(rt_auxvar)
  nullify(nwt_auxvar)

  allocate(material_auxvar)
  call MaterialAuxVarInit(material_auxvar,option)
  material_auxvar%porosity = option%flow%reference_porosity

  select type(r=>reaction_base)
    class is(reaction_rt_type)
      reaction => r
    class is(reaction_nw_type)
      reaction_nw => r
  end select

  cur_coupler => coupler_list%first
  do
    if (.not.associated(cur_coupler)) exit

    if (.not.associated(cur_coupler%tran_condition)) then
      option%io_buffer = 'Null transport condition found in coupler'
      if (len_trim(cur_coupler%name) > 1) then
        option%io_buffer = trim(option%io_buffer) // &
                           ' "' // trim(cur_coupler%name) // '"'
      endif
      call PrintErrMsg(option)
    endif

    cur_constraint_coupler => &
        cur_coupler%tran_condition%constraint_coupler_list
    do
      if (.not.associated(cur_constraint_coupler)) exit
      global_auxvar => cur_constraint_coupler%global_auxvar
      if (associated(cur_coupler%flow_condition)) then
        if (associated(cur_coupler%flow_condition%pressure)) then
          if (associated(cur_coupler%flow_condition%pressure%dataset)) then
            ! only use dataset value if the dataset is of type ascii
            select type(dataset=>cur_coupler%flow_condition%pressure%dataset)
              class is(dataset_ascii_type)
                global_auxvar%pres = dataset%rarray(1)
              class default
                ! otherwise, we don't know which pressure to use at this point,
                ! but we need to re-equilibrate at each cell
                cur_constraint_coupler%equilibrate_at_each_cell = PETSC_TRUE
                global_auxvar%pres = option%flow%reference_pressure
            end select
          else
            global_auxvar%pres = option%flow%reference_pressure
          endif
        else
          global_auxvar%pres = option%flow%reference_pressure
        endif
        if (associated(cur_coupler%flow_condition%temperature)) then
          if (associated(cur_coupler%flow_condition%temperature%dataset)) then
            ! only use dataset value if the dataset is of type ascii
            select type(dataset=>cur_coupler%flow_condition%temperature%dataset)
              class is(dataset_ascii_type)
                global_auxvar%temp = dataset%rarray(1)
              class default
                ! otherwise, we don't know which temperature to use at this
                ! point, but we need to re-equilibrate at each cell
                cur_constraint_coupler%equilibrate_at_each_cell = PETSC_TRUE
                global_auxvar%temp = option%flow%reference_temperature
            end select
          else
            global_auxvar%temp = option%flow%reference_temperature
          endif
        else
          global_auxvar%temp = option%flow%reference_temperature
        endif

        call EOSWaterDensity(global_auxvar%temp, &
                             global_auxvar%pres(1), &
                             global_auxvar%den_kg(1), &
                             dum1,ierr)
      else
        global_auxvar%pres = option%flow%reference_pressure
        global_auxvar%temp = option%flow%reference_temperature
        global_auxvar%den_kg(option%liquid_phase) = &
          option%flow%reference_density(option%liquid_phase)
      endif
      global_auxvar%sat = option%flow%reference_saturation

      if (option%transport%nphase > 1) then
        ! gas phase not considered explicitly on flow side
        global_auxvar%den_kg(option%gas_phase) = &
          option%flow%reference_density(option%gas_phase)
        global_auxvar%sat(option%gas_phase) = &
          1.d0 - global_auxvar%sat(option%liquid_phase)
      endif

      select type(constraint_coupler=>cur_constraint_coupler)
        class is(tran_constraint_coupler_rt_type)
          ! set this pointer for use below in CO2
          rt_auxvar => constraint_coupler%rt_auxvar
          call ReactionEquilibrateConstraint(rt_auxvar, &
                                  global_auxvar,material_auxvar,reaction, &
                         TranConstraintRTCast(constraint_coupler%constraint), &
                                  constraint_coupler%num_iterations, &
                                  PETSC_FALSE,option)
        class is(tran_constraint_coupler_nwt_type)
          nwt_auxvar => constraint_coupler%nwt_auxvar
          call NWTEquilibrateConstraint(reaction_nw, &
                        TranConstraintNWTCast(constraint_coupler%constraint), &
                                      nwt_auxvar,global_auxvar, &
                                      material_auxvar,option)
      end select
      ! update CO2 mole fraction for CO2 modes
      select case(option%iflowmode)
      ! TODO(jenn) Add error message saying you can't use NW Transport with MPH_MODE, etc. Here is not the best place, do it some where sooner in the set up.
        case(MPH_MODE)
          if ( (cur_coupler%flow_condition%iphase == 1) .and. &
               (associated(reaction)) ) then
            dum1 = RCO2MoleFraction(rt_auxvar,global_auxvar,reaction,option)
            cur_coupler%flow_condition%concentration%dataset%rarray(1) = dum1
            if (associated(cur_coupler%flow_aux_real_var)) then
              cur_coupler%flow_aux_real_var(MPH_CONCENTRATION_DOF,:) = dum1
            endif
          endif
      end select
      cur_constraint_coupler => cur_constraint_coupler%next
    enddo
    cur_coupler => cur_coupler%next
  enddo

  call MaterialAuxVarStrip(material_auxvar)
  deallocate(material_auxvar)

end subroutine PatchInitCouplerConstraints

! ************************************************************************** !

subroutine PatchUpdateUniformVelocity(patch,velocity,option)
  !
  ! Assigns uniform velocity in connection list
  ! darcy velocities
  !
  ! Author: Glenn Hammond
  ! Date: 02/20/08
  !

  use Option_module
  use Coupler_module
  use Condition_module
  use Connection_module

  implicit none

  type(patch_type), pointer :: patch
  PetscReal :: velocity(:)
  type(option_type), pointer :: option

  type(grid_type), pointer :: grid
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn, sum_connection, iphase
  PetscReal :: phase_velocity(3,option%transport%nphase)
  PetscReal :: vdarcy

  grid => patch%grid

  do iphase = 0, option%transport%nphase-1
    phase_velocity(1:3,iphase+1) = velocity(1+iphase*3:3+iphase*3)
  enddo

  ! Internal Flux Terms -----------------------------------
  cur_connection_set => grid%internal_connection_set_list%first
  sum_connection = 0
  do
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      do iphase = 1, option%transport%nphase
        vdarcy = dot_product(phase_velocity(:,iphase), &
                             cur_connection_set%dist(1:3,iconn))
        patch%internal_velocities(iphase,sum_connection) = vdarcy
      enddo
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
      do iphase = 1, option%transport%nphase
        vdarcy = dot_product(phase_velocity(:,iphase), &
                             cur_connection_set%dist(1:3,iconn))
        patch%boundary_velocities(iphase,sum_connection) = vdarcy
      enddo
    enddo
    boundary_condition => boundary_condition%next
  enddo

end subroutine PatchUpdateUniformVelocity

! ************************************************************************** !

subroutine PatchGetVariable1(patch,field,reaction_base,option, &
                             output_option,vec,ivar,isubvar,isubvar2)
  !
  ! PatchGetVariable: Extracts variables indexed by ivar and isubvar from a patch
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/08
  !

  use Grid_module
  use Option_module
  use Field_module

  use ERT_Aux_module
  use Mphase_Aux_module
  use TH_Aux_module
  use Richards_Aux_module
  use Reaction_Mineral_module
  use Reaction_Redox_module
  use Reaction_module
  use Reactive_Transport_Aux_module
  use Reaction_Surface_Complexation_Aux_module
  use General_Aux_module, only : general_fmw => fmw_comp, &
                                 GAS_STATE, LIQUID_STATE
  use WIPP_Flow_Aux_module, only : WIPPFloScalePerm
  use ZFlow_Aux_module
  use Output_Aux_module
  use Variables_module
  use Material_Aux_module
  use Reaction_Base_module

  implicit none

  type(option_type), pointer :: option
  class(reaction_base_type), pointer :: reaction_base
  type(output_option_type), pointer :: output_option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: isubvar2
  PetscInt :: iphase

  PetscInt :: local_id, ghosted_id
  type(grid_type), pointer :: grid
  type(material_auxvar_type), pointer :: material_auxvars(:)
  class(reaction_rt_type), pointer :: reaction
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscReal :: xmass, lnQKgas, ehfac, eh0, pe0, ph0, tk
  PetscReal :: tempreal
  PetscInt :: tempint, tempint2
  PetscInt :: irate, istate, irxn, ifo2, jcomp, comp_id
  PetscInt :: ivar_temp
  PetscErrorCode :: ierr


  grid => patch%grid
  material_auxvars => patch%aux%Material%auxvars
  reaction => ReactionCast(reaction_base)

  call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
  vec_ptr(:) = UNINITIALIZED_DOUBLE

  iphase = 1
  select case(ivar)
    case(TEMPERATURE,LIQUID_PRESSURE,GAS_PRESSURE,AIR_PRESSURE, &
         LIQUID_SATURATION,GAS_SATURATION,HYDRATE_SATURATION,ICE_SATURATION, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY, &
         LIQUID_DENSITY,GAS_DENSITY,GAS_DENSITY_MOL,LIQUID_VISCOSITY, &
         GAS_VISCOSITY,CAPILLARY_PRESSURE,LIQUID_DENSITY_MOL, &
         LIQUID_MOBILITY,GAS_MOBILITY,SC_FUGA_COEFF,ICE_DENSITY, &
         LIQUID_HEAD,VAPOR_PRESSURE,SATURATION_PRESSURE,DERIVATIVE, &
         MAXIMUM_PRESSURE,LIQUID_MASS_FRACTION,GAS_MASS_FRACTION, &
         SOLUTE_CONCENTRATION)

      if (associated(patch%aux%TH)) then
        select case(ivar)
          case(CAPILLARY_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%TH%auxvars(grid%nL2G(local_id))%pc
            enddo
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%auxvars(grid%nL2G(local_id))%temp
            enddo
          case(LIQUID_PRESSURE,MAXIMUM_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%auxvars(grid%nL2G(local_id))%pres(1)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%auxvars(grid%nL2G(local_id))%sat(1)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%auxvars(grid%nL2G(local_id))%den_kg(1)
            enddo
          case(GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY,GAS_VISCOSITY)
            call PatchUnsupportedVariable('TH','for gas phase',option)
          case(GAS_SATURATION)
            do local_id=1,grid%nlmax
              if (option%flow%th_freezing) then
                vec_ptr(local_id) = &
                  patch%aux%TH%auxvars(grid%nL2G(local_id))%ice%sat_gas
              else
                vec_ptr(local_id) = &
                  1.d0 - patch%aux%Global%auxvars(grid%nL2G(local_id))%sat(1)
              endif
            enddo
          case(ICE_SATURATION)
            if (option%flow%th_freezing) then
              do local_id=1,grid%nlmax
                vec_ptr(local_id) = &
                  patch%aux%TH%auxvars(grid%nL2G(local_id))%ice%sat_ice
              enddo
            else
              call PrintErrMsg(option,'ICE_SATURATION not supported by &
                                      &without freezing option TH')
            endif
          case(ICE_DENSITY)
            if (option%flow%th_freezing) then
              do local_id=1,grid%nlmax
                vec_ptr(local_id) = &
                  patch%aux%TH%auxvars(grid%nL2G(local_id))%ice%den_ice*FMWH2O
              enddo
            else
              call PrintErrMsg(option,'ICE_DENSITY not supported without &
                                      &freezing option in TH')
            endif
          case(LIQUID_VISCOSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                  patch%aux%TH%auxvars(grid%nL2G(local_id))%vis
            enddo
          case(LIQUID_MOBILITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                  patch%aux%TH%auxvars(grid%nL2G(local_id))%kvr
            enddo
          case(LIQUID_MOLE_FRACTION)
            call PatchUnsupportedVariable('TH','LIQUID_MOLE_FRACTION',option)
          case(LIQUID_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                  patch%aux%TH%auxvars(grid%nL2G(local_id))%u
            enddo
          case default
            call PatchUnsupportedVariable('TH',ivar,option)
        end select

      else if (associated(patch%aux%Richards)) then

        select case(ivar)
          case(TEMPERATURE)
            call PatchUnsupportedVariable('RICHARDS','TEMPERATURE',option)
          case(GAS_SATURATION)
            if (option%transport%nphase == 1) then
              call PatchUnsupportedVariable('RICHARDS','GAS_SATURATION',option)
            endif
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%auxvars(grid%nL2G(local_id))%sat(2)
            enddo
          case(ICE_SATURATION)
            call PatchUnsupportedVariable('RICHARDS','ICE_SATURATION',option)
          case(ICE_DENSITY)
            call PatchUnsupportedVariable('RICHARDS','ICE_DENSITY',option)
          case(GAS_DENSITY)
            call PatchUnsupportedVariable('RICHARDS','GAS_DENSITY',option)
          case(LIQUID_MOLE_FRACTION)
            call PatchUnsupportedVariable('RICHARDS','LIQUID_MOLE_FRACTION', &
                                          option)
          case(GAS_MOLE_FRACTION)
            call PatchUnsupportedVariable('RICHARDS','GAS_MOLE_FRACTION',option)
          case(LIQUID_ENERGY)
            call PatchUnsupportedVariable('RICHARDS','LIQUID_ENERGY',option)
          case(GAS_ENERGY)
            call PatchUnsupportedVariable('RICHARDS','GAS_ENERGY',option)
          case(LIQUID_VISCOSITY)
            do local_id = 1, grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = &
                  patch%aux%Richards%auxvars(ghosted_id)%kr / &
                  patch%aux%Richards%auxvars(ghosted_id)%kvr
            enddo
          case(GAS_VISCOSITY)
            call PatchUnsupportedVariable('RICHARDS','GAS_VISCOSITY',option)
          case(GAS_MOBILITY)
            call PatchUnsupportedVariable('RICHARDS','GAS_MOBILITY',option)
          case(LIQUID_PRESSURE,MAXIMUM_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%auxvars(grid%nL2G(local_id))%pres(1)
            enddo
          case(CAPILLARY_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Richards%auxvars(grid%nL2G(local_id))%pc
            enddo
          case(LIQUID_HEAD)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%auxvars(grid%nL2G(local_id))%pres(1)/ &
                EARTH_GRAVITY/ &
                patch%aux%Global%auxvars(grid%nL2G(local_id))%den_kg(1)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%auxvars(grid%nL2G(local_id))%sat(1)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%auxvars(grid%nL2G(local_id))%den_kg(1)
            enddo
          case(LIQUID_MOBILITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                  patch%aux%Richards%auxvars(grid%nL2G(local_id))%kvr
            enddo
          case default
            call PatchUnsupportedVariable('RICHARDS',ivar,option)
        end select

      else if (associated(patch%aux%PNF)) then

        select case(ivar)
          case(LIQUID_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%PNF%auxvars(grid%nL2G(local_id))%head
            enddo
          case default
            call PatchUnsupportedVariable('PNF',ivar,option)
          end select

      else if (associated(patch%aux%ZFlow)) then

        select case(ivar)
          case(LIQUID_PRESSURE,MAXIMUM_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%ZFlow%auxvars(ZERO_INTEGER,grid%nL2G(local_id))%pres
            enddo
          case(CAPILLARY_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%ZFlow%auxvars(ZERO_INTEGER,grid%nL2G(local_id))%pc
            enddo
          case(LIQUID_HEAD)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%ZFlow%auxvars(ZERO_INTEGER,grid%nL2G(local_id))%pres/ &
                EARTH_GRAVITY/ &
                zflow_density_kg
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%ZFlow%auxvars(ZERO_INTEGER,grid%nL2G(local_id))%sat
            enddo
          case(SOLUTE_CONCENTRATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%ZFlow%auxvars(ZERO_INTEGER,grid%nL2G(local_id))%conc
            enddo
          case(DERIVATIVE)
            select case(isubvar)
              case(ZFLOW_LIQ_SAT_WRT_LIQ_PRES)
                do local_id=1,grid%nlmax
                  vec_ptr(local_id) = &
                    patch%aux%ZFlow%auxvars(ZERO_INTEGER, &
                                            grid%nL2G(local_id))%dsat_dp
                enddo
              case default
                call PatchUnsupportedVariable('ZFLOW','DERIVATIVE', &
                                              isubvar,option)
            end select
          case default
            call PatchUnsupportedVariable('ZFLOW',ivar,option)
        end select

      else if (associated(patch%aux%Mphase)) then

        select case(ivar)

          case(MAXIMUM_PRESSURE)
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = &
                  maxval(patch%aux%Global%auxvars(ghosted_id)%pres(1:2))
            enddo
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global% &
                auxvars(grid%nL2G(local_id))%temp
            enddo
          case(LIQUID_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global% &
                auxvars(grid%nL2G(local_id))%pres(1)
            enddo
          case(GAS_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global% &
                auxvars(grid%nL2G(local_id))%pres(2)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global% &
                auxvars(grid%nL2G(local_id))%sat(1)
            enddo
          case(LIQUID_DENSITY)
!geh: CO2 Mass Balance Fix (change to #if 0 to scale the mixture density by the water mole fraction)
#if 1
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global% &
                auxvars(grid%nL2G(local_id))%den_kg(1)
            enddo
#else
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%auxvars(grid%nL2G(local_id))%den(1) * &
                patch%aux%Mphase%auxvars(grid%nL2G(local_id))% &
                  auxvar_elem(0)%xmol(1) * FMWH2O
            enddo
#endif
          case(LIQUID_VISCOSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase% &
                auxvars(grid%nL2G(local_id))%auxvar_elem(0)%vis(1)
            enddo
          case(LIQUID_MOBILITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase% &
                auxvars(grid%nL2G(local_id))%auxvar_elem(0)%kvr(1)
            enddo
          case(GAS_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global% &
                auxvars(grid%nL2G(local_id))%sat(2)
            enddo
          case(GAS_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase% &
                auxvars(grid%nL2G(local_id))%auxvar_elem(0)%xmol(2+isubvar)
            enddo
          case(GAS_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase% &
                auxvars(grid%nL2G(local_id))%auxvar_elem(0)%u(2)
            enddo
          case(GAS_VISCOSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase% &
                auxvars(grid%nL2G(local_id))%auxvar_elem(0)%vis(2)
            enddo
          case(GAS_MOBILITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase% &
                auxvars(grid%nL2G(local_id))%auxvar_elem(0)%kvr(2)
            enddo
          case(GAS_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global% &
                auxvars(grid%nL2G(local_id))%den_kg(2)
            enddo
          case(GAS_DENSITY_MOL)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%auxvars(grid%nL2G(local_id))%den(2)
            enddo
          case(SC_FUGA_COEFF)
            if (.not.associated(patch%aux%Global%auxvars(1)%fugacoeff) .and. &
                OptionPrintToScreen(option))then
               print *,'ERROR: fugacoeff not allocated for ', &
                       option%iflowmode, 1
            endif
            do local_id=1,grid%nlmax
             vec_ptr(local_id) = patch%aux%Global%&
                auxvars(grid%nL2G(local_id))%fugacoeff(1)
            enddo
          case(LIQUID_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%&
                auxvars(grid%nL2G(local_id))%auxvar_elem(0)%xmol(isubvar)
            enddo
          case(LIQUID_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Mphase%&
                auxvars(grid%nL2G(local_id))%auxvar_elem(0)%u(1)
            enddo
          case default
            call PatchUnsupportedVariable('MPHASE',ivar,option)
        end select

      else if (associated(patch%aux%General)) then

        select case(ivar)
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%General%auxvars(ZERO_INTEGER, &
                                          grid%nL2G(local_id))%temp
            enddo
          case(MAXIMUM_PRESSURE)
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = &
                  maxval(patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                           pres(option%liquid_phase:option%gas_phase))
            enddo
          case(LIQUID_PRESSURE)
            if (output_option%filter_non_state_variables) then
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                if (patch%aux%Global%auxvars(ghosted_id)%istate /= &
                    GAS_STATE) then
                  vec_ptr(local_id) = &
                    patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      pres(option%liquid_phase)
                else
                  vec_ptr(local_id) = 0.d0
                endif
              enddo
            else
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                vec_ptr(local_id) = &
                  patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                    pres(option%liquid_phase)
              enddo
            endif
          case(GAS_PRESSURE)
            if (output_option%filter_non_state_variables) then
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                if (patch%aux%Global%auxvars(ghosted_id)%istate /= &
                    LIQUID_STATE) then
                  vec_ptr(local_id) = &
                    patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      pres(option%gas_phase)
                else
                  vec_ptr(local_id) = 0.d0
                endif
              enddo
            else
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                vec_ptr(local_id) = &
                  patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                    pres(option%gas_phase)
              enddo
            endif
          case(AIR_PRESSURE)
            if (output_option%filter_non_state_variables) then
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                if (patch%aux%Global%auxvars(ghosted_id)%istate /= &
                    LIQUID_STATE) then
                  vec_ptr(local_id) = &
                    patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      pres(option%air_pressure_id)
                else
                  vec_ptr(local_id) = 0.d0
                endif
              enddo
            else
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                vec_ptr(local_id) = &
                  patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                    pres(option%air_pressure_id)
              enddo
            endif
          case(CAPILLARY_PRESSURE)
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = &
                patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                  pres(option%capillary_pressure_id)
            enddo
          case(VAPOR_PRESSURE)
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = &
                patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                  pres(option%vapor_pressure_id)
            enddo
          case(SATURATION_PRESSURE)
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = &
                patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                  pres(option%saturation_pressure_id)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%sat(option%liquid_phase)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%den_kg(option%liquid_phase)
            enddo
          case(LIQUID_DENSITY_MOL)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%den(option%liquid_phase)
            enddo
          case(LIQUID_ENERGY)
            if (isubvar == ZERO_INTEGER) then
              do local_id=1,grid%nlmax
                vec_ptr(local_id) = patch%aux%General%auxvars(ZERO_INTEGER, &
                    grid%nL2G(local_id))%U(option%liquid_phase)
              enddo
            else
              do local_id=1,grid%nlmax
                vec_ptr(local_id) = patch%aux%General%auxvars(ZERO_INTEGER, &
                      grid%nL2G(local_id))%U(option%liquid_phase) * &
                    patch%aux%General%auxvars(ZERO_INTEGER, &
                      grid%nL2G(local_id))%den(option%liquid_phase)
              enddo
            endif
          case(LIQUID_MOLE_FRACTION,LIQUID_MASS_FRACTION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%xmol(isubvar,option%liquid_phase)
            enddo
            if (ivar == LIQUID_MASS_FRACTION) then
              tempint = isubvar
              tempint2 = tempint+1
              if (tempint2 > 2) tempint2 = 1
              vec_ptr(:) = vec_ptr(:)*general_fmw(tempint) / &
                           (vec_ptr(:)*general_fmw(tempint) + &
                            (1.d0-vec_ptr(:))*general_fmw(tempint2))
            endif
          case(LIQUID_MOBILITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%mobility(option%liquid_phase)
            enddo
          case(LIQUID_VISCOSITY)
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = &
                patch%aux%General%auxvars(ZERO_INTEGER, &
                  ghosted_id)%kr(option%liquid_phase) / &
                patch%aux%General%auxvars(ZERO_INTEGER, &
                  ghosted_id)%mobility(option%liquid_phase)
            enddo
          case(GAS_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%sat(option%gas_phase)
            enddo
          case(GAS_ENERGY)
            if (isubvar == ZERO_INTEGER) then
              do local_id=1,grid%nlmax
                vec_ptr(local_id) = patch%aux%General%auxvars(ZERO_INTEGER, &
                    grid%nL2G(local_id))%U(option%gas_phase)
              enddo
            else
              do local_id=1,grid%nlmax
                vec_ptr(local_id) = patch%aux%General%auxvars(ZERO_INTEGER, &
                      grid%nL2G(local_id))%U(option%gas_phase) * &
                    patch%aux%General%auxvars(ZERO_INTEGER, &
                      grid%nL2G(local_id))%den(option%gas_phase)
              enddo
            endif
          case(GAS_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%den_kg(option%gas_phase)
            enddo
          case(GAS_DENSITY_MOL)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%den(option%gas_phase)
            enddo
          case(GAS_MOLE_FRACTION,GAS_MASS_FRACTION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%xmol(isubvar,option%gas_phase)
            enddo
            if (ivar == GAS_MASS_FRACTION) then
              tempint = isubvar
              tempint2 = tempint+1
              if (tempint2 > 2) tempint2 = 1
              vec_ptr(:) = vec_ptr(:)*general_fmw(tempint) / &
                           (vec_ptr(:)*general_fmw(tempint) + &
                            (1.d0-vec_ptr(:))*general_fmw(tempint2))
            endif
          case(GAS_MOBILITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%General%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%mobility(option%gas_phase)
            enddo
          case(GAS_VISCOSITY)
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = &
                patch%aux%General%auxvars(ZERO_INTEGER, &
                  ghosted_id)%kr(option%gas_phase) / &
                patch%aux%General%auxvars(ZERO_INTEGER, &
                  ghosted_id)%mobility(option%gas_phase)
            enddo
          case default
            call PatchUnsupportedVariable('GENERAL',ivar,option)
        end select

      else if (associated(patch%aux%Hydrate)) then
        select case(ivar)
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                                          grid%nL2G(local_id))%temp
            enddo
          case(MAXIMUM_PRESSURE)
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = &
                  maxval(patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                           pres(option%liquid_phase:option%gas_phase))
            enddo
          case(LIQUID_PRESSURE)
            if (output_option%filter_non_state_variables) then
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                if (patch%aux%Global%auxvars(ghosted_id)%istate /= &
                    GAS_STATE) then
                  vec_ptr(local_id) = &
                    patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      pres(option%liquid_phase)
                else
                  vec_ptr(local_id) = 0.d0
                endif
              enddo
            else
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                vec_ptr(local_id) = &
                  patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                    pres(option%liquid_phase)
              enddo
            endif
          case(GAS_PRESSURE)
            if (output_option%filter_non_state_variables) then
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                if (patch%aux%Global%auxvars(ghosted_id)%istate /= &
                    LIQUID_STATE) then
                  vec_ptr(local_id) = &
                    patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      pres(option%gas_phase)
                else
                  vec_ptr(local_id) = 0.d0
                endif
              enddo
            else
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                vec_ptr(local_id) = &
                  patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                    pres(option%gas_phase)
              enddo
            endif
          case(AIR_PRESSURE)
            if (output_option%filter_non_state_variables) then
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                if (patch%aux%Global%auxvars(ghosted_id)%istate /= &
                    LIQUID_STATE) then
                  vec_ptr(local_id) = &
                    patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      pres(option%air_pressure_id)
                else
                  vec_ptr(local_id) = 0.d0
                endif
              enddo
            else
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                vec_ptr(local_id) = &
                  patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                    pres(option%air_pressure_id)
              enddo
            endif
          case(CAPILLARY_PRESSURE)
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = &
                patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                  pres(option%capillary_pressure_id)
            enddo
          case(VAPOR_PRESSURE)
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = &
                patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                  pres(option%vapor_pressure_id)
            enddo
          case(SATURATION_PRESSURE)
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = &
                patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                  pres(option%saturation_pressure_id)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%sat(option%liquid_phase)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%den_kg(option%liquid_phase)
            enddo
          case(LIQUID_DENSITY_MOL)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%den(option%liquid_phase)
            enddo
          case(LIQUID_ENERGY)
            if (isubvar == ZERO_INTEGER) then
              do local_id=1,grid%nlmax
                vec_ptr(local_id) = patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                    grid%nL2G(local_id))%U(option%liquid_phase)
              enddo
            else
              do local_id=1,grid%nlmax
                vec_ptr(local_id) = patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                      grid%nL2G(local_id))%U(option%liquid_phase) * &
                    patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                      grid%nL2G(local_id))%den(option%liquid_phase)
              enddo
            endif
          case(LIQUID_MOLE_FRACTION,LIQUID_MASS_FRACTION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%xmol(isubvar,option%liquid_phase)
            enddo
            if (ivar == LIQUID_MASS_FRACTION) then
              tempint = isubvar
              tempint2 = tempint+1
              if (tempint2 > 2) tempint2 = 1
              !MAN: need to put in proper conversion here
              vec_ptr(:) = vec_ptr(:)*general_fmw(tempint) / &
                           (vec_ptr(:)*general_fmw(tempint) + &
                            (1.d0-vec_ptr(:))*general_fmw(tempint2))
            endif
          case(LIQUID_MOBILITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%mobility(option%liquid_phase)
            enddo
          case(LIQUID_VISCOSITY)
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = &
                patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                  ghosted_id)%kr(option%liquid_phase) / &
                patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                  ghosted_id)%mobility(option%liquid_phase)
            enddo
          case(GAS_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%sat(option%gas_phase)
            enddo
          case(GAS_ENERGY)
            if (isubvar == ZERO_INTEGER) then
              do local_id=1,grid%nlmax
                vec_ptr(local_id) = patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                    grid%nL2G(local_id))%U(option%gas_phase)
              enddo
            else
              do local_id=1,grid%nlmax
                vec_ptr(local_id) = patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                      grid%nL2G(local_id))%U(option%gas_phase) * &
                    patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                      grid%nL2G(local_id))%den(option%gas_phase)
              enddo
            endif
          case(GAS_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%den_kg(option%gas_phase)
            enddo
          case(GAS_DENSITY_MOL)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%den(option%gas_phase)
            enddo
          case(GAS_MOLE_FRACTION,GAS_MASS_FRACTION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%xmol(isubvar,option%gas_phase)
            enddo
            if (ivar == GAS_MASS_FRACTION) then
              tempint = isubvar
              tempint2 = tempint+1
              if (tempint2 > 2) tempint2 = 1
              !MAN: need to put in proper conversion here
              vec_ptr(:) = vec_ptr(:)*general_fmw(tempint) / &
                           (vec_ptr(:)*general_fmw(tempint) + &
                            (1.d0-vec_ptr(:))*general_fmw(tempint2))
            endif
          case(GAS_MOBILITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%mobility(option%gas_phase)
            enddo
          case(HYDRATE_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%sat(option%hydrate_phase)
            enddo
          case(ICE_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%sat(option%ice_phase)
            enddo
          case(GAS_VISCOSITY)
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = &
                patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                  ghosted_id)%kr(option%gas_phase) / &
                patch%aux%Hydrate%auxvars(ZERO_INTEGER, &
                  ghosted_id)%mobility(option%gas_phase)
            enddo
          case default
            call PatchUnsupportedVariable('HYDRATE',ivar,option)
        end select

      else if (associated(patch%aux%WIPPFlo)) then
        select case(ivar)
          case(MAXIMUM_PRESSURE)
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = &
                  maxval(patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                           pres(option%liquid_phase:option%gas_phase))
            enddo
          case(LIQUID_PRESSURE)
            if (output_option%filter_non_state_variables) then
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                if (patch%aux%Global%auxvars(ghosted_id)%istate /= &
                    GAS_STATE) then
                  vec_ptr(local_id) = &
                    patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                      pres(option%liquid_phase)
                else
                  vec_ptr(local_id) = 0.d0
                endif
              enddo
            else
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                vec_ptr(local_id) = &
                  patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                    pres(option%liquid_phase)
              enddo
            endif
          case(GAS_PRESSURE)
            if (output_option%filter_non_state_variables) then
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                if (patch%aux%Global%auxvars(ghosted_id)%istate /= &
                    LIQUID_STATE) then
                  vec_ptr(local_id) = &
                    patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                      pres(option%gas_phase)
                else
                  vec_ptr(local_id) = 0.d0
                endif
              enddo
            else
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                vec_ptr(local_id) = &
                  patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                    pres(option%gas_phase)
              enddo
            endif
          case(CAPILLARY_PRESSURE)
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = &
                patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                  pres(option%capillary_pressure_id)
            enddo
          case(SATURATION_PRESSURE)
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = &
                patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                  pres(option%saturation_pressure_id)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%sat(option%liquid_phase)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%den_kg(option%liquid_phase)
            enddo
          case(LIQUID_DENSITY_MOL)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%den(option%liquid_phase)
            enddo
          case(LIQUID_MOBILITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%mobility(option%liquid_phase)
            enddo
          case(LIQUID_VISCOSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%mu(option%liquid_phase)
            enddo
          case(GAS_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%sat(option%gas_phase)
            enddo
          case(GAS_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%den_kg(option%gas_phase)
            enddo
          case(GAS_DENSITY_MOL)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%den(option%gas_phase)
            enddo
          case(GAS_MOBILITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%mobility(option%gas_phase)
            enddo
          case(GAS_VISCOSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER, &
                  grid%nL2G(local_id))%mu(option%gas_phase)
            enddo
          case default
            call PatchUnsupportedVariable('WIPP_FLOW',ivar,option)
        end select
      else ! null flow mode
        select case(ivar)
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Global%auxvars( &
                  grid%nL2G(local_id))%sat(option%liquid_phase)
            enddo
          case default
            call PatchUnsupportedVariable('NULL_FLOW',ivar,option)
        end select
      endif

    ! NUCLEAR_WASTE_TRANSPORT:
    case(TOTAL_BULK_CONC,AQUEOUS_EQ_CONC,MNRL_EQ_CONC,SORB_EQ_CONC, &
         MNRL_VOLUME_FRACTION)
      select case(ivar)
        case(TOTAL_BULK_CONC)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
            patch%aux%NWT%auxvars(grid%nL2G(local_id))%total_bulk_conc(isubvar)
            if ( (patch%aux%NWT%truncate_output) .and. &
                 ((vec_ptr(local_id) < 1.d-99) .and. (vec_ptr(local_id) > 0.d0)) ) then
              vec_ptr(local_id) = 1.d-99
            endif
          enddo
        case(AQUEOUS_EQ_CONC)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%NWT%auxvars(grid%nL2G(local_id))%aqueous_eq_conc(isubvar)
            if ( (patch%aux%NWT%truncate_output) .and. &
                 ((vec_ptr(local_id) < 1.d-99) .and. (vec_ptr(local_id) > 0.d0)) ) then
              vec_ptr(local_id) = 1.d-99
            endif
          enddo
        case(MNRL_EQ_CONC)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%NWT%auxvars(grid%nL2G(local_id))%mnrl_eq_conc(isubvar)
            if ( (patch%aux%NWT%truncate_output) .and. &
                 ((vec_ptr(local_id) < 1.d-99) .and. (vec_ptr(local_id) > 0.d0)) ) then
              vec_ptr(local_id) = 1.d-99
            endif
          enddo
        case(SORB_EQ_CONC)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%NWT%auxvars(grid%nL2G(local_id))%sorb_eq_conc(isubvar)
            if ( (patch%aux%NWT%truncate_output) .and. &
                 ((vec_ptr(local_id) < 1.d-99) .and. (vec_ptr(local_id) > 0.d0)) ) then
              vec_ptr(local_id) = 1.d-99
            endif
          enddo
        case(MNRL_VOLUME_FRACTION)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%NWT%auxvars(grid%nL2G(local_id))%mnrl_vol_frac(isubvar)
            if ( (patch%aux%NWT%truncate_output) .and. &
                 ((vec_ptr(local_id) < 1.d-99) .and. (vec_ptr(local_id) > 0.d0)) ) then
              vec_ptr(local_id) = 1.d-99
            endif
          enddo
      end select


    case(PH,PE,EH,O2,PRIMARY_MOLALITY,PRIMARY_MOLARITY,SECONDARY_MOLALITY, &
         SECONDARY_MOLARITY,TOTAL_MOLALITY,TOTAL_MOLARITY, &
         MINERAL_RATE,MINERAL_VOLUME_FRACTION,MINERAL_SATURATION_INDEX, &
         MINERAL_SURFACE_AREA, &
         SURFACE_CMPLX,SURFACE_CMPLX_FREE,SURFACE_SITE_DENSITY, &
         KIN_SURFACE_CMPLX,KIN_SURFACE_CMPLX_FREE, PRIMARY_ACTIVITY_COEF, &
         SECONDARY_ACTIVITY_COEF,PRIMARY_KD,TOTAL_SORBED,TOTAL_SORBED_MOBILE, &
         COLLOID_MOBILE,COLLOID_IMMOBILE,AGE,TOTAL_BULK,IMMOBILE_SPECIES, &
         GAS_CONCENTRATION,REACTION_AUXILIARY)

      select case(ivar)
        case(PH)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            call RRedoxCalcpH(patch%aux%RT%auxvars(ghosted_id), &
                              patch%aux%Global%auxvars(ghosted_id), &
                              reaction,ph0,option)
            vec_ptr(local_id) = ph0
          enddo
        case(EH)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            call RRedoxCalcEhpe(patch%aux%RT%auxvars(ghosted_id), &
                                  patch%aux%Global%auxvars(ghosted_id), &
                                  reaction,eh0,pe0,option)
            vec_ptr(local_id) = eh0
          enddo
        case(PE)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            call RRedoxCalcEhpe(patch%aux%RT%auxvars(ghosted_id), &
                                  patch%aux%Global%auxvars(ghosted_id), &
                                  reaction,eh0,pe0,option)
            vec_ptr(local_id) = pe0
          enddo
        case(O2)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            call RRedoxCalcLnFO2(patch%aux%RT%auxvars(ghosted_id), &
                                 patch%aux%Global%auxvars(ghosted_id), &
                                 reaction,lnQKgas,option)
            vec_ptr(local_id) = lnQKgas
          enddo
        case(PRIMARY_MOLALITY)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(grid%nL2G(local_id))%pri_molal(isubvar)
          enddo
        case(PRIMARY_MOLARITY)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            if (associated(patch%aux%Global%auxvars(ghosted_id)%xmass)) then
              xmass = patch%aux%Global%auxvars(ghosted_id)%xmass(iphase)
            else
              xmass = 1.d0
            endif
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar) * xmass * &
              (patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase)/1000.d0)
          enddo
        case(SECONDARY_MOLALITY)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(ghosted_id)%sec_molal(isubvar)
          enddo
        case(SECONDARY_MOLARITY)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            if (associated(patch%aux%Global%auxvars(ghosted_id)%xmass)) then
              xmass = patch%aux%Global%auxvars(ghosted_id)%xmass(iphase)
            else
              xmass = 1.d0
            endif
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(ghosted_id)%sec_molal(isubvar) * xmass * &
              (patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase)/1000.d0)
          enddo
        case(TOTAL_MOLALITY)
          do local_id=1,grid%nlmax
            ghosted_id =grid%nL2G(local_id)
            if (associated(patch%aux%Global%auxvars(ghosted_id)%xmass)) then
              xmass = patch%aux%Global%auxvars(ghosted_id)%xmass(iphase)
            else
              xmass = 1.d0
            endif
            if (patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase) > 0.d0) then
              vec_ptr(local_id) = &
                patch%aux%RT%auxvars(ghosted_id)%total(isubvar,iphase) / &
                xmass / &
                (patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase)/1000.d0)
            else
              vec_ptr(local_id) = 0.d0
            endif
          enddo
        case(TOTAL_MOLARITY)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(grid%nL2G(local_id))%total(isubvar,iphase)
          enddo
        case(TOTAL_BULK) ! mol/m^3 bulk
          ! add in total molarity and convert to mol/m^3 bulk
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(ghosted_id)%total(isubvar,iphase) * &
              patch%aux%Material%auxvars(ghosted_id)%porosity * &
                                                             ! mol/L -> mol/m^3
              patch%aux%Global%auxvars(ghosted_id)%sat(iphase) * 1.d-3
          enddo
          ! add in total sorbed.  already in mol/m^3 bulk
          if (patch%reaction%nsorb > 0) then
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              if (patch%reaction%surface_complexation%neqsrfcplxrxn > 0) then
                vec_ptr(local_id) = vec_ptr(local_id) + &
                  patch%aux%RT%auxvars(ghosted_id)%total_sorb_eq(isubvar)
              endif
              if (patch%reaction%surface_complexation%nkinmrsrfcplxrxn > 0) then
                do irxn = 1, &
                   patch%reaction%surface_complexation%nkinmrsrfcplxrxn
                  do irate = 1, &
                     patch%reaction%surface_complexation%kinmr_nrate(irxn)
                    vec_ptr(local_id) = vec_ptr(local_id) + &
                      patch%aux%RT%auxvars(ghosted_id)% &
                        kinmr_total_sorb(isubvar,irate,irxn)
                  enddo
                enddo
              endif
            enddo
          endif
        case(GAS_CONCENTRATION)
          iphase = 2
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            if (patch%aux%Global%auxvars(ghosted_id)%sat(iphase) > 0.d0) then
              vec_ptr(local_id) = &
                patch%aux%RT%auxvars(ghosted_id)%gas_pp(isubvar)
            else
              vec_ptr(local_id) = 0.d0
            endif
          enddo
        case(MINERAL_VOLUME_FRACTION)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(grid%nL2G(local_id))%mnrl_volfrac(isubvar)
          enddo
        case(MINERAL_SURFACE_AREA)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(grid%nL2G(local_id))%mnrl_area(isubvar)
          enddo
        case(MINERAL_RATE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(grid%nL2G(local_id))%mnrl_rate(isubvar)
          enddo
        case(MINERAL_SATURATION_INDEX)
          do local_id = 1, grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            vec_ptr(local_id) = &
              RMineralSaturationIndex(isubvar, &
                                      patch%aux%RT%auxvars(ghosted_id), &
                                      patch%aux%Global%auxvars(ghosted_id), &
                                      reaction,option)
          enddo
        case(IMMOBILE_SPECIES)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(grid%nL2G(local_id))%immobile(isubvar)
          enddo
        case(SURFACE_CMPLX)
          if (associated(patch%aux%RT%auxvars(1)%eqsrfcplx_conc)) then
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%RT%auxvars(grid%nL2G(local_id))% &
                                    eqsrfcplx_conc(isubvar)
            enddo
          else
            vec_ptr = UNINITIALIZED_DOUBLE
          endif
        case(SURFACE_SITE_DENSITY)
          tempreal = &
            reaction%surface_complexation%srfcplxrxn_site_density(isubvar)
          select case(reaction%surface_complexation% &
                        srfcplxrxn_surf_type(isubvar))
            case(ROCK_SURFACE)
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                vec_ptr(local_id) = tempreal* &
                        material_auxvars(ghosted_id)%soil_particle_density * &
                        (1.d0-material_auxvars(ghosted_id)%porosity)
              enddo
            case(MINERAL_SURFACE)
              tempint = &
                reaction%surface_complexation%srfcplxrxn_to_surf(isubvar)
              do local_id=1,grid%nlmax
                vec_ptr(local_id) = tempreal* &
                                    patch%aux%RT%auxvars(grid%nL2G(local_id))% &
                                      mnrl_volfrac(tempint)
              enddo
            case(COLLOID_SURFACE)
                option%io_buffer = 'Printing of surface site density for ' // &
                                     'colloidal surfaces not implemented.'
                call PrintErrMsg(option)
            case(NULL_SURFACE)
              do local_id=1,grid%nlmax
                vec_ptr(local_id) = tempreal
              enddo
          end select
        case(SURFACE_CMPLX_FREE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%aux%RT%auxvars(grid%nL2G(local_id))% &
              srfcplxrxn_free_site_conc(isubvar)
          enddo
        case(KIN_SURFACE_CMPLX)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%aux%RT%auxvars(grid%nL2G(local_id))% &
              kinsrfcplx_conc(isubvar,1)
          enddo
        case(KIN_SURFACE_CMPLX_FREE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = patch%aux%RT%auxvars(grid%nL2G(local_id))% &
              kinsrfcplx_free_site_conc(isubvar)
          enddo
        case(PRIMARY_ACTIVITY_COEF)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(ghosted_id)%pri_act_coef(isubvar)
          enddo
        case(SECONDARY_ACTIVITY_COEF)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(ghosted_id)%sec_act_coef(isubvar)
          enddo
        case(PRIMARY_KD)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            call ReactionComputeKd(isubvar,vec_ptr(local_id), &
                                   patch%aux%RT%auxvars(ghosted_id), &
                                   patch%aux%Global%auxvars(ghosted_id), &
                                   patch%aux%Material%auxvars(ghosted_id), &
                                   patch%reaction,option)
          enddo
        case(TOTAL_SORBED)
          if (patch%reaction%nsorb > 0) then
            if (patch%reaction%neqsorb > 0) then
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                vec_ptr(local_id) = &
                  patch%aux%RT%auxvars(ghosted_id)%total_sorb_eq(isubvar)
              enddo
            endif
            if (patch%reaction%surface_complexation%nkinmrsrfcplxrxn > 0) then
              do local_id=1,grid%nlmax
                ghosted_id = grid%nL2G(local_id)
                vec_ptr(local_id) = 0.d0
                do irxn = 1, &
                  patch%reaction%surface_complexation%nkinmrsrfcplxrxn
                  do irate = 1, &
                    patch%reaction%surface_complexation%kinmr_nrate(irxn)
                    vec_ptr(local_id) = vec_ptr(local_id) + &
                      patch%aux%RT%auxvars(ghosted_id)% &
                        kinmr_total_sorb(isubvar,irate,irxn)
                  enddo
                enddo
              enddo
            endif
          endif
        case(TOTAL_SORBED_MOBILE)
          if (patch%reaction%nsorb > 0 .and. patch%reaction%ncollcomp > 0) then
            do local_id=1,grid%nlmax
              ghosted_id = grid%nL2G(local_id)
              vec_ptr(local_id) = patch%aux%RT%auxvars(ghosted_id)%colloid% &
                total_eq_mob(isubvar)
            enddo
          endif
        case(COLLOID_MOBILE)
          if (patch%reaction%print_tot_conc_type == TOTAL_MOLALITY) then
            do local_id=1,grid%nlmax
              ghosted_id =grid%nL2G(local_id)
              if (patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase) > &
                  0.d0) then
                vec_ptr(local_id) = &
                  patch%aux%RT%auxvars(ghosted_id)% &
                    colloid%conc_mob(isubvar) / &
                  (patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase)/1000.d0)
              else
                vec_ptr(local_id) = 0.d0
              endif
            enddo
          else
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%RT%auxvars(grid%nL2G(local_id))% &
                                    colloid%conc_mob(isubvar)
            enddo
          endif
        case(COLLOID_IMMOBILE)
          if (patch%reaction%print_tot_conc_type == TOTAL_MOLALITY) then
            do local_id=1,grid%nlmax
              ghosted_id =grid%nL2G(local_id)
              if (patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase) > &
                  0.d0) then
                vec_ptr(local_id) = &
                  patch%aux%RT%auxvars(ghosted_id)% &
                    colloid%conc_imb(isubvar) / &
                  (patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase)/1000.d0)
              else
                vec_ptr(local_id) = 0.d0
              endif
            enddo
          else
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%RT%auxvars(grid%nL2G(local_id))% &
                                    colloid%conc_imb(isubvar)
            enddo
          endif
        case(AGE)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            if (patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar) > &
                0.d0) then
              vec_ptr(local_id) = &
                patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar) / &
                patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar2) / &
                output_option%tconv
            endif
          enddo
        case(REACTION_AUXILIARY)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%RT%auxvars(grid%nL2G(local_id))%auxiliary_data(isubvar)
          enddo
        case default
          call PatchUnsupportedVariable('REACTIVE_TRANSPORT',ivar,option)
      end select
    case(STATE,PHASE)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = &
          patch%aux%Global%auxvars(grid%nL2G(local_id))%istate
      enddo
    case(POROSITY,BASE_POROSITY,INITIAL_POROSITY, &
         VOLUME,TORTUOSITY,SOIL_COMPRESSIBILITY, &
         SOIL_REFERENCE_PRESSURE,ELECTRICAL_CONDUCTIVITY)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = &
          MaterialAuxVarGetValue(material_auxvars(grid%nL2G(local_id)),ivar)
      enddo
    case(PERMEABILITY,PERMEABILITY_X,PERMEABILITY_Y, PERMEABILITY_Z, &
         PERMEABILITY_XY,PERMEABILITY_XZ,PERMEABILITY_YZ, &
         GAS_PERMEABILITY,GAS_PERMEABILITY_X,GAS_PERMEABILITY_Y, &
         GAS_PERMEABILITY_Z)
      ivar_temp = ivar
      ! only liquid permeabilities in x, y, z are stored.
      select case(ivar)
        case(PERMEABILITY,GAS_PERMEABILITY,GAS_PERMEABILITY_X)
          ivar_temp = PERMEABILITY_X
        case(GAS_PERMEABILITY_Y)
          ivar_temp = PERMEABILITY_Y
        case(GAS_PERMEABILITY_Z)
          ivar_temp = PERMEABILITY_Z
      end select
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = &
          MaterialAuxVarGetValue(material_auxvars(grid%nL2G(local_id)), &
                                 ivar_temp)
      enddo
      select case(option%iflowmode)
        case(WF_MODE)
          do local_id=1,grid%nlmax
            ghosted_id = grid%nL2G(local_id)
            call WIPPFloScalePerm(patch%aux%WIPPFlo%auxvars(ZERO_INTEGER, &
                                                            ghosted_id), &
                                  material_auxvars(ghosted_id), &
                                  vec_ptr(local_id),ivar)
          enddo
      end select
    case(LIQUID_RELATIVE_PERMEABILITY)
      select case(option%iflowmode)
        case(RICHARDS_MODE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%Richards%auxvars(grid%nL2G(local_id))%kr
          enddo
        case(ZFLOW_MODE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%ZFlow%auxvars(ZERO_INTEGER,grid%nL2G(local_id))%kr
          enddo
        case(TH_MODE,TH_TS_MODE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%TH%auxvars(grid%nL2G(local_id))%kvr * &
              patch%aux%TH%auxvars(grid%nL2G(local_id))%vis
          enddo
        case(G_MODE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%General%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                kr(option%liquid_phase)
          enddo
        case(H_MODE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%Hydrate%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                kr(option%liquid_phase)
          enddo
        case(WF_MODE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                kr(option%liquid_phase)
          enddo
        case default
          option%io_buffer = 'Output of liquid phase relative permeability &
            &not supported for current flow mode.'
      end select
    case(GAS_RELATIVE_PERMEABILITY)
      select case(option%iflowmode)
        case(G_MODE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%General%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                kr(option%gas_phase)
          enddo
        case(H_MODE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%Hydrate%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                kr(option%gas_phase)
          enddo
        case(WF_MODE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                kr(option%gas_phase)
          enddo
        case default
          option%io_buffer = 'Output of gas phase relative permeability &
            &not supported for current flow mode.'
      end select
    case(MATERIAL_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = &
          patch%imat_internal_to_external(abs(patch%imat(grid%nL2G(local_id))))
      enddo
    case(FRACTURE)
      vec_ptr = 0.d0
      do local_id=1,grid%nlmax
        ghosted_id = grid%nL2G(local_id)
        if (associated(material_auxvars(ghosted_id)%fracture)) then
          if (material_auxvars(ghosted_id)%fracture%fracture_is_on) then
            vec_ptr(local_id) = 1.d0
          endif
        endif
      enddo
    case(ELECTRICAL_POTENTIAL)
      call ERTAuxCheckElectrodeBounds(size(patch%aux%ERT%auxvars(1)% &
                                      potential),isubvar,isubvar2,option)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = &
          patch%aux%ERT%auxvars(grid%nL2G(local_id))%potential(isubvar)
      enddo
    case(ELECTRICAL_POTENTIAL_DIPOLE)
      call ERTAuxCheckElectrodeBounds(size(patch%aux%ERT%auxvars(1)% &
                                      potential),isubvar,isubvar2,option)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = &
          patch%aux%ERT%auxvars(grid%nL2G(local_id))%potential(isubvar) - &
          patch%aux%ERT%auxvars(grid%nL2G(local_id))%potential(isubvar2)
      enddo
    case(ELECTRICAL_JACOBIAN)
      call ERTAuxCheckElectrodeBounds(size(patch%aux%ERT%auxvars(1)% &
                                      jacobian),isubvar,isubvar2,option)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = &
          patch%aux%ERT%auxvars(grid%nL2G(local_id))%jacobian(isubvar)
      enddo
    case(PROCESS_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = option%myrank
      enddo
    case(NATURAL_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = grid%nG2A(grid%nL2G(local_id))
      enddo
    case(SALINITY)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = &
          patch%aux%Global%auxvars(grid%nL2G(local_id))%m_nacl(ONE_INTEGER)
      enddo
    case(RESIDUAL)
      call VecRestoreArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
      call VecStrideGather(field%flow_r,isubvar-1,vec,INSERT_VALUES,ierr)
      call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
    case(X_COORDINATE)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = grid%x(grid%nL2G(local_id))
      enddo
    case(Y_COORDINATE)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = grid%y(grid%nL2G(local_id))
      enddo
    case(Z_COORDINATE)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = grid%z(grid%nL2G(local_id))
      enddo
    case(K_ORTHOGONALITY_ERROR)
      if (grid%itype == IMPLICIT_UNSTRUCTURED_GRID) then
        call PatchGetKOrthogonalityError(grid, material_auxvars, vec_ptr)
        !vec_ptr(:) = UNINITIALIZED_DOUBLE
      endif
    ! PM Material Transform
    case(SMECTITE)
      if (associated(patch%aux%MTransform)) then
        select case(ivar)
          case(SMECTITE)
            do local_id=1,grid%nlmax
              if (associated(patch%aux%MTransform% &
                  auxvars(grid%nL2G(local_id))%il_aux)) then
                vec_ptr(local_id) = &
                  patch%aux%MTransform%auxvars(grid%nL2G(local_id))%il_aux%fs
              endif
            enddo
        end select
      endif
    case default
      call PatchUnsupportedVariable(ivar,option)
  end select

  call VecRestoreArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine PatchGetVariable1

! ************************************************************************** !

function PatchGetVariableValueAtCell(patch,field,reaction_base,option, &
                                     output_option,ghosted_id, &
                                     ivar,isubvar,isubvar2)
  !
  ! Returns variables indexed by ivar,
  ! isubvar, local id from Reactive Transport type
  !
  ! Author: Glenn Hammond
  ! Date: 02/11/08
  !

  use Grid_module
  use Option_module
  use Field_module

  use Mphase_Aux_module
  use TH_Aux_module
  use Richards_Aux_module
  use Reactive_Transport_Aux_module
  use Reaction_Mineral_module
  use Reaction_Redox_module
  use Reaction_module
  use Reaction_Mineral_Aux_module
  use Reaction_Surface_Complexation_Aux_module
  use Output_Aux_module
  use ERT_Aux_module
  use Variables_module
  use General_Aux_module, only : general_fmw => fmw_comp, &
                                 GAS_STATE, LIQUID_STATE
  use WIPP_Flow_Aux_module, only : WIPPFloScalePerm
  use ZFlow_Aux_module
  use Material_Aux_module

  implicit none

  PetscReal :: PatchGetVariableValueAtCell
  type(option_type), pointer :: option
  class(reaction_base_type), pointer :: reaction_base
  type(output_option_type), pointer :: output_option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(material_auxvar_type), pointer :: material_auxvars(:)
  class(reaction_rt_type), pointer :: reaction
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: tempint, tempint2
  PetscInt :: isubvar2
  PetscInt :: iphase
  PetscInt :: ghosted_id
  PetscInt :: local_id
  PetscInt :: ivar_temp

  PetscReal :: value, xmass, lnQKgas, tk, ehfac, eh0, pe0, ph0
  PetscInt :: irate, istate, irxn, ifo2, jcomp, comp_id
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr2(:)
  PetscErrorCode :: ierr

  grid => patch%grid
  material_auxvars => patch%aux%Material%auxvars
  reaction => ReactionCast(reaction_base)

  value = UNINITIALIZED_DOUBLE

  ! inactive grid cell
  if (patch%imat(ghosted_id) <= 0) then
    PatchGetVariableValueAtCell = 0.d0
    return
  endif

  iphase = 1
  xmass = 1.d0
  if (associated(patch%aux%Global%auxvars(ghosted_id)%xmass)) &
    xmass = patch%aux%Global%auxvars(ghosted_id)%xmass(iphase)

  select case(ivar)
    case(TEMPERATURE,LIQUID_PRESSURE,GAS_PRESSURE, &
         LIQUID_SATURATION,GAS_SATURATION,ICE_SATURATION,HYDRATE_SATURATION, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY, &
         LIQUID_DENSITY,GAS_DENSITY,GAS_DENSITY_MOL,LIQUID_VISCOSITY, &
         GAS_VISCOSITY,AIR_PRESSURE,CAPILLARY_PRESSURE, &
         LIQUID_MOBILITY,GAS_MOBILITY,SC_FUGA_COEFF,ICE_DENSITY, &
         SECONDARY_TEMPERATURE,LIQUID_DENSITY_MOL, &
         LIQUID_HEAD,VAPOR_PRESSURE,SATURATION_PRESSURE,MAXIMUM_PRESSURE, &
         LIQUID_MASS_FRACTION,GAS_MASS_FRACTION,SOLUTE_CONCENTRATION)

      if (associated(patch%aux%TH)) then
        select case(ivar)
          case(TEMPERATURE)
            value = patch%aux%Global%auxvars(ghosted_id)%temp
          case(LIQUID_PRESSURE,MAXIMUM_PRESSURE)
            value = patch%aux%Global%auxvars(ghosted_id)%pres(1)
          case(LIQUID_SATURATION)
            value = patch%aux%Global%auxvars(ghosted_id)%sat(1)
          case(LIQUID_DENSITY)
            value = patch%aux%Global%auxvars(ghosted_id)%den_kg(1)
          case(LIQUID_VISCOSITY)
            value = patch%aux%TH%auxvars(ghosted_id)%vis
          case(LIQUID_MOBILITY)
            value = patch%aux%TH%auxvars(ghosted_id)%kvr
          case(GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY)
            call PatchUnsupportedVariable('TH','GAS_MOLE_FRACTION',option)
          case(GAS_SATURATION)
            if (option%flow%th_freezing) then
              value = patch%aux%TH%auxvars(ghosted_id)%ice%sat_gas
            else
              value = 0.d0
            endif
          case(CAPILLARY_PRESSURE)
            value = patch%aux%TH%auxvars(ghosted_id)%pc
          case(ICE_SATURATION)
            if (option%flow%th_freezing) then
              value = patch%aux%TH%auxvars(ghosted_id)%ice%sat_ice
            endif
          case(ICE_DENSITY)
            if (option%flow%th_freezing) then
              value = patch%aux%TH%auxvars(ghosted_id)%ice%den_ice*FMWH2O
            endif
          case(LIQUID_MOLE_FRACTION)
            call PatchUnsupportedVariable('TH','LIQUID_MOLE_FRACTION',option)
          case(LIQUID_ENERGY)
            value = patch%aux%TH%auxvars(ghosted_id)%u
          case(SECONDARY_TEMPERATURE)
            local_id = grid%nG2L(ghosted_id)
            value = patch%aux%SC_heat%sec_heat_vars(local_id)%sec_temp(isubvar)
          case default
            call PatchUnsupportedVariable('TH',ivar,option)
        end select
      else if (associated(patch%aux%Richards)) then
        select case(ivar)
          case(TEMPERATURE)
            call PatchUnsupportedVariable('RICHARDS','TEMPERATURE',option)
          case(GAS_SATURATION)
            if (option%transport%nphase == 1) then
              call PatchUnsupportedVariable('RICHARDS','GAS_SATURATION',option)
            endif
            value = patch%aux%Global%auxvars(ghosted_id)%sat(2)
          case(GAS_DENSITY)
            call PatchUnsupportedVariable('RICHARDS','GAS_DENSITY',option)
          case(LIQUID_MOLE_FRACTION)
            call PatchUnsupportedVariable('RICHARDS','LIQUID_MOLE_FRACTION', &
                                          option)
          case(GAS_MOLE_FRACTION)
            call PatchUnsupportedVariable('RICHARDS','GAS_MOLE_FRACTION',option)
          case(LIQUID_ENERGY)
            call PatchUnsupportedVariable('RICHARDS','LIQUID_ENERGY',option)
          case(GAS_ENERGY)
            call PatchUnsupportedVariable('RICHARDS','GAS_ENERGY',option)
          case(LIQUID_PRESSURE,MAXIMUM_PRESSURE)
            value = patch%aux%Global%auxvars(ghosted_id)%pres(1)
          case(CAPILLARY_PRESSURE)
            value = patch%aux%Richards%auxvars(ghosted_id)%pc
          case(LIQUID_HEAD)
            value = patch%aux%Global%auxvars(ghosted_id)%pres(1)/ &
                    EARTH_GRAVITY/ &
                    patch%aux%Global%auxvars(ghosted_id)%den_kg(1)
          case(LIQUID_SATURATION)
            value = patch%aux%Global%auxvars(ghosted_id)%sat(1)
          case(LIQUID_DENSITY)
            value = patch%aux%Global%auxvars(ghosted_id)%den_kg(1)
          case(LIQUID_MOBILITY)
            value = patch%aux%Richards%auxvars(ghosted_id)%kvr
          case(LIQUID_VISCOSITY)
            value = patch%aux%Richards%auxvars(ghosted_id)%kr / &
                    patch%aux%Richards%auxvars(ghosted_id)%kvr
          case default
            call PatchUnsupportedVariable('RICHARDS',ivar,option)
        end select
      else if (associated(patch%aux%PNF)) then
        select case(ivar)
          case(LIQUID_PRESSURE)
            value = patch%aux%PNF%auxvars(ghosted_id)%head
          case default
            call PatchUnsupportedVariable('PNF',ivar,option)
        end select
      else if (associated(patch%aux%ZFlow)) then
        select case(ivar)
          case(LIQUID_PRESSURE,MAXIMUM_PRESSURE)
            value = patch%aux%ZFlow%auxvars(ZERO_INTEGER,ghosted_id)%pres
          case(CAPILLARY_PRESSURE)
            value = patch%aux%ZFlow%auxvars(ZERO_INTEGER,ghosted_id)%pc
          case(LIQUID_HEAD)
            value = patch%aux%ZFlow%auxvars(ZERO_INTEGER,ghosted_id)%pres/ &
                    EARTH_GRAVITY/zflow_density_kg
          case(LIQUID_SATURATION)
            value = patch%aux%ZFlow%auxvars(ZERO_INTEGER,ghosted_id)%sat
          case(SOLUTE_CONCENTRATION)
            value = patch%aux%ZFlow%auxvars(ZERO_INTEGER,ghosted_id)%conc
          case default
            call PatchUnsupportedVariable('ZFLOW',ivar,option)
        end select
      else if (associated(patch%aux%Mphase)) then
        select case(ivar)
          case(MAXIMUM_PRESSURE)
            value = maxval(patch%aux%Global%auxvars(ghosted_id)%pres(1:2))
          case(TEMPERATURE)
            value = patch%aux%Global%auxvars(ghosted_id)%temp
          case(LIQUID_PRESSURE)
            value = patch%aux%Global%auxvars(ghosted_id)%pres(1)
          case(GAS_PRESSURE)
            value = patch%aux%Global%auxvars(ghosted_id)%pres(2)
          case(LIQUID_SATURATION)
            value = patch%aux%Global%auxvars(ghosted_id)%sat(1)
          case(LIQUID_MOLE_FRACTION)
            if (patch%aux%Global%auxvars(ghosted_id)%sat(1) > 0.d0) then
              value = patch%aux%Mphase%auxvars(ghosted_id)% &
                        auxvar_elem(0)%xmol(isubvar)
            else
              value = 0.d0
            endif
          case(LIQUID_ENERGY)
            value = patch%aux%Mphase%auxvars(ghosted_id)%auxvar_elem(0)%u(1)
          case(LIQUID_DENSITY)
            value = patch%aux%Global%auxvars(ghosted_id)%den_kg(1)
          case(LIQUID_VISCOSITY)
            value = patch%aux%Mphase%auxvars(ghosted_id)%auxvar_elem(0)%vis(1)
          case(LIQUID_MOBILITY)
            value = patch%aux%Mphase%auxvars(ghosted_id)%auxvar_elem(0)%kvr(1)
          case(GAS_SATURATION)
            value = patch%aux%Global%auxvars(ghosted_id)%sat(2)
          case(GAS_MOLE_FRACTION)
            if (patch%aux%Global%auxvars(ghosted_id)%sat(2) > 0.d0) then
              value = patch%aux%Mphase%auxvars(ghosted_id)% &
                        auxvar_elem(0)%xmol(2+isubvar)
            else
              value = 0.d0
            endif
          case(GAS_ENERGY)
            value = patch%aux%Mphase%auxvars(ghosted_id)%auxvar_elem(0)%u(2)
          case(GAS_DENSITY)
            value = patch%aux%Global%auxvars(ghosted_id)%den_kg(2)
          case(GAS_VISCOSITY)
            value = patch%aux%Mphase%auxvars(ghosted_id)%auxvar_elem(0)%vis(2)
          case(GAS_MOBILITY)
            value = patch%aux%Mphase%auxvars(ghosted_id)%auxvar_elem(0)%kvr(2)
          case(GAS_DENSITY_MOL)
            value = patch%aux%Global%auxvars(ghosted_id)%den(2)
          case(SC_FUGA_COEFF)
            value = patch%aux%Global%auxvars(ghosted_id)%fugacoeff(1)
          case(SECONDARY_TEMPERATURE)
            local_id = grid%nG2L(ghosted_id)
            value = patch%aux%SC_heat%sec_heat_vars(local_id)%sec_temp(isubvar)
          case default
            call PatchUnsupportedVariable('MPHASE',ivar,option)
        end select
      else if (associated(patch%aux%General)) then
        select case(ivar)
          case(TEMPERATURE)
            value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)%temp
          case(MAXIMUM_PRESSURE)
            value = maxval(patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                           pres(option%liquid_phase:option%gas_phase))
          case(LIQUID_PRESSURE)
            if (output_option%filter_non_state_variables) then
              if (patch%aux%Global%auxvars(ghosted_id)%istate /= GAS_STATE) then
                value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                          pres(option%liquid_phase)
              else
                value = 0.d0
              endif
            else
              value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                        pres(option%liquid_phase)
            endif
          case(GAS_PRESSURE)
            if (output_option%filter_non_state_variables) then
              if (patch%aux%Global%auxvars(ghosted_id)%istate /= &
                  LIQUID_STATE) then
                value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                          pres(option%gas_phase)
              else
                value = 0.d0
              endif
            else
              value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                        pres(option%gas_phase)
            endif
          case(AIR_PRESSURE)
            if (output_option%filter_non_state_variables) then
              if (patch%aux%Global%auxvars(ghosted_id)%istate /= &
                  LIQUID_STATE) then
                value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                          pres(option%air_pressure_id)
              else
                value = 0.d0
              endif
            else
              value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                        pres(option%air_pressure_id)
            endif
          case(CAPILLARY_PRESSURE)
            value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      pres(option%capillary_pressure_id)
          case(VAPOR_PRESSURE)
            value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      pres(option%vapor_pressure_id)
          case(SATURATION_PRESSURE)
            value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      pres(option%saturation_pressure_id)
          case(LIQUID_SATURATION)
            value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      sat(option%liquid_phase)
          case(LIQUID_DENSITY)
            value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      den_kg(option%liquid_phase)
          case(LIQUID_DENSITY_MOL)
            value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      den(option%liquid_phase)
          case(LIQUID_ENERGY)
            if (isubvar == ZERO_INTEGER) then
              value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                        U(option%liquid_phase)
            else
              value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                        U(option%liquid_phase) * &
                      patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                        den(option%liquid_phase)
            endif
          case(LIQUID_MOLE_FRACTION,LIQUID_MASS_FRACTION)
            value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      xmol(isubvar,option%liquid_phase)
            if (ivar == LIQUID_MASS_FRACTION) then
              tempint = isubvar
              tempint2 = tempint+1
              if (tempint2 > 2) tempint2 = 1
              value = value*general_fmw(tempint) / &
                      (value*general_fmw(tempint) + &
                       (1.d0-value)*general_fmw(tempint2))
            endif
          case(LIQUID_MOBILITY)
            value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      mobility(option%liquid_phase)
          case(LIQUID_VISCOSITY)
            value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      kr(option%liquid_phase) / &
                    patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      mobility(option%liquid_phase)
          case(GAS_SATURATION)
            value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      sat(option%gas_phase)
          case(GAS_DENSITY)
            value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      den_kg(option%gas_phase)
          case(GAS_DENSITY_MOL)
            value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      den(option%gas_phase)
          case(GAS_ENERGY)
            if (isubvar == ZERO_INTEGER) then
              value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                        U(option%gas_phase)
            else
              value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                        U(option%gas_phase) * &
                      patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                        den(option%gas_phase)
            endif
          case(GAS_MOLE_FRACTION,GAS_MASS_FRACTION)
            value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      xmol(isubvar,option%gas_phase)
            if (ivar == GAS_MASS_FRACTION) then
              tempint = isubvar
              tempint2 = tempint+1
              if (tempint2 > 2) tempint2 = 1
              value = value*general_fmw(tempint) / &
                      (value*general_fmw(tempint) + &
                       (1.d0-value)*general_fmw(tempint2))
            endif
          case(GAS_MOBILITY)
            value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      mobility(option%gas_phase)
          case(GAS_VISCOSITY)
            value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      kr(option%gas_phase) / &
                    patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      mobility(option%gas_phase)
          case(ICE_SATURATION)
            value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      sat(option%ice_phase)
          case(HYDRATE_SATURATION)
            value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                      sat(option%hydrate_phase)
          case default
            call PatchUnsupportedVariable('GENERAL',ivar,option)
        end select

      else if (associated(patch%aux%Hydrate)) then
        select case(ivar)
          case(TEMPERATURE)
            value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)%temp
          case(MAXIMUM_PRESSURE)
            value = maxval(patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                           pres(option%liquid_phase:option%gas_phase))
          case(LIQUID_PRESSURE)
            if (output_option%filter_non_state_variables) then
              if (patch%aux%Global%auxvars(ghosted_id)%istate /= GAS_STATE) then
                value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                          pres(option%liquid_phase)
              else
                value = 0.d0
              endif
            else
              value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                        pres(option%liquid_phase)
            endif
          case(GAS_PRESSURE)
            if (output_option%filter_non_state_variables) then
              if (patch%aux%Global%auxvars(ghosted_id)%istate /= &
                  LIQUID_STATE) then
                value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                          pres(option%gas_phase)
              else
                value = 0.d0
              endif
            else
              value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                        pres(option%gas_phase)
            endif
          case(AIR_PRESSURE)
            if (output_option%filter_non_state_variables) then
              if (patch%aux%Global%auxvars(ghosted_id)%istate /= &
                  LIQUID_STATE) then
                value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                          pres(option%air_pressure_id)
              else
                value = 0.d0
              endif
            else
              value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                        pres(option%air_pressure_id)
            endif
          case(CAPILLARY_PRESSURE)
            value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      pres(option%capillary_pressure_id)
          case(VAPOR_PRESSURE)
            value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      pres(option%vapor_pressure_id)
          case(SATURATION_PRESSURE)
            value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      pres(option%saturation_pressure_id)
          case(LIQUID_SATURATION)
            value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      sat(option%liquid_phase)
          case(LIQUID_DENSITY)
            value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      den_kg(option%liquid_phase)
          case(LIQUID_DENSITY_MOL)
            value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      den(option%liquid_phase)
          case(LIQUID_ENERGY)
            if (isubvar == ZERO_INTEGER) then
              value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                        U(option%liquid_phase)
            else
              value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                        U(option%liquid_phase) * &
                      patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                        den(option%liquid_phase)
            endif
          case(LIQUID_MOLE_FRACTION,LIQUID_MASS_FRACTION)
            value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      xmol(isubvar,option%liquid_phase)
            if (ivar == LIQUID_MASS_FRACTION) then
              tempint = isubvar
              tempint2 = tempint+1
              !MAN: correct this
              if (tempint2 > 2) tempint2 = 1
              value = value*general_fmw(tempint) / &
                      (value*general_fmw(tempint) + &
                       (1.d0-value)*general_fmw(tempint2))
            endif
          case(LIQUID_MOBILITY)
            value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      mobility(option%liquid_phase)
          case(LIQUID_VISCOSITY)
            value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      kr(option%liquid_phase) / &
                    patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      mobility(option%liquid_phase)
          case(GAS_SATURATION)
            value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      sat(option%gas_phase)
          case(GAS_DENSITY)
            value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      den_kg(option%gas_phase)
          case(GAS_DENSITY_MOL)
            value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      den(option%gas_phase)
          case(GAS_ENERGY)
            if (isubvar == ZERO_INTEGER) then
              value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                        U(option%gas_phase)
            else
              value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                        U(option%gas_phase) * &
                      patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                        den(option%gas_phase)
            endif
          case(GAS_MOLE_FRACTION,GAS_MASS_FRACTION)
            value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      xmol(isubvar,option%gas_phase)
            if (ivar == GAS_MASS_FRACTION) then
              tempint = isubvar
              tempint2 = tempint+1
              if (tempint2 > 2) tempint2 = 1
              !MAN: correct this
              value = value*general_fmw(tempint) / &
                      (value*general_fmw(tempint) + &
                       (1.d0-value)*general_fmw(tempint2))
            endif
          case(GAS_MOBILITY)
            value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      mobility(option%gas_phase)
          case(GAS_VISCOSITY)
            value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      kr(option%gas_phase) / &
                    patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      mobility(option%gas_phase)
          case(ICE_SATURATION)
            value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      sat(option%ice_phase)
          case(HYDRATE_SATURATION)
            value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                      sat(option%hydrate_phase)
          case default
            call PatchUnsupportedVariable('HYDRATE',ivar,option)
        end select

      else if (associated(patch%aux%WIPPFlo)) then
        select case(ivar)
          case(MAXIMUM_PRESSURE)
            value = maxval(patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                           pres(option%liquid_phase:option%gas_phase))
          case(LIQUID_PRESSURE)
            if (output_option%filter_non_state_variables) then
              if (patch%aux%Global%auxvars(ghosted_id)%istate /= GAS_STATE) then
                value = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                          pres(option%liquid_phase)
              else
                value = 0.d0
              endif
            else
              value = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                        pres(option%liquid_phase)
            endif
          case(GAS_PRESSURE)
            if (output_option%filter_non_state_variables) then
              if (patch%aux%Global%auxvars(ghosted_id)%istate /= &
                  LIQUID_STATE) then
                value = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                          pres(option%gas_phase)
              else
                value = 0.d0
              endif
            else
              value = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                        pres(option%gas_phase)
            endif
          case(CAPILLARY_PRESSURE)
            value = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                      pres(option%capillary_pressure_id)
          case(SATURATION_PRESSURE)
            value = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                      pres(option%saturation_pressure_id)
          case(LIQUID_SATURATION)
            value = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                      sat(option%liquid_phase)
          case(LIQUID_DENSITY)
            value = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                      den_kg(option%liquid_phase)
          case(LIQUID_DENSITY_MOL)
            value = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                      den(option%liquid_phase)
          case(LIQUID_MOBILITY)
            value = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                      mobility(option%liquid_phase)
          case(LIQUID_VISCOSITY)
            value = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                      mu(option%liquid_phase)
          case(GAS_SATURATION)
            value = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                      sat(option%gas_phase)
          case(GAS_DENSITY)
            value = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                      den_kg(option%gas_phase)
          case(GAS_DENSITY_MOL)
            value = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                      den(option%gas_phase)
          case(GAS_MOBILITY)
            value = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                      mobility(option%gas_phase)
          case(GAS_VISCOSITY)
            value = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                      mu(option%gas_phase)
          case default
            call PatchUnsupportedVariable('WIPP_FLOW',ivar,option)
        end select
      else ! null flow mode
        select case(ivar)
          case(LIQUID_SATURATION)
            value = patch%aux%Global%auxvars(ghosted_id)% &
                      sat(option%liquid_phase)
          case default
            call PatchUnsupportedVariable('NULL_FLOW',ivar,option)
        end select
      endif

    ! NUCLEAR_WASTE_TRANSPORT:
    case(TOTAL_BULK_CONC,AQUEOUS_EQ_CONC,MNRL_EQ_CONC,SORB_EQ_CONC, &
         MNRL_VOLUME_FRACTION,NWT_AUXILIARY)
      select case(ivar)
        case(TOTAL_BULK_CONC)
          value = patch%aux%NWT%auxvars(ghosted_id)%total_bulk_conc(isubvar)
        case(AQUEOUS_EQ_CONC)
          value = patch%aux%NWT%auxvars(ghosted_id)%aqueous_eq_conc(isubvar)
        case(MNRL_EQ_CONC)
          value = patch%aux%NWT%auxvars(ghosted_id)%mnrl_eq_conc(isubvar)
        case(SORB_EQ_CONC)
          value = patch%aux%NWT%auxvars(ghosted_id)%sorb_eq_conc(isubvar)
        case(MNRL_VOLUME_FRACTION)
          value = patch%aux%NWT%auxvars(ghosted_id)%mnrl_vol_frac(isubvar)
        case(NWT_AUXILIARY)
          value = patch%aux%NWT%auxvars(ghosted_id)%auxiliary_data(isubvar)
      end select
      if ( (patch%aux%NWT%truncate_output) .and. &
           ((value < 1.d-99) .and. (value > 0.d0)) ) then
        value = 1.d-99
      endif

    case(PH,PE,EH,O2,PRIMARY_MOLALITY,PRIMARY_MOLARITY,SECONDARY_MOLALITY, &
         SECONDARY_MOLARITY, TOTAL_MOLALITY,TOTAL_MOLARITY, &
         MINERAL_VOLUME_FRACTION,MINERAL_RATE,MINERAL_SATURATION_INDEX, &
         MINERAL_SURFACE_AREA, &
         SURFACE_CMPLX,SURFACE_CMPLX_FREE,SURFACE_SITE_DENSITY, &
         KIN_SURFACE_CMPLX,KIN_SURFACE_CMPLX_FREE, PRIMARY_ACTIVITY_COEF, &
         SECONDARY_ACTIVITY_COEF,PRIMARY_KD, TOTAL_SORBED, &
         TOTAL_SORBED_MOBILE,COLLOID_MOBILE,COLLOID_IMMOBILE,AGE,TOTAL_BULK, &
         IMMOBILE_SPECIES,GAS_CONCENTRATION,REACTION_AUXILIARY)

      select case(ivar)
        case(PH)
          call RRedoxCalcpH(patch%aux%RT%auxvars(ghosted_id), &
                            patch%aux%Global%auxvars(ghosted_id), &
                            reaction,ph0,option)
          value = ph0
        case(EH)
          call RRedoxCalcEhpe(patch%aux%RT%auxvars(ghosted_id), &
                                patch%aux%Global%auxvars(ghosted_id), &
                                reaction,eh0,pe0,option)
          value = eh0
        case(PE)
          call RRedoxCalcEhpe(patch%aux%RT%auxvars(ghosted_id), &
                                patch%aux%Global%auxvars(ghosted_id), &
                                reaction,eh0,pe0,option)
          value = pe0
        case(O2)
          call RRedoxCalcLnFO2(patch%aux%RT%auxvars(ghosted_id), &
                               patch%aux%Global%auxvars(ghosted_id), &
                               reaction,lnQKgas,option)
          value = lnQKgas * LN_TO_LOG
        case(PRIMARY_MOLALITY)
          value = patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar)
        case(PRIMARY_MOLARITY)
          value = patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar)*xmass * &
                  patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase) / 1000.d0
        case(SECONDARY_MOLALITY)
          value = patch%aux%RT%auxvars(ghosted_id)%sec_molal(isubvar)
        case(SECONDARY_MOLARITY)
          value = patch%aux%RT%auxvars(ghosted_id)%sec_molal(isubvar)*xmass * &
                  patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase) / 1000.d0
        case(TOTAL_MOLALITY)
          value = patch%aux%RT%auxvars(ghosted_id)%total(isubvar,iphase) / &
                  xmass / &
                  patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase)*1000.d0
        case(TOTAL_MOLARITY)
          value = patch%aux%RT%auxvars(ghosted_id)%total(isubvar,iphase)
        case(TOTAL_BULK) ! mol/m^3 bulk
          ! add in total molarity and convert to mol/m^3 bulk
          value = &
              patch%aux%RT%auxvars(ghosted_id)%total(isubvar,iphase) * &
              patch%aux%Material%auxvars(ghosted_id)%porosity * &
                                                              ! mol/L -> mol/m^3
              patch%aux%Global%auxvars(ghosted_id)%sat(iphase) * 1.d-3
          ! add in total sorbed.  already in mol/m^3 bulk
          if (patch%reaction%nsorb > 0) then
            if (patch%reaction%surface_complexation%neqsrfcplxrxn > 0) then
              value = value + &
                patch%aux%RT%auxvars(ghosted_id)%total_sorb_eq(isubvar)
            endif
            if (patch%reaction%surface_complexation%nkinmrsrfcplxrxn > 0) then
              do irxn = 1, patch%reaction%surface_complexation%nkinmrsrfcplxrxn
                do irate = 1, &
                  patch%reaction%surface_complexation%kinmr_nrate(irxn)
                  value = value + &
                      patch%aux%RT%auxvars(ghosted_id)% &
                        kinmr_total_sorb(isubvar,irate,irxn)
                enddo
              enddo
            endif
          endif
        case(GAS_CONCENTRATION)
          iphase = 2
          if (patch%aux%Global%auxvars(ghosted_id)%sat(iphase) > 0.d0) then
            value = patch%aux%RT%auxvars(ghosted_id)%gas_pp(isubvar)
          else
            value = 0.d0
          endif
        case(MINERAL_VOLUME_FRACTION)
          value = patch%aux%RT%auxvars(ghosted_id)%mnrl_volfrac(isubvar)
        case(MINERAL_SURFACE_AREA)
          value = patch%aux%RT%auxvars(ghosted_id)%mnrl_area(isubvar)
        case(MINERAL_RATE)
          value = patch%aux%RT%auxvars(ghosted_id)%mnrl_rate(isubvar)
        case(MINERAL_SATURATION_INDEX)
          value = RMineralSaturationIndex(isubvar, &
                                         patch%aux%RT%auxvars(ghosted_id), &
                                         patch%aux%Global%auxvars(ghosted_id), &
                                         reaction,option)
        case(IMMOBILE_SPECIES)
          value = patch%aux%RT%auxvars(ghosted_id)%immobile(isubvar)
        case(SURFACE_CMPLX)
          if (associated(patch%aux%RT%auxvars(ghosted_id)%eqsrfcplx_conc)) then
            value = patch%aux%RT%auxvars(ghosted_id)%eqsrfcplx_conc(isubvar)
          else
            value = UNINITIALIZED_DOUBLE
          endif
        case(SURFACE_CMPLX_FREE)
          value = &
            patch%aux%RT%auxvars(ghosted_id)%srfcplxrxn_free_site_conc(isubvar)
        case(SURFACE_SITE_DENSITY)
          select case(reaction%surface_complexation% &
                        srfcplxrxn_surf_type(isubvar))
            case(ROCK_SURFACE)
              value = reaction%surface_complexation% &
                        srfcplxrxn_site_density(isubvar)* &
                      material_auxvars(ghosted_id)%soil_particle_density * &
                      (1.d0-material_auxvars(ghosted_id)%porosity)
            case(MINERAL_SURFACE)
              value = reaction%surface_complexation% &
                        srfcplxrxn_site_density(isubvar)* &
                      patch%aux%RT%auxvars(ghosted_id)% &
                        mnrl_volfrac(reaction%surface_complexation% &
                                       srfcplxrxn_to_surf(isubvar))
            case(COLLOID_SURFACE)
                option%io_buffer = 'Printing of surface site density for ' // &
                  'colloidal surfaces not implemented.'
                call PrintErrMsg(option)
            case(NULL_SURFACE)
              value = reaction%surface_complexation% &
                        srfcplxrxn_site_density(isubvar)
          end select
        case(KIN_SURFACE_CMPLX)
          value = patch%aux%RT%auxvars(ghosted_id)%kinsrfcplx_conc(isubvar,1)
        case(KIN_SURFACE_CMPLX_FREE)
          value = &
            patch%aux%RT%auxvars(ghosted_id)%kinsrfcplx_free_site_conc(isubvar)
        case(PRIMARY_ACTIVITY_COEF)
          value = patch%aux%RT%auxvars(ghosted_id)%pri_act_coef(isubvar)
        case(SECONDARY_ACTIVITY_COEF)
          value = patch%aux%RT%auxvars(ghosted_id)%sec_act_coef(isubvar)
        case(PRIMARY_KD)
          call ReactionComputeKd(isubvar,value, &
                                 patch%aux%RT%auxvars(ghosted_id), &
                                 patch%aux%Global%auxvars(ghosted_id), &
                                 material_auxvars(ghosted_id), &
                                 patch%reaction,option)
        case(TOTAL_SORBED)
          if (patch%reaction%nsorb > 0) then
            if (patch%reaction%neqsorb > 0) then
              value = patch%aux%RT%auxvars(ghosted_id)%total_sorb_eq(isubvar)
            endif
            if (patch%reaction%surface_complexation%nkinmrsrfcplxrxn > 0) then
              value = 0.d0
              do irxn = 1, patch%reaction%surface_complexation%nkinmrsrfcplxrxn
                do irate = 1, &
                  patch%reaction%surface_complexation%kinmr_nrate(irxn)
                  value = value + &
                    patch%aux%RT%auxvars(ghosted_id)% &
                      kinmr_total_sorb(isubvar,irate,irxn)
                enddo
              enddo
            endif
          endif
        case(TOTAL_SORBED_MOBILE)
          if (patch%reaction%nsorb > 0 .and. patch%reaction%ncollcomp > 0) then
            value = &
              patch%aux%RT%auxvars(ghosted_id)%colloid%total_eq_mob(isubvar)
          endif
        case(COLLOID_MOBILE)
          if (patch%reaction%print_tot_conc_type == TOTAL_MOLALITY) then
            value = patch%aux%RT%auxvars(ghosted_id)% &
                      colloid%conc_mob(isubvar) / &
                    patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase)*1000.d0
          else
            value = patch%aux%RT%auxvars(ghosted_id)%colloid%conc_mob(isubvar)
          endif
        case(COLLOID_IMMOBILE)
          if (patch%reaction%print_tot_conc_type == TOTAL_MOLALITY) then
            value = patch%aux%RT%auxvars(ghosted_id)% &
                      colloid%conc_imb(isubvar) / &
                    patch%aux%Global%auxvars(ghosted_id)%den_kg(iphase)*1000.d0
          else
            value = patch%aux%RT%auxvars(ghosted_id)%colloid%conc_imb(isubvar)
          endif
        case(AGE)
          if (patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar) > &
              0.d0) then
            value = patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar) / &
            patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar2) / &
            output_option%tconv
          endif
        case(REACTION_AUXILIARY)
          value = patch%aux%RT%auxvars(ghosted_id)%auxiliary_data(isubvar)
        case default
          call PatchUnsupportedVariable('REACTIVE_TRANSPORT',ivar,option)
      end select
    case(STATE,PHASE)
      value = patch%aux%Global%auxvars(ghosted_id)%istate
    case(POROSITY,BASE_POROSITY,INITIAL_POROSITY, &
         VOLUME,TORTUOSITY,SOIL_COMPRESSIBILITY,SOIL_REFERENCE_PRESSURE, &
         ELECTRICAL_CONDUCTIVITY)
      value = MaterialAuxVarGetValue(material_auxvars(ghosted_id),ivar)
    case(PERMEABILITY,PERMEABILITY_X,PERMEABILITY_Y, PERMEABILITY_Z, &
         PERMEABILITY_XY,PERMEABILITY_XZ,PERMEABILITY_YZ, &
         GAS_PERMEABILITY,GAS_PERMEABILITY_X,GAS_PERMEABILITY_Y, &
         GAS_PERMEABILITY_Z)
      ivar_temp = ivar
      ! only liquid permeabilities in x, y, z are stored.
      select case(ivar)
        case(PERMEABILITY,GAS_PERMEABILITY,GAS_PERMEABILITY_X)
          ivar_temp = PERMEABILITY_X
        case(GAS_PERMEABILITY_Y)
          ivar_temp = PERMEABILITY_Y
        case(GAS_PERMEABILITY_Z)
          ivar_temp = PERMEABILITY_Z
      end select
      value = MaterialAuxVarGetValue(material_auxvars(ghosted_id),ivar_temp)
      select case(option%iflowmode)
        case(WF_MODE)
          call WIPPFloScalePerm(patch%aux%WIPPFlo%auxvars(ZERO_INTEGER, &
                                                          ghosted_id), &
                                material_auxvars(ghosted_id), &
                                value,ivar)
      end select
    case(LIQUID_RELATIVE_PERMEABILITY)
      select case(option%iflowmode)
        case(RICHARDS_MODE)
          value = patch%aux%Richards%auxvars(ghosted_id)%kr
        case(ZFLOW_MODE)
          value = patch%aux%ZFlow%auxvars(ZERO_INTEGER,ghosted_id)%kr
        case(TH_MODE,TH_TS_MODE)
          value = patch%aux%TH%auxvars(ghosted_id)%kvr / &
                  patch%aux%TH%auxvars(ghosted_id)%vis
        case(G_MODE)
          value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                    kr(option%liquid_phase)
        case(H_MODE)
          value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                    kr(option%liquid_phase)
        case(WF_MODE)
          value = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                    kr(option%liquid_phase)
        case default
          option%io_buffer = 'Output of liquid phase relative permeability &
            &not supported for current flow mode.'
      end select
    case(GAS_RELATIVE_PERMEABILITY)
      select case(option%iflowmode)
        case(G_MODE)
          value = patch%aux%General%auxvars(ZERO_INTEGER,ghosted_id)% &
                    kr(option%gas_phase)
        case(H_MODE)
          value = patch%aux%Hydrate%auxvars(ZERO_INTEGER,ghosted_id)% &
                    kr(option%gas_phase)
        case(WF_MODE)
          value = patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,ghosted_id)% &
                    kr(option%gas_phase)
        case default
          option%io_buffer = 'Output of gas phase relative permeability &
            &not supported for current flow mode.'
      end select
    case(MATERIAL_ID)
      value = patch%imat_internal_to_external(abs(patch%imat(ghosted_id)))
    case(FRACTURE)
      value = 0.d0
      if (associated(material_auxvars(ghosted_id)%fracture)) then
        if (material_auxvars(ghosted_id)%fracture%fracture_is_on) then
          value = 1.d0
        endif
      endif
    case(ELECTRICAL_POTENTIAL_DIPOLE)
      call ERTAuxCheckElectrodeBounds(size(patch%aux%ERT%auxvars(1)% &
                                      potential),isubvar,isubvar2,option)
      value = patch%aux%ERT%auxvars(ghosted_id)%potential(isubvar) - &
              patch%aux%ERT%auxvars(ghosted_id)%potential(isubvar2)
    case(ELECTRICAL_POTENTIAL)
      call ERTAuxCheckElectrodeBounds(size(patch%aux%ERT%auxvars(1)% &
                                      potential),isubvar,isubvar2,option)
      value = patch%aux%ERT%auxvars(ghosted_id)%potential(isubvar)
    case(ELECTRICAL_JACOBIAN)
      call ERTAuxCheckElectrodeBounds(size(patch%aux%ERT%auxvars(1)% &
                                      jacobian),isubvar,isubvar2,option)
      value = patch%aux%ERT%auxvars(ghosted_id)%jacobian(isubvar)
    case(PROCESS_ID)
      value = grid%nG2A(ghosted_id)
    case(NATURAL_ID)
      value = option%myrank
    ! Need to fix the below two cases (they assume only one component) -- SK 02/06/13
    case(SECONDARY_CONCENTRATION)
      ! Note that the units are in mol/kg
      local_id = grid%nG2L(ghosted_id)
      value = patch%aux%SC_RT%sec_transport_vars(local_id)% &
                                   sec_rt_auxvar(isubvar)%total(isubvar2,1)
    case(SECONDARY_CONCENTRATION_GAS)
      local_id = grid%nG2L(ghosted_id)
      value = patch%aux%SC_RT%sec_transport_vars(local_id)% &
                                   sec_rt_auxvar(isubvar)%gas_pp(isubvar2)
    case(SEC_MIN_VOLFRAC)
      local_id = grid%nG2L(ghosted_id)
      value = patch%aux%SC_RT%sec_transport_vars(local_id)% &
              sec_rt_auxvar(isubvar)%mnrl_volfrac(isubvar2)
    case(SEC_MIN_RATE)
      local_id = grid%nG2L(ghosted_id)
      value = patch%aux%SC_RT%sec_transport_vars(local_id)% &
              sec_rt_auxvar(isubvar)%mnrl_rate(isubvar2)
    case(SEC_MIN_SI)
      local_id = grid%nG2L(ghosted_id)
      value = RMineralSaturationIndex(isubvar2,&
                                      patch%aux%SC_RT% &
                                      sec_transport_vars(local_id)% &
                                      sec_rt_auxvar(isubvar), &
                                      patch%aux%Global%auxvars(ghosted_id),&
                                      reaction,option)
    case(SALINITY)
      value = patch%aux%Global%auxvars(ghosted_id)%m_nacl(ONE_INTEGER)
    case(RESIDUAL)
      local_id = grid%nG2L(ghosted_id)
      call VecGetArrayF90(field%flow_r,vec_ptr2,ierr);CHKERRQ(ierr)
      value = vec_ptr2((local_id-1)*option%nflowdof+isubvar)
      call VecRestoreArrayF90(field%flow_r,vec_ptr2,ierr);CHKERRQ(ierr)
    case(X_COORDINATE)
      value = grid%x(ghosted_id)
    case(Y_COORDINATE)
      value = grid%y(ghosted_id)
    case(Z_COORDINATE)
      value = grid%z(ghosted_id)
    case default
      call PatchUnsupportedVariable(ivar,option)
  end select

  PatchGetVariableValueAtCell = value

end function PatchGetVariableValueAtCell

! ************************************************************************** !

subroutine PatchSetVariable(patch,field,option,vec,vec_format,ivar,isubvar)
  !
  ! Sets variables indexed by ivar and isubvar within a patch
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/08
  !

  use Grid_module
  use Option_module
  use Field_module
  use Variables_module
  use General_Aux_module
  use WIPP_Flow_Aux_module
  use Material_Aux_module

  implicit none

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  Vec :: vec
  PetscInt :: vec_format
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: iphase, istate

  PetscInt :: local_id, ghosted_id
  type(grid_type), pointer :: grid
  type(material_auxvar_type), pointer :: material_auxvars(:)
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscErrorCode :: ierr

  grid => patch%grid
  material_auxvars => patch%aux%Material%auxvars

  call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)

  if (vec_format == NATURAL) then
    call PrintErrMsg(option,&
                     'NATURAL vector format not supported by PatchSetVariable')
  endif

  iphase = 1
  select case(ivar)
    case(TEMPERATURE,LIQUID_PRESSURE,GAS_PRESSURE,LIQUID_SATURATION, &
         GAS_SATURATION,AIR_PRESSURE,CAPILLARY_PRESSURE, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY, &
         LIQUID_DENSITY,GAS_DENSITY,GAS_DENSITY_MOL,LIQUID_VISCOSITY, &
         GAS_VISCOSITY,HYDRATE_SATURATION,ICE_SATURATION, &
         LIQUID_MOBILITY,GAS_MOBILITY)

      if (associated(patch%aux%TH)) then
        select case(ivar)
          case(TEMPERATURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%auxvars(grid%nL2G(local_id))%temp = &
                  vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%auxvars(ghosted_id)%temp = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_PRESSURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%auxvars(grid%nL2G(local_id))%pres = &
                  vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%auxvars(ghosted_id)%pres = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_SATURATION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%auxvars(grid%nL2G(local_id))%sat = &
                  vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%auxvars(ghosted_id)%sat = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_DENSITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%auxvars(grid%nL2G(local_id))%den_kg = &
                  vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%auxvars(ghosted_id)%den_kg = &
                  vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY)
            call PrintErrMsg(option,'GAS_MOLE_FRACTION not supported by TH')
          case(GAS_SATURATION)
            if (option%flow%th_freezing) then
              if (vec_format == GLOBAL) then
                do local_id=1,grid%nlmax
                  patch%aux%TH%auxvars(grid%nL2G(local_id))%ice%sat_gas = &
                    vec_ptr(local_id)
                enddo
              else if (vec_format == LOCAL) then
                do ghosted_id=1,grid%ngmax
                  patch%aux%TH%auxvars(ghosted_id)%ice%sat_gas = vec_ptr(ghosted_id)
                enddo
              endif
            endif
          case(ICE_SATURATION)
            if (option%flow%th_freezing) then
              if (vec_format == GLOBAL) then
                do local_id=1,grid%nlmax
                  patch%aux%TH%auxvars(grid%nL2G(local_id))%ice%sat_ice = &
                    vec_ptr(local_id)
                enddo
              else if (vec_format == LOCAL) then
                do ghosted_id=1,grid%ngmax
                  patch%aux%TH%auxvars(ghosted_id)%ice%sat_ice = vec_ptr(ghosted_id)
                enddo
              endif
            endif
          case(ICE_DENSITY)
            if (option%flow%th_freezing) then
              if (vec_format == GLOBAL) then
                do local_id=1,grid%nlmax
                  patch%aux%TH%auxvars(grid%nL2G(local_id))%ice%den_ice = &
                    vec_ptr(local_id)
                enddo
              else if (vec_format == LOCAL) then
                do ghosted_id=1,grid%ngmax
                  patch%aux%TH%auxvars(ghosted_id)%ice%den_ice = vec_ptr(ghosted_id)
                enddo
              endif
            endif
          case(LIQUID_VISCOSITY)
          case(GAS_VISCOSITY)
          case(LIQUID_MOLE_FRACTION)
            call PrintErrMsg(option,'LIQUID_MOLE_FRACTION not supported by TH')
          case(LIQUID_ENERGY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%TH%auxvars(grid%nL2G(local_id))%u = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%TH%auxvars(ghosted_id)%u = vec_ptr(ghosted_id)
              enddo
            endif
        end select
      else if (associated(patch%aux%Richards)) then
        select case(ivar)
          case(TEMPERATURE)
            call PrintErrMsg(option,'TEMPERATURE not supported by Richards')
          case(GAS_SATURATION)
            call PrintErrMsg(option,'GAS_SATURATION not supported by Richards')
          case(GAS_DENSITY)
            call PrintErrMsg(option,'GAS_DENSITY not supported by Richards')
          case(LIQUID_MOLE_FRACTION)
            call PrintErrMsg(option,'LIQUID_MOLE_FRACTION not supported by Richards')
          case(GAS_MOLE_FRACTION)
            call PrintErrMsg(option,'GAS_MOLE_FRACTION not supported by Richards')
          case(LIQUID_VISCOSITY)
            call PrintErrMsg(option,'LIQUID_VISCOSITY not supported by Richards')
          case(GAS_VISCOSITY)
            call PrintErrMsg(option,'GAS_VISCOSITY not supported by Richards')
          case(GAS_MOBILITY)
            call PrintErrMsg(option,'GAS_MOBILITY not supported by Richards')
          case(LIQUID_ENERGY)
            call PrintErrMsg(option,'LIQUID_ENERGY not supported by Richards')
          case(GAS_ENERGY)
            call PrintErrMsg(option,'GAS_ENERGY not supported by Richards')
          case(LIQUID_PRESSURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%auxvars(grid%nL2G(local_id))%pres(1) = &
                  vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%auxvars(ghosted_id)%pres(1) = &
                  vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_SATURATION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%auxvars(grid%nL2G(local_id))%sat(1) = &
                  vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%auxvars(ghosted_id)%sat(1) = &
                  vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_DENSITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%auxvars(grid%nL2G(local_id))%den_kg(1) = &
                  vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%auxvars(ghosted_id)%den_kg(1) = &
                  vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_MOBILITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Richards%auxvars(grid%nL2G(local_id))%kvr = &
                  vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Richards%auxvars(ghosted_id)%kvr = vec_ptr(ghosted_id)
              enddo
            endif
        end select
      else if (associated(patch%aux%Mphase)) then
        select case(ivar)
          case(TEMPERATURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%auxvars(grid%nL2G(local_id))% &
                  auxvar_elem(0)%temp = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%auxvars(ghosted_id)% &
                  auxvar_elem(0)%temp = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_PRESSURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%auxvars(grid%nL2G(local_id))%pres(1) = &
                  vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%auxvars(ghosted_id)%pres(1) = &
                  vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_PRESSURE)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%auxvars(grid%nL2G(local_id))%pres(2) = &
                  vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%auxvars(ghosted_id)%pres(2) = &
                  vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_SATURATION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%auxvars(grid%nL2G(local_id))% &
                  auxvar_elem(0)%sat(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%auxvars(ghosted_id)% &
                  auxvar_elem(0)%sat(1) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_DENSITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%auxvars(grid%nL2G(local_id))% &
                  auxvar_elem(0)%den(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%auxvars(ghosted_id)% &
                  auxvar_elem(0)%den(1) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_VISCOSITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%auxvars(grid%nL2G(local_id))% &
                  auxvar_elem(0)%vis(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%auxvars(ghosted_id)% &
                  auxvar_elem(0)%vis(1) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_MOBILITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%auxvars(grid%nL2G(local_id))% &
                  auxvar_elem(0)%kvr(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%auxvars(ghosted_id)% &
                  auxvar_elem(0)%kvr(1) = vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_SATURATION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%auxvars(grid%nL2G(local_id))% &
                  auxvar_elem(0)%sat(2) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%auxvars(ghosted_id)% &
                  auxvar_elem(0)%sat(2) = vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_MOLE_FRACTION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%auxvars(grid%nL2G(local_id))% &
                  auxvar_elem(0)%xmol(2+isubvar) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%auxvars(ghosted_id)% &
                  auxvar_elem(0)%xmol(2+isubvar) = vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_ENERGY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%auxvars(grid%nL2G(local_id))% &
                  auxvar_elem(0)%u(2) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%auxvars(ghosted_id)% &
                  auxvar_elem(0)%u(2) = vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_DENSITY, GAS_DENSITY_MOL)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%auxvars(grid%nL2G(local_id))% &
                  auxvar_elem(0)%den(2) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%auxvars(ghosted_id)% &
                  auxvar_elem(0)%den(2) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_MOLE_FRACTION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%auxvars(grid%nL2G(local_id))% &
                  auxvar_elem(0)%xmol(isubvar) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%auxvars(ghosted_id)% &
                  auxvar_elem(0)%xmol(isubvar) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_ENERGY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Mphase%auxvars(grid%nL2G(local_id))% &
                  auxvar_elem(0)%u(1) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Mphase%auxvars(ghosted_id)% &
                  auxvar_elem(0)%u(1) = vec_ptr(ghosted_id)
              enddo
            endif
        end select
      else if (associated(patch%aux%General)) then
        select case(ivar)
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              patch%aux%General%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                temp = vec_ptr(local_id)
            enddo
          case(LIQUID_PRESSURE)
            do local_id=1,grid%nlmax
              patch%aux%General%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                pres(option%liquid_phase) = vec_ptr(local_id)
            enddo
          case(GAS_PRESSURE)
            do local_id=1,grid%nlmax
              patch%aux%General%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                pres(option%gas_phase) = vec_ptr(local_id)
            enddo
          case(AIR_PRESSURE)
            do local_id=1,grid%nlmax
              patch%aux%General%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                pres(option%air_pressure_id) = vec_ptr(local_id)
            enddo
          case(CAPILLARY_PRESSURE)
            do local_id=1,grid%nlmax
              patch%aux%General%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                pres(option%capillary_pressure_id) = vec_ptr(local_id)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              patch%aux%General%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                sat(option%liquid_phase) = vec_ptr(local_id)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              patch%aux%General%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
               den_kg(option%liquid_phase) = vec_ptr(local_id)
            enddo
          case(LIQUID_ENERGY)
            do local_id=1,grid%nlmax
              patch%aux%General%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                U(option%liquid_phase) = vec_ptr(local_id)
            enddo
          case(LIQUID_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              patch%aux%General%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                xmol(isubvar,option%liquid_phase) = vec_ptr(local_id)
            enddo
          case(GAS_SATURATION)
            do local_id=1,grid%nlmax
              patch%aux%General%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                sat(option%gas_phase) = vec_ptr(local_id)
            enddo
          case(GAS_DENSITY)
            do local_id=1,grid%nlmax
              patch%aux%General%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                den_kg(option%gas_phase) = vec_ptr(local_id)
            enddo
          case(GAS_ENERGY)
            do local_id=1,grid%nlmax
              patch%aux%General%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                U(option%gas_phase) = vec_ptr(local_id)
            enddo
          case(GAS_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              patch%aux%General%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                xmol(isubvar,option%gas_phase) = vec_ptr(local_id)
            enddo
          case(HYDRATE_SATURATION)
            do local_id=1,grid%nlmax
              patch%aux%General%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                sat(option%hydrate_phase) = vec_ptr(local_id)
            enddo
          case(ICE_SATURATION)
            do local_id=1,grid%nlmax
              patch%aux%General%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                sat(option%ice_phase) = vec_ptr(local_id)
            enddo
        end select
      else if (associated(patch%aux%Hydrate)) then
        select case(ivar)
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              patch%aux%Hydrate%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                temp = vec_ptr(local_id)
            enddo
          case(LIQUID_PRESSURE)
            do local_id=1,grid%nlmax
              patch%aux%Hydrate%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                pres(option%liquid_phase) = vec_ptr(local_id)
            enddo
          case(GAS_PRESSURE)
            do local_id=1,grid%nlmax
              patch%aux%Hydrate%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                pres(option%gas_phase) = vec_ptr(local_id)
            enddo
          case(AIR_PRESSURE)
            do local_id=1,grid%nlmax
              patch%aux%Hydrate%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                pres(option%air_pressure_id) = vec_ptr(local_id)
            enddo
          case(CAPILLARY_PRESSURE)
            do local_id=1,grid%nlmax
              patch%aux%Hydrate%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                pres(option%capillary_pressure_id) = vec_ptr(local_id)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              patch%aux%Hydrate%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                sat(option%liquid_phase) = vec_ptr(local_id)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              patch%aux%Hydrate%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
               den_kg(option%liquid_phase) = vec_ptr(local_id)
            enddo
          case(LIQUID_ENERGY)
            do local_id=1,grid%nlmax
              patch%aux%Hydrate%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                U(option%liquid_phase) = vec_ptr(local_id)
            enddo
          case(LIQUID_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              patch%aux%Hydrate%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                xmol(isubvar,option%liquid_phase) = vec_ptr(local_id)
            enddo
          case(GAS_SATURATION)
            do local_id=1,grid%nlmax
              patch%aux%Hydrate%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                sat(option%gas_phase) = vec_ptr(local_id)
            enddo
          case(GAS_DENSITY)
            do local_id=1,grid%nlmax
              patch%aux%Hydrate%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                den_kg(option%gas_phase) = vec_ptr(local_id)
            enddo
          case(GAS_ENERGY)
            do local_id=1,grid%nlmax
              patch%aux%Hydrate%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                U(option%gas_phase) = vec_ptr(local_id)
            enddo
          case(GAS_MOLE_FRACTION)
            do local_id=1,grid%nlmax
              patch%aux%Hydrate%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                xmol(isubvar,option%gas_phase) = vec_ptr(local_id)
            enddo
          case(HYDRATE_SATURATION)
            do local_id=1,grid%nlmax
              patch%aux%Hydrate%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                sat(option%hydrate_phase) = vec_ptr(local_id)
            enddo
          case(ICE_SATURATION)
            do local_id=1,grid%nlmax
              patch%aux%Hydrate%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                sat(option%ice_phase) = vec_ptr(local_id)
            enddo
        end select
      else if (associated(patch%aux%WIPPFlo)) then
        select case(ivar)
          case(LIQUID_PRESSURE)
            do local_id=1,grid%nlmax
              patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                pres(option%liquid_phase) = vec_ptr(local_id)
            enddo
          case(GAS_PRESSURE)
            do local_id=1,grid%nlmax
              patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                pres(option%gas_phase) = vec_ptr(local_id)
            enddo
          case(CAPILLARY_PRESSURE)
            do local_id=1,grid%nlmax
              patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                pres(option%capillary_pressure_id) = vec_ptr(local_id)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                sat(option%liquid_phase) = vec_ptr(local_id)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
               den_kg(option%liquid_phase) = vec_ptr(local_id)
            enddo
          case(GAS_SATURATION)
            do local_id=1,grid%nlmax
              patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                sat(option%gas_phase) = vec_ptr(local_id)
            enddo
          case(GAS_DENSITY)
            do local_id=1,grid%nlmax
              patch%aux%WIPPFlo%auxvars(ZERO_INTEGER,grid%nL2G(local_id))% &
                den_kg(option%gas_phase) = vec_ptr(local_id)
            enddo
        end select
      endif

    ! NUCLEAR_WASTE_TRANSPORT:
    case(TOTAL_BULK_CONC,AQUEOUS_EQ_CONC,MNRL_EQ_CONC,SORB_EQ_CONC, &
         MNRL_VOLUME_FRACTION)
      select case(ivar)
        case(TOTAL_BULK_CONC)
          call PrintErrMsg(option,'Setting of TOTAL_BULK_CONC at grid cell &
                           &not yet supported.')
        case(AQUEOUS_EQ_CONC)
          call PrintErrMsg(option,'Setting of AQUEOUS_EQ_CONC at grid cell &
                           &not yet supported.')
        case(MNRL_EQ_CONC)
          call PrintErrMsg(option,'Setting of MNRL_EQ_CONC at grid cell &
                           &not yet supported.')
        case(SORB_EQ_CONC)
          call PrintErrMsg(option,'Setting of SORB_EQ_CONC at grid cell &
                           &not yet supported.')
        case(MNRL_VOLUME_FRACTION)
          call PrintErrMsg(option,'Setting of MNRL_VOLUME_FRACTION at grid &
                           &cell not yet supported.')
      end select

    case(PRIMARY_MOLALITY,TOTAL_MOLARITY,MINERAL_VOLUME_FRACTION, &
         PRIMARY_ACTIVITY_COEF,SECONDARY_ACTIVITY_COEF,IMMOBILE_SPECIES, &
         GAS_CONCENTRATION,REACTION_AUXILIARY,MINERAL_SURFACE_AREA)
      select case(ivar)
        case(PRIMARY_MOLALITY)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%auxvars(grid%nL2G(local_id))%pri_molal(isubvar) = &
                vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%auxvars(ghosted_id)%pri_molal(isubvar) = &
                vec_ptr(ghosted_id)
            enddo
          endif
        case(TOTAL_MOLARITY)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%auxvars(grid%nL2G(local_id))% &
                total(isubvar,iphase) = vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%auxvars(ghosted_id)% &
                total(isubvar,iphase) = vec_ptr(ghosted_id)
            enddo
          endif
        case(GAS_CONCENTRATION)
          option%io_buffer = 'Active gas concentrations cannot be set in &
            &PatchSetVariable.'
          call PrintErrMsg(option)
        case(MINERAL_VOLUME_FRACTION)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%auxvars(grid%nL2G(local_id))% &
                mnrl_volfrac(isubvar) = vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%auxvars(ghosted_id)% &
                mnrl_volfrac(isubvar) = vec_ptr(ghosted_id)
            enddo
          endif
        case(MINERAL_SURFACE_AREA)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%auxvars(grid%nL2G(local_id))% &
                mnrl_area(isubvar) = vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%auxvars(ghosted_id)% &
                mnrl_area(isubvar) = vec_ptr(ghosted_id)
            enddo
          endif
        case(IMMOBILE_SPECIES)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%auxvars(grid%nL2G(local_id))% &
                immobile(isubvar) = vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%auxvars(ghosted_id)% &
                immobile(isubvar) = vec_ptr(ghosted_id)
            enddo
          endif
        case(PRIMARY_ACTIVITY_COEF)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%auxvars(grid%nL2G(local_id))% &
                pri_act_coef(isubvar) = vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%auxvars(ghosted_id)% &
                pri_act_coef(isubvar) = vec_ptr(ghosted_id)
            enddo
          endif
        case(SECONDARY_ACTIVITY_COEF)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%auxvars(grid%nL2G(local_id))% &
                sec_act_coef(isubvar) = vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%auxvars(ghosted_id)% &
                sec_act_coef(isubvar) = vec_ptr(ghosted_id)
            enddo
          endif
        case(REACTION_AUXILIARY)
          if (vec_format == GLOBAL) then
            do local_id=1,grid%nlmax
              patch%aux%RT%auxvars(grid%nL2G(local_id))% &
                auxiliary_data(isubvar) = vec_ptr(local_id)
            enddo
          else if (vec_format == LOCAL) then
            do ghosted_id=1,grid%ngmax
              patch%aux%RT%auxvars(ghosted_id)% &
                auxiliary_data(isubvar) = vec_ptr(ghosted_id)
            enddo
          endif
      end select
    case(PRIMARY_MOLARITY,TOTAL_MOLALITY, &
         SECONDARY_MOLALITY,SECONDARY_MOLARITY, &
         COLLOID_MOBILE,COLLOID_IMMOBILE)
      select case(ivar)
        case(PRIMARY_MOLARITY)
          call PrintErrMsg(option,'Setting of primary molarity at grid cell not supported.')
        case(SECONDARY_MOLALITY)
          call PrintErrMsg(option,'Setting of secondary molality at grid cell not supported.')
        case(SECONDARY_MOLARITY)
          call PrintErrMsg(option,'Setting of secondary molarity at grid cell not supported.')
        case(TOTAL_MOLALITY)
          call PrintErrMsg(option,'Setting of total molality at grid cell not supported.')
        case(COLLOID_MOBILE)
          call PrintErrMsg(option,'Setting of mobile colloid concentration at grid cell not supported.')
        case(COLLOID_IMMOBILE)
          call PrintErrMsg(option,'Setting of immobile colloid concentration at grid cell not supported.')
      end select
    case(POROSITY,BASE_POROSITY,INITIAL_POROSITY,ELECTRICAL_CONDUCTIVITY)
      if (vec_format == GLOBAL) then
        do local_id=1,grid%nlmax
          call MaterialAuxVarSetValue(material_auxvars(grid%nL2G(local_id)), &
                                      ivar,vec_ptr(local_id))
        enddo
      else if (vec_format == LOCAL) then
        do ghosted_id=1,grid%ngmax
          call MaterialAuxVarSetValue(material_auxvars(ghosted_id), &
                                      ivar,vec_ptr(ghosted_id))
        enddo
      endif
    case(STATE,PHASE)
      do local_id=1,grid%nlmax
        patch%aux%Global%auxvars(grid%nL2G(local_id))%istate = &
          int(vec_ptr(local_id)+1.d-10)
      enddo
    case(VOLUME,TORTUOSITY,SOIL_COMPRESSIBILITY,SOIL_REFERENCE_PRESSURE)
      option%io_buffer = 'Setting of volume, tortuosity, ' // &
        'soil compressibility or soil reference pressure in ' // &
        '"PatchSetVariable" not supported.'
      call PrintErrMsg(option)
    case(PERMEABILITY,PERMEABILITY_X,PERMEABILITY_Y,PERMEABILITY_Z, &
         GAS_PERMEABILITY,GAS_PERMEABILITY_X,GAS_PERMEABILITY_Y, &
         GAS_PERMEABILITY_Z)
      option%io_buffer = 'Setting of permeability in "PatchSetVariable"' // &
        ' not supported.'
      call PrintErrMsg(option)
    case(MATERIAL_ID)
      !geh: this would require the creation of a permanent mapping between
      !     external and internal material ids, which we want to avoid.
      call PrintErrMsg(option, &
                       'Cannot set MATERIAL_ID through PatchSetVariable()')
      if (vec_format == GLOBAL) then
        do local_id=1,grid%nlmax
          patch%imat(grid%nL2G(local_id)) = int(vec_ptr(local_id))
        enddo
      else if (vec_format == LOCAL) then
        patch%imat(1:grid%ngmax) = int(vec_ptr(1:grid%ngmax))
      endif
    case(ELECTRICAL_POTENTIAL)
      do local_id=1,grid%nlmax
        patch%aux%ERT%auxvars(grid%nL2G(local_id))%potential(isubvar) = &
          vec_ptr(local_id)
      enddo
    case(ELECTRICAL_JACOBIAN)
      do local_id=1,grid%nlmax
        patch%aux%ERT%auxvars(grid%nL2G(local_id))%jacobian(isubvar) = &
          vec_ptr(local_id)
      enddo
    case(PROCESS_ID)
      call PrintErrMsg(option, &
                       'Cannot set PROCESS_ID through PatchSetVariable()')
    case(NATURAL_ID)
      call PrintErrMsg(option, &
                       'Cannot set NATURAL_ID through PatchSetVariable()')
    case default
      write(option%io_buffer, &
            '(''IVAR ('',i3,'') not found in PatchSetVariable'')') ivar
      call PrintErrMsg(option)
  end select

  call VecRestoreArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine PatchSetVariable

! ************************************************************************** !

subroutine PatchCountCells(patch,total_count,active_count)
  !
  ! Counts # of active and inactive grid cells
  !
  ! Author: Glenn Hammond
  ! Date: 06/01/10
  !

  use Option_module

  implicit none

  type(patch_type) :: patch
  PetscInt :: total_count
  PetscInt :: active_count

  type(grid_type), pointer :: grid
  PetscInt :: local_id

  grid => patch%grid

  total_count = grid%nlmax

  active_count = 0
  do local_id = 1, grid%nlmax
    if (patch%imat(grid%nL2G(local_id)) <= 0) cycle
    active_count = active_count + 1
  enddo

end subroutine PatchCountCells

! ************************************************************************** !

subroutine PatchCalculateCFL1Timestep(patch,option,max_dt_cfl_1, &
                                      max_pore_velocity)
  !
  ! Calculates largest time step to preserves a
  ! CFL # of 1 in a patch
  !
  ! Author: Glenn Hammond
  ! Date: 10/06/11
  !

  use Option_module
  use Connection_module
  use Coupler_module
  use Field_module
  use Global_Aux_module
  use Material_Aux_module

  implicit none

  type(patch_type) :: patch
  type(option_type) :: option
  PetscReal :: max_dt_cfl_1
  PetscReal :: max_pore_velocity

  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: por_sat_ave, por_sat_min, v_darcy, v_pore_ave, v_pore_max
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: iphase
  PetscReal :: tempreal(2)
  PetscReal :: dt_cfl_1
  PetscErrorCode :: ierr

  field => patch%field
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  grid => patch%grid

  max_dt_cfl_1 = MAX_DOUBLE
  max_pore_velocity = 0.d0

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
      distance = cur_connection_set%dist(0,iconn)
      fraction_upwind = cur_connection_set%dist(-1,iconn)
      do iphase = 1, option%nphase
        ! if the phase is not present in either cell, skip the connection
        if (.not.(global_auxvars(ghosted_id_up)%sat(iphase) > 0.d0 .and. &
                  global_auxvars(ghosted_id_dn)%sat(iphase) > 0.d0)) cycle
        por_sat_min = min(material_auxvars(ghosted_id_up)%porosity* &
                          global_auxvars(ghosted_id_up)%sat(iphase), &
                          material_auxvars(ghosted_id_dn)%porosity* &
                          global_auxvars(ghosted_id_dn)%sat(iphase))
        por_sat_ave = (fraction_upwind* &
                       material_auxvars(ghosted_id_up)%porosity* &
                       global_auxvars(ghosted_id_up)%sat(iphase) + &
                      (1.d0-fraction_upwind)* &
                      material_auxvars(ghosted_id_dn)%porosity* &
                      global_auxvars(ghosted_id_dn)%sat(iphase))
        v_darcy = patch%internal_velocities(iphase,sum_connection)
        v_pore_max = v_darcy / por_sat_min
        v_pore_ave = v_darcy / por_sat_ave
        !geh: I use v_pore_max to ensure that we limit the cfl based on the
        !     highest velocity through the face.  If porosity*saturation
        !     varies, the pore water velocity will be highest on the side
        !     of the face with the smalled value of porosity*saturation.
        dt_cfl_1 = distance / dabs(v_pore_max)
        max_dt_cfl_1 = min(dt_cfl_1,max_dt_cfl_1)
        max_pore_velocity = max(v_pore_max,max_pore_velocity)
      enddo
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id_dn = cur_connection_set%id_dn(iconn)
      ghosted_id_dn = grid%nL2G(local_id_dn)
      if (patch%imat(ghosted_id_dn) <= 0) cycle
      !geh: since on boundary, dist must be scaled by 2.d0
      distance = 2.d0*cur_connection_set%dist(0,iconn)
      do iphase = 1, option%nphase
        ! the _ave variable is being reused. it is actually, max
        por_sat_ave = material_auxvars(ghosted_id_dn)%porosity* &
                      global_auxvars(ghosted_id_dn)%sat(iphase)
        v_darcy = patch%boundary_velocities(iphase,sum_connection)
        v_pore_ave = v_darcy / por_sat_ave
        dt_cfl_1 = distance / dabs(v_pore_ave)
        max_dt_cfl_1 = min(dt_cfl_1,max_dt_cfl_1)
        max_pore_velocity = max(v_pore_ave,max_pore_velocity)
      enddo
    enddo
    boundary_condition => boundary_condition%next
  enddo

  tempreal(1) = max_dt_cfl_1
  tempreal(2) = -max_pore_velocity
  call MPI_Allreduce(MPI_IN_PLACE,tempreal,TWO_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_MIN, &
                     option%mycomm,ierr)

  max_dt_cfl_1 = tempreal(1)
  max_pore_velocity = -tempreal(2)

end subroutine PatchCalculateCFL1Timestep

! ************************************************************************** !

function PatchGetVarNameFromKeyword(keyword,option)
  !
  ! Returns the name of variable defined by keyword
  !
  ! Author: Glenn Hammond
  ! Date: 07/28/11
  !

  use Option_module

  implicit none

  character(len=MAXWORDLENGTH) :: keyword
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: PatchGetVarNameFromKeyword
  character(len=MAXSTRINGLENGTH) :: var_name

  select case(keyword)
    case('PROCESS_ID')
      var_name = 'Processor ID'
    case('NATURAL_ID')
      var_name = 'Natural ID'
    case default
      option%io_buffer = 'Keyword "' // trim(keyword) // '" not ' // &
                         'recognized in PatchGetIvarsFromKeyword()'
      call PrintErrMsg(option)
  end select

  PatchGetVarNameFromKeyword = var_name

end function PatchGetVarNameFromKeyword

! ************************************************************************** !

subroutine PatchGetIvarsFromKeyword(keyword,ivar,isubvar,var_type,option)
  !
  ! Returns the ivar and isubvars for extracting
  ! datasets using PatchGet/PatchSet routines
  !
  ! Author: Glenn Hammond
  ! Date: 07/28/11
  !

  use Option_module
  use Variables_module

  implicit none

  character(len=MAXWORDLENGTH) :: keyword
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: var_type
  type(option_type) :: option

  select case(keyword)
    case('PROCESS_ID')
      ivar = PROCESS_ID
      isubvar = ZERO_INTEGER
      var_type = INT_VAR
    case('NATURAL_ID')
      ivar = NATURAL_ID
      isubvar = ZERO_INTEGER
      var_type = INT_VAR
    case default
      option%io_buffer = 'Keyword "' // trim(keyword) // '" not ' // &
                         'recognized in PatchGetIvarsFromKeyword()'
      call PrintErrMsg(option)
  end select

end subroutine

! ************************************************************************** !

subroutine PatchGetVariable2(patch,option,output_option,vec, &
                             ivar,isubvar,isubvar2)
  !
  ! PatchGetVariable: Extracts variables indexed by ivar and isubvar from a patch
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/08
  !

  use Grid_module
  use Option_module
  use Output_Aux_module
  use Variables_module

  implicit none

  type(option_type), pointer :: option
  !class(reaction_rt_type), pointer :: reaction
  type(output_option_type), pointer :: output_option
  type(patch_type), pointer :: patch
  Vec :: vec
  PetscInt :: ivar
  PetscInt :: isubvar
  PetscInt :: isubvar2
  PetscInt :: iphase

  PetscInt :: local_id, ghosted_id
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscReal :: xmass
  PetscReal :: tempreal
  PetscInt :: tempint
  PetscInt :: irate, istate, irxn
  PetscErrorCode :: ierr

  grid => patch%grid

  call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)

  iphase = 1

  select case(ivar)
    case(MATERIAL_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = &
          patch%imat_internal_to_external(abs(patch%imat(grid%nL2G(local_id))))
      enddo
    case(PROCESS_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = option%myrank
      enddo
    case(NATURAL_ID)
      do local_id=1,grid%nlmax
        vec_ptr(local_id) = grid%nG2A(grid%nL2G(local_id))
      enddo
    case default
      write(option%io_buffer, &
            '(''IVAR ('',i3,'') not found in PatchGetVariable'')') ivar
      call PrintErrMsg(option)
  end select

end subroutine PatchGetVariable2

! ************************************************************************** !

subroutine PatchGetCellCenteredVelocities(patch,iphase,velocities)
  !
  ! Calculates the Darcy velocity at the center of all cells in a patch
  !
  ! Author: Glenn Hammond
  ! Date: 01/31/14
  !
  use Connection_module
  use Coupler_module

  implicit none

  type(patch_type), pointer :: patch
  PetscInt :: iphase
  PetscReal, intent(out) :: velocities(:,:)

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: sum_connection, iconn, num_connections
  PetscReal, allocatable :: sum_area(:,:), sum_velocity(:,:)
  PetscReal :: area(3), velocity(3)
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: local_id
  PetscInt :: i

  grid => patch%grid

  allocate(sum_velocity(3,grid%nlmax))
  allocate(sum_area(3,grid%nlmax))
  sum_velocity(:,:) = 0.d0
  sum_area(:,:) = 0.d0

  ! interior velocities
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
      local_id_dn = grid%nG2L(ghosted_id_dn) ! = zero for ghost nodes
      ! velocities are stored as the downwind face of the upwind cell
      area = cur_connection_set%area(iconn)* &
             cur_connection_set%dist(1:3,iconn)
      velocity = patch%internal_velocities(iphase,sum_connection)*area
      if (local_id_up > 0) then
        sum_velocity(:,local_id_up) = sum_velocity(:,local_id_up) + velocity
        sum_area(:,local_id_up) = sum_area(:,local_id_up) + dabs(area)
      endif
      if (local_id_dn > 0) then
        sum_velocity(:,local_id_dn) = sum_velocity(:,local_id_dn) + velocity
        sum_area(:,local_id_dn) = sum_area(:,local_id_dn) + dabs(area)
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! boundary velocities
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      area = cur_connection_set%area(iconn)* &
             cur_connection_set%dist(1:3,iconn)
      velocity = patch%boundary_velocities(iphase,sum_connection)*area
      sum_velocity(:,local_id) = sum_velocity(:,local_id) + velocity
      sum_area(:,local_id) = sum_area(:,local_id) + dabs(area)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! divide by total area
  do local_id=1,grid%nlmax
    do i=1,3
      if (sum_area(i,local_id) > 0.d0) then
        velocities(i,local_id) = sum_velocity(i,local_id) / &
                                 sum_area(i,local_id)
      else
        velocities(i,local_id) = 0.d0
      endif
    enddo
  enddo

  deallocate(sum_velocity)
  deallocate(sum_area)

end subroutine PatchGetCellCenteredVelocities

! ************************************************************************** !

subroutine PatchGetKOrthogonalityError(grid, material_auxvars, error)
  !
  ! Calculates the K orthogonality error at all cells in a patch
  ! Error from https://www.onepetro.org/conference-paper/SPE-37998-MS
  ! Work only for unstructured implicit grid type
  !
  ! Author: Moise Rousseau
  ! Date: 06/28/20

  use Option_module
  use Connection_module
  use Geometry_module
  use Material_Aux_module
  use Variables_module, only : PERMEABILITY_X, PERMEABILITY_Y, &
                               PERMEABILITY_Z, PERMEABILITY_XY, &
                               PERMEABILITY_XZ, PERMEABILITY_YZ
  use Utility_module, only : DotProduct, CrossProduct
  !use Grid_Unstructured_Aux_module

  implicit none

  type(grid_type), pointer :: grid
  type(material_auxvar_type), pointer :: material_auxvars(:)
  PetscReal, pointer :: error(:)
  class (connection_set_type), pointer :: cur_connection_set

  PetscInt :: cell_id, face_id, i, iconn
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
  PetscReal :: v1(3),v2(3),normal(3), num(3), den(3)
  PetscReal :: face_error, normal_magnitude
  type(point3d_type) :: point1, point2, point3
  PetscReal, allocatable :: array_normal(:), array_error(:)
  PetscErrorCode :: ierr
  PetscReal :: kx, kxy, kxz, ky, kyz, kz

  type(option_type) :: option

  !unstructured_grid => patch%grid%unstructured_grid
  cur_connection_set => grid%internal_connection_set_list%first

  allocate(array_normal(grid%nlmax))
  allocate(array_error(grid%nlmax))
  array_normal(:) = 0.0
  array_error(:) = 0.0

  do iconn = 1, cur_connection_set%num_connections
    ! For each face we need face normal and the vecteur between the two center
    ! shared by the face
    face_id = cur_connection_set%face_id(iconn)

    ! Normal computation
    point1 = grid%unstructured_grid%vertices( &
                          grid%unstructured_grid%face_to_vertex(1,face_id))
    point2 = grid%unstructured_grid%vertices( &
                          grid%unstructured_grid%face_to_vertex(2,face_id))
    point3 = grid%unstructured_grid%vertices( &
                          grid%unstructured_grid%face_to_vertex(3,face_id))
    v1(1) = point1%x-point2%x
    v1(2) = point1%y-point2%y
    v1(3) = point1%z-point2%z
    v2(1) = point1%x-point3%x
    v2(2) = point1%y-point3%y
    v2(3) = point1%z-point3%z
    normal = CrossProduct(v1,v2)
    normal_magnitude = sqrt(DotProduct(normal,normal))
    normal = normal / normal_magnitude * grid%unstructured_grid%face_area(face_id)

    !Cell ids from both side of the face
    ghosted_id_up = cur_connection_set%id_up(iconn)
    ghosted_id_dn = cur_connection_set%id_dn(iconn)
    local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
    local_id_dn = grid%nG2L(ghosted_id_dn) ! = zero for ghost nodes

    ! Error for upwind cell
    if (local_id_up > 0) then
      ! get permeability
      kx = MaterialAuxVarGetValue(material_auxvars(local_id_up),PERMEABILITY_X)
      ky = MaterialAuxVarGetValue(material_auxvars(local_id_up),PERMEABILITY_Y)
      kz = MaterialAuxVarGetValue(material_auxvars(local_id_up),PERMEABILITY_Z)
      kxy = MaterialAuxVarGetValue(material_auxvars(local_id_up),PERMEABILITY_XY)
      kxz = MaterialAuxVarGetValue(material_auxvars(local_id_up),PERMEABILITY_XZ)
      kyz = MaterialAuxVarGetValue(material_auxvars(local_id_up),PERMEABILITY_YZ)
      ! error computation
      ! denominator
      den = 0.0
      den(1) = normal(1)*kx + normal(2)*kxy + normal(3)*kxz
      den(2) = normal(1)*kxy + normal(2)*ky + normal(3)*kyz
      den(3) = normal(1)*kxz + normal(2)*kyz + normal(3)*kz
      ! numerator
      num = den - DotProduct(den,cur_connection_set%dist(1:3,iconn)) * &
                                      cur_connection_set%dist(1:3,iconn)
      ! face error
      array_error(local_id_up) = array_error(local_id_up) + sqrt(DotProduct(num,num))
      array_normal(local_id_up) = array_normal(local_id_up) + &
                                      sqrt(DotProduct(den,den))
    endif

    if (local_id_dn > 0) then
      ! get permeability
      kx = MaterialAuxVarGetValue(material_auxvars(local_id_dn),PERMEABILITY_X)
      ky = MaterialAuxVarGetValue(material_auxvars(local_id_dn),PERMEABILITY_Y)
      kz = MaterialAuxVarGetValue(material_auxvars(local_id_dn),PERMEABILITY_Z)
      kxy = MaterialAuxVarGetValue(material_auxvars(local_id_dn),PERMEABILITY_XY)
      kxz = MaterialAuxVarGetValue(material_auxvars(local_id_dn),PERMEABILITY_XZ)
      kyz = MaterialAuxVarGetValue(material_auxvars(local_id_dn),PERMEABILITY_YZ)
      ! error computation
      ! denominator
      den = 0.0
      den(1) = normal(1)*kx + normal(2)*kxy + normal(3)*kxz
      den(2) = normal(1)*kxy + normal(2)*ky + normal(3)*kyz
      den(3) = normal(1)*kxz + normal(2)*kyz + normal(3)*kz
      ! numerator
      num = den - DotProduct(den,cur_connection_set%dist(1:3,iconn)) * &
                                      cur_connection_set%dist(1:3,iconn)
      ! face error
      array_error(local_id_dn) = array_error(local_id_dn) + sqrt(DotProduct(num,num))
      array_normal(local_id_dn) = array_normal(local_id_dn) + &
                                      sqrt(DotProduct(den,den))
    endif
  enddo
  !Compute final error
  do i = 1, grid%nlmax
    error(i) = array_error(i) / array_normal(i)
  enddo

  deallocate(array_normal)
  deallocate(array_error)

end subroutine PatchGetKOrthogonalityError

! ************************************************************************** !

subroutine PatchGetIntegralFluxConnections(patch,integral_flux,option)
  !
  !
  ! Returns a list of internal and boundary connection ids for cell
  ! interfaces defined by an integral flux object.
  !
  ! Author: Glenn Hammond
  ! Date: 10/20/14, 01/31/18
  !
  use Option_module
  use String_module
  use Integral_Flux_module
  use Geometry_module
  use Utility_module
  use Connection_module
  use Coupler_module
  use Grid_Unstructured_Cell_module, only : MAX_VERT_PER_FACE


  implicit none

  type(patch_type) :: patch
  type(integral_flux_type) :: integral_flux
  type(option_type) :: option

  type(point3d_type), pointer :: polygon(:)
  type(plane_type), pointer :: plane
  PetscReal,pointer :: coordinates_and_directions(:,:)
  PetscInt,pointer :: vertices(:,:)
  PetscInt,pointer :: by_cell_ids(:,:)
  type(grid_type), pointer :: grid
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(coupler_type), pointer :: boundary_condition

  PetscInt, pointer :: connections(:)
  PetscInt :: idir
  PetscInt :: icount
  PetscInt :: array_size
  PetscInt :: sum_connection
  PetscInt :: iconn
  PetscInt :: ivert
  PetscInt :: i, ii
  PetscInt :: face_id
  PetscInt :: local_id
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: natural_id_up, natural_id_dn
  PetscReal :: fraction_upwind
  PetscReal :: magnitude
  PetscReal :: v1(3), v2(3), cp(3)
  PetscReal :: x, y, z
  PetscReal :: coord(3)
  PetscReal :: unit_direction(3)
  PetscReal, parameter :: absolute_tolerance = 1.d-10
  PetscReal, parameter :: relative_tolerance = 1.d-6
  PetscBool :: found
  PetscBool :: found2
  PetscBool :: reverse_direction
  PetscBool :: same_direction
  PetscBool, pointer :: yet_to_be_found(:)
  PetscInt :: ipass
  PetscInt :: face_vertices_natural(MAX_VERT_PER_FACE)
  PetscInt :: ivert1, ivert2
  PetscInt :: iv1, iv2
  PetscInt :: num_vertices1, num_vertices2
  PetscInt :: num_to_be_found
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscErrorCode :: ierr

  nullify(yet_to_be_found)
  nullify(connections)

  grid => patch%grid
  polygon => integral_flux%polygon
  plane => integral_flux%plane
  coordinates_and_directions => integral_flux%coordinates_and_directions
  vertices => integral_flux%vertices
  by_cell_ids => integral_flux%cell_ids
  num_to_be_found = 0

  error_string = 'error string missing in PatchGetIntegralFluxConnections'

  if (associated(plane)) then
    error_string = 'plane coincides with an internal cell boundary &
      &or a boundary condition'
  endif

  if (associated(polygon)) then
    error_string = 'polygon coincides with an internal cell boundary &
      &or a boundary condition'
    ! determine orientation of polygon
    allocate(plane)
    if (size(polygon) > 2) then
      call GeometryComputePlaneWithPoints(plane, &
                                   polygon(1)%x,polygon(1)%y,polygon(1)%z, &
                                   polygon(2)%x,polygon(2)%y,polygon(2)%z, &
                                   polygon(3)%x,polygon(3)%y,polygon(3)%z)
    else
      if (Equal(polygon(1)%x,polygon(2)%x) .and. &
          .not. Equal(polygon(1)%y,polygon(2)%y) .and. &
          .not. Equal(polygon(1)%z,polygon(2)%z)) then
        call GeometryComputePlaneWithPoints(plane, &
                                     polygon(1)%x,polygon(1)%y,polygon(1)%z, &
                                     polygon(1)%x,polygon(1)%y,polygon(2)%z, &
                                     polygon(1)%x,polygon(2)%y,polygon(2)%z)
      else if (.not. Equal(polygon(1)%x,polygon(2)%x) .and. &
               Equal(polygon(1)%y,polygon(2)%y) .and. &
               .not. Equal(polygon(1)%z,polygon(2)%z)) then
        call GeometryComputePlaneWithPoints(plane, &
                                     polygon(1)%x,polygon(1)%y,polygon(1)%z, &
                                     polygon(1)%x,polygon(1)%y,polygon(2)%z, &
                                     polygon(2)%x,polygon(1)%y,polygon(2)%z)
      else if (.not. Equal(polygon(1)%x,polygon(2)%x) .and. &
               .not. Equal(polygon(1)%y,polygon(2)%y) .and. &
               Equal(polygon(1)%z,polygon(2)%z)) then
        call GeometryComputePlaneWithPoints(plane, &
                                     polygon(1)%x,polygon(1)%y,polygon(1)%z, &
                                     polygon(1)%x,polygon(2)%y,polygon(1)%z, &
                                     polygon(2)%x,polygon(2)%y,polygon(1)%z)
      else
        if (OptionPrintToScreen(option)) write(*,*) 'pt1: ', &
          polygon(1)%x, polygon(1)%y, polygon(1)%z
        if (OptionPrintToScreen(option)) write(*,*) 'pt2: ', &
          polygon(2)%x, polygon(2)%y, polygon(2)%z
        option%io_buffer = 'An integral flux polygon defined by 2 points must &
          & be a plane.'
        call PrintErrMsg(option)
      endif
    endif
    icount = 0
    if (dabs(plane%A) > absolute_tolerance) then
      idir = 1
      icount = icount + 1
    endif
    if (dabs(plane%B) > absolute_tolerance) then
      idir = 2
      icount = icount + 1
    endif
    if (dabs(plane%C) > absolute_tolerance) then
      idir = 3
      icount = icount + 1
    endif
    if (icount /= 1) then
      option%io_buffer = 'Polygon defined in integral flux "' // &
        trim(adjustl(integral_flux%name)) // &
        '" must be aligned with grid axes.'
      call PrintErrMsg(option)
    endif
    integral_flux%plane => plane
  endif

  if (associated(coordinates_and_directions)) then
    error_string = 'coordinates and directions coincide with an internal &
      &cell boundary or a boundary condition'
    num_to_be_found = size(coordinates_and_directions,2)
    allocate(yet_to_be_found(num_to_be_found))
    yet_to_be_found = PETSC_TRUE
  endif

  if (associated(vertices)) then
    error_string = 'vertices match the face and are ordered properly &
      &(clockwise or counter-clockwise)'
    if (grid%itype /= IMPLICIT_UNSTRUCTURED_GRID) then
      option%io_buffer = 'INTEGRAL_FLUXES defined by VERTICES are only &
        &supported for implicit unstructured grids: ' // &
        trim(integral_flux%name) // '.'
      call PrintErrMsg(option)
    endif
    num_to_be_found = size(vertices,2)
    allocate(yet_to_be_found(num_to_be_found))
    yet_to_be_found = PETSC_TRUE
  endif

  if (associated(by_cell_ids)) then
    error_string = 'cell ids match the an actual face between these cells'
    num_to_be_found = size(by_cell_ids,2)
    allocate(yet_to_be_found(num_to_be_found))
    yet_to_be_found = PETSC_TRUE
  endif

  array_size = 100
  allocate(connections(array_size))

  ipass = 1
  nullify(boundary_condition)
  nullify(cur_connection_set)
  do
    select case(ipass)
      case(1) ! internal connections
        cur_connection_set => grid%internal_connection_set_list%first
        sum_connection = 0
        icount = 0
      case(2) ! boundary connections
        ! sets up first boundary condition in list
        if (.not.associated(boundary_condition)) then
          boundary_condition => patch%boundary_condition_list%first
          sum_connection = 0
          icount = 0
        endif
        if (associated(boundary_condition)) then
          cur_connection_set => boundary_condition%connection_set
        else
          nullify(cur_connection_set)
        endif
    end select
    do
      if (.not.associated(cur_connection_set)) exit
      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        magnitude = cur_connection_set%dist(0,iconn)
        unit_direction(:) = &
          cur_connection_set%dist(X_DIRECTION:Z_DIRECTION,iconn)
        select case(ipass)
          case(1) ! internal connections
            ghosted_id_up = cur_connection_set%id_up(iconn)
            ghosted_id_dn = cur_connection_set%id_dn(iconn)
            local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
            local_id_dn = grid%nG2L(ghosted_id_dn)
            ! natural ids are used to ensure that a connection is only mapped
            ! on one process if running in parallel
            natural_id_up = grid%nG2A(ghosted_id_up)
            natural_id_dn = grid%nG2A(ghosted_id_dn)
            ! if one of the cells is ghosted, the process stores the flux only
            ! when the upwind cell is non-ghosted.
            if (local_id_up <= 0 .or. local_id_dn <= 0) then
              ! we are on a ghost boundary. we cannot rely on local or ghost
              ! numbering to determine on which process to store the
              ! connection as these numbering schemes vary across
              ! processes. use natural numbering and only store if the
              ! on process cell has the lower number
              if ((local_id_up > 0 .and. natural_id_up > natural_id_dn) .or. &
                  (local_id_dn > 0 .and. natural_id_dn > natural_id_up)) then
                ! skip this connection on this process
                cycle
              endif
            endif
            fraction_upwind = cur_connection_set%dist(-1,iconn)
            x = grid%x(ghosted_id_up) + fraction_upwind * magnitude * &
                                        unit_direction(X_DIRECTION)
            y = grid%y(ghosted_id_up) + fraction_upwind * magnitude * &
                                        unit_direction(Y_DIRECTION)
            z = grid%z(ghosted_id_up) + fraction_upwind * magnitude * &
                                        unit_direction(Z_DIRECTION)
          case(2) ! boundary connections
            local_id = cur_connection_set%id_dn(iconn)
            ghosted_id = grid%nL2G(local_id)
            fraction_upwind = 1.d0
            x = grid%x(ghosted_id) - fraction_upwind * magnitude * &
                                     unit_direction(X_DIRECTION)
            y = grid%y(ghosted_id) - fraction_upwind * magnitude * &
                                     unit_direction(Y_DIRECTION)
            z = grid%z(ghosted_id) - fraction_upwind * magnitude * &
                                     unit_direction(Z_DIRECTION)
        end select
        found = PETSC_FALSE
        same_direction = PETSC_TRUE
        if (associated(plane)) then
          if (minval(unit_direction) < 0.d0) then
            same_direction = PETSC_FALSE
          endif
          found = dabs(GeomComputeDistanceFromPlane(plane,x,y,z)) < &
                  relative_tolerance
          if (found .and. associated(polygon)) then
            found = GeometryPointInPolygon(x,y,z,iabs(idir),polygon)
          endif
        endif
        if (.not.found .and. associated(coordinates_and_directions)) then
          do i = 1, num_to_be_found
            if (.not.yet_to_be_found(i)) cycle
            v1(1) = coordinates_and_directions(1,i) - x
            v1(2) = coordinates_and_directions(2,i) - y
            v1(3) = coordinates_and_directions(3,i) - z
            ! same point?
            if (DotProduct(v1,v1)/magnitude < relative_tolerance) then
              v2(1) = coordinates_and_directions(4,i)
              v2(2) = coordinates_and_directions(5,i)
              v2(3) = coordinates_and_directions(6,i)
              ! same direction?
              cp = CrossProduct(v2,unit_direction)
              if (DotProduct(cp,cp)/magnitude < relative_tolerance) then
                found = PETSC_TRUE
                yet_to_be_found(i) = PETSC_FALSE
                ! could be opposite direction
                if (minval(v2*unit_direction) < 0.d0) then
                  same_direction = PETSC_FALSE
                endif
                exit
              endif
            endif
          enddo
        endif
        if (.not.found .and. associated(vertices)) then
          select case(ipass)
            case(1) ! internal connections
              face_id = grid%unstructured_grid%connection_to_face(iconn)
            case(2) ! boundary connections
              face_id = boundary_condition%connection_set%face_id(iconn)
          end select
          do ivert = 1, MAX_VERT_PER_FACE
            face_vertices_natural(ivert) = &
              grid%unstructured_grid% &
                vertex_ids_natural(grid%unstructured_grid% &
                                     face_to_vertex(ivert,face_id))
          enddo
          num_vertices1 = MAX_VERT_PER_FACE
          do
            if (face_vertices_natural(num_vertices1) > 0) exit
            num_vertices1 = num_vertices1 - 1
          enddo
          do i = 1, num_to_be_found
            if (.not.yet_to_be_found(i)) cycle
            num_vertices2 = size(vertices,1)
            do
              if (Initialized(vertices(num_vertices2,i))) exit
              num_vertices2 = num_vertices2 - 1
            enddo
            if (num_vertices1 /= num_vertices2) cycle
            found2 = PETSC_FALSE
            ! find the first vertex
            do ivert1 = 1, num_vertices1
              do ivert2 = 1, num_vertices2
                if (face_vertices_natural(ivert1) == vertices(ivert2,i)) then
                  found2 = PETSC_TRUE
                  exit
                endif
              enddo
              if (found2) exit
            enddo
            if (found2) then
              reverse_direction = PETSC_FALSE
              ! search forward direction
              iv1 = ivert1
              iv2 = ivert2
              do ii = 1, num_vertices1
                if (face_vertices_natural(iv1) /= vertices(iv2,i)) then
                  found2 = PETSC_FALSE
                  exit
                endif
                iv1 = iv1 + 1
                if (iv1 > num_vertices1) iv1 = 1
                iv2 = iv2 + 1
                if (iv2 > num_vertices1) iv2 = 1
              enddo
              ! search backward direction
              if (.not.found2) then
                reverse_direction = PETSC_TRUE
                found2 = PETSC_TRUE
                iv1 = ivert1
                iv2 = ivert2
                do ii = 1, num_vertices1
                  if (face_vertices_natural(iv1) /= vertices(iv2,i)) then
                    found2 = PETSC_FALSE
                    exit
                  endif
                  iv1 = iv1 + 1
                  if (iv1 > num_vertices1) iv1 = 1
                  iv2 = iv2 - 1
                  if (iv2 < 1) iv2 = num_vertices1
                enddo
              endif
            endif
            if (found2) then
              yet_to_be_found(i) = PETSC_FALSE
              found = PETSC_TRUE
              if (reverse_direction) same_direction = PETSC_FALSE
              exit
            endif
          enddo
        endif
        if (.not.found .and. associated(by_cell_ids)) then
          select case(ipass)
            case(1) ! internal connections
              do i = 1, num_to_be_found
                if (natural_id_dn == by_cell_ids(1,i) .and. &
                    natural_id_up == by_cell_ids(2,i)) then
                  yet_to_be_found(i) = PETSC_FALSE
                  same_direction = PETSC_FALSE
                  found = PETSC_TRUE
                  exit
                elseif (natural_id_up == by_cell_ids(1,i) .and. &
                        natural_id_dn == by_cell_ids(2,i)) then
                  yet_to_be_found(i) = PETSC_FALSE
                  found = PETSC_TRUE
                  exit
                endif
              enddo
            case(2) ! boundary connections
              ! not yet supported
          end select
        endif
        if (found) then
          icount = icount + 1
          if (icount > size(connections)) then
            call ReallocateArray(connections,array_size)
          endif
          if (same_direction) then
            connections(icount) = sum_connection
          else
            connections(icount) = -sum_connection
          endif
        endif
      enddo
      cur_connection_set => cur_connection_set%next
    enddo
    select case(ipass)
      case(1) ! internal connections
        if (icount > 0) then
          allocate(integral_flux%internal_connections(icount))
          integral_flux%internal_connections = connections(1:icount)
        endif
        icount = 0
        ipass = ipass + 1
      case(2) ! boundary connections
        if (associated(boundary_condition)) then
          boundary_condition => boundary_condition%next
        endif
        if (.not.associated(boundary_condition)) then
          if (icount > 0) then
            allocate(integral_flux%boundary_connections(icount))
            integral_flux%boundary_connections = connections(1:icount)
          endif
          exit
        endif
    end select
  enddo

  call DeallocateArray(yet_to_be_found)
  call DeallocateArray(connections)

  icount = 0
  if (associated(integral_flux%internal_connections)) then
    icount = icount + size(integral_flux%internal_connections)
  endif
  if (associated(integral_flux%boundary_connections)) then
    icount = icount + size(integral_flux%boundary_connections)
  endif
  call MPI_Allreduce(MPI_IN_PLACE,icount,ONE_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_SUM,option%mycomm,ierr)
  if (icount == 0) then
    option%io_buffer = 'Zero connections found for INTEGRAL_FLUX "' // &
      trim(adjustl(integral_flux%name)) // &
      '".  Please ensure that the ' // trim(error_string) // '.'
    call PrintErrMsg(option)
  else if (num_to_be_found > 0 .and. icount /= num_to_be_found) then
    if (icount > num_to_be_found) then
      option%io_buffer = trim(StringWrite(icount-num_to_be_found)) // &
        ' extra face(s) [beyond ' // trim(StringWrite(num_to_be_found)) // &
        '] found for INTEGRAL_FLUX "' // &
        trim(adjustl(integral_flux%name)) // &
        '".  There is likely an issue with the parallel implementation. &
        &Please email pflotran-dev@googlegroups.com.'
    else
      option%io_buffer = trim(StringWrite(num_to_be_found-icount)) // &
        ' face(s) of ' // trim(StringWrite(num_to_be_found)) // &
        ' missed for INTEGRAL_FLUX "' // &
        trim(adjustl(integral_flux%name)) // &
        '".  Please ensure that the ' // trim(error_string) // '.'
    endif
    call PrintErrMsg(option)
  else
    write(option%io_buffer,*) icount
    option%io_buffer = trim(StringWrite(icount)) // ' connections found &
      &for integral flux "' // trim(adjustl(integral_flux%name)) // '".'
    call PrintMsg(option)
  endif

end subroutine PatchGetIntegralFluxConnections

! **************************************************************************** !

subroutine PatchCouplerInputRecord(patch)
  !
  ! Prints ingested coupler information to the input record file.
  !
  ! Author: Jenn Frederick
  ! Date: 04/18/2016
  !
  use Coupler_module

  implicit none

  type(patch_type), pointer :: patch

  type(coupler_type), pointer :: cur_coupler
  character(len=MAXWORDLENGTH) :: word1, word2
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: k
  PetscInt :: id = INPUT_RECORD_UNIT

  k = 0

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'INITIAL CONDITIONS'

  ! Initial conditions
  cur_coupler => patch%initial_condition_list%first
  do
    if (.not.associated(cur_coupler)) exit
    k = k + 1
    write(id,'(a29)',advance='no') 'initial condition listed: '
    write(word1,*) k
    write(id,'(a)') '#' // adjustl(trim(word1))
    write(id,'(a29)',advance='no') 'applies to region: '
    write(id,'(a)') adjustl(trim(cur_coupler%region_name))
    if (len_trim(cur_coupler%flow_condition_name) > 0) then
      write(id,'(a29)',advance='no') 'flow condition name: '
      write(id,'(a)') adjustl(trim(cur_coupler%flow_condition_name))
    endif
    if (len_trim(cur_coupler%tran_condition_name) > 0) then
      write(id,'(a29)',advance='no') 'transport condition name: '
      write(id,'(a)') adjustl(trim(cur_coupler%tran_condition_name))
    endif
    write(id,'(a29)') '---------------------------: '
    cur_coupler => cur_coupler%next
  enddo

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'BOUNDARY CONDITIONS'

  ! Boundary conditions
  cur_coupler => patch%boundary_condition_list%first
  do
    if (.not.associated(cur_coupler)) exit
    write(id,'(a29)',advance='no') 'boundary condition name: '
    write(id,'(a)') adjustl(trim(cur_coupler%name))
    write(id,'(a29)',advance='no') 'applies to region: '
    write(id,'(a)') adjustl(trim(cur_coupler%region_name))
    if (len_trim(cur_coupler%flow_condition_name) > 0) then
      write(id,'(a29)',advance='no') 'flow condition name: '
      write(id,'(a)') adjustl(trim(cur_coupler%flow_condition_name))
    endif
    if (len_trim(cur_coupler%tran_condition_name) > 0) then
      write(id,'(a29)',advance='no') 'transport condition name: '
      write(id,'(a)') adjustl(trim(cur_coupler%tran_condition_name))
    endif
    write(id,'(a29)') '---------------------------: '
    cur_coupler => cur_coupler%next
  enddo

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'SOURCE-SINKS'

  ! Source-Sink conditions
  cur_coupler => patch%source_sink_list%first
  do
    if (.not.associated(cur_coupler)) exit
    write(id,'(a29)',advance='no') 'source-sink name: '
    write(id,'(a)') adjustl(trim(cur_coupler%name))
    write(id,'(a29)',advance='no') 'applies to region: '
    write(id,'(a)') adjustl(trim(cur_coupler%region_name))
    if (len_trim(cur_coupler%flow_condition_name) > 0) then
      write(id,'(a29)',advance='no') 'flow condition name: '
      write(id,'(a)') adjustl(trim(cur_coupler%flow_condition_name))
    endif
    if (len_trim(cur_coupler%tran_condition_name) > 0) then
      write(id,'(a29)',advance='no') 'transport condition name: '
      write(id,'(a)') adjustl(trim(cur_coupler%tran_condition_name))
    endif
    write(id,'(a29)') '---------------------------: '
    cur_coupler => cur_coupler%next
  enddo

end subroutine PatchCouplerInputRecord

! **************************************************************************** !

subroutine PatchGetWaterMassInRegion(cell_ids,num_cells,patch,option, &
                                     global_water_mass)
  !
  ! Calculates the water mass in a region in kg
  !
  ! Author: Satish Karra, Heeho Park
  ! Date: 09/20/2016
  !
  use Global_Aux_module
  use Material_Aux_module
  use Grid_module
  use Option_module

  implicit none

  PetscInt, pointer :: cell_ids(:)
  PetscInt :: num_cells
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  PetscReal :: global_water_mass

  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  PetscReal :: m3_water, kg_water
  PetscInt :: k, j, m
  PetscInt :: local_id, ghosted_id
  PetscErrorCode :: ierr
  PetscReal :: local_water_mass

  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  local_water_mass = 0.d0
  global_water_mass = 0.d0

  ! Loop through all cells in the region:
  do k = 1,num_cells
    local_id = cell_ids(k)
    ghosted_id = patch%grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    m3_water = material_auxvars(ghosted_id)%porosity * &         ! [-]
               global_auxvars(ghosted_id)%sat(LIQUID_PHASE) * &  ! [water %]
               material_auxvars(ghosted_id)%volume               ! [m^3-bulk]
    kg_water = m3_water*global_auxvars(ghosted_id)% &            ! [m^3-water]
               den_kg(LIQUID_PHASE)                              ! [kg/m^3-water]
    local_water_mass = local_water_mass + kg_water
  enddo ! Cell loop

  ! Sum the local_water_mass across all processes that own the region:
  call MPI_Allreduce(local_water_mass,global_water_mass,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)

end subroutine PatchGetWaterMassInRegion

! **************************************************************************** !

subroutine PatchGetCompMassInRegionAssign(region_list, &
           mass_balance_region_list,option)
  !
  ! Assigns patch%region information to the mass balance region object
  !
  ! Author: Jenn Frederick
  ! Date: 04/26/2016
  !
  use Output_Aux_module
  use Region_module
  use String_module

  implicit none

  type(region_list_type), pointer :: region_list
  type(mass_balance_region_type), pointer :: mass_balance_region_list
  type(option_type), pointer :: option

  type(region_type), pointer :: cur_region
  type(mass_balance_region_type), pointer :: cur_mbr
  PetscBool :: success

  cur_mbr => mass_balance_region_list
  do
    if (.not.associated(cur_mbr)) exit
    ! Loop through patch%region_list to find wanted region:
    cur_region => region_list%first
    do
      if (.not.associated(cur_region)) exit
      success = PETSC_TRUE
      if (StringCompareIgnoreCase(cur_region%name,cur_mbr%region_name)) exit
      success = PETSC_FALSE
      cur_region => cur_region%next
    enddo
    ! If the wanted region was not found, throw an error msg:
    if (.not.success) then
      option%io_buffer = 'Region ' // trim(cur_mbr%region_name) // ' not &
                          &found among listed regions.'
      call PrintErrMsg(option)
    endif
    ! Assign the mass balance region the wanted region's info:
    cur_mbr%num_cells = cur_region%num_cells
    cur_mbr%region_cell_ids => cur_region%cell_ids
    ! Go to next mass balance region
    cur_mbr => cur_mbr%next
  enddo

end subroutine PatchGetCompMassInRegionAssign

! ************************************************************************** !

subroutine PatchVerifyDatasetGriddedForFlux(dataset,coupler,option)
  !
  ! Verifies that a dataset being used to define fluxes adheres to
  ! several rules that attempt to minimize mass balance error.
  !
  ! Author: Jennifer Frederick
  ! Date: 11/04/2016
  !

  use Option_module
  use Coupler_module
  use Dataset_Gridded_HDF5_class

  implicit none

  class(dataset_gridded_hdf5_type) :: dataset
  type(coupler_type) :: coupler
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscInt :: i, dataset_size

  ! check if dataset is cell-centered:
  if (.not.dataset%is_cell_centered) then
    option%io_buffer = 'Dataset ' // trim(dataset%hdf5_dataset_name) // &
      " must be cell-centered for fluxes. You must set attribute: &
      &h5grp.attrs['Cell Centered'] = True."
    call PrintErrMsg(option)
  endif
  ! check if the dimensions match:
  dataset_size = 1
  do i=1,size(dataset%dims)
    dataset_size = dataset%dims(i)*dataset_size
  enddo
  if (coupler%connection_set%num_connections /= dataset_size) then
    write(string,*) dataset%dims
    write(string2,*) coupler%connection_set%num_connections
    option%io_buffer = 'Dataset ' // trim(dataset%hdf5_dataset_name) // &
      " must have a value for each cell on the boundary defined by &
      &REGION " // trim(coupler%region%name) // '. The dataset dimension &
      &is ' // adjustl(trim(string)) // ' but the number of boundary &
      &connections is ' // adjustl(trim(string2)) // '.'
    call PrintErrMsg(option)
  endif
  ! check if the interpolation method is STEP:
  if (.not.dataset%interpolation_method == INTERPOLATION_STEP) then
    option%io_buffer = 'Dataset ' // trim(dataset%hdf5_dataset_name) // &
      " must be assigned the STEP interpolation method for fluxes. You &
      &must set attribute: h5grp.attrs['Interpolation Method'] = &
      &np.string_('STEP')."
    call PrintErrMsg(option)
  endif

end subroutine PatchVerifyDatasetGriddedForFlux

! ************************************************************************** !

subroutine PatchSetupUpwindDirection(patch,option)
  !
  ! Sets up upwind direction arrays
  !
  ! Author: Glenn Hammond
  ! Date: 10/22/18
  !
  use Grid_module
  use Coupler_module
  use Connection_module
  use Option_module
  use Upwind_Direction_module

  implicit none

  type(patch_type) :: patch
  type(option_type) :: option

  type(grid_type), pointer :: grid
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(coupler_type), pointer :: boundary_condition
  PetscInt, pointer :: upwind_direction(:,:)
  PetscInt, pointer :: upwind_direction_bc(:,:)
  PetscInt :: sum_connection
  PetscInt :: i, iconn
  PetscReal :: dist(3)

  grid => patch%grid

  call UpwindDirectionInit()

  ! allocate arrays for upwind directions
  ! internal connections
  connection_set_list => grid%internal_connection_set_list
  sum_connection = ConnectionGetNumberInList(connection_set_list)
  allocate(upwind_direction(option%nphase,sum_connection))
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
      upwind_direction(:,sum_connection) = i
    enddo
    cur_connection_set => cur_connection_set%next
  enddo
  ! boundary connections
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_condition_list)
  boundary_condition => patch%boundary_condition_list%first
  allocate(upwind_direction_bc(option%nphase,sum_connection))
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
      upwind_direction_bc(:,sum_connection) = i
    enddo
    boundary_condition => boundary_condition%next
  enddo

  patch%flow_upwind_direction => upwind_direction
  patch%flow_upwind_direction_bc => upwind_direction_bc

end subroutine PatchSetupUpwindDirection

! ************************************************************************** !

subroutine PatchUnsupportedVariable1(process_name,variable_name,ivar,option)
  !
  ! Prints an error message when an unsupported output variable is requested.
  !
  ! Author: Glenn Hammond
  ! Date: 04/22/19

  implicit none

  character(len=*) :: process_name
  character(len=*) :: variable_name
  PetscInt :: ivar
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word, word2

  if (len_trim(variable_name) > 1 .and. Initialized(ivar)) then
    option%io_buffer = 'Both Variable name and ID passed as initialized in &
      &PatchUnsupportedVariable.'
    call PrintErrMsg(option)
  else if (len_trim(variable_name) > 1) then
    word = trim(variable_name)
    word2= 'variable name'
  else
    write(word,*) ivar
    word = 'ID ' // trim(adjustl(word))
    word2= 'integer ID'
  endif
  if (len_trim(process_name) > 1) then
    string = 'Variable ' // trim(adjustl(word)) // &
      ' not supported for the process model: ' // &
      trim(adjustl(process_name)) // '.'
  else
    string = 'Variable ' // trim(adjustl(word)) // &
      ' not supported for the employed process models.'
  endif
  option%io_buffer = trim(string) // &
    ' Please look at variables.F90 to match this ' // trim(word2) // &
    ' with a variable parameter constant.  If you feel this is in &
    &error, please email this message to pflotran-dev@googlegroups.com.'
  call PrintErrMsg(option)

end subroutine PatchUnsupportedVariable1

! ************************************************************************** !

subroutine PatchUnsupportedVariable2(process_name,ivar,option)
  !
  ! Author: Glenn Hammond
  ! Date: 04/22/19

  implicit none

  character(len=*) :: process_name
  PetscInt :: ivar
  type(option_type) :: option

  call PatchUnsupportedVariable(process_name,'',ivar,option)

end subroutine PatchUnsupportedVariable2

! ************************************************************************** !

subroutine PatchUnsupportedVariable3(process_name,variable_name,option)
  !
  ! Author: Glenn Hammond
  ! Date: 04/22/19

  implicit none

  character(len=*) :: process_name
  character(len=*) :: variable_name
  type(option_type) :: option

  call PatchUnsupportedVariable(process_name,variable_name, &
                                UNINITIALIZED_INTEGER,option)

end subroutine PatchUnsupportedVariable3

! ************************************************************************** !

subroutine PatchUnsupportedVariable4(ivar,option)
  !
  ! Author: Glenn Hammond
  ! Date: 04/22/19

  implicit none

  PetscInt :: ivar
  type(option_type) :: option

  call PatchUnsupportedVariable('','',ivar,option)

end subroutine PatchUnsupportedVariable4

! ************************************************************************** !

subroutine PatchDestroyList(patch_list)
  !
  ! Deallocates a patch list and array of patches
  !
  ! Author: Glenn Hammond
  ! Date: 10/15/07
  !

  implicit none

  type(patch_list_type), pointer :: patch_list

  type(patch_type), pointer :: cur_patch, prev_patch

  if (.not.associated(patch_list)) return

  if (associated(patch_list%array)) deallocate(patch_list%array)
  nullify(patch_list%array)

  cur_patch => patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    prev_patch => cur_patch
    cur_patch => cur_patch%next
    call PatchDestroy(prev_patch)
  enddo

  nullify(patch_list%first)
  nullify(patch_list%last)
  patch_list%num_patch_objects = 0

  deallocate(patch_list)
  nullify(patch_list)

end subroutine PatchDestroyList

! ************************************************************************** !

subroutine PatchDestroy(patch)
  !
  ! Deallocates a patch object
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  use Utility_module, only : DeallocateArray

  implicit none

  type(patch_type), pointer :: patch

  call DeallocateArray(patch%imat)
  call DeallocateArray(patch%imat_internal_to_external)
  call DeallocateArray(patch%cc_id)
  call DeallocateArray(patch%cct_id)
  call DeallocateArray(patch%mtf_id)
  call DeallocateArray(patch%internal_velocities)
  call DeallocateArray(patch%boundary_velocities)
  call DeallocateArray(patch%internal_tran_coefs)
  call DeallocateArray(patch%boundary_tran_coefs)
  call DeallocateArray(patch%internal_flow_fluxes)
  call DeallocateArray(patch%boundary_flow_fluxes)
  call DeallocateArray(patch%ss_flow_fluxes)
  call DeallocateArray(patch%internal_tran_fluxes)
  call DeallocateArray(patch%boundary_tran_fluxes)
  call DeallocateArray(patch%ss_tran_fluxes)
  call DeallocateArray(patch%ss_flow_vol_fluxes)

  call DeallocateArray(patch%flow_upwind_direction)
  call DeallocateArray(patch%flow_upwind_direction_bc)


  if (associated(patch%material_property_array)) &
    deallocate(patch%material_property_array)
  nullify(patch%material_property_array)
  ! Since this linked list will be destroyed by realization, just nullify here
  nullify(patch%material_properties)
  if (associated(patch%saturation_function_array)) &
    deallocate(patch%saturation_function_array)
  nullify(patch%saturation_function_array)
  ! Since this linked list will be destroyed by realization, just nullify here
  nullify(patch%saturation_functions)
  if (associated(patch%characteristic_curves_array)) &
    deallocate(patch%characteristic_curves_array)
  nullify(patch%characteristic_curves_array)
  ! Since this linked list will be destroyed by realization, just nullify here
  nullify(patch%characteristic_curves)

  if (associated(patch%char_curves_thermal_array)) &
       deallocate(patch%char_curves_thermal_array)
  nullify(patch%char_curves_thermal_array)
  nullify(patch%characteristic_curves_thermal)

  if (associated(patch%material_transform_array)) &
    deallocate(patch%material_transform_array)
  nullify(patch%material_transform_array)
  nullify(patch%material_transform)
  
  ! solely nullify grid since destroyed in discretization
  nullify(patch%grid)
  call RegionDestroyList(patch%region_list)
  call CouplerDestroyList(patch%boundary_condition_list)
  call CouplerDestroyList(patch%initial_condition_list)
  call CouplerDestroyList(patch%source_sink_list)

  call ObservationDestroyList(patch%observation_list)
  call IntegralFluxDestroyList(patch%integral_flux_list)
  call StrataDestroyList(patch%strata_list)

  call AuxDestroy(patch%aux)

  ! these are solely pointers, must not destroy.
  nullify(patch%reaction_base)
  nullify(patch%reaction)
  nullify(patch%reaction_nw)
  nullify(patch%datasets)
  nullify(patch%field)

  deallocate(patch)
  nullify(patch)

end subroutine PatchDestroy

end module Patch_module
