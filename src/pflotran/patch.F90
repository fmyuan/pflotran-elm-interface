module Patch_module

#include "petsc/finclude/petscvec.h"
  use petscvec

  use Option_module
  use Grid_module
  use Coupler_module
  use Observation_module
  use Integral_Flux_module
  use Strata_module
  use Region_module
  use Dataset_Base_class
  use Material_module
  use Field_module
  use Characteristic_Curves_module

  use Auxiliary_module

  use PFLOTRAN_Constants_module

  use MpFlow_Aux_module

  implicit none

  private

  type, public :: patch_type

    PetscInt :: id

    ! These arrays will be used by all modes, mode-specific arrays should
    ! go in the auxiliary data stucture for that mode
    PetscInt, pointer :: imat(:)
    PetscInt, pointer :: imat_internal_to_external(:)
    PetscInt, pointer :: sat_func_id(:)

    PetscReal, pointer :: internal_velocities(:,:)
    PetscReal, pointer :: boundary_velocities(:,:)
    PetscReal, pointer :: internal_flow_fluxes(:,:)
    PetscReal, pointer :: boundary_flow_fluxes(:,:)
    ! fluid fluxes in moles/sec
    PetscReal, pointer :: ss_flow_fluxes(:,:)
    ! volumetric flux (m^3/sec) for liquid phase needed for transport
    PetscReal, pointer :: ss_flow_vol_fluxes(:,:)

    ! for TH surface/subsurface
    PetscReal, pointer :: boundary_energy_flux(:,:)

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
    class(characteristic_curves_type), pointer :: characteristic_curves
    type(characteristic_curves_ptr_type), pointer :: characteristic_curves_array(:)

    type(strata_list_type), pointer :: strata_list
    type(observation_list_type), pointer :: observation_list
    type(integral_flux_list_type), pointer :: integral_flux_list

    ! Pointers to objects in mother realization object
    type(field_type), pointer :: field
    class(dataset_base_type), pointer :: datasets

    type(auxiliary_type) :: aux

    type(patch_type), pointer :: next

    PetscInt :: surf_or_subsurf_flag  ! Flag to identify if the current patch
                                      ! is a surface or subsurface (default)

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
            PatchLocalizeRegions, &
            PatchGetVariable, PatchGetVariableValueAtCell, &
            PatchSetVariable, &
            PatchCountCells, PatchGetIvarsFromKeyword, &
            PatchGetVarNameFromKeyword, &
            PatchCalculateCFL1Timestep, &
            PatchGetCellCenteredVelocities, &
            PatchGetCompMassInRegion, &
            PatchGetWaterMassInRegion, &
            PatchGetCompMassInRegionAssign, &
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
  patch%surf_or_subsurf_flag = SUBSURFACE
  nullify(patch%imat)
  nullify(patch%imat_internal_to_external)
  nullify(patch%sat_func_id)
  nullify(patch%internal_velocities)
  nullify(patch%boundary_velocities)
  nullify(patch%internal_flow_fluxes)
  nullify(patch%boundary_flow_fluxes)
  nullify(patch%ss_flow_fluxes)
  nullify(patch%ss_flow_vol_fluxes)
  nullify(patch%boundary_energy_flux)

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
  nullify(patch%characteristic_curves)
  nullify(patch%characteristic_curves_array)

  allocate(patch%observation_list)
  call ObservationInitList(patch%observation_list)
  allocate(patch%integral_flux_list)
  call IntegralFluxInitList(patch%integral_flux_list)
  allocate(patch%strata_list)
  call StrataInitList(patch%strata_list)

  call AuxInit(patch%aux)

  nullify(patch%field)
  nullify(patch%datasets)

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

subroutine PatchProcessCouplers(patch,flow_conditions, &
                                option)

  !
  ! Assigns conditions and regions to couplers
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  use Option_module
  use Material_module
  use Condition_module

  use Connection_module

  implicit none

  type(patch_type) :: patch
  type(condition_list_type) :: flow_conditions
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
        if (.not.associated(coupler%flow_condition%liq_rate)) then
          temp_int = 0
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
        ! gb: Depending on a surface/subsurface patch, use corresponding
        !     material properties
        if (patch%surf_or_subsurf_flag == SUBSURFACE) then
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
      case(OBSERVATION_SCALAR)
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
        if (observation%region%num_cells == 0) then
          ! remove the observation object
          call ObservationRemoveFromList(observation,patch%observation_list)
        endif
      case(OBSERVATION_FLUX)
        coupler => CouplerGetPtrFromList(observation%linkage_name, &
                                         patch%boundary_condition_list)
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
  enddo

  temp_int = ConnectionGetNumberInList(patch%grid%internal_connection_set_list)
  temp_int = max(temp_int,1)
  nphase = max(option%nfluids,1)

  ! all simulations
  allocate(patch%internal_velocities(nphase,temp_int))
  patch%internal_velocities = 0.d0

  ! flow
  if (option%nflowdof > 0) then
    if (option%flow%store_fluxes .or. &
        (patch%surf_or_subsurf_flag == SURFACE)) then
      allocate(patch%internal_flow_fluxes(option%nflowdof,temp_int))
      patch%internal_flow_fluxes = 0.d0
    endif
  endif


  temp_int = CouplerGetNumConnectionsInList(patch%boundary_condition_list)

  if (temp_int > 0) then
    ! all simulations
    allocate(patch%boundary_velocities(nphase,temp_int))
    patch%boundary_velocities = 0.d0
    ! flow
    if (option%nflowdof > 0) then
      if (option%flow%store_fluxes .or. &
          (patch%surf_or_subsurf_flag == SURFACE)) then
        allocate(patch%boundary_flow_fluxes(option%nflowdof,temp_int))
        patch%boundary_flow_fluxes = 0.d0
      endif
      ! surface/subsurface storage
      if (option%iflowmode == MPFLOW_MODE) then
        allocate(patch%boundary_energy_flux(2,temp_int))  ! why 2?
        patch%boundary_energy_flux = 0.d0
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

  use Global_Aux_module
  use Condition_module
  use MpFlow_Aux_module

  implicit none

  type(coupler_list_type), pointer :: coupler_list
  type(patch_type), pointer :: patch
  type(option_type) :: option

  PetscInt :: num_connections
  PetscBool :: force_update_flag

  type(coupler_type), pointer :: coupler

  PetscInt :: idof
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

          if (associated(coupler%flow_condition%liq_pressure) .or. &
              associated(coupler%flow_condition%liq_molarity) .or. &
              associated(coupler%flow_condition%liq_saturation) .or. &
              associated(coupler%flow_condition%temperature)) then

            ! allocate arrays that match the number of connections
            select case(option%iflowmode)

               case(MPFLOW_MODE)
                temp_int = 2
                select case(coupler%flow_condition%liq_pressure%itype)
                  case(HYDROSTATIC_CONDUCTANCE_BC, &
                       DIRICHLET_CONDUCTANCE_BC, &
                       HET_HYDROSTATIC_CONDUCTANCE_BC)
                    temp_int = temp_int + 1
                end select
                allocate(coupler%flow_aux_real_var(temp_int,num_connections))
                allocate(coupler%flow_aux_int_var(1,num_connections))
                coupler%flow_aux_real_var = 0.d0
                coupler%flow_aux_int_var = 0

              case default
                option%io_buffer = 'Failed allocation for flow condition "' // &
                  trim(coupler%flow_condition%name)
                call PrintErrMsg(option)
            end select

          else if (associated(coupler%flow_condition%liq_rate)) then
            option%io_buffer = 'Flow condition "' // &
              trim(coupler%flow_condition%name) // '" can only be used in a &
              &SOURCE_SINK since a rate is prescribed.'
            call PrintErrMsg(option)
          endif ! associated(coupler%flow_condition%liq_pressure)

        else if (coupler%itype == SRC_SINK_COUPLER_TYPE) then

          if (associated(coupler%flow_condition%liq_rate)) then

            select case(coupler%flow_condition%liq_rate%itype)
              case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS, &
                   VOLUMETRIC_RATE_SS,MASS_RATE_SS, &
                   HET_VOL_RATE_SS,HET_MASS_RATE_SS)
                select case(option%iflowmode)
                  case(MPFLOW_MODE)
                    allocate(coupler%flow_aux_real_var(option%nflowdof,num_connections))
                    coupler%flow_aux_real_var = 0.d0
                  case default
                    string = GetSubConditionName(coupler%flow_condition%liq_rate%itype)
                    option%io_buffer='Source/Sink of rate%itype = "' // &
                      trim(adjustl(string)) // '", not implemented in this mode.'
                    call PrintErrMsg(option)
                end select
              case default
                string = GetSubConditionName(coupler%flow_condition%liq_rate%itype)
                option%io_buffer = &
                  FlowConditionUnknownItype(coupler%flow_condition,'rate', &
                                            string)
                call PrintErrMsg(option)
            end select
          endif ! associated(coupler%flow_condition%liq_rate)
        endif ! coupler%itype == SRC_SINK_COUPLER_TYPE
      endif ! associated(coupler%flow_condition)
    endif ! associated(coupler%connection_set)

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
  use Saturation_module


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
          case(MPFLOW_MODE)
            call PatchUpdateCouplerAuxVarsTH(patch,coupler,option)
        end select
      endif
    endif

    coupler => coupler%next
  enddo

end subroutine PatchUpdateCouplerAuxVars

! ************************************************************************** !

! ************************************************************************** !

subroutine PatchUpdateCouplerAuxVarsTH(patch,coupler,option)
  !
  ! Updates flow auxiliary variables associated
  ! with a coupler for MPFLOW_MODE
  !
  ! Author: Glenn Hammond
  ! Date: 11/26/07
  !

  use Option_module
  use Condition_module
  use Hydrostatic_module
  use Saturation_module


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

  !type(flow_condition_type), pointer :: general
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
  !PetscInt :: iphase

  num_connections = coupler%connection_set%num_connections

  flow_condition => coupler%flow_condition

  if (associated(flow_condition%liq_pressure)) then
    !geh: this is a fix for an Intel compiler bug. Not sure why Intel cannot
    !     access flow_condition%iphase directly....
    !iphase = flow_condition%iphase
    !coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = iphase
!    coupler%flow_aux_int_var(COUPLER_IPHASE_INDEX,1:num_connections) = &
!                                                        flow_condition%iphase
    select case(flow_condition%liq_pressure%itype)
      case(DIRICHLET_BC,NEUMANN_BC,ZERO_GRADIENT_BC,SPILLOVER_BC, &
           DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC)
        select type(selector =>flow_condition%liq_pressure%dataset)
          class is(dataset_ascii_type)
            coupler%flow_aux_real_var(MPFLOW_PRESSURE_DOF,1:num_connections) = &
              selector%rarray(1)
          class is(dataset_gridded_hdf5_type)
            call PatchUpdateCouplerGridDataset(coupler,option, &
                                               patch%grid,selector, &
                                               MPFLOW_PRESSURE_DOF)
          class is(dataset_common_hdf5_type)
            ! skip cell indexed datasets used in initial conditions
          class default
            call PrintMsg(option,'th%pressure%itype,DIRICHLET_BC')
            call DatasetUnknownClass(selector,option, &
                                     'PatchUpdateCouplerAuxVarsTH')
        end select
        select case(flow_condition%liq_pressure%itype)
          case(DIRICHLET_CONDUCTANCE_BC)
            coupler%flow_aux_real_var(MPFLOW_CONDUCTANCE_DOF, &
                                      1:num_connections) = &
                                           flow_condition%liq_pressure%aux_real(1)
        end select
      case(HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC)
        call HydrostaticUpdateCoupler(coupler,option,patch%grid)
      case(HET_DIRICHLET_BC,HET_HYDROSTATIC_SEEPAGE_BC, &
           HET_HYDROSTATIC_CONDUCTANCE_BC)
        call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                flow_condition%liq_pressure%dataset,MPFLOW_PRESSURE_DOF,option)
        if (flow_condition%liq_pressure%itype == &
            HET_HYDROSTATIC_CONDUCTANCE_BC) then
          coupler%flow_aux_real_var(MPFLOW_CONDUCTANCE_DOF,1:num_connections) = &
            flow_condition%liq_pressure%aux_real(1)
        endif
      case(HET_SURF_HYDROSTATIC_SEEPAGE_BC)
        ! Do nothing, since this BC type is only used for coupling of
        ! surface-subsurface model
      case default
        string = &
          GetSubConditionName(flow_condition%liq_pressure%itype)
        option%io_buffer = &
          FlowConditionUnknownItype(flow_condition,'TH pressure',string)
        call PrintErrMsg(option)
    end select
    if (associated(flow_condition%temperature)) then
      select case(flow_condition%temperature%itype)
        case(DIRICHLET_BC,ZERO_GRADIENT_BC)
          select type(selector =>flow_condition%temperature%dataset)
            class is(dataset_ascii_type)
              if (flow_condition%liq_pressure%itype /= HYDROSTATIC_BC .or. &
                 (flow_condition%liq_pressure%itype == HYDROSTATIC_BC .and. &
                 flow_condition%temperature%itype /= DIRICHLET_BC)) then
                coupler%flow_aux_real_var(MPFLOW_TEMPERATURE_DOF, &
                                          1:num_connections) = &
                  selector%rarray(1)
              endif
            class is(dataset_gridded_hdf5_type)
              call PatchUpdateCouplerGridDataset(coupler,option, &
                                                 patch%grid,selector, &
                                                 MPFLOW_TEMPERATURE_DOF)
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
                  MPFLOW_TEMPERATURE_DOF,option)
        case default
          string = &
            GetSubConditionName(flow_condition%temperature%itype)
          option%io_buffer = &
            FlowConditionUnknownItype(flow_condition,'TH temperature',string)
          call PrintErrMsg(option)
      end select
    endif
    if (associated(flow_condition%energy_flux)) then
      coupler%flow_aux_real_var(MPFLOW_TEMPERATURE_DOF,1:num_connections) = &
        flow_condition%energy_flux%dataset%rarray(1)
    endif
  endif

  apply_temp_cond = PETSC_FALSE
  if (associated(flow_condition%temperature) .and. associated(flow_condition%liq_pressure)) then
    if (flow_condition%liq_pressure%itype /= HYDROSTATIC_BC) then
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
            coupler%flow_aux_real_var(MPFLOW_TEMPERATURE_DOF, &
                                      1:num_connections) = &
              selector%rarray(1)
          class is(dataset_gridded_hdf5_type)
            call PatchUpdateCouplerGridDataset(coupler,option, &
                                               patch%grid,selector, &
                                               MPFLOW_TEMPERATURE_DOF)
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
                MPFLOW_TEMPERATURE_DOF,option)
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
            coupler%flow_aux_real_var(MPFLOW_TEMPERATURE_DOF, &
                                      1:num_connections) = &
              selector%rarray(1)
          class is(dataset_gridded_hdf5_type)
            call PatchUpdateCouplerGridDataset(coupler,option, &
                                               patch%grid,selector, &
                                               MPFLOW_TEMPERATURE_DOF)
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
  if (associated(flow_condition%liq_rate)) then
    select case(flow_condition%liq_rate%itype)
      case (HET_MASS_RATE_SS,HET_VOL_RATE_SS)
        call PatchUpdateHetroCouplerAuxVars(patch,coupler, &
                                            flow_condition%liq_rate%dataset, &
                                            MPFLOW_PRESSURE_DOF,option)
      case(SCALED_MASS_RATE_SS,SCALED_VOLUMETRIC_RATE_SS)
        call PatchScaleSourceSink(patch,coupler,flow_condition%liq_rate%isubtype, &
                                  option)
        rate_scale_type = flow_condition%liq_rate%isubtype
      case(MASS_RATE_SS,VOLUMETRIC_RATE_SS)
      ! do nothing here
      case default
        string = &
          GetSubConditionName(flow_condition%liq_rate%itype)
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
        coupler%flow_aux_real_var(MPFLOW_TEMPERATURE_DOF,1:num_connections) = &
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
                MPFLOW_TEMPERATURE_DOF,option)
      case default
        string = &
          GetSubConditionName(flow_condition%energy_rate%itype)
        option%io_buffer = &
          FlowConditionUnknownItype(flow_condition,'TH energy rate',string)
        call PrintErrMsg(option)
    end select
  endif
  if (associated(flow_condition%liq_saturation)) then
    call SaturationUpdateCoupler(coupler,option,patch%grid, &
                                 patch%characteristic_curves_array, &
                                 patch%sat_func_id)
  endif

end subroutine PatchUpdateCouplerAuxVarsTH

! ************************************************************************** !

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
  use Material_Aux_class
  use Variables_module, only : PERMEABILITY_X

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
  class(material_auxvar_type), pointer :: material_auxvars(:)

  field => patch%field
  grid => patch%grid
  material_auxvars => patch%aux%Material%auxvars

  grid => patch%grid

  call VecZeroEntries(field%work,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)

  cur_connection_set => source_sink%connection_set

  select case(iscale_type)
    case(SCALE_BY_VOLUME)
      do iconn = 1, cur_connection_set%num_connections
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        vec_ptr(local_id) = vec_ptr(local_id) + &
          material_auxvars(ghosted_id)%volume
      enddo
    case(SCALE_BY_PERM)
      do iconn = 1, cur_connection_set%num_connections
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
        vec_ptr(local_id) = vec_ptr(local_id) + &
          ! this function protects from error in gfortran compiler when indexing
          ! the permeability array
          MaterialAuxVarGetValue(material_auxvars(ghosted_id), &
                                 PERMEABILITY_X) * &
          material_auxvars(ghosted_id)%volume
      enddo
    case(SCALE_BY_NEIGHBOR_PERM)
      do iconn = 1, cur_connection_set%num_connections
        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)
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
          sum = sum + &
                MaterialAuxVarGetValue(material_auxvars(neighbor_ghosted_id), &
                                       PERMEABILITY_X) * &
                grid%structured_grid%dx(neighbor_ghosted_id)* &
                grid%structured_grid%dz(neighbor_ghosted_id)
        enddo
        ! z-direction
        do while (icount < x_count + y_count + z_count)
          icount = icount + 1
          neighbor_ghosted_id = ghosted_neighbors(icount)
          sum = sum + &
                MaterialAuxVarGetValue(material_auxvars(neighbor_ghosted_id), &
                                       PERMEABILITY_X) * &
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
      case(MPFLOW_MODE)
        source_sink%flow_aux_real_var(ONE_INTEGER,iconn) = &
          vec_ptr(local_id)
    end select
  enddo
  call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)

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

  if (isub_condition>option%nflowdof*option%nfluids) then
    option%io_buffer='ERROR: PatchUpdateHetroCouplerAuxVars  '// &
      'isub_condition > option%nflowdof*(option%nfluids.'
    call PrintErrMsg(option)
  endif

  if (option%iflowmode/=MPFLOW_MODE) then
    option%io_buffer='PatchUpdateHetroCouplerAuxVars only implemented '// &
      ' for TH mode.'
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

! ************************************************************************** !

! ************************************************************************** !

subroutine PatchGetVariable1(patch,field,option, &
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

  use MpFlow_Aux_module
  use Output_Aux_module
  use Variables_module
  use Material_Aux_class

  implicit none

  type(option_type), pointer :: option
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
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscReal :: xmass, lnQKgas, ehfac, eh0, pe0, ph0, tk
  PetscReal :: tempreal
  PetscInt :: tempint, tempint2
  PetscInt :: irate, istate, irxn, ifo2, jcomp, comp_id
  PetscInt :: ivar_temp
  PetscErrorCode :: ierr


  grid => patch%grid
  material_auxvars => patch%aux%Material%auxvars


  call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)
  vec_ptr(:) = UNINITIALIZED_DOUBLE

  iphase = 1
  select case(ivar)
    case(TEMPERATURE,LIQUID_PRESSURE,GAS_PRESSURE,AIR_PRESSURE, &
         LIQUID_SATURATION,GAS_SATURATION,ICE_SATURATION, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY, &
         LIQUID_DENSITY,GAS_DENSITY,GAS_DENSITY_MOL,LIQUID_VISCOSITY, &
         GAS_VISCOSITY,CAPILLARY_PRESSURE,LIQUID_DENSITY_MOL, &
         LIQUID_MOBILITY,GAS_MOBILITY,SC_FUGA_COEFF,ICE_DENSITY, &
         LIQUID_HEAD,VAPOR_PRESSURE,SATURATION_PRESSURE, &
         MAXIMUM_PRESSURE,LIQUID_MASS_FRACTION,GAS_MASS_FRACTION, &
         OIL_PRESSURE,OIL_SATURATION,OIL_DENSITY,OIL_DENSITY_MOL,OIL_ENERGY, &
         OIL_MOBILITY,OIL_VISCOSITY,BUBBLE_POINT, &
         SOLVENT_PRESSURE,SOLVENT_SATURATION,SOLVENT_DENSITY, &
         SOLVENT_DENSITY_MOL,SOLVENT_ENERGY,SOLVENT_MOBILITY )

      if (associated(patch%aux%Flow)) then
        select case(ivar)
          case(CAPILLARY_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = patch%aux%Flow%auxvars(grid%nL2G(local_id))%pc
            enddo
          case(TEMPERATURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%auxvars(grid%nL2G(local_id))%temp
            enddo
          case(LIQUID_PRESSURE,MAXIMUM_PRESSURE)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%auxvars(grid%nL2G(local_id))%pres(LIQUID_PHASE)
            enddo
          case(LIQUID_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%auxvars(grid%nL2G(local_id))%sat(LIQUID_PHASE)
            enddo
          case(LIQUID_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Global%auxvars(grid%nL2G(local_id))%den_kg(LIQUID_PHASE)
            enddo
          case(GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY,GAS_VISCOSITY)
            call PatchUnsupportedVariable('MPFLOW','for gas phase',option)
          case(GAS_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Flow%auxvars(grid%nL2G(local_id))%sat(GAS_PHASE)
            enddo
          case(ICE_SATURATION)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Flow%auxvars(grid%nL2G(local_id))%sat(SOLID_PHASE)
            enddo
          case(ICE_DENSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                patch%aux%Flow%auxvars(grid%nL2G(local_id))%den(SOLID_PHASE)*FMWH2O
            enddo
          case(LIQUID_VISCOSITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                  patch%aux%Flow%auxvars(grid%nL2G(local_id))%vis
            enddo
          case(LIQUID_MOBILITY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                  patch%aux%Flow%auxvars(grid%nL2G(local_id))%kvr
            enddo
          case(LIQUID_MOLE_FRACTION)
            call PatchUnsupportedVariable('MPFLOW','LIQUID_MOLE_FRACTION',option)
          case(LIQUID_ENERGY)
            do local_id=1,grid%nlmax
              vec_ptr(local_id) = &
                  patch%aux%Flow%auxvars(grid%nL2G(local_id))%u(LIQUID_PHASE)
            enddo
          case default
            call PatchUnsupportedVariable('MPFLOW',ivar,option)
        end select

      endif

    case(POROSITY,BASE_POROSITY,INITIAL_POROSITY, &
         VOLUME,TORTUOSITY,SOIL_COMPRESSIBILITY, &
         SOIL_REFERENCE_PRESSURE)
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

    case(LIQUID_RELATIVE_PERMEABILITY)
      select case(option%iflowmode)
        case(MPFLOW_MODE)
          do local_id=1,grid%nlmax
            vec_ptr(local_id) = &
              patch%aux%Flow%auxvars(grid%nL2G(local_id))%kvr * &
              patch%aux%Flow%auxvars(grid%nL2G(local_id))%vis
          enddo
        case default
          option%io_buffer = 'Output of liquid phase relative permeability &
            &not supported for current flow mode.'
      end select
    case(GAS_RELATIVE_PERMEABILITY)
      select case(option%iflowmode)
        case(MPFLOW_MODE)
 !         do local_id=1,grid%nlmax  ! TODO
 !           vec_ptr(local_id) = &
 !             patch%aux%Flow%auxvars(grid%nL2G(local_id))% &
 !               kr(GAS_PHASE)
 !         enddo
        case default
          option%io_buffer = 'Output of gas phase relative permeability &
            &not supported for current flow mode.'
      end select
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
    case default
      call PatchUnsupportedVariable(ivar,option)
  end select

  call VecRestoreArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)

end subroutine PatchGetVariable1

! ************************************************************************** !

function PatchGetVariableValueAtCell(patch,field,option, &
                                     output_option,ghosted_id, &
                                     ivar,isubvar,isubvar2)
  !
  ! Returns variables indexed by ivar,
  ! isubvar, local id
  !
  ! Author: Glenn Hammond
  ! Date: 02/11/08
  !

  use Grid_module
  use Option_module
  use Field_module


  use MpFlow_Aux_module
  use Output_Aux_module
  use Variables_module
  use Material_Aux_class

  implicit none

  PetscReal :: PatchGetVariableValueAtCell
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  class(material_auxvar_type), pointer :: material_auxvars(:)
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
         LIQUID_SATURATION,GAS_SATURATION,ICE_SATURATION, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY, &
         LIQUID_DENSITY,GAS_DENSITY,GAS_DENSITY_MOL,LIQUID_VISCOSITY, &
         GAS_VISCOSITY,AIR_PRESSURE,CAPILLARY_PRESSURE, &
         LIQUID_MOBILITY,GAS_MOBILITY,SC_FUGA_COEFF,ICE_DENSITY, &
         SECONDARY_TEMPERATURE,LIQUID_DENSITY_MOL, &
         LIQUID_HEAD,VAPOR_PRESSURE,SATURATION_PRESSURE,MAXIMUM_PRESSURE, &
         LIQUID_MASS_FRACTION,GAS_MASS_FRACTION, &
         OIL_PRESSURE,OIL_SATURATION,OIL_DENSITY,OIL_DENSITY_MOL,OIL_ENERGY, &
         OIL_MOBILITY,OIL_VISCOSITY,BUBBLE_POINT, &
         SOLVENT_PRESSURE,SOLVENT_SATURATION,SOLVENT_DENSITY, &
         SOLVENT_DENSITY_MOL,SOLVENT_ENERGY,SOLVENT_MOBILITY)

      if (associated(patch%aux%Flow)) then
        select case(ivar)
          case(TEMPERATURE)
            value = patch%aux%Global%auxvars(ghosted_id)%temp
          case(LIQUID_PRESSURE,MAXIMUM_PRESSURE)
            value = patch%aux%Global%auxvars(ghosted_id)%pres(LIQUID_PHASE)
          case(CAPILLARY_PRESSURE)
            value = patch%aux%Flow%auxvars(ghosted_id)%pc
          case(LIQUID_SATURATION)
            value = patch%aux%Global%auxvars(ghosted_id)%sat(LIQUID_PHASE)
          case(LIQUID_DENSITY)
            value = patch%aux%Global%auxvars(ghosted_id)%den_kg(LIQUID_PHASE)
          case(LIQUID_VISCOSITY)
            value = patch%aux%Flow%auxvars(ghosted_id)%vis
          case(LIQUID_MOBILITY)
            value = patch%aux%Flow%auxvars(ghosted_id)%kvr
          case(LIQUID_MOLE_FRACTION)
            call PatchUnsupportedVariable('MPFLOW','LIQUID_MOLE_FRACTION',option)
          case(LIQUID_ENERGY)
            value = patch%aux%Flow%auxvars(ghosted_id)%u(LIQUID_PHASE)
          case(GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY)
            call PatchUnsupportedVariable('MPFLOW','GAS_MOLE_FRACTION',option)
          case(GAS_SATURATION)
            value = patch%aux%Flow%auxvars(ghosted_id)%sat(GAS_PHASE)
          case(ICE_SATURATION)
            value = patch%aux%Flow%auxvars(ghosted_id)%sat(SOLID_PHASE)
          case(ICE_DENSITY)
            value = patch%aux%Flow%auxvars(ghosted_id)%den(SOLID_PHASE)*FMWH2O
          case default
            call PatchUnsupportedVariable('MPFLOW',ivar,option)
        end select

      endif

    !case(STATE,PHASE)
    !  value = patch%aux%Global%auxvars(ghosted_id)%istate

    case(POROSITY,BASE_POROSITY,INITIAL_POROSITY, &
         VOLUME,TORTUOSITY,SOIL_COMPRESSIBILITY,SOIL_REFERENCE_PRESSURE)
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

    case(LIQUID_RELATIVE_PERMEABILITY)
      select case(option%iflowmode)
        case(MPFLOW_MODE)
          value = patch%aux%Flow%auxvars(ghosted_id)%kvr / &
                  patch%aux%Flow%auxvars(ghosted_id)%vis
        case default
          option%io_buffer = 'Output of liquid phase relative permeability &
            &not supported for current flow mode.'
      end select
    case(GAS_RELATIVE_PERMEABILITY)
      select case(option%iflowmode)
        case(MPFLOW_MODE)
!          value = patch%aux%Flow%auxvars(ghosted_id)% &
!                    kr(GAS_PHASE)
        case default
          option%io_buffer = 'Output of gas phase relative permeability &
            &not supported for current flow mode.'
      end select
    case(MATERIAL_ID)
      value = patch%imat_internal_to_external(abs(patch%imat(ghosted_id)))
    case(PROCESS_ID)
      value = grid%nG2A(ghosted_id)
    case(NATURAL_ID)
      value = option%myrank
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
  use Material_Aux_class

  implicit none

  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  Vec :: vec
  PetscInt :: vec_format
  PetscInt :: ivar
  PetscInt :: isubvar

  PetscInt :: local_id, ghosted_id
  type(grid_type), pointer :: grid
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscErrorCode :: ierr

  grid => patch%grid
  material_auxvars => patch%aux%Material%auxvars

  call VecGetArrayF90(vec,vec_ptr,ierr);CHKERRQ(ierr)

  if (vec_format == NATURAL) then
    call PrintErrMsg(option,&
                     'NATURAL vector format not supported by PatchSetVariable')
  endif

  select case(ivar)
    case(TEMPERATURE,LIQUID_PRESSURE,GAS_PRESSURE,LIQUID_SATURATION, &
         GAS_SATURATION,AIR_PRESSURE,CAPILLARY_PRESSURE, &
         LIQUID_MOLE_FRACTION,GAS_MOLE_FRACTION,LIQUID_ENERGY,GAS_ENERGY, &
         LIQUID_DENSITY,GAS_DENSITY,GAS_DENSITY_MOL,LIQUID_VISCOSITY, &
         GAS_VISCOSITY,ICE_SATURATION, &
         LIQUID_MOBILITY,GAS_MOBILITY)

      if (associated(patch%aux%Flow)) then
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
                patch%aux%Global%auxvars(grid%nL2G(local_id))%pres(LIQUID_PHASE) = &
                  vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%auxvars(ghosted_id)%pres(LIQUID_PHASE) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_SATURATION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%auxvars(grid%nL2G(local_id))%sat(LIQUID_PHASE) = &
                  vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%auxvars(ghosted_id)%sat(LIQUID_PHASE) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_DENSITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Global%auxvars(grid%nL2G(local_id))%den(LIQUID_PHASE) = &
                  vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Global%auxvars(ghosted_id)%den_kg(LIQUID_PHASE) = &
                  vec_ptr(ghosted_id)
              enddo
            endif
          case(GAS_MOLE_FRACTION,GAS_ENERGY,GAS_DENSITY)
            call PrintErrMsg(option,'GAS_MOLE_FRACTION not supported by TH')
          case(GAS_SATURATION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Flow%auxvars(grid%nL2G(local_id))%sat(GAS_PHASE) = &
                  vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Flow%auxvars(ghosted_id)%sat(GAS_PHASE) = vec_ptr(ghosted_id)
              enddo
            endif
          case(ICE_SATURATION)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Flow%auxvars(grid%nL2G(local_id))%sat(SOLID_PHASE) = &
                  vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Flow%auxvars(ghosted_id)%sat(SOLID_PHASE) = vec_ptr(ghosted_id)
              enddo
            endif
          case(ICE_DENSITY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Flow%auxvars(grid%nL2G(local_id))%den(SOLID_PHASE) = &
                  vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Flow%auxvars(ghosted_id)%den(SOLID_PHASE) = vec_ptr(ghosted_id)
              enddo
            endif
          case(LIQUID_VISCOSITY)
          case(GAS_VISCOSITY)
          case(LIQUID_MOLE_FRACTION)
            call PrintErrMsg(option,'LIQUID_MOLE_FRACTION not supported by MPFLOW mode')
          case(LIQUID_ENERGY)
            if (vec_format == GLOBAL) then
              do local_id=1,grid%nlmax
                patch%aux%Flow%auxvars(grid%nL2G(local_id))%u(LIQUID_PHASE) = vec_ptr(local_id)
              enddo
            else if (vec_format == LOCAL) then
              do ghosted_id=1,grid%ngmax
                patch%aux%Flow%auxvars(ghosted_id)%u(LIQUID_PHASE) = vec_ptr(ghosted_id)
              enddo
            endif
        end select
      endif
      
    case(PRIMARY_MOLARITY,TOTAL_MOLALITY, &
         COLLOID_MOBILE,COLLOID_IMMOBILE)
      select case(ivar)
        case(PRIMARY_MOLARITY)
          call PrintErrMsg(option,'Setting of primary molarity at grid cell not supported.')
        case(TOTAL_MOLALITY)
          call PrintErrMsg(option,'Setting of total molality at grid cell not supported.')
        case(COLLOID_MOBILE)
          call PrintErrMsg(option,'Setting of mobile colloid concentration at grid cell not supported.')
        case(COLLOID_IMMOBILE)
          call PrintErrMsg(option,'Setting of immobile colloid concentration at grid cell not supported.')
      end select
    case(POROSITY,BASE_POROSITY,INITIAL_POROSITY)
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

subroutine PatchCalculateCFL1Timestep(patch,option,max_dt_cfl_1)
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
  use Material_Aux_class

  implicit none

  type(patch_type) :: patch
  type(option_type) :: option
  PetscReal :: max_dt_cfl_1

  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  PetscInt :: iconn
  PetscInt :: sum_connection
  PetscReal :: distance, fraction_upwind
  PetscReal :: por_sat_ave, por_sat_min, v_darcy, v_pore_ave, v_pore_max
  PetscInt :: local_id_up, local_id_dn
  PetscInt :: ghosted_id_up, ghosted_id_dn
  PetscInt :: ifluid

  PetscReal :: dt_cfl_1
  PetscErrorCode :: ierr

  field => patch%field
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  grid => patch%grid

  max_dt_cfl_1 = 1.d20

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
      do ifluid = 1, option%nfluids
        ! if the phase is not present in either cell, skip the connection
        if (.not.(global_auxvars(ghosted_id_up)%sat(ifluid) > 0.d0 .and. &
                  global_auxvars(ghosted_id_dn)%sat(ifluid) > 0.d0)) cycle
        por_sat_min = min(material_auxvars(ghosted_id_up)%porosity* &
                          global_auxvars(ghosted_id_up)%sat(ifluid), &
                          material_auxvars(ghosted_id_dn)%porosity* &
                          global_auxvars(ghosted_id_dn)%sat(ifluid))
        por_sat_ave = (fraction_upwind* &
                       material_auxvars(ghosted_id_up)%porosity* &
                       global_auxvars(ghosted_id_up)%sat(ifluid) + &
                      (1.d0-fraction_upwind)* &
                      material_auxvars(ghosted_id_dn)%porosity* &
                      global_auxvars(ghosted_id_dn)%sat(ifluid))
        v_darcy = patch%internal_velocities(ifluid,sum_connection)
        v_pore_max = v_darcy / por_sat_min
        v_pore_ave = v_darcy / por_sat_ave
        !geh: I use v_por_max to ensure that we limit the cfl based on the
        !     highest velocity through the face.  If porosity*saturation
        !     varies, the pore water velocity will be highest on the side
        !     of the face with the smalled value of porosity*saturation.
        dt_cfl_1 = distance / dabs(v_pore_max)
        max_dt_cfl_1 = min(dt_cfl_1,max_dt_cfl_1)
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
      do ifluid = 1, option%nfluids
        por_sat_ave = material_auxvars(ghosted_id_dn)%porosity* &
                      global_auxvars(ghosted_id_dn)%sat(ifluid)
        v_darcy = patch%boundary_velocities(ifluid,sum_connection)
        v_pore_ave = v_darcy / por_sat_ave
        dt_cfl_1 = distance / dabs(v_pore_ave)
        max_dt_cfl_1 = min(dt_cfl_1,max_dt_cfl_1)
      enddo
    enddo
    boundary_condition => boundary_condition%next
  enddo

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
  type(grid_type), pointer :: grid
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  type(coupler_type), pointer :: boundary_condition

  character(len=MAXWORDLENGTH) :: word
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
  PetscInt :: local_id_up
  PetscInt :: ghosted_id
  PetscInt :: ghosted_id_up
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
  PetscErrorCode :: ierr

  nullify(yet_to_be_found)
  nullify(connections)

  grid => patch%grid
  polygon => integral_flux%polygon
  plane => integral_flux%plane
  coordinates_and_directions => integral_flux%coordinates_and_directions
  vertices => integral_flux%vertices
  num_to_be_found = 0

  if (associated(polygon)) then
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
    num_to_be_found = size(coordinates_and_directions,2)
    allocate(yet_to_be_found(num_to_be_found))
    yet_to_be_found = PETSC_TRUE
  endif

  if (associated(vertices)) then
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
        ! sets up first boundayr condition in list
        if (.not.associated(boundary_condition)) then
          boundary_condition => patch%boundary_condition_list%first
          sum_connection = 0
          icount = 0
        endif
        cur_connection_set => boundary_condition%connection_set
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
            local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
            ! if one of the cells is ghosted, the process stores the flux only
            ! when the upwind cell is non-ghosted.
            if (local_id_up <= 0) cycle
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
          face_id = grid%unstructured_grid%connection_to_face(iconn)
          face_vertices_natural = &
            grid%unstructured_grid%face_to_vertex_natural(:,face_id)
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
            do ivert1 = 1, num_vertices1
              do ivert2 = 1, num_vertices2
                if (face_vertices_natural(ivert1) == vertices(ivert2,i)) then
                  found2 = PETSC_TRUE
                  exit
                endif
              enddo
            enddo
            if (found2) then
              reverse_direction = PETSC_FALSE
              ! search forward direction
              iv1 = ivert1
              iv2 = ivert2
              do ii = 1, num_vertices1
                if (face_vertices_natural(iv1) /= vertices(iv2,i)) then
                  found2 = PETSC_FALSE
                endif
                iv1 = iv1 + 1
                if (iv1 > num_vertices1) iv1 = 1
                iv2 = iv2 + 1
                if (iv2 > num_vertices1) iv2 = 1
              enddo
              ! search backward direction
              if (.not.found) then
                reverse_direction = PETSC_TRUE
                found2 = PETSC_TRUE
                iv1 = ivert1
                iv2 = ivert2
                do ii = 1, num_vertices1
                  if (face_vertices_natural(iv1) /= vertices(iv2,i)) then
                    found2 = PETSC_FALSE
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
        boundary_condition => boundary_condition%next
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
      '".  Please ensure that the polygon coincides with an internal &
      &cell boundary or a boundary condition.'
    call PrintErrMsg(option)
  else if (num_to_be_found > 0 .and. icount /= num_to_be_found) then
    write(word,*) num_to_be_found - icount
    option%io_buffer = trim(adjustl(word)) // &
      ' face(s) missed for INTEGRAL_FLUX "' // &
      trim(adjustl(integral_flux%name)) // &
      '".  Please ensure that the polygon coincides with an internal &
      &cell boundary or a boundary condition.'
    call PrintErrMsg(option)
  else
    write(option%io_buffer,*) icount
    option%io_buffer = trim(adjustl(option%io_buffer)) // ' connections found &
      &for integral flux "' // trim(adjustl(integral_flux%name)) // '".'
    call PrintMsg(option)
  endif

end subroutine PatchGetIntegralFluxConnections
! **************************************************************************** !

! **************************************************************************** !

subroutine PatchGetCompMassInRegion(cell_ids,num_cells,patch,option, &
                                    global_total_mass)
  !
  ! Calculates the total mass (aqueous, sorbed, and precipitated) in a region
  ! in units of mol.
  !
  ! Author: Jenn Frederick
  ! Date: 04/25/2016
  !
  use Global_Aux_module
  use Material_Aux_class

  use Grid_module
  use Option_module

  implicit none

  PetscInt, pointer :: cell_ids(:)
  PetscInt :: num_cells
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  PetscReal :: global_total_mass  ! [mol]

  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
  PetscReal :: aq_species_mass    ! [mol]
  PetscReal :: sorb_species_mass  ! [mol]
  PetscReal :: ppt_species_mass   ! [mol]
  PetscReal :: m3_water           ! [m^3-water]
  PetscReal :: m3_bulk            ! [m^3-bulk]
  PetscInt :: k, j, m
  PetscInt :: local_id, ghosted_id
  PetscErrorCode :: ierr
  PetscReal :: local_total_mass

  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  local_total_mass = 0.d0
  global_total_mass = 0.d0

  ! Loop through all cells in the region:
  do k = 1,num_cells
    local_id = cell_ids(k)
    ghosted_id = patch%grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    m3_water = material_auxvars(ghosted_id)%porosity * &         ! [-]
               global_auxvars(ghosted_id)%sat(LIQUID_PHASE) * &  ! [water]
               material_auxvars(ghosted_id)%volume               ! [m^3-bulk]
    m3_bulk = material_auxvars(ghosted_id)%volume                ! [m^3-bulk]
  enddo ! Cell loop

  ! Sum the local_total_mass across all processes that own the region:
  call MPI_Allreduce(local_total_mass,global_total_mass,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_SUM,option%mycomm,ierr)

end subroutine PatchGetCompMassInRegion

! **************************************************************************** !

subroutine PatchGetWaterMassInRegion(cell_ids,num_cells,patch,option, &
                                     global_water_mass)
  !
  ! Calculates the water mass in a region in kg
  !
  ! Author: Satish Karra
  ! Date: 09/20/2016
  !
  use Global_Aux_module
  use Material_Aux_class
  use Grid_module
  use Option_module

  implicit none

  PetscInt, pointer :: cell_ids(:)
  PetscInt :: num_cells
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  PetscReal :: global_water_mass

  type(global_auxvar_type), pointer :: global_auxvars(:)
  class(material_auxvar_type), pointer :: material_auxvars(:)
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
               global_auxvars(ghosted_id)%sat(LIQUID_PHASE) * &  ! [water]
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
  allocate(upwind_direction(option%nfluids,sum_connection))
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
  allocate(upwind_direction_bc(option%nfluids,sum_connection))
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
  call DeallocateArray(patch%sat_func_id)
  call DeallocateArray(patch%internal_velocities)
  call DeallocateArray(patch%boundary_velocities)
  call DeallocateArray(patch%internal_flow_fluxes)
  call DeallocateArray(patch%boundary_flow_fluxes)
  call DeallocateArray(patch%ss_flow_fluxes)
  call DeallocateArray(patch%ss_flow_vol_fluxes)

  call DeallocateArray(patch%boundary_energy_flux)

  call DeallocateArray(patch%flow_upwind_direction)
  call DeallocateArray(patch%flow_upwind_direction_bc)


  if (associated(patch%material_property_array)) &
    deallocate(patch%material_property_array)
  nullify(patch%material_property_array)
  ! Since this linked list will be destroyed by realization, just nullify here
  nullify(patch%material_properties)

  if (associated(patch%characteristic_curves_array)) &
    deallocate(patch%characteristic_curves_array)
  nullify(patch%characteristic_curves_array)
  ! Since this linked list will be destroyed by realization, just nullify here
  nullify(patch%characteristic_curves)


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

  nullify(patch%datasets)
  nullify(patch%field)

  deallocate(patch)
  nullify(patch)

end subroutine PatchDestroy

end module Patch_module
