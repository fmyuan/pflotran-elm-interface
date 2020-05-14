module Realization_Subsurface_class
  
#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Realization_Base_class
  use Option_module
  use Input_Aux_module
  use Region_module
  use Condition_module
  use Well_Data_class
  use Transport_Constraint_Base_module
  use Transport_Constraint_module
  use Material_module
  use Saturation_Function_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Thermal_module
  use Dataset_Base_class
  use Fluid_module
  use Patch_module
  use Reaction_Aux_module
  use NW_Transport_Aux_module
  

  implicit none

private

  type, public, extends(realization_base_type) :: realization_subsurface_type

    type(region_list_type), pointer :: region_list
    type(condition_list_type), pointer :: flow_conditions
    type(well_data_list_type), pointer :: well_data=>null()
    type(tran_condition_list_type), pointer :: transport_conditions
    type(tran_constraint_list_type), pointer :: transport_constraints
    
    class(tran_constraint_base_type), pointer :: sec_transport_constraint
    type(material_property_type), pointer :: material_properties
    type(fluid_property_type), pointer :: fluid_properties
    type(fluid_property_type), pointer :: fluid_property_array(:)
    type(saturation_function_type), pointer :: saturation_functions
    class(characteristic_curves_type), pointer :: characteristic_curves
    class(cc_thermal_type), pointer :: thermal_characteristic_curves
    class(dataset_base_type), pointer :: datasets
    
    class(dataset_base_type), pointer :: uniform_velocity_dataset
    character(len=MAXSTRINGLENGTH) :: nonuniform_velocity_filename

    class(reaction_rt_type), pointer :: reaction
    class(reaction_nw_type), pointer :: reaction_nw
    
  end type realization_subsurface_type

  interface RealizationCreate
    module procedure RealizationCreate1
    module procedure RealizationCreate2
  end interface
  
  public :: RealizationCreate, &
            RealizationStrip, &
            RealizationProcessCouplers, &
            RealizationInitAllCouplerAuxVars, &
            RealizationProcessConditions, &
            RealizationProcessDatasets, &
            RealizationAddWaypointsToList, &
            RealizationCreateDiscretization, &
            RealizationLocalizeRegions, &
            RealizationAddCoupler, &
            RealizationAddStrata, &
            RealizUpdateUniformVelocity, &
            RealizationRevertFlowParameters, &
            RealizStoreRestartFlowParams, &
!            RealizationGetVariable, &
!            RealizGetVariableValueAtCell, &
!            RealizationSetVariable, &
            RealizationPrintCouplers, &
            RealizationInitConstraints, &
            RealProcessMatPropAndSatFunc, &
            RealProcessFluidProperties, &
            RealizationUpdatePropertiesTS, &
            RealizationUpdatePropertiesNI, &
            RealizationCalcMineralPorosity, &
            RealizationCountCells, &
            RealizationPrintGridStatistics, &
            RealizationPassPtrsToPatches, &
            RealLocalToLocalWithArray, &
            RealizationCalculateCFL1Timestep, &
            RealizUpdateAllCouplerAuxVars, &
            RealizUnInitializedVarsFlow, &
            RealizUnInitializedVarsTran, &
            RealizationLimitDTByCFL

  !TODO(intel)
  ! public from Realization_Base_class
  !public :: RealizationGetVariable

contains

! ************************************************************************** !

function RealizationCreate1()
  ! 
  ! Allocates and initializes a new Realization object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  implicit none
  
  class(realization_subsurface_type), pointer :: RealizationCreate1
  
  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  
  nullify(option)
  RealizationCreate1 => RealizationCreate2(option)
  
end function RealizationCreate1  

! ************************************************************************** !

function RealizationCreate2(option)
  ! 
  ! Allocates and initializes a new Realization object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 
  
  implicit none
  
  type(option_type), pointer :: option
  
  class(realization_subsurface_type), pointer :: RealizationCreate2
  
  class(realization_subsurface_type), pointer :: realization
  
  allocate(realization)
  call RealizationBaseInit(realization,option)

  allocate(realization%region_list)
  call RegionInitList(realization%region_list)

  allocate(realization%flow_conditions)
  call FlowConditionInitList(realization%flow_conditions)
! Allocate well_data and create its list of wells
  allocate(realization%well_data)
  call WellDataInitList(realization%well_data,option%nphase)
  allocate(realization%transport_conditions)
  call TranConditionInitList(realization%transport_conditions)
  allocate(realization%transport_constraints)
  call TranConstraintInitList(realization%transport_constraints)

  nullify(realization%material_properties)
  nullify(realization%fluid_properties)
  nullify(realization%fluid_property_array)
  nullify(realization%saturation_functions)
  nullify(realization%characteristic_curves)
  if (option%iflowmode == G_MODE) then
    nullify(realization%thermal_characteristic_curves)
  end if
  nullify(realization%datasets)
  nullify(realization%uniform_velocity_dataset)
  nullify(realization%sec_transport_constraint)
  nullify(realization%reaction)
  nullify(realization%reaction_nw)
  realization%nonuniform_velocity_filename = ''

  RealizationCreate2 => realization
  
end function RealizationCreate2 

! ************************************************************************** !

subroutine RealizationCreateDiscretization(realization)
  ! 
  ! Creates grid
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Field_module
  use Grid_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_module, only : UGridEnsureRightHandRule
  use Grid_Structured_module, only : StructGridCreateTVDGhosts
  use Coupler_module
  use Discretization_module
  use Grid_Unstructured_Cell_module
  use DM_Kludge_module
  use Variables_module, only : VOLUME
  use Communicator_Structured_class, only : StructuredCommunicatorCreate
  use Communicator_Unstructured_class, only : UnstructuredCommunicatorCreate
  
  implicit none

  class(realization_subsurface_type) :: realization
  
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(coupler_type), pointer :: boundary_condition
  PetscErrorCode :: ierr
  PetscInt :: ivar

  option => realization%option
  field => realization%field
 
  discretization => realization%discretization
  
  call DiscretizationCreateDMs(discretization, option%nflowdof, &
                               option%ntrandof, option%nphase, &
                               option%ngeomechdof, option%n_stress_strain_dof, &
                               option)

  ! 1 degree of freedom, global
  call DiscretizationCreateVector(discretization,ONEDOF,field%work, &
                                  GLOBAL,option)
  call DiscretizationDuplicateVector(discretization,field%work, &
                                     field%porosity0)
  call DiscretizationDuplicateVector(discretization,field%work, &
                                     field%tortuosity0)
  call DiscretizationDuplicateVector(discretization,field%work, &
                                     field%volume0)
!geh: this is now allocated in 
!     init_subsurface.F90:SubsurfAllocMatPropDataStructs()
!  call DiscretizationDuplicateVector(discretization,field%work, &
!                                     field%compressibility0)
  if (option%flow%transient_porosity) then
    call DiscretizationDuplicateVector(discretization,field%work, &
                                       field%porosity_base_store)
    call DiscretizationDuplicateVector(discretization,field%work, &
                                       field%porosity_t)
    call DiscretizationDuplicateVector(discretization,field%work, &
                                       field%porosity_tpdt)
  endif

  if (option%geomech_on) then
    call DiscretizationDuplicateVector(discretization,field%work, &
                                       field%porosity_base_store)
    call DiscretizationDuplicateVector(discretization,field%work, &
                                       field%porosity_geomech_store)
  endif

  ! 1 degree of freedom, local
  call DiscretizationCreateVector(discretization,ONEDOF,field%work_loc, &
                                  LOCAL,option)
  
  if (option%nflowdof > 0) then

    ! 1-dof global  
    call DiscretizationDuplicateVector(discretization,field%work, &
                                       field%perm0_xx)
    call DiscretizationDuplicateVector(discretization,field%work, &
                                       field%perm0_yy)
    call DiscretizationDuplicateVector(discretization,field%work, &
                                       field%perm0_zz)

    ! 1-dof local
    call DiscretizationDuplicateVector(discretization,field%work_loc, &
                                       field%ithrm_loc)
    call DiscretizationDuplicateVector(discretization,field%work_loc, &
                                       field%icap_loc)
    
    ! ndof degrees of freedom, global
    call DiscretizationCreateVector(discretization,NFLOWDOF,field%flow_xx, &
                                    GLOBAL,option)
    call DiscretizationDuplicateVector(discretization,field%flow_xx, &
                                       field%flow_yy)
    call DiscretizationDuplicateVector(discretization,field%flow_xx, &
                                       field%flow_dxx)
    call DiscretizationDuplicateVector(discretization,field%flow_xx, &
                                       field%flow_r)
    call DiscretizationDuplicateVector(discretization,field%flow_xx, &
                                       field%flow_accum)
    call DiscretizationDuplicateVector(discretization,field%flow_xx, &
                                       field%flow_accum2)

    ! ndof degrees of freedom, local
    call DiscretizationCreateVector(discretization,NFLOWDOF,field%flow_xx_loc, &
                                    LOCAL,option)

    if ((option%iflowmode == RICHARDS_TS_MODE) .or. &
        (option%iflowmode == TH_TS_MODE)) then
      call DiscretizationCreateVector(discretization,NFLOWDOF,field%flow_xxdot, &
                                      GLOBAL,option)

      call DiscretizationCreateVector(discretization,NFLOWDOF,field%flow_xxdot_loc, &
                                      LOCAL,option)

    endif

  endif

  if (option%ntrandof > 0) then
    if (option%transport%reactive_transport_coupling == GLOBAL_IMPLICIT) then
      ! ndof degrees of freedom, global
      call DiscretizationCreateVector(discretization,NTRANDOF,field%tran_xx, &
                                      GLOBAL,option)
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_yy)
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_dxx)
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_r)
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_accum)

      ! ndof degrees of freedom, local
      call DiscretizationCreateVector(discretization,NTRANDOF,field%tran_xx_loc, &
                                      LOCAL,option)
      
      if (associated(realization%reaction_base)) then
        if (realization%reaction_base%use_log_formulation) then
          call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                             field%tran_log_xx)
          call DiscretizationDuplicateVector(discretization,field%tran_xx_loc, &
                                             field%tran_work_loc)
        endif
      endif
 
    else ! operator splitting
      ! ndof degrees of freedom, global
      ! create the 1 dof vector for solving the individual linear systems
      call DiscretizationCreateVector(discretization,ONEDOF,field%tran_rhs_coef, &
                                      GLOBAL,option)
      ! create the ntran dof vector for storage of the solution
      call DiscretizationCreateVector(discretization,NTRANDOF,field%tran_xx, &
                                      GLOBAL,option)
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_yy)
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_dxx)
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_rhs)

      ! ndof degrees of freedom, local
      ! again, just for storage of the current colution
      call DiscretizationCreateVector(discretization,NTRANDOF,field%tran_xx_loc, &
                                      LOCAL,option)

    endif
    
  endif

  grid => discretization%grid
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      ! set up nG2L, nL2G, etc.
      call GridMapIndices(grid, &
                          discretization%dm_1dof, &
                          discretization%stencil_type,&
                          option)
      if (option%itranmode == EXPLICIT_ADVECTION) then
        call StructGridCreateTVDGhosts(grid%structured_grid, &
                                       realization%reaction%naqcomp, &
                                       field%tran_xx, &
                                       discretization%dm_1dof%dm, &
                                       field%tvd_ghosts, &
                                       discretization%tvd_ghost_scatter, &
                                       option)
      endif
      call GridComputeSpacing(grid,discretization%origin_global,option)
      call GridComputeCoordinates(grid,discretization%origin_global,option)
      call GridComputeVolumes(grid,field%volume0,option)
      ! set up internal connectivity, distance, etc.
      call GridComputeInternalConnect(grid,option)
    case(UNSTRUCTURED_GRID)
      ! set up nG2L, NL2G, etc.
      call GridMapIndices(grid, &
                          discretization%dm_1dof, &
                          discretization%stencil_type,&
                          option)
      call GridComputeCoordinates(grid,discretization%origin_global,option, & 
                                    discretization%dm_1dof%ugdm) 
      if (grid%itype == IMPLICIT_UNSTRUCTURED_GRID) then
        call UGridEnsureRightHandRule(grid%unstructured_grid,grid%x, &
                                      grid%y,grid%z,grid%nG2A,grid%nL2G, &
                                      option)
      endif
      ! set up internal connectivity, distance, etc.
      call GridComputeInternalConnect(grid,option, &
                                      discretization%dm_1dof%ugdm) 
      call GridComputeVolumes(grid,field%volume0,option)
  end select
  call GridPrintExtents(grid,option)
 
  ! initialize to UNINITIALIZED_DOUBLE for check later that verifies all values 
  ! have been set
  call VecSet(field%porosity0,UNINITIALIZED_DOUBLE,ierr);CHKERRQ(ierr)

  ! Allocate vectors to hold temporally average output quantites
  if (realization%output_option%aveg_output_variable_list%nvars>0) then

    field%nvars = realization%output_option%aveg_output_variable_list%nvars
    allocate(field%avg_vars_vec(field%nvars))

    do ivar=1,field%nvars
      call DiscretizationDuplicateVector(discretization,field%work, &
                                         field%avg_vars_vec(ivar))
      call VecSet(field%avg_vars_vec(ivar),0.d0,ierr);CHKERRQ(ierr)
    enddo
  endif
       
  ! Allocate vectors to hold flowrate quantities
  if (realization%output_option%print_hdf5_mass_flowrate.or. &
     realization%output_option%print_hdf5_energy_flowrate.or. &
     realization%output_option%print_hdf5_aveg_mass_flowrate.or. &
     realization%output_option%print_hdf5_aveg_energy_flowrate) then
    call VecCreateMPI(option%mycomm, &
        (option%nflowdof*MAX_FACE_PER_CELL+1)*realization%patch%grid%nlmax, &
        PETSC_DETERMINE,field%flowrate_inst,ierr);CHKERRQ(ierr)
    call VecSet(field%flowrate_inst,0.d0,ierr);CHKERRQ(ierr)
  endif

  ! Allocate vectors to hold velocity at face
  if (realization%output_option%print_hdf5_vel_face) then

    ! vx
    call VecCreateMPI(option%mycomm, &
        (option%nflowspec*MAX_FACE_PER_CELL+1)*realization%patch%grid%nlmax, &
        PETSC_DETERMINE,field%vx_face_inst,ierr);CHKERRQ(ierr)
    call VecSet(field%vx_face_inst,0.d0,ierr);CHKERRQ(ierr)

    ! vy and vz
    call VecDuplicate(field%vx_face_inst,field%vy_face_inst, &
                      ierr);CHKERRQ(ierr)
    call VecDuplicate(field%vx_face_inst,field%vz_face_inst, &
                      ierr);CHKERRQ(ierr)
  endif

  if (realization%output_option%print_explicit_flowrate) then
    call VecCreateMPI(option%mycomm, &
         size(grid%unstructured_grid%explicit_grid%connections,2), &
         PETSC_DETERMINE,field%flowrate_inst,ierr);CHKERRQ(ierr)
    call VecSet(field%flowrate_inst,0.d0,ierr);CHKERRQ(ierr)
  endif
    
  ! If average flowrate has to be saved, create a vector for it
  if (realization%output_option%print_hdf5_aveg_mass_flowrate.or. &
      realization%output_option%print_hdf5_aveg_energy_flowrate) then
    call VecCreateMPI(option%mycomm, &
        (option%nflowdof*MAX_FACE_PER_CELL+1)*realization%patch%grid%nlmax, &
        PETSC_DETERMINE,field%flowrate_aveg,ierr);CHKERRQ(ierr)
    call VecSet(field%flowrate_aveg,0.d0,ierr);CHKERRQ(ierr)
  endif

  select case(realization%discretization%itype)
    case(STRUCTURED_GRID)
      realization%comm1 => StructuredCommunicatorCreate()
    case(UNSTRUCTURED_GRID)
      realization%comm1 => UnstructuredCommunicatorCreate()
  end select
  call realization%comm1%SetDM(discretization%dm_1dof)

  if (option%flow%quasi_3d) then
    call RealizCreateFlowMassTransferVec(realization)
  endif

end subroutine RealizationCreateDiscretization

! ************************************************************************** !

subroutine RealizationLocalizeRegions(realization)
  ! 
  ! Localizes regions within each patch
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  use Option_module
  use String_module
  use Grid_module

  implicit none
  
  class(realization_subsurface_type) :: realization
  
  type (region_type), pointer :: cur_region, cur_region2
  type(option_type), pointer :: option
  type(region_type), pointer :: region

  option => realization%option

  ! check to ensure that region names are not duplicated
  cur_region => realization%region_list%first
  do
    if (.not.associated(cur_region)) exit
    cur_region2 => cur_region%next
    do
      if (.not.associated(cur_region2)) exit
      if (StringCompare(cur_region%name,cur_region2%name,MAXWORDLENGTH)) then
        option%io_buffer = 'Duplicate region names: ' // trim(cur_region%name)
        call PrintErrMsg(option)
      endif
      cur_region2 => cur_region2%next
    enddo
    cur_region => cur_region%next
  enddo

  call PatchLocalizeRegions(realization%patch,realization%region_list, &
                            realization%option)
  ! destroy realization's copy of region list as it can be confused with the
  ! localized patch regions later in teh simulation.
  call RegionDestroyList(realization%region_list)

  ! compute regional connections for inline surface flow
  if (option%inline_surface_flow) then
     region => RegionGetPtrFromList(option%inline_surface_region_name, &
          realization%patch%region_list)
     if (.not.associated(region)) then
        option%io_buffer = 'realization_subsurface.F90:RealizationLocalize&
             &Regions() --> Could not find a required region named "' // &
             trim(option%inline_surface_region_name) // &
             '" from the list of regions.'
        call PrintErrMsg(option)
     endif
     call GridRestrictRegionalConnect(realization%patch%grid,region)
   endif
   
end subroutine RealizationLocalizeRegions

! ************************************************************************** !

subroutine RealizationPassPtrsToPatches(realization)
  ! 
  ! Sets patch%field => realization%field
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/11
  ! 

  use Option_module

  implicit none
  
  class(realization_subsurface_type) :: realization
  
  realization%patch%field => realization%field
  realization%patch%datasets => realization%datasets
  realization%patch%reaction => realization%reaction
  realization%patch%reaction_nw => realization%reaction_nw
  realization%patch%reaction_base => realization%reaction_base
  
end subroutine RealizationPassPtrsToPatches

! ************************************************************************** !

subroutine RealizationAddCoupler(realization,coupler)
  ! 
  ! Adds a copy of a coupler to a list
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  use Coupler_module

  implicit none
  
  class(realization_subsurface_type) :: realization
  type(coupler_type), pointer :: coupler
  
  type(patch_type), pointer :: patch
  
  type(coupler_type), pointer :: new_coupler
  
  patch => realization%patch
  
  ! only add to flow list for now, since they will be split out later
  new_coupler => CouplerCreate(coupler)
  select case(coupler%itype)
    case(BOUNDARY_COUPLER_TYPE)
      call CouplerAddToList(new_coupler,patch%boundary_condition_list)
    case(INITIAL_COUPLER_TYPE)
      call CouplerAddToList(new_coupler,patch%initial_condition_list)
    case(SRC_SINK_COUPLER_TYPE)
      call CouplerAddToList(new_coupler,patch%source_sink_list)
  end select
  nullify(new_coupler)

  call CouplerDestroy(coupler)
 
end subroutine RealizationAddCoupler

! ************************************************************************** !

subroutine RealizationAddStrata(realization,strata)
  ! 
  ! Adds a copy of a strata to a list
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  use Strata_module

  implicit none
  
  class(realization_subsurface_type) :: realization
  type(strata_type), pointer :: strata
  
  type(strata_type), pointer :: new_strata
  
  new_strata => StrataCreate(strata)
  call StrataAddToList(new_strata,realization%patch%strata_list)
  nullify(new_strata)
  
  call StrataDestroy(strata)
 
end subroutine RealizationAddStrata


! ************************************************************************** !

subroutine RealizationProcessCouplers(realization)
  ! 
  ! Sets connectivity and pointers for couplers
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  use Option_module

  implicit none
  
  class(realization_subsurface_type) :: realization
  
  call PatchProcessCouplers( realization%patch,realization%flow_conditions, &
                             realization%transport_conditions, &
                             realization%option)
  
end subroutine RealizationProcessCouplers

! ************************************************************************** !

subroutine RealizationProcessDatasets(realization)
  ! 
  ! Processes datasets before they are linked to anything else
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/20/16
  ! 
  use Dataset_module
  
  implicit none
  
  class(realization_subsurface_type) :: realization

  call DatasetScreenForNonCellIndexed(realization%datasets,realization%option)
  
end subroutine RealizationProcessDatasets

! ************************************************************************** !

subroutine RealizationProcessConditions(realization)
  ! 
  ! Sets up auxiliary data associated with
  ! conditions
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  ! 
  use Data_Mediator_Base_class
  use Data_Mediator_Dataset_class
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  class(data_mediator_base_type), pointer :: cur_data_mediator

  if (realization%option%nflowdof > 0) then
    call RealProcessFlowConditions(realization)
  endif
  if (realization%option%ntrandof > 0) then
    call RealProcessTranConditions(realization)
  endif
  
  ! update data mediators
  cur_data_mediator => realization%flow_data_mediator_list
  do
    if (.not.associated(cur_data_mediator)) exit
    call RealizCreateFlowMassTransferVec(realization)
    select type(cur_data_mediator)
      class is(data_mediator_dataset_type)
        call DataMediatorDatasetInit(cur_data_mediator, &
                                     realization%discretization, &
                                     realization%datasets, &
                                     realization%option)
        call cur_data_mediator%Update(realization%field%flow_mass_transfer, &
                                      realization%option)
      class default
    end select
    cur_data_mediator => cur_data_mediator%next
  enddo

  cur_data_mediator => realization%tran_data_mediator_list
  do
    if (.not.associated(cur_data_mediator)) exit
    call RealizCreateTranMassTransferVec(realization)
    select type(cur_data_mediator)
      class is(data_mediator_dataset_type)
        call DataMediatorDatasetInit(cur_data_mediator, &
                                     realization%discretization, &
                                     realization%datasets, &
                                     realization%option)
        call cur_data_mediator%Update(realization%field%tran_mass_transfer, &
                                      realization%option)
      class default
    end select
    cur_data_mediator => cur_data_mediator%next
  enddo

end subroutine RealizationProcessConditions

! ************************************************************************** !

subroutine RealProcessMatPropAndSatFunc(realization)
  ! 
  ! Sets up linkeage between material properties
  ! and saturation function, auxiliary arrays
  ! and datasets
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09, 01/12/11
  ! 

  use String_module
  use Dataset_Common_HDF5_class
  use Dataset_module
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  
  PetscBool :: found
  PetscInt :: i
  type(option_type), pointer :: option
  type(material_property_type), pointer :: cur_material_property
  type(patch_type), pointer :: patch
  character(len=MAXSTRINGLENGTH) :: string
  class(dataset_base_type), pointer :: dataset

  option => realization%option
  patch => realization%patch
  
  ! set up mirrored pointer arrays within patches to saturation functions
  ! and material properties
  patch%material_properties => realization%material_properties
  call MaterialPropConvertListToArray(patch%material_properties, &
                                      patch%material_property_array, &
                                      option)
  if (associated(realization%saturation_functions)) then
    patch%saturation_functions => realization%saturation_functions
    call SaturatFuncConvertListToArray(patch%saturation_functions, &
                                       patch%saturation_function_array, &
                                       option)
  endif
  if (associated(realization%characteristic_curves)) then
    patch%characteristic_curves => realization%characteristic_curves
    call CharCurvesConvertListToArray(patch%characteristic_curves, &
                                      patch%characteristic_curves_array, &
                                      option)
  endif

  if (option%iflowmode == G_MODE) then
  ! set up analogous mapping to thermal characteristic curves, if used    
    if (associated(realization%thermal_characteristic_curves)) then
      patch%thermal_characteristic_curves => &
           realization%thermal_characteristic_curves
      call CharCurvesThermalConvertListToArray( &
           patch%thermal_characteristic_curves, &
           patch%thermal_characteristic_curves_array, option)
    end if
  end if
  
  ! create mapping of internal to external material id
  call MaterialCreateIntToExtMapping(patch%material_property_array, &
                                     patch%imat_internal_to_external)

  cur_material_property => realization%material_properties
  do
    if (.not.associated(cur_material_property)) exit

    ! obtain saturation function id
    if (option%iflowmode /= NULL_MODE) then
      if (associated(patch%saturation_function_array)) then
        cur_material_property%saturation_function_id = &
          SaturationFunctionGetID(patch%saturation_functions, &
                             cur_material_property%saturation_function_name, &
                             cur_material_property%name,option)
      endif
      if (associated(patch%characteristic_curves_array)) then
        cur_material_property%saturation_function_id = &
          CharacteristicCurvesGetID(patch%characteristic_curves_array, &
                             cur_material_property%saturation_function_name, &
                             cur_material_property%name,option)
      endif
      if (cur_material_property%saturation_function_id == 0) then
        option%io_buffer = 'Characteristic curve "' // &
          trim(cur_material_property%saturation_function_name) // &
          '" not found.'
        call PrintErrMsg(option)
      else
        if (associated(patch%characteristic_curves_array)) then
          call CharCurvesProcessTables(patch%characteristic_curves_array(  &
                        cur_material_property%saturation_function_id)%ptr,option)
        end if
      end if
    endif

    ! thermal conducitivity function id 
    if (option%iflowmode == G_MODE) then
      if (associated(patch%thermal_characteristic_curves_array)) then
        cur_material_property%thermal_conductivity_function_id = &
             CharacteristicCurvesThermalGetID( &
             patch%thermal_characteristic_curves_array, &
             cur_material_property%thermal_conductivity_function_name, &
             cur_material_property%name,option)
      end if
      if (cur_material_property%thermal_conductivity_function_id == 0) then
        option%io_buffer = 'Thermal characteristic curve "' // &
          trim(cur_material_property%thermal_conductivity_function_name) // &
          '" not found.'
        call PrintErrMsg(option)
      end if
    end if

    ! if named, link dataset to property
    if (associated(cur_material_property%porosity_dataset)) then
      string = 'MATERIAL_PROPERTY(' // trim(cur_material_property%name) // &
               '),POROSITY'
      dataset => &
        DatasetBaseGetPointer(realization%datasets, &
                              cur_material_property%porosity_dataset%name, &
                              string,option)
      call DatasetDestroy(cur_material_property%porosity_dataset)
      select type(dataset)
        class is (dataset_common_hdf5_type)
          cur_material_property%porosity_dataset => dataset
        class default
          option%io_buffer = 'Incorrect dataset type for porosity.'
          call PrintErrMsg(option)
      end select
    endif
    if (associated(cur_material_property%tortuosity_dataset)) then
      string = 'MATERIAL_PROPERTY(' // trim(cur_material_property%name) // &
               '),TORTUOSITY'
      dataset => &
        DatasetBaseGetPointer(realization%datasets, &
                              cur_material_property%tortuosity_dataset%name, &
                              string,option)
      call DatasetDestroy(cur_material_property%tortuosity_dataset)
      select type(dataset)
        class is (dataset_common_hdf5_type)
          cur_material_property%tortuosity_dataset => dataset
        class default
          option%io_buffer = 'Incorrect dataset type for tortuosity.'
          call PrintErrMsg(option)
      end select
    endif
    if (associated(cur_material_property%permeability_dataset)) then
      string = 'MATERIAL_PROPERTY(' // trim(cur_material_property%name) // &
               '),PERMEABILITY or PERMEABILITY X'
      dataset => &
        DatasetBaseGetPointer(realization%datasets, &
                            cur_material_property%permeability_dataset%name, &
                            string,option)
      call DatasetDestroy(cur_material_property%permeability_dataset)
      select type(dataset)
        class is (dataset_common_hdf5_type)
          cur_material_property%permeability_dataset => dataset
        class default
          option%io_buffer = 'Incorrect dataset type for permeability or &
                             &permeability X.'
          call PrintErrMsg(option)
      end select      
    endif
    if (associated(cur_material_property%permeability_dataset_y)) then
      string = 'MATERIAL_PROPERTY(' // trim(cur_material_property%name) // &
               '),PERMEABILITY Y'
      dataset => &
        DatasetBaseGetPointer(realization%datasets, &
                            cur_material_property%permeability_dataset_y%name, &
                            string,option)
      call DatasetDestroy(cur_material_property%permeability_dataset_y)
      select type(dataset)
        class is (dataset_common_hdf5_type)
          cur_material_property%permeability_dataset_y => dataset
        class default
          option%io_buffer = 'Incorrect dataset type for permeability Y.'
          call PrintErrMsg(option)
      end select      
    endif
    if (associated(cur_material_property%permeability_dataset_z)) then
      string = 'MATERIAL_PROPERTY(' // trim(cur_material_property%name) // &
               '),PERMEABILITY Z'
      dataset => &
        DatasetBaseGetPointer(realization%datasets, &
                            cur_material_property%permeability_dataset_z%name, &
                            string,option)
      call DatasetDestroy(cur_material_property%permeability_dataset_z)
      select type(dataset)
        class is (dataset_common_hdf5_type)
          cur_material_property%permeability_dataset_z => dataset
        class default
          option%io_buffer = 'Incorrect dataset type for permeability Z.'
          call PrintErrMsg(option)
      end select      
    endif
    if (associated(cur_material_property%soil_reference_pressure_dataset)) then
      string = 'MATERIAL_PROPERTY(' // trim(cur_material_property%name) // &
               '),SOIL_REFERENCE_PRESSURE'
      dataset => &
        DatasetBaseGetPointer(realization%datasets, &
                 cur_material_property%soil_reference_pressure_dataset%name, &
                 string,option)
      call DatasetDestroy(cur_material_property%soil_reference_pressure_dataset)
      select type(dataset)
        class is (dataset_common_hdf5_type)
          cur_material_property%soil_reference_pressure_dataset => dataset
        class default
          option%io_buffer = 'Incorrect dataset type for soil reference &
                              &pressure.'
          call PrintErrMsg(option)
      end select
    endif
    if (associated(cur_material_property%compressibility_dataset)) then
      string = 'MATERIAL_PROPERTY(' // trim(cur_material_property%name) // &
               '),SOIL_COMPRESSIBILITY'
      dataset => &
        DatasetBaseGetPointer(realization%datasets, &
                         cur_material_property%compressibility_dataset%name, &
                         string,option)
      call DatasetDestroy(cur_material_property%compressibility_dataset)
      select type(dataset)
        class is (dataset_common_hdf5_type)
          cur_material_property%compressibility_dataset => dataset
        class default
          option%io_buffer = 'Incorrect dataset type for soil_compressibility.'
          call PrintErrMsg(option)
      end select      
    endif
    
    cur_material_property => cur_material_property%next
  enddo
  
  
end subroutine RealProcessMatPropAndSatFunc

! ************************************************************************** !

subroutine RealProcessFluidProperties(realization)
  ! 
  ! Sets up linkeage with fluid properties
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09
  ! 

  use Grid_Grdecl_module, only : GetSatnumSet

  implicit none

  class(realization_subsurface_type) :: realization

  PetscBool :: found
  type(option_type), pointer :: option
  type(fluid_property_type), pointer :: cur_fluid_property
  PetscInt :: icc, ncc, maxsatn
  PetscBool :: satnum_set, ccset

  option => realization%option

  found = PETSC_FALSE
  cur_fluid_property => realization%fluid_properties                            
  do                                      
    if (.not.associated(cur_fluid_property)) exit
    found = PETSC_TRUE
    select case(trim(cur_fluid_property%phase_name))
      case('LIQUID')
        cur_fluid_property%phase_id = LIQUID_PHASE
      case('GAS')
        cur_fluid_property%phase_id = GAS_PHASE
      case default
        cur_fluid_property%phase_id = LIQUID_PHASE
    end select
    cur_fluid_property => cur_fluid_property%next
  enddo

  if (option%ntrandof > 0 .and. .not.found) then
    option%io_buffer = 'A fluid property must be present in input file' // &
                       ' for solute transport'
  endif

  ! If saturation table numbers set,
  ! check that matches characteristic curves count

  satnum_set = GetSatnumSet(maxsatn)
  if( satnum_set ) then
    ccset = associated(realization%patch%characteristic_curves_array)
    if (ccset) then
      ncc = size(realization%patch%characteristic_curves_array(:))
      if( maxsatn > ncc ) then
        option%io_buffer = &
         'SATNUM data does not match CHARACTERISTIC CURVES count'
        call PrintErrMsg(option)
      endif
      do icc = 1, ncc
        call CharCurvesProcessTables( &
          realization%patch%characteristic_curves_array(icc)%ptr,option)
      end do
    else
      option%io_buffer = 'SATNUM data but no CHARACTERISTIC CURVES'
      call PrintErrMsg(option)
    end if
  endif

end subroutine RealProcessFluidProperties

! ************************************************************************** !

subroutine RealProcessFlowConditions(realization)
  ! 
  ! Sets linkage of flow conditions to dataset
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/11
  ! 

  use Dataset_Base_class
  use Dataset_module

  implicit none

  class(realization_subsurface_type) :: realization
  
  type(flow_condition_type), pointer :: cur_flow_condition
  type(flow_sub_condition_type), pointer :: cur_flow_sub_condition
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i
  
  option => realization%option
  
  ! loop over flow conditions looking for linkage to datasets
  cur_flow_condition => realization%flow_conditions%first
  do
    if (.not.associated(cur_flow_condition)) exit
    string = 'flow_condition ' // trim(cur_flow_condition%name)
    ! find datum dataset
    call DatasetFindInList(realization%datasets,cur_flow_condition%datum, &
                           cur_flow_condition%default_time_storage, &
                           string,option)
    select case(option%iflowmode)
      case default
        do i = 1, size(cur_flow_condition%sub_condition_ptr)
          ! find dataset
          call DatasetFindInList(realization%datasets, &
                 cur_flow_condition%sub_condition_ptr(i)%ptr%dataset, &
                 cur_flow_condition%default_time_storage, &
                 string,option)
          ! find gradient dataset
          call DatasetFindInList(realization%datasets, &
                 cur_flow_condition%sub_condition_ptr(i)%ptr%gradient, &
                 cur_flow_condition%default_time_storage, &
                 string,option)
        enddo
    end select
    cur_flow_condition => cur_flow_condition%next
  enddo

end subroutine RealProcessFlowConditions

! ************************************************************************** !

subroutine RealProcessTranConditions(realization)
  ! 
  ! Sets up auxiliary data associated with
  ! transport conditions
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/14/08
  ! 

  use Reaction_module
  use String_module
  use Transport_Constraint_Base_module
  use Transport_Constraint_NWT_module
  use Transport_Constraint_RT_module
  use Transport_Constraint_module
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  
  PetscBool :: found, coupling_needed
  type(option_type), pointer :: option
  type(tran_condition_type), pointer :: cur_condition
  class(tran_constraint_coupler_base_type), pointer :: cur_constraint_coupler
  class(tran_constraint_base_type), pointer :: cur_constraint, &
                                              another_constraint
  
  option => realization%option
  coupling_needed = PETSC_FALSE
  
  ! check for duplicate constraint names
  cur_constraint => realization%transport_constraints%first
  do
    if (.not.associated(cur_constraint)) exit
      another_constraint => cur_constraint%next
      ! now compare names
      found = PETSC_FALSE
      do
        if (.not.associated(another_constraint)) exit
        if (StringCompare(cur_constraint%name,another_constraint%name, &
            MAXWORDLENGTH)) then
          found = PETSC_TRUE
        endif
        another_constraint => another_constraint%next
      enddo
      if (found) then
        option%io_buffer = 'Duplicate transport constraints named "' // &
                 trim(cur_constraint%name) // '"'
        call PrintErrMsg(realization%option)
      endif
    cur_constraint => cur_constraint%next
  enddo
  
  ! initialize constraints
  cur_constraint => realization%transport_constraints%first
  do
    if (.not.associated(cur_constraint)) exit
    select type(constraint=>cur_constraint)
      class is (tran_constraint_rt_type)
        call ReactionProcessConstraint(realization%reaction, &
                                       constraint,realization%option)
      class is (tran_constraint_nwt_type)
        call NWTConstraintProcess(realization%reaction_nw, &
                                  constraint,realization%option)
    end select
    cur_constraint => cur_constraint%next
  enddo
  
  if (option%use_mc) then
    select type(constraint=>realization%sec_transport_constraint)
      class is (tran_constraint_rt_type)
        call ReactionProcessConstraint(realization%reaction, &
                                       constraint,realization%option)
    end select
  endif
  
  ! tie constraints to couplers, if not already associated
  cur_condition => realization%transport_conditions%first
  do
    if (.not.associated(cur_condition)) exit
    cur_constraint_coupler => cur_condition%constraint_coupler_list
    do
      if (.not.associated(cur_constraint_coupler)) exit
      ! if constraint exists, it was coupled during the embedded read.
      if (.not.associated(cur_constraint_coupler%constraint)) then
        cur_constraint => realization%transport_constraints%first
        do
          if (.not.associated(cur_constraint)) exit
          if (StringCompare(cur_constraint%name, &
                            cur_constraint_coupler%constraint_name, &
                            MAXWORDLENGTH)) then
            cur_constraint_coupler%constraint => cur_constraint
            cur_constraint_coupler%equilibrate_at_each_cell = &
              (realization%reaction_base%equilibrate_at_each_cell .or. &
               cur_constraint%equilibrate_at_each_cell)
            exit
          endif
          cur_constraint => cur_constraint%next
        enddo
        if (.not.associated(cur_constraint_coupler%constraint)) then
          option%io_buffer = 'Transport constraint "' // &
                   trim(cur_constraint_coupler%constraint_name) // &
                   '" not found in input file constraints.'
          call PrintErrMsg(realization%option)
        endif
      endif
      cur_constraint_coupler => cur_constraint_coupler%next
    enddo
!TODO(geh) remove this?
    if (associated(cur_condition%constraint_coupler_list%next)) then 
      ! there are more than one
      cur_condition%is_transient = PETSC_TRUE
    else
      cur_condition%is_transient = PETSC_FALSE
    endif
    cur_condition => cur_condition%next
  enddo
 
  ! final details for setup
  cur_condition => realization%transport_conditions%first
  do
    if (.not.associated(cur_condition)) exit
    ! is the condition transient?
    if (associated(cur_condition%constraint_coupler_list%next)) then 
      ! there are more than one
      cur_condition%is_transient = PETSC_TRUE
    else
      cur_condition%is_transient = PETSC_FALSE
    endif
    ! set pointer to first constraint coupler
    cur_condition%cur_constraint_coupler => &
                         cur_condition%constraint_coupler_list
    cur_condition => cur_condition%next
  enddo

end subroutine RealProcessTranConditions

! ************************************************************************** !

subroutine RealizationInitConstraints(realization)
  ! 
  ! Initializes constraint concentrations
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/04/08
  ! 

  implicit none

  class(realization_subsurface_type) :: realization
  
  type(patch_type), pointer :: cur_patch
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit
    call PatchInitConstraints(cur_patch,realization%reaction_base, &
                              realization%option)
    cur_patch => cur_patch%next
  enddo
 
end subroutine RealizationInitConstraints

! ************************************************************************** !

subroutine RealizationPrintCouplers(realization)
  ! 
  ! Print boundary and initial condition data
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/28/08
  ! 

  use Coupler_module
  use Reaction_Aux_module
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  
  type(patch_type), pointer :: cur_patch
  type(coupler_type), pointer :: cur_coupler
  type(option_type), pointer :: option
  class(reaction_rt_type), pointer :: reaction
 
  option => realization%option
  reaction => realization%reaction
 
  if (.not.OptionPrintToFile(option)) return
  
  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit

    cur_coupler => cur_patch%initial_condition_list%first
    do
      if (.not.associated(cur_coupler)) exit
      call RealizationPrintCoupler(cur_coupler,reaction,option)    
      cur_coupler => cur_coupler%next
    enddo
     
    cur_coupler => cur_patch%boundary_condition_list%first
    do
      if (.not.associated(cur_coupler)) exit
      call RealizationPrintCoupler(cur_coupler,reaction,option)    
      cur_coupler => cur_coupler%next
    enddo
     
    cur_coupler => cur_patch%source_sink_list%first
    do
      if (.not.associated(cur_coupler)) exit
      call RealizationPrintCoupler(cur_coupler,reaction,option)    
      cur_coupler => cur_coupler%next
    enddo

    cur_patch => cur_patch%next
  enddo
    
end subroutine RealizationPrintCouplers

! ************************************************************************** !

subroutine RealizationPrintCoupler(coupler,reaction,option)
  ! 
  ! Prints boundary and initial condition coupler
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/28/08
  ! 
  use Coupler_module
  use Reaction_module
  use Reaction_Aux_module
  use Transport_Constraint_Base_module
  use Transport_Constraint_RT_module
  
  implicit none
  
  type(coupler_type) :: coupler
  type(option_type) :: option
  class(reaction_rt_type), pointer :: reaction
  
  character(len=MAXSTRINGLENGTH) :: string
  
  type(flow_condition_type), pointer :: flow_condition
  type(tran_condition_type), pointer :: tran_condition
  type(region_type), pointer :: region
  class(tran_constraint_coupler_base_type), pointer :: constraint_coupler
   
98 format(40('=+'))
99 format(80('-'))
  
  flow_condition => coupler%flow_condition
  tran_condition => coupler%tran_condition
  region => coupler%region

  write(option%fid_out,*)
  write(option%fid_out,98)


  select case(coupler%itype)
    case(INITIAL_COUPLER_TYPE)
      string = 'Initial Condition'
    case(BOUNDARY_COUPLER_TYPE)
      string = 'Boundary Condition'
    case(SRC_SINK_COUPLER_TYPE)
      string = 'Source Sink'
  end select
  write(option%fid_out,'(/,2x,a,/)') trim(string)

  write(option%fid_out,99)
101 format(5x,'     Flow Condition: ',2x,a)
  if (associated(flow_condition)) &
    write(option%fid_out,101) trim(flow_condition%name)
102 format(5x,'Transport Condition: ',2x,a)
  if (associated(tran_condition)) &
    write(option%fid_out,102) trim(tran_condition%name)
103 format(5x,'             Region: ',2x,a)
  if (associated(region)) &
    write(option%fid_out,103) trim(region%name)
  write(option%fid_out,99)
  
  if (associated(flow_condition)) then
    call FlowConditionPrint(flow_condition,option)
  endif
  if (associated(tran_condition)) then
    constraint_coupler => tran_condition%cur_constraint_coupler
    write(option%fid_out,'(/,2x,''Transport Condition: '',a)') &
      trim(tran_condition%name)
    select type(c=>constraint_coupler)
      class is (tran_constraint_coupler_rt_type)
        call ReactionPrintConstraint(c,reaction,option)
        write(option%fid_out,'(/)')
        write(option%fid_out,99)
    end select
  endif
 
end subroutine RealizationPrintCoupler

! ************************************************************************** !

subroutine RealizationInitAllCouplerAuxVars(realization)
  ! 
  ! RealizationInitCouplerAuxVars: Initializes coupler auxillary variables
  ! within list
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  use Option_module

  implicit none
  
  class(realization_subsurface_type) :: realization
  
  !geh: Must update conditions prior to initializing the aux vars.  
  !     Otherwise, datasets will not have been read for routines such as
  !     hydrostatic and auxvars will be initialized to garbage.
  call FlowConditionUpdate(realization%flow_conditions,realization%option)
  call TranConditionUpdate(realization%transport_conditions, &
                           realization%option)
  call PatchInitAllCouplerAuxVars(realization%patch,realization%option)
   
end subroutine RealizationInitAllCouplerAuxVars

! ************************************************************************** !

subroutine RealizUpdateAllCouplerAuxVars(realization,force_update_flag)
  ! 
  ! Updates auxiliary variables associated
  ! with couplers in lis
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  use Option_module

  implicit none
  
  class(realization_subsurface_type) :: realization
  PetscBool :: force_update_flag

  !TODO(geh): separate flow from transport in these calls
  call PatchUpdateAllCouplerAuxVars(realization%patch,force_update_flag, &
                                    realization%option)

end subroutine RealizUpdateAllCouplerAuxVars

! ************************************************************************** !

subroutine RealizationRevertFlowParameters(realization)
  ! 
  ! Overwrites porosity/permeability in materials_auxvars with values stored in 
  ! Vecs
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/09/08
  ! 

  use Option_module
  use Field_module
  use Discretization_module
  use Material_Aux_class, only : material_type, &
                              POROSITY_CURRENT, POROSITY_BASE, POROSITY_INITIAL
  use Variables_module, only : PERMEABILITY_X, PERMEABILITY_Y, PERMEABILITY_Z, &
                               POROSITY

  implicit none
  
  class(realization_subsurface_type) :: realization
  
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(material_type), pointer :: Material
  
  option => realization%option
  field => realization%field
  discretization => realization%discretization
  Material => realization%patch%aux%Material

  if (option%nflowdof > 0) then
    call DiscretizationGlobalToLocal(discretization,field%perm0_xx, &
                                     field%work_loc,ONEDOF)  
    call MaterialSetAuxVarVecLoc(Material,field%work_loc,PERMEABILITY_X, &
                                 ZERO_INTEGER)
    call DiscretizationGlobalToLocal(discretization,field%perm0_yy, &
                                     field%work_loc,ONEDOF)  
    call MaterialSetAuxVarVecLoc(Material,field%work_loc,PERMEABILITY_Y, &
                                 ZERO_INTEGER)
    call DiscretizationGlobalToLocal(discretization,field%perm0_zz, &
                                     field%work_loc,ONEDOF)  
    call MaterialSetAuxVarVecLoc(Material,field%work_loc,PERMEABILITY_Z, &
                                 ZERO_INTEGER)
  endif   
  call DiscretizationGlobalToLocal(discretization,field%porosity0, &
                                   field%work_loc,ONEDOF)  
  call MaterialSetAuxVarVecLoc(Material,field%work_loc,POROSITY, &
                               POROSITY_INITIAL)
  call MaterialSetAuxVarVecLoc(Material,field%work_loc,POROSITY, &
                               POROSITY_BASE)
  call MaterialSetAuxVarVecLoc(Material,field%work_loc,POROSITY, &
                               POROSITY_CURRENT)
  ! tortuosity is not currently checkpointed
!  call DiscretizationGlobalToLocal(discretization,field%tortuosity0, &
!                                   field%work_loc,ONEDOF)  
!  call MaterialSetAuxVarVecLoc(Material,field%work_loc,TORTUOSITY, &
!                               ZERO_INTEGER)

end subroutine RealizationRevertFlowParameters

! ************************************************************************** !

subroutine RealizStoreRestartFlowParams(realization)
  ! 
  ! Overwrites porosity/permeability Vecs with restart values stored in 
  ! material_auxvars
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/09/08
  ! 

  use Option_module
  use Field_module
  use Discretization_module
  use Material_Aux_class
  use Variables_module

  implicit none
  
  class(realization_subsurface_type) :: realization
  
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  type(material_type), pointer :: Material
  
  option => realization%option
  field => realization%field
  discretization => realization%discretization
  Material => realization%patch%aux%Material

  if (option%nflowdof > 0) then
    call MaterialGetAuxVarVecLoc(Material,field%work_loc,PERMEABILITY_X, &
                                 ZERO_INTEGER)
    call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                     field%perm0_xx,ONEDOF)
    call MaterialGetAuxVarVecLoc(Material,field%work_loc,PERMEABILITY_Y, &
                                 ZERO_INTEGER)
    call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                     field%perm0_yy,ONEDOF)
    call MaterialGetAuxVarVecLoc(Material,field%work_loc,PERMEABILITY_Z, &
                                 ZERO_INTEGER)
    call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                     field%perm0_zz,ONEDOF)
  endif   
  call MaterialGetAuxVarVecLoc(Material,field%work_loc,POROSITY, &
                               POROSITY_BASE)
  ! might as well update initial and base at the same time
  call MaterialSetAuxVarVecLoc(Material,field%work_loc,POROSITY, &
                               POROSITY_INITIAL)
  call MaterialSetAuxVarVecLoc(Material,field%work_loc,POROSITY, &
                               POROSITY_CURRENT)
  call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                   field%porosity0,ONEDOF)
  ! tortuosity is not currently checkpointed
!  call DiscretizationLocalToGlobal(discretization,field%work_loc, &
!                                   field%tortuosity0,ONEDOF)
!  call MaterialGetAuxVarVecLoc(Material,field%work_loc,TORTUOSITY, &
!                               ZERO_INTEGER)

end subroutine RealizStoreRestartFlowParams

! ************************************************************************** !

subroutine RealizUpdateUniformVelocity(realization)
  ! 
  ! Assigns uniform velocity for transport
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  use Option_module
  use Dataset_module

  implicit none
  
  class(realization_subsurface_type) :: realization
  
  call DatasetUpdate(realization%uniform_velocity_dataset, &
                     realization%option)
  call PatchUpdateUniformVelocity(realization%patch, &
                            realization%uniform_velocity_dataset%rarray, &
                            realization%option)
 
end subroutine RealizUpdateUniformVelocity

! ************************************************************************** !

subroutine RealizationAddWaypointsToList(realization,waypoint_list)
  ! 
  ! Creates waypoints associated with source/sinks
  ! boundary conditions, etc. and add to list
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  use Option_module
  use Waypoint_module
  use Time_Storage_module
  use Data_Mediator_Base_class
  use Data_Mediator_Dataset_class
  use Strata_module

  implicit none
  
  class(realization_subsurface_type) :: realization
  type(waypoint_list_type) :: waypoint_list

  type(flow_condition_type), pointer :: cur_flow_condition
  type(tran_condition_type), pointer :: cur_tran_condition
  type(flow_sub_condition_type), pointer :: sub_condition
  class(tran_constraint_coupler_base_type), pointer :: cur_constraint_coupler
  class(data_mediator_base_type), pointer :: cur_data_mediator
  type(waypoint_type), pointer :: waypoint, cur_waypoint
  type(option_type), pointer :: option
  type(strata_type), pointer :: cur_strata
  type(time_storage_type), pointer :: time_storage_ptr
  PetscInt :: itime, isub_condition
  PetscReal :: temp_real, final_time
  PetscReal, pointer :: times(:)

  option => realization%option
  nullify(times)
  
  ! set flag for final output
  cur_waypoint => waypoint_list%first
  do
    if (.not.associated(cur_waypoint)) exit
    if (cur_waypoint%final) then
      cur_waypoint%print_snap_output = &
        realization%output_option%print_final_snap
      exit
    endif
    cur_waypoint => cur_waypoint%next
  enddo
  ! use final time in conditional below
  if (associated(cur_waypoint)) then
    final_time = cur_waypoint%time
  else
    option%io_buffer = 'Final time not found in RealizationAddWaypointsToList'
    call PrintErrMsg(option)
  endif

  ! add update of flow conditions
  cur_flow_condition => realization%flow_conditions%first
  do
    if (.not.associated(cur_flow_condition)) exit
    if (cur_flow_condition%sync_time_with_update) then
      do isub_condition = 1, cur_flow_condition%num_sub_conditions
        sub_condition => cur_flow_condition%sub_condition_ptr(isub_condition)%ptr
        !TODO(geh): check if this updated more than simply the flow_dataset (i.e. datum and gradient)
        !geh: followup - no, datum/gradient are not considered.  Should they be considered?
        call TimeStorageGetTimes(sub_condition%dataset%time_storage, option, &
                                final_time, times)
        if (associated(times)) then
          if (size(times) > 1000) then
            option%io_buffer = 'For flow condition "' // &
              trim(cur_flow_condition%name) // &
              '" dataset "' // trim(sub_condition%name) // &
              '", the number of times is excessive for synchronization ' // &
              'with waypoints.'
            call PrintErrMsg(option)
          endif
          do itime = 1, size(times)
            waypoint => WaypointCreate()
            waypoint%time = times(itime)
            waypoint%update_conditions = PETSC_TRUE
            call WaypointInsertInList(waypoint,waypoint_list)
          enddo
          deallocate(times)
          nullify(times)
        endif
      enddo
    endif
    cur_flow_condition => cur_flow_condition%next
  enddo
      
  ! add update of transport conditions
  cur_tran_condition => realization%transport_conditions%first
  do
    if (.not.associated(cur_tran_condition)) exit
    if (cur_tran_condition%is_transient) then
      cur_constraint_coupler => cur_tran_condition%constraint_coupler_list
      do
        if (.not.associated(cur_constraint_coupler)) exit
        if (cur_constraint_coupler%time > 1.d-40) then
          waypoint => WaypointCreate()
          waypoint%time = cur_constraint_coupler%time
          waypoint%update_conditions = PETSC_TRUE
          call WaypointInsertInList(waypoint,waypoint_list)
        endif
        cur_constraint_coupler => cur_constraint_coupler%next
      enddo
    endif
    cur_tran_condition => cur_tran_condition%next
  enddo

  ! add update of velocity fields
  if (associated(realization%uniform_velocity_dataset)) then
    time_storage_ptr => realization%uniform_velocity_dataset%time_storage
    if (associated(time_storage_ptr)) then
      if (time_storage_ptr%times(1) > 1.d-40 .or. &
          time_storage_ptr%max_time_index > 1) then
        do itime = 1, size(time_storage_ptr%times)
          waypoint => WaypointCreate()
          waypoint%time = time_storage_ptr%times(itime)
          waypoint%update_conditions = PETSC_TRUE
          call WaypointInsertInList(waypoint,waypoint_list)
        enddo
      endif
    endif
  endif
  
  ! add waypoints for flow mass transfer
  if (associated(realization%flow_data_mediator_list)) then
    cur_data_mediator => realization%flow_data_mediator_list
    do
      if (.not.associated(cur_data_mediator)) exit
      select type(cur_data_mediator)
        class is(data_mediator_dataset_type)
          time_storage_ptr => cur_data_mediator%dataset%time_storage
          if (associated(time_storage_ptr)) then
            do itime = 1, time_storage_ptr%max_time_index
              waypoint => WaypointCreate()
              waypoint%time = time_storage_ptr%times(itime)
              waypoint%update_conditions = PETSC_TRUE
              call WaypointInsertInList(waypoint,waypoint_list)
            enddo
          endif
        class default
      end select 
      cur_data_mediator => cur_data_mediator%next
    enddo
  endif  

  ! add waypoints for rt mass transfer
  if (associated(realization%tran_data_mediator_list)) then
    cur_data_mediator => realization%tran_data_mediator_list
    do
      if (.not.associated(cur_data_mediator)) exit
      select type(cur_data_mediator)
        class is(data_mediator_dataset_type)
          time_storage_ptr => cur_data_mediator%dataset%time_storage
          if (associated(time_storage_ptr)) then
            do itime = 1, time_storage_ptr%max_time_index
              waypoint => WaypointCreate()
              waypoint%time = time_storage_ptr%times(itime)
              waypoint%update_conditions = PETSC_TRUE
              call WaypointInsertInList(waypoint,waypoint_list)
            enddo
          endif
        class default
      end select           
      cur_data_mediator => cur_data_mediator%next
    enddo
  endif 

  ! add in strata that change over time
  cur_strata => realization%patch%strata_list%first
  do
    if (.not.associated(cur_strata)) exit
    if (Initialized(cur_strata%start_time)) then
      waypoint => WaypointCreate()
      waypoint%time = cur_strata%start_time
      waypoint%sync = PETSC_TRUE
      call WaypointInsertInList(waypoint,waypoint_list)
    endif
    if (Initialized(cur_strata%final_time)) then
      waypoint => WaypointCreate()
      waypoint%time = cur_strata%final_time
      waypoint%sync = PETSC_TRUE
      call WaypointInsertInList(waypoint,waypoint_list)
    endif
    cur_strata => cur_strata%next
  enddo

end subroutine RealizationAddWaypointsToList

! ************************************************************************** !

subroutine RealizationUpdatePropertiesTS(realization)
  ! 
  ! Updates coupled properties at each grid cell
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/05/09
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Discretization_module
  use Field_module
  use Grid_module
  use Material_Aux_class
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Variables_module, only : POROSITY, TORTUOSITY, PERMEABILITY_X, &
                               PERMEABILITY_Y, PERMEABILITY_Z
 
  implicit none
  
  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  class(reaction_rt_type), pointer :: reaction
  type(grid_type), pointer :: grid
  type(material_property_ptr_type), pointer :: material_property_array(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:) 
  type(discretization_type), pointer :: discretization
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: local_id, ghosted_id
  PetscInt :: imnrl, imnrl1, imnrl_armor, imat
  PetscReal :: sum_volfrac
  PetscReal :: scale, porosity_scale, volfrac_scale
  PetscBool :: porosity_updated
  PetscReal, pointer :: vec_p(:)
  PetscReal, pointer :: porosity0_p(:)
  PetscReal, pointer :: tortuosity0_p(:)
  PetscReal, pointer :: perm0_xx_p(:), perm0_yy_p(:), perm0_zz_p(:)
  PetscReal, pointer :: perm_ptr(:)
  PetscReal :: min_value
  PetscReal :: critical_porosity
  PetscReal :: porosity_base_
  PetscInt :: ivalue
  PetscErrorCode :: ierr

  option => realization%option
  discretization => realization%discretization
  patch => realization%patch
  field => realization%field
  reaction => realization%reaction
  grid => patch%grid
  material_property_array => patch%material_property_array
  rt_auxvars => patch%aux%RT%auxvars
  material_auxvars => patch%aux%Material%auxvars

  porosity_updated = PETSC_FALSE
  if (reaction%update_porosity) then
    porosity_updated = PETSC_TRUE
    call RealizationCalcMineralPorosity(realization)
  endif
  
  if ((porosity_updated .and. &
       (reaction%update_tortuosity .or. &
        reaction%update_permeability)) .or. &
      ! if porosity ratio is used in mineral surface area update, we must
      ! recalculate it every time.
      (reaction%update_mineral_surface_area .and. &
       reaction%update_mnrl_surf_with_porosity)) then
    call VecGetArrayF90(field%porosity0,porosity0_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      vec_p(local_id) = material_auxvars(ghosted_id)%porosity_base / &
                        porosity0_p(local_id)
    enddo
    call VecRestoreArrayF90(field%porosity0,porosity0_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
  endif      

  if (reaction%update_mineral_surface_area) then

    if (reaction%update_mnrl_surf_with_porosity) then
      ! placing the get/restore array calls within the condition will
      ! avoid improper access.
      call VecGetArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
    endif

    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      do imnrl = 1, reaction%mineral%nkinmnrl

        porosity_scale = 1.d0
        if (reaction%update_mnrl_surf_with_porosity) then
          porosity_scale = vec_p(local_id)** &
             reaction%mineral%kinmnrl_surf_area_porosity_pwr(imnrl)
!       geh: srf_area_vol_frac_pwr must be defined on a per mineral basis, not
!       solely material type.
!       material_property_array(patch%imat(ghosted_id))%ptr%mnrl_surf_area_porosity_pwr
        endif

        volfrac_scale = 1.d0
        if (rt_auxvars(ghosted_id)%mnrl_volfrac0(imnrl) > 0.d0) then
          volfrac_scale = (rt_auxvars(ghosted_id)%mnrl_volfrac(imnrl)/ &
                         rt_auxvars(ghosted_id)%mnrl_volfrac0(imnrl))** &
             reaction%mineral%kinmnrl_surf_area_vol_frac_pwr(imnrl)
!       geh: srf_area_vol_frac_pwr must be defined on a per mineral basis, not
!       solely material type.
!       material_property_array(patch%imat(ghosted_id))%ptr%mnrl_surf_area_volfrac_pwr
!         rt_auxvars(ghosted_id)%mnrl_area(imnrl) = &
!           rt_auxvars(ghosted_id)%mnrl_area0(imnrl)*porosity_scale*volfrac_scale
!       else
!         rt_auxvars(ghosted_id)%mnrl_area(imnrl) = &
!           rt_auxvars(ghosted_id)%mnrl_area0(imnrl)
        endif

        rt_auxvars(ghosted_id)%mnrl_area(imnrl) = &
            rt_auxvars(ghosted_id)%mnrl_area0(imnrl)*porosity_scale*volfrac_scale

        if (reaction%update_armor_mineral_surface .and. &
            reaction%mineral%kinmnrl_armor_crit_vol_frac(imnrl) > 0.d0) then
          imnrl_armor = imnrl
          do imnrl1 = 1, reaction%mineral%nkinmnrl
            if (reaction%mineral%kinmnrl_armor_min_names(imnrl) == &
                reaction%mineral%kinmnrl_names(imnrl1)) then
              imnrl_armor = imnrl1
              exit
            endif
          enddo

!         print *,'update-armor: ',imnrl,imnrl_armor, &
!         reaction%mineral%kinmnrl_armor_min_names(imnrl_armor)

!       check for negative surface area armoring correction
          if (reaction%mineral%kinmnrl_armor_crit_vol_frac(imnrl) > &
              rt_auxvars(ghosted_id)%mnrl_volfrac(imnrl_armor)) then

            if (reaction%update_armor_mineral_surface_flag == 0) then ! surface unarmored
              rt_auxvars(ghosted_id)%mnrl_area(imnrl) = &
                rt_auxvars(ghosted_id)%mnrl_area(imnrl) * &
                ((reaction%mineral%kinmnrl_armor_crit_vol_frac(imnrl) &
                - rt_auxvars(ghosted_id)%mnrl_volfrac(imnrl_armor))/ &
                reaction%mineral%kinmnrl_armor_crit_vol_frac(imnrl))** &
                reaction%mineral%kinmnrl_surf_area_vol_frac_pwr(imnrl)
            else
              rt_auxvars(ghosted_id)%mnrl_area(imnrl) = rt_auxvars(ghosted_id)%mnrl_area0(imnrl)
              reaction%update_armor_mineral_surface_flag = 0
            endif
          else
            rt_auxvars(ghosted_id)%mnrl_area(imnrl) = 0.d0
            reaction%update_armor_mineral_surface_flag = 1 ! surface armored
          endif
        endif

!       print *,'update min srf: ',imnrl,local_id,reaction%mineral%kinmnrl_names(imnrl), &
!       reaction%mineral%kinmnrl_armor_min_names(imnrl), &
!       reaction%update_armor_mineral_surface, &
!       rt_auxvars(ghosted_id)%mnrl_area(imnrl), &
!       reaction%mineral%kinmnrl_armor_pwr(imnrl), &
!       reaction%mineral%kinmnrl_armor_crit_vol_frac(imnrl), &
!       rt_auxvars(ghosted_id)%mnrl_volfrac(imnrl_armor), &
!       rt_auxvars(ghosted_id)%mnrl_volfrac(imnrl)
      enddo
    enddo

    if (reaction%update_mnrl_surf_with_porosity) then
      call VecRestoreArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
    endif
!geh:remove
    call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 TORTUOSITY,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 TORTUOSITY,ZERO_INTEGER)
  endif
      
  if (reaction%update_tortuosity) then
    call VecGetArrayF90(field%tortuosity0,tortuosity0_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      scale = vec_p(local_id)** &
        material_property_array(patch%imat(ghosted_id))%ptr%tortuosity_pwr
      material_auxvars(ghosted_id)%tortuosity = &
        tortuosity0_p(local_id)*scale
    enddo
    call VecRestoreArrayF90(field%tortuosity0,tortuosity0_p, &
                            ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
    call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 TORTUOSITY,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 TORTUOSITY,ZERO_INTEGER)
  endif
      
  if (reaction%update_permeability) then
    call VecGetArrayF90(field%perm0_xx,perm0_xx_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%perm0_zz,perm0_zz_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%perm0_yy,perm0_yy_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%porosity0,porosity0_p,ierr);CHKERRQ(ierr)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      imat = patch%imat(ghosted_id)
      critical_porosity = material_property_array(imat)%ptr% &
                            permeability_crit_por
      porosity_base_ = material_auxvars(ghosted_id)%porosity_base
      scale = 0.d0
      if (porosity_base_ > critical_porosity .and. &
          porosity0_p(local_id) > critical_porosity) then
        scale = ((porosity_base_ - critical_porosity) / &
                 (porosity0_p(local_id) - critical_porosity)) ** &
                material_property_array(imat)%ptr%permeability_pwr
      endif
      scale = max(material_property_array(imat)%ptr% &
                    permeability_min_scale_fac,scale)
      !geh: this is a kludge for gfortran.  the code reports errors when 
      !     material_auxvars(ghosted_id)%permeability is used.
      ! This is not an issue with Intel
#if 1
      perm_ptr => material_auxvars(ghosted_id)%permeability
      perm_ptr(perm_xx_index) = perm0_xx_p(local_id)*scale
      perm_ptr(perm_yy_index) = perm0_yy_p(local_id)*scale
      perm_ptr(perm_zz_index) = perm0_zz_p(local_id)*scale
#else
      material_auxvars(ghosted_id)%permeability(perm_xx_index) = &
        perm0_xx_p(local_id)*scale
      material_auxvars(ghosted_id)%permeability(perm_yy_index) = &
        perm0_yy_p(local_id)*scale
      material_auxvars(ghosted_id)%permeability(perm_zz_index) = &
        perm0_zz_p(local_id)*scale
#endif
    enddo
    call VecRestoreArrayF90(field%perm0_xx,perm0_xx_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%perm0_zz,perm0_zz_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%perm0_yy,perm0_yy_p,ierr);CHKERRQ(ierr)

    call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_X,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_X,ZERO_INTEGER)
    call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_Y,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_Y,ZERO_INTEGER)
    call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_Z,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 PERMEABILITY_Z,ZERO_INTEGER)
  endif  
  
  ! perform check to ensure that porosity is bounded between 0 and 1
  ! since it is calculated as 1.d-sum_volfrac, it cannot be > 1
  call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               POROSITY,POROSITY_BASE)
  call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                  field%work,ONEDOF)
  call VecMin(field%work,ivalue,min_value,ierr);CHKERRQ(ierr)
  if (min_value < 0.d0) then
    write(option%io_buffer,*) 'Sum of mineral volume fractions has ' // &
      'exceeded 1.d0 at cell (note PETSc numbering): ', ivalue
    call PrintErrMsg(option)
  endif
   
end subroutine RealizationUpdatePropertiesTS

! ************************************************************************** !

subroutine RealizationUpdatePropertiesNI(realization)
  ! 
  ! Updates coupled properties at each grid cell
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/05/09
  ! 

  use Discretization_module
  use Field_module
  use Grid_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Material_Aux_class
  use Variables_module, only : POROSITY, TORTUOSITY, PERMEABILITY_X, &
                               PERMEABILITY_Y, PERMEABILITY_Z
 
  implicit none
  
  class(realization_subsurface_type) :: realization

#if 0
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  class(reaction_rt_type), pointer :: reaction
  type(grid_type), pointer :: grid
  type(material_property_ptr_type), pointer :: material_property_array(:)
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:) 
  type(discretization_type), pointer :: discretization
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: local_id, ghosted_id
  PetscInt :: imnrl, imnrl1, imnrl_armor, imat
  PetscReal :: sum_volfrac
  PetscReal :: scale, porosity_scale, volfrac_scale
  PetscBool :: porosity_updated
  PetscReal, pointer :: vec_p(:)
  PetscReal, pointer :: porosity0_p(:)
  PetscReal, pointer :: porosity_mnrl_loc_p(:)
  PetscReal, pointer :: tortuosity0_p(:)
  PetscReal, pointer :: perm0_xx_p(:), perm0_yy_p(:), perm0_zz_p(:)
  PetscReal :: min_value  
  PetscInt :: ivalue
  PetscErrorCode :: ierr

  option => realization%option
  discretization => realization%discretization
  patch => realization%patch
  field => realization%field
  reaction => realization%reaction
  grid => patch%grid
  material_property_array => patch%material_property_array
  rt_auxvars => patch%aux%RT%auxvars
  material_auxvars => patch%aux%Material%auxvars
#endif

end subroutine RealizationUpdatePropertiesNI

! ************************************************************************** !

subroutine RealizationCalcMineralPorosity(realization)
  ! 
  ! Calculates porosity based on the sum of mineral volume fractions
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/03/14
  !

  use Discretization_module
  use Field_module
  use Grid_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Material_Aux_class
  use Variables_module, only : POROSITY
 
  implicit none
  
  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  class(reaction_rt_type), pointer :: reaction
  type(grid_type), pointer :: grid
  type(reactive_transport_auxvar_type), pointer :: rt_auxvars(:) 
  type(discretization_type), pointer :: discretization
  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: local_id, ghosted_id
  PetscInt :: imnrl 
  PetscReal :: sum_volfrac
  PetscErrorCode :: ierr

  option => realization%option
  discretization => realization%discretization
  patch => realization%patch
  field => realization%field
  reaction => realization%reaction
  grid => patch%grid
  rt_auxvars => patch%aux%RT%auxvars
  material_auxvars => patch%aux%Material%auxvars

  if (reaction%mineral%nkinmnrl > 0) then
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      ! Go ahead and compute for inactive cells since their porosity does
      ! not matter (avoid check on active/inactive)
      sum_volfrac = 0.d0
      do imnrl = 1, reaction%mineral%nkinmnrl
        sum_volfrac = sum_volfrac + &
                      rt_auxvars(ghosted_id)%mnrl_volfrac(imnrl)
      enddo 
      ! the adjusted porosity becomes:
      ! 1 - sum(mineral volume fractions), but is truncated.
      material_auxvars(ghosted_id)%porosity_base = &
        max(1.d0-sum_volfrac,reaction%minimum_porosity)
    enddo
  endif
  ! update ghosted porosities
  call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               POROSITY,POROSITY_BASE)
  call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)
  call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               POROSITY,POROSITY_BASE)
  call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                               POROSITY,POROSITY_CURRENT)
!  call MaterialSetAuxVarScalar(patch%aux%Material,UNINITIALIZED_DOUBLE, &
!                               POROSITY,POROSITY_CURRENT)

end subroutine RealizationCalcMineralPorosity

! ************************************************************************** !

subroutine RealLocalToLocalWithArray(realization,array_id)
  ! 
  ! Takes an F90 array that is ghosted
  ! and updates the ghosted values
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/09/11
  ! 

  use Discretization_module
  use Field_module
  use Grid_module

  implicit none

  class(realization_subsurface_type) :: realization
  PetscInt :: array_id
  
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field

  field => realization%field
  patch => realization%patch

  grid => patch%grid
  select case(array_id)
    case(MATERIAL_ID_ARRAY)
      call GridCopyIntegerArrayToVec(grid,patch%imat,field%work_loc, &
                                     grid%ngmax)
    case(SATURATION_FUNCTION_ID_ARRAY)
      call GridCopyIntegerArrayToVec(grid,patch%sat_func_id, &
                                     field%work_loc, grid%ngmax)
  end select

  call DiscretizationLocalToLocal(realization%discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)

  select case(array_id)
    case(MATERIAL_ID_ARRAY)
      call GridCopyVecToIntegerArray(grid,patch%imat,field%work_loc, &
                                      grid%ngmax)
    case(SATURATION_FUNCTION_ID_ARRAY)
      call GridCopyVecToIntegerArray(grid,patch%sat_func_id, &
                                      field%work_loc, grid%ngmax)
  end select

end subroutine RealLocalToLocalWithArray

! ************************************************************************** !

subroutine RealizationCountCells(realization,global_total_count, &
                                 global_active_count,total_count,active_count)
  ! 
  ! Counts # of active and inactive grid cells
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/01/10
  ! 

  use Option_module

  implicit none
  
  class(realization_subsurface_type) :: realization
  PetscInt :: global_total_count
  PetscInt :: global_active_count
  PetscInt :: total_count
  PetscInt :: active_count
  
  PetscInt :: patch_total_count
  PetscInt :: patch_active_count
  PetscInt :: temp_int_in(2), temp_int_out(2)
  PetscErrorCode :: ierr
  
  type(patch_type), pointer :: patch
  
  total_count = 0
  active_count = 0
    
  patch => realization%patch
  call PatchCountCells(patch,patch_total_count,patch_active_count)
  total_count = total_count + patch_total_count
  active_count = active_count + patch_active_count
  
  temp_int_in(1) = total_count
  temp_int_in(2) = active_count
  call MPI_Allreduce(temp_int_in,temp_int_out,TWO_INTEGER_MPI,MPIU_INTEGER, &
                     MPI_SUM,realization%option%mycomm,ierr)
  global_total_count = temp_int_out(1)
  global_active_count = temp_int_out(2)

end subroutine RealizationCountCells

! ************************************************************************** !

subroutine RealizationPrintGridStatistics(realization)
  ! 
  ! Prints statistics regarding the numerical
  ! discretization
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/01/10
  ! 

  use Grid_module

  implicit none

  class(realization_subsurface_type) :: realization
  
  type(option_type), pointer :: option
  type(grid_type), pointer :: grid

  PetscInt :: i1, i2, i3
  PetscReal :: r1, r2, r3
  PetscInt :: global_total_count, global_active_count
  PetscInt :: total_count, active_count
  PetscReal :: total_min, total_max, total_mean, total_variance
  PetscReal :: active_min, active_max, active_mean, active_variance
  PetscInt :: inactive_histogram(12), temp_int_out(12)
  PetscReal :: inactive_percentages(12)
  PetscErrorCode :: ierr

  option => realization%option
  grid => realization%patch%grid

  ! print # of active and inactive grid cells
  call RealizationCountCells(realization,global_total_count, &
                             global_active_count,total_count,active_count)
  r1 = dble(total_count)
  call OptionMaxMinMeanVariance(r1,total_max, &
                                total_min,total_mean, &
                                total_variance,PETSC_TRUE,option)
  r1 = dble(active_count)
  call OptionMaxMinMeanVariance(r1,active_max, &
                                active_min,active_mean, &
                                active_variance,PETSC_TRUE,option)
                  
  r1 = dble(active_count) / dble(total_count)    
  inactive_histogram = 0                          
  if (r1 >= (1.d0-1.d-8)) then
    inactive_histogram(12) = 1
  else if (r1 >= .9d0 .and. r1 < (1.d0-1.d-8)) then
    inactive_histogram(11) = 1
  else if (r1 >= .8d0 .and. r1 < .9d0) then
    inactive_histogram(10) = 1
  else if (r1 >= .7d0 .and. r1 < .8d0) then
    inactive_histogram(9) = 1
  else if (r1 >= .6d0 .and. r1 < .7d0) then
    inactive_histogram(8) = 1
  else if (r1 >= .5d0 .and. r1 < .6d0) then
    inactive_histogram(7) = 1
  else if (r1 >= .4d0 .and. r1 < .5d0) then
    inactive_histogram(6) = 1
  else if (r1 >= .3d0 .and. r1 < .4d0) then
    inactive_histogram(5) = 1
  else if (r1 >= .2d0 .and. r1 < .3d0) then
    inactive_histogram(4) = 1
  else if (r1 >= .1d0 .and. r1 < .2d0) then
    inactive_histogram(3) = 1
  else if (r1 > 1.d-20 .and. r1 < .1d0) then
    inactive_histogram(2) = 1
  else if (r1 < 1.d-20) then
    inactive_histogram(1) = 1
  endif
  
  call MPI_Allreduce(inactive_histogram,temp_int_out,TWELVE_INTEGER_MPI, &
                     MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)

  ! why I cannot use *100, I do not know....geh
  inactive_percentages = dble(temp_int_out)/dble(option%mycommsize)*10.d0
  inactive_percentages = inactive_percentages+1.d-8

  r1 = 0.d0
  do i1 = 1, 12
    r1 = r1 + inactive_percentages(i1)
  enddo
                                
  i1 = UNINITIALIZED_INTEGER
  i2 = UNINITIALIZED_INTEGER
  i3 = UNINITIALIZED_INTEGER
  if (associated(grid%structured_grid)) then
    i1 = grid%structured_grid%npx_final
    i2 = grid%structured_grid%npy_final
    i3 = grid%structured_grid%npz_final
  endif
  if (OptionPrintToScreen(option)) then
    write(*,'(/," Grid Stats:",/, &
              & "                       Global # cells: ",i12,/, &
              & "                Global # active cells: ",i12,/, &
              & "                              # cores: ",i12,/, &
              & "         Processor core decomposition: ",3i6,/, &
              & "               Maximum # cells / core: ",i12,/, &
              & "               Minimum # cells / core: ",i12,/, &
              & "               Average # cells / core: ",1pe12.4,/, &
              & "               Std Dev # cells / core: ",1pe12.4,/, &
              & "        Maximum # active cells / core: ",i12,/, &
              & "        Minimum # active cells / core: ",i12,/, &
              & "        Average # active cells / core: ",1pe12.4,/, &
              & "        Std Dev # active cells / core: ",1pe12.4,/,/, &
              & "        % cores with % active cells =       0%: ",1f7.2,/, &
              & "        % cores with % active cells =  0.1-10%: ",1f7.2,/, &
              & "        % cores with % active cells =   10-20%: ",1f7.2,/, &
              & "        % cores with % active cells =   20-30%: ",1f7.2,/, &
              & "        % cores with % active cells =   30-40%: ",1f7.2,/, &
              & "        % cores with % active cells =   40-50%: ",1f7.2,/, &
              & "        % cores with % active cells =   50-60%: ",1f7.2,/, &
              & "        % cores with % active cells =   60-70%: ",1f7.2,/, &
              & "        % cores with % active cells =   70-80%: ",1f7.2,/, &
              & "        % cores with % active cells =   80-90%: ",1f7.2,/, &
              & "        % cores with % active cells = 90-99.9%: ",1f7.2,/, &
              & "        % cores with % active cells =     100%: ",1f7.2,/, &
              & "                                        Check : ",1f7.2,/)') &
           global_total_count, &
           global_active_count, &
           option%mycommsize, &
           i1,i2,i3, &
           int(total_max+1.d-4), &
           int(total_min+1.d-4), &
           total_mean, sqrt(total_variance), &
           int(active_max+1.d-4), &
           int(active_min+1.d-4), &
           active_mean, sqrt(active_variance), &
           inactive_percentages(1), &
           inactive_percentages(2), &
           inactive_percentages(3), &
           inactive_percentages(4), &
           inactive_percentages(5), &
           inactive_percentages(6), &
           inactive_percentages(7), &
           inactive_percentages(8), &
           inactive_percentages(9), &
           inactive_percentages(10), &
           inactive_percentages(11), &
           inactive_percentages(12), &
           r1
  endif
  if (OptionPrintToFile(option)) then
    write(option%fid_out,'(/," Grid Stats:",/, &
               & "                       Global # cells: ",i12,/, &
               & "                Global # active cells: ",i12,/, &
               & "                              # cores: ",i12,/, &
               & "         Processor core decomposition: ",3i6,/, &
               & "               Maximum # cells / core: ",i12,/, &
               & "               Minimum # cells / core: ",i12,/, &
               & "               Average # cells / core: ",1pe12.4,/, &
               & "               Std Dev # cells / core: ",1pe12.4,/, &
               & "        Maximum # active cells / core: ",i12,/, &
               & "        Minimum # active cells / core: ",i12,/, &
               & "        Average # active cells / core: ",1pe12.4,/, &
               & "        Std Dev # active cells / core: ",1pe12.4,/,/, &
               & "        % cores with % active cells =       0%: ",1f7.2,/, &
               & "        % cores with % active cells =  0.1-10%: ",1f7.2,/, &
               & "        % cores with % active cells =   10-20%: ",1f7.2,/, &
               & "        % cores with % active cells =   20-30%: ",1f7.2,/, &
               & "        % cores with % active cells =   30-40%: ",1f7.2,/, &
               & "        % cores with % active cells =   40-50%: ",1f7.2,/, &
               & "        % cores with % active cells =   50-60%: ",1f7.2,/, &
               & "        % cores with % active cells =   60-70%: ",1f7.2,/, &
               & "        % cores with % active cells =   70-80%: ",1f7.2,/, &
               & "        % cores with % active cells =   80-90%: ",1f7.2,/, &
               & "        % cores with % active cells = 90-99.9%: ",1f7.2,/, &
               & "        % cores with % active cells =     100%: ",1f7.2,/, &
               & "                                        Check : ",1f7.2,/)') &
           global_total_count, &
           global_active_count, &
           option%mycommsize, &
           i1,i2,i3, &
           int(total_max+1.d-4), &
           int(total_min+1.d-4), &
           total_mean, sqrt(total_variance), &
           int(active_max+1.d-4), &
           int(active_min+1.d-4), &
           active_mean, sqrt(active_variance), &
           inactive_percentages(1), &
           inactive_percentages(2), &
           inactive_percentages(3), &
           inactive_percentages(4), &
           inactive_percentages(5), &
           inactive_percentages(6), &
           inactive_percentages(7), &
           inactive_percentages(8), &
           inactive_percentages(9), &
           inactive_percentages(10), &
           inactive_percentages(11), &
           inactive_percentages(12), &
           r1
  endif

end subroutine RealizationPrintGridStatistics

! ************************************************************************** !

subroutine RealizationCalculateCFL1Timestep(realization,max_dt_cfl_1)
  ! 
  ! Calculates largest time step that
  ! preserves a CFL # of 1 in a realization
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/07/11
  ! 

  implicit none

  class(realization_subsurface_type) realization
  PetscReal :: max_dt_cfl_1
  
  type(patch_type), pointer :: patch
  PetscReal :: max_dt_cfl_1_patch
  PetscErrorCode :: ierr
  
  max_dt_cfl_1 = 1.d20
  patch => realization%patch
  call PatchCalculateCFL1Timestep(patch,realization%option, &
                                  max_dt_cfl_1_patch)
  max_dt_cfl_1 = min(max_dt_cfl_1,max_dt_cfl_1_patch)

  ! get the minimum across all cores
  call MPI_Allreduce(MPI_IN_PLACE,max_dt_cfl_1,ONE_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_MIN, &
                     realization%option%mycomm,ierr)

end subroutine RealizationCalculateCFL1Timestep

! ************************************************************************** !

subroutine RealizUnInitializedVarsFlow(realization)
  ! 
  ! Checks for uninitialized flow variables
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/06/16
  ! 
  use Option_module
  use Material_Aux_class
  use Variables_module, only : VOLUME, BASE_POROSITY, PERMEABILITY_X, &
                               PERMEABILITY_Y, PERMEABILITY_Z

  implicit none
  
  class(realization_subsurface_type) :: realization

  character(len=MAXWORDLENGTH) :: var_name
  PetscInt :: i

  call RealizUnInitializedVar1(realization,VOLUME,'volume')
  ! mineral porosity is the base, unmodified porosity
  call RealizUnInitializedVar1(realization,BASE_POROSITY,'porosity')
  call RealizUnInitializedVar1(realization,PERMEABILITY_X,'permeability X')
  call RealizUnInitializedVar1(realization,PERMEABILITY_Y,'permeability Y')
  call RealizUnInitializedVar1(realization,PERMEABILITY_Z,'permeability Z')
  do i = 1, max_material_index
    var_name = MaterialAuxIndexToPropertyName(i)
    call RealizUnInitializedVar1(realization,i,var_name)
  enddo

end subroutine RealizUnInitializedVarsFlow

! ************************************************************************** !

subroutine RealizUnInitializedVarsTran(realization)
  ! 
  ! Checks for uninitialized transport variables
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/06/16
  ! 

  use Grid_module
  use Patch_module
  use Option_module
  use Material_module
  use Material_Aux_class
  use Variables_module, only : VOLUME, BASE_POROSITY, TORTUOSITY

  implicit none
  
  class(realization_subsurface_type) :: realization

  call RealizUnInitializedVar1(realization,VOLUME,'volume')
  ! mineral porosity is the base, unmodified porosity
  call RealizUnInitializedVar1(realization,BASE_POROSITY,'porosity')
  call RealizUnInitializedVar1(realization,TORTUOSITY,'tortuosity')

end subroutine RealizUnInitializedVarsTran

! ************************************************************************** !

subroutine RealizUnInitializedVar1(realization,ivar,var_name)
  ! 
  ! Checks whether a variable is initialized at all active grid cells
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/06/16
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  use Field_module
  use Patch_module
  use Grid_module

  implicit none
  
  class(realization_subsurface_type) :: realization
  PetscInt :: ivar
  character(len=*) :: var_name
  
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  PetscReal, pointer :: vec_p(:)
  PetscInt :: local_id
  PetscInt :: imin
  PetscReal :: rmin
  character(len=MAXWORDLENGTH) :: word
  PetscErrorCode :: ierr

  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  call RealizationGetVariable(realization,realization%field%work, &
                              ivar,ZERO_INTEGER)
  ! apply mask to filter inactive cells
  call VecGetArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax
    ! if inactive, set to 1.d-40
    if (patch%imat(grid%nL2G(local_id)) <= 0) vec_p(local_id) = 1.d-40
  enddo
  call VecRestoreArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
  call VecMin(field%work,imin,rmin,ierr);CHKERRQ(ierr)
  if (Uninitialized(rmin)) then
    write(word,*) imin+1 ! zero to one based indexing
    option%io_buffer = 'Incorrect assignment of variable (' &
      // trim(var_name) // ',cell=' // trim(adjustl(word)) // ').'
    call PrintErrMsgToDev(option,'send your input deck.')
  endif

end subroutine RealizUnInitializedVar1

! ************************************************************************** !

subroutine RealizationLimitDTByCFL(realization,cfl_governor,dt,dt_max)
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/09/16 
  !
  use Option_module
  use Output_Aux_module

  implicit none

  class(realization_subsurface_type) :: realization
  PetscReal :: cfl_governor
  PetscReal :: dt
  PetscReal :: dt_max

  PetscReal :: max_dt_cfl_1
  PetscReal :: prev_dt
  type(output_option_type), pointer :: output_option

  if (Initialized(cfl_governor)) then
    call RealizationCalculateCFL1Timestep(realization,max_dt_cfl_1)
    if (dt/cfl_governor > max_dt_cfl_1) then
      prev_dt = dt
      dt = max_dt_cfl_1*cfl_governor
      ! have to set dt_max here so that timestepper%dt_max is truncated
      ! for timestepper_base%revert_dt in TimestepperBaseSetTargetTime
      dt_max = dt
      output_option => realization%output_option
      if (OptionPrintToScreen(realization%option)) then
        write(*, &
          '(" CFL Limiting (",f4.1,"): ",1pe12.4," -> ",1pe12.4," [",a,"]")') &
              cfl_governor,prev_dt/output_option%tconv, &
              dt/output_option%tconv,trim(output_option%tunit)
      endif
      if (OptionPrintToFile(realization%option)) then
        write(realization%option%fid_out, &
          '(" CFL Limiting (",f4.1,"): ",1pe12.4," -> ",1pe12.4," [",a,"]")') &
              cfl_governor,prev_dt/output_option%tconv, &
              dt/output_option%tconv,trim(output_option%tunit)
      endif
    endif
  endif

end subroutine RealizationLimitDTByCFL

! ************************************************************************** !

subroutine RealizationStrip(this)
  ! 
  ! Deallocates a realization
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  use Dataset_module
  use Output_Eclipse_module, only : ReleaseEwriterBuffers

  implicit none
  
  class(realization_subsurface_type) :: this
  
  call RealizationBaseStrip(this)
  call RegionDestroyList(this%region_list)
  
  call FlowConditionDestroyList(this%flow_conditions)

  !  Destroy the list of wells held by well_data
  call WellDataDestroyList(this%well_data, this%option)
  !  Release output buffers held by Output_Eclipse_module
  call ReleaseEwriterBuffers()

  call TranConditionDestroyList(this%transport_conditions)
  call TranConstraintListDestroy(this%transport_constraints)

  if (associated(this%fluid_property_array)) &
    deallocate(this%fluid_property_array)
  nullify(this%fluid_property_array)
  call FluidPropertyDestroy(this%fluid_properties)
  
  call MaterialPropertyDestroy(this%material_properties)

  call SaturationFunctionDestroy(this%saturation_functions)
  call CharacteristicCurvesDestroy(this%characteristic_curves)  

  if (this%option%iflowmode == G_MODE) then
    if (associated(this%thermal_characteristic_curves)) then
      call ThermalCharacteristicCurvesDestroy(this%thermal_characteristic_curves)
    end if
  end if
  
  call DatasetDestroy(this%datasets)
  
  call DatasetDestroy(this%uniform_velocity_dataset)

  ! nullify since they are pointers to reaction_base in realization_base
  nullify(this%reaction)
  nullify(this%reaction_nw)
  
  call TranConstraintDestroy(this%sec_transport_constraint)
  
end subroutine RealizationStrip

end module Realization_Subsurface_class
