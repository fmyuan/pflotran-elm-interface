module Realization_Subsurface_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Realization_Base_class
  use Option_module
  use Input_Aux_module
  use Region_module
  use Condition_module
  use Transport_Constraint_Base_module
  use Transport_Constraint_module
  use Material_module
  use Saturation_Function_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Thermal_module
  use Material_Transform_module
  use Dataset_Base_class
  use Fluid_module
  use Patch_module
  use Reaction_Aux_module
  use NW_Transport_Aux_module
  use Survey_module


  implicit none

private

  type, public, extends(realization_base_type) :: realization_subsurface_type

    type(region_list_type), pointer :: region_list
    type(condition_list_type), pointer :: flow_conditions
    type(tran_condition_list_type), pointer :: transport_conditions
    type(tran_constraint_list_type), pointer :: transport_constraints
    type(geop_condition_list_type), pointer :: geophysics_conditions

    type(material_property_type), pointer :: material_properties
    type(fluid_property_type), pointer :: fluid_properties
    type(fluid_property_type), pointer :: fluid_property_array(:)
    type(saturation_function_type), pointer :: saturation_functions
    class(characteristic_curves_type), pointer :: characteristic_curves
    class(cc_thermal_type), pointer :: characteristic_curves_thermal
    type(material_transform_type), pointer :: material_transform
    class(dataset_base_type), pointer :: datasets

    class(dataset_base_type), pointer :: uniform_velocity_dataset
    character(len=MAXSTRINGLENGTH) :: nonuniform_velocity_filename

    class(reaction_rt_type), pointer :: reaction
    class(reaction_nw_type), pointer :: reaction_nw

    type(survey_type), pointer :: survey

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
            RealizationLimitDTByCFL, &
            RealizationReadGeopSurveyFile, &
            RealizationCheckConsistency, &
            RealizationPrintStateAtCells

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
  allocate(realization%transport_conditions)
  call TranConditionInitList(realization%transport_conditions)
  allocate(realization%transport_constraints)
  call TranConstraintInitList(realization%transport_constraints)
  allocate(realization%geophysics_conditions)
  call GeopConditionInitList(realization%geophysics_conditions)

  nullify(realization%material_properties)
  nullify(realization%fluid_properties)
  nullify(realization%fluid_property_array)
  nullify(realization%saturation_functions)
  nullify(realization%characteristic_curves)
  nullify(realization%characteristic_curves_thermal)
  nullify(realization%material_transform)
  nullify(realization%datasets)
  nullify(realization%uniform_velocity_dataset)
  nullify(realization%reaction)
  nullify(realization%reaction_nw)
  nullify(realization%survey)
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
  use Communicator_Structured_class, only : StructuredCommunicatorCreate
  use Communicator_Unstructured_class, only : UnstructuredCommunicatorCreate

  implicit none

  class(realization_subsurface_type) :: realization

  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(option_type), pointer :: option
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
    if (option%flow%full_perm_tensor) then
      call DiscretizationDuplicateVector(discretization,field%work, &
                                         field%perm0_xy)
      call DiscretizationDuplicateVector(discretization,field%work, &
                                         field%perm0_xz)
      call DiscretizationDuplicateVector(discretization,field%work, &
                                         field%perm0_yz)
    endif

    ! 1-dof local
    !call DiscretizationDuplicateVector(discretization,field%work_loc, &
    !                                   field%xyz)

    ! ndof degrees of freedom, global
    call DiscretizationCreateVector(discretization,NFLOWDOF,field%flow_xx, &
                                    GLOBAL,option)
    ! for scaling pressure in the range of saturation
    if (option%flow%scale_all_pressure) then
    call DiscretizationDuplicateVector(discretization,field%flow_xx, &
                                       field%flow_scaled_xx)
    call DiscretizationDuplicateVector(discretization,field%flow_xx, &
                                       field%flow_work_loc)
    endif
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
    call DiscretizationCreateVector(discretization,NFLOWDOF, &
                                    field%flow_xx_loc,LOCAL,option)

    if ((option%iflowmode == RICHARDS_TS_MODE) .or. &
        (option%iflowmode == TH_TS_MODE)) then
      call DiscretizationCreateVector(discretization,NFLOWDOF, &
                                      field%flow_xxdot, &
                                      GLOBAL,option)

      call DiscretizationCreateVector(discretization,NFLOWDOF, &
                                      field%flow_xxdot_loc, &
                                      LOCAL,option)
    endif

    if (option%iflowmode == PNF_MODE) then
      call DiscretizationDuplicateVector(discretization,field%flow_xx, &
                                         field%flow_rhs)
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
      ! create the ntran dof vector for storage of the solution
      call DiscretizationCreateVector(discretization,NTRANDOF,field%tran_xx, &
                                      GLOBAL,option)
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_yy)
      call DiscretizationDuplicateVector(discretization,field%tran_xx, &
                                         field%tran_dxx)

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
                      (option%nflowdof*MAX_FACE_PER_CELL+1)*realization% &
                        patch%grid%nlmax, &
                      PETSC_DETERMINE,field%flowrate_inst,ierr);CHKERRQ(ierr)
    call VecSet(field%flowrate_inst,0.d0,ierr);CHKERRQ(ierr)
  endif

  ! Allocate vectors to hold velocity at face
  if (realization%output_option%print_hdf5_vel_face) then

    ! vx
    call VecCreateMPI(option%mycomm, &
                      (option%nflowspec*MAX_FACE_PER_CELL+1)*realization% &
                        patch%grid%nlmax, &
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
                      size(grid%unstructured_grid%explicit_grid% &
                        connections,2), &
                      PETSC_DETERMINE,field%flowrate_inst,ierr);CHKERRQ(ierr)
    call VecSet(field%flowrate_inst,0.d0,ierr);CHKERRQ(ierr)
  endif

  ! If average flowrate has to be saved, create a vector for it
  if (realization%output_option%print_hdf5_aveg_mass_flowrate.or. &
      realization%output_option%print_hdf5_aveg_energy_flowrate) then
    call VecCreateMPI(option%mycomm, &
                      (option%nflowdof*MAX_FACE_PER_CELL+1)*realization% &
                        patch%grid%nlmax, &
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

  call PatchProcessCouplers(realization%patch,realization%flow_conditions, &
                            realization%transport_conditions, &
                            realization%geophysics_conditions,realization%option)

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
  use Characteristic_Curves_Thermal_module
  use TH_Aux_module, only : th_ice_model


  implicit none

  class(realization_subsurface_type) :: realization

  PetscInt :: i, num_mat_prop
  type(option_type), pointer :: option
  type(material_property_type), pointer :: cur_material_property
  PetscReal, allocatable :: check_thermal_conductivity(:,:)
  type(patch_type), pointer :: patch
  character(len=MAXSTRINGLENGTH) :: string, verify_string, mat_string
  class(dataset_base_type), pointer :: dataset
  class(cc_thermal_type), pointer :: thermal_cc

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

  ! set up analogous mapping to thermal characteristic curves, if used
  num_mat_prop = size(patch%material_property_array)
  allocate(check_thermal_conductivity(2,num_mat_prop))
  do i = 1, num_mat_prop
    if (associated(patch%material_property_array(i)%ptr)) then
      check_thermal_conductivity(1,i) = &
        patch%material_property_array(i)%ptr%thermal_conductivity_dry
      check_thermal_conductivity(2,i) = &
        patch%material_property_array(i)%ptr%thermal_conductivity_wet
    endif
  enddo
  if (associated(realization%characteristic_curves_thermal)) then
    if (maxval(check_thermal_conductivity(:,:)) >= 0.d0) then
      option%io_buffer = 'Cannot combine material-based thermal conductivity'//&
                         ' input format with thermal characteristic curves. '//&
                         'Use TCC with "DEFAULT" specification instead.'
      call PrintErrMsg(option)
    endif
    patch%characteristic_curves_thermal => &
         realization%characteristic_curves_thermal
    call CharCurvesThermalConvertListToArray( &
         patch%characteristic_curves_thermal, &
         patch%char_curves_thermal_array, option)
    do i = 1, size(patch%char_curves_thermal_array)
      select type(tcf => patch%char_curves_thermal_array(i)%ptr% &
                  thermal_conductivity_function)
      class is(kT_composite_type)
        call CompositeTCCList(patch%characteristic_curves_thermal, &
                              tcf,option)
      end select
    enddo
  else if (maxval(check_thermal_conductivity(:,:)) >= 0.d0) then
    ! use default tcc curve for legacy thermal conductivity input by material
    do i = 1, num_mat_prop
      if (.not. option%iflowmode == G_MODE) then
        ! some modes outside of general will only use one thermal conductivity
        ! if that is the case, use default values as fallback options
        if (patch%material_property_array(i)%ptr%thermal_conductivity_wet == &
          UNINITIALIZED_DOUBLE) then
          patch%material_property_array(i)%ptr%thermal_conductivity_wet = 2.d0
        endif
        if (patch%material_property_array(i)%ptr%thermal_conductivity_dry == &
          UNINITIALIZED_DOUBLE) then
          patch%material_property_array(i)%ptr%thermal_conductivity_dry = 5.d-1
        endif
      endif
      thermal_cc => CharCurvesThermalCreate()
      thermal_cc%name = patch%material_property_array(i)%ptr% &
                                thermal_conductivity_func_name
      if (option%iflowmode == TH_MODE .or. option%iflowmode == TH_TS_MODE) then
        thermal_cc%thermal_conductivity_function => TCFFrozenCreate()
        call TCFAssignFrozen(thermal_cc%thermal_conductivity_function,&
          patch%material_property_array(i)%ptr%thermal_conductivity_wet, &
          patch%material_property_array(i)%ptr%thermal_conductivity_dry, &
          patch%material_property_array(i)%ptr%thermal_conductivity_frozen, &
          patch%material_property_array(i)%ptr%alpha, &
          patch%material_property_array(i)%ptr%alpha_fr, &
          th_ice_model, &
          option)
      else
        thermal_cc%thermal_conductivity_function => TCFDefaultCreate()
        call TCFAssignDefault(thermal_cc%thermal_conductivity_function,&
          patch%material_property_array(i)%ptr%thermal_conductivity_wet, &
          patch%material_property_array(i)%ptr%thermal_conductivity_dry, &
          patch%material_property_array(i)%ptr%alpha, &
          option)
      endif
      if (associated(thermal_cc%thermal_conductivity_function)) then
        write (mat_string,*) patch%material_property_array(i)%ptr%external_id
        verify_string = 'THERMAL_CHARACTERISTIC_CURVES(' // &
          trim(thermal_cc%name) // ') used for material ID #'// &
          trim(adjustl(mat_string)) // '. '
        call thermal_cc%thermal_conductivity_function% &
          Verify(verify_string,option)
      else
        write (mat_string,*) patch%material_property_array(i)%ptr%external_id
        option%io_buffer = 'A thermal conductivity function has &
          &not been set under THERMAL_CHARACTERISTIC_CURVES "' // &
          trim(thermal_cc%name) // '" intended for material ID #'// &
          trim(adjustl(mat_string)) // '. '
      endif
      call CharCurvesThermalAddToList(thermal_cc, &
        realization%characteristic_curves_thermal)
      nullify(thermal_cc)
    enddo
    ! afterwards, proceed with normal TCC procedure
    if (associated(realization%characteristic_curves_thermal)) then
      patch%characteristic_curves_thermal => &
           realization%characteristic_curves_thermal
      call CharCurvesThermalConvertListToArray( &
         patch%characteristic_curves_thermal, &
         patch%char_curves_thermal_array, option)
      do i = 1, size(patch%char_curves_thermal_array)
         select type(tcf => patch%char_curves_thermal_array(i)%ptr% &
                     thermal_conductivity_function)
         class is(kT_composite_type)
           call CompositeTCCList(patch%characteristic_curves_thermal, &
                                 tcf,option)
         end select
      enddo
    else
      option%io_buffer = 'Manual assignments of DEFAULT thermal '//&
                         'characteristic curve failed!'
      call PrintErrMsg(option)
    endif
  endif
  deallocate(check_thermal_conductivity)

  ! create mapping of internal to external material id
  call MaterialCreateIntToExtMapping(patch%material_property_array, &
                                     patch%imat_internal_to_external)

  cur_material_property => realization%material_properties
  do
    if (.not.associated(cur_material_property)) exit

    ! obtain saturation function id
    select case(option%iflowmode)
      case(NULL_MODE,PNF_MODE)
      case default
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
        endif
    end select

    ! thermal conducitivity function id
    if (associated(patch%char_curves_thermal_array)) then
      if (cur_material_property%thermal_conductivity_function_id < 1) then
        cur_material_property%thermal_conductivity_function_id = &
           CharCurvesThermalGetID( &
           patch%char_curves_thermal_array, &
           cur_material_property%thermal_conductivity_func_name, &
           cur_material_property%name,option)
      endif
    endif
    if (cur_material_property%thermal_conductivity_function_id == 0) then
      option%io_buffer = 'Thermal characteristic curve "' // &
        trim(cur_material_property%thermal_conductivity_func_name) // &
        '" not found.'
        call PrintErrMsg(option)
    endif

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
    if (associated(cur_material_property%electrical_conductivity_dataset)) then
      string = 'MATERIAL_PROPERTY(' // trim(cur_material_property%name) // &
               '),ELECTRICAL_CONDUCTIVITY'
      dataset => &
        DatasetBaseGetPointer(realization%datasets, &
                              cur_material_property% &
                                electrical_conductivity_dataset%name, &
                              string,option)
      call DatasetDestroy(cur_material_property%electrical_conductivity_dataset)
      select type(dataset)
        class is (dataset_common_hdf5_type)
          cur_material_property%electrical_conductivity_dataset => dataset
        class default
          option%io_buffer = 'Incorrect dataset type for electrical &
                             &conductivity.'
          call PrintErrMsg(option)
      end select
    endif
    if (associated(cur_material_property%multicontinuum)) then
      if (associated(cur_material_property%multicontinuum%epsilon_dataset)) then
        string = 'MATERIAL_PROPERTY(' // trim(cur_material_property%name) // &
                 '),EPSILON'
        dataset => &
          DatasetBaseGetPointer(realization%datasets, &
                                cur_material_property% &
                                  multicontinuum%epsilon_dataset%name, &
                                string,option)
        call DatasetDestroy(cur_material_property% &
                              multicontinuum%epsilon_dataset)
        select type(dataset)
          class is (dataset_common_hdf5_type)
            cur_material_property%multicontinuum%epsilon_dataset => dataset
          class default
            option%io_buffer = 'Incorrect dataset type for epsilon.'
            call PrintErrMsg(option)
        end select
      endif
      if (associated(cur_material_property%multicontinuum%half_matrix_width_dataset)) then
        string = 'MATERIAL_PROPERTY(' // trim(cur_material_property%name) // &
                 '),LENGTH'
        dataset => &
          DatasetBaseGetPointer(realization%datasets, &
                                cur_material_property% &
                                  multicontinuum%half_matrix_width_dataset%name, &
                                string,option)
        call DatasetDestroy(cur_material_property% &
                              multicontinuum%half_matrix_width_dataset)
        select type(dataset)
          class is (dataset_common_hdf5_type)
            cur_material_property%multicontinuum%half_matrix_width_dataset => dataset
          class default
            option%io_buffer = 'Incorrect dataset type for length.'
            call PrintErrMsg(option)
        end select
      endif
      if (associated(cur_material_property%multicontinuum%ncells_dataset)) then
        string = 'MATERIAL_PROPERTY(' // trim(cur_material_property%name) // &
                 '),NUM_CELLS'
        dataset => &
          DatasetBaseGetPointer(realization%datasets, &
                                cur_material_property% &
                                  multicontinuum%ncells_dataset%name, &
                                string,option)
        call DatasetDestroy(cur_material_property% &
                              multicontinuum%ncells_dataset)
        select type(dataset)
          class is (dataset_common_hdf5_type)
            cur_material_property%multicontinuum%ncells_dataset => dataset
          class default
            option%io_buffer = 'Incorrect dataset type for number of secondary cells.'
            call PrintErrMsg(option)
        end select
      endif
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
    if (associated(cur_material_property%permeability_dataset_xy)) then
      string = 'MATERIAL_PROPERTY(' // trim(cur_material_property%name) // &
               '),PERMEABILITY XY'
      dataset => &
        DatasetBaseGetPointer(realization%datasets, &
                          cur_material_property%permeability_dataset_xy%name, &
                          string,option)
      call DatasetDestroy(cur_material_property%permeability_dataset_xy)
      select type(dataset)
        class is (dataset_common_hdf5_type)
          cur_material_property%permeability_dataset_xy => dataset
        class default
          option%io_buffer = 'Incorrect dataset type for permeability XY.'
          call PrintErrMsg(option)
      end select
    endif
    if (associated(cur_material_property%permeability_dataset_xz)) then
      string = 'MATERIAL_PROPERTY(' // trim(cur_material_property%name) // &
               '),PERMEABILITY XZ'
      dataset => &
        DatasetBaseGetPointer(realization%datasets, &
                          cur_material_property%permeability_dataset_xz%name, &
                          string,option)
      call DatasetDestroy(cur_material_property%permeability_dataset_xz)
      select type(dataset)
        class is (dataset_common_hdf5_type)
          cur_material_property%permeability_dataset_xz => dataset
        class default
          option%io_buffer = 'Incorrect dataset type for permeability XZ.'
          call PrintErrMsg(option)
      end select
    endif
    if (associated(cur_material_property%permeability_dataset_yz)) then
      string = 'MATERIAL_PROPERTY(' // trim(cur_material_property%name) // &
               '),PERMEABILITY YZ'
      dataset => &
        DatasetBaseGetPointer(realization%datasets, &
                          cur_material_property%permeability_dataset_yz%name, &
                          string,option)
      call DatasetDestroy(cur_material_property%permeability_dataset_yz)
      select type(dataset)
        class is (dataset_common_hdf5_type)
          cur_material_property%permeability_dataset_yz => dataset
        class default
          option%io_buffer = 'Incorrect dataset type for permeability YZ.'
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

  implicit none

  class(realization_subsurface_type) :: realization

  PetscBool :: found
  type(option_type), pointer :: option
  type(fluid_property_type), pointer :: cur_fluid_property

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
    if (option%use_sc) then
      cur_constraint_coupler => cur_condition%sec_constraint_coupler
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
              exit
            endif
            cur_constraint => cur_constraint%next
          enddo
          if (.not.associated(cur_constraint_coupler%constraint)) then
            option%io_buffer = 'Secondary constraint "' // &
                     trim(cur_constraint_coupler%constraint_name) // &
                     '" not found in input file constraints.'
            call PrintErrMsg(realization%option)
          endif
        endif
        cur_constraint_coupler => cur_constraint_coupler%next
      enddo
    endif

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
  class(tran_constraint_coupler_rt_type), pointer :: constraint_rt_coupler

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
        ! the following three lines are a work around for an intel compiler bug
        ! claiming that c is not a pointer, though it points to a pointer
        !call ReactionPrintConstraint(c%global_auxvar,c%rt_auxvar,c, &
        constraint_rt_coupler => c
        call ReactionPrintConstraint(c%global_auxvar,c%rt_auxvar, &
                                     constraint_rt_coupler,reaction,option)
        write(option%fid_out,'(/)')
        write(option%fid_out,99)
    end select
  endif

end subroutine RealizationPrintCoupler

! ************************************************************************** !

subroutine RealizationPrintStateAtCells(realization)
  !
  ! loops over cells and prints their reactive transport states
  !
  ! Author: Glenn Hammond
  ! Date: 04/21/23

  use String_module

  implicit none

  class(realization_subsurface_type) :: realization

  PetscInt :: i
  PetscInt :: icell

  if (.not.associated(realization%reaction%print_cells)) return

  if (realization%option%comm%size > 1) then
    call PrintErrMsg(realization%option,'Printing of cell states for &
           &reactive transport not supported in parallel.')
  endif
  i = maxval(realization%reaction%print_cells)
  if (i > realization%discretization%grid%nmax) then
    realization%option%io_buffer = 'A cell id (' // &
      StringWrite(i) // ') specified under CHEMISTRY,&
      &OUTPUT,PRINT_CELLS is larger than the maximum cell id (' // &
      StringWrite(realization%discretization%grid%nmax) // ').'
    call PrintErrMsg(realization%option)
  endif

  write(realization%option%fid_out,'(/,40("+="),//,&
        &"  States at select cells")')
  do i = 1, size(realization%reaction%print_cells)
    icell = realization%reaction%print_cells(i)
    write(realization%option%fid_out,'(/,80("-"),//,"  State at cell ",a,/)') &
      StringWrite(icell)
    call RealizationPrintStateAtCell( &
           realization%patch%aux%Global%auxvars(icell), &
           realization%patch%aux%RT%auxvars(icell), &
           realization%reaction,realization%option)
  enddo
  if (i > 1) write(realization%option%fid_out,'(/,40("+="),/)')

end subroutine RealizationPrintStateAtCells

! ************************************************************************** !

subroutine RealizationPrintStateAtCell(global_auxvar,rt_auxvar, &
                                       reaction,option)
  !
  ! Prints reactive transport state variabiables for a given grid cell
  !
  ! Author: Glenn Hammond
  ! Date: 04/21/23
  !
  use Global_Aux_module
  use Option_module
  use Reaction_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Transport_Constraint_RT_module

  implicit none

  type(global_auxvar_type) :: global_auxvar
  type(reactive_transport_auxvar_type) :: rt_auxvar
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  class(tran_constraint_coupler_rt_type), pointer :: null_constraint_coupler

  nullify(null_constraint_coupler)

  call ReactionPrintConstraint(global_auxvar,rt_auxvar, &
                               null_constraint_coupler, &
                               reaction,option)

end subroutine RealizationPrintStateAtCell

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
  use Material_Aux_module, only : material_type, &
                              POROSITY_CURRENT, POROSITY_BASE, POROSITY_INITIAL
  use Variables_module, only : PERMEABILITY_X, PERMEABILITY_Y, PERMEABILITY_Z, &
                               PERMEABILITY_XY, PERMEABILITY_XZ, &
                               PERMEABILITY_YZ, POROSITY

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
    if (option%flow%full_perm_tensor) then
      call DiscretizationGlobalToLocal(discretization,field%perm0_xy, &
                                       field%work_loc,ONEDOF)
      call MaterialSetAuxVarVecLoc(Material,field%work_loc,PERMEABILITY_XY, &
                                   ZERO_INTEGER)
      call DiscretizationGlobalToLocal(discretization,field%perm0_xz, &
                                       field%work_loc,ONEDOF)
      call MaterialSetAuxVarVecLoc(Material,field%work_loc,PERMEABILITY_XZ, &
                                   ZERO_INTEGER)
      call DiscretizationGlobalToLocal(discretization,field%perm0_yz, &
                                       field%work_loc,ONEDOF)
      call MaterialSetAuxVarVecLoc(Material,field%work_loc,PERMEABILITY_YZ, &
                                   ZERO_INTEGER)
     endif
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
  use Material_Aux_module
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
    if (option%flow%full_perm_tensor) then
      call MaterialGetAuxVarVecLoc(Material,field%work_loc,PERMEABILITY_XY, &
                                   ZERO_INTEGER)
      call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                       field%perm0_xy,ONEDOF)
      call MaterialGetAuxVarVecLoc(Material,field%work_loc,PERMEABILITY_XZ, &
                                   ZERO_INTEGER)
      call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                       field%perm0_xz,ONEDOF)
      call MaterialGetAuxVarVecLoc(Material,field%work_loc,PERMEABILITY_YZ, &
                                   ZERO_INTEGER)
      call DiscretizationLocalToGlobal(discretization,field%work_loc, &
                                       field%perm0_yz,ONEDOF)
    endif
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
  PetscReal :: final_time
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
    if (FlowConditionHasRateOrFlux(cur_flow_condition) .or. &
        cur_flow_condition%sync_time_with_update) then
      do isub_condition = 1, cur_flow_condition%num_sub_conditions
        sub_condition => cur_flow_condition%sub_condition_ptr(isub_condition)%ptr
        !TODO(geh): check if this updated more than simply the flow_dataset (i.e. datum and gradient)
        !geh: followup - no, datum/gradient are not considered.  Should they be considered?
        call TimeStorageGetTimes(sub_condition%dataset%time_storage, option, &
                                final_time, times)
        if (associated(times)) then
          if (size(times) > 20000) then
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
  use Material_Aux_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Reaction_Mineral_module
  use Variables_module, only : POROSITY, TORTUOSITY, PERMEABILITY_X, &
                               PERMEABILITY_Y, PERMEABILITY_Z, &
                               PERMEABILITY_XY, PERMEABILITY_XZ, &
                               PERMEABILITY_YZ

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
  type(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: local_id, ghosted_id
  PetscInt :: imnrl, imat
  PetscReal :: scale
  PetscBool :: porosity_updated
  PetscReal, pointer :: porosity0_p(:)
  PetscReal, pointer :: tortuosity0_p(:)
  PetscReal, pointer :: perm0_xx_p(:), perm0_yy_p(:), perm0_zz_p(:)
  PetscReal, pointer :: perm0_xy_p(:), perm0_xz_p(:), perm0_yz_p(:)
  PetscReal :: min_value
  PetscReal :: critical_porosity
  PetscReal :: porosity_base_
  PetscReal :: temp_porosity
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

  if (reaction%update_mineral_surface_area) then

    nullify(porosity0_p)
    if (reaction%update_mnrl_surf_with_porosity) then
      ! placing the get/restore array calls within the condition will
      ! avoid improper access.
      call VecGetArrayReadF90(field%porosity0,porosity0_p,ierr);CHKERRQ(ierr)
    endif

    temp_porosity = UNINITIALIZED_DOUBLE
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      do imnrl = 1, reaction%mineral%nkinmnrl
        if (associated(porosity0_p)) then
          temp_porosity = porosity0_p(local_id)
        endif
        call MineralUpdateSpecSurfaceArea(reaction,rt_auxvars(ghosted_id), &
                                          material_auxvars(ghosted_id), &
                                          temp_porosity,option)
      enddo
    enddo

    if (reaction%update_mnrl_surf_with_porosity) then
      call VecRestoreArrayReadF90(field%porosity0,porosity0_p, &
                                  ierr);CHKERRQ(ierr)
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
    call VecGetArrayReadF90(field%tortuosity0,tortuosity0_p, &
                            ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(field%porosity0,porosity0_p,ierr);CHKERRQ(ierr)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      scale = (material_auxvars(ghosted_id)%porosity_base / &
               porosity0_p(local_id))** &
        material_property_array(patch%imat(ghosted_id))%ptr%tortuosity_pwr
      material_auxvars(ghosted_id)%tortuosity = &
        tortuosity0_p(local_id)*scale
    enddo
    call VecRestoreArrayReadF90(field%tortuosity0,tortuosity0_p, &
                                ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(field%porosity0,porosity0_p, &
                                ierr);CHKERRQ(ierr)
    call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 TORTUOSITY,ZERO_INTEGER)
    call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                    field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                 TORTUOSITY,ZERO_INTEGER)
  endif

  if (reaction%update_permeability) then
    call VecGetArrayReadF90(field%perm0_xx,perm0_xx_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(field%perm0_zz,perm0_zz_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(field%perm0_yy,perm0_yy_p,ierr);CHKERRQ(ierr)
    if (option%flow%full_perm_tensor) then
      call VecGetArrayReadF90(field%perm0_xy,perm0_xy_p,ierr);CHKERRQ(ierr)
      call VecGetArrayReadF90(field%perm0_xz,perm0_xz_p,ierr);CHKERRQ(ierr)
      call VecGetArrayReadF90(field%perm0_yz,perm0_yz_p,ierr);CHKERRQ(ierr)
    endif
    call VecGetArrayReadF90(field%porosity0,porosity0_p,ierr);CHKERRQ(ierr)
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
      material_auxvars(ghosted_id)%permeability(perm_xx_index) = &
          perm0_xx_p(local_id)*scale
      material_auxvars(ghosted_id)%permeability(perm_yy_index) = &
          perm0_yy_p(local_id)*scale
      material_auxvars(ghosted_id)%permeability(perm_zz_index) = &
          perm0_zz_p(local_id)*scale
      if (option%flow%full_perm_tensor) then
        material_auxvars(ghosted_id)%permeability(perm_xy_index) = &
          perm0_xy_p(local_id)*scale
        material_auxvars(ghosted_id)%permeability(perm_xz_index) = &
          perm0_xz_p(local_id)*scale
        material_auxvars(ghosted_id)%permeability(perm_yz_index) = &
          perm0_yz_p(local_id)*scale
      endif
    enddo
    call VecRestoreArrayReadF90(field%perm0_xx,perm0_xx_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(field%perm0_zz,perm0_zz_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(field%perm0_yy,perm0_yy_p,ierr);CHKERRQ(ierr)
    if (option%flow%full_perm_tensor) then
      call VecRestoreArrayReadF90(field%perm0_xy,perm0_xy_p, &
                                  ierr);CHKERRQ(ierr)
      call VecRestoreArrayReadF90(field%perm0_xz,perm0_xz_p, &
                                  ierr);CHKERRQ(ierr)
      call VecRestoreArrayReadF90(field%perm0_yz,perm0_yz_p, &
                                  ierr);CHKERRQ(ierr)
    endif
    call VecRestoreArrayReadF90(field%porosity0,porosity0_p, &
                                ierr);CHKERRQ(ierr)

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
    if (option%flow%full_perm_tensor) then
      call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                   PERMEABILITY_XY,ZERO_INTEGER)
      call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                      field%work_loc,ONEDOF)
      call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                   PERMEABILITY_XY,ZERO_INTEGER)
      call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                   PERMEABILITY_XZ,ZERO_INTEGER)
      call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                      field%work_loc,ONEDOF)
      call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                   PERMEABILITY_XZ,ZERO_INTEGER)
      call MaterialGetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                   PERMEABILITY_YZ,ZERO_INTEGER)
      call DiscretizationLocalToLocal(discretization,field%work_loc, &
                                      field%work_loc,ONEDOF)
      call MaterialSetAuxVarVecLoc(patch%aux%Material,field%work_loc, &
                                   PERMEABILITY_YZ,ZERO_INTEGER)
    endif
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

  implicit none

  class(realization_subsurface_type) :: realization

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
  use Material_Aux_module
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
  type(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: local_id, ghosted_id
  PetscInt :: imnrl
  PetscReal :: sum_volfrac

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
    case(CC_ID_ARRAY)
      call GridCopyIntegerArrayToVec(grid,patch%cc_id, &
                                     field%work_loc, grid%ngmax)
    case(CCT_ID_ARRAY)
      call GridCopyIntegerArrayToVec(grid,patch%cct_id, &
                                     field%work_loc, grid%ngmax)
    case(MTF_ID_ARRAY)
      call GridCopyIntegerArrayToVec(grid,patch%mtf_id, &
                                     field%work_loc, grid%ngmax)
  end select

  call DiscretizationLocalToLocal(realization%discretization,field%work_loc, &
                                  field%work_loc,ONEDOF)

  select case(array_id)
    case(MATERIAL_ID_ARRAY)
      call GridCopyVecToIntegerArray(grid,patch%imat,field%work_loc, &
                                      grid%ngmax)
    case(CC_ID_ARRAY)
      call GridCopyVecToIntegerArray(grid,patch%cc_id, &
                                      field%work_loc, grid%ngmax)
    case(CCT_ID_ARRAY)
      call GridCopyVecToIntegerArray(grid,patch%cct_id, &
                                      field%work_loc, grid%ngmax)
    case(MTF_ID_ARRAY)
      call GridCopyVecToIntegerArray(grid,patch%mtf_id, &
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
                     MPI_SUM,realization%option%mycomm,ierr);CHKERRQ(ierr)
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
  PetscReal :: r1
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
                     MPIU_INTEGER,MPI_SUM,option%mycomm,ierr);CHKERRQ(ierr)

  ! why I cannot use *100, I do not know....geh
  inactive_percentages = dble(temp_int_out)/dble(option%comm%size)*10.d0
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
           option%comm%size, &
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
           option%comm%size, &
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

subroutine RealizationCalculateCFL1Timestep(realization,max_dt_cfl_1, &
                                            max_pore_velocity)
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
  PetscReal :: max_pore_velocity

  max_dt_cfl_1 = MAX_DOUBLE
  max_pore_velocity = 0.d0
  call PatchCalculateCFL1Timestep(realization%patch,realization%option, &
                                  max_dt_cfl_1, &
                                  max_pore_velocity)

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
  use Material_Aux_module
  use Variables_module, only : VOLUME, BASE_POROSITY, PERMEABILITY_X, &
                               PERMEABILITY_Y, PERMEABILITY_Z, &
                               PERMEABILITY_XY, PERMEABILITY_XZ, &
                               PERMEABILITY_YZ

  implicit none

  class(realization_subsurface_type) :: realization

  PetscInt :: i

  call RealizUnInitializedVar1(realization,VOLUME,'volume')
  ! mineral porosity is the base, unmodified porosity
  call RealizUnInitializedVar1(realization,BASE_POROSITY,'porosity')
  call RealizUnInitializedVar1(realization,PERMEABILITY_X,'permeability X')
  call RealizUnInitializedVar1(realization,PERMEABILITY_Y,'permeability Y')
  call RealizUnInitializedVar1(realization,PERMEABILITY_Z,'permeability Z')
  if (realization%option%flow%full_perm_tensor) then
    call RealizUnInitializedVar1(realization,PERMEABILITY_XY,'permeability XY')
    call RealizUnInitializedVar1(realization,PERMEABILITY_XZ,'permeability XZ')
    call RealizUnInitializedVar1(realization,PERMEABILITY_YZ,'permeability YZ')
  endif
  do i = 1, max_material_index
    call RealizUnInitializedVar1(realization, &
                   realization%patch%aux%Material%soil_properties_ivar(i), &
                   realization%patch%aux%Material%soil_properties_name(i))
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
  use Material_Aux_module
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
    call PrintErrMsgToDev(option,'send your input deck')
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
  PetscReal :: max_pore_velocity
  PetscReal :: prev_dt
  PetscBool :: print_to_screen, print_to_file
  character(len=MAXSTRINGLENGTH) :: string
  type(output_option_type), pointer :: output_option

  if (Initialized(cfl_governor)) then
    call RealizationCalculateCFL1Timestep(realization,max_dt_cfl_1, &
                                          max_pore_velocity)
    print_to_screen = OptionPrintToScreen(realization%option)
    print_to_file = OptionPrintToFile(realization%option)
    if (print_to_screen .or. print_to_file) then
      output_option => realization%output_option
      write(string,'(" Maximum Pore Velocity: ",1pe12.4," [m/",a,"]")') &
        max_pore_velocity*output_option%tconv,trim(output_option%tunit)
      if (print_to_screen) write(STDOUT_UNIT,'(a)') trim(string)
      if (print_to_file) write(realization%option%fid_out,'(a)') trim(string)
    endif
    if (dt/cfl_governor > max_dt_cfl_1) then
      prev_dt = dt
      dt = max_dt_cfl_1*cfl_governor
      ! have to set dt_max here so that timestepper%dt_max is truncated
      ! for timestepper_base%revert_dt in TimestepperBaseSetTargetTime
      dt_max = dt
      if (print_to_screen .or. print_to_file) then
        write(string,'(" CFL Limiting (",f4.1,"): ",1pe12.4," -> ",1pe12.4, &
              &" [",a,"]")') &
              cfl_governor,prev_dt/output_option%tconv, &
              dt/output_option%tconv,trim(output_option%tunit)
        if (print_to_screen) write(STDOUT_UNIT,'(a)') trim(string)
        if (print_to_file) write(realization%option%fid_out,'(a)') trim(string)
      endif
    endif
  endif

end subroutine RealizationLimitDTByCFL

! ************************************************************************** !

subroutine RealizationReadGeopSurveyFile(realization)
  !
  ! Read geophysics survey file
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 02/01/21
  !

  use Input_Aux_module
  use Option_module
  use Grid_module
  use Survey_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(survey_type), pointer :: survey

  type(input_type), pointer :: input_tmp

  option => realization%option
  grid => realization%patch%grid
  survey => realization%survey

  if (.not.associated(survey)) then
    option%io_buffer = 'There should be a SURVEY card in input file for &
                        &geophysics process models.'
    call PrintErrMsg(option)
  endif

  input_tmp => InputCreate(IUNIT_TEMP,survey%filename,option)

  if (option%igeopmode == ERT_MODE) then
    call SurveyReadERT(survey,grid,input_tmp,option)
  else
    option%io_buffer = 'Only ERT mode is supported for &
        &RealizationReadGeopSurveyFile.'
    call PrintErrMsg(option)
  endif
  call InputDestroy(input_tmp)

end subroutine RealizationReadGeopSurveyFile

! ************************************************************************** !

subroutine RealizationCheckConsistency(r1,r2)
  !
  ! Checks to ensure that two separate realizations have the same
  ! sized objects (used in inversion calculations to ensure that the
  ! prerequisite have the same sized data structures as the main
  ! forward run)
  !
  ! Author: Glenn Hammond
  ! Date: 12/02/22
  !
  use String_module

  implicit none

  class(realization_subsurface_type) :: r1, r2

  type(option_type), pointer :: o1, o2
  PetscBool :: error_found
  character(len=MAXSTRINGLENGTH) :: s
  PetscInt :: i
  PetscErrorCode :: ierr

  o1 => r1%option
  o2 => r2%option

  s = new_line('a') // &
      'Checking consistency of outer and prerequisite realizations:' // &
      new_line('a')
  call PrintMsg(o1,s)

  error_found = PETSC_FALSE
  s = 'Grid type'
  if (IntegersDiffer(r1%patch%grid%itype,r2%patch%grid%itype,s)) then
    error_found = PETSC_TRUE
    call PrintMsg(o1,s)
  endif
  s = 'Number of grid cells'
  if (IntegersDiffer(r1%patch%grid%nmax,r2%patch%grid%nmax,s)) then
    error_found = PETSC_TRUE
    call PrintMsg(o1,s)
  endif
  s = 'Material IDs differ'
  i = min(size(r1%patch%imat),size(r2%patch%imat))
  i = maxval(abs(r1%patch%imat(1:i)-r2%patch%imat(1:i)))
  call MPI_Allreduce(MPI_IN_PLACE,i,ONE_INTEGER_MPI,MPI_INTEGER, &
                     MPI_MAX,o1%mycomm,ierr);CHKERRQ(ierr)
  if (i > 0) then
    error_found = PETSC_TRUE
    call PrintMsg(o1,s)
  endif
  s = 'Number of material properties'
  if (IntegersDiffer(size(r1%patch%material_property_array), &
                     size(r2%patch%material_property_array),s)) then
    error_found = PETSC_TRUE
    call PrintMsg(o1,s)
  endif
  s = 'Number of characteristic curves'
  if (IntegersDiffer(size(r1%patch%characteristic_curves_array), &
                     size(r2%patch%characteristic_curves_array),s)) then
    error_found = PETSC_TRUE
    call PrintMsg(o1,s)
  endif
  s = 'Number of flow degrees of freedom'
  if (IntegersDiffer(o1%nflowdof,o2%nflowdof,s)) then
    error_found = PETSC_TRUE
    call PrintMsg(o1,s)
  endif
  s = 'Number of transport degrees of freedom'
  if (IntegersDiffer(o1%ntrandof,o2%ntrandof,s)) then
    error_found = PETSC_TRUE
    call PrintMsg(o1,s)
  endif

  if (error_found) then
    call PrintErrMsg(o1,'Inconsistent realizations.')
  endif
  call PrintMsg(o1,'Realizations are consistent.')

contains

  function IntegersDiffer(i1,i2,error_string)
    PetscInt :: i1, i2
    character(*) :: error_string
    PetscBool :: IntegersDiffer
    IntegersDiffer = (i1 /= i2)
    if (IntegersDiffer) then
      error_string = trim(error_string) // ' differ: ' // &
                     trim(StringWrite(i1)) // ' vs ' // &
                     trim(StringWrite(i2))
    endif
  end function IntegersDiffer

end subroutine RealizationCheckConsistency

! ************************************************************************** !

subroutine RealizationStrip(this)
  !
  ! Deallocates a realization
  !
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  !

  use Dataset_module

  implicit none

  class(realization_subsurface_type) :: this

  call RealizationBaseStrip(this)
  call RegionDestroyList(this%region_list)

  call FlowConditionDestroyList(this%flow_conditions)
  call TranConditionDestroyList(this%transport_conditions)
  call TranConstraintListDestroy(this%transport_constraints)
  call GeopConditionDestroyList(this%geophysics_conditions)

  if (associated(this%fluid_property_array)) &
    deallocate(this%fluid_property_array)
  nullify(this%fluid_property_array)
  call FluidPropertyDestroy(this%fluid_properties)

  call MaterialPropertyDestroy(this%material_properties)

  call SaturationFunctionDestroy(this%saturation_functions)
  call CharacteristicCurvesDestroy(this%characteristic_curves)

  if (associated(this%characteristic_curves_thermal)) then
    call CharCurvesThermalDestroy(this%characteristic_curves_thermal)
  endif

  nullify(this%material_transform)

  call DatasetDestroy(this%datasets)

  call DatasetDestroy(this%uniform_velocity_dataset)

  ! nullify since they are pointers to reaction_base in realization_base
  nullify(this%reaction)
  nullify(this%reaction_nw)

  if (associated(this%survey)) then
    call SurveyDestroy(this%survey)
  endif

end subroutine RealizationStrip

end module Realization_Subsurface_class
