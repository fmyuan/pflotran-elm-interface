
module Factory_Subsurface_Geomechanics_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: FactorySubsurfGeomechInitSimulation!, & ! used in factory_subsurf
contains

! ************************************************************************** !

subroutine FactorySubsurfGeomechInitSimulation(simulation, pm_geomech)

  use Simulation_Subsurface_class
  use Init_Common_module
  use Option_module
  use PM_Base_class
  use PM_Base_Pointer_module
  use PM_Geomechanics_Force_class
  use PMC_Base_class
  use PMC_Geomechanics_class
  use PFLOTRAN_Constants_module
  use Geomechanics_Discretization_module
  use Geomechanics_Force_module
  use Geomechanics_Realization_class
  use Geomechanics_Regression_module
  use Simulation_Aux_module
  use Realization_Subsurface_class
  use Realization_Base_class
  use Timestepper_Steady_class
  use Input_Aux_module
  use Logging_module
  use Output_Aux_module
  use Waypoint_module
  use Init_Subsurface_Geomech_module

  implicit none

  class(simulation_subsurface_type) :: simulation
  class(pm_geomech_force_type), pointer :: pm_geomech

  type(option_type), pointer :: option
  class(realization_subsurface_type), pointer :: subsurf_realization
  class(realization_geomech_type), pointer :: geomech_realization
  class(pmc_base_type), pointer :: cur_process_model_coupler
  type(gmdm_ptr_type), pointer :: dm_ptr
  class(pmc_geomechanics_type), pointer :: pmc_geomech
  class(timestepper_steady_type), pointer :: timestepper
  type(geomechanics_regression_type), pointer :: geomech_regression
  PetscErrorCode :: ierr

  if (.not. associated(pm_geomech)) return

  option => simulation%option
  geomech_realization => simulation%geomech%realization
  subsurf_realization => simulation%realization

  geomech_regression => simulation%geomech%regression

  ! initialize geomech realization
  !call SubsurfGeomechInitSetupRealization(simulation)
  !call InitSubsurfGeomechSetupRealization(simulation%realization,&
  !                                        simulation%geomech%realization)
  !call InitSubsurfGeomechSetupRealization(simulation)
  call SubsurfGeomechInitSetupRealization(simulation)

  call pm_geomech%PMGeomechForceSetRealization(geomech_realization)
  call pm_geomech%Setup()

  !pmc_geomech => GeomechPMC(simulation)
  pmc_geomech => simulation%geomech%process_model_coupler
  timestepper => TimestepperSteadyCast(pmc_geomech%timestepper)
  call pmc_geomech%SetupSolvers()

  ! Here I first calculate the linear part of the jacobian and store it
  ! since the jacobian is always linear with geomech (even when coupled with
  ! flow since we are performing sequential coupling). Although
  ! SNESSetJacobian is called, nothing is done there and PETSc just
  ! re-uses the linear Jacobian at all iterations and times
  call MatSetOption(timestepper%solver%M,MAT_NEW_NONZERO_ALLOCATION_ERR, &
                    PETSC_FALSE,ierr);CHKERRQ(ierr)
  call GeomechForceJacobianLinearPart(timestepper%solver%M, &
                                      geomech_realization)
  call MatSetOption(timestepper%solver%M,MAT_NEW_NONZERO_ALLOCATION_ERR, &
                    PETSC_TRUE,ierr);CHKERRQ(ierr)
  nullify(simulation%process_model_coupler_list)

  ! sim_aux: Create PETSc Vectors and VectorScatters
  if (option%ngeomechdof > 0) then

    call GeomechCreateGeomechSubsurfVec(subsurf_realization, &
                                        geomech_realization)
    call SimAuxCopySubsurfVec(simulation%sim_aux,subsurf_realization%field%work)

    call GeomechCreateSubsurfStressStrainVec(subsurf_realization, &
                                             geomech_realization)
    call SimAuxCopySubsurfGeomechVec(simulation%sim_aux, &
          geomech_realization%geomech_field%strain_subsurf)

    call GeomechRealizMapSubsurfGeomechGrid(subsurf_realization, &
                                            geomech_realization, &
                                            option)

    dm_ptr => GeomechDiscretizationGetDMPtrFromIndex( &
                geomech_realization%geomech_discretization, ONEDOF)

    call SimAuxCopyVecScatter(simulation%sim_aux, &
                              dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                              SUBSURF_TO_GEOMECHANICS)
    call SimAuxCopyVecScatter(simulation%sim_aux, &
                              dm_ptr%gmdm%scatter_geomech_to_subsurf_ndof, &
                              GEOMECHANICS_TO_SUBSURF)
  endif

  call GeomechanicsRegressionCreateMapping(geomech_regression, &
                                           geomech_realization)

  ! sim_aux: Set pointer
  simulation%flow_process_model_coupler%sim_aux => simulation%sim_aux
  if (associated(simulation%tran_process_model_coupler)) &
    simulation%tran_process_model_coupler%sim_aux => simulation%sim_aux
  if (option%ngeomechdof>0 .and. associated(pmc_geomech)) &
    pmc_geomech%sim_aux => simulation%sim_aux

  ! set geomech as not master
  pmc_geomech%is_master = PETSC_FALSE
  ! link geomech and master
  !simulation%process_model_coupler_list => &
  !  simulation%geomech_process_model_coupler_new
  ! jaa testing.. set flow as the master
  simulation%process_model_coupler_list => &
    simulation%flow_process_model_coupler
  ! link subsurface flow as peer
  !simulation%process_model_coupler_list%peer => &
  !  simulation%flow_process_model_coupler
  ! jaa testing.. set geomech as a child
  simulation%process_model_coupler_list%child => &
    pmc_geomech
  ! Set data in sim_aux
  cur_process_model_coupler => simulation%process_model_coupler_list
  call cur_process_model_coupler%SetAuxData()
  !if (associated(cur_process_model_coupler%peer)) then
  !  cur_process_model_coupler => cur_process_model_coupler%peer
  if (associated(cur_process_model_coupler%child)) then
    cur_process_model_coupler => cur_process_model_coupler%child
    call cur_process_model_coupler%GetAuxData()
    call cur_process_model_coupler%SetAuxData()
    select type(pmc => cur_process_model_coupler)
      class is(pmc_geomechanics_type)
        call GeomechStoreInitialPressTemp(pmc%geomech_realization)
    end select
  endif

  call GeomechanicsJumpStart(simulation%geomech)

end subroutine FactorySubsurfGeomechInitSimulation

! ************************************************************************** !

!subroutine FactorySubsurfGeomechReadSimBlock(input,pm)
!  !
!  ! Author: Piyoosh Jaysaval
!  ! Date: 01/25/21
!  !
!  use Input_Aux_module
!  use Option_module
!  use String_module
!
!  use PM_Base_class
!  use PM_ERT_class
!
!  implicit none
!
!  type(input_type), pointer :: input
!  class(pm_base_type), pointer :: pm
!
!  type(option_type), pointer :: option
!  character(len=MAXWORDLENGTH) :: word
!  character(len=MAXSTRINGLENGTH) :: error_string
!
!  option => pm%option
!
!  error_string = 'SIMULATION,PROCESS_MODELS,SUBSURFACE_GEOMECHANICS'
!
!  call InputPushBlock(input,option)
!  do
!    call InputReadPflotranString(input,option)
!    if (InputCheckExit(input,option)) exit
!    call InputReadCard(input,option,word,PETSC_FALSE)
!    call StringToUpper(word)
!    select case(word)
!      case('OPTIONS')
!        call pm%ReadSimulationOptionsBlock(input)
!      case default
!        call InputKeywordUnrecognized(input,word,error_string,option)
!    end select
!  enddo
!  call InputPopBlock(input,option)
!
!end subroutine FactorySubsurfGeomechReadSimBlock

! ************************************************************************** !

!subroutine SubsurfGeomechicsInitReadRequiredCards(geomech_realization,input)
!  !
!  ! Reads the required input file cards
!  ! related to geomechanics
!  !
!  ! Author: Satish Karra, LANL
!  ! Date: 05/23/13
!  !
!
!  use Geomechanics_Discretization_module
!  use Geomechanics_Realization_class
!  use Geomechanics_Patch_module
!  use Geomechanics_Grid_module
!  use Input_Aux_module
!  use String_module
!  use Patch_module
!  use Option_module
!
!  implicit none
!
!  class(realization_geomech_type) :: geomech_realization
!  type(input_type), pointer :: input
!
!  character(len=MAXSTRINGLENGTH) :: string
!  type(option_type), pointer :: option
!
!  option => geomech_realization%option
!
!! Read in select required cards
!!.........................................................................
!
!  ! GEOMECHANICS information
!  string = "GEOMECHANICS"
!  call InputFindStringInFile(input,option,string)
!  if (InputError(input)) return
!  option%ngeomechdof = 3  ! displacements in x, y, z directions
!  option%n_stress_strain_dof = 6
!
!  string = "GEOMECHANICS_GRID"
!  call InputFindStringInFile(input,option,string)
!  call GeomechanicsInit(geomech_realization,input,option)
!
!
!end subroutine SubsurfGeomechicsInitReadRequiredCards

! ************************************************************************** !

!subroutine GeomechanicsInit(geomech_realization,input,option)
!  !
!  ! Reads the required geomechanics data from input file
!  !
!  ! Author: Satish Karra, LANL
!  ! Date: 05/23/13
!  !
!
!  use Option_module
!  use Input_Aux_module
!  use String_module
!  use Geomechanics_Grid_module
!  use Geomechanics_Grid_Aux_module
!  use Geomechanics_Discretization_module
!  use Geomechanics_Realization_class
!  use Geomechanics_Patch_module
!  use Grid_Unstructured_Aux_module
!  use Grid_Unstructured_module
!
!  implicit none
!
!  class(realization_geomech_type) :: geomech_realization
!  type(geomech_discretization_type), pointer :: geomech_discretization
!  type(geomech_patch_type), pointer :: patch
!  type(input_type), pointer :: input
!  type(option_type), pointer :: option
!  character(len=MAXWORDLENGTH) :: word
!  type(grid_unstructured_type), pointer :: ugrid
!  character(len=MAXWORDLENGTH) :: card
!
!  geomech_discretization       => geomech_realization%geomech_discretization
!
!  input%ierr = 0
!  ! we initialize the word to blanks to avoid error reported by valgrind
!  word = ''
!
!  call InputPushBlock(input,option)
!  do
!    call InputReadPflotranString(input,option)
!    call InputReadStringErrorMsg(input,option,card)
!    if (InputCheckExit(input,option)) exit
!    call InputReadCard(input,option,word)
!    call InputErrorMsg(input,option,'keyword','GEOMECHANICS')
!    call StringToUpper(word)
!
!    select case(trim(word))
!      case ('TYPE')
!        call InputReadCard(input,option,word)
!        call InputErrorMsg(input,option,'keyword','TYPE')
!        call StringToUpper(word)
!
!        select case(trim(word))
!          case ('UNSTRUCTURED')
!            geomech_discretization%itype = UNSTRUCTURED_GRID
!            call InputReadFilename(input,option,geomech_discretization%filename)
!            call InputErrorMsg(input,option,'keyword','filename')
!
!            geomech_discretization%grid  => GMGridCreate()
!            ugrid => UGridCreate()
!            call UGridRead(ugrid,geomech_discretization%filename,option)
!            call UGridDecompose(ugrid,option)
!            call CopySubsurfaceGridtoGeomechGrid(ugrid, &
!                                                 geomech_discretization%grid, &
!                                                 option)
!            patch => GeomechanicsPatchCreate()
!            patch%geomech_grid => geomech_discretization%grid
!            geomech_realization%geomech_patch => patch
!          case default
!            option%io_buffer = 'Geomechanics supports only unstructured grid'
!            call PrintErrMsg(option)
!        end select
!      case ('GRAVITY')
!        call InputReadDouble(input,option,option%geomech_gravity(X_DIRECTION))
!        call InputErrorMsg(input,option,'x-direction','GEOMECH GRAVITY')
!        call InputReadDouble(input,option,option%geomech_gravity(Y_DIRECTION))
!        call InputErrorMsg(input,option,'y-direction','GEOMECH GRAVITY')
!        call InputReadDouble(input,option,option%geomech_gravity(Z_DIRECTION))
!        call InputErrorMsg(input,option,'z-direction','GEOMECH GRAVITY')
!        if (OptionIsIORank(option) .and. OptionPrintToScreen(option)) &
!            write(option%fid_out,'(/," *GEOMECH_GRAV",/, &
!            & "  gravity    = "," [m/s^2]",3x,1p3e12.4 &
!            & )') option%geomech_gravity(1:3)
!      case default
!        call InputKeywordUnrecognized(input,word,'GEOMECHANICS_GRID',option)
!    end select
!  enddo
!  call InputPopBlock(input,option)
!
!end subroutine GeomechanicsInit

! ************************************************************************** !

subroutine SubsurfGeomechInitMatPropToGeomechRegions(geomech_realization)
  !
  ! This routine assigns geomech material
  ! properties to associated regions
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  !
  use Geomechanics_Realization_class
  use Geomechanics_Discretization_module
  use Geomechanics_Strata_module
  use Geomechanics_Region_module
  use Geomechanics_Material_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Field_module
  use Geomechanics_Patch_module
  use Realization_Subsurface_class, only : MATERIAL_ID_ARRAY
  use Option_module

  implicit none

  class(realization_geomech_type) :: geomech_realization

  PetscInt :: ivertex, local_id, ghosted_id, geomech_material_id
  PetscInt :: istart, iend
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscErrorCode :: ierr

  type(option_type), pointer :: option
  type(geomech_grid_type), pointer :: grid
  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_field_type), pointer :: field
  type(geomech_strata_type), pointer :: strata
  type(geomech_patch_type), pointer :: patch

  type(geomech_material_property_type), pointer :: geomech_material_property
  type(geomech_material_property_type), pointer :: null_geomech_material_property
  type(gm_region_type), pointer :: region
  PetscBool :: update_ghosted_material_ids
  PetscReal, pointer :: imech_loc_p(:)

  option => geomech_realization%option
  geomech_discretization => geomech_realization%geomech_discretization
  field => geomech_realization%geomech_field
  patch => geomech_realization%geomech_patch

  ! loop over all patches and allocation material id arrays
  if (.not.associated(patch%imat)) then
    allocate(patch%imat(patch%geomech_grid%ngmax_node))
    ! initialize to "unset"
    patch%imat = UNINITIALIZED_INTEGER
  endif

  ! if material ids are set based on region, as opposed to being read in
  ! we must communicate the ghosted ids.  This flag toggles this operation.
  update_ghosted_material_ids = PETSC_FALSE
  grid => patch%geomech_grid
  strata => patch%geomech_strata_list%first
  do
    if (.not.associated(strata)) exit
    ! Read in cell by cell material ids if they exist
    if (.not.associated(strata%region) .and. strata%active) then
      option%io_buffer = 'Reading of material prop from file for' // &
        ' geomech is not implemented.'
      call PrintErrMsgByRank(option)
    ! Otherwise, set based on region
    else if (strata%active) then
      update_ghosted_material_ids = PETSC_TRUE
      region => strata%region
      geomech_material_property => strata%material_property
      if (associated(region)) then
        istart = 1
        iend = region%num_verts
      else
        istart = 1
        iend = grid%nlmax_node
      endif
      do ivertex = istart, iend
        if (associated(region)) then
          local_id = region%vertex_ids(ivertex)
        else
          local_id = ivertex
        endif
        ghosted_id = grid%nL2G(local_id)
        patch%imat(ghosted_id) = geomech_material_property%id
      enddo
    endif
    strata => strata%next
  enddo

  if (update_ghosted_material_ids) then
    ! update ghosted material ids
    call GeomechRealizLocalToLocalWithArray(geomech_realization, &
                                            MATERIAL_ID_ARRAY)
  endif

  ! set cell by cell material properties
  ! create null material property for inactive cells
  null_geomech_material_property => GeomechanicsMaterialPropertyCreate()
  call VecGetArrayF90(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)
  do local_id = 1, grid%nlmax_node
    ghosted_id = grid%nL2G(local_id)
    geomech_material_id = patch%imat(ghosted_id)
    if (geomech_material_id == 0) then ! accomodate inactive cells
      geomech_material_property = null_geomech_material_property
    else if ( geomech_material_id > 0 .and. &
              geomech_material_id <= &
              size(geomech_realization%geomech_material_property_array)) then
      geomech_material_property => &
         geomech_realization% &
           geomech_material_property_array(geomech_material_id)%ptr
      if (.not.associated(geomech_material_property)) then
        write(dataset_name,*) geomech_material_id
        option%io_buffer = 'No material property for geomech material id ' // &
                            trim(adjustl(dataset_name)) &
                            //  ' defined in input file.'
        call PrintErrMsgByRank(option)
      endif
    else if (Uninitialized(geomech_material_id)) then
      write(dataset_name,*) grid%nG2A(ghosted_id)
      option%io_buffer = 'Uninitialized geomech material id in patch at cell ' // &
                         trim(adjustl(dataset_name))
      call PrintErrMsgByRank(option)
    else if (geomech_material_id > size(geomech_realization% &
      geomech_material_property_array)) then
      write(option%io_buffer,*) geomech_material_id
      option%io_buffer = 'Unmatched geomech material id in patch:' // &
        adjustl(trim(option%io_buffer))
      call PrintErrMsgByRank(option)
    else
      option%io_buffer = 'Something messed up with geomech material ids. ' // &
        ' Possibly material ids not assigned to all grid cells. ' // &
        ' Contact Glenn/Satish!'
      call PrintErrMsgByRank(option)
    endif
    imech_loc_p(ghosted_id) = geomech_material_property%id
  enddo ! local_id - loop
  call VecRestoreArrayF90(field%imech_loc,imech_loc_p,ierr);CHKERRQ(ierr)

  call GeomechanicsMaterialPropertyDestroy(null_geomech_material_property)
  nullify(null_geomech_material_property)

  call GeomechDiscretizationLocalToLocal(geomech_discretization,field%imech_loc, &
                                         field%imech_loc,ONEDOF)

end subroutine SubsurfGeomechInitMatPropToGeomechRegions

! ************************************************************************** !

subroutine SubsurfGeomechInitSetupRealization(simulation)
  !
  ! Initializes material property data structres and assign them to the domain.
  !
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  !
  use Simulation_Geomechanics_class
  use Geomechanics_Realization_class
  use Geomechanics_Global_module
  use Geomechanics_Force_module
  use Realization_Subsurface_class
  use Simulation_Subsurface_class

  use Option_module
  use Waypoint_module

  implicit none

  class(simulation_subsurface_type) :: simulation

  class(realization_subsurface_type), pointer :: subsurf_realization
  class(realization_geomech_type), pointer :: geomech_realization
  type(option_type), pointer :: option

  subsurf_realization => simulation%realization
  geomech_realization => simulation%geomech%realization
  option => subsurf_realization%option

  call GeomechRealizCreateDiscretization(geomech_realization)

  if (option%geomech_subsurf_coupling /= 0) then
    call GeomechCreateGeomechSubsurfVec(subsurf_realization, &
                                        geomech_realization)
    call GeomechCreateSubsurfStressStrainVec(subsurf_realization, &
                                              geomech_realization)

    call GeomechRealizMapSubsurfGeomechGrid(subsurf_realization, &
                                            geomech_realization, &
                                            option)
  endif
  call GeomechRealizLocalizeRegions(geomech_realization)
  call GeomechRealizPassFieldPtrToPatch(geomech_realization)
  call GeomechRealizProcessMatProp(geomech_realization)
  call GeomechRealizProcessGeomechCouplers(geomech_realization)
  call GeomechRealizProcessGeomechConditions(geomech_realization)
  call SubsurfGeomechInitMatPropToGeomechRegions(geomech_realization)
  call GeomechRealizInitAllCouplerAuxVars(geomech_realization)
  call GeomechRealizPrintCouplers(geomech_realization)
  call GeomechGridElemSharedByNodes(geomech_realization,option)
  call GeomechForceSetup(geomech_realization)
  call GeomechGlobalSetup(geomech_realization)

  ! SK: We are solving quasi-steady state solution for geomechanics.
  ! Initial condition is not needed, hence CondControlAssignFlowInitCondGeomech
  ! is not needed, at this point.
  call GeomechForceUpdateAuxVars(geomech_realization)

end subroutine SubsurfGeomechInitSetupRealization

end module Factory_Subsurface_Geomechanics_module

