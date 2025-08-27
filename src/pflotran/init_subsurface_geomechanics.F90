module Init_Subsurface_Geomech_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: InitSubsurfGeomechReadRequiredCards, &
            InitSubsurfGeomechReadInput, &
            InitSubsurfGeomechJumpStart, & ! remove later
            InitSubsurfGeomechSetupRealization, &
            InitSubsurfGeomechInitSimulation, &
            InitSubsurfGeomechSetGeomechMode, &
            InitSubsurfGeomechChkInactiveCells, &
            InitSubsurfGeomechSetupPMC, &
            InitSubsurfGeomechReadSimBlock
contains

! ************************************************************************** !

subroutine InitSubsurfGeomechReadRequiredCards(geomech_realization,input)
  !
  ! Reads the required input file cards
  ! related to geomechanics
  !
  ! Author: Satish Karra, LANL
  ! Date: 05/23/13
  !
  ! jaa: moved from factory_geomechanics.F90 on 1/28/25

  use Geomechanics_Discretization_module
  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Geomechanics_Grid_module
  use Input_Aux_module
  use String_module
  use Patch_module
  use Option_module

  implicit none

  class(realization_geomech_type) :: geomech_realization
  type(input_type), pointer :: input

  character(len=MAXSTRINGLENGTH) :: string
  type(option_type), pointer :: option

  option => geomech_realization%option

! Read in select required cards
!.........................................................................

  ! GEOMECHANICS information
  string = "GEOMECHANICS"
  call InputFindStringInFile(input,option,string)
  if (InputError(input)) return

  string = "GEOMECHANICS_GRID"
  call InputFindStringInFile(input,option,string)
  call InitSubsurfGeomechReadGridBlock(geomech_realization,input,option)

end subroutine InitSubsurfGeomechReadRequiredCards

! ************************************************************************** !

subroutine InitSubsurfGeomechReadInput(geomech,geomech_solver, &
                                     input,option,output_option)
  !
  ! Reads the geomechanics input data
  !
  ! Author: Satish Karra, LANL
  ! Date: 05/23/13
  !
  ! jaa: moved from factory_geomechanics.F90 on 1/28/25

  use Option_module
  use Input_Aux_module
  use String_module
  use Geomechanics_Discretization_module
  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Material_module
  use Geomechanics_Region_module
  use Geomechanics_Debug_module
  use Geomechanics_Strata_module
  use Geomechanics_Condition_module
  use Geomechanics_Coupler_module
  use Geomechanics_Regression_module
  use Output_Aux_module
  use Output_Tecplot_module
  use Realization_Base_class
  use Solver_module
  use Units_module
  use Waypoint_module
  use Dataset_Base_class
  use Dataset_module
  use Dataset_Common_HDF5_class
  use Utility_module, only : DeallocateArray, UtilityReadArray
  use Geomechanics_Attr_module

  ! Still need to add other geomech modules for output, etc once created

  implicit none

  type(solver_type), pointer :: geomech_solver
  type(input_type), pointer :: input
  type(geomechanics_attr_type), pointer:: geomech
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option

  class(realization_geomech_type), pointer :: geomech_realization
  class(dataset_base_type), pointer :: dataset

  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_material_property_type),pointer :: geomech_material_property
  type(geomech_grid_type), pointer :: grid
  type(gm_region_type), pointer :: region
  type(geomech_strata_type), pointer :: strata
  type(geomech_condition_type), pointer :: condition
  type(geomech_coupler_type), pointer :: coupler
  type(waypoint_list_type), pointer :: waypoint_list

  character(len=MAXWORDLENGTH) :: word, internal_units
  character(len=MAXWORDLENGTH) :: card
  character(len=1) :: backslash

  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
  input%ierr = 0
! we initialize the word to blanks to avoid error reported by valgrind
  word = ''

  waypoint_list => geomech%waypoint_list
  geomech_realization => geomech%realization
  geomech_discretization => geomech_realization%geomech_discretization

  if (associated(geomech_realization%geomech_patch)) grid => &
    geomech_realization%geomech_patch%geomech_grid

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword','GEOMECHANICS')
    call StringToUpper(word)
    option%io_buffer = 'word :: ' // trim(word)
    call PrintMsg(option)

    select case(trim(word))

      !.........................................................................
      ! Read geomechanics grid information
      case ('GEOMECHANICS_GRID')
        call InputSkipToEND(input,option,trim(word))

      !.........................................................................
      ! Read geomechanics material information
      case ('GEOMECHANICS_MATERIAL_PROPERTY')
        geomech_material_property => GeomechanicsMaterialPropertyCreate()

        call InputReadWord(input,option,geomech_material_property%name, &
                           PETSC_TRUE)

        call InputErrorMsg(input,option,'name','GEOMECHANICS_MATERIAL_PROPERTY')
        call GeomechanicsMaterialPropertyRead(geomech_material_property,input, &
                                              option)
        call GeomechanicsMaterialPropertyAddToList(geomech_material_property, &
                                geomech_realization%geomech_material_properties)
        nullify(geomech_material_property)

      case ('GEOMECHANICS_SET_REF_P_T_TO_IC')
        option%geomechanics%set_ref_pres_and_temp_to_IC = PETSC_TRUE

      !.........................................................................
      ! Read geomechanics datasets
      case ('GEOMECHANICS_DATASET')
          nullify(dataset)
          call DatasetRead(input,dataset,option)
          call DatasetBaseAddToList(dataset,geomech_realization%geomech_datasets)
          nullify(dataset)

      !.........................................................................
      case ('GEOMECHANICS_REGION')
        region => GeomechRegionCreate()
        call InputReadWord(input,option,region%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'name','GEOMECHANICS_REGION')
        call PrintMsg(option,region%name)
        call GeomechRegionRead(region,input,option)
        ! we don't copy regions down to patches quite yet, since we
        ! don't want to duplicate IO in reading the regions
        call GeomechRegionAddToList(region,geomech_realization%geomech_region_list)
        nullify(region)

      !.........................................................................
      case ('GEOMECHANICS_CONDITION')
        condition => GeomechConditionCreate(option)
        call InputReadWord(input,option,condition%name,PETSC_TRUE)
        call InputErrorMsg(input,option,'GEOMECHANICS_CONDITION','name')
        call PrintMsg(option,condition%name)
        call GeomechConditionRead(condition,input,option)
        call GeomechConditionAddToList(condition,geomech_realization%geomech_conditions)
        nullify(condition)

     !.........................................................................
      case ('GEOMECHANICS_BOUNDARY_CONDITION')
        coupler =>  GeomechCouplerCreate(GM_BOUNDARY_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Geomech Boundary Condition name')
        call GeomechCouplerRead(coupler,input,option)
        call GeomechRealizAddGeomechCoupler(geomech_realization,coupler)
        nullify(coupler)

      !.........................................................................
      case ('GEOMECHANICS_SRC_SINK')
        coupler => GeomechCouplerCreate(GM_SRC_SINK_COUPLER_TYPE)
        call InputReadWord(input,option,coupler%name,PETSC_TRUE)
        call InputDefaultMsg(input,option,'Source Sink name')
        call GeomechCouplerRead(coupler,input,option)
        call GeomechRealizAddGeomechCoupler(geomech_realization,coupler)
        nullify(coupler)

     !....................
      case ('GEOMECHANICS_LINEAR_SOLVER')
        call SolverReadLinear(geomech_solver,input,option)

      !.....................
      case ('GEOMECHANICS_REGRESSION')
        call GeomechanicsRegressionRead(geomech%regression,input,option)

      !.........................................................................
      case ('GEOMECHANICS_TIME')
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'word','GEOMECHANICS_TIME')
          select case(trim(word))
            case('COUPLING_TIMESTEP_SIZE')
              call InputReadDouble(input,option,geomech_realization%dt_coupling)
              call InputErrorMsg(input,option, &
                                 'Coupling Timestep Size','GEOMECHANICS_TIME')
              internal_units = 'sec'
              call InputReadAndConvertUnits(input, &
                                            geomech_realization%dt_coupling, &
                                            internal_units,'GEOMECHANICS_TIME,&
                                            &COUPLING_TIMESTEP_SIZE',option)
            case default
              call InputKeywordUnrecognized(input,word, &
                                            'GEOMECHANICS_TIME',option)
            end select
        enddo
        call InputPopBlock(input,option)

      !.........................................................................
      case ('GEOMECHANICS_DEBUG')
        call GeomechDebugRead(geomech_realization%geomech_debug,input,option)

      !.........................................................................
      case ('GEOMECHANICS_MAPPING_FILE')
        call InputReadFilename(input,option,grid%mapping_filename)
        call InputErrorMsg(input,option,'keyword','mapping_file')
        call GeomechSubsurfMapFromFilename(grid,grid%mapping_filename,option)

      !.........................................................................
      case ('GEOMECHANICS_OUTPUT')
        call InputPushBlock(input,option)
        do
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,card)
          if (InputCheckExit(input,option)) exit
          call InputReadCard(input,option,word)
          call InputErrorMsg(input,option,'keyword','GEOMECHANICS_OUTPUT')
          call StringToUpper(word)
          select case(trim(word))
            case('TIMES')
              option%io_buffer = 'Subsurface times are now used for ' // &
              'geomechanics as well. No need for TIMES keyword under ' // &
              'GEOMECHANICS_OUTPUT.'
              call PrintWrnMsg(option)
            case('FORMAT')
              call InputReadCard(input,option,word)
              call InputErrorMsg(input,option,'keyword','GEOMECHANICS_OUTPUT,&
                                                         &FORMAT')
              call StringToUpper(word)
              select case(trim(word))
                case ('HDF5')
                  output_option%print_hdf5 = PETSC_TRUE
                  call InputReadCard(input,option,word)
                  call InputDefaultMsg(input,option, &
                                       'GEOMECHANICS_OUTPUT,FORMAT,HDF5,&
                                        &# FILES')
                  if (len_trim(word) > 1) then
                    call StringToUpper(word)
                    select case(trim(word))
                      case('SINGLE_FILE')
                        output_option%print_single_h5_file = PETSC_TRUE
                      case('MULTIPLE_FILES')
                        output_option%print_single_h5_file = PETSC_FALSE
                      case default
                        option%io_buffer = 'HDF5 keyword (' // trim(word) // &
                          ') not recongnized.  Use "SINGLE_FILE" or ' // &
                          '"MULTIPLE_FILES".'
                        call PrintErrMsg(option)
                    end select
                  endif
                case ('TECPLOT')
                  output_option%print_tecplot = PETSC_TRUE
                  call InputReadCard(input,option,word)
                  call InputErrorMsg(input,option,'TECPLOT','GEOMECHANICS_OUTPUT,FORMAT')
                  call StringToUpper(word)
                  output_option%tecplot_format = TECPLOT_FEQUADRILATERAL_FORMAT ! By default it is unstructured
                case ('VTK')
                  output_option%print_vtk = PETSC_TRUE
                case default
                  call InputKeywordUnrecognized(input,word, &
                                 'GEOMECHANICS_OUTPUT,FORMAT',option)
              end select
            case default
              call InputKeywordUnrecognized(input,word, &
                             'GEOMECHANICS_OUTPUT',option)
          end select
        enddo
        call InputPopBlock(input,option)

      !.........................................................................
      case ('GEOMECHANICS_STRATIGRAPHY','GEOMECHANICS_STRATA')
        strata => GeomechStrataCreate()
        call GeomechStrataRead(strata,input,option)
        call GeomechRealizAddStrata(geomech_realization,strata)
        nullify(strata)

      !.........................................................................
      case ('GEOMECHANICS_IMPROVE_TET_WEIGHT')
        option%geomechanics%improve_tet_weighting = PETSC_TRUE

      !.........................................................................
      case ('END_GEOMECHANICS')
        exit

      !.........................................................................
      case default
        call InputKeywordUnrecognized(input,word, &
                                 'GeomechanicsInitReadInput',option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine InitSubsurfGeomechReadInput

! ************************************************************************** !

subroutine InitSubsurfGeomechJumpStart(geomech)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  !
  ! jaa: moved from factory_geomechanics.F90 on 1/28/25

  use Geomechanics_Realization_class
  use Option_module
  use Timestepper_KSP_class
  use Output_Aux_module
  use Output_module, only : Output, OutputPrintCouplers
  use Output_Geomechanics_module
  use Logging_module
  use Condition_Control_module
  use Simulation_Subsurface_class
  use PMC_Geomechanics_class
  use Geomechanics_Attr_module

  implicit none

  type(geomechanics_attr_type) :: geomech

  class(realization_geomech_type), pointer :: geomech_realization
  class(pmc_geomechanics_type), pointer :: geomech_pmc
  class(timestepper_ksp_type), pointer :: geomech_timestepper

  PetscBool :: snapshot_plot_flag,observation_plot_flag,massbal_plot_flag
  PetscBool :: geomech_read
  PetscBool :: failure
  PetscErrorCode :: ierr

  geomech_realization => geomech%realization
  geomech_pmc => geomech%process_model_coupler
  geomech_timestepper => TimestepperKSPCast(geomech_pmc%timestepper)

  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           "-vecload_block_size",failure,ierr);CHKERRQ(ierr)

  geomech_timestepper%name = 'GEOMECHANICS'

  snapshot_plot_flag = PETSC_FALSE
  observation_plot_flag = PETSC_FALSE
  massbal_plot_flag = PETSC_FALSE
  geomech_read = PETSC_FALSE
  failure = PETSC_FALSE

  call OutputGeomechInit(geomech_timestepper%steps)

  ! pushed in INIT_STAGE()
  call PetscLogStagePop(ierr);CHKERRQ(ierr)

  ! popped in TS_STAGE()
  call PetscLogStagePush(logging%stage(TS_STAGE),ierr);CHKERRQ(ierr)

end subroutine InitSubsurfGeomechJumpStart

! ************************************************************************** !

subroutine InitSubsurfGeomechReadGridBlock(geomech_realization,input,option)
  !
  ! Reads the required geomechanics data from input file
  !
  ! Author: Satish Karra, LANL
  ! Date: 05/23/13
  !
  ! jaa: moved from factory_geomechanics.F90 on 1/28/25

  use Option_module
  use Input_Aux_module
  use String_module
  use Geomechanics_Grid_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Discretization_module
  use Geomechanics_Realization_class
  use Geomechanics_Patch_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_module

  implicit none

  class(realization_geomech_type) :: geomech_realization
  type(input_type), pointer :: input
  type(option_type), pointer :: option

  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_patch_type), pointer :: patch
  character(len=MAXWORDLENGTH) :: word
  type(grid_unstructured_type), pointer :: ugrid
  character(len=MAXWORDLENGTH) :: card

  geomech_discretization => geomech_realization%geomech_discretization

  input%ierr = 0
  ! we initialize the word to blanks to avoid error reported by valgrind
  word = ''

  geomech_discretization%grid  => GMGridCreate()
  ugrid => UGridCreate()

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    call InputReadStringErrorMsg(input,option,card)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword','GEOMECHANICS')
    call StringToUpper(word)

    select case(trim(word))
      case ('TYPE')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'keyword','TYPE')
        call StringToUpper(word)

        select case(trim(word))
          case ('UNSTRUCTURED')
            geomech_discretization%itype = UNSTRUCTURED_GRID
            call InputReadFilename(input,option,geomech_discretization%filename)
            call InputErrorMsg(input,option,'keyword','filename')
          case default
            option%io_buffer = 'Geomechanics supports only unstructured grid'
            call PrintErrMsg(option)
        end select
      case ('GRAVITY')
        call InputReadDouble(input,option,option%geomechanics% &
                             gravity(X_DIRECTION))
        call InputErrorMsg(input,option,'x-direction','GEOMECH GRAVITY')
        call InputReadDouble(input,option,option%geomechanics% &
                             gravity(Y_DIRECTION))
        call InputErrorMsg(input,option,'y-direction','GEOMECH GRAVITY')
        call InputReadDouble(input,option,option%geomechanics% &
                             gravity(Z_DIRECTION))
        call InputErrorMsg(input,option,'z-direction','GEOMECH GRAVITY')
        if (OptionIsIORank(option) .and. OptionPrintToScreen(option)) &
            write(option%fid_out,'(/," *GEOMECH_GRAV",/, &
            & "  gravity    = "," [m/s^2]",3x,1p3e12.4 &
            & )') option%geomechanics%gravity(1:3)
      case ('MAX_CELLS_SHARING_A_VERTEX')
        call InputReadInt(input,option,ugrid%max_cells_sharing_a_vertex)
        call InputErrorMsg(input,option,'max_cells_sharing_a_vertex', &
                           'GEOMECHANICS_GRID')
      case default
        call InputKeywordUnrecognized(input,word,'GEOMECHANICS_GRID',option)
    end select
  enddo
  call InputPopBlock(input,option)

  call UGridRead(ugrid,geomech_discretization%filename,option)
  call UGridDecompose(ugrid,option)
  call CopySubsurfaceGridtoGeomechGrid(ugrid, &
                                       geomech_discretization%grid, &
                                       option)
  patch => GeomechanicsPatchCreate()
  patch%geomech_grid => geomech_discretization%grid
  geomech_realization%geomech_patch => patch

end subroutine InitSubsurfGeomechReadGridBlock

! ************************************************************************** !

subroutine InitSubsurfGeomechSetupRealization(subsurf_realization, &
                                              geomech_realization)
  !
  ! Initializes material property data structres and assign them to the domain.
  !
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  !
  ! jaa: moved from factory_geomechanics.F90 on 1/28/25

  use Geomechanics_Realization_class
  use Geomechanics_Global_module
  use Geomechanics_Force_module
  use Realization_Subsurface_class
  use Simulation_Subsurface_class

  use Option_module

  implicit none

  class(realization_subsurface_type), pointer :: subsurf_realization
  class(realization_geomech_type), pointer :: geomech_realization

  type(option_type), pointer :: option

  option => subsurf_realization%option

  call GeomechRealizCreateDiscretization(geomech_realization)

  if (option%geomechanics%flow_coupling /= 0 .or. &
      option%geomechanics%geophysics_coupling /= 0) then
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
  call InitMatPropToGeomechRegions(geomech_realization)
  call GeomechRealizInitAllCouplerAuxVars(geomech_realization)
  call GeomechRealizPrintCouplers(geomech_realization)
  call GeomechGridElemSharedByNodes(geomech_realization,option)
  call GeomechForceSetup(geomech_realization)
  call GeomechGlobalSetup(geomech_realization)

  ! SK: We are solving quasi-steady state solution for geomechanics.
  ! Initial condition is not needed, hence CondControlAssignFlowInitCondGeomech
  ! is not needed, at this point.
  call GeomechForceUpdateAuxVars(geomech_realization)

end subroutine InitSubsurfGeomechSetupRealization

! ************************************************************************** !

subroutine InitMatPropToGeomechRegions(geomech_realization)
  !
  ! This routine assigns geomech material
  ! properties to associated regions
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/17/13
  !
  ! jaa: moved from factory_geomechanics.F90 on 1/28/25

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

  Vec :: temp_vec
  PetscReal, pointer :: temp_vec_p(:)

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

  ! read in any user-defined geomech property fields
  call GeomechDiscretizationDuplicateVector(geomech_discretization, &
                                            field%press_loc, &
                                            temp_vec)
  do geomech_material_id = 1, size(patch%geomech_material_property_array)
    geomech_material_property => &
            patch%geomech_material_property_array(geomech_material_id)%ptr
    if (.not.associated(geomech_material_property)) cycle
    ! Young's modulus
    if (associated(geomech_material_property%youngs_modulus_dataset)) then
      ! Set the value of field%youngs_modulus for this material to the dataset values
      call GeomechReadDatasetToVecWithMask(geomech_realization, &
            geomech_material_property%youngs_modulus_dataset, &
            geomech_material_property%id,PETSC_FALSE,field%youngs_modulus,temp_vec)
    else
      ! Set the value of field%youngs_modulus for this material to the constant value
      call VecGetArrayF90(field%youngs_modulus,temp_vec_p,ierr);CHKERRQ(ierr)
      do local_id = 1, grid%nlmax_node
        if (patch%imat(grid%nL2G(local_id)) == geomech_material_property%id) then
          temp_vec_p(local_id) = geomech_material_property%youngs_modulus
        endif
      enddo
      call VecRestoreArrayF90(field%youngs_modulus,temp_vec_p,ierr);CHKERRQ(ierr)
    endif
    ! Poisson's ratio
    if (associated(geomech_material_property%poissons_ratio_dataset)) then
      ! Set the value of field%poissons_ratio for this material to the dataset values
      call GeomechReadDatasetToVecWithMask(geomech_realization, &
            geomech_material_property%poissons_ratio_dataset, &
            geomech_material_property%id,PETSC_FALSE,field%poissons_ratio,temp_vec)
    else
      ! Set the value of field%poissons_ratio for this material to the constant value
      call VecGetArrayF90(field%poissons_ratio,temp_vec_p,ierr);CHKERRQ(ierr)
      do local_id = 1, grid%nlmax_node
        if (patch%imat(grid%nL2G(local_id)) == geomech_material_property%id) then
          temp_vec_p(local_id) = geomech_material_property%poissons_ratio
        endif
      enddo
      call VecRestoreArrayF90(field%poissons_ratio,temp_vec_p,ierr);CHKERRQ(ierr)
    endif
    ! Density
    if (associated(geomech_material_property%density_dataset)) then
      ! Set the value of field%density for this material to the dataset values
      call GeomechReadDatasetToVecWithMask(geomech_realization, &
            geomech_material_property%density_dataset, &
            geomech_material_property%id,PETSC_FALSE,field%density,temp_vec)
    else
      ! Set the value of field%density for this material to the constant value
      call VecGetArrayF90(field%density,temp_vec_p,ierr);CHKERRQ(ierr)
      do local_id = 1, grid%nlmax_node
        if (patch%imat(grid%nL2G(local_id)) == geomech_material_property%id) then
          temp_vec_p(local_id) = geomech_material_property%density
        endif
      enddo
      call VecRestoreArrayF90(field%density,temp_vec_p,ierr);CHKERRQ(ierr)
    endif
    ! Biot's coefficient
    if (associated(geomech_material_property%biot_coeff_dataset)) then
      ! Set the value of field%biot_coeff for this material to the dataset values
      call GeomechReadDatasetToVecWithMask(geomech_realization, &
            geomech_material_property%biot_coeff_dataset, &
            geomech_material_property%id,PETSC_FALSE,field%biot_coeff,temp_vec)
    else
      ! Set the value of field%biot_coeff for this material to the constant value
      call VecGetArrayF90(field%biot_coeff,temp_vec_p,ierr);CHKERRQ(ierr)
      do local_id = 1, grid%nlmax_node
        if (patch%imat(grid%nL2G(local_id)) == geomech_material_property%id) then
          temp_vec_p(local_id) = geomech_material_property%biot_coeff
        endif
      enddo
      call VecRestoreArrayF90(field%biot_coeff,temp_vec_p,ierr);CHKERRQ(ierr)
    endif
    ! Thermal expansion coefficient
    if (associated(geomech_material_property%thermal_exp_coeff_dataset)) then
      ! Set the value of field%thermal_exp_coeff for this material to the dataset values
      call GeomechReadDatasetToVecWithMask(geomech_realization, &
            geomech_material_property%thermal_exp_coeff_dataset, &
            geomech_material_property%id,PETSC_FALSE,field%thermal_exp_coeff,temp_vec)
    else
      ! Set the value of field%thermal_exp_coeff for this material to the constant value
      call VecGetArrayF90(field%thermal_exp_coeff,temp_vec_p,ierr);CHKERRQ(ierr)
      do local_id = 1, grid%nlmax_node
        if (patch%imat(grid%nL2G(local_id)) == geomech_material_property%id) then
          temp_vec_p(local_id) = geomech_material_property%thermal_exp_coeff
        endif
      enddo
      call VecRestoreArrayF90(field%thermal_exp_coeff,temp_vec_p,ierr);CHKERRQ(ierr)
    endif
  enddo
  call VecDestroy(temp_vec,ierr);CHKERRQ(ierr)

  call GeomechanicsMaterialPropertyDestroy(null_geomech_material_property)
  nullify(null_geomech_material_property)

  call GeomechDiscretizationLocalToLocal(geomech_discretization,field%imech_loc, &
                                         field%imech_loc,ONEDOF)

end subroutine InitMatPropToGeomechRegions

! ************************************************************************** !

subroutine InitSubsurfGeomechInitSimulation(simulation, pm_geomech)
  !
  ! This routine initializes geomechanics process
  ! model components for factory subsurface
  ! after linkages have been established
  !
  ! Author: Jumanah Al Kubaisy
  ! Date: 2/10/25
  !
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
  use Timestepper_KSP_class
  use Logging_module
  use Output_Aux_module
  use Waypoint_module

  implicit none

  class(simulation_subsurface_type) :: simulation
  class(pm_geomech_force_type), pointer :: pm_geomech

  type(option_type), pointer :: option
  class(realization_subsurface_type), pointer :: subsurf_realization
  class(realization_geomech_type), pointer :: geomech_realization
  class(pmc_base_type), pointer :: cur_process_model_coupler
  class(pmc_base_type), pointer :: pmc_dummy
  type(gmdm_ptr_type), pointer :: dm_ptr
  class(pmc_geomechanics_type), pointer :: pmc_geomech
  class(timestepper_ksp_type), pointer :: timestepper
  type(geomechanics_regression_type), pointer :: geomech_regression
  PetscErrorCode :: ierr

  if (.not. associated(pm_geomech)) return

  nullify(pmc_dummy)

  option => simulation%option
  geomech_realization => simulation%geomech%realization
  subsurf_realization => simulation%realization
  pmc_geomech => simulation%geomech%process_model_coupler
  timestepper => TimestepperKSPCast(pmc_geomech%timestepper)

  geomech_regression => simulation%geomech%regression

  ! initialize geomech realization
  call InitSubsurfGeomechSetupRealization(simulation%realization,&
                                          simulation%geomech%realization)

  call pm_geomech%PMGeomechForceSetRealization(geomech_realization, &
                                               subsurf_realization)
  call pm_geomech%Setup()

  call pmc_geomech%SetupSolvers()

  ! Here I first calculate the linear part of the jacobian and store it
  ! since the jacobian is always linear with geomech (even when coupled with
  ! flow since we are performing sequential coupling). Although
  ! SNESSetJacobian is called, nothing is done there and PETSc just
  ! re-uses the linear Jacobian at all iterations and times
  call MatSetOption(timestepper%solver%M,MAT_NEW_NONZERO_ALLOCATION_ERR, &
                    PETSC_FALSE,ierr);CHKERRQ(ierr)
  call GeomechForceAssembleCoeffMatrix(timestepper%solver%M, &
                                      geomech_realization)
  call MatSetOption(timestepper%solver%M,MAT_NEW_NONZERO_ALLOCATION_ERR, &
                    PETSC_TRUE,ierr);CHKERRQ(ierr)
  nullify(simulation%process_model_coupler_list)

  ! sim_aux: Create PETSc Vectors and VectorScatters
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
  ! jaa: set flow as the master
  simulation%process_model_coupler_list => &
    simulation%flow_process_model_coupler
  ! link subsurface flow as peer
  ! jaa: set geomech as a child
!  simulation%process_model_coupler_list%child => &
!    pmc_geomech
  call PMCBaseSetChildPeerPtr(pmc_geomech%CastToBase(),PM_CHILD, &
                    simulation%flow_process_model_coupler%CastToBase(), &
                    pmc_dummy,PM_APPEND)

  call InitSubsurfGeomechChkInactiveCells(geomech_realization, &
                                          subsurf_realization)

  ! Set data in sim_aux
  cur_process_model_coupler => simulation%process_model_coupler_list
  call simulation%flow_process_model_coupler%SetAuxData()
  call simulation%geomech%process_model_coupler%GetAuxData()
  call simulation%geomech%process_model_coupler%SetAuxData()
  ! this is solely for casting to pmc geomech
  select type(pmc => simulation%geomech%process_model_coupler)
    class is(pmc_geomechanics_type)
      call GeomechStoreInitialPressTemp(pmc%geomech_realization)
  end select

  call InitSubsurfGeomechJumpStart(simulation%geomech)

end subroutine InitSubsurfGeomechInitSimulation

! ************************************************************************** !

subroutine InitSubsurfGeomechSetGeomechMode(pm_geomech,option)
  !
  ! sets geomech options
  !
  ! Author: Jumanah Al Kubaisy
  ! Date: 2/10/25
  !
  use Option_module
  use PM_Geomechanics_Force_class

  implicit none

  type(option_type) :: option
  class(pm_geomech_force_type), pointer :: pm_geomech

  if (.not.associated(pm_geomech)) then
    return
  endif

  select type(pm_geomech)
    class is (pm_geomech_force_type)
      option%igeommode = LINEAR_ELASTICITY_MODE
      option%geommode = "GEOMECHANICS"
      option%ngeomechdof = 3 ! displacements in x, y, z directions
      option%n_stress_strain_dof = 6
    class default
      option%io_buffer = 'Unrecognized geomechanics class in '// &
                          'InitSubsurfGeomechSetGeomechMode'
      call PrintErrMsg(option)
  end select

end subroutine InitSubsurfGeomechSetGeomechMode

! ************************************************************************** !

subroutine InitSubsurfGeomechChkInactiveCells(geomech_realization, &
                                             subsurf_realization)
  !
  ! checks if geomech nodes were mapped to inactive flow cells
  !
  ! Author: Glenn, Jumanah
  ! Date: 2/10/25
  !
  use Realization_Subsurface_class
  use Geomechanics_Realization_class
  use Geomechanics_Discretization_module
  use Option_module

  implicit none

  class(realization_subsurface_type) :: subsurf_realization
  class(realization_geomech_type) :: geomech_realization

  type(option_type), pointer :: option
  type(gmdm_ptr_type), pointer :: dm_ptr

  PetscErrorCode :: ierr
  PetscBool :: error_found
  PetscInt :: geomech_local_id, subsurf_local_id, geomech_ghosted_id
  PetscInt :: subsurf_ghosted_id
  PetscReal, pointer :: subsurf_vec_1dof(:)

  error_found = PETSC_FALSE

  option => geomech_realization%option
  dm_ptr => GeomechDiscretizationGetDMPtrFromIndex( &
            geomech_realization%geomech_discretization, ONEDOF)

  call VecSet(geomech_realization%geomech_field%subsurf_vec_1dof,-777.d0, &
              ierr);CHKERRQ(ierr)
  call VecSet(geomech_realization%geomech_field%press,-888.d0, &
              ierr);CHKERRQ(ierr)
  call VecGetArrayF90(geomech_realization%geomech_field%subsurf_vec_1dof, &
              subsurf_vec_1dof,ierr);CHKERRQ(ierr)
  do subsurf_local_id = 1, subsurf_realization%patch%grid%nlmax
    subsurf_ghosted_id = subsurf_realization%patch%grid%nL2G(subsurf_local_id)
    subsurf_vec_1dof(subsurf_local_id) = subsurf_realization%patch%imat( &
                                                  subsurf_ghosted_id)
  enddo
  call VecRestoreArrayF90(geomech_realization%geomech_field%subsurf_vec_1dof, &
                          subsurf_vec_1dof,ierr);CHKERRQ(ierr)
  ! Scatter the data
  call VecScatterBegin(dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                       geomech_realization%geomech_field%subsurf_vec_1dof, &
                       geomech_realization%geomech_field%press, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecScatterEnd(dm_ptr%gmdm%scatter_subsurf_to_geomech_ndof, &
                     geomech_realization%geomech_field%subsurf_vec_1dof, &
                     geomech_realization%geomech_field%press, &
                     INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(geomech_realization%geomech_field%press, &
                      subsurf_vec_1dof, ierr);CHKERRQ(ierr)
  do geomech_local_id = 1, geomech_realization%geomech_patch%geomech_grid% &
                           nlmax_node
    geomech_ghosted_id = geomech_realization%geomech_patch%geomech_grid% &
                         nL2G(geomech_local_id)
    if (nint(subsurf_vec_1dof(geomech_local_id)) <= 0) error_found = PETSC_TRUE
  enddo
  call MPI_Allreduce(MPI_IN_PLACE,error_found,ONE_INTEGER_MPI,MPI_LOGICAL, &
                     MPI_LOR,option%mycomm,ierr);CHKERRQ(ierr)
  if (error_found)then
    option%io_buffer = 'Cannot map inactive flow cell to geomechanics '//&
                       'node in the GEOMECHANICS_MAPPING_FILE! '
    call PrintErrMsg(option)
  endif

end subroutine InitSubsurfGeomechChkInactiveCells

! ************************************************************************** !

subroutine InitSubsurfGeomechSetupPMC(simulation,pm_geomech, &
                                     pmc_name,input)
  !
  ! refactored from factory_geomechanics.F90
  !
  ! Author: Jumanah Al Kubaisy
  ! Date: 2/10/25
  !
  use Realization_Subsurface_class
  use Option_module
  use Logging_module
  use Input_Aux_module
  use PM_Geomechanics_Force_class
  use Geomechanics_Realization_class
  use Timestepper_KSP_class
  use PMC_Geomechanics_class
  use Output_Aux_module
  use Waypoint_module
  use Simulation_Subsurface_class

  implicit none

  class(simulation_subsurface_type) :: simulation
  class(pm_geomech_force_type), pointer :: pm_geomech
  character(len=*) :: pmc_name
  type(input_type), pointer :: input

  character(len=MAXSTRINGLENGTH) :: string
  class(realization_subsurface_type), pointer :: subsurf_realization
  type(option_type), pointer :: option

  class(pmc_geomechanics_type), pointer :: pmc_geomech
  class(realization_geomech_type), pointer :: geomech_realization
  class(timestepper_ksp_type), pointer :: timestepper

  subsurf_realization => simulation%realization
  option => subsurf_realization%option
  subsurf_realization%output_option => simulation%output_option

  geomech_realization => GeomechRealizCreate(option)
  simulation%geomech%realization => geomech_realization

  input => InputCreate(IN_UNIT,option%input_filename,option)
  call InitSubsurfGeomechReadRequiredCards(geomech_realization,input)
  pmc_geomech => PMCGeomechanicsCreate()

  call pmc_geomech%SetName(pmc_name)
  call pmc_geomech%SetOption(option)
  simulation%geomech%process_model_coupler => pmc_geomech
  pmc_geomech%waypoint_list => simulation%waypoint_list_subsurface
  pmc_geomech%pm_list => pm_geomech
  pmc_geomech%pm_ptr%pm => pm_geomech
  pmc_geomech%geomech_realization => geomech_realization
  pm_geomech%geomech_realization => geomech_realization
  pmc_geomech%subsurf_realization => simulation%realization
  pm_geomech%subsurf_realization => simulation%realization

  ! add time integrator
  timestepper => TimestepperKSPCreate()
  pmc_geomech%timestepper => timestepper

  ! add solver
  call pm_geomech%InitializeSolver()
  timestepper%solver => pm_geomech%solver

  ! set up logging stage
  string = trim(pmc_geomech%name) // 'Geomechanics'
  call LoggingCreateStage(string,pmc_geomech%stage)

  string = 'GEOMECHANICS'
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)
  geomech_realization%output_option => &
    OutputOptionDuplicate(simulation%output_option)
  call OutputVariableListDestroy(geomech_realization%output_option% &
                                   output_snap_variable_list)
  call OutputVariableListDestroy(geomech_realization%output_option% &
                                   output_obs_variable_list)
  geomech_realization%output_option%output_snap_variable_list => &
    OutputVariableListCreate()
  geomech_realization%output_option%output_obs_variable_list => &
    OutputVariableListCreate()
  call InitSubsurfGeomechReadInput(simulation%geomech, &
                                   timestepper%solver, &
                                   input,option, &
                                   geomech_realization%output_option)
  pm_geomech%output_option => geomech_realization%output_option

  ! Hijack subsurface waypoint to geomechanics waypoint
  ! Subsurface controls the output now
  ! Always have snapshot on at t=0
  pmc_geomech%waypoint_list%first%print_snap_output = PETSC_TRUE

  ! link geomech and flow timestepper waypoints to geomech way point list
  if (associated(pmc_geomech)) then
    call pmc_geomech%SetWaypointPtr(pmc_geomech%waypoint_list)
    if (associated(simulation%flow_process_model_coupler)) then
      call simulation%flow_process_model_coupler% &
             SetWaypointPtr(pmc_geomech%waypoint_list)
    endif
  endif

  ! print the waypoints when debug flag is on
  if (geomech_realization%geomech_debug%print_waypoints) then
    call WaypointListPrint(pmc_geomech%waypoint_list,option, &
                           geomech_realization%output_option)
  endif

end subroutine InitSubsurfGeomechSetupPMC

! ************************************************************************** !

subroutine InitSubsurfGeomechReadSimBlock(input,pm)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/25/21
  !
  ! jaa: moved from factory_geomechanics.F90
  !
  use Input_Aux_module
  use Option_module
  use String_module

  use PM_Base_class
  use PM_ERT_class

  implicit none

  type(input_type), pointer :: input
  class(pm_base_type), pointer :: pm

  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  option => pm%option

  error_string = 'SIMULATION,PROCESS_MODELS,SUBSURFACE_GEOMECHANICS'

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    select case(word)
      case('OPTIONS')
        call pm%ReadSimulationOptionsBlock(input)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine InitSubsurfGeomechReadSimBlock

! ************************************************************************** !

subroutine GeomechReadDatasetToVecWithMask(geomech_realization,dataset, &
                                           geomech_material_id,read_all_values,vec,temp_vec)
  !
  ! Reads a geomechanics dataset into a PETSc Vec
  ! (based on SubsurfReadDatasetToVecWithMask)
  !
  ! Author: Kyle Mosley, WSP
  ! Date: 07/2025
  !
  use Geomechanics_Realization_class
  use Geomechanics_Field_module
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Patch_module
  use Option_module
  use Input_Aux_module
  Use HDF5_module
  Use Dataset_Base_class
  use Dataset_Common_HDF5_class
  use Dataset_Gridded_HDF5_class

  implicit none

  ! input/output declarations
  class(realization_geomech_type) :: geomech_realization
  class(dataset_base_type) :: dataset
  PetscInt :: geomech_material_id
  PetscBool :: read_all_values
  Vec :: vec
  Vec :: temp_vec

  type(geomech_field_type), pointer :: field
  type(geomech_patch_type), pointer :: patch
  type(geomech_grid_type), pointer :: grid
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscInt :: local_id
  PetscErrorCode :: ierr
  PetscReal, pointer :: vec_p(:)
  PetscReal, pointer :: work_p(:)

  field => geomech_realization%geomech_field
  patch => geomech_realization%geomech_patch
  grid => patch%geomech_grid
  option => geomech_realization%option

  call VecGetArrayF90(vec,vec_p,ierr);CHKERRQ(ierr)
  if (index(dataset%filename,'.h5') > 0) then
    group_name = ''
    dataset_name = dataset%name
    select type(dataset)
    class is(dataset_common_hdf5_type)


      dataset_name = dataset%hdf5_dataset_name
      call HDF5ReadCellIndexedRealArrayGM(geomech_realization,temp_vec, &
                                        dataset%filename, &
                                        group_name,dataset_name, &
                                        dataset%realization_dependent)
      call VecGetArrayF90(temp_vec,work_p,ierr);CHKERRQ(ierr)
      if (read_all_values) then
        do local_id = 1, grid%nlmax_node
          vec_p(local_id) = work_p(local_id)
        enddo
      else
        do local_id = 1, grid%nlmax_node
          if (patch%imat(grid%nL2G(local_id)) == geomech_material_id) then
            vec_p(local_id) = work_p(local_id)
          endif
        enddo
      endif
      call VecRestoreArrayF90(temp_vec,work_p,ierr);CHKERRQ(ierr)
    class default
        option%io_buffer = 'Dataset "' // trim(dataset%name) // '" is of the &
          &wrong type for GeomechReadDatasetToVecWithMask()'
        call PrintErrMsg(option)
    end select
  endif
  call VecRestoreArrayF90(vec,vec_p,ierr);CHKERRQ(ierr)

  end subroutine GeomechReadDatasetToVecWithMask

! ************************************************************************** !

end module Init_Subsurface_Geomech_module
