module Init_Subsurface_Geomech_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: InitSubsurfGeomechReadRequiredCards, &
            InitSubsurfGeomechReadInput 

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
  option%ngeomechdof = 3  ! displacements in x, y, z directions
  option%n_stress_strain_dof = 6

  string = "GEOMECHANICS_GRID"
  call InputFindStringInFile(input,option,string)
  call GeomechanicsInit(geomech_realization,input,option)


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
  ! comments (jaa) moved from factory_geomechanics

  !use Simulation_Geomechanics_class
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
  use Utility_module, only : DeallocateArray, UtilityReadArray
  use Geomechanics_Attr_module

  ! Still need to add other geomech modules for output, etc once created

  implicit none

  !class(simulation_geomechanics_type) :: simulation
  type(solver_type), pointer :: geomech_solver
  type(input_type), pointer :: input
  type(geomechanics_attr_type), pointer:: geomech
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option

  class(realization_geomech_type), pointer :: geomech_realization
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

      !.........................................................................
      case('NEWTON_SOLVER')
        call InputReadCard(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('GEOMECHANICS')
            call SolverReadNewton(geomech_solver,input,option)
        end select

     !....................
      case ('LINEAR_SOLVER')
        call InputReadCard(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case(word)
          case('GEOMECHANICS')
            call SolverReadLinear(geomech_solver,input,option)
        end select

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

subroutine GeomechanicsInit(geomech_realization,input,option)
  !
  ! Reads the required geomechanics data from input file
  !
  ! Author: Satish Karra, LANL
  ! Date: 05/23/13
  !

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
  type(geomech_discretization_type), pointer :: geomech_discretization
  type(geomech_patch_type), pointer :: patch
  type(input_type), pointer :: input
  type(option_type), pointer :: option
  character(len=MAXWORDLENGTH) :: word
  type(grid_unstructured_type), pointer :: ugrid
  character(len=MAXWORDLENGTH) :: card

  geomech_discretization       => geomech_realization%geomech_discretization

  input%ierr = 0
  ! we initialize the word to blanks to avoid error reported by valgrind
  word = ''

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

            geomech_discretization%grid  => GMGridCreate()
            ugrid => UGridCreate()
            call UGridRead(ugrid,geomech_discretization%filename,option)
            call UGridDecompose(ugrid,option)
            call CopySubsurfaceGridtoGeomechGrid(ugrid, &
                                                 geomech_discretization%grid, &
                                                 option)
            patch => GeomechanicsPatchCreate()
            patch%geomech_grid => geomech_discretization%grid
            geomech_realization%geomech_patch => patch
          case default
            option%io_buffer = 'Geomechanics supports only unstructured grid'
            call PrintErrMsg(option)
        end select
      case ('GRAVITY')
        call InputReadDouble(input,option,option%geomech_gravity(X_DIRECTION))
        call InputErrorMsg(input,option,'x-direction','GEOMECH GRAVITY')
        call InputReadDouble(input,option,option%geomech_gravity(Y_DIRECTION))
        call InputErrorMsg(input,option,'y-direction','GEOMECH GRAVITY')
        call InputReadDouble(input,option,option%geomech_gravity(Z_DIRECTION))
        call InputErrorMsg(input,option,'z-direction','GEOMECH GRAVITY')
        if (OptionIsIORank(option) .and. OptionPrintToScreen(option)) &
            write(option%fid_out,'(/," *GEOMECH_GRAV",/, &
            & "  gravity    = "," [m/s^2]",3x,1p3e12.4 &
            & )') option%geomech_gravity(1:3)
      case default
        call InputKeywordUnrecognized(input,word,'GEOMECHANICS_GRID',option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine GeomechanicsInit

! ************************************************************************** !

end module Init_Subsurface_Geomech_module
