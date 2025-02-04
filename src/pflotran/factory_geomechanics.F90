
module Factory_Geomechanics_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes

  use Simulation_Geomechanics_class
  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: FactoryGeomechanicsInitialize, &
            FactoryGeomechReadSimBlock

contains

! ************************************************************************** !

subroutine FactoryGeomechanicsInitialize(simulation)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  !

  implicit none

  class(simulation_geomechanics_type) :: simulation

  ! NOTE: PETSc must already have been initialized here!
  call GeomechanicsInitializePostPETSc(simulation)

end subroutine FactoryGeomechanicsInitialize

! ************************************************************************** !

subroutine GeomechanicsInitializePostPETSc(simulation)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL and Satish Karra, LANL
  ! Date: 01/01/14, 02/10/15
  !
  use Simulation_Subsurface_class
  use Factory_Subsurface_module
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
  use Simulation_Base_class

  implicit none

  class(simulation_geomechanics_type) :: simulation

  type(option_type), pointer :: option
  class(realization_subsurface_type), pointer :: subsurf_realization
  class(realization_geomech_type), pointer :: geomech_realization
  class(pmc_base_type), pointer :: cur_process_model_coupler
  type(gmdm_ptr_type), pointer :: dm_ptr
  class(pm_base_type), pointer :: cur_pm, prev_pm
  class(pm_geomech_force_type), pointer :: pm_geomech
  class(pmc_geomechanics_type), pointer :: pmc_geomech
  class(timestepper_steady_type), pointer :: timestepper
  character(len=MAXSTRINGLENGTH) :: string
  type(input_type), pointer :: input
  PetscErrorCode :: ierr

  option => simulation%option
  nullify(timestepper)

  nullify(prev_pm)
  cur_pm => simulation%process_model_list
  do
    if (.not.associated(cur_pm)) exit
    select type(cur_pm)
      class is(pm_geomech_force_type)
        pm_geomech => cur_pm
        if (associated(prev_pm)) then
          prev_pm%next => cur_pm%next
        else
          simulation%process_model_list => cur_pm%next
        endif
        exit
      class default
    end select
    prev_pm => cur_pm
    cur_pm => cur_pm%next
  enddo

  call InitSubsurfGeomechSetGeomechMode(pm_geomech,option)

  call FactorySubsurfaceInitPostPetsc(simulation)
  simulation%process_model_coupler_list%is_master = PETSC_TRUE

  subsurf_realization => simulation%realization
  subsurf_realization%output_option => OutputOptionDuplicate(simulation%output_option)

  geomech_realization => GeomechRealizCreate(option)
  simulation%geomech%realization => geomech_realization
  input => InputCreate(IN_UNIT,option%input_filename,option)
  call InitSubsurfGeomechReadRequiredCards(geomech_realization,input)

  pmc_geomech => PMCGeomechanicsCreate()
  simulation%geomech%process_model_coupler => pmc_geomech
  pmc_geomech%name = 'PMCGeomech'
  pmc_geomech%option => option
  pmc_geomech%waypoint_list => simulation%waypoint_list_subsurface
  pmc_geomech%pm_list => pm_geomech
  pmc_geomech%pm_ptr%pm => pm_geomech
  pmc_geomech%geomech_realization => geomech_realization
  pmc_geomech%subsurf_realization => simulation%realization
  pm_geomech%subsurf_realization => simulation%realization

  ! add time integrator
  timestepper => TimestepperSteadyCreate()
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
  nullify(geomech_realization%output_option%output_snap_variable_list)
  nullify(geomech_realization%output_option%output_obs_variable_list)
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
  if (associated(simulation%geomech%process_model_coupler)) then
    call simulation%geomech%process_model_coupler% &
           SetWaypointPtr(pmc_geomech%waypoint_list)
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

  ! initialize geomech realization
  call InitSubsurfGeomechSetupRealization(subsurf_realization, &
                                          geomech_realization)

  call pm_geomech%PMGeomechForceSetRealization(geomech_realization)
  call pm_geomech%Setup()

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

  call GeomechanicsRegressionCreateMapping(simulation%geomech%regression, &
                                           geomech_realization)

  ! sim_aux: Set pointer
  simulation%flow_process_model_coupler%sim_aux => simulation%sim_aux
  if (associated(simulation%tran_process_model_coupler)) &
    simulation%tran_process_model_coupler%sim_aux => simulation%sim_aux
  if (option%ngeomechdof>0 .and. &
     associated(simulation%geomech%process_model_coupler)) &
    simulation%geomech%process_model_coupler%sim_aux => simulation%sim_aux

  ! set geomech as not master
  simulation%geomech%process_model_coupler%is_master = PETSC_FALSE
  ! link geomech and master
  simulation%process_model_coupler_list => &
    simulation%geomech%process_model_coupler
  ! link subsurface flow as peer
  simulation%process_model_coupler_list%peer => &
    simulation%flow_process_model_coupler

  call InitSubsurfGeomechChkInactiveCells(geomech_realization, &
                                          subsurf_realization)

  ! Set data in sim_aux
  cur_process_model_coupler => simulation%process_model_coupler_list
  call cur_process_model_coupler%SetAuxData()
  if (associated(cur_process_model_coupler%peer)) then
    cur_process_model_coupler => cur_process_model_coupler%peer
    call cur_process_model_coupler%GetAuxData()
    call cur_process_model_coupler%SetAuxData()
    select type(pmc => cur_process_model_coupler)
      class is(pmc_geomechanics_type)
        call GeomechStoreInitialPressTemp(pmc%geomech_realization)
    end select
  endif

  call InitSubsurfGeomechJumpStart(simulation%geomech)
  call InputDestroy(input)

end subroutine GeomechanicsInitializePostPETSc

! ************************************************************************** !

subroutine FactoryGeomechReadSimBlock(input,pm)
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/25/21
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

end subroutine FactoryGeomechReadSimBlock

! ************************************************************************** !

end module Factory_Geomechanics_module
