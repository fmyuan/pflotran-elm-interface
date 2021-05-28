module Factory_PFLOTRAN_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: FactoryPFLOTRANInitPrePetsc, &
            FactoryPFLOTRANInitPostPetsc, &
            FactoryPFLOTRANFinalize

contains

! ************************************************************************** !

subroutine FactoryPFLOTRANInitPrePetsc(option)
!
! Sets up PFLOTRAN subsurface simulation framework prior to PETSc
!   initialization
! Author: Glenn Hammond
! Date: 06/07/13
!
  use Option_module
  use Input_Aux_module
  use HDF5_Aux_module
  use Logging_module

  implicit none

  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: bool_flag
  PetscBool :: option_found
  PetscErrorCode :: ierr

  ! NOTE: Cannot add anything that requires PETSc in this routine as PETSc
  !       has not yet been initialized.

  call FactoryPFLOTRANReadCommandLine(option)

  if (option%verbosity > 0) then
    call PetscLogDefaultBegin(ierr);CHKERRQ(ierr)
    string = '-log_view'
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, &
                                  string, ierr);CHKERRQ(ierr)
  endif

  call OptionBeginTiming(option)

end subroutine FactoryPFLOTRANInitPrePetsc

! ************************************************************************** !

subroutine FactoryPFLOTRANInitPostPetsc(simulation,driver,option)
!
! Sets up PFLOTRAN subsurface simulation framework after PETSc initialization
! Author: Glenn Hammond
! Date: 06/17/13
!
  use Driver_module
  use Option_module
  use Simulation_Base_class
  use Output_Aux_module
  use Logging_module
  use EOS_module
  use Input_Aux_module

  implicit none

  class(simulation_base_type), pointer :: simulation
  class(driver_type), pointer :: driver
  type(option_type), pointer :: option

  character(len=MAXSTRINGLENGTH) :: filename
  PetscErrorCode :: ierr

  ! popped in SimulationBaseInitializeRun()
  call PetscLogStagePush(logging%stage(INIT_STAGE),ierr);CHKERRQ(ierr)
  call PetscLogEventBegin(logging%event_init,ierr);CHKERRQ(ierr)

  call EOSInit()
  filename = trim(option%global_prefix) // trim(option%group_prefix) // &
             '.out'
  if (option%myrank == option%io_rank .and. option%print_to_file) then
    open(option%fid_out, file=filename, action="write", status="unknown")
  endif

  call OptionPrintPFLOTRANHeader(option)
  call FactoryPFLOTRANReadSimulationBlk(simulation,driver,option)
  call InputCheckKeywordBlockCount(option)
  ! Must come after simulation is initialized so that proper stages are setup
  ! for process models.  This call sets flag that disables the creation of
  ! new stages, which is necessary for multisimulation
  call LoggingSetupComplete()

end subroutine FactoryPFLOTRANInitPostPetsc

! ************************************************************************** !

subroutine FactoryPFLOTRANReadSimulationBlk(simulation,driver,option)
!
! Sets up PFLOTRAN subsurface simulation framework after PETSc initialization
! Author: Glenn Hammond
! Date: 06/17/13
!
  use Driver_module
  use Option_module
  use Input_Aux_module
  use String_module

  use Simulation_Base_class
  use Simulation_Subsurface_class
  use Simulation_Geomechanics_class
  use PM_Base_class
  use PMC_Base_class
  use Checkpoint_module
  use Output_Aux_module
  use Waypoint_module
  use Units_module

  use Factory_Subsurface_module
  use Factory_Geomechanics_module

  implicit none

  class(simulation_base_type), pointer :: simulation
  class(driver_type), pointer :: driver
  type(option_type), pointer :: option

  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: simulation_type

  class(pm_base_type), pointer :: pm_master
  class(pm_base_type), pointer :: cur_pm
  type(checkpoint_option_type), pointer :: checkpoint_option
  type(waypoint_list_type), pointer :: checkpoint_waypoint_list

  class(pmc_base_type), pointer :: pmc_master

  PetscBool :: print_ekg

  nullify(pm_master)
  nullify(cur_pm)

  nullify(pmc_master)
  nullify(checkpoint_option)
  nullify(checkpoint_waypoint_list)
  print_ekg = PETSC_FALSE

  input => InputCreate(IN_UNIT,option%input_filename,option)

  simulation_type = ''
  string = 'SIMULATION'
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)
  word = ''
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'PROCESS_MODEL','SIMULATION')

    call StringToUpper(word)
    select case(trim(word))
      case('SIMULATION_TYPE')
          call InputReadCard(input,option,simulation_type,PETSC_FALSE)
          call InputErrorMsg(input,option,'simulation_type', &
                             'SIMULATION')
      case('PROCESS_MODELS')
        call FactoryPFLOTRANReadSimProcessModels(input,pm_master,option)
      case('MASTER')
        call FactoryPFLOTRANSetupPMCHierarchy(input,option,pmc_master)
      case('PRINT_EKG')
        option%print_ekg = PETSC_TRUE
      case('CHECKPOINT')
        checkpoint_option => CheckpointOptionCreate()
        checkpoint_waypoint_list => WaypointListCreate()
        call CheckpointRead(input,option,checkpoint_option, &
                            checkpoint_waypoint_list)
      case ('RESTART')
        call FactoryPFLOTRANReadRestart(input,option)
      case('INPUT_RECORD_FILE')
        option%input_record = PETSC_TRUE
        call OpenAndWriteInputRecord(option)
      case default
        call InputKeywordUnrecognized(input,word,'SIMULATION',option)
    end select
  enddo
  call InputPopBlock(input,option)
  call InputDestroy(input)

  if (.not.associated(pm_master)) then
    option%io_buffer = 'No process models defined in SIMULATION block.'
    call PrintErrMsg(option)
  endif

  if (option%print_ekg) then
    cur_pm => pm_master
    do
      if (.not.associated(cur_pm)) exit
      cur_pm%print_ekg = PETSC_TRUE
      cur_pm => cur_pm%next
    enddo
  endif

  if (.not.associated(simulation)) then
    ! create the simulation objects
    select case(simulation_type)
      case('SUBSURFACE')
        simulation => SimSubsurfCreate(driver,option)
      case('GEOMECHANICS_SUBSURFACE')
        simulation => GeomechanicsSimulationCreate(driver,option)
      case default
        if (len_trim(simulation_type) == 0) then
          option%io_buffer = 'A SIMULATION_TYPE (e.g. "SIMULATION_TYPE &
            &SUBSURFACE") must be specified within the SIMULATION block.'
          call PrintErrMsg(option)
        endif
        call InputKeywordUnrecognized(input,simulation_type, &
                      'SIMULATION,SIMULATION_TYPE',option)
    end select
  endif

  select type(simulation)
    class is(simulation_subsurface_type)
      call WaypointListMerge(simulation%waypoint_list_outer, &
                             checkpoint_waypoint_list,option)
  end select

  select type(simulation)
    class is(simulation_subsurface_type)
      simulation%process_model_list => pm_master
      simulation%checkpoint_option => checkpoint_option
  end select

  select type(simulation)
    class is(simulation_subsurface_type)
      call FactorySubsurfaceInitialize(simulation)
    class is(simulation_geomechanics_type)
      call FactoryGeomechanicsInitialize(simulation)
  end select

end subroutine FactoryPFLOTRANReadSimulationBlk

! ************************************************************************** !

subroutine FactoryPFLOTRANReadSimProcessModels(input,pm_master,option)
!
! Reads in the process models listed in simulation block
!
! Author: Glenn Hammond
! Date: 02/05/21
!
  use Option_module
  use Input_Aux_module
  use String_module

  use PM_Base_class
  use PM_Geomechanics_Force_class
  use PM_Auxiliary_class

  use Factory_Subsurface_module
  use Factory_Geomechanics_module

  implicit none

  class(pm_base_type), pointer :: pm_master
  type(input_type), pointer :: input
  type(option_type), pointer :: option

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: pm_name
  class(pm_base_type), pointer :: cur_pm
  class(pm_base_type), pointer :: new_pm

  nullify(cur_pm)
  nullify(new_pm)

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'process_model', &
                       'SIMULATION,PROCESS_MODELS')
    call InputReadWord(input,option,pm_name,PETSC_TRUE)
    if (InputError(input)) then
      input%err_buf = 'Process Model Name'
      call InputDefaultMsg(input,option)
      pm_name = ''
    endif
    call StringToUpper(word)
    select case(trim(word))
      case('SUBSURFACE_FLOW')
        call FactorySubsurfaceReadFlowPM(input,option,new_pm)
      case('SUBSURFACE_TRANSPORT')
        call FactorySubsurfaceReadTransportPM(input,option,new_pm)
      case('NUCLEAR_WASTE_TRANSPORT')
        if (OptionPrintToScreen(option)) then
          print *
          print *, 'SIMULATION'
          print *, '  SIMULATION_TYPE SUBSURFACE'
          print *, '  PROCESS_MODELS'
          print *, '    SUBSURFACE_TRANSPORT'
          print *, '      MODE NWT'
          print *, '      OPTIONS'
          print *, '      /'
          print *, '    /'
          print *, '  /'
          print *, 'END'
          print *
        endif
        option%io_buffer = "PFLOTRAN's NUCLEAR_WASTE_TRANSPORT &
          &process model has been refactored to use the &
          &combination of the SUBSURFACE_TRANSPORT and 'MODE &
          &NWT' keywords and an (optional) OPTIONS block. &
          &Please use the keywords above in reformatting the &
          &SIMULATION block."
        call PrintErrMsg(option)
      case('WASTE_FORM')
        call FactorySubsurfaceReadWasteFormPM(input,option,new_pm)
      case('UFD_DECAY')
        call FactorySubsurfaceReadUFDDecayPM(input,option,new_pm)
      case('UFD_BIOSPHERE')
        call FactorySubsurfReadUFDBiospherePM(input,option,new_pm)
      case('WIPP_SOURCE_SINK')
        option%io_buffer = 'Do not include the WIPP_SOURCE_SINK block &
          &unless you are running in WIPP_FLOW mode and intend to &
          &include gas generation.'
        call PrintErrMsg(option)
      case('GEOMECHANICS_SUBSURFACE')
        option%geomech_on = PETSC_TRUE
        new_pm => PMGeomechForceCreate()
      case('SUBSURFACE_GEOPHYSICS')
        call FactorySubsurfReadGeophysicsPM(input,option,new_pm)
      case('AUXILIARY')
        if (len_trim(pm_name) < 1) then
          option%io_buffer = 'AUXILIARY process models must have a name.'
          call PrintErrMsg(option)
        endif
        new_pm => PMAuxiliaryCreate()
        input%buf = pm_name
        call PMAuxiliaryRead(input,option,PMAuxiliaryCast(new_pm))
      case default
        call InputKeywordUnrecognized(input,word, &
               'SIMULATION,PROCESS_MODELS',option)
    end select
    if (.not.associated(new_pm%option)) new_pm%option => option
    if (len_trim(pm_name) > 0) then
      new_pm%name = pm_name
    endif
    if (associated(cur_pm)) then
      cur_pm%next => new_pm
    else
      cur_pm => new_pm
    endif
    if (.not.associated(pm_master)) then
      pm_master => new_pm
    endif
    cur_pm => new_pm
    nullify(new_pm)
  enddo
  call InputPopBlock(input,option)

end subroutine FactoryPFLOTRANReadSimProcessModels

! ************************************************************************** !

subroutine FactoryPFLOTRANReadRestart(input,option)
!
! Read the restart block within the simulation block
! Author: Glenn Hammond
! Date: 02/05/21
!
  use Option_module
  use Input_Aux_module
  use PMC_Base_class
  use String_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: word
  PetscBool :: realization_dependent_restart
  PetscInt :: i

  option%restart_flag = PETSC_TRUE
  realization_dependent_restart = PETSC_FALSE
  ! this section preserves the legacy implementation
  call InputReadFilename(input,option,option%restart_filename)
  if (input%ierr == 0) then
    call InputReadWord(input,option,word,PETSC_TRUE)
    if (input%ierr == 0) then
      option%restart_time = 0.d0
    endif
    return
    ! end legacy implementation
  endif
  input%ierr = 0
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call StringToUpper(word)
    select case(word)
      case('FILENAME')
        call InputReadFilename(input,option,option%restart_filename)
        call InputErrorMsg(input,option,'RESTART','filename')
      case('RESET_TO_TIME_ZERO')
        ! any value but UNINITIALIZED_DOUBLE will set back to zero.
        option%restart_time = 0.d0
      case('REALIZATION_DEPENDENT')
        realization_dependent_restart = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(input,word, &
                                      'SIMULATION,RESTART',option)
    end select
  enddo
  call InputPopBlock(input,option)
  if (realization_dependent_restart) then
    ! insert realization id
    i = index(option%restart_filename,'-restart')
    if (i > 0) then
      string = option%restart_filename(1:i-1) // &
               trim(option%group_prefix) // &
             option%restart_filename(i:len_trim(option%restart_filename))
    else
      option%io_buffer = 'Realization-dependent restart requires that &
        &"-restart" be present in the restart file name so that the &
        &realization id can be inserted prior to -restart.  E.g. &
        &pflotran-restart.h5 -> pflotranR1-restart.h5'
      call PrintErrMsg(option)
    endif
    option%restart_filename = trim(string)
  endif

end subroutine FactoryPFLOTRANReadRestart

! ************************************************************************** !

recursive subroutine FactoryPFLOTRANSetupPMCHierarchy(input,option,pmc)
!
! Forms a linked list of named dummy pmcs as placeholders
! Author: Glenn Hammond
! Date: 12/10/14
!
  use Option_module
  use Input_Aux_module
  use PMC_Base_class
  use String_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option
  class(pmc_base_type), pointer :: pmc

  character(len=MAXWORDLENGTH) :: word

  call InputReadWord(input,option,word,PETSC_TRUE)
  call InputErrorMsg(input,option,'PMC name','SIMULATION')
    ! at this point, we are creating a
  pmc => PMCBaseCreate()
  pmc%name = word

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'CHILD or PEER','SIMULATION')
    call StringToUpper(word)
    select case(trim(word))
      case('PEER')
        call FactoryPFLOTRANSetupPMCHierarchy(input,option,pmc%peer)
      case('CHILD')
        call FactoryPFLOTRANSetupPMCHierarchy(input,option,pmc%child)
      case default
        call InputKeywordUnrecognized(input,word, &
                                      'FactoryPFLOTRANSetupPMCHierarchy',option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine FactoryPFLOTRANSetupPMCHierarchy

! ************************************************************************** !

recursive subroutine FactoryPFLOTRANLinkPMToPMC(input,option,pmc,pm)
!
! Forms a linked list of named dummy pmcs as placeholders
! Author: Glenn Hammond
! Date: 12/10/14
!
  use Option_module
  use Input_Aux_module
  use String_module
  use PM_Base_class
  use PMC_Base_class

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option
  class(pmc_base_type), pointer :: pmc
  class(pm_base_type), pointer :: pm

  if (.not.associated(pmc)) return

  print *, pmc%name, pm%name
  if (StringCompareIgnoreCase(pmc%name,pm%name)) then
    pmc%pm_list => pm
    return
  endif

  call FactoryPFLOTRANLinkPMToPMC(input,option,pmc%peer,pm)
  call FactoryPFLOTRANLinkPMToPMC(input,option,pmc%child,pm)

end subroutine FactoryPFLOTRANLinkPMToPMC

! ************************************************************************** !

subroutine FactoryPFLOTRANFinalize(option)
!
! Destroys PFLOTRAN subsurface simulation framework
! Author: Glenn Hammond
! Date: 06/07/13
!
  use Option_module
  use HDF5_Aux_module

  implicit none

  type(option_type) :: option
  PetscErrorCode :: ierr

  call MPI_Barrier(option%driver%comm%global_comm,ierr)

end subroutine FactoryPFLOTRANFinalize

! ************************************************************************** !

subroutine FactoryPFLOTRANReadCommandLine(option)
  !
  ! Initializes PFLOTRAN output filenames, etc.
  !
  ! Author: Glenn Hammond
  ! Date: 06/06/13
  !

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string, string2
  PetscBool :: option_found
  PetscBool :: bool_flag
  PetscBool :: pflotranin_option_found
  PetscBool :: input_prefix_option_found
  PetscBool :: output_dir_found
  PetscBool :: output_file_prefix_found
  character(len=MAXSTRINGLENGTH), pointer :: strings(:)
  PetscInt :: i
  PetscErrorCode :: ierr

  ! check for non-default input filename
  option%input_filename = 'pflotran.in'
  string = '-pflotranin'
  call InputGetCommandLineString(string,option%input_filename, &
                                 pflotranin_option_found,option)
  string = '-input_prefix'
  call InputGetCommandLineString(string,option%input_prefix, &
                                 input_prefix_option_found,option)

  if (pflotranin_option_found .and. input_prefix_option_found) then
    option%io_buffer = 'Cannot specify both "-pflotranin" and ' // &
      '"-input_prefix" on the command lines.'
    call PrintErrMsg(option)
  else if (pflotranin_option_found) then
    strings => StringSplit(option%input_filename,'.')
    if (size(strings) > 1) then
      option%input_prefix = StringsMerge(strings(1:size(strings)-1),'.')
    else
      option%input_prefix = strings(1)
    endif
    deallocate(strings)
    nullify(strings)
  else if (input_prefix_option_found) then
    option%input_filename = trim(option%input_prefix) // '.in'
  endif

  string = '-output_prefix'
  call InputGetCommandLineString(string,option%global_prefix,option_found,option)

  if (.not.option_found) then
    option%global_prefix = option%input_prefix
    option%output_file_name_prefix = option%input_prefix
  else
    call InputReadFileDirNamePrefix(option%global_prefix, &
                                    option%output_file_name_prefix, &
                                    option%output_dir)
  end if

  string = '-screen_output'
  call InputGetCommandLineTruth(string,option%print_to_screen,option_found,option)

  string = '-file_output'
  call InputGetCommandLineTruth(string,option%print_to_file,option_found,option)

  string = '-v'
  call InputGetCommandLineInt(string,i,option_found,option)
  if (option_found) option%verbosity = i

  string = '-successful_exit_code'
  call InputGetCommandLineInt(string,i,option_found,option)
  if (option_found) option%exit_code = i

  string = '-keyword_screen_output'
  call InputGetCommandLineTruth(string,option%keyword_logging_screen_output, &
                                option_found,option)

  ! this will get overwritten later if stochastic
  string = '-realization_id'
  call InputGetCommandLineInt(string,i,option_found,option)
  if (option_found) then
    if (i < 1) then
      option%io_buffer = 'realization_id must be greater than zero.'
      call PrintErrMsg(option)
    endif
    option%id = i
  endif

end subroutine FactoryPFLOTRANReadCommandLine

end module Factory_PFLOTRAN_module
