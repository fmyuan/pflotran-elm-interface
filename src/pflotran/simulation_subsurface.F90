module Simulation_Subsurface_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Simulation_Base_class
  use Simulation_Aux_module
  use Option_module
  use Output_Aux_module
  use PM_Base_class
  use PMC_Base_class
  use PMC_Geophysics_class
  use PMC_Subsurface_class
  use Realization_Subsurface_class
  use Waypoint_module
  use Regression_module

  implicit none

  private

  type, public, extends(simulation_base_type) :: simulation_subsurface_type
    type(option_type), pointer :: option
    type(output_option_type), pointer :: output_option
    PetscInt :: stop_flag
    class(pmc_base_type), pointer :: process_model_coupler_list
    class(pm_base_type), pointer :: process_model_list
    type(simulation_aux_type), pointer :: sim_aux
    ! pointer to flow process model coupler
    class(pmc_subsurface_type), pointer :: flow_process_model_coupler
    ! pointer to transport process model coupler
    class(pmc_subsurface_type), pointer :: tran_process_model_coupler
    ! pointer to geophysics process model coupler
    class(pmc_geophysics_type), pointer :: geop_process_model_coupler
    ! pointer to realization object shared by flow and reactive transport
    class(realization_subsurface_type), pointer :: realization
    ! regression object
    type(regression_type), pointer :: regression
    type(waypoint_list_type), pointer :: waypoint_list_subsurface
    type(waypoint_list_type), pointer :: waypoint_list_outer ! outer sync loop
  contains
    procedure, public :: JumpStart => SimSubsurfJumpStart
    procedure, public :: InitializeRun => SimSubsurfInitializeRun
    procedure, public :: InputRecord => SimSubsurfInputRecord
    procedure, public :: ExecuteRun => SimSubsurfExecuteRun
    procedure, public :: RunToTime => SimSubsurfRunToTime
    procedure, public :: FinalizeRun => SimSubsurfFinalizeRun
    procedure, public :: Strip => SimSubsurfStrip
  end type simulation_subsurface_type


  public :: SimSubsurfCreate, &
            SimSubsurfInit, &
            SimSubsurfCast, &
            SimSubsurfInitializeRun, &
            SimSubsurfFinalizeRun, &
            SimSubsurfStrip, &
            SimSubsurfDestroy

  public :: SimSubsurfGetFinalWaypointTime

contains

! ************************************************************************** !

function SimSubsurfCreate(driver,option)
  !
  ! Allocates and initializes a new simulation object
  !
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  !
  use Driver_module
  use Option_module

  implicit none

  class(driver_type), pointer :: driver
  type(option_type), pointer :: option

  class(simulation_subsurface_type), pointer :: SimSubsurfCreate

#ifdef DEBUG
  print *, 'SimSubsurfCreate'
#endif

  allocate(SimSubsurfCreate)
  call SimSubsurfInit(SimSubsurfCreate,driver,option)

end function SimSubsurfCreate

! ************************************************************************** !

subroutine SimSubsurfInit(this,driver,option)
  !
  ! Initializes simulation values
  !
  ! Author: Glenn Hammond
  ! Date: 04/22/13
  !
  use Timestepper_Base_class, only : TS_CONTINUE
  use Output_Aux_module
  use Waypoint_module
  use Driver_module
  use Option_module

  implicit none

  class(simulation_subsurface_type) :: this
  class(driver_type), pointer :: driver
  type(option_type), pointer :: option

#ifdef DEBUG
  call PrintMsg(this%option,'SimSubsurfInit()')
#endif

  call SimulationBaseInit(this,driver)
  this%option => option
  this%output_option => OutputOptionCreate()
  nullify(this%process_model_coupler_list)
  nullify(this%process_model_list)
  this%sim_aux => SimAuxCreate()
  this%stop_flag = TS_CONTINUE
  nullify(this%flow_process_model_coupler)
  nullify(this%tran_process_model_coupler)
  nullify(this%geop_process_model_coupler)
  nullify(this%realization)
  nullify(this%regression)
  this%waypoint_list_subsurface => WaypointListCreate()
  this%waypoint_list_outer => WaypointListCreate()

end subroutine SimSubsurfInit

! ************************************************************************** !

function SimSubsurfCast(simulation)
  !
  ! Casts any simulation_type to simulation_subsurface_type if of that type
  !
  ! Author: Glenn Hammond
  ! Date: 02/23/21
  !
  implicit none

  class(simulation_base_type), pointer :: simulation

  class(simulation_subsurface_type), pointer :: SimSubsurfCast

  nullify(SimSubsurfCast)
  select type(simulation)
    class is(simulation_subsurface_type)
      SimSubsurfCast=> simulation
  end select

end function SimSubsurfCast

! ************************************************************************** !

subroutine SimSubsurfInitializeRun(this)
  !
  ! Initializes simulation
  !
  ! Author: Glenn Hammond
  ! Date: 02/15/21
  !

  use Logging_module
  use Option_module
  use Option_Checkpoint_module
  use Output_module
  use hdf5

  implicit none


  class(simulation_subsurface_type) :: this

  integer(HID_T) :: h5_chk_grp_id
  PetscViewer :: viewer
  PetscBool :: flag
  PetscErrorCode :: ierr

#ifdef DEBUG
  call PrintMsg(this%option,'SimSubsurfInitializeRun()')
#endif

  call SimulationBaseInitializeRun(this)
  ! print error message if binary checkpoint/restart is used in
  ! combination with unstructured gridding
  flag = PETSC_FALSE
  if (associated(this%option%checkpoint)) then
    if (this%option%checkpoint%format == CHECKPOINT_BINARY) then
      flag = PETSC_TRUE
    endif
  endif
  if (index(this%option%restart_filename,'.chk') > 0) then
    flag = PETSC_TRUE
  endif
  if (flag) then
    flag = PETSC_FALSE
    select type(s=>this)
      class is(simulation_subsurface_type)
        ! also covers simulation_geomechanics_type
        if (s%realization%patch%grid%itype /= STRUCTURED_GRID) then
          flag = PETSC_TRUE
        endif
      class default
        this%option%io_buffer = 'Unknown simulation class in &
          &SimSubsurfInitializeRun'
        call PrintErrMsg(this%option)
    end select
    if (flag) then
        this%option%io_buffer = 'Binary Checkpoint/Restart (.chk format) &
          &is not supported for unstructured grids.  Please use HDF5 &
          &(.h5 format).'
        call PrintErrMsg(this%option)
    endif
  endif

  ! the user may request output of variable that do not exist for the
  ! the requested process models; this routine should catch such issues.
  call OutputEnsureVariablesExist(this%output_option,this%option)
  call SimSubsurfForbiddenCombinations(this)

  if (associated(this%process_model_coupler_list)) then
    if (this%option%restart_flag) then
      if (index(this%option%restart_filename,'.chk') > 0) then
        call this%process_model_coupler_list%RestartBinary(viewer)
      elseif (index(this%option%restart_filename,'.h5') > 0) then
        call this%process_model_coupler_list%RestartHDF5(h5_chk_grp_id)
      else
        this%option%io_buffer = 'Unknown restart filename format. ' // &
        'Only *.chk and *.h5 supported.'
        call PrintErrMsg(this%option)
      endif
    endif

    ! initialize performs overwrite of restart, if applicable
    call this%process_model_coupler_list%InitializeRun()
    call this%JumpStart()
  endif

  call SimulationBaseInputRecordPrint(this,this%option)
  call PrintMsg(this%option,'')
  call PrintMsg(this%option,' Finished Initialization')
  call PrintMsg(this%option,'')
  if (OptionPrintToFile(this%option)) then
    inquire(this%option%fid_inputrecord,opened=flag)
    if (flag) close(this%option%fid_inputrecord)
  endif
  call PetscLogEventEnd(logging%event_init,ierr);CHKERRQ(ierr)
  ! pushed in PFLOTRANInitializePostPetsc()
  call PetscLogStagePop(ierr);CHKERRQ(ierr)

  ! popped in FinalizeRun()
  call PetscLogStagePush(logging%stage(TS_STAGE),ierr);CHKERRQ(ierr)

end subroutine SimSubsurfInitializeRun

! ************************************************************************** !

subroutine SimSubsurfInputRecord(this)
  !
  ! Writes ingested information to the input record file.
  !
  ! Author: Jenn Frederick, SNL
  ! Date: 03/17/2016
  !
  use Option_module
  use Output_module
  use Checkpoint_module
  use Discretization_module
  use Reaction_Aux_module
  use Region_module
  use Strata_module
  use Material_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Thermal_module
  use Patch_module
  use Condition_module
  use EOS_module
  use Waypoint_module

  implicit none

  class(simulation_subsurface_type) :: this

  PetscInt :: id = INPUT_RECORD_UNIT

  ! print checkpoint information
  call CheckpointInputRecord(this%option%checkpoint,this%waypoint_list_outer)

  write(id,'(a)') ' '
  ! print process model coupler and process model information
  call this%process_model_coupler_list%InputRecord()

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') 'simulation type: '
  write(id,'(a)') 'subsurface'
  write(id,'(a29)',advance='no') 'flow mode: '
  select case(this%realization%option%iflowmode)
    case(MPH_MODE)
      write(id,'(a)') 'multi-phase'
    case(RICHARDS_MODE)
      write(id,'(a)') 'richards'
    case(ZFLOW_MODE)
      write(id,'(a)') 'zflow'
    case(PNF_MODE)
      write(id,'(a)') 'pore network flow'
    case(RICHARDS_TS_MODE)
      write(id,'(a)') 'richards_ts'
    case(G_MODE)
      write(id,'(a)') 'general'
    case(H_MODE)
      write(id,'(a)') 'hydrate'
    case(WF_MODE)
      write(id,'(a)') 'wipp flow'
    case(TH_MODE)
      write(id,'(a)') 'thermo-hydro'
    case(TH_TS_MODE)
      write(id,'(a)') 'thermo-hydro_ts'
  end select

  ! print time information
  call WaypointInputRecord(this%output_option,this%waypoint_list_subsurface)

  ! print output file information
  call OutputInputRecord(this%output_option,this%waypoint_list_subsurface)

  ! print grid/discretization information
  call DiscretizationInputRecord(this%realization%discretization)

  ! print region information
  call RegionInputRecord(this%realization%patch%region_list)

  ! print strata information
  call StrataInputRecord(this%realization%patch%strata_list)

  ! print material property information
  call MaterialPropInputRecord(this%realization%material_properties)

  ! print characteristic curves information
  call CharCurvesInputRecord(this%realization%patch%characteristic_curves)

  ! print thermal characteristic curve info
  call CharCurvesThermalInputRecord( &
       this%realization%patch%characteristic_curves_thermal)

  ! print chemistry and reactive transport information
  call ReactionInputRecord(this%realization%reaction)

  ! print coupler information (ICs, BCs, SSs)
  call PatchCouplerInputRecord(this%realization%patch)

  ! print flow and trans condition information
  call FlowCondInputRecord(this%realization%flow_conditions, &
                           this%realization%option)
  call TranCondInputRecord(this%realization%transport_conditions, &
                           this%realization%option)

  ! print equation of state (eos) information
  call EOSInputRecord()

end subroutine SimSubsurfInputRecord

! ************************************************************************** !

subroutine SimSubsurfJumpStart(this)
  !
  ! Initializes simulation
  !
  ! Author: Glenn Hammond
  ! Date: 08/11/14
  !
  use Logging_module
  use Output_module
  use Option_module
  use Output_Aux_module
  use String_module
  use Timestepper_Base_class

  implicit none

  class(simulation_subsurface_type) :: this

  class(timestepper_base_type), pointer :: master_timestepper
  class(timestepper_base_type), pointer :: flow_timestepper
  class(timestepper_base_type), pointer :: tran_timestepper
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  PetscBool :: snapshot_plot_flag, observation_plot_flag, massbal_plot_flag
  PetscErrorCode :: ierr
  PetscBool :: bypass_final_time_check

  bypass_final_time_check = PETSC_FALSE

  call PetscOptionsHasName(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER, &
                           "-bypass_final_time_check",bypass_final_time_check, &
                           ierr);CHKERRQ(ierr)
#ifdef DEBUG
  call PrintMsg(this%option,'SimSubsurfJumpStart()')
#endif

  nullify(master_timestepper)
  nullify(flow_timestepper)
  nullify(tran_timestepper)
  snapshot_plot_flag = PETSC_FALSE
  observation_plot_flag = PETSC_FALSE
  massbal_plot_flag = PETSC_FALSE

  option => this%option
  output_option => this%output_option

  ! first time stepper is master
  master_timestepper => this%process_model_coupler_list%timestepper
  if (associated(this%flow_process_model_coupler)) then
    flow_timestepper => this%flow_process_model_coupler%timestepper
  endif
  if (associated(this%tran_process_model_coupler)) then
    tran_timestepper => this%tran_process_model_coupler%timestepper
  endif

  !if TIMESTEPPER->MAX_STEPS < 0, print out solution composition only
  if (master_timestepper%max_time_step < 0) then
    call PrintMsg(option,'')
    write(option%io_buffer,*) master_timestepper%max_time_step
    option%io_buffer = 'The maximum # of time steps (' // &
      trim(StringWrite(master_timestepper%max_time_step)) // &
      '), specified by TIMESTEPPER->MAX_STEPS, has been met.  Stopping....'
    call PrintMsg(option)
    call PrintMsg(option,'')
    option%driver%status = DONE
    this%stop_flag = TS_STOP_MAX_TIME_STEP
    return
  endif

  ! print initial condition output if not a restarted sim
  call OutputInit(option,master_timestepper%steps)
  if (output_option%plot_number == 0 .and. &
      master_timestepper%max_time_step >= 0) then
    snapshot_plot_flag = output_option%print_initial_snap
    observation_plot_flag = output_option%print_initial_obs
    massbal_plot_flag = output_option%print_initial_massbal
    call Output(this%realization,snapshot_plot_flag,observation_plot_flag, &
                massbal_plot_flag)
  endif

  !if TIMESTEPPER->MAX_STEPS < 1, print out initial condition only
  if (master_timestepper%max_time_step < 1) then
    call PrintMsg(option,'')
    option%io_buffer = 'The maximum # of time steps (' // &
      trim(StringWrite(master_timestepper%max_time_step)) // &
      '), specified by TIMESTEPPER->MAX_STEPS, has been met.  Stopping....'
    call PrintMsg(option)
    call PrintMsg(option,'')
    option%driver%status = DONE
    this%stop_flag = TS_STOP_MAX_TIME_STEP
    return
  endif

  ! increment plot number so that 000 is always the initial condition,
  ! and nothing else
  if (output_option%plot_number == 0) output_option%plot_number = 1

  if (associated(flow_timestepper)) then
    if (.not.associated(flow_timestepper%cur_waypoint)) then
      if (.not. bypass_final_time_check) then
        option%io_buffer = &
          'Null flow waypoint list; final time likely equal to start time &
          &or final simulation time needs to be extended on a restart.'
        call PrintMsg(option)
        option%driver%status = FAIL
        return
      endif
    else
      flow_timestepper%dt_max = flow_timestepper%cur_waypoint%dt_max
    endif
  endif
  if (associated(tran_timestepper)) then
    if (.not.associated(tran_timestepper%cur_waypoint)) then
      if (.not. bypass_final_time_check) then
        option%io_buffer = &
          'Null transport waypoint list; final time likely equal to start &
          &or final simulation time needs to be extended on a restart.'
        call PrintMsg(option)
        option%driver%status = FAIL
        return
      endif
    else
      tran_timestepper%dt_max = tran_timestepper%cur_waypoint%dt_max
    endif
  endif

  if (associated(flow_timestepper)) &
    flow_timestepper%start_time_step = flow_timestepper%steps + 1
  if (associated(tran_timestepper)) &
    tran_timestepper%start_time_step = tran_timestepper%steps + 1

  if (this%realization%debug%print_regions) then
    call OutputPrintRegions(this%realization)
    if (this%realization%discretization%itype == UNSTRUCTURED_GRID .and. &
        this%realization%patch%grid%itype == IMPLICIT_UNSTRUCTURED_GRID) then
      call OutputPrintRegionsH5(this%realization)
    endif
  endif

  if (this%realization%debug%print_couplers) then
    call OutputPrintCouplers(this%realization,ZERO_INTEGER)
    if (this%realization%discretization%itype == UNSTRUCTURED_GRID .and. &
        this%realization%patch%grid%itype == IMPLICIT_UNSTRUCTURED_GRID) then
      call OutputPrintCouplersH5(this%realization,ZERO_INTEGER)
    endif
  endif

end subroutine SimSubsurfJumpStart

! ************************************************************************** !

subroutine SimSubsurfExecuteRun(this)
  !
  ! Execute a simulation
  !
  ! Author: Glenn Hammond
  ! Date: 02/15/21
  !
  use Waypoint_module
  use Timestepper_Base_class, only : TS_CONTINUE
  use Checkpoint_module

  implicit none

  class(simulation_subsurface_type) :: this

  PetscReal :: final_time
  type(waypoint_type), pointer :: cur_waypoint
  character(len=MAXSTRINGLENGTH) :: append_name

#ifdef DEBUG
  call PrintMsg(this%option,'SimSubsurfExecuteRun()')
#endif

  if (.not.associated(this%process_model_coupler_list)) then
    return
  endif

  final_time = SimSubsurfGetFinalWaypointTime(this)
  cur_waypoint => this%waypoint_list_outer%first
  if (cur_waypoint%print_checkpoint) then
    append_name = &
         CheckpointAppendNameAtTime(this%process_model_coupler_list% &
                                        option%time, &
                                    this%process_model_coupler_list%option)
    call this%process_model_coupler_list%Checkpoint(append_name)
  endif
  call WaypointSkipToTime(cur_waypoint,this%option%time)
  do
    if (this%stop_flag /= TS_CONTINUE) exit ! end simulation
    if (.not.associated(cur_waypoint)) exit
    call this%RunToTime(min(final_time,cur_waypoint%time))
    cur_waypoint => cur_waypoint%next
  enddo
  append_name = '-restart'
  if (associated(this%option%checkpoint)) then
    call this%process_model_coupler_list%Checkpoint(append_name)
  endif

end subroutine SimSubsurfExecuteRun

! ************************************************************************** !

subroutine SimSubsurfRunToTime(this,target_time)
  !
  ! Executes simulation
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  use Option_module
  use Simulation_Aux_module

  implicit none

  class(simulation_subsurface_type) :: this
  PetscReal :: target_time

#ifdef DEBUG
  call PrintMsg(this%option,'SimSubsurfRunToTime()')
#endif

  call this%process_model_coupler_list%RunToTime(target_time,this%stop_flag)

end subroutine SimSubsurfRunToTime


! ************************************************************************** !

function SimSubsurfGetFinalWaypointTime(this)
  !
  ! Returns the earliest final waypoint time from the top layer of process
  ! model couplers.
  !
  ! Author: Glenn Hammond
  ! Date: 06/12/13
  !
  use Waypoint_module

  implicit none

  class(simulation_subsurface_type) :: this

  PetscReal :: SimSubsurfGetFinalWaypointTime

  class(pmc_base_type), pointer :: cur_process_model_coupler
  PetscReal :: final_time

  SimSubsurfGetFinalWaypointTime = 0.d0

  cur_process_model_coupler => this%process_model_coupler_list
  do
    if (.not.associated(cur_process_model_coupler)) exit
    final_time = WaypointListGetFinalTime(cur_process_model_coupler% &
                                            waypoint_list)
    if (SimSubsurfGetFinalWaypointTime < 1.d-40 .or. &
        final_time < SimSubsurfGetFinalWaypointTime) then
      SimSubsurfGetFinalWaypointTime = final_time
    endif
    cur_process_model_coupler => cur_process_model_coupler%peer
  enddo

end function SimSubsurfGetFinalWaypointTime

! ************************************************************************** !

subroutine SimSubsurfForbiddenCombinations(this)
  !
  ! Throws error messages when forbidden combinations of processes/process
  ! models are requested.
  !
  ! Author: Glenn Hammond
  ! Date: 05/31/22

  use Strata_module

  implicit none

  class(simulation_subsurface_type) :: this

  ! cannot update porosity based on mineral volume fractions and evolve
  ! strata at the same time.
  if (associated(this%realization%reaction)) then
    if (StrataEvolves(this%realization%patch%strata_list) .and. &
        this%realization%reaction%update_porosity) then
      call PrintErrMsg(this%option,'Time dependent STRATA and the update of &
            &porosity based on mineral volume fractions cannot be used &
            &simultaneously.')
    endif
  endif

end subroutine SimSubsurfForbiddenCombinations

! ************************************************************************** !

subroutine SimSubsurfFinalizeRun(this)
  !
  ! Finalizes simulation
  !
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  !

  use Logging_module
  use Timestepper_Base_class
  use String_module, only : StringWrite
  use Reaction_Sandbox_module, only : RSandboxDestroy
  use SrcSink_Sandbox_module, only : SSSandboxDestroyList
  use WIPP_module, only : WIPPDestroy
  use Klinkenberg_module, only : KlinkenbergDestroy
  use CLM_Rxn_module, only : RCLMRxnDestroy
  use Output_EKG_module

  implicit none

  class(simulation_subsurface_type) :: this

  character(MAXSTRINGLENGTH) :: string
  class(timestepper_base_type), pointer :: flow_timestepper
  class(timestepper_base_type), pointer :: tran_timestepper
  PetscErrorCode :: ierr

  if (this%stop_flag /= TS_STOP_END_SIMULATION) then
    select case(this%stop_flag)
      case(TS_STOP_WALLCLOCK_EXCEEDED)
        string = 'Wallclock stop time exceeded.  Exiting!'
      case(TS_STOP_MAX_TIME_STEP)
        string = 'Maximum timestep exceeded.  Exiting!'
      case(TS_STOP_FAILURE)
        string = 'Simulation failed.  Exiting!'
        this%driver%exit_code = EXIT_FAILURE
      case default
        string = 'Simulation stopped for unknown reason (' // &
                trim(StringWrite(this%stop_flag)) // ').'
    end select
    if (OptionPrintToScreen(this%option)) write(*,'(/,a,/)') trim(string)
  endif

  ! pushed in InitializeRun()
  call PetscLogStagePop(ierr);CHKERRQ(ierr)
  ! popped in below
  call PetscLogStagePush(logging%stage(FINAL_STAGE),ierr);CHKERRQ(ierr)

  if (associated(this%process_model_coupler_list)) then
    call this%process_model_coupler_list%FinalizeRun()
  endif

  nullify(flow_timestepper)
  nullify(tran_timestepper)
  if (associated(this%flow_process_model_coupler)) then
    flow_timestepper => this%flow_process_model_coupler%timestepper
    call SSSandboxDestroyList()
    call WIPPDestroy()
    call KlinkenbergDestroy()
  endif
  if (associated(this%tran_process_model_coupler)) then
    tran_timestepper => this%tran_process_model_coupler%timestepper
    if (this%option%itranmode == RT_MODE) then
      call RSandboxDestroy()
      call RCLMRxnDestroy()
    endif
  endif

  select case(this%stop_flag)
    case(TS_STOP_END_SIMULATION,TS_STOP_MAX_TIME_STEP)
      call RegressionOutput(this%regression,this%realization, &
                            flow_timestepper,tran_timestepper)
  end select

  call SimulationBaseFinalizeRun(this)
  if (OptionIsIORank(this%option)) then
    call SimulationBaseWriteTimes(this,this%option%fid_out)
  endif

  ! close output files
  call PetscLogStagePop(ierr);CHKERRQ(ierr)
  if (OptionPrintToFile(this%option)) then
    close(this%option%fid_out)
    call OutputEKGFinalize()
  endif

end subroutine SimSubsurfFinalizeRun

! ************************************************************************** !

subroutine SimSubsurfStrip(this)
  !
  ! Deallocates members of subsurface simulation
  !
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  use Waypoint_module
  use Input_Aux_module
  use EOS_module

  implicit none

  class(simulation_subsurface_type) :: this

#ifdef DEBUG
  call PrintMsg(this%option,'SimSubsurfStrip()')
#endif

  call SimulationBaseStrip(this)

  call SimAuxDestroy(this%sim_aux)
  call OutputOptionDestroy(this%output_option)
  if (associated(this%process_model_coupler_list)) then
call this%process_model_coupler_list%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%process_model_coupler_list)
    nullify(this%process_model_coupler_list)
  endif
  call InputDbaseDestroy()
  call AllEOSDBaseDestroy()
  call RealizationStrip(this%realization)
  deallocate(this%realization)
  nullify(this%realization)
  call RegressionDestroy(this%regression)
  call WaypointListDestroy(this%waypoint_list_subsurface)
  call WaypointListDestroy(this%waypoint_list_outer)
  call OptionDestroy(this%option)

end subroutine SimSubsurfStrip

! ************************************************************************** !

subroutine SimSubsurfDestroy(simulation)
  !
  ! Deallocates a simulation
  !
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  !

  implicit none

  class(simulation_subsurface_type), pointer :: simulation

#ifdef DEBUG
  call PrintMsg(simulation%option,'SimulationDestroy()')
#endif

  if (.not.associated(simulation)) return

  call simulation%Strip()
  deallocate(simulation)
  nullify(simulation)

end subroutine SimSubsurfDestroy

end module Simulation_Subsurface_class
