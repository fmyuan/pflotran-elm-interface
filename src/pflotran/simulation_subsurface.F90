module Simulation_Subsurface_class
  
#include "petsc/finclude/petscsys.h"
  use petscsys  
  use Simulation_Base_class
  use Regression_module
  use Option_module
  use PMC_Subsurface_class
  use PMC_Third_Party_class
  use PMC_Base_class
  use Realization_Subsurface_class
  use Waypoint_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(simulation_base_type) :: simulation_subsurface_type
    ! pointer to flow process model coupler
    class(pmc_subsurface_type), pointer :: flow_process_model_coupler
    ! pointer to transport process model coupler
    class(pmc_subsurface_type), pointer :: tran_process_model_coupler
    ! pointer to realization object shared by flow and reactive transport
    class(realization_subsurface_type), pointer :: realization 
    ! regression object
    type(regression_type), pointer :: regression
    type(waypoint_list_type), pointer :: waypoint_list_subsurface
  contains
    procedure, public :: Init => SimSubsurfInit
    procedure, public :: JumpStart => SimSubsurfJumpStart
    procedure, public :: InitializeRun => SimSubsurfInitializeRun
    procedure, public :: InputRecord => SimSubsurfInputRecord
!    procedure, public :: ExecuteRun
!    procedure, public :: RunToTime
    procedure, public :: FinalizeRun => SimSubsurfFinalizeRun
    procedure, public :: Strip => SimSubsurfStrip
  end type simulation_subsurface_type
  
  public :: SimSubsurfCreate, &
            SimSubsurfInit, &
            SimSubsurfFinalizeRun, &
            SimSubsurfStrip, &
            SimSubsurfDestroy
  
contains

! ************************************************************************** !

function SimSubsurfCreate(option)
  ! 
  ! Allocates and initializes a new simulation object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  use Option_module
  
  implicit none
  
  type(option_type), pointer :: option

  class(simulation_subsurface_type), pointer :: SimSubsurfCreate
  
#ifdef DEBUG
  print *, 'SimSubsurfCreate'
#endif
  
  allocate(SimSubsurfCreate)
  call SimSubsurfCreate%Init(option)
  
end function SimSubsurfCreate

! ************************************************************************** !

subroutine SimSubsurfInit(this,option)
  ! 
  ! Initializes simulation values
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/22/13
  ! 
  use Waypoint_module
  use Option_module
  
  implicit none
  
  class(simulation_subsurface_type) :: this
  type(option_type), pointer :: option
  
  call SimulationBaseInit(this,option)
  nullify(this%flow_process_model_coupler)
  nullify(this%tran_process_model_coupler)
  nullify(this%realization)
  nullify(this%regression)
  this%waypoint_list_subsurface => WaypointListCreate()
  
end subroutine SimSubsurfInit

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
  use Output_Aux_module
  use hdf5

  implicit none
  

  class(simulation_subsurface_type) :: this

  PetscBool :: flag
  PetscErrorCode :: ierr
  
#ifdef DEBUG
  call PrintMsg(this%option,'SimSubsurfInitializeRun()')
#endif
  
  ! print error message if binary checkpoint/restart is used in
  ! combination with unstructured gridding
  flag = PETSC_FALSE
  if (associated(this%checkpoint_option)) then
    if (this%checkpoint_option%format == CHECKPOINT_BINARY) then
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

  call SimulationBaseInitializeRun(this)
  
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

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id = INPUT_RECORD_UNIT
  
  if (OptionPrintToScreen(this%option)) then
    write (*,*) 'Printing input record file.'
  endif
  
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
  use Timestepper_Base_class

  implicit none
  
  class(simulation_subsurface_type) :: this

  class(timestepper_base_type), pointer :: master_timestepper
  class(timestepper_base_type), pointer :: flow_timestepper
  class(timestepper_base_type), pointer :: tran_timestepper
  type(option_type), pointer :: option
  type(output_option_type), pointer :: output_option
  PetscBool :: snapshot_plot_flag, observation_plot_flag, massbal_plot_flag
  
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
                       trim(adjustl(option%io_buffer)) // &
                       '), specified by TIMESTEPPER->MAX_STEPS, ' // &
                       'has been met.  Stopping....'  
    call PrintMsg(option)
    call PrintMsg(option,'')
    option%status = DONE
    this%stop_flag = TS_STOP_MAX_TIME_STEP
    return
  endif

  ! print initial condition output if not a restarted sim
  call OutputInit(option,master_timestepper%steps)
  if (output_option%plot_number == 0 .and. &
      master_timestepper%max_time_step >= 0) then
    if (output_option%print_initial_snap) snapshot_plot_flag = PETSC_TRUE
    if (output_option%print_initial_obs) observation_plot_flag = PETSC_TRUE
    if (output_option%print_initial_massbal) massbal_plot_flag = PETSC_FALSE
    call Output(this%realization,snapshot_plot_flag,observation_plot_flag, &
                massbal_plot_flag)
  endif
  
  !if TIMESTEPPER->MAX_STEPS < 1, print out initial condition only
  if (master_timestepper%max_time_step < 1) then
    call PrintMsg(option,'')
    write(option%io_buffer,*) master_timestepper%max_time_step
    option%io_buffer = 'The maximum # of time steps (' // &
                       trim(adjustl(option%io_buffer)) // &
                       '), specified by TIMESTEPPER->MAX_STEPS, ' // &
                       'has been met.  Stopping....'  
    call PrintMsg(option)
    call PrintMsg(option,'')
    option%status = DONE
    this%stop_flag = TS_STOP_MAX_TIME_STEP
    return
  endif

  ! increment plot number so that 000 is always the initial condition, 
  ! and nothing else
  if (output_option%plot_number == 0) output_option%plot_number = 1

  if (associated(flow_timestepper)) then
    if (.not.associated(flow_timestepper%cur_waypoint)) then
      option%io_buffer = &
        'Null flow waypoint list; final time likely equal to start time.&
        &time or simulation time needs to be extended on a restart.'
      call PrintMsg(option)
      option%status = FAIL
      return
    else
      flow_timestepper%dt_max = flow_timestepper%cur_waypoint%dt_max
    endif
  endif  
  if (associated(tran_timestepper)) then
    if (.not.associated(tran_timestepper%cur_waypoint)) then
      option%io_buffer = &
        'Null transport waypoint list; final time likely equal to start &
        &time or simulation time needs to be extended on a restart.'
      call PrintMsg(option)
      option%status = FAIL
      return
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

subroutine SimSubsurfFinalizeRun(this)
  ! 
  ! Finalizes simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  use Timestepper_Base_class
  use Reaction_Sandbox_module, only : RSandboxDestroy
  use SrcSink_Sandbox_module, only : SSSandboxDestroyList
  use WIPP_module, only : WIPPDestroy
  use Klinkenberg_module, only : KlinkenbergDestroy
  use CLM_Rxn_module, only : RCLMRxnDestroy

  implicit none
  
  class(simulation_subsurface_type) :: this
  
  PetscErrorCode :: ierr
  
  class(timestepper_base_type), pointer :: flow_timestepper
  class(timestepper_base_type), pointer :: tran_timestepper

#ifdef DEBUG
  call PrintMsg(this%option,'SimSubsurfFinalizeRun()')
#endif
  
  call SimulationBaseFinalizeRun(this)
  
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
  
end subroutine SimSubsurfFinalizeRun

! ************************************************************************** !

subroutine SimSubsurfStrip(this)
  ! 
  ! Deallocates members of subsurface simulation
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  ! 

  implicit none
  
  class(simulation_subsurface_type) :: this
  
#ifdef DEBUG
  call PrintMsg(this%option,'SimSubsurfStrip()')
#endif
  
  call SimulationBaseStrip(this)
  call RealizationStrip(this%realization)
  deallocate(this%realization)
  nullify(this%realization)
  call RegressionDestroy(this%regression)
  call WaypointListDestroy(this%waypoint_list_subsurface)
  
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
