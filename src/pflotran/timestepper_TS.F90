module Timestepper_TS_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Timestepper_Base_class
  use Solver_module
  use Waypoint_module

  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal

  implicit none

  private

  type, public, extends(timestepper_base_type) :: timestepper_TS_type
    PetscReal :: dt_max_allowable

    PetscInt :: num_newton_iterations ! number of Newton iterations in a time step
    PetscInt :: num_linear_iterations ! number of linear solver iterations in a time step
    PetscInt :: cumulative_newton_iterations       ! Total number of Newton iterations
    PetscInt :: cumulative_linear_iterations     ! Total number of linear iterations

    PetscInt :: iaccel        ! Accelerator index
    ! An array of multiplicative factors that specify how to increase time step.
    PetscReal, pointer :: tfac(:)
    PetscInt :: ntfac             ! size of tfac

  contains
    procedure, public :: CheckpointBinary => TimestepperTSCheckpointBinary
    procedure, public :: Init => TimestepperTSInit
    procedure, public :: ReadSelectCase => TimestepperTSReadSelectCase
    procedure, public :: RestartBinary => TimestepperTSRestartBinary
    procedure, public :: Reset => TimestepperTSReset
    procedure, public :: InputRecord => TimestepperSurfInputRecord
    procedure, public :: Strip => TimestepperTSStrip
    procedure, public :: StepDT => TimestepperTSStepDT
    procedure, public :: UpdateDT => TimestepperTSUpdateDT
    procedure, public :: FinalizeRun => TimestepperTSFinalizeRun
  end type timestepper_TS_type

  ! For checkpointing
  type, public, extends(stepper_base_header_type) :: timestepper_TS_header_type
    real*8 :: dt_max_allowable
  end type timestepper_TS_header_type
  PetscSizeT, parameter, private :: bagsize = 80 ! 64 (base) + 16 (SNES)

  interface PetscBagGetData
    subroutine PetscBagGetData(bag,header,ierr)
      import :: timestepper_TS_header_type
      implicit none
      PetscBag :: bag
      class(timestepper_TS_header_type), pointer :: header
      PetscErrorCode :: ierr
    end subroutine
  end interface PetscBagGetData

  public :: TimestepperTSCreate

contains

! ************************************************************************** !

function TimestepperTSCreate()
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/19/18
  !

  implicit none

  class(timestepper_TS_type), pointer :: TimestepperTSCreate

  class(timestepper_TS_type), pointer :: ts_timestepper

  allocate(ts_timestepper)

  call ts_timestepper%Init()

  TimestepperTSCreate => ts_timestepper

end function TimestepperTSCreate

! ************************************************************************** !

subroutine TimestepperTSInit(this)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/19/18
  !

  implicit none

  class (timestepper_TS_type) :: this

  call TimestepperBaseInit(this)

  this%num_newton_iterations = 0
  this%num_linear_iterations = 0

  this%cumulative_newton_iterations = 0
  this%cumulative_linear_iterations = 0

  this%dt_max_allowable = 0.d0
  this%iaccel = 5
  this%ntfac = 13
  allocate(this%tfac(13))
  this%tfac(1)  = 2.0d0; this%tfac(2)  = 2.0d0
  this%tfac(3)  = 2.0d0; this%tfac(4)  = 2.0d0
  this%tfac(5)  = 2.0d0; this%tfac(6)  = 1.8d0
  this%tfac(7)  = 1.6d0; this%tfac(8)  = 1.4d0
  this%tfac(9)  = 1.2d0; this%tfac(10) = 1.0d0
  this%tfac(11) = 1.0d0; this%tfac(12) = 1.0d0
  this%tfac(13) = 1.0d0

end subroutine TimestepperTSInit

! ************************************************************************** !

subroutine TimestepperTSReadSelectCase(this,input,keyword,found, &
                                       error_string,option)
  !
  ! Reads select case statement for TS
  !
  ! Author: Glenn Hammond
  ! Date: 03/16/20
  !

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module

  implicit none

  class(timestepper_TS_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  found = PETSC_TRUE
  call TimestepperBaseReadSelectCase(this,input,keyword,found, &
                                     error_string,option)
  if (found) return

  found = PETSC_TRUE
  select case(trim(keyword))

    case('TS_ACCELERATION')
      call InputReadInt(input,option,this%iaccel)
      call InputDefaultMsg(input,option,'iaccel')

    case('DT_FACTOR')
      string='time_step_factor'
      call UtilityReadArray(this%tfac,NEG_ONE_INTEGER,string,input, &
          option)
      this%ntfac = size(this%tfac)

    case default
      found = PETSC_FALSE
  end select

end subroutine TimestepperTSReadSelectCase

! ************************************************************************** !
! ************************************************************************** !

subroutine TimestepperTSStepDT(this,process_model,stop_flag)
  !
  ! This is a dummy routine added to be extended in timestepper_TS_type
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/19/18
  !

#include "petsc/finclude/petscts.h"
  use petscts
  use PM_Base_class
  use PM_Subsurface_Flow_class
  use Option_module
  use Output_module, only : Output
  use String_module

  implicit none

  class(timestepper_TS_type) :: this
  class(pm_base_type) :: process_model
  PetscInt :: stop_flag

  PetscReal :: time
  PetscReal :: dtime
  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  Vec :: residual_vec
  PetscReal :: fnorm, inorm, scaled_fnorm
  TSConvergedReason :: ts_reason
  PetscInt :: num_linear_iterations
  PetscInt :: num_newton_iterations
  PetscLogDouble :: log_start_time
  PetscLogDouble :: log_end_time
  PetscErrorCode :: ierr

  solver => process_model%solver
  option => process_model%option

  option%dt = this%dt

  call process_model%InitializeTimestep()

  call process_model%PreSolve()

  call TSSetTime(solver%ts,0.d0,ierr);CHKERRQ(ierr)
  call TSSetMaxTime(solver%ts,option%flow_dt,ierr);CHKERRQ(ierr)
  call TSSetTimeStep(solver%ts,option%flow_dt,ierr);CHKERRQ(ierr)

  call TSSetStepNumber(solver%ts,ZERO_INTEGER,ierr);CHKERRQ(ierr)
  call TSSetExactFinalTime(solver%ts,TS_EXACTFINALTIME_MATCHSTEP, &
                           ierr);CHKERRQ(ierr)

  call PetscTime(log_start_time,ierr);CHKERRQ(ierr)
  call TSSolve(solver%ts,process_model%solution_vec,ierr);CHKERRQ(ierr)
  call PetscTime(log_end_time,ierr);CHKERRQ(ierr)

  this%cumulative_solver_time = &
    this%cumulative_solver_time + &
    (log_end_time - log_start_time)

  call TSGetConvergedReason(solver%ts,ts_reason,ierr);CHKERRQ(ierr)
  call TSGetTime(solver%ts,time,ierr);CHKERRQ(ierr)
  call TSGetTimeStep(solver%ts,dtime,ierr);CHKERRQ(ierr)
  call TSGetSNESIterations(solver%ts,num_newton_iterations, &
                           ierr);CHKERRQ(ierr)
  call TSGetKSPIterations(solver%ts,num_linear_iterations,ierr);CHKERRQ(ierr)

  this%num_newton_iterations = num_newton_iterations
  this%cumulative_newton_iterations = this%cumulative_newton_iterations + this%num_newton_iterations

  this%num_linear_iterations = num_linear_iterations
  this%cumulative_linear_iterations = this%cumulative_linear_iterations + this%num_linear_iterations

  if (ts_reason<0) then
     write(*,*)'TS failed to converge. Stopping execution!'
     stop
  endif

  call TSGetRHSFunction(solver%ts,residual_vec,PETSC_NULL_FUNCTION, &
                        PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
  call VecNorm(residual_vec,NORM_2,fnorm,ierr);CHKERRQ(ierr)
  call VecNorm(residual_vec,NORM_INFINITY,inorm,ierr);CHKERRQ(ierr)

  this%steps = this%steps + 1

  call TimestepperBasePrintStepInfo(this,process_model%output_option, &
                                    ts_reason,option)
  write(option%io_buffer,'("  newton = ",i3," [",i8,"]", " linear = ",i5, &
                         &" [",i10,"]")') &
           num_newton_iterations,this%cumulative_newton_iterations, &
           num_linear_iterations,this%cumulative_linear_iterations
  call PrintMsg(option)

  select type(process_model)
    class is (pm_subsurface_flow_type)
      scaled_fnorm = fnorm/process_model%realization_base% &
                     discretization%grid%nmax
    class default
      scaled_fnorm = fnorm
  end select

  option%io_buffer = '  --> TS SNES Linear/Non-Linear Iterations = ' // &
    trim(StringWrite(this%num_linear_iterations)) // ' / ' // &
    trim(StringWrite(this%num_newton_iterations))
  call PrintMsg(option)
  option%io_buffer = '  --> TS SNES Residual: ' // &
    trim(StringWrite('(e14.6)',fnorm)) // ' ' // &
    trim(StringWrite('(e14.6)',scaled_fnorm)) // ' ' // &
    trim(StringWrite('(e14.6)',inorm))
  call PrintMsg(option)

  call process_model%FinalizeTimestep()

  call process_model%PostSolve()

end subroutine TimestepperTSStepDT

! ************************************************************************** !

subroutine TimestepperTSCheckpointBinary(this,viewer,option)
  !
  ! This checkpoints parameters/variables associated with surface-timestepper
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/19/18
  !

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_TS_type) :: this
  PetscViewer :: viewer
  type(option_type) :: option

  class(timestepper_TS_header_type), pointer :: header
  PetscBag :: bag
  PetscErrorCode :: ierr

  call PetscBagCreate(option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
  call TimestepperTSRegisterHeader(this,bag,header)
  call TimestepperTSSetHeader(this,bag,header)
  call PetscBagView(bag,viewer,ierr);CHKERRQ(ierr)
  call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)

end subroutine TimestepperTSCheckpointBinary

! ************************************************************************** !

subroutine TimestepperTSRestartBinary(this,viewer,option)
  !
  ! This checkpoints parameters/variables associated with surface-timestepper
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/19/18
  !

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_TS_type) :: this
  PetscViewer :: viewer
  type(option_type) :: option

  class(timestepper_TS_header_type), pointer :: header
  PetscBag :: bag
  PetscErrorCode :: ierr

  call PetscBagCreate(option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
  call TimestepperTSRegisterHeader(this,bag,header)
  call PetscBagLoad(viewer,bag,ierr);CHKERRQ(ierr)
  call TimestepperTSGetHeader(this,header)
  call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)

end subroutine TimestepperTSRestartBinary

! ************************************************************************** !

subroutine TimestepperTSRegisterHeader(this,bag,header)
  !
  ! This subroutine register header entries for surface-flow.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/19/18
  !

  use Option_module

  implicit none

  class(timestepper_TS_type) :: this
  class(timestepper_TS_header_type) :: header
  PetscBag :: bag

  PetscErrorCode :: ierr

  ! bagsize = 2 * 8 bytes = 16 bytes
  call PetscBagRegisterReal(bag,header%dt_max_allowable,0.d0, &
                            "dt_max_allowable","",ierr);CHKERRQ(ierr)

  call TimestepperBaseRegisterHeader(this,bag,header)

end subroutine TimestepperTSRegisterHeader

! ************************************************************************** !

subroutine TimestepperTSSetHeader(this,bag,header)
  !
  ! This subroutine sets values in checkpoint header.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/19/18
  !

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_TS_type) :: this
  class(timestepper_TS_header_type) :: header
  PetscBag :: bag

  header%dt_max_allowable = this%dt_max_allowable

  call TimestepperBaseSetHeader(this,bag,header)

end subroutine TimestepperTSSetHeader

! ************************************************************************** !

subroutine TimestepperTSGetHeader(this,header)
  !
  ! This subroutine gets values in checkpoint header.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/19/18
  !

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"

  class(timestepper_TS_type) :: this
  class(timestepper_TS_header_type) :: header

  PetscErrorCode :: ierr

  this%dt_max_allowable = header%dt_max_allowable

  call TimestepperBaseGetHeader(this,header)

  call TSSetTime(this%solver%ts,this%target_time,ierr);CHKERRQ(ierr)

end subroutine TimestepperTSGetHeader

! ************************************************************************** !

subroutine TimestepperTSReset(this)

  implicit none

  class(timestepper_TS_type) :: this

end subroutine TimestepperTSReset

! ************************************************************************** !

subroutine TimestepperTSPrintInfo(this,aux_string,option)
  !
  ! Prints settings for base timestepper.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/19/18
  !
  use Option_module

  implicit none

#include "petsc/finclude/petscts.h"

  class(timestepper_TS_type) :: this
  character(len=*) :: aux_string
  type(option_type) :: option

  PetscErrorCode :: ierr

  call TimestepperBasePrintInfo(this,'TS',option)
  if (OptionPrintToScreen(option)) then
    write(*,*) ' '
    write(*,*) 'TS Solver:'
    call TSView(this%solver%ts,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
  endif
  call SolverPrintNewtonInfo(this%solver,this%name,option)
  call SolverPrintLinearInfo(this%solver,this%name,option)

end subroutine TimestepperTSPrintInfo

! ************************************************************************** !

subroutine TimestepperSurfInputRecord(this)
  !
  ! Prints information about the time stepper to the input record.
  ! To get a## format, must match that in simulation types.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/19/18
  !

  implicit none

  class(timestepper_TS_type) :: this

  PetscInt :: id
  character(len=MAXWORDLENGTH) :: word

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pmc timestepper: '
  write(id,'(a)') this%name

  write(id,'(a29)',advance='no') 'max timestep size: '
  write(word,*) this%dt_max_allowable
  write(id,'(a)') trim(adjustl(word)) // ' sec'

end subroutine TimestepperSurfInputRecord

! ************************************************************************** !

subroutine TimestepperTSStrip(this)
  !
  ! Deallocates members of a surface time stepper
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/19/18
  !

  implicit none

  class(timestepper_TS_type) :: this

  call TimestepperBaseStrip(this)

end subroutine TimestepperTSStrip

! ************************************************************************** !

subroutine TimestepperTSDestroy(this)
  !
  ! Deallocates a surface time stepper
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/19/18
  !

  implicit none

  class(timestepper_TS_type) :: this

  call TimestepperTSStrip(this)

end subroutine TimestepperTSDestroy

! ************************************************************************** !
subroutine TimestepperTSUpdateDT(this,process_model)
  !
  ! Updates time step
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/19/18
  !

  use PM_Base_class

  implicit none

  class(timestepper_TS_type) :: this
  class(pm_base_type) :: process_model

  PetscBool :: update_time_step

  update_time_step = PETSC_TRUE

  if (this%time_step_cut_flag) then
    this%num_constant_time_steps = 1
  else if (this%num_constant_time_steps > 0) then
    ! otherwise, only increment if the constant time step counter was
    ! initialized to 1
    this%num_constant_time_steps = &
      this%num_constant_time_steps + 1
  endif

  ! num_constant_time_steps = 0: normal time stepping with growing steps
  ! num_constant_time_steps > 0: restriction of constant time steps until
  !                              constant_time_step_threshold is met
  if (this%num_constant_time_steps > &
      this%constant_time_step_threshold) then
    this%num_constant_time_steps = 0
  else if (this%num_constant_time_steps > 0) then
    ! do not increase time step size
    update_time_step = PETSC_FALSE
  endif

  call process_model%UpdateTimestep(update_time_step, &
                                    this%dt, &
                                    this%dt_min, &
                                    this%dt_max, &
                                    this%iaccel, &
                                    this%num_newton_iterations, &
                                    this%tfac, &
                                    this%time_step_max_growth_factor)

end subroutine TimestepperTSUpdateDT

! ************************************************************************** !

recursive subroutine TimestepperTSFinalizeRun(this,option)
  !
  ! Finalizes the time stepping
  !
  ! Author: Gautam Bisht
  ! Date: 07/03/18

  use Option_module
  use String_module

  implicit none

  class(timestepper_TS_type) :: this
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

#ifdef DEBUG
  call PrintMsg(option,'TimestepperSNESFinalizeRun()')
#endif

  string = ' ' // trim(this%name) // &
    ' PETSc TS steps = ' // trim(StringWrite(this%steps)) // &
    ' newton = ' // trim(StringWrite(this%cumulative_newton_iterations)) // &
    ' linear = ' // trim(StringWrite(this%cumulative_linear_iterations)) // &
    ' cuts = ' // trim(StringWrite(this%cumulative_time_step_cuts))
  call PrintMsg(option,string)
  string = ' ' // trim(this%name) // &
    ' PETSc TS time = ' // &
    trim(StringWrite('(f12.1)',this%cumulative_solver_time)) // ' seconds'
  call PrintMsg(option,string)

end subroutine TimestepperTSFinalizeRun

end module Timestepper_TS_class
