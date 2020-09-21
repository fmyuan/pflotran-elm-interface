module Timestepper_SNES_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Solver_module
  use Convergence_module
  use Timestepper_Base_class

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(timestepper_base_type) :: timestepper_SNES_type

    PetscInt :: num_newton_iterations ! number of Newton iterations in a time step
    PetscInt :: num_linear_iterations ! number of linear solver iterations in a time step
    PetscInt :: cumulative_newton_iterations       ! Total number of Newton iterations
    PetscInt :: cumulative_linear_iterations     ! Total number of linear iterations
    PetscInt :: cumulative_wasted_linear_iterations
    PetscInt :: cumulative_wasted_newton_iterations

    PetscInt :: iaccel        ! Accelerator index
    ! An array of multiplicative factors that specify how to increase time step.
    PetscReal, pointer :: tfac(:)
    PetscInt :: ntfac             ! size of tfac

    ! rescue mode related parameters - heeho
    PetscBool :: rescue_mode  ! increase dt when a simulation is in 
                              ! a hole or oscillation
    PetscInt  :: rescue_frequency       ! how often do we rescue?
    PetscReal :: rescue_factor          ! what do we increase dt by?
    PetscReal :: rescue_step_threshold  ! how small should dt be compared to t
    PetscInt  :: rescue_step_counter    ! counter to trigger rescue

  contains

    procedure, public :: ReadSelectCase => TimestepperSNESReadSelectCase
    procedure, public :: Init => TimestepperSNESInit
!    procedure, public :: SetTargetTime => TimestepperBaseSetTargetTime
    procedure, public :: StepDT => TimestepperSNESStepDT
    procedure, public :: UpdateDT => TimestepperSNESUpdateDT
    procedure, public :: CheckpointBinary => TimestepperSNESCheckpointBinary
    procedure, public :: RestartBinary => TimestepperSNESRestartBinary
    procedure, public :: CheckpointHDF5 => TimestepperSNESCheckpointHDF5
    procedure, public :: RestartHDF5 => TimestepperSNESRestartHDF5
    procedure, public :: Reset => TimestepperSNESReset
    procedure, public :: PrintInfo => TimestepperSNESPrintInfo
    procedure, public :: InputRecord => TimestepperSNESInputRecord
    procedure, public :: FinalizeRun => TimestepperSNESFinalizeRun
    procedure, public :: Strip => TimestepperSNESStrip
    procedure, public :: Destroy => TimestepperSNESDestroy

  end type timestepper_SNES_type

  ! For checkpointing
  type, public, extends(stepper_base_header_type) :: stepper_SNES_header_type
    PetscInt :: cumulative_newton_iterations
    PetscInt :: cumulative_linear_iterations
    PetscInt :: num_newton_iterations
  end type stepper_SNES_header_type

  interface PetscBagGetData
    subroutine PetscBagGetData(bag,header,ierr)
      import :: stepper_SNES_header_type
      implicit none
#include "petsc/finclude/petscbag.h"
      PetscBag :: bag
      class(stepper_SNES_header_type), pointer :: header
      PetscErrorCode :: ierr
    end subroutine
  end interface PetscBagGetData

  public :: TimestepperSNESCreate, &
            TimestepperSNESPrintInfo, &
            TimestepperSNESInit, &
            TimestepperSNESReadSelectCase

contains

! ************************************************************************** !

function TimestepperSNESCreate()
  !
  ! Allocates and initializes a new Timestepper object
  !
  ! Author: Glenn Hammond
  ! Date: 07/22/13
  !

  implicit none

  class(timestepper_SNES_type), pointer :: TimestepperSNESCreate

  class(timestepper_SNES_type), pointer :: stepper

  allocate(stepper)
  call stepper%Init()

  TimestepperSNESCreate => stepper

end function TimestepperSNESCreate

! ************************************************************************** !

subroutine TimestepperSNESInit(this)
  !
  ! Allocates and initializes a new Timestepper object
  !
  ! Author: Glenn Hammond
  ! Date: 07/22/13
  !

  implicit none

  class(timestepper_SNES_type) :: this

  call TimestepperBaseInit(this)

  this%num_newton_iterations = 0
  this%num_linear_iterations = 0

  this%cumulative_newton_iterations = 0
  this%cumulative_linear_iterations = 0
  this%cumulative_wasted_linear_iterations = 0
  this%cumulative_wasted_newton_iterations = 0

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

  ! rescue mode defaults - heeho
  this%rescue_mode = PETSC_FALSE
  this%rescue_frequency = 100
  this%rescue_factor = 1.0d3
  this%rescue_step_threshold = 1.0d-5
  this%rescue_step_counter = 0

end subroutine TimestepperSNESInit

! ************************************************************************** !

subroutine TimestepperSNESReadSelectCase(this,input,keyword,found, &
                                       error_string,option)
  !
  ! Reads select case statement for SNES
  !
  ! Author: Glenn Hammond
  ! Date: 03/16/20
  !

  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module

  implicit none

  class(timestepper_SNES_type) :: this
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

    case ('RESCUE_MODE')  ! heeho
      this%rescue_mode = PETSC_TRUE
      input%ierr = 0
      call InputPushBlock(input,option)
      do
        call InputReadPflotranString(input,option)
        if (InputCheckExit(input,option)) exit

        call InputReadCard(input,option,keyword)
        call InputErrorMsg(input,option,'keyword','RESCUE MODE')
        call StringToUpper(keyword)

        select case(trim(keyword))
          case('RESCUE_FACTOR','FACTOR','RFAC')
            call InputReadDouble(input,option,this%rescue_factor)
            call InputErrorMsg(input,option, &
                               'rescue factor ', &
                               'RESCUE mode options')
          case('RESCUE_FREQUENCY','FREQUENCY','RFREQ')
            call InputReadInt(input,option,this%rescue_frequency)
            call InputErrorMsg(input,option, &
                               'rescue frequency ', &
                               'RESCUE mode options')
          case('RESCUE_STEP_THRESHOLD','THRESHOLD','STEP_THRESHOLD','RTHRESH')
            call InputReadDouble(input,option,this%rescue_step_threshold)
            call InputErrorMsg(input,option, &
                               'rescue threshold ', &
                               'RESCUE mode options')
          case default
            option%io_buffer  = 'Timestepper RESCUE MODE option: ' // trim(keyword) // &
                              ' unknown.'
            call PrintErrMsg(option)
        end select
      enddo
      call InputPopBlock(input,option)

    case default
      found = PETSC_FALSE
  end select

end subroutine TimestepperSNESReadSelectCase

! ************************************************************************** !

subroutine TimestepperSNESUpdateDT(this,process_model)
  !
  ! Updates time step
  !
  ! Author: Glenn Hammond
  ! Date: 07/22/13
  !

  use Option_module
  use PM_Base_class

  implicit none

  class(timestepper_SNES_type) :: this
  class(pm_base_type) :: process_model
  type(option_type), pointer :: option
  
  PetscBool :: update_time_step
  
  option => process_model%option

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

  if (update_time_step .and. this%iaccel /= 0) then

    call process_model%UpdateTimestep(this%dt, &
                                      this%dt_min, &
                                      this%dt_max, &
                                      this%iaccel, &
                                      this%num_newton_iterations, &
                                      this%tfac, &
                                      this%time_step_max_growth_factor)

  endif

  ! rescue mode - heeho
  if (this%rescue_mode) then
    if (this%target_time*this%rescue_step_threshold > this%dt) then
      this%rescue_step_counter = this%rescue_step_counter + 1
    else
      if (this%rescue_step_counter > 0) then
        ! subtract from the counter if it recovers itself
        this%rescue_step_counter = this%rescue_step_counter - 1
      endif
    endif
    if (this%rescue_step_counter > this%rescue_frequency) then
      this%rescue_step_counter = 0
#if 0
      if (2**this%max_time_step_cuts < this%rescue_factor) then
        ! can't jump timestep more than the max time step cut is allowed
        if (this%max_time_step_cuts < 2) then
          option%io_buffer  = 'max time step cut too small for rescue mode. exiting.'
          call PrintErrMsg(option)
        endif
        this%rescue_factor = 2**(this%max_time_step_cuts-2)
        option%io_buffer = 'rescue factor too big. automatically adjusted.'
        call PrintMsg(option)
      endif
#endif
      this%dt = this%dt * this%rescue_factor
      option%io_buffer = 'rescue mode activated. jumping time step size.'
      call PrintMsg(option)
    endif
  endif

end subroutine TimestepperSNESUpdateDT

! ************************************************************************** !

subroutine TimestepperSNESStepDT(this,process_model,stop_flag)
  !
  ! Steps forward one step in time
  !
  ! Author: Glenn Hammond
  ! Date: 07/22/13
  !

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use Option_module
  use Output_module, only : Output, OutputFindNaNOrInfInVec
  use Output_EKG_module, only : IUNIT_EKG

  implicit none

  class(timestepper_SNES_type) :: this
  class(pm_base_type) :: process_model
  PetscInt :: stop_flag

  SNESConvergedReason :: snes_reason
  PetscInt :: icut

  type(solver_type), pointer :: solver
  type(option_type), pointer :: option

  PetscLogDouble :: log_start_time
  PetscLogDouble :: log_end_time
  PetscInt :: num_newton_iterations
  PetscInt :: num_linear_iterations
  PetscInt :: sum_newton_iterations
  PetscInt :: sum_linear_iterations
  PetscInt :: sum_wasted_linear_iterations
  PetscInt :: sum_wasted_newton_iterations
  character(len=MAXWORDLENGTH) :: tunit

  PetscReal :: tconv
  PetscReal :: fnorm, inorm, scaled_fnorm
  PetscBool :: snapshot_plot_flag, observation_plot_flag, massbal_plot_flag
  Vec :: residual_vec
  PetscErrorCode :: ierr

  solver => process_model%solver
  option => process_model%option

!geh: for debugging
#ifdef DEBUG
  write(process_model%option%io_buffer,'(es12.5)') this%dt
  process_model%option%io_buffer = 'StepperStepDT(' // &
    trim(adjustl(process_model%option%io_buffer)) // ')'
  call PrintMsg(process_model%option)
#endif

  tconv = process_model%output_option%tconv
  tunit = process_model%output_option%tunit
  sum_linear_iterations = 0
  sum_wasted_linear_iterations = 0
  sum_wasted_newton_iterations = 0
  sum_newton_iterations = 0
  icut = 0

  option%dt = this%dt
  option%time = this%target_time-this%dt

  call process_model%InitializeTimestep()

  do

    call process_model%PreSolve()

    call PetscTime(log_start_time,ierr);CHKERRQ(ierr)

    call SNESSolve(solver%snes,PETSC_NULL_VEC,process_model%solution_vec, &
                   ierr);CHKERRQ(ierr)

    call PetscTime(log_end_time,ierr);CHKERRQ(ierr)

    this%cumulative_solver_time = &
      this%cumulative_solver_time + &
      (log_end_time - log_start_time)

    call SNESGetIterationNumber(solver%snes,num_newton_iterations, &
                                ierr);CHKERRQ(ierr)
    call SNESGetLinearSolveIterations(solver%snes,num_linear_iterations, &
                                      ierr);CHKERRQ(ierr)
    call SNESGetConvergedReason(solver%snes,snes_reason,ierr);CHKERRQ(ierr)

    sum_newton_iterations = sum_newton_iterations + num_newton_iterations
    sum_linear_iterations = sum_linear_iterations + num_linear_iterations

    if (snes_reason <= 0 .or. .not. process_model%AcceptSolution()) then
      sum_wasted_linear_iterations = sum_wasted_linear_iterations + &
           num_linear_iterations
      sum_wasted_newton_iterations = sum_wasted_newton_iterations + &
           num_newton_iterations
      ! The Newton solver diverged, so try reducing the time step.
      icut = icut + 1
      this%time_step_cut_flag = PETSC_TRUE
      ! if a cut occurs on the last time step, the stop_flag will have been
      ! set to TS_STOP_END_SIMULATION.  Set back to TS_CONTINUE to prevent
      ! premature ending of simulation.
      if (stop_flag /= TS_STOP_MAX_TIME_STEP) stop_flag = TS_CONTINUE

      if (icut > this%max_time_step_cuts .or. this%dt < this%dt_min) then

        if (icut > this%max_time_step_cuts) then
          option%io_buffer = ' Stopping: Time step cut criteria exceeded.'
          call PrintMsg(option)
          write(option%io_buffer, &
                '("    icut =",i3,", max_time_step_cuts=",i3)') &
                icut,this%max_time_step_cuts
          call PrintMsg(option)
        endif
        if (this%dt < this%dt_min) then
          option%io_buffer = ' Stopping: Time step size is less than the &
                             &minimum allowable time step.'
          call PrintMsg(option)
          write(option%io_buffer, &
                '("    dt   =",es15.7,", dt_min=",es15.7," [",a,"]")') &
               this%dt/tconv,this%dt_min/tconv,trim(tunit)
          call PrintMsg(option)
       endif

        process_model%output_option%plot_name = trim(process_model%name)// &
          '_cut_to_failure'
        snapshot_plot_flag = PETSC_TRUE
        observation_plot_flag = PETSC_FALSE
        massbal_plot_flag = PETSC_FALSE
        call Output(process_model%realization_base,snapshot_plot_flag, &
                    observation_plot_flag,massbal_plot_flag)
        stop_flag = TS_STOP_FAILURE
        return
      endif

      this%target_time = this%target_time - this%dt

      this%dt = this%time_step_reduction_factor * this%dt

      write(option%io_buffer,'(''-> Cut time step: snes='',i3, &
           &   '' icut= '',i2,''['',i3,'']'','' t= '',1pe12.5, '' dt= '', &
           &   1pe12.5)')  snes_reason,icut, &
           this%cumulative_time_step_cuts+icut, &
           option%time/tconv, &
           this%dt/tconv
      call PrintMsg(option)
      if (snes_reason < SNES_CONVERGED_ITERATING) then
        call SolverNewtonPrintFailedReason(solver,option)
        if (solver%verbose_logging) then
          select case(snes_reason)
            case(SNES_DIVERGED_FNORM_NAN)
              ! attempt to find cells with NaNs.
              call SNESGetFunction(solver%snes,residual_vec, &
                                   PETSC_NULL_FUNCTION,PETSC_NULL_INTEGER, &
                                   ierr);CHKERRQ(ierr)
              call OutputFindNaNOrInfInVec(residual_vec, &
                                           process_model%realization_base% &
                                             discretization%grid,option)
          end select
        endif
      endif

      this%target_time = this%target_time + this%dt
      option%dt = this%dt
      call process_model%TimeCut()

    else
      ! The Newton solver converged, so we can exit.
      exit
    endif
  enddo

  this%steps = this%steps + 1
  this%cumulative_newton_iterations = &
    this%cumulative_newton_iterations + sum_newton_iterations
  this%cumulative_linear_iterations = &
    this%cumulative_linear_iterations + sum_linear_iterations
  this%cumulative_wasted_linear_iterations = &
       this%cumulative_wasted_linear_iterations + sum_wasted_linear_iterations
  this%cumulative_wasted_newton_iterations = &
       this%cumulative_wasted_newton_iterations + sum_wasted_newton_iterations
  this%cumulative_time_step_cuts = &
    this%cumulative_time_step_cuts + icut

  this%num_newton_iterations = num_newton_iterations
  this%num_linear_iterations = num_linear_iterations

  ! print screen output
  call SNESGetFunction(solver%snes,residual_vec,PETSC_NULL_FUNCTION, &
                       PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
  call VecNorm(residual_vec,NORM_2,fnorm,ierr);CHKERRQ(ierr)
  call VecNorm(residual_vec,NORM_INFINITY,inorm,ierr);CHKERRQ(ierr)
  if (option%print_screen_flag) then
      write(*, '(/," Step ",i6," Time= ",1pe12.5," Dt= ",1pe12.5, &
           & " [",a,"]", " snes_conv_reason: ",i4,/,"  newton = ",i3, &
           & " [",i8,"]", " linear = ",i5," [",i10,"]"," cuts = ",i2, &
           & " [",i4,"]")') &
           this%steps, &
           this%target_time/tconv, &
           this%dt/tconv, &
           trim(tunit),snes_reason,sum_newton_iterations, &
           this%cumulative_newton_iterations,sum_linear_iterations, &
           this%cumulative_linear_iterations,icut, &
           this%cumulative_time_step_cuts


    if (associated(process_model%realization_base%discretization%grid)) then
       scaled_fnorm = fnorm/process_model%realization_base% &
                        discretization%grid%nmax
    else
       scaled_fnorm = fnorm
    endif

    print *,' --> SNES Linear/Non-Linear Iterations = ', &
             num_linear_iterations,' / ',num_newton_iterations
    write(*,'("  --> SNES Residual: ",1p3e14.6)') fnorm, scaled_fnorm, inorm
  endif

  if (option%print_file_flag) then
    write(option%fid_out, '(" Step ",i6," Time= ",1pe12.5," Dt= ",1pe12.5, &
      & " [",a,"]"," snes_conv_reason: ",i4,/,"  newton = ",i3, &
      & " [",i8,"]", " linear = ",i5," [",i10,"]"," cuts = ",i2," [",i4,"]")') &
      this%steps, &
      this%target_time/tconv, &
      this%dt/tconv, &
      trim(tunit),snes_reason,sum_newton_iterations, &
      this%cumulative_newton_iterations,sum_linear_iterations, &
      this%cumulative_linear_iterations,icut, &
      this%cumulative_time_step_cuts
  endif

  option%time = this%target_time
  call process_model%FinalizeTimestep()

  if (this%print_ekg .and. OptionPrintToFile(option)) then
100 format(a32," TIMESTEP ",i10,2es16.8,a,i3,i5,i3,i5,i5,i10)
    write(IUNIT_EKG,100) trim(this%name), this%steps, this%target_time/tconv, &
      this%dt/tconv, trim(tunit), &
      icut, this%cumulative_time_step_cuts, &
      sum_newton_iterations, this%cumulative_newton_iterations, &
      sum_linear_iterations, this%cumulative_linear_iterations
  endif

  if (option%print_screen_flag) print *, ""

end subroutine TimestepperSNESStepDT

! ************************************************************************** !

subroutine TimestepperSNESCheckpointBinary(this,viewer,option)
  !
  ! Checkpoints parameters/variables associated with
  ! a time stepper.
  !
  ! Author: Glenn Hammond
  ! Date: 07/25/13
  !
  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_SNES_type) :: this
  PetscViewer :: viewer
  type(option_type) :: option

  class(stepper_SNES_header_type), pointer :: header
  type(stepper_SNES_header_type) :: dummy_header
  character(len=1),pointer :: dummy_char(:)
  PetscBag :: bag
  PetscSizeT :: bagsize
  PetscErrorCode :: ierr

  bagsize = size(transfer(dummy_header,dummy_char))

  call PetscBagCreate(option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
  call TimestepperSNESRegisterHeader(this,bag,header)
  call TimestepperSNESSetHeader(this,bag,header)
  call PetscBagView(bag,viewer,ierr);CHKERRQ(ierr)
  call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)

end subroutine TimestepperSNESCheckpointBinary

! ************************************************************************** !

subroutine TimestepperSNESRegisterHeader(this,bag,header)
  !
  ! Register header entries.
  !
  ! Author: Glenn Hammond
  ! Date: 07/30/13
  !

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_SNES_type) :: this
  class(stepper_SNES_header_type) :: header
  PetscBag :: bag

  PetscErrorCode :: ierr

  call PetscBagRegisterInt(bag,header%cumulative_newton_iterations,0, &
                           "cumulative_newton_iterations","", &
                           ierr);CHKERRQ(ierr)
  call PetscBagRegisterInt(bag,header%cumulative_linear_iterations,0, &
                           "cumulative_linear_iterations","", &
                           ierr);CHKERRQ(ierr)
! need to add cumulative wasted linear iterations
  call PetscBagRegisterInt(bag,header%num_newton_iterations,0, &
                           "num_newton_iterations","",ierr);CHKERRQ(ierr)

  call TimestepperBaseRegisterHeader(this,bag,header)

end subroutine TimestepperSNESRegisterHeader

! ************************************************************************** !

subroutine TimestepperSNESSetHeader(this,bag,header)
  !
  ! Sets values in checkpoint header.
  !
  ! Author: Glenn Hammond
  ! Date: 07/25/13
  !

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_SNES_type) :: this
  class(stepper_SNES_header_type) :: header
  PetscBag :: bag

  PetscErrorCode :: ierr

  header%cumulative_newton_iterations = this%cumulative_newton_iterations
  header%cumulative_linear_iterations = this%cumulative_linear_iterations
  header%num_newton_iterations = this%num_newton_iterations

  call TimestepperBaseSetHeader(this,bag,header)

end subroutine TimestepperSNESSetHeader

! ************************************************************************** !

subroutine TimestepperSNESRestartBinary(this,viewer,option)
  !
  ! Checkpoints parameters/variables associated with
  ! a time stepper.
  !
  ! Author: Glenn Hammond
  ! Date: 07/25/13
  !

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_SNES_type) :: this
  PetscViewer :: viewer
  type(option_type) :: option

  class(stepper_SNES_header_type), pointer :: header
  type(stepper_SNES_header_type) :: dummy_header
  character(len=1),pointer :: dummy_char(:)
  PetscBag :: bag
  PetscSizeT :: bagsize
  PetscErrorCode :: ierr

  bagsize = size(transfer(dummy_header,dummy_char))

  call PetscBagCreate(option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
  call TimestepperSNESRegisterHeader(this,bag,header)
  call PetscBagLoad(viewer,bag,ierr);CHKERRQ(ierr)
  call TimestepperSNESGetHeader(this,header)
  call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)

end subroutine TimestepperSNESRestartBinary

! ************************************************************************** !

subroutine TimestepperSNESCheckpointHDF5(this, h5_chk_grp_id, option)
  !
  ! Checkpoints parameters/variables associated with
  ! a time stepper.
  !
  ! Author: Gautam Bisht
  ! Date: 07/30/15
  !
  use Option_module
  use hdf5
  use Checkpoint_module, only : CheckPointWriteIntDatasetHDF5
  use Checkpoint_module, only : CheckPointWriteRealDatasetHDF5

  implicit none

  class(timestepper_SNES_type) :: this
  integer(HID_T) :: h5_chk_grp_id
  type(option_type) :: option

  integer(HSIZE_T), pointer :: dims(:)
  integer(HSIZE_T), pointer :: start(:)
  integer(HSIZE_T), pointer :: stride(:)
  integer(HSIZE_T), pointer :: length(:)
  integer(HID_T) :: timestepper_grp_id

  PetscMPIInt :: dataset_rank
  character(len=MAXSTRINGLENGTH) :: dataset_name
  character(len=MAXSTRINGLENGTH) :: string
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_array(:)
  PetscReal, pointer :: real_array(:)
  PetscMPIInt :: hdf5_err

  string = "Timestepper"
  call h5gcreate_f(h5_chk_grp_id, string, timestepper_grp_id, &
                   hdf5_err, OBJECT_NAMELEN_DEFAULT_F)

  allocate(start(1))
  allocate(dims(1))
  allocate(length(1))
  allocate(stride(1))
  allocate(int_array(1))
  allocate(real_array(1))

  dataset_rank = 1
  dims(1) = ONE_INTEGER
  start(1) = 0
  length(1) = ONE_INTEGER
  stride(1) = ONE_INTEGER

  dataset_name = "Cumulative_newton_iterations" // CHAR(0)
  int_array(1) = this%cumulative_newton_iterations
  call CheckPointWriteIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)

  dataset_name = "Cumulative_linear_iterations" // CHAR(0)
  int_array(1) = this%cumulative_linear_iterations
  call CheckPointWriteIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)

  dataset_name = "Num_newton_iterations" // CHAR(0)
  int_array(1) = this%num_newton_iterations
  call CheckPointWriteIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)

  dataset_name = "Time" // CHAR(0)
  real_array(1) = this%target_time
  call CheckPointWriteRealDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, real_array, option)

  dataset_name = "Dt" // CHAR(0)
  real_array(1) = this%dt
  call CheckPointWriteRealDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, real_array, option)

  dataset_name = "Prev_dt" // CHAR(0)
  real_array(1) = this%prev_dt
  call CheckPointWriteRealDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, real_array, option)

  dataset_name = "Num_steps" // CHAR(0)
  int_array(1) = this%steps
  call CheckPointWriteIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)

  dataset_name = "Cumulative_time_step_cuts" // CHAR(0)
  int_array(1) = this%cumulative_time_step_cuts
  call CheckPointWriteIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)

  dataset_name = "Num_constant_time_steps" // CHAR(0)
  int_array(1) = this%num_constant_time_steps
  call CheckPointWriteIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)

  dataset_name = "Num_contig_revert_due_to_sync" // CHAR(0)
  int_array(1) = this%num_contig_revert_due_to_sync
  call CheckPointWriteIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)

  dataset_name = "Revert_dt" // CHAR(0)
  int_array(1) = ZERO_INTEGER
  if (this%revert_dt) int_array(1) = ONE_INTEGER
  call CheckPointWriteIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)

  call h5gclose_f(timestepper_grp_id, hdf5_err)

  deallocate(start)
  deallocate(dims)
  deallocate(length)
  deallocate(stride)
  deallocate(int_array)
  deallocate(real_array)

end subroutine TimestepperSNESCheckpointHDF5

! ************************************************************************** !

subroutine TimestepperSNESRestartHDF5(this, h5_chk_grp_id, option)
  !
  ! Restarts parameters/variables associated with
  ! a time stepper.
  !
  ! Author: Gautam Bisht
  ! Date: 08/16/15
  !
  use Option_module
  use hdf5
  use HDF5_Aux_module
  use Checkpoint_module, only : CheckPointReadIntDatasetHDF5
  use Checkpoint_module, only : CheckPointReadRealDatasetHDF5

  implicit none

  class(timestepper_SNES_type) :: this
  integer(HID_T) :: h5_chk_grp_id
  type(option_type) :: option

  integer(HSIZE_T), pointer :: dims(:)
  integer(HSIZE_T), pointer :: start(:)
  integer(HSIZE_T), pointer :: stride(:)
  integer(HSIZE_T), pointer :: length(:)
  integer(HID_T) :: timestepper_grp_id

  PetscMPIInt :: dataset_rank
  character(len=MAXSTRINGLENGTH) :: dataset_name
  character(len=MAXSTRINGLENGTH) :: string
  ! must be 'integer' so that ibuffer does not switch to 64-bit integers
  ! when PETSc is configured with --with-64-bit-indices=yes.
  integer, pointer :: int_array(:)
  PetscReal, pointer :: real_array(:)
  PetscMPIInt :: hdf5_err

  string = "Timestepper"
  call HDF5GroupOpen(h5_chk_grp_id,string,timestepper_grp_id,option)

  allocate(start(1))
  allocate(dims(1))
  allocate(length(1))
  allocate(stride(1))
  allocate(int_array(1))
  allocate(real_array(1))

  dataset_rank = 1
  dims(1) = ONE_INTEGER
  start(1) = 0
  length(1) = ONE_INTEGER
  stride(1) = ONE_INTEGER

  dataset_name = "Cumulative_newton_iterations" // CHAR(0)
  call CheckPointReadIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                    dataset_rank, dims, start, length, &
                                    stride, int_array, option)
  this%cumulative_newton_iterations = int_array(1)

  dataset_name = "Cumulative_linear_iterations" // CHAR(0)
  call CheckPointReadIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                    dataset_rank, dims, start, length, &
                                    stride, int_array, option)
  this%cumulative_linear_iterations = int_array(1)

  dataset_name = "Num_newton_iterations" // CHAR(0)
  call CheckPointReadIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                    dataset_rank, dims, start, length, &
                                    stride, int_array, option)
  this%num_newton_iterations = int_array(1)

  dataset_name = "Time" // CHAR(0)
  call CheckPointReadRealDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, real_array, option)
  this%target_time = real_array(1)

  dataset_name = "Dt" // CHAR(0)
  call CheckPointReadRealDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, real_array, option)
  this%dt = real_array(1)

  dataset_name = "Prev_dt" // CHAR(0)
  call CheckPointReadRealDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, real_array, option)
  this%prev_dt = real_array(1)

  dataset_name = "Num_steps" // CHAR(0)
  call CheckPointReadIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)
  this%steps = int_array(1)

  dataset_name = "Cumulative_time_step_cuts" // CHAR(0)
  call CheckPointReadIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)
  this%cumulative_time_step_cuts = int_array(1)

  dataset_name = "Num_constant_time_steps" // CHAR(0)
  call CheckPointReadIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)
  this%num_constant_time_steps = int_array(1)

  dataset_name = "Num_contig_revert_due_to_sync" // CHAR(0)
  call CheckPointReadIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)
  this%num_contig_revert_due_to_sync = int_array(1)

  dataset_name = "Revert_dt" // CHAR(0)
  call CheckPointReadIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                     dataset_rank, dims, start, length, &
                                     stride, int_array, option)
  this%revert_dt = (int_array(1) == ONE_INTEGER)

  call h5gclose_f(timestepper_grp_id, hdf5_err)

  deallocate(start)
  deallocate(dims)
  deallocate(length)
  deallocate(stride)
  deallocate(int_array)
  deallocate(real_array)

end subroutine TimestepperSNESRestartHDF5

! ************************************************************************** !

subroutine TimestepperSNESGetHeader(this,header)
  !
  ! Gets values in checkpoint header.
  !
  ! Author: Glenn Hammond
  ! Date: 07/25/13
  !

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_SNES_type) :: this
  class(stepper_SNES_header_type) :: header

  this%cumulative_newton_iterations = header%cumulative_newton_iterations
  this%cumulative_linear_iterations = header%cumulative_linear_iterations
  this%num_newton_iterations = header%num_newton_iterations

  call TimestepperBaseGetHeader(this,header)

end subroutine TimestepperSNESGetHeader

! ************************************************************************** !

subroutine TimestepperSNESReset(this)
  !
  ! Zeros timestepper object members.
  !
  ! Author: Glenn Hammond
  ! Date: 01/20/14
  !

  implicit none

  class(timestepper_SNES_type) :: this

  this%cumulative_newton_iterations = 0
  this%cumulative_linear_iterations = 0
  this%num_newton_iterations = 0

  call TimestepperBaseReset(this)

end subroutine TimestepperSNESReset

! ************************************************************************** !

subroutine TimestepperSNESPrintInfo(this,option)
  !
  ! Prints settings for base timestepper.
  !
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  !
  use Option_module
  use String_module

  implicit none

  class(timestepper_SNES_type) :: this
  type(option_type) :: option

  PetscInt :: fids(2)
  PetscInt :: i
  character(len=MAXSTRINGLENGTH), allocatable :: strings(:)

  fids = OptionGetFIDs(option)
  call StringWriteToUnits(fids,' ')
  call StringWriteToUnits(fids,trim(this%name) // ' Time Stepper')

  ! have to allocate since ntfac can be infinite
  allocate(strings(this%ntfac+20))
  strings = ''
  strings(1) = 'acceleration: ' // &
                           StringWrite(String1Or2(this%iaccel>0,'on','off'))
  if (this%iaccel > 0) then
    strings(2) = 'acceleration threshold: ' // StringWrite(this%iaccel)
  endif
  strings(3) = 'number of ramp entries: ' // StringWrite(this%iaccel)
  do i = 1, this%ntfac
    strings(i+3) = 'ramp entry #' // trim(StringWrite(i)) // ': ' // &
                   StringWriteF(this%tfac(i))
  enddo
  call StringsCenter(strings,30,':')
  call StringWriteToUnits(fids,strings)
  deallocate(strings)

  call TimestepperBasePrintInfo(this,option)
  call SolverPrintNewtonInfo(this%solver,this%name,option)
  call SolverPrintLinearInfo(this%solver,this%name,option)

end subroutine TimestepperSNESPrintInfo

! ************************************************************************** !

subroutine TimestepperSNESInputRecord(this)
  !
  ! Prints information about the time stepper to the input record.
  ! To get a## format, must match that in simulation types.
  !
  ! Author: Jenn Frederick, SNL
  ! Date: 03/17/2016
  !

  implicit none

  class(timestepper_SNES_type) :: this

  PetscInt :: id
  character(len=MAXWORDLENGTH) :: word

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pmc timestepper: '
  write(id,'(a)') this%name

  write(id,'(a29)',advance='no') 'initial timestep size: '
  write(word,*) this%dt_init
  write(id,'(a)') trim(adjustl(word)) // ' sec'

end subroutine TimestepperSNESInputRecord

! ************************************************************************** !

recursive subroutine TimestepperSNESFinalizeRun(this,option)
  !
  ! Finalizes the time stepping
  !
  ! Author: Glenn Hammond
  ! Date: 07/22/13
  !

  use Option_module

  implicit none

  class(timestepper_SNES_type) :: this
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

#ifdef DEBUG
  call PrintMsg(option,'TimestepperSNESFinalizeRun()')
#endif

  if (OptionPrintToScreen(option)) then
    write(*,'(/,x,a," TS SNES steps = ",i6," newton = ",i8," linear = ",i10, &
            & " cuts = ",i6)') &
            trim(this%name), &
            this%steps, &
            this%cumulative_newton_iterations, &
            this%cumulative_linear_iterations, &
            this%cumulative_time_step_cuts
    write(string,'(i12)') this%cumulative_wasted_linear_iterations

    write(*,'(x,a)') trim(this%name) // &
      ' TS SNES Wasted Linear Iterations = ' // trim(adjustl(string))
    write(string,'(i12)') this%cumulative_wasted_newton_iterations
    write(*,'(x,a)') trim(this%name) // &
      ' TS SNES Wasted Newton Iterations = ' // trim(adjustl(string))

    write(string,'(f12.1)') this%cumulative_solver_time
    write(*,'(x,a)') trim(this%name) // ' TS SNES SNES time = ' // &
      trim(adjustl(string)) // ' seconds'
  endif

end subroutine TimestepperSNESFinalizeRun

! ************************************************************************** !

subroutine TimestepperSNESStrip(this)
  !
  ! Deallocates members of a time stepper
  !
  ! Author: Glenn Hammond
  ! Date: 07/22/13
  !

  implicit none

  class(timestepper_SNES_type) :: this

  call TimestepperBaseStrip(this)

  if (associated(this%tfac)) deallocate(this%tfac)
  nullify(this%tfac)

end subroutine TimestepperSNESStrip

! ************************************************************************** !

subroutine TimestepperSNESDestroy(this)
  !
  ! Deallocates a time stepper
  !
  ! Author: Glenn Hammond
  ! Date: 07/22/13
  !

  implicit none

  class(timestepper_SNES_type) :: this

  call TimestepperSNESStrip(this)

!  deallocate(this)
!  nullify(this)

end subroutine TimestepperSNESDestroy

end module Timestepper_SNES_class
