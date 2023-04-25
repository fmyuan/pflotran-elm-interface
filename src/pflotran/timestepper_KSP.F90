module Timestepper_KSP_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Solver_module
  use Convergence_module
  use Timestepper_Base_class

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(timestepper_base_type) :: timestepper_KSP_type

    PetscInt :: num_linear_iterations ! number of linear solver iterations in a time step
    PetscInt :: cumulative_linear_iterations     ! Total number of linear iterations
    PetscInt :: cumulative_wasted_linear_iterations
  contains

    procedure, public :: ReadSelectCase => TimestepperKSPReadSelectCase
    procedure, public :: Init => TimestepperKSPInit
    procedure, public :: UpdateDT => TimestepperKSPUpdateDT
    procedure, public :: StepDT => TimestepperKSPStepDT
    procedure, public :: CheckpointBinary => TimestepperKSPCheckpointBinary
    procedure, public :: RestartBinary => TimestepperKSPRestartBinary
    procedure, public :: CheckpointHDF5 => TimestepperKSPCheckpointHDF5
    procedure, public :: RestartHDF5 => TimestepperKSPRestartHDF5
    procedure, public :: Reset => TimestepperKSPReset
    procedure, public :: PrintInfo => TimestepperKSPPrintInfo
    procedure, public :: InputRecord => TimestepperKSPInputRecord
    procedure, public :: FinalizeRun => TimestepperKSPFinalizeRun
    procedure, public :: Strip => TimestepperKSPStrip
    procedure, public :: Destroy => TimestepperKSPDestroy

  end type timestepper_KSP_type

  ! For checkpointing
  type, public, extends(stepper_base_header_type) :: stepper_KSP_header_type
    PetscInt :: cumulative_linear_iterations
  end type stepper_KSP_header_type

  interface PetscBagGetData
    subroutine PetscBagGetData(bag,header,ierr)
      import :: stepper_KSP_header_type
      implicit none
#include "petsc/finclude/petscbag.h"
      PetscBag :: bag
      class(stepper_KSP_header_type), pointer :: header
      PetscErrorCode :: ierr
    end subroutine
  end interface PetscBagGetData

  public :: TimestepperKSPCreate, &
            TimestepperKSPPrintInfo, &
            TimestepperKSPInit

contains

! ************************************************************************** !

function TimestepperKSPCreate()
  !
  ! Allocates and initializes a new Timestepper object
  !
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  !

  implicit none

  class(timestepper_KSP_type), pointer :: TimestepperKSPCreate

  class(timestepper_KSP_type), pointer :: stepper

  allocate(stepper)
  call stepper%Init()

  stepper%solver => SolverCreate()

  TimestepperKSPCreate => stepper

end function TimestepperKSPCreate

! ************************************************************************** !

subroutine TimestepperKSPInit(this)
  !
  ! Allocates and initializes a new Timestepper object
  !
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  !

  implicit none

  class(timestepper_KSP_type) :: this

  call TimestepperBaseInit(this)

  this%num_linear_iterations = 0

  this%cumulative_linear_iterations = 0
  this%cumulative_wasted_linear_iterations = 0

  nullify(this%solver)

end subroutine TimestepperKSPInit

! ************************************************************************** !

subroutine TimestepperKSPReadSelectCase(this,input,keyword,found, &
                                        error_string,option)
  !
  ! Updates time step
  !
  ! Author: Glenn Hammond
  ! Date: 03/20/13
  !

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module

  implicit none

  class(timestepper_KSP_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option

  found = PETSC_TRUE
  call TimestepperBaseReadSelectCase(this,input,keyword,found, &
                                     error_string,option)
  if (found) return

  found = PETSC_TRUE
  select case(trim(keyword))
    case default
      found = PETSC_FALSE
  end select

end subroutine TimestepperKSPReadSelectCase

! ************************************************************************** !

subroutine TimestepperKSPUpdateDT(this,process_model)
  !
  ! Updates time step
  !
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  !

  use PM_Base_class

  implicit none

  class(timestepper_KSP_type) :: this
  class(pm_base_type) :: process_model

  PetscBool :: update_time_step
  PetscInt :: dummy_int
  PetscReal :: dummy_array(1)

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
                                    dummy_int, dummy_int, dummy_array, &
                                    this%time_step_max_growth_factor)

end subroutine TimestepperKSPUpdateDT

! ************************************************************************** !

subroutine TimestepperKSPStepDT(this,process_model,stop_flag)
  !
  ! Steps forward one step in time
  !
  ! Author: Glenn Hammond
  ! Date: 09/03/21
  !

#include "petsc/finclude/petscksp.h"
  use petscksp
  use PM_Base_class
  use Option_module
  use Output_EKG_module, only : IUNIT_EKG

  implicit none

  class(timestepper_KSP_type) :: this
  class(pm_base_type) :: process_model
  PetscInt :: stop_flag


  type(solver_type), pointer :: solver
  type(option_type), pointer :: option

  PetscLogDouble :: log_start_time
  PetscLogDouble :: log_end_time
  PetscInt :: num_linear_iterations
  PetscInt :: sum_linear_iterations
  PetscInt :: sum_wasted_linear_iterations
  PetscInt :: icut
  character(len=MAXWORDLENGTH) :: tunit
  PetscReal :: tconv

  Vec :: rhs
  Vec :: solution
  KSPConvergedReason :: ksp_reason

  PetscErrorCode :: ierr

  solver => process_model%solver
  option => process_model%option

  tconv = process_model%output_option%tconv
  tunit = process_model%output_option%tunit
  sum_linear_iterations = 0
  sum_wasted_linear_iterations = 0
  icut = 0

  option%dt = this%dt
  option%time = this%target_time-this%dt

  call process_model%InitializeTimestep()

  do

    ! these pointers are set within SetupLinearSystem()
    solution = PETSC_NULL_VEC
    rhs = PETSC_NULL_VEC
    call process_model%SetupLinearSystem(solver%M,solution,rhs,ierr)

    call PetscTime(log_start_time,ierr);CHKERRQ(ierr)

    call KSPSetOperators(solver%ksp,solver%M,solver%Mpre,ierr);CHKERRQ(ierr)
    call KSPSolve(solver%ksp,rhs,solution,ierr);CHKERRQ(ierr)

    call PetscTime(log_end_time,ierr);CHKERRQ(ierr)

    this%cumulative_solver_time = &
      this%cumulative_solver_time + &
      (log_end_time - log_start_time)

    call KSPGetIterationNumber(solver%ksp,num_linear_iterations, &
                               ierr);CHKERRQ(ierr)
    call KSPGetConvergedReason(solver%ksp,ksp_reason,ierr);CHKERRQ(ierr)

    sum_linear_iterations = sum_linear_iterations + num_linear_iterations

    if (ksp_reason <= 0 .or. .not. process_model%AcceptSolution()) then
      sum_wasted_linear_iterations = sum_wasted_linear_iterations + &
           num_linear_iterations
      ! The linear solver diverged, so try reducing the time step.
      call this%CutDT(process_model,icut,stop_flag,'ksp', &
                      ksp_reason,option)
      if (ksp_reason < 0) then
        call SolverLinearPrintFailedReason(solver,option)
        if (solver%verbose_logging) then
          ! add any verbose logging (see timestepper_BE)
        endif
      endif
      if (stop_flag == TS_STOP_FAILURE) return

      this%target_time = this%target_time + this%dt
      option%dt = this%dt
      call process_model%TimeCut()

    else
      ! The Newton solver converged, so we can exit.
      exit
    endif
  enddo

  this%steps = this%steps + 1
  this%cumulative_linear_iterations = &
    this%cumulative_linear_iterations + sum_linear_iterations
  this%cumulative_wasted_linear_iterations = &
       this%cumulative_wasted_linear_iterations + sum_wasted_linear_iterations
  this%cumulative_time_step_cuts = &
    this%cumulative_time_step_cuts + icut

  this%num_linear_iterations = num_linear_iterations

  call process_model%PostSolve()

  call TimestepperBasePrintStepInfo(this,process_model%output_option, &
                                    ksp_reason,option)
  write(option%io_buffer,'("  linear = ",i5, &
                         &" [",i10,"]"," cuts = ",i2," [",i4,"]")') &
           sum_linear_iterations,this%cumulative_linear_iterations, &
           icut,this%cumulative_time_step_cuts
  call PrintMsg(option)

  option%time = this%target_time
  call process_model%FinalizeTimestep()

  if (this%print_ekg .and. OptionPrintToFile(option)) then
100 format(a32," TIMESTEP ",i10,2es16.8,a,i3,i5,i5,i10)
    write(IUNIT_EKG,100) trim(this%name), this%steps, &
      this%target_time/process_model%output_option%tconv, &
      this%dt/process_model%output_option%tconv, &
      trim(process_model%output_option%tunit), &
      icut, this%cumulative_time_step_cuts, &
      sum_linear_iterations, this%cumulative_linear_iterations
  endif

end subroutine TimestepperKSPStepDT

! ************************************************************************** !

subroutine TimestepperKSPCheckpointBinary(this,viewer,option)
  !
  ! Checkpoints parameters/variables associated with
  ! a time stepper.
  !
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  !
  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_KSP_type) :: this
  PetscViewer :: viewer
  type(option_type) :: option

  class(stepper_KSP_header_type), pointer :: header
  type(stepper_KSP_header_type) :: dummy_header
  character(len=1),pointer :: dummy_char(:)
  PetscBag :: bag
  PetscSizeT :: bagsize
  PetscErrorCode :: ierr

  bagsize = size(transfer(dummy_header,dummy_char))

  call PetscBagCreate(option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
  call TimestepperKSPRegisterHeader(this,bag,header)
  call TimestepperKSPSetHeader(this,bag,header)
  call PetscBagView(bag,viewer,ierr);CHKERRQ(ierr)
  call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)

end subroutine TimestepperKSPCheckpointBinary

! ************************************************************************** !

subroutine TimestepperKSPRegisterHeader(this,bag,header)
  !
  ! Register header entries.
  !
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  !

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_KSP_type) :: this
  class(stepper_KSP_header_type) :: header
  PetscBag :: bag

  PetscErrorCode :: ierr

  call PetscBagRegisterInt(bag,header%cumulative_linear_iterations,0, &
                           "cumulative_linear_iterations","", &
                           ierr);CHKERRQ(ierr)
! need to add cumulative wasted linear iterations

  call TimestepperBaseRegisterHeader(this,bag,header)

end subroutine TimestepperKSPRegisterHeader

! ************************************************************************** !

subroutine TimestepperKSPSetHeader(this,bag,header)
  !
  ! Sets values in checkpoint header.
  !
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  !

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_KSP_type) :: this
  class(stepper_KSP_header_type) :: header
  PetscBag :: bag

  header%cumulative_linear_iterations = this%cumulative_linear_iterations

  call TimestepperBaseSetHeader(this,bag,header)

end subroutine TimestepperKSPSetHeader

! ************************************************************************** !

subroutine TimestepperKSPRestartBinary(this,viewer,option)
  !
  ! Checkpoints parameters/variables associated with
  ! a time stepper.
  !
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  !

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_KSP_type) :: this
  PetscViewer :: viewer
  type(option_type) :: option

  class(stepper_KSP_header_type), pointer :: header
  type(stepper_KSP_header_type) :: dummy_header
  character(len=1),pointer :: dummy_char(:)
  PetscBag :: bag
  PetscSizeT :: bagsize
  PetscErrorCode :: ierr

  bagsize = size(transfer(dummy_header,dummy_char))

  call PetscBagCreate(option%mycomm,bagsize,bag,ierr);CHKERRQ(ierr)
  call PetscBagGetData(bag,header,ierr);CHKERRQ(ierr)
  call TimestepperKSPRegisterHeader(this,bag,header)
  call PetscBagLoad(viewer,bag,ierr);CHKERRQ(ierr)
  call TimestepperKSPGetHeader(this,header)
  call PetscBagDestroy(bag,ierr);CHKERRQ(ierr)

end subroutine TimestepperKSPRestartBinary

! ************************************************************************** !

subroutine TimestepperKSPCheckpointHDF5(this, h5_chk_grp_id, option)
  !
  ! Checkpoints parameters/variables associated with
  ! a time stepper.
  !
  ! Author: Gautam Bisht
  ! Date: 12/06/19
  !
  use Option_module
  use hdf5
  use HDF5_Aux_module
  use Checkpoint_module, only : CheckPointWriteIntDatasetHDF5
  use Checkpoint_module, only : CheckPointWriteRealDatasetHDF5

  implicit none

  class(timestepper_KSP_type) :: this
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

  string = "Timestepper"
  call HDF5GroupCreate(h5_chk_grp_id, string, timestepper_grp_id, option)

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

  dataset_name = "Cumulative_linear_iterations" // CHAR(0)
  int_array(1) = this%cumulative_linear_iterations
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

  call HDF5GroupClose(timestepper_grp_id,option)

  deallocate(start)
  deallocate(dims)
  deallocate(length)
  deallocate(stride)
  deallocate(int_array)
  deallocate(real_array)

end subroutine TimestepperKSPCheckpointHDF5

! ************************************************************************** !

subroutine TimestepperKSPRestartHDF5(this, h5_chk_grp_id, option)
  !
  ! Restarts parameters/variables associated with
  ! a time stepper.
  !
  ! Author: Gautam Bisht
  ! Date: 12/06/19
  !
  use Option_module
  use hdf5
  use HDF5_Aux_module
  use Checkpoint_module, only : CheckPointReadIntDatasetHDF5
  use Checkpoint_module, only : CheckPointReadRealDatasetHDF5

  implicit none

  class(timestepper_KSP_type) :: this
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

  dataset_name = "Cumulative_linear_iterations" // CHAR(0)
  call CheckPointReadIntDatasetHDF5(timestepper_grp_id, dataset_name, &
                                    dataset_rank, dims, start, length, &
                                    stride, int_array, option)
  this%cumulative_linear_iterations = int_array(1)

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

  call HDF5GroupClose(timestepper_grp_id,option)

  deallocate(start)
  deallocate(dims)
  deallocate(length)
  deallocate(stride)
  deallocate(int_array)
  deallocate(real_array)

end subroutine TimestepperKSPRestartHDF5

! ************************************************************************** !

subroutine TimestepperKSPGetHeader(this,header)
  !
  ! Gets values in checkpoint header.
  !
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  !

  use Option_module

  implicit none

#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscbag.h"

  class(timestepper_KSP_type) :: this
  class(stepper_KSP_header_type) :: header

  this%cumulative_linear_iterations = header%cumulative_linear_iterations

  call TimestepperBaseGetHeader(this,header)

end subroutine TimestepperKSPGetHeader

! ************************************************************************** !

subroutine TimestepperKSPReset(this)
  !
  ! Zeros timestepper object members.
  !
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  !

  implicit none

  class(timestepper_KSP_type) :: this

  this%cumulative_linear_iterations = 0

  call TimestepperBaseReset(this)

end subroutine TimestepperKSPReset

! ************************************************************************** !

subroutine TimestepperKSPPrintInfo(this,aux_string,option)
  !
  ! Prints settings for base timestepper.
  !
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  !
  use Option_module
  use String_module

  implicit none

  class(timestepper_KSP_type) :: this
  character(len=*) :: aux_string
  type(option_type) :: option

  call TimestepperBasePrintInfo(this,'KSP',option)
  call SolverPrintNewtonInfo(this%solver,this%name,option)
  call SolverPrintLinearInfo(this%solver,this%name,option)

end subroutine TimestepperKSPPrintInfo

! ************************************************************************** !

subroutine TimestepperKSPInputRecord(this)
  !
  ! Prints information about the time stepper to the input record.
  ! To get a## format, must match that in simulation types.
  !
  ! Author: Jenn Frederick, SNL
  ! Date: 12/06/19
  !

  implicit none

  class(timestepper_KSP_type) :: this

  PetscInt :: id
  character(len=MAXWORDLENGTH) :: word

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pmc timestepper: '
  write(id,'(a)') this%name

  write(id,'(a29)',advance='no') 'initial timestep size: '
  write(word,*) this%dt_init
  write(id,'(a)') trim(adjustl(word)) // ' sec'

end subroutine TimestepperKSPInputRecord

! ************************************************************************** !

recursive subroutine TimestepperKSPFinalizeRun(this,option)
  !
  ! Finalizes the time stepping
  !
  ! Author: Glenn Hammond
  ! Date: 12/06/19

  use Option_module
  use String_module

  implicit none

  class(timestepper_KSP_type) :: this
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

#ifdef DEBUG
  call PrintMsg(option,'TimestepperKSPFinalizeRun()')
#endif

  call PrintMsg(option,'')
  string = ' ' // trim(this%name) // &
    ' TS KSP steps = ' // trim(StringWrite(this%steps)) // &
    ' linear = ' // trim(StringWrite(this%cumulative_linear_iterations)) // &
    ' cuts = ' // trim(StringWrite(this%cumulative_time_step_cuts))
  call PrintMsg(option,string)
  string = ' ' // trim(this%name) // &
    ' TS KSP Wasted Linear Iterations = ' // &
    trim(StringWrite(this%cumulative_wasted_linear_iterations))
  call PrintMsg(option,string)
  string = ' ' // trim(this%name) // &
    ' TS KSP time = ' // &
    trim(StringWrite('(f12.1)',this%cumulative_solver_time)) // ' seconds'
  call PrintMsg(option,string)

end subroutine TimestepperKSPFinalizeRun

! ************************************************************************** !

subroutine TimestepperKSPStrip(this)
  !
  ! Deallocates members of a time stepper
  !
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  !

  implicit none

  class(timestepper_KSP_type) :: this

  call TimestepperBaseStrip(this)

end subroutine TimestepperKSPStrip

! ************************************************************************** !

subroutine TimestepperKSPDestroy(this)
  !
  ! Deallocates a time stepper
  !
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  !

  implicit none

  class(timestepper_KSP_type) :: this

  call TimestepperKSPStrip(this)

end subroutine TimestepperKSPDestroy

end module Timestepper_KSP_class
