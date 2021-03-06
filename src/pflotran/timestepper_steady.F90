module Timestepper_Steady_class

#include "petsc/finclude/petscsys.h"
  use petscsys
 
  use Timestepper_BE_class
  use Timestepper_Base_class
  use Convergence_module
  use Solver_module
  use Waypoint_module
  use PFLOTRAN_Constants_module

  implicit none

  type, public, extends(timestepper_BE_type) :: timestepper_steady_type
  
  contains

 !   procedure, public :: Init => TimestepperSteadyInit
    procedure, public :: StepDT => TimestepperSteadyStepDT
    procedure, public :: UpdateDT => TimestepperSteadyUpdateDT
    procedure, public :: InputRecord => TimestepperSteadyInputRecord

  end type timestepper_steady_type

  public :: TimestepperSteadyCreate, &
            TimestepperSteadyCreateFromBE, &
            TimestepperSteadyCreateFromBase

contains

! ************************************************************************** !

function TimestepperSteadyCreate()
  ! 
  ! This routine creates timestepper for steady solve
  ! 
  ! Author: Gautam Bisht, LBNL and Satish Karra, LANL
  ! Date: 01/01/14, 04/07/2015
  ! 

  implicit none

  class(timestepper_steady_type), pointer :: TimestepperSteadyCreate

  class(timestepper_steady_type), pointer :: stepper

  allocate(stepper)
  call stepper%Init()

  TimestepperSteadyCreate => stepper

end function TimestepperSteadyCreate

! ************************************************************************** !

function TimestepperSteadyCreateFromBase(timestepper_base)
  ! 
  ! This routine creates timestepper for steady solve
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 06/29/18
  ! 

  implicit none

  class(timestepper_base_type) :: timestepper_base

  class(timestepper_steady_type), pointer :: TimestepperSteadyCreateFromBase

  class(timestepper_steady_type), pointer :: stepper

  allocate(stepper)
  call stepper%Init()

  stepper%name = timestepper_base%name
  stepper%steps = timestepper_base%steps
  stepper%num_constant_time_steps = timestepper_base%num_constant_time_steps

  stepper%max_time_step = timestepper_base%max_time_step
  stepper%max_time_step_cuts = timestepper_base%max_time_step_cuts
  stepper%constant_time_step_threshold = timestepper_base%constant_time_step_threshold

  stepper%cumulative_time_step_cuts = timestepper_base%cumulative_time_step_cuts
  stepper%cumulative_solver_time = timestepper_base%cumulative_solver_time

  stepper%dt = timestepper_base%dt
  stepper%prev_dt = timestepper_base%prev_dt
  stepper%dt_init = timestepper_base%dt_init
  stepper%dt_max = timestepper_base%dt_max
  stepper%revert_dt = timestepper_base%revert_dt
  stepper%num_contig_revert_due_to_sync = timestepper_base%num_contig_revert_due_to_sync

  stepper%init_to_steady_state = timestepper_base%init_to_steady_state
  stepper%run_as_steady_state = timestepper_base%run_as_steady_state
  stepper%steady_state_rel_tol = timestepper_base%steady_state_rel_tol

  stepper%time_step_cut_flag = timestepper_base%time_step_cut_flag

  stepper%start_time = timestepper_base%start_time
  stepper%start_time_step = timestepper_base%start_time_step
  stepper%time_step_tolerance = timestepper_base%time_step_tolerance
  stepper%target_time = timestepper_base%target_time

  stepper%cur_waypoint => timestepper_base%cur_waypoint
  stepper%prev_waypoint => timestepper_base%prev_waypoint

  stepper%solver => timestepper_base%solver

  TimestepperSteadyCreateFromBase => stepper

end function TimestepperSteadyCreateFromBase

! ************************************************************************** !

function TimestepperSteadyCreateFromBE(timestepper_BE)
  ! 
  ! This routine creates timestepper for steady solve
  ! 
  ! Author: Gautam Bisht, LBNL and Satish Karra, LANL
  ! Date: 01/01/14, 04/07/2015
  ! 

  implicit none

  class(timestepper_BE_type) :: timestepper_BE

  class(timestepper_steady_type), pointer :: TimestepperSteadyCreateFromBE

  class(timestepper_steady_type), pointer :: stepper
  PetscInt :: i

  allocate(stepper)
  call stepper%Init()
  
  stepper%name = timestepper_BE%name
  stepper%steps = timestepper_BE%steps
  stepper%num_constant_time_steps = timestepper_BE%num_constant_time_steps

  stepper%max_time_step = timestepper_BE%max_time_step
  stepper%max_time_step_cuts = timestepper_BE%max_time_step_cuts
  stepper%constant_time_step_threshold = timestepper_BE%constant_time_step_threshold

  stepper%cumulative_time_step_cuts = timestepper_BE%cumulative_time_step_cuts    
  stepper%cumulative_solver_time = timestepper_BE%cumulative_solver_time

  stepper%start_time = timestepper_BE%start_time
  stepper%start_time_step = timestepper_BE%start_time_step
  stepper%time_step_tolerance = timestepper_BE%time_step_tolerance
  stepper%target_time = timestepper_BE%target_time
  
  stepper%prev_dt = timestepper_BE%prev_dt
  stepper%dt = timestepper_BE%dt
  stepper%dt_init = timestepper_BE%dt_init
  stepper%dt_max = timestepper_BE%dt_max
  
  stepper%time_step_cut_flag = timestepper_BE%time_step_cut_flag
  
  stepper%init_to_steady_state = timestepper_BE%init_to_steady_state
  stepper%steady_state_rel_tol = timestepper_BE%steady_state_rel_tol
  stepper%run_as_steady_state = timestepper_BE%run_as_steady_state
  
  stepper%cur_waypoint => timestepper_BE%cur_waypoint
  stepper%prev_waypoint => timestepper_BE%prev_waypoint
  stepper%revert_dt = timestepper_BE%revert_dt
  stepper%num_contig_revert_due_to_sync = timestepper_BE%num_contig_revert_due_to_sync
      
  stepper%num_newton_iterations = timestepper_BE%num_newton_iterations
  stepper%num_linear_iterations = timestepper_BE%num_linear_iterations

  stepper%cumulative_newton_iterations = timestepper_BE%cumulative_newton_iterations
  stepper%cumulative_linear_iterations = timestepper_BE%cumulative_linear_iterations

  stepper%iaccel = timestepper_BE%iaccel
  stepper%ntfac = timestepper_BE%ntfac
  allocate(stepper%tfac(timestepper_BE%ntfac))
  do i = 1, timestepper_BE%ntfac
    stepper%tfac(i) = timestepper_BE%tfac(i)
  enddo

  stepper%solver => timestepper_BE%solver

  TimestepperSteadyCreateFromBE => stepper

end function TimestepperSteadyCreateFromBE

! ************************************************************************** !

subroutine TimestepperSteadyInit(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL and Satish Karra, LANL
  ! Date: 01/01/14, 04/07/2015
  ! 

  implicit none
  
  class(timestepper_steady_type) :: this

  call TimestepperBaseInit(this)

end subroutine TimestepperSteadyInit


! ************************************************************************** !

subroutine TimestepperSteadyUpdateDT(this,process_model)
  ! 
  ! Updates time step
  ! 
  ! Author: Gautam Bisht, LBNL and Satish Karra, LANL
  ! Date: 01/01/14, 04/07/2015
  ! 

  use PM_Base_class
  use Option_module
  
  implicit none

  class(timestepper_steady_type) :: this
  class(pm_base_type) :: process_model
  

end subroutine TimestepperSteadyUpdateDT

! ************************************************************************** !

subroutine TimestepperSteadyStepDT(this, process_model, stop_flag)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL and Satish Karra, LANL
  ! Date: 01/01/14, 04/07/2015
  ! 

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use Option_module
  use PM_Geomechanics_Force_class
  use Output_module, only : Output

  use Solver_module

  implicit none

  class(timestepper_steady_type) :: this
  class(pm_base_type) :: process_model
  PetscInt :: stop_flag

  PetscBool :: failure
  PetscLogDouble :: log_start_time
  PetscLogDouble :: log_end_time
  PetscErrorCode :: ierr
  PetscInt :: sum_newton_iterations, sum_linear_iterations
  PetscInt :: num_newton_iterations, num_linear_iterations
  PetscInt :: snes_reason
  PetscInt :: icut
  PetscReal :: fnorm
  PetscReal :: inorm
  PetscReal :: scaled_fnorm
  Vec :: residual_vec

  type(option_type), pointer :: option
  type(solver_type), pointer :: solver

  solver => process_model%solver
  option => process_model%option

  sum_newton_iterations = 0
  sum_linear_iterations = 0
  icut = 0

  call process_model%InitializeTimestep()

  call process_model%PreSolve()

  call PetscTime(log_start_time, ierr);CHKERRQ(ierr)

  call SNESSolve(solver%snes, PETSC_NULL_VEC, &
                 process_model%solution_vec, ierr);CHKERRQ(ierr)

  call PetscTime(log_end_time, ierr);CHKERRQ(ierr)

  this%cumulative_solver_time = &
      this%cumulative_solver_time + &
      (log_end_time - log_start_time)
     
  call SNESGetIterationNumber(solver%snes, num_newton_iterations,  &
                              ierr);CHKERRQ(ierr)
  call SNESGetLinearSolveIterations(solver%snes, num_linear_iterations,  &
                                    ierr);CHKERRQ(ierr)
  call SNESGetConvergedReason(solver%snes, snes_reason, ierr);CHKERRQ(ierr)

  sum_newton_iterations = sum_newton_iterations + num_newton_iterations
  sum_linear_iterations = sum_linear_iterations + num_linear_iterations

  if (snes_reason <= 0) then
    if (option%print_screen_flag) then
      print *, 'Newton solver failed to converge in steady-solve, reason: ', &
                snes_reason
    endif
    failure = PETSC_TRUE
    stop_flag = TS_STOP_END_SIMULATION
    return
  endif
  
  this%steps = this%steps + 1
  this%cumulative_newton_iterations = &
    this%cumulative_newton_iterations + sum_newton_iterations
  this%cumulative_linear_iterations = &
    this%cumulative_linear_iterations + sum_linear_iterations
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
    select type(pm => process_model)
      class is(pm_geomech_force_type)
        if (associated(pm%geomech_realization%geomech_discretization%grid)) then
          scaled_fnorm = fnorm/pm%geomech_realization%geomech_discretization% &
                          grid%nmax_node
        else
          scaled_fnorm = fnorm
        endif
    end select
    write(*,*) ''
    print *,' --> SNES Linear/Non-Linear Iterations = ', &
             num_linear_iterations,' / ',num_newton_iterations
    write(*,'(" --> SNES Residual: ",1p3e14.6)') fnorm, scaled_fnorm, inorm 
  endif
  
  if (option%print_screen_flag) print *, ""
  
  if (option%print_file_flag) then
    write(option%fid_out, '(" STEADY-SOLVE ",i6," snes_conv_reason: ",i4,/, &
      &"  newton = ",i3," [",i8,"]", &
      & " linear = ",i5," [",i10,"]")') &
      this%steps, &
      snes_reason,sum_newton_iterations, &
      this%cumulative_newton_iterations,sum_linear_iterations, &
      this%cumulative_linear_iterations
  endif  

  option%time = this%target_time
  call process_model%FinalizeTimestep()
  
  if (option%print_screen_flag) print *, ""
  ! check if the steady state option is selected in input deck
  ! if yes then end simulation
  if (option%steady_state) stop_flag = TS_STOP_END_SIMULATION

end subroutine TimestepperSteadyStepDT

! ************************************************************************** !

subroutine TimestepperSteadyInputRecord(this)
  ! 
  ! Prints information about the time stepper to the input record.
  ! To get a## format, must match that in simulation types.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/17/2016
  ! 
  
  implicit none
  
  class(timestepper_steady_type) :: this

  PetscInt :: id
  character(len=MAXWORDLENGTH) :: word
   
  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pmc timestepper: '
  write(id,'(a)') this%name

end subroutine TimestepperSteadyInputRecord

end module Timestepper_Steady_class
