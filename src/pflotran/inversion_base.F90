module Inversion_Base_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Driver_class
  use Timer_class
  use Option_Inversion_module

  implicit none

  private

  type, public :: inversion_base_type
    class(driver_type), pointer :: driver
    class(timer_type), pointer :: timer
    type(inversion_option_type), pointer :: inversion_option
    PetscInt :: iteration        ! iteration number
    PetscInt :: maximum_iteration        ! Maximum iteration number
    PetscBool :: converged               ! convergence flag
  contains
    procedure, public :: Init => InversionBaseInit
    procedure, public :: ReadBlock => InversionBaseReadBlock
    procedure, public :: Step => InversionBaseStep
    procedure, public :: InitializeForwardRun => InversionBaseThisAndOption
    procedure, public :: SetupForwardRunLinkage => InversionBaseThisOnlyError
    procedure, public :: ConnectToForwardRun => InversionBaseThisOnlyError
    procedure, public :: ExecuteForwardRun => InversionBaseThisOnlyError
    procedure, public :: DestroyForwardRun => InversionBaseThisOnlyError
    procedure, public :: Checkpoint => InversionBaseThisOnly
    procedure, public :: RestartReadData => InversionBaseThisOnly
    procedure, public :: InitializeIterationNumber => &
                           InversionBaseInitIterationNum
    procedure, public :: IncrementIteration => InversionBaseIncrementIteration
    procedure, public :: EvaluateCostFunction => InversionBaseThisOnlyError
    procedure, public :: CheckConvergence => InversionBaseThisOnlyError
    procedure, public :: WriteIterationInfo => InversionBaseThisOnlyError
    procedure, public :: CalculateSensitivity => InversionBaseThisOnlyError
    procedure, public :: ScaleSensitivity => InversionBaseThisOnlyError
    procedure, public :: CalculateUpdate => InversionBaseThisOnlyError
    procedure, public :: UpdateParameters => InversionBaseThisOnlyError
    procedure, public :: UpdateRegularizationParameters => &
                           InversionBaseThisOnlyError
    procedure, public :: Finalize => InversionBaseFinalize
    procedure, public :: Strip => InversionBaseStrip
  end type inversion_base_type

  public :: InversionBaseInit, &
            InversionBaseReadSelectCase, &
            InversionBaseFinalize, &
            InversionBaseStrip, &
            InversionBaseDestroy

contains

! ************************************************************************** !

subroutine InversionBaseInit(this,driver)
  !
  ! Initializes a new inversion object
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_base_type) :: this
  class(driver_type), pointer :: driver

  this%driver => driver
  this%timer => TimerCreate()
  call this%timer%Start()
  this%inversion_option => OptionInversionCreate()

  this%iteration = 0
  this%maximum_iteration = UNINITIALIZED_INTEGER
  this%converged = PETSC_FALSE

end subroutine InversionBaseInit

! ************************************************************************** !

subroutine InversionBaseReadBlock(this,input,option)

  use Input_Aux_module
  use Option_module
  use String_module

  class(inversion_base_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  call this%driver%PrintErrMsg('InversionBaseReadBlock must be extended.')

end subroutine InversionBaseReadBlock

! ************************************************************************** !

subroutine InversionBaseReadSelectCase(this,input,keyword,found, &
                                       error_string,option)

  use Input_Aux_module
  use Option_module

  class(inversion_base_type) :: this
  type(input_type) :: input

  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option

  found = PETSC_TRUE
  select case(trim(keyword))
    case('MAX_INVERSION_ITERATION','MAXIMUM_NUMBER_OF_ITERATIONS')
      call InputReadInt(input,option,this%maximum_iteration)
      call InputErrorMsg(input,option,'MAXIMUM_NUMBER_OF_ITERATIONS', &
                         error_string)
    case default
      found = PETSC_FALSE
  end select

end subroutine InversionBaseReadSelectCase

! ************************************************************************** !

subroutine InversionBaseStep(this)
  !
  ! Performes a single inversion iteration (forward runs, Jacobians, update)
  !
  ! Author: Glenn Hammond
  ! Date: 03/21/22

  use Option_module

  class(inversion_base_type) :: this

  type(option_type), pointer :: option

  nullify(option)
  call this%InitializeForwardRun(option)
  call this%SetupForwardRunLinkage()
  call this%ConnectToForwardRun()
  call this%ExecuteForwardRun()
  call this%CheckConvergence()
  call this%WriteIterationInfo()
  call this%Checkpoint()
  if (.not.this%converged) then
    call this%CalculateSensitivity()
    call this%ScaleSensitivity()
    call this%CalculateUpdate()
    call this%UpdateRegularizationParameters()
  endif
  call this%DestroyForwardRun()

  this%converged = PETSC_FALSE
  if (this%iteration > this%maximum_iteration) this%converged = PETSC_TRUE

end subroutine InversionBaseStep

! ************************************************************************** !

subroutine InversionBaseThisAndOption(this,option)

  use Option_module

  class(inversion_base_type) :: this
  type(option_type), pointer :: option

end subroutine InversionBaseThisAndOption

! ************************************************************************** !

subroutine InversionBaseThisOnly(this)

  class(inversion_base_type) :: this

end subroutine InversionBaseThisOnly

! ************************************************************************** !

subroutine InversionBaseThisOnlyError(this)

  class(inversion_base_type) :: this

  call this%driver%PrintErrMsg('An inversion routine that passes only "this" &
    &must be extended from the Base implementation.')

end subroutine InversionBaseThisOnlyError

! ************************************************************************** !

subroutine InversionBaseInitIterationNum(this)
  !
  ! Sets starting iteration number
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 07/09/21

  class(inversion_base_type) :: this

  this%iteration = 1

end subroutine InversionBaseInitIterationNum

! ************************************************************************** !

subroutine InversionBaseIncrementIteration(this)
  !
  ! Sets starting iteration number
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 07/09/21

  class(inversion_base_type) :: this

  this%iteration = this%iteration + 1

end subroutine InversionBaseIncrementIteration

! ************************************************************************** !

subroutine InversionBaseOutputSensitivity(this,suffix)

  class(inversion_base_type) :: this
  character(len=*) :: suffix

  call this%driver%PrintErrMsg('InversionBaseOutputSensitivity &
    &must be extended from the Base implementation.')

end subroutine InversionBaseOutputSensitivity

! ************************************************************************** !

subroutine InversionBaseFinalize(this)
  !
  ! Finalizes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_base_type) :: this

  PetscLogDouble :: total_time

  call this%timer%Stop()
  total_time = this%timer%GetCumulativeTime()

end subroutine InversionBaseFinalize

! ************************************************************************** !

subroutine InversionBaseStrip(this)
  !
  ! Deallocates members of inversion base
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_base_type) :: this

  nullify(this%driver)
  call TimerDestroy(this%timer)
  call OptionInversionDestroy(this%inversion_option)

end subroutine InversionBaseStrip

! ************************************************************************** !

subroutine InversionBaseDestroy(inversion)
  !
  ! Deallocates a inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_base_type), pointer :: inversion

  if (.not.associated(inversion)) return

  call inversion%Strip()
  deallocate(inversion)
  nullify(inversion)

end subroutine InversionBaseDestroy

end module Inversion_Base_class
