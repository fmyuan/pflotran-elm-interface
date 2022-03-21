module Inversion_Base_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Driver_module
  use Timer_class
  use Option_Inversion_module

  implicit none

  private

  type, public :: inversion_base_type
    class(driver_type), pointer :: driver
    class(timer_type), pointer :: timer
    type(inversion_option_type), pointer :: inversion_option
    PetscInt :: iteration                ! iteration number
    PetscInt :: maximum_iteration        ! Maximum iteration number
    PetscBool :: converg_flag            ! convergence flag
  contains
    procedure, public :: Init => InversionBaseInit
    procedure, public :: Initialize => InversionBaseThisOnly
    procedure, public :: ReadBlock => InversionBaseReadBlock
    procedure, public :: Step => InversionBaseThisOnly
    procedure, public :: InitializeForwardRun => InversionBaseThisAndOption
    procedure, public :: ConnectToFowardRun => InversionBaseThisOnly
    procedure, public :: ExecuteForwardRun => InversionBaseThisOnly
    procedure, public :: DestroyForwardRun => InversionBaseThisOnly
    procedure, public :: InitializeIterationNumber => &
                           InversionBaseInitIterationNum
    procedure, public :: IncrementIteration => InversionBaseIncrementIteration
    procedure, public :: UpdateParameters => InversionBaseThisOnly
    procedure, public :: CalculateUpdate => InversionBaseThisOnly
    procedure, public :: CalculateSensitivity => InversionBaseThisOnly
    procedure, public :: OutputSensitivity => InversionBaseOutputSensitivity
    procedure, public :: Invert => InversionBaseThisOnly
    procedure, public :: CheckConvergence => InversionBaseThisOnly
    procedure, public :: EvaluateCostFunction => InversionBaseThisOnly
    procedure, public :: UpdateRegularizParameters => InversionBaseThisOnly
    procedure, public :: WriteIterationInfo => InversionBaseThisOnly
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
  this%converg_flag = PETSC_FALSE

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

subroutine InversionBaseThisAndOption(this,option)

  use Option_module

  class(inversion_base_type) :: this
  type(option_type), pointer :: option

end subroutine InversionBaseThisAndOption

! ************************************************************************** !

subroutine InversionBaseThisOnly(this)

  class(inversion_base_type) :: this

  call this%driver%PrintErrMsg('An inversion routine that passes only "this" &
    &must be extended from the Base implementation.')

end subroutine InversionBaseThisOnly

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
