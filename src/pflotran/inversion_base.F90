module Inversion_Base_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Driver_module
  use Timer_class

  implicit none

  private

  type, public :: inversion_base_type
    class(driver_type), pointer :: driver
    class(timer_type), pointer :: timer
    PetscInt :: iteration                ! iteration number
    PetscBool :: converg_flag            ! convergence flag
  contains
    procedure, public :: Init => InversionBaseInit
    procedure, public :: Initialize => InversionBaseInitialize
    procedure, public :: ReadBlock => InversionBaseReadBlock
    procedure, public :: Step => InversionBaseStep
    procedure, public :: SetIteration => InversionBaseSetIteration
    procedure, public :: IncrementIteration => InversionBaseIncrementIteration
    procedure, public :: UpdateParameters => InversionBaseUpdateParameters
    procedure, public :: CalculateUpdate => InversionBaseCalculateUpdate
    procedure, public :: CheckConvergence => InversionBaseCheckConvergence
    procedure, public :: EvaluateCostFunction => InvBaseEvaluateCostFunction
    procedure, public :: UpdateRegularizationParameters => &
                           InvBaseUpdateRegularizParams
    procedure, public :: WriteIterationInfo => InversionBaseWriteIterationInfo
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

  this%iteration = 1
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

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  error_string = 'Base Inversion'

  input%ierr = 0
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    found = PETSC_TRUE
    call InversionBaseReadSelectCase(this,input,keyword,found, &
                                     error_string,option)
    if (found) cycle

    select case(trim(keyword))
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select

  enddo
  call InputPopBlock(input,option)

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
    case default
      found = PETSC_FALSE
  end select

end subroutine InversionBaseReadSelectCase

! ************************************************************************** !

subroutine InversionBaseInitialize(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_base_type) :: this

end subroutine InversionBaseInitialize


! ************************************************************************** !

subroutine InversionBaseSetIteration(this,i)
  !
  ! Sets starting iteration number
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 07/09/21

  class(inversion_base_type) :: this
  PetscInt :: i

  this%iteration = i

end subroutine InversionBaseSetIteration

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

subroutine InversionBaseStep(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_base_type) :: this

  call this%driver%PrintErrMsg('InversionBaseStep must be extended.')

end subroutine InversionBaseStep

! ************************************************************************** !

subroutine InversionBaseUpdateParameters(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_base_type) :: this

end subroutine InversionBaseUpdateParameters

! ************************************************************************** !

subroutine InversionBaseCalculateUpdate(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_base_type) :: this

end subroutine InversionBaseCalculateUpdate

! ************************************************************************** !

subroutine InversionBaseCheckConvergence(this)
  !
  ! Check inversion convergence
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
  !
  class(inversion_base_type) :: this

end subroutine InversionBaseCheckConvergence

! ************************************************************************** !

subroutine InvBaseEvaluateCostFunction(this)
  !
  ! Computes cost functions
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
  !
  class(inversion_base_type) :: this

end subroutine InvBaseEvaluateCostFunction

! ************************************************************************** !

subroutine InvBaseUpdateRegularizParams(this)
  !
  ! Computes cost functions
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/18/21
  !
  class(inversion_base_type) :: this

end subroutine InvBaseUpdateRegularizParams

! ************************************************************************** !

subroutine InversionBaseWriteIterationInfo(this)
  !
  ! Writes inversion run info
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 07/01/21
  !
  class(inversion_base_type) :: this

end subroutine InversionBaseWriteIterationInfo

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
