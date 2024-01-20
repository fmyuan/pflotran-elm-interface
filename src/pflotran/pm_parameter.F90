module PM_Parameter_class

#include "petsc/finclude/petscvec.h"
  use petscvec
  use PM_Base_class
  use Realization_Subsurface_class
  use Communicator_Base_class
  use Parameter_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscInt, parameter, public :: UPDATE_AFTER_LAST_PM = 0
  PetscInt, parameter, public :: UPDATE_AFTER_FLOW = 1

  PetscInt, parameter :: PARAMETER_BY_CELL = 1
  PetscInt, parameter :: PARAMETER_BY_MATERIAL = 2

  type, public, extends(pm_base_type) :: pm_parameter_type
    class(realization_subsurface_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    type(parameter_type), pointer :: parameter
    PetscInt :: when_to_update
    procedure(PMParameterUpdate), pointer :: Update => null()
  contains
    procedure, public :: Setup => PMParameterSetup
    procedure, public :: SetRealization => PMParameterSetRealization
    procedure, public :: InitializeRun => PMParameterInitializeRun
    procedure, public :: FinalizeRun => PMParameterFinalizeRun
    procedure, public :: Destroy => PMParameterDestroy
  end type pm_parameter_type

  ! interface blocks
  interface
    subroutine PMParameterUpdate(this,time,ierr)
      import :: pm_parameter_type
      implicit none
      class(pm_parameter_type) :: this
      PetscReal :: time
      PetscErrorCode :: ierr
    end subroutine PMParameterUpdate
  end interface

  public :: PMParameterCreate, &
            PMParameterInit, &
            PMParameterCast, &
            PMParameterRead, &
            PMParameterSortPMs

contains

! ************************************************************************** !

function PMParameterCreate()
  !
  ! Creates an parameter process model
  !
  ! Author: Glenn Hammond
  ! Date: 12/22/23
  !

  class(pm_parameter_type), pointer :: PMParameterCreate

  class(pm_parameter_type), pointer :: pm

  allocate(pm)
  call PMParameterInit(pm)

  PMParameterCreate => pm

end function PMParameterCreate

! ************************************************************************** !

subroutine PMParameterInit(this)
  !
  ! Initializes parameter process model
  !
  ! Author: Glenn Hammond
  ! Date: 12/22/23

  class(pm_parameter_type) :: this

  call PMBaseInit(this)
  this%when_to_update = UNINITIALIZED_INTEGER
  nullify(this%parameter)
  nullify(this%realization)
  nullify(this%comm1)
  this%Update => PMParameterUpdateDoNothing

end subroutine PMParameterInit

! ************************************************************************** !

subroutine PMParameterSetup(this)
  !
  ! Sets up parameter process model
  !
  ! Author: Glenn Hammond
  ! Date: 12/22/23

  use Option_module
  use String_module

  class(pm_parameter_type) :: this

  type(parameter_type), pointer :: cur_parameter
  character(len=MAXWORDLENGTH) :: parameter_name

  parameter_name = this%parameter%name
  ! destroy the placeholder providing name
  call ParameterDestroy(this%parameter)
  if (len_trim(parameter_name)== 0) then
    this%option%io_buffer = 'No PARAMETER_NAME specified for Parameter &
      &process model "' // trim(this%name) // '".'
    call PrintErrMsg(this%option)
  endif

  cur_parameter => this%realization%parameter_list
  do
    if (.not.associated(cur_parameter)) exit
    if (StringCompare(cur_parameter%name,parameter_name)) then
      this%parameter => cur_parameter
      exit
    endif
    cur_parameter => cur_parameter%next
  enddo

  if (.not.associated(this%parameter)) then
    this%option%io_buffer = 'Parameter "' // trim(parameter_name) // &
      '" not found among available parameters.'
    call PrintErrMsg(this%option)
  endif

  this%header = 'PARAMETER (' // trim(this%name) // '->' // &
                trim(parameter_name) // ')'

end subroutine PMParameterSetup

! ************************************************************************** !

function PMParameterCast(this)
  !
  ! Casts a base process model to parameter
  !
  ! Author: Glenn Hammond
  ! Date: 12/22/23

  use Option_module

  class(pm_base_type), pointer :: this

  class(pm_parameter_type), pointer :: PMParameterCast

  nullify(PMParameterCast)
  if (.not.associated(this)) return
  select type (this)
    class is (pm_parameter_type)
      PMParameterCast => this
    class default
      this%option%io_buffer = 'Cannot cast pm_base_type to pm_parameter_type &
        &in PMParameterCast.'
      call PrintErrMsg(this%option)
  end select

end function PMParameterCast

! ************************************************************************** !

subroutine PMParameterRead(input,option,this)
  !
  ! Author: Glenn Hammond
  ! Date: 12/22/23
  !
  use Input_Aux_module
  use Option_module
  use String_module

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_parameter_type), pointer :: this

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_str

  error_str = 'SIMULATION,PROCESS_MODELS,PARAMETER'

  input%ierr = 0
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_str)
    call StringToUpper(keyword)

    select case(trim(keyword))
      case('PARAMETER_NAME')
        ! at this point, the parameter object solely stores the parameter name
        this%parameter => ParameterCreate()
        call InputReadWord(input,option,this%parameter%name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_str)
      case('WHEN_TO_UPDATE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_str)
        call StringToUpper(word)
        select case(trim(word))
          case('AFTER_LAST_PM')
            this%when_to_update = UPDATE_AFTER_LAST_PM
          case('AFTER_FLOW')
            this%when_to_update = UPDATE_AFTER_FLOW
          case default
            call InputKeywordUnrecognized(input,word,trim(error_str)//word, &
                                          option)
        end select
      case default
        call InputKeywordUnrecognized(input,keyword,error_str,option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine PMParameterRead

! ************************************************************************** !

subroutine PMParameterSetRealization(this,realization)
  !
  ! Sets the realization pointer
  !
  ! Author: Glenn Hammond
  ! Date: 01/15/23

  use Realization_Subsurface_class

  class(pm_parameter_type) :: this
  class(realization_subsurface_type), pointer :: realization

  this%realization => realization
  this%realization_base => realization

end subroutine PMParameterSetRealization

! ************************************************************************** !

recursive subroutine PMParameterInitializeRun(this)
  !
  ! Initializes the process model for the simulation
  !
  ! Author: Glenn Hammond
  ! Date: 12/22/23

  class(pm_parameter_type) :: this

  PetscErrorCode :: ierr

  call this%Update(0.d0,ierr)

end subroutine PMParameterInitializeRun

! ************************************************************************** !

subroutine PMParameterSortPMs(pm_list)
  !
  ! Reorder the list of PMs and ensure that they are aligned with desired
  ! time of execution
  !
  ! Author: Glenn Hammond
  ! Date: 12/22/23

  class(pm_parameter_type), pointer :: pm_list

  class(pm_parameter_type), pointer :: cur_pm, prev_pm, next_pm
  PetscBool :: swapped

  swapped = PETSC_TRUE
  do
    if (.not.swapped) exit
    swapped = PETSC_FALSE
    nullify(prev_pm)
    cur_pm => pm_list
    do
      if (.not.associated(cur_pm)) exit
      next_pm => PMParameterCast(cur_pm%next)
      select case(cur_pm%when_to_update)
        case(UPDATE_AFTER_LAST_PM)
          if (associated(next_pm)) then
            if (next_pm%when_to_update /= UPDATE_AFTER_LAST_PM) then
              call PMParameterSwapPM(pm_list,prev_pm)
              swapped = PETSC_TRUE
            endif
          endif
      end select
      prev_pm => cur_pm
      cur_pm => PMParameterCast(cur_pm%next)
    enddo
  enddo

end subroutine PMParameterSortPMs

! ************************************************************************** !

subroutine PMParameterSwapPM(pm_list,pm0)
  !
  ! Swaps two process models (pm1 and pm2) in a list
  !
  ! Author: Glenn Hammond
  ! Date: 01/18/24

  class(pm_parameter_type), pointer :: pm_list, pm0

  class(pm_parameter_type), pointer :: pm1, pm2, pm3

  ! list
  !     \
  !      -> pm1 -> pm2 -> pm3
  !     /
  !  pm0

  if (associated(pm0)) then ! beginning of list (pm0 does not exist)
    pm1 => PMParameterCast(pm0%next)
    pm2 => PMParameterCast(pm0%next%next)
    pm3 => PMParameterCast(pm0%next%next%next)
    pm0%next => pm2
  else
    pm1 => pm_list
    pm2 => PMParameterCast(pm_list%next)
    pm3 => PMParameterCast(pm_list%next%next)
    pm_list => pm2
  endif
  pm2%next => pm1
  pm1%next => pm3

end subroutine PMParameterSwapPM

! ************************************************************************** !

subroutine PMParameterUpdateDoNothing(this,time,ierr)
  !
  ! Updates the parameter
  !
  ! Author: Glenn Hammond
  ! Date: 01/18/24

  use Option_module

  class(pm_parameter_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  this%option%io_buffer = 'Doing nothing for PARAMETER "' // &
    trim(this%name) // '".'
  call PrintMsg(this%option)
  ierr = 0

end subroutine PMParameterUpdateDoNothing

! ************************************************************************** !

subroutine PMParameterUpdateReadDataset(this,time,ierr)
  !
  ! Updates the parameter
  !
  ! Author: Glenn Hammond
  ! Date: 01/18/24

  class(pm_parameter_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  ierr = 0

end subroutine PMParameterUpdateReadDataset

! ************************************************************************** !

recursive subroutine PMParameterFinalizeRun(this)
  !
  ! Finalizes the run
  !
  ! Author: Glenn Hammond
  ! Date: 12/22/23
  !

  class(pm_parameter_type) :: this

  ! do something here

  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif

end subroutine PMParameterFinalizeRun

! ************************************************************************** !

subroutine PMParameterDestroy(this)
  !
  ! Destroys parameter process model
  !
  ! Author: Glenn Hammond
  ! Date: 12/22/23

  class(pm_parameter_type) :: this

  call PMBaseDestroy(this)

  ! destroyed in realization
  nullify(this%realization)
  nullify(this%comm1)
  nullify(this%option)

end subroutine PMParameterDestroy

end module PM_Parameter_class
