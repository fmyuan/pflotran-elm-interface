module PM_Parameter_class

#include "petsc/finclude/petscvec.h"
  use petscvec
  use PM_Base_class
  use Realization_Subsurface_class
  use Communicator_Base_class

  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscInt, parameter :: PARAMETER_BY_CELL = 1
  PetscInt, parameter :: PARAMETER_BY_MATERIAL = 2

  type, public, extends(pm_base_type) :: pm_parameter_type
    class(realization_subsurface_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    type(parameter_type), pointer :: parameter
  contains
    procedure, public :: Setup => PMParameterSetup
    procedure, public :: InitializeRun => PMParameterInitializeRun
    procedure, public :: FinalizeRun => PMParameterFinalizeRun
    procedure, public :: Destroy => PMParameterDestroy
  end type pm_parameter_type

  type, public :: parameter_type
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXWORDLENGTH) :: dataset_name
    PetscInt :: index_in_global_parameter

  end type parameter_type

  public :: PMParameterCreate, &
            PMParameterInit, &
            PMParameterCast, &
            PMParameterRead

contains

! ************************************************************************** !

function PMParameterCreate()
  !
  ! Creates an parameter process model
  !
  ! Author: Glenn Hammond
  ! Date: 12/22/23
  !

  implicit none

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

  implicit none

  class(pm_parameter_type) :: this

  call PMBaseInit(this)
  nullify(this%realization)
  nullify(this%comm1)
  nullify(this%parameter)

end subroutine PMParameterInit

! ************************************************************************** !

subroutine PMParameterSetup(this)
  !
  ! Sets up parameter process model
  !
  ! Author: Glenn Hammond
  ! Date: 12/22/23

  implicit none

  class(pm_parameter_type) :: this

end subroutine PMParameterSetup

! ************************************************************************** !

function PMParameterCast(this)
  !
  ! Casts a base process model to parameter
  !
  ! Author: Glenn Hammond
  ! Date: 12/22/23

  use Option_module

  implicit none

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

subroutine PMParameterRead(this,input)
  !
  ! Author: Glenn Hammond
  ! Date: 12/22/23
  !
  use Input_Aux_module
  use Option_module
  use String_module

  implicit none

  class(pm_parameter_type), pointer :: this
  type(input_type), pointer :: input

  type(option_type), pointer :: option
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
      case('')
      case default
        call InputKeywordUnrecognized(input,keyword,error_str,option)
    end select
  enddo

end subroutine PMParameterRead

! ************************************************************************** !

recursive subroutine PMParameterInitializeRun(this)
  !
  ! Initializes the time stepping
  !
  ! Author: Glenn Hammond
  ! Date: 12/22/23

  use Condition_module
  use Coupler_module
  use Option_module
  use Reaction_Aux_module

  implicit none

  class(pm_parameter_type) :: this

end subroutine PMParameterInitializeRun

! ************************************************************************** !

recursive subroutine PMParameterFinalizeRun(this)
  !
  ! Finalizes the run
  !
  ! Author: Glenn Hammond
  ! Date: 12/22/23
  !

  implicit none

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

  implicit none

  class(pm_parameter_type) :: this

  call PMBaseDestroy(this)

  ! destroyed in realization
  nullify(this%realization)
  nullify(this%comm1)
  nullify(this%option)

end subroutine PMParameterDestroy

end module PM_Parameter_class
