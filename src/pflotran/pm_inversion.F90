module PM_Inversion_class

#include "petsc/finclude/petscvec.h"
  use petscvec
  use PM_Base_class
  use Realization_Subsurface_class
  use Communicator_Base_class

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_base_type) :: pm_inversion_type
    class(realization_subsurface_type), pointer :: realization
    character(len=MAXWORDLENGTH) :: ctype
    procedure(PMInversionEvaluate), pointer :: Evaluate => null()
  contains
    procedure, public :: Setup => PMInversionSetup
    procedure, public :: InitializeRun => PMInversionInitializeRun
    procedure, public :: FinalizeRun => PMInversionFinalizeRun
    procedure, public :: Destroy => PMInversionDestroy
  end type pm_inversion_type

  ! interface blocks
  interface
    subroutine PMInversionEvaluate(this,time,ierr)
      import :: pm_inversion_type
      implicit none
      class(pm_inversion_type) :: this
      PetscReal :: time
      PetscErrorCode :: ierr
    end subroutine PMInversionEvaluate
  end interface

  public :: PMInversionCreate, &
            PMInversionInit, &
            PMInversionSetFunctionPointer

contains

! ************************************************************************** !

function PMInversionCreate()
  !
  ! Creates an inversion process model
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22
  !

  implicit none

  class(pm_inversion_type), pointer :: PMInversionCreate

  class(pm_inversion_type), pointer :: pm

  allocate(pm)
  call PMInversionInit(pm)

  PMInversionCreate => pm

end function PMInversionCreate

! ************************************************************************** !

subroutine PMInversionInit(this)
  !
  ! Initializes inversion process model
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22

  implicit none

  class(pm_inversion_type) :: this

  nullify(this%realization)
  this%ctype = ''
  this%name = ''

  call PMBaseInit(this)
  ! restart not currently supported for inversion pm's, and not needed.
  this%skip_restart = PETSC_TRUE

end subroutine PMInversionInit

! ************************************************************************** !

subroutine PMInversionSetup(this)
  !
  ! Dummy setup routine
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22

  implicit none

  class(pm_inversion_type) :: this

end subroutine PMInversionSetup

! ************************************************************************** !

subroutine PMInversionSetFunctionPointer(this,string)
  !
  ! Initializes inversion process model
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22

  use Option_module

  implicit none

  class(pm_inversion_type) :: this
  character(len=*) :: string

  this%ctype = trim(string)
  select case(string)
    case('INVERSION_MEASUREMENT')
      this%Evaluate => PMInversionInversionMeasurement
      this%header = 'INVERSION MEASUREMENT'
    case('INVERSION_ADJOINT')
      this%Evaluate => PMInversionInversionAdjoint
      this%header = 'INVERSION ADJOINT'
    case default
      this%option%io_buffer = 'Function pointer "' // trim(string) // '" not &
        &found among available functions in PMInversionSetFunctionPointer.'
      call PrintErrMsg(this%option)
  end select

end subroutine PMInversionSetFunctionPointer

! ************************************************************************** !

recursive subroutine PMInversionInitializeRun(this)
  !
  ! Initializes the time stepping
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22

  use Condition_module
  use Coupler_module
  use Option_module

  implicit none

  class(pm_inversion_type) :: this

end subroutine PMInversionInitializeRun

! ************************************************************************** !

subroutine PMInversionInversionMeasurement(this,time,ierr)
  !
  ! Takes measurements for inversion calculation
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22

  use Inversion_TS_Aux_module
  use Option_module
  use Realization_Subsurface_class
  use Utility_module

  implicit none

  class(pm_inversion_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  type(inversion_forward_aux_type), pointer :: inversion_forward_aux

  ierr = 0
  inversion_forward_aux => this%realization%patch%aux%inversion_forward_aux
  if (associated(inversion_forward_aux)) then
    if (Equal(inversion_forward_aux%sync_times( &
                inversion_forward_aux%isync_time),time)) then
      call RealizationGetObservedVariables(this%realization)
      call InversionForwardAuxMeasure(inversion_forward_aux,time,this%option)
    else
      call PrintMsg(this%option,'  No measurement requested at this time.')
    endif
  endif

end subroutine PMInversionInversionMeasurement

! ************************************************************************** !

subroutine PMInversionInversionAdjoint(this,time,ierr)
  !
  ! Initializes inversion process model
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22

  use Inversion_TS_Aux_module
  use Realization_Subsurface_class

  implicit none

  class(pm_inversion_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  ierr = 0
  if (associated(this%realization%patch%aux%inversion_forward_aux)) then
    call InversionForwardAuxStep(this%realization%patch%aux% &
                                   inversion_forward_aux,time)
  endif

end subroutine PMInversionInversionAdjoint

! ************************************************************************** !

recursive subroutine PMInversionFinalizeRun(this)
  !
  ! Finalizes the run
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22
  !

  implicit none

  class(pm_inversion_type) :: this

  ! do something here

  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif

end subroutine PMInversionFinalizeRun

! ************************************************************************** !

subroutine PMInversionDestroy(this)
  !
  ! Destroys inversion process model
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/22

  implicit none

  class(pm_inversion_type) :: this

  call PMBaseDestroy(this)

  ! destroyed in realization
  nullify(this%realization)
  nullify(this%option)

end subroutine PMInversionDestroy

end module PM_Inversion_class
