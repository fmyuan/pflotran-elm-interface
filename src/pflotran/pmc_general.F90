module PMC_General_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PMC_Base_class
  use PM_Base_class
  use PFLOTRAN_Constants_module

  implicit none


  private
  type, public, extends(pmc_base_type) :: pmc_general_type
    class(pm_base_type), pointer :: pm
    PetscBool :: evaluate_at_end_of_simulation
  contains
    procedure, public :: Init => PMCGeneralInit
    procedure, public :: RunToTime => PMCGeneralRunToTime
    procedure, public :: Destroy => PMCGeneralDestroy
  end type pmc_general_type

  public :: PMCGeneralCreate

contains

! ************************************************************************** !

function PMCGeneralCreate(name_,pm_base)
  !
  ! Allocates and initializes a new process_model_coupler
  ! object.
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  implicit none

  character(len=*) :: name_
  class(pm_base_type), pointer :: pm_base

  class(pmc_general_type), pointer :: PMCGeneralCreate

  class(pmc_general_type), pointer :: pmc

  allocate(pmc)
  call pmc%Init()

  if (len_trim(name_) > 0) then
    call pmc%SetName(name_)
  endif
  call pmc%SetOption(pm_base%option)
  pmc%pm_list => pm_base
  pmc%pm => pm_base

  PMCGeneralCreate => pmc

end function PMCGeneralCreate

! ************************************************************************** !

subroutine PMCGeneralInit(this)
  !
  ! Initializes a new process model coupler object.
  !
  ! Author: Glenn Hammond
  ! Date: 06/10/13
  !

  implicit none

  class(pmc_general_type) :: this

  call PMCBaseInit(this)
  this%name = 'PMCGeneral'
  this%evaluate_at_end_of_simulation = PETSC_TRUE

end subroutine PMCGeneralInit

! ************************************************************************** !

recursive subroutine PMCGeneralRunToTime(this,sync_time,stop_flag)
  !
  ! Runs the actual simulation.
  !
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  !
  use Option_module
  use PM_Auxiliary_class
  use PM_Inversion_class
  use PM_Parameter_class
  use Timestepper_Base_class

  implicit none

#include "petsc/finclude/petscviewer.h"

  class(pmc_general_type), target :: this
  PetscReal :: sync_time
  PetscInt :: stop_flag

  PetscInt :: local_stop_flag
  PetscErrorCode :: ierr

  if (stop_flag == TS_STOP_FAILURE) return

  if (this%stage /= 0) then
    call PetscLogStagePush(this%stage,ierr);CHKERRQ(ierr)
  endif
  this%option%io_buffer = trim(this%name)
  call PrintVerboseMsg(this%option)

  ! Get data of other process-model
  call this%GetAuxData()

  local_stop_flag = TS_CONTINUE
  if (stop_flag /= TS_STOP_END_SIMULATION .or. &
      this%evaluate_at_end_of_simulation) then
    call this%PrintHeader()
    ! must use ierr here due to 32-/64-bit integer issues
    select type(pm_=>this%pm)
      class is(pm_auxiliary_type)
        call pm_%Evaluate(sync_time,ierr)
      class is(pm_inversion_type)
        call pm_%Evaluate(sync_time,ierr)
      class is(pm_parameter_type)
        call pm_%Update(sync_time,ierr)
      class default
        this%option%io_buffer = 'PMC General not configured for PM ' // &
          trim(pm_%name) // '.'
        call PrintErrMsg(this%option)
    end select
    local_stop_flag = ierr
  endif

  ! Run underlying process model couplers
  if (associated(this%child)) then
    ! Set data needed by process-models
    call this%SetAuxData()
    call this%child%RunToTime(this%timestepper%target_time,local_stop_flag)
    ! Get data from other process-models
    call this%GetAuxData()
  endif

  ! Set data needed by process-model
  call this%SetAuxData()

  ! Run neighboring process model couplers
  if (associated(this%peer)) then
    call this%peer%RunToTime(sync_time,local_stop_flag)
  endif

  stop_flag = max(stop_flag,local_stop_flag)

  if (this%stage /= 0) then
    call PetscLogStagePop(ierr);CHKERRQ(ierr)
  endif

end subroutine PMCGeneralRunToTime

! ************************************************************************** !
!
! PMCGeneralFinalizeRun: Finalizes the time stepping
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine PMCGeneralFinalizeRun(this)
  !
  ! Finalizes the time stepping
  !
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  !

  use Option_module

  implicit none

  class(pmc_general_type) :: this

end subroutine PMCGeneralFinalizeRun

! ************************************************************************** !

subroutine PMCGeneralStrip(this)
  !
  ! Deallocates members of PMC General.
  !
  ! Author: Glenn Hammond
  ! Date: 01/13/14

  implicit none

  class(pmc_general_type) :: this

  call PMCBaseStrip(this)
  nullify(this%pm)

end subroutine PMCGeneralStrip

! ************************************************************************** !

recursive subroutine PMCGeneralDestroy(this)
  !
  ! ProcessModelCouplerDestroy: Deallocates a process_model_coupler object
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  !

  use Option_module

  implicit none

  class(pmc_general_type) :: this

  call PMCGeneralStrip(this)

  if (associated(this%child)) then
    call this%child%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%child)
    nullify(this%child)
  endif

  if (associated(this%peer)) then
    call this%peer%Destroy()
    ! destroy does not currently destroy; it strips
    deallocate(this%peer)
    nullify(this%peer)
  endif

end subroutine PMCGeneralDestroy

end module PMC_General_class
