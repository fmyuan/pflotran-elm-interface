module PM_Unit_Test_class

#include "petsc/finclude/petscvec.h"
  use petscvec
  use PM_Base_class
  use Realization_Subsurface_class
  use Communicator_Base_class
  use Option_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_base_type) :: pm_unittest_type
    class(realization_subsurface_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    character(len=MAXWORDLENGTH) :: filename
    PetscReal :: tolerance
    type(pm_unittest_type), pointer :: next_unittest
    procedure(PMUnitTestRunTest), pointer :: RunUnitTest => null()
  contains
    procedure, public :: Setup => PMUnitTestSetup
    procedure, public :: SetRealization => PMUnitTestSetRealization
    procedure, public :: InitializeRun => PMUnitTestInitializeRun
    procedure, public :: RunUnitTests => PMUnitTestRunTests
    procedure, public :: FinalizeRun => PMUnitTestFinalizeRun
    procedure, public :: Destroy => PMUnitTestDestroy
  end type pm_unittest_type

  ! interface blocks
  interface
    subroutine PMUnitTestRunTest(this,realization)
      import :: pm_unittest_type
      import :: realization_subsurface_type
      implicit none
      class(pm_unittest_type) :: this
      type(realization_subsurface_type) :: realization
    end subroutine PMUnitTestRunTest
  end interface

  public :: PMUnitTestCreate, &
            PMUnitTestInit, &
            PMUnitTestCast

contains

! ************************************************************************** !

function PMUnitTestCast(this)
  !
  ! Casts a base process model to unittest
  !
  ! Author: David Fukuyama, SNL
  ! Date: 03/06/25

  use Option_module

  class(pm_base_type), pointer :: this

  class(pm_unittest_type), pointer :: PMUnitTestCast

  nullify(PMUnitTestCast)
  if (.not.associated(this)) return
  select type (this)
    class is (pm_unittest_type)
      PMUnitTestCast => this
    class default
      this%option%io_buffer = 'Cannot cast pm_base_type to pm_unittest_type &
        &in PMUnitTestCast.'
      call PrintErrMsg(this%option)
  end select

end function PMUnitTestCast

! ************************************************************************** !

function PMUnitTestCreate()
  !
  ! Creates an unittest process model
  !
  ! Author: David Fukuyama, SNL
  ! Date: 03/06/25
  !

  class(pm_unittest_type), pointer :: PMUnitTestCreate

  class(pm_unittest_type), pointer :: pm

  allocate(pm)
  call PMUnitTestInit(pm)

  PMUnitTestCreate => pm

end function PMUnitTestCreate

! ************************************************************************** !

subroutine PMUnitTestInit(this)
  !
  ! Initializes unittest process model
  !
  ! Author: David Fukuyama, SNL
  ! Date: 03/06/25

  class(pm_unittest_type) :: this

  call PMBaseInit(this)
  nullify(this%realization)
  nullify(this%comm1)
  nullify(this%next_unittest)
  this%tolerance = 1.d-10

end subroutine PMUnitTestInit

! ************************************************************************** !

subroutine PMUnitTestSetup(this)
  !
  ! Sets up unittest process model
  !
  ! Author: David Fukuyama, SNL
  ! Date: 03/06/25

  use Option_module
  use String_module

  class(pm_unittest_type) :: this

  call this%SetRealization()

  this%header = 'UNITTEST'

end subroutine PMUnitTestSetup

! ************************************************************************** !

subroutine PMUnitTestSetRealization(this)
  !
  ! Sets the realization pointer
  !
  ! Author: David Fukuyama, SNL
  ! Date: 01/15/23

  use Realization_Subsurface_class

  class(pm_unittest_type) :: this

  this%realization => RealizationCast(this%realization_base)

end subroutine PMUnitTestSetRealization

! ************************************************************************** !

recursive subroutine PMUnitTestInitializeRun(this)
  !
  ! Initializes the process model for the simulation
  !
  ! Author: David Fukuyama, SNL
  ! Date: 03/06/25

  class(pm_unittest_type) :: this

  ! do nothing

end subroutine PMUnitTestInitializeRun

! ************************************************************************** !

subroutine PMUnitTestRunTests(this,sync_time,ierr)
  !
  ! Author: David Fukuyama, SNL
  ! Date: 03/06/25
  !
  implicit none

  class(pm_unittest_type) :: this
  class(pm_unittest_type), pointer :: next_unittest
  PetscReal :: sync_time
  PetscErrorCode :: ierr

  next_unittest => this%next_unittest
  call this%RunUnitTest(this%realization)
  do
    call next_unittest%RunUnitTest(this%realization)
    next_unittest => next_unittest%next_unittest
    if (.not.associated(next_unittest%RunUnitTest))return
  end do

end subroutine PMUnitTestRunTests

! ************************************************************************** !

recursive subroutine PMUnitTestFinalizeRun(this)
  !
  ! Finalizes the run
  !
  ! Author: David Fukuyama, SNL
  ! Date: 03/06/25
  !

  class(pm_unittest_type) :: this

  ! do something here

  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif

end subroutine PMUnitTestFinalizeRun

! ************************************************************************** !

subroutine PMUnitTestDestroy(this)
  !
  ! Destroys unittest process model
  !
  ! Author: David Fukuyama, SNL
  ! Date: 03/06/25

  class(pm_unittest_type) :: this

  call PMBaseDestroy(this)

  ! destroyed in realization
  nullify(this%realization)
  nullify(this%comm1)
  nullify(this%option)

end subroutine PMUnitTestDestroy

end module PM_Unit_Test_class
