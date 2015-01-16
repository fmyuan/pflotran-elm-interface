module PMC_Third_Party_class

  use PMC_Base_class
  use Realization_class
  use Option_module

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
  
  type, public, extends(pmc_base_type) :: pmc_third_party_type
    class(realization_type), pointer :: realization
  contains
    procedure, public :: Init => PMCThirdPartyInit
!    procedure, public :: InitializeRun => PMCThirdPartyInitializeRun
    procedure, public :: RunToTime => PMCThirdPartyRunToTime
    procedure, public :: FinalizeRun => PMCThirdPartyFinalizeRun
    procedure, public :: Destroy => PMCThirdPartyDestroy
    procedure, public :: GetAuxData => PMCThirdPartyGetAuxData
  end type pmc_third_party_type
  
  public :: PMCThirdPartyCreate
  
contains

! ************************************************************************** !

function PMCThirdPartyCreate()
  ! 
  ! Allocates and initializes a new
  ! process_model_coupler object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 

  implicit none
  
  class(pmc_third_party_type), pointer :: PMCThirdPartyCreate
  
  class(pmc_third_party_type), pointer :: pmc

  allocate(pmc)
  call pmc%Init()
  
  PMCThirdPartyCreate => pmc  
  
end function PMCThirdPartyCreate

! ************************************************************************** !

subroutine PMCThirdPartyInit(this)
  ! 
  ! Initializes a new process model coupler object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 

  implicit none
  
  class(pmc_third_party_type) :: this
  
  call PMCBaseInit(this)
  this%name = 'PMCThirdParty'
  nullify(this%realization) 
  
end subroutine PMCThirdPartyInit

! ************************************************************************** !

recursive subroutine PMCThirdPartyRunToTime(this,sync_time,stop_flag)
  ! 
  ! Runs the actual simulation.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 

  use Timestepper_Base_class, only : TS_CONTINUE

  implicit none
  
  class(pmc_third_party_type), target :: this
  PetscReal :: sync_time
  PetscInt :: stop_flag
  
  class(pmc_base_type), pointer :: pmc_base
  PetscInt :: local_stop_flag
  
  this%option%io_buffer = trim(this%name)
  call printVerboseMsg(this%option)
  
  call this%GetAuxData()
  
  local_stop_flag = TS_CONTINUE

  ! add time stepping here.

  ! Run neighboring process model couplers
  if (associated(this%child)) then
    call this%child%RunToTime(sync_time,local_stop_flag)
  endif

  ! Run neighboring process model couplers
  if (associated(this%peer)) then
    call this%peer%RunToTime(sync_time,local_stop_flag)
  endif

  stop_flag = max(stop_flag,local_stop_flag)  
  
end subroutine PMCThirdPartyRunToTime

! ************************************************************************** !

subroutine PMCThirdPartyGetAuxData(this)
  ! 
  ! Runs the actual simulation.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 

  implicit none
  
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscviewer.h"

  class(pmc_third_party_type) :: this

  PetscErrorCode :: ierr
  PetscViewer :: viewer

end subroutine PMCThirdPartyGetAuxData

! ************************************************************************** !

recursive subroutine PMCThirdPartyFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 

  implicit none
  
  class(pmc_third_party_type) :: this
  
  call printMsg(this%option,'PMCThirdParty%FinalizeRun()')
  
end subroutine PMCThirdPartyFinalizeRun

! ************************************************************************** !

subroutine PMCThirdPartyStrip(this)
  !
  ! Deallocates members of PMC Subsurface.
  !
  ! Author: Glenn Hammond
  ! Date: 01/13/14
  
  implicit none
  
  class(pmc_third_party_type) :: this
  
  PetscErrorCode :: ierr

  call PMCBaseStrip(this)
  nullify(this%realization)
  
end subroutine PMCThirdPartyStrip

! ************************************************************************** !

recursive subroutine PMCThirdPartyDestroy(this)
  ! 
  ! Deallocates a process_model_coupler object
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/02/13
  ! 

  use Utility_module, only: DeallocateArray 

  implicit none
  
  class(pmc_third_party_type) :: this

  PetscErrorCode :: ierr
  
  call printMsg(this%option,'PMCThirdParty%Destroy()')
  
  call PMCThirdPartyStrip(this)
  
  if (associated(this%child)) then
    call this%child%Destroy()
  endif 
  
  if (associated(this%peer)) then
    call this%peer%Destroy()
  endif  

end subroutine PMCThirdPartyDestroy
  
end module PMC_Third_Party_class