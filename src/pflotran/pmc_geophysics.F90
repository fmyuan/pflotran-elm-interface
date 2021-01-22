module PMC_Geophysics_class

  use PMC_Base_class
  use Realization_Subsurface_class

  use PFLOTRAN_Constants_module

#include "petsc/finclude/petscts.h"
  use petscts

  implicit none

  
  private

  type, public, extends(pmc_base_type) :: pmc_geophysics_type
    class(realization_subsurface_type), pointer :: realization
  contains
    procedure, public :: Init => PMCGeophysicsInit
    procedure, public :: Destroy => PMCGeophysicsDestroy
  end type pmc_geophysics_type
  
  public :: PMCGeophysicsCreate, &
            PMCGeophysicsInit, &
            PMCGeophysicsStrip
  
contains

! ************************************************************************** !

function PMCGeophysicsCreate()
  ! 
  ! Allocates and initializes a new process_model_coupler
  ! object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pmc_geophysics_type), pointer :: PMCGeophysicsCreate
  
  class(pmc_geophysics_type), pointer :: pmc

  allocate(pmc)
  call pmc%Init()
  
  PMCGeophysicsCreate => pmc  
  
end function PMCGeophysicsCreate

! ************************************************************************** !

subroutine PMCGeophysicsInit(this)
  ! 
  ! Initializes a new process model coupler object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/10/13
  ! 
  ! for some reason, Intel with VS want this explicitly specified.
  use PMC_Base_class, only : PMCBaseInit
  
  implicit none
  
  class(pmc_geophysics_type) :: this
  
  call PMCBaseInit(this)
  this%name = 'PMCGeophysics'
  nullify(this%realization)

end subroutine PMCGeophysicsInit

! ************************************************************************** !

recursive subroutine PMCGeophysicsFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  use Option_module
  
  implicit none
  
  class(pmc_geophysics_type) :: this
  
#ifdef DEBUG
  call PrintMsg(this%option,'PMCGeophysics%FinalizeRun()')
#endif
  
  nullify(this%realization)
  
end subroutine PMCGeophysicsFinalizeRun

! ************************************************************************** !

subroutine PMCGeophysicsStrip(this)
  !
  ! Deallocates members of PMC Geophysics.
  !
  ! Author: Glenn Hammond
  ! Date: 01/13/14
  
  implicit none
  
  class(pmc_geophysics_type) :: this

  call PMCBaseStrip(this)
  nullify(this%realization)

end subroutine PMCGeophysicsStrip

! ************************************************************************** !

recursive subroutine PMCGeophysicsDestroy(this)
  ! 
  ! ProcessModelCouplerDestroy: Deallocates a process_model_coupler object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Option_module

  implicit none
  
  class(pmc_geophysics_type) :: this
  
#ifdef DEBUG
  call PrintMsg(this%option,'PMCGeophysics%Destroy()')
#endif

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
  
  call PMCGeophysicsStrip(this)
  
end subroutine PMCGeophysicsDestroy

! ************************************************************************** !
  
end module PMC_Geophysics_class
