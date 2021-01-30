module PMC_Geophysics_class

  use PMC_Base_class
  use PM_ERT_class
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
    procedure, public :: SetupSolvers => PMCGeophysicsSetupSolvers    
    procedure, public :: Destroy => PMCGeophysicsDestroy
  end type pmc_geophysics_type

  public :: PMCGeophysicsCreate, &
            PMCGeophysicsInit

contains

! ************************************************************************** !

function PMCGeophysicsCreate()
  !
  ! Allocates and initializes a new process_model_coupler
  ! object for geophysics.
  !
  ! Author: Glenn Hammond
  ! Date: 01/29/21
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
  ! Date: 01/29/21
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

subroutine PMCGeophysicsSetupSolvers(this)
  !
  ! Author: Glenn Hammond & Piyoosh Jaysaval
  ! Date: 01/29/21
  !
  use Option_module
  use Timestepper_KSP_class

  implicit none

  class(pmc_geophysics_type) :: this

  type(option_type), pointer :: option
  class(pm_ert_type), pointer :: pm_ert
  PetscErrorCode :: ierr

  option => this%option

  select type(ts=>this%timestepper)
    class is(timestepper_KSP_type)
    class default
      option%io_buffer = 'A KSP timestepper must be used for geophysics.'
      call PrintErrMsg(option)
  end select

  select type(pm=>this%pm_ptr%pm)
    class is(pm_ert_type)
      pm_ert => pm
      ! Sets up solver for ERT
      call PMERTSetupSolver(pm_ert)
    class default
      option%io_buffer = 'Solver setup implemented only for ERT &
                          &geophysics process model.'
      call PrintErrMsg(option)  
  end select

end subroutine PMCGeophysicsSetupSolvers

! ************************************************************************** !

subroutine PMCGeophysicsStepDT(this,stop_flag)
  !
  ! Solves a round of steady-state geophysics solutions.
  !
  ! Author: Glenn Hammond & Piyoosh Jaysaval
  ! Date: 01/29/21
  !
  use Option_module
  
  implicit none

  class(pmc_geophysics_type) :: this
  PetscInt :: stop_flag

  class(pm_ert_type), pointer :: pm_ert
  type(option_type), pointer :: option

  PetscLogDouble :: log_outer_start_time
  PetscLogDouble :: log_end_time
  PetscErrorCode :: ierr

  call PetscTime(log_outer_start_time,ierr);CHKERRQ(ierr)

  option => this%option

  select type(pm=>this%pm_ptr%pm)
    class is(pm_ert_type)
      pm_ert => pm
      call PMERTSolve(pm_ert) 
    class default
      option%io_buffer = 'StepDT implemented only for ERT &
                          &geophysics process model.'
      call PrintErrMsg(option)   
  end select

  call PetscTime(log_end_time,ierr);CHKERRQ(ierr)
  this%cumulative_time = this%cumulative_time + &
    log_end_time - log_outer_start_time
  ! call this%timer%Start()
  ! do ielectrode = 1, this%num_electrodes
  ! enddo
  ! call this%timer%Stop()

end subroutine PMCGeophysicsStepDT

! ************************************************************************** !

recursive subroutine PMCGeophysicsFinalizeRun(this)
  !
  ! Finalizes the time stepping
  !
  ! Author: Glenn Hammond
  ! Date: 01/29/21
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
  ! Date: 01/29/21

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
  ! Date: 01/29/21
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
