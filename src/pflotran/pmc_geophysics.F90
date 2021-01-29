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
  ! Author: Glenn Hammond
  ! Date: 01/29/21
  !
  use Option_module
  use Solver_module
  use Discretization_module
  use PM_ERT_class
  use Timestepper_KSP_class

  implicit none

  class(pmc_geophysics_type) :: this

  type(option_type), pointer :: option
  type(solver_type), pointer :: solver
  character(len=MAXSTRINGLENGTH) :: string
  class(pm_ert_type), pointer :: pm_ert
  PetscErrorCode :: ierr

  option => this%option
  solver => this%timestepper%solver

  select type(ts=>this%timestepper)
    class is(timestepper_KSP_type)
    class default
      option%io_buffer = 'A KSP timestepper must be used for geophysics.'
      call PrintErrMsg(option)
  end select

  select type(pm=>this%pm_ptr%pm)
    class is(pm_ert_type)
      pm_ert => pm
  end select

  call SolverCreateKSP(solver,option%mycomm)

  call PrintMsg(option,"  Beginning setup of TRAN KSP")
  call KSPSetOptionsPrefix(solver%ksp, "tran_",ierr);CHKERRQ(ierr)
  call SolverCheckCommandLine(solver)

  solver%J_mat_type = MATAIJ
  solver%Jpre_mat_type = MATAIJ
  !TODO(geh): solver%J -> solver%M and XXXCreateJacobian -> XXXCreateMatrix
  call DiscretizationCreateJacobian(pm_ert%realization%discretization, &
                                    ONEDOF, &
                                    solver%Jpre_mat_type, &
                                    solver%Jpre,option)
  call MatSetOptionsPrefix(solver%Jpre,"tran_",ierr);CHKERRQ(ierr)
  solver%J = solver%Jpre

  ! Have PETSc do a KSP_View() at the end of each solve if
  ! verbosity > 0.
  if (option%verbosity >= 2) then
    string = '-tran_ksp_view'
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, &
                                  string, ierr);CHKERRQ(ierr)
  endif

  call SolverSetKSPOptions(solver,option)

end subroutine PMCGeophysicsSetupSolvers

! ************************************************************************** !

subroutine PMCGeophysicsStepDT(this,stop_flag)
  !
  ! Solves a round of steady-state geophysics solutions for each electrode.
  !
  ! Author: Glenn Hammond
  ! Date: 01/29/21
  !
  use Option_module
  
  implicit none

  class(pmc_geophysics_type) :: this
  PetscInt :: stop_flag

  class(pm_ert_type), pointer :: pm_ert
  type(option_type), pointer :: option

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: tunit
  
  KSPConvergedReason :: ksp_reason
  PetscLogDouble :: log_outer_start_time
  PetscLogDouble :: log_start_time
  PetscLogDouble :: log_ksp_start_time
  PetscLogDouble :: log_end_time
  PetscErrorCode :: ierr

  call PetscTime(log_outer_start_time,ierr);CHKERRQ(ierr)

  option => this%option

  select type(pm=>this%pm_ptr%pm)
    class is(pm_ert_type)
      pm_ert => pm
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
