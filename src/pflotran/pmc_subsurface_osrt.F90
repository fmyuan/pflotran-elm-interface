module PMC_Subsurface_OSRT_class

#include "petsc/finclude/petscts.h"
  use petscts

  use PMC_Subsurface_class
  use Realization_Subsurface_class

  use PFLOTRAN_Constants_module

  implicit none

  
  private

  type, public, extends(pmc_subsurface_type) :: pmc_subsurface_osrt_type
  contains
    procedure, public :: Init => PMCSubsurfaceOSRTInit
    procedure, public :: SetupSolvers => PMCSubsurfaceOSRTSetupSolvers
    procedure, public :: StepDT => PMCSubsurfaceOSRTStepDT
    procedure, public :: Destroy => PMCSubsurfaceOSRTDestroy
  end type pmc_subsurface_osrt_type
  
  public :: PMCSubsurfaceOSRTCreate, &
            PMCSubsurfaceOSRTInit
  
contains

! ************************************************************************** !

function PMCSubsurfaceOSRTCreate()
  ! 
  ! Allocates and initializes a new process_model_coupler
  ! object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  ! 

  implicit none
  
  class(pmc_subsurface_osrt_type), pointer :: PMCSubsurfaceOSRTCreate
  
  class(pmc_subsurface_osrt_type), pointer :: pmc

  allocate(pmc)
  call pmc%Init()
  
  PMCSubsurfaceOSRTCreate => pmc  
  
end function PMCSubsurfaceOSRTCreate

! ************************************************************************** !

subroutine PMCSubsurfaceOSRTInit(this)
  ! 
  ! Initializes a new process model coupler object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  ! 
  
  implicit none
  
  class(pmc_subsurface_osrt_type) :: this
  
  call PMCSubsurfaceInit(this)
  this%name = 'PMCSubsurfaceOSRT'

end subroutine PMCSubsurfaceOSRTInit

! ************************************************************************** !

subroutine PMCSubsurfaceOSRTSetupSolvers(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  ! 
  use Option_module
  use Solver_module
  use Discretization_module
  use PM_RT_class

  implicit none

  class(pmc_subsurface_osrt_type) :: this

  type(option_type), pointer :: option
  type(solver_type), pointer :: solver
  character(len=MAXSTRINGLENGTH) :: string
  class(pm_rt_type), pointer :: pm_rt
  PetscErrorCode :: ierr

  option => this%option
  solver => this%timestepper%solver

  select type(pm=>this%pm_ptr%pm)
    class is(pm_rt_type)
      pm_rt => pm
  end select

  call SolverCreateKSP(solver,option%mycomm)
  ! set solver pointer within pm for convergence purposes
  call this%pm_ptr%pm%SetSolver(solver)

  call PrintMsg(option,"  Beginning setup of TRAN KSP")
  call KSPSetOptionsPrefix(solver%ksp, "tran_",ierr);CHKERRQ(ierr)
  call SolverCheckCommandLine(solver)

  solver%J_mat_type = MATAIJ
  solver%Jpre_mat_type = MATAIJ
  !TODO(geh): solver%J -> solver%M and XXXCreateJacobian -> XXXCreateMatrix
  call DiscretizationCreateJacobian(pm_rt%realization%discretization, &
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

end subroutine PMCSubsurfaceOSRTSetupSolvers

! ************************************************************************** !

subroutine PMCSubsurfaceOSRTStepDT(this,local_stop_flag)
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  ! 
  implicit none

  class(pmc_subsurface_osrt_type) :: this
  PetscInt :: local_stop_flag

  call this%timestepper%StepDT(this%pm_list,local_stop_flag)

end subroutine PMCSubsurfaceOSRTStepDT

! ************************************************************************** !

recursive subroutine PMCSubsurfaceOSRTFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  ! 

  use Option_module
  
  implicit none
  
  class(pmc_subsurface_osrt_type) :: this
  
  nullify(this%realization)
  
end subroutine PMCSubsurfaceOSRTFinalizeRun

! ************************************************************************** !

subroutine PMCSubsurfaceOSRTStrip(this)
  !
  ! Deallocates members of PMC Subsurface.
  !
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  
  implicit none
  
  class(pmc_subsurface_osrt_type) :: this

  call PMCSubsurfaceStrip(this)
  nullify(this%realization)

end subroutine PMCSubsurfaceOSRTStrip

! ************************************************************************** !

recursive subroutine PMCSubsurfaceOSRTDestroy(this)
  ! 
  ! ProcessModelCouplerDestroy: Deallocates a process_model_coupler object
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/06/19
  ! 

  use Option_module

  implicit none
  
  class(pmc_subsurface_osrt_type) :: this
  
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
  
  !TODO(geh): place this routine in PMC_Base_class and redirect Strip() to
  !           avoid creating all these Destroy routines
  call PMCSubsurfaceOSRTStrip(this)
  
end subroutine PMCSubsurfaceOSRTDestroy

! ************************************************************************** !
  
end module PMC_Subsurface_OSRT_class
