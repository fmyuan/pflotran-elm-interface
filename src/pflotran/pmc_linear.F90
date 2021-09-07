module PMC_Linear_class

  use PMC_Subsurface_class
  use Realization_Subsurface_class

  use PFLOTRAN_Constants_module

#include "petsc/finclude/petscts.h"
  use petscts

  implicit none


  private

  type, public, extends(pmc_subsurface_type) :: pmc_linear_type
  contains
    procedure, public :: Init => PMCLinearInit
    procedure, public :: SetupSolvers => PMCLinearSetupSolvers
!    procedure, public :: FinalizeRun => PMCLinearFinalizeRun
    procedure, public :: Destroy => PMCLinearDestroy
  end type pmc_linear_type

  public :: PMCLinearCreate, &
            PMCLinearInit, &
            PMCLinearStrip

contains

! ************************************************************************** !

function PMCLinearCreate()
  !
  ! Allocates and initializes a new process_model_coupler object.
  !
  ! Author: Glenn Hammond
  ! Date: 08/30/21
  !

  implicit none

  class(pmc_linear_type), pointer :: PMCLinearCreate

  class(pmc_linear_type), pointer :: pmc

  allocate(pmc)
  call pmc%Init()

  PMCLinearCreate => pmc

end function PMCLinearCreate

! ************************************************************************** !

subroutine PMCLinearInit(this)
  !
  ! Initializes a new process model coupler object.
  !
  ! Author: Glenn Hammond
  ! Date: 08/30/21
  !
  implicit none

  class(pmc_linear_type) :: this

  call PMCSubsurfaceInit(this)
  this%name = 'PMCLinear'
  nullify(this%realization)

end subroutine PMCLinearInit

! ************************************************************************** !

subroutine PMCLinearSetupSolvers(this)
  !
  ! Author: Glenn Hammond
  ! Date: 08/30/21
  !
  use Discretization_module
  use Option_module
  use Solver_module
  use PM_Subsurface_Flow_class

  implicit none

  class(pmc_linear_type) :: this

  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  class(realization_subsurface_type), pointer :: realization
  character(len=MAXSTRINGLENGTH) :: string
  PetscErrorCode :: ierr

  option => this%option
  solver => this%timestepper%solver

  call SolverCreateKSP(solver,option%mycomm)

  select type(pm => this%pm_ptr%pm)
  ! ----- subsurface flow
    class is(pm_subsurface_flow_type)
      call PrintMsg(option,"  Beginning setup of FLOW KSP ")
      if (OptionPrintToScreen(option)) then
        write(*,'(" number of dofs = ",i3)') option%nflowdof
        select case(option%iflowmode)
          case(PNF_MODE)
            write(*,'(" mode = PNF: h")')
          case default
            option%io_buffer = 'Flow mode ' // trim(option%flowmode) // &
              ' not supported by PMC Linear.'
            call PrintErrMsg(option)
        end select
      endif

      call KSPSetOptionsPrefix(solver%ksp,"flow_",ierr);CHKERRQ(ierr)
      call SolverCheckCommandLine(solver)

      if (Uninitialized(solver%Mpre_mat_type) .and. &
          Uninitialized(solver%M_mat_type)) then
        solver%Mpre_mat_type = MATAIJ
        solver%M_mat_type = solver%Mpre_mat_type
      else if (Uninitialized(solver%M_mat_type)) then
        solver%M_mat_type = solver%Mpre_mat_type
      endif

      call DiscretizationCreateMatrix(pm%realization%discretization, &
                                      NFLOWDOF, &
                                      solver%Mpre_mat_type, &
                                      solver%Mpre, &
                                      option)

      call MatSetOptionsPrefix(solver%Mpre,"flow_",ierr);CHKERRQ(ierr)

      if (solver%Mpre_mat_type == solver%M_mat_type) then
        solver%M = solver%Mpre
      else
        call DiscretizationCreateMatrix(pm%realization%discretization, &
                                        NFLOWDOF, &
                                        solver%M_mat_type, &
                                        solver%M, &
                                        option)

        call MatSetOptionsPrefix(solver%M,"flow_",ierr);CHKERRQ(ierr)
      endif

      if (option%verbosity >= 2) then
        string = '-flow_ksp_view'
        call PetscOptionsInsertString(PETSC_NULL_OPTIONS,string, &
                                      ierr);CHKERRQ(ierr)
      endif

      call PrintMsg(option,"  Finished setting up FLOW SNES ")
    class default
      option%io_buffer = 'Process model not supported by PMCLinearSetupSolvers'
      call PrintErrMsg(option)
  end select
  call SolverSetKSPOptions(solver,option)

end subroutine PMCLinearSetupSolvers

! ************************************************************************** !

recursive subroutine PMCLinearFinalizeRun(this)
  !
  ! Finalizes the time stepping
  !
  ! Author: Glenn Hammond
  ! Date: 08/30/21
  !

  ! remove if PMC_Base_class is included at op
  use PMC_Base_class

  implicit none

  class(pmc_linear_type) :: this

  nullify(this%realization)
  call PMCBaseFinalizeRun(this)

end subroutine PMCLinearFinalizeRun

! ************************************************************************** !

subroutine PMCLinearStrip(this)
  !
  ! Deallocates members of PMC Subsurface.
  !
  ! Author: Glenn Hammond
  ! Date: 08/30/21

  implicit none

  class(pmc_linear_type) :: this

  call PMCSubsurfaceStrip(this)
  nullify(this%realization)

end subroutine PMCLinearStrip

! ************************************************************************** !

recursive subroutine PMCLinearDestroy(this)
  !
  ! ProcessModelCouplerDestroy: Deallocates a process_model_coupler object
  !
  ! Author: Glenn Hammond
  ! Date: 08/30/21
  !
  implicit none

  class(pmc_linear_type) :: this

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

  call PMCLinearStrip(this)

end subroutine PMCLinearDestroy

! ************************************************************************** !

end module PMC_Linear_class
