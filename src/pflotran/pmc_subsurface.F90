module PMC_Subsurface_class

  use PMC_Base_class
  use Realization_Subsurface_class

  use PFLOTRAN_Constants_module

#include "petsc/finclude/petscts.h"
  use petscts

  implicit none

  
  private

  type, public, extends(pmc_base_type) :: pmc_subsurface_type
    class(realization_subsurface_type), pointer :: realization
  contains
    procedure, public :: Init => PMCSubsurfaceInit
    procedure, public :: SetupSolvers => PMCSubsurfaceSetupSolvers
    procedure, public :: GetAuxData => PMCSubsurfaceGetAuxData
    procedure, public :: SetAuxData => PMCSubsurfaceSetAuxData
    procedure, public :: Destroy => PMCSubsurfaceDestroy
  end type pmc_subsurface_type
  
  public :: PMCSubsurfaceCreate
  
contains

! ************************************************************************** !

function PMCSubsurfaceCreate()
  ! 
  ! Allocates and initializes a new process_model_coupler
  ! object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pmc_subsurface_type), pointer :: PMCSubsurfaceCreate
  
  class(pmc_subsurface_type), pointer :: pmc

#ifdef DEBUG
  print *, 'PMCSubsurface%Create()'
#endif
  
  allocate(pmc)
  call pmc%Init()
  
  PMCSubsurfaceCreate => pmc  
  
end function PMCSubsurfaceCreate

! ************************************************************************** !

subroutine PMCSubsurfaceInit(this)
  ! 
  ! Initializes a new process model coupler object.
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/10/13
  ! 
  ! for some reason, Intel with VS want this explicitly specified.
  use PMC_Base_class, only : PMCBaseInit
  
  implicit none
  
  class(pmc_subsurface_type) :: this
  
#ifdef DEBUG
  print *, 'PMCSubsurface%Init()'
#endif
  
  call PMCBaseInit(this)
  this%name = 'PMCSubsurface'
  nullify(this%realization)

end subroutine PMCSubsurfaceInit

! ************************************************************************** !

subroutine PMCSubsurfaceSetupSolvers(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 06/22/18
  ! 
  use Option_module
  use Timestepper_BE_class


  implicit none

  class(pmc_subsurface_type) :: this

  type(option_type), pointer :: option

  option => this%option

  if (associated(this%timestepper)) then
    select type(ts => this%timestepper)
      class is(timestepper_BE_type)
        call PMCSubsurfaceSetupSolvers_TimestepperBE(this)
      class default
        option%io_buffer = &
          'Unknown timestepper found in PMCSubsurfaceSetupSolvers '
        call PrintErrMsg(option)
    end select
  endif

end subroutine PMCSubsurfaceSetupSolvers

! ************************************************************************** !

subroutine PMCSubsurfaceSetupSolvers_TimestepperBE(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 
  use Convergence_module
  use Discretization_module
  use Realization_Subsurface_class
  use Option_module
  use Flowmode_Aux_module
  use PMC_Base_class
  use PM_Base_Pointer_module
  use PM_Base_class
  use PM_Subsurface_Flow_class
  use PM_TH_class
  use Solver_module
  use Timestepper_Base_class
  use Timestepper_BE_class

  implicit none

  class(pmc_subsurface_type) :: this

  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  class(realization_subsurface_type), pointer :: realization
  SNESLineSearch :: linesearch  
  character(len=MAXSTRINGLENGTH) :: string  
  PetscBool :: add_pre_check, check_update, check_post_convergence
  PetscInt :: itransport
  PetscInt :: trans_coupling
  SNESType :: snes_type
  PetscErrorCode :: ierr

#ifdef DEBUG
  call PrintMsg(this%option,'PMCSubsurface%SetupSolvers()')
#endif

  option => this%option
  solver => this%timestepper%solver
  
  check_update = PETSC_FALSE
  itransport = 0
  trans_coupling = 10000 ! set to high value (not 0 or 1)

  call SolverCreateSNES(solver,option%mycomm)
  call SNESGetLineSearch(this%timestepper%solver%snes,linesearch, &
                         ierr);CHKERRQ(ierr)
  ! set solver pointer within pm for convergence purposes
  call this%pm_ptr%pm%SetSolver(solver)
  select type(pm => this%pm_ptr%pm)
  ! ----- subsurface flow
    class is(pm_subsurface_flow_type)
      call PrintMsg(option,"  Beginning setup of FLOW SNES ")
      if (solver%J_mat_type == MATAIJ .and. &
          option%iflowmode /= TH_MODE) then

        option%io_buffer = 'AIJ matrix not supported for current &
          &mode: '// option%flowmode
        call PrintErrMsg(option)
      endif
      if (OptionPrintToScreen(option)) then
        write(*,'(" number of dofs = ",i3,", number of &
                  &fluids = ",i3,i2)') option%nflowdof, option%nfluids
        select case(option%iflowmode)
          case(TH_MODE)
            write(*,'(" mode = TH: p, sg/X, T")')

        end select
      endif

      select case(option%iflowmode)
        case(TH_MODE)
          call SNESGetType(solver%snes,snes_type,ierr);CHKERRQ(ierr)
          if (trim(snes_type) == 'newtontr') then
            flow_using_newtontr = PETSC_TRUE
          endif
      end select

      call SNESSetOptionsPrefix(solver%snes, "flow_",ierr);CHKERRQ(ierr)
      call SolverCheckCommandLine(solver)

! ----- Set up the J and Jpre matrices -----
! 1) If neither J_mat_type or Jpre_mat_type are specified, set to default.
! 2) If only one of J_mat_type and Jpre_mat_type are specified, then default 
!    to setting the other to the same value (except for MATMFFD case).
! 3) Once J_mat_type and Jpre_mat_type are set appropriately, then 
!    * If J_mat_type == Jpre_mat_type, then set solver%J = solver%Jpre
!    * Otherwise 
!      - Create different matrices for each.
!      - Inside Jacobian routines, will need to check for 
!        solver%J != solver%Jpre, and populate two matrices if so.

      if (Uninitialized(solver%Jpre_mat_type) .and. &
          Uninitialized(solver%J_mat_type)) then
        ! Matrix types not specified, so set to default.
        solver%Jpre_mat_type = MATBAIJ
        solver%J_mat_type = solver%Jpre_mat_type
      else if (Uninitialized(solver%Jpre_mat_type)) then
        if (solver%J_mat_type == MATMFFD) then
          solver%Jpre_mat_type = MATBAIJ
        else
          solver%Jpre_mat_type = solver%J_mat_type
        endif
      else if (Uninitialized(solver%J_mat_type)) then
        solver%J_mat_type = solver%Jpre_mat_type
      endif

      if (associated(solver%cprstash)) then
        call CPRWorkersCreate(pm, solver, option)
      endif

      call DiscretizationCreateJacobian(pm%realization%discretization, &
                                        NFLOWDOF, &
                                        solver%Jpre_mat_type, &
                                        solver%Jpre, &
                                        option)

      call MatSetOptionsPrefix(solver%Jpre,"flow_",ierr);CHKERRQ(ierr)

            if (solver%Jpre_mat_type == solver%J_mat_type) then
        solver%J = solver%Jpre
            else
              call DiscretizationCreateJacobian(pm%realization%discretization, &
                                                NFLOWDOF, &
                                                solver%J_mat_type, &
                                                solver%J, &
                                                option)

              call MatSetOptionsPrefix(solver%J,"flow_",ierr);CHKERRQ(ierr)
      endif

      if (solver%use_galerkin_mg) then
        call DiscretizationCreateInterpolation( &
                       pm%realization%discretization,NFLOWDOF, &
                       solver%interpolation, &
                       solver%galerkin_mg_levels_x, &
                       solver%galerkin_mg_levels_y, &
                       solver%galerkin_mg_levels_z, &
                       option)
      endif
    
      if (solver%J_mat_type == MATMFFD) then
        call MatCreateSNESMF(solver%snes,solver%J,ierr);CHKERRQ(ierr)
      endif

      ! by default turn off line search
      call SNESGetLineSearch(solver%snes, linesearch, ierr);CHKERRQ(ierr)
      call SNESLineSearchSetType(linesearch, SNESLINESEARCHBASIC,  &
                                  ierr);CHKERRQ(ierr)
      ! Have PETSc do a SNES_View() at the end of each solve if 
      ! verbosity > 0.
      if (option%verbosity >= 2) then
        string = '-flow_snes_view'
        call PetscOptionsInsertString(PETSC_NULL_OPTIONS,string, &
                                      ierr);CHKERRQ(ierr)
      endif

      ! If we are using a structured grid, set the corresponding flow 
      ! DA as the DA for the PCEXOTIC preconditioner, in case we 
      ! choose to use it. The PCSetDA() call is ignored if the 
      ! PCEXOTIC preconditioner is no used.  We need to put this call 
      ! after SolverCreateSNES() so that KSPSetFromOptions() will 
      ! already have been called.  I also note that this 
      ! preconditioner is intended only for the flow 
      ! solver.  --RTM
      if (pm%realization%discretization%itype == STRUCTURED_GRID) then
        call PCSetDM(solver%pc, &
                     pm%realization%discretization%dm_nflowdof%dm, &
                     ierr);CHKERRQ(ierr)
      endif

      call SNESSetConvergenceTest(solver%snes, &
#if defined(USE_PM_AS_PETSC_CONTEXT)
                                  PMCheckConvergence, &
                                  this%pm_ptr%pm, &
#else
                                        PMCheckConvergencePtr, &
                                        this%pm_ptr, &
#endif
                                  PETSC_NULL_FUNCTION,ierr);CHKERRQ(ierr)

      if (pm%check_post_convergence) then
        call SNESLineSearchSetPostCheck(linesearch, &
#if defined(USE_PM_AS_PETSC_CONTEXT)
                                        PMCheckUpdatePost, &
                                        this%pm_ptr%pm, &
#else
                                        PMCheckUpdatePostPtr, &
                                        this%pm_ptr, &
#endif
                                        ierr);CHKERRQ(ierr)
        !geh: it is possible that the other side has not been set
        pm%check_post_convergence = PETSC_TRUE
      endif
                                  
      add_pre_check = PETSC_FALSE
      select type(pm)
        class is(pm_th_type)
          add_pre_check = PETSC_TRUE
      end select

        if (add_pre_check) then
#if defined(USE_PM_AS_PETSC_CONTEXT)
            call SNESLineSearchSetPreCheck(linesearch, &
                                           PMCheckUpdatePre, &
                                           this%pm_ptr%pm, &
                                           ierr);CHKERRQ(ierr)
#else
            call SNESLineSearchSetPreCheck(linesearch, &
                                           PMCheckUpdatePrePtr, &
                                           this%pm_ptr, &
                                           ierr);CHKERRQ(ierr)
#endif
        endif

      call PrintMsg(option,"  Finished setting up FLOW SNES ")
  end select
  
  call SNESSetFunction(this%timestepper%solver%snes, &
                       this%pm_ptr%pm%residual_vec, &
#if defined(USE_PM_AS_PETSC_CONTEXT)
                       PMResidual, &
                       this%pm_ptr%pm, &
#else
                       PMResidualPtr, &
                       this%pm_ptr, &
#endif
                       ierr);CHKERRQ(ierr)
  call SNESSetJacobian(this%timestepper%solver%snes, &
                       solver%J, &
                       solver%Jpre, &
#if defined(USE_PM_AS_PETSC_CONTEXT)
                       PMJacobian, &
                       this%pm_ptr%pm, &
#else
                       PMJacobianPtr, &
                       this%pm_ptr, &
#endif
                       ierr);CHKERRQ(ierr)
  call SolverSetSNESOptions(solver,option)

end subroutine PMCSubsurfaceSetupSolvers_TimestepperBE

! ************************************************************************** !

subroutine PMCSubsurfaceGetAuxData(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 10/24/13
  ! 

  implicit none

  class(pmc_subsurface_type) :: this

end subroutine PMCSubsurfaceGetAuxData

! ************************************************************************** !

subroutine PMCSubsurfaceSetAuxData(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 10/24/13
  ! 

  implicit none

  class(pmc_subsurface_type) :: this

end subroutine PMCSubsurfaceSetAuxData

! ************************************************************************** !
!
! PMCSubsurfaceFinalizeRun: Finalizes the time stepping
! author: Glenn Hammond
! date: 03/18/13
!
! ************************************************************************** !
recursive subroutine PMCSubsurfaceFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  ! 

  use Option_module
  
  implicit none
  
  class(pmc_subsurface_type) :: this
  
#ifdef DEBUG
  call PrintMsg(this%option,'PMCSubsurface%FinalizeRun()')
#endif
  
  nullify(this%realization)
  
end subroutine PMCSubsurfaceFinalizeRun

! ************************************************************************** !

subroutine CPRWorkersCreate(pm, solver, option)
  ! 
  ! create all the worker/storage matrices/vectors that will be needed for the
  ! cpr preconditioner
  !
  ! Author: Daniel Stone
  ! Date: Oct 2017 - March 2018
  ! 

  use PM_Subsurface_Flow_class
  use Solver_module
  use Option_module
  use Discretization_module
  
  implicit none

  class(pm_subsurface_flow_type) :: pm
  class(solver_type) :: solver
  class(option_type) :: option
  MatType :: cpr_ap_mat_type



  cpr_ap_mat_type =  MATAIJ
  call DiscretizationCreateJacobian(pm%realization%discretization, &
                                    ONEDOF, &
                                    cpr_ap_mat_type, &
                                    solver%cprstash%Ap, &
                                    option)

  call DiscretizationCreateVector(pm%realization%discretization, &
                                  NFLOWDOF, solver%cprstash%T1r, &
                                  GLOBAL, option)
  call DiscretizationCreateVector(pm%realization%discretization, &
                                  NFLOWDOF, solver%cprstash%r2, &
                                  GLOBAL, option)

  call DiscretizationCreateVector(pm%realization%discretization, &
                                  ONEDOF, solver%cprstash%s,  &
                                  GLOBAL, option)
  call DiscretizationCreateVector(pm%realization%discretization, &
                                  ONEDOF, solver%cprstash%z, &
                                  GLOBAL, option)

  call DiscretizationCreateVector(pm%realization%discretization, &
                                  NFLOWDOF, solver%cprstash%factors1vec, &
                                  GLOBAL, option)
  call DiscretizationCreateVector(pm%realization%discretization, &
                                  NFLOWDOF, solver%cprstash%factors2vec, &
                                  GLOBAL, option)
end subroutine CPRWorkersCreate

! ************************************************************************** !

subroutine PMCSubsurfaceStrip(this)
  !
  ! Deallocates members of PMC Subsurface.
  !
  ! Author: Glenn Hammond
  ! Date: 01/13/14
  
  implicit none
  
  class(pmc_subsurface_type) :: this

  call PMCBaseStrip(this)
  nullify(this%realization)

end subroutine PMCSubsurfaceStrip

! ************************************************************************** !

recursive subroutine PMCSubsurfaceDestroy(this)
  ! 
  ! ProcessModelCouplerDestroy: Deallocates a process_model_coupler object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Option_module

  implicit none
  
  class(pmc_subsurface_type) :: this
  
#ifdef DEBUG
  call PrintMsg(this%option,'PMCSubsurface%Destroy()')
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
  
  call PMCSubsurfaceStrip(this)
  
end subroutine PMCSubsurfaceDestroy

! ************************************************************************** !
  
end module PMC_Subsurface_class
