module PMC_Subsurface_class

  use PMC_Base_class
  use Realization_Subsurface_class

  use PFLOTRAN_Constants_module

#include "petsc/finclude/petscmat.h"
  use petscmat
#include "petsc/finclude/petscsys.h"
  use petscsys

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
  use Timestepper_TS_class

  implicit none

  class(pmc_subsurface_type) :: this

  type(option_type), pointer :: option

  option => this%option

  if (associated(this%timestepper)) then
    select type(ts => this%timestepper)
      class is(timestepper_TS_type)
        call PMCSubsurfaceSetupSolvers_TS(this)
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
#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Convergence_module
  use Discretization_module
  use Realization_Subsurface_class
  use Option_module
  use PMC_Base_class
  use PM_Base_Pointer_module
  use PM_Base_class
  use PM_Subsurface_Flow_class
  use PM_General_class
  use PM_WIPP_Flow_class
  use PM_Richards_class
  use PM_TH_class
  use PM_RT_class
  use PM_NWT_class
  use PM_Waste_Form_class
  use PM_UFD_Decay_class
  use PM_TOilIms_class
  use PM_TOWG_class
  use Secondary_Continuum_module, only : SecondaryRTUpdateIterate  
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
  PetscInt :: itransport, RT, NWT
  PetscInt :: trans_coupling
  PetscErrorCode :: ierr

#ifdef DEBUG
  call PrintMsg(this%option,'PMCSubsurface%SetupSolvers()')
#endif

  option => this%option
  solver => this%timestepper%solver
  
  check_update = PETSC_FALSE
  itransport = 0
  trans_coupling = 10000 ! set to high value (not 0 or 1)
  RT = 1
  NWT = 2

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
          option%iflowmode /= RICHARDS_MODE) then

        option%io_buffer = 'AIJ matrix not supported for current &
          &mode: '// option%flowmode
        call PrintErrMsg(option)
      endif
      if (OptionPrintToScreen(option)) then
        write(*,'(" number of dofs = ",i3,", number of &
                  &phases = ",i3,i2)') option%nflowdof,option%nphase
        select case(option%iflowmode)
          case(FLASH2_MODE)
            write(*,'(" mode = FLASH2: p, T, s/X")')
          case(MPH_MODE)
            write(*,'(" mode = MPH: p, T, s/X")')
          case(IMS_MODE)
            write(*,'(" mode = IMS: p, T, s")')
          case(MIS_MODE)
            write(*,'(" mode = MIS: p, Xs")')
          case(TH_MODE)
            write(*,'(" mode = TH: p, T")')
          case(RICHARDS_MODE)
            write(*,'(" mode = Richards: p")')
          case(G_MODE) 
            write(*,'(" mode = General: p, sg/X, T")')
          case(WF_MODE) 
            write(*,'(" mode = WIPP Flow: p, sg")')
          case(TOIL_IMS_MODE)   
        end select
      endif

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
        class is(pm_richards_type)
          if (Initialized(pm%pressure_dampening_factor) .or. &
              Initialized(pm%saturation_change_limit)) then
              add_pre_check = PETSC_TRUE
          endif
        class is(pm_general_type)
              add_pre_check = PETSC_TRUE
        class is(pm_wippflo_type)
              add_pre_check = PETSC_TRUE
        class is(pm_toil_ims_type)
              add_pre_check = PETSC_TRUE
        class is(pm_towg_type)
              add_pre_check = PETSC_TRUE
        class is(pm_th_type)
          if (Initialized(pm%pressure_dampening_factor) .or. &
              Initialized(pm%pressure_change_limit) .or. &
              Initialized(pm%temperature_change_limit)) then
              add_pre_check = PETSC_TRUE
          endif
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
  ! ----- subsurface reactive transport
    class is(pm_rt_type)
      itransport = RT
      check_post_convergence = pm%check_post_convergence
      trans_coupling = option%transport%reactive_transport_coupling
      check_update = pm%realization%reaction%check_update
      discretization => pm%realization%discretization
      realization => pm%realization
      if (OptionPrintToScreen(option)) then
        write(*,'(" mode = Reactive Transport")')
      endif
  ! ----- nuclear waste transport
    class is(pm_nwt_type)
      itransport = NWT
      check_post_convergence = pm%controls%check_post_convergence
      trans_coupling = option%transport%nw_transport_coupling
      check_update = pm%controls%check_update
      discretization => pm%realization%discretization
      realization => pm%realization
      if (OptionPrintToScreen(option)) then
        write(*,'(" mode = Nuclear Waste Transport")')
      endif
  end select
  
  if ( (itransport == RT) .or. (itransport == NWT) ) then
    call PrintMsg(option,"  Beginning setup of TRAN SNES ")
    call SNESSetOptionsPrefix(solver%snes, "tran_",ierr);CHKERRQ(ierr)
    call SolverCheckCommandLine(solver)
    
    if (trans_coupling == GLOBAL_IMPLICIT) then
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

      call DiscretizationCreateJacobian(discretization, &
                                        NTRANDOF, &
                                        solver%Jpre_mat_type, &
                                        solver%Jpre,option)
    else
      solver%J_mat_type = MATAIJ
      solver%Jpre_mat_type = MATAIJ
      call DiscretizationCreateJacobian(discretization, &
                                        ONEDOF, &
                                        solver%Jpre_mat_type, &
                                        solver%Jpre,option)
    endif

    call MatSetOptionsPrefix(solver%Jpre,"tran_",ierr);CHKERRQ(ierr)

    if (solver%Jpre_mat_type == solver%J_mat_type) then
      solver%J = solver%Jpre
    else
      call DiscretizationCreateJacobian(discretization, &
                                        NTRANDOF, &
                                        solver%J_mat_type, &
                                        solver%J, &
                                        option)

      call MatSetOptionsPrefix(solver%J,"tran_",ierr);CHKERRQ(ierr)
    endif

    if (solver%use_galerkin_mg) then
      call DiscretizationCreateInterpolation( &
                       discretization,NTRANDOF, &
                       solver%interpolation, &
                       solver%galerkin_mg_levels_x, &
                       solver%galerkin_mg_levels_y, &
                       solver%galerkin_mg_levels_z, &
                       option)
    endif

    if (trans_coupling == GLOBAL_IMPLICIT) then

      if (solver%J_mat_type == MATMFFD) then
        call MatCreateSNESMF(solver%snes,solver%J, &
                             ierr);CHKERRQ(ierr)
      endif
      
      ! this could be changed in the future if there is a way to 
      ! ensure that the linesearch update does not perturb 
      ! concentrations negative.
      call SNESGetLineSearch(solver%snes, linesearch, &
                             ierr);CHKERRQ(ierr)
      call SNESLineSearchSetType(linesearch, SNESLINESEARCHBASIC,  &
                                  ierr);CHKERRQ(ierr)
      
      if (option%use_mc .and. (itransport == RT)) then
        call SNESLineSearchSetPostCheck(linesearch, &
                              SecondaryRTUpdateIterate, &
                              realization,ierr);CHKERRQ(ierr)
      endif
      
      ! Have PETSc do a SNES_View() at the end of each solve if 
      ! verbosity > 0.
      if (option%verbosity >= 2) then
        string = '-tran_snes_view'
        call PetscOptionsInsertString(PETSC_NULL_OPTIONS, &
                                      string, ierr);CHKERRQ(ierr)
      endif

    endif

    if (trans_coupling == GLOBAL_IMPLICIT) then
      call SNESSetConvergenceTest(solver%snes, &
#if defined(USE_PM_AS_PETSC_CONTEXT)
                                  PMCheckConvergence, &
                                  this%pm_ptr%pm, &
#else
                                  PMCheckConvergencePtr, &
                                  this%pm_ptr, &
#endif
                                  PETSC_NULL_FUNCTION,ierr);CHKERRQ(ierr)
    endif
    if (this%pm_ptr%pm%print_EKG .or. option%use_mc .or. &
        check_post_convergence) then
      call SNESLineSearchSetPostCheck(linesearch, &
#if defined(USE_PM_AS_PETSC_CONTEXT)
                                      PMCheckUpdatePost, &
                                      this%pm_ptr%pm, &
#else
                                      PMCheckUpdatePostPtr, &
                                      this%pm_ptr, &
#endif
                                      ierr);CHKERRQ(ierr)
      if (this%pm_ptr%pm%print_EKG) then
        check_post_convergence = PETSC_TRUE
      endif
    endif
    if (check_update) then
      call SNESLineSearchSetPreCheck(linesearch, &
#if defined(USE_PM_AS_PETSC_CONTEXT)
                                     PMCheckUpdatePre, &
                                     this%pm_ptr%pm, &
#else
                                     PMCheckUpdatePrePtr, &
                                     this%pm_ptr, &
#endif
                                     ierr);CHKERRQ(ierr)
    endif
    call PrintMsg(option,"  Finished setting up TRAN SNES ")
      
  endif ! RT or NWT
      

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

subroutine PMCSubsurfaceSetupSolvers_TS(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 06/20/2018
  ! 
#include "petsc/finclude/petscts.h"
#include "petsc/finclude/petscsnes.h"
  use petscts
  use petscsnes
  use Convergence_module
  use Discretization_module
  use Option_module
  use PMC_Base_class
  use PM_Base_Pointer_module
  use PM_Base_class
  use PM_Subsurface_Flow_class
  use PM_Richards_class
  use Solver_module
  use Timestepper_Base_class
  use Timestepper_TS_class

  implicit none

  class(pmc_subsurface_type) :: this

  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  SNESLineSearch :: linesearch  
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: add_pre_check
  SNES :: snes
  PetscErrorCode :: ierr

  option => this%option

  select type(ts => this%timestepper)
    class is(timestepper_TS_type)
      solver => ts%solver
    class default
      call PrintErrMsg(option,"Attempting to set up PETSc TS when" // &
       " timestepper is not of TS type")
  end select

  call SolverCreateTS(solver,option%mycomm)

  call TSSetEquationType(solver%ts,TS_EQ_IMPLICIT, &
                        ierr);CHKERRQ(ierr)

  call TSSetType(solver%ts, TSBEULER, ierr); CHKERRQ(ierr)

  ! set solver pointer within pm for convergence purposes
  call this%pm_ptr%pm%SetSolver(solver)

  select type(pm => this%pm_ptr%pm)

    class is (pm_subsurface_flow_type)
      call PrintMsg(option,"  Beginning setup of FLOW SNES ")

      select case(option%iflowmode)
        case(RICHARDS_TS_MODE,TH_TS_MODE)
        case default
          option%io_buffer = 'Timestepper TS unsupported for mode: '// option%flowmode
          call PrintErrMsg(option)
        end select

        if (OptionPrintToScreen(option)) then
          write(*,'(" number of dofs = ",i3,", number of &
                    &phases = ",i3,i2)') option%nflowdof,option%nphase
          select case(option%iflowmode)
            case(RICHARDS_TS_MODE)
              write(*,'(" mode = Richards: p")')
            case(TH_TS_MODE)
              write(*,'(" mode = TH: p, T")')
          end select
        endif

        call TSSetOptionsPrefix(solver%ts, "flow_",ierr);CHKERRQ(ierr)
        call TSSetFromOptions(solver%ts,ierr);CHKERRQ(ierr)

        call SolverCheckCommandLine(solver)

        solver%Jpre_mat_type = MATBAIJ

        call DiscretizationCreateJacobian(pm%realization%discretization, &
                                          NFLOWDOF, &
                                          solver%Jpre_mat_type, &
                                          solver%Jpre, &
                                          option)

        call MatSetOptionsPrefix(solver%Jpre,"flow_",ierr);CHKERRQ(ierr)

        if (solver%J_mat_type /= MATMFFD) then
          solver%J = solver%Jpre
        endif

      call TSSetIFunction(solver%ts, &
                          this%pm_ptr%pm%residual_vec, &
                          PMIFunctionPtr, &
                          this%pm_ptr, &
                          ierr);CHKERRQ(ierr)

      call TSSetIJacobian(solver%ts, &
                          solver%J, &
                          solver%Jpre, &
                          PMIJacobianPtr, &
                          this%pm_ptr, &
                          ierr);CHKERRQ(ierr)

      call TSGetSNES(solver%ts,snes,ierr); CHKERRQ(ierr)

      call SNESSetConvergenceTest(snes, &
                                  PMCheckConvergencePtr, &
                                  this%pm_ptr, &
                                  PETSC_NULL_FUNCTION,ierr);CHKERRQ(ierr)

      call PrintMsg(option,"  Finished setting up FLOW SNES ")
  end select

end subroutine PMCSubsurfaceSetupSolvers_TS

! ************************************************************************** !

subroutine PMCSubsurfaceGetAuxData(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 10/24/13
  ! 

  implicit none

  class(pmc_subsurface_type) :: this

  if (this%option%surf_flow_on) call PMCSubsurfaceGetAuxDataFromSurf(this)
  if (this%option%ngeomechdof > 0) call PMCSubsurfaceGetAuxDataFromGeomech(this)

end subroutine PMCSubsurfaceGetAuxData

! ************************************************************************** !

subroutine PMCSubsurfaceSetAuxData(this)
  ! 
  ! Author: Gautam Bisht
  ! Date: 10/24/13
  ! 

  implicit none

  class(pmc_subsurface_type) :: this

  if (this%option%surf_flow_on) call PMCSubsurfaceSetAuxDataForSurf(this)
  if (this%option%ngeomechdof > 0) call PMCSubsurfaceSetAuxDataForGeomech(this)

end subroutine PMCSubsurfaceSetAuxData

! ************************************************************************** !

subroutine PMCSubsurfaceGetAuxDataFromSurf(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 08/22/13
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Connection_module
  use Coupler_module
  use Field_module
  use Grid_module
  use Option_module
  use Patch_module
!  use Realization_Base_class
  use Realization_Subsurface_class
  use String_module
  use EOS_Water_module

  implicit none

  class(pmc_subsurface_type) :: this
  
  class(realization_subsurface_type), pointer :: realization
  type (patch_type),pointer :: patch
  type (grid_type),pointer :: grid
  type (coupler_list_type), pointer :: coupler_list
  type (coupler_type), pointer :: coupler
  type (option_type), pointer :: option
  type (field_type),pointer :: field
  type (connection_set_type), pointer :: cur_connection_set
  PetscBool :: coupler_found
  PetscInt :: iconn
  PetscReal :: den
  PetscReal :: dt
  PetscReal :: surfpress
  PetscReal :: dum1
  PetscReal, pointer :: mflux_p(:)
  PetscReal, pointer :: hflux_p(:)
  PetscReal, pointer :: head_p(:)
  PetscReal, pointer :: temp_p(:)
  PetscErrorCode :: ierr

#ifdef DEBUG
  print *, 'PMCSubsurfaceGetAuxData()'
#endif

  dt = this%option%surf_subsurf_coupling_flow_dt

  if (associated(this%sim_aux)) then

    select type (pmc => this)
      class is (pmc_subsurface_type)

      if (this%sim_aux%subsurf_mflux_exchange_with_surf /= PETSC_NULL_VEC) then
        ! PETSc Vector to store relevant mass-flux data between
        ! surface-subsurface model exists

        patch      => pmc%realization%patch
        grid       => pmc%realization%discretization%grid
        field      => pmc%realization%field
        option     => pmc%realization%option

        select case(this%option%iflowmode)
          case (RICHARDS_MODE,RICHARDS_TS_MODE)
            call VecScatterBegin(pmc%sim_aux%surf_to_subsurf, &
                                 pmc%sim_aux%surf_mflux_exchange_with_subsurf, &
                                 pmc%sim_aux%subsurf_mflux_exchange_with_surf, &
                                 INSERT_VALUES,SCATTER_FORWARD, &
                                 ierr);CHKERRQ(ierr)
            call VecScatterEnd(pmc%sim_aux%surf_to_subsurf, &
                               pmc%sim_aux%surf_mflux_exchange_with_subsurf, &
                               pmc%sim_aux%subsurf_mflux_exchange_with_surf, &
                               INSERT_VALUES,SCATTER_FORWARD, &
                               ierr);CHKERRQ(ierr)

            call VecScatterBegin(pmc%sim_aux%surf_to_subsurf, &
                                 pmc%sim_aux%surf_head, &
                                 pmc%sim_aux%subsurf_pres_top_bc, &
                                 INSERT_VALUES,SCATTER_FORWARD, &
                                 ierr);CHKERRQ(ierr)
            call VecScatterEnd(pmc%sim_aux%surf_to_subsurf, &
                               pmc%sim_aux%surf_head, &
                               pmc%sim_aux%subsurf_pres_top_bc, &
                               INSERT_VALUES,SCATTER_FORWARD, &
                               ierr);CHKERRQ(ierr)
            call EOSWaterdensity(option%reference_temperature, &
                                 option%reference_pressure,den,dum1,ierr)

#if 0
            coupler_list => patch%source_sink_list
            coupler => coupler_list%first
            do
              if (.not.associated(coupler)) exit

              ! FLOW
              if (associated(coupler%flow_aux_real_var)) then

                ! Find the BC from the list of BCs
                if (StringCompare(coupler%name,'from_surface_ss')) then
                  coupler_found = PETSC_TRUE
                  
                  call VecGetArrayF90(pmc%sim_aux%subsurf_mflux_exchange_with_surf, &
                                      mflux_p,ierr);CHKERRQ(ierr)
                  do iconn = 1,coupler%connection_set%num_connections
                    !coupler%flow_aux_real_var(ONE_INTEGER,iconn) = -mflux_p(iconn)/dt*den
                  enddo
                  call VecRestoreArrayF90(pmc%sim_aux%subsurf_mflux_exchange_with_surf, &
                                          mflux_p,ierr);CHKERRQ(ierr)

                  call VecSet(pmc%sim_aux%surf_mflux_exchange_with_subsurf,0.d0, &
                              ierr);CHKERRQ(ierr)
                endif
              endif

              coupler => coupler%next
            enddo
#endif

            coupler_list => patch%boundary_condition_list
            coupler => coupler_list%first
            do
              if (.not.associated(coupler)) exit

              ! FLOW
              if (associated(coupler%flow_aux_real_var)) then
                ! Find the BC from the list of BCs
                if (StringCompare(coupler%name,'from_surface_bc')) then
                  coupler_found = PETSC_TRUE
                  call VecGetArrayF90(pmc%sim_aux%subsurf_pres_top_bc, &
                                      head_p,ierr);CHKERRQ(ierr)
                  do iconn = 1,coupler%connection_set%num_connections
                    surfpress = head_p(iconn)*(abs(option%gravity(3)))*den + &
                                option%reference_pressure
                    coupler%flow_aux_real_var(RICHARDS_PRESSURE_DOF,iconn) = &
                    surfpress
                  enddo
                  call VecRestoreArrayF90(pmc%sim_aux%subsurf_pres_top_bc, &
                                          head_p,ierr);CHKERRQ(ierr)
                endif
              endif
              coupler => coupler%next
            enddo

          case (TH_MODE,TH_TS_MODE)
            call VecScatterBegin(pmc%sim_aux%surf_to_subsurf, &
                                 pmc%sim_aux%surf_head, &
                                 pmc%sim_aux%subsurf_pres_top_bc, &
                                 INSERT_VALUES,SCATTER_FORWARD, &
                                 ierr);CHKERRQ(ierr)
            call VecScatterEnd(pmc%sim_aux%surf_to_subsurf, &
                               pmc%sim_aux%surf_head, &
                               pmc%sim_aux%subsurf_pres_top_bc, &
                               INSERT_VALUES,SCATTER_FORWARD, &
                               ierr);CHKERRQ(ierr)

            call VecScatterBegin(pmc%sim_aux%surf_to_subsurf, &
                                 pmc%sim_aux%surf_temp, &
                                 pmc%sim_aux%subsurf_temp_top_bc, &
                                 INSERT_VALUES,SCATTER_FORWARD, &
                                 ierr);CHKERRQ(ierr)
            call VecScatterEnd(pmc%sim_aux%surf_to_subsurf, &
                               pmc%sim_aux%surf_temp, &
                               pmc%sim_aux%subsurf_temp_top_bc, &
                               INSERT_VALUES,SCATTER_FORWARD, &
                               ierr);CHKERRQ(ierr)

            call VecScatterBegin(pmc%sim_aux%surf_to_subsurf, &
                                 pmc%sim_aux%surf_hflux_exchange_with_subsurf, &
                                 pmc%sim_aux%subsurf_mflux_exchange_with_surf, &
                                 INSERT_VALUES,SCATTER_FORWARD, &
                                 ierr);CHKERRQ(ierr)
            call VecScatterEnd(pmc%sim_aux%surf_to_subsurf, &
                               pmc%sim_aux%surf_hflux_exchange_with_subsurf, &
                               pmc%sim_aux%subsurf_mflux_exchange_with_surf, &
                               INSERT_VALUES,SCATTER_FORWARD, &
                               ierr);CHKERRQ(ierr)

            coupler_list => patch%boundary_condition_list
            coupler => coupler_list%first
            do
              if (.not.associated(coupler)) exit

              ! FLOW
              if (associated(coupler%flow_aux_real_var)) then
                ! Find the BC from the list of BCs
                if (StringCompare(coupler%name,'from_surface_bc')) then
                  coupler_found = PETSC_TRUE

                  call VecGetArrayF90(pmc%sim_aux%subsurf_pres_top_bc, &
                                      head_p,ierr);CHKERRQ(ierr)
                  call VecGetArrayF90(pmc%sim_aux%subsurf_temp_top_bc, &
                                      temp_p,ierr);CHKERRQ(ierr)

                  do iconn = 1,coupler%connection_set%num_connections

                    ! The pressure value needed to computed density should
                    ! be surf_press and not reference_pressure. But,
                    ! surf_pressure depends on density.
                    call EOSWaterdensity(temp_p(iconn), option%reference_pressure, &
                                         den,dum1,ierr)

                    surfpress = head_p(iconn)*(abs(option%gravity(3)))*den + &
                                option%reference_pressure
                    coupler%flow_aux_real_var(TH_PRESSURE_DOF,iconn) = &
                      surfpress
                    coupler%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
                      temp_p(iconn)
                  enddo

                  call VecRestoreArrayF90(pmc%sim_aux%subsurf_pres_top_bc, &
                                          head_p,ierr);CHKERRQ(ierr)
                  call VecRestoreArrayF90(pmc%sim_aux%subsurf_temp_top_bc, &
                                      temp_p,ierr);CHKERRQ(ierr)
                endif
              endif

              if (StringCompare(coupler%name,'from_atm_subsurface_bc')) then
                coupler_found = PETSC_TRUE

                call VecGetArrayF90(pmc%sim_aux%subsurf_mflux_exchange_with_surf, &
                                    mflux_p,ierr);CHKERRQ(ierr)

                do iconn = 1,coupler%connection_set%num_connections
                  coupler%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn) = &
                    mflux_p(iconn)
                enddo

                call VecRestoreArrayF90(pmc%sim_aux%subsurf_mflux_exchange_with_surf, &
                                    mflux_p,ierr);CHKERRQ(ierr)
              endif

              coupler => coupler%next
            enddo

          case default
            this%option%io_buffer='PMCSubsurfaceGetAuxData() not supported for this mode.'
            call PrintErrMsg(this%option)

        end select

        if ( .not. coupler_found) then
          option%io_buffer = 'Coupler not found in PMCSubsurfaceGetAuxData()'
          call PrintErrMsg(option)
        endif
      endif

    end select

  endif ! if (associated(this%sim_aux))

end subroutine PMCSubsurfaceGetAuxDataFromSurf

! ************************************************************************** !

subroutine PMCSubsurfaceSetAuxDataForSurf(this)
  ! 
  ! This routine sets auxiliary to be exchanged between process-models.
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 08/21/13
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Grid_module
  use String_module
  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Coupler_module
  use Field_module
  use Connection_module
  use Realization_Base_class
  use EOS_Water_module

  implicit none

  class(pmc_subsurface_type) :: this
  
  class(realization_subsurface_type), pointer :: realization
  type (patch_type),pointer :: patch
  type (grid_type),pointer :: grid
  type (coupler_list_type), pointer :: coupler_list
  type (coupler_type), pointer :: coupler
  type (option_type), pointer :: option
  type (field_type),pointer :: field
  type (connection_set_type), pointer :: cur_connection_set
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iconn
  PetscInt :: istart
  PetscInt :: iend
  PetscReal :: den
  PetscReal :: dum1
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal, pointer :: pres_top_bc_p(:)
  PetscReal, pointer :: temp_top_bc_p(:)
  PetscReal, pointer :: head_p(:)
  PetscErrorCode :: ierr

#ifdef DEBUG
  print *, 'PMCSubsurfaceSetAuxData()'
#endif

  if (associated(this%sim_aux)) then

    select type (pmc => this)
      class is (pmc_subsurface_type)

        if (this%sim_aux%subsurf_pres_top_bc/= PETSC_NULL_VEC) then
          ! PETSc Vector to store relevant subsurface-flow data for
          ! surface-flow model exists

          patch      => pmc%realization%patch
          grid       => pmc%realization%discretization%grid
          field      => pmc%realization%field
          option     => pmc%realization%option

          call EOSWaterdensity(option%reference_temperature, option%reference_pressure, &
                               den,dum1,ierr)
          coupler_list => patch%boundary_condition_list
          coupler => coupler_list%first
          do
            if (.not.associated(coupler)) exit

            ! FLOW
            if (associated(coupler%flow_aux_real_var)) then

              ! Find the BC from the list of BCs
              if (StringCompare(coupler%name,'from_surface_bc')) then
                select case(this%option%iflowmode)
                  case (RICHARDS_MODE,RICHARDS_TS_MODE)
                    call VecGetArrayF90(this%sim_aux%subsurf_pres_top_bc, &
                                        pres_top_bc_p,ierr);CHKERRQ(ierr)
                    do iconn = 1,coupler%connection_set%num_connections
                      pres_top_bc_p(iconn) = &
                        coupler%flow_aux_real_var(RICHARDS_PRESSURE_DOF,iconn)
                    enddo
                    call VecRestoreArrayF90(this%sim_aux%subsurf_pres_top_bc, &
                                            pres_top_bc_p,ierr);CHKERRQ(ierr)
                  case (TH_MODE,TH_TS_MODE)
                    call VecGetArrayF90(this%sim_aux%subsurf_pres_top_bc, &
                                        pres_top_bc_p,ierr);CHKERRQ(ierr)
                    call VecGetArrayF90(this%sim_aux%subsurf_temp_top_bc, &
                                        temp_top_bc_p,ierr);CHKERRQ(ierr)

                    do iconn = 1,coupler%connection_set%num_connections
                      pres_top_bc_p(iconn) = &
                        coupler%flow_aux_real_var(TH_PRESSURE_DOF,iconn)
                      temp_top_bc_p(iconn) = &
                        coupler%flow_aux_real_var(TH_TEMPERATURE_DOF,iconn)
                    enddo

                    call VecRestoreArrayF90(this%sim_aux%subsurf_pres_top_bc, &
                                            pres_top_bc_p,ierr);CHKERRQ(ierr)
                    call VecRestoreArrayF90(this%sim_aux%subsurf_temp_top_bc, &
                                            temp_top_bc_p,ierr);CHKERRQ(ierr)
                    case default
                      option%io_buffer = 'PMCSubsurfaceGetAuxData() not ' // &
                        'supported in this FLOW_MODE'
                      call PrintErrMsg(option)
                end select
              endif
            endif

            coupler => coupler%next
          enddo

        endif
    end select

  endif

end subroutine PMCSubsurfaceSetAuxDataForSurf

! ************************************************************************** !

subroutine PMCSubsurfaceGetAuxDataFromGeomech(this)
  !
  ! This routine updates subsurface data from geomechanics process model.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/04/14

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Discretization_module, only : DiscretizationLocalToGlobal
  use Field_module
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use PFLOTRAN_Constants_module
  use Material_Aux_class
  use Material_module
  use Variables_module, only : POROSITY

  implicit none

  class (pmc_subsurface_type) :: this

  type(grid_type), pointer :: subsurf_grid
  type(option_type), pointer :: option
  type(field_type), pointer :: subsurf_field

  PetscScalar, pointer :: sim_por_p(:)
  PetscScalar, pointer :: work_loc_p(:)
  class(material_auxvar_type), pointer :: subsurf_material_auxvars(:)

  PetscInt :: local_id
  PetscInt :: ghosted_id

  PetscErrorCode :: ierr
  PetscViewer :: viewer

#ifdef GEOMECH_DEBUG
  print *, 'PMCSubsurfaceGetAuxDataFromGeomech()'
#endif

  if (associated(this%sim_aux)) then
    select type (pmc => this)
      class is (pmc_subsurface_type)
        option        => pmc%option
        subsurf_grid  => pmc%realization%discretization%grid
        subsurf_field => pmc%realization%field
        subsurf_material_auxvars => pmc%realization%patch%aux%Material%auxvars

        if (option%geomech_subsurf_coupling == GEOMECH_TWO_WAY_COUPLED) then

          call VecGetArrayF90(pmc%sim_aux%subsurf_por,sim_por_p,  &
                              ierr);CHKERRQ(ierr)
          call VecGetArrayF90(subsurf_field%work_loc,work_loc_p,  &
                              ierr);CHKERRQ(ierr)

          do local_id = 1, subsurf_grid%nlmax
            ghosted_id = subsurf_grid%nL2G(local_id)
            work_loc_p(ghosted_id) = sim_por_p(local_id)
          enddo
            
          call VecRestoreArrayF90(pmc%sim_aux%subsurf_por,sim_por_p,  &
                                  ierr);CHKERRQ(ierr)
          call VecRestoreArrayF90(subsurf_field%work_loc,work_loc_p,  &
                              ierr);CHKERRQ(ierr)

          call DiscretizationLocalToGlobal(pmc%realization%discretization, &
                                          subsurf_field%work_loc, &
                                          subsurf_field%porosity_geomech_store, &
                                          ONEDOF)
#ifdef GEOMECH_DEBUG
          call PetscViewerASCIIOpen(pmc%realization%option%mycomm, &
                                    'porosity_geomech_store.out', &
                                    viewer,ierr);CHKERRQ(ierr)
          call VecView(subsurf_field%porosity_geomech_store,viewer,ierr);CHKERRQ(ierr)
          call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif


        endif
    end select
  endif

end subroutine PMCSubsurfaceGetAuxDataFromGeomech

! ************************************************************************** !

subroutine PMCSubsurfaceSetAuxDataForGeomech(this)
  !
  ! This routine sets auxiliary needed by geomechanics process model.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/04/14

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Material_Aux_class
  use PFLOTRAN_Constants_module

  implicit none

  class (pmc_subsurface_type) :: this

  type(grid_type), pointer :: subsurf_grid
  type(option_type), pointer :: option
  type(field_type), pointer :: subsurf_field

  PetscScalar, pointer :: xx_loc_p(:)
  PetscScalar, pointer :: pres_p(:)
  PetscScalar, pointer :: temp_p(:)
  PetscScalar, pointer :: sub_por_loc_p(:)
  PetscScalar, pointer :: sim_por0_p(:)
  PetscScalar, pointer :: sim_perm0_p(:) !DANNY - added this 11/7/16
  
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: pres_dof
  PetscInt :: temp_dof

  class(material_auxvar_type), pointer :: material_auxvars(:)

  PetscErrorCode :: ierr

#ifdef GEOMECH_DEBUG
  print *, 'PMCSubsurfaceSetAuxDataForGeomech()'
#endif


  select case(this%option%iflowmode)
    case (TH_MODE,TH_TS_MODE)
      pres_dof = TH_PRESSURE_DOF
      temp_dof = TH_TEMPERATURE_DOF
    case (MPH_MODE)
      pres_dof = MPH_PRESSURE_DOF
      temp_dof = MPH_TEMPERATURE_DOF
    case(RICHARDS_MODE,RICHARDS_TS_MODE)
      pres_dof = RICHARDS_PRESSURE_DOF
    case default
      this%option%io_buffer = 'PMCSubsurfaceSetAuxDataForGeomech() not ' // &
        'supported for ' // trim(this%option%flowmode)
      call PrintErrMsg(this%option)
  end select

  if (associated(this%sim_aux)) then

    select type (pmc => this)
      class is (pmc_subsurface_type)

        option        => pmc%option
        subsurf_grid  => pmc%realization%discretization%grid
        subsurf_field => pmc%realization%field


        ! Extract pressure, temperature and porosity from subsurface realization
        call VecGetArrayF90(subsurf_field%flow_xx_loc, xx_loc_p,  &
                            ierr);CHKERRQ(ierr)
        call VecGetArrayF90(pmc%sim_aux%subsurf_pres, pres_p,  &
                            ierr);CHKERRQ(ierr)
        call VecGetArrayF90(pmc%sim_aux%subsurf_temp, temp_p,  &
                            ierr);CHKERRQ(ierr)

        do local_id = 1, subsurf_grid%nlmax
          ghosted_id = subsurf_grid%nL2G(local_id)
          pres_p(local_id) = xx_loc_p(option%nflowdof*(ghosted_id - 1) + &
                                      pres_dof)
          if (this%option%iflowmode == RICHARDS_MODE .or. &
              this%option%iflowmode == RICHARDS_TS_MODE) then
            temp_p(local_id) = this%option%reference_temperature
          else
            temp_p(local_id) = xx_loc_p(option%nflowdof*(ghosted_id - 1) + &
                                        temp_dof)
          endif
        enddo

        call VecRestoreArrayF90(subsurf_field%flow_xx_loc, xx_loc_p,  &
                                ierr);CHKERRQ(ierr)
        call VecRestoreArrayF90(pmc%sim_aux%subsurf_pres, pres_p,  &
                                ierr);CHKERRQ(ierr)
        call VecRestoreArrayF90(pmc%sim_aux%subsurf_temp, temp_p,  &
                                ierr);CHKERRQ(ierr)

        if (pmc%timestepper%steps == 0) then
          material_auxvars => pmc%realization%patch%aux%Material%auxvars
          call VecGetArrayF90(pmc%sim_aux%subsurf_por0, sim_por0_p,  &
                              ierr);CHKERRQ(ierr)
          call VecGetArrayF90(pmc%sim_aux%subsurf_perm0, sim_perm0_p,  &
                              ierr);CHKERRQ(ierr)

          ghosted_id = subsurf_grid%nL2G(1)
          
          do local_id = 1, subsurf_grid%nlmax
            ghosted_id = subsurf_grid%nL2G(local_id)
            sim_por0_p(local_id) = material_auxvars(ghosted_id)%porosity ! Set here.  
            sim_perm0_p(local_id) = material_auxvars(ghosted_id)%permeability(perm_xx_index) ! assuming isotropic perm
          enddo
          call VecRestoreArrayF90(pmc%sim_aux%subsurf_por0, sim_por0_p,  &
                                  ierr);CHKERRQ(ierr)
          call VecRestoreArrayF90(pmc%sim_aux%subsurf_perm0, sim_perm0_p,  &
                                  ierr);CHKERRQ(ierr)
        endif
    end select
  endif

end subroutine PMCSubsurfaceSetAuxDataForGeomech

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
