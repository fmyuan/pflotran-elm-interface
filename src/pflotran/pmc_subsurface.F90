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
    procedure, public :: FinalizeRun => PMCSubsurfaceFinalizeRun
    procedure, public :: Destroy => PMCSubsurfaceDestroy
  end type pmc_subsurface_type

  public :: PMCSubsurfaceCreate, &
            PMCSubsurfaceInit, &
            PMCSubsurfaceStrip

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
  use Timestepper_SNES_class
  use Timestepper_TS_class

  implicit none

  class(pmc_subsurface_type) :: this

  type(option_type), pointer :: option

  option => this%option

  if (associated(this%timestepper)) then
    select type(ts => this%timestepper)
      class is(timestepper_TS_type)
        call PMCSubsurfaceSetupSolvers_TS(this)
      class is(timestepper_SNES_type)
        call PMCSubsurfaceSetupSolvers_TimestepperSNES(this)
      class default
        option%io_buffer = &
          'Unknown timestepper found in PMCSubsurfaceSetupSolvers '
        call PrintErrMsg(option)
    end select
  endif

end subroutine PMCSubsurfaceSetupSolvers

! ************************************************************************** !

subroutine PMCSubsurfaceSetupSolvers_TimestepperSNES(this)
  !
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  !
  use Convergence_module
  use Discretization_module
  use Realization_Subsurface_class
  use Option_module
  use General_Aux_module
  use PMC_Base_class
  use PM_Base_Pointer_module
  use PM_Base_class
  use PM_Subsurface_Flow_class
  use PM_General_class
  use PM_WIPP_Flow_class
  use PM_Richards_class
  use PM_TH_class
  use PM_ZFlow_class
  use PM_RT_class
  use PM_NWT_class
  use PM_Waste_Form_class
  use PM_UFD_Decay_class
  use PM_Well_class
  use Solver_module
  use Timestepper_Base_class
  use Timestepper_SNES_class

  implicit none

  class(pmc_subsurface_type) :: this

  type(solver_type), pointer :: solver
  type(option_type), pointer :: option
  type(discretization_type), pointer :: discretization
  class(realization_subsurface_type), pointer :: realization
  SNESLineSearch :: linesearch
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: add_pre_check, check_update, check_post_convergence
  PetscInt :: itranmode
  PetscInt :: transport_coupling
  PetscErrorCode :: ierr

  option => this%option
  solver => this%timestepper%solver

  check_update = PETSC_FALSE
  itranmode = NULL_MODE
  transport_coupling = UNINITIALIZED_INTEGER

  select type(pm => this%pm_ptr%pm)
  ! ----- subsurface flow
    class is(pm_subsurface_flow_type)
      call PrintMsg(option,"  Beginning setup of FLOW SNES ")
      call SolverCreateSNES(solver,option%mycomm,'flow_',option)
      if (solver%M_mat_type == MATAIJ .and. &
          option%iflowmode /= RICHARDS_MODE) then
        option%io_buffer = 'AIJ matrix not supported for current &
          &mode: '// option%flowmode
        call PrintErrMsg(option)
      endif
      write(option%io_buffer,'(" number of dofs = ",i3,", number of &
                &phases = ",i3,i2)') option%nflowdof,option%nphase
      call PrintMsg(option)
      select case(option%iflowmode)
        case(MPH_MODE)
          string = " mode = MPH: p, T, s/X"
        case(TH_MODE)
          string = " mode = TH: p, T"
        case(RICHARDS_MODE)
          string = " mode = Richards: p"
        case(PNF_MODE)
          string = " mode = PNF: h"
        case(ZFLOW_MODE)
          string = " mode = ZFlow: p"
        case(G_MODE)
          string = " mode = General: p, sg/X, T"
        case(H_MODE)
          string = " mode = Hydrate: p, sg/sh/si/X, T"
        case(WF_MODE)
          string = " mode = WIPP Flow: p, sg"
        case default
          string = "mode unknown"
      end select
      call PrintMsg(option,string)

! ----- Set up the J and Jpre matrices -----
! 1) If neither J_mat_type or Jpre_mat_type are specified, set to default.
! 2) If only one of J_mat_type and Jpre_mat_type are specified, then default
!    to setting the other to the same value (except for MATMFFD case).
! 3) Once J_mat_type and Jpre_mat_type are set appropriately, then
!    * If J_mat_type == Jpre_mat_type, then set solver%M = solver%Mpre
!    * Otherwise
!      - Create different matrices for each.
!      - Inside Jacobian routines, will need to check for
!        solver%M != solver%Mpre, and populate two matrices if so.

      if (Uninitialized(solver%Mpre_mat_type) .and. &
          Uninitialized(solver%M_mat_type)) then
        ! Matrix types not specified, so set to default.
        solver%Mpre_mat_type = MATBAIJ
        solver%M_mat_type = solver%Mpre_mat_type
      else if (Uninitialized(solver%Mpre_mat_type)) then
        if (solver%M_mat_type == MATMFFD) then
          solver%Mpre_mat_type = MATBAIJ
        else
          solver%Mpre_mat_type = solver%M_mat_type
        endif
      else if (Uninitialized(solver%M_mat_type)) then
        solver%M_mat_type = solver%Mpre_mat_type
      endif

      if (associated(solver%cprstash)) then
        call CPRWorkersCreate(pm, solver, option)
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

      if (solver%use_galerkin_mg) then
        call DiscretizationCreateInterpolation( &
                       pm%realization%discretization,NFLOWDOF, &
                       solver%interpolation, &
                       solver%galerkin_mg_levels_x, &
                       solver%galerkin_mg_levels_y, &
                       solver%galerkin_mg_levels_z, &
                       option)
      endif

      if (solver%M_mat_type == MATMFFD) then
        call MatCreateSNESMF(solver%snes,solver%M,ierr);CHKERRQ(ierr)
      endif

      if (solver%snes_type == SNESNEWTONLS) then
        ! by default turn off line search
        call SNESGetLineSearch(solver%snes,linesearch,ierr);CHKERRQ(ierr)
        call SNESLineSearchSetType(linesearch,SNESLINESEARCHBASIC, &
                                   ierr);CHKERRQ(ierr)
      endif

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
        call PCSetDM(solver%pc,pm%realization%discretization%dm_nflowdof%dm, &
                     ierr);CHKERRQ(ierr)
      endif

      call SNESSetConvergenceTest(solver%snes,PMCheckConvergencePtr, &
                                  this%pm_ptr,PETSC_NULL_FUNCTION, &
                                  ierr);CHKERRQ(ierr)

      if (pm%check_post_convergence) then
        select case(solver%snes_type)
          case(SNESNEWTONTR)
            call SNESNewtonTRSetPostCheck(solver%snes,PMCheckUpdatePostTRPtr, &
                                          this%pm_ptr,ierr);CHKERRQ(ierr)
          case(SNESNEWTONLS)
            call SNESLineSearchSetPostCheck(linesearch,PMCheckUpdatePostPtr, &
                                            this%pm_ptr,ierr);CHKERRQ(ierr)
          case(SNESNEWTONTRDC)
            call SNESNewtonTRDCSetPostCheck(solver%snes, &
                                            PMCheckUpdatePostTRPtr, &
                                            this%pm_ptr,ierr);CHKERRQ(ierr)
        end select
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
        class is(pm_zflow_type)
              add_pre_check = PETSC_TRUE
        class is(pm_th_type)
          if (Initialized(pm%pressure_dampening_factor) .or. &
              Initialized(pm%pressure_change_limit) .or. &
              Initialized(pm%temperature_change_limit)) then
              add_pre_check = PETSC_TRUE
          endif
      end select

      if (add_pre_check) then
        select case(solver%snes_type)
          case(SNESNEWTONTR)
            call SNESNewtonTRSetPreCheck(solver%snes,PMCheckUpdatePreTRPtr, &
                                         this%pm_ptr,ierr);CHKERRQ(ierr)
          case(SNESNEWTONLS)
            call SNESLineSearchSetPreCheck(linesearch,PMCheckUpdatePrePtr, &
                                           this%pm_ptr,ierr);CHKERRQ(ierr)
          case(SNESNEWTONTRDC)
            call SNESNewtonTRDCSetPreCheck(solver%snes, &
                                           PMCheckUpdatePreTRPtr, &
                                           this%pm_ptr,ierr);CHKERRQ(ierr)
        end select
      endif

      call PrintMsg(option,"  Finished setting up FLOW SNES ")
  ! ----- subsurface reactive transport
    class is(pm_rt_type)
      check_post_convergence = pm%check_post_convergence
      itranmode = option%itranmode
      transport_coupling = option%transport%reactive_transport_coupling
      check_update = pm%realization%reaction%check_update
      discretization => pm%realization%discretization
      realization => pm%realization
      option%io_buffer = "mode = Reactive Transport"
      call PrintMsg(option)
  ! ----- nuclear waste transport
    class is(pm_nwt_type)
      check_post_convergence = pm%controls%check_post_convergence
      itranmode = option%itranmode
      transport_coupling = option%transport%nw_transport_coupling
      check_update = pm%controls%check_update
      discretization => pm%realization%discretization
      realization => pm%realization
      option%io_buffer = "mode = Nuclear Waste Transport"
      call PrintMsg(option)
  end select

  if (itranmode == RT_MODE .or. itranmode == NWT_MODE) then
    call PrintMsg(option,"  Beginning setup of TRAN SNES ")
    call SolverCreateSNES(solver,option%mycomm,'tran_',option)

    if (transport_coupling == GLOBAL_IMPLICIT) then
      if (Uninitialized(solver%Mpre_mat_type) .and. &
          Uninitialized(solver%M_mat_type)) then
        ! Matrix types not specified, so set to default.
        solver%Mpre_mat_type = MATBAIJ
        solver%M_mat_type = solver%Mpre_mat_type
      else if (Uninitialized(solver%Mpre_mat_type)) then
        if (solver%M_mat_type == MATMFFD) then
          solver%Mpre_mat_type = MATBAIJ
        else
          solver%Mpre_mat_type = solver%M_mat_type
        endif
      else if (Uninitialized(solver%M_mat_type)) then
        solver%M_mat_type = solver%Mpre_mat_type
      endif

      call DiscretizationCreateMatrix(discretization, &
                                      NTRANDOF, &
                                      solver%Mpre_mat_type, &
                                      solver%Mpre,option)
    else
      solver%M_mat_type = MATAIJ
      solver%Mpre_mat_type = MATAIJ
      call DiscretizationCreateMatrix(discretization, &
                                      ONEDOF, &
                                      solver%Mpre_mat_type, &
                                      solver%Mpre,option)
    endif

    call MatSetOptionsPrefix(solver%Mpre,"tran_",ierr);CHKERRQ(ierr)

    if (solver%Mpre_mat_type == solver%M_mat_type) then
      solver%M = solver%Mpre
    else
      call DiscretizationCreateMatrix(discretization, &
                                      NTRANDOF, &
                                      solver%M_mat_type, &
                                      solver%M, &
                                      option)

      call MatSetOptionsPrefix(solver%M,"tran_",ierr);CHKERRQ(ierr)
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

    if (transport_coupling == GLOBAL_IMPLICIT) then

      if (solver%M_mat_type == MATMFFD) then
        call MatCreateSNESMF(solver%snes,solver%M,ierr);CHKERRQ(ierr)
      endif

      ! this could be changed in the future if there is a way to
      ! ensure that the linesearch update does not perturb
      ! concentrations negative.
      if (solver%snes_type == SNESNEWTONLS) then
        call SNESGetLineSearch(solver%snes,linesearch,ierr);CHKERRQ(ierr)
        call SNESLineSearchSetType(linesearch,SNESLINESEARCHBASIC, &
                                   ierr);CHKERRQ(ierr)
      endif

      ! Have PETSc do a SNES_View() at the end of each solve if
      ! verbosity > 0.
      if (option%verbosity >= 2) then
        string = '-tran_snes_view'
        call PetscOptionsInsertString(PETSC_NULL_OPTIONS,string, &
                                      ierr);CHKERRQ(ierr)
      endif

    endif

    if (transport_coupling == GLOBAL_IMPLICIT) then
      call SNESSetConvergenceTest(solver%snes,PMCheckConvergencePtr, &
                                  this%pm_ptr,PETSC_NULL_FUNCTION, &
                                  ierr);CHKERRQ(ierr)
    endif
    if (this%pm_ptr%pm%print_EKG .or. option%use_sc .or. &
        check_post_convergence) then
      select case(solver%snes_type)
        case(SNESNEWTONLS)
          call SNESLineSearchSetPostCheck(linesearch,PMCheckUpdatePostPtr, &
                                          this%pm_ptr,ierr);CHKERRQ(ierr)
        case(SNESNEWTONTRDC)
          call SNESNewtonTRDCSetPostCheck(solver%snes, &
                                          PMCheckUpdatePostTRPtr, &
                                          this%pm_ptr,ierr);CHKERRQ(ierr)
      end select
      if (this%pm_ptr%pm%print_EKG) then
        check_post_convergence = PETSC_TRUE
      endif
    endif
    if (check_update) then
      select case(solver%snes_type)
        case(SNESNEWTONLS)
          call SNESLineSearchSetPreCheck(linesearch,PMCheckUpdatePrePtr, &
                                         this%pm_ptr,ierr);CHKERRQ(ierr)
        case(SNESNEWTONTRDC)
          call SNESNewtonTRDCSetPreCheck(solver%snes, &
                                         PMCheckUpdatePreTRPtr, &
                                         this%pm_ptr,ierr);CHKERRQ(ierr)
      end select
    endif
    call PrintMsg(option,"  Finished setting up TRAN SNES ")

  endif


  call SNESSetFunction(solver%snes,this%pm_ptr%pm%residual_vec,PMResidualPtr, &
                       this%pm_ptr,ierr);CHKERRQ(ierr)
  call SNESSetJacobian(solver%snes,solver%M,solver%Mpre,PMJacobianPtr, &
                       this%pm_ptr,ierr);CHKERRQ(ierr)
  call SolverSetSNESOptions(solver,option)

end subroutine PMCSubsurfaceSetupSolvers_TimestepperSNES

! ************************************************************************** !

subroutine PMCSubsurfaceSetupSolvers_TS(this)
  !
  ! Author: Gautam Bisht
  ! Date: 06/20/2018
  !
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

  call TSSetEquationType(solver%ts,TS_EQ_IMPLICIT,ierr);CHKERRQ(ierr)

  call TSSetType(solver%ts,TSBEULER,ierr);CHKERRQ(ierr)

  select type(pm => this%pm_ptr%pm)

    class is (pm_subsurface_flow_type)
      call PrintMsg(option,"  Beginning setup of FLOW SNES ")

      select case(option%iflowmode)
        case(RICHARDS_TS_MODE,TH_TS_MODE)
        case default
          option%io_buffer = 'Timestepper TS unsupported for mode: '// option%flowmode
          call PrintErrMsg(option)
        end select

        write(option%io_buffer,'(" number of dofs = ",i3,", number of &
                  &phases = ",i3,i2)') option%nflowdof,option%nphase
        call PrintMsg(option)
        select case(option%iflowmode)
          case(RICHARDS_TS_MODE)
            option%io_buffer = " mode = Richards: p"
          case(TH_TS_MODE)
            option%io_buffer = " mode = TH: p, T"
        end select
        call PrintMsg(option)

        call TSSetOptionsPrefix(solver%ts,"flow_",ierr);CHKERRQ(ierr)
        call TSSetFromOptions(solver%ts,ierr);CHKERRQ(ierr)

        call SolverCheckCommandLine(solver)

        solver%Mpre_mat_type = MATBAIJ

        call DiscretizationCreateMatrix(pm%realization%discretization, &
                                        NFLOWDOF, &
                                        solver%Mpre_mat_type, &
                                        solver%Mpre, &
                                        option)

        call MatSetOptionsPrefix(solver%Mpre,"flow_",ierr);CHKERRQ(ierr)

        if (solver%M_mat_type /= MATMFFD) then
          solver%M = solver%Mpre
        endif

      call TSSetIFunction(solver%ts,this%pm_ptr%pm%residual_vec, &
                          PMIFunctionPtr,this%pm_ptr,ierr);CHKERRQ(ierr)

      call TSSetIJacobian(solver%ts,solver%M,solver%Mpre,PMIJacobianPtr, &
                          this%pm_ptr,ierr);CHKERRQ(ierr)

      call TSGetSNES(solver%ts,solver%snes,ierr);CHKERRQ(ierr)
      call SNESGetKSP(solver%snes,solver%ksp,ierr);CHKERRQ(ierr)
      call KSPGetPC(solver%ksp,solver%pc,ierr);CHKERRQ(ierr)

      call SNESSetConvergenceTest(solver%snes,PMCheckConvergencePtr, &
                                  this%pm_ptr,PETSC_NULL_FUNCTION, &
                                  ierr);CHKERRQ(ierr)

      call PrintMsg(option,"  Finished setting up FLOW SNES ")
      call SolverSetSNESOptions(solver,option)

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

  !TODO(geh): create a class for get/set and add it to a dynamic linked list
  !           based on the process models involved. there is no need to
  !           extend this class just to override this function as there
  !           may be more than one PM for which to apply get/set. the
  !           proposed accomplishes the task.
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

  if (this%option%ngeomechdof > 0) call PMCSubsurfaceSetAuxDataForGeomech(this)

end subroutine PMCSubsurfaceSetAuxData

! ************************************************************************** !

subroutine PMCSubsurfaceGetAuxDataFromGeomech(this)
  !
  ! This routine updates subsurface data from geomechanics process model.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/04/14

  use Discretization_module, only : DiscretizationLocalToGlobal
  use Field_module
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use PFLOTRAN_Constants_module
  use Material_Aux_module
  use Material_module
  use Variables_module, only : POROSITY

  implicit none

  class (pmc_subsurface_type) :: this

  type(grid_type), pointer :: subsurf_grid
  type(option_type), pointer :: option
  type(field_type), pointer :: subsurf_field

  PetscScalar, pointer :: sim_por_p(:)
  PetscScalar, pointer :: work_loc_p(:)
  type(material_auxvar_type), pointer :: subsurf_material_auxvars(:)

  PetscInt :: local_id
  PetscInt :: ghosted_id

  PetscErrorCode :: ierr

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

          call VecGetArrayF90(pmc%sim_aux%subsurf_por,sim_por_p, &
                              ierr);CHKERRQ(ierr)
          call VecGetArrayF90(subsurf_field%work_loc,work_loc_p, &
                              ierr);CHKERRQ(ierr)

          do local_id = 1, subsurf_grid%nlmax
            ghosted_id = subsurf_grid%nL2G(local_id)
            work_loc_p(ghosted_id) = sim_por_p(local_id)
          enddo

          call VecRestoreArrayF90(pmc%sim_aux%subsurf_por,sim_por_p, &
                                  ierr);CHKERRQ(ierr)
          call VecRestoreArrayF90(subsurf_field%work_loc,work_loc_p, &
                                  ierr);CHKERRQ(ierr)

          call DiscretizationLocalToGlobal(pmc%realization%discretization, &
                                          subsurf_field%work_loc, &
                                          subsurf_field%porosity_geomech_store, &
                                          ONEDOF)
#ifdef GEOMECH_DEBUG
          call PetscViewerASCIIOpen(pmc%realization%option%mycomm, &
                                    'porosity_geomech_store.out',viewer, &
                                    ierr);CHKERRQ(ierr)
          call VecView(subsurf_field%porosity_geomech_store,viewer, &
                       ierr);CHKERRQ(ierr)
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

  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Material_Aux_module
  use ZFlow_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  class (pmc_subsurface_type) :: this

  type(grid_type), pointer :: subsurf_grid
  type(option_type), pointer :: option
  type(field_type), pointer :: subsurf_field

  PetscScalar, pointer :: xx_loc_p(:)
  PetscScalar, pointer :: pres_p(:)
  PetscScalar, pointer :: temp_p(:)
  PetscScalar, pointer :: sim_por0_p(:)
  PetscScalar, pointer :: sim_perm0_p(:) !DANNY - added this 11/7/16

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: pres_dof
  PetscInt :: temp_dof

  type(material_auxvar_type), pointer :: material_auxvars(:)

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
    case(ZFLOW_MODE)
      if (zflow_liq_flow_eq > 0) then
        pres_dof = zflow_liq_flow_eq
      else
        this%option%io_buffer = 'Geomechanics cannot be run without water &
          &mass conservation in ZFLOW.'
        call PrintErrMsg(this%option)
      endif
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
        call VecGetArrayF90(subsurf_field%flow_xx_loc,xx_loc_p, &
                            ierr);CHKERRQ(ierr)
        call VecGetArrayF90(pmc%sim_aux%subsurf_pres,pres_p, &
                            ierr);CHKERRQ(ierr)
        call VecGetArrayF90(pmc%sim_aux%subsurf_temp,temp_p, &
                            ierr);CHKERRQ(ierr)

        do local_id = 1, subsurf_grid%nlmax
          ghosted_id = subsurf_grid%nL2G(local_id)
          pres_p(local_id) = xx_loc_p(option%nflowdof*(ghosted_id - 1) + &
                                      pres_dof)
          if (this%option%iflowmode == RICHARDS_MODE .or. &
              this%option%iflowmode == RICHARDS_TS_MODE .or. &
              this%option%iflowmode == ZFLOW_MODE) then
            temp_p(local_id) = this%option%flow%reference_temperature
          else
            temp_p(local_id) = xx_loc_p(option%nflowdof*(ghosted_id - 1) + &
                                        temp_dof)
          endif
        enddo

        call VecRestoreArrayF90(subsurf_field%flow_xx_loc,xx_loc_p, &
                                ierr);CHKERRQ(ierr)
        call VecRestoreArrayF90(pmc%sim_aux%subsurf_pres,pres_p, &
                                ierr);CHKERRQ(ierr)
        call VecRestoreArrayF90(pmc%sim_aux%subsurf_temp,temp_p, &
                                ierr);CHKERRQ(ierr)

        if (pmc%timestepper%steps == 0) then
          material_auxvars => pmc%realization%patch%aux%Material%auxvars
          call VecGetArrayF90(pmc%sim_aux%subsurf_por0,sim_por0_p, &
                              ierr);CHKERRQ(ierr)
          call VecGetArrayF90(pmc%sim_aux%subsurf_perm0,sim_perm0_p, &
                              ierr);CHKERRQ(ierr)

          ghosted_id = subsurf_grid%nL2G(1)

          do local_id = 1, subsurf_grid%nlmax
            ghosted_id = subsurf_grid%nL2G(local_id)
            sim_por0_p(local_id) = &
              material_auxvars(ghosted_id)%porosity ! Set here.
            ! assuming isotropic perm
            sim_perm0_p(local_id) = &
              material_auxvars(ghosted_id)%permeability(perm_xx_index)
          enddo
          call VecRestoreArrayF90(pmc%sim_aux%subsurf_por0,sim_por0_p, &
                                  ierr);CHKERRQ(ierr)
          call VecRestoreArrayF90(pmc%sim_aux%subsurf_perm0,sim_perm0_p, &
                                  ierr);CHKERRQ(ierr)
        endif
    end select
  endif

end subroutine PMCSubsurfaceSetAuxDataForGeomech

! ************************************************************************** !

recursive subroutine PMCSubsurfaceFinalizeRun(this)
  !
  ! Finalizes the time stepping
  !
  ! Author: Glenn Hammond
  ! Date: 03/18/13
  !

  use Option_module

  ! needed for Intel compilers
  use PMC_Base_class, only : PMCBaseFinalizeRun

  implicit none

  class(pmc_subsurface_type) :: this

#ifdef DEBUG
  call PrintMsg(this%option,'PMCSubsurface%FinalizeRun()')
#endif

  nullify(this%realization)
  call PMCBaseFinalizeRun(this)

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
  call DiscretizationCreateMatrix(pm%realization%discretization, &
                                  ONEDOF, &
                                  cpr_ap_mat_type, &
                                  solver%cprstash%Ap, &
                                  option)
  call DiscretizationCreateMatrix(pm%realization%discretization, &
                                  ONEDOF, &
                                  cpr_ap_mat_type, &
                                  solver%cprstash%As, &
                                  option)
  call DiscretizationCreateVector(pm%realization%discretization, &
                                  NFLOWDOF, solver%cprstash%T1r, &
                                  GLOBAL, option)
  call DiscretizationCreateVector(pm%realization%discretization, &
                                  NFLOWDOF, solver%cprstash%T3r, &
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
  call DiscretizationCreateVector(pm%realization%discretization, &
                                  NFLOWDOF, solver%cprstash%factors3vec, &
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
