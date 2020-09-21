module PM_TH_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class
!geh: using TH_module here fails with gfortran (internal compiler error)
!  use TH_module
  use Realization_Subsurface_class
  use Communicator_Base_module
  use Option_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  PetscInt, parameter :: ABS_UPDATE_INDEX = 1
  PetscInt, parameter :: REL_UPDATE_INDEX = 2
  PetscInt, parameter :: RESIDUAL_INDEX = 3
  PetscInt, parameter :: SCALED_RESIDUAL_INDEX = 4
  PetscInt, parameter :: MAX_INDEX = SCALED_RESIDUAL_INDEX

  type, public, extends(pm_subsurface_flow_type) :: pm_th_type
    class(communicator_type), pointer :: commN
    PetscBool :: converged_flag(2,MAX_INDEX)
    PetscInt :: converged_cell(2,MAX_INDEX)
    PetscReal :: converged_real(2,MAX_INDEX)
    PetscReal :: residual_abs_inf_tol(2)
    PetscReal :: residual_scaled_inf_tol(2)
    PetscReal :: abs_update_inf_tol(2)
    PetscReal :: rel_update_inf_tol(2)
  contains
    procedure, public :: Setup => PMTHSetup
    procedure, public :: ReadSimulationOptionsBlock => PMTHReadSimOptionsBlock
    procedure, public :: ReadNewtonBlock => PMTHReadNewtonSelectCase
    procedure, public :: InitializeTimestep => PMTHInitializeTimestep
    procedure, public :: Residual => PMTHResidual
    procedure, public :: Jacobian => PMTHJacobian
    procedure, public :: UpdateTimestep => PMTHUpdateTimestep
    procedure, public :: PreSolve => PMTHPreSolve
    procedure, public :: PostSolve => PMTHPostSolve
    procedure, public :: CheckUpdatePre => PMTHCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMTHCheckUpdatePost
    procedure, public :: CheckConvergence => PMTHCheckConvergence
    procedure, public :: TimeCut => PMTHTimeCut
    procedure, public :: UpdateSolution => PMTHUpdateSolution
    procedure, public :: UpdateAuxVars => PMTHUpdateAuxVars
    procedure, public :: MaxChange => PMTHMaxChange
    procedure, public :: ComputeMassBalance => PMTHComputeMassBalance
    procedure, public :: InputRecord => PMTHInputRecord
    procedure, public :: Destroy => PMTHDestroy
  end type pm_th_type
  
  public :: PMTHCreate, &
            PMTHDestroy, &
            PMTHCheckConvergence
  
contains

! ************************************************************************** !

function PMTHCreate()
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  implicit none
  
  class(pm_th_type), pointer :: PMTHCreate

  class(pm_th_type), pointer :: this

  PetscReal, parameter :: pres_abs_inf_tol = 1.d0
  PetscReal, parameter :: temp_abs_inf_tol = 1.d-5
  PetscReal, parameter :: abs_update_inf_tol(2) = &
                            [pres_abs_inf_tol,temp_abs_inf_tol]
  PetscReal, parameter :: pres_rel_inf_tol = 1.d-5
  PetscReal, parameter :: temp_rel_inf_tol = 1.d-5
  PetscReal, parameter :: rel_update_inf_tol(2) = &
                            [pres_rel_inf_tol,temp_rel_inf_tol]
  PetscReal, parameter :: residual_abs_inf_tol(2) = 1.d-5
  PetscReal, parameter :: residual_scaled_inf_tol(2) = 1.d-5
  
#ifdef PM_TH_DEBUG
  print *, 'PMTHCreate()'
#endif  

  allocate(this)

  nullify(this%commN)

  call PMSubsurfaceFlowInit(this)
  this%name = 'TH Flow'
  this%header = 'TH FLOW'

  this%residual_abs_inf_tol = residual_abs_inf_tol
  this%residual_scaled_inf_tol = residual_scaled_inf_tol
  this%abs_update_inf_tol = abs_update_inf_tol
  this%rel_update_inf_tol = rel_update_inf_tol

  PMTHCreate => this
  
end function PMTHCreate

! ************************************************************************** !

subroutine PMTHReadSimOptionsBlock(this,input)
  ! 
  ! Reads input file parameters associated with the TH process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/29/15

  use Input_Aux_module
  use String_module
  use Utility_module
  use EOS_Water_module  
  use Option_module
  use TH_Aux_module
 
  implicit none
  
  class(pm_th_type) :: this
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option
  PetscBool :: found
  PetscReal :: tempreal

  option => this%option

  error_string = 'TH Options'
  
  input%ierr = 0
  call InputPushBlock(input,option)
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    found = PETSC_FALSE
    call PMSubsurfFlowReadSimOptionsSC(this,input,keyword,found, &
                                       error_string,option)
    if (found) cycle
    
    select case(trim(keyword))
      case('FREEZING')
        option%freezing = PETSC_TRUE
        option%io_buffer = ' TH: using FREEZING submode!'
        call PrintMsg(option)
        ! Override the default setting for TH-mode with freezing
        call EOSWaterSetDensity('PAINTER')
        call EOSWaterSetEnthalpy('PAINTER')
      case('ICE_MODEL')
        call InputReadCard(input,option,keyword,PETSC_FALSE)
        call StringToUpper(keyword)
        select case (trim(keyword))
          case ('PAINTER_EXPLICIT')
            th_ice_model = PAINTER_EXPLICIT
          case ('PAINTER_KARRA_IMPLICIT')
            th_ice_model = PAINTER_KARRA_IMPLICIT
          case ('PAINTER_KARRA_EXPLICIT')
            th_ice_model = PAINTER_KARRA_EXPLICIT
          case ('PAINTER_KARRA_EXPLICIT_NOCRYO')
            th_ice_model = PAINTER_KARRA_EXPLICIT_NOCRYO
          case ('DALL_AMICO')
            th_ice_model = DALL_AMICO
          case default
            option%io_buffer = 'Cannot identify the specificed ice model. &
             &Specify PAINTER_EXPLICIT or PAINTER_KARRA_IMPLICIT &
             &or PAINTER_KARRA_EXPLICIT or PAINTER_KARRA_EXPLICIT_NOCRYO &
             &or DALL_AMICO.'
            call PrintErrMsg(option)
          end select
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)
  
end subroutine PMTHReadSimOptionsBlock

! ************************************************************************** !

subroutine PMTHReadNewtonSelectCase(this,input,keyword,found, &
                                    error_string,option)
  ! 
  ! Reads input file parameters associated with the TH process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/16/20

  use Input_Aux_module
  use String_module
  use Utility_module
  use Option_module
  use TH_Aux_module
 
  implicit none
  
  class(pm_th_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option

  PetscBool :: found
  PetscReal :: tempreal
  PetscInt :: lid, eid

  option => this%option

  lid = option%liquid_phase
  eid = option%energy_id
  
  error_string = 'TH Newton Solver'
  
  found = PETSC_FALSE
  call PMSubsurfaceFlowReadNewtonSelectCase(this,input,keyword,found, &
                                            error_string,option)
  if (found) return
    
  found = PETSC_TRUE
  select case(trim(keyword))

    ! Tolerances

    ! All Residual
    case('RESIDUAL_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%residual_abs_inf_tol(:) = tempreal
      this%residual_scaled_inf_tol(:) = tempreal

    ! Absolute Residual
    case('RESIDUAL_ABS_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%residual_abs_inf_tol(:) = tempreal
    case('LIQUID_RESIDUAL_ABS_INF_TOL')
      call InputReadDouble(input,option,this%residual_abs_inf_tol(lid))
      call InputErrorMsg(input,option,keyword,error_string)
    case('ENERGY_RESIDUAL_ABS_INF_TOL')
      call InputReadDouble(input,option,this%residual_abs_inf_tol(eid))
      call InputErrorMsg(input,option,keyword,error_string)

    ! Scaled Residual
    case('RESIDUAL_SCALED_INF_TOL','ITOL_SCALED_RESIDUAL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%residual_scaled_inf_tol(:) = tempreal
    case('LIQUID_RESIDUAL_SCALED_INF_TOL')
      call InputReadDouble(input,option,this%residual_scaled_inf_tol(lid))
      call InputErrorMsg(input,option,keyword,error_string)
    case('ENERGY_RESIDUAL_SCALED_INF_TOL')
      call InputReadDouble(input,option,this%residual_scaled_inf_tol(eid))
      call InputErrorMsg(input,option,keyword,error_string)

    ! All Updates
    case('UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%abs_update_inf_tol(:) = tempreal
      this%rel_update_inf_tol(:) = tempreal

    ! Absolute Updates
    case('ABS_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%abs_update_inf_tol(:) = tempreal
    case('PRES_ABS_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%abs_update_inf_tol(1) = tempreal
    case('TEMP_ABS_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%abs_update_inf_tol(2) = tempreal

    ! Relative Updates
    case('REL_UPDATE_INF_TOL','ITOL_RELATIVE_UPDATE')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%rel_update_inf_tol(:) = tempreal
    case('PRES_REL_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%rel_update_inf_tol(1) = tempreal
    case('TEMP_REL_UPDATE_INF_TOL')
      call InputReadDouble(input,option,tempreal)
      call InputErrorMsg(input,option,keyword,error_string)
      this%rel_update_inf_tol(2) = tempreal

    case default
      found = PETSC_FALSE

  end select
  
end subroutine PMTHReadNewtonSelectCase

! ************************************************************************** !

subroutine PMTHSetup(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use Discretization_module
  use Communicator_Structured_class
  use Communicator_Unstructured_class
  use Grid_module 

  implicit none
  
  class(pm_th_type) :: this

  call PMSubsurfaceFlowSetup(this)
  
  ! set up communicator
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID)
      this%commN => StructuredCommunicatorCreate()
    case(UNSTRUCTURED_GRID)
      this%commN => UnstructuredCommunicatorCreate()
  end select
  call this%commN%SetDM(this%realization%discretization%dm_nflowdof)

end subroutine PMTHSetup

! ************************************************************************** !

subroutine PMTHInitializeTimestep(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THInitializeTimestep
  use Option_module

  implicit none
  
  class(pm_th_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)

  ! update porosity

  call THInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)
  
end subroutine PMTHInitializeTimestep

! ************************************************************************** !

subroutine PMTHPreSolve(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use Global_module

  implicit none
  
  class(pm_th_type) :: this

  call PMSubsurfaceFlowPreSolve(this)

end subroutine PMTHPreSolve

! ************************************************************************** !

subroutine PMTHPostSolve(this)
  ! 
  ! Date: 03/14/13
  ! 

  use Global_module

  implicit none
  
  class(pm_th_type) :: this
  
end subroutine PMTHPostSolve

! ************************************************************************** !

subroutine PMTHUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                              num_newton_iterations,tfac, &
                              time_step_max_growth_factor)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL

  implicit none
  
  class(pm_th_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor
  
  PetscReal :: fac
  PetscReal :: ut
  PetscReal :: up
  PetscReal :: utmp
  PetscReal :: dtt
  PetscReal :: dt_u
  PetscReal :: dt_tfac
  PetscInt :: ifac
  
#ifdef PM_TH_DEBUG
  call PrintMsg(this%option,'PMTH%UpdateTimestep()')
#endif
  
  if (iacceleration > 0) then
    fac = 0.5d0
    if (num_newton_iterations >= iacceleration) then
      fac = 0.33d0
      ut = 0.d0
    else
      up = this%pressure_change_governor/(this%max_pressure_change+0.1)
      utmp = this%temperature_change_governor/ &
             (this%max_temperature_change+1.d-5)
      ut = min(up,utmp)
    endif
    dtt = fac * dt * (1.d0 + ut)
  else
    ifac = max(min(num_newton_iterations,size(tfac)),1)
    dt_tfac = tfac(ifac) * dt

    fac = 0.5d0
    up = this%pressure_change_governor/(this%max_pressure_change+0.1)
    utmp = this%temperature_change_governor/ &
           (this%max_temperature_change+1.d-5)
    ut = min(up,utmp)
    dt_u = fac * dt * (1.d0 + ut)

    dtt = min(dt_tfac,dt_u)
  endif

  dtt = min(time_step_max_growth_factor*dt,dtt)
  if (dtt > dt_max) dtt = dt_max
  ! geh: There used to be code here that cut the time step if it is too
  !      large relative to the simulation time.  This has been removed.
  dtt = max(dtt,dt_min)
  dt = dtt

  call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt,dt_max)
  
end subroutine PMTHUpdateTimestep

! ************************************************************************** !

subroutine PMTHResidual(this,snes,xx,r,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THResidual

  implicit none
  
  class(pm_th_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  call THResidual(snes,xx,r,this%realization,ierr)

end subroutine PMTHResidual

! ************************************************************************** !

subroutine PMTHJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THJacobian

  implicit none
  
  class(pm_th_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  call THJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMTHJacobian

! ************************************************************************** !

subroutine PMTHCheckUpdatePre(this,snes,X,dX,changed,ierr)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Option_module
  use Saturation_Function_module
  use Patch_module
  use TH_Aux_module
  use Global_Aux_module

  implicit none
  
  class(pm_th_type) :: this
  SNES :: snes
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: X_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscReal, pointer :: r_p(:)
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(TH_auxvar_type), pointer :: TH_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  PetscInt :: local_id, ghosted_id
  PetscReal :: P0, P1, P_R, delP, delP_old
  PetscReal :: scale, press_limit, temp_limit
  PetscInt :: iend, istart
  
  patch => this%realization%patch
  grid => patch%grid
  option => this%realization%option
  field => this%realization%field
  TH_auxvars => patch%aux%TH%auxvars
  global_auxvars => patch%aux%Global%auxvars

  if (Initialized(this%pressure_change_limit)) then

    call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(X,X_p,ierr);CHKERRQ(ierr)

    press_limit = dabs(this%pressure_change_limit)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      P0 = X_p(istart)
      delP = dX_p(istart)
      if (press_limit < dabs(delP)) then
        write(option%io_buffer,'("dP_trunc:",1i7,2es15.7)') &         
          grid%nG2A(grid%nL2G(local_id)),press_limit,dabs(delP)
        call PrintMsgAnyRank(option)
      endif
      delP = sign(min(dabs(delP),press_limit),delP)
      dX_p(istart) = delP
    enddo
    
    call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(X,X_p,ierr);CHKERRQ(ierr)

  endif
  
  if (dabs(this%temperature_change_limit) > 0.d0) then
      
    call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(X,X_p,ierr);CHKERRQ(ierr)

    temp_limit = dabs(this%temperature_change_limit)
    do local_id = 1, grid%nlmax
      ghosted_id = grid%nL2G(local_id)
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      P0 = X_p(iend)
      delP = dX_p(iend)
      if (abs(delP) > abs(temp_limit)) then
        write(option%io_buffer,'("dT_trunc:",1i7,2es15.7)') &
          grid%nG2A(grid%nL2G(local_id)),temp_limit,dabs(delP)
        call PrintMsgAnyRank(option)
      endif
      delP = sign(min(dabs(delP),temp_limit),delP)
      dX_p(iend) = delP
    enddo
    
    call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(X,X_p,ierr);CHKERRQ(ierr)
    
  endif


  if (Initialized(this%pressure_dampening_factor)) then
    ! P^p+1 = P^p - dP^p
    P_R = option%reference_pressure
    scale = this%pressure_dampening_factor

    call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(X,X_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    do local_id = 1, grid%nlmax
      iend = local_id*option%nflowdof
      istart = iend-option%nflowdof+1
      delP = dX_p(istart)
      P0 = X_p(istart)
      P1 = P0 - delP
      if (P0 < P_R .and. P1 > P_R) then
        write(option%io_buffer,'("U -> S:",1i7,2f12.1)') &
          grid%nG2A(grid%nL2G(local_id)),P0,P1 
        call PrintMsgAnyRank(option)
#if 0
        ghosted_id = grid%nL2G(local_id)
        call RichardsPrintAuxVars(rich_auxvars(ghosted_id), &
                                  global_auxvars(ghosted_id),ghosted_id)
        write(option%io_buffer,'("Residual:",es15.7)') r_p(istart)
        call PrintMsgAnyRank(option)
#endif
      else if (P1 < P_R .and. P0 > P_R) then
        write(option%io_buffer,'("S -> U:",1i7,2f12.1)') &
          grid%nG2A(grid%nL2G(local_id)),P0,P1
        call PrintMsgAnyRank(option)
#if 0
        ghosted_id = grid%nL2G(local_id)
        call RichardsPrintAuxVars(rich_auxvars(ghosted_id), &
                                  global_auxvars(ghosted_id),ghosted_id)
        write(option%io_buffer,'("Residual:",es15.7)') r_p(istart)
        call PrintMsgAnyRank(option)
#endif
      endif
      ! transition from unsaturated to saturated
      if (P0 < P_R .and. P1 > P_R) then
        dX_p(istart) = scale*delP
      endif
    enddo
    call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(X,X_p,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
  endif

end subroutine PMTHCheckUpdatePre

! ************************************************************************** !

subroutine PMTHCheckUpdatePost(this,snes,X0,dX,X1,dX_changed, &
                                  X1_changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/21/18
  ! 

  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Option_module
  use Patch_module

  implicit none
  
  class(pm_th_type) :: this
  SNES :: snes
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: X0_p(:)
  PetscReal, pointer :: dX_p(:)
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  PetscInt :: local_id, ghosted_id, natural_id
  PetscInt :: offset, ival, idof
  PetscReal :: dX_, X0_, dX_X0
  PetscBool :: converged_abs_update_flag(2)
  PetscBool :: converged_rel_update_flag(2)
  PetscInt :: converged_abs_update_cell(2)
  PetscInt :: converged_rel_update_cell(2)
  PetscReal :: converged_abs_update_real(2)
  PetscReal :: converged_rel_update_real(2)
  PetscBool :: converged_absolute
  PetscBool :: converged_relative

  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch

  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE

  call VecGetArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)
  converged_abs_update_flag = PETSC_TRUE
  converged_rel_update_flag = PETSC_TRUE
  converged_abs_update_cell = ZERO_INTEGER
  converged_rel_update_cell = ZERO_INTEGER
  converged_abs_update_real = 0.d0
  converged_rel_update_real = 0.d0
  do local_id = 1, grid%nlmax
    offset = (local_id-1)*option%nflowdof
    ghosted_id = grid%nL2G(local_id)
    natural_id = grid%nG2A(ghosted_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    do idof = 1, option%nflowdof
      ival = offset+idof
      ! infinity norms on update
      converged_absolute = PETSC_TRUE
      converged_relative = PETSC_TRUE
      dX_ = dabs(dX_p(ival))
      X0_ = dabs(X0_p(ival))
      X0_ = max(X0_,1.d-40)
      dX_X0 = dX_/X0_
      if (dX_ > this%abs_update_inf_tol(idof)) then
        converged_absolute = PETSC_FALSE
      endif
      if (converged_abs_update_real(idof) < dX_) then
        converged_abs_update_real(idof) = dX_
        converged_abs_update_cell(idof) = natural_id
      endif
      if (dX_X0 > this%rel_update_inf_tol(idof)) then
        converged_relative = PETSC_FALSE
      endif
      if (converged_rel_update_real(idof) < dX_X0) then
        converged_rel_update_real(idof) = dX_X0
        converged_rel_update_cell(idof) = natural_id
      endif
      ! only enter this condition if both are not converged
      if (.not.(converged_absolute .or. converged_relative)) then
        if (.not.converged_absolute) then
          converged_abs_update_flag(idof) = PETSC_FALSE
        endif
        if (.not.converged_relative) then
          converged_rel_update_flag(idof) = PETSC_FALSE
        endif
      endif
    enddo
  enddo
  call VecRestoreArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)

  this%converged_flag(:,ABS_UPDATE_INDEX) = converged_abs_update_flag(:)
  this%converged_flag(:,REL_UPDATE_INDEX) = converged_rel_update_flag(:)
  this%converged_real(:,ABS_UPDATE_INDEX) = converged_abs_update_real(:)
  this%converged_real(:,REL_UPDATE_INDEX) = converged_rel_update_real(:)
  this%converged_cell(:,ABS_UPDATE_INDEX) = converged_abs_update_cell(:)
  this%converged_cell(:,REL_UPDATE_INDEX) = converged_rel_update_cell(:)

end subroutine PMTHCheckUpdatePost

! ************************************************************************** !

subroutine PMTHCheckConvergence(this,snes,it,xnorm,unorm,fnorm,reason,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 11/20/18
  !
  use Convergence_module
  use General_Aux_module
  use Global_Aux_module
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  use String_module

  implicit none

  class(pm_th_type) :: this
  SNES :: snes
  PetscInt :: it
  PetscReal :: xnorm
  PetscReal :: unorm
  PetscReal :: fnorm
  SNESConvergedReason :: reason
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum2_p(:)
  PetscInt :: local_id, ghosted_id, natural_id
  PetscInt :: offset, ival, idof, itol
  PetscReal :: R, A, R_A
  PetscReal, parameter :: A_zero = 1.d-15
  PetscBool :: converged_abs_residual_flag(2)
  PetscReal :: converged_abs_residual_real(2)
  PetscInt :: converged_abs_residual_cell(2)
  PetscBool :: converged_scaled_residual_flag(2)
  PetscReal :: converged_scaled_residual_real(2)
  PetscInt :: converged_scaled_residual_cell(2)
  PetscBool :: converged_absolute
  PetscBool :: converged_scaled
  PetscMPIInt :: mpi_int
  character(len=MAXSTRINGLENGTH) :: string
  character(len=17), parameter :: dof_string(2) = &
                                   ['Liquid Pressure  ','Temperature      ']
  character(len=15), parameter :: tol_string(MAX_INDEX) = &
    ['Absolute Update','Relative Update','Residual       ','Scaled Residual']

  patch => this%realization%patch
  option => this%realization%option
  field => this%realization%field
  grid => patch%grid

  if (this%check_post_convergence) then
    call VecGetArrayReadF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    call VecGetArrayReadF90(field%flow_accum2,accum2_p,ierr);CHKERRQ(ierr)
    converged_abs_residual_flag = PETSC_TRUE
    converged_abs_residual_real = 0.d0
    converged_abs_residual_cell = ZERO_INTEGER
    converged_scaled_residual_flag = PETSC_TRUE
    converged_scaled_residual_real = 0.d0
    converged_scaled_residual_cell = ZERO_INTEGER
    do local_id = 1, grid%nlmax
      offset = (local_id-1)*option%nflowdof
      ghosted_id = grid%nL2G(local_id)
      natural_id = grid%nG2A(ghosted_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      do idof = 1, option%nflowdof
        ival = offset+idof
        converged_absolute = PETSC_TRUE
        converged_scaled = PETSC_TRUE
        ! infinity norms on residual
        R = dabs(r_p(ival))
        A = dabs(accum2_p(ival))
        R_A = R/A
        if (R > this%residual_abs_inf_tol(idof)) then
          converged_absolute = PETSC_FALSE
        endif
        ! find max value regardless of convergence
        if (converged_abs_residual_real(idof) < R) then
          converged_abs_residual_real(idof) = R
          converged_abs_residual_cell(idof) = natural_id
        endif
        if (A > A_zero) then
          if (R_A > this%residual_scaled_inf_tol(idof)) then
            converged_scaled = PETSC_FALSE
          endif
          ! find max value regardless of convergence
          if (converged_scaled_residual_real(idof) < R_A) then
            converged_scaled_residual_real(idof) = R_A
            converged_scaled_residual_cell(idof) = natural_id
          endif
        endif
        ! only enter this condition if both are not converged
        if (.not.(converged_absolute .or. converged_scaled)) then
          if (.not.converged_absolute) then
            converged_abs_residual_flag(idof) = PETSC_FALSE
          endif
          if (.not.converged_scaled) then
            converged_scaled_residual_flag(idof) = PETSC_FALSE
          endif
        endif
      enddo
    enddo
    call VecRestoreArrayReadF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
    call VecRestoreArrayReadF90(field%flow_accum2,accum2_p,ierr);CHKERRQ(ierr)
  
    this%converged_flag(:,RESIDUAL_INDEX) = converged_abs_residual_flag(:)
    this%converged_real(:,RESIDUAL_INDEX) = converged_abs_residual_real(:)
    this%converged_cell(:,RESIDUAL_INDEX) = converged_abs_residual_cell(:)
    this%converged_flag(:,SCALED_RESIDUAL_INDEX) = &
                                         converged_scaled_residual_flag(:)
    this%converged_real(:,SCALED_RESIDUAL_INDEX) = &
                                         converged_scaled_residual_real(:)
    this%converged_cell(:,SCALED_RESIDUAL_INDEX) = &
                                         converged_scaled_residual_cell(:)
    mpi_int = 2*MAX_INDEX
    ! do not perform an all reduce on cell id as this info is not printed 
    ! in parallel
    call MPI_Allreduce(MPI_IN_PLACE,this%converged_flag,mpi_int, &
                       MPI_INTEGER,MPI_LAND,option%mycomm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,this%converged_real,mpi_int, &
                       MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)

    option%convergence = CONVERGENCE_CONVERGED
    do itol = 1, MAX_INDEX
      do idof = 1, option%nflowdof
        if (.not.this%converged_flag(idof,itol)) then
          option%convergence = CONVERGENCE_KEEP_ITERATING
          if (this%logging_verbosity > 0) then
            string = '   ' // trim(tol_string(itol)) // ', ' // &
              dof_string(idof)
            if (option%mycommsize == 1) then
              string = trim(string) // ' (' // &
                trim(StringFormatInt(this%converged_cell(idof,itol))) &
                // ')'
            endif
            string = trim(string) // ' : ' // &
              StringFormatDouble(this%converged_real(idof,itol))
            call OptionPrint(string,option)
          endif
        endif
      enddo
    enddo
    if (this%logging_verbosity > 0 .and. it > 0 .and. &
        option%convergence == CONVERGENCE_CONVERGED) then
      string = '   Converged'
      call OptionPrint(string,option)
      write(string,'(4x," R:",2es8.1)') this%converged_real(:,RESIDUAL_INDEX)
      call OptionPrint(string,option)
      write(string,'(4x,"SR:",2es8.1)') &
        this%converged_real(:,SCALED_RESIDUAL_INDEX)
      call OptionPrint(string,option)
      write(string,'(4x,"AU:",2es8.1)') this%converged_real(:,ABS_UPDATE_INDEX)
      call OptionPrint(string,option)
      write(string,'(4x,"RU:",2es8.1)') this%converged_real(:,REL_UPDATE_INDEX)
      call OptionPrint(string,option)
    endif
  endif

  call PMSubsurfaceFlowCheckConvergence(this,snes,it,xnorm,unorm,fnorm, &
                                        reason,ierr)

end subroutine PMTHCheckConvergence

! ************************************************************************** !

subroutine PMTHTimeCut(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THTimeCut

  implicit none
  
  class(pm_th_type) :: this
  
  call PMSubsurfaceFlowTimeCut(this)
  call THTimeCut(this%realization)

end subroutine PMTHTimeCut

! ************************************************************************** !

subroutine PMTHUpdateSolution(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THUpdateSolution, THUpdateSurfaceBC

  implicit none
  
  class(pm_th_type) :: this
  
  call PMSubsurfaceFlowUpdateSolution(this)
  call THUpdateSolution(this%realization)
  if (this%option%surf_flow_on) &
    call THUpdateSurfaceBC(this%realization)

end subroutine PMTHUpdateSolution     

! ************************************************************************** !

subroutine PMTHUpdateAuxVars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use TH_module, only : THUpdateAuxVars
  
  implicit none
  
  class(pm_th_type) :: this

  call THUpdateAuxVars(this%realization)

end subroutine PMTHUpdateAuxVars   

! ************************************************************************** !

subroutine PMTHMaxChange(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THMaxChange
  use Option_module

  implicit none
  
  class(pm_th_type) :: this
  character(len=MAXSTRINGLENGTH) :: string
  
  call THMaxChange(this%realization,this%max_pressure_change, &
                   this%max_temperature_change)
  write(string,'("  --> max chng: dpmx= ",1pe12.4," dtmpmx= ",1pe12.4)') &
      this%max_pressure_change,this%max_temperature_change
  call OptionPrint(string,this%option)

end subroutine PMTHMaxChange

! ************************************************************************** !

subroutine PMTHComputeMassBalance(this,mass_balance_array)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THComputeMassBalance

  implicit none
  
  class(pm_th_type) :: this
  PetscReal :: mass_balance_array(:)
  
  call THComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMTHComputeMassBalance

! ************************************************************************** !

subroutine PMTHInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  ! 
  
  implicit none
  
  class(pm_th_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'thermo-hydro'

end subroutine PMTHInputRecord

! ************************************************************************** !

subroutine PMTHDestroy(this)
  ! 
  ! This routine
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 03/90/13
  ! 

  use TH_module, only : THDestroy

  implicit none
  
  class(pm_th_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  call this%commN%Destroy()

  ! preserve this ordering
  call THDestroy(this%realization%patch)
  call PMSubsurfaceFlowDestroy(this)

end subroutine PMTHDestroy

end module PM_TH_class
