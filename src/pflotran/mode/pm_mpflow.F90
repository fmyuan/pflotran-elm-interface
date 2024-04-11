module PM_MpFlow_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class
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

  type, public, extends(pm_subsurface_flow_type) :: pm_mpflow_type
    class(communicator_type), pointer :: commN
    PetscBool :: converged_flag(2,MAX_INDEX)
    PetscInt  :: converged_cell(2,MAX_INDEX)
    PetscReal :: converged_real(2,MAX_INDEX)
    PetscReal :: residual_abs_inf_tol(2)
    PetscReal :: residual_scaled_inf_tol(2)
    PetscReal :: abs_update_inf_tol(2)
    PetscReal :: rel_update_inf_tol(2)
  contains
    procedure, public :: Setup => PMMpFlowSetup
    procedure, public :: ReadSimulationBlock => PMMpFlowRead
    procedure, public :: InitializeTimestep => PMMpFlowInitializeTimestep
    procedure, public :: Residual => PMMpFlowResidual
    procedure, public :: Jacobian => PMMpFlowJacobian
    procedure, public :: UpdateTimestep => PMMpFlowUpdateTimestep
    procedure, public :: PreSolve => PMMpFlowPreSolve
    procedure, public :: PostSolve => PMMpFlowPostSolve
    procedure, public :: CheckUpdatePre => PMMpFlowCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMMpFlowCheckUpdatePost
    procedure, public :: CheckConvergence => PMMpFlowCheckConvergence
    procedure, public :: TimeCut => PMMpFlowTimeCut
    procedure, public :: UpdateSolution => PMMpFlowUpdateSolution
    procedure, public :: UpdateAuxVars => PMMpFlowUpdateAuxVars
    procedure, public :: MaxChange => PMMpFlowMaxChange
    procedure, public :: ComputeMassBalance => PMMpFlowComputeMassBalance
    procedure, public :: Destroy => PMMpFlowDestroy
  end type pm_mpflow_type
  
  public :: PMMpFlowCreate, &
            PMMpFlowDestroy, &
            PMMpFlowCheckConvergence
  
contains

! ************************************************************************** !

function PMMpFlowCreate()
  ! 
  ! This routine
  ! 
  ! Author: F-M Yuan, ESD/CCSI-ORNL
  ! Date: 12/23/2019
  ! 

  implicit none
  
  class(pm_mpflow_type), pointer :: PMMpFlowCreate

  class(pm_mpflow_type), pointer :: this

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
  print *, 'PMMpFlowCreate()'
#endif  

  allocate(this)

  nullify(this%commN)

  call PMSubsurfaceFlowCreate(this)
  this%name = 'MpFlow'
  this%header = 'MpFlow MODE'

  this%residual_abs_inf_tol = residual_abs_inf_tol
  this%residual_scaled_inf_tol = residual_scaled_inf_tol
  this%abs_update_inf_tol = abs_update_inf_tol
  this%rel_update_inf_tol = rel_update_inf_tol

  PMMpFlowCreate => this
  
end function PMMpFlowCreate

! ************************************************************************** !

subroutine PMMpFlowRead(this,input)
  ! 
  ! Reads input file parameters associated with the MpFlow process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/29/15
  use Input_Util_module
  use String_module
  use Utility_module
  use EOS_Water_module  
  use Option_module
 
  implicit none
  
  class(pm_mpflow_type) :: this
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option
  PetscBool :: found
  PetscInt :: lid, eid
  PetscReal :: tempreal

  option => this%option
  lid = PRESSURE_DOF
  eid = TEMPERATURE_DOF
  
  error_string = 'MPFLOW Options'
  
  input%ierr = 0
  call InputPushBlock(input,option)
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)

    found = PETSC_FALSE
    call PMSubsurfaceFlowReadSelectCase(this,input,word,found, &
                                        error_string,option)
    if (found) cycle
    
    select case(trim(word))
      ! Tolerances

      ! All Residual
      case('RESIDUAL_INF_TOL')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,word,error_string)
        this%residual_abs_inf_tol(:) = tempreal
        this%residual_scaled_inf_tol(:) = tempreal

      ! Absolute Residual
      case('RESIDUAL_ABS_INF_TOL')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,word,error_string)
        this%residual_abs_inf_tol(:) = tempreal
      case('LIQUID_RESIDUAL_ABS_INF_TOL')
        call InputReadDouble(input,option,this%residual_abs_inf_tol(lid))
        call InputErrorMsg(input,option,word,error_string)
      case('ENERGY_RESIDUAL_ABS_INF_TOL')
        call InputReadDouble(input,option,this%residual_abs_inf_tol(eid))
        call InputErrorMsg(input,option,word,error_string)

      ! Scaled Residual
      case('RESIDUAL_SCALED_INF_TOL','ITOL_SCALED_RESIDUAL')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,word,error_string)
        this%residual_scaled_inf_tol(:) = tempreal
      case('LIQUID_RESIDUAL_SCALED_INF_TOL')
        call InputReadDouble(input,option,this%residual_scaled_inf_tol(lid))
        call InputErrorMsg(input,option,word,error_string)
      case('ENERGY_RESIDUAL_SCALED_INF_TOL')
        call InputReadDouble(input,option,this%residual_scaled_inf_tol(eid))
        call InputErrorMsg(input,option,word,error_string)

      ! All Updates
      case('UPDATE_INF_TOL')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,word,error_string)
        this%abs_update_inf_tol(:) = tempreal
        this%rel_update_inf_tol(:) = tempreal

      ! Absolute Updates
      case('ABS_UPDATE_INF_TOL')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,word,error_string)
        this%abs_update_inf_tol(:) = tempreal
      case('PRES_ABS_UPDATE_INF_TOL')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,word,error_string)
        this%abs_update_inf_tol(1) = tempreal
      case('TEMP_ABS_UPDATE_INF_TOL')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,word,error_string)
        this%abs_update_inf_tol(2) = tempreal

      ! Relative Updates
      case('REL_UPDATE_INF_TOL','ITOL_RELATIVE_UPDATE')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,word,error_string)
        this%rel_update_inf_tol(:) = tempreal
      case('PRES_REL_UPDATE_INF_TOL')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,word,error_string)
        this%rel_update_inf_tol(1) = tempreal
      case('TEMP_REL_UPDATE_INF_TOL')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,word,error_string)
        this%rel_update_inf_tol(2) = tempreal

!....................
      case ('ISOTHERMAL')
        option_flow%isothermal = PETSC_TRUE
        option%io_buffer = 'ISOTHERMAL option is ON in MPFLOW mode.'
        call printMsg(option)

!....................
      case ('ONLY_THERMAL','ISOBARIC')
        option_flow%isobaric = PETSC_TRUE
        option%io_buffer = 'ONLY_THERMAL or ISOBARIC is ON in MPFLOW mode.'
        call printMsg(option)

!....................
      case ('ONLY_VERTICAL_FLOW')
        option_flow%only_vertical_flow = PETSC_TRUE
        option%io_buffer = 'ONLY_VERTICAL_FLOW implemented in MPFLOW mode.'
        call printMsg(option)

!....................
      case('ICE_MODEL')
        option%io_buffer = ' MPFLOW: using ICE submode option ON.'
        call printMsg(option)
        ! Override the default setting for MPFLOW-mode with freezing
        call EOSWaterSetDensity('PAINTER')
        call EOSWaterSetEnthalpy('PAINTER')

        call InputReadWord(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case (trim(word))
          case ('PAINTER_EXPLICIT')
            option_flow%ice_model = PAINTER_EXPLICIT
          case ('PAINTER_KARRA_IMPLICIT')
            option_flow%ice_model = PAINTER_KARRA_IMPLICIT
          case ('PAINTER_KARRA_EXPLICIT')
            option_flow%ice_model = PAINTER_KARRA_EXPLICIT
          case ('PAINTER_KARRA_EXPLICIT_NOCRYO')
            option_flow%ice_model = PAINTER_KARRA_EXPLICIT_NOCRYO
          case ('PAINTER_KARRA_EXPLICIT_SMOOTH')
            option_flow%ice_model = PAINTER_KARRA_EXPLICIT_SMOOTH
            call InputReadDouble(input,option,tempreal)
            call InputDefaultMsg(input,option,'freezing-thawing smoothing')
            if(tempreal > 1.d-10) option_flow%frzthw_halfwidth = tempreal
          case ('DALL_AMICO')
            option_flow%ice_model = DALL_AMICO
          case default
            option%io_buffer = 'Cannot identify the specificed ice model.' // &
             'Specify PAINTER_EXPLICIT or PAINTER_KARRA_IMPLICIT' // &
             ' or PAINTER_KARRA_EXPLICIT or PAINTER_KARRA_EXPLICIT_NOCRYO ' // &
             ' or PAINTER_KARRA_EXPLICIT_SMOOTH ' // &
             ' or DALL_AMICO.' // &
             ' MPFLOW MODE with isothermal option is ON'

            option_flow%isothermal = PETSC_TRUE

            call printMsg(option)
          end select
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)
  
end subroutine PMMpFlowRead

! ************************************************************************** !

subroutine PMMpFlowSetup(this)
  ! 
  ! This routine
  ! 
  ! Author: F-M Yuan, ESD/CCSI-ORNL
  ! Date: 12/23/2019
  ! 

  use Discretization_module
  use Communicator_Structured_class
  use Communicator_Unstructured_class
  use Grid_module 

  implicit none
  
  class(pm_mpflow_type) :: this

  call PMSubsurfaceFlowSetup(this)
  
  ! set up communicator
  select case(this%realization%discretization%itype)
    case(STRUCTURED_GRID)
      this%commN => StructuredCommunicatorCreate()
    case(UNSTRUCTURED_GRID)
      this%commN => UnstructuredCommunicatorCreate()
  end select
  call this%commN%SetDM(this%realization%discretization%dm_nflowdof)

end subroutine PMMpFlowSetup

! ************************************************************************** !

subroutine PMMpFlowInitializeTimestep(this)
  ! 
  ! This routine
  ! 
  ! Author: F-M Yuan, ESD/CCSI-ORNL
  ! Date: 12/23/2019
  ! 

  use MpFlow_module, only : MpFlowInitializeTimestep
  use Option_module

  implicit none
  
  class(pm_mpflow_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)

  ! update porosity
  call this%comm1%LocalToLocal(this%realization%field%icap_loc, &
                               this%realization%field%icap_loc)
  call this%comm1%LocalToLocal(this%realization%field%ithrm_loc, &
                               this%realization%field%ithrm_loc)

  call MpFlowInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)
  
end subroutine PMMpFlowInitializeTimestep

! ************************************************************************** !

subroutine PMMpFlowPreSolve(this)
  ! 
  ! This routine
  ! 
  ! Author: F-M Yuan, ESD/CCSI-ORNL
  ! Date: 12/23/2019
  ! 

  use Global_module

  implicit none
  
  class(pm_mpflow_type) :: this

  call PMSubsurfaceFlowPreSolve(this)

end subroutine PMMpFlowPreSolve

! ************************************************************************** !

subroutine PMMpFlowPostSolve(this)
  ! 
  ! Date: 03/14/13
  ! 

  use Global_module

  implicit none
  
  class(pm_mpflow_type) :: this

  ! nothing ?

end subroutine PMMpFlowPostSolve

! ************************************************************************** !

subroutine PMMpFlowUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                              num_newton_iterations,tfac, &
                              time_step_max_growth_factor)
  ! 
  ! This routine
  ! 
  ! Author: F-M Yuan, ESD/CCSI-ORNL
  ! Date: 12/23/2019
  ! 
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL

  implicit none
  
  class(pm_mpflow_type) :: this
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

  call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt)
  
end subroutine PMMpFlowUpdateTimestep

! ************************************************************************** !

subroutine PMMpFlowResidual(this,snes,xx,r,ierr)
  ! 
  ! This routine
  ! 
  ! Author: F-M Yuan, ESD/CCSI-ORNL
  ! Date: 12/23/2019
  ! 

  use MpFlow_module, only : MpFlowResidual

  implicit none
  
  class(pm_mpflow_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr

  call MpFlowResidual(snes,xx,r,this%realization,ierr)

end subroutine PMMpFlowResidual

! ************************************************************************** !

subroutine PMMpFlowJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! This routine
  ! 
  ! Author: F-M Yuan, ESD/CCSI-ORNL
  ! Date: 12/23/2019
  ! 

  use MpFlow_module, only : MpFlowJacobian

  implicit none
  
  class(pm_mpflow_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  call MpFlowJacobian(A,B,this%realization,ierr)

end subroutine PMMpFlowJacobian

! ************************************************************************** !

subroutine PMMpFlowCheckUpdatePre(this,line_search,X,dX,changed,ierr)
  ! 
  ! This routine
  ! 
  ! Author: F-M Yuan, ESD/CCSI-ORNL
  ! Date: 12/23/2019
  ! 

  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Option_module
  use Patch_module

  implicit none
  
  class(pm_mpflow_type) :: this
  SNESLineSearch :: line_search
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
  PetscInt :: local_id, ghosted_id
  PetscReal :: P0, P1, P_R, delP
  PetscReal :: scale, press_limit, temp_limit
  PetscInt :: iend, istart

  patch => this%realization%patch
  grid => patch%grid
  option => this%realization%option
  field => this%realization%field

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
      else if (P1 < P_R .and. P0 > P_R) then
        write(option%io_buffer,'("S -> U:",1i7,2f12.1)') &
          grid%nG2A(grid%nL2G(local_id)),P0,P1
        call PrintMsgAnyRank(option)
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

end subroutine PMMpFlowCheckUpdatePre

! ************************************************************************** !

subroutine PMMpFlowCheckUpdatePost(this,line_search,X0,dX,X1,dX_changed, &
                                  X1_changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/21/18
  ! 
  ! A NOTE (fmyuan, 2018-10-17): this option may have mass-balance issue if not carefully
  !    setting the 'flow_itol_scaled_res' values. Temperally off now.

  use Option_module
  use Patch_module
  use Grid_module

  implicit none
  
  class(pm_mpflow_type) :: this
  SNESLineSearch :: line_search
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: X0_p(:)
  PetscReal, pointer :: dX_p(:)
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
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

  option => this%realization%option
  patch => this%realization%patch
  grid  => patch%grid

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

end subroutine PMMpFlowCheckUpdatePost

! ************************************************************************** !

subroutine PMMpFlowCheckConvergence(this,snes,it,xnorm,unorm,fnorm,reason,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 11/20/18
  !
  use Convergence_module
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

  class(pm_mpflow_type) :: this
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
      ghosted_id = grid%nL2G(local_id)
      natural_id = grid%nG2A(ghosted_id)
      if (patch%imat(ghosted_id) <= 0) cycle

      offset = (local_id-1)*option%nflowdof
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

end subroutine PMMpFlowCheckConvergence

! ************************************************************************** !

subroutine PMMpFlowTimeCut(this)
  ! 
  ! This routine
  ! 
  ! Author: F-M Yuan, ESD/CCSI-ORNL
  ! Date: 12/23/2019
  ! 

  use MpFlow_module, only : MpFlowTimeCut

  implicit none
  
  class(pm_mpflow_type) :: this
  
  call PMSubsurfaceFlowTimeCut(this)
  call MpFlowTimeCut(this%realization)

end subroutine PMMpFlowTimeCut

! ************************************************************************** !

subroutine PMMpFlowUpdateSolution(this)
  ! 
  ! This routine
  ! 
  ! Author: F-M Yuan, ESD/CCSI-ORNL
  ! Date: 12/23/2019
  ! 

  use MpFlow_module, only : MpFlowUpdateSolution

  implicit none
  
  class(pm_mpflow_type) :: this

  call PMSubsurfaceFlowUpdateSolution(this)
  call MpFlowUpdateSolution(this%realization)

end subroutine PMMpFlowUpdateSolution

! ************************************************************************** !

subroutine PMMpFlowUpdateAuxVars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14

  use MpFlow_module, only : MpFlowUpdateAuxVars
  
  implicit none
  
  class(pm_mpflow_type) :: this

  call MpFlowUpdateAuxVars(this%realization)

end subroutine PMMpFlowUpdateAuxVars

! ************************************************************************** !

subroutine PMMpFlowMaxChange(this)
  ! 
  ! This routine
  ! 
  ! Author: F-M Yuan, ESD/CCSI-ORNL
  ! Date: 12/23/2019
  ! 

  use MpFlow_module, only : MpFlowMaxChange
  use Option_module

  implicit none
  
  class(pm_mpflow_type) :: this
  character(len=MAXSTRINGLENGTH) :: string

  call MpFlowMaxChange(this%realization,this%max_pressure_change, &
                   this%max_temperature_change)

#ifndef CLM_PFLOTRAN
  write(string,'("  --> max chng: dpmx= ",1pe12.4," dtmpmx= ",1pe12.4)') &
      this%max_pressure_change,this%max_temperature_change
  call OptionPrint(string,this%option)
#endif

end subroutine PMMpFlowMaxChange

! ************************************************************************** !

subroutine PMMpFlowComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: F-M Yuan, ESD/CCSI-ORNL
  ! Date: 12/23/2019
  ! 

  use MpFlow_module, only : MpFlowComputeMassBalance

  implicit none
  
  class(pm_mpflow_type) :: this
  PetscReal,pointer :: mass_balance_array(:,:)
  
  call MpFlowComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMMpFlowComputeMassBalance

! ************************************************************************** !

! ************************************************************************** !

subroutine PMMpFlowDestroy(this)
  ! 
  ! Author: F-M Yuan, ESD/CCSI-ORNL
  ! Date: 12/23/2019
  ! 

  use MpFlow_module, only : MpFlowDestroy

  implicit none
  
  class(pm_mpflow_type) :: this

  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  call this%commN%Destroy()

  ! preserve this ordering
  call MpFlowDestroy(this%realization%patch)
  call PMSubsurfaceFlowDestroy(this)

end subroutine PMMpFlowDestroy

end module PM_MpFlow_class
