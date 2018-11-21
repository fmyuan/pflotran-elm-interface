module PM_General_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_subsurface_flow_type) :: pm_general_type
    PetscInt, pointer :: max_change_ivar(:)
    PetscInt, pointer :: max_change_isubvar(:)
    PetscInt :: converged_flag(3,3,3)
    PetscReal :: converged_real(3,3,3)
  contains
    procedure, public :: Read => PMGeneralRead
    procedure, public :: InitializeRun => PMGeneralInitializeRun
    procedure, public :: InitializeTimestep => PMGeneralInitializeTimestep
    procedure, public :: Residual => PMGeneralResidual
    procedure, public :: Jacobian => PMGeneralJacobian
    procedure, public :: UpdateTimestep => PMGeneralUpdateTimestep
    procedure, public :: PreSolve => PMGeneralPreSolve
    procedure, public :: PostSolve => PMGeneralPostSolve
    procedure, public :: CheckUpdatePre => PMGeneralCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMGeneralCheckUpdatePost
    procedure, public :: CheckConvergence => PMGeneralCheckConvergence
    procedure, public :: TimeCut => PMGeneralTimeCut
    procedure, public :: UpdateSolution => PMGeneralUpdateSolution
    procedure, public :: UpdateAuxVars => PMGeneralUpdateAuxVars
    procedure, public :: MaxChange => PMGeneralMaxChange
    procedure, public :: ComputeMassBalance => PMGeneralComputeMassBalance
    procedure, public :: InputRecord => PMGeneralInputRecord
    procedure, public :: CheckpointBinary => PMGeneralCheckpointBinary
    procedure, public :: RestartBinary => PMGeneralRestartBinary
    procedure, public :: Destroy => PMGeneralDestroy
  end type pm_general_type
  
  public :: PMGeneralCreate
  
contains

! ************************************************************************** !

function PMGeneralCreate()
  ! 
  ! Creates General process models shell
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 
  use Variables_module, only : LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
                               LIQUID_MOLE_FRACTION, TEMPERATURE, &
                               GAS_SATURATION
  use Upwind_Direction_module

  implicit none
  
  class(pm_general_type), pointer :: PMGeneralCreate

  class(pm_general_type), pointer :: general_pm
  
#ifdef PM_GENERAL_DEBUG  
  print *, 'PMGeneralCreate()'
#endif  

  allocate(general_pm)
  allocate(general_pm%max_change_ivar(6))
  general_pm%max_change_ivar = [LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
                                LIQUID_MOLE_FRACTION, TEMPERATURE, &
                                GAS_SATURATION]
  allocate(general_pm%max_change_isubvar(6))
                                   ! UNINITIALIZED_INTEGER avoids zeroing of 
                                   ! pressures not represented in phase
                                       ! 2 = air in xmol(air,liquid)
  general_pm%max_change_isubvar = [0,0,0,2,0,0]
  
  call PMSubsurfaceFlowCreate(general_pm)
  general_pm%name = 'General Multiphase Flow'
  general_pm%header = 'GENERAL MULTIPHASE FLOW'

  ! turn off default upwinding which is set to PETSC_TRUE in
  !  upwind_direction.F90
  fix_upwind_direction = PETSC_FALSE
  general_pm%check_post_convergence = PETSC_TRUE

  PMGeneralCreate => general_pm
  
end function PMGeneralCreate

! ************************************************************************** !

subroutine PMGeneralRead(this,input)
  ! 
  ! Sets up SNES solvers.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/04/15
  !
  use General_module
  use General_Aux_module
  use Input_Aux_module
  use String_module
  use Option_module

  implicit none
  
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: keyword, word
  class(pm_general_type) :: this
  type(option_type), pointer :: option
  PetscReal :: tempreal
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  option => this%option

  error_string = 'General Options'
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)
    
    found = PETSC_FALSE
    call PMSubsurfaceFlowReadSelectCase(this,input,keyword,found, &
                                        error_string,option)    
    if (found) cycle
    
    select case(trim(keyword))
      case('ITOL_SCALED_RESIDUAL')
        call InputReadDouble(input,option,general_itol_scaled_res)
        call InputDefaultMsg(input,option,'general_itol_scaled_res')
        this%check_post_convergence = PETSC_TRUE
      case('ITOL_RELATIVE_UPDATE')
        call InputReadDouble(input,option,general_itol_rel_update)
        call InputDefaultMsg(input,option,'general_itol_rel_update')
        this%check_post_convergence = PETSC_TRUE        
      case('TOUGH2_ITOL_SCALED_RESIDUAL')
        ! since general_tough2_itol_scaled_res_e1 is an array, we must read
        ! the tolerance into a double and copy it to the array.
        tempreal = UNINITIALIZED_DOUBLE
        call InputReadDouble(input,option,tempreal)
        ! tempreal will remain uninitialized if the read fails.
        call InputDefaultMsg(input,option,'tough_itol_scaled_residual_e1')
        if (Initialized(tempreal)) then
          general_tough2_itol_scaled_res_e1 = tempreal
        endif
        call InputReadDouble(input,option,general_tough2_itol_scaled_res_e2)
        call InputDefaultMsg(input,option,'tough_itol_scaled_residual_e2')
        general_tough2_conv_criteria = PETSC_TRUE
        this%check_post_convergence = PETSC_TRUE
      case('T2_ITOL_SCALED_RESIDUAL_TEMP')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option, &
                           'tough_itol_scaled_residual_e1 for temperature', &
                           error_string)
        general_tough2_itol_scaled_res_e1(3,:) = tempreal
      case('WINDOW_EPSILON') 
        call InputReadDouble(input,option,window_epsilon)
        call InputErrorMsg(input,option,'window epsilon',error_string)
      case('GAS_COMPONENT_FORMULA_WEIGHT')
        !geh: assuming gas component is index 2
        call InputReadDouble(input,option,fmw_comp(2))
        call InputErrorMsg(input,option,'gas component formula wt.', &
                           error_string)
      case('TWO_PHASE_ENERGY_DOF')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'two_phase_energy_dof',error_string)
        call GeneralAuxSetEnergyDOF(word,option)
      case('ISOTHERMAL')
        general_isothermal = PETSC_TRUE
      case('NO_AIR')
        general_no_air = PETSC_TRUE
      case('MAXIMUM_PRESSURE_CHANGE')
        call InputReadDouble(input,option,general_max_pressure_change)
        call InputErrorMsg(input,option,'maximum pressure change', &
                           error_string)
      case('MAX_ITERATION_BEFORE_DAMPING')
        call InputReadInt(input,option,general_max_it_before_damping)
        call InputErrorMsg(input,option,'maximum iteration before damping', &
                           error_string)
      case('DAMPING_FACTOR')
        call InputReadDouble(input,option,general_damping_factor)
        call InputErrorMsg(input,option,'damping factor',error_string)
#if 0        
      case('GOVERN_MAXIMUM_PRESSURE_CHANGE')
        call InputReadDouble(input,option,this%dPmax_allowable)
        call InputErrorMsg(input,option,'maximum allowable pressure change', &
                           error_string)
      case('GOVERN_MAXIMUM_TEMPERATURE_CHANGE')
        call InputReadDouble(input,option,this%dTmax_allowable)
        call InputErrorMsg(input,option, &
                           'maximum allowable temperature change', &
                           error_string)
      case('GOVERN_MAXIMUM_SATURATION_CHANGE')
        call InputReadDouble(input,option,this%dSmax_allowable)
        call InputErrorMsg(input,option,'maximum allowable saturation change', &
                           error_string)
      case('GOVERN_MAXIMUM_MOLE_FRACTION_CHANGE')
        call InputReadDouble(input,option,this%dXmax_allowable)
        call InputErrorMsg(input,option, &
                           'maximum allowable mole fraction change', &
                           error_string)
#endif
      case('DEBUG_CELL')
        call InputReadInt(input,option,general_debug_cell_id)
        call InputErrorMsg(input,option,'debug cell id',error_string)
      case('NO_TEMP_DEPENDENT_DIFFUSION')
        general_temp_dep_gas_air_diff = PETSC_FALSE
      case('DIFFUSE_XMASS')
        general_diffuse_xmol = PETSC_FALSE
      case('HARMONIC_GAS_DIFFUSIVE_DENSITY')
        general_harmonic_diff_density = PETSC_TRUE
      case('ARITHMETIC_GAS_DIFFUSIVE_DENSITY')
        general_harmonic_diff_density = PETSC_FALSE
      case('IMMISCIBLE')
        general_immiscible = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(keyword,'GENERAL Mode',option)
    end select
    
  enddo  

  if (general_isothermal .and. &
      general_2ph_energy_dof == GENERAL_AIR_PRESSURE_INDEX) then
    option%io_buffer = 'Isothermal GENERAL mode may only be run with ' // &
                       'temperature as the two phase energy dof.'
    call printErrMsg(option)
  endif

end subroutine PMGeneralRead

! ************************************************************************** !

recursive subroutine PMGeneralInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14 

  use Realization_Base_class
  
  implicit none
  
  class(pm_general_type) :: this
  
  PetscInt :: i
  PetscErrorCode :: ierr

  ! need to allocate vectors for max change
  call VecDuplicateVecsF90(this%realization%field%work,SIX_INTEGER, &
                           this%realization%field%max_change_vecs, &
                           ierr);CHKERRQ(ierr)
  ! set initial values
  do i = 1, 6
    call RealizationGetVariable(this%realization, &
                                this%realization%field%max_change_vecs(i), &
                                this%max_change_ivar(i), &
                                this%max_change_isubvar(i))
  enddo

  ! call parent implementation
  call PMSubsurfaceFlowInitializeRun(this)

end subroutine PMGeneralInitializeRun

! ************************************************************************** !

subroutine PMGeneralInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralInitializeTimestep
  use Global_module
  use Variables_module, only : TORTUOSITY
  use Material_module, only : MaterialAuxVarCommunicate
  use Option_module
  
  implicit none
  
  class(pm_general_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)                                 
!geh:remove   everywhere                                
  call MaterialAuxVarCommunicate(this%comm1, &
                                 this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc,TORTUOSITY, &
                                 ZERO_INTEGER)
                                 
  call GeneralInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)                                 
  
end subroutine PMGeneralInitializeTimestep

! ************************************************************************** !

subroutine PMGeneralPreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13

  implicit none

  class(pm_general_type) :: this

  call PMSubsurfaceFlowPreSolve(this)

end subroutine PMGeneralPreSolve

! ************************************************************************** !

subroutine PMGeneralPostSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13

  implicit none

  class(pm_general_type) :: this

end subroutine PMGeneralPostSolve

! ************************************************************************** !

subroutine PMGeneralUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                                   num_newton_iterations,tfac, &
                                   time_step_max_growth_factor)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Realization_Base_class, only : RealizationGetVariable
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL
  use Field_module
  use Global_module, only : GlobalSetAuxVarVecLoc
  use Variables_module, only : LIQUID_SATURATION, GAS_SATURATION

  implicit none
  
  class(pm_general_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor
  
  PetscReal :: fac
  PetscInt :: ifac
  PetscReal :: up, ut, ux, us, umin
  PetscReal :: dtt
  type(field_type), pointer :: field
  
#ifdef PM_GENERAL_DEBUG  
  call printMsg(this%option,'PMGeneral%UpdateTimestep()')
#endif
  
  fac = 0.5d0
  if (num_newton_iterations >= iacceleration) then
    fac = 0.33d0
    umin = 0.d0
  else
    up = this%pressure_change_governor/(this%max_pressure_change+0.1)
    ut = this%temperature_change_governor/(this%max_temperature_change+1.d-5)
    ux = this%xmol_change_governor/(this%max_xmol_change+1.d-5)
    us = this%saturation_change_governor/(this%max_saturation_change+1.d-5)
    umin = min(up,ut,ux,us)
  endif
  ifac = max(min(num_newton_iterations,size(tfac)),1)
  dtt = fac * dt * (1.d0 + umin)
  dtt = min(time_step_max_growth_factor*dt,dtt)
  dt = min(dtt,tfac(ifac)*dt,dt_max)
  dt = max(dt,dt_min)

  if (Initialized(this%cfl_governor)) then
    ! Since saturations are not stored in global_auxvar for general mode, we
    ! must copy them over for the CFL check
    ! liquid saturation
    field => this%realization%field
    call RealizationGetVariable(this%realization,field%work, &
                                LIQUID_SATURATION,ZERO_INTEGER)
    call this%realization%comm1%GlobalToLocal(field%work,field%work_loc)
    call GlobalSetAuxVarVecLoc(this%realization,field%work_loc, &
                               LIQUID_SATURATION,TIME_NULL)
    call RealizationGetVariable(this%realization,field%work, &
                                GAS_SATURATION,ZERO_INTEGER)
    call this%realization%comm1%GlobalToLocal(field%work,field%work_loc)
    call GlobalSetAuxVarVecLoc(this%realization,field%work_loc, &
                               GAS_SATURATION,TIME_NULL)
    call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt)
  endif

end subroutine PMGeneralUpdateTimestep

! ************************************************************************** !

subroutine PMGeneralResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralResidual

  implicit none
  
  class(pm_general_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  call PMSubsurfaceFlowUpdatePropertiesNI(this)
  call GeneralResidual(snes,xx,r,this%realization,ierr)

end subroutine PMGeneralResidual

! ************************************************************************** !

subroutine PMGeneralJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralJacobian

  implicit none
  
  class(pm_general_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  call GeneralJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMGeneralJacobian

! ************************************************************************** !

subroutine PMGeneralCheckUpdatePre(this,line_search,X,dX,changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Option_module
  use Saturation_Function_module
  use Patch_module
  use General_Aux_module
  use Global_Aux_module
  
  implicit none
  
  class(pm_general_type) :: this
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
  type(general_auxvar_type), pointer :: gen_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
#ifdef DEBUG_GENERAL_INFO
  character(len=MAXSTRINGLENGTH) :: string, string2, string3
  character(len=MAXWORDLENGTH) :: cell_id_word
  PetscInt, parameter :: max_cell_id = 10
  PetscInt :: cell_id, cell_locator(0:max_cell_id)
  PetscInt :: i, ii
#endif
  PetscInt :: local_id, ghosted_id
  PetscInt :: offset
  PetscInt :: liquid_pressure_index, gas_pressure_index, air_pressure_index
  PetscInt :: saturation_index, temperature_index, xmol_index
  PetscInt :: lid, gid, apid, cpid, vpid, spid
  PetscInt :: pgas_index
  PetscReal :: liquid_pressure0, liquid_pressure1, del_liquid_pressure
  PetscReal :: gas_pressure0, gas_pressure1, del_gas_pressure
  PetscReal :: air_pressure0, air_pressure1, del_air_pressure
  PetscReal :: temperature0, temperature1, del_temperature
  PetscReal :: saturation0, saturation1, del_saturation
  PetscReal :: xmol0, xmol1, del_xmol
  PetscReal :: max_saturation_change = 0.125d0
  PetscReal :: max_temperature_change = 10.d0
  PetscReal :: min_pressure
  PetscReal :: scale, temp_scale, temp_real
  PetscReal, parameter :: tolerance = 0.99d0
  PetscReal, parameter :: initial_scale = 1.d0
  PetscReal, parameter :: ALMOST_ZERO = 1.d-10
  PetscReal, parameter :: ALMOST_ONE = 1.d0-ALMOST_ZERO
  SNES :: snes
  PetscInt :: newton_iteration
  
  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  gen_auxvars => this%realization%patch%aux%General%auxvars
  global_auxvars => this%realization%patch%aux%Global%auxvars

  patch => this%realization%patch

  spid = option%saturation_pressure_id
  apid = option%air_pressure_id

  call SNESLineSearchGetSNES(line_search,snes,ierr)
  call SNESGetIterationNumber(snes,newton_iteration,ierr)

  call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

  changed = PETSC_TRUE

  call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

end subroutine PMGeneralCheckUpdatePre

! ************************************************************************** !

subroutine PMGeneralCheckUpdatePost(this,line_search,X0,dX,X1,dX_changed, &
                                    X1_changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/20/18
  ! 
  use General_Aux_module
  use Global_Aux_module
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  
  implicit none
  
  class(pm_general_type) :: this
  SNESLineSearch :: line_search
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
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  PetscInt :: local_id, ghosted_id, natural_id
  PetscInt :: offset, ival, idof
  PetscReal :: dX_, dX_X0
  PetscReal, parameter :: inf_pres_tol = 1.d0
  PetscReal, parameter :: inf_temp_tol = 1.d-5
  PetscReal, parameter :: inf_sat_tol = 1.d-6
  PetscReal, parameter :: inf_xmol_tol = 1.d-6
  PetscReal, parameter :: inf_norm_update_tol(3,3) = &
    reshape([inf_pres_tol,inf_xmol_tol,inf_temp_tol, &
             inf_pres_tol,inf_pres_tol,inf_temp_tol, &
             inf_pres_tol,inf_pres_tol,inf_sat_tol], &
            shape(inf_norm_update_tol)) * &
            1.d0 ! change to 0.d0 to zero tolerances
  PetscReal, parameter :: inf_rel_pres_tol = 1.d-3
  PetscReal, parameter :: inf_rel_temp_tol = 1.d-3
  PetscReal, parameter :: inf_rel_sat_tol = 1.d-3
  PetscReal, parameter :: inf_rel_xmol_tol = 1.d-3
  PetscReal, parameter :: inf_norm_rel_update_tol(3,3) = &
    reshape([inf_rel_pres_tol,inf_rel_xmol_tol,inf_rel_temp_tol, &
             inf_rel_pres_tol,inf_rel_pres_tol,inf_rel_temp_tol, &
             inf_rel_pres_tol,inf_rel_pres_tol,inf_rel_sat_tol], &
            shape(inf_norm_rel_update_tol)) * &
            1.d0 ! change to 0.d0 to zero tolerances
  PetscInt :: converged_update_flag(3,3)
  PetscInt :: converged_rel_update_flag(3,3)
  PetscReal :: converged_update_real(3,3)
  PetscReal :: converged_rel_update_real(3,3)
  PetscInt :: istate
  PetscBool :: converged_absolute
  PetscBool :: converged_relative
  
  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch
  global_auxvars => patch%aux%Global%auxvars
  
  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE

  call VecGetArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)
  converged_update_flag = ZERO_INTEGER
  converged_rel_update_flag = ZERO_INTEGER
  converged_update_real = 0.d0
  converged_rel_update_real = 0.d0
  do local_id = 1, grid%nlmax
    offset = (local_id-1)*option%nflowdof
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    istate = global_auxvars(ghosted_id)%istate
    do idof = 1, option%nflowdof
      ival = offset+idof
      ! infinity norms on update
      converged_absolute = PETSC_TRUE
      converged_relative = PETSC_TRUE
      dX_ = dabs(dX_p(ival))
      dX_X0 = dabs(dX_/X0_p(ival))
      if (dX_ > inf_norm_update_tol(idof,istate)) then
        converged_absolute = PETSC_FALSE
      endif
      if (dX_X0 > inf_norm_rel_update_tol(idof,istate)) then
        converged_relative = PETSC_FALSE
      endif
      if (.not.(converged_absolute .or. converged_relative)) then
        natural_id = grid%nG2A(ghosted_id)
        if (.not.converged_absolute .and. &
            converged_update_real(idof,istate) < dX_) then
          converged_update_real(idof,istate) = dX_
          converged_update_flag(idof,istate) = natural_id
        endif
        if (.not.converged_relative .and. &
            converged_rel_update_real(idof,istate) < dX_) then
          converged_rel_update_real(idof,istate) = dX_X0
          converged_rel_update_flag(idof,istate) = natural_id
        endif
      endif
    enddo
  enddo
  call VecRestoreArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)

  this%converged_flag(:,:,1) = converged_update_flag(:,:)
  this%converged_flag(:,:,2) = converged_rel_update_flag(:,:)
  this%converged_real(:,:,1) = converged_update_real(:,:)
  this%converged_real(:,:,2) = converged_rel_update_real(:,:)

end subroutine PMGeneralCheckUpdatePost

! ************************************************************************** !

subroutine PMGeneralCheckConvergence(this,snes,it,xnorm,unorm,fnorm, &
                                     reason,ierr)
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

  class(pm_general_type) :: this
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
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum2_p(:)
  PetscInt :: local_id, ghosted_id
  PetscInt :: offset, ival, idof, itol
  PetscReal :: R, A, R_A
  PetscReal, parameter :: A_zero = 1.d-15
  PetscReal, parameter :: inf_tol(3) = 1.d-5
  PetscInt :: converged_residual_flag(3,3)
  PetscReal :: converged_residual_real(3,3)
  PetscInt :: istate
  PetscMPIInt :: mpi_int
  character(len=MAXSTRINGLENGTH) :: string
  character(len=12), parameter :: state_string(3) = &
    ['Liquid State','Gas State   ','2Phase State']
  character(len=17), parameter :: dof_string(3,3) = &
    reshape(['Liquid Pressure  ','Air Mole Fraction','Temperature      ', &
             'Gas Pressure     ','Gas Saturation   ','Temperature      ', &
             'Gas Pressure     ','Air Pressure     ','Temperature      '], &
             shape(dof_string))
  character(len=15), parameter :: tol_string(3) = &
    ['Absolute Update','Relative Update','Residual       ']
  
  patch => this%realization%patch
  option => this%realization%option
  field => this%realization%field
  grid => patch%grid
  global_auxvars => patch%aux%Global%auxvars

  call VecGetArrayReadF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%flow_accum2,accum2_p,ierr);CHKERRQ(ierr)
  converged_residual_flag = ZERO_INTEGER
  converged_residual_real = 0.d0
  do local_id = 1, grid%nlmax
    offset = (local_id-1)*option%nflowdof
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    istate = global_auxvars(ghosted_id)%istate
    do idof = 1, option%nflowdof
      ival = offset+idof
      ! infinity norms on residual
      R = dabs(r_p(ival))
      A = dabs(accum2_p(ival))
      R_A = R/A
      if (A > A_zero) then
        if (R > inf_tol(idof)) then
          if (R_A > inf_tol(idof)) then
            if (converged_residual_real(idof,istate) <  R) then
              converged_residual_real(idof,istate) = R
              converged_residual_flag(idof,istate) = grid%nG2A(ghosted_id)
            endif
          endif
        endif
      endif
    enddo
  enddo
  call VecRestoreArrayReadF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%flow_accum2,accum2_p,ierr);CHKERRQ(ierr)

  this%converged_flag(:,:,3) = converged_residual_flag(:,:)
  this%converged_real(:,:,3) = converged_residual_real(:,:)
  mpi_int = 27
  call MPI_Allreduce(MPI_IN_PLACE,this%converged_flag,mpi_int, &
                     MPI_INTEGER,MPI_MAX,option%mycomm,ierr)
  call MPI_Allreduce(MPI_IN_PLACE,this%converged_real,mpi_int, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)

  option%convergence = CONVERGENCE_CONVERGED
  if (maxval(this%converged_flag) > ZERO_INTEGER) then
    option%convergence = CONVERGENCE_KEEP_ITERATING
  endif

  do itol = 1, 3 
    do istate = 1, 3
      do idof = 1, 3
        if (this%converged_flag(idof,istate,itol) /= ZERO_INTEGER) then
          string = '   ' // trim(tol_string(itol)) // ', ' // &
            trim(state_string(istate)) // ', ' // &
            trim(dof_string(idof,istate)) // ' (' // &
            trim(StringFormatInt(this%converged_flag(idof,istate,itol))) // &
            ') : ' // &
            StringFormatDouble(this%converged_real(idof,istate,itol))
          call OptionPrint(string,option)
        endif
      enddo
    enddo
  enddo

  call PMSubsurfaceFlowCheckConvergence(this,snes,it,xnorm,unorm,fnorm, &
                                        reason,ierr)

end subroutine PMGeneralCheckConvergence

! ************************************************************************** !

subroutine PMGeneralTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralTimeCut

  implicit none
  
  class(pm_general_type) :: this
  
  call PMSubsurfaceFlowTimeCut(this)
  call GeneralTimeCut(this%realization)

end subroutine PMGeneralTimeCut

! ************************************************************************** !

subroutine PMGeneralUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralUpdateSolution, &
                             GeneralMapBCAuxVarsToGlobal

  implicit none
  
  class(pm_general_type) :: this
  
  call PMSubsurfaceFlowUpdateSolution(this)
  call GeneralUpdateSolution(this%realization)
  call GeneralMapBCAuxVarsToGlobal(this%realization)

end subroutine PMGeneralUpdateSolution     

! ************************************************************************** !

subroutine PMGeneralUpdateAuxVars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/21/14
  use General_module, only : GeneralUpdateAuxVars

  implicit none
  
  class(pm_general_type) :: this

  call GeneralUpdateAuxVars(this%realization,PETSC_FALSE)

end subroutine PMGeneralUpdateAuxVars   

! ************************************************************************** !

subroutine PMGeneralMaxChange(this)
  ! 
  ! Not needed given GeneralMaxChange is called in PostSolve
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use Realization_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use Global_Aux_module
  use General_Aux_module
  use Variables_module, only : LIQUID_PRESSURE, LIQUID_MOLE_FRACTION, &
                               TEMPERATURE, GAS_PRESSURE, AIR_PRESSURE, &
                               GAS_SATURATION
  implicit none
  
  class(pm_general_type) :: this
  
  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_ptr(:), vec_ptr2(:)
  PetscReal :: max_change_local(6)
  PetscReal :: max_change_global(6)
  PetscReal :: max_change
  PetscInt :: i, j
  PetscInt :: local_id, ghosted_id

  
  PetscErrorCode :: ierr
  
  realization => this%realization
  option => realization%option
  field => realization%field
  grid => realization%patch%grid

  max_change_global = 0.d0
  max_change_local = 0.d0
  
  ! max change variables: [LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
  !                        LIQUID_MOLE_FRACTION, TEMPERATURE, GAS_SATURATION]
  do i = 1, 6
    call RealizationGetVariable(realization,field%work, &
                                this%max_change_ivar(i), &
                                this%max_change_isubvar(i))
    ! yes, we could use VecWAXPY and a norm here, but we need the ability
    ! to customize
    call VecGetArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(field%max_change_vecs(i),vec_ptr2,ierr);CHKERRQ(ierr)
    max_change = 0.d0
    do j = 1, grid%nlmax
      ! have to weed out cells that changed state
      if (dabs(vec_ptr(j)) > 1.d-40 .and. dabs(vec_ptr2(j)) > 1.d-40) then
        max_change = max(max_change,dabs(vec_ptr(j)-vec_ptr2(j)))
      endif
    enddo
    max_change_local(i) = max_change
    call VecRestoreArrayF90(field%work,vec_ptr,ierr);CHKERRQ(ierr)
    call VecRestoreArrayF90(field%max_change_vecs(i),vec_ptr2, &
                            ierr);CHKERRQ(ierr)
    call VecCopy(field%work,field%max_change_vecs(i),ierr);CHKERRQ(ierr)
  enddo
  call MPI_Allreduce(max_change_local,max_change_global,SIX_INTEGER, &
                      MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
  ! print them out
  if (OptionPrintToScreen(option)) then
    write(*,'("  --> max chng: dpl= ",1pe12.4, " dpg= ",1pe12.4,&
      & " dpa= ",1pe12.4,/,15x," dxa= ",1pe12.4,"  dt= ",1pe12.4,&
      & " dsg= ",1pe12.4)') &
      max_change_global(1:6)
  endif
  if (OptionPrintToFile(option)) then
    write(option%fid_out,'("  --> max chng: dpl= ",1pe12.4, " dpg= ",1pe12.4,&
      & " dpa= ",1pe12.4,/,15x," dxa= ",1pe12.4,"  dt= ",1pe12.4, &
      & " dsg= ",1pe12.4)') &
      max_change_global(1:6)
  endif

  ! max change variables: [LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
  !                        LIQUID_MOLE_FRACTION, TEMPERATURE, GAS_SATURATION]
  ! ignore air pressure as it jumps during phase change
  this%max_pressure_change = maxval(max_change_global(1:2))
  this%max_xmol_change = max_change_global(4)
  this%max_temperature_change = max_change_global(5)
  this%max_saturation_change = max_change_global(6)
  
end subroutine PMGeneralMaxChange

! ************************************************************************** !

subroutine PMGeneralComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralComputeMassBalance

  implicit none
  
  class(pm_general_type) :: this
  PetscReal :: mass_balance_array(:)
  
  call GeneralComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMGeneralComputeMassBalance

! ************************************************************************** !

subroutine PMGeneralInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  ! 
  
  implicit none
  
  class(pm_general_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'general'
  if (this%check_post_convergence) then
    write(id,'(a29)',advance='no') 'ITOL_SCALED_RESIDUAL: '
    write(id,'(a)') 'ON'
    write(id,'(a29)',advance='no') 'ITOL_RELATIVE_RESIDUAL: '
    write(id,'(a)') 'ON'
  endif

end subroutine PMGeneralInputRecord

! ************************************************************************** !

subroutine PMGeneralCheckpointBinary(this,viewer)
  ! 
  ! Checkpoints data associated with General PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/18/15

  use Checkpoint_module
  use Global_module
  use Variables_module, only : STATE

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_general_type) :: this
  PetscViewer :: viewer
  
  call PMSubsurfaceFlowCheckpointBinary(this,viewer)
  
end subroutine PMGeneralCheckpointBinary

! ************************************************************************** !

subroutine PMGeneralRestartBinary(this,viewer)
  ! 
  ! Restarts data associated with General PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/18/15

  use Checkpoint_module
  use Global_module
  use Variables_module, only : STATE

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_general_type) :: this
  PetscViewer :: viewer
  
  call PMSubsurfaceFlowRestartBinary(this,viewer)
  
end subroutine PMGeneralRestartBinary
! ************************************************************************** !

subroutine PMGeneralDestroy(this)
  ! 
  ! Destroys General process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  use General_module, only : GeneralDestroy

  implicit none
  
  class(pm_general_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  deallocate(this%max_change_ivar)
  nullify(this%max_change_ivar)
  deallocate(this%max_change_isubvar)
  nullify(this%max_change_isubvar)

  ! preserve this ordering
  call GeneralDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)
  
end subroutine PMGeneralDestroy
  
end module PM_General_class
