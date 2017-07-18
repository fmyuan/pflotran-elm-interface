module PM_WIPP_Flow_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use PM_Base_class
  use PM_Subsurface_Flow_class
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_subsurface_flow_type) :: pm_wippflo_type
    PetscInt, pointer :: max_change_ivar(:)
    PetscReal :: liquid_equation_tolerance
    PetscReal :: gas_equation_tolerance
    PetscReal :: liquid_pressure_tolerance
    PetscReal :: gas_saturation_tolerance
  contains
    procedure, public :: Read => PMWIPPFloRead
    procedure, public :: InitializeRun => PMWIPPFloInitializeRun
    procedure, public :: InitializeTimestep => PMWIPPFloInitializeTimestep
    procedure, public :: Residual => PMWIPPFloResidual
    procedure, public :: Jacobian => PMWIPPFloJacobian
    procedure, public :: UpdateTimestep => PMWIPPFloUpdateTimestep
    procedure, public :: PreSolve => PMWIPPFloPreSolve
    procedure, public :: PostSolve => PMWIPPFloPostSolve
    procedure, public :: CheckUpdatePre => PMWIPPFloCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMWIPPFloCheckUpdatePost
    procedure, public :: TimeCut => PMWIPPFloTimeCut
    procedure, public :: UpdateSolution => PMWIPPFloUpdateSolution
    procedure, public :: UpdateAuxVars => PMWIPPFloUpdateAuxVars
    procedure, public :: MaxChange => PMWIPPFloMaxChange
    procedure, public :: ComputeMassBalance => PMWIPPFloComputeMassBalance
    procedure, public :: InputRecord => PMWIPPFloInputRecord
    procedure, public :: CheckpointBinary => PMWIPPFloCheckpointBinary
    procedure, public :: RestartBinary => PMWIPPFloRestartBinary
    procedure, public :: Destroy => PMWIPPFloDestroy
  end type pm_wippflo_type
  
  public :: PMWIPPFloCreate
  
contains

! ************************************************************************** !

function PMWIPPFloCreate()
  ! 
  ! Creates WIPPFlo process models shell
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
  use Variables_module, only : LIQUID_PRESSURE, GAS_PRESSURE, AIR_PRESSURE, &
                               LIQUID_MOLE_FRACTION, TEMPERATURE, &
                               GAS_SATURATION
  implicit none
  
  class(pm_wippflo_type), pointer :: PMWIPPFloCreate

  class(pm_wippflo_type), pointer :: wippflo_pm
  
#ifdef PM_WIPPFLO_DEBUG  
  print *, 'PMWIPPFloCreate()'
#endif  

  allocate(wippflo_pm)
  allocate(wippflo_pm%max_change_ivar(3))
  wippflo_pm%max_change_ivar = [LIQUID_PRESSURE, GAS_PRESSURE, GAS_SATURATION]
  
  call PMSubsurfaceFlowCreate(wippflo_pm)
  wippflo_pm%name = 'WIPP Immiscible Multiphase Flow'

  wippflo_pm%check_post_convergence = PETSC_TRUE

  ! defaults from BRAGFLO input deck
  wippflo_pm%liquid_equation_tolerance = 1.d-2
  wippflo_pm%gas_equation_tolerance = 1.d-2
  wippflo_pm%liquid_pressure_tolerance = 1.d-2
  wippflo_pm%gas_saturation_tolerance = 1.d-3

  PMWIPPFloCreate => wippflo_pm
  
end function PMWIPPFloCreate

! ************************************************************************** !

subroutine PMWIPPFloRead(this,input)
  ! 
  ! Sets up SNES solvers.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  !
  use WIPP_Flow_module
  use WIPP_Flow_Aux_module
  use Input_Aux_module
  use String_module
  use Option_module

  implicit none
  
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: keyword, word
  class(pm_wippflo_type) :: this
  type(option_type), pointer :: option
  PetscReal :: tempreal
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  option => this%option

  error_string = 'WIPP Flow Options'
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)
    
    found = PETSC_FALSE
    call PMSubsurfaceFlowReadSelectCase(this,input,keyword,found,option)    
    if (found) cycle
    
    select case(trim(keyword))
      case('LIQUID_EQUATION_TOLERANCE')
        call InputReadDouble(input,option,this%liquid_equation_tolerance)
        call InputDefaultMsg(input,option,'LIQUID_EQUATION_TOLERANCE')
      case('GAS_EQUATION_TOLERANCE')
        call InputReadDouble(input,option,this%gas_equation_tolerance)
        call InputDefaultMsg(input,option,'GAS_EQUATION_TOLERANCE')
      case('LIQUID_PRESSURE_TOLERANCE')
        call InputReadDouble(input,option,this%liquid_pressure_tolerance)
        call InputDefaultMsg(input,option,'LIQUID_PRESSURE_TOLERANCE')
      case('GAS_SATURATION_TOLERANCE')
        call InputReadDouble(input,option,this%gas_saturation_tolerance)
        call InputDefaultMsg(input,option,'GAS_SATURATION_TOLERANCE')
      case('GAS_COMPONENT_FORMULA_WEIGHT')
        call InputReadDouble(input,option,fmw_comp(2))
        call InputErrorMsg(input,option,'gas component formula wt.', &
                           error_string)
      case('MAXIMUM_PRESSURE_CHANGE')
        call InputReadDouble(input,option,wippflo_max_pressure_change)
        call InputErrorMsg(input,option,'maximum pressure change', &
                           error_string)
      case('MAX_ITERATION_BEFORE_DAMPING')
        call InputReadInt(input,option,wippflo_max_it_before_damping)
        call InputErrorMsg(input,option,'maximum iteration before damping', &
                           error_string)
      case('DAMPING_FACTOR')
        call InputReadDouble(input,option,wippflo_damping_factor)
        call InputErrorMsg(input,option,'damping factor',error_string)
      case('FIX_UPWIND_DIRECTION')
        wippflo_fix_upwind_direction = PETSC_TRUE
      case('UNFIX_UPWIND_DIRECTION')
        wippflo_fix_upwind_direction = PETSC_FALSE
      case('COUNT_UPWIND_DIRECTION_FLIP')
        wippflo_count_upwind_dir_flip = PETSC_TRUE
      case('NO_FRACTURE')
        wippflo_use_fracture = PETSC_FALSE
      case('NO_CREEP_CLOSURE')
        wippflo_use_creep_closure = PETSC_FALSE
      case default
        call InputKeywordUnrecognized(keyword,'WIPP Flow Mode',option)
    end select
    
  enddo  
  
end subroutine PMWIPPFloRead

! ************************************************************************** !

recursive subroutine PMWIPPFloInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17

  use Realization_Base_class
  
  implicit none
  
  class(pm_wippflo_type) :: this
  
  PetscInt :: i
  PetscErrorCode :: ierr

  ! need to allocate vectors for max change
  call VecDuplicateVecsF90(this%realization%field%work,SIX_INTEGER, &
                           this%realization%field%max_change_vecs, &
                           ierr);CHKERRQ(ierr)
  ! set initial values
  do i = 1, 3
    call RealizationGetVariable(this%realization, &
                                this%realization%field%max_change_vecs(i), &
                                this%max_change_ivar(i),ZERO_INTEGER)
  enddo

  ! call parent implementation
  call PMSubsurfaceFlowInitializeRun(this)

end subroutine PMWIPPFloInitializeRun

! ************************************************************************** !

subroutine PMWIPPFloInitializeTimestep(this)
  ! 
  ! Should not need this as it is called in PreSolve.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use WIPP_Flow_module, only : WIPPFloInitializeTimestep
  use Global_module
  use Variables_module, only : TORTUOSITY
  use Material_module, only : MaterialAuxVarCommunicate
  
  implicit none
  
  class(pm_wippflo_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)                                 
  if (this%option%print_screen_flag) then
    write(*,'(/,2("=")," WIPP FLOW MODE ",62("="))')
  endif
  
  call WIPPFloInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)                                 
  
end subroutine PMWIPPFloInitializeTimestep

! ************************************************************************** !

subroutine PMWIPPFloPreSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17

  implicit none

  class(pm_wippflo_type) :: this

  call PMSubsurfaceFlowPreSolve(this)

end subroutine PMWIPPFloPreSolve

! ************************************************************************** !

subroutine PMWIPPFloPostSolve(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  
  use WIPP_Flow_Common_module
  use WIPP_Flow_Aux_module
  use Option_module

  implicit none

  class(pm_wippflo_type) :: this

  PetscInt, save :: lr = 0, gr = 0, lbr = 0, gbr = 0
  PetscInt, save :: lj = 0, gj = 0, lbj = 0, gbj = 0

  if (wippflo_fix_upwind_direction .and. &
      wippflo_count_upwind_dir_flip .and. &
      OptionPrintToScreen(this%realization%option)) then
    write(*,'(6x,"Res: ",4i5," : ",4i7)') &
      liq_upwind_flip_count_by_res-lr, &
      gas_upwind_flip_count_by_res-gr, &
      liq_bc_upwind_flip_count_by_res-lbr, &
      gas_bc_upwind_flip_count_by_res-gbr, &
      liq_upwind_flip_count_by_res, &
      gas_upwind_flip_count_by_res, &
      liq_bc_upwind_flip_count_by_res, &
      gas_bc_upwind_flip_count_by_res
    write(*,'(6x,"Jac: ",4i5," : ",4i7)') &
      liq_upwind_flip_count_by_jac-lj, &
      gas_upwind_flip_count_by_jac-gj, &
      liq_bc_upwind_flip_count_by_jac-lbj, &
      gas_bc_upwind_flip_count_by_jac-gbj, &
      liq_upwind_flip_count_by_jac, &
      gas_upwind_flip_count_by_jac, &
      liq_bc_upwind_flip_count_by_jac, &
      gas_bc_upwind_flip_count_by_jac
  endif

  lr = liq_upwind_flip_count_by_res
  gr = gas_upwind_flip_count_by_res
  lbr = liq_bc_upwind_flip_count_by_res
  gbr = gas_bc_upwind_flip_count_by_res
  lj = liq_upwind_flip_count_by_jac
  gj = gas_upwind_flip_count_by_jac
  lbj = liq_bc_upwind_flip_count_by_jac
  gbj = gas_bc_upwind_flip_count_by_jac

end subroutine PMWIPPFloPostSolve

! ************************************************************************** !

subroutine PMWIPPFloUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                                    num_newton_iterations,tfac)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Realization_Base_class, only : RealizationGetVariable
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL
  use Field_module
  use Global_module, only : GlobalSetAuxVarVecLoc
  use Variables_module, only : LIQUID_SATURATION, GAS_SATURATION

  implicit none
  
  class(pm_wippflo_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  
  PetscReal :: fac
  PetscInt :: ifac
  PetscReal :: up, us, umin
  PetscReal :: dtt
  type(field_type), pointer :: field
  
  fac = 0.5d0
  if (num_newton_iterations >= iacceleration) then
    fac = 0.33d0
    umin = 0.d0
  else
    up = this%pressure_change_governor/(this%max_pressure_change+0.1)
    us = this%saturation_change_governor/(this%max_saturation_change+1.d-5)
    umin = min(up,us)
  endif
  ifac = max(min(num_newton_iterations,size(tfac)),1)
  dtt = fac * dt * (1.d0 + umin)
  dt = min(dtt,tfac(ifac)*dt,dt_max)
  dt = max(dt,dt_min)

  if (Initialized(this%cfl_governor)) then
    ! Since saturations are not stored in global_auxvar for wipp flow mode, we
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

end subroutine PMWIPPFloUpdateTimestep

! ************************************************************************** !

subroutine PMWIPPFloResidual(this,snes,xx,r,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use WIPP_Flow_module, only : WIPPFloResidual

  implicit none
  
  class(pm_wippflo_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  
  call PMSubsurfaceFlowUpdatePropertiesNI(this)
  call WIPPFloResidual(snes,xx,r,this%realization,ierr)

  call this%PostSolve()

end subroutine PMWIPPFloResidual

! ************************************************************************** !

subroutine PMWIPPFloJacobian(this,snes,xx,A,B,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use WIPP_Flow_module, only : WIPPFloJacobian

  implicit none
  
  class(pm_wippflo_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  
  call WIPPFloJacobian(snes,xx,A,B,this%realization,ierr)

end subroutine PMWIPPFloJacobian

! ************************************************************************** !

subroutine PMWIPPFloCheckUpdatePre(this,line_search,X,dX,changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Option_module
  use Saturation_Function_module
  use Patch_module
  use WIPP_Flow_Aux_module
  use Global_Aux_module
  
  implicit none
  
  class(pm_wippflo_type) :: this
  SNESLineSearch :: line_search
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr
  
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  PetscReal, pointer :: X_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: offset
  PetscInt :: saturation_index
  PetscInt :: temp_int(2)
  PetscBool :: cut_timestep
  PetscBool :: force_another_iteration
  PetscReal :: saturation0, saturation1, del_saturation

  SNES :: snes
  
  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  wippflo_auxvars => this%realization%patch%aux%WIPPFlo%auxvars
  global_auxvars => this%realization%patch%aux%Global%auxvars

  patch => this%realization%patch

  call SNESLineSearchGetSNES(line_search,snes,ierr)

  call VecGetArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

  changed = PETSC_TRUE
  cut_timestep = PETSC_FALSE
  force_another_iteration = PETSC_FALSE

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    offset = (local_id-1)*option%nflowdof
    saturation_index = offset + WIPPFLO_GAS_SATURATION_DOF
    del_saturation = dX_p(saturation_index)
    saturation0 = X_p(saturation_index)
    saturation1 = saturation0 - del_saturation
    if (saturation1 < 0.d0) then
      if (saturation1 < -0.2d0) then
        cut_timestep = PETSC_TRUE
      else 
        if (saturation1 < -1.d-3) then
          force_another_iteration = PETSC_TRUE
        endif
        ! set saturation to zero
        dX_p(saturation_index) = saturation0
      endif
    else if (saturation1 > 1.d0) then
      if (saturation1 > 1.d0 + 0.2d0) then
        cut_timestep = PETSC_TRUE
      else 
        if (saturation1 > 1.d0 + 1.d-3) then
          force_another_iteration = PETSC_TRUE
        endif
        ! set saturation to one
        dX_p(saturation_index) = saturation0 - 1.d0
      endif
    endif
  enddo

  call VecRestoreArrayF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X,X_p,ierr);CHKERRQ(ierr)

  temp_int(:) = 0
  if (force_another_iteration) temp_int(1) = 1
  if (cut_timestep) temp_int(2) = 1
  call MPI_Allreduce(MPI_IN_PLACE,temp_int,TWO_INTEGER_MPI, &
                     MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
  option%convergence = CONVERGENCE_KEEP_ITERATING
  if (temp_int(1) > 0) option%convergence = CONVERGENCE_FORCE_ITERATION
  if (temp_int(2) > 0) option%convergence = CONVERGENCE_CUT_TIMESTEP

end subroutine PMWIPPFloCheckUpdatePre

! ************************************************************************** !

subroutine PMWIPPFloCheckUpdatePost(this,line_search,X0,dX,X1,dX_changed, &
                                    X1_changed,ierr)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 
  use WIPP_Flow_Aux_module
  use Global_Aux_module
  use Grid_module
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Field_module
  use Patch_module
  use Option_module
  use Material_Aux_class  
  use Output_EKG_module
  
  implicit none
  
  class(pm_wippflo_type) :: this
  SNESLineSearch :: line_search
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(wippflo_auxvar_type), pointer :: wippflo_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)  
  class(material_auxvar_type), pointer :: material_auxvars(:)  
  type(material_parameter_type), pointer :: material_parameter

  PetscReal, pointer :: X0_p(:)
  PetscReal, pointer :: X1_p(:)
  PetscReal, pointer :: dX_p(:)
  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum_p(:), accum_p2(:)

  PetscInt :: local_id, ghosted_id
  PetscInt :: offset , idof
  PetscInt :: liquid_equation_index
  PetscInt :: gas_equation_index
  PetscInt :: saturation_index
  PetscInt :: pressure_index
  PetscInt :: temp_int(5)
  PetscReal :: temp_real(4)

  PetscReal, parameter :: zero_saturation = 1.d-15
  PetscReal, parameter :: zero_accumulation = 1.d-15

  PetscReal :: residual
  PetscReal :: accumulation
  PetscReal :: residual_over_accumulation
  PetscReal :: abs_X
  PetscReal :: abs_dX
  PetscBool :: converged_liquid_equation
  PetscBool :: converged_gas_equation
  PetscBool :: converged_liquid_pressure
  PetscBool :: converged_gas_saturation
  PetscReal :: max_liq_eq
  PetscReal :: max_gas_eq
  PetscReal :: max_liq_pres_rel_change
  PetscReal :: max_gas_sat_change
  PetscReal :: min_liq_pressure
  PetscReal :: min_gas_pressure
  PetscInt :: max_liq_eq_cell
  PetscInt :: max_gas_eq_cell
  PetscInt :: max_liq_pres_rel_change_cell
  PetscInt :: max_gas_sat_change_cell
  PetscInt :: min_liq_pressure_cell
  PetscInt :: min_gas_pressure_cell
  PetscReal :: abs_dX_over_absX
  character(len=4) :: reason
  PetscInt :: istart, iend
  
  grid => this%realization%patch%grid
  option => this%realization%option
  field => this%realization%field
  patch => this%realization%patch
  wippflo_auxvars => patch%aux%WIPPFlo%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  material_parameter => patch%aux%Material%material_parameter

  dX_changed = PETSC_FALSE
  X1_changed = PETSC_FALSE
  
  call VecGetArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(X1,X1_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)
  call VecGetArrayReadF90(field%flow_accum2,accum_p2,ierr);CHKERRQ(ierr)
  converged_liquid_equation = PETSC_TRUE
  converged_gas_equation = PETSC_TRUE
  converged_liquid_pressure = PETSC_TRUE
  converged_gas_saturation = PETSC_TRUE
  max_liq_eq = 0.d0
  max_gas_eq = 0.d0
  max_liq_pres_rel_change = 0.d0
  max_gas_sat_change = 0.d0
  min_liq_pressure = 1.d20
  min_gas_pressure = 1.d20
  max_liq_eq_cell = 0
  max_gas_eq_cell = 0
  max_liq_pres_rel_change_cell = 0
  max_gas_sat_change_cell = 0
  min_liq_pressure_cell = 0
  min_gas_pressure_cell = 0
  do local_id = 1, grid%nlmax
    offset = (local_id-1)*option%nflowdof
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    liquid_equation_index = offset + WIPPFLO_LIQUID_EQUATION_INDEX
    gas_equation_index = offset + WIPPFLO_GAS_EQUATION_INDEX
    pressure_index = offset + WIPPFLO_LIQUID_PRESSURE_DOF
    saturation_index = offset + WIPPFLO_GAS_SATURATION_DOF
    
    ! liquid component equation
    if (accum_p2(liquid_equation_index) > zero_accumulation) then 
      residual = r_p(liquid_equation_index)
      accumulation = accum_p2(liquid_equation_index)
      residual_over_accumulation = residual / accumulation
      if (max_liq_eq < dabs(residual_over_accumulation)) then
        max_liq_eq_cell = local_id
        max_liq_eq = residual_over_accumulation
      endif
      if (dabs(residual) > this%liquid_equation_tolerance) then
        if (dabs(residual_over_accumulation) > &
            this%liquid_equation_tolerance) then
          converged_liquid_equation = PETSC_FALSE
        endif
      endif
    endif

    ! gas component equation
    if (accum_p2(gas_equation_index) > zero_accumulation .and. &
        X1_p(gas_equation_index) > zero_saturation) then
      residual = r_p(gas_equation_index)
      accumulation = accum_p2(gas_equation_index)
      residual_over_accumulation = residual / accumulation
      if (max_gas_eq < dabs(residual_over_accumulation)) then
        max_gas_eq_cell = local_id
        max_gas_eq = residual_over_accumulation
      endif
      if (dabs(residual) > this%gas_equation_tolerance) then 
        if (dabs(residual_over_accumulation) > &
            this%gas_equation_tolerance) then
          converged_gas_equation = PETSC_FALSE
        endif
      endif
    endif

    ! maximum relative change in liquid pressure
    abs_X = dabs(X1_p(pressure_index))
    if (abs_X > 0.d0) then
      abs_dX_over_absX = dabs(dX_p(pressure_index))/abs_X
      if (dabs(max_liq_pres_rel_change) < abs_dX_over_absX) then
        max_liq_pres_rel_change_cell = local_id
        max_liq_pres_rel_change = dX_p(pressure_index)/abs_X
      endif
      if (abs_dX_over_absX >= this%liquid_pressure_tolerance) then
        converged_liquid_pressure = PETSC_FALSE
      endif
    endif
    
    ! gas saturation
    abs_dX = dabs(dX_p(saturation_index))
    if (abs_dX > 0.d0) then
      ! BRAGFLO uses -log10(abs_dX), why not check if abs_dX > tol?
!      if (-log10(abs_dX) <= &
      if (dabs(max_gas_sat_change) < abs_dX) then
        max_gas_sat_change_cell = local_id
        max_gas_sat_change = dX_p(saturation_index)
      endif
      if (abs_dX > this%gas_saturation_tolerance) then
        converged_gas_saturation = PETSC_FALSE
      endif
    endif

    ! liquid pressure
    if (X1_p(pressure_index) < min_liq_pressure) then
      min_liq_pressure = X1_p(pressure_index)
      min_liq_pressure_cell = local_id
    endif

    ! gas pressure
    if (wippflo_auxvars(0,ghosted_id)%pres(2) < min_gas_pressure) then
      min_gas_pressure = wippflo_auxvars(0,ghosted_id)%pres(2)
      min_gas_pressure_cell = local_id
    endif
  enddo

  temp_int = 0
  if (.not.converged_liquid_equation) temp_int(1) = 1
  if (.not.converged_gas_equation) temp_int(2) = 1
  if (.not.converged_liquid_pressure) temp_int(3) = 1
  if (.not.converged_gas_saturation) temp_int(4) = 1
  if (min_gas_pressure < 0.d0) temp_int(5) = 1
  call MPI_Allreduce(MPI_IN_PLACE,temp_int,FIVE_INTEGER_MPI, &
                     MPIU_INTEGER,MPI_MAX,option%mycomm,ierr)
  temp_real(1) = max_liq_eq
  temp_real(2) = max_gas_eq
  temp_real(3) = max_liq_pres_rel_change
  temp_real(4) = max_gas_sat_change
  call MPI_Allreduce(MPI_IN_PLACE,temp_real,FOUR_INTEGER_MPI, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
  max_liq_eq = temp_real(1)
  max_gas_eq = temp_real(2)
  max_liq_pres_rel_change = temp_real(3)
  max_gas_sat_change = temp_real(4)

  if (option%convergence == CONVERGENCE_CUT_TIMESTEP) then
    reason = '!cut'
  else if (temp_int(5) > 0) then
    reason = 'negP'
    option%convergence = CONVERGENCE_CUT_TIMESTEP
  else if (option%convergence == CONVERGENCE_FORCE_ITERATION) then
    reason = '!it '
  else if (maxval(temp_int) > 0) then
    option%convergence = CONVERGENCE_KEEP_ITERATING
    reason = '    '
    if (temp_int(1) > 0) reason(1:1) = 'L'
    if (temp_int(2) > 0) reason(2:2) = 'G'
    if (temp_int(3) > 0) reason(3:3) = 'P'
    if (temp_int(4) > 0) reason(4:4) = 'S'
  else
    reason = '----'
    option%convergence = CONVERGENCE_CONVERGED
  endif
  if (OptionPrintToScreen(option)) then
    if (option%mycommsize > 1) then
      write(*,'(4x,"Reason: ",a4,4es10.2)') reason, &
        max_liq_eq, max_gas_eq, &
        max_liq_pres_rel_change, max_gas_sat_change
    else
      write(*,'(4x,"Reason: ",a4,4(i5,es10.2))') reason, &
        max_liq_eq_cell, max_liq_eq, &
        max_gas_eq_cell, max_gas_eq, &
        max_liq_pres_rel_change_cell, max_liq_pres_rel_change, &
        max_gas_sat_change_cell, max_gas_sat_change
! for debugging
#if 0
      local_id = min_gas_pressure_cell
      offset = (local_id-1)*option%nflowdof
      istart = offset + 1
      iend = offset + option%nflowdof
      ghosted_id = grid%nL2G(local_id)
      write(*,'(4x,i5,3es11.3)') local_id, &
        wippflo_auxvars(0,ghosted_id)%pres(1:2), &
        wippflo_auxvars(0,ghosted_id)%pres(option%capillary_pressure_id)
#endif
    endif
  endif

  call VecRestoreArrayReadF90(dX,dX_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X0,X0_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(X1,X1_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%flow_r,r_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%flow_accum,accum_p,ierr);CHKERRQ(ierr)
  call VecRestoreArrayReadF90(field%flow_accum2,accum_p2,ierr);CHKERRQ(ierr)
                               
end subroutine PMWIPPFloCheckUpdatePost

! ************************************************************************** !

subroutine PMWIPPFloTimeCut(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use WIPP_Flow_module, only : WIPPFloTimeCut

  implicit none
  
  class(pm_wippflo_type) :: this
  
  call PMSubsurfaceFlowTimeCut(this)
  call WIPPFloTimeCut(this%realization)

end subroutine PMWIPPFloTimeCut

! ************************************************************************** !

subroutine PMWIPPFloUpdateSolution(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use WIPP_Flow_module, only : WIPPFloUpdateSolution, &
                             WIPPFloMapBCAuxVarsToGlobal

  implicit none
  
  class(pm_wippflo_type) :: this
  
  call PMSubsurfaceFlowUpdateSolution(this)
  call WIPPFloUpdateSolution(this%realization)
  call WIPPFloMapBCAuxVarsToGlobal(this%realization)

end subroutine PMWIPPFloUpdateSolution     

! ************************************************************************** !

subroutine PMWIPPFloUpdateAuxVars(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  use WIPP_Flow_module, only : WIPPFloUpdateAuxVars

  implicit none
  
  class(pm_wippflo_type) :: this

  call WIPPFloUpdateAuxVars(this%realization)

end subroutine PMWIPPFloUpdateAuxVars   

! ************************************************************************** !

subroutine PMWIPPFloMaxChange(this)
  ! 
  ! Not needed given WIPPFloMaxChange is called in PostSolve
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use Realization_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use WIPP_Flow_Aux_module
  use Variables_module, only : LIQUID_PRESSURE, LIQUID_MOLE_FRACTION, &
                               TEMPERATURE, GAS_PRESSURE, AIR_PRESSURE, &
                               GAS_SATURATION
  implicit none
  
  class(pm_wippflo_type) :: this
  
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
  
  ! max change variables: [LIQUID_PRESSURE, GAS_SATURATION]
  do i = 1, 3
    call RealizationGetVariable(realization,field%work, &
                                this%max_change_ivar(i),ZERO_INTEGER)
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
    write(*,'("  --> max chng: dpl= ",1pe12.4, " dpg= ",1pe12.4, &
             &" dsg= ",1pe12.4)') &
      max_change_global(1:3)
  endif
  if (OptionPrintToFile(option)) then
    write(option%fid_out,'("  --> max chng: dpl= ",1pe12.4, " dpg= ", &
                          &1pe12.4, " dsg= ",1pe12.4)') &
      max_change_global(1:3)
  endif

  ! max change variables: [LIQUID_PRESSURE, GAS_PRESSURE, GAS_SATURATION]
  this%max_pressure_change = maxval(max_change_global(1:2))
  this%max_saturation_change = max_change_global(3)
  
end subroutine PMWIPPFloMaxChange

! ************************************************************************** !

subroutine PMWIPPFloComputeMassBalance(this,mass_balance_array)
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use WIPP_Flow_module, only : WIPPFloComputeMassBalance

  implicit none
  
  class(pm_wippflo_type) :: this
  PetscReal :: mass_balance_array(:)
  
  call WIPPFloComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMWIPPFloComputeMassBalance

! ************************************************************************** !

subroutine PMWIPPFloInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 07/11/17
  ! 
  
  implicit none
  
  class(pm_wippflo_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'wipp flow'

end subroutine PMWIPPFloInputRecord

! ************************************************************************** !

subroutine PMWIPPFloCheckpointBinary(this,viewer)
  ! 
  ! Checkpoints data associated with WIPPFlo PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17

  use Checkpoint_module
  use Global_module
  use Variables_module, only : STATE

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_wippflo_type) :: this
  PetscViewer :: viewer
  
  call PMSubsurfaceFlowCheckpointBinary(this,viewer)
  
end subroutine PMWIPPFloCheckpointBinary

! ************************************************************************** !

subroutine PMWIPPFloRestartBinary(this,viewer)
  ! 
  ! Restarts data associated with WIPPFlo PM
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17

  use Checkpoint_module
  use Global_module
  use Variables_module, only : STATE

  implicit none
#include "petsc/finclude/petscviewer.h"      

  class(pm_wippflo_type) :: this
  PetscViewer :: viewer
  
  call PMSubsurfaceFlowRestartBinary(this,viewer)
  
end subroutine PMWIPPFloRestartBinary

! ************************************************************************** !

subroutine PMWIPPFloDestroy(this)
  ! 
  ! Destroys WIPPFlo process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/11/17
  ! 

  use WIPP_Flow_module, only : WIPPFloDestroy

  implicit none
  
  class(pm_wippflo_type) :: this
  
  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  deallocate(this%max_change_ivar)
  nullify(this%max_change_ivar)

  ! preserve this ordering
  call WIPPFloDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)
  
end subroutine PMWIPPFloDestroy
  
end module PM_WIPP_Flow_class
