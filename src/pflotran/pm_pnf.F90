module PM_PNF_class

#include "petsc/finclude/petscsnes.h"
  use petscsnes

  use PM_Subsurface_Flow_class
  use Communicator_Base_class
  use Realization_Subsurface_class
  use PFLOTRAN_Constants_module
  use PNF_Aux_module

  implicit none

  private

  type, public, extends(pm_subsurface_flow_type) :: pm_pnf_type
    Vec :: max_pressure_change_vec
    PetscBool :: use_darcy
  contains
    procedure, public :: ReadSimulationOptionsBlock => &
                           PMPNFReadSimOptionsBlock
    procedure, public :: ReadTSBlock => PMPNFReadTSSelectCase
    procedure, public :: InitializeRun => PMPNFInitializeRun
    procedure, public :: InitializeTimestep => PMPNFInitializeTimestep
    procedure, public :: UpdateTimestep => PMPNFUpdateTimestep
    procedure, public :: FinalizeTimestep => PMPNFFinalizeTimestep
    procedure, public :: PreSolve => PMPNFPreSolve
    procedure, public :: SetupLinearSystem => PMPNFSetupLinearSystem
    procedure, public :: PostSolve => PMPNFPostSolve
    procedure, public :: TimeCut => PMPNFTimeCut
    procedure, public :: UpdateSolution => PMPNFUpdateSolution
    procedure, public :: UpdateAuxVars => PMPNFUpdateAuxVars
    procedure, public :: MaxChange => PMPNFMaxChange
    procedure, public :: ComputeMassBalance => PMPNFComputeMassBalance
    procedure, public :: InputRecord => PMPNFInputRecord
    procedure, public :: CheckpointBinary => PMPNFCheckpointBinary
    procedure, public :: RestartBinary => PMPNFRestartBinary
    procedure, public :: Destroy => PMPNFDestroy
  end type pm_pnf_type

  public :: PMPNFCreate, &
            PMPNFInitObject, &
            PMPNFInitializeRun, &
            PMPNFFinalizeTimestep, &
            PMPNFDestroy

contains

! ************************************************************************** !

function PMPNFCreate()
  !
  ! Creates PNF process models shell
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !
  implicit none

  class(pm_pnf_type), pointer :: PMPNFCreate

  class(pm_pnf_type), pointer :: pnf_pm

  allocate(pnf_pm)
  call PMPNFInitObject(pnf_pm)

  PMPNFCreate => pnf_pm

end function PMPNFCreate

! ************************************************************************** !

subroutine PMPNFInitObject(this)
  !
  ! Creates PNF process models shell
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !
  use EOS_Water_module, only : EOSWaterSetDensity

  implicit none

  class(pm_pnf_type) :: this

  PetscReal :: array(1)

  call PMSubsurfaceFlowInit(this)
  this%name = 'PN Flow'
  this%header = 'PN FLOW'

  this%max_pressure_change_vec = PETSC_NULL_VEC
  this%use_darcy = PETSC_FALSE

  array(1) = pnf_density_kg ! dist is the aux array
  call EOSWaterSetDensity('CONSTANT',array)

end subroutine PMPNFInitObject

! ************************************************************************** !

subroutine PMPNFReadSimOptionsBlock(this,input)
  !
  ! Read PNF options input block
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !
  use PNF_module
  use PNF_Aux_module
  use Input_Aux_module
  use String_module
  use Option_module
  use Utility_module

  implicit none

  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: keyword
  class(pm_pnf_type) :: this
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  option => this%option

  error_string = 'PNF Options'

  input%ierr = 0
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    found = PETSC_FALSE
    call PMSubsurfFlowReadSimOptionsSC(this,input,keyword,found, &
                                       error_string,option)
    if (found) cycle

    select case(trim(keyword))
      case('USE_DARCY')
        this%use_darcy = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(input,keyword,'PNF Mode',option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine PMPNFReadSimOptionsBlock

! ************************************************************************** !

subroutine PMPNFReadTSSelectCase(this,input,keyword,found, &
                                   error_string,option)
  !
  ! Read timestepper settings specific to the PNF process model
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21

  use Input_Aux_module
  use String_module
  use Option_module

  implicit none

  class(pm_pnf_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option

  found = PETSC_TRUE
  call PMSubsurfaceFlowReadTSSelectCase(this,input,keyword,found, &
                                        error_string,option)
  if (found) return

  found = PETSC_TRUE
  select case(trim(keyword))
    case default
      found = PETSC_FALSE
  end select

end subroutine PMPNFReadTSSelectCase

! ************************************************************************** !

recursive subroutine PMPNFInitializeRun(this)
  !
  ! Initializes the PNF mode run.
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21

  use Discretization_module
  use Realization_Base_class
  use Variables_module, only : LIQUID_PRESSURE

  implicit none

  class(pm_pnf_type) :: this

  call DiscretizationDuplicateVector(this%realization%discretization, &
                                     this%realization%field%work, &
                                     this%max_pressure_change_vec)
  call RealizationGetVariable(this%realization, &
                              this%max_pressure_change_vec, &
                              LIQUID_PRESSURE,ZERO_INTEGER)

  ! call parent implementation
  call PMSubsurfaceFlowInitializeRun(this)

end subroutine PMPNFInitializeRun

! ************************************************************************** !

subroutine PMPNFInitializeTimestep(this)
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !

  use PNF_module, only : PNFInitializeTimestep
  use PNF_Aux_module
  use Global_module
  use Option_module

  implicit none

  class(pm_pnf_type) :: this

  call PMSubsurfaceFlowInitializeTimestepA(this)
  call PNFInitializeTimestep(this%realization)
  call PMSubsurfaceFlowInitializeTimestepB(this)

end subroutine PMPNFInitializeTimestep

! ************************************************************************** !

subroutine PMPNFFinalizeTimestep(this)
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !
  implicit none

  class(pm_pnf_type) :: this

  call PMSubsurfaceFlowFinalizeTimestep(this)

end subroutine PMPNFFinalizeTimestep

! ************************************************************************** !

subroutine PMPNFPreSolve(this)
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21

  implicit none

  class(pm_pnf_type) :: this

  call PMSubsurfaceFlowPreSolve(this)

end subroutine PMPNFPreSolve

! ************************************************************************** !

subroutine PMPNFPostSolve(this)
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21

  implicit none

  class(pm_pnf_type) :: this

end subroutine PMPNFPostSolve

! ************************************************************************** !

subroutine PMPNFSetupLinearSystem(this,A,solution,right_hand_side,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21

  use Grid_module
  use Patch_module
  use Material_Aux_module
  use Coupler_module
  use Connection_module
  use Option_module

  implicit none

  class(pm_pnf_type) :: this
  Vec :: right_hand_side
  Vec :: solution
  Mat :: A
  PetscErrorCode :: ierr

  PetscViewer :: viewer

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(connection_set_list_type), pointer :: connection_set_list

  PetscInt :: local_id, local_id_up, local_id_dn
  PetscInt :: ghosted_id, ghosted_id_up, ghosted_id_dn
  PetscInt :: iconn, sum_connection
  PetscInt :: itype
  PetscReal :: area
  PetscReal :: rate
  PetscReal :: rvalue
  PetscReal :: tempreal(1,1)
  PetscReal, parameter :: g_sup_h_constant = 0.4217d0 &
!                          * pnf_density_kg * EARTH_GRAVITY &
                          / (12.d0 * pnf_viscosity)
  PetscReal, pointer :: rhs_ptr(:)

  solution = this%realization%field%flow_xx
  right_hand_side = this%realization%field%flow_rhs

  option => this%realization%option
  patch => this%realization%patch
  grid => patch%grid

  material_auxvars => patch%aux%Material%auxvars

  call MatZeroEntries(A,ierr);CHKERRQ(ierr)
  call VecZeroEntries(right_hand_side,ierr);CHKERRQ(ierr)

  call VecGetArrayF90(right_hand_side,rhs_ptr,ierr);CHKERRQ(ierr)
#if 0
  do local_id = 1, grid%nlmax  ! For each local node do...
    ghosted_id = grid%nL2G(local_id)
    tempreal = pnf_density_kg * material_auxvars(ghosted_id)%volume / &
               option%flow_dt
    call MatSetValuesLocal(A,1,ghosted_id-1,1,ghosted_id-1,tempreal, &
                           ADD_VALUES,ierr);CHKERRQ(ierr)
    rhs_ptr(local_id) = tempreal(1,1)
  enddo
#endif

  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0
  do
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      local_id_up = grid%nG2L(ghosted_id_up)
      local_id_dn = grid%nG2L(ghosted_id_dn)
      if (this%use_darcy) then
        tempreal = pnf_density_kg * &
                  material_auxvars(ghosted_id_up)%permeability(1) * &
                  cur_connection_set%area(iconn) / &
                  (pnf_viscosity * cur_connection_set%dist(0,iconn))
      else
        tempreal = g_sup_h_constant * cur_connection_set%area(iconn)**2 / &
                  cur_connection_set%dist(0,iconn)
      endif
      if (local_id_up > 0) then
        call MatSetValuesLocal(A,1,ghosted_id_up-1,1,ghosted_id_up-1,tempreal, &
                               ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesLocal(A,1,ghosted_id_up-1,1,ghosted_id_dn-1, &
                               -tempreal,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
      if (local_id_dn > 0) then
        call MatSetValuesLocal(A,1,ghosted_id_dn-1,1,ghosted_id_dn-1,tempreal, &
                               ADD_VALUES,ierr);CHKERRQ(ierr)
        call MatSetValuesLocal(A,1,ghosted_id_dn-1,1,ghosted_id_up-1, &
                               -tempreal,ADD_VALUES,ierr);CHKERRQ(ierr)
      endif
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    itype = boundary_condition%flow_condition%itype(PNF_LIQUID_PRESSURE_DOF)
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      area = cur_connection_set%area(iconn)
      rvalue = boundary_condition%flow_aux_real_var(1,iconn)
      select case(itype)
        case(DIRICHLET_BC)
          if (this%use_darcy) then
            tempreal = pnf_density_kg * &
                      material_auxvars(ghosted_id)%permeability(1) * &
                      cur_connection_set%area(iconn) / &
                      (pnf_viscosity * cur_connection_set%dist(0,iconn))
          else
            tempreal = g_sup_h_constant * area**2 / &  ! w^3*b
                      cur_connection_set%dist(0,iconn)
          endif
          call MatSetValuesLocal(A,1,ghosted_id-1,1,ghosted_id-1,tempreal, &
                                 ADD_VALUES,ierr);CHKERRQ(ierr)
          tempreal = tempreal * rvalue
        case(NEUMANN_BC)
          if (this%use_darcy) then
            tempreal = rvalue * area * pnf_density_kg
          else
            tempreal = rvalue * area
          endif
      end select
      rhs_ptr(local_id) = tempreal(1,1)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit

    cur_connection_set => source_sink%connection_set
    itype = source_sink%flow_condition%rate%itype
    rate = source_sink%flow_condition%rate%dataset%rarray(1)
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      tempreal = source_sink%flow_condition%rate%dataset%rarray(1)
      select case(itype)
      case(MASS_RATE_SS)
        if (this%use_darcy) then
          tempreal = rate
        else
          tempreal = rate / pnf_density_kg
        endif
      case(VOLUMETRIC_RATE_SS)
        if (this%use_darcy) then
          tempreal = rate * pnf_density_kg
        else
          tempreal = rate
        endif
        case default
          option%io_buffer = 'src_sink_type not supported in PNFSrcSink'
          call PrintErrMsg(option)
      end select
      rhs_ptr(local_id) = tempreal(1,1)
    enddo
    source_sink => source_sink%next
  enddo

  call VecRestoreArrayF90(right_hand_side,rhs_ptr,ierr);CHKERRQ(ierr)
  call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
  call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

  if (this%realization%debug%matview_Matrix) then
    call PetscViewerASCIIOpen(option%mycomm,'PNFmatrix.mat',viewer, &
                              ierr);CHKERRQ(ierr)
    call MatView(A,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (this%realization%debug%vecview_residual) then
    call PetscViewerASCIIOpen(option%mycomm,'PNFrhs.vec',viewer, &
                              ierr);CHKERRQ(ierr)
    call VecView(right_hand_side,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif
  if (this%realization%debug%vecview_solution) then
    call PetscViewerASCIIOpen(option%mycomm,'PNFsolution.vec',viewer, &
                              ierr);CHKERRQ(ierr)
    call VecView(solution,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
  endif

end subroutine PMPNFSetupLinearSystem

! ************************************************************************** !

subroutine PMPNFCalculateVelocities(this)
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21

  use Option_module
  use Grid_module
  use Patch_module
  use Coupler_module
  use Connection_module
  use Material_Aux_module

  implicit none

  class(pm_pnf_type) :: this

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(connection_set_list_type), pointer :: connection_set_list
  PetscReal, pointer :: vec_loc_ptr(:)

  PetscInt :: local_id, local_id_up, local_id_dn
  PetscInt :: ghosted_id, ghosted_id_up, ghosted_id_dn
  PetscInt :: iconn, sum_connection
  PetscInt :: itype
  PetscReal :: area
  PetscReal :: rvalue
  PetscReal :: velocity
  PetscReal :: tempreal
  PetscReal, parameter :: g_sup_h_constant = 0.4217d0 &
!                          * pnf_density_kg * EARTH_GRAVITY &
                          / (12.d0 * pnf_viscosity)
  PetscErrorCode :: ierr
  option => this%realization%option
  patch => this%realization%patch
  grid => patch%grid

  material_auxvars => patch%aux%Material%auxvars

  call VecGetArrayReadF90(this%realization%field%flow_xx_loc,vec_loc_ptr, &
                          ierr);CHKERRQ(ierr)

  ! Interior Flux Terms -----------------------------------
  connection_set_list => grid%internal_connection_set_list
  cur_connection_set => connection_set_list%first
  sum_connection = 0
  do
    if (.not.associated(cur_connection_set)) exit
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      ghosted_id_up = cur_connection_set%id_up(iconn)
      ghosted_id_dn = cur_connection_set%id_dn(iconn)
      local_id_up = grid%nG2L(ghosted_id_up)
      local_id_dn = grid%nG2L(ghosted_id_dn)
      if (this%use_darcy) then
        tempreal = material_auxvars(ghosted_id_up)%permeability(1) / &
                  (pnf_viscosity * cur_connection_set%dist(0,iconn))
      else
        tempreal = g_sup_h_constant * cur_connection_set%area(iconn) / & ! **2 -> **1
                  cur_connection_set%dist(0,iconn)
      endif
      velocity = tempreal * &
                 (vec_loc_ptr(ghosted_id_up) - vec_loc_ptr(ghosted_id_dn))
      patch%internal_velocities(1,sum_connection) = velocity
    enddo
    cur_connection_set => cur_connection_set%next
  enddo

  ! Boundary Flux Terms -----------------------------------
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    itype = boundary_condition%flow_condition%itype(PNF_LIQUID_PRESSURE_DOF)
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      area = cur_connection_set%area(iconn)
      rvalue = boundary_condition%flow_aux_real_var(1,iconn)
      select case(itype)
        case(DIRICHLET_BC)
          if (this%use_darcy) then
            tempreal = material_auxvars(ghosted_id)%permeability(1) / &
                      (pnf_viscosity * cur_connection_set%dist(0,iconn))
          else
            tempreal = g_sup_h_constant * area / &   ! **2 -> **1
                      cur_connection_set%dist(0,iconn)
          endif
          velocity = tempreal * (rvalue - vec_loc_ptr(ghosted_id))
        case(NEUMANN_BC)
          velocity = rvalue
      end select
      patch%boundary_velocities(1,sum_connection) = velocity
    enddo
    boundary_condition => boundary_condition%next
  enddo

  call VecRestoreArrayReadF90(this%realization%field%flow_xx_loc,vec_loc_ptr, &
                              ierr);CHKERRQ(ierr)

end subroutine PMPNFCalculateVelocities

! ************************************************************************** !

subroutine PMPNFUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                                 num_newton_iterations,tfac, &
                                 time_step_max_growth_factor)
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !
  use Option_module
  use Realization_Base_class, only : RealizationGetVariable
  use Realization_Subsurface_class, only : RealizationLimitDTByCFL
  use Field_module
  use Global_module, only : GlobalSetAuxVarVecLoc
  use Utility_module, only : Equal

  implicit none

  class(pm_pnf_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min ! DO NOT USE (see comment below)
  PetscReal :: dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor

  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: pres_ratio
  PetscReal :: dt_prev

  dt_prev = dt

  ! calculate the time step ramping factor
  pres_ratio = (2.d0*this%pressure_change_governor)/ &
               (this%pressure_change_governor+this%max_pressure_change)
  ! pick minimum time step from calc'd ramping factor or maximum ramping factor
  dt = min(pres_ratio*dt,time_step_max_growth_factor*dt)
  ! make sure time step is within bounds given in the input deck
  dt = min(dt,dt_max)
  if (this%logging_verbosity > 0) then
    if (Equal(dt,dt_max)) then
      string = 'maximum time step size'
    else if (pres_ratio > time_step_max_growth_factor) then
      string = 'maximum time step growth factor'
    else
      string = 'liquid pressure governor'
    endif
    string = 'TS update: ' // trim(string)
    call PrintMsg(this%option,string)
  endif

  if (Initialized(this%cfl_governor)) then
    call RealizationLimitDTByCFL(this%realization,this%cfl_governor,dt,dt_max)
  endif

end subroutine PMPNFUpdateTimestep

! ************************************************************************** !

subroutine PMPNFTimeCut(this)
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !

  use PNF_module, only : PNFTimeCut
  use Discretization_module, only : DiscretizationGlobalToLocal

  implicit none

  class(pm_pnf_type) :: this

  call PMSubsurfaceFlowTimeCut(this)
  call DiscretizationGlobalToLocal(this%realization%discretization, &
                                   this%realization%field%flow_xx, &
                                   this%realization%field%flow_xx_loc,NFLOWDOF)
  call PNFTimeCut(this%realization)

end subroutine PMPNFTimeCut

! ************************************************************************** !

subroutine PMPNFUpdateSolution(this)
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !

  use PNF_module, only : PNFUpdateSolution

  implicit none

  class(pm_pnf_type) :: this

  call PMSubsurfaceFlowUpdateSolution(this)
  call PNFUpdateSolution(this%realization)
  call PMPNFCalculateVelocities(this)
  !call PNFMapBCAuxVarsToGlobal(this%realization)

end subroutine PMPNFUpdateSolution

! ************************************************************************** !

subroutine PMPNFUpdateAuxVars(this)
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  use PNF_module, only : PNFUpdateAuxVars

  implicit none

  class(pm_pnf_type) :: this

  call PNFUpdateAuxVars(this%realization)

end subroutine PMPNFUpdateAuxVars

! ************************************************************************** !

subroutine PMPNFMaxChange(this)
  !
  ! Not needed given PNFMaxChange is called in PostSolve
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !

  use Realization_Base_class
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Grid_module
  use PNF_Aux_module
  use Variables_module, only : LIQUID_PRESSURE, LIQUID_SATURATION

  implicit none

  class(pm_pnf_type) :: this

  class(realization_subsurface_type), pointer :: realization
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  PetscReal, pointer :: vec_old_ptr(:), vec_new_ptr(:)
  PetscReal :: max_change_local(1)
  PetscReal :: max_change_global(1)
  PetscReal :: max_change, change
  PetscInt :: j

  PetscErrorCode :: ierr

  realization => this%realization
  option => realization%option
  field => realization%field
  grid => realization%patch%grid

  max_change_global = 0.d0
  max_change_local = 0.d0

  call RealizationGetVariable(realization,field%work, &
                              LIQUID_PRESSURE,ZERO_INTEGER)
  ! yes, we could use VecWAXPY and a norm here, but we need the ability
  ! to customize
  call VecGetArrayF90(field%work,vec_new_ptr,ierr);CHKERRQ(ierr)
  call VecGetArrayF90(this%max_pressure_change_vec,vec_old_ptr, &
                      ierr);CHKERRQ(ierr)
  max_change = 0.d0
  do j = 1, grid%nlmax
    change = dabs(vec_new_ptr(j)-vec_old_ptr(j))
    max_change = max(max_change,change)
  enddo
  max_change_local(1) = max_change
  call VecRestoreArrayF90(field%work,vec_new_ptr,ierr);CHKERRQ(ierr)
  call VecRestoreArrayF90(this%max_pressure_change_vec,vec_old_ptr, &
                          ierr);CHKERRQ(ierr)
  call VecCopy(field%work,this%max_pressure_change_vec,ierr);CHKERRQ(ierr)
  call MPI_Allreduce(max_change_local,max_change_global,ONE_INTEGER, &
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm, &
                     ierr);CHKERRQ(ierr)
  ! print them out
  write(option%io_buffer,'("  --> max change: dpl= ",1pe12.4)') &
        max_change_global(1)
  call PrintMsg(option)

  this%max_pressure_change = max_change_global(1)

end subroutine PMPNFMaxChange

! ************************************************************************** !

subroutine PMPNFComputeMassBalance(this,mass_balance_array)
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !

  use PNF_module, only : PNFComputeMassBalance

  implicit none

  class(pm_pnf_type) :: this
  PetscReal :: mass_balance_array(:)

  call PNFComputeMassBalance(this%realization,mass_balance_array)

end subroutine PMPNFComputeMassBalance

! ************************************************************************** !

subroutine PMPNFInputRecord(this)
  !
  ! Writes ingested information to the input record file.
  !
  ! Author: Jenn Frederick, SNL
  ! Date: 08/27/21
  !

  implicit none

  class(pm_pnf_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name
  write(id,'(a29)',advance='no') 'mode: '
  write(id,'(a)') 'pnf'

end subroutine PMPNFInputRecord

! ************************************************************************** !

subroutine PMPNFCheckpointBinary(this,viewer)
  !
  ! Checkpoints data associated with PNF PM
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21

  use Checkpoint_module
  use Global_module
  use Variables_module, only : STATE

  implicit none
#include "petsc/finclude/petscviewer.h"

  class(pm_pnf_type) :: this
  PetscViewer :: viewer

  call PMSubsurfaceFlowCheckpointBinary(this,viewer)

end subroutine PMPNFCheckpointBinary

! ************************************************************************** !

subroutine PMPNFRestartBinary(this,viewer)
  !
  ! Restarts data associated with PNF PM
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21

  use Checkpoint_module
  use Global_module
  use Variables_module, only : STATE

  implicit none
#include "petsc/finclude/petscviewer.h"

  class(pm_pnf_type) :: this
  PetscViewer :: viewer

  call PMSubsurfaceFlowRestartBinary(this,viewer)

end subroutine PMPNFRestartBinary

! ************************************************************************** !

subroutine PMPNFDestroy(this)
  !
  ! Destroys PNF process model
  !
  ! Author: Glenn Hammond
  ! Date: 08/27/21
  !

  use PNF_module, only : PNFDestroy
  use Utility_module, only : DeallocateArray

  implicit none

  class(pm_pnf_type) :: this

  PetscErrorCode :: ierr

  if (associated(this%next)) then
    call this%next%Destroy()
  endif

  call VecDestroy(this%max_pressure_change_vec,ierr);CHKERRQ(ierr)
  call PNFDestroy(this%realization)
  call PMSubsurfaceFlowDestroy(this)

end subroutine PMPNFDestroy

end module PM_PNF_class
