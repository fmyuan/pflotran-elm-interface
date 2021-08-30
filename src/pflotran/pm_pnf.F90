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
    PetscReal :: liq_pres_change_ts_governor
    Vec :: max_pressure_change_vec
  contains
    procedure, public :: ReadSimulationOptionsBlock => &
                           PMPNFReadSimOptionsBlock
    procedure, public :: ReadTSBlock => PMPNFReadTSSelectCase
    procedure, public :: InitializeRun => PMPNFInitializeRun
    procedure, public :: InitializeTimestep => PMPNFInitializeTimestep
    procedure, public :: UpdateTimestep => PMPNFUpdateTimestep
    procedure, public :: FinalizeTimestep => PMPNFFinalizeTimestep
    procedure, public :: PreSolve => PMPNFPreSolve
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
      case('NO_ACCUMULATION')
        pnf_calc_accum = PETSC_FALSE
      case('NO_FLUX')
        pnf_calc_flux = PETSC_FALSE
      case('NO_BCFLUX')
        pnf_calc_bcflux = PETSC_FALSE
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
    case('LIQ_PRES_CHANGE_TS_GOVERNOR')
      call InputReadDouble(input,option,this%liq_pres_change_ts_governor)
      call InputErrorMsg(input,option,keyword,error_string)
      ! units conversion since it is absolute
      call InputReadAndConvertUnits(input,this%liq_pres_change_ts_governor, &
                                    'Pa',keyword,option)
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

  implicit none

  class(pm_pnf_type) :: this

  call DiscretizationDuplicateVector(this%realization%discretization, &
                                     this%realization%field%work, &
                                     this%max_pressure_change_vec)

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

  use Upwind_Direction_module
  use Option_module

  implicit none

  class(pm_pnf_type) :: this

end subroutine PMPNFPostSolve

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
  pres_ratio = (2.d0*this%liq_pres_change_ts_governor)/ &
               (this%liq_pres_change_ts_governor+this%max_pressure_change)
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
    call OptionPrint(string,this%option)
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

  implicit none

  class(pm_pnf_type) :: this

  call PMSubsurfaceFlowTimeCut(this)
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
                     MPI_DOUBLE_PRECISION,MPI_MAX,option%mycomm,ierr)
  ! print them out
  if (OptionPrintToScreen(option)) then
    write(*,'("  --> max chng: dpl= ",1pe12.4)') &
      max_change_global(1)
  endif
  if (OptionPrintToFile(option)) then
    write(option%fid_out,'("  --> max chng: dpl= ",1pe12.4)') &
      max_change_global(1)
  endif

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
