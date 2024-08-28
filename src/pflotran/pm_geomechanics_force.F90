module PM_Geomechanics_Force_class

#include "petsc/finclude/petscts.h"
  use petscts
  use PM_Base_class
  use Geomechanics_Realization_class
  use Realization_Subsurface_class
  use Communicator_Base_class
  use Option_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_base_type) :: pm_geomech_force_type
    class(realization_geomech_type), pointer :: geomech_realization
    class(realization_subsurface_type), pointer :: subsurf_realization
    class(communicator_type), pointer :: comm1
  contains
    procedure, public :: ReadSimulationOptionsBlock => &
                           PMGeomechReadSimOptionsBlock
    procedure, public :: Setup => PMGeomechForceSetup
    procedure, public :: PMGeomechForceSetRealization
    procedure, public :: InitializeRun => PMGeomechForceInitializeRun
    procedure, public :: FinalizeRun => PMGeomechForceFinalizeRun
    procedure, public :: InitializeTimestep => PMGeomechForceInitializeTimestep
    procedure, public :: CheckConvergence => PMGeomechCheckConvergence
    procedure, public :: AcceptSolution => PMGeomechAcceptSolution
    procedure, public :: Residual => PMGeomechForceResidual
    procedure, public :: Jacobian => PMGeomechForceJacobian
    procedure, public :: PreSolve => PMGeomechForcePreSolve
    procedure, public :: UpdateSolution => PMGeomechForceUpdateSolution
    procedure, public :: CheckpointBinary => PMGeomechForceCheckpointBinary
    procedure, public :: RestartBinary => PMGeomechForceRestartBinary
    procedure, public :: InputRecord => PMGeomechForceInputRecord
    procedure, public :: Destroy => PMGeomechForceDestroy
    procedure, public :: FinalizeTimestep => PMGeomechForceFinalizeTimestep
  end type pm_geomech_force_type

  public :: PMGeomechForceCreate

contains

! ************************************************************************** !

function PMGeomechForceCreate()
  !
  ! This routine creates
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  !

  implicit none

  class(pm_geomech_force_type), pointer :: PMGeomechForceCreate

  class(pm_geomech_force_type), pointer :: geomech_force_pm

  allocate(geomech_force_pm)
  nullify(geomech_force_pm%option)
  nullify(geomech_force_pm%output_option)
  nullify(geomech_force_pm%geomech_realization)
  nullify(geomech_force_pm%subsurf_realization)
  nullify(geomech_force_pm%comm1)

  call PMBaseInit(geomech_force_pm)
  geomech_force_pm%header = 'GEOMECHANICS'

  PMGeomechForceCreate => geomech_force_pm

end function PMGeomechForceCreate

! ************************************************************************** !

subroutine PMGeomechReadSimOptionsBlock(this,input)
  !
  ! Reads input file parameters associated with the geomechanics process model
  !
  ! Author: Glenn Hammond
  ! Date: 03/16/20

  use Input_Aux_module
  use String_module
  use Utility_module
  use Option_module

  implicit none

  class(pm_geomech_force_type) :: this
  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: keyword, word
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option

  option => this%option

  error_string = 'SIMULATION,PROCESS_MODELS,SUBSURFACE_GEOMECHANICS,OPTIONS'

  input%ierr = 0
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)
    select case(trim(keyword))
      case('COUPLING')
        call InputReadCard(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case (word)
          case ('ONE_WAY_COUPLED')
            option%geomech_subsurf_coupling = GEOMECH_ONE_WAY_COUPLED
          case ('TWO_WAY_COUPLED')
            option%geomech_subsurf_coupling = GEOMECH_TWO_WAY_COUPLED
          case ('COUPLE_ERT')
            option%geomech_subsurf_coupling = GEOMECH_ERT_COUPLING
          case default
            call InputKeywordUnrecognized(input,word, &
                                trim(error_string)//',COUPLING',option)
        end select
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine PMGeomechReadSimOptionsBlock

! ************************************************************************** !

subroutine PMGeomechForceSetup(this)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  !
  use Geomechanics_Discretization_module
  use Communicator_Structured_class
  use Communicator_Unstructured_class
  use Grid_module

  implicit none

  class(pm_geomech_force_type) :: this

  ! set up communicator
  select case(this%geomech_realization%geomech_discretization%itype)
    case(STRUCTURED_GRID)
      this%comm1 => StructuredCommunicatorCreate()
    case(UNSTRUCTURED_GRID)
      this%comm1 => UnstructuredCommunicatorCreate()
  end select

  !call this%comm1%SetDM(this%geomech_realization%geomech_discretization%dm_1dof)

end subroutine PMGeomechForceSetup

! ************************************************************************** !

recursive subroutine PMGeomechForceInitializeRun(this)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  !

  use Geomechanics_Force_module, only : GeomechUpdateSolution

  implicit none

  class(pm_geomech_force_type) :: this

end subroutine PMGeomechForceInitializeRun

! ************************************************************************** !

recursive subroutine PMGeomechForceFinalizeRun(this)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  !

  implicit none

  class(pm_geomech_force_type) :: this

#ifdef PM_GEOMECH_FORCE_DEBUG
  call PrintMsg(this%option,'PMGeomechForce%FinalizeRun()')
#endif

  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif

end subroutine PMGeomechForceFinalizeRun

! ************************************************************************** !

subroutine PMGeomechForceSetRealization(this, geomech_realization)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  !

  use Grid_module

  implicit none

  class(pm_geomech_force_type) :: this
  class(realization_geomech_type), pointer :: geomech_realization

  this%geomech_realization => geomech_realization
  this%realization_base => geomech_realization

  this%solution_vec = geomech_realization%geomech_field%disp_xx
  this%residual_vec = geomech_realization%geomech_field%disp_r

end subroutine PMGeomechForceSetRealization

! ************************************************************************** !

subroutine PMGeomechForceInitializeTimestep(this)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  !

  use Geomechanics_Force_module, only : GeomechanicsForceInitialGuess
  use Global_module

  implicit none

  class(pm_geomech_force_type) :: this

#ifdef PM_GEOMECH_FORCE_DEBUG
  call PrintMsg(this%option,'PMGeomechForce%InitializeTimestep()')
#endif

  call GeomechanicsForceInitialGuess(this%geomech_realization)

end subroutine PMGeomechForceInitializeTimestep

! ************************************************************************** !

subroutine PMGeomechForceResidual(this,snes,xx,r,ierr)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  !

  use Geomechanics_Force_module, only : GeomechForceResidual

  implicit none

  class(pm_geomech_force_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr

#ifdef PM_GEOMECH_FORCE_DEBUG
  call PrintMsg(this%option,'PMGeomechForce%Residual()')
#endif

  call GeomechForceResidual(snes,xx,r,this%geomech_realization,ierr)

end subroutine PMGeomechForceResidual

! ************************************************************************** !

subroutine PMGeomechForceJacobian(this,snes,xx,A,B,ierr)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  !

  use Geomechanics_Force_module, only : GeomechForceJacobian

  implicit none

  class(pm_geomech_force_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr

#ifdef PM_GEOMECH_FORCE_DEBUG
  call PrintMsg(this%option,'PMGeomechForce%Jacobian()')
#endif

  call GeomechForceJacobian(snes,xx,A,B,this%geomech_realization,ierr)

end subroutine PMGeomechForceJacobian

! ************************************************************************** !

subroutine PMGeomechForcePreSolve(this)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  !

  implicit none

  class(pm_geomech_force_type) :: this

end subroutine PMGeomechForcePreSolve

! ************************************************************************** !

subroutine PMGeomechForceUpdateSolution(this)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  !

  use Geomechanics_Force_module, only : GeomechUpdateSolution, &
                                        GeomechStoreInitialDisp, &
                                        GeomechForceUpdateAuxVars
  use Geomechanics_Condition_module
  use Geomechanics_Realization_class, only : &
                                 GeomechRealizUpdateAllCouplerAuxVars


  implicit none

  class(pm_geomech_force_type) :: this

  PetscBool :: force_update_flag = PETSC_FALSE

#ifdef PM_GEOMECH_FORCE_DEBUG
  call PrintMsg(this%option,'PMGeomechForce%UpdateSolution()')
#endif

  ! begin from RealizationUpdate()
  call GeomechConditionUpdate(this%geomech_realization%geomech_conditions, &
                              this%geomech_realization%option)

  call GeomechUpdateSolution(this%geomech_realization)
  call GeomechRealizUpdateAllCouplerAuxVars(this%geomech_realization, &
                                            force_update_flag)
  if (this%option%geomech_initial) then
    call GeomechStoreInitialDisp(this%geomech_realization)
    this%option%geomech_initial = PETSC_FALSE
  endif
  call GeomechForceUpdateAuxVars(this%geomech_realization)

end subroutine PMGeomechForceUpdateSolution

! ************************************************************************** !

function PMGeomechAcceptSolution(this)
  !
  ! Author: Glenn Hammond
  ! Date: 03/19/21

  implicit none

  class(pm_geomech_force_type) :: this

  PetscBool :: PMGeomechAcceptSolution

  ! do nothing
  PMGeomechAcceptSolution = PETSC_TRUE

end function PMGeomechAcceptSolution

! ************************************************************************** !

subroutine PMGeomechCheckConvergence(this,snes,it,xnorm,unorm,fnorm, &
                                     reason,ierr)
  !
  ! Author: Glenn Hammond
  ! Date: 11/15/17
  !
  use Convergence_module
  use Grid_module

  implicit none

  class(pm_geomech_force_type) :: this
  SNES :: snes
  PetscInt :: it
  PetscReal :: xnorm
  PetscReal :: unorm
  PetscReal :: fnorm
  SNESConvergedReason :: reason
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid

  nullify(grid)

  call ConvergenceTest(snes,it,xnorm,unorm,fnorm,reason, &
                       grid,this%option,this%solver,ierr)

end subroutine PMGeomechCheckConvergence

! ************************************************************************** !

subroutine PMGeomechForceFinalizeTimestep(this)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  !

  use Global_module

  implicit none

  class(pm_geomech_force_type) :: this

#ifdef PM_GEOMECH_FORCE_DEBUG
  call PrintMsg(this%option,'PMGeomechForce%FinalizeTimestep()')
#endif

end subroutine PMGeomechForceFinalizeTimestep

! ************************************************************************** !

subroutine PMGeomechForceCheckpointBinary(this,viewer)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  !

  use Checkpoint_module

  implicit none
#include "petsc/finclude/petscviewer.h"

  class(pm_geomech_force_type) :: this
  PetscViewer :: viewer

  call PrintErrMsg(this%option,'add code for checkpointing Geomech in PM approach')

end subroutine PMGeomechForceCheckpointBinary

! ************************************************************************** !

subroutine PMGeomechForceRestartBinary(this,viewer)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  !

  use Checkpoint_module

  implicit none
#include "petsc/finclude/petscviewer.h"

  class(pm_geomech_force_type) :: this
  PetscViewer :: viewer

  call PrintErrMsg(this%option,'add code for restarting Geomech in PM approach')

end subroutine PMGeomechForceRestartBinary

! ************************************************************************** !

subroutine PMGeomechForceInputRecord(this)
  !
  ! Writes ingested information to the input record file.
  !
  ! Author: Jenn Frederick, SNL
  ! Date: 03/21/2016
  !

  implicit none

  class(pm_geomech_force_type) :: this

  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

end subroutine PMGeomechForceInputRecord

! ************************************************************************** !

subroutine PMGeomechForceDestroy(this)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  !

  use Geomechanics_Realization_class, only : GeomechRealizDestroy

  implicit none

  class(pm_geomech_force_type) :: this

  if (associated(this%next)) then
    call this%next%Destroy()
  endif
  call PMBaseDestroy(this)

#ifdef PM_GEOMECH_FORCE_DEBUG
  call PrintMsg(this%option,'PMGeomechForce%Destroy()')
#endif

  call GeomechRealizDestroy(this%geomech_realization)
  nullify(this%subsurf_realization)

  call this%comm1%Destroy()

end subroutine PMGeomechForceDestroy

end module PM_Geomechanics_Force_class
