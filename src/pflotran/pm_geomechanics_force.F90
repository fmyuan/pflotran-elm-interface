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
    procedure, public :: AcceptSolution => PMGeomechAcceptSolution
    procedure, public :: SetupLinearSystem => PMGeomechForceSetupLinearSystem
    procedure, public :: PreSolve => PMGeomechForcePreSolve
    procedure, public :: PostSolve => PMGeomechForcePostSolve
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
  PetscBool :: found

  option => this%option

  error_string = 'SIMULATION,PROCESS_MODELS,SUBSURFACE_GEOMECHANICS,OPTIONS'

  input%ierr = INPUT_ERROR_NONE
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)

    found = PETSC_TRUE
    call PMBaseReadSimOptionsSelectCase(this,input,keyword,found, &
                                        error_string,option)
    if (found) cycle

    select case(trim(keyword))
      case('FLOW_COUPLING')
        call InputReadCard(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case (word)
          case ('ONE_WAY_COUPLED')
            option%geomechanics%flow_coupling = GEOMECH_ONE_WAY_COUPLED
          case ('TWO_WAY_COUPLED')
            option%geomechanics%flow_coupling = GEOMECH_TWO_WAY_COUPLED
          case default
            call InputKeywordUnrecognized(input,word, &
                                trim(error_string)//','//trim(keyword),option)
        end select
      case('GEOPHYSICS_COUPLING')
        call InputReadCard(input,option,word,PETSC_FALSE)
        call StringToUpper(word)
        select case (word)
          case ('COUPLE_ERT')
            option%geomechanics%geophysics_coupling = GEOMECH_ERT_COUPLING
          case default
            call InputKeywordUnrecognized(input,word, &
                                trim(error_string)//','//trim(keyword),option)
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

  use Geomechanics_Force_module

  implicit none

  class(pm_geomech_force_type) :: this

  call GeomechStoreInitialPressTemp(this%geomech_realization)

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

subroutine PMGeomechForceSetRealization(this, geomech_realization, &
                                        subsurf_realization)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  !

  use Grid_module
  use Realization_Subsurface_class

  implicit none

  class(pm_geomech_force_type) :: this
  class(realization_geomech_type), pointer :: geomech_realization
  class(realization_subsurface_type), pointer :: subsurf_realization

  this%geomech_realization => geomech_realization
  this%realization_base => subsurf_realization

  this%solution_vec = geomech_realization%geomech_field%disp_xx

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

  ! call GeomechanicsForceInitialGuess(this%geomech_realization)

end subroutine PMGeomechForceInitializeTimestep

! ************************************************************************** !

subroutine PMGeomechForceSetupLinearSystem(this,A,solution,right_hand_side,ierr)
  !
  ! This routine calculates the linear system
  !
  ! Author: Satish Karra, PNNL
  ! Date: 06/04/2025
  !

  use Geomechanics_Force_module, only : GeomechForceSetupLinearSystem

  implicit none

  class(pm_geomech_force_type) :: this
  Vec :: right_hand_side
  Vec :: solution
  Mat :: A
  PetscErrorCode :: ierr

#ifdef PM_GEOMECH_FORCE_DEBUG
  call PrintMsg(this%option,'PMGeomechForce%SetupLinearSystem()')
#endif

  call GeomechForceSetupLinearSystem(A,solution,right_hand_side, &
                                     this%geomech_realization,ierr)

end subroutine PMGeomechForceSetupLinearSystem


! ************************************************************************** !

subroutine PMGeomechForcePreSolve(this)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 12/31/13
  !

  use Geomechanics_Force_module, only : GeomechUpdateSolution, &
                                        GeomechStoreInitialDisp, &
                                        GeomechForceUpdateAuxVars

  implicit none

  class(pm_geomech_force_type) :: this

  PetscBool :: force_update_flag = PETSC_TRUE

  call GeomechRealizUpdateAllCouplerAuxVars(this%geomech_realization, &
                                            force_update_flag)

end subroutine PMGeomechForcePreSolve

! ************************************************************************** !

subroutine PMGeomechForcePostSolve(this)
  !
  ! Author: Satish Karra, PNNL
  ! Date: 06/04/2025
  !

  implicit none

  class(pm_geomech_force_type) :: this

end subroutine PMGeomechForcePostSolve

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
  if (this%option%geomechanics%initial_flag) then
    call GeomechStoreInitialDisp(this%geomech_realization)
    this%option%geomechanics%initial_flag = PETSC_FALSE
  endif

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
