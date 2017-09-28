module PM_Base_class

#include "petsc/finclude/petscts.h"
  use petscts
  use Option_module
  use Output_Aux_module
  use Realization_Base_class
  use Solver_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: pm_base_type
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXWORDLENGTH) :: header
    type(option_type), pointer :: option
    type(output_option_type), pointer :: output_option
    Vec :: solution_vec
    Vec :: residual_vec
    PetscBool :: print_ekg
    PetscBool :: skip_restart
    PetscBool :: steady_state
    ! solver now has to originate in the pm to support pm-dependent defaults
    type(solver_type), pointer :: solver
    class(realization_base_type), pointer :: realization_base
    class(pm_base_type), pointer :: next
  contains
    procedure, public :: Setup => PMBaseSetup
    procedure, public :: ReadSimulationOptionsBlock => PMBaseReadSimOptionsBlock
    procedure, public :: ReadNewtonBlock => PMBaseReadSelectCaseStop
    procedure, public :: ReadTSBlock => PMBaseReadSelectCaseStop
    procedure, public :: ReadPMBlock => PMBaseReadPMBlock
    procedure, public :: InitializeRun => PMBaseThisOnly
    procedure, public :: InputRecord => PMBaseInputRecord
    procedure, public :: InitializeSolver => PMBaseInitializeSolver
    procedure, public :: FinalizeRun => PMBaseThisOnly
    procedure, public :: Residual => PMBaseResidual
    procedure, public :: Jacobian => PMBaseJacobian
    procedure, public :: UpdateTimestep => PMBaseUpdateTimestep
    procedure, public :: InitializeTimestep => PMBaseThisOnly
    procedure, public :: SetupSolvers => PMBaseThisOnly
    procedure, public :: PreSolve => PMBaseThisOnly
    procedure, public :: Solve => PMBaseThisTimeError
    procedure, public :: PostSolve => PMBaseThisOnly
    procedure, public :: FinalizeTimestep => PMBaseThisOnly
    procedure, public :: AcceptSolution => PMBaseFunctionThisOnly
    procedure, public :: CheckUpdatePre => PMBaseCheckUpdatePre
    procedure, public :: CheckUpdatePost => PMBaseCheckUpdatePost
    procedure, public :: CheckConvergence => PMBaseCheckConvergence
    procedure, public :: TimeCut => PMBaseThisOnly
    procedure, public :: UpdateSolution => PMBaseThisOnly
    procedure, public :: UpdateAuxVars => PMBaseThisOnly
    procedure, public :: MaxChange => PMBaseThisOnly
    procedure, public :: ComputeMassBalance => PMBaseComputeMassBalance
    procedure, public :: Destroy => PMBaseDestroy
    procedure, public :: RHSFunction => PMBaseRHSFunction
    procedure, public :: IFunction => PMBaseIFunction
    procedure, public :: IJacobian => PMBaseIJacobian
    procedure, public :: CheckpointBinary => PMBaseCheckpointBinary
    procedure, public :: RestartBinary => PMBaseCheckpointBinary
    procedure, public :: CheckpointHDF5 => PMBaseCheckpointHDF5
    procedure, public :: RestartHDF5 => PMBaseCheckpointHDF5
    procedure, public :: PrintErrMsg => PMBasePrintErrMsg
  end type pm_base_type
  
  type, public :: pm_base_header_type
    PetscInt :: ndof
  end type pm_base_header_type
    
  public :: PMBaseInit, &
            PMBaseInputRecord, &
            PMBaseInitializeSolver, &
            PMBaseReadSimOptionsSelectCase, &
            PMBasePrintHeader, &
            PMBaseResidual, &
            PMBaseJacobian, &
            PMBaseRHSFunction, &
            PMBaseDestroy
  
contains

! ************************************************************************** !

subroutine PMBaseInit(this)

  implicit none
  
  class(pm_base_type) :: this  

  ! Cannot allocate here.  Allocation takes place in daughter class
  this%name = ''
  this%header = ''
  nullify(this%option)
  nullify(this%output_option)
  nullify(this%realization_base)
  nullify(this%solver)
  this%solution_vec = PETSC_NULL_VEC
  this%residual_vec = PETSC_NULL_VEC
  this%print_ekg = PETSC_FALSE
  this%steady_state = PETSC_FALSE
  this%skip_restart = PETSC_FALSE
  nullify(this%next)
  
end subroutine PMBaseInit

! ************************************************************************** !

subroutine PMBaseReadSimOptionsBlock(this,input)
  use Input_Aux_module
  use String_module
  implicit none
  class(pm_base_type) :: this
  type(input_type), pointer :: input

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option
  PetscBool :: found

  option => this%option

  error_string = 'Base Options'

  input%ierr = 0
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
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select

  enddo
  call InputPopBlock(input,option)
  
end subroutine PMBaseReadSimOptionsBlock

! ************************************************************************** !

subroutine PMBaseReadPMBlock(this,input)
  use Input_Aux_module
  implicit none
  class(pm_base_type) :: this
  type(input_type), pointer :: input
  this%option%exit_code = EXIT_FAILURE
  this%option%io_buffer = 'A member routine PMBaseReadPMBlock must &
             &extend for: ' // trim(this%name)
  call PrintErrMsg(this%option)
end subroutine PMBaseReadPMBlock

! ************************************************************************** !

subroutine PMBaseReadSimOptionsSelectCase(this,input,keyword,found, &
                                          error_string,option)

  use Input_Aux_module

  implicit none
  class(pm_base_type) :: this
  type(input_type) :: input

  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type) :: option

  found = PETSC_TRUE
  select case(trim(keyword))
    case('STEADY_STATE')
      this%steady_state = PETSC_TRUE
    case('SKIP_RESTART')
      this%skip_restart = PETSC_TRUE
    case default
      found = PETSC_FALSE
  end select

end subroutine PMBaseReadSimOptionsSelectCase

! ************************************************************************** !

subroutine PMBaseReadSelectCaseStop(this,input,keyword,found, &
                                    error_string,option)
  use Input_Aux_module
  implicit none
  class(pm_base_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option
!  call this%PrintErrMsg('PMBaseReadSelectCaseStop')
  found = PETSC_FALSE
end subroutine PMBaseReadSelectCaseStop

! ************************************************************************** !

subroutine PMBaseSetup(this)
  implicit none
  class(pm_base_type) :: this
  call this%PrintErrMsg('PMBaseSetup')
end subroutine PMBaseSetup

! ************************************************************************** !

subroutine PMBaseInputRecord(this)
  implicit none
  class(pm_base_type) :: this
  call this%PrintErrMsg('PMBaseInputRecord')
end subroutine PMBaseInputRecord

! ************************************************************************** !

subroutine PMBaseResidual(this,snes,xx,r,ierr)
  implicit none
  class(pm_base_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  call this%PrintErrMsg('PMBaseResidual')
end subroutine PMBaseResidual

! ************************************************************************** !

subroutine PMBaseJacobian(this,snes,xx,A,B,ierr)
  implicit none
  class(pm_base_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  call this%PrintErrMsg('PMBaseJacobian')
end subroutine PMBaseJacobian

! ************************************************************************** !

!TODO(geh): replace anything TS BE-related with an array that can be
!           packed/unpacked on either side.
subroutine PMBaseUpdateTimestep(this,dt,dt_min,dt_max,iacceleration, &
                                num_newton_iterations,tfac, &
                                time_step_max_growth_factor)
  implicit none
  class(pm_base_type) :: this
  PetscReal :: dt
  PetscReal :: dt_min,dt_max
  PetscInt :: iacceleration
  PetscInt :: num_newton_iterations
  PetscReal :: tfac(:)
  PetscReal :: time_step_max_growth_factor
  call this%PrintErrMsg('PMBaseUpdateTimestep')
end subroutine PMBaseUpdateTimestep

! ************************************************************************** !

subroutine PMBaseCheckUpdatePre(this,snes,X,dX,changed,ierr)
  implicit none
  class(pm_base_type) :: this
  SNES :: snes
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr
  call this%PrintErrMsg('PMBaseCheckUpdatePre')
end subroutine PMBaseCheckUpdatePre

! ************************************************************************** !

subroutine PMBaseCheckUpdatePost(this,snes,X0,dX,X1,dX_changed, &
                                 X1_changed,ierr)
  implicit none
  class(pm_base_type) :: this
  SNES :: snes
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr
  call this%PrintErrMsg('PMBaseCheckUpdatePost')
end subroutine PMBaseCheckUpdatePost

! ************************************************************************** !

subroutine PMBaseCheckConvergence(this,snes,it,xnorm,unorm,fnorm,reason,ierr)
  implicit none
  class(pm_base_type) :: this
  SNES :: snes
  PetscInt :: it
  PetscReal :: xnorm
  PetscReal :: unorm
  PetscReal :: fnorm
  SNESConvergedReason :: reason
  PetscErrorCode :: ierr
  call this%PrintErrMsg('PMBaseCheckConvergence')
end subroutine PMBaseCheckConvergence

! ************************************************************************** !

subroutine PMBaseThisOnly(this)
  implicit none
  class(pm_base_type) :: this
  call this%PrintErrMsg('PMBaseThisOnly')
end subroutine PMBaseThisOnly

! ************************************************************************** !

subroutine PMBaseThisTime(this,time)
  implicit none
  class(pm_base_type) :: this
  PetscReal :: time
  call this%PrintErrMsg('PMBaseThisTime')
end subroutine PMBaseThisTime

! ************************************************************************** !

subroutine PMBaseThisTimeError(this,time,ierr)
  implicit none
  class(pm_base_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
  call this%PrintErrMsg('PMBaseThisTimeError')
end subroutine PMBaseThisTimeError

! ************************************************************************** !

function PMBaseFunctionThisOnly(this)
  implicit none
  class(pm_base_type) :: this
  PetscBool ::  PMBaseFunctionThisOnly
  PMBaseFunctionThisOnly = PETSC_TRUE
  call this%PrintErrMsg('PMBaseFunctionThisOnly')
end function PMBaseFunctionThisOnly

! ************************************************************************** !

subroutine PMBaseComputeMassBalance(this,mass_balance_array)
  implicit none
  class(pm_base_type) :: this
  PetscReal :: mass_balance_array(:)
  call this%PrintErrMsg('PMBaseComputeMassBalance')
end subroutine PMBaseComputeMassBalance


! ************************************************************************** !

subroutine PMBaseInitializeSolver(this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/15/17

  use Solver_module

  implicit none

  class(pm_base_type) :: this

  this%solver => SolverCreate()

end subroutine PMBaseInitializeSolver

! ************************************************************************** !

subroutine PMBaseRHSFunction(this,ts,time,xx,ff,ierr)
  implicit none
  class(pm_base_type) :: this
  TS :: ts
  PetscReal :: time
  Vec :: xx
  Vec :: ff
  PetscErrorCode :: ierr
  call this%PrintErrMsg('PMBaseRHSFunction')
end subroutine PMBaseRHSFunction

! ************************************************************************** !

subroutine PMBaseIFunction(this,ts,time,U,Udot,F,ierr)
  implicit none
  class(pm_base_type) :: this
  TS :: ts
  PetscReal :: time
  Vec :: U, Udot
  Vec :: F
  PetscErrorCode :: ierr
  call this%PrintErrMsg('PMBaseIFunction')
end subroutine PMBaseIFunction

! ************************************************************************** !

subroutine PMBaseIJacobian(this,ts,time,U,Udot,shift,A,B,ierr)
  implicit none
  class(pm_base_type) :: this
  TS :: ts
  PetscReal :: time
  Vec :: U, Udot
  PetscReal :: shift
  Mat :: A, B
  PetscErrorCode :: ierr
  call this%PrintErrMsg('PMBaseIJacobian')
end subroutine PMBaseIJacobian

! ************************************************************************** !

subroutine PMBaseCheckpointBinary(this,viewer)
  implicit none
#include "petsc/finclude/petscviewer.h"      
  class(pm_base_type) :: this
  PetscViewer :: viewer
!  call this%PrintErrMsg('PMBaseCheckpointBinary')
end subroutine PMBaseCheckpointBinary

! ************************************************************************** !

subroutine PMBaseCheckpointHDF5(this, pm_grp_id)

  use hdf5

  implicit none

  class(pm_base_type) :: this
  integer(HID_T) :: pm_grp_id
!  call this%PrintErrMsg('PMBaseReadSelectCaseStop')

end subroutine PMBaseCheckpointHDF5

! ************************************************************************** !

subroutine PMBasePrintHeader(this)
  !
  ! Prints PM header to screen and file
  !
  ! Author: Glenn Hammond
  ! Date: 08/06/18
  !
  use Option_module
  use String_module

  implicit none

  class(pm_base_type) :: this

  character(len=MAXSTRINGLENGTH) :: string

  if (this%option%print_screen_flag .or. this%option%print_file_flag) then

  if (len_trim(this%header) == 0) then
    this%option%io_buffer = &
      'header name needs to be set for PMBaseInitializeTimestep'
    call PrintErrMsg(this%option)
  endif
  string = '(2("=")," ' // trim(this%header) // ' ",' // &
           trim(StringWrite(80-len_trim(this%header)-4)) // '("="))'
  write(string,string)
  call OptionPrint('',this%option)
  call OptionPrint(string,this%option)

  endif

end subroutine PMBasePrintHeader

! ************************************************************************** !

subroutine PMBasePrintErrMsg(this,subroutine_name)
  implicit none
  class(pm_base_type) :: this
  character(len=*) :: subroutine_name
  this%option%exit_code = EXIT_FAILURE
  this%option%io_buffer = 'A member routine ' // trim(subroutine_name) // &
         ' must extend for: ' //  trim(this%name)
  call PrintErrMsg(this%option)
end subroutine PMBasePrintErrMsg

! ************************************************************************** !

subroutine PMBaseDestroy(this)

  implicit none
  class(pm_base_type) :: this

  nullify(this%option)
  nullify(this%output_option)
  nullify(this%realization_base)
  nullify(this%next)
  call SolverDestroy(this%solver)

end subroutine PMBaseDestroy

end module PM_Base_class
