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
    !TODO(geh): remove solver and place required convergence settings elsewhere
    type(solver_type), pointer :: solver
    class(realization_base_type), pointer :: realization_base
    class(pm_base_type), pointer :: next
  contains
    procedure, public :: Setup => PMBaseSetup
    procedure, public :: ReadSimulationBlock => PMBaseRead
    procedure, public :: ReadPMBlock => PMBaseRead
    procedure, public :: InitializeRun => PMBaseThisOnly
    procedure, public :: InputRecord => PMBaseInputRecord
    procedure, public :: SetSolver => PMBaseSetSolver
    procedure, public :: FinalizeRun => PMBaseThisOnly
    procedure, public :: Residual => PMBaseResidual
    procedure, public :: Jacobian => PMBaseJacobian
    procedure, public :: UpdateTimestep => PMBaseUpdateTimestep
    procedure, public :: InitializeTimestep => PMBaseThisOnly
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
    procedure, public :: Destroy => PMBaseThisOnly
    procedure, public :: RHSFunction => PMBaseRHSFunction
    procedure, public :: IFunction => PMBaseIFunction
    procedure, public :: IJacobian => PMBaseIJacobian
    procedure, public :: CheckpointBinary => PMBaseCheckpointBinary
    procedure, public :: RestartBinary => PMBaseCheckpointBinary
    procedure, public :: CheckpointHDF5 => PMBaseCheckpointHDF5
    procedure, public :: RestartHDF5 => PMBaseCheckpointHDF5
  end type pm_base_type
  
  type, public :: pm_base_header_type
    PetscInt :: ndof
  end type pm_base_header_type
    
  public :: PMBaseInit, &
            PMBaseInputRecord, &
            PMBaseReadSelectCase, &
            PMBasePrintHeader, &
            PMBaseResidual, &
            PMBaseJacobian, &
            PMBaseRHSFunction
  
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
  this%skip_restart = PETSC_FALSE
  nullify(this%next)
  
end subroutine PMBaseInit

! ************************************************************************** !

subroutine PMBaseRead(this,input)
  use Input_Aux_module
  implicit none
  class(pm_base_type) :: this
  type(input_type), pointer :: input
  print *, 'Must extend PMBaseRead for: ' // trim(this%name)
  stop
end subroutine PMBaseRead

! ************************************************************************** !

subroutine PMBaseReadSelectCase(this,input,keyword,found,error_string,option)

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
    case('SKIP_RESTART')
      this%skip_restart = PETSC_TRUE
    case default
      found = PETSC_FALSE
  end select

end subroutine PMBaseReadSelectCase

! ************************************************************************** !

subroutine PMBaseSetup(this)
  implicit none
  class(pm_base_type) :: this
  print *, 'Must extend PMBaseSetup for: ' // trim(this%name)
  stop
end subroutine PMBaseSetup

! ************************************************************************** !

subroutine PMBaseInputRecord(this)
  implicit none
  class(pm_base_type) :: this
  print *, 'Must extend PMBaseInputRecord for: ' // trim(this%name)
  stop
end subroutine PMBaseInputRecord

! ************************************************************************** !

subroutine PMBaseResidual(this,snes,xx,r,ierr)
  implicit none
  class(pm_base_type) :: this
  SNES :: snes
  Vec :: xx
  Vec :: r
  PetscErrorCode :: ierr
  print *, 'Must extend PMBaseResidual for: ' // trim(this%name)
  stop
end subroutine PMBaseResidual

! ************************************************************************** !

subroutine PMBaseJacobian(this,snes,xx,A,B,ierr)
  implicit none
  class(pm_base_type) :: this
  SNES :: snes
  Vec :: xx
  Mat :: A, B
  PetscErrorCode :: ierr
  print *, 'Must extend PMBaseJacobian for: ' // trim(this%name)
  stop
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
  print *, 'Must extend PMBaseUpdateTimestep for: ' // trim(this%name)
  stop
end subroutine PMBaseUpdateTimestep

! ************************************************************************** !

subroutine PMBaseCheckUpdatePre(this,line_search,X,dX,changed,ierr)
  implicit none
  class(pm_base_type) :: this
  SNESLineSearch :: line_search
  Vec :: X
  Vec :: dX
  PetscBool :: changed
  PetscErrorCode :: ierr
  print *, 'Must extend PMBaseCheckUpdatePre for: ' // trim(this%name)
  stop
end subroutine PMBaseCheckUpdatePre

! ************************************************************************** !

subroutine PMBaseCheckUpdatePost(this,line_search,X0,dX,X1,dX_changed, &
                                 X1_changed,ierr)
  implicit none
  class(pm_base_type) :: this
  SNESLineSearch :: line_search
  Vec :: X0
  Vec :: dX
  Vec :: X1
  PetscBool :: dX_changed
  PetscBool :: X1_changed
  PetscErrorCode :: ierr
  print *, 'Must extend PMBaseCheckUpdatePost for: ' // trim(this%name)
  stop
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
  print *, 'Must extend PMBaseCheckConvergence for: ' // trim(this%name)
  stop
end subroutine PMBaseCheckConvergence

! ************************************************************************** !

subroutine PMBaseThisOnly(this)
  implicit none
  class(pm_base_type) :: this
  print *, 'Must extend PMBaseThisOnly for: ' // trim(this%name)
  stop
end subroutine PMBaseThisOnly

! ************************************************************************** !

subroutine PMBaseThisTime(this,time)
  implicit none
  class(pm_base_type) :: this
  PetscReal :: time
  print *, 'Must extend PMBaseThisTime for: ' // trim(this%name)
  stop
end subroutine PMBaseThisTime

! ************************************************************************** !

subroutine PMBaseThisTimeError(this,time,ierr)
  implicit none
  class(pm_base_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
  print *, 'Must extend PMBaseThisTimeError for: ' // trim(this%name)
  stop
end subroutine PMBaseThisTimeError

! ************************************************************************** !

function PMBaseFunctionThisOnly(this)
  implicit none
  class(pm_base_type) :: this
  PetscBool ::  PMBaseFunctionThisOnly
  PMBaseFunctionThisOnly = PETSC_TRUE
  print *, 'Must extend PMBaseFunctionThisOnly for: ' // trim(this%name)
  stop
end function PMBaseFunctionThisOnly

! ************************************************************************** !

subroutine PMBaseComputeMassBalance(this,mass_balance_array)
  implicit none
  class(pm_base_type) :: this
  PetscReal :: mass_balance_array(:)
  print *, 'Must extend PMBaseComputeMassBalance for: ' // trim(this%name)
  stop
end subroutine PMBaseComputeMassBalance


! ************************************************************************** !

subroutine PMBaseSetSolver(this,solver)
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/15/17

  use Solver_module

  implicit none

  class(pm_base_type) :: this
  type(solver_type), pointer :: solver

  this%solver => solver

end subroutine PMBaseSetSolver

! ************************************************************************** !

subroutine PMBaseRHSFunction(this,ts,time,xx,ff,ierr)
  implicit none
  class(pm_base_type) :: this
  TS :: ts
  PetscReal :: time
  Vec :: xx
  Vec :: ff
  PetscErrorCode :: ierr
  print *, 'Must extend PMBaseRHSFunction for: ' // trim(this%name)
  stop
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
  print *, 'Must extend PMBaseIFunction for: ' // trim(this%name)
  stop
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
  print *, 'Must extend PMBaseIJacobian for: ' // trim(this%name)
  stop
end subroutine PMBaseIJacobian

! ************************************************************************** !

subroutine PMBaseCheckpointBinary(this,viewer)
  implicit none
#include "petsc/finclude/petscviewer.h"      
  class(pm_base_type) :: this
  PetscViewer :: viewer
!  print *, 'Must extend PMBaseCheckpointBinary/RestartBinary.'
!  stop
end subroutine PMBaseCheckpointBinary

! ************************************************************************** !

subroutine PMBaseCheckpointHDF5(this, pm_grp_id)

  use hdf5

  implicit none

  class(pm_base_type) :: this
  integer(HID_T) :: pm_grp_id
!  print *, 'Must extend PMBaseCheckpointHDF5/RestartHDF5.'
!  stop

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

  if (len_trim(this%header) == 0) then
    print *, 'header name needs to be set for PMBaseInitializeTimestep'
    stop
  endif
  string = '(2("=")," ' // trim(this%header) // ' ",' // &
           trim(StringWrite(80-len_trim(this%header)-4)) // '("="))'
  write(string,string)
  call OptionPrint('',this%option)
  call OptionPrint(string,this%option)

end subroutine PMBasePrintHeader

end module PM_Base_class
