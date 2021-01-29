module PM_ERT_class

#include "petsc/finclude/petscksp.h"
  use petscksp

  use PM_Base_class
  use Realization_Subsurface_class
  use Communicator_Base_module  
  use Option_module
  use ERT_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_base_type) :: pm_ert_type
    class(realization_subsurface_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    Vec :: rhs
  contains
    procedure, public :: Setup => PMERTSetup
    procedure, public :: ReadSimulationOptionsBlock => PMERTReadSimOptionsBlock
    procedure, public :: SetRealization => PMERTSetRealization
    procedure, public :: InitializeRun => PMERTInitializeRun
    procedure, public :: FinalizeRun => PMERTFinalizeRun
    procedure, public :: AcceptSolution => PMERTAcceptSolution
    procedure, public :: UpdateSolution => PMERTUpdateSolution
    procedure, public :: UpdateAuxVars => PMERTUpdateAuxVars
    procedure, public :: InputRecord => PMERTInputRecord
    procedure, public :: Destroy => PMERTDestroy
  end type pm_ert_type
  
  public :: PMERTCreate, &
            PMERTInit, &
            PMERTInitializeRun, &
            PMERTStrip

contains

! ************************************************************************** !

function PMERTCreate()
  ! 
  ! Creates ert process model
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  ! 
  implicit none
  
  class(pm_ert_type), pointer :: PMERTCreate

  class(pm_ert_type), pointer :: pm_ert
  
  allocate(pm_ert)
  call PMERTInit(pm_ert)
  pm_ert%name = 'Electrical Resistivity Tomography'
  pm_ert%header = 'ERT'
  
  PMERTCreate => pm_ert
  
end function PMERTCreate

! ************************************************************************** !

subroutine PMERTInit(pm_ert)
  ! 
  ! Initializes ert process model
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  ! 
  implicit none
  
  class(pm_ert_type) :: pm_ert
  
  call PMBaseInit(pm_ert)
  nullify(pm_ert%realization)
  nullify(pm_ert%comm1)

  pm_ert%rhs = PETSC_NULL_VEC

end subroutine PMERTInit

! ************************************************************************** !

subroutine PMERTReadSimOptionsBlock(this,input)
  ! 
  ! Reads input file parameters associated with the ert process model
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  !
  use Input_Aux_module
  use String_module
  use Option_module
 
  implicit none
  
  class(pm_ert_type) :: this
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  type(option_type), pointer :: option
  PetscBool :: found

  option => this%option
  
  error_string = 'ERT Options'
  
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
    call PMBaseReadSimOptionsSelectCase(this,input,keyword,found, &
                                        error_string,option)
    if (found) cycle

    select case(trim(keyword))
      ! Add various options for ERT if needed here
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)
  
end subroutine PMERTReadSimOptionsBlock

! ************************************************************************** !

subroutine PMERTSetup(this)
  ! 
  ! Initializes variables associated with ert
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  ! 

  implicit none
  
  class(pm_ert_type) :: this

  type(ert_type), pointer :: ert

  ! set the communicator
  this%comm1 => this%realization%comm1

end subroutine PMERTSetup

! ************************************************************************** !

subroutine PMERTSetRealization(this,realization)
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  ! 

  use Realization_Subsurface_class  

  implicit none
  
  class(pm_ert_type) :: this
  class(realization_subsurface_type), pointer :: realization

  this%realization => realization
  this%realization_base => realization
  
end subroutine PMERTSetRealization

! ************************************************************************** !

recursive subroutine PMERTInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  ! 
  use Discretization_module

  implicit none
  
  class(pm_ert_type) :: this
  PetscErrorCode :: ierr

  call DiscretizationDuplicateVector(this%realization%discretization, &
                                     this%realization%field%work,this%rhs)

  ! Initialize to zeros
  call VecZeroEntries(this%rhs,ierr); CHKERRQ(ierr)                             
  
end subroutine PMERTInitializeRun

! ************************************************************************** !

subroutine PMERTSolve(this)

  ! Solves the linear systsem for ERT
  ! for each electrodes
  ! Author: Piyoosh Jaysaval
  ! Date: 01/27/2021
  !
  use Patch_module
  use Grid_module
  use Solver_module
  use Field_module  
  use ERT_module

  implicit none

  class(pm_ert_type) :: this

  class(realization_subsurface_type), pointer :: realization
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(solver_type), pointer :: solver
  type(field_type), pointer :: field
  type(ert_auxvar_type), pointer :: ert_auxvars(:)

  PetscInt :: ielec,nelec
  PetscInt :: elec_id, local_elec_id
  PetscInt :: local_id
  PetscInt :: ghosted_id 
  PetscInt :: num_linear_iterations
  PetscInt :: sum_linear_iterations  
  PetscReal :: val 
  PetscReal, pointer :: vec_ptr(:)

  !PetscLogDouble :: log_start_time
  PetscLogDouble :: log_ksp_start_time
  PetscLogDouble :: log_end_time
  PetscErrorCode :: ierr
  KSPConvergedReason :: ksp_reason  

  ! Forward solve start
  call PetscTime(log_ksp_start_time,ierr); CHKERRQ(ierr) 

  solver => this%solver 
  realization => this%realization 
  field => realization%field  
  patch => realization%patch
  grid => patch%grid

  ert_auxvars => patch%aux%ERT%auxvars

  sum_linear_iterations = 0
  num_linear_iterations = 0

  ! Build System matrix
  call ERTCalculateMatrix(realization,solver%M)
  call KSPSetOperators(solver%ksp,solver%M,solver%M,ierr);CHKERRQ(ierr)
  
  ! TODO: Get nelec from survey file
  nelec = 1

  do ielec=1,nelec
    
    ! NB. solution is stored in field%work -> this can be an initial guess
    call VecZeroEntries(field%work,ierr); CHKERRQ(ierr)
 
    ! RHS
    call VecZeroEntries(this%rhs,ierr); CHKERRQ(ierr)
    call VecGetArrayF90(this%rhs,vec_ptr,ierr); CHKERRQ(ierr)    
    ! TODO: Find cell id where the electrode is located
    ! set first cell for now
    ! Only for a local processor
    elec_id = 1

    if (elec_id > 0) then
      ! it should qualify on only proc
      val = 1.0
      vec_ptr(elec_id) = val
    endif
    call VecRestoreArrayF90(this%rhs,vec_ptr,ierr); CHKERRQ(ierr)

    ! Solve system
    call PetscTime(log_ksp_start_time,ierr); CHKERRQ(ierr) 
    call KSPSolve(solver%ksp,this%rhs,field%work,ierr); CHKERRQ(ierr)
    call PetscTime(log_end_time,ierr); CHKERRQ(ierr)
    
    call VecGetArrayF90(field%work,vec_ptr,ierr); CHKERRQ(ierr)
    ! store potentials for each electrode 
    do local_id=1,grid%nlmax
       ghosted_id = grid%nL2G(local_id)
       if (patch%imat(ghosted_id) <= 0) cycle
       ! TODO -> store in auxvars(i) for ERT?
       !ert_auxvars(ghosted_id)%potential(ielec) = vec_ptr(local_id)
       ert_auxvars(ghosted_id)%potential(1) = vec_ptr(local_id)
    enddo 
    call VecRestoreArrayF90(field%work,vec_ptr,ierr); CHKERRQ(ierr)
    call KSPGetIterationNumber(solver%ksp,num_linear_iterations,ierr)
    call KSPGetConvergedReason(solver%ksp,ksp_reason,ierr)
    sum_linear_iterations = sum_linear_iterations + num_linear_iterations
  enddo

  call PetscTime(log_end_time,ierr); CHKERRQ(ierr)

end subroutine PMERTSolve
! ************************************************************************** !

function PMERTAcceptSolution(this)
  ! 
  ! PMERTAcceptSolution:
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  ! 

  implicit none
  
  class(pm_ert_type) :: this
  
  PetscBool :: PMERTAcceptSolution
  
  ! do nothing
  PMERTAcceptSolution = PETSC_TRUE
  
end function PMERTAcceptSolution

! ************************************************************************** !

recursive subroutine PMERTFinalizeRun(this)
  ! 
  ! Finalizes the time stepping
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  ! 

  implicit none
  
  class(pm_ert_type) :: this
  
  ! do something here
  
  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif  
  
end subroutine PMERTFinalizeRun

! ************************************************************************** !

subroutine PMERTUpdateSolution(this)
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  ! 

  implicit none
  
  class(pm_ert_type) :: this
  
end subroutine PMERTUpdateSolution

! ************************************************************************** !

subroutine PMERTUpdateAuxVars(this)
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21

  
  implicit none
  
  class(pm_ert_type) :: this

end subroutine PMERTUpdateAuxVars  

! ************************************************************************** !

subroutine PMERTInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  ! 
  
  implicit none
  
  class(pm_ert_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

end subroutine PMERTInputRecord

! ************************************************************************** !

subroutine PMERTStrip(this)
  ! 
  ! Strips members of ERT process model
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  ! 

  implicit none
  
  class(pm_ert_type) :: this

  PetscErrorCode :: ierr

  call PMBaseDestroy(this)
  ! destroyed in realization
  nullify(this%realization)
  nullify(this%comm1)

  if (this%rhs /= PETSC_NULL_VEC) then
    call VecDestroy(this%rhs,ierr);CHKERRQ(ierr)
  endif  

end subroutine PMERTStrip
  
! ************************************************************************** !

subroutine PMERTDestroy(this)
  ! 
  ! Destroys ERT process model
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  ! 
  implicit none
  
  class(pm_ert_type) :: this

  call PMERTStrip(this)

end subroutine PMERTDestroy
  
end module PM_ERT_class
