module PM_ERT_class

#include "petsc/finclude/petscmat.h"
  use petscmat

  use PM_Base_class
  use ERT_Aux_module
  use Realization_Subsurface_class
  use Communicator_Base_module  
  use Option_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_base_type) :: pm_ert_type
    class(realization_subsurface_type), pointer :: realization
    class(communicator_type), pointer :: comm1
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

  implicit none
  
  class(pm_ert_type) :: this
  PetscErrorCode :: ierr
  
end subroutine PMERTInitializeRun

! ************************************************************************** !

function PMERTAcceptSolution(this)
  ! 
  ! PMRichardsAcceptSolution:
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

  use Reactive_Transport_module, only : RTUpdateAuxVars
  
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
  ! Strips members of RT process model
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  ! 
  use Utility_module, only : DeallocateArray

  implicit none
  
  class(pm_ert_type) :: this

  call PMBaseDestroy(this)
  ! destroyed in realization
  nullify(this%realization)
  nullify(this%comm1)

end subroutine PMERTStrip
  
! ************************************************************************** !

subroutine PMERTDestroy(this)
  ! 
  ! Destroys RT process model
  ! 
  ! Author: Piyoosh Jaysaval
  ! Date: 01/22/21
  ! 
  implicit none
  
  class(pm_ert_type) :: this

  call PMERTStrip(this)

end subroutine PMERTDestroy
  
end module PM_ERT_class
