module PM_Auxiliary_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PM_Base_class
  use Realization_Subsurface_class
  use Communicator_Base_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public, extends(pm_base_type) :: pm_auxiliary_type
    class(realization_subsurface_type), pointer :: realization
    class(communicator_type), pointer :: comm1
    character(len=MAXWORDLENGTH) :: ctype
    type(pm_auxiliary_salinity_type), pointer :: salinity
    procedure(PMAuxliaryEvaluate), pointer :: Evaluate => null()
  contains
    procedure, public :: Setup => PMAuxiliarySetup
    procedure, public :: InitializeRun => PMAuxiliaryInitializeRun
    procedure, public :: FinalizeRun => PMAuxiliaryFinalizeRun
    procedure, public :: InputRecord => PMAuxiliaryInputRecord
    procedure, public :: Destroy => PMAuxiliaryDestroy
  end type pm_auxiliary_type

  type :: pm_auxiliary_salinity_type
    PetscInt :: nspecies
    character(len=MAXWORDLENGTH) :: species_names(6)
    PetscInt :: ispecies(6)
    PetscReal :: molecular_weights(6)
  end type pm_auxiliary_salinity_type
  
  ! interface blocks
  interface
    subroutine PMAuxliaryEvaluate(this,time,ierr)
      import :: pm_auxiliary_type
      implicit none
      class(pm_auxiliary_type) :: this
      PetscReal :: time
      PetscErrorCode :: ierr      
    end subroutine PMAuxliaryEvaluate
  end interface
  
  public :: PMAuxiliaryCreate, &
            PMAuxiliaryInit, &
            PMAuxiliaryCast, &
            PMAuxiliaryRead, &
            PMAuxiliaryInputRecord, &
            PMAuxiliarySetFunctionPointer
  
contains

! ************************************************************************** !

function PMAuxiliaryCreate()
  ! 
  ! Creates reactive transport process models shell
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/14/13
  ! 

  implicit none
  
  class(pm_auxiliary_type), pointer :: PMAuxiliaryCreate

  class(pm_auxiliary_type), pointer :: pm

  allocate(pm)
  call PMAuxiliaryInit(pm)
  
  PMAuxiliaryCreate => pm

end function PMAuxiliaryCreate
  
! ************************************************************************** !

subroutine PMAuxiliaryInit(this)
  ! 
  ! Initializes auxiliary process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  implicit none
  
  class(pm_auxiliary_type) :: this  

  nullify(this%realization)
  nullify(this%comm1)
  nullify(this%salinity)
  this%ctype = ''
  this%name = ''
  
  call PMBaseInit(this)
  ! restart not currently supported for auxiliary pm's, and not needed.
  this%skip_restart = PETSC_TRUE
  
end subroutine PMAuxiliaryInit

! ************************************************************************** !

subroutine PMAuxiliarySetup(this)
  ! 
  ! Sets up auxiliary process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/02/16

  implicit none
  
  class(pm_auxiliary_type) :: this  

end subroutine PMAuxiliarySetup

! ************************************************************************** !

function PMAuxiliaryCast(this)
  ! 
  ! Initializes auxiliary process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  implicit none
  
  class(pm_base_type), pointer :: this  

  class(pm_auxiliary_type), pointer :: PMAuxiliaryCast  
  
  nullify(PMAuxiliaryCast)
  if (.not.associated(this)) return
  select type (this)
    class is (pm_auxiliary_type)
      PMAuxiliaryCast => this
    class default
      !geh: have default here to pass a null pointer if not of type ascii
  end select
  
end function PMAuxiliaryCast

! ************************************************************************** !

subroutine PMAuxiliaryRead(input, option, this)
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/11/13
  !
  use Input_Aux_module
  use Option_module
  use String_module

  implicit none

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  class(pm_auxiliary_type), pointer :: this

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscInt :: i
  PetscReal :: tempreal

  error_string = 'SIMULATION,PROCESS_MODELS,AUXILIARY'
  call InputReadCard(input,option,word,PETSC_FALSE)
  call InputErrorMsg(input,option,'type',error_string)
  call StringToUpper(word)
  error_string = trim(error_string) // ',' // trim(word)

  this%ctype = word
  select case(word)
    case('SALINITY')
      option%flow%density_depends_on_salinity = PETSC_TRUE
      allocate(this%salinity)
      this%salinity%nspecies = 0
      this%salinity%species_names = ''
      this%salinity%ispecies = UNINITIALIZED_INTEGER
      this%salinity%molecular_weights = UNINITIALIZED_DOUBLE
      i = 0
      word = ''
      call InputPushBlock(input,option)
      do
        call InputReadPflotranString(input,option)
        if (InputCheckExit(input,option)) exit
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(word)
        select case(word)
          case('SPECIES')
            i = i + 1
            if (i > 6) then
              option%io_buffer = 'Number of salinity species exceeded.'
              call PrintErrMsgToDev(option, &
                                    'ask for the maximum number of salinity &
                                    &species to be increased')
            endif
            call InputReadWord(input,option,this%salinity% &
                                 species_names(i),PETSC_TRUE)
            call InputErrorMsg(input,option,'species_name',error_string)
            call InputReadDouble(input,option,tempreal)
            if (input%ierr == 0) then
              this%salinity%molecular_weights(i) = tempreal
            else
              ! for now let's print an error message.  Decide on whether to
              ! read from database later.
              call InputErrorMsg(input,option,'molecular weight',error_string)
            endif
          case default
            error_string = trim(error_string) // 'SALINITY'
            call InputKeywordUnrecognized(input,word,error_string,option)
        end select
      enddo
      call InputPopBlock(input,option)
      this%salinity%nspecies = i
    case default
      call InputKeywordUnrecognized(input,word,error_string,option)
  end select

  call PMAuxiliarySetFunctionPointer(this,this%ctype)

end subroutine PMAuxiliaryRead

! ************************************************************************** !

subroutine PMAuxiliarySetFunctionPointer(this,string)
  ! 
  ! Initializes auxiliary process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  use Option_module
  
  implicit none
  
  class(pm_auxiliary_type) :: this
  character(len=*) :: string
  
  this%ctype = trim(string)
  select case(string)
    case('EVOLVING_STRATA')
      this%Evaluate => PMAuxiliaryEvolvingStrata
      this%name = 'auxiliary evolving strata'
      this%header = 'AUXILIARY EVOLVING STRATA'
    case default
      this%option%io_buffer = 'Function pointer "' // trim(string) // '" not &
        &found among available functions in PMAuxiliarySetFunctionPointer.'
      call PrintErrMsg(this%option)
  end select
  
end subroutine PMAuxiliarySetFunctionPointer

! ************************************************************************** !

recursive subroutine PMAuxiliaryInitializeRun(this)
  ! 
  ! Initializes the time stepping
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/10/16 

  implicit none

  class(pm_auxiliary_type) :: this
  
  PetscReal :: time
  PetscInt :: i
  PetscErrorCode :: ierr
  
  ierr = 0
  time = 0.d0
  select case(this%ctype)
    case('EVOLVING_STRATA')
!      call MatSetOption(Jacobian,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_FALSE, &
!                        ierr);CHKERRQ(ierr)
!      call MatSetOption(Jacobian,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE, &
!                        ierr);CHKERRQ(ierr)
   end select

end subroutine PMAuxiliaryInitializeRun

! ************************************************************************** !

recursive subroutine PMAuxiliaryFinalizeRun(this)
  ! 
  ! Finalizes the run
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/22/18
  ! 

  implicit none

  class(pm_auxiliary_type) :: this

  ! do something here

  if (associated(this%next)) then
    call this%next%FinalizeRun()
  endif

end subroutine PMAuxiliaryFinalizeRun

! ************************************************************************** !

subroutine PMAuxiliaryEvolvingStrata(this,time,ierr)
  ! 
  ! Initializes auxiliary process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  use Init_Subsurface_module

  implicit none
  
  class(pm_auxiliary_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr

  PetscInt :: ndof

  call PMBasePrintHeader(this)

  ierr = 0
  call InitSubsurfAssignMatIDsToRegns(this%realization)
  call InitSubsurfAssignMatProperties(this%realization)
  call InitSubsurfaceSetupZeroArrays(this%realization)
  
end subroutine PMAuxiliaryEvolvingStrata

! ************************************************************************** !

subroutine PMAuxiliaryInputRecord(this)
  ! 
  ! Writes ingested information to the input record file.
  ! 
  ! Author: Jenn Frederick, SNL
  ! Date: 04/21/2016
  ! 
  
  implicit none
  
  class(pm_auxiliary_type) :: this

  character(len=MAXWORDLENGTH) :: word
  PetscInt :: id

  id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'pm: '
  write(id,'(a)') this%name

end subroutine PMAuxiliaryInputRecord

! ************************************************************************** !

subroutine PMAuxiliaryDestroy(this)
  ! 
  ! Destroys auxiliary process model
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/10/16

  implicit none
  
  class(pm_auxiliary_type) :: this
  
  ! destroyed in realization
  nullify(this%realization)
  nullify(this%comm1)
  nullify(this%option)
  nullify(this%output_option)

  if (associated(this%salinity)) then
    deallocate(this%salinity)
    nullify(this%salinity)
  endif
  
end subroutine PMAuxiliaryDestroy

end module PM_Auxiliary_class
