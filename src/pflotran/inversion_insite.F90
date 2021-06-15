module Inversion_INSITE_class

#include "petsc/finclude/petscvec.h"
  use petscvec

  use PFLOTRAN_Constants_module
  use Inversion_Base_class
  use Realization_Subsurface_class

  implicit none

  private

  type, public, extends(inversion_base_type) :: inversion_insite_type
    class(realization_subsurface_type), pointer :: realization
    Vec :: quantity_of_interest
    PetscInt :: iqoi
  contains
    procedure, public :: Init => InversionINSITEInit
    procedure, public :: Initialize => InversionINSITEInitialize
    procedure, public :: ReadBlock => InversionINSITEReadBlock
    procedure, public :: UpdateParameters => InversionINSITEUpdateParameters
    procedure, public :: CalculateInverse => InversionINSITECalculateInverse
    procedure, public :: CheckConvergence => InversionINSITECheckConvergence
    procedure, public :: CostFunctions => InversionINSITECostFunctions
    procedure, public :: Finalize => InversionINSITEFinalize
    procedure, public :: Strip => InversionINSITEStrip
  end type inversion_insite_type

  public :: InversionINSITECreate, &
            InversionINSITEStrip, &
            InversionINSITEDestroy

contains

! ************************************************************************** !

function InversionINSITECreate(driver)
  !
  ! Allocates and initializes a new inversion object
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  use Driver_module

  class(driver_type), pointer :: driver

  class(inversion_insite_type), pointer :: InversionINSITECreate

  allocate(InversionINSITECreate)
  call InversionINSITECreate%Init(driver)

end function InversionINSITECreate

! ************************************************************************** !

subroutine InversionINSITEInit(this,driver)
  !
  ! Initializes a new inversion object
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  use Variables_module, only : ELECTRICAL_CONDUCTIVITY, PERMEABILITY
  use Driver_module

  class(inversion_insite_type) :: this
  class(driver_type), pointer :: driver

  this%quantity_of_interest = PETSC_NULL_VEC
  this%iqoi = UNINITIALIZED_INTEGER
  nullify(this%realization)
  call InversionBaseInit(this,driver)

end subroutine InversionINSITEInit

! ************************************************************************** !

subroutine InversionINSITEReadBlock(this,input,option)
  ! 
  ! Reads input file parameters associated an INSITE inversion
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/25/16
  !
  use Input_Aux_module
  use String_module
  use Option_module
  use Variables_module, only : PERMEABILITY, ELECTRICAL_CONDUCTIVITY
 
  class(inversion_insite_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  error_string = 'INSITE Inversion'
  
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
    call InversionBaseReadSelectCase(this,input,keyword,found, &
                                     error_string,option)
    if (found) cycle

    select case(trim(keyword))
      case('QUANTITY_OF_INTEREST')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,error_string)
        call StringToUpper(word)
        select case(word)
          case('PERMEABILITY')
            this%iqoi = PERMEABILITY
          case('ELECTRICAL_CONDUCTIVITY')
            this%iqoi = ELECTRICAL_CONDUCTIVITY
          case default
            call InputKeywordUnrecognized(input,keyword,error_string,option)
        end select
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)
  
end subroutine InversionINSITEReadBlock

! ************************************************************************** !

subroutine InversionINSITEInitialize(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  use Discretization_module
  use Material_module
  use Option_module
  use Variables_module, only : PERMEABILITY, ELECTRICAL_CONDUCTIVITY

  class(inversion_insite_type) :: this

  PetscBool :: exists
  character(len=MAXWORDLENGTH) :: word
  PetscErrorCode :: ierr

  ! theck to ensure that quantity of interest exists
  exists = PETSC_FALSE
  select case(this%iqoi)
    case(PERMEABILITY)
      if (this%realization%option%iflowmode /= NULL_MODE) exists = PETSC_TRUE
      word = 'PERMEABILITY'
    case(ELECTRICAL_CONDUCTIVITY)
      if (this%realization%option%igeopmode /= NULL_MODE) exists = PETSC_TRUE
      word = 'ELECTRICAL_CONDUCTIVITY'
    case default
  end select
  if (.not.exists) then
    this%realization%option%io_buffer = 'Inversion for ' // trim(word) // &
      ' cannot be performed with the specified process models.'
    call PrintErrMsg(this%realization%option)
  endif

  ! non-ghosted Vec
  call VecDuplicate(this%realization%field%work, &
                    this%quantity_of_interest,ierr);CHKERRQ(ierr)
  ! ghosted Vec
  call MaterialGetAuxVarVecLoc(this%realization%patch%aux%Material, &
                               this%realization%field%work_loc, &
                               this%iqoi,ZERO_INTEGER)
  call DiscretizationLocalToGlobal(this%realization%discretization, &
                                   this%realization%field%work_loc, &
                                   this%quantity_of_interest,ONEDOF)

end subroutine InversionINSITEInitialize

! ************************************************************************** !

subroutine InversionINSITEUpdateParameters(this)
  !
  ! Updates input parameters
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  use Discretization_module
  use Material_module

  class(inversion_insite_type) :: this

  if (this%quantity_of_interest == PETSC_NULL_VEC) then
    call this%Initialize()
  else
    call DiscretizationGlobalToLocal(this%realization%discretization, &
                                     this%quantity_of_interest, &
                                     this%realization%field%work_loc,ONEDOF)
    call MaterialSetAuxVarVecLoc(this%realization%patch%aux%Material, &
                                 this%realization%field%work_loc, &
                                 this%iqoi, &
                                 ZERO_INTEGER)
  endif

end subroutine InversionINSITEUpdateParameters

! ************************************************************************** !

subroutine InversionINSITECalculateInverse(this)
  !
  ! Calculates update to input parameters
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_insite_type) :: this

  PetscErrorCode :: ierr

  if (this%quantity_of_interest /= PETSC_NULL_VEC) then
    call VecScale(this%quantity_of_interest,2.d0,ierr);CHKERRQ(ierr)
  endif

end subroutine InversionINSITECalculateInverse

! ************************************************************************** !

subroutine InversionINSITECheckConvergence(this)
  !
  ! Check inversion convergence
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
  !
  class(inversion_insite_type) :: this

end subroutine InversionINSITECheckConvergence

! ************************************************************************** !

subroutine InversionINSITECostFunctions(this)
  !
  ! Computes cost functions
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/14/21
  !
  class(inversion_insite_type) :: this

end subroutine InversionINSITECostFunctions

! ************************************************************************** !

subroutine InversionINSITEFinalize(this)
  !
  ! Finalizes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_insite_type) :: this

  call InversionBaseFinalize(this)

end subroutine InversionINSITEFinalize

! ************************************************************************** !

subroutine InversionINSITEStrip(this)
  !
  ! Deallocates members of inversion insite
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_insite_type) :: this

  PetscErrorCode :: ierr

  call InversionBaseStrip(this)

  nullify(this%realization)
  if (this%quantity_of_interest /= PETSC_NULL_VEC) then
    call VecDestroy(this%quantity_of_interest,ierr);CHKERRQ(ierr)
  endif

end subroutine InversionINSITEStrip

! ************************************************************************** !

subroutine InversionINSITEDestroy(inversion)
  !
  ! Deallocates a inversion
  !
  ! Author: Glenn Hammond
  ! Date: 06/04/21
  !
  class(inversion_insite_type), pointer :: inversion

  if (.not.associated(inversion)) return

  call inversion%Strip()
  deallocate(inversion)
  nullify(inversion)

end subroutine InversionINSITEDestroy

end module Inversion_INSITE_class
