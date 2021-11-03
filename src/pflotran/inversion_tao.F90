module Inversion_Tao_class

#include "petsc/finclude/petscksp.h"
  use petscksp

  use PFLOTRAN_Constants_module
  use Inversion_Base_class
  use Inversion_Subsurface_class

  implicit none

  private

  type, public, extends(inversion_subsurface_type) :: inversion_tao_type
    Mat :: Jsensitivity
  contains
    procedure, public :: Init => InversionTaoInit
    procedure, public :: ReadBlock => InversionTaoReadBlock
    procedure, public :: Initialize => InversionTaoInitialize
    procedure, public :: Step => InversionTaoStep
    procedure, public :: Finalize => InversionTaoFinalize
    procedure, public :: Strip => InversionTaoStrip
  end type inversion_tao_type

  public :: InversionTaoCreate, &
            InversionTaoFinalize, &
            InversionTaoStrip, &
            InversionTaoDestroy

contains

! ************************************************************************** !

function InversionTaoCreate(driver)
  !
  ! Allocates and initializes a new tao inversion object
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/21
  !
  use Driver_module

  class(driver_type), pointer :: driver

  class(inversion_tao_type), pointer :: InversionTaoCreate

  allocate(InversionTaoCreate)
  call InversionTaoCreate%Init(driver)

end function InversionTaoCreate

! ************************************************************************** !

subroutine InversionTaoInit(this,driver)
  !
  ! Initializes a new inversion object
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/21
  !
  use Driver_module

  class(inversion_tao_type) :: this
  class(driver_type), pointer :: driver

  call InversionSubsurfaceInit(this,driver)

  this%Jsensitivity = PETSC_NULL_MAT

end subroutine InversionTaoInit

! ************************************************************************** !

subroutine InversionTaoReadBlock(this,input,option)

  use Input_Aux_module
  use Option_module
  use String_module

  class(inversion_tao_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found

  error_string = 'Tao Inversion'

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
    call InversionSubsurfReadSelectCase(this,input,keyword,found, &
                                        error_string,option)
    if (found) cycle

    select case(trim(keyword))
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select

  enddo
  call InputPopBlock(input,option)

end subroutine InversionTaoReadBlock

! ************************************************************************** !

subroutine InversionTaoInitialize(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/21
  !
  class(inversion_tao_type) :: this

  call InversionBaseInitialize(this)

end subroutine InversionTaoInitialize

! ************************************************************************** !

subroutine InversionTaoStep(this)
  !
  ! Execute a simulation
  !
  ! Author: Glenn Hammond
  ! Date: 05/27/21

  use Option_module
  use Factory_Forward_module

  class(inversion_tao_type) :: this

  type(option_type), pointer :: option

  option => OptionCreate()
  write(option%group_prefix,'(i6)') this%iteration+1
  option%group_prefix = 'Run' // trim(adjustl(option%group_prefix))
  call OptionSetDriver(option,this%driver)
  call FactoryForwardInitialize(this%forward_simulation, &
                                this%forward_simulation_filename,option)
  call this%forward_simulation%InitializeRun()
  if (option%status == PROCEED) then
    call this%forward_simulation%ExecuteRun()
  endif
  call this%forward_simulation%FinalizeRun()
  call this%forward_simulation%Strip()
  deallocate(this%forward_simulation)
  nullify(this%forward_simulation)

  this%converg_flag = PETSC_FALSE
  if (this%iteration > this%maximum_iteration) this%converg_flag = PETSC_TRUE

end subroutine InversionTaoStep

! ************************************************************************** !

subroutine InversionTaoUpdateParameters(this)
  !
  ! Initializes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/21
  !
  class(inversion_tao_type) :: this

  PetscInt :: num_measurements
  PetscInt :: num_constraints
  PetscInt :: num_rows
  PetscInt, allocatable :: d_nnz(:)
  PetscInt, allocatable :: o_nnz(:)
  PetscErrorCode :: ierr

  num_measurements = 5
  num_constraints = 4
  num_rows = num_measurements + num_constraints

  call MatCreateDense(this%driver%comm%mycomm, &
                      num_measurements, &
                      this%realization%patch%grid%nlmax, &
                      num_measurements, &
                      this%realization%patch%grid%nmax, &
                      PETSC_NULL_SCALAR, &
                      this%Jsensitivity,ierr)


end subroutine InversionTaoUpdateParameters

! ************************************************************************** !

subroutine InversionTaoFinalize(this)
  !
  ! Finalizes inversion
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/21
  !
  class(inversion_tao_type) :: this

  call InversionBaseStrip(this)

end subroutine InversionTaoFinalize

! ************************************************************************** !

subroutine InversionTaoStrip(this)
  !
  ! Deallocates members of inversion tao
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/21
  !
  class(inversion_tao_type) :: this

   call InversionSubsurfaceStrip(this)

end subroutine InversionTaoStrip

! ************************************************************************** !

subroutine InversionTaoDestroy(inversion)
  !
  ! Deallocates a inversion
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/21
  !
  class(inversion_tao_type), pointer :: inversion

  if (.not.associated(inversion)) return

  call inversion%Strip()
  deallocate(inversion)
  nullify(inversion)

end subroutine InversionTaoDestroy

end module Inversion_Tao_class
