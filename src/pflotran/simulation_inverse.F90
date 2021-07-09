module Simulation_Inverse_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Simulation_Base_class
  use Simulation_Subsurface_class
  use Inversion_Base_class

  implicit none

  private

  type, public, extends(simulation_base_type) :: &
                                      simulation_inverse_type
    class(simulation_subsurface_type), pointer :: forward_simulation
    class(inversion_base_type), pointer :: inversion
    character(len=MAXSTRINGLENGTH) :: forward_simulation_filename
  contains
    procedure, public :: Init => SimulationInverseInit
    procedure, public :: InitializeRun => SimulationInverseInitializeRun
    procedure, public :: UpdateParameters => SimulationInvUpdateParameters
    procedure, public :: CalculateUpdate => SimulationInvCalculateUpdate
    procedure, public :: CheckBeta => SimulationInvCheckBeta
    procedure, public :: CheckConvergence => SimulationInvCheckConvergence
    procedure, public :: WriteIterationInfo => SimulationInvWriteIterationInfo
    procedure, public :: ExecuteRun => SimulationInverseExecuteRun
    procedure, public :: FinalizeRun => SimulationInverseFinalizeRun
    procedure, public :: Strip => SimulationInverseStrip
  end type simulation_inverse_type

  public :: SimulationInverseCreate, &
            SimulationInverseInit, &
            SimulationInverseRead, &
            SimulationInverseInitializeRun, &
            SimulationInverseFinalizeRun, &
            SimulationInverseStrip, &
            SimulationInverseDestroy

contains

! ************************************************************************** !

function SimulationInverseCreate(driver)
  !
  ! Allocates and initializes a new simulation object
  !
  ! Author: Glenn Hammond
  ! Date: 05/27/21

   use Driver_module

  class(simulation_inverse_type), pointer :: SimulationInverseCreate
  class(driver_type), pointer :: driver

  allocate(SimulationInverseCreate)
  call SimulationInverseInit(SimulationInverseCreate,driver)

end function SimulationInverseCreate

! ************************************************************************** !

subroutine SimulationInverseInit(this,driver)
  !
  ! Initializes simulation values
  !
  ! Author: Glenn Hammond
  ! Date: 05/27/21

  use Driver_module

  class(simulation_inverse_type) :: this
  class(driver_type), pointer :: driver

  call SimulationBaseInit(this,driver)
  nullify(this%forward_simulation)
  nullify(this%inversion)
  this%forward_simulation_filename = ''

end subroutine SimulationInverseInit

! ************************************************************************** !

subroutine SimulationInverseRead(this,option)
  !
  ! Initializes simulation values
  !
  ! Author: Glenn Hammond
  ! Date: 05/27/21

  use Option_module
  use Input_Aux_module
  use String_module
  use Utility_module
  use Inversion_ERT_class
  use Inversion_ERT_class
  use Inversion_INSITE_class

  class(simulation_inverse_type) :: this
  type(option_type), pointer :: option

  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: error_string
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXWORDLENGTH) :: word

  error_string = 'SIMULATION,INVERSION'

  input => InputCreate(IN_UNIT,this%driver%input_filename,option)

  string = 'SIMULATION'
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)
  keyword = ''
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'',error_string)

    call StringToUpper(keyword)
    select case(trim(keyword))
      case('SIMULATION_TYPE')
      case('INVERSION')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'inversion type', &
                           trim(error_string)//','//keyword)
        call StringToUpper(word)
        select case(word)
          case('INSITE')
            this%inversion => InversionINSITECreate(this%driver)
          case('ERT')
            this%inversion => InversionERTCreate(this%driver)
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        end select
        call this%inversion%ReadBlock(input,option)
      case('FORWARD_SIMULATION_FILENAME')
        call InputReadFilename(input,option,this%forward_simulation_filename)
        call InputErrorMsg(input,option,keyword,error_string)
      case default
        call InputKeywordUnrecognized(input,keyword,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)
  call InputDestroy(input)

end subroutine SimulationInverseRead

! ************************************************************************** !

subroutine SimulationInverseInitializeRun(this)
  !
  ! Initializes simulation
  !
  ! Author: Glenn Hammond
  ! Date: 05/27/21

  use Option_module
  use Input_Aux_module
  use Communicator_Aux_module
  use Inversion_ERT_class
  use Inversion_INSITE_class

  class(simulation_inverse_type) :: this

  type(option_type), pointer :: option
  PetscInt :: i
  PetscInt :: offset, delta, remainder
  PetscInt :: realization_id
  character(len=MAXSTRINGLENGTH) :: string
  PetscBool :: option_found
  PetscInt, pointer :: realization_ids_from_file(:)
  character(len=MAXSTRINGLENGTH) :: filename
  type(input_type), pointer :: input
  PetscErrorCode :: ierr

  option => OptionCreate()
  call OptionSetDriver(option,this%driver)

  call SimulationBaseInitializeRun(this)

  call OptionDestroy(option)

end subroutine SimulationInverseInitializeRun

! ************************************************************************** !

subroutine SimulationInverseExecuteRun(this)
  !
  ! Execute a simulation
  !
  ! Author: Glenn Hammond
  ! Date: 05/27/21

  use Option_module
  use Factory_Forward_module
  use Simulation_Subsurface_class
  use Inversion_ERT_class
  use Inversion_INSITE_class

  class(simulation_inverse_type) :: this

  type(option_type), pointer :: option

  PetscInt :: iteration

  iteration = 0
  do
    if (this%inversion%converg_flag) exit
    option => OptionCreate()
    write(option%group_prefix,'(i6)') iteration+1
    option%group_prefix = 'Run' // trim(adjustl(option%group_prefix))
    call OptionSetDriver(option,this%driver)
    call FactoryForwardInitialize(this%forward_simulation, &
                                  this%forward_simulation_filename,option)
    select type(i=>this%inversion)
      class is(inversion_insite_type)
        i%realization => this%forward_simulation%realization
      class is(inversion_ert_type)
      i%realization => this%forward_simulation%realization
    end select
    call this%UpdateParameters()
    call this%forward_simulation%InitializeRun()
    if (option%status == PROCEED) then
      call this%forward_simulation%ExecuteRun()
    endif
    call this%CheckConvergence()
    call this%CalculateUpdate()
    call this%WriteIterationInfo()
    call this%CheckBeta()
    call this%forward_simulation%FinalizeRun()
    call this%forward_simulation%Strip()
    deallocate(this%forward_simulation)
    nullify(this%forward_simulation)
    iteration = iteration + 1
  enddo

end subroutine SimulationInverseExecuteRun

! ************************************************************************** !

subroutine SimulationInvUpdateParameters(this)
  !
  ! Alters forward run parameters based on inverse calculation
  !
  ! Author: Glenn Hammond
  ! Date: 06/02/21

  class(simulation_inverse_type) :: this

  call this%inversion%UpdateParameters()

end subroutine SimulationInvUpdateParameters

! ************************************************************************** !

subroutine SimulationInvCalculateUpdate(this)
  !
  ! Calculates updated parameters
  !
  ! Author: Glenn Hammond
  ! Date: 06/02/21

  class(simulation_inverse_type) :: this

  call this%inversion%CalculateUpdate()

end subroutine SimulationInvCalculateUpdate

! ************************************************************************** !

subroutine SimulationInvCheckConvergence(this)
  !
  ! Calculates updated parameters
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/18/21

  class(simulation_inverse_type) :: this

  call this%inversion%CheckConvergence()

end subroutine SimulationInvCheckConvergence

! ************************************************************************** !

subroutine SimulationInvCheckBeta(this)
  !
  ! Checks beta if it needs cooling or not
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 06/21/21

  class(simulation_inverse_type) :: this

  call this%inversion%CheckBeta()

end subroutine SimulationInvCheckBeta

! ************************************************************************** !

subroutine SimulationInvWriteIterationInfo(this)
  !
  ! Writes information after each iteration
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 07/07/21

  class(simulation_inverse_type) :: this

  PetscInt :: fid
  PetscBool :: print_to_file,print_to_screen

  fid = this%driver%fid_out
  print_to_file = this%driver%print_to_file
  print_to_screen = this%driver%print_to_screen

  if (this%driver%comm%global_rank == this%driver%io_rank) then
    call this%inversion%WriteIterationInfo(fid,print_to_file,print_to_screen)
  endif

end subroutine SimulationInvWriteIterationInfo

! ************************************************************************** !

subroutine SimulationInverseFinalizeRun(this)
  !
  ! Finalizes simulation
  !
  ! Author: Glenn Hammond
  ! Date: 05/27/21

  class(simulation_inverse_type) :: this

  call this%inversion%Finalize()
  call SimulationBaseFinalizeRun(this)
  if (this%driver%comm%global_rank == this%driver%io_rank) then
    call SimulationBaseWriteTimes(this,this%driver%fid_out)
  endif

end subroutine SimulationInverseFinalizeRun

! ************************************************************************** !

subroutine SimulationInverseStrip(this)

  ! Deallocates members of inverse simulation

  ! Author: Glenn Hammond
  ! Date: 05/27/21

  class(simulation_inverse_type) :: this

  call SimulationBaseStrip(this)
  if (associated(this%forward_simulation)) then
    print *, 'Why is forward simulation still associated in &
             &SimulationInverseStrip?'
    stop
  endif
  nullify(this%forward_simulation)
  call this%inversion%Strip()
  deallocate(this%inversion)
  nullify(this%inversion)

end subroutine SimulationInverseStrip

! ************************************************************************** !

subroutine SimulationInverseDestroy(simulation)
  !
  ! Deallocates a simulation
  !
  ! Author: Glenn Hammond
  ! Date: 05/27/21

  class(simulation_inverse_type), pointer :: simulation

  call simulation%Strip()
  deallocate(simulation)
  nullify(simulation)

end subroutine SimulationInverseDestroy

end module Simulation_Inverse_class
