module Simulation_Inverse_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Simulation_Base_class
  use Inversion_Base_class

  implicit none

  private

  type, public, extends(simulation_base_type) :: &
                                      simulation_inverse_type
    class(inversion_base_type), pointer :: inversion
  contains
    procedure, public :: Init => SimulationInverseInit
    procedure, public :: InitializeRun => SimulationInverseInitializeRun
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

   use Driver_class

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

  use Driver_class

  class(simulation_inverse_type) :: this
  class(driver_type), pointer :: driver

  call SimulationBaseInit(this,driver)
  nullify(this%inversion)

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
  use Inversion_Subsurface_class
  use Inversion_ZFlow_class

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
          case('ERT')
            this%inversion => InversionERTCreate(this%driver)
          case('TEST_SENSITIVITY_JACOBIAN')
            this%inversion => InversionSubsurfaceCreate(this%driver)
          case('ZFLOW')
            this%inversion => InversionZFlowCreate(this%driver)
          case default
            call InputKeywordUnrecognized(input,word,error_string,option)
        end select
        call this%inversion%ReadBlock(input,option)
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

  class(simulation_inverse_type) :: this

  type(option_type), pointer :: option
  type(comm_type), pointer :: newcomm
  PetscInt :: num_groups
  PetscErrorCode :: ierr

  ! create process groups here
  nullify(newcomm)
  num_groups = 1
  if (this%inversion%inversion_option%split_comms .and. &
      this%driver%comm%size > 1) then
    num_groups = 2
  endif
  call CommCreateProcessGroups(this%driver%comm,num_groups,PETSC_TRUE, &
                               this%inversion%inversion_option%invcomm,ierr)
  if (ierr /= 0) then
    call this%driver%PrintErrMsg('Unevenly sized MPI comm groups.')
  endif
  if (this%inversion%inversion_option%invcomm%group_id > 1) then
    call CommDestroy(this%inversion%inversion_option%invcomm)
  endif
  call CommCreateProcessGroups(this%driver%comm,num_groups,PETSC_TRUE, &
                               this%inversion%inversion_option%forcomm,ierr)
  if (this%inversion%inversion_option%forcomm%group_id > 1) then
!    call CommDestroy(this%inversion%inversion_option%forcomm)
  endif

  option => OptionCreate()
  call OptionSetDriver(option,this%driver)
  call OptionSetComm(option,this%driver%comm) ! doesn't matter which comm
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

  class(simulation_inverse_type) :: this

  call this%inversion%InitializeIterationNumber()
  do
    if (this%inversion%converged) exit
    call this%inversion%Step()
    call this%inversion%IncrementIteration()
  enddo

end subroutine SimulationInverseExecuteRun

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
  if (this%driver%IsIORank()) then
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
