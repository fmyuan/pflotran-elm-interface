module Simulation_MultiRealization_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Simulation_Base_class

  implicit none

  private

  type, public, extends(simulation_base_type) :: &
                                      simulation_multirealization_type
    PetscInt :: num_groups
    PetscInt :: num_realizations
    PetscInt :: num_local_realizations
    PetscInt :: cur_realization
    PetscInt, pointer :: realization_ids(:)
  contains
    procedure, public :: Init => SimulationMRInit
    procedure, public :: InitializeRun => SimulationMRInitializeRun
    procedure, public :: ExecuteRun => SimulationMRExecuteRun
    procedure, public :: FinalizeRun => SimulationMRFinalizeRun
    procedure, public :: Strip => SimulationMRStrip
  end type simulation_multirealization_type

  public :: SimulationMRCreate, &
            SimulationMRInit, &
            SimulationMRInitializeRun, &
            SimulationMRFinalizeRun, &
            SimulationMRStrip, &
            SimulationMRDestroy

contains

! ************************************************************************** !

function SimulationMRCreate(driver,option)
  !
  ! Allocates and initializes a new simulation object
  !
  ! Author: Glenn Hammond
  ! Date: 05/27/21

   use Driver_module
   use Option_module

  class(simulation_multirealization_type), pointer :: SimulationMRCreate
  class(driver_type), pointer :: driver
  type(option_type), pointer :: option

  allocate(SimulationMRCreate)
  call SimulationMRCreate%Init(driver,option)

end function SimulationMRCreate

! ************************************************************************** !

subroutine SimulationMRInit(this,driver,option)
  !
  ! Initializes simulation values
  !
  ! Author: Glenn Hammond
  ! Date: 05/27/21

  use Driver_module
  use Option_module

  class(simulation_multirealization_type) :: this
  class(driver_type), pointer :: driver
  type(option_type), pointer :: option

  call SimulationBaseInit(this,driver)
  this%num_groups = 0
  this%num_realizations = 0
  this%num_local_realizations = 0
  this%cur_realization = 0
  nullify(this%realization_ids)

end subroutine SimulationMRInit

! ************************************************************************** !

subroutine SimulationMRInitializeRun(this)
  !
  ! Initializes simulation
  !
  ! Author: Glenn Hammond
  ! Date: 05/27/21

  use Option_module
  use Input_Aux_module
  use Communicator_Aux_module

  class(simulation_multirealization_type) :: this

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

  ! query user for number of communicator groups and realizations
  string = '-num_groups'
  call InputGetCommandLineInt(string,this%num_groups, &
                              option_found,option)

  string = '-num_realizations'
  call InputGetCommandLineInt(string,this%num_realizations, &
                              option_found,option)

  ! read realization ids from a file - contributed by Xingyuan
  string = '-realization_ids_file'
  call InputGetCommandLineString(string,filename,option_found,option)
  if (option_found) then
    input => InputCreate(IUNIT_TEMP,filename,option)
    if (this%num_realizations == 0) then
      option%io_buffer = '"-num_realizations <int>" must be specified ' // &
        'with an integer value matching the number of ids in ' // &
        '"-realization_ids_file <string>".'
      call PrintErrMsg(option)
    endif
    allocate(realization_ids_from_file(this%num_realizations))
    realization_ids_from_file = 0
    string = &
      '# of realization ids read from file may be too few in StochasticInit()'
    do i = 1, this%num_realizations
      call InputReadPflotranString(input,option)
      call InputReadStringErrorMsg(input,option,string)
      call InputReadInt(input,option,realization_ids_from_file(i))
      call InputErrorMsg(input,option,'realization id', &
                         'MultiSimulationInitialize')
    enddo
    call InputDestroy(input)
  else
    nullify(realization_ids_from_file)
  endif

  ! Realization offset contributed by Xingyuan.  This allows one to specify the
  ! smallest/lowest realization id (other than zero) in a stochastic simulation
  string = '-realization_offset'
  call InputGetCommandLineInt(string,offset,option_found,option)
  if (.not.option_found) then
    offset = 0
  endif

  ! error checking
  if (this%num_groups == 0) then
    option%io_buffer = 'Number of stochastic processor groups not ' // &
                       'initialized. Setting to 1.'
    call PrintWrnMsg(option)
    this%num_groups = 1
  endif
  if (this%num_realizations == 0) then
    option%io_buffer = 'Number of stochastic realizations not ' // &
                       'initialized. Setting to 1.'
    call PrintWrnMsg(option)
    this%num_realizations = 1
  endif
  if (this%num_groups > this%num_realizations) then
    write(string,*) this%num_realizations
    option%io_buffer = 'Number of stochastic realizations (' // adjustl(string)
    write(string,*) this%num_groups
    option%io_buffer = trim(option%io_buffer) // ') must be equal to &
         &or greater than number of processor groups (' // adjustl(string)
    option%io_buffer = trim(option%io_buffer) // ').'
    call PrintErrMsg(option)
  endif

  call CommCreateProcessorGroups(option%driver%comm,this%num_groups)
  call OptionUpdateComm(option)

  ! divvy up the realizations
  this%num_local_realizations = this%num_realizations / &
                                           this%num_groups
  remainder = this%num_realizations - this%num_groups * &
                                         this%num_local_realizations

  ! offset is initialized above after check for '-realization_offset'
  do i = 1, option%driver%comm%mygroup_id-1
    delta = this%num_local_realizations
    if (i < remainder) delta = delta + 1
    offset = offset + delta
  enddo

  if (option%driver%comm%mygroup_id < remainder) &
    this%num_local_realizations = this%num_local_realizations + 1
  allocate(this%realization_ids(this%num_local_realizations))
  this%realization_ids = 0
  do i = 1, this%num_local_realizations
    this%realization_ids(i) = offset + i
  enddo

  ! map ids from file - contributed by Xingyuan
  if (associated(realization_ids_from_file)) then
    do i = 1, this%num_local_realizations
      this%realization_ids(i) = &
        realization_ids_from_file(this%realization_ids(i))
    enddo
  endif

  call OptionDestroy(option)

end subroutine SimulationMRInitializeRun

! ************************************************************************** !

subroutine SimulationMRExecuteRun(this)
  !
  ! Execute a simulation
  !
  ! Author: Glenn Hammond
  ! Date: 05/27/21

  use Option_module
  use Factory_Forward_module
  use Simulation_Subsurface_class

  class(simulation_multirealization_type) :: this

  type(option_type), pointer :: option
  class(simulation_subsurface_type), pointer :: forward_simulation

  nullify(forward_simulation)
  do
    option => OptionCreate()
    call OptionSetDriver(option,this%driver)
    call SimulationMRIncrement(this,option)
    call FactoryForwardInitialize(forward_simulation,option)
    call forward_simulation%InitializeRun()
    if (option%status == PROCEED) then
      call forward_simulation%ExecuteRun()
    endif
    call forward_simulation%FinalizeRun()
    call forward_simulation%Strip()
    deallocate(forward_simulation)
    nullify(forward_simulation)
    if (this%cur_realization >= this%num_local_realizations) exit
  enddo

end subroutine SimulationMRExecuteRun

! ************************************************************************** !

subroutine SimulationMRIncrement(this,option)
  ! Author: Glenn Hammond
  ! Date: 02/04/09, 01/06/14

  use Option_module

  class(simulation_multirealization_type) :: this
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  call OptionInitRealization(option)

  this%cur_realization = this%cur_realization + 1
  ! Set group prefix based on id
  option%id = &
    this%realization_ids(this%cur_realization)
  write(string,'(i6)') option%id
  option%group_prefix = 'R' // trim(adjustl(string))

end subroutine SimulationMRIncrement

! ************************************************************************** !

subroutine SimulationMRFinalizeRun(this)
  !
  ! Finalizes simulation
  !
  ! Author: Glenn Hammond
  ! Date: 05/27/21

  class(simulation_multirealization_type) :: this

  call SimulationBaseFinalizeRun(this)

end subroutine SimulationMRFinalizeRun

! ************************************************************************** !

subroutine SimulationMRStrip(this)

  ! Deallocates members of multirealization simulation

  ! Author: Glenn Hammond
  ! Date: 05/27/21

  class(simulation_multirealization_type) :: this

  call SimulationBaseStrip(this)

end subroutine SimulationMRStrip

! ************************************************************************** !

subroutine SimulationMRDestroy(simulation)
  !
  ! Deallocates a simulation
  !
  ! Author: Glenn Hammond
  ! Date: 05/27/21

  class(simulation_multirealization_type), pointer :: simulation

  call simulation%Strip()
  deallocate(simulation)
  nullify(simulation)

end subroutine SimulationMRDestroy

end module Simulation_MultiRealization_class
