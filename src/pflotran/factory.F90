module Factory_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: Initialize, &
            CreateSimulation, &
            Finalize

contains

! ************************************************************************** !

subroutine Initialize(driver)
  !
  ! Initializes libraries for simulation
  !
  ! Author: Glenn Hammond
  ! Date: 05/10/21

  use Driver_module
  use HDF5_Aux_module
  use Input_Aux_module
  use Option_module
  use Logging_module
  use String_module

  class(driver_type), pointer :: driver

  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH), pointer :: strings(:)
  PetscBool :: pflotranin_option_found
  PetscBool :: input_prefix_option_found
  PetscBool :: option_found

  option => OptionCreate()
  call OptionSetDriver(option,driver)

  ! check for non-default input filename
  driver%input_filename = 'pflotran.in'
  string = '-pflotranin'
  call InputGetCommandLineString(string,driver%input_filename, &
                                 pflotranin_option_found,option)
  string = '-input_prefix'
  call InputGetCommandLineString(string,driver%input_prefix, &
                                 input_prefix_option_found,option)

  if (pflotranin_option_found .and. input_prefix_option_found) then
    call driver%PrintErrMsg('Cannot specify both "-pflotranin" and &
      &"-input_prefix" on the command lines.')
  else if (pflotranin_option_found) then
    strings => StringSplit(driver%input_filename,'.')
    if (size(strings) > 1) then
      driver%input_prefix = StringsMerge(strings(1:size(strings)-1),'.')
    else
      driver%input_prefix = strings(1)
    endif
    deallocate(strings)
    nullify(strings)
  else if (input_prefix_option_found) then
    driver%input_filename = trim(driver%input_prefix) // '.in'
  endif

  string = '-output_prefix'
  call InputGetCommandLineString(string,driver%global_prefix,option_found, &
                                 option)

  call HDF5Init()
  call LoggingCreate()
  call OptionDestroy(option)

end subroutine Initialize

! ************************************************************************** !

function CreateSimulation(driver)
  !
  ! Opens the input file, reads and creates a simulation object matching the
  ! type requested
  !
  ! Author: Glenn Hammond
  ! Date: 05/25/21

  use Driver_module
  use Input_Aux_module
  use Option_module
  use Simulation_Base_class
  use Simulation_Geomechanics_class
  use Simulation_Subsurface_class
  use Simulation_MultiRealization_class

  class(driver_type), pointer :: driver

  class(simulation_base_type), pointer :: CreateSimulation

  type(input_type), pointer :: input
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: simulation_type
  PetscBool :: option_found, bool_flag

  option => OptionCreate()
  call OptionSetDriver(option,driver)
  input => InputCreate(IN_UNIT,driver%input_filename,option)
  string = 'SIMULATION_TYPE'
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)
  call InputReadWord(input,option,simulation_type,PETSC_TRUE)

  string = '-stochastic'
  call InputGetCommandLineTruth(string,bool_flag,option_found,option)
  if (option_found) simulation_type = 'MULTIREALIZATION'

  select case(simulation_type)
    case('MULTIREALIZATION')
      CreateSimulation => SimulationMRCreate(driver,option)
    case('SUBSURFACE')
      CreateSimulation => SimSubsurfCreate(driver,option)
    case('GEOMECHANICS_SUBSURFACE')
      CreateSimulation => GeomechanicsSimulationCreate(driver,option)
    case('INVERSE')
    case default
      call driver%PrintErrMsg('Unrecognized SIMULATION_TYPE ' // &
        trim(simulation_type))
  end select
  call InputDestroy(input)
  call OptionDestroy(option)

end function CreateSimulation

! ************************************************************************** !

subroutine Finalize()
  !
  ! Initializes libraries for simulation
  !
  ! Author: Glenn Hammond
  ! Date: 05/10/21
  !
  use HDF5_Aux_module
  use Logging_module


  PetscErrorCode :: ierr

  call LoggingDestroy()
  call HDF5Finalize()
  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                            '-options_left','no',ierr);CHKERRQ(ierr)
  ! list any PETSc objects that have not been freed - for debugging
  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                            '-objects_left','yes',ierr);CHKERRQ(ierr)
  call PetscFinalize(ierr);CHKERRQ(ierr)

end subroutine Finalize

end module Factory_module
