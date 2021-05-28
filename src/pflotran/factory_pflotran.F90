module Factory_PFLOTRAN_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: FactoryPFLOTRANInitialize, &
            FactoryPFLOTRANCreateSimulation, &
            FactoryPFLOTRANFinalize

contains

! ************************************************************************** !

subroutine FactoryPFLOTRANInitialize(driver)
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

end subroutine FactoryPFLOTRANInitialize

! ************************************************************************** !

function FactoryPFLOTRANCreateSimulation(driver)
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

  class(simulation_base_type), pointer :: FactoryPFLOTRANCreateSimulation

  class(simulation_base_type), pointer :: simulation
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
      simulation => SimulationMRCreate(driver,option)
    case('SUBSURFACE')
      simulation => SimSubsurfCreate(driver,option)
    case('GEOMECHANICS_SUBSURFACE')
      simulation => GeomechanicsSimulationCreate(driver,option)
    case('INVERSE')
    case default
      call driver%PrintErrMsg('Unrecognized SIMULATION_TYPE ' // &
        trim(simulation_type))
  end select
  call InputDestroy(input)
  call OptionDestroy(option)

  FactoryPFLOTRANCreateSimulation => simulation

end function FactoryPFLOTRANCreateSimulation

! ************************************************************************** !

subroutine FactoryPFLOTRANFinalize(driver)
  !
  ! Initializes libraries for simulation
  !
  ! Author: Glenn Hammond
  ! Date: 05/10/21
  !
  use Driver_module
  use HDF5_Aux_module
  use Logging_module

  class(driver_type) :: driver

  PetscErrorCode :: ierr

  call MPI_Barrier(driver%comm%global_comm,ierr)
  call LoggingDestroy() ! can only be called once, even for mult-realization
  call HDF5Finalize()
  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                            '-options_left','no',ierr);CHKERRQ(ierr)
  ! list any PETSc objects that have not been freed - for debugging
  call PetscOptionsSetValue(PETSC_NULL_OPTIONS, &
                            '-objects_left','yes',ierr);CHKERRQ(ierr)
  call PetscFinalize(ierr);CHKERRQ(ierr)

end subroutine FactoryPFLOTRANFinalize

end module Factory_PFLOTRAN_module
