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

subroutine FactoryPFLOTRANInitialize(driver,simulation)
  !
  ! Initializes libraries for simulation
  !
  ! Author: Glenn Hammond
  ! Date: 05/10/21

  use Communicator_Aux_module
  use Driver_module
  use HDF5_Aux_module
  use Input_Aux_module
  use Option_module
  use Logging_module
  use String_module
  use Simulation_Base_class

  class(driver_type), pointer :: driver
  class(simulation_base_type), pointer :: simulation

  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH), pointer :: strings(:)
  PetscBool :: pflotranin_option_found
  PetscBool :: input_prefix_option_found
  PetscBool :: option_found

  call CommInitPetsc(driver%comm)
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

  string = '-screen_output'
  call InputGetCommandLineTruth(string,driver%print_to_screen, &
                                option_found,option)

  string = '-file_output'
  call InputGetCommandLineTruth(string,driver%print_to_file, &
                                option_found,option)

  if (pflotranin_option_found .and. input_prefix_option_found) then
    call driver%PrintErrMsg('Cannot specify both "-pflotranin" and &
      &"-input_prefix" on the command lines.')
  else if (pflotranin_option_found) then
    driver%input_prefix = StringStripFilenameSuffix(driver%input_filename)
  else if (input_prefix_option_found) then
    driver%input_filename = trim(driver%input_prefix) // '.in'
  endif

  driver%global_prefix = driver%input_prefix

  call HDF5Init()
  call LoggingCreate()
  call OptionDestroy(option)

  simulation => FactoryPFLOTRANCreateSimulation(driver)

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
  use Factory_Forward_module
  use Simulation_Base_class
  use Simulation_Geomechanics_class
  use Simulation_Inverse_class
  use Simulation_MultiRealization_class
  use Simulation_Subsurface_class

  class(driver_type), pointer :: driver

  class(simulation_base_type), pointer :: FactoryPFLOTRANCreateSimulation

  class(simulation_base_type), pointer :: simulation
  type(input_type), pointer :: input
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXWORDLENGTH) :: simulation_type
  PetscBool :: stochastic_option_found, bool_flag

  option => OptionCreate()
  call OptionSetDriver(option,driver)
  input => InputCreate(IN_UNIT,driver%input_filename,option)
  string = 'SIMULATION_TYPE'
  call InputFindStringInFile(input,option,string)
  call InputFindStringErrorMsg(input,option,string)
  call InputReadWord(input,option,simulation_type,PETSC_TRUE)
  call InputDestroy(input)

  string = '-stochastic'
  call InputGetCommandLineTruth(string,bool_flag,stochastic_option_found,option)
  if (stochastic_option_found) simulation_type = 'MULTIREALIZATION'

  select case(simulation_type)
    case('MULTIREALIZATION','INVERSE')
      driver%fid_out = DRIVER_OUT_UNIT
      string = trim(driver%global_prefix) // '.out'
      if (driver%comm%myrank == driver%io_rank .and. driver%print_to_file) then
        open(driver%fid_out, file=string, action="write", status="unknown")
      endif
  end select

  select case(simulation_type)
    case('MULTIREALIZATION')
      simulation => SimulationMRCreate(driver)
    case('INVERSE')
      simulation => SimulationInverseCreate(driver)
    case('SUBSURFACE','GEOMECHANICS_SUBSURFACE')
      select case(simulation_type)
        case('SUBSURFACE')
          simulation => SimSubsurfCreate(driver,option)
        case('GEOMECHANICS_SUBSURFACE')
          simulation => GeomechanicsSimulationCreate(driver,option)
      end select
    case default
      call driver%PrintErrMsg('Unrecognized SIMULATION_TYPE ' // &
        trim(simulation_type))
  end select

  select type(simulation)
    class is(simulation_subsurface_type)
      call FactoryForwardInitialize(simulation,driver%input_filename,option)
    class is(simulation_inverse_type)
      call SimulationInverseRead(simulation,option)
    class is(simulation_multirealization_type)
      if (stochastic_option_found) then
        simulation%forward_simulation_filename = driver%input_filename
      else
        call SimulationMRRead(simulation,option)
      endif
  end select

  select type(simulation)
    class is(simulation_subsurface_type)
    class default
      ! only destroy option if not a forward simulation
      call OptionDestroy(option)
  end select

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
