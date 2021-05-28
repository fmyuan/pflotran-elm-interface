!=======================================================================
! Please see LICENSE and COPYRIGHT files at top of repository.
!=======================================================================
program pflotran

  use Communicator_Aux_module
  use Driver_module
  use Option_module
  use Simulation_Base_class
  use Factory_PFLOTRAN_module
  use Factory_Forward_module
  use Factory_Subsurface_module
  use PFLOTRAN_Constants_module
  use PFLOTRAN_Provenance_module, only : PrintProvenanceToScreen

  implicit none

#include "petsc/finclude/petscsys.h"

  class(simulation_base_type), pointer :: simulation
  class(driver_type), pointer :: driver
  type(option_type), pointer :: option
  PetscInt :: iflag
  PetscErrorCode :: ierr

  driver => DriverCreate()
  call CommInitPetsc(driver%comm)
  call FactoryPFLOTRANInitialize(driver)
  simulation => FactoryPFLOTRANCreateSimulation(driver)
  option => OptionCreate()
  call OptionSetDriver(option,driver)
  call FactoryForwardInitialize(simulation,option)
  call simulation%InitializeRun()
  if (driver%status == PROCEED) then
    call simulation%ExecuteRun()
  endif
  call simulation%FinalizeRun(option)
  call SimulationBaseDestroy(simulation)
  call FactoryForwardFinalize(option)
  driver%exit_code = option%exit_code
  call OptionFinalize(option)
  call FactoryPFLOTRANFinalize(driver)
  iflag = driver%exit_code
  call DriverDestroy(driver)
  call exit(iflag)

end program pflotran
