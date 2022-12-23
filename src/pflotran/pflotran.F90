!=======================================================================
! Please see LICENSE and COPYRIGHT files at top of repository.
!=======================================================================
program pflotran

  use Driver_class
  use Simulation_Base_class
  use Factory_PFLOTRAN_module
  use PFLOTRAN_Constants_module

  implicit none

#include "petsc/finclude/petscsys.h"

  class(simulation_base_type), pointer :: simulation
  class(driver_type), pointer :: driver
  PetscInt :: iflag

  driver => DriverCreate()
  call FactoryPFLOTRANInitialize(driver,simulation)
  call simulation%InitializeRun()
  if (driver%status == PROCEED) then
    call simulation%ExecuteRun()
  endif
  call simulation%FinalizeRun()
  call SimulationBaseDestroy(simulation)
  call FactoryPFLOTRANFinalize(driver)
  iflag = driver%exit_code
  call DriverDestroy(driver)
  call exit(iflag)

end program pflotran
