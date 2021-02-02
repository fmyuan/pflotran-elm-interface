!=======================================================================
! Please see LICENSE and COPYRIGHT files at top of repository.
!=======================================================================
program pflotran
  
  use Option_module
  use Simulation_Base_class
  use Multi_Simulation_module
  use Factory_PFLOTRAN_module
  use Factory_Subsurface_module
  use PFLOTRAN_Constants_module
  use PFLOTRAN_Provenance_module, only : PrintProvenanceToScreen
  
  implicit none

#include "petsc/finclude/petscsys.h"

  class(simulation_base_type), pointer :: simulation
  ! multisimulation enables multiple simulations to be run concurrently
  ! and/or one after another until a specified set of simulations has 
  ! completed.
  type(multi_simulation_type), pointer :: multisimulation
  type(option_type), pointer :: option
  
  nullify(simulation)
  nullify(multisimulation)
  option => OptionCreate()
  call OptionInitMPI(option)
  call PFLOTRANInitializePrePetsc(multisimulation,option)
  call OptionInitPetsc(option)
  if (option%myrank == option%io_rank .and. option%print_to_screen) then
    !call PrintProvenanceToScreen()
  endif

  do ! multi-simulation loop
    call PFLOTRANInitializePostPetsc(simulation,multisimulation,option)

    call simulation%InitializeRun()

    if (option%status == PROCEED) then
      call simulation%ExecuteRun()
    endif
stop
    call simulation%FinalizeRun()
    call simulation%Strip()
    deallocate(simulation)
    nullify(simulation)
    call PFLOTRANFinalize(option)
    if (MultiSimulationDone(multisimulation)) exit
  enddo ! multi-simulation loop
  call OptionFinalize(option)

end program pflotran
