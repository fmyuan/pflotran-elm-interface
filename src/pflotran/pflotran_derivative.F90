program pflotran_derivative
#include "petsc/finclude/petscsys.h"

  use Option_module
  use General_Derivative_module
  use EOS_module

  implicit none

  type(option_type), pointer :: option

  option => OptionCreate()
  call OptionInitMPI(option)
  call OptionInitPetsc(option)
  call EOSInit()
  call GeneralDerivativeDriver(option)
  call OptionFinalize(option)

end program pflotran_derivative
