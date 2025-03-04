module Init_Subsurface_Tran_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private


  public :: InitFlowGlobalAuxVar

contains

! ************************************************************************** !

subroutine InitFlowGlobalAuxVar(realization,option)
  !
  ! Initializes flow global aux var after transport setup routines.
  !
  ! Author: Glenn Hammond
  ! Date: 12/04/14
  !
  use Realization_Subsurface_class
  use Option_module
  use Global_module
  use Variables_module

  implicit none

  class(realization_subsurface_type) :: realization
  type(option_type), pointer :: option

  ! initialize densities and saturations
  if (option%nflowdof == 0) then
    call GlobalSetAuxVarScalar(realization,option%flow%reference_pressure, &
                               LIQUID_PRESSURE)
    call GlobalSetAuxVarScalar(realization,option%flow%reference_temperature, &
                               TEMPERATURE)
    call GlobalSetAuxVarScalar(realization,option%flow%reference_saturation, &
                               LIQUID_SATURATION)
    call GlobalSetAuxVarScalar(realization, &
                           option%flow%reference_density(option%liquid_phase), &
                               LIQUID_DENSITY)
    if (option%transport%nphase > 1) then
      call GlobalSetAuxVarScalar(realization, &
                                 1.d0-option%flow%reference_saturation, &
                                 GAS_SATURATION)
      call GlobalSetAuxVarScalar(realization, &
                              option%flow%reference_density(option%gas_phase), &
                                 GAS_DENSITY)
    endif
  else
    call GlobalUpdateAuxVars(realization,TIME_T,0.d0)
    call GlobalWeightAuxVars(realization,0.d0)
  endif

end subroutine InitFlowGlobalAuxVar

! ************************************************************************** !

end module Init_Subsurface_Tran_module
