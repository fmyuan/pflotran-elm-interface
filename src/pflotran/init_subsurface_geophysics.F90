module Init_Subsurface_Geop_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: InitSubsurfGeopSetupRealization

contains

! ************************************************************************** !

subroutine InitSubsurfGeopSetupRealization(realization)
  !
  ! Initializes geophsyics data structres and assign them to the domain.
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 01/29/21

  use Realization_Subsurface_class
  use Option_module
  use ERT_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option

  option => realization%option

  select case(option%igeopmode)
    case(ERT_MODE)
      call ERTSetup(realization)
    case default
      option%io_buffer = 'Unknown igeopmode found during <Mode>Setup'
      call PrintErrMsg(option)
    end select

end subroutine

end module Init_Subsurface_Geop_module
