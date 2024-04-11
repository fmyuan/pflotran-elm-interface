module Mode_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Global_Aux_module
  use MpFlow_Aux_module
  use Material_Aux_class
  use InlineSurface_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: auxiliary_type 
    type(global_type), pointer :: Global
    type(MpFlow_type), pointer :: MpFlow
    type(material_type), pointer :: Material
    type(inlinesurface_type), pointer :: InlineSurface
  end type auxiliary_type
  
  public :: AuxInit, &
            AuxDestroy

contains

! ************************************************************************** !

subroutine AuxInit(aux)
  ! 
  ! Nullifies pointers in auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/09/08
  ! 

  implicit none
  
  type(auxiliary_type) :: aux
  
  nullify(aux%Global)
  nullify(aux%MpFlow)

  nullify(aux%Material)
  nullify(aux%InlineSurface)

end subroutine AuxInit

! ************************************************************************** !

subroutine AuxDestroy(aux)
  ! 
  ! Deallocates any allocated pointers in auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/09/08
  ! 

  implicit none
  
  type(auxiliary_type) :: aux
  
  call GlobalAuxDestroy(aux%Global)

  call MpFlowAuxDestroy(aux%MpFlow)
  call MaterialAuxDestroy(aux%Material)
  call InlineSurfaceAuxDestroy(aux%InlineSurface)
  
  nullify(aux%Global)
  nullify(aux%MpFlow)
  nullify(aux%Material)
  nullify(aux%InlineSurface)

end subroutine AuxDestroy

end module Mode_Aux_module
