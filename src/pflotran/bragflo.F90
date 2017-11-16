module Bragflo_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes

  implicit none
  
  private 

  public :: BragfloResidual, &
            BragfloJacobian

contains

! ************************************************************************** !

subroutine BragfloResidual(snes,xx,r,realization,pmwss_ptr,ierr)
  ! 
  ! Computes the residual equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/16/17
  ! 
  use Realization_Subsurface_class
  use WIPP_Flow_module
  use PM_WIPP_SrcSink_class

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  type(realization_subsurface_type) :: realization
  class(pm_wipp_srcsink_type), pointer :: pmwss_ptr
  PetscErrorCode :: ierr
  
  call WIPPFloResidual(snes,xx,r,realization,pmwss_ptr,ierr)

end subroutine BragfloResidual

! ************************************************************************** !

subroutine BragfloJacobian(snes,xx,A,B,realization,pmwss_ptr,ierr)
  ! 
  ! Computes the Jacobian
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/16/17
  ! 
  use Realization_Subsurface_class
  use WIPP_Flow_module
  use PM_WIPP_SrcSink_class

  implicit none

  SNES :: snes
  Vec :: xx
  Mat :: A, B
  type(realization_subsurface_type) :: realization
  class(pm_wipp_srcsink_type), pointer :: pmwss_ptr
  PetscErrorCode :: ierr

  call WIPPFloJacobian(snes,xx,A,B,realization,pmwss_ptr,ierr)

end subroutine BragfloJacobian

end module Bragflo_module
