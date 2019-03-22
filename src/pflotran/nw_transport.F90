module NW_Transport_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Global_Aux_module
  use Material_Aux_class
  use PM_NWT_class
  use NW_Transport_Aux_module
  
  use PFLOTRAN_Constants_module
  
  implicit none
  
  private
  
  PetscReal, parameter :: perturbation_tolerance = 1.d-5
  
  public :: NWTTimeCut, &
            NWTSetup
            
contains

! ************************************************************************** !

subroutine NWTTimeCut(realization)
  ! 
  ! Resets arrays for a time step cut.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/12/2019
  ! 
 
  use Realization_Subsurface_class
  use Option_module
  use Field_module
  use Global_module
  !use Secondary_Continuum_module, only : SecondaryRTTimeCut
 
  implicit none
  
  type(realization_subsurface_type) :: realization
  type(field_type), pointer :: field
  type(option_type), pointer :: option
  
  PetscErrorCode :: ierr

  field => realization%field
  option => realization%option
 
  ! copy previous solution back to current solution
  call VecCopy(field%tran_yy,field%tran_xx,ierr);CHKERRQ(ierr)
  
  ! set densities and saturations to t+dt
  if (realization%option%nflowdof > 0) then
    call GlobalWeightAuxVars(realization, &
                             realization%option%transport%tran_weight_t1)
  endif

  !if (option%use_mc) then
  !  call SecondaryRTTimeCut(realization)
  !endif
 
end subroutine NWTTimeCut

! ************************************************************************** !

subroutine NWTSetup(realization)
  ! 
  ! Sets up the nuclear waste transport realization.
  ! Author: Jenn Frederick
  ! Date: 03/12/2019
  ! 
  
  use Realization_Subsurface_class
  
  implicit none

  type(realization_subsurface_type) :: realization
  
end subroutine NWTSetup



end module NW_Transport_module