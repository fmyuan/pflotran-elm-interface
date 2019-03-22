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
            NWTSetup, &
            NWTProcessConstraint
            
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

! ************************************************************************** !

subroutine NWTProcessConstraint(nw_trans,constraint_name, &
                                nwt_species_constraint,option)
  ! 
  ! Ensures ordering of species is consistant between the nw_trans object
  ! and the constraint object. I don't know what the point of this is.
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/22/2019
  ! 
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Option_module
  use String_module
  use Utility_module
  use Transport_Constraint_module
  
  implicit none
  
  type(nw_trans_realization_type), pointer :: nw_trans
  character(len=MAXWORDLENGTH) :: constraint_name
  type(nwt_species_constraint_type), pointer :: nwt_species_constraint
  type(option_type) :: option
  
  PetscBool :: found
  PetscInt :: icomp, jcomp
  PetscReal :: constraint_conc(nw_trans%params%ncomp)
  character(len=MAXWORDLENGTH) :: constraint_species_names( &
                                                         nw_trans%params%ncomp)
  
  constraint_conc = 0.d0
  constraint_species_names = ''
  
  do icomp = 1, nw_trans%params%ncomp
    found = PETSC_FALSE
    do jcomp = 1, nw_trans%params%ncomp
      if (StringCompare(nwt_species_constraint%names(icomp), &
                        nw_trans%species_names(jcomp),MAXWORDLENGTH)) then
        found = PETSC_TRUE
        exit
      endif
    enddo
    if (.not.found) then
      option%io_buffer = &
               'Species ' // trim(nwt_species_constraint%names(icomp)) // &
               ' from CONSTRAINT ' // trim(constraint_name) // &
               ' not found among species.'
      call printErrMsg(option)
    else
      constraint_conc(jcomp) = nwt_species_constraint%constraint_conc(icomp)
      constraint_species_names(jcomp) = nwt_species_constraint%names(icomp)
    endif
  enddo
  
  ! place ordered constraint parameters back in original arrays
  nwt_species_constraint%constraint_conc = constraint_conc
  nwt_species_constraint%names = constraint_species_names

  
end subroutine NWTProcessConstraint

end module NW_Transport_module