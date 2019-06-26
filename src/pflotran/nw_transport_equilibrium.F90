module NWT_Equilibrium_module

#include "petsc/finclude/petscsnes.h"
use PFLOTRAN_Constants_module
use NW_Transport_Aux_module

implicit none

public :: NWTEquilibrateConstraint

contains

! ************************************************************************** !

subroutine NWTEquilibrateConstraint(nw_trans,nwt_species_constraint, &
                                    nwt_auxvar,global_auxvar,material_auxvar, &
                                    option)
  ! 
  ! Calculates the transport constraints based on equilibrium conditions.
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/30/2019
  ! 
  
  use Option_module
  use Global_Aux_module
  use Material_Aux_class
  use NWT_Constraint_module
  
  implicit none
  
  type(nw_trans_realization_type), pointer :: nw_trans
  type(nwt_species_constraint_type), pointer :: nwt_species_constraint
  type(nw_transport_auxvar_type) :: nwt_auxvar
  type(global_auxvar_type) :: global_auxvar
  class(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  
  PetscInt :: ispecies
  PetscInt :: c_type
  
  do ispecies = 1,nw_trans%params%nspecies
  
    c_type = nwt_species_constraint%constraint_type(ispecies)
    ! jenn:todo Each of these needs to be modified still:
    select case(c_type)
      case(CONSTRAINT_T_EQUILIBRIUM)
        nwt_auxvar%total_bulk_conc(ispecies) = &
                              nwt_species_constraint%constraint_conc(ispecies)
        nwt_auxvar%aqueous_eq_conc(ispecies) = &
                              nwt_auxvar%total_bulk_conc(ispecies) / &
                              global_auxvar%sat(LIQUID_PHASE) / &
                              material_auxvar%porosity
        nwt_auxvar%sorb_eq_conc(ispecies) = 1.d-40  ! placeholder
        nwt_auxvar%mnrl_eq_conc(ispecies) = 1.d-40  ! placeholder
      case(CONSTRAINT_AQ_EQUILIBRIUM)
        nwt_auxvar%aqueous_eq_conc(ispecies) = &
                              nwt_species_constraint%constraint_conc(ispecies)
        nwt_auxvar%total_bulk_conc(ispecies) = &
                              nwt_auxvar%aqueous_eq_conc(ispecies) * &
                              global_auxvar%sat(LIQUID_PHASE) * &
                              material_auxvar%porosity
        nwt_auxvar%sorb_eq_conc(ispecies) = 1.d-40  ! placeholder
        nwt_auxvar%mnrl_eq_conc(ispecies) = 1.d-40  ! placeholder
      case(CONSTRAINT_PPT_EQUILIBRIUM)
        nwt_auxvar%mnrl_eq_conc(ispecies) = &
                              nwt_species_constraint%constraint_conc(ispecies)
        nwt_auxvar%aqueous_eq_conc(ispecies) = 1.d-40  ! placeholder
        nwt_auxvar%sorb_eq_conc(ispecies) = 1.d-40  ! placeholder
        nwt_auxvar%total_bulk_conc(ispecies) = 1.d-40  ! placeholder
      case(CONSTRAINT_SB_EQUILIBRIUM)
        nwt_auxvar%sorb_eq_conc(ispecies) = &
                              nwt_species_constraint%constraint_conc(ispecies)
        nwt_auxvar%aqueous_eq_conc(ispecies) = 1.d-40  ! placeholder
        nwt_auxvar%mnrl_eq_conc(ispecies) = 1.d-40  ! placeholder
        nwt_auxvar%total_bulk_conc(ispecies) = 1.d-40  ! placeholder
    end select
  
  enddo

end subroutine NWTEquilibrateConstraint

! ************************************************************************** !

end module NWT_Equilibrium_module
