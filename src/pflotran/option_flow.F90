module Option_Flow_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: flow_option_type 

    PetscReal :: reference_temperature
    PetscReal :: reference_pressure
    PetscReal :: reference_density(MAX_PHASE)
    PetscReal :: reference_porosity
    PetscReal :: reference_saturation
  
    PetscBool :: store_fluxes
    PetscBool :: transient_porosity
    PetscBool :: creep_closure_on
    PetscBool :: fracture_on
    PetscBool :: only_vertical_flow
    PetscBool :: density_depends_on_salinity
    PetscBool :: quasi_3d
    PetscBool :: numerical_derivatives
    PetscBool :: numerical_derivatives_compare
    PetscBool :: num_as_alyt_derivs
    PetscBool :: only_energy_eq
    PetscBool :: full_perm_tensor
    PetscBool :: steady_state

    ! If true, permeability changes due to pressure
    PetscBool :: update_flow_perm 
    ! Type of averaging scheme for relative permeability
    PetscInt :: rel_perm_aveg
    ! For WIPP_type pc-sat characteristic curves that use Pct
    PetscBool :: pct_updated
    ! flag to use inline surface flow in Richards mode
    PetscBool :: inline_surface_flow
    PetscReal :: inline_surface_Mannings_coeff
    character(len=MAXSTRINGLENGTH) :: inline_surface_region_name
    ! flag to use freezing model in TH mode
    PetscBool :: th_freezing
    ! If true, then secondary init temp is different from prim. init temp
    PetscBool :: set_secondary_init_temp

    PetscReal :: minimum_hydrostatic_pressure

  end type flow_option_type
  
  public :: OptionFlowCreate, &
            OptionFlowInitAll, &
            OptionFlowInitRealization, &
            OptionFlowDestroy

contains

! ************************************************************************** !

function OptionFlowCreate()
  ! 
  ! Allocates and initializes a new Option object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  implicit none
  
  type(flow_option_type), pointer :: OptionFlowCreate
  
  type(flow_option_type), pointer :: option
  
  allocate(option)

  ! DO NOT initialize members of the option type here.  One must decide 
  ! whether the member needs initialization once for all stochastic 
  ! simulations or initialization for every realization (e.g. within multiple 
  ! stochastic simulations).  This is done in OptionInitAll() and
  ! OptionInitRealization()
  call OptionFlowInitAll(option)
  OptionFlowCreate => option
  
end function OptionFlowCreate

! ************************************************************************** !

subroutine OptionFlowInitAll(option)
  ! 
  ! Initializes all option variables
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  implicit none
  
  type(flow_option_type) :: option
  
  ! These variables should only be initialized once at the beginning of a
  ! PFLOTRAN run (regardless of whether stochastic)
  
  call OptionFlowInitRealization(option)

end subroutine OptionFlowInitAll

! ************************************************************************** !

subroutine OptionFlowInitRealization(option)
  ! 
  ! Initializes option variables specific to a single
  ! realization
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

  implicit none
  
  type(flow_option_type) :: option
  
  ! These variables should be initialized once at the beginning of every 
  ! PFLOTRAN realization or simulation of a single realization

  option%reference_pressure = 101325.d0
  option%reference_temperature = 25.d0
  option%reference_density = 0.d0
  option%reference_porosity = 0.25d0
  option%reference_saturation = 1.d0
    
  option%store_fluxes = PETSC_FALSE
  option%transient_porosity = PETSC_FALSE
  option%creep_closure_on = PETSC_FALSE
  option%fracture_on = PETSC_FALSE
  option%only_vertical_flow = PETSC_FALSE
  option%density_depends_on_salinity = PETSC_FALSE
  option%quasi_3d = PETSC_FALSE
  option%numerical_derivatives = PETSC_FALSE
  option%numerical_derivatives_compare = petsc_false
  option%num_as_alyt_derivs= petsc_false
  option%only_energy_eq = PETSC_FALSE
  option%full_perm_tensor = PETSC_FALSE

  option%set_secondary_init_temp = PETSC_FALSE
  option%update_flow_perm = PETSC_FALSE
  option%rel_perm_aveg = UPWIND
  option%pct_updated = PETSC_FALSE
  option%inline_surface_flow           = PETSC_FALSE
  option%inline_surface_Mannings_coeff = 0.02d0
  option%inline_surface_region_name    = ""
  option%set_secondary_init_temp = PETSC_FALSE
  option%minimum_hydrostatic_pressure = -1.d20
  option%th_freezing = PETSC_FALSE
  option%steady_state = PETSC_FALSE

end subroutine OptionFlowInitRealization

! ************************************************************************** !

subroutine OptionFlowDestroy(option)
  ! 
  ! Deallocates an option
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  ! 

  implicit none
  
  type(flow_option_type), pointer :: option
  
  if (.not.associated(option)) return
  ! all kinds of stuff needs to be added here.

  ! all the below should be placed somewhere other than option.F90
  deallocate(option)
  nullify(option)
  
end subroutine OptionFlowDestroy

end module Option_Flow_module
