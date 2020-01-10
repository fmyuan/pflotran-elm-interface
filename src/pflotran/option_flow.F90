module Option_Flow_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: flow_option_type 
    ! fluid(s): it's NOT phase(s), rather moveable (mobile) matter or mixture by e.g. pressure
    PetscInt :: nfluid      ! 1: liq_fluid, 2: +air_fluid, 3:+oil_fluid
    PetscInt :: nspecliq    ! species no. in xxx_fluid, with ZERO(0) for non-existing of such a fluid
    PetscInt :: nspecgas    !
    PetscInt :: nspecoil    !
    PetscInt :: nspecflow   !
  
    PetscReal :: dt    ! The size of the time step for flow mode.

    PetscBool :: store_fluxes
    PetscBool :: transient_porosity
    PetscBool :: density_depends_on_salinity
    PetscBool :: numerical_derivatives
    PetscBool :: numerical_derivatives_compare
    PetscBool :: num_as_alyt_derivs
    PetscBool :: resdef
    PetscBool :: flowSolverLinearDone
    PetscBool :: flowTimestepperDone


    PetscBool :: only_vertical_flow
    PetscBool :: isothermal
    PetscBool :: isobaric
    PetscInt  :: ice_model            ! specify water/ice/vapor phase partitioning model
    PetscReal :: frzthw_halfwidth     ! freezing-thawing smoothing half-width (oC)

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

  option%nfluid    = 0
  option%nspecliq  = UNINITIALIZED_INTEGER
  option%nspecgas  = UNINITIALIZED_INTEGER
  option%nspecoil  = UNINITIALIZED_INTEGER
  option%nspecflow = UNINITIALIZED_INTEGER

  option%dt = 0.d0
    
  option%store_fluxes = PETSC_FALSE
  option%transient_porosity = PETSC_FALSE
  option%only_vertical_flow = PETSC_FALSE
  option%density_depends_on_salinity = PETSC_FALSE
  option%numerical_derivatives = PETSC_FALSE
  option%numerical_derivatives_compare = PETSC_FALSE
  option%num_as_alyt_derivs = PETSC_FALSE
  option%resdef = PETSC_FALSE
  option%flowSolverLinearDone = PETSC_FALSE
  option%flowTimestepperDone = PETSC_FALSE

  option%isothermal  = PETSC_FALSE
  option%isobaric    = PETSC_FALSE
  option%ice_model        = PAINTER_EXPLICIT
  option%frzthw_halfwidth = UNINITIALIZED_DOUBLE

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
