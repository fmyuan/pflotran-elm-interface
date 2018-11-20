module Criticality_module


! MODULE DESCRIPTION:
! ===========================================================================
! This module calculates the radionuclide and energy source terms due 
! to waste form criticality.
! ===========================================================================
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none
  
  ! Stores variables relevant to criticality calculations
  type, public :: criticality_mechanism_base_type
    PetscReal :: heat_released
    PetscReal :: sw
    PetscReal :: rho_w
    PetscReal :: temperature
    PetscReal :: k_effective
  end type criticality_mechanism_base_type
  
  ! Stores information regarding the criticality event
  type, public :: criticality_event_type
    PetscBool :: steady_state
    PetscReal :: criticality_time
    PetscBool :: crit_flag
  end type criticality_event_type

  ! Criticality process model object. Extends the pm_waste_form type to
  ! include relevant variables for criticality consequence calculations.
  type, public :: criticality_type
    type(criticality_mechanism_base_type), pointer :: crit_mech
    type(criticality_event_type), pointer :: crit_event
    type(criticality_type), pointer :: next
    
  end type criticality_type
  
! -------------------------------------------------------------------
  public :: CriticalityInit, &
            ReadCriticality, &
            CriticalityCalc

! -------------------------------------------------------------------

contains

subroutine CriticalityInit(this)

  implicit none
  
  type(criticality_type) :: this
  
  allocate(this%crit_mech)
  allocate(this%crit_event)
  nullify(this%next)

  
  this%crit_mech%heat_released=UNINITIALIZED_DOUBLE
  this%crit_mech%sw=UNINITIALIZED_DOUBLE
  this%crit_mech%rho_w=UNINITIALIZED_DOUBLE
  this%crit_mech%temperature=UNINITIALIZED_DOUBLE
  this%crit_mech%k_effective=UNINITIALIZED_DOUBLE
  
  
  this%crit_event%steady_state = PETSC_FALSE
  this%crit_event%crit_flag = PETSC_FALSE
  this%crit_event%criticality_time = UNINITIALIZED_DOUBLE
  
  
end subroutine CriticalityInit

! ************************************************************************** !

function CriticalityCreate()
  
  implicit none
  
  type(criticality_type), pointer :: CriticalityCreate
  type(criticality_type), pointer :: crit
  
  allocate(crit)
  call CriticalityInit(crit)
  
  CriticalityCreate => crit
  
end function CriticalityCreate

! ************************************************************************** !

subroutine ReadCriticality(this,input,option,keyword,error_string,found)
  
  use Input_Aux_module
  use Reaction_Aux_module, only: GetPrimarySpeciesIDFromName
  use Option_module
  use Condition_module, only : ConditionReadValues
  use Dataset_Ascii_class
  use String_module
  use Units_module
  use Region_module
  use Coupler_module
!   use Simulation_Subsurface_class
  
  implicit none
  
  type(criticality_type), pointer :: this
!   class(simulation_subsurface_type) :: simulation
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
  
  character(len=MAXWORDLENGTH) :: word
  type(criticality_type), pointer :: new_criticality, cur_criticality
  type(coupler_type), pointer :: coupler
!   class(realization_subsurface_type), pointer :: realization
  
!   realization => simulation%realization
  
  error_string = trim(error_string) // ',WASTE_FORM'
  found = PETSC_TRUE
  
  select case(trim(keyword))
    case('CRITICALITY')
      allocate(new_criticality)
      new_criticality => CriticalityCreate()
!       coupler => CouplerCreate(SRC_SINK_COUPLER_TYPE)
!       
!       call InputReadWord(input,option,coupler%name,PETSC_TRUE)
!       call InputDefaultMsg(input,option,'Source Sink name')
!       call CouplerRead(coupler,input,option)
!       call RealizationAddCoupler(realization,coupler)
      
      do
        call InputReadPflotranString(input, option)
        if (InputError(input)) exit
        if (InputCheckExit(input,option)) exit
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'keyword',error_string)
        call StringToUpper(word)
        select case (trim(word))
          case('CRITICALITY_TIME')
            call InputReadDouble(input,option,new_criticality% &
                                 crit_event%criticality_time)
        end select
      enddo
      this => new_criticality
  end select    
  
  
  found=PETSC_FALSE

end subroutine ReadCriticality

! ************************************************************************** !

subroutine CriticalityCalc(this,time,ierr)
  
  implicit none
  
  type(criticality_type) :: this
  PetscReal :: time
  PetscErrorCode :: ierr
  
  this%crit_mech%heat_released = time/1000
  
  
end subroutine CriticalityCalc

! ************************************************************************** !



end module Criticality_module
