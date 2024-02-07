module Reaction_Inhibition_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Reaction_Inhibition_Aux_module
  use Reaction_Aux_module

  implicit none

  private

  public :: ReactionInhibitionRead

contains

! ************************************************************************** !

subroutine ReactionInhibitionRead(inhibition,input,option,reaction_name, &
                                  error_string)
  !
  ! Reads an inhibition factor for a reaction
  !
  ! Author: Glenn Hammond
  ! Date: 11/27/23
  !
  use Option_module
  use String_module
  use Input_Aux_module

  implicit none

  type(inhibition_type) :: inhibition
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=*) :: reaction_name
  character(len=*) :: error_string

  character(len=MAXWORDLENGTH) :: word
  PetscInt, parameter :: INHIBIT_ABOVE = 1
  PetscInt, parameter :: INHIBIT_BELOW = 2
  PetscInt :: inhibit_above_or_below
  PetscBool :: require_above_or_below

  inhibit_above_or_below = UNINITIALIZED_INTEGER
  require_above_or_below = PETSC_FALSE
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)
    select case(trim(word))
      case('SPECIES_NAME')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,word,error_string)
        inhibition%species_name = word
      case('TYPE')
        call InputReadCard(input,option,word)
        call InputErrorMsg(input,option,word,error_string)
        call StringToUpper(word)
        select case(word)
          case('MONOD')
            inhibition%itype = INHIBITION_MONOD
            require_above_or_below = PETSC_TRUE
          case('INVERSE_MONOD')
            call InputKeywordDeprecated(word,'a combination of &
              &TYPE MONOD and INHIBIT_BELOW_THRESHOLD',option)
          case('THRESHOLD')
            inhibition%itype = INHIBITION_THRESHOLD
            require_above_or_below = PETSC_TRUE
            call InputReadDouble(input,option,inhibition%inhibition_constant)
            if (.not.InputError(input)) then
              option%io_buffer = 'The INHIBITION THRESHOLD scaling factor &
                &must now be specified under a separate keyword "SCALING_&
                &FACTOR <float>".'
              call PrintErrMsg(option)
            endif
          case('SMOOTHSTEP')
            inhibition%itype = INHIBITION_SMOOTHSTEP
            require_above_or_below = PETSC_TRUE
            call InputReadDouble(input,option,inhibition%inhibition_constant)
            if (.not.InputError(input)) then
              option%io_buffer = 'The INHIBITION SMOOTHSTEP interval &
                &must now be specified under a separate keyword "SMOOTHSTEP_&
                &INTERVAL <float>".'
              call PrintErrMsg(option)
            endif
          case default
            call InputKeywordUnrecognized(input,word, &
                                          trim(error_string)//'TYPE',option)
        end select
      case('THRESHOLD_CONCENTRATION')
        call InputReadDouble(input,option,inhibition%inhibition_constant)
        call InputErrorMsg(input,option,word,error_string)
      case('SCALING_FACTOR','SMOOTHSTEP_INTERVAL')
        call InputReadDouble(input,option,inhibition%inhibition_constant2)
        call InputErrorMsg(input,option,word,error_string)
      case('INHIBIT_BELOW_THRESHOLD')
        inhibit_above_or_below = INHIBIT_BELOW
      case('INHIBIT_ABOVE_THRESHOLD')
        inhibit_above_or_below = INHIBIT_ABOVE
      case('INHIBITION_CONSTANT')
        call InputKeywordDeprecated(word,'a combination of &
          &THRESHOLD_CONCENTRATION and (INHIBIT_ABOVE_THRESHOLD or &
          &INHIBIT_BELOW_THRESHOLD)',option)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (Uninitialized(inhibit_above_or_below) .and. require_above_or_below) then
    option%io_buffer = 'Please specify whether to INHIBIT_ABOVE_THRESHOLD or &
      &INHIBIT_BELOW_THRESHOLD concentration for MONOD, THRESHOLD or SMOOTHSTEP &
      &inhibition in "' // &
      trim(error_string) // ' : ' // trim(reaction_name) // '".'
    call PrintErrMsg(option)
  endif

  select case(inhibit_above_or_below)
    case(INHIBIT_ABOVE)
      inhibition%inhibition_constant = -1.d0 * &
                                       dabs(inhibition%inhibition_constant)
    case(INHIBIT_BELOW)
      inhibition%inhibition_constant = dabs(inhibition%inhibition_constant)
  end select

  if (len_trim(inhibition%species_name) < 2 .or. &
      inhibition%itype == 0 .or. &
      Uninitialized(inhibition%inhibition_constant)) then
    option%io_buffer = 'A SPECIES_NAME, TYPE, and INHIBITION_CON' // &
      'STANT must be defined for INHIBITION in REACTION "' // &
      trim(reaction_name) // '".'
    call PrintErrMsg(option)
  endif

  ! set defaults
  if (Uninitialized(inhibition%inhibition_constant2)) then
    select case(inhibition%itype)
      case(INHIBITION_THRESHOLD) ! scaling factor
        inhibition%inhibition_constant2 = 1.d5 / &
                                          dabs(inhibition%inhibition_constant)
      case(INHIBITION_SMOOTHSTEP) ! interval
        inhibition%inhibition_constant2 = 3.d0
    end select
  endif

end subroutine ReactionInhibitionRead

end module Reaction_Inhibition_module
