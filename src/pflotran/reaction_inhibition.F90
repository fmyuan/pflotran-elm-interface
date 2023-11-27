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
          case('INVERSE_MONOD')
            inhibition%itype = INHIBITION_INVERSE_MONOD
          case('THRESHOLD')
            inhibition%itype = INHIBITION_THRESHOLD
            call InputReadDouble(input,option, &
                                 inhibition%inhibition_constant2)
            call InputErrorMsg(input,option,'scaling factor', &
                               trim(error_string)//',TYPE,THRESHOLD')
          case('SMOOTHSTEP')
            inhibition%itype = INHIBITION_SMOOTHSTEP
            call InputReadDouble(input,option, &
                                 inhibition%inhibition_constant2)
            call InputErrorMsg(input,option,'interval', &
                               trim(error_string)//',TYPE,SMOOTHSTEP')
          case default
            call InputKeywordUnrecognized(input,word, &
                                          trim(error_string)//'TYPE',option)
        end select
      case('INHIBITION_CONSTANT')
        call InputReadDouble(input,option,inhibition%inhibition_constant)
        call InputErrorMsg(input,option,word,error_string)
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (len_trim(inhibition%species_name) < 2 .or. &
      inhibition%itype == 0 .or. &
      Uninitialized(inhibition%inhibition_constant)) then
    option%io_buffer = 'A SPECIES_NAME, TYPE, and INHIBITION_CON' // &
      'STANT must be defined for INHIBITION in REACTION "' // &
      trim(reaction_name) // '".'
    call PrintErrMsg(option)
  endif

end subroutine ReactionInhibitionRead

end module Reaction_Inhibition_module
