module Carbon_Sandbox_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Carbon_Sandbox_Base_class
  use Carbon_Sandbox_MEND_class
  use Carbon_Sandbox_MEND20_class

  use PFLOTRAN_Constants_module

  implicit none

  private

  class(carbon_sandbox_base_type), pointer, public :: carbon_sandbox_list

  interface CarbonSandboxRead
    module procedure CarbonSandboxRead1
    module procedure CarbonSandboxRead2
  end interface

  interface CarbonSandboxDestroy
    module procedure CarbonSandboxDestroy1
    module procedure CarbonSandboxDestroy2
  end interface

  public :: CarbonSandboxInit, &
            CarbonSandboxRead, &
            CarbonSandboxSkipInput, &
            CarbonSandboxSetup, &
            CarbonSandboxEvaluate, &
            CarbonSandboxDestroy

contains

! ************************************************************************** !

subroutine CarbonSandboxInit(option)
  !
  ! Initializes the carbon sandbox list
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  use Option_module
  implicit none
  type(option_type) :: option

  if (associated(carbon_sandbox_list)) then
    call CarbonSandboxDestroy()
  endif
  nullify(carbon_sandbox_list)

end subroutine CarbonSandboxInit

! ************************************************************************** !

subroutine CarbonSandboxRead1(input,option)
  !
  ! Reads input deck for reaction sandbox parameters
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option

  call CarbonSandboxRead(carbon_sandbox_list,input,option)

end subroutine CarbonSandboxRead1

! ************************************************************************** !

subroutine CarbonSandboxRead2(local_sandbox_list,input,option)
  !
  ! Reads input deck for reaction sandbox parameters
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module

  implicit none

  class(carbon_sandbox_base_type), pointer :: local_sandbox_list
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: err_string
  class(carbon_sandbox_base_type), pointer :: new_sandbox, cur_sandbox

  err_string = 'CHEMISTRY,CARBON_SANDBOX'

  call InputReadCard(input,option,keyword)
  call InputErrorMsg(input,option,'keyword',err_string)
  call StringToUpper(keyword)
  nullify(new_sandbox)
  call InputPushBlock(input,option)
  select case(trim(keyword))
    case('FIRST_ORDER')
      new_sandbox => CarbonBaseCreate()
    case('MEND')
      new_sandbox => CarbonMENDCreate()
    case('MEND20')
      new_sandbox => CarbonMEND20Create()
    case default
      call InputKeywordUnrecognized(input,keyword,err_string,option)
  end select
  call new_sandbox%ReadInput(input,option)
  if (.not.associated(local_sandbox_list)) then
    local_sandbox_list => new_sandbox
  else
    cur_sandbox => local_sandbox_list
    do
      if (.not.associated(cur_sandbox%next)) exit
      cur_sandbox => cur_sandbox%next
    enddo
    cur_sandbox%next => new_sandbox
  endif
  nullify(new_sandbox)
  call InputPopBlock(input,option)

end subroutine CarbonSandboxRead2

! ************************************************************************** !

subroutine CarbonSandboxSkipInput(input,option)
  !
  ! Intelligently skips over REACTION_SANDBOX block
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  use Option_module
  use String_module
  use Input_Aux_module
  use Utility_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option

  class(carbon_sandbox_base_type), pointer :: dummy_list

  nullify(dummy_list)
  call CarbonSandboxRead(dummy_list,input,option)
  call CarbonSandboxDestroy(dummy_list)

end subroutine CarbonSandboxSkipInput

! ************************************************************************** !

subroutine CarbonSandboxSetup(reaction,option)
  !
  ! Configures the reactions and associated data structures in entire list
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  use Option_module
  use Reaction_Aux_module, only : reaction_rt_type

  implicit none

  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  class(carbon_sandbox_base_type), pointer :: cur_sandbox

  ! sandbox reactions
  cur_sandbox => carbon_sandbox_list
  do
    if (.not.associated(cur_sandbox)) exit
    call cur_sandbox%Setup(reaction,option)
    cur_sandbox => cur_sandbox%next
  enddo

end subroutine CarbonSandboxSetup

! ************************************************************************** !

subroutine CarbonSandboxEvaluate(Residual,Jacobian,compute_derivative, &
                                 rt_auxvar,global_auxvar,material_auxvar, &
                                 reaction,option)
  !
  ! Evaluates reaction rate expression storing residual and/or Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_module, only: material_auxvar_type

  implicit none

  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  class(carbon_sandbox_base_type), pointer :: cur_sandbox

  cur_sandbox => carbon_sandbox_list
  do
    if (.not.associated(cur_sandbox)) exit
    call cur_sandbox%Evaluate(Residual,Jacobian,compute_derivative, &
                              rt_auxvar,global_auxvar,material_auxvar, &
                              reaction,option)
    cur_sandbox => cur_sandbox%next
  enddo

end subroutine CarbonSandboxEvaluate

! ************************************************************************** !

subroutine CarbonSandboxDestroy1()
  !
  ! Destroys master sandbox list
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !

  implicit none

  call CarbonSandboxDestroy(carbon_sandbox_list)

end subroutine CarbonSandboxDestroy1

! ************************************************************************** !

subroutine CarbonSandboxDestroy2(local_sandbox_list)
  !
  ! Destroys arbitrary sandbox list
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !

  implicit none

  class(carbon_sandbox_base_type), pointer :: local_sandbox_list

  class(carbon_sandbox_base_type), pointer :: cur_sandbox, next_sandbox

  ! sandbox reactions
  cur_sandbox => local_sandbox_list
  do
    if (.not.associated(cur_sandbox)) exit
    next_sandbox => cur_sandbox%next
    call cur_sandbox%Strip()
    deallocate(cur_sandbox)
    cur_sandbox => next_sandbox
  enddo
  nullify(local_sandbox_list)

end subroutine CarbonSandboxDestroy2

end module Carbon_Sandbox_module
