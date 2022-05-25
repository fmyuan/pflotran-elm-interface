module Reaction_Sandbox_Lambda_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Global_Aux_module
  use PFLOTRAN_Constants_module
  use Reaction_Sandbox_Base_class
  use Reactive_Transport_Aux_module

  implicit none

  private

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_lambda_type
    character(len=MAXSTRINGLENGTH) :: reaction_network_filename
    PetscInt, pointer :: species_ids(:,:)
    PetscReal, pointer :: stoich(:,:)
  contains
    procedure, public :: ReadInput => LambdaRead
    procedure, public :: Setup => LambdaSetup
    procedure, public :: Evaluate => LambdaEvaluate
    procedure, public :: Destroy => LambdaDestroy
  end type reaction_sandbox_lambda_type

  public :: LambdaCreate

contains

! ************************************************************************** !

function LambdaCreate()
  !
  ! Allocates lambda reaction object.
  !
  ! Author: Katie Muller
  ! Date: 05/25/22

  implicit none

  class(reaction_sandbox_lambda_type), pointer :: LambdaCreate

  allocate(LambdaCreate)
  LambdaCreate%reaction_network_filename = ''
  nullify(LambdaCreate%species_ids)
  nullify(LambdaCreate%stoich)
  nullify(LambdaCreate%next)

end function LambdaCreate

! ************************************************************************** !

subroutine LambdaRead(this,input,option)
  !
  ! Reads input deck for lambda reaction parameters
  !
  ! Author: Katie Muller
  ! Date: 05/25/22

  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(reaction_sandbox_lambda_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  PetscInt :: i
  character(len=MAXWORDLENGTH) :: word, internal_units

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,LAMBDA')
    call StringToUpper(word)

    select case(trim(word))
      case('REACTION_NETWORK')
        call InputReadFilename(input,option,this%reaction_network_filename)
        call InputErrorMsg(input,option,word, &
                           'CHEMISTRY,REACTION_SANDBOX,LAMBDA')
      case default
        call InputKeywordUnrecognized(input,word, &
                     'CHEMISTRY,REACTION_SANDBOX,LAMBDA',option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine LambdaRead

! ************************************************************************** !

subroutine LambdaSetup(this,reaction,option)
  !
  ! Sets up the lambda reaction with parameters either read from the
  ! input deck or hardwired.
  !
  ! Author: Katie Muller
  ! Date: 05/25/22

  use Option_module
  use Reaction_Aux_module

  implicit none

  class(reaction_sandbox_lambda_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  call ReactionNetworkToStoich(reaction,this%reaction_network_filename, &
                               this%species_ids,this%stoich,option)

end subroutine LambdaSetup

! ************************************************************************** !

subroutine LambdaEvaluate(this,Residual,Jacobian,compute_derivative, &
                          rt_auxvar,global_auxvar,material_auxvar, &
                          reaction,option)
  !
  ! Evaluates the reaction storing the Residual and/or Jacobian
  !
  ! Author: Katie Muller
  ! Date: 05/25/22

  use Option_module
  use Reaction_Aux_module
  use Material_Aux_module

  implicit none

  class(reaction_sandbox_lambda_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscInt, parameter :: iphase = 1
  PetscReal :: L_water

  ! just a check; can be removed.
  PetscInt :: icomp, irxn
  do irxn = 1, size(this%species_ids,2)
    print *, 'Reaction: ', irxn
    do icomp = 1, this%species_ids(0,irxn)
      print *, '  ', this%stoich(icomp,irxn), ' ', &
               reaction%primary_species_names(this%species_ids(icomp,irxn))
    enddo
  enddo

  option%io_buffer = ''
  call PrintMsg(option)
  option%io_buffer = 'LambdaEvaluate is currently set up to only print out &
    &the reaction network and stop.'
  call PrintMsg(option)
  stop

end subroutine LambdaEvaluate

! ************************************************************************** !

subroutine LambdaDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this module
  !
  ! Author: Katie Muller
  ! Date: 05/25/22

  use Utility_module

  implicit none

  class(reaction_sandbox_lambda_type) :: this

  call DeallocateArray(this%species_ids)
  call DeallocateArray(this%stoich)

end subroutine LambdaDestroy

end module Reaction_Sandbox_Lambda_class
