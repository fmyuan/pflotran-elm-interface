module Reaction_Sandbox_Radon_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use Global_Aux_module
  use Reaction_Sandbox_Base_class
  use Reactive_Transport_Aux_module

  implicit none

  private

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_radon_type
    character(len=MAXWORDLENGTH) :: species_name
    character(len=MAXWORDLENGTH) :: mineral_name
    PetscInt :: species_id
    PetscInt :: mineral_id
    PetscReal :: radon_generation_rate ! mol/m^3-sec (mineral)
  contains
    procedure, public :: ReadInput => RadonRead
    procedure, public :: Setup => RadonSetup
    procedure, public :: Evaluate => RadonEvaluate
    procedure, public :: Destroy => RadonDestroy
  end type reaction_sandbox_radon_type

  public :: RadonCreate

contains

! ************************************************************************** !

function RadonCreate()
  !
  ! Allocates radon reaction object.
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/22
  !

  implicit none

  class(reaction_sandbox_radon_type), pointer :: RadonCreate

  allocate(RadonCreate)
  RadonCreate%species_name = ''
  RadonCreate%mineral_name = ''
  RadonCreate%species_id = UNINITIALIZED_INTEGER
  RadonCreate%mineral_id = UNINITIALIZED_INTEGER
  RadonCreate%radon_generation_rate = UNINITIALIZED_DOUBLE
  nullify(RadonCreate%next)

end function RadonCreate

! ************************************************************************** !

subroutine RadonRead(this,input,option)
  !
  ! Reads input deck for radon reaction parameters (if any)
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/22
  !
  use Option_module
  use String_module
  use Input_Aux_module
  use Units_module, only : UnitsConvertToInternal

  implicit none

  class(reaction_sandbox_radon_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word, internal_units

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', &
                       'CHEMISTRY,REACTION_SANDBOX,RADON')
    call StringToUpper(word)

    select case(trim(word))
      case('SPECIES_NAME')
        call InputReadWord(input,option,this%species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,word, &
                           'CHEMISTRY,REACTION_SANDBOX,RADON')
      case('MINERAL_NAME')
        call InputReadWord(input,option,this%mineral_name,PETSC_TRUE)
        call InputErrorMsg(input,option,word, &
                           'CHEMISTRY,REACTION_SANDBOX,RADON')
      case('RADON_GENERATION_RATE')
        call InputReadDouble(input,option,this%radon_generation_rate)
        call InputErrorMsg(input,option,word, &
                           'CHEMISTRY,REACTION_SANDBOX,RADON')
        internal_units = 'mol/m^3-sec'
        call InputReadAndConvertUnits(input,this%radon_generation_rate, &
                                internal_units,'CHEMISTRY,REACTION_SANDBOX,&
                                &RADON,'//word,option)
      case default
        call InputKeywordUnrecognized(input,word, &
                     'CHEMISTRY,REACTION_SANDBOX,RADON',option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine RadonRead

! ************************************************************************** !

subroutine RadonSetup(this,reaction,option)
  !
  ! Sets up the radon reaction with parameters either read from the
  ! input deck or hardwired.
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/22
  !
  use Reaction_Aux_module, only : reaction_rt_type, GetPrimarySpeciesIDFromName
  use Reaction_Mineral_Aux_module, only : GetMineralIDFromName
  use Option_module

  implicit none

  class(reaction_sandbox_radon_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  this%species_id = &
    GetPrimarySpeciesIDFromName(this%species_name,reaction,option)
  this%mineral_id = &
    GetMineralIDFromName(this%mineral_name,reaction%mineral,option)
  if (Uninitialized(this%radon_generation_rate)) then
    option%io_buffer = 'A "RADON_GENERATION_RATE" must be defined in the &
      &RADON REACTION_SANDBOX.'
    call PrintErrMsg(option)
  endif

end subroutine RadonSetup

! ************************************************************************** !

subroutine RadonEvaluate(this,Residual,Jacobian,compute_derivative, &
                         rt_auxvar,global_auxvar,material_auxvar, &
                         reaction,option)
  !
  ! Evaluates the reaction storing the Residual and/or Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/22
  !
  use Option_module
  use Reaction_Aux_module
  use Material_Aux_module

  implicit none

  class(reaction_sandbox_radon_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp)
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  Residual(this%species_id) = Residual(this%species_id) - &
    (1.d0) * & ! positive stoichiometry for generation
    this%radon_generation_rate * &  ! mole/m^3-sec
    rt_auxvar%mnrl_volfrac(this%mineral_id) * & ! m^3 mnrl/m^3 bulk
    material_auxvar%volume ! m^3 bulk

  if (compute_derivative) then

    ! no derivative as generation rate is solely a function of mineral
    ! volume fraction which is not a primary dependent variable.

  endif

end subroutine RadonEvaluate

! ************************************************************************** !

subroutine RadonDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this module
  !
  ! Author: Glenn Hammond
  ! Date: 03/14/22
  !
  implicit none

  class(reaction_sandbox_radon_type) :: this

end subroutine RadonDestroy

end module Reaction_Sandbox_Radon_class
