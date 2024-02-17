module Carbon_Sandbox_MEND_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Carbon_Sandbox_Base_class
  use Option_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, extends(carbon_sandbox_base_type), public :: carbon_sandbox_mend_type
  contains
    procedure, public :: ReadInput => CarbonMENDReadInput
    procedure, public :: Setup => CarbonMENDSetup
!    procedure, public :: Evaluate => CarbonMENDEvaluate
    procedure, public :: Strip => CarbonMENDStrip
  end type carbon_sandbox_mend_type

  type, extends(carbon_sandbox_rxn_base_type), public :: &
                                                 carbon_sandbox_rxn_mend_type
  contains
    procedure, public :: ReadInput => CarbonRxnMENDReadInput
    procedure, public :: Setup => CarbonRxnMENDSetup
    procedure, public :: Evaluate => CarbonRxnMENDEvaluate
    procedure, public :: Strip => CarbonRxnMENDStrip
  end type carbon_sandbox_rxn_mend_type

  public :: CarbonMENDCreate, &
            CarbonRxnMENDCreate

contains

! ************************************************************************** !

function CarbonMENDCreate()
  !
  ! Allocates and initializes the carbon sandbox object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  class(carbon_sandbox_mend_type), pointer :: CarbonMENDCreate

  class(carbon_sandbox_mend_type), pointer :: this

  allocate(this)
  call CarbonBaseInit(this)

  CarbonMENDCreate => this

end function CarbonMENDCreate

! ************************************************************************** !

subroutine CarbonMENDReadInput(this,input,option)
  !
  ! Reads parameters from input deck
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Input_Aux_module
  use Option_module
  use String_module

  class(carbon_sandbox_mend_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: err_string
  PetscBool :: found

  err_string = 'CHEMISTRY,CARBON_SANDBOX,MEND'
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',err_string)
    call StringToUpper(keyword)

    found = PETSC_FALSE
    call CarbonBaseReadSelectCase(this,input,keyword,found,err_string,option)
    if (found) cycle

    select case(keyword)
      case default
        call InputKeywordUnrecognized(input,keyword,err_string,option)
    end select
  enddo

end subroutine CarbonMENDReadInput

! ************************************************************************** !

subroutine CarbonMENDSetup(this,reaction,option)
  !
  ! Configures the reaction and associated data structures
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Option_module
  use Reaction_Aux_module

  class(carbon_sandbox_mend_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  call CarbonBaseSetup(this,reaction,option)

end subroutine CarbonMENDSetup

! ************************************************************************** !

subroutine CarbonMENDStrip(this)
  !
  ! Destroys members of the carbon sandbox MEND object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Utility_module

  class(carbon_sandbox_mend_type) :: this

  call CarbonBaseStrip(this)

end subroutine CarbonMENDStrip

! ************************************************************************** !

function CarbonRxnMENDCreate()
  !
  ! Allocates and initializes the carbon sandbox rxn object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  class(carbon_sandbox_rxn_mend_type), pointer :: CarbonRxnMENDCreate

  class(carbon_sandbox_rxn_mend_type), pointer :: this

  allocate(this)
  this%rate_constant = UNINITIALIZED_DOUBLE
  this%reaction_string = ''
  nullify(this%aux)
  nullify(this%reaction_equation)
  nullify(this%next)

  CarbonRxnMENDCreate => this

end function CarbonRxnMENDCreate

! ************************************************************************** !

subroutine CarbonRxnMENDReadInput(this,input,option)
  !
  ! Reads reaction parameters from input deck
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Input_Aux_module
  use Option_module
  use String_module

  class(carbon_sandbox_rxn_mend_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: err_string

  err_string = 'CHEMISTRY,CARBON_SANDBOX,REACTION_NETWORK'

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',err_string)

    call StringToUpper(keyword)
    select case(trim(keyword))
      case('RATE_CONSTANT')
        call InputReadDouble(input,option,this%rate_constant)
        call InputErrorMsg(input,option,keyword,err_string)
      case('REACTION_EQUATION')
        this%reaction_string = trim(input%buf)
      case default
        call InputKeywordUnrecognized(input,keyword,err_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (Uninitialized(this%rate_constant)) then
    option%io_buffer = 'RATE_CONSTANT not defined in ' // trim(err_string)
    call PrintErrMsg(option)
  endif
  if (len_trim(this%reaction_string) <= 1) then
    option%io_buffer = 'A REACTION_EQUATION is not defined in ' // &
      trim(err_string)
    call PrintErrMsg(option)
  endif

end subroutine CarbonRxnMENDReadInput

! ************************************************************************** !

subroutine CarbonRxnMENDSetup(this,aux,reaction,option)
  !
  ! Configures the reaction and associated data structures
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Reaction_Aux_module

  class(carbon_sandbox_rxn_mend_type) :: this
  type(carbon_sandbox_aux_type), pointer :: aux
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  call CarbonRxnBaseSetup(this,aux,reaction,option)

end subroutine CarbonRxnMENDSetup

! ************************************************************************** !

subroutine CarbonRxnMENDEvaluate(this,Residual,Jacobian,option)
  !
  ! Evaluates the rate expression
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Option_module
  use Reaction_Inhibition_Aux_module

  class(carbon_sandbox_rxn_mend_type) :: this
  PetscReal :: Residual(:)
  PetscReal :: Jacobian(:,:)
  type(option_type) :: option

  PetscInt, pointer :: specid(:)
  PetscReal, pointer :: stoich(:)
  PetscReal :: effective_rate
  PetscReal :: conc,inhibition,inhibition_factor
  PetscReal :: dummy
  PetscInt :: i, icomp, ncomp

  ncomp = this%reaction_equation%nspec
  specid => this%reaction_equation%specid
  stoich => this%reaction_equation%stoich

  effective_rate = log(this%rate_constant)
  inhibition = 0.d0
  do i = 1, ncomp
    if (stoich(i) > 0.d0) cycle
    icomp = specid(i)
    ! subtract due to negative stoichiometry
    effective_rate = effective_rate - stoich(i) * this%aux%ln_conc(icomp)
    conc = this%aux%conc(icomp) / this%aux%liter_water
    call ReactionInhibitionSmoothStep(conc,1.d-20,inhibition_factor,dummy)
    inhibition = inhibition + log(inhibition_factor)
  enddo
  effective_rate = exp(effective_rate+inhibition)
  do i = 1, ncomp
    icomp = specid(i)
    Residual(icomp) = Residual(icomp) - stoich(i)*effective_rate
  enddo

end subroutine CarbonRxnMENDEvaluate

! ************************************************************************** !

recursive subroutine CarbonRxnMENDStrip(this)
  !
  ! Destroys the carbon sandbox object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  class(carbon_sandbox_rxn_mend_type) :: this

  call CarbonRxnBaseStrip(this)

end subroutine CarbonRxnMENDStrip

end module Carbon_Sandbox_MEND_class
