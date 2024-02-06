module Carbon_Sandbox_Base_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Option_module
  use PFLOTRAN_Constants_module
  use Reaction_Equation_module

  implicit none

  private

  type, public :: carbon_sandbox_base_type
    type(carbon_sandbox_aux_type), pointer :: aux
    class(carbon_sandbox_rxn_base_type), pointer :: rxn_list
    class(carbon_sandbox_base_type), pointer :: next
  contains
!    procedure, public :: ReadInput => CarbonBaseReadInput
    procedure, public :: Setup => CarbonBaseSetup
    procedure, public :: Evaluate => CarbonBaseEvaluate
    procedure, public :: Destroy => CarbonBaseDestroy
  end type carbon_sandbox_base_type

  type, public :: carbon_sandbox_rxn_base_type
    type(carbon_sandbox_aux_type), pointer :: aux
    PetscReal :: rate_constant
    type(reaction_equation_type), pointer :: reaction_equation
    character(len=MAXSTRINGLENGTH) :: reaction_string
    class(carbon_sandbox_rxn_base_type), pointer :: next
  contains
!    procedure, public :: Setup => CarbonRxnBaseSetup
!    procedure, public :: Evaluate => CarbonRxnBaseEvaluate
  end type carbon_sandbox_rxn_base_type

  type, public :: carbon_sandbox_aux_type
    type(option_type), pointer :: option
    PetscReal :: liquid_saturation
    PetscReal :: liquid_density
    PetscReal :: porosity
    PetscReal :: cell_volume
    PetscReal, pointer :: conc(:)
    PetscReal, pointer :: ln_conc(:)
    PetscReal, pointer :: inhibition_conc(:)
  end type carbon_sandbox_aux_type

  public :: CarbonBaseCreate, &
            CarbonRxnBaseCreate

contains

! ************************************************************************** !

function CarbonBaseCreate(option)
  !
  ! Allocates and initializes the carbon sandbox object
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  type(option_type), pointer :: option

  class(carbon_sandbox_base_type), pointer :: CarbonBaseCreate

  class(carbon_sandbox_base_type), pointer :: this

  allocate(this)
  nullify(this%aux)
  nullify(this%rxn_list)
  nullify(this%next)

  CarbonBaseCreate => this

end function CarbonBaseCreate

! ************************************************************************** !

function CarbonRxnBaseCreate(aux)
  !
  ! Allocates and initializes the carbon sandbox rxn object
  !
  ! Author: Glenn Hammond
  ! Date: 02/05/24
  !
  class(carbon_sandbox_aux_type), pointer :: aux

  class(carbon_sandbox_rxn_base_type), pointer :: CarbonRxnBaseCreate

  class(carbon_sandbox_rxn_base_type), pointer :: this

  allocate(this)
  this%rate_constant = UNINITIALIZED_DOUBLE
  this%reaction_string = ''
  nullify(this%aux)
  nullify(this%reaction_equation)
  nullify(this%next)

  CarbonRxnBaseCreate => this

end function CarbonRxnBaseCreate

! ************************************************************************** !

subroutine CarbonBaseReadInput(this,input,option)
  !
  ! Reads parameters from input deck
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  use Input_Aux_module
  use Option_module
  use String_module

  class(carbon_sandbox_rxn_base_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: err_string

  err_string = 'CHEMISTRY,CARBON_SANDBOX,BASE'

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
      case('REACTION')
        this%reaction_string = trim(input%buf)
      case default
        call InputKeywordUnrecognized(input,keyword,err_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (Uninitialized(this%rate_constant)) then
    option%io_buffer = 'Uninitialized RATE_CONSTANT in ' // trim(err_string)
    call PrintErrMsg(option)
  endif
  if (len_trim(this%reaction_string) <= 1) then
    option%io_buffer = 'A REACTION is undefined in ' // &
      trim(err_string)
    call PrintErrMsg(option)
  endif

end subroutine CarbonBaseReadInput

! ************************************************************************** !

subroutine CarbonBaseSetup(this,reaction,option)
  !
  ! Configures the reaction and associated data structures
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  use Option_module
  use Reaction_Aux_module

  class(carbon_sandbox_base_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type), pointer :: option

  allocate(this%aux)
  this%aux%option => option
  this%aux%liquid_saturation = UNINITIALIZED_DOUBLE
  this%aux%liquid_density = UNINITIALIZED_DOUBLE
  this%aux%porosity = UNINITIALIZED_DOUBLE
  this%aux%cell_volume = UNINITIALIZED_DOUBLE
  allocate(this%aux%conc(reaction%ncomp))
  this%aux%conc = UNINITIALIZED_DOUBLE
  allocate(this%aux%ln_conc(reaction%ncomp))
  this%aux%ln_conc = UNINITIALIZED_DOUBLE
  allocate(this%aux%inhibition_conc(reaction%ncomp))
  this%aux%inhibition_conc = UNINITIALIZED_DOUBLE

end subroutine CarbonBaseSetup

! ************************************************************************** !

subroutine CarbonBaseRxnSetup(this,reaction,option)
  !
  ! Configures the reaction and associated data structures
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  use Reaction_Aux_module

  class(carbon_sandbox_rxn_base_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  this%reaction_equation => &
      ReactionEquationCreateFromString(this%reaction_string, &
                                       reaction%naqcomp, &
                                       reaction%offset_aqueous, &
                                       reaction%primary_species_names, &
                                       reaction%nimcomp, &
                                       reaction%offset_immobile, &
                                       reaction%immobile%names, &
                                       PETSC_FALSE,option)

end subroutine CarbonBaseRxnSetup

! ************************************************************************** !

subroutine CarbonBaseEvaluate(this,Residual,Jacobian,compute_derivative, &
                              rt_auxvar,global_auxvar,material_auxvar, &
                              reaction,option)
  !
  ! Evaluates the rate expression
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  use Global_Aux_module
  use Material_Aux_module
  use Option_module
  use Reaction_Aux_module
  use Reaction_Inhibition_Aux_module
  use Reactive_Transport_Aux_module

  class(carbon_sandbox_base_type) :: this
  class(reaction_rt_type) :: reaction
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(:)
  PetscReal :: Jacobian(:,:)
  PetscBool :: compute_derivative
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
#if 0
  PetscInt, pointer :: specid(:)
  PetscReal, pointer :: stoich(:)
  PetscReal :: effective_rate
  PetscReal :: mol_spec(reaction%ncomp)
  PetscReal :: ln_mol_spec(reaction%ncomp)
  PetscReal :: conc,inhibition,inhibition_factor
  PetscReal :: liter_water
  PetscReal :: dummy
  PetscInt :: i, icomp, ncomp

  ncomp = this%reaction_equation%nspec
  specid => this%reaction_equation%specid
  stoich => this%reaction_equation%stoich

  liter_water = material_auxvar%volume* &
                material_auxvar%porosity* &
                global_auxvar%sat(1)*1.d3

  do i = 1, reaction%naqcomp
    mol_spec(i) = rt_auxvar%pri_molal(i)*liter_water
  enddo
  do i = 1, reaction%offset_immobile+reaction%immobile%nimmobile
    mol_spec(reaction%offset_immobile+i) = &
      rt_auxvar%immobile(i)*material_auxvar%volume
  enddo
  ln_mol_spec = log(mol_spec)

  effective_rate = log(this%rate_constant)
  inhibition = 0.d0
  do i = 1, ncomp
    if (stoich(i) > 0.d0) cycle
    icomp = specid(i)
    ! subtract due to negative stoichiometry
    effective_rate = effective_rate - stoich(i) * ln_mol_spec(icomp)
    conc = mol_spec(icomp) / liter_water
    call ReactionInhibitionSmoothStep(conc,1.d-20,inhibition_factor,dummy)
    inhibition = inhibition + log(inhibition_factor)
  enddo
  effective_rate = exp(effective_rate+inhibition)
  do i = 1, ncomp
    icomp = specid(i)
    Residual(icomp) = Residual(icomp) - stoich(i)*effective_rate
  enddo
#endif
end subroutine CarbonBaseEvaluate

! ************************************************************************** !

subroutine CarbonBaseDestroy(this)
  !
  ! Destroys the carbon sandbox object
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  use Utility_module

  class(carbon_sandbox_base_type) :: this

  call DeallocateArray(this%aux%conc)
  call DeallocateArray(this%aux%ln_conc)
  call DeallocateArray(this%aux%conc)
  nullify(this%aux%option)
  deallocate(this%aux)
  nullify(this%aux)
  call CarbonRxnBaseDestroy(this%rxn_list)

end subroutine CarbonBaseDestroy

! ************************************************************************** !

recursive subroutine CarbonRxnBaseDestroy(this)
  !
  ! Destroys the carbon sandbox object
  !
  ! Author: Glenn Hammond
  ! Date: 02/02/24
  !
  class(carbon_sandbox_rxn_base_type), pointer :: this

  if (.not.associated(this)) return

  call ReactionEquationDestroy(this%reaction_equation)
  nullify(this%aux)

  deallocate(this)
  nullify(this)

end subroutine CarbonRxnBaseDestroy

end module Carbon_Sandbox_Base_class
