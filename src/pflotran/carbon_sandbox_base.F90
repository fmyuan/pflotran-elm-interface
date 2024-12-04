module Carbon_Sandbox_Base_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Option_module
  use PFLOTRAN_Constants_module
  use Reaction_Equation_module

  implicit none

  private

  PetscInt, parameter, public :: CARBON_UNITS_MOLALITY = 1
  PetscInt, parameter, public :: CARBON_UNITS_MOLARITY = 2
  PetscInt, parameter, public :: CARBON_UNITS_MOLE_PER_KG_SOIL = 3
  PetscInt, parameter, public :: CARBON_UNITS_MOLE_PER_M3_BULK = 4
  PetscInt, parameter, public :: CARBON_UNITS_MOLES = 5

  type, public :: carbon_sandbox_base_type
    PetscInt :: concentration_units
    type(carbon_sandbox_aux_type), pointer :: aux
    class(carbon_sandbox_rxn_base_type), pointer :: rxn_list
    class(carbon_sandbox_base_type), pointer :: next
  contains
    procedure, public :: ReadInput => CarbonBaseReadInput
    procedure, public :: Setup => CarbonBaseSetup
    procedure, public :: GetConcentrationUnits => &
                           CarbonBaseGetConcentrationUnits
    procedure, public :: EnforceConcentrationUnits => &
                           CarbonBaseEnforceConcUnits
    procedure, public :: MapStateVariables => CarbonBaseMapStateVariables
    procedure, public :: Evaluate => CarbonBaseEvaluate
    procedure, public :: Strip => CarbonBaseStrip
  end type carbon_sandbox_base_type

  type, public :: carbon_sandbox_rxn_base_type
    type(carbon_sandbox_aux_type), pointer :: aux
    PetscReal :: rate_constant
    type(reaction_equation_type), pointer :: reaction_equation
    character(len=MAXSTRINGLENGTH) :: reaction_string
    class(carbon_sandbox_rxn_base_type), pointer :: next
  contains
    procedure, public :: ReadInput => CarbonRxnBaseReadInput
    procedure, public :: Setup => CarbonRxnBaseSetup
    procedure, public :: Evaluate => CarbonRxnBaseEvaluate
    procedure, public :: Strip => CarbonRxnBaseStrip
  end type carbon_sandbox_rxn_base_type

  type, public :: carbon_sandbox_aux_type
    PetscReal :: liquid_saturation
    PetscReal :: liquid_density
    PetscReal :: porosity
    PetscReal :: cell_volume
    PetscReal :: liter_water
    PetscReal :: temperature
    PetscReal, pointer :: conc(:)
    PetscReal, pointer :: ln_conc(:)
    PetscReal, pointer :: inhibition_conc(:)
  end type carbon_sandbox_aux_type

  public :: CarbonBaseCreate, &
            CarbonRxnBaseCreate, &
            CarbonBaseInit, &
            CarbonBaseReadSelectCase, &
            CarbonBaseSetup, &
            CarbonBaseStrip, &
            CarbonRxnBaseSetup, &
            CarbonRxnBaseStrip

contains

! ************************************************************************** !

function CarbonBaseCreate()
  !
  ! Allocates and initializes the carbon sandbox object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  class(carbon_sandbox_base_type), pointer :: CarbonBaseCreate

  class(carbon_sandbox_base_type), pointer :: this

  allocate(this)
  call CarbonBaseInit(this)

  CarbonBaseCreate => this

end function CarbonBaseCreate

! ************************************************************************** !

subroutine CarbonBaseInit(this)
  !
  ! Initializes members of carbon sandbox base class
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  class(carbon_sandbox_base_type) :: this

  this%concentration_units = UNINITIALIZED_INTEGER
  nullify(this%aux)
  nullify(this%rxn_list)
  nullify(this%next)

end subroutine CarbonBaseInit

! ************************************************************************** !

subroutine CarbonBaseReadInput(this,input,option)
  !
  ! Reads parameters from input deck
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Input_Aux_module
  use Option_module
  use String_module

  class(carbon_sandbox_base_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: err_string
  PetscBool :: found

  err_string = 'CHEMISTRY,CARBON_SANDBOX,BASE'
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword',err_string)
    call StringToUpper(keyword)

    found = PETSC_FALSE
    call CarbonBaseReadSelectCase(this,input,keyword,found,err_string,option)
    if (.not.found) then
      call InputKeywordUnrecognized(input,keyword,err_string,option)
    endif
  enddo

end subroutine CarbonBaseReadInput

! ************************************************************************** !

subroutine CarbonBaseReadSelectCase(this,input,keyword,found,err_string,option)
  !
  ! Reads parameters from input deck
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Input_Aux_module
  use Option_module
  use String_module

  class(carbon_sandbox_base_type) :: this
  type(input_type), pointer :: input
  character(len=MAXWORDLENGTH) :: keyword
  PetscBool :: found
  character(len=MAXSTRINGLENGTH) :: err_string
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  class(carbon_sandbox_rxn_base_type), pointer :: new_rxn

  found = PETSC_TRUE
  select case(trim(keyword))
    case('CONCENTRATION_UNITS')
      call InputReadWord(input,option,word,PETSC_TRUE)
      call InputErrorMsg(input,option,keyword,err_string)
      call StringToUpper(word)
      select case(word)
        case('MOLALITY')
          this%concentration_units = CARBON_UNITS_MOLALITY
        case('MOLARITY')
          this%concentration_units = CARBON_UNITS_MOLARITY
        case('MOLE_PER_KG_SOIL')
          this%concentration_units = CARBON_UNITS_MOLE_PER_KG_SOIL
        case('MOLE_PER_CUBIC_METER_BULK')
          this%concentration_units = CARBON_UNITS_MOLE_PER_M3_BULK
        case('MOLES')
          this%concentration_units = CARBON_UNITS_MOLES
        case default
          err_string = trim(err_string) // ',CONCENTRATION_UNITS'
          call InputKeywordUnrecognized(input,word,err_string,option)
      end select
    case('REACTION_NETWORK')
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
          case('REACTION')
            new_rxn => CarbonRxnBaseCreate()
            call new_rxn%ReadInput(input,option)
            call CarbonRxnBaseAppend(new_rxn,this%rxn_list)
          case default
            call InputKeywordUnrecognized(input,keyword,err_string,option)
        end select
      enddo
      call InputPopBlock(input,option)
    case default
      found = PETSC_FALSE
  end select

end subroutine CarbonBaseReadSelectCase

! ************************************************************************** !

subroutine CarbonBaseSetup(this,reaction,option)
  !
  ! Configures the reaction and associated data structures
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Option_module
  use Reaction_Aux_module

  class(carbon_sandbox_base_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  class(carbon_sandbox_rxn_base_type), pointer :: cur_rxn

  if (Uninitialized(this%concentration_units)) then
    option%io_buffer = 'CONCENTRATION_UNITS not defined in CarbonBaseSetup.'
    call PrintErrMsg(option)
  endif

  allocate(this%aux)
  this%aux%liquid_saturation = UNINITIALIZED_DOUBLE
  this%aux%liquid_density = UNINITIALIZED_DOUBLE
  this%aux%porosity = UNINITIALIZED_DOUBLE
  this%aux%cell_volume = UNINITIALIZED_DOUBLE
  this%aux%liter_water = UNINITIALIZED_DOUBLE
  this%aux%temperature = UNINITIALIZED_DOUBLE
  allocate(this%aux%conc(reaction%ncomp))
  this%aux%conc = UNINITIALIZED_DOUBLE
  allocate(this%aux%ln_conc(reaction%ncomp))
  this%aux%ln_conc = UNINITIALIZED_DOUBLE
  allocate(this%aux%inhibition_conc(reaction%ncomp))
  this%aux%inhibition_conc = UNINITIALIZED_DOUBLE

  cur_rxn => this%rxn_list
  do
    if (.not.associated(cur_rxn)) exit
    call cur_rxn%Setup(this%aux,reaction,option)
    cur_rxn => cur_rxn%next
  enddo

end subroutine CarbonBaseSetup

! ************************************************************************** !

subroutine CarbonBaseEnforceConcUnits(this,iunits,option)
  !
  ! Checks whether appropriate carbon units have been set
  !
  ! Author: Glenn Hammond
  ! Date: 02/19/24
  !
  use Option_module

  class(carbon_sandbox_base_type) :: this
  PetscInt :: iunits
  type(option_type) :: option

  if (this%concentration_units /= iunits) then
    option%io_buffer = 'Only concentration units of "' // &
      this%GetConcentrationUnits() // &
      '" are supported by the carbon sandbox.'
    call PrintErrMsg(option)
  endif

end subroutine CarbonBaseEnforceConcUnits

! ************************************************************************** !

function CarbonBaseGetConcentrationUnits(this)
  !
  ! Checks whether appropriate carbon units have been set
  !
  ! Author: Glenn Hammond
  ! Date: 02/19/24
  !
  class(carbon_sandbox_base_type) :: this

  character(len=:), allocatable :: CarbonBaseGetConcentrationUnits

  CarbonBaseGetConcentrationUnits = &
    CarbonBaseConcUnitsIntToString(this%concentration_units)

end function CarbonBaseGetConcentrationUnits

! ************************************************************************** !

function CarbonBaseConcUnitsIntToString(iunits)
  !
  ! Maps the integer representing carbon units to the corresponding string
  !
  ! Author: Glenn Hammond
  ! Date: 02/19/24
  !
  PetscInt :: iunits

  character(len=:), allocatable :: CarbonBaseConcUnitsIntToString

  CarbonBaseConcUnitsIntToString = ''
  select case(iunits)
    case(CARBON_UNITS_MOLALITY)
      CarbonBaseConcUnitsIntToString = 'MOLALITY'
    case(CARBON_UNITS_MOLARITY)
      CarbonBaseConcUnitsIntToString = 'MOLARITY'
    case(CARBON_UNITS_MOLE_PER_KG_SOIL)
      CarbonBaseConcUnitsIntToString = 'MOLE_PER_KG_SOIL'
    case(CARBON_UNITS_MOLE_PER_M3_BULK)
      CarbonBaseConcUnitsIntToString = 'MOLE_PER_CUBIC_METER_BULK'
    case(CARBON_UNITS_MOLES)
      CarbonBaseConcUnitsIntToString = 'MOLES'
  end select

end function CarbonBaseConcUnitsIntToString

! ************************************************************************** !

function CarbonBaseConcUnitsStringToInt(string)
  !
  ! Maps the string representing carbon units to the corresponding integer
  !
  ! Author: Glenn Hammond
  ! Date: 02/19/24
  !
  character(len=*) :: string

  PetscInt :: CarbonBaseConcUnitsStringToInt

  CarbonBaseConcUnitsStringToInt = UNINITIALIZED_INTEGER
  select case(string)
    case('MOLALITY')
      CarbonBaseConcUnitsStringToInt = CARBON_UNITS_MOLALITY
    case('MOLARITY')
      CarbonBaseConcUnitsStringToInt = CARBON_UNITS_MOLARITY
    case('MOLE_PER_KG_SOIL')
      CarbonBaseConcUnitsStringToInt = CARBON_UNITS_MOLE_PER_KG_SOIL
    case('MOLE_PER_CUBIC_METER_BULK')
      CarbonBaseConcUnitsStringToInt = CARBON_UNITS_MOLE_PER_M3_BULK
    case('MOLES')
      CarbonBaseConcUnitsStringToInt = CARBON_UNITS_MOLES
  end select

end function CarbonBaseConcUnitsStringToInt

! ************************************************************************** !

subroutine CarbonBaseMapStateVariables(this,rt_auxvar,global_auxvar, &
                                       material_auxvar,reaction,option)
  !
  ! Maps external state variables to aux object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Global_Aux_module
  use Material_Aux_module
  use Option_module
  use Reaction_Aux_module
  use Reaction_Inhibition_Aux_module
  use Reactive_Transport_Aux_module
  use String_module

  class(carbon_sandbox_base_type) :: this
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  PetscInt :: i
  PetscReal :: conc_volume

  this%aux%liquid_saturation = global_auxvar%sat(1)
  this%aux%liquid_density = global_auxvar%den_kg(1)
  this%aux%porosity = material_auxvar%porosity
  this%aux%cell_volume = material_auxvar%volume
  this%aux%liter_water = this%aux%liquid_saturation * &
                         this%aux%porosity * &
                         this%aux%cell_volume * 1.d3
  this%aux%temperature = global_auxvar%temp

  do i = 1, reaction%naqcomp
    this%aux%conc(i) = rt_auxvar%pri_molal(i)*this%aux%liter_water
  enddo
  do i = 1, reaction%immobile%nimmobile
    this%aux%conc(reaction%offset_immobile+i) = &
      rt_auxvar%immobile(i)*this%aux%cell_volume
  enddo
  conc_volume = UNINITIALIZED_DOUBLE
  select case(this%concentration_units)
    case(CARBON_UNITS_MOLALITY)
      conc_volume = this%aux%liter_water*this%aux%liquid_density*1.d-3
    case(CARBON_UNITS_MOLARITY)
      conc_volume = this%aux%liter_water
    case(CARBON_UNITS_MOLE_PER_KG_SOIL)
      if (Uninitialized(material_auxvar%soil_particle_density)) then
        option%io_buffer = 'ROCK_DENSITY must be defined for each material &
          &to use units of MOLE_PER_KG_SOIL.'
        call PrintErrMsg(option)
      endif
      conc_volume = (1.d0-this%aux%porosity) * &
                    material_auxvar%soil_particle_density * &
                    this%aux%cell_volume
    case(CARBON_UNITS_MOLE_PER_M3_BULK)
      conc_volume = this%aux%cell_volume
    case(CARBON_UNITS_MOLES)
      conc_volume = 1.d0
    case default
      option%io_buffer = 'Unrecognized concentration units in &
        &CarbonBaseMapStateVariables(): ' // &
        StringWrite(this%concentration_units)
      call PrintErrMsg(option)
  end select
  this%aux%conc = this%aux%conc/conc_volume
  this%aux%ln_conc = log(this%aux%conc)
  !this%aux%inhibition_conc =

end subroutine CarbonBaseMapStateVariables

! ************************************************************************** !

subroutine CarbonBaseEvaluate(this,Residual,Jacobian,compute_derivative, &
                              rt_auxvar,global_auxvar,material_auxvar, &
                              reaction,option)
  !
  ! Evaluates the rate expression
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
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

  class(carbon_sandbox_rxn_base_type), pointer :: cur_rxn

  call this%MapStateVariables(rt_auxvar,global_auxvar,material_auxvar, &
                              reaction,option)
  cur_rxn => this%rxn_list
  do
    if (.not.associated(cur_rxn)) exit
    call cur_rxn%Evaluate(Residual,Jacobian,option)
    cur_rxn => cur_rxn%next
  enddo

end subroutine CarbonBaseEvaluate

! ************************************************************************** !

subroutine CarbonBaseDestroyAux(aux)
  !
  ! Destroys the carbon sandbox aux object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Utility_module

  type(carbon_sandbox_aux_type), pointer :: aux

  if (.not.associated(aux)) return

  call DeallocateArray(aux%conc)
  call DeallocateArray(aux%ln_conc)
  call DeallocateArray(aux%conc)
  deallocate(aux)
  nullify(aux)

end subroutine CarbonBaseDestroyAux

! ************************************************************************** !

subroutine CarbonBaseStrip(this)
  !
  ! Destroys the carbon sandbox object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  class(carbon_sandbox_base_type) :: this

  call CarbonBaseDestroyAux(this%aux)
  call CarbonRxnBaseDestroyList(this%rxn_list)

end subroutine CarbonBaseStrip

! ************************************************************************** !

subroutine CarbonBaseDestroy(this)
  !
  ! Destroys the carbon sandbox object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  class(carbon_sandbox_base_type), pointer :: this

  call CarbonBaseStrip(this)
  deallocate(this)
  nullify(this)

end subroutine CarbonBaseDestroy

! ************************************************************************** !

function CarbonRxnBaseCreate()
  !
  ! Allocates and initializes the carbon sandbox rxn object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
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

recursive subroutine CarbonRxnBaseAppend(new_rxn,rxn_list)
  !
  ! Appends a new reaction to a reaction list
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  class(carbon_sandbox_rxn_base_type), pointer :: new_rxn
  class(carbon_sandbox_rxn_base_type), pointer :: rxn_list

  if (associated(rxn_list)) then
    call CarbonRxnBaseAppend(new_rxn,rxn_list%next)
    return
  endif
  rxn_list => new_rxn

end subroutine CarbonRxnBaseAppend

! ************************************************************************** !

subroutine CarbonRxnBaseReadInput(this,input,option)
  !
  ! Reads reaction parameters from input deck
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Input_Aux_module
  use Option_module
  use String_module

  class(carbon_sandbox_rxn_base_type) :: this
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

end subroutine CarbonRxnBaseReadInput

! ************************************************************************** !

subroutine CarbonRxnBaseSetup(this,aux,reaction,option)
  !
  ! Configures the reaction and associated data structures
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Reaction_Aux_module
  use Reaction_Equation_module

  class(carbon_sandbox_rxn_base_type) :: this
  type(carbon_sandbox_aux_type), pointer :: aux
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  if (Uninitialized(this%rate_constant)) then
    option%io_buffer = 'RATE_CONSTANT not defined in CarbonRxnBaseSetup.'
    call PrintErrMsg(option)
  endif
  if (len_trim(this%reaction_string) <= 1) then
    option%io_buffer = 'A REACTION_EQUATION is not defined in &
      &CarbonRxnBaseSetup.'
    call PrintErrMsg(option)
  endif

  this%aux => aux
  this%reaction_equation => &
      ReactionEquationCreateFromString(this%reaction_string,option)
  call ReactionEquationMapSpeciesNames(this%reaction_equation, &
                                       reaction%naqcomp, &
                                       reaction%offset_aqueous, &
                                       reaction%primary_species_names, &
                                       reaction%nimcomp, &
                                       reaction%offset_immobile, &
                                       reaction%immobile%names, &
                                       PETSC_FALSE,option)
  call ReactionEquationRemoveSpecies(this%reaction_equation,h2oname,option)

end subroutine CarbonRxnBaseSetup

! ************************************************************************** !

subroutine CarbonRxnBaseEvaluate(this,Residual,Jacobian,option)
  !
  ! Evaluates the rate expression
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  use Option_module
  use Reaction_Inhibition_Aux_module

  class(carbon_sandbox_rxn_base_type) :: this
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

end subroutine CarbonRxnBaseEvaluate

! ************************************************************************** !

subroutine CarbonRxnBaseStrip(this)
  !
  ! Destroys members of carbon sandbox object but not the object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  class(carbon_sandbox_rxn_base_type) :: this

  call ReactionEquationDestroy(this%reaction_equation)
  nullify(this%aux)
  nullify(this%next)

end subroutine CarbonRxnBaseStrip

! ************************************************************************** !

subroutine CarbonRxnBaseDestroyList(rxn_list)
  !
  ! Destroys the carbon sandbox object
  !
  ! Author: Glenn Hammond
  ! Date: 02/16/24
  !
  class(carbon_sandbox_rxn_base_type), pointer :: rxn_list

  class(carbon_sandbox_rxn_base_type), pointer :: cur_rxn, next_rxn

  cur_rxn => rxn_list
  do
    if (.not.associated(cur_rxn)) exit
    next_rxn => cur_rxn%next
    call cur_rxn%Strip()
    deallocate(cur_rxn)
    nullify(cur_rxn)
    cur_rxn => next_rxn
  enddo

end subroutine CarbonRxnBaseDestroyList

end module Carbon_Sandbox_Base_class
