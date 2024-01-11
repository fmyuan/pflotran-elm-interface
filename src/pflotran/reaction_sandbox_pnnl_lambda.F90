module Reaction_Sandbox_Lambda_class

#include "petsc/finclude/petscsys.h"
  use petscsys
  use Global_Aux_module
  use PFLOTRAN_Constants_module
  use Reaction_Sandbox_Base_class
  use Reactive_Transport_Aux_module

  implicit none

  private

  PetscInt, parameter :: LAMBDA_THRESHOLD_INHIBITION = 1
  PetscInt, parameter :: LAMBDA_SMOOTHSTEP_INHIBITION = 2

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_lambda_type
    character(len=MAXSTRINGLENGTH) :: reaction_network_filename
    character(len=MAXWORDLENGTH) :: scaling_mineral_name
    PetscReal, pointer :: stoich(:,:)
    PetscInt, pointer :: i_donor(:)
    PetscInt :: inhibition_type

    !Number of species, reactions and carbon sources in problem
    PetscInt :: n_species
    PetscInt :: n_rxn
    PetscInt :: n_donor

    !Lambda Parameters
    PetscReal :: mu_max
    PetscReal :: vh
    PetscReal :: k_deg
    PetscReal :: cc
    PetscReal :: nh4_inhibit
    PetscInt :: i_o2
    PetscInt :: i_biomass
    PetscInt :: i_nh4
    PetscInt :: i_scaling_mineral

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
  ! Author: Katie Muller
  ! Date: 05/25/22

  implicit none

  class(reaction_sandbox_lambda_type), pointer :: LambdaCreate

  allocate(LambdaCreate)
  LambdaCreate%reaction_network_filename = ''
  LambdaCreate%scaling_mineral_name = ''

  LambdaCreate%n_species = UNINITIALIZED_INTEGER
  LambdaCreate%n_rxn = UNINITIALIZED_INTEGER
  LambdaCreate%n_donor = UNINITIALIZED_INTEGER
  LambdaCreate%i_o2 = UNINITIALIZED_INTEGER
  LambdaCreate%i_biomass = UNINITIALIZED_INTEGER
  LambdaCreate%i_nh4 = UNINITIALIZED_INTEGER
  LambdaCreate%i_scaling_mineral = UNINITIALIZED_INTEGER
  LambdaCreate%inhibition_type = LAMBDA_THRESHOLD_INHIBITION

  LambdaCreate%mu_max = UNINITIALIZED_DOUBLE
  LambdaCreate%vh = 1.d0 ! m^3
  LambdaCreate%k_deg = UNINITIALIZED_DOUBLE
  LambdaCreate%cc = UNINITIALIZED_DOUBLE
  LambdaCreate%nh4_inhibit = UNINITIALIZED_DOUBLE

  nullify(LambdaCreate%stoich)
  nullify(LambdaCreate%i_donor)
  nullify(LambdaCreate%next)

end function LambdaCreate

! ************************************************************************** !

subroutine LambdaRead(this,input,option)
  !
  ! Reads input deck for lambda reaction parameters
  ! Author: Katie Muller
  ! Date: 07/05/22

  use Option_module
  use String_module
  use Input_Aux_module

  implicit none

  class(reaction_sandbox_lambda_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: error_string
  error_string = 'CHEMISTRY,RXN_SANDBOX,LAMBDA'

  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword', &
                       trim(error_string))
    call StringToUpper(word)

    select case(trim(word))
      case('REACTION_NETWORK')
        call InputReadFilename(input,option,this%reaction_network_filename)
        call InputErrorMsg(input,option,word, &
                        trim(error_string))

      case('MU_MAX')
        call InputReadDouble(input,option,this%mu_max)
        call InputErrorMsg(input,option,'mu_max',error_string)
        call InputReadAndConvertUnits(input,this%mu_max,'1/sec|mol/mol-sec',&
                        trim(error_string)//',mu_max',option)

      case('VH')
        call InputReadDouble(input,option,this%vh)
        call InputErrorMsg(input,option,'vh',error_string)
        call InputReadAndConvertUnits(input,this%vh,'m^3',&
                         trim(error_string)//',vh',option)

      case('CC')
        call InputReadDouble(input,option,this%cc)
        call InputErrorMsg(input,option,'cc',error_string)
        call InputReadAndConvertUnits(input,this%cc,'M',&
                          trim(error_string)//',cc',option)

      case('K_DEG')
        call InputReadDouble(input,option,this%k_deg)
        call InputErrorMsg(input,option,'k_deg',error_string)
        call InputReadAndConvertUnits(input,this%k_deg,'1/sec',&
                          trim(error_string)//',k_deg',option)

      case('NH4_INHIBIT')
        call InputReadDouble(input,option,this%nh4_inhibit)
        call InputErrorMsg(input,option,'nh4_inhibit',error_string)
        call InputReadAndConvertUnits(input,this%nh4_inhibit,'M',&
                          trim(error_string)//',nh4_inhibit',option)

      case('INHIBITION_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,word,trim(error_string))
        call STringToUpper(word)
        select case(word)
          case('THRESHOLD')
            this%inhibition_type = LAMBDA_THRESHOLD_INHIBITION
          case('SMOOTHSTEP')
            this%inhibition_type = LAMBDA_SMOOTHSTEP_INHIBITION
          case default
            error_string = trim(error_string) // ',INHIBITION_TYPE'
            call InputKeywordUnrecognized(input,word,error_string ,option)
        end select

      case('SCALING_MINERAL')
        call InputReadWord(input,option,this%scaling_mineral_name,PETSC_TRUE)
        call InputErrorMsg(input,option,word,trim(error_string))

      case default
        call InputKeywordUnrecognized(input,word,error_string ,option)
    end select
  enddo
  call InputPopBlock(input,option)

end subroutine LambdaRead

! ************************************************************************** !

subroutine LambdaSetup(this,reaction,option)
  !
  ! Sets up the lambda reaction with parameters either read from the
  ! input deck or wired.
  !
  ! Author: Katie Muller
  ! Date: 07/05/22

  use Option_module
  use Utility_module
  use Reaction_Aux_module
  use Reaction_Mineral_Aux_module

  implicit none

  class(reaction_sandbox_lambda_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: word
  PetscInt :: i, irxn
  PetscReal, pointer :: stoich(:,:)
  PetscInt, pointer :: species_ids(:,:)

  call ReactionNetworkToStoich(reaction,this%reaction_network_filename, &
                               species_ids,stoich,option)

  ! Determines the number of rxns and species in the problem
  ! and allocates arrays

  this%n_rxn = size(species_ids,2)
  this%n_species = size(reaction%primary_species_names)

  allocate(this%stoich(this%n_species,this%n_rxn))
  allocate(this%i_donor(this%n_rxn))

  this%stoich = 0.d0
  this%i_donor = UNINITIALIZED_INTEGER

  do irxn = 1, this%n_rxn
    do i = 1, species_ids(0,irxn)
      this%stoich(species_ids(i,irxn),irxn) = stoich(i,irxn)
    enddo
  enddo

  call DeallocateArray(stoich)

  word = 'O2(aq)'
  this%i_o2 = &
    ReactionGetPriSpeciesIDFromName(word,reaction,option)

  word = 'NH4+'
  this%i_nh4 = &
    ReactionGetPriSpeciesIDFromName(word,reaction,option)

  word = 'BIOMASS'
  this%i_biomass = &
    ReactionGetPriSpeciesIDFromName(word,reaction,option)

  if (len_trim(this%scaling_mineral_name) > 0) then
    this%i_scaling_mineral = &
      ReactionMnrlGetKinMnrlIDFromName(this%scaling_mineral_name, &
                                       reaction%mineral,option)
  endif

  ! Input file must have carbon species indicated by the word "DONOR"
  ! to be used as donor

  do irxn = 1, this%n_rxn
    do i = 1, species_ids(0,irxn)
      if (index(reaction%primary_species_names(species_ids(i,irxn)),'DONOR') &
        > 0) then
        this%i_donor(irxn) = species_ids(i,irxn)
        exit
      endif
    enddo
    if (UnInitialized(this%i_donor(irxn))) then
      option%io_buffer = 'No DONOR is specified for the reaction. &
        &Please ensure DONOR is used to specify the carbon donor(s) &
        &for each reaction.'
      call PrintErrMsg(option)
    endif
  enddo

  call DeallocateArray(species_ids)
end subroutine LambdaSetup

! ************************************************************************** !

subroutine LambdaEvaluate(this,Residual,Jacobian,compute_derivative, &
                          rt_auxvar,global_auxvar,material_auxvar, &
                          reaction,option)
  !
  ! Evaluates the reaction storing the Residual and/or Jacobian
  !
  ! Author: Katie Muller
  ! Date: 07/05/22

  use Material_Aux_module
  use Option_module
  use Reaction_Aux_module
  use Reaction_Inhibition_Aux_module

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
  PetscReal :: molality_to_molarity
  PetscReal :: vh_L, sumkin, Biomass_mod
  PetscReal :: nh4_inhibition, tempreal
  PetscReal :: C_reactant_inhibit

  PetscReal :: Reactant_inhibition(this%n_species)
  PetscReal :: C_aq(this%n_species)
  PetscReal :: rkin(this%n_rxn)
  PetscReal :: R(this%n_rxn)
  PetscReal :: Rate(this%n_species)
  PetscReal :: u(this%n_rxn)
  PetscReal :: mu_max

  PetscInt :: icomp, irxn, i_carbon

  L_water = material_auxvar%porosity*global_auxvar%sat(iphase)* &
            material_auxvar%volume*1.d3 ! m^3 -> L

  molality_to_molarity = global_auxvar%den_kg(iphase)*1.d-3
  ! concentrations are molarities [M]

  mu_max = this%mu_max
  if (Initialized(this%i_scaling_mineral)) then
    mu_max = mu_max * rt_auxvar%mnrl_volfrac(this%i_scaling_mineral)
  endif

  do icomp = 1, this%n_species
    C_aq(icomp) = rt_auxvar%pri_molal(icomp)* &
                  rt_auxvar%pri_act_coef(icomp)*molality_to_molarity
  enddo

  vh_L = 1000*this%vh  !converting vh from m^3 as input to L

  do irxn = 1, this%n_rxn
    i_carbon = this%i_donor(irxn)

    rkin(irxn) = mu_max * &
                 exp(-(ABS(this%stoich(i_carbon,irxn)) / &
                   (vh_L * C_aq(i_carbon)))) * &
                 exp(-(ABS(this%stoich(this%i_o2,irxn)) / &
                   (vh_L * C_aq(this%i_o2))))  ![1/sec]
  enddo
  ! Cybernetic Formulation
  ! Relative contribution, ui (unitless)
  sumkin = 0.d0
  do irxn = 1, this%n_rxn
    sumkin = sumkin + rkin(irxn)
  enddo

  u = rkin / sumkin
  R = u * rkin ![1/sec]


  C_reactant_inhibit = 1.d-18
  select case(this%inhibition_type)
    case(LAMBDA_SMOOTHSTEP_INHIBITION)
!      call ReactionThreshInhibitSmoothstep(C_aq(this%i_nh4), &
      call ReactionInhibitionSmoothstep(C_aq(this%i_nh4), &
                                             dabs(this%nh4_inhibit),1.d-3, &
                                             nh4_inhibition,tempreal)
      do icomp = 1, this%n_species
        call ReactionInhibitionSmoothstep(C_aq(icomp),C_reactant_inhibit, &
                                          1.d-3, &
                                          Reactant_inhibition(icomp), &
                                          tempreal)
      enddo
    case(LAMBDA_THRESHOLD_INHIBITION)
      call ReactionInhibitionThreshold(C_aq(this%i_nh4), &
                                       dabs(this%nh4_inhibit), &
                                       nh4_inhibition,tempreal)
      do icomp = 1, this%n_species
        call ReactionInhibitionThreshold(C_aq(icomp),C_reactant_inhibit, &
                                         Reactant_inhibition(icomp),tempreal)
      enddo
  end select

  ! Reactions are modulated by biomass concentration
  ! Biomass is moduluated by a carrying capacity (CC)
  Biomass_mod = C_aq(this%i_biomass) * (1 - C_aq(this%i_biomass) / this%cc)
  Biomass_mod = max(Biomass_mod, 0.d0)

  Rate = 0.d0

  do irxn = 1, this%n_rxn
    do icomp = 1, this%n_species
      if (this%stoich(icomp,irxn) < 0.d0) then
        if (icomp == this%i_nh4) then
          R(irxn) = R(irxn) * nh4_inhibition
        else
          R(irxn) = R(irxn) * Reactant_inhibition(icomp)
        endif
      endif
    enddo
    if (C_aq(this%i_biomass) > this%cc) then
      R(irxn) = 0
    endif
    Rate(:) = Rate(:) + this%stoich(:,irxn) * R(irxn)
  enddo
  Rate(:) = Rate(:) * Biomass_mod * L_water
  Rate(this%i_biomass) = Rate(this%i_biomass) - &
                         this%k_deg * C_aq(this%i_biomass) * L_water

  ! Residuals
  Residual(:) = Residual(:) - Rate(:)

  if (compute_derivative) then
    option%io_buffer = 'REACTION_SANDBOX LAMBDA must be run with &
      &NUMERICAL_JACOBIAN listed in the NUMERICAL_METHODS TRANSPORT &
      &NEWTON_SOLVER block as analytical derivatives are not calculated &
      &in the sandbox evaluate routine.'
    call PrintErrMsg(option)
  endif

end subroutine LambdaEvaluate

! ************************************************************************** !

subroutine LambdaDestroy(this)
  !
  ! Destroys allocatable or pointer objects created in this module
  !
  ! Author: Katie Muller
  ! Date: 07/05/22

  use Utility_module

  implicit none
  class(reaction_sandbox_lambda_type) :: this

  call DeallocateArray(this%stoich)

end subroutine LambdaDestroy

end module Reaction_Sandbox_Lambda_class
