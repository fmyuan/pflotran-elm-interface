module Carbon_Sandbox_MEND_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Carbon_Sandbox_Base_class
  use Option_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, extends(carbon_sandbox_base_type), public :: carbon_sandbox_mend_type
    character(len=MAXWORDLENGTH) :: B_species_name
    character(len=MAXWORDLENGTH) :: D_species_name
    character(len=MAXWORDLENGTH) :: EM_species_name
    character(len=MAXWORDLENGTH) :: EP_species_name
    character(len=MAXWORDLENGTH) :: IC_species_name
    character(len=MAXWORDLENGTH) :: M_species_name
    character(len=MAXWORDLENGTH) :: P_species_name
    character(len=MAXWORDLENGTH) :: Q_species_name
    PetscInt :: B_species_index
    PetscInt :: D_species_index
    PetscInt :: EM_species_index
    PetscInt :: EP_species_index
    PetscInt :: IC_species_index
    PetscInt :: M_species_index
    PetscInt :: P_species_index
    PetscInt :: Q_species_index
    PetscReal :: V_P
    PetscReal :: K_P
    PetscReal :: V_M
    PetscReal :: K_M
    PetscReal :: V_D
    PetscReal :: K_D
    PetscReal :: m_R
    PetscReal :: E_C
    PetscReal :: f_D
    PetscReal :: g_D
    PetscReal :: P_EP
    PetscReal :: P_EM
    PetscReal :: r_EP
    PetscReal :: r_EM
    PetscReal :: Q_max
    PetscReal :: K_ads
    PetscReal :: K_des
    PetscReal :: K_BA
    PetscReal :: I_P
    PetscReal :: I_D
    !PetscReal :: fI_D ! ratio I_D/I_P (currently unsupported)
    PetscReal :: reference_temperature
  contains
    procedure, public :: ReadInput => CarbonMENDReadInput
    procedure, public :: Setup => CarbonMENDSetup
    procedure, public :: Evaluate => CarbonMENDEvaluate
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
  this%B_species_name = 'biomass'
  this%D_species_name = 'doc'
  this%EM_species_name = 'mineral_enzyme'
  this%EP_species_name = 'particulate_enzyme'
  this%IC_species_name = 'inorganic_carbon'
  this%M_species_name = 'mineral_associated_carbon'
  this%P_species_name = 'particulate_carbon'
  this%Q_species_name = 'adsorbed_carbon'
  this%B_species_index = UNINITIALIZED_INTEGER
  this%D_species_index = UNINITIALIZED_INTEGER
  this%EM_species_index = UNINITIALIZED_INTEGER
  this%EP_species_index = UNINITIALIZED_INTEGER
  this%IC_species_index = UNINITIALIZED_INTEGER
  this%M_species_index = UNINITIALIZED_INTEGER
  this%P_species_index = UNINITIALIZED_INTEGER
  this%Q_species_index = UNINITIALIZED_INTEGER
  this%V_P = UNINITIALIZED_DOUBLE
  this%K_P = UNINITIALIZED_DOUBLE
  this%V_M = UNINITIALIZED_DOUBLE
  this%K_M = UNINITIALIZED_DOUBLE
  this%V_D = UNINITIALIZED_DOUBLE
  this%K_D = UNINITIALIZED_DOUBLE
  this%m_R = UNINITIALIZED_DOUBLE
  this%E_C = UNINITIALIZED_DOUBLE
  this%f_D = UNINITIALIZED_DOUBLE
  this%g_D = UNINITIALIZED_DOUBLE
  this%P_EP = UNINITIALIZED_DOUBLE
  this%P_EM = UNINITIALIZED_DOUBLE
  this%r_EP = UNINITIALIZED_DOUBLE
  this%r_EM = UNINITIALIZED_DOUBLE
  this%Q_max = UNINITIALIZED_DOUBLE
  this%K_ads = UNINITIALIZED_DOUBLE
  this%K_des = UNINITIALIZED_DOUBLE
  this%K_BA = UNINITIALIZED_DOUBLE
  this%I_P = UNINITIALIZED_DOUBLE
  this%I_D = UNINITIALIZED_DOUBLE
  !this%fI_D = UNINITIALIZED_DOUBLE
  this%reference_temperature = UNINITIALIZED_DOUBLE

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
  character(len=MAXWORDLENGTH) :: internal_units
  character(len=MAXSTRINGLENGTH) :: err_string
  PetscReal :: tempreal
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
      case('B_SPECIES_NAME')
        call InputReadWord(input,option,this%B_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('D_SPECIES_NAME')
        call InputReadWord(input,option,this%D_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('EM_SPECIES_NAME')
        call InputReadWord(input,option,this%EM_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('EP_SPECIES_NAME')
        call InputReadWord(input,option,this%EP_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('IC_SPECIES_NAME')
        call InputReadWord(input,option,this%IC_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('M_SPECIES_NAME')
        call InputReadWord(input,option,this%M_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('P_SPECIES_NAME')
        call InputReadWord(input,option,this%P_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('Q_SPECIES_NAME')
        call InputReadWord(input,option,this%Q_species_name,PETSC_TRUE)
        call InputErrorMsg(input,option,keyword,err_string)
      case('V_P','V_M','V_D','M_R','R_EP','R_EM','K_ADS','K_DES', &
           'I_P','I_D')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,err_string)
        internal_units = 'kg/kg-s'
        call InputReadAndConvertUnits(input,tempreal,internal_units, &
                                      keyword,option)
        select case(keyword)
          case('V_P')
            this%V_P = tempreal
          case('V_M')
            this%V_M = tempreal
          case('V_D')
            this%V_D = tempreal
          case('K_D')
            this%K_D = tempreal
          case('M_R')
            this%M_R = tempreal
          case('R_EP')
            this%r_EP = tempreal
          case('R_EM')
            this%r_EM = tempreal
          case('K_ADS')
            this%K_ads = tempreal
          case('K_DES')
            this%K_des = tempreal
          case('I_P')
            this%I_P = tempreal
          case('I_D')
            this%I_D = tempreal
        end select
      case('K_P','K_M','K_D','Q_MAX')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,err_string)
        internal_units = 'kg/kg'
        call InputReadAndConvertUnits(input,tempreal,internal_units, &
                                      keyword,option)
        select case(keyword)
          case('K_P')
            this%K_P = tempreal
          case('K_M')
            this%K_M = tempreal
          case('K_D')
            this%K_D = tempreal
          case('Q_MAX')
            this%Q_max = tempreal
        end select
      case('K_BA')
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,err_string)
        internal_units = 'kg/kg'  ! these are inverted units
        call InputReadAndConvertUnits(input,tempreal,internal_units, &
                                      keyword,option)
        this%K_BA = tempreal
      case('E_C','F_D','G_D','P_EP','P_EM') ! unitless
        call InputReadDouble(input,option,tempreal)
        call InputErrorMsg(input,option,keyword,err_string)
        select case(keyword)
          case('E_C')
            this%E_C = tempreal
          case('F_D')
            this%f_D = tempreal
          case('G_D')
            this%g_D = tempreal
          case('P_EP')
            this%p_EP = tempreal
          case('P_EM')
            this%p_EM = tempreal
        end select
      case('REFERENCE_TEMPERATURE')
        call InputReadDouble(input,option,this%reference_temperature)
        call InputErrorMsg(input,option,keyword,err_string)
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
  use Material_Aux_module
  use Option_module
  use Reaction_Aux_module
  use Reaction_Immobile_Aux_module

  class(carbon_sandbox_mend_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  call CarbonBaseSetup(this,reaction,option)

  this%B_species_index = &
    ReactionAuxGetPriSpecIDFromName(this%B_species_name,reaction,option)
  this%D_species_index = &
    ReactionAuxGetPriSpecIDFromName(this%D_species_name,reaction,option)
  this%EM_species_index = &
    ReactionAuxGetPriSpecIDFromName(this%EM_species_name,reaction,option)
  this%EP_species_index = &
    ReactionAuxGetPriSpecIDFromName(this%EP_species_name,reaction,option)
  this%IC_species_index = &
    ReactionAuxGetPriSpecIDFromName(this%IC_species_name,reaction,option)
  this%M_species_index = &
    ReactionImGetSpeciesIDFromName(this%M_species_name,reaction%immobile, &
                                   PETSC_FALSE,option) + &
    reaction%offset_immobile
  if (this%M_species_index <= 0) then
    option%io_buffer = 'Species "' // trim(this%M_species_name) // &
      '" not found among the available immobile species. Mineral associated &
      &biomass must be an immobile species.'
    call PrintErrMsg(option)
  endif
  this%P_species_index = &
    ReactionAuxGetPriSpecIDFromName(this%P_species_name,reaction,option)
  this%Q_species_index = &
    ReactionImGetSpeciesIDFromName(this%Q_species_name,reaction%immobile, &
                                   PETSC_FALSE,option) + &
    reaction%offset_immobile
  if (this%M_species_index <= 0) then
    option%io_buffer = 'Species "' // trim(this%M_species_name) // &
      '" not found among the available immobile species. Adsorbed biomass &
      &must be an immobile species.'
    call PrintErrMsg(option)
  endif

  if (Uninitialized(this%V_P)) then
    option%io_buffer = 'V_P not defined in MEND carbon sandbox.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%K_P)) then
    option%io_buffer = 'K_P not defined in MEND carbon sandbox.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%V_M)) then
    option%io_buffer = 'V_M not defined in MEND carbon sandbox.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%K_M)) then
    option%io_buffer = 'K_M not defined in MEND carbon sandbox.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%V_D)) then
    option%io_buffer = 'V_D not defined in MEND carbon sandbox.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%K_D)) then
    option%io_buffer = 'K_D not defined in MEND carbon sandbox.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%m_R)) then
    option%io_buffer = 'm_R not defined in MEND carbon sandbox.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%E_C)) then
    option%io_buffer = 'E_C not defined in MEND carbon sandbox.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%f_D)) then
    option%io_buffer = 'f_D not defined in MEND carbon sandbox.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%g_D)) then
    option%io_buffer = 'g_D not defined in MEND carbon sandbox.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%P_EP)) then
    option%io_buffer = 'P_EP not defined in MEND carbon sandbox.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%P_EM)) then
    option%io_buffer = 'P_EM not defined in MEND carbon sandbox.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%r_EP)) then
    option%io_buffer = 'r_EP not defined in MEND carbon sandbox.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%r_EM)) then
    option%io_buffer = 'r_EM not defined in MEND carbon sandbox.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%Q_max)) then
    option%io_buffer = 'Q_max not defined in MEND carbon sandbox.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%K_ads)) then
    option%io_buffer = 'K_ads not defined in MEND carbon sandbox.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%K_des)) then
    option%io_buffer = 'K_des not defined in MEND carbon sandbox.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%K_BA)) then
    option%io_buffer = 'K_BA not defined in MEND carbon sandbox.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%K_ads)) then
    option%io_buffer = 'K_ads not defined in MEND carbon sandbox.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%I_D)) then
    ! these source terms are not required
    this%I_D = 0.d0
  endif
  if (Uninitialized(this%I_P)) then
    ! these source terms are not required
    this%I_P = 0.d0
  endif
  if (Uninitialized(this%reference_temperature)) then
    this%reference_temperature = option%flow%reference_temperature
  endif

  call this%EnforceConcentrationUnits(CARBON_UNITS_MOLE_PER_KG_SOIL,option)

end subroutine CarbonMENDSetup

! ************************************************************************** !

subroutine CarbonMENDEvaluate(this,Residual,Jacobian,compute_derivative, &
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
  use Utility_module

  class(carbon_sandbox_mend_type) :: this
  class(reaction_rt_type) :: reaction
  ! the following arrays must be declared after reaction
  PetscReal :: Residual(:)
  PetscReal :: Jacobian(:,:)
  PetscBool :: compute_derivative
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(option_type) :: option

  PetscReal :: B, D, EM, EP, IC, M, P, Q
  PetscReal :: D_monod, M_monod, P_monod
  PetscReal :: one_over_E_C
  PetscReal :: F1, F2, F3, F4, F5, F6, F7, F8, F9EP, F9EM, F10EP, F10EM
  PetscReal :: dB, dD, dEM, dEP, dIC, dM, dP, dQ
  PetscReal :: conc_units_conversion
  PetscReal :: V_P_adj, V_M_adj, V_D_adj, K_P_adj, K_M_adj, K_D_adj
  PetscReal :: m_R_adj, E_C_adj
  PetscReal :: K_ads_adj, K_des_adj
  PetscReal :: t
  PetscReal :: tempreal

  call this%MapStateVariables(rt_auxvar,global_auxvar,material_auxvar, &
                              reaction,option)

  t = this%aux%temperature
  ! mgC/gsoil = molC/kgsoil*gC/molC*mgC/gC*kgsoil/gsoil
  conc_units_conversion = 30.d0*1000.d0*1.d-3
  B = this%aux%conc(this%B_species_index)*conc_units_conversion
  D = this%aux%conc(this%D_species_index)*conc_units_conversion
  EM = this%aux%conc(this%EM_species_index)*conc_units_conversion
  EP = this%aux%conc(this%EP_species_index)*conc_units_conversion
  IC = this%aux%conc(this%IC_species_index)*conc_units_conversion
  M = this%aux%conc(this%M_species_index)*conc_units_conversion
  P = this%aux%conc(this%P_species_index)*conc_units_conversion
  Q = this%aux%conc(this%Q_species_index)*conc_units_conversion

  V_P_adj = this%V_P*Arrhenius(53.d0,t,this%reference_temperature)
  tempreal = Arrhenius(47.d0,t,this%reference_temperature)
  V_M_adj = this%V_M*tempreal
  V_D_adj = this%V_D*tempreal
  tempreal = Arrhenius(30.d0,t,this%reference_temperature)
  K_P_adj = this%K_P*tempreal
  K_M_adj = this%K_M*tempreal
  K_D_adj = this%K_D*tempreal
  m_R_adj = this%m_R*Arrhenius(20.d0,t,this%reference_temperature)
  K_ads_adj = this%K_ads*Arrhenius(5.d0,t,this%reference_temperature)
  K_des_adj = this%K_des*Arrhenius(20.d0,t,this%reference_temperature)
  E_C_adj = &
    max(min(this%E_C - 0.012d0*(t-this%reference_temperature),1.d0),0.d0)

  D_monod = D/(this%K_D+D)
  M_monod = M/(this%K_M+M)
  P_monod = P/(this%K_P+P)
  one_over_E_C = 1.d0 / E_C_adj
  F1 = one_over_E_C*(V_D_adj+m_R_adj)*B*D_monod
  F2 = V_P_adj*EP*P_monod
  F3 = V_M_adj*EM*M_monod
  F4 = (one_over_E_C-1.d0)*V_D_adj*B
  F5 = (one_over_E_C-1.d0)*m_R_adj*B*D_monod
  F6 = K_ads_adj*(1.d0-Q/this%Q_max)*D
  F7 = K_des_adj*Q/this%Q_max
  F8 = (1.d0-this%p_EP-this%p_EM)*m_R_adj*B
  F9EM = this%p_EM*m_R_adj*B
  F9EP = this%p_EP*m_R_adj*B
  F10EM = this%r_EM*EM
  F10EP = this%r_EP*EP

  dB = F1 - (F4 + F5) - F8 - (F9EP + F9EM)
  dD = this%I_D + this%f_D*F2 + this%g_D*F8 + F3 + (F10EP + F10EM) - &
       F1 - (F6 - F7)
  dEM = F9EM - F10EM
  dEP = F9EP - F10EP
  dM = (1.d0-this%f_D)*F2 - F3
  dP = this%I_P + (1.d0-this%g_D)*F8 - F2
  dQ = F6 - F7
  dIC = F4 + F5

  Residual(this%B_species_index) = Residual(this%B_species_index) - dB
  Residual(this%D_species_index) = Residual(this%D_species_index) - dD
  Residual(this%EM_species_index) = Residual(this%EM_species_index) - dEM
  Residual(this%EP_species_index) = Residual(this%EP_species_index) - dEP
  Residual(this%IC_species_index) = Residual(this%IC_species_index) - dIC
  Residual(this%M_species_index) = Residual(this%M_species_index) - dM
  Residual(this%P_species_index) = Residual(this%P_species_index) - dP
  Residual(this%Q_species_index) = Residual(this%Q_species_index) - dQ

end subroutine CarbonMENDEvaluate

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
